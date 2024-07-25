#include <map>
#include <queue>
#include <mutex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "Common.h"
#include "RankBitvector.h"
#include "UnionFind.h"
#include "EdlibWrapper.h"
#include "CompressedStringIndex.h"
#include "TwobitString.h"
#include "ChunkUnitigGraph.h"
#include "ChunkGraphWriter.h"
#include "KmerIterator.h"
#include "ConsensusMaker.h"
#include "SequenceHelper.h"
#include "GraphCleaner.h"
#include "ChunkExtractor.h"
#include "SequenceIdentitySplitter.h"
#include "ChunkHelper.h"
#include "TransitiveClosure.h"
#include "ChunkPhasing.h"
#include "ChunkResolution.h"

double mismatchFraction;

void countReadRepetitiveUnitigs(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "counting read-repetitive unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	std::vector<uint8_t> chunkRepeatCount;
	chunkRepeatCount.resize(chunkBelongsToUnitig.size(), 0);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		phmap::flat_hash_set<size_t> foundHere;
		phmap::flat_hash_set<size_t> repetitiveHere;
		for (auto t : chunksPerRead[i])
		{
			size_t chunk = std::get<2>(t) & maskUint64_t;
			if (foundHere.count(chunk) == 1) repetitiveHere.insert(chunk);
			foundHere.insert(chunk);
		}
		for (size_t chunk : repetitiveHere)
		{
			if (chunk >= chunkRepeatCount.size()) continue;
			if (chunkRepeatCount[chunk] >= 2) continue;
			chunkRepeatCount[chunk] += 1;
		}
	}
	std::vector<bool> unitigRepetitive;
	unitigRepetitive.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < chunkRepeatCount.size(); i++)
	{
		if (chunkRepeatCount[i] < 2) continue;
		if (chunkBelongsToUnitig[i] == std::numeric_limits<size_t>::max()) continue;
		unitigRepetitive[chunkBelongsToUnitig[i]] = true;
	}
	size_t countRepetitive = 0;
	std::ofstream file { "repetitive.txt" };
	for (size_t i = 0; i < unitigRepetitive.size(); i++)
	{
		if (!unitigRepetitive[i]) continue;
		countRepetitive += 1;
		file << i << std::endl;
	}
	std::cerr << countRepetitive << " read-repetitive unitigs" << std::endl;
}

std::vector<std::pair<size_t, size_t>> getUniqueRanges(std::vector<std::pair<size_t, size_t>> raw)
{
	assert(raw.size() >= 1);
	std::sort(raw.begin(), raw.end());
	for (size_t j = raw.size()-1; j < raw.size() && j > 0; j--)
	{
		if (raw[j] != raw[j-1]) continue;
		std::swap(raw[j], raw.back());
		raw.pop_back();
	}
	assert(raw.size() >= 1);
	std::sort(raw.begin(), raw.end());
	return raw;
}

void expandChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t expandedSize)
{
	std::vector<size_t> chunkLengths;
	{
		std::vector<std::vector<size_t>> rawChunkLengths;
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
				size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
				size_t length = std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) + 1;
				while (chunk >= rawChunkLengths.size()) rawChunkLengths.emplace_back();
				rawChunkLengths[chunk].emplace_back(length);
			}
		}
		chunkLengths.resize(rawChunkLengths.size());
		for (size_t i = 0; i < rawChunkLengths.size(); i++)
		{
			assert(rawChunkLengths[i].size() >= 1);
			std::sort(rawChunkLengths[i].begin(), rawChunkLengths[i].end());
			chunkLengths[i] = rawChunkLengths[i][rawChunkLengths[i].size()/2];
		}
	}
	std::map<std::vector<uint64_t>, size_t> chunkmerToNewChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<std::pair<size_t, size_t>> newChunkmerRanges;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) newChunkmerRanges.emplace_back(j, j);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t end = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			while (end < chunksPerRead[i].size() && length < expandedSize)
			{
				end += 1;
				if (end == chunksPerRead[i].size()) break;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end]))) break;
				length += chunkLengths[std::get<2>(chunksPerRead[i][end]) & maskUint64_t];
				length -= kmerSize;
			}
			if (length < expandedSize) continue;
			if (end == chunksPerRead[i].size()) continue;
			newChunkmerRanges.emplace_back(j, end);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t start = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			while (start > 0 && length < expandedSize)
			{
				start -= 1;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][start]))) break;
				length += chunkLengths[std::get<2>(chunksPerRead[i][start]) & maskUint64_t];
				length -= kmerSize;
			}
			if (length < expandedSize) continue;
			newChunkmerRanges.emplace_back(start, j);
		}
		if (newChunkmerRanges.size() == 0)
		{
			chunksPerRead[i].clear();
			continue;
		}
		newChunkmerRanges = getUniqueRanges(newChunkmerRanges);
		assert(newChunkmerRanges.size() >= 1);
		for (size_t j = 1; j < newChunkmerRanges.size(); j++)
		{
			assert(newChunkmerRanges[j].first >= newChunkmerRanges[j-1].first);
			assert(newChunkmerRanges[j].second >= newChunkmerRanges[j-1].second);
		}
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunks;
		for (size_t j = 0; j < newChunkmerRanges.size(); j++)
		{
			size_t newchunk = std::numeric_limits<size_t>::max();
			std::vector<uint64_t> chunkmer;
			bool valid = true;
			for (size_t k = newChunkmerRanges[j].first; k <= newChunkmerRanges[j].second; k++)
			{
				chunkmer.emplace_back(std::get<2>(chunksPerRead[i][k]));
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][k]))) valid = false;
			}
			if (valid)
			{
				if (chunkmerToNewChunk.count(chunkmer) == 1)
				{
					newchunk = chunkmerToNewChunk.at(chunkmer) + firstBitUint64_t;
				}
				else
				{
					newchunk = chunkmerToNewChunk.size();
					chunkmerToNewChunk[chunkmer] = newchunk;
					newchunk += firstBitUint64_t;
				}
			}
			else
			{
				newchunk = std::numeric_limits<size_t>::max();
			}
			newChunks.emplace_back(std::get<0>(chunksPerRead[i][newChunkmerRanges[j].first]), std::get<1>(chunksPerRead[i][newChunkmerRanges[j].second]), newchunk);
		}
		std::swap(chunksPerRead[i], newChunks);
	}
}

void resegmentChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::vector<bool> canMergeRight;
	std::vector<bool> canMergeLeft;
	std::vector<size_t> uniqueRight;
	std::vector<size_t> uniqueLeft;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			while (canMergeLeft.size() <= chunk) canMergeLeft.emplace_back(true);
		}
	}
	canMergeRight.resize(canMergeLeft.size(), true);
	uniqueRight.resize(canMergeLeft.size(), std::numeric_limits<size_t>::max());
	uniqueLeft.resize(canMergeLeft.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i][0])))
		{
			if (std::get<2>(chunksPerRead[i][0]) & firstBitUint64_t)
			{
				canMergeLeft[std::get<2>(chunksPerRead[i][0]) & maskUint64_t] = false;
			}
			else
			{
				canMergeRight[std::get<2>(chunksPerRead[i][0]) & maskUint64_t] = false;
			}
		}
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i].back())))
		{
			if (std::get<2>(chunksPerRead[i].back()) & firstBitUint64_t)
			{
				canMergeRight[std::get<2>(chunksPerRead[i].back()) & maskUint64_t] = false;
			}
			else
			{
				canMergeLeft[std::get<2>(chunksPerRead[i].back()) & maskUint64_t] = false;
			}
		}
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])))
			{
				if (!NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
				{
					if (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t)
					{
						canMergeLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
					}
					else
					{
						canMergeRight[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
					}
				}
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				assert(!NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])));
				if (std::get<2>(chunksPerRead[i][j-1]) & firstBitUint64_t)
				{
					canMergeRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;
				}
				else
				{
					canMergeLeft[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;
				}
				continue;
			}
			assert(!NonexistantChunk(std::get<2>(chunksPerRead[i][j])));
			assert(!NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])));
			if (std::get<2>(chunksPerRead[i][j-1]))
			{
				if (uniqueRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] == std::numeric_limits<size_t>::max()) uniqueRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = std::get<2>(chunksPerRead[i][j]);
				if (uniqueRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] != std::get<2>(chunksPerRead[i][j])) canMergeRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;;
			}
			else
			{
				if (uniqueLeft[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] == std::numeric_limits<size_t>::max()) uniqueLeft[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t;
				if (uniqueLeft[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] != (std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t)) canMergeLeft[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;;
			}
			if (std::get<2>(chunksPerRead[i][j]))
			{
				if (uniqueLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] == std::numeric_limits<size_t>::max()) uniqueLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = std::get<2>(chunksPerRead[i][j-1]);
				if (uniqueLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] != std::get<2>(chunksPerRead[i][j-1])) canMergeLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;;
			}
			else
			{
				if (uniqueRight[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] == std::numeric_limits<size_t>::max()) uniqueRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = std::get<2>(chunksPerRead[i][j-1]) ^ firstBitUint64_t;
				if (uniqueRight[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] != (std::get<2>(chunksPerRead[i][j-1]) ^ firstBitUint64_t)) canMergeRight[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;;
			}
		}
	}
	std::map<std::vector<uint64_t>, size_t> chunkmerToNewChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunksHere;
		size_t currentStart = 0;
		size_t currentEnd = 0;
		std::vector<uint64_t> currentMer;
		currentMer.emplace_back(std::get<2>(chunksPerRead[i][0]));
		for (size_t j = 1; j <= chunksPerRead[i].size(); j++)
		{
			assert(currentMer.size() >= 1);
			assert(currentEnd == j-1);
			bool canContinueMerge = true;
			if (j == chunksPerRead[i].size() || NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || NonexistantChunk(currentMer.back()))
			{
				canContinueMerge = false;
			}
			else
			{
				if (currentMer.back() & firstBitUint64_t)
				{
					if (!canMergeRight[currentMer.back() & maskUint64_t])
					{
						canContinueMerge = false;
					}
				}
				else
				{
					if (!canMergeLeft[currentMer.back() & maskUint64_t])
					{
						canContinueMerge = false;
					}
				}
				if (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t)
				{
					if (!canMergeLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
					{
						canContinueMerge = false;
					}
				}
				else
				{
					if (!canMergeRight[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
					{
						canContinueMerge = false;
					}
				}
			}
			if (canContinueMerge)
			{
				currentEnd = j;
				currentMer.emplace_back(std::get<2>(chunksPerRead[i][j]));
				continue;
			}
			size_t newMer = std::numeric_limits<size_t>::max();
			bool validmer = true;
			for (size_t k = currentStart; k <= currentEnd; k++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][k]))) validmer = false;
			}
			if (validmer)
			{
				if (chunkmerToNewChunk.count(currentMer) == 1)
				{
					newMer = chunkmerToNewChunk.at(currentMer);
				}
				else
				{
					newMer = chunkmerToNewChunk.size();
					chunkmerToNewChunk[currentMer] = newMer;
				}
			}
			else
			{
				newMer = std::numeric_limits<size_t>::max() ^ firstBitUint64_t;
			}
			newChunksHere.emplace_back(std::get<0>(chunksPerRead[i][currentStart]), std::get<1>(chunksPerRead[i][currentEnd]), newMer + firstBitUint64_t);
			currentMer.clear();
			currentStart = j;
			currentEnd = j;
			if (j < chunksPerRead[i].size()) currentMer.emplace_back(std::get<2>(chunksPerRead[i][j]));
		}
		assert(currentEnd == chunksPerRead[i].size());
		assert(currentMer.size() == 0);
		std::swap(chunksPerRead[i], newChunksHere);
	}
}

std::pair<std::pair<std::vector<uint64_t>, std::vector<uint64_t>>, bool> getChunkContext(const std::vector<std::tuple<size_t, size_t, uint64_t>>& chunks, const std::vector<size_t>& chunkLengths, const size_t j, const size_t expandedSize)
{
	std::pair<std::vector<uint64_t>, std::vector<uint64_t>> context;
	bool valid = true;
	assert(j < chunks.size());
	if (NonexistantChunk(std::get<2>(chunks[j]))) return std::make_pair(context, false);
	size_t leftoverPerSide = 0;
	if (chunkLengths[std::get<2>(chunks[j]) & maskUint64_t] < expandedSize)
	{
		leftoverPerSide = (expandedSize - chunkLengths[std::get<2>(chunks[j]) & maskUint64_t]) / 2;
	}
	if (leftoverPerSide >= 1)
	{
		size_t leftoverLeftside = leftoverPerSide;
		size_t leftAdd = 0;
		while (leftoverLeftside >= 1)
		{
			if (j > leftAdd)
			{
				leftAdd += 1;
				if (NonexistantChunk(std::get<2>(chunks[j-leftAdd])))
				{
					valid = false;
					break;
				}
				context.first.emplace_back(std::get<2>(chunks[j-leftAdd]));
				if (leftoverLeftside < chunkLengths[std::get<2>(chunks[j-leftAdd]) & maskUint64_t])
				{
					leftoverLeftside = 0;
					break;
				}
				leftoverLeftside -= chunkLengths[std::get<2>(chunks[j-leftAdd]) & maskUint64_t];
			}
			else
			{
				valid = false;
				break;
			}
		}
		if (leftoverLeftside > 0) valid = false;
		size_t leftoverRightside = leftoverPerSide;
		size_t rightAdd = 0;
		while (leftoverRightside >= 1)
		{
			if (j+rightAdd+1 < chunks.size())
			{
				rightAdd += 1;
				if (NonexistantChunk(std::get<2>(chunks[j+rightAdd])))
				{
					valid = false;
					break;
				}
				context.second.emplace_back(std::get<2>(chunks[j+rightAdd]));
				if (leftoverRightside < chunkLengths[std::get<2>(chunks[j+rightAdd]) & maskUint64_t])
				{
					leftoverRightside = 0;
					break;
				}
				leftoverRightside -= chunkLengths[std::get<2>(chunks[j+rightAdd]) & maskUint64_t];
			}
			else
			{
				valid = false;
				break;
			}
		}
		if (leftoverRightside > 0) valid = false;
	}
	return std::make_pair(context, valid);
}

bool contextMatches(const std::pair<std::vector<uint64_t>, std::vector<uint64_t>>& subcontext, const std::pair<std::vector<uint64_t>, std::vector<uint64_t>>& supercontext)
{
	if (subcontext.first.size() > supercontext.first.size()) return false;
	if (subcontext.second.size() > supercontext.second.size()) return false;
	for (size_t i = 0; i < subcontext.first.size(); i++)
	{
		if (subcontext.first[i] != supercontext.first[i]) return false;
	}
	for (size_t i = 0; i < subcontext.second.size(); i++)
	{
		if (subcontext.second[i] != supercontext.second[i]) return false;
	}
	return true;
}

void contextResolve(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t expandedSize)
{
	std::vector<size_t> chunkLengths;
	{
		std::vector<std::vector<size_t>> rawChunkLengths;
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
				size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
				size_t length = std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) + 1;
				while (chunk >= rawChunkLengths.size()) rawChunkLengths.emplace_back();
				rawChunkLengths[chunk].emplace_back(length);
			}
		}
		chunkLengths.resize(rawChunkLengths.size());
		for (size_t i = 0; i < rawChunkLengths.size(); i++)
		{
			assert(rawChunkLengths[i].size() >= 1);
			std::sort(rawChunkLengths[i].begin(), rawChunkLengths[i].end());
			chunkLengths[i] = rawChunkLengths[i][rawChunkLengths[i].size()/2];
		}
	}
	std::vector<std::map<std::pair<std::vector<uint64_t>, std::vector<uint64_t>>, size_t>> contextToNewChunk;
	contextToNewChunk.resize(chunkLengths.size());
	size_t nextNum = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				continue;
			}
			std::pair<std::vector<uint64_t>, std::vector<uint64_t>> context;
			bool valid = true;
			std::tie(context, valid) = getChunkContext(chunksPerRead[i], chunkLengths, j, expandedSize);
			if (!valid) continue;
			assert(std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t);
			if (contextToNewChunk[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].count(context) == 0)
			{
				contextToNewChunk[std::get<2>(chunksPerRead[i][j]) & maskUint64_t][context] = nextNum;
				nextNum += 1;
			}
		}
	}
	phmap::flat_hash_map<uint64_t, size_t> unresolvedChunkMapping;
	for (size_t i = 0; i < contextToNewChunk.size(); i++)
	{
		if (contextToNewChunk[i].size() != 0) continue;
		unresolvedChunkMapping[i] = nextNum;
		nextNum += 1;
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunksHere;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				newChunksHere.emplace_back(std::get<0>(chunksPerRead[i][j]), std::get<1>(chunksPerRead[i][j]), std::numeric_limits<size_t>::max());
				continue;
			}
			if (unresolvedChunkMapping.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 1)
			{
				newChunksHere.emplace_back(std::get<0>(chunksPerRead[i][j]), std::get<1>(chunksPerRead[i][j]), unresolvedChunkMapping.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) + firstBitUint64_t);
				continue;
			}
			std::pair<std::vector<uint64_t>, std::vector<uint64_t>> context;
			bool valid = true;
			std::tie(context, valid) = getChunkContext(chunksPerRead[i], chunkLengths, j, expandedSize);
			if (valid)
			{
				assert(contextToNewChunk.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t).count(context) == 1);
				newChunksHere.emplace_back(std::get<0>(chunksPerRead[i][j]), std::get<1>(chunksPerRead[i][j]), contextToNewChunk.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t).at(context) + firstBitUint64_t);
			}
			else
			{
				size_t uniqueMatchingContext = std::numeric_limits<size_t>::max();
				for (auto pair : contextToNewChunk.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t))
				{
					if (contextMatches(context, pair.first))
					{
						if (uniqueMatchingContext == std::numeric_limits<size_t>::max())
						{
							uniqueMatchingContext = pair.second;
						}
						else
						{
							uniqueMatchingContext = std::numeric_limits<size_t>::max();
							break;
						}
					}
				}
				newChunksHere.emplace_back(std::get<0>(chunksPerRead[i][j]), std::get<1>(chunksPerRead[i][j]), uniqueMatchingContext | firstBitUint64_t);
			}
		}
		std::swap(chunksPerRead[i], newChunksHere);
	}
}

void mergeFakeBubbles(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
	std::vector<bool> allowedNode;
	SparseEdgeContainer allowedEdges;
	std::tie(allowedNode, allowedEdges) = getAllowedNodesAndEdges(coverages, edgeCoverage, chunksPerRead);
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> merges;
	phmap::flat_hash_set<size_t> removeNodes;
	for (size_t i = 0; i < allowedNode.size(); i++)
	{
		if (!allowedNode[i]) continue;
		if (allowedEdges.getEdges(std::make_pair(i, true)).size() != 1) continue;
		if (allowedEdges.getEdges(std::make_pair(i, false)).size() != 1) continue;
		std::pair<size_t, bool> bubbleStartPair = allowedEdges.getEdges(std::make_pair(i, false))[0];
		assert(!bubbleStartPair.second);
		std::pair<size_t, bool> bubbleEndPair = allowedEdges.getEdges(std::make_pair(i, true))[0];
		assert(bubbleEndPair.second);
		size_t bubbleStart = bubbleStartPair.first;
		size_t bubbleEnd = bubbleEndPair.first;
		assert(allowedNode[bubbleStart]);
		assert(allowedNode[bubbleEnd]);
		if (allowedEdges.getEdges(std::make_pair(bubbleStart, true)).size() != 2) continue;
		if (allowedEdges.getEdges(std::make_pair(bubbleEnd, false)).size() != 2) continue;
		size_t otherAlleleFirst = std::numeric_limits<size_t>::max();;
		size_t otherAlleleSecond = std::numeric_limits<size_t>::max();
		for (auto edge : allowedEdges.getEdges(std::make_pair(bubbleStart, true)))
		{
			assert(edge.second);
			if (edge.first == i) continue;
			otherAlleleFirst = edge.first;
		}
		for (auto edge : allowedEdges.getEdges(std::make_pair(bubbleEnd, false)))
		{
			assert(!edge.second);
			if (edge.first == i) continue;
			otherAlleleSecond = edge.first;
		}
		if (allowedEdges.getEdges(std::make_pair(otherAlleleFirst, false)).size() != 1) continue;
		if (allowedEdges.getEdges(std::make_pair(otherAlleleFirst, true)).size() != 1) continue;
		if (allowedEdges.getEdges(std::make_pair(otherAlleleSecond, false)).size() != 1) continue;
		if (allowedEdges.getEdges(std::make_pair(otherAlleleSecond, true)).size() != 1) continue;
		if (allowedEdges.getEdges(std::make_pair(otherAlleleFirst, true))[0] != std::make_pair(otherAlleleSecond, true)) continue;
		size_t alleleLength = lengths[i][lengths[i].size()/2];
		size_t otherAlleleLength = lengths[otherAlleleFirst][lengths[otherAlleleFirst].size()/2] + lengths[otherAlleleSecond][lengths[otherAlleleSecond].size()/2] - kmerSize;
		if (alleleLength + 50 < otherAlleleLength) continue;
		if (alleleLength > otherAlleleLength + 50) continue;
//		std::cerr << "merge " << otherAlleleFirst << " " << otherAlleleSecond << " into " << i << std::endl;
		merges[std::make_pair(otherAlleleFirst, otherAlleleSecond)] = i;
		removeNodes.insert(otherAlleleFirst);
		removeNodes.insert(otherAlleleSecond);
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		for (size_t j = chunksPerRead[i].size()-1; j > 0; j--)
		{
			if (merges.count(std::make_pair(std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t, std::get<2>(chunksPerRead[i][j]) & maskUint64_t)) == 0) continue;
			chunksPerRead[i].emplace_back(std::get<0>(chunksPerRead[i][j-1]), std::get<1>(chunksPerRead[i][j]), merges.at(std::make_pair(std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t, std::get<2>(chunksPerRead[i][j]) & maskUint64_t)) + firstBitUint64_t);
		}
		for (size_t j = chunksPerRead[i].size()-1; j < chunksPerRead[i].size(); j--)
		{
			if (removeNodes.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			std::swap(chunksPerRead[i][j], chunksPerRead[i].back());
			chunksPerRead[i].pop_back();
		}
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
	}
}

void makeGraph(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize, const size_t windowSize, const size_t middleSkip, const size_t startStage)
{
	std::cerr << "start at stage " << startStage << std::endl;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead;
	if (startStage > 0)
	{
		chunksPerRead = readChunksFromFakePathsFile("fakepaths" + std::to_string(startStage) + ".txt");
		while (chunksPerRead.size() < sequenceIndex.size()*2) chunksPerRead.emplace_back();
		std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
	switch(startStage)
	{
		case 0:
			std::cerr << "getting chunks from reads" << std::endl;
			chunksPerRead = getMinimizerBoundedChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize);
			//addMissingPiecesBetweenChunks(chunksPerRead, k);
			{
				size_t numChunks = 0;
				for (size_t i = 0; i < chunksPerRead.size(); i++)
				{
					numChunks += chunksPerRead[i].size();
				}
				std::cerr << numChunks << " chunks" << std::endl;
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
			splitPerFirstLastKmers(sequenceIndex, chunksPerRead, kmerSize);
			splitPerLength(chunksPerRead, numThreads);
			writeStage(0, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			mergeFakeBubbles(chunksPerRead, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(1, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 1:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerBaseCounts(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(2, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 2:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerMinHashes(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(3, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 3:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(4, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 4:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentityRoughly(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentity(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(5, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 5:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			contextResolve(chunksPerRead, kmerSize, 1000);
//			expandChunks(chunksPerRead, kmerSize, 5000);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			writeStage(51, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 51:
//			resegmentChunks(chunksPerRead, rawReadLengths, approxOneHapCoverage, kmerSize);
//			writeStage(52, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 52:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(6, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 6:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(7, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 7:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerNearestNeighborPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(8, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 8:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerCorrectedKmerClustering(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(81, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 81:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerCorrectedKmerPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(82, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 82:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSNPTransitiveClosureClustering(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(9, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 9:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			contextResolve(chunksPerRead, kmerSize, 2000);
			expandChunks(chunksPerRead, kmerSize, 10000);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(90, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 90:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resegmentChunks(chunksPerRead, rawReadLengths, approxOneHapCoverage, kmerSize);
			writeStage(91, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 91:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSNPTransitiveClosureClustering(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(92, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 92:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(10, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 10:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(11, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 11:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(12, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 12:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(13, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 13:
//			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(14, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 14:
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, 2, kmerSize);
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.25, kmerSize);
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.5, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(15, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 15:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(16, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 16:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(17, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 17:
//			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			writeStage(18, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 18:
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			countReadRepetitiveUnitigs(chunksPerRead, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(19, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadChunkSequences("sequences-chunk19.txt", rawReadLengths, chunksPerRead, sequenceIndex);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeBidirectedUnitigGraph("graph-dbg-final.gfa", "paths-dbg-final.gaf", "unitigs-dbg-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadUnitigSequences("sequences-dbg-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 19:
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 50000, 100000);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(20, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 20:
			writeBidirectedUnitigGraph("graph-final.gfa", "paths-final.gaf", "unitigs-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadUnitigSequences("sequences-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
}

int main(int argc, char** argv)
{
	const size_t numThreads = std::stoi(argv[1]);
	const size_t k = std::stoull(argv[2]);
	const size_t windowSize = std::stoull(argv[3]);
	const size_t middleSkip = std::stoull(argv[4]);
	const double approxOneHapCoverage = std::stod(argv[5]);
	mismatchFraction = std::stod(argv[6]); // try 2-3x average error rate
	const size_t startStage = std::stoull(argv[7]);
	std::vector<std::string> readFiles;
	for (size_t i = 8; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	FastaCompressor::CompressedStringIndex sequenceIndex { 5, 100 };
	std::vector<size_t> readBasepairLengths;
	for (auto file : readFiles)
	{
		std::cerr << "reading from file " << file << std::endl;
		readFilesAndAddToSequenceIndex(file, sequenceIndex, readBasepairLengths, numThreads);
	}
	sequenceIndex.removeConstructionVariables();
	size_t lastReal = readBasepairLengths.size();
	for (size_t i = 0; i < lastReal; i++)
	{
		readBasepairLengths.emplace_back(readBasepairLengths[i]);
	}
	std::cerr << sequenceIndex.size() << " reads" << std::endl;
	makeGraph(sequenceIndex, readBasepairLengths, numThreads, approxOneHapCoverage, k, windowSize, middleSkip, startStage);
}
