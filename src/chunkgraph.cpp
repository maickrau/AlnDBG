#include <map>
#include <queue>
#include <mutex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cxxopts.hpp>
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
#include "OverlapMatcher.h"
#include "PathWalker.h"
#include "TrioKmerCounter.h"

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
				if (end+1 == chunksPerRead[i].size()) break;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end+1]))) break;
				end += 1;
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
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][start-1]))) break;
				start -= 1;
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
	assert(chunksPerRead.size() % 2 == 0);
	assert(rawReadLengths.size() == chunksPerRead.size());
	std::vector<bool> canMergeRight;
	std::vector<bool> canMergeLeft;
	std::vector<size_t> uniqueRight;
	std::vector<size_t> uniqueLeft;
	std::vector<size_t> complementChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			while (canMergeLeft.size() <= chunk) canMergeLeft.emplace_back(true);
		}
	}
	complementChunk.resize(canMergeLeft.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		size_t otheri = i;
		if (i >= chunksPerRead.size()/2)
		{
			otheri = i - chunksPerRead.size()/2;
		}
		else
		{
			otheri = i + chunksPerRead.size()/2;
		}
		assert(chunksPerRead[i].size() == chunksPerRead[otheri].size());
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t otherj = chunksPerRead[otheri].size() - 1 - j;
			assert(otherj < chunksPerRead[otheri].size());
			size_t chunkHere = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			size_t otherChunk = std::get<2>(chunksPerRead[otheri][otherj]) & maskUint64_t;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				assert(NonexistantChunk(std::get<2>(chunksPerRead[otheri][otherj])));
				continue;
			}
			assert(!NonexistantChunk(std::get<2>(chunksPerRead[otheri][otherj])));
			assert(chunkHere < complementChunk.size());
			assert(otherChunk < complementChunk.size());
			assert(complementChunk[chunkHere] == std::numeric_limits<size_t>::max() || complementChunk[chunkHere] == otherChunk);
			assert(complementChunk[otherChunk] == std::numeric_limits<size_t>::max() || complementChunk[otherChunk] == chunkHere);
			complementChunk[chunkHere] = otherChunk;
			complementChunk[otherChunk] = chunkHere;
		}
	}
	for (size_t i = 0; i < complementChunk.size(); i++)
	{
		assert(complementChunk[i] == std::numeric_limits<size_t>::max() || complementChunk[complementChunk[i]] == i);
	}
	canMergeRight.resize(canMergeLeft.size(), true);
	uniqueRight.resize(canMergeLeft.size(), std::numeric_limits<size_t>::max());
	uniqueLeft.resize(canMergeLeft.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i][0])))
		{
			assert(std::get<2>(chunksPerRead[i][0]) & firstBitUint64_t);
			canMergeLeft[std::get<2>(chunksPerRead[i][0]) & maskUint64_t] = false;
		}
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i].back())))
		{
			assert(std::get<2>(chunksPerRead[i].back()) & firstBitUint64_t);
			canMergeRight[std::get<2>(chunksPerRead[i].back()) & maskUint64_t] = false;
		}
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])))
			{
				if (!NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
				{
					assert(std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t);
					canMergeLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				}
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				assert(!NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])));
				assert(std::get<2>(chunksPerRead[i][j-1]) & firstBitUint64_t);
				canMergeRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;
				continue;
			}
			if (std::get<0>(chunksPerRead[i][j-1]) != std::get<0>(chunksPerRead[i][j]))
			{
				canMergeRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;
			}
			else
			{
				assert(std::get<1>(chunksPerRead[i][j-1]) < std::get<1>(chunksPerRead[i][j]));
			}
			if (std::get<1>(chunksPerRead[i][j-1]) != std::get<1>(chunksPerRead[i][j]))
			{
				canMergeLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
			}
			else
			{
				assert(std::get<0>(chunksPerRead[i][j-1]) < std::get<0>(chunksPerRead[i][j]));
			}
			assert(!NonexistantChunk(std::get<2>(chunksPerRead[i][j])));
			assert(!NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])));
			assert(std::get<2>(chunksPerRead[i][j-1]) & firstBitUint64_t);
			if (uniqueRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] == std::numeric_limits<size_t>::max()) uniqueRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = std::get<2>(chunksPerRead[i][j]);
			if (uniqueRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] != std::get<2>(chunksPerRead[i][j])) canMergeRight[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;
			assert(std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t);
			if (uniqueLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] == std::numeric_limits<size_t>::max()) uniqueLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = std::get<2>(chunksPerRead[i][j-1]);
			if (uniqueLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] != std::get<2>(chunksPerRead[i][j-1])) canMergeLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (complementChunk[i] == std::numeric_limits<size_t>::max())
		{
			assert(uniqueRight[i] == std::numeric_limits<size_t>::max());
			assert(uniqueLeft[i] == std::numeric_limits<size_t>::max());
			continue;
		}
		size_t other = complementChunk[i];
		if (other == i) continue; // palindrome!
		assert(canMergeLeft[i] == canMergeRight[other]);
		assert(canMergeRight[i] == canMergeLeft[other]);
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunksHere;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				newChunksHere.emplace_back(chunksPerRead[i][j]);
				continue;
			}
			if (canMergeLeft[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
			{
				assert(j >= 1);
				assert(std::get<0>(chunksPerRead[i][j-1]) < std::get<0>(chunksPerRead[i][j]));
				assert(std::get<1>(chunksPerRead[i][j-1]) == std::get<1>(chunksPerRead[i][j]));
				continue;
			}
			if (canMergeRight[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
			{
				assert(j+1 < chunksPerRead[i].size());
				assert(std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j+1]));
				assert(std::get<1>(chunksPerRead[i][j]) < std::get<1>(chunksPerRead[i][j+1]));
				continue;
			}
			newChunksHere.emplace_back(chunksPerRead[i][j]);
			continue;
		}
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

int approxChunkMatchGoodDiagonal(const size_t leftLength, const size_t rightLength, const phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, size_t>>>& positionsInLeft, const std::vector<std::tuple<size_t, size_t, uint64_t>>& chunks)
{
	// diagonal, leftstart, leftend, rightstart, rightend
	std::vector<std::tuple<int, size_t, size_t, size_t, size_t>> matches;
	for (auto t : chunks)
	{
		if (positionsInLeft.count(std::get<2>(t)) == 0) continue;
		for (auto pos : positionsInLeft.at(std::get<2>(t)))
		{
			matches.emplace_back((int)pos.first - (int)std::get<0>(t), pos.first, pos.second, std::get<0>(t), std::get<1>(t));
		}
	}
	assert(matches.size() >= 1);
	std::sort(matches.begin(), matches.end());
	size_t currentClusterLeftMin = std::get<1>(matches[0]);
	size_t currentClusterLeftMax = std::get<2>(matches[0]);
	size_t currentClusterRightMin = std::get<3>(matches[0]);
	size_t currentClusterRightMax = std::get<4>(matches[0]);
	size_t currentClusterStart = 0;
	int lastDiagonal = std::get<0>(matches[0]);
	int bestDiagonal = std::numeric_limits<int>::max();
	size_t bestDiagonalMatchSize = 0;
	for (size_t i = 0; i < matches.size(); i++)
	{
		auto t = matches[i];
		assert(std::get<0>(t) >= lastDiagonal);
		if (std::get<0>(t) > lastDiagonal + 100)
		{
			if (currentClusterLeftMax - currentClusterLeftMin > 5000 && currentClusterRightMax - currentClusterRightMin > 5000 && (currentClusterLeftMin < 10000 || currentClusterRightMin < 10000) && (currentClusterLeftMax + 10000 > leftLength || currentClusterRightMax + 10000 > rightLength))
			{
				size_t matchLength = std::min(currentClusterLeftMax - currentClusterLeftMin, currentClusterRightMax - currentClusterRightMin);
				if (matchLength > bestDiagonalMatchSize)
				{
					bestDiagonalMatchSize = matchLength;
					bestDiagonal = std::get<0>(matches[currentClusterStart + (i-currentClusterStart)/2]);
				}
			}
			currentClusterLeftMin = std::get<1>(t);
			currentClusterLeftMax = std::get<2>(t);
			currentClusterRightMin = std::get<3>(t);
			currentClusterRightMax = std::get<4>(t);
			lastDiagonal = std::get<0>(t);
			currentClusterStart = i;
		}
		currentClusterLeftMin = std::min(currentClusterLeftMin, std::get<1>(t));
		currentClusterLeftMax = std::max(currentClusterLeftMax, std::get<2>(t));
		currentClusterRightMin = std::min(currentClusterRightMin, std::get<3>(t));
		currentClusterRightMax = std::max(currentClusterRightMax, std::get<4>(t));
		lastDiagonal = std::get<0>(t);
	}
	if (currentClusterLeftMax - currentClusterLeftMin > 5000 && currentClusterRightMax - currentClusterRightMin > 5000 && (currentClusterLeftMin < 10000 || currentClusterRightMin < 10000) && (currentClusterLeftMax + 10000 > leftLength || currentClusterRightMax + 10000 > rightLength))
	{
		size_t matchLength = std::min(currentClusterLeftMax - currentClusterLeftMin, currentClusterRightMax - currentClusterRightMin);
		if (matchLength > bestDiagonalMatchSize)
		{
			bestDiagonalMatchSize = matchLength;
			bestDiagonal = std::get<0>(matches[currentClusterStart + (matches.size()-currentClusterStart)/2]);
		}
	}
	return bestDiagonal;
}

std::pair<std::vector<size_t>, std::vector<int>> filterGoodChunkMatches(const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunks, const phmap::flat_hash_set<size_t>& unfiltered, const size_t refRead)
{
	phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, size_t>>> positionsInRef;
	for (size_t i = 0; i < chunks[refRead].size(); i++)
	{
		positionsInRef[std::get<2>(chunks[refRead][i])].emplace_back(std::get<0>(chunks[refRead][i]), std::get<1>(chunks[refRead][i]));
	}
	std::pair<std::vector<size_t>, std::vector<int>> result;
	for (auto read : unfiltered)
	{
		int diagonal = approxChunkMatchGoodDiagonal(rawReadLengths[refRead], rawReadLengths[read], positionsInRef, chunks[read]);
		if (diagonal == std::numeric_limits<int>::max()) continue;
		result.first.emplace_back(read);
		result.second.emplace_back(diagonal);
	}
	return result;
}

phmap::flat_hash_set<size_t> getRoughMatchingReads(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunks, const std::vector<std::vector<size_t>>& readsPerChunk, const size_t refRead, const size_t reverseReadOffset)
{
	phmap::flat_hash_set<size_t> result;
	for (size_t j = 0; j < chunks[refRead].size(); j++)
	{
		size_t chunk = std::get<2>(chunks[refRead][j]) & maskUint64_t;
		assert(chunk < readsPerChunk.size());
		if (readsPerChunk[chunk].size() < 1000)
		{
			result.insert(readsPerChunk[chunk].begin(), readsPerChunk[chunk].end());
		}
	}
	if (result.count(refRead) == 1) result.erase(refRead);
	if (result.count(refRead + reverseReadOffset) == 1) result.erase(refRead + reverseReadOffset);
	return result;
}

std::vector<std::vector<bool>> getGoodKmersFromAlignments(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const size_t kmerSize, const size_t windowSize, const size_t middleSkip)
{
	const size_t blockSize = 2000;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "getting good kmers" << std::endl;
	std::vector<std::vector<bool>> kmerIsGood;
	kmerIsGood.resize(sequenceIndex.size());
	for (size_t i = 0; i < kmerIsGood.size(); i++)
	{
		kmerIsGood[i].resize(rawReadLengths[i], true);
	}
	auto chunks = getGapEnclosingChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, 1000, 3000);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	splitPerFirstLastKmers(sequenceIndex, chunks, kmerSize, numThreads);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	splitPerLength(chunks, numThreads);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	splitPerBaseCounts(sequenceIndex, rawReadLengths, chunks, numThreads);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//	splitPerSequenceIdentityRoughly(sequenceIndex, rawReadLengths, chunks, numThreads);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//	splitPerSequenceIdentity(sequenceIndex, rawReadLengths, chunks, numThreads);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	assert(chunks.size() == sequenceIndex.size()*2);
	assert(rawReadLengths.size() == sequenceIndex.size()*2);
	assert(chunks.size() == rawReadLengths.size());
	std::vector<std::vector<size_t>> readsPerChunk;
	for (size_t i = 0; i < chunks.size(); i++)
	{
		for (size_t j = 0; j < chunks[i].size(); j++)
		{
			size_t chunk = std::get<2>(chunks[i][j]) & maskUint64_t;
			while (readsPerChunk.size() <= chunk) readsPerChunk.emplace_back();
			readsPerChunk[chunk].emplace_back(i);
		}
	}
	for (size_t i = 0; i < readsPerChunk.size(); i++)
	{
		phmap::flat_hash_set<size_t> uniqs { readsPerChunk[i].begin(), readsPerChunk[i].end() };
		readsPerChunk[i].clear();
		readsPerChunk[i].insert(readsPerChunk[i].end(), uniqs.begin(), uniqs.end());
	}
	std::vector<size_t> readBelongsToBlock;
	std::vector<std::vector<size_t>> readsWithinBlock;
	{
		std::vector<size_t> smallestNonsingularChunkInRead;
		for (size_t i = 0; i < sequenceIndex.size(); i++)
		{
			assert(chunks[i].size() == chunks[i+sequenceIndex.size()].size());
			size_t smallestHere = std::numeric_limits<size_t>::max();
			for (size_t j = 0; j < chunks[i].size(); j++)
			{
				size_t chunk = std::get<2>(chunks[i][j]) & maskUint64_t;
				if (readsPerChunk[chunk].size() < 2) continue;
				smallestHere = std::min(smallestHere, chunk);
			}
			for (size_t j = 0; j < chunks[i + sequenceIndex.size()].size(); j++)
			{
				size_t chunk = std::get<2>(chunks[i + sequenceIndex.size()][j]) & maskUint64_t;
				if (readsPerChunk[chunk].size() < 2) continue;
				smallestHere = std::min(smallestHere, chunk);
			}
			smallestNonsingularChunkInRead.emplace_back(smallestHere);
		}
		std::vector<size_t> readOrder;
		for (size_t i = 0; i < sequenceIndex.size(); i++)
		{
			readOrder.emplace_back(i);
		}
		std::sort(readOrder.begin(), readOrder.end(), [&smallestNonsingularChunkInRead](const size_t left, const size_t right) { return smallestNonsingularChunkInRead[left] < smallestNonsingularChunkInRead[right]; });
		readsWithinBlock.resize((sequenceIndex.size() + blockSize - 1) / blockSize);
		readBelongsToBlock.resize(sequenceIndex.size(), std::numeric_limits<size_t>::max());
		for (size_t i = 0; i < sequenceIndex.size(); i++)
		{
			assert(readBelongsToBlock[readOrder[i]] == std::numeric_limits<size_t>::max());
			readBelongsToBlock[readOrder[i]] = i / blockSize;
			readsWithinBlock[i / blockSize].emplace_back(readOrder[i]);
		}
	}
	assert(readsWithinBlock.size() >= 1);
	std::mutex printMutex;
	for (size_t block = 0; block < readsWithinBlock.size(); block++)
	{
		std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
		std::cerr << "begin block " << block << " / " << readsWithinBlock.size() << " reads " << readsWithinBlock[block].size() << std::endl;
		auto startTime = getTime();
		std::vector<ReadOverlapInformation> readInformationHere;
		readInformationHere.resize(readsWithinBlock[block].size());
		std::vector<ReadOverlapInformation> reverseReadInformationHere;
		reverseReadInformationHere.resize(readsWithinBlock[block].size());
		phmap::flat_hash_map<size_t, size_t> readIndexWithinBlock;
		for (size_t i = 0; i < readsWithinBlock[block].size(); i++)
		{
			assert(readIndexWithinBlock.count(readsWithinBlock[block][i]) == 0);
			readIndexWithinBlock[readsWithinBlock[block][i]] = i;
		}
		auto getReadInfos = getTime();
		iterateMultithreaded(0, readsWithinBlock[block].size(), numThreads, [&sequenceIndex, &readsWithinBlock, &reverseReadInformationHere, &readInformationHere, block, kmerSize](const size_t index)
		{
			readInformationHere[index] = getReadOverlapInformation(sequenceIndex, readsWithinBlock[block][index], kmerSize);
			reverseReadInformationHere[index] = getReadOverlapInformation(sequenceIndex, readsWithinBlock[block][index] + sequenceIndex.size(), kmerSize);
		});
		auto afterReadInfos = getTime();
		std::atomic<size_t> totalMatchesWithinBlock;
		totalMatchesWithinBlock = 0;
		iterateMultithreaded(0, readsWithinBlock[block].size(), numThreads, [&sequenceIndex, &readsPerChunk, &chunks, &kmerIsGood, &rawReadLengths, &totalMatchesWithinBlock, &readIndexWithinBlock, &readsWithinBlock, &printMutex, &readInformationHere, &reverseReadInformationHere, &readBelongsToBlock, block](const size_t readIndex)
		{
			std::string readName = sequenceIndex.getName(readsWithinBlock[block][readIndex]);
			{
				std::lock_guard<std::mutex> lock { printMutex };
				std::cerr << "block " << block << " begin overlapping read " << readIndex << "/" << sequenceIndex.size() << " id " << readsWithinBlock[block][readIndex] << " name " << readName << std::endl;
			}
			auto readStartTime = getTime();
			phmap::flat_hash_set<size_t> fwmatchingReadsBeforeFilter = getRoughMatchingReads(chunks, readsPerChunk, readsWithinBlock[block][readIndex], sequenceIndex.size());
			phmap::flat_hash_set<size_t> fwMatchingReads;
			for (auto t : fwmatchingReadsBeforeFilter)
			{
				size_t check = t < sequenceIndex.size() ? t : t - sequenceIndex.size();
				if (readBelongsToBlock[check] != block) continue;
				fwMatchingReads.emplace(t);
			}
			std::vector<size_t> filteredFwMatchingReads;
			std::vector<int> filteredFwMatchingDiagonals;
			std::tie(filteredFwMatchingReads, filteredFwMatchingDiagonals) = filterGoodChunkMatches(rawReadLengths, chunks, fwMatchingReads, readsWithinBlock[block][readIndex]);
			size_t countFwMatchingReads = 0;
			for (size_t k = 0; k < filteredFwMatchingReads.size(); k++)
			{
				const size_t otherRead = filteredFwMatchingReads[k];
				size_t lastPos = std::numeric_limits<size_t>::max();
				ReadOverlapInformation& otherReadInfo = (otherRead < sequenceIndex.size()) ? (readInformationHere[readIndexWithinBlock[otherRead]]) : (reverseReadInformationHere[readIndexWithinBlock[otherRead - sequenceIndex.size()]]);
				iterateUniqueKmerMatches(readInformationHere[readIndex], otherReadInfo, filteredFwMatchingDiagonals[k], [&kmerIsGood, &lastPos, &countFwMatchingReads, &readsWithinBlock, block, readIndex](const size_t otherReadIndex, const size_t thisReadPos, const size_t otherReadPos)
				{
					if (lastPos != std::numeric_limits<size_t>::max())
					{
						assert(thisReadPos > lastPos);
						for (size_t j = lastPos+1; j < thisReadPos; j++)
						{
							kmerIsGood[readsWithinBlock[block][readIndex]][j] = false;
						}
					}
					else
					{
						countFwMatchingReads += 1;
					}
					lastPos = thisReadPos;
				});
			}
			auto readEndTime = getTime();
			totalMatchesWithinBlock += countFwMatchingReads;
			{
				std::lock_guard<std::mutex> lock { printMutex };
				std::cerr << "block " << block << " read " << readIndex << "/" << readsWithinBlock[block].size() << " id " << readsWithinBlock[block][readIndex] << " name " << readName << " length " << rawReadLengths[readsWithinBlock[block][readIndex]] << " has " << countFwMatchingReads << " matches, time " << formatTime(readStartTime, readEndTime) << ", pre-filter overlaps " << fwMatchingReads.size() << ", pre-kmer overlaps " << filteredFwMatchingReads.size() << std::endl;
			}
		});
		auto endTime = getTime();
		std::cerr << "finished block " << block << " / " << readsWithinBlock.size() << ", total matches " << totalMatchesWithinBlock << ", total time: " << formatTime(startTime, endTime) << std::endl;
		std::cerr << "times: " << formatTime(startTime, getReadInfos) << " " << formatTime(getReadInfos, afterReadInfos) << " " << formatTime(afterReadInfos, endTime) << std::endl;
	}
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "begin overlapping" << std::endl;
	iterateMultithreaded(0, sequenceIndex.size(), numThreads, [&sequenceIndex, &readsPerChunk, &chunks, &kmerIsGood, &rawReadLengths, &printMutex, &readBelongsToBlock, kmerSize](const size_t i)
	{
		auto startTime = getTime();
		std::string readName = sequenceIndex.getName(i);
		{
			std::lock_guard<std::mutex> lock { printMutex };
			std::cerr << "begin overlapping read " << i << "/" << sequenceIndex.size() << " name " << readName << std::endl;
		}
		phmap::flat_hash_set<size_t> fwmatchingReadsBeforeFilter = getRoughMatchingReads(chunks, readsPerChunk, i, sequenceIndex.size());
		phmap::flat_hash_set<size_t> fwmatchingReads;
		for (auto t : fwmatchingReadsBeforeFilter)
		{
			size_t check = t < sequenceIndex.size() ? t : t - sequenceIndex.size();
			if (readBelongsToBlock[i] == readBelongsToBlock[check]) continue;
			fwmatchingReads.emplace(t);
		}
		std::vector<size_t> filteredFwMatchingReads;
		std::vector<int> filteredFwMatchingDiagonals;
		std::tie(filteredFwMatchingReads, filteredFwMatchingDiagonals) = filterGoodChunkMatches(rawReadLengths, chunks, fwmatchingReads, i);
		size_t lastRead = std::numeric_limits<size_t>::max();
		size_t lastPos = 0;
		size_t countFwMatchingReads = 0;
		iterateUniqueKmerMatches(sequenceIndex, i, filteredFwMatchingReads, filteredFwMatchingDiagonals, kmerSize, [&kmerIsGood, &lastRead, &lastPos, &countFwMatchingReads, i](const size_t otherReadIndex, const size_t thisReadPos, const size_t otherReadPos)
		{
			if (lastRead == otherReadIndex)
			{
				assert(thisReadPos > lastPos);
				for (size_t j = lastPos+1; j < thisReadPos; j++)
				{
					kmerIsGood[i][j] = false;
				}
			}
			else
			{
				countFwMatchingReads += 1;
			}
			lastRead = otherReadIndex;
			lastPos = thisReadPos;
		});
		auto endTime = getTime();
		{
			std::lock_guard<std::mutex> lock { printMutex };
			std::cerr << "read " << i << " name " << readName << " length " << rawReadLengths[i] << " has " << countFwMatchingReads << " matches, time " << formatTime(startTime, endTime) << ", pre-filter overlaps " << fwmatchingReads.size() << ", pre-kmer overlaps " << filteredFwMatchingReads.size() << std::endl;
		}
	});
	return kmerIsGood;
}

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getCorrectnessWeightedChunksPerRead(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const size_t kmerSize, const size_t windowSize, const size_t middleSkip)
{
	std::vector<std::vector<bool>> kmerIsGood = getGoodKmersFromAlignments(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize, middleSkip);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "getting chunks with good kmers" << std::endl;
	auto result = getWeightedMinimizerBoundedChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize, kmerIsGood);
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			assert((std::get<2>(result[i][j]) & firstBitUint64_t) == 0);
			std::get<2>(result[i][j]) |= firstBitUint64_t;
		}
	}
	return result;
}

void makeGraph(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const TrioKmerCounter& trioHapmers, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize, const size_t windowSize, const size_t startStage)
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
			//chunksPerRead = getCorrectnessWeightedChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize, middleSkip);
			chunksPerRead = getMinimizerBoundedChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize);
			for (size_t i = 0; i < chunksPerRead.size(); i++)
			{
				for (size_t j = 0; j < chunksPerRead[i].size(); j++)
				{
					assert((std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t) == 0);
					std::get<2>(chunksPerRead[i][j]) |= firstBitUint64_t;
				}
			}
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
			writeStage(1, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 1:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerFirstLastKmers(sequenceIndex, chunksPerRead, kmerSize, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerLength(chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			mergeFakeBubbles(chunksPerRead, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(2, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 2:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerBaseCounts(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(3, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 3:
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			writeStage(4, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 4:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentityRoughly(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			// splitPerSequenceIdentity(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			// std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(5, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 5:
			if (trioHapmers.notEmpty())
			{
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				splitPerTrioPhasing(sequenceIndex, rawReadLengths, chunksPerRead, trioHapmers, numThreads);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				writeStage(6, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
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
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			writeStage(6, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 6:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(7, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 7:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(8, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 8:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerCorrectedKmerClustering(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			writeStage(81, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 81:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerCorrectedKmerPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			writeStage(82, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 82:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSNPTransitiveClosureClustering(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(83, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 83:
			if (trioHapmers.notEmpty())
			{
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				splitPerTrioPhasing(sequenceIndex, rawReadLengths, chunksPerRead, trioHapmers, numThreads);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				writeStage(9, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
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
			if (trioHapmers.notEmpty())
			{
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				splitPerTrioPhasing(sequenceIndex, rawReadLengths, chunksPerRead, trioHapmers, numThreads);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				writeStage(92, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
			[[fallthrough]];
		case 92:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSNPTransitiveClosureClustering(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(93, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 93:
			if (trioHapmers.notEmpty())
			{
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				splitPerTrioPhasing(sequenceIndex, rawReadLengths, chunksPerRead, trioHapmers, numThreads);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				writeStage(94, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
			[[fallthrough]];
		case 94:
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
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			getContigPathsAndConsensuses(chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize, "paths-final.txt", "contigs-final.fa", numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 50000, 100000);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			writeStage(20, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
//			[[fallthrough]];
//		case 20:
//			writeBidirectedUnitigGraph("graph-final.gfa", "paths-final.gaf", "unitigs-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			writeReadUnitigSequences("sequences-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
}

int main(int argc, char** argv)
{
	std::cerr << "chunkgraph " << VERSION << std::endl;
	cxxopts::Options options { "chunkgraph" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input reads. Multiple files can be input with -i file1.fa -i file2.fa etc (required)", cxxopts::value<std::vector<std::string>>())
		("t,threads", "Number of threads", cxxopts::value<size_t>()->default_value("1"))
		("k", "K-mer size", cxxopts::value<size_t>()->default_value("11"))
		("w", "Window size", cxxopts::value<size_t>()->default_value("5000"))
		("avg-hap-coverage", "Average single haplotype coverage", cxxopts::value<double>())
		("max-error-rate", "Maximum error rate. Try 2-3x average read error rate", cxxopts::value<double>()->default_value("0.03"))
		("restart", "Restart from stage", cxxopts::value<size_t>())
		("parent-1-reads", "Trio reads from parent 1", cxxopts::value<std::string>())
		("parent-2-reads", "Trio reads from parent 2", cxxopts::value<std::string>())
		;
	auto params = options.parse(argc, argv);
	if (params.count("v") == 1)
	{
		std::cerr << "Version: " << VERSION << std::endl;
		exit(0);
	}
	if (params.count("h") == 1)
	{
		std::cerr << options.help();
		exit(0);
	}
	bool paramError = false;
	if (params.count("avg-hap-coverage") == 0)
	{
		std::cerr << "Average single haplotype coverage (--avg-hap-coverage) is required" << std::endl;
		paramError = true;
	}
	if (params.count("i") == 0)
	{
		std::cerr << "Input reads are required" << std::endl;
		paramError = true;
	}
	if (params.count("parent-1-reads") == 1 && params.count("parent-2-reads") == 0)
	{
		std::cerr << "Parent 2 reads missing" << std::endl;
		paramError = true;
	}
	if (params.count("parent-1-reads") == 0 && params.count("parent-2-reads") == 1)
	{
		std::cerr << "Parent 1 reads missing" << std::endl;
		paramError = true;
	}
	if (paramError) std::abort();
	const size_t numThreads = params["t"].as<size_t>();
	const size_t k = params["k"].as<size_t>();
	const size_t windowSize = params["w"].as<size_t>();
	const double approxOneHapCoverage = params["avg-hap-coverage"].as<double>();
	mismatchFraction = params["max-error-rate"].as<double>();
	size_t startStage = 0;
	if (params.count("restart") == 1) startStage = params["restart"].as<size_t>();
	std::vector<std::string> readFiles = params["i"].as<std::vector<std::string>>();
	TrioKmerCounter trioHapmers;
	if (params.count("parent-1-reads") == 1)
	{
		assert(params.count("parent-2-reads") == 1);
		std::cerr << "getting parent-specific kmers" << std::endl;
		std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
		trioHapmers.initialize(params["parent-1-reads"].as<std::string>(), params["parent-2-reads"].as<std::string>(), 31, 20);
		std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
		std::cerr << trioHapmers.hap1KmerCount() << " parent 1 specific kmers" << std::endl;
		std::cerr << trioHapmers.hap2KmerCount() << " parent 2 specific kmers" << std::endl;
	}
	FastaCompressor::CompressedStringIndex sequenceIndex { 5, 100 };
	std::vector<size_t> readBasepairLengths;
	for (auto file : readFiles)
	{
		std::cerr << "reading from file " << file << std::endl;
		readFilesAndAddToSequenceIndex(file, sequenceIndex, readBasepairLengths, numThreads);
	}
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "postprocessing sequence index" << std::endl;
	sequenceIndex.removeConstructionVariables(numThreads);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	size_t lastReal = readBasepairLengths.size();
	for (size_t i = 0; i < lastReal; i++)
	{
		readBasepairLengths.emplace_back(readBasepairLengths[i]);
	}
	std::cerr << sequenceIndex.size() << " reads" << std::endl;
	makeGraph(sequenceIndex, readBasepairLengths, trioHapmers, numThreads, approxOneHapCoverage, k, windowSize, startStage);
}
