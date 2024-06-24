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

double mismatchFraction;

size_t findBiggestSmallerIndex(const std::vector<size_t>& nums, const size_t value)
{
	if (nums.size() < 5)
	{
		assert(nums[0] < value);
		for (size_t i = 1; i < nums.size(); i++)
		{
			if (nums[i] > value) return i-1;
		}
		return nums.size()-1;
	}
	size_t low = 0;
	size_t high = nums.size();
	while (high > low)
	{
		size_t mid = (low + high) / 2;
		assert(nums[mid] != value);
		if (nums[mid] < value)
		{
			low = mid+1;
		}
		else if (nums[mid] > value)
		{
			high = mid;
		}
	}
	assert(low == high);
	assert(low <= nums.size());
	assert(low == nums.size() || nums[low] > value);
	assert(low > 0);
	assert(nums[low-1] < value);
	return low-1;
}

void splitPerLength(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by length" << std::endl;
	const double differenceFraction = 0.02;
	const size_t differenceConstant = 50;
	std::vector<std::vector<size_t>> lengthsPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			while ((std::get<2>(t) & maskUint64_t) >= lengthsPerChunk.size()) lengthsPerChunk.emplace_back();
			lengthsPerChunk[std::get<2>(t) & maskUint64_t].emplace_back(std::get<1>(t) - std::get<0>(t));
		}
	}
	std::vector<std::vector<size_t>> splitters;
	splitters.resize(lengthsPerChunk.size());
	iterateMultithreaded(0, lengthsPerChunk.size(), numThreads, [&lengthsPerChunk, &splitters, differenceFraction, differenceConstant](const size_t i)
	{
		std::sort(lengthsPerChunk[i].begin(), lengthsPerChunk[i].end());
		splitters[i].emplace_back(0);
		for (size_t j = 1; j < lengthsPerChunk[i].size(); j++)
		{
			size_t distance = lengthsPerChunk[i][j] - lengthsPerChunk[i][j-1];
			if (distance < std::max((size_t)(lengthsPerChunk[i][j] * differenceFraction), (size_t)differenceConstant)) continue;
			splitters[i].emplace_back((lengthsPerChunk[i][j] + lengthsPerChunk[i][j-1])/2);
		}
	});
	std::vector<std::vector<size_t>> splitterToNode;
	splitterToNode.resize(lengthsPerChunk.size());
	size_t nextNum = 0;
	for (size_t i = 0; i < splitters.size(); i++)
	{
		for (size_t j = 0; j < splitters[i].size(); j++)
		{
			splitterToNode[i].emplace_back(nextNum);
			nextNum += 1;
		}
	}
	iterateMultithreaded(0, chunksPerRead.size(), numThreads, [&splitters, &splitterToNode, &chunksPerRead](const size_t i)
	{
		for (auto& t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			assert(std::get<1>(t) > std::get<0>(t));
			size_t distance = std::get<1>(t) - std::get<0>(t);
			assert((std::get<2>(t) & maskUint64_t) < splitters.size());
			size_t index = findBiggestSmallerIndex(splitters[std::get<2>(t) & maskUint64_t], distance);
			assert(index < splitterToNode[std::get<2>(t) & maskUint64_t].size());
			std::get<2>(t) = (std::get<2>(t) & firstBitUint64_t) + splitterToNode[std::get<2>(t) & maskUint64_t][index];
		}
	});
}

std::vector<size_t> getMinHashes(const std::string& sequence, const size_t k, const size_t count)
{
	assert(sequence.size() > k + count);
	uint64_t kmer = 0;
	uint64_t mask = (1ull << 2ull*k) - 1;
	for (size_t i = 0; i < k; i++)
	{
		kmer <<= 2;
		switch(sequence[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
		}
	}
	std::priority_queue<size_t> queue;
	queue.emplace(hashFwAndBw(kmer, k));
	for (size_t i = k; i < sequence.size(); i++)
	{
		kmer <<= 2;
		switch(sequence[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
		}
		kmer &= mask;
		size_t h = hashFwAndBw(kmer, k);
		if (queue.size() < count)
		{
			queue.emplace(h);
		}
		else if (h < queue.top())
		{
			queue.pop();
			queue.emplace(h);
		}
	}
	std::vector<size_t> result;
	while (queue.size() > 0)
	{
		result.emplace_back(queue.top());
		queue.pop();
	}
	return result;
}

void splitPerFirstLastKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	assert(kmerSize <= 31);
	size_t nextNum = 0;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> endKmersToNumber;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		std::string readseq = getReadSequence(sequenceIndex, i);
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<0>(chunksPerRead[i][j]) < std::get<1>(chunksPerRead[i][j]));
			assert(std::get<1>(chunksPerRead[i][j]) < readseq.size());
			uint64_t firstKmer = 0;
			uint64_t lastKmer = 0;
			for (size_t k = 0; k < kmerSize; k++)
			{
				firstKmer <<= 2;
				lastKmer <<= 2;
				switch(readseq[std::get<0>(chunksPerRead[i][j])+k])
				{
					case 'A':
						firstKmer += 0;
						break;
					case 'C':
						firstKmer += 1;
						break;
					case 'G':
						firstKmer += 2;
						break;
					case 'T':
						firstKmer += 3;
						break;
					default:
						assert(false);
				}
				switch(readseq[std::get<1>(chunksPerRead[i][j])-k])
				{
					case 'A':
						lastKmer += 3;
						break;
					case 'C':
						lastKmer += 2;
						break;
					case 'G':
						lastKmer += 1;
						break;
					case 'T':
						lastKmer += 0;
						break;
					default:
						assert(false);
				}
			}
			auto key = std::make_pair(firstKmer, lastKmer);
			if (endKmersToNumber.count(key) == 0)
			{
				endKmersToNumber[key] = nextNum;
				nextNum += 1;
			}
			std::get<2>(chunksPerRead[i][j]) = endKmersToNumber.at(key) + firstBitUint64_t;
		}
	}
}

void splitPerBaseCounts(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t mismatchFloor = 10;
	std::cerr << "splitting by base counts" << std::endl;
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, mismatchFloor](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		std::vector<std::vector<size_t>> countsPerOccurrence;
		countsPerOccurrence.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			std::string chunkSequence = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
			assert(chunkSequence.size() == std::get<1>(t)-std::get<0>(t)+1);
			countsPerOccurrence[j].resize(4);
			for (size_t k = 0; k < chunkSequence.size(); k++)
			{
				switch(chunkSequence[k])
				{
				case 'A':
					countsPerOccurrence[j][0] += 1;
					break;
				case 'C':
					countsPerOccurrence[j][1] += 1;
					break;
				case 'G':
					countsPerOccurrence[j][2] += 1;
					break;
				case 'T':
					countsPerOccurrence[j][3] += 1;
					break;
				default:
					assert(false);
				}
			}
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				size_t maxEdits = std::min(std::get<1>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) - std::get<0>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]), std::get<1>(chunksPerRead[occurrencesPerChunk[i][k].first][occurrencesPerChunk[i][k].second]) - std::get<0>(chunksPerRead[occurrencesPerChunk[i][k].first][occurrencesPerChunk[i][k].second]));
				maxEdits *= mismatchFraction;
				maxEdits = std::max(maxEdits, mismatchFloor);
				size_t edits = 0;
				for (size_t c = 0; c < 4; c++)
				{
					if (countsPerOccurrence[j][c] > countsPerOccurrence[k][c])
					{
						edits += countsPerOccurrence[j][c] - countsPerOccurrence[k][c];
					}
					else
					{
						edits += countsPerOccurrence[k][c] - countsPerOccurrence[j][c];
					}
				}
				if (edits < maxEdits) merge(parent, j, k);
			}
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "base count splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
}

void splitPerMinHashes(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t kmerSize = 11;
	const size_t hashCount = 10;
	std::cerr << "splitting by minhash" << std::endl;
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, kmerSize, hashCount](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		phmap::flat_hash_map<size_t, size_t> parent;
		std::vector<size_t> oneHashPerLocation;
		size_t smallestSequence = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			smallestSequence = std::min(smallestSequence, std::get<1>(t)-std::get<0>(t)+1);
		}
		if (smallestSequence <= kmerSize + hashCount + 1)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			std::vector<uint64_t> minHashes;
			auto chunkSequence = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
			minHashes = getMinHashes(chunkSequence, kmerSize, hashCount);
			assert(minHashes.size() >= 1);
			for (auto hash : minHashes)
			{
				if (parent.count(hash) == 0) parent[hash] = hash;
			}
			for (auto hash : minHashes)
			{
				merge(parent, hash, minHashes[0]);
			}
			oneHashPerLocation.emplace_back(minHashes[0]);
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> clusterToNode;
			for (auto hash : oneHashPerLocation)
			{
				if (clusterToNode.count(find(parent, hash)) == 1) continue;
				clusterToNode[find(parent, hash)] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "minhash splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, oneHashPerLocation[j])) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
}

void resolveVerySmallNodes(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const size_t maxResolveSize, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving very small nodes" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	phmap::flat_hash_set<size_t> resolvableTinies;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.chunksInUnitig[i].size() != 1) continue;
		if (graph.unitigLengths[i] >= maxResolveSize) continue;
		if (graph.edges.getEdges(std::make_pair(i, true)).size() < 2 && graph.edges.getEdges(std::make_pair(i, false)).size() < 2) continue;
		resolvableTinies.insert(graph.chunksInUnitig[i][0] & maskUint64_t);
	}
	size_t nextNum = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (!NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) nextNum = std::max(nextNum, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			if (resolvableTinies.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			if (std::get<0>(chunksPerRead[i][j]) == 0)
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (std::get<1>(chunksPerRead[i][j])+1 == rawReadLengths[i])
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1])+1)
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (j+1 < chunksPerRead[i].size() && std::get<1>(chunksPerRead[i][j])+1 == std::get<1>(chunksPerRead[i][j-1]))
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
		}
	}
	for (size_t chunk : resolvableTinies)
	{
		std::cerr << "resolve tiny chunk " << chunk << std::endl;
	}
	nextNum += 1;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeToNumber;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<size_t> eraseIndices;
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunks;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (resolvableTinies.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			eraseIndices.emplace_back(j);
			if (j > 0)
			{
				uint64_t here = std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t;
				uint64_t neighbor = std::get<2>(chunksPerRead[i][j-1]) ^ firstBitUint64_t;
				if (!NonexistantChunk(neighbor))
				{
					std::pair<uint64_t, uint64_t> key { here, neighbor };
					if (edgeToNumber.count(key) == 0)
					{
						edgeToNumber[key] = nextNum;
						nextNum += 1;
					}
					uint64_t node = edgeToNumber.at(key);
					newChunks.emplace_back(std::get<0>(chunksPerRead[i][j])-1, std::get<1>(chunksPerRead[i][j]), node);
				}
			}
			if (j+1 < chunksPerRead[i].size())
			{
				uint64_t here = std::get<2>(chunksPerRead[i][j]);
				uint64_t neighbor = std::get<2>(chunksPerRead[i][j+1]);
				if (!NonexistantChunk(neighbor))
				{
					std::pair<uint64_t, uint64_t> key { here, neighbor };
					if (edgeToNumber.count(key) == 0)
					{
						edgeToNumber[key] = nextNum;
						nextNum += 1;
					}
					uint64_t node = edgeToNumber.at(key) + firstBitUint64_t;
					newChunks.emplace_back(std::get<0>(chunksPerRead[i][j]), std::get<1>(chunksPerRead[i][j])+1, node);
				}
			}
		}
		if (eraseIndices.size() == 0 && newChunks.size() == 0) continue;
		for (size_t j = eraseIndices.size()-1; j < eraseIndices.size(); j--)
		{
			size_t index = eraseIndices[j];
			chunksPerRead[i].erase(chunksPerRead[i].begin()+index);
		}
		chunksPerRead[i].insert(chunksPerRead[i].end(), newChunks.begin(), newChunks.end());
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
	}
}

void resolveSemiAmbiguousUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving semiambiguous resolvable unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> tripletsPerUnitig;
	tripletsPerUnitig.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				uint64_t prev, after;
				if (mid & firstBitUint64_t)
				{
					prev = readPaths[i][j].path[k-1];
					after = readPaths[i][j].path[k+1];
				}
				else
				{
					prev = readPaths[i][j].path[k+1] ^ firstBitUint64_t;
					after = readPaths[i][j].path[k-1] ^ firstBitUint64_t;
				}
				tripletsPerUnitig[mid & maskUint64_t][std::make_pair(prev, after)] += 1;
			}
		}
	}
	std::vector<bool> canResolveUnitig;
	canResolveUnitig.resize(graph.unitigLengths.size(), false);
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> prevToAllele;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> afterToAllele;
	prevToAllele.resize(graph.unitigLengths.size());
	afterToAllele.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (tripletsPerUnitig[i].size() < 2) continue;
		phmap::flat_hash_set<uint64_t> coveredPrev;
		phmap::flat_hash_set<uint64_t> coveredAfter;
		phmap::flat_hash_map<uint64_t, uint64_t> parent;
		bool valid = true;
		for (auto triplet : tripletsPerUnitig[i])
		{
			if (triplet.second < 2)
			{
				valid = false;
				break;
			}
			coveredPrev.insert(triplet.first.first);
			coveredAfter.insert(triplet.first.second);
			assert((triplet.first.first & (firstBitUint64_t >> 1)) == 0);
			uint64_t prev = triplet.first.first + (firstBitUint64_t >> 1);
			uint64_t after = triplet.first.second;
			if (parent.count(prev) == 0) parent[prev] = prev;
			if (parent.count(after) == 0) parent[after] = after;
			merge(parent, prev, after);
		}
		if (!valid) continue;
		assert(coveredPrev.size() <= graph.edges.getEdges(std::make_pair(i, false)).size());
		assert(coveredAfter.size() <= graph.edges.getEdges(std::make_pair(i, true)).size());
		if (coveredPrev.size() < graph.edges.getEdges(std::make_pair(i, false)).size()) continue;
		if (coveredAfter.size() < graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		phmap::flat_hash_map<uint64_t, size_t> clusterToAllele;
		size_t alleleNum = 0;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t after = triplet.first.second;
			uint64_t cluster = find(parent, after);
			if (clusterToAllele.count(cluster) == 0)
			{
				clusterToAllele[cluster] = alleleNum;
				alleleNum += 1;
			}
		}
		assert(clusterToAllele.size() == alleleNum);
		assert(clusterToAllele.size() >= 1);
		if (clusterToAllele.size() == 1) continue;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t prev = triplet.first.first;
			uint64_t after = triplet.first.second;
			prevToAllele[i][prev] = clusterToAllele.at(find(parent, after));
			afterToAllele[i][after] = clusterToAllele.at(find(parent, after));
		}
		canResolveUnitig[i] = true;
	}
	phmap::flat_hash_set<size_t> readsAllowedToDiscard;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				if (!canResolveUnitig[mid & maskUint64_t]) continue;
				uint64_t prev = std::numeric_limits<size_t>::max();
				uint64_t after = std::numeric_limits<size_t>::max();
				if (k > 0)
				{
					prev = readPaths[i][j].path[k-1];
				}
				else if (j > 0)
				{
					prev = readPaths[i][j-1].path.back();
				}
				if (k+1 < readPaths[i][j].path.size())
				{
					after = readPaths[i][j].path[k+1];
				}
				else if (j+1 < readPaths[i].size())
				{
					after = readPaths[i][j+1].path[0];
				}
				if ((mid ^ firstBitUint64_t) & firstBitUint64_t)
				{
					std::swap(prev, after);
					if (prev != std::numeric_limits<size_t>::max()) prev ^= firstBitUint64_t;
					if (after != std::numeric_limits<size_t>::max()) after ^= firstBitUint64_t;
				}
				if (prev == std::numeric_limits<size_t>::max() && after == std::numeric_limits<size_t>::max())
				{
					if (readPaths[i].size() == 1 && readPaths[i][0].path.size() == 1)
					{
						readsAllowedToDiscard.insert(i);
						continue;
					}
					else
					{
						canResolveUnitig[mid & maskUint64_t] = false;
						continue;
					}
				}
				size_t prevAllele = std::numeric_limits<size_t>::max();
				size_t afterAllele = std::numeric_limits<size_t>::max();
				if (prev != std::numeric_limits<size_t>::max() && prevToAllele[mid & maskUint64_t].count(prev) == 1)
				{
					prevAllele = prevToAllele[mid & maskUint64_t].at(prev);
				}
				if (after != std::numeric_limits<size_t>::max() && afterToAllele[mid & maskUint64_t].count(after) == 1)
				{
					afterAllele = afterToAllele[mid & maskUint64_t].at(after);
				}
				if (prevAllele == std::numeric_limits<size_t>::max() && afterAllele == std::numeric_limits<size_t>::max())
				{
					canResolveUnitig[mid & maskUint64_t] = false;
					continue;
				}
				if (prevAllele != std::numeric_limits<size_t>::max() && afterAllele != std::numeric_limits<size_t>::max() && prevAllele != afterAllele)
				{
					canResolveUnitig[mid & maskUint64_t] = false;
					continue;
				}
			}
		}
	}
	assert(readPaths.size() == chunksPerRead.size());
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			assert(chunkBelongsToUnitig[chunk] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, size_t>> chunkAlleleReplacement;
	size_t nextNum = 0;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> unitigAllelesInThisRead;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (!canResolveUnitig[readPaths[i][j].path[k] & maskUint64_t]) continue;
				uint64_t prev = std::numeric_limits<size_t>::max();
				uint64_t after = std::numeric_limits<size_t>::max();
				if (k > 0)
				{
					prev = readPaths[i][j].path[k-1];
				}
				else if (j > 0)
				{
					prev = readPaths[i][j-1].path.back();
				}
				if (k+1 < readPaths[i][j].path.size())
				{
					after = readPaths[i][j].path[k+1];
				}
				else if (j+1 < readPaths[i].size())
				{
					after = readPaths[i][j+1].path[0];
				}
				if ((readPaths[i][j].path[k] ^ firstBitUint64_t) & firstBitUint64_t)
				{
					std::swap(prev, after);
					if (prev != std::numeric_limits<size_t>::max()) prev ^= firstBitUint64_t;
					if (after != std::numeric_limits<size_t>::max()) after ^= firstBitUint64_t;
				}
				if (prev == std::numeric_limits<size_t>::max() && after == std::numeric_limits<size_t>::max())
				{
					assert(readsAllowedToDiscard.count(i) == 1);
					assert(readPaths[i].size() == 1);
					assert(readPaths[i][0].path.size() == 1);
					continue;
				}
				assert(prev != std::numeric_limits<size_t>::max() || after != std::numeric_limits<size_t>::max());
				size_t alleleBefore = std::numeric_limits<size_t>::max();
				size_t alleleAfter = std::numeric_limits<size_t>::max();
				if (prev != std::numeric_limits<size_t>::max() && prevToAllele[readPaths[i][j].path[k] & maskUint64_t].count(prev) == 1) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(prev);
				if (after != std::numeric_limits<size_t>::max() && afterToAllele[readPaths[i][j].path[k] & maskUint64_t].count(after) == 1) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(after);
				assert(alleleBefore != std::numeric_limits<size_t>::max() || alleleAfter != std::numeric_limits<size_t>::max());
				assert(alleleBefore == alleleAfter || alleleBefore == std::numeric_limits<size_t>::max() || alleleAfter == std::numeric_limits<size_t>::max());
				size_t allele = alleleBefore;
				if (alleleBefore == std::numeric_limits<size_t>::max()) allele = alleleAfter;
				unitigAllelesInThisRead.emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k] & maskUint64_t, allele);
			}
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			uint64_t node = std::get<2>(chunksPerRead[i][j]);
			if (NonexistantChunk(node)) continue;
			size_t allele = std::numeric_limits<size_t>::max();
			if ((node & maskUint64_t) < chunkBelongsToUnitig.size() && chunkBelongsToUnitig[node & maskUint64_t] != std::numeric_limits<size_t>::max())
			{
				size_t unitig = chunkBelongsToUnitig[node & maskUint64_t];
				if (canResolveUnitig[unitig])
				{
					if (readsAllowedToDiscard.count(i) == 1)
					{
						assert(readPaths[i].size() == 1);
						assert(readPaths[i][0].path.size() == 1);
						std::get<2>(chunksPerRead[i][j]) = std::numeric_limits<size_t>::max();
						continue;
					}
					for (auto t : unitigAllelesInThisRead)
					{
						if (std::get<2>(t) != unitig) continue;
						if (std::get<0>(t) > std::get<0>(chunksPerRead[i][j])) continue;
						if (std::get<1>(t) < std::get<1>(chunksPerRead[i][j])) continue;
						assert(allele == std::numeric_limits<size_t>::max());
						allele = std::get<3>(t);
					}
					assert(allele != std::numeric_limits<size_t>::max());
				}
			}
			if (chunkAlleleReplacement.count(node & maskUint64_t) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t];
			}
			if (chunkAlleleReplacement.at(node & maskUint64_t).count(allele) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t][allele] = nextNum;
				nextNum += 1;
			}
			uint64_t replacement = chunkAlleleReplacement.at(node & maskUint64_t).at(allele);
			assert((replacement & firstBitUint64_t) == 0);
			replacement += (node & firstBitUint64_t);
			std::get<2>(chunksPerRead[i][j]) = replacement;
		}
	}
	size_t countUnitigsReplaced = 0;
	size_t countChunksReplaced = 0;
	for (size_t i = 0; i < canResolveUnitig.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		countUnitigsReplaced += 1;
		countChunksReplaced += graph.chunksInUnitig[i].size();
	}
	std::cerr << "semiambiguously resolvable unitig resolution resolved " << countUnitigsReplaced << " unitigs, " << countChunksReplaced << " chunks" << std::endl;
}

void resolveBetweenLongUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t longNodeThreshold, const size_t maxDistance, const size_t kmerSize)
{
	const uint64_t noNodeWithinDistance = (std::numeric_limits<size_t>::max() ^ (firstBitUint64_t>>1)) - 1;
	std::cerr << "resolving between long unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	phmap::flat_hash_set<size_t> chunkInsideLongUnitig;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < longNodeThreshold) continue;
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			chunkInsideLongUnitig.emplace(node & maskUint64_t);
		}
	}
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> longNodesInReads;
	longNodesInReads.resize(readPaths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (graph.unitigLengths[readPaths[i][j].path[k] & maskUint64_t] < longNodeThreshold) continue;
				longNodesInReads[i].emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k]);
			}
		}
		for (size_t j = 1; j < longNodesInReads[i].size(); j++)
		{
			assert(std::get<0>(longNodesInReads[i][j]) > std::get<0>(longNodesInReads[i][j-1]));
			assert(std::get<1>(longNodesInReads[i][j]) > std::get<1>(longNodesInReads[i][j-1]));
		}
	}
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	size_t countSplitted = 0;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&chunksPerRead, &occurrencesPerChunk, &longNodesInReads, &rawReadLengths, &nextNum, &resultMutex, &chunkInsideLongUnitig, &countSplitted, maxDistance, noNodeWithinDistance](const size_t chunki)
	{
		if (chunkInsideLongUnitig.count(chunki) == 1 || occurrencesPerChunk[chunki].size() < 4)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		std::vector<std::tuple<uint64_t, uint64_t, size_t, size_t>> predecessorsAndSuccessors;
		for (auto t : occurrencesPerChunk[chunki])
		{
			size_t read = t.first;
			size_t offset = t.second;
			bool fw = std::get<2>(chunksPerRead[read][offset]) & firstBitUint64_t;
			size_t chunkstartpos = std::get<0>(chunksPerRead[read][offset]);
			size_t chunkendpos = std::get<1>(chunksPerRead[read][offset]);
			uint64_t predecessor = std::numeric_limits<size_t>::max();
			uint64_t successor = std::numeric_limits<size_t>::max();
			for (size_t j = 0; j < longNodesInReads[read].size(); j++)
			{
				if (std::get<1>(longNodesInReads[read][j]) < chunkendpos && std::get<1>(longNodesInReads[read][j])+maxDistance > chunkstartpos)
				{
					predecessor = std::get<2>(longNodesInReads[read][j]);
				}
				if (std::get<0>(longNodesInReads[read][j]) > chunkstartpos && std::get<0>(longNodesInReads[read][j]) < chunkendpos+maxDistance && successor == std::numeric_limits<size_t>::max())
				{
					successor = std::get<2>(longNodesInReads[read][j]);
				}
			}
			if (chunkstartpos > maxDistance && predecessor == std::numeric_limits<size_t>::max()) predecessor = noNodeWithinDistance;
			if (chunkendpos + maxDistance < rawReadLengths[read] && successor == std::numeric_limits<size_t>::max()) successor = noNodeWithinDistance;
			if (!fw)
			{
				std::swap(predecessor, successor);
				if (predecessor != std::numeric_limits<size_t>::max() && predecessor != noNodeWithinDistance) predecessor ^= firstBitUint64_t;
				if (successor != std::numeric_limits<size_t>::max() && successor != noNodeWithinDistance) successor ^= firstBitUint64_t;
			}
			predecessorsAndSuccessors.emplace_back(predecessor, successor, read, offset);
		}
		assert(predecessorsAndSuccessors.size() == occurrencesPerChunk[chunki].size());
		bool hasAny = false;
		phmap::flat_hash_set<uint64_t> mergedSomewhere;
		phmap::flat_hash_map<uint64_t, uint64_t> parent;
		for (auto t : predecessorsAndSuccessors)
		{
			assert(std::get<1>(t) == std::numeric_limits<size_t>::max() || ((std::get<1>(t) & (firstBitUint64_t>>1)) == 0));
			if (std::get<0>(t) != std::numeric_limits<size_t>::max() && parent.count(std::get<0>(t)) == 0)
			{
				hasAny = true;
				parent[std::get<0>(t)] = std::get<0>(t);
			}
			if (std::get<1>(t) != std::numeric_limits<size_t>::max() && parent.count(std::get<1>(t) + (firstBitUint64_t>>1)) == 0)
			{
				hasAny = true;
				parent[std::get<1>(t) + (firstBitUint64_t>>1)] = std::get<1>(t) + (firstBitUint64_t>>1);
			}
			if (std::get<0>(t) != std::numeric_limits<size_t>::max() && std::get<1>(t) != std::numeric_limits<size_t>::max())
			{
				mergedSomewhere.insert(std::get<0>(t));
				mergedSomewhere.insert(std::get<1>(t) + (firstBitUint64_t>>1));
				merge(parent, std::get<0>(t), std::get<1>(t) + (firstBitUint64_t>>1));
			}
		}
		if (!hasAny)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "chunk " << chunki << " doesn't have any" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		bool canSplit = true;
		for (auto pair : parent)
		{
			if (mergedSomewhere.count(pair.first) == 0)
			{
				canSplit = false;
				break;
			}
		}
		if (!canSplit)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "can't split chunk " << chunki << " because of missing merge" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		phmap::flat_hash_map<uint64_t, size_t> clusterToNode;
		size_t nextId = 0;
		for (uint64_t node : mergedSomewhere)
		{
			uint64_t found = find(parent, node);
			if (clusterToNode.count(found) == 1)
			{
				continue;
			}
			clusterToNode[found] = nextId;
			nextId += 1;
		}
		phmap::flat_hash_map<uint64_t, size_t> clusterCoverage;
		for (auto t : predecessorsAndSuccessors)
		{
			if (std::get<0>(t) != std::numeric_limits<size_t>::max())
			{
				clusterCoverage[clusterToNode.at(find(parent, std::get<0>(t)))] += 1;
			}
			else if (std::get<1>(t) != std::numeric_limits<size_t>::max())
			{
				clusterCoverage[clusterToNode.at(find(parent, std::get<1>(t) + (firstBitUint64_t>>1)))] += 1;
			}
		}
		for (auto pair : clusterCoverage)
		{
			if (pair.second < 2)
			{
				canSplit = false;
			}
		}
		if (!canSplit)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "can't split chunk " << chunki << " because of low coverage" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		assert(clusterToNode.size() >= 1);
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			if (clusterToNode.size() >= 2) countSplitted += 1;
//			std::cerr << "chunk " << chunki << " split to " << clusterToNode.size() << " chunks" << std::endl;
			for (auto t : predecessorsAndSuccessors)
			{
				if (std::get<0>(t) == std::numeric_limits<size_t>::max() && std::get<1>(t) == std::numeric_limits<size_t>::max())
				{
					std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) = std::numeric_limits<size_t>::max();
					continue;
				}
				size_t cluster;
				if (std::get<0>(t) != std::numeric_limits<size_t>::max())
				{
					cluster = clusterToNode.at(find(parent, std::get<0>(t)));
				}
				else
				{
					assert(std::get<1>(t) != std::numeric_limits<size_t>::max());
					cluster = clusterToNode.at(find(parent, std::get<1>(t) + (firstBitUint64_t>>1)));
				}
				std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) = nextNum + cluster + (std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) & firstBitUint64_t);
			}
			nextNum += clusterToNode.size();
			return;
		}
	});
	std::cerr << "between long unitigs resolved " << countSplitted << " chunks" << std::endl;
}

void resolveUnambiguouslyResolvableUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving unambiguously resolvable unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<bool> edgesBalanced;
	std::vector<bool> unitigHasContainedPath;
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> tripletsPerUnitig;
	edgesBalanced.resize(graph.unitigLengths.size(), false);
	unitigHasContainedPath.resize(graph.unitigLengths.size(), false);
	tripletsPerUnitig.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.edges.getEdges(std::make_pair(i, true)).size() == graph.edges.getEdges(std::make_pair(i, false)).size())
		{
			if (graph.edges.getEdges(std::make_pair(i, true)).size() >= 2)
			{
				edgesBalanced[i] = true;
			}
		}
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			if (readPaths[i][j].path.size() == 1)
			{
				unitigHasContainedPath[readPaths[i][j].path[0] & maskUint64_t] = true;
				continue;
			}
			assert(readPaths[i][j].path.size() >= 2);
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				if (!edgesBalanced[mid & maskUint64_t]) continue;
				if (unitigHasContainedPath[mid & maskUint64_t]) continue;
				uint64_t prev, after;
				if (mid & firstBitUint64_t)
				{
					prev = readPaths[i][j].path[k-1];
					after = readPaths[i][j].path[k+1];
				}
				else
				{
					prev = readPaths[i][j].path[k+1] ^ firstBitUint64_t;
					after = readPaths[i][j].path[k-1] ^ firstBitUint64_t;
				}
				tripletsPerUnitig[mid & maskUint64_t][std::make_pair(prev, after)] += 1;
			}
		}
	}
	std::vector<bool> canResolveUnitig;
	canResolveUnitig.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (!edgesBalanced[i]) continue;
		if (unitigHasContainedPath[i]) continue;
		if (tripletsPerUnitig[i].size() != graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		phmap::flat_hash_set<uint64_t> nodesPrev;
		phmap::flat_hash_set<uint64_t> nodesAfter;
		bool valid = true;
		for (auto triplet : tripletsPerUnitig[i])
		{
			if (triplet.second < 2)
			{
				valid = false;
				break;
			}
			nodesPrev.insert(triplet.first.first);
			nodesAfter.insert(triplet.first.second);
		}
		if (!valid) continue;
		if (nodesPrev.size() != nodesAfter.size()) continue;
		if (nodesPrev.size() != graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		canResolveUnitig[i] = true;
	}
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> prevToAllele;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> afterToAllele;
	prevToAllele.resize(graph.unitigLengths.size());
	afterToAllele.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		size_t alleleNum = 0;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t prev = triplet.first.first;
			uint64_t after = triplet.first.second;
			assert(prevToAllele[i].count(prev) == 0);
			assert(afterToAllele[i].count(after) == 0);
			prevToAllele[i][prev] = alleleNum;
			afterToAllele[i][after] = alleleNum;
			alleleNum += 1;
		}
	}
	assert(readPaths.size() == chunksPerRead.size());
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			assert(chunkBelongsToUnitig[chunk] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, size_t>> chunkAlleleReplacement;
	size_t nextNum = 0;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> unitigAllelesInThisRead;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (!canResolveUnitig[readPaths[i][j].path[k] & maskUint64_t]) continue;
				size_t alleleBefore = std::numeric_limits<size_t>::max();
				size_t alleleAfter = std::numeric_limits<size_t>::max();
				if (readPaths[i][j].path[k] & firstBitUint64_t)
				{
					if (k > 0) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k-1]);
					if (k+1 < readPaths[i][j].path.size()) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k+1]);
				}
				else
				{
					if (k > 0) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k-1] ^ firstBitUint64_t);
					if (k+1 < readPaths[i][j].path.size()) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k+1] ^ firstBitUint64_t);
				}
				assert(alleleBefore != std::numeric_limits<size_t>::max() || alleleAfter != std::numeric_limits<size_t>::max());
				assert(alleleBefore == alleleAfter || alleleBefore == std::numeric_limits<size_t>::max() || alleleAfter == std::numeric_limits<size_t>::max());
				size_t allele = alleleBefore;
				if (alleleBefore == std::numeric_limits<size_t>::max()) allele = alleleAfter;
				unitigAllelesInThisRead.emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k] & maskUint64_t, allele);
			}
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			uint64_t node = std::get<2>(chunksPerRead[i][j]);
			if (NonexistantChunk(node)) continue;
			size_t allele = std::numeric_limits<size_t>::max();
			if ((node & maskUint64_t) < chunkBelongsToUnitig.size() && chunkBelongsToUnitig[node & maskUint64_t] != std::numeric_limits<size_t>::max())
			{
				size_t unitig = chunkBelongsToUnitig[node & maskUint64_t];
				if (canResolveUnitig[unitig])
				{
					for (auto t : unitigAllelesInThisRead)
					{
						if (std::get<2>(t) != unitig) continue;
						if (std::get<0>(t) > std::get<0>(chunksPerRead[i][j])) continue;
						if (std::get<1>(t) < std::get<1>(chunksPerRead[i][j])) continue;
						assert(allele == std::numeric_limits<size_t>::max());
						allele = std::get<3>(t);
					}
					assert(allele != std::numeric_limits<size_t>::max());
				}
			}
			if (chunkAlleleReplacement.count(node & maskUint64_t) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t];
			}
			if (chunkAlleleReplacement.at(node & maskUint64_t).count(allele) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t][allele] = nextNum;
				nextNum += 1;
			}
			uint64_t replacement = chunkAlleleReplacement.at(node & maskUint64_t).at(allele);
			assert((replacement & firstBitUint64_t) == 0);
			replacement += (node & firstBitUint64_t);
			std::get<2>(chunksPerRead[i][j]) = replacement;
		}
	}
	size_t countUnitigsReplaced = 0;
	size_t countChunksReplaced = 0;
	for (size_t i = 0; i < canResolveUnitig.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		countUnitigsReplaced += 1;
		countChunksReplaced += graph.chunksInUnitig[i].size();
	}
	std::cerr << "unambiguously resolvable unitig resolution resolved " << countUnitigsReplaced << " unitigs, " << countChunksReplaced << " chunks" << std::endl;
}

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

std::vector<std::pair<size_t, size_t>> getNoncontainedUniqueRanges(std::vector<std::pair<size_t, size_t>> raw)
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
	size_t longestExtent = raw[0].second;
	std::vector<bool> contained;
	contained.resize(raw.size(), false);
	for (size_t j = 1; j < raw.size(); j++)
	{
		assert(raw[j].first >= raw[j-1].first);
		if (raw[j].second <= longestExtent) contained[j] = true;
		if (raw[j].first == raw[j-1].first)
		{
			assert(raw[j].second > raw[j-1].second);
			contained[j-1] = true;
		}
		longestExtent = std::max(longestExtent, raw[j].second);
	}
	for (size_t j = raw.size()-1; j < raw.size(); j--)
	{
		if (!contained[j]) continue;
		std::swap(raw[j], raw.back());
		raw.pop_back();
	}
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
			size_t end = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			while (end < chunksPerRead[i].size() && length < expandedSize)
			{
				end += 1;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end]))) break;
				if (end == chunksPerRead[i].size()) break;
				length += chunkLengths[std::get<2>(chunksPerRead[i][end]) & maskUint64_t];
				length -= kmerSize;
			}
			if (length < expandedSize) continue;
			if (end == chunksPerRead[i].size()) continue;
			newChunkmerRanges.emplace_back(j, end);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t start = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
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
		newChunkmerRanges = getNoncontainedUniqueRanges(newChunkmerRanges);
		assert(newChunkmerRanges.size() >= 1);
		for (size_t j = 1; j < newChunkmerRanges.size(); j++)
		{
			assert(newChunkmerRanges[j].first > newChunkmerRanges[j-1].first);
			assert(newChunkmerRanges[j].second > newChunkmerRanges[j-1].second);
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

void resolveTinyNodesRecklessly(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize)
{
	size_t nextNum = 0;
	std::mutex resultMutex;
	size_t countResolved = 0;
	size_t countResolvedTo = 0;
	std::vector<std::tuple<size_t, size_t, size_t>> replacements;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &countResolved, &countResolvedTo, &replacements, kmerSize](const size_t chunki, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		if (occurrencesPerChunk[chunki].size() <= 2)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				replacements.emplace_back(occurrencesPerChunk[chunki][j].first, occurrencesPerChunk[chunki][j].second, nextNum);
			}
			nextNum += 1;
			return;
		}
		bool isSmall = true;
		for (auto pair : occurrencesPerChunk[chunki])
		{
			auto t = chunksPerRead[pair.first][pair.second];
			if (std::get<1>(t) - std::get<0>(t) >= kmerSize*2)
			{
				isSmall = false;
				break;
			}
		}
		if (!isSmall)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				replacements.emplace_back(occurrencesPerChunk[chunki][j].first, occurrencesPerChunk[chunki][j].second, nextNum);
			}
			nextNum += 1;
			return;
		}
		phmap::flat_hash_map<size_t, size_t> parent;
		for (auto pair : occurrencesPerChunk[chunki])
		{
			size_t prev = std::numeric_limits<size_t>::max();
			size_t next = std::numeric_limits<size_t>::max();
			if (pair.second > 0 && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second-1])))
			{
				prev = std::get<2>(chunksPerRead[pair.first][pair.second-1]);
				if (parent.count(prev) == 0) parent[prev] = prev;
			}
			if (pair.second+1 < chunksPerRead[pair.first].size() && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second+1])))
			{
				next = std::get<2>(chunksPerRead[pair.first][pair.second+1]) + (firstBitUint64_t >> 1);
				if (parent.count(next) == 0) parent[next] = next;
			}
			if (prev != std::numeric_limits<size_t>::max() && next != std::numeric_limits<size_t>::max())
			{
				merge(parent, prev, next);
			}
		}
		phmap::flat_hash_map<size_t, size_t> clusterToNode;
		size_t nextClusterNum = 1;
		for (auto pair : parent)
		{
			if (pair.second != pair.first) continue;
			clusterToNode[pair.second] = nextClusterNum;
			nextClusterNum += 1;
		}
		size_t zeroOffset = 0;
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			zeroOffset = nextNum;
			nextNum += nextClusterNum;
			if (nextClusterNum > 2)
			{
				countResolved += 1;
				countResolvedTo += nextClusterNum-1;
			}
			for (auto pair : occurrencesPerChunk[chunki])
			{
				size_t clusterNum = 0;
				if (pair.second > 0 && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second-1])))
				{
					clusterNum = clusterToNode.at(find(parent, std::get<2>(chunksPerRead[pair.first][pair.second-1])));
				}
				if (pair.second+1 < chunksPerRead[pair.first].size() && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second+1])))
				{
					size_t foundClusterNum = clusterToNode.at(find(parent, std::get<2>(chunksPerRead[pair.first][pair.second+1]) + (firstBitUint64_t >> 1)));
					assert(clusterNum == 0 || foundClusterNum == clusterNum);
					clusterNum = foundClusterNum;
				}
				replacements.emplace_back(pair.first, pair.second, zeroOffset + clusterNum);
			}
		}
	});
	for (auto t : replacements)
	{
		std::get<2>(chunksPerRead[std::get<0>(t)][std::get<1>(t)]) = std::get<2>(t) + (std::get<2>(chunksPerRead[std::get<0>(t)][std::get<1>(t)]) & firstBitUint64_t);
	}
	std::cerr << "tiny node resolution resolved " << countResolved << " chunks to " << countResolvedTo << " chunks" << std::endl;
}

void removeTinyProblemNodes(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			assert(std::get<2>(t) & firstBitUint64_t);
			maxChunk = std::max(maxChunk, std::get<2>(t) & maskUint64_t);
		}
	}
	std::vector<size_t> uniquePredecessor;
	std::vector<size_t> uniqueSuccessor;
	uniquePredecessor.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	uniqueSuccessor.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1]))) continue;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t prev = std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t;
			size_t curr = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			assert(prev < uniqueSuccessor.size());
			assert(curr < uniquePredecessor.size());
			if (uniqueSuccessor[prev] == std::numeric_limits<size_t>::max())
			{
				uniqueSuccessor[prev] = curr;
			}
			else if (uniqueSuccessor[prev] != curr)
			{
				uniqueSuccessor[prev] = std::numeric_limits<size_t>::max()-1;
			}
			if (uniquePredecessor[curr] == std::numeric_limits<size_t>::max())
			{
				uniquePredecessor[curr] = prev;
			}
			else if (uniquePredecessor[curr] != prev)
			{
				uniquePredecessor[curr] = std::numeric_limits<size_t>::max()-1;
			}
		}
	}
	std::vector<bool> canBeRemoved;
	canBeRemoved.resize(maxChunk+1, true);
	for (size_t i = 0; i < uniquePredecessor.size(); i++)
	{
		if (uniquePredecessor[i] != std::numeric_limits<size_t>::max()-1 && uniqueSuccessor[i] != std::numeric_limits<size_t>::max()-1)
		{
			canBeRemoved[i] = false;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i][0])) && std::get<1>(chunksPerRead[i][0]) - std::get<0>(chunksPerRead[i][0]) >= 2*kmerSize) canBeRemoved[std::get<2>(chunksPerRead[i][0]) & maskUint64_t] = false;
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i].back())) && std::get<1>(chunksPerRead[i].back()) - std::get<0>(chunksPerRead[i].back()) >= 2*kmerSize) canBeRemoved[std::get<2>(chunksPerRead[i].back()) & maskUint64_t] = false;
		for (size_t j = 1; j+1 < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			if (std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) >= 2*kmerSize)
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j])] = false;
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j+1])))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j])] = false;
				continue;
			}
			if (canBeRemoved[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (canBeRemoved[std::get<2>(chunksPerRead[i][j+1]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (std::get<1>(chunksPerRead[i][j-1])+1 <= std::get<0>(chunksPerRead[i][j+1]))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
		}
	}
	for (size_t i = 0; i < chunksPerRead.size()/2; i++)
	{
		size_t other = i + chunksPerRead.size()/2;
		assert(chunksPerRead[i].size() == chunksPerRead[other].size());
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) - std::get<0>(chunksPerRead[other][chunksPerRead[other].size()-1-j]));
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				assert(NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])));
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])))
			{
				std::cerr << i << " " << other << std::endl;
				std::cerr << j << " " << chunksPerRead[i].size() << std::endl;
				std::cerr << std::get<2>(chunksPerRead[i][j]) << " " << ((std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(chunksPerRead[i][j]) & maskUint64_t) << std::endl;
				std::cerr << std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) << " " << ((std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t) << std::endl;
			}
			assert(!NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])));
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t] = false;
			}
			if (!canBeRemoved[std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
			}
		}
	}
	size_t countRemovableChunks = 0;
	for (size_t i = 0; i < canBeRemoved.size(); i++)
	{
		if (!canBeRemoved[i]) continue;
		countRemovableChunks += 1;
	}
	size_t countRemoved = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = chunksPerRead[i].size()-1; j < chunksPerRead[i].size(); j--)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			std::swap(chunksPerRead[i][j], chunksPerRead[i].back());
			chunksPerRead[i].pop_back();
			countRemoved += 1;
		}
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<0>(chunksPerRead[i][j]) > std::get<0>(chunksPerRead[i][j-1]));
			assert(std::get<1>(chunksPerRead[i][j]) > std::get<1>(chunksPerRead[i][j-1]));
		}
	}
	std::cerr << countRemovableChunks << " distinct removable chunks, removed " << countRemoved << " total chunks" << std::endl;
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
			expandChunks(chunksPerRead, kmerSize, 5000);
			writeStage(51, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 51:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(6, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 6:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(7, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 7:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerNearestNeighborPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(8, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 8:
			splitPerCorrectedKmerPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(9, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 9:
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
			writeReadChunkSequences("sequences-chunk19.txt", rawReadLengths, chunksPerRead, sequenceIndex);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(19, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
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
