#include <queue>
#include <iostream>
#include "ChunkResolution.h"
#include "Common.h"
#include "ChunkHelper.h"
#include "UnionFind.h"
#include "ChunkUnitigGraph.h"
#include "SequenceHelper.h"
#include "KmerIterator.h"
#include "CanonHelper.h"
#include "TransitiveClosure.h"

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

void splitPerLength(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by length" << std::endl;
	const double differenceFraction = 0.02;
	const size_t differenceConstant = 50;
	size_t nextNum = 0;
	std::mutex resultMutex;
	auto oldChunks = chunksPerRead;
	iterateCanonicalChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, differenceFraction, differenceConstant, kmerSize](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		if (occurrencesPerChunk[i].size() == 0) return;
		std::vector<std::pair<size_t, size_t>> lengthsAndIndicesHere;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			lengthsAndIndicesHere.emplace_back(std::get<1>(t) - std::get<0>(t), j);
		}
		std::sort(lengthsAndIndicesHere.begin(), lengthsAndIndicesHere.end());
		size_t clusterStart = 0;
		for (size_t j = 1; j <= lengthsAndIndicesHere.size(); j++)
		{
			if (j < lengthsAndIndicesHere.size())
			{
				if (lengthsAndIndicesHere[j].first - lengthsAndIndicesHere[j-1].first < std::max((size_t)(lengthsAndIndicesHere[j-1].first * differenceFraction), (size_t)differenceConstant))
				{
					if (lengthsAndIndicesHere[j-1].first > 2*kmerSize)
					{
						continue;
					}
					else if (lengthsAndIndicesHere[j].first - lengthsAndIndicesHere[j-1].first < 2)
					{
						continue;
					}
				}
			}
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t k = clusterStart; k < j; k++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][lengthsAndIndicesHere[k].second].first][occurrencesPerChunk[i][lengthsAndIndicesHere[k].second].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][lengthsAndIndicesHere[k].second].first][occurrencesPerChunk[i][lengthsAndIndicesHere[k].second].second]) & firstBitUint64_t) + nextNum;
			}
			nextNum += 1;
			clusterStart = j;
		}
	});
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
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

void splitPerFirstLastKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	assert(kmerSize <= 31);
	size_t nextNum = 0;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> endKmersToNumber;
	std::mutex resultMutex;
	iterateMultithreaded(0, chunksPerRead.size(), numThreads, [&sequenceIndex, &nextNum, &resultMutex, &endKmersToNumber, &chunksPerRead, kmerSize](const size_t i)
	{
		if (chunksPerRead[i].size() == 0) return;
		std::string readseq = getReadSequence(sequenceIndex, i);
		std::vector<std::pair<size_t, size_t>> keys;
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
			keys.emplace_back(firstKmer, lastKmer);
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				auto key = keys[j];
				if (endKmersToNumber.count(key) == 0)
				{
					endKmersToNumber[key] = nextNum;
					nextNum += 1;
				}
				std::get<2>(chunksPerRead[i][j]) = endKmersToNumber.at(key) + firstBitUint64_t;
			}
		}
	});
}

void splitPerBaseCounts(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	const size_t mismatchFloor = 10;
	std::cerr << "splitting by base counts" << std::endl;
	size_t nextNum = 0;
	std::mutex resultMutex;
	auto oldChunks = chunksPerRead;
	size_t countSplitted = 0;
	size_t countSplittedTo = 0;
	iterateCanonicalChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, &countSplitted, &countSplittedTo, kmerSize, mismatchFloor](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		if (occurrencesPerChunk[i].size() < 10)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + nextNum;
			}
			nextNum += 1;
			return;
		}
		bool allSmall = true;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			if (std::get<1>(chunksPerRead[occurrencesPerChunk[i][0].first][occurrencesPerChunk[i][0].second]) - std::get<0>(chunksPerRead[occurrencesPerChunk[i][0].first][occurrencesPerChunk[i][0].second]) < 2*kmerSize) continue;
			allSmall = false;
			break;
		}
		if (allSmall)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + nextNum;
			}
			nextNum += 1;
			return;
		}
		std::vector<std::vector<size_t>> countsPerOccurrence;
		countsPerOccurrence.resize(occurrencesPerChunk[i].size());
		size_t maxDistance = 0;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			std::string chunkSequence = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
			maxDistance = std::max(maxDistance, chunkSequence.size());
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
		maxDistance *= mismatchFraction;
		maxDistance = std::max(maxDistance, mismatchFloor);
		std::vector<size_t> parent = getFastTransitiveClosure(countsPerOccurrence.size(), maxDistance, [&countsPerOccurrence](const size_t left, const size_t right, const size_t maxDist)
		{
			size_t result = 0;
			for (size_t c = 0; c < 4; c++)
			{
				if (countsPerOccurrence[left][c] > countsPerOccurrence[right][c])
				{
					result += countsPerOccurrence[left][c] - countsPerOccurrence[right][c];
				}
				else
				{
					result += countsPerOccurrence[right][c] - countsPerOccurrence[left][c];
				}
			}
			return result;
		});
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
			if (keyToNode.size() >= 2)
			{
				countSplitted += 1;
				countSplittedTo += keyToNode.size();
			}
//			std::cerr << "base count splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
	std::cerr << "base count splitted " << countSplitted << " chunks to " << countSplittedTo << " chunks" << std::endl;
}

void splitPerMinHashes(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t kmerSize = 11;
	const size_t hashCount = 10;
	std::cerr << "splitting by minhash" << std::endl;
	size_t nextNum = 0;
	std::mutex resultMutex;
	auto oldChunks = chunksPerRead;
	iterateCanonicalChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, kmerSize, hashCount](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
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
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
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
				unitigAllelesInThisRead.emplace_back(readPaths[i][j].chunksInPathnode[k].first, readPaths[i][j].chunksInPathnode[k].second, readPaths[i][j].path[k] & maskUint64_t, allele);
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
						if (std::get<0>(t) > j) continue;
						if (std::get<1>(t) < j) continue;
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
				unitigAllelesInThisRead.emplace_back(readPaths[i][j].chunksInPathnode[k].first, readPaths[i][j].chunksInPathnode[k].second, readPaths[i][j].path[k] & maskUint64_t, allele);
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
						if (std::get<0>(t) > j) continue;
						if (std::get<1>(t) < j) continue;
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

void resolveTinyNodesRecklessly(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize)
{
	size_t nextNum = 0;
	std::mutex resultMutex;
	size_t countResolved = 0;
	size_t countResolvedTo = 0;
	std::vector<std::tuple<size_t, size_t, size_t>> replacements;
	auto oldChunks = chunksPerRead;
	iterateCanonicalChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &countResolved, &countResolvedTo, &replacements, kmerSize](const size_t chunki, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
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
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
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

void resolveBetweenTangles(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t longUnitigThreshold, const size_t kmerSize)
{
	std::cerr << "resolving between tangles" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<bool> unitigIsLong;
	unitigIsLong.resize(graph.unitigLengths.size(), false);
	size_t countLongUnitigs = 0;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] >= longUnitigThreshold)
		{
			unitigIsLong[i] = true;
			countLongUnitigs += 1;
		}
	}
	std::cerr << "unitigs " << graph.unitigLengths.size() << " of which long " << countLongUnitigs << std::endl;
	if (countLongUnitigs < 2) return;
	std::vector<size_t> parent;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		parent.emplace_back(i*2);
		parent.emplace_back(i*2+1);
	}
	for (size_t i = 0; i < unitigIsLong.size(); i++)
	{
		if (unitigIsLong[i]) continue;
		merge(parent, i*2, i*2+1);
	}
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		for (auto edge : graph.edges.getEdges(std::make_pair(i, true)))
		{
			merge(parent, i*2+1, edge.first*2 + (edge.second ? 0 : 1));
		}
		for (auto edge : graph.edges.getEdges(std::make_pair(i, false)))
		{
			merge(parent, i*2, edge.first*2 + (edge.second ? 0 : 1));
		}
	}
	phmap::flat_hash_set<size_t> tangles;
	for (size_t i = 0; i < parent.size(); i++)
	{
		tangles.emplace(find(parent, i));
	}
	for (size_t i = 0; i < parent.size(); i++)
	{
		std::cerr << "tip " << ((i%2) ? ">" : "<") << (i/2) << " in tangle " << find(parent, i) << std::endl;
	}
	phmap::flat_hash_map<size_t, std::vector<uint64_t>> longUnitigsPerTangle;
	for (size_t i = 0; i < unitigIsLong.size(); i++)
	{
		if (!unitigIsLong[i]) continue;
		longUnitigsPerTangle[find(parent, i*2+1)].emplace_back(i + firstBitUint64_t);
		longUnitigsPerTangle[find(parent, i*2)].emplace_back(i);
	}
	for (const auto& pair : longUnitigsPerTangle)
	{
		if (pair.second.size() < 4)
		{
			assert(tangles.count(pair.first) == 1);
			tangles.erase(pair.first);
		}
	}
	std::cerr << tangles.size() << " tangles after filtering for long unitig touchers" << std::endl;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, size_t>> connectionCounts;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		uint64_t lastLongUnitig = std::numeric_limits<uint64_t>::max();
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				size_t unitig = readPaths[i][j].path[k] & maskUint64_t;
				assert(unitig < unitigIsLong.size());
				if (!unitigIsLong[unitig]) continue;
				uint64_t longUnitigHere = readPaths[i][j].path[k];
				if (lastLongUnitig == std::numeric_limits<size_t>::max())
				{
					lastLongUnitig = longUnitigHere;
					continue;
				}
				if (find(parent, 2 * (lastLongUnitig & maskUint64_t) + ((lastLongUnitig & firstBitUint64_t) ? 1 : 0)) != find(parent, 2 * (longUnitigHere & maskUint64_t) + ((longUnitigHere & firstBitUint64_t) ? 0 : 1)))
				{
					lastLongUnitig = longUnitigHere;
					continue;
				}
				connectionCounts[lastLongUnitig][longUnitigHere ^ firstBitUint64_t] += 1;
				connectionCounts[longUnitigHere ^ firstBitUint64_t][lastLongUnitig] += 1;
				lastLongUnitig = longUnitigHere;
			}
		}
	}
	phmap::flat_hash_set<size_t> forbiddenTangles;
	for (size_t i = 0; i < unitigIsLong.size(); i++)
	{
		if (!unitigIsLong[i]) continue;
		if (connectionCounts.count(i) == 0)
		{
			forbiddenTangles.emplace(find(parent, 2*i));
		}
		else
		{
			bool good = true;
			size_t minCoverage = 4;
			if (connectionCounts.at(i).size() == 1)
			{
				minCoverage = 2;
			}
			for (auto pair : connectionCounts.at(i))
			{
				if (pair.second < minCoverage) good = false;
			}
			if (!good) forbiddenTangles.emplace(find(parent, 2*i));
		}
		if (connectionCounts.count(i+firstBitUint64_t) == 0)
		{
			forbiddenTangles.emplace(find(parent, 2*i+1));
		}
		else
		{
			bool good = true;
			size_t minCoverage = 4;
			if (connectionCounts.at(i+firstBitUint64_t).size() == 1)
			{
				minCoverage = 2;
			}
			for (auto pair : connectionCounts.at(i+firstBitUint64_t))
			{
				if (pair.second < minCoverage) good = false;
			}
			if (!good) forbiddenTangles.emplace(find(parent, 2*i+1));
		}
	}
	for (auto tangle : forbiddenTangles)
	{
		if (tangles.count(tangle) == 1) tangles.erase(tangle);
	}
	std::cerr << tangles.size() << " tangles after filtering for connections" << std::endl;
	phmap::flat_hash_map<uint64_t, uint64_t> longUnitigTipParent;
	for (const auto& pair : connectionCounts)
	{
		size_t tangle = find(parent, 2 * (pair.first & maskUint64_t) + ((pair.first & firstBitUint64_t) ? 1 : 0));
		if (tangles.count(tangle) == 0) continue;
		longUnitigTipParent[pair.first] = pair.first;
	}
	for (const auto& pair : connectionCounts)
	{
		size_t tangle = find(parent, 2 * (pair.first & maskUint64_t) + ((pair.first & firstBitUint64_t) ? 1 : 0));
		if (tangles.count(tangle) == 0) continue;
		assert(longUnitigTipParent.count(pair.first) == 1);
		for (auto pair2 : pair.second)
		{
			assert(longUnitigTipParent.count(pair2.first) == 1);
			merge(longUnitigTipParent, pair.first, pair2.first);
		}
	}
	phmap::flat_hash_set<size_t> tanglesWhichHaveShortNodes;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (unitigIsLong[i]) continue;
		size_t tangle = find(parent, 2*i);
		if (tangles.count(tangle) == 0) continue;
		tanglesWhichHaveShortNodes.insert(tangle);
	}
	phmap::flat_hash_set<size_t> tanglesForbiddenByLackOfShortNodes;
	for (size_t tangle : tangles)
	{
		if (tanglesWhichHaveShortNodes.count(tangle) == 1) continue;
		tanglesForbiddenByLackOfShortNodes.emplace(tangle);
	}
	for (size_t tangle : tanglesForbiddenByLackOfShortNodes)
	{
		assert(tangles.count(tangle) == 1);
		tangles.erase(tangle);
	}
	std::cerr << tangles.size() << " tangles after filtering for short nodes" << std::endl;
	phmap::flat_hash_set<size_t> tanglesForbiddenByLackOfResolutions;
	for (size_t tangle : tangles)
	{
		phmap::flat_hash_set<uint64_t> longUnitigTipClustersHere;
		assert(longUnitigsPerTangle.count(tangle) == 1);
		for (uint64_t longUnitigTip : longUnitigsPerTangle.at(tangle))
		{
			assert(longUnitigTipParent.count(longUnitigTip) == 1);
			longUnitigTipClustersHere.emplace(find(longUnitigTipParent, longUnitigTip));
		}
		if (longUnitigTipClustersHere.size() < 2) tanglesForbiddenByLackOfResolutions.insert(tangle);
	}
	for (size_t tangle : tanglesForbiddenByLackOfResolutions)
	{
		assert(tangles.count(tangle) == 1);
		tangles.erase(tangle);
	}
	std::cerr << tangles.size() << " tangles after filtering for resolutions" << std::endl;
	{
		phmap::flat_hash_map<size_t, std::vector<uint64_t>> nodesPerTangle;
		for (size_t i = 0; i < parent.size(); i++)
		{
			nodesPerTangle[find(parent, i)].emplace_back(i);
		}
		for (size_t tangle : tangles)
		{
			std::cerr << "resolve tangle " << tangle << " with tips:";
			for (uint64_t node : nodesPerTangle.at(tangle))
			{
				std::cerr << " " << (node % 2 ? ">" : "<") << (node/2);
			}
			std::cerr << std::endl;
		}
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, size_t>>>> readUnitigBordersPerTangle; // tangle -> read -> pathindex, nodeindex
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				size_t unitig = readPaths[i][j].path[k] & maskUint64_t;
				if (!unitigIsLong[unitig]) continue;
				size_t prevTangle = find(parent, 2*unitig);
				if (tangles.count(prevTangle) == 1) readUnitigBordersPerTangle[prevTangle][i].emplace_back(j, k);
				size_t nextTangle = find(parent, 2*unitig+1);
				if (tangles.count(nextTangle) == 1) readUnitigBordersPerTangle[nextTangle][i].emplace_back(j, k);
			}
		}
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, std::vector<std::tuple<size_t, size_t, size_t, size_t>>>> readPartsPerTangle; // tangle -> read -> startpathindex, startnodeindex, endpathindex, endnodeindex (INCLUSIVE)
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		size_t currentTangle = std::numeric_limits<size_t>::max();
		size_t currentStartPath = std::numeric_limits<size_t>::max();
		size_t currentStartNode = std::numeric_limits<size_t>::max();
		size_t currentEndPath = std::numeric_limits<size_t>::max();
		size_t currentEndNode = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				size_t unitig = readPaths[i][j].path[k] & maskUint64_t;
				size_t tangleHere = find(parent, 2*unitig);
				if (unitigIsLong[unitig] || tangles.count(tangleHere) == 0 || tangleHere != currentTangle)
				{
					if (currentTangle != std::numeric_limits<size_t>::max())
					{
						readPartsPerTangle[currentTangle][i].emplace_back(currentStartPath, currentStartNode, currentEndPath, currentEndNode);
					}
					currentStartPath = std::numeric_limits<size_t>::max();
					currentStartNode = std::numeric_limits<size_t>::max();
					currentTangle = std::numeric_limits<size_t>::max();
					if (unitigIsLong[unitig] || tangles.count(tangleHere) == 0)
					{
						continue;
					}
				}
				if (currentTangle == std::numeric_limits<size_t>::max())
				{
					currentTangle = tangleHere;
					currentStartPath = j;
					currentStartNode = k;
				}
				assert(currentTangle == tangleHere);
				currentEndPath = j;
				currentEndNode = k;
			}
		}
		if (currentTangle != std::numeric_limits<size_t>::max())
		{
			readPartsPerTangle[currentTangle][i].emplace_back(currentStartPath, currentStartNode, currentEndPath, currentEndNode);
		}
	}
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t>>> readTanglePartAssignments; // startpathindex, startnodeindex, endpathindex, endnodeindex, allele
	readTanglePartAssignments.resize(readPaths.size());
	for (size_t tangle : tangles)
	{
		assert(readUnitigBordersPerTangle.count(tangle) == 1);
		assert(readPartsPerTangle.count(tangle) == 1);
		phmap::flat_hash_map<size_t, uint64_t> innerNodesBelongToResolution;
		for (const auto& pair : readPartsPerTangle.at(tangle))
		{
			size_t read = pair.first;
			if (readUnitigBordersPerTangle.at(tangle).count(read) == 0) continue;
			for (const auto& partt : readPartsPerTangle.at(tangle).at(read))
			{
				size_t tangle = find(parent, 2*readPaths[read][std::get<0>(partt)].path[std::get<1>(partt)]);
				size_t alleleIndexBefore = std::numeric_limits<size_t>::max();
				size_t alleleIndexAfter = std::numeric_limits<size_t>::max();
				for (const auto& bordert : readUnitigBordersPerTangle.at(tangle).at(read))
				{
					if (std::get<0>(bordert) < std::get<0>(partt) || (std::get<0>(bordert) == std::get<0>(partt) && std::get<1>(bordert) < std::get<1>(partt)))
					{
						uint64_t node = readPaths[read][std::get<0>(bordert)].path[std::get<1>(bordert)];
						size_t foundTangle = find(parent, 2*(node & maskUint64_t) + ((node & firstBitUint64_t) ? 1 : 0));
						if (foundTangle != tangle) continue;
						alleleIndexBefore = find(longUnitigTipParent, node);
					}
					if (alleleIndexAfter == std::numeric_limits<size_t>::max())
					{
						if (std::get<0>(bordert) > std::get<0>(partt) || (std::get<0>(bordert) == std::get<0>(partt) && std::get<1>(bordert) > std::get<1>(partt)))
						{
							uint64_t node = readPaths[read][std::get<0>(bordert)].path[std::get<1>(bordert)] ^ firstBitUint64_t;
							size_t foundTangle = find(parent, 2*(node & maskUint64_t) + ((node & firstBitUint64_t) ? 1 : 0));
							if (foundTangle != tangle) continue;
							alleleIndexAfter = find(longUnitigTipParent, node);
						}
					}
				}
				assert(alleleIndexAfter == alleleIndexBefore || alleleIndexAfter == std::numeric_limits<size_t>::max() || alleleIndexBefore == std::numeric_limits<size_t>::max());
				if (alleleIndexAfter == std::numeric_limits<size_t>::max() && alleleIndexBefore == std::numeric_limits<size_t>::max()) continue;
				size_t allele = alleleIndexBefore;
				if (alleleIndexBefore == std::numeric_limits<size_t>::max()) allele = alleleIndexAfter;
				for (size_t j = std::get<0>(partt); j <= std::get<2>(partt); j++)
				{
					for (size_t k = (j == std::get<0>(partt) ? std::get<1>(partt) : 0); k <= (j == std::get<2>(partt) ? std::get<3>(partt) : readPaths[read][j].path.size()-1); k++)
					{
						size_t node = readPaths[read][j].path[k] & maskUint64_t;
						if (innerNodesBelongToResolution.count(node) == 1 && innerNodesBelongToResolution.at(node) != allele)
						{
							innerNodesBelongToResolution[node] = std::numeric_limits<size_t>::max();
						}
						else
						{
							innerNodesBelongToResolution[node] = allele;
						}
					}
				}
				readTanglePartAssignments[read].emplace_back(std::get<0>(partt), std::get<1>(partt), std::get<2>(partt), std::get<3>(partt), allele);
			}
		}
		for (const auto& pair : readPartsPerTangle.at(tangle))
		{
			size_t read = pair.first;
			for (const auto& partt : readPartsPerTangle.at(tangle).at(read))
			{
				bool alreadyAssigned = false;
				for (auto t : readTanglePartAssignments[read])
				{
					if (std::get<0>(t) == std::get<0>(partt) && std::get<1>(t) == std::get<1>(partt) && std::get<2>(t) == std::get<2>(partt) && std::get<3>(t) == std::get<3>(partt))
					{
						alreadyAssigned = true;
						break;
					}
				}
				if (alreadyAssigned) continue;
				phmap::flat_hash_set<size_t> potentialAlleles;
				for (size_t j = std::get<0>(partt); j <= std::get<2>(partt); j++)
				{
					for (size_t k = (j == std::get<0>(partt) ? std::get<1>(partt) : 0); k <= (j == std::get<2>(partt) ? std::get<3>(partt) : readPaths[read][j].path.size()-1); k++)
					{
						size_t node = readPaths[read][j].path[k] & maskUint64_t;
						if (innerNodesBelongToResolution.count(node) == 0) continue;
						if (innerNodesBelongToResolution.at(node) == std::numeric_limits<size_t>::max()) continue;
						potentialAlleles.emplace(innerNodesBelongToResolution.at(node));
					}
				}
				if (potentialAlleles.size() != 1) continue;
				size_t allele = *potentialAlleles.begin();
				readTanglePartAssignments[read].emplace_back(std::get<0>(partt), std::get<1>(partt), std::get<2>(partt), std::get<3>(partt), allele);
			}
		}
	}
	size_t nextNum = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				continue;
			}
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			nextNum = std::max(nextNum, chunk);
		}
	}
	nextNum += 1;
	phmap::flat_hash_set<size_t> chunksWillBeReplaced;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (unitigIsLong[i]) continue;
		assert(find(parent, i*2) == find(parent, i*2+1));
		if (tangles.count(find(parent, i*2)) == 0) continue;
		for (uint64_t chunk : graph.chunksInUnitig[i])
		{
			chunksWillBeReplaced.emplace(chunk & maskUint64_t);
		}
	}
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> chunkPlusAlleleToNewChunk;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t>> replacementSpans;
		for (auto t : readTanglePartAssignments[i])
		{
			size_t chunkStart = readPaths[i][std::get<0>(t)].chunksInPathnode[std::get<1>(t)].first;
			size_t chunkEnd = readPaths[i][std::get<2>(t)].chunksInPathnode[std::get<3>(t)].second;
			size_t allele = std::get<4>(t);
			assert(chunkEnd >= chunkStart);
			assert(chunkEnd < chunksPerRead[i].size());
			for (size_t j = chunkStart; j <= chunkEnd; j++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
				if (chunksWillBeReplaced.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
				std::pair<size_t, size_t> replacement { std::get<2>(chunksPerRead[i][j]) & maskUint64_t, allele };
				size_t newChunk = nextNum;
				if (chunkPlusAlleleToNewChunk.count(replacement) == 0)
				{
					chunkPlusAlleleToNewChunk[replacement] = newChunk;
					nextNum += 1;
				}
				else
				{
					newChunk = chunkPlusAlleleToNewChunk.at(replacement);
				}
				std::get<2>(chunksPerRead[i][j]) = newChunk + (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t);
			}
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (chunksWillBeReplaced.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			std::get<2>(chunksPerRead[i][j]) = std::numeric_limits<size_t>::max();
		}
	}
}
