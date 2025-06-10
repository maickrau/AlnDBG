#include <map>
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

void splitPerLength(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double differenceFraction, const size_t differenceConstant, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by length" << std::endl;
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
		std::vector<size_t> orderedOccurrences;
		for (size_t j = 0; j < countsPerOccurrence.size(); j++)
		{
			orderedOccurrences.emplace_back(j);
		}
		std::sort(orderedOccurrences.begin(), orderedOccurrences.end(), [&countsPerOccurrence](const size_t left, const size_t right) { return countsPerOccurrence[left] < countsPerOccurrence[right]; });
		std::vector<size_t> parent;
		for (size_t j = 0; j < countsPerOccurrence.size(); j++)
		{
			parent.emplace_back(j);
		}
		std::vector<size_t> uniqueIndices;
		uniqueIndices.emplace_back(orderedOccurrences[0]);
		for (size_t j = 1; j < orderedOccurrences.size(); j++)
		{
			if (countsPerOccurrence[orderedOccurrences[j-1]] != countsPerOccurrence[orderedOccurrences[j]])
			{
				uniqueIndices.emplace_back(orderedOccurrences[j]);
				continue;
			}
			merge(parent, orderedOccurrences[j-1], orderedOccurrences[j]);
		}
		std::string info = "chunk with coverage " + std::to_string(parent.size()) + " has " + std::to_string(uniqueIndices.size()) + " uniques\n";
		std::cerr << info;
		std::vector<size_t> uniqueParent = getFastTransitiveClosure(uniqueIndices.size(), maxDistance, [&countsPerOccurrence, &uniqueIndices](const size_t left, const size_t right, const size_t maxDist)
		{
			size_t result = 0;
			for (size_t c = 0; c < 4; c++)
			{
				if (countsPerOccurrence[uniqueIndices[left]][c] > countsPerOccurrence[uniqueIndices[right]][c])
				{
					result += countsPerOccurrence[uniqueIndices[left]][c] - countsPerOccurrence[uniqueIndices[right]][c];
				}
				else
				{
					result += countsPerOccurrence[uniqueIndices[right]][c] - countsPerOccurrence[uniqueIndices[left]][c];
				}
			}
			return result;
		});
		for (size_t i = 0; i < uniqueParent.size(); i++)
		{
			merge(parent, uniqueIndices[i], uniqueIndices[find(uniqueParent, i)]);
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
	assert(chunksPerRead.size()%2 == 0);
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
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
	assert(chunksPerRead.size()%2 == 0);
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
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
	assert(chunksPerRead.size()%2 == 0);
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
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
	assert(chunksPerRead.size()%2 == 0);
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
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

void resolveBetweenTanglesInner(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t longUnitigThreshold, const size_t kmerSize, const bool allowGaps)
{
	std::cerr << "resolving between tangles" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	assert(chunksPerRead.size()%2 == 0);
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
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
	phmap::flat_hash_set<size_t> tipsConnectAcrossTangles;
	if (allowGaps)
	{
		std::vector<std::pair<uint64_t, uint64_t>> merges;
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
					if (lastLongUnitig != std::numeric_limits<size_t>::max())
					{
						if (k >= 1 || readPaths[i][j].pathLeftClipChunks == 0)
						{
							if (find(parent, 2 * (lastLongUnitig & maskUint64_t) + ((lastLongUnitig & firstBitUint64_t) ? 1 : 0)) != find(parent, 2 * (longUnitigHere & maskUint64_t) + ((longUnitigHere & firstBitUint64_t) ? 0 : 1)))
							{
								merges.emplace_back(2 * (lastLongUnitig & maskUint64_t) + ((lastLongUnitig & firstBitUint64_t) ? 1 : 0), 2 * (longUnitigHere & maskUint64_t) + ((longUnitigHere & firstBitUint64_t) ? 0 : 1));
								tipsConnectAcrossTangles.emplace(2 * (lastLongUnitig & maskUint64_t) + ((lastLongUnitig & firstBitUint64_t) ? 1 : 0));
								tipsConnectAcrossTangles.emplace(2 * (longUnitigHere & maskUint64_t) + ((longUnitigHere & firstBitUint64_t) ? 0 : 1));
							}
						}
					}
					lastLongUnitig = std::numeric_limits<size_t>::max();
					if (k+1 < readPaths[i][j].path.size() || readPaths[i][j].pathRightClipChunks == 0)
					{
						lastLongUnitig = longUnitigHere;
					}
				}
			}
		}
		for (auto pair : merges)
		{
			merge(parent, pair.first, pair.second);
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
	phmap::flat_hash_set<size_t> tanglesWithoutAnyLongNodes;
	for (size_t tangle : tangles)
	{
		if (longUnitigsPerTangle.count(tangle) == 0) tanglesWithoutAnyLongNodes.insert(tangle);
	}
	for (size_t tangle : tanglesWithoutAnyLongNodes)
	{
		assert(tangles.count(tangle) == 1);
		tangles.erase(tangle);
	}
	std::cerr << tangles.size() << " tangles after filtering for long unitig touchers" << std::endl;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, size_t>> connectionCounts;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		uint64_t lastLongUnitig = std::numeric_limits<uint64_t>::max();
		size_t lastLongUnitigChunkPos = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				size_t unitig = readPaths[i][j].path[k] & maskUint64_t;
				assert(unitig < unitigIsLong.size());
				if (!unitigIsLong[unitig]) continue;
				uint64_t longUnitigHere = readPaths[i][j].path[k];
				size_t chunkPosHere = graph.chunksInUnitig[longUnitigHere & maskUint64_t].size();
				if (k+1 == readPaths[i][j].path.size())
				{
					assert(chunkPosHere >= readPaths[i][j].pathRightClipChunks);
					chunkPosHere -= readPaths[i][j].pathRightClipChunks;
				}
				if (lastLongUnitig == std::numeric_limits<size_t>::max())
				{
					lastLongUnitig = longUnitigHere;
					lastLongUnitigChunkPos = chunkPosHere;
					continue;
				}
				if (lastLongUnitig == longUnitigHere)
				{
					if (chunkPosHere > lastLongUnitigChunkPos)
					{
						lastLongUnitig = longUnitigHere;
						lastLongUnitigChunkPos = chunkPosHere;
						continue;
					}
				}
				if (find(parent, 2 * (lastLongUnitig & maskUint64_t) + ((lastLongUnitig & firstBitUint64_t) ? 1 : 0)) != find(parent, 2 * (longUnitigHere & maskUint64_t) + ((longUnitigHere & firstBitUint64_t) ? 0 : 1)))
				{
					lastLongUnitig = longUnitigHere;
					lastLongUnitigChunkPos = chunkPosHere;
					continue;
				}
				connectionCounts[lastLongUnitig][longUnitigHere ^ firstBitUint64_t] += 1;
				connectionCounts[longUnitigHere ^ firstBitUint64_t][lastLongUnitig] += 1;
				lastLongUnitig = longUnitigHere;
				lastLongUnitigChunkPos = chunkPosHere;
			}
		}
	}
	phmap::flat_hash_set<size_t> forbiddenTangles;
	for (size_t i = 0; i < unitigIsLong.size(); i++)
	{
		if (!unitigIsLong[i]) continue;
		if (connectionCounts.count(i + firstBitUint64_t) == 0)
		{
			forbiddenTangles.emplace(find(parent, 2*i+1));
		}
		else
		{
			bool good = true;
			size_t minCoverage = 4;
			if (connectionCounts.at(i + firstBitUint64_t).size() == 1)
			{
				minCoverage = 2;
			}
			for (auto pair : connectionCounts.at(i + firstBitUint64_t))
			{
				if (pair.second < minCoverage) good = false;
			}
			if (!good) forbiddenTangles.emplace(find(parent, 2*i+1));
		}
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
	}
	for (auto tangle : forbiddenTangles)
	{
		if (tangles.count(tangle) == 1) tangles.erase(tangle);
	}
	std::cerr << tangles.size() << " tangles after filtering for connections" << std::endl;
	phmap::flat_hash_set<size_t> tanglesWhichHaveShortNodes;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (unitigIsLong[i]) continue;
		assert(find(parent, 2*i) == find(parent, 2*i+1));
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
				assert(find(parent, 2*readPaths[read][std::get<0>(partt)].path[std::get<1>(partt)]) == tangle);
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
				if (alleleIndexAfter == std::numeric_limits<size_t>::max() && alleleIndexBefore == std::numeric_limits<size_t>::max()) continue;
				if (alleleIndexAfter != std::numeric_limits<size_t>::max() && alleleIndexBefore != std::numeric_limits<size_t>::max() && alleleIndexAfter != alleleIndexBefore)
				{
					readTanglePartAssignments[read].emplace_back(std::get<0>(partt), std::get<1>(partt), std::get<2>(partt), std::get<3>(partt), std::numeric_limits<size_t>::max());
					continue;
				}
				assert(alleleIndexAfter == alleleIndexBefore || alleleIndexAfter == std::numeric_limits<size_t>::max() || alleleIndexBefore == std::numeric_limits<size_t>::max());
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
				if (allele == std::numeric_limits<size_t>::max())
				{
					std::get<2>(chunksPerRead[i][j]) = std::numeric_limits<size_t>::max();
					continue;
				}
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

void resolveBetweenTangles(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t longUnitigThreshold, const size_t kmerSize)
{
	resolveBetweenTanglesInner(chunksPerRead, approxOneHapCoverage, longUnitigThreshold, kmerSize, false);
}

void resolveBetweenTanglesAllowGaps(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t longUnitigThreshold, const size_t kmerSize)
{
	resolveBetweenTanglesInner(chunksPerRead, approxOneHapCoverage, longUnitigThreshold, kmerSize, true);
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

bool expandResolvableChunksOnce(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize, const size_t expandedSize)
{
	const size_t minSolidCoverage = 4;
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			maxChunk = std::max(maxChunk, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
		}
	}
	maxChunk += 1;
	std::vector<std::pair<size_t, size_t>> chunkPositionInUnitig;
	std::vector<size_t> unitigChunkCounts;
	chunkPositionInUnitig.resize(maxChunk, std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	std::vector<bool> chunkCanExpand;
	chunkCanExpand.resize(maxChunk, false);
	{
		ChunkUnitigGraph graph;
		std::vector<std::vector<UnitigPath>> readPaths;
		assert(chunksPerRead.size()%2 == 0);
		std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
		unitigChunkCounts.resize(graph.unitigLengths.size(), 0);
		std::vector<bool> canExpandUnitig;
		canExpandUnitig.resize(graph.unitigLengths.size(), false);
		for (size_t i = 0; i < graph.unitigLengths.size(); i++)
		{
			unitigChunkCounts[i] = graph.chunksInUnitig[i].size();
			if (graph.unitigLengths[i] > expandedSize) continue;
			if (graph.edges.getEdges(std::make_pair(i, true)).size() == 1 || graph.edges.getEdges(std::make_pair(i, false)).size() == 1) continue;
			if (graph.edges.getEdges(std::make_pair(i, true)).size() == 0 && graph.edges.getEdges(std::make_pair(i, false)).size() == 0) continue;
			canExpandUnitig[i] = true;
		}
		phmap::flat_hash_map<size_t, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> tripletCoverages;
		for (size_t i = 0; i < readPaths.size(); i++)
		{
			for (size_t j = 0; j < readPaths[i].size(); j++)
			{
				for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
				{
					assert(readPaths[i][j].path[k] & firstBitUint64_t);
					size_t unitig = readPaths[i][j].path[k] & maskUint64_t;
					if (!canExpandUnitig[unitig]) continue;
					uint64_t prev = readPaths[i][j].path[k-1];
					uint64_t next = readPaths[i][j].path[k+1];
					if ((readPaths[i][j].path[k] ^ firstBitUint64_t) & firstBitUint64_t)
					{
						std::swap(prev, next);
						prev ^= firstBitUint64_t;
						next ^= firstBitUint64_t;
					}
					tripletCoverages[unitig][std::make_pair(prev, next)] += 1;
				}
			}
		}
		for (const auto& pair : tripletCoverages)
		{
			assert(canExpandUnitig[pair.first]);
			canExpandUnitig[pair.first] = false;
			if (tripletCoverages[pair.first].size() < graph.edges.getEdges(std::make_pair(pair.first, true)).size()) continue;
			if (tripletCoverages[pair.first].size() < graph.edges.getEdges(std::make_pair(pair.first, false)).size()) continue;
			phmap::flat_hash_set<uint64_t> solidPredecessors;
			phmap::flat_hash_set<uint64_t> solidSuccessors;
			bool valid = true;
			for (auto pair2 : pair.second)
			{
				if (pair2.second == 1) continue;
				if (pair2.second < minSolidCoverage)
				{
					valid = false;
					break;
				}
				assert(pair2.second >= minSolidCoverage);
				solidPredecessors.insert(pair2.first.first);
				solidSuccessors.insert(pair2.first.second);
			}
			assert(solidPredecessors.size() <= graph.edges.getEdges(std::make_pair(pair.first, false)).size());
			assert(solidSuccessors.size() <= graph.edges.getEdges(std::make_pair(pair.first, true)).size());
			if (solidPredecessors.size() != graph.edges.getEdges(std::make_pair(pair.first, false)).size()) continue;
			if (solidSuccessors.size() != graph.edges.getEdges(std::make_pair(pair.first, true)).size()) continue;
			if (!valid) continue;
			canExpandUnitig[pair.first] = true;
		}
		for (size_t i = 0; i < graph.unitigLengths.size(); i++)
		{
			if (!canExpandUnitig[i]) continue;
			for (size_t j = 0; j < graph.chunksInUnitig[i].size(); j++)
			{
				assert(graph.chunksInUnitig[i][j] & firstBitUint64_t);
				size_t chunk = graph.chunksInUnitig[i][j] & maskUint64_t;
				assert(chunk < chunkCanExpand.size());
				assert(!chunkCanExpand[chunk]);
				chunkCanExpand[chunk] = true;
				chunkPositionInUnitig[chunk] = std::make_pair(i, j);
			}
		}
	}
	std::map<std::vector<uint64_t>, size_t> chunkmerToNewChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<std::pair<size_t, size_t>> newChunkmerRanges;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || !chunkCanExpand[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) newChunkmerRanges.emplace_back(j, j);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || !chunkCanExpand[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			size_t end = j;
			while (end < chunksPerRead[i].size())
			{
				end += 1;
				if (end == chunksPerRead[i].size()) break;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end]))) break;
				if (!chunkCanExpand[std::get<2>(chunksPerRead[i][end]) & maskUint64_t]) break;
				if (chunkPositionInUnitig[std::get<2>(chunksPerRead[i][end]) & maskUint64_t].first != chunkPositionInUnitig[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].first) break;
				if (chunkPositionInUnitig[std::get<2>(chunksPerRead[i][end]) & maskUint64_t].second != chunkPositionInUnitig[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].second + (end - j)) break;
			}
			if (end == chunksPerRead[i].size()) continue;
			assert(end > j);
			if (chunkPositionInUnitig[std::get<2>(chunksPerRead[i][end]) & maskUint64_t].first == chunkPositionInUnitig[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].first)
			{
				if (chunkPositionInUnitig[std::get<2>(chunksPerRead[i][end]) & maskUint64_t].second == chunkPositionInUnitig[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].second + (end - j))
				{
					continue;
				}
			}
			newChunkmerRanges.emplace_back(j, end);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || !chunkCanExpand[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			size_t start = j;
			while (true)
			{
				if (start == 0)
				{
					start = std::numeric_limits<size_t>::max();
					break;
				}
				start -= 1;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][start]))) break;
				if (!chunkCanExpand[std::get<2>(chunksPerRead[i][start]) & maskUint64_t]) break;
				if (chunkPositionInUnitig[std::get<2>(chunksPerRead[i][start]) & maskUint64_t].first != chunkPositionInUnitig[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].first) break;
				if (chunkPositionInUnitig[std::get<2>(chunksPerRead[i][start]) & maskUint64_t].second + (j - start) != chunkPositionInUnitig[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].second) break;
			}
			if (start == std::numeric_limits<size_t>::max()) continue;
			assert(start < chunksPerRead[i].size());
			assert(start < j);
			if (chunkPositionInUnitig[std::get<2>(chunksPerRead[i][start]) & maskUint64_t].first == chunkPositionInUnitig[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].first)
			{
				if (chunkPositionInUnitig[std::get<2>(chunksPerRead[i][start]) & maskUint64_t].second + (j - start) == chunkPositionInUnitig[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].second)
				{
					continue;
				}
			}
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
			if (!(newChunkmerRanges[j].first <= newChunkmerRanges[j-1].second+1))
			{
				std::cerr << i << " " << j << std::endl;
				for (size_t k = 0; k < newChunkmerRanges.size(); k++)
				{
					std::cerr << newChunkmerRanges[k].first << "-" << newChunkmerRanges[k].second << " ";
				}
				std::cerr << std::endl;
			}
			assert(newChunkmerRanges[j].first <= newChunkmerRanges[j-1].second+1);
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
	size_t countNews = 0;
	for (const auto& pair : chunkmerToNewChunk)
	{
		assert(pair.first.size() >= 1);
		if (pair.first.size() >= 2)
		{
			countNews += 1;
		}
	}
	std::cerr << "resolved " << countNews << " new unitigs" << std::endl;
	return (countNews >= 1);
}

void expandResolvableChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize, const size_t expandedSize)
{
	std::cerr << "resolve resolvable chunks" << std::endl;
	while (true)
	{
		bool resolvedAnything = expandResolvableChunksOnce(chunksPerRead, approxOneHapCoverage, kmerSize, expandedSize);
		if (!resolvedAnything) break;
	}
}

void expandChunksUntilSolids(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize, const size_t expandedSize)
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
	std::vector<bool> chunkIsAnchor;
	chunkIsAnchor.resize(chunkLengths.size(), true);
	{
		ChunkUnitigGraph graph;
		std::vector<std::vector<UnitigPath>> readPaths;
		assert(chunksPerRead.size()%2 == 0);
		std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
		for (size_t i = 0; i < graph.unitigLengths.size(); i++)
		{
			if (graph.unitigLengths[i] > expandedSize) continue;
			for (size_t j = 0; j < graph.chunksInUnitig[i].size(); j++)
			{
				size_t chunk = graph.chunksInUnitig[i][j] & maskUint64_t;
				assert(chunk < chunkIsAnchor.size());
				assert(chunkIsAnchor[chunk]);
				chunkIsAnchor[chunk] = false;
			}
		}
	}
	std::map<std::vector<uint64_t>, size_t> chunkmerToNewChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<std::pair<size_t, size_t>> newChunkmerRanges;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || chunkIsAnchor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) newChunkmerRanges.emplace_back(j, j);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t end = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			while (end < chunksPerRead[i].size() && length < expandedSize)
			{
				if (end+1 == chunksPerRead[i].size()) break;
				end += 1;
				length += chunkLengths[std::get<2>(chunksPerRead[i][end]) & maskUint64_t];
				length -= kmerSize;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end]))) break;
				if (chunkIsAnchor[std::get<2>(chunksPerRead[i][end]) & maskUint64_t]) break;
			}
			assert(end < chunksPerRead[i].size());
			bool startAnchor = false;
			bool endAnchor = false;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || chunkIsAnchor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) startAnchor = true;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][end])) || chunkIsAnchor[std::get<2>(chunksPerRead[i][end]) & maskUint64_t]) endAnchor = true;
			if (startAnchor && endAnchor && end == j+1) continue;
			if (length < expandedSize && (!startAnchor || !endAnchor)) continue;
			if (end == chunksPerRead[i].size()) continue;
			newChunkmerRanges.emplace_back(j, end);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t start = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			while (start > 0 && length < expandedSize)
			{
				start -= 1;
				length += chunkLengths[std::get<2>(chunksPerRead[i][start]) & maskUint64_t];
				length -= kmerSize;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][start]))) break;
				if (chunkIsAnchor[std::get<2>(chunksPerRead[i][start]) & maskUint64_t]) break;
			}
			assert(start < chunksPerRead[i].size());
			bool startAnchor = false;
			bool endAnchor = false;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][start])) || chunkIsAnchor[std::get<2>(chunksPerRead[i][start]) & maskUint64_t]) startAnchor = true;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || chunkIsAnchor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) endAnchor = true;
			if (startAnchor && endAnchor && j == start+1) continue;
			if (length < expandedSize && (!startAnchor || !endAnchor)) continue;
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
			if (!(newChunkmerRanges[j].first <= newChunkmerRanges[j-1].second+1))
			{
				std::cerr << i << " " << j << std::endl;
				for (size_t k = 0; k < newChunkmerRanges.size(); k++)
				{
					std::cerr << newChunkmerRanges[k].first << "-" << newChunkmerRanges[k].second << " ";
				}
				std::cerr << std::endl;
			}
			assert(newChunkmerRanges[j].first <= newChunkmerRanges[j-1].second+1);
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

void addGapChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() < 2) continue;
		bool addedAny = false;
		for (size_t j = chunksPerRead[i].size()-1; j > 0; j--)
		{
			if (std::get<1>(chunksPerRead[i][j-1]) >= std::get<0>(chunksPerRead[i][j])+kmerSize-1) continue;
			addedAny = true;
			chunksPerRead[i].emplace_back(std::get<1>(chunksPerRead[i][j-1])-kmerSize+1, std::get<0>(chunksPerRead[i][j])+kmerSize-1, std::numeric_limits<size_t>::max());
		}
		if (addedAny)
		{
			std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
		}
	}
}

void expandChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t expandedSize)
{
	std::vector<size_t> chunkLengths;
	std::vector<bool> keepChunkAsIs;
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
		keepChunkAsIs.resize(rawChunkLengths.size(), false);
		for (size_t i = 0; i < rawChunkLengths.size(); i++)
		{
			assert(rawChunkLengths[i].size() >= 1);
			std::sort(rawChunkLengths[i].begin(), rawChunkLengths[i].end());
			chunkLengths[i] = rawChunkLengths[i][rawChunkLengths[i].size()/2];
			if (rawChunkLengths[i].size() == 1) keepChunkAsIs[i] = true;
		}
	}
	std::map<std::vector<uint64_t>, size_t> chunkmerToNewChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<std::pair<size_t, size_t>> newChunkmerRanges;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || keepChunkAsIs[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) newChunkmerRanges.emplace_back(j, j);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || keepChunkAsIs[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			size_t end = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			while (end < chunksPerRead[i].size() && length < expandedSize)
			{
				if (end+1 == chunksPerRead[i].size()) break;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end+1])) || keepChunkAsIs[std::get<2>(chunksPerRead[i][end+1]) & maskUint64_t]) break;
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
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || keepChunkAsIs[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			size_t start = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			while (start > 0 && length < expandedSize)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][start-1])) || keepChunkAsIs[std::get<2>(chunksPerRead[i][start-1]) & maskUint64_t]) break;
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
	for (size_t i = 0; i < complementChunk.size(); i++)
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
	std::vector<bool> chunkIsAnchor;
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
		chunkIsAnchor.resize(rawChunkLengths.size(), false);
		for (size_t i = 0; i < rawChunkLengths.size(); i++)
		{
			assert(rawChunkLengths[i].size() >= 1);
			if (rawChunkLengths[i].size() == 1) chunkIsAnchor[i] = true;
			std::sort(rawChunkLengths[i].begin(), rawChunkLengths[i].end());
			chunkLengths[i] = rawChunkLengths[i][rawChunkLengths[i].size()/2];
		}
	}
	std::map<std::vector<uint64_t>, size_t> chunkmerToNewChunk;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunkParent;
	std::vector<std::vector<std::pair<size_t, size_t>>> newChunkParent;
	chunkParent.resize(chunksPerRead.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		chunkParent[i].resize(chunksPerRead[i].size());
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			chunkParent[i][j] = std::make_pair(i, j);
		}
	}
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
			if (chunkIsAnchor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			size_t end = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			while (end < chunksPerRead[i].size() && length < expandedSize)
			{
				if (end+1 == chunksPerRead[i].size()) break;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end+1]))) break;
				if (chunkIsAnchor[std::get<2>(chunksPerRead[i][end+1]) & maskUint64_t]) break;
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
			if (chunkIsAnchor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			size_t start = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			while (start > 0 && length < expandedSize)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][start-1]))) break;
				if (chunkIsAnchor[std::get<2>(chunksPerRead[i][start-1]) & maskUint64_t]) break;
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
					size_t newChunkIndex = chunkmerToNewChunk.at(chunkmer);
					newchunk = newChunkIndex;
					for (size_t k = 0; k < chunkmer.size(); k++)
					{
						merge(chunkParent, std::make_pair(i, newChunkmerRanges[j].first+k), newChunkParent[newchunk][k]);
					}
				}
				else
				{
					newchunk = chunkmerToNewChunk.size();
					chunkmerToNewChunk[chunkmer] = newchunk;
					newChunkParent.emplace_back();
					for (size_t k = 0; k < chunkmer.size(); k++)
					{
						newChunkParent.back().emplace_back(i, newChunkmerRanges[j].first+k);
					}
				}
			}
			else
			{
				newchunk = std::numeric_limits<size_t>::max();
			}
		}
	}
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> newChunkClusterToIndex;
	size_t nextNum = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (chunkParent[i][j].first != i || chunkParent[i][j].second != j) continue;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			newChunkClusterToIndex[chunkParent[i][j]] = nextNum;
			nextNum += 1;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			std::get<2>(chunksPerRead[i][j]) = newChunkClusterToIndex.at(find(chunkParent, chunkParent[i][j])) + firstBitUint64_t;
		}
	}
}

void fragmentChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const size_t kmerSize)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> parent;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunkParent;
	parent.resize(chunksPerRead.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (minimizerPositionsPerRead[i].size() < 2)
		{
			assert(chunksPerRead[i].size() == 0);
			continue;
		}
		parent[i].resize(minimizerPositionsPerRead[i].size()-1);
		for (size_t j = 0; j < minimizerPositionsPerRead[i].size()-1; j++)
		{
			parent[i][j] = std::make_pair(i, j);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			while (chunkParent.size() <= chunk) chunkParent.emplace_back();
			if (chunkParent[chunk].size() == 0)
			{
				size_t firstMinimizer = std::numeric_limits<size_t>::max();
				size_t lastMinimizer = std::numeric_limits<size_t>::max();
				for (size_t k = 0; k < minimizerPositionsPerRead[i].size(); k++)
				{
					if (minimizerPositionsPerRead[i][k] == std::get<0>(chunksPerRead[i][j]))
					{
						firstMinimizer = k;
					}
					if (minimizerPositionsPerRead[i][k] + kmerSize - 1 == std::get<1>(chunksPerRead[i][j]))
					{
						lastMinimizer = k;
					}
				}
				assert(firstMinimizer != std::numeric_limits<size_t>::max());
				assert(lastMinimizer != std::numeric_limits<size_t>::max());
				assert(lastMinimizer > firstMinimizer);
				for (size_t k = firstMinimizer; k < lastMinimizer; k++)
				{
					chunkParent[chunk].emplace_back(std::make_pair(i, k));
				}
			}
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			assert(chunkParent.size() > chunk);
			assert(chunkParent[chunk].size() >= 1);
			size_t firstMinimizer = std::numeric_limits<size_t>::max();
			size_t lastMinimizer = std::numeric_limits<size_t>::max();
			for (size_t k = 0; k < minimizerPositionsPerRead[i].size(); k++)
			{
				if (minimizerPositionsPerRead[i][k] == std::get<0>(chunksPerRead[i][j]))
				{
					firstMinimizer = k;
				}
				if (minimizerPositionsPerRead[i][k] + kmerSize - 1 == std::get<1>(chunksPerRead[i][j]))
				{
					lastMinimizer = k;
				}
			}
			assert(firstMinimizer != std::numeric_limits<size_t>::max());
			assert(lastMinimizer != std::numeric_limits<size_t>::max());
			assert(lastMinimizer > firstMinimizer);
			assert(lastMinimizer - firstMinimizer == chunkParent[chunk].size());
			for (size_t k = firstMinimizer; k < lastMinimizer; k++)
			{
				merge(parent, std::make_pair(i, k), chunkParent[chunk][k-firstMinimizer]);
			}
		}
	}
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> newChunks;
	size_t nextNum = 0;
	for (size_t i = 0; i < parent.size(); i++)
	{
		for (size_t j = 0; j < parent[i].size(); j++)
		{
			if (parent[i][j].first != i || parent[i][j].second != j) continue;
			newChunks[std::make_pair(i, j)] = nextNum;
			nextNum += 1;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (parent[i].size() == 0)
		{
			assert(chunksPerRead[i].size() == 0);
			assert(minimizerPositionsPerRead[i].size() < 2);
			continue;
		}
		assert(parent[i].size()+1 == minimizerPositionsPerRead[i].size());
		chunksPerRead[i].clear();
		for (size_t j = 0; j < parent[i].size(); j++)
		{
			chunksPerRead[i].emplace_back(minimizerPositionsPerRead[i][j], minimizerPositionsPerRead[i][j+1] + kmerSize - 1, newChunks.at(find(parent, std::make_pair(i, j))) + firstBitUint64_t);
		}
	}
}

void resplitFalselyMergedChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& oldChunks)
{
	phmap::flat_hash_map<size_t, size_t> oneMatch;
	phmap::flat_hash_map<size_t, phmap::flat_hash_set<size_t>> extraMatches;
	assert(chunksPerRead.size() == oldChunks.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		size_t newIndex = 0;
		size_t oldIndex = 0;
		while (newIndex < chunksPerRead[i].size() && oldIndex < oldChunks[i].size())
		{
			if (std::get<0>(chunksPerRead[i][newIndex]) < std::get<0>(oldChunks[i][oldIndex]))
			{
				newIndex += 1;
				continue;
			}
			if (std::get<0>(chunksPerRead[i][newIndex]) > std::get<0>(oldChunks[i][oldIndex]))
			{
				oldIndex += 1;
				continue;
			}
			assert(std::get<0>(chunksPerRead[i][newIndex]) == std::get<0>(oldChunks[i][oldIndex]));
			if (std::get<1>(chunksPerRead[i][newIndex]) < std::get<1>(oldChunks[i][oldIndex]))
			{
				newIndex += 1;
				continue;
			}
			if (std::get<1>(chunksPerRead[i][newIndex]) > std::get<1>(oldChunks[i][oldIndex]))
			{
				oldIndex += 1;
				continue;
			}
			assert(std::get<1>(chunksPerRead[i][newIndex]) == std::get<1>(oldChunks[i][oldIndex]));
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][newIndex])))
			{
				newIndex += 1;
				continue;
			}
			if (NonexistantChunk(std::get<2>(oldChunks[i][oldIndex])))
			{
				oldIndex += 1;
				continue;
			}
			size_t newChunk = std::get<2>(chunksPerRead[i][newIndex]) & maskUint64_t;
			size_t oldChunk = std::get<2>(oldChunks[i][oldIndex]) & maskUint64_t;
			newIndex += 1;
			oldIndex += 1;
			if (oneMatch.count(newChunk) == 0)
			{
				oneMatch[newChunk] = oldChunk;
				continue;
			}
			if (oneMatch.at(newChunk) == oldChunk)
			{
				continue;
			}
			extraMatches[newChunk].insert(oldChunk);
		}
	}
	std::cerr << "re-split " << extraMatches.size() << " chunks" << std::endl;
	if (extraMatches.size() == 0) return;
	size_t nextNum = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			nextNum = std::max(nextNum, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
		}
	}
	nextNum += 1;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> remapping;
	for (const auto& pair : extraMatches)
	{
		assert(oneMatch.count(pair.first) == 1);
		std::cerr << "re-split chunk " << pair.first << " to " << (pair.second.size()+1) << " chunks" << std::endl;
		for (size_t old : pair.second)
		{
			remapping[std::make_pair(pair.first, old)] = nextNum;
			nextNum += 1;
		}
		remapping[std::make_pair(pair.first, oneMatch.at(pair.first))] = nextNum;
		nextNum += 1;
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		size_t newIndex = 0;
		size_t oldIndex = 0;
		while (newIndex < chunksPerRead[i].size() && oldIndex < oldChunks[i].size())
		{
			if (std::get<0>(chunksPerRead[i][newIndex]) < std::get<0>(oldChunks[i][oldIndex]))
			{
				newIndex += 1;
				continue;
			}
			if (std::get<0>(chunksPerRead[i][newIndex]) > std::get<0>(oldChunks[i][oldIndex]))
			{
				oldIndex += 1;
				continue;
			}
			assert(std::get<0>(chunksPerRead[i][newIndex]) == std::get<0>(oldChunks[i][oldIndex]));
			if (std::get<1>(chunksPerRead[i][newIndex]) < std::get<1>(oldChunks[i][oldIndex]))
			{
				newIndex += 1;
				continue;
			}
			if (std::get<1>(chunksPerRead[i][newIndex]) > std::get<1>(oldChunks[i][oldIndex]))
			{
				oldIndex += 1;
				continue;
			}
			assert(std::get<1>(chunksPerRead[i][newIndex]) == std::get<1>(oldChunks[i][oldIndex]));
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][newIndex])))
			{
				newIndex += 1;
				continue;
			}
			if (NonexistantChunk(std::get<2>(oldChunks[i][oldIndex])))
			{
				oldIndex += 1;
				continue;
			}
			size_t newChunk = std::get<2>(chunksPerRead[i][newIndex]) & maskUint64_t;
			size_t oldChunk = std::get<2>(oldChunks[i][oldIndex]) & maskUint64_t;
			assert(extraMatches.count(newChunk) == remapping.count(std::make_pair(newChunk, oldChunk)));
			if (remapping.count(std::make_pair(newChunk, oldChunk)) == 1)
			{
				std::get<2>(chunksPerRead[i][newIndex]) = remapping.at(std::make_pair(newChunk, oldChunk)) + firstBitUint64_t;
			}
			newIndex += 1;
			oldIndex += 1;
		}
	}
}

void removeBadShortHighCoverageChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	std::vector<size_t> shortCoverage;
	std::vector<size_t> longCoverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			while (shortCoverage.size() <= chunk) shortCoverage.emplace_back(0);
			while (longCoverage.size() <= chunk) longCoverage.emplace_back(0);
			if (std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) >= 2*kmerSize)
			{
				longCoverage[chunk] += 1;
			}
			else
			{
				shortCoverage[chunk] += 1;
			}
		}
	}
	std::vector<bool> remove;
	remove.resize(shortCoverage.size(), false);
	for (size_t i = 0; i < shortCoverage.size(); i++)
	{
		if (shortCoverage[i] < 1000000) continue;
		if (longCoverage[i] > 100000) continue;
		std::cerr << "remove short chunk " << i << " with coverages " << shortCoverage[i] << " " << longCoverage[i] << std::endl;
		remove[i] = true;
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) >= 2*kmerSize) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			if (remove[chunk])
			{
				std::get<2>(chunksPerRead[i][j]) = std::numeric_limits<size_t>::max();
			}
		}
	}
}

void mergeNonexistentChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() < 2) continue;
		bool removedAny = false;
		for (size_t j = chunksPerRead[i].size()-1; j > 0; j--)
		{
			if (!NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!NonexistantChunk(std::get<2>(chunksPerRead[i][j-1]))) continue;
			assert(std::get<0>(chunksPerRead[i][j-1]) <= std::get<0>(chunksPerRead[i][j]));
			assert(std::get<1>(chunksPerRead[i][j-1]) <= std::get<1>(chunksPerRead[i][j]));
			std::get<1>(chunksPerRead[i][j-1]) = std::get<1>(chunksPerRead[i][j]);
			std::swap(chunksPerRead[i][j], chunksPerRead[i].back());
			chunksPerRead[i].pop_back();
			removedAny = true;
		}
		if (removedAny) std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
	}
}

std::tuple<uint64_t, uint64_t, uint64_t> checkYForkEdge(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double averageCoverage, const size_t minLongUnitigLength, const std::pair<size_t, bool> startNode)
{
	std::tuple<uint64_t, uint64_t, uint64_t> notFound { std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max() };
	if (graph.edges.getEdges(startNode).size() != 2) return notFound;
	uint64_t firstEdge = std::numeric_limits<uint64_t>::max();
	uint64_t secondEdge = std::numeric_limits<uint64_t>::max();
	for (auto edge : graph.edges.getEdges(startNode))
	{
		if (graph.edges.getEdges(reverse(edge)).size() != 1) return notFound;
		if (graph.unitigLengths[edge.first] < minLongUnitigLength) return notFound;
		if (graph.coverages[edge.first] < 0.5 * averageCoverage) return notFound;
		if (graph.coverages[edge.first] > 1.5 * averageCoverage) return notFound;
		if (firstEdge == std::numeric_limits<uint64_t>::max())
		{
			firstEdge = edge.first + (edge.second ? 0 : firstBitUint64_t); // reversed edge, ie points towards y fork
		}
		else
		{
			secondEdge = edge.first + (edge.second ? 0 : firstBitUint64_t); // reversed edge, ie points towards y fork
		}
	}
	assert(firstEdge != std::numeric_limits<uint64_t>::max());
	assert(secondEdge != std::numeric_limits<uint64_t>::max());
	return std::make_tuple(startNode.first + (startNode.second ? firstBitUint64_t : 0), firstEdge, secondEdge);
}

void fixYForks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize)
{
	const size_t minLongUnitigLength = 100000;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	assert(chunksPerRead.size()%2 == 0);
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
	double coverageSum = 0;
	double coverageDivisor = 0;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < minLongUnitigLength) continue;
		coverageSum += graph.unitigLengths[i] * graph.coverages[i];
		coverageDivisor += graph.unitigLengths[i];
	}
	if (coverageDivisor == 0) return;
	double averageCoverage = coverageSum/coverageDivisor;
	std::cerr << "average coverage " << averageCoverage << std::endl;
	std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> yForks; // stem, first branch, second branch. all three point towards y fork (middle)
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < minLongUnitigLength) continue;
		if (graph.coverages[i] < averageCoverage * 0.5) continue;
		if (graph.coverages[i] > averageCoverage * 1.5) continue;
		auto found = checkYForkEdge(chunksPerRead, graph, readPaths, averageCoverage, minLongUnitigLength, std::make_pair(i, true));
		if (std::get<0>(found) != std::numeric_limits<uint64_t>::max())
		{
			assert(std::get<1>(found) != std::numeric_limits<uint64_t>::max());
			assert(std::get<2>(found) != std::numeric_limits<uint64_t>::max());
			yForks.emplace_back(found);
		}
		found = checkYForkEdge(chunksPerRead, graph, readPaths, averageCoverage, minLongUnitigLength, std::make_pair(i, false));
		if (std::get<0>(found) != std::numeric_limits<uint64_t>::max())
		{
			assert(std::get<1>(found) != std::numeric_limits<uint64_t>::max());
			assert(std::get<2>(found) != std::numeric_limits<uint64_t>::max());
			yForks.emplace_back(found);
		}
	}
	if (yForks.size() == 0) return;
	phmap::flat_hash_map<uint64_t, size_t> unitigBelongsToForkBranch;
	phmap::flat_hash_map<uint64_t, size_t> unitigIsForkStem;
	for (size_t i = 0; i < yForks.size(); i++)
	{
		auto t = yForks[i];
		assert(unitigIsForkStem.count(std::get<0>(t)) == 0);
		unitigIsForkStem[std::get<0>(t)] = i;
		assert(unitigBelongsToForkBranch.count(std::get<1>(t)) == 0);
		unitigBelongsToForkBranch[std::get<1>(t)] = i;
		assert(unitigBelongsToForkBranch.count(std::get<2>(t)) == 0);
		unitigBelongsToForkBranch[std::get<2>(t)] = i;
	}
	phmap::flat_hash_set<uint64_t> longnodeTips;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < minLongUnitigLength) continue;
		if (graph.coverages[i] < 0.5 * averageCoverage) continue;
		if (graph.coverages[i] > 1.5 * averageCoverage) continue;
		if (graph.edges.getEdges(std::make_pair(i, true)).size() == 0)
		{
			longnodeTips.insert(i + firstBitUint64_t);
		}
		if (graph.edges.getEdges(std::make_pair(i, false)).size() == 0)
		{
			longnodeTips.insert(i);
		}
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, size_t>> forkBranchLongNodeConnections;
	phmap::flat_hash_map<uint64_t, size_t> forkBranchLongestHomozygousChunks;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<uint64_t> longNodes;
		std::vector<std::pair<size_t, size_t>> unitigChunksInLongnode;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				assert(readPaths[i][j].path[k] & firstBitUint64_t);
				size_t unitig = readPaths[i][j].path[k] & maskUint64_t;
				assert(unitig < graph.unitigLengths.size());
				if (graph.unitigLengths[unitig] < minLongUnitigLength) continue;
				if (graph.coverages[unitig] < 0.5 * averageCoverage) continue;
				if (graph.coverages[unitig] > 1.5 * averageCoverage) continue;
				longNodes.emplace_back(readPaths[i][j].path[k]);
				unitigChunksInLongnode.emplace_back(0, graph.chunksInUnitig[unitig].size()-1);
				if (k == 0 && readPaths[i][j].pathLeftClipChunks > 0)
				{
					unitigChunksInLongnode.back().first = readPaths[i][j].pathLeftClipChunks;
				}
				if (k+1 == readPaths[i][j].path.size() && readPaths[i][j].pathRightClipChunks > 0)
				{
					unitigChunksInLongnode.back().second -= readPaths[i][j].pathRightClipChunks;
				}
				assert(unitigChunksInLongnode.back().second >= unitigChunksInLongnode.back().first);
				assert(unitigChunksInLongnode.back().second < graph.chunksInUnitig[unitig].size());
			}
		}
		for (size_t j = 0; j+1 < longNodes.size(); j++)
		{
			if (longnodeTips.count(longNodes[j]) == 0) continue;
			if (unitigBelongsToForkBranch.count(longNodes[j+1] ^ firstBitUint64_t) == 1)
			{
				forkBranchLongNodeConnections[longNodes[j+1] ^ firstBitUint64_t][longNodes[j] ^ firstBitUint64_t] += 1;
			}
			else if (j+2 < longNodes.size() && unitigIsForkStem.count(longNodes[j+1]) == 1 && unitigBelongsToForkBranch.count(longNodes[j+2] ^ firstBitUint64_t) == 1 && unitigIsForkStem.at(longNodes[j+1]) == unitigBelongsToForkBranch.at(longNodes[j+2] ^ firstBitUint64_t))
			{
				forkBranchLongNodeConnections[longNodes[j+2] ^ firstBitUint64_t][longNodes[j] ^ firstBitUint64_t] += 1;
			}
		}
		for (size_t j = 1; j < longNodes.size(); j++)
		{
			if (longnodeTips.count(longNodes[j] ^ firstBitUint64_t) == 0) continue;
			if (unitigBelongsToForkBranch.count(longNodes[j-1]) == 1)
			{
				forkBranchLongNodeConnections[longNodes[j-1]][longNodes[j]] += 1;
			}
			else if (j >= 2 && unitigIsForkStem.count(longNodes[j-1] ^ firstBitUint64_t) == 1 && unitigBelongsToForkBranch.count(longNodes[j-2]) == 1 && unitigIsForkStem.at(longNodes[j-1] ^ firstBitUint64_t) == unitigBelongsToForkBranch.at(longNodes[j-2]))
			{
				forkBranchLongNodeConnections[longNodes[j-2]][longNodes[j]] += 1;
			}
		}
		for (size_t j = 0; j+1 < longNodes.size(); j++)
		{
			if (unitigBelongsToForkBranch.count(longNodes[j]) == 0) continue;
			if (unitigIsForkStem.count(longNodes[j+1] ^ firstBitUint64_t) == 0) continue;
			if (unitigBelongsToForkBranch.at(longNodes[j]) != unitigIsForkStem.at(longNodes[j+1] ^ firstBitUint64_t)) continue;
			forkBranchLongestHomozygousChunks[longNodes[j]] = std::max(forkBranchLongestHomozygousChunks[longNodes[j]], unitigChunksInLongnode[j+1].second + 1);
		}
		for (size_t j = 1; j < longNodes.size(); j++)
		{
			if (unitigBelongsToForkBranch.count(longNodes[j] ^ firstBitUint64_t) == 0) continue;
			if (unitigIsForkStem.count(longNodes[j-1]) == 0) continue;
			if (unitigBelongsToForkBranch.at(longNodes[j] ^ firstBitUint64_t) != unitigIsForkStem.at(longNodes[j-1])) continue;
			forkBranchLongestHomozygousChunks[longNodes[j] ^ firstBitUint64_t] = std::max(forkBranchLongestHomozygousChunks[longNodes[j] ^ firstBitUint64_t], graph.chunksInUnitig[longNodes[j-1] & maskUint64_t].size() - unitigChunksInLongnode[j-1].first);
		}
	}
	phmap::flat_hash_map<size_t, std::pair<size_t, size_t>> chunkImpliesForkAllele;
	phmap::flat_hash_map<size_t, size_t> chunkIsForkHomozygous;
	phmap::flat_hash_map<size_t, size_t> chunkNewAlleleNum;
	size_t nextNum = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			nextNum = std::max(nextNum, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
		}
	}
	nextNum += 1;
	for (size_t i = 0; i < yForks.size(); i++)
	{
		uint64_t intactBranch = std::numeric_limits<uint64_t>::max();
		uint64_t brokenBranch = std::numeric_limits<uint64_t>::max();
		assert(forkBranchLongestHomozygousChunks.count(std::get<1>(yForks[i])) == 1);
		assert(forkBranchLongestHomozygousChunks.count(std::get<2>(yForks[i])) == 1);
		if (forkBranchLongestHomozygousChunks.at(std::get<1>(yForks[i])) > forkBranchLongestHomozygousChunks.at(std::get<2>(yForks[i])))
		{
			intactBranch = std::get<1>(yForks[i]);
			brokenBranch = std::get<2>(yForks[i]);
		}
		if (forkBranchLongestHomozygousChunks.at(std::get<1>(yForks[i])) < forkBranchLongestHomozygousChunks.at(std::get<2>(yForks[i])))
		{
			intactBranch = std::get<2>(yForks[i]);
			brokenBranch = std::get<1>(yForks[i]);
		}
		if (intactBranch == std::numeric_limits<uint64_t>::max())
		{
			// both branches extend equally deep into stem node, can't pick which one was broken
			continue;
		}
		phmap::flat_hash_set<uint64_t> otherEndCandidates;
		if (forkBranchLongNodeConnections.count(brokenBranch) == 0) continue;
		if (forkBranchLongNodeConnections.count(intactBranch) == 1)
		{
			bool allIntactConnectionsSmall = false;
			for (auto pair : forkBranchLongNodeConnections.at(intactBranch))
			{
				if (pair.second != 1)
				{
					allIntactConnectionsSmall = false;
					break;
				}
			}
			if (allIntactConnectionsSmall) continue;
		}
		for (auto pair : forkBranchLongNodeConnections.at(brokenBranch))
		{
			if (pair.second == 1) continue;
			otherEndCandidates.insert(pair.first);
		}
		if (otherEndCandidates.size() != 1) continue;
		size_t otherEnd = *otherEndCandidates.begin();
		size_t homozygousCount = forkBranchLongestHomozygousChunks.at(brokenBranch);
		assert(homozygousCount >= 1);
		assert(homozygousCount < graph.chunksInUnitig[std::get<0>(yForks[i])].size());
		std::cerr << "mend stem " << ((std::get<0>(yForks[i]) & firstBitUint64_t) ? ">" : "<") << (std::get<0>(yForks[i]) & maskUint64_t) << " intact branch " << ((intactBranch & firstBitUint64_t) ? ">" : "<") << (intactBranch & maskUint64_t) << " broken branch " << ((brokenBranch & firstBitUint64_t) ? ">" : "<") << (brokenBranch & maskUint64_t) << " to " << ((otherEnd & firstBitUint64_t) ? ">" : "<") << (otherEnd & maskUint64_t) << " homozygous chunk count " << homozygousCount << std::endl;
		if (std::get<0>(yForks[i]) & firstBitUint64_t)
		{
			for (size_t j = 0; j+homozygousCount < graph.chunksInUnitig[std::get<0>(yForks[i]) & maskUint64_t].size(); j++)
			{
				chunkImpliesForkAllele[graph.chunksInUnitig[std::get<0>(yForks[i])][j] & maskUint64_t] = std::make_pair(i, 0);
			}
			for (size_t j = graph.chunksInUnitig[std::get<0>(yForks[i]) & maskUint64_t].size() - homozygousCount; j < graph.chunksInUnitig[std::get<0>(yForks[i]) & maskUint64_t].size(); j++)
			{
				chunkIsForkHomozygous[graph.chunksInUnitig[std::get<0>(yForks[i])][j] & maskUint64_t] = i;
				chunkNewAlleleNum[graph.chunksInUnitig[std::get<0>(yForks[i])][j] & maskUint64_t] = nextNum;
				nextNum += 2;
			}
		}
		else
		{
			for (size_t j = homozygousCount; j < graph.chunksInUnitig[std::get<0>(yForks[i]) & maskUint64_t].size(); j++)
			{
				chunkImpliesForkAllele[graph.chunksInUnitig[std::get<0>(yForks[i])][j] & maskUint64_t] = std::make_pair(i, 0);
			}
			for (size_t j = 0; j < homozygousCount; j++)
			{
				chunkIsForkHomozygous[graph.chunksInUnitig[std::get<0>(yForks[i])][j] & maskUint64_t] = i;
				chunkNewAlleleNum[graph.chunksInUnitig[std::get<0>(yForks[i])][j] & maskUint64_t] = nextNum;
				nextNum += 2;
			}
		}
		for (uint64_t chunk : graph.chunksInUnitig[intactBranch & maskUint64_t])
		{
			chunkImpliesForkAllele[chunk & maskUint64_t] = std::make_pair(i, 0);
		}
		for (uint64_t chunk : graph.chunksInUnitig[brokenBranch & maskUint64_t])
		{
			chunkImpliesForkAllele[chunk & maskUint64_t] = std::make_pair(i, 1);
		}
		for (uint64_t chunk : graph.chunksInUnitig[otherEnd & maskUint64_t])
		{
			chunkImpliesForkAllele[chunk & maskUint64_t] = std::make_pair(i, 1);
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		phmap::flat_hash_map<size_t, std::pair<size_t, size_t>> forkAlleleCounts;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (chunkImpliesForkAllele.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			std::pair<size_t, size_t> forkAndAllele = chunkImpliesForkAllele.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			assert(forkAndAllele.second == 0 || forkAndAllele.second == 1);
			if (forkAndAllele.second == 0)
			{
				forkAlleleCounts[forkAndAllele.first].first += 1;
			}
			if (forkAndAllele.second == 1)
			{
				forkAlleleCounts[forkAndAllele.first].second += 1;
			}
		}
		phmap::flat_hash_map<size_t, size_t> forkAlleleAssignments;
		for (auto pair : forkAlleleCounts)
		{
			if (pair.second.first > pair.second.second)
			{
				forkAlleleAssignments[pair.first] = 0;
			}
			if (pair.second.first < pair.second.second)
			{
				forkAlleleAssignments[pair.first] = 1;
			}
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (chunkIsForkHomozygous.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			size_t fork = chunkIsForkHomozygous.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			if (forkAlleleAssignments.count(fork) == 0)
			{
				std::get<2>(chunksPerRead[i][j]) = std::numeric_limits<uint64_t>::max();
				continue;
			}
			size_t allele = forkAlleleAssignments.at(fork);
			assert(allele == 0 || allele == 1);
			std::get<2>(chunksPerRead[i][j]) = chunkNewAlleleNum[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] + allele + firstBitUint64_t;
		}
	}
}
