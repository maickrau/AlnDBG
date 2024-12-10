#include <iostream>
#include <cassert>
#include <thread>
#include "SequenceIdentitySplitter.h"
#include "ChunkHelper.h"
#include "Common.h"
#include "UnionFind.h"
#include "EdlibWrapper.h"
#include "SequenceHelper.h"
#include "TransitiveClosure.h"
#include "CanonHelper.h"

void splitPerSequenceIdentityRoughly(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by sequence identity" << std::endl;
	const size_t mismatchFloor = 10;
	size_t nextNum = 0;
	size_t totalSplitted = 0;
	size_t totalSplittedTo = 0;
	std::mutex resultMutex;
	auto oldChunks = chunksPerRead;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	std::vector<size_t> iterationOrder;
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		if (!canonical[i]) continue;
		iterationOrder.emplace_back(i);
	}
	std::sort(iterationOrder.begin(), iterationOrder.end(), [&occurrencesPerChunk](size_t left, size_t right) { return occurrencesPerChunk[left].size() > occurrencesPerChunk[right].size(); });
	size_t firstSmall = iterationOrder.size();
	for (size_t i = 0; i < iterationOrder.size(); i++)
	{
		if (occurrencesPerChunk[iterationOrder[i]].size() > 1000) continue;
		firstSmall = i;
		break;
	}
	for (size_t iterationIndex = 0; iterationIndex < firstSmall; iterationIndex++)
	{
		const size_t i = iterationOrder[iterationIndex];
		std::cerr << "begin big chunk with coverage " << occurrencesPerChunk[i].size() << std::endl;
		auto startTime = getTime();
		bool allSmall = true;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto chunk = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			if (std::get<1>(chunk) - std::get<0>(chunk) < 2*kmerSize) continue;
			allSmall = false;
			break;
		}
		if (allSmall)
		{
			auto endTime = getTime();
			std::cerr << "skip short big chunk with coverage " << occurrencesPerChunk[i].size() << " time " << formatTime(startTime, endTime) << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			continue;
		}
		std::sort(occurrencesPerChunk[i].begin(), occurrencesPerChunk[i].end(), [&chunksPerRead](const auto left, const auto right)
		{
			auto tleft = chunksPerRead[left.first][left.second];
			auto tright = chunksPerRead[right.first][right.second];
			return (std::get<1>(tleft) - std::get<0>(tleft)) < (std::get<1>(tright) - std::get<0>(tright));
		});
		std::vector<std::string> sequences;
		sequences.reserve(occurrencesPerChunk[i].size());
		size_t longestLength = 0;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
			longestLength = std::max(longestLength, sequences.back().size());
		}
		std::vector<size_t> parent = getFastTransitiveClosureMultithread(sequences.size(), std::max((size_t)(longestLength * mismatchFraction), mismatchFloor), numThreads, [&sequences](const size_t i, const size_t j, const size_t maxDist) { return getNumMismatches(sequences[i], sequences[j], maxDist); }, [&sequences, mismatchFloor](const size_t i, const size_t j) { return std::max((size_t)(std::min(sequences[i].size(), sequences[j].size()) * mismatchFraction), mismatchFloor); });
		phmap::flat_hash_map<size_t, size_t> parentToCluster;
		for (size_t j = 0; j < parent.size(); j++)
		{
			find(parent, j);
		}
		size_t nextCluster = 0;
		for (size_t j = 0; j < parent.size(); j++)
		{
			if (parent[j] != j) continue;
			parentToCluster[j] = nextCluster;
			nextCluster += 1;
		}
		auto endTime = getTime();
		assert(nextCluster >= 1);
		{
			if (nextCluster > 1)
			{
				totalSplitted += 1;
				totalSplittedTo += nextCluster;
			}
			std::cerr << "split big chunk with coverage " << occurrencesPerChunk[i].size() << " into " << nextCluster << " chunks time " << formatTime(startTime, endTime) << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + parentToCluster.at(parent[j]) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += nextCluster;
		}
	}
	iterateMultithreaded(firstSmall, iterationOrder.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, &totalSplitted, &totalSplittedTo, &iterationOrder, &occurrencesPerChunk, mismatchFloor](const size_t iterationIndex)
	{
		const size_t i = iterationOrder[iterationIndex];
		if (occurrencesPerChunk[i].size() < 2)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			// std::cerr << "skip chunk with coverage " << occurrencesPerChunk[i].size() << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			std::cerr << "begin chunk with coverage " << occurrencesPerChunk[i].size() << std::endl;
		}
		auto startTime = getTime();
		std::sort(occurrencesPerChunk[i].begin(), occurrencesPerChunk[i].end(), [&chunksPerRead](const auto left, const auto right)
		{
			auto tleft = chunksPerRead[left.first][left.second];
			auto tright = chunksPerRead[right.first][right.second];
			return (std::get<1>(tleft) - std::get<0>(tleft)) < (std::get<1>(tright) - std::get<0>(tright));
		});
		std::vector<std::string> sequences;
		sequences.reserve(occurrencesPerChunk[i].size());
		size_t longestLength = 0;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
			longestLength = std::max(longestLength, sequences.back().size());
		}
		std::vector<size_t> parent = getFastTransitiveClosure(sequences.size(), std::max((size_t)(longestLength * mismatchFraction), mismatchFloor), [&sequences](const size_t i, const size_t j, const size_t maxDist) { return getNumMismatches(sequences[i], sequences[j], maxDist); }, [&sequences, mismatchFloor](const size_t i, const size_t j) { return std::max((size_t)(std::min(sequences[i].size(), sequences[j].size()) * mismatchFraction), mismatchFloor); });
		phmap::flat_hash_map<size_t, size_t> parentToCluster;
		for (size_t j = 0; j < parent.size(); j++)
		{
			find(parent, j);
		}
		size_t nextCluster = 0;
		for (size_t j = 0; j < parent.size(); j++)
		{
			if (parent[j] != j) continue;
			parentToCluster[j] = nextCluster;
			nextCluster += 1;
		}
		auto endTime = getTime();
		assert(nextCluster >= 1);
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			if (nextCluster > 1)
			{
				totalSplitted += 1;
				totalSplittedTo += nextCluster;
			}
			std::cerr << "split chunk with coverage " << occurrencesPerChunk[i].size() << " into " << nextCluster << " chunks time " << formatTime(startTime, endTime) << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + parentToCluster.at(parent[j]) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += nextCluster;
		}
	});
	std::cerr << "rough identity splitting splitted " << totalSplitted << " chunks to " << totalSplittedTo << " chunks" << std::endl;
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerSequenceIdentity(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by sequence identity" << std::endl;
	const size_t mismatchFloor = 10;
	auto oldChunks = chunksPerRead;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	std::vector<size_t> iterationOrder;
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		if (!canonical[i]) continue;
		iterationOrder.emplace_back(i);
	}
	std::sort(iterationOrder.begin(), iterationOrder.end(), [&occurrencesPerChunk](size_t left, size_t right) { return occurrencesPerChunk[left].size() > occurrencesPerChunk[right].size(); });
	size_t firstSmall = iterationOrder.size();
	for (size_t i = 0; i < iterationOrder.size(); i++)
	{
		if (occurrencesPerChunk[iterationOrder[i]].size() > 200 || occurrencesPerChunk[iterationOrder[i]].size() > 200*numThreads) continue;
		firstSmall = i;
		break;
	}
	size_t nextNum = 0;
	for (size_t iterationIndex = 0; iterationIndex < firstSmall; iterationIndex++)
	{
		const size_t i = iterationOrder[iterationIndex];
		// std::cerr << "try split chunk " << i << std::endl;
		// auto startTime = getTime();
		std::vector<std::pair<std::string, size_t>> sequencesPerOccurrence;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sequencesPerOccurrence.emplace_back();
			sequencesPerOccurrence.back().second = j;
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequencesPerOccurrence.back().first = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
		}
		std::sort(sequencesPerOccurrence.begin(), sequencesPerOccurrence.end(), [](const auto& left, const auto& right) {
			if (left.first.size() < right.first.size()) return true;
			if (left.first.size() > right.first.size()) return false;
			return left.first < right.first;
		});
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		const size_t blockSize = 50;
		const size_t numBlocks = (sequencesPerOccurrence.size()+blockSize-1)/blockSize;
		// std::cerr << "split part 1" << std::endl;
		iterateMultithreaded(0, numBlocks, numThreads, [&sequencesPerOccurrence, &parent, blockSize, i, mismatchFloor](size_t index)
		{
			for (size_t j = index * blockSize+1; j < (index+1) * blockSize && j < sequencesPerOccurrence.size(); j++)
			{
				for (size_t k = j-1; k < j && k >= index*blockSize; k--)
				{
					if (find(parent, k) == find(parent, j)) continue;
					assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
					size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
					if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches) break;
					size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
					if (mismatches > maxMismatches) continue;
					merge(parent, j, k);
				}
			}
		});
		for (size_t index = 1; index < numBlocks; index++)
		{
			size_t j = index * blockSize;
			size_t k = j-1;
			assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
			size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
			size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
			if (mismatches > maxMismatches) continue;
			merge(parent, j, k);
		}
		for (size_t i = 0; i < parent.size(); i++)
		{
			find(parent, i);
		}
		std::vector<size_t> activeBlocks;
		for (size_t i = 1; i < numBlocks; i++)
		{
			activeBlocks.emplace_back(i);
		}
		for (size_t blockDistance = 1; blockDistance < numBlocks; blockDistance++)
		{
			// std::cerr << "split part 2 dist " << blockDistance << std::endl;
			std::vector<std::pair<size_t, size_t>> merges;
			std::mutex mergeMutex;
			phmap::flat_hash_set<size_t> removeBlocks;
			iterateMultithreaded(0, activeBlocks.size(), numThreads, [&sequencesPerOccurrence, &parent, &merges, &mergeMutex, &activeBlocks, &removeBlocks, blockSize, blockDistance, i, mismatchFloor](size_t blockindex)
			{
				size_t index = activeBlocks[blockindex];
				phmap::flat_hash_map<size_t, size_t> extraParent;
				std::vector<std::pair<size_t, size_t>> mergesHere;
				for (size_t j = index * blockSize; j < (index+1) * blockSize && j < sequencesPerOccurrence.size(); j++)
				{
					for (size_t k = (index-blockDistance+1) * blockSize - 1; k < j && k >= (index-blockDistance) * blockSize; k--)
					{
						size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
						if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches)
						{
							if (j == index * blockSize && k == (index-blockDistance+1) * blockSize - 1)
							{
								assert(mergesHere.size() == 0);
								std::lock_guard<std::mutex> lock { mergeMutex };
								removeBlocks.emplace(index);
								return;
							}
							break;
						}
						if (parent[k] == parent[j]) continue;
						if (extraParent.count(parent[k]) == 0) extraParent[parent[k]] = parent[k];
						if (extraParent.count(parent[j]) == 0) extraParent[parent[j]] = parent[j];
						if (find(extraParent, parent[k]) == find(extraParent, parent[j])) continue;
						assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
						size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
						if (mismatches > maxMismatches) continue;
						merge(extraParent, parent[j], parent[k]);
						mergesHere.emplace_back(parent[j], parent[k]);
					}
				}
				{
					std::lock_guard<std::mutex> lock { mergeMutex };
					merges.insert(merges.end(), mergesHere.begin(), mergesHere.end());
				}
			});
			// std::cerr << merges.size() << " merges" << std::endl;
			for (auto pair : merges)
			{
				merge(parent, pair.first, pair.second);
			}
			if (merges.size() >= 1)
			{
				for (size_t i = 0; i < parent.size(); i++)
				{
					find(parent, i);
				}
			}
			removeBlocks.emplace(blockDistance);
			for (size_t i = activeBlocks.size()-1; i < activeBlocks.size(); i--)
			{
				if (removeBlocks.count(activeBlocks[i]) == 0) continue;
				std::swap(activeBlocks[i], activeBlocks.back());
				activeBlocks.pop_back();
			}
			// std::cerr << "remove " << removeBlocks.size() << ", " << activeBlocks.size() << " remaining" << std::endl;
			if (activeBlocks.size() == 0) break;
		}
		// auto endTime = getTime();
		phmap::flat_hash_map<size_t, size_t> keyToNode;
		for (size_t j = 0; j < parent.size(); j++)
		{
			size_t key = find(parent, j);
			if (keyToNode.count(key) == 1) continue;
			keyToNode[key] = nextNum;
			nextNum += 1;
		}
		// std::cerr << "sequence identity big chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " time " << formatTime(startTime, endTime) << " splitted to " << keyToNode.size() << std::endl;
		// std::cerr << "sequence identity splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			size_t index = sequencesPerOccurrence[j].second;
			std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
		}
	}
	std::mutex resultMutex;
	iterateMultithreaded(firstSmall, iterationOrder.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &sequenceIndex, &iterationOrder, &rawReadLengths, &canonical, mismatchFloor](const size_t iterationIndex)
	{
		const size_t i = iterationOrder[iterationIndex];
		if (!canonical[i]) return;
		// {
		// 	std::lock_guard<std::mutex> lock { resultMutex };
		// 	std::cerr << "sequence identity split chunk " << i << " coverage " << occurrencesPerChunk[i].size() << std::endl;
		// }
		// size_t coverage = occurrencesPerChunk[i].size();
		// size_t totalCheckablePairs = 0;
		// size_t totalMergedPairs = 0;
		// size_t totalNotMergedPairs = 0;
		// size_t totalDistinct = 0;
		// auto startTime = getTime();
		std::vector<std::pair<std::string, size_t>> sequencesPerOccurrence;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sequencesPerOccurrence.emplace_back();
			sequencesPerOccurrence.back().second = j;
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequencesPerOccurrence.back().first = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
		}
		std::sort(sequencesPerOccurrence.begin(), sequencesPerOccurrence.end(), [](const auto& left, const auto& right) {
			if (left.first.size() < right.first.size()) return true;
			if (left.first.size() > right.first.size()) return false;
			return left.first < right.first;
		});
		std::vector<bool> differentFromPredecessor;
		differentFromPredecessor.resize(occurrencesPerChunk[i].size(), false);
		differentFromPredecessor[0] = true;
		// totalDistinct += 1;
		for (size_t j = 1; j < differentFromPredecessor.size(); j++)
		{
			if (sequencesPerOccurrence[j].first == sequencesPerOccurrence[j-1].first) continue;
			differentFromPredecessor[j] = true;
			// totalDistinct += 1;
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
		{
			if (!differentFromPredecessor[j])
			{
				merge(parent, j, j-1);
				continue;
			}
			for (size_t k = j-1; k < occurrencesPerChunk[i].size(); k--)
			{
				if (!differentFromPredecessor[k]) continue;
				// totalCheckablePairs += 1;
				if (find(parent, k) == find(parent, j)) continue;
				assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
				size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
				if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches) break;
				size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
				if (mismatches > maxMismatches)
				{
					// totalNotMergedPairs += 1;
					continue;
				}
				// totalMergedPairs += 1;
				merge(parent, j, k);
			}
		}
		// auto endTime = getTime();
		phmap::flat_hash_map<size_t, size_t> keyToNode;
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
			// std::cerr << "sequence identity chunk " << i << " coverage " << coverage << " distinct " << totalDistinct << " checkable pairs " << totalCheckablePairs << " merged " << totalMergedPairs << " not merged " << totalNotMergedPairs << " time " << formatTime(startTime, endTime) << " splitted to " << keyToNode.size() << std::endl;
			// std::cerr << "sequence identity splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t index = sequencesPerOccurrence[j].second;
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}
