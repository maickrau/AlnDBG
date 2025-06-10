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
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		std::vector<size_t> lengthCounts;
		phmap::flat_hash_map<std::string, size_t> distinctSequences;
		std::vector<size_t> indicesWithDistinctSequences;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
			assert(j == 0 || sequences[j].size() >= sequences[j-1].size());
			if (j > 0 && sequences.back().size() > sequences[j-1].size())
			{
				distinctSequences.clear();
			}
			if (distinctSequences.count(sequences.back()) == 1)
			{
				merge(parent, j, distinctSequences.at(sequences.back()));
				continue;
			}
			while (lengthCounts.size() <= sequences.back().size()) lengthCounts.emplace_back(0);
			lengthCounts[sequences.back().size()] += 1;
			indicesWithDistinctSequences.emplace_back(j);
			distinctSequences[sequences.back()] = j;
		}
		{
			decltype(distinctSequences) tmp;
			std::swap(tmp, distinctSequences);
		}
		size_t maxRangeSize = 0;
		for (size_t i = 21; i < lengthCounts.size(); i++)
		{
			size_t counts = 0;
			for (size_t j = i-1; j >= i-std::max<size_t>(i*mismatchFraction, mismatchFloor) && j < lengthCounts.size(); j--)
			{
				counts += lengthCounts[j];
			}
			maxRangeSize = std::max(maxRangeSize, counts);
		}
		std::vector<std::pair<size_t, size_t>> ranges;
		if (maxRangeSize > 1000)
		{
			std::vector<std::pair<size_t, size_t>> lengthRanges;
			size_t currentRangeStart = 0;
			size_t currentBlockMinimum = std::numeric_limits<size_t>::max();
			size_t currentBlockMinimumIndex = 0;
			bool hadAnySolidsSoFar = false;
			for (size_t i = 21; i < lengthCounts.size(); i++)
			{
				size_t counts = 0;
				for (size_t j = i-1; j >= i-std::max<size_t>(i*mismatchFraction, mismatchFloor) && j < lengthCounts.size(); j--)
				{
					counts += lengthCounts[j];
				}
				if (counts >= maxRangeSize/1000)
				{
					hadAnySolidsSoFar = true;
				}
				if (!hadAnySolidsSoFar) continue;
				if (counts >= maxRangeSize/1000)
				{
					if (currentBlockMinimum != std::numeric_limits<size_t>::max())
					{
						lengthRanges.emplace_back(currentRangeStart, currentBlockMinimumIndex);
						currentRangeStart = currentBlockMinimumIndex - std::max<size_t>(currentBlockMinimumIndex * mismatchFraction, mismatchFloor);
					}
					currentBlockMinimum = std::numeric_limits<size_t>::max();
				}
				else
				{
					if (counts < currentBlockMinimum)
					{
						currentBlockMinimum = counts;
						currentBlockMinimumIndex = i;
					}
				}
			}
			lengthRanges.emplace_back(currentRangeStart, sequences.back().size());
			ranges.resize(lengthRanges.size(), std::make_pair(std::numeric_limits<size_t>::max(), 0));
			for (size_t i = 0; i < indicesWithDistinctSequences.size(); i++)
			{
				for (size_t j = 0; j < lengthRanges.size(); j++)
				{
					if (sequences[indicesWithDistinctSequences[i]].size() < lengthRanges[j].first) continue;
					if (sequences[indicesWithDistinctSequences[i]].size() > lengthRanges[j].second) continue;
					ranges[j].first = std::min<size_t>(ranges[j].first, i);
					ranges[j].second = std::max<size_t>(ranges[j].second, i+1);
				}
			}
		}
		else
		{
			ranges.emplace_back(0, indicesWithDistinctSequences.size());
		}
		for (auto range : ranges)
		{
			assert(range.second > range.first);
			longestLength = sequences[indicesWithDistinctSequences[range.second-1]].size();
			std::vector<size_t> parentOfDistinctSequences = getFastTransitiveClosureMultithread(range.second-range.first, std::max((size_t)(longestLength * mismatchFraction), mismatchFloor), longestLength, numThreads, [&sequences, &indicesWithDistinctSequences, range](const size_t i, const size_t j, const size_t maxDist) { return getNumMismatches(sequences[indicesWithDistinctSequences[range.first+i]], sequences[indicesWithDistinctSequences[range.first+j]], maxDist); }, [&sequences, mismatchFloor, &indicesWithDistinctSequences, range](const size_t i, const size_t j) { return std::max((size_t)(std::min(sequences[indicesWithDistinctSequences[range.first+i]].size(), sequences[indicesWithDistinctSequences[range.first+j]].size()) * mismatchFraction), mismatchFloor); });
			for (size_t j = 0; j < parentOfDistinctSequences.size(); j++)
			{
				merge(parent, indicesWithDistinctSequences[range.first+j], indicesWithDistinctSequences[range.first+find(parentOfDistinctSequences, j)]);
			}
		}
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
