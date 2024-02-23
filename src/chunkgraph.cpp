#include <mutex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "UnionFind.h"
#include "CommonUtils.h"
#include "ReadStorage.h"
#include "MatchIndex.h"
#include "KmerGraph.h"
#include "KmerMatcher.h"
#include "MatchGroup.h"
#include "UnitigGraph.h"
#include "GraphCleaner.h"
#include "AlnHaploFilter.h"
#include "GraphPhaser.h"
#include "GraphResolver.h"
#include "AnchorFinder.h"
#include "ChunkmerFilter.h"

void writeGraph(const std::string& filename, const std::vector<std::vector<size_t>>& lengths, const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeCoverage, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps)
{
	std::ofstream graph { filename };
	for (size_t i = 0; i < coverages.size(); i++)
	{
		assert(lengths[i].size() >= coverages[i]);
		if (coverages[i] < 2) continue;
		size_t length = lengths[i][lengths[i].size()/2];
		size_t coverage = coverages[i];
		graph << "S\t" << i << "\t*\tLN:i:" << length << "\tll:i:" << coverage << "\tFC:i:" << length*coverage << std::endl;
	}
	for (auto pair : edgeCoverage)
	{
		if (pair.second < 2) continue;
		assert(edgeOverlaps.count(pair.first) == 1);
		graph << "L\t" << (pair.first.first & maskUint64_t) << "\t" << ((pair.first.first & firstBitUint64_t) ? "+" : "-") << "\t" << (pair.first.second & maskUint64_t) << "\t" << ((pair.first.second & firstBitUint64_t) ? "+" : "-") << "\t" << edgeOverlaps.at(pair.first) << "M\tec:i:" << pair.second << std::endl;
	}
}

std::pair<std::vector<std::vector<size_t>>, std::vector<size_t>> getLengthsAndCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if ((std::get<2>(chunksPerRead[i][j]) & maskUint64_t) >= coverages.size())
			{
				coverages.resize((std::get<2>(chunksPerRead[i][j]) & maskUint64_t)+1, 0);
				lengths.resize((std::get<2>(chunksPerRead[i][j]) & maskUint64_t)+1);
			}
			if (j > 1 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			lengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].emplace_back(std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]));
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (j > 1 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			coverages[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] += 1;
		}
	}
	return std::make_pair(lengths, coverages);
}

phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getEdgeOverlaps(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<size_t>> overlaps;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			auto prev = std::get<2>(chunksPerRead[i][j-1]);
			auto curr = std::get<2>(chunksPerRead[i][j]);
			auto pairkey = MBG::canon(std::make_pair(prev & maskUint64_t, prev & firstBitUint64_t), std::make_pair(curr & maskUint64_t, curr & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
			size_t overlap = 0;
			if (std::get<1>(chunksPerRead[i][j-1]) > std::get<0>(chunksPerRead[i][j])) overlap = std::get<1>(chunksPerRead[i][j-1]) - std::get<0>(chunksPerRead[i][j]);
			overlaps[key].push_back(overlap);
		}
	}
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> result;
	for (auto& pair : overlaps)
	{
		std::sort(pair.second.begin(), pair.second.end());
		result[pair.first] = pair.second[pair.second.size()/2];
	}
	return result;
}

phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getEdgeCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			auto prev = std::get<2>(chunksPerRead[i][j-1]);
			auto curr = std::get<2>(chunksPerRead[i][j]);
			auto pairkey = MBG::canon(std::make_pair(prev & maskUint64_t, prev & firstBitUint64_t), std::make_pair(curr & maskUint64_t, curr & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
			edgeCoverage[key] += 1;
		}
	}
	return edgeCoverage;
}

void removeContainedChunks(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, std::vector<bool>& useTheseChunks)
{
	phmap::flat_hash_set<size_t> contained;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		size_t coveredUntil = 0;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (std::get<1>(chunksPerRead[i][j]) < coveredUntil)
			{
				contained.insert(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			else
			{
				for (size_t k = 0; k < j; k++)
				{
					if (std::get<0>(chunksPerRead[i][k]) == std::get<0>(chunksPerRead[i][j])) break;
					if (std::get<1>(chunksPerRead[i][k]) >= std::get<1>(chunksPerRead[i][j]))
					{
						contained.insert(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
						break;
					}
				}
			}
			coveredUntil = std::max(coveredUntil, std::get<1>(chunksPerRead[i][j]));
		}
	}
	for (size_t i = 0; i < useTheseChunks.size(); i++)
	{
		if (!useTheseChunks[i]) continue;
		if (contained.count(i) == 1) useTheseChunks[i] = false;
	}
}

void mergePerLocation(const MatchIndex& matchIndex, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::cerr << "merging per location" << std::endl;
	std::vector<std::pair<size_t, bool>> parent = getParent(chunksPerRead);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto found = find(parent, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			uint64_t result = found.first + (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t);
			if (!found.second) result ^= firstBitUint64_t;
			std::get<2>(chunksPerRead[i][j]) = result;
		}
	}
}

void splitPerLength(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::cerr << "splitting per length" << std::endl;
	const double differenceFraction = 0.05;
	const size_t differenceConstant = 50;
	phmap::flat_hash_map<size_t, std::vector<size_t>> lengthsPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			lengthsPerChunk[std::get<2>(t) & maskUint64_t].emplace_back(std::get<1>(t) - std::get<0>(t));
		}
	}
	phmap::flat_hash_map<size_t, std::vector<size_t>> splitters;
	for (auto& pair : lengthsPerChunk)
	{
		std::sort(pair.second.begin(), pair.second.end());
		splitters[pair.first].emplace_back(0);
		for (size_t i = 1; i < pair.second.size(); i++)
		{
			size_t distance = pair.second[i] - pair.second[i-1];
			if (distance < std::max((size_t)(pair.second[i] * differenceFraction), (size_t)differenceConstant)) continue;
			splitters[pair.first].emplace_back((pair.second[i] + pair.second[i-1])/2);
		}
	}
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> splitterToNode;
	size_t nextNum = 0;
	for (const auto& pair : splitters)
	{
		for (const size_t len : pair.second)
		{
			auto key = std::make_pair(pair.first, len);
			assert(splitterToNode.count(key) == 0);
			assert(nextNum == splitterToNode.size());
			splitterToNode[key] = nextNum;
			nextNum += 1;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto& t : chunksPerRead[i])
		{
			assert(std::get<1>(t) > std::get<0>(t));
			size_t distance = std::get<1>(t) - std::get<0>(t);
			size_t newNode = std::numeric_limits<size_t>::max();
			size_t dist = std::numeric_limits<size_t>::max();
			for (size_t len : splitters[std::get<2>(t) & maskUint64_t])
			{
				if (len >= distance) break;
				dist = len;
			}
			assert(dist != std::numeric_limits<size_t>::max());
			assert(dist < distance);
			auto key = std::make_pair(std::get<2>(t) & maskUint64_t, dist);
			assert(splitterToNode.count(key) == 1);
			std::get<2>(t) = (std::get<2>(t) & firstBitUint64_t) + splitterToNode.at(key);
		}
	}
}

template <typename F1, typename F2>
size_t getNumMismatches(F1 leftSequenceGetter, F2 rightSequenceGetter, const size_t leftSequenceLength, const size_t rightSequenceLength, const size_t maxMismatches)
{
	assert(leftSequenceLength >= rightSequenceLength);
	assert(rightSequenceLength + maxMismatches >= leftSequenceLength);
	// left sequence: left to right
	// right sequence: up to down
	// offset: diagonally down-right
	// diagonal: left-right
	std::vector<size_t> offset;
	size_t maxBefore = (leftSequenceLength - rightSequenceLength) + (maxMismatches - (leftSequenceLength - rightSequenceLength))/2;
	size_t maxAfter = (maxMismatches - (leftSequenceLength - rightSequenceLength))/2;
	assert(maxBefore < maxMismatches);
	assert(maxAfter < maxMismatches);
	offset.resize(maxBefore + maxAfter + 1, std::numeric_limits<size_t>::max());
	offset[maxBefore] = 0;
	while (offset[maxBefore] < leftSequenceLength && offset[maxBefore] < rightSequenceLength && leftSequenceGetter(offset[maxBefore]) == rightSequenceGetter(offset[maxBefore])) offset[maxBefore] += 1;
	if (leftSequenceLength == rightSequenceLength && offset[maxBefore] == leftSequenceLength) return 0;
	// zero diagonal (top left corner) at maxBefore
	for (size_t score = 1; score < maxMismatches; score++)
	{
		if (score <= maxBefore) offset[maxBefore-score] = score;
		if (maxBefore+score < offset.size()) offset[maxBefore+score] = score;
		for (size_t diagonal = (score < maxBefore ? (maxBefore - score) : 0); diagonal <= maxBefore + score && diagonal < offset.size(); diagonal++)
		{
			assert(offset[diagonal] != std::numeric_limits<size_t>::max());
			offset[diagonal] += 1;
			if (diagonal > 0 && offset[diagonal-1] != std::numeric_limits<size_t>::max())
			{
				offset[diagonal] = std::max(offset[diagonal], offset[diagonal-1]);
			}
			if (diagonal+1 < offset.size() && offset[diagonal+1] != std::numeric_limits<size_t>::max())
			{
				offset[diagonal] = std::max(offset[diagonal], offset[diagonal+1]+1);
			}
			offset[diagonal] = std::min(offset[diagonal], std::min(leftSequenceLength, rightSequenceLength+maxBefore-diagonal));
		}
		for (size_t diagonal = (score < maxBefore ? (maxBefore - score) : 0); diagonal <= maxBefore + score && diagonal < offset.size(); diagonal++)
		{
			assert(offset[diagonal] <= leftSequenceLength);
			assert(offset[diagonal]+diagonal-maxBefore <= rightSequenceLength);
			while (offset[diagonal] < leftSequenceLength && offset[diagonal]+diagonal-maxBefore < rightSequenceLength && leftSequenceGetter(offset[diagonal]) == rightSequenceGetter(offset[diagonal]+diagonal-maxBefore)) offset[diagonal] += 1;
			if (offset[diagonal] == leftSequenceLength && offset[diagonal]+diagonal-maxBefore == rightSequenceLength) return score;
		}
	}
	return maxMismatches+1;
}

size_t getNumMismatches(const std::vector<TwobitString>& readSequences, const size_t leftRead, const size_t rightRead, const std::tuple<size_t, size_t, uint64_t> left, const std::tuple<size_t, size_t, uint64_t> right, const size_t maxMismatches)
{
	assert((std::get<2>(left) & maskUint64_t) == (std::get<2>(right) & maskUint64_t));
	size_t leftLen = std::get<1>(left) - std::get<0>(left);
	size_t rightLen = std::get<1>(right) - std::get<0>(right);
	if (rightLen > leftLen) return getNumMismatches(readSequences, rightRead, leftRead, right, left, maxMismatches);
	assert(leftLen >= rightLen);
	if (leftLen > rightLen + maxMismatches) return maxMismatches+1;
	if (std::get<2>(left) & firstBitUint64_t)
	{
		if (std::get<2>(right) & firstBitUint64_t)
		{
			return getNumMismatches([&readSequences, left, leftRead](size_t pos){return readSequences[leftRead].get(std::get<0>(left)+pos);}, [&readSequences, right, rightRead](size_t pos){return readSequences[rightRead].get(std::get<0>(right)+pos);}, std::get<1>(left)-std::get<0>(left)+1, std::get<1>(right)-std::get<0>(right)+1, maxMismatches);
		}
		else
		{
			return getNumMismatches([&readSequences, left, leftRead](size_t pos){return readSequences[leftRead].get(std::get<0>(left)+pos);}, [&readSequences, right, rightRead](size_t pos){return 3-readSequences[rightRead].get(std::get<1>(right)-pos);}, std::get<1>(left)-std::get<0>(left)+1, std::get<1>(right)-std::get<0>(right)+1, maxMismatches);
		}
	}
	else
	{
		if (std::get<2>(right) & firstBitUint64_t)
		{
			return getNumMismatches([&readSequences, left, leftRead](size_t pos){return 3-readSequences[leftRead].get(std::get<1>(left)-pos);}, [&readSequences, right, rightRead](size_t pos){return readSequences[rightRead].get(std::get<0>(right)+pos);}, std::get<1>(left)-std::get<0>(left)+1, std::get<1>(right)-std::get<0>(right)+1, maxMismatches);
		}
		else
		{
			return getNumMismatches([&readSequences, left, leftRead](size_t pos){return 3-readSequences[leftRead].get(std::get<1>(left)-pos);}, [&readSequences, right, rightRead](size_t pos){return 3-readSequences[rightRead].get(std::get<1>(right)-pos);}, std::get<1>(left)-std::get<0>(left)+1, std::get<1>(right)-std::get<0>(right)+1, maxMismatches);
		}
	}
}

void splitPerSequenceIdentity(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::cerr << "splitting per sequence identity" << std::endl;
	const double mismatchFraction = 0.05;
	const size_t mismatchFloor = 10;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		std::cerr << "split chunk " << i << " coverage " << occurrencesPerChunk[i].size() << std::endl;
		std::sort(occurrencesPerChunk[i].begin(), occurrencesPerChunk[i].end(), [&chunksPerRead](auto left, auto right) { return std::get<1>(chunksPerRead[left.first][left.second])-std::get<0>(chunksPerRead[left.first][left.second]) < std::get<1>(chunksPerRead[right.first][right.second])-std::get<0>(chunksPerRead[right.first][right.second]); });
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
		{
			for (size_t k = j-1; k < occurrencesPerChunk[i].size(); k--)
			{
				if (find(parent, k) == find(parent, j)) continue;
				auto left = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				auto right = chunksPerRead[occurrencesPerChunk[i][k].first][occurrencesPerChunk[i][k].second];
				assert(std::get<1>(left) - std::get<0>(left) >= std::get<1>(right) - std::get<0>(right));
				size_t maxMismatches = std::max(mismatchFloor, (size_t)((std::get<1>(left) - std::get<0>(left))*mismatchFraction));
				if ((std::get<1>(left) - std::get<0>(left)) - (std::get<1>(right) - std::get<0>(right)) >= maxMismatches) break;
				size_t mismatches = getNumMismatches(readSequences, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][k].first, left, right, maxMismatches);
				if (mismatches > maxMismatches) continue;
				merge(parent, j, k);
			}
		}
		phmap::flat_hash_map<size_t, size_t> keyToNode;
		for (size_t j = 0; j < parent.size(); j++)
		{
			size_t key = find(parent, j);
			if (keyToNode.count(key) == 1) continue;
			keyToNode[key] = nextNum;
			nextNum += 1;
		}
		std::cerr << "splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
		}
	}
}

void makeGraph(const MatchIndex& matchIndex, const std::vector<size_t>& rawReadLengths, const std::vector<TwobitString>& readSequences)
{
	std::vector<bool> useTheseChunks;
	useTheseChunks.resize(matchIndex.numWindowChunks() - matchIndex.numUniqueChunks(), true);
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead = getChunksPerRead(matchIndex, rawReadLengths, useTheseChunks);
	mergePerLocation(matchIndex, chunksPerRead);
	splitPerLength(chunksPerRead);
	{
		std::cerr << "writing first graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round1.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
	splitPerSequenceIdentity(readSequences, chunksPerRead);
	{
		std::cerr << "writing second graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round2.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths2.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
	/*
	for (size_t i = 0; i < useTheseChunks.size(); i++)
	{
		if (coverages[find(parent, i).first] >= 40)
		{
			useTheseChunks[i] = false;
		}
	}
	phmap::flat_hash_map<size_t, size_t> countRepetitive;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		phmap::flat_hash_map<uint64_t, size_t> lastPos;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			std::pair<size_t, bool> pairkey = find(parent, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			if ((std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t) ^ firstBitUint64_t) pairkey.second = !pairkey.second;
			uint64_t key = pairkey.first + (pairkey.second ? firstBitUint64_t : 0);
			if (lastPos.count(key) == 1)
			{
				if (lastPos.at(key) + 20000 > std::get<0>(chunksPerRead[i][j]) && lastPos.at(key) != std::get<1>(chunksPerRead[i][j]))
				{
					countRepetitive[key & maskUint64_t] += 1;
				}
			}
			lastPos[key] = std::get<1>(chunksPerRead[i][j]);
		}
	}
	for (size_t i = 0; i < useTheseChunks.size(); i++)
	{
		if (countRepetitive[find(parent, i).first] >= 5)
		{
			useTheseChunks[i] = false;
		}
	}
	chunksPerRead = getChunksPerRead(matchIndex, rawReadLengths, useTheseChunks);
	parent = getParent(chunksPerRead);
	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead, parent);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead, parent);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead, parent);
	writeGraph("graph.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);*/
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	std::vector<std::string> readFiles;
	for (size_t i = 5; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	std::vector<TwobitString> readSequences;
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		std::mutex sequenceMutex;
		storage.iterateReadsFromFile(file, numThreads, false, [&matchIndex, &indexMutex, &sequenceMutex, &readSequences](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
			{
				std::lock_guard<std::mutex> lock { sequenceMutex };
				while (readSequences.size() <= readName) readSequences.emplace_back();
				readSequences[readName] = sequence;
			}
		});
	}
	std::cerr << readSequences.size() << " reads" << std::endl;
	std::cerr << matchIndex.numWindowChunks() << " distinct windowchunks" << std::endl;
	std::cerr << matchIndex.numUniqueChunks() << " windowchunks have only one read" << std::endl;
	matchIndex.clearConstructionVariablesAndCompact();
	const std::vector<size_t>& readBasepairLengths = storage.getRawReadLengths();
	const std::vector<std::string>& readNames = storage.getNames();
	for (size_t i = 0; i < readNames.size(); i++)
	{
		std::cerr << "readname " << i << " " << readNames[i] << std::endl;
	}
	makeGraph(matchIndex, readBasepairLengths, readSequences);
}
