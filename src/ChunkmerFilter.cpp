#include <tuple>
#include "ChunkmerFilter.h"
#include "Common.h"
#include "UnionFind.h"

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getChunksPerRead(const MatchIndex& matchIndex, const std::vector<size_t>& rawReadLengths, const std::vector<bool>& useTheseChunks)
{
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead;
	chunksPerRead.resize(rawReadLengths.size());
	matchIndex.iterateChunks([&chunksPerRead, &rawReadLengths, &useTheseChunks](const size_t i, const ReadMatchposStorage& reads)
	{
		if (!useTheseChunks[i]) return;
		for (std::tuple<uint32_t, uint32_t, uint32_t> t : reads)
		{
			chunksPerRead[std::get<0>(t)].emplace_back(std::get<1>(t), std::get<2>(t), i);
		}
	});
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto& t : chunksPerRead[i])
		{
			if (std::get<0>(t) & 0x80000000)
			{
				assert(std::get<1>(t) & 0x80000000);
				std::swap(std::get<0>(t), std::get<1>(t));
				std::get<0>(t) ^= 0x80000000;
				std::get<1>(t) ^= 0x80000000;
				std::get<2>(t) ^= firstBitUint64_t;
				std::get<0>(t) = rawReadLengths[i] - 1 - std::get<0>(t);
				std::get<1>(t) = rawReadLengths[i] - 1 - std::get<1>(t);
				assert(std::get<0>(t) < std::get<1>(t));
				assert(std::get<1>(t) < rawReadLengths[i]);
			}
			else
			{
				assert(std::get<0>(t) < std::get<1>(t));
				assert(std::get<1>(t) < rawReadLengths[i]);
			}
		}
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
	}
	return chunksPerRead;
}

std::vector<std::pair<size_t, bool>> getParent(const MatchIndex& matchIndex, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<std::pair<size_t, bool>> parent;
	const size_t count = matchIndex.numWindowChunks() - matchIndex.numUniqueChunks();
	for (size_t i = 0; i < count; i++)
	{
		parent.emplace_back(i, true);
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (std::get<0>(chunksPerRead[i][j]) != std::get<0>(chunksPerRead[i][j-1])) continue;
			if (std::get<1>(chunksPerRead[i][j]) != std::get<1>(chunksPerRead[i][j-1])) continue;
			std::cerr << (std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t) << " " << (std::get<2>(chunksPerRead[i][j]) & maskUint64_t) << " " << (((std::get<2>(chunksPerRead[i][j-1]) & firstBitUint64_t) == (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t)) ? "fw" : "bw") << std::endl;
			merge(parent, std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t, std::get<2>(chunksPerRead[i][j]) & maskUint64_t, (std::get<2>(chunksPerRead[i][j-1]) & firstBitUint64_t) == (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t));
		}
	}
	return parent;
}

std::vector<bool> getFilteredValidChunks(const MatchIndex& matchIndex, const std::vector<size_t>& rawReadLengths, const size_t maxCoverage, const size_t repetitiveLength, const size_t repetitiveCoverage)
{
	std::vector<bool> useTheseChunks;
	useTheseChunks.resize(matchIndex.numWindowChunks() - matchIndex.numUniqueChunks(), true);
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead = getChunksPerRead(matchIndex, rawReadLengths, useTheseChunks);
	std::vector<std::pair<size_t, bool>> parent = getParent(matchIndex, chunksPerRead);
	std::vector<size_t> coverages;
	coverages.resize(parent.size(), 0);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (j > 1 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			auto key = find(parent, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			coverages[key.first] += 1;
		}
	}
	for (size_t i = 0; i < useTheseChunks.size(); i++)
	{
		if (coverages[find(parent, i).first] >= maxCoverage)
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
			if (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t) pairkey.second = !pairkey.second;
			uint64_t key = pairkey.first + (pairkey.second ? firstBitUint64_t : 0);
			if (lastPos.count(key) == 1)
			{
				if (lastPos.at(key) + repetitiveLength > std::get<0>(chunksPerRead[i][j]) && lastPos.at(key) != std::get<1>(chunksPerRead[i][j]))
				{
					countRepetitive[key & maskUint64_t] += 1;
				}
			}
			lastPos[key] = std::get<1>(chunksPerRead[i][j]);
		}
	}
	for (size_t i = 0; i < useTheseChunks.size(); i++)
	{
		if (countRepetitive[find(parent, i).first] >= repetitiveCoverage)
		{
			useTheseChunks[i] = false;
		}
	}
	return useTheseChunks;
}

