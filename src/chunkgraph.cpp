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
			if (std::get<0>(t) & 0x70000000)
			{
				assert(std::get<1>(t) & 0x70000000);
				std::swap(std::get<0>(t), std::get<1>(t));
				std::get<0>(t) ^= 0x70000000;
				std::get<1>(t) ^= 0x70000000;
				std::get<2>(t) ^= firstBitUint64_t;
				std::get<0>(t) = rawReadLengths[i] - std::get<0>(t);
				std::get<1>(t) = rawReadLengths[i] - std::get<1>(t);
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
			merge(parent, std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t, std::get<2>(chunksPerRead[i][j]) & maskUint64_t, (std::get<2>(chunksPerRead[i][j-1]) & firstBitUint64_t) == (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t));
		}
	}
	return parent;
}

void writeGraph(const std::string& filename, const std::vector<std::vector<size_t>>& lengths, const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeCoverage)
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
		graph << "L\t" << (pair.first.first & maskUint64_t) << "\t" << ((pair.first.first & firstBitUint64_t) ? "+" : "-") << "\t" << (pair.first.second & maskUint64_t) << "\t" << ((pair.first.second & firstBitUint64_t) ? "+" : "-") << "\t0M\tec:i:" << pair.second << std::endl;
	}
}

std::pair<std::vector<std::vector<size_t>>, std::vector<size_t>> getLengthsAndCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, std::vector<std::pair<size_t, bool>>& parent)
{
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	coverages.resize(parent.size(), 0);
	lengths.resize(parent.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto key = find(parent, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			lengths[key.first].emplace_back(std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]));
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (j > 1 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			auto key = find(parent, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			coverages[key.first] += 1;
		}
	}
	return std::make_pair(lengths, coverages);
}

phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getEdgeCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, std::vector<std::pair<size_t, bool>>& parent)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			auto prev = find(parent, std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t);
			auto curr = find(parent, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			if (std::get<2>(chunksPerRead[i][j-1]) & firstBitUint64_t) prev.second = !prev.second;
			if (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t) curr.second = !curr.second;
			auto pairkey = canon(prev, curr);
			std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
			edgeCoverage[key] += 1;
		}
	}
	return edgeCoverage;
}

void makeGraph(const MatchIndex& matchIndex, const std::vector<size_t>& rawReadLengths)
{
	std::vector<bool> useTheseChunks;
	useTheseChunks.resize(matchIndex.numWindowChunks() - matchIndex.numUniqueChunks(), true);
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead = getChunksPerRead(matchIndex, rawReadLengths, useTheseChunks);
	std::vector<std::pair<size_t, bool>> parent = getParent(matchIndex, chunksPerRead);
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead, parent);
	{
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead, parent);
		writeGraph("graph-round1.gfa", lengths, coverages, edgeCoverage);
	}
	for (size_t i = 0; i < useTheseChunks.size(); i++)
	{
		if (coverages[find(parent, i).first] >= 40)
		{
			useTheseChunks[i] = false;
		}
	}
	// chunksPerRead = getChunksPerRead(matchIndex, rawReadLengths, useTheseChunks);
	// parent = getParent(matchIndex, chunksPerRead);
	// {
	// 	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead, parent);
	// 	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead, parent);
	// 	writeGraph("graph-round2.gfa", lengths, coverages, edgeCoverage);
	// }
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
	parent = getParent(matchIndex, chunksPerRead);
	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead, parent);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead, parent);
	writeGraph("graph.gfa", lengths, coverages, edgeCoverage);
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
	auto rawReadLengths = storage.getRawReadLengths();
	makeGraph(matchIndex, rawReadLengths);
}
