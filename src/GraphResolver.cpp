#include <iostream>
#include <phmap.h>
#include "MBGCommon.h"
#include "GraphResolver.h"

phmap::flat_hash_set<size_t> getPossiblyResolvableNodes(const UnitigGraph& unitigGraph, const SparseEdgeContainer& activeEdges, const double averageOneHaplotypeCoverage)
{
	phmap::flat_hash_set<size_t> result;
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (activeEdges.getEdges(std::make_pair(i, true)).size() < 2) continue;
		if (activeEdges.getEdges(std::make_pair(i, true)).size() != activeEdges.getEdges(std::make_pair(i, false)).size()) continue;
		if (unitigGraph.coverages[i] < averageOneHaplotypeCoverage * 1.5) continue;
		bool valid = true;
		double totalFwCoverage = 0;
		for (auto edge : activeEdges.getEdges(std::make_pair(i, true)))
		{
			auto key = canon(std::make_pair(i, true), edge);
			size_t edgeCoverage = unitigGraph.edgeCoverages.get(key.first, key.second);
			totalFwCoverage += edgeCoverage;
			if (edgeCoverage < averageOneHaplotypeCoverage * 0.5) valid = false;
		}
		if (!valid) continue;
		if (unitigGraph.coverages[i] + averageOneHaplotypeCoverage < totalFwCoverage) continue;
		if (unitigGraph.coverages[i] > totalFwCoverage + averageOneHaplotypeCoverage) continue;
		double totalBwCoverage = 0;
		for (auto edge : activeEdges.getEdges(std::make_pair(i, false)))
		{
			auto key = canon(std::make_pair(i, false), edge);
			size_t edgeCoverage = unitigGraph.edgeCoverages.get(key.first, key.second);
			totalBwCoverage += edgeCoverage;
			if (edgeCoverage < averageOneHaplotypeCoverage * 0.5) valid = false;
		}
		if (valid) result.insert(i);
		if (unitigGraph.coverages[i] + averageOneHaplotypeCoverage < totalBwCoverage) continue;
		if (unitigGraph.coverages[i] > totalBwCoverage + averageOneHaplotypeCoverage) continue;
	}
	return result;
}

phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, size_t>>> getValidResolutions(const UnitigGraph& unitigGraph, const SparseEdgeContainer& activeEdges, const std::vector<ReadPathBundle>& readPaths, const phmap::flat_hash_set<size_t>& possiblyResolvableNodes, const double averageOneHaplotypeCoverage)
{
	std::vector<phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>> tripletCoverages;
	tripletCoverages.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 1; k+1 < readPaths[i].paths[j].path.size(); k++)
			{
				if (possiblyResolvableNodes.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 0) continue;
				uint64_t prev = readPaths[i].paths[j].path[k-1];
				uint64_t next = readPaths[i].paths[j].path[k+1];
				if ((readPaths[i].paths[j].path[k] & firstBitUint64_t) == 0)
				{
					std::swap(prev, next);
					prev ^= firstBitUint64_t;
					next ^= firstBitUint64_t;
				}
				tripletCoverages[readPaths[i].paths[j].path[k] & maskUint64_t][std::make_pair(prev, next)] += 1;
			}
		}
	}
	std::vector<phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>> filteredTripletCoverages;
	filteredTripletCoverages.resize(tripletCoverages.size());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		for (auto pair : tripletCoverages[i])
		{
			if (pair.second < 2) continue;
			filteredTripletCoverages[i].emplace(pair);
		}
	}
	phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, size_t>>> result;
	for (const size_t i : possiblyResolvableNodes)
	{
		if (filteredTripletCoverages[i].size() != activeEdges.getEdges(std::make_pair(i, true)).size()) continue;
		size_t totalTripletCoverage = 0;
		phmap::flat_hash_set<size_t> allPrevs;
		phmap::flat_hash_set<size_t> allNexts;
		bool valid = true;
		for (auto pair : filteredTripletCoverages[i])
		{
			allPrevs.insert(pair.first.first);
			allNexts.insert(pair.first.second);
			totalTripletCoverage += pair.second;
			if ((double)pair.second * 1.5 < (double)unitigGraph.coverages[pair.first.first & maskUint64_t]) valid = false;
			if ((double)pair.second * 1.5 < (double)unitigGraph.coverages[pair.first.second & maskUint64_t]) valid = false;
			if ((double)pair.second > (double)unitigGraph.coverages[pair.first.first & maskUint64_t] * 1.5) valid = false;
			if ((double)pair.second > (double)unitigGraph.coverages[pair.first.second & maskUint64_t] * 1.5) valid = false;
		}
		if (!valid) continue;
		if (allPrevs.size() != allNexts.size()) continue;
		if (allNexts.size() != activeEdges.getEdges(std::make_pair(i, true)).size()) continue;
		if ((double)totalTripletCoverage * 1.5 < (double)unitigGraph.coverages[i]) continue;
		if ((double)totalTripletCoverage > (double)unitigGraph.coverages[i] * 1.5) continue;
		for (auto pair : filteredTripletCoverages[i])
		{
			result[i].push_back(pair.first);
		}
	}
	for (auto& pair : result)
	{
		std::sort(pair.second.begin(), pair.second.end());
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolve(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, size_t>>>& resolutions)
{
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<uint64_t, size_t>> prevNodeToResolution;
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<uint64_t, size_t>> nextNodeToResolution;
	size_t nextNewNode = unitigGraph.nodeCount();
	for (const auto& pair : resolutions)
	{
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			prevNodeToResolution[pair.first][pair.second[i].first] = nextNewNode+i;
			nextNodeToResolution[pair.first][pair.second[i].second] = nextNewNode+i;
			resultGraph.lengths.push_back(resultGraph.lengths[pair.first]);
			resultGraph.coverages.push_back(0);
		}
		nextNewNode += pair.second.size();
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				if (resolutions.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 0) continue;
				size_t resultNode = readPaths[i].paths[j].path[k] & maskUint64_t;
				uint64_t prev = std::numeric_limits<uint64_t>::max();
				if (k > 0) prev = readPaths[i].paths[j].path[k-1];
				uint64_t next = std::numeric_limits<uint64_t>::max();
				if (k < readPaths[i].paths[j].path.size()-1) next = readPaths[i].paths[j].path[k+1];
				if ((readPaths[i].paths[j].path[k] & firstBitUint64_t) == 0)
				{
					std::swap(prev, next);
					if (prev != std::numeric_limits<uint64_t>::max()) prev ^= firstBitUint64_t;
					if (next != std::numeric_limits<uint64_t>::max()) next ^= firstBitUint64_t;
				}
				if (prev != std::numeric_limits<uint64_t>::max() && next != std::numeric_limits<uint64_t>::max())
				{
					assert(prevNodeToResolution.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 1);
					assert(nextNodeToResolution.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 1);
					assert(prevNodeToResolution.at(readPaths[i].paths[j].path[k] & maskUint64_t).count(prev) == 1);
					assert(nextNodeToResolution.at(readPaths[i].paths[j].path[k] & maskUint64_t).count(next) == 1);
					size_t prevResult = prevNodeToResolution.at(readPaths[i].paths[j].path[k] & maskUint64_t).at(prev);
					size_t nextResult = nextNodeToResolution.at(readPaths[i].paths[j].path[k] & maskUint64_t).at(next);
					if (prevResult == nextResult) resultNode = prevResult;
				}
				else if (prev != std::numeric_limits<uint64_t>::max())
				{
					resultNode = prevNodeToResolution.at(readPaths[i].paths[j].path[k] & maskUint64_t).at(prev);
				}
				else if (next != std::numeric_limits<uint64_t>::max())
				{
					resultNode = nextNodeToResolution.at(readPaths[i].paths[j].path[k] & maskUint64_t).at(next);
				}
				else
				{
					assert(readPaths[i].paths[j].path.size() == 1);
				}
				resultPaths[i].paths[j].path[k] = resultNode + (readPaths[i].paths[j].path[k] & firstBitUint64_t);
			}
		}
	}
	return std::make_pair(std::move(resultGraph), std::move(resultPaths));
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolveSimpleStructures(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double averageOneHaplotypeCoverage)
{
	auto edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	phmap::flat_hash_set<size_t> possiblyResolvableNodes = getPossiblyResolvableNodes(unitigGraph, edges, averageOneHaplotypeCoverage);
	auto resolutions = getValidResolutions(unitigGraph, edges, readPaths, possiblyResolvableNodes, averageOneHaplotypeCoverage);
	std::cerr << resolutions.size() << " simple repeat nodes resolved" << std::endl;
	UnitigGraph newGraph;
	std::vector<ReadPathBundle> newPaths;
	std::tie(newGraph, newPaths) = resolve(unitigGraph, readPaths, resolutions);
	RankBitvector kept;
	kept.resize(newGraph.nodeCount());
	for (size_t i = 0; i < kept.size(); i++)
	{
		kept.set(i, true);
	}
	for (const auto& pair : resolutions)
	{
		kept.set(pair.first, false);
	}
	kept.buildRanks();
	return filterUnitigGraph(newGraph, newPaths, kept);
}
