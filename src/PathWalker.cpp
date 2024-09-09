#include <map>
#include <limits>
#include <iostream>
#include "ChunkGraphWriter.h"
#include "PathWalker.h"
#include "Common.h"
#include "phmap.h"
#include "UnionFind.h"

double getAverageLongNodeCoverage(const ChunkUnitigGraph& graph, const size_t minLength)
{
	double totalLength = 0;
	double totalCoverage = 0;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < minLength) continue;
		totalLength += graph.unitigLengths[i];
		totalCoverage += graph.unitigLengths[i] * graph.coverages[i];
	}
	return totalCoverage / totalLength;
}

phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>> getNodeSpanningTriplets(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double estimatedSingleCopyCoverage)
{
	std::vector<bool> canBeSpanned;
	canBeSpanned.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.edges.getEdges(std::make_pair(i, true)).size() < 2) continue;
		if (graph.edges.getEdges(std::make_pair(i, false)).size() != graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		bool valid = true;
		for (auto edge : graph.edges.getEdges(std::make_pair(i, true)))
		{
			if (graph.edges.getEdges(reverse(edge)).size() != 1)
			{
				valid = false;
				break;
			}
		}
		for (auto edge : graph.edges.getEdges(std::make_pair(i, false)))
		{
			if (graph.edges.getEdges(reverse(edge)).size() != 1)
			{
				valid = false;
				break;
			}
		}
		if (!valid) continue;
		canBeSpanned[i] = true;
	}
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> spanningTripletCoverage;
	spanningTripletCoverage.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				if (!canBeSpanned[readPaths[i][j].path[k] & maskUint64_t]) continue;
				uint64_t prev = readPaths[i][j].path[k-1];
				uint64_t after = readPaths[i][j].path[k+1];
				if (readPaths[i][j].path[k] & firstBitUint64_t)
				{
				}
				else
				{
					std::swap(prev, after);
					prev ^= firstBitUint64_t;
					after ^= firstBitUint64_t;
				}
				spanningTripletCoverage[readPaths[i][j].path[k] & maskUint64_t][std::make_pair(prev, after)] += 1;
			}
		}
	}
	phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>> result;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (!canBeSpanned[i]) continue;
		size_t unitigCopyCount = graph.edges.getEdges(std::make_pair(i, true)).size();
		assert(unitigCopyCount >= 2);
		phmap::flat_hash_set<uint64_t> coveredPredecessors;
		phmap::flat_hash_set<uint64_t> coveredSuccessors;
		phmap::flat_hash_set<std::pair<size_t, size_t>> solids;
		for (auto pair : spanningTripletCoverage[i])
		{
			if (pair.second < 2) continue;
			coveredPredecessors.insert(pair.first.first);
			coveredSuccessors.insert(pair.first.second);
			solids.insert(pair.first);
		}
		if (solids.size() > unitigCopyCount) continue;
		if (spanningTripletCoverage[i].size() > unitigCopyCount+1) continue;
		if (solids.size() == unitigCopyCount-1 && spanningTripletCoverage[i].size() == unitigCopyCount)
		{
			for (auto pair : spanningTripletCoverage[i])
			{
				coveredPredecessors.insert(pair.first.first);
				coveredSuccessors.insert(pair.first.second);
				solids.insert(pair.first);
			}
		}
		if (solids.size() < 2 && unitigCopyCount == 2)
		{
			if (solids.size() == 1 && spanningTripletCoverage[i].size() == 1)
			{
				phmap::flat_hash_set<uint64_t> uncoveredPredecessor;
				phmap::flat_hash_set<uint64_t> uncoveredSuccessor;
				for (auto edge : graph.edges.getEdges(std::make_pair(i, true)))
				{
					uncoveredSuccessor.insert(edge.first + (edge.second ? firstBitUint64_t : 0));
				}
				for (auto edge : graph.edges.getEdges(std::make_pair(i, false)))
				{
					uncoveredPredecessor.insert(edge.first + (edge.second ? 0 : firstBitUint64_t));
				}
				for (auto pair : solids)
				{
					assert(uncoveredPredecessor.count(pair.first) == 1);
					assert(uncoveredSuccessor.count(pair.second) == 1);
					uncoveredPredecessor.erase(pair.first);
					uncoveredSuccessor.erase(pair.second);
				}
				if (uncoveredSuccessor.size() == 1 && uncoveredPredecessor.size() == 1)
				{
					std::pair<uint64_t, uint64_t> newUncoveredTriplet { *uncoveredPredecessor.begin(), *uncoveredSuccessor.begin() };
					solids.insert(newUncoveredTriplet);
					coveredPredecessors.insert(newUncoveredTriplet.first);
					coveredSuccessors.insert(newUncoveredTriplet.second);
				}
			}
			if (solids.size() == 0 && spanningTripletCoverage[i].size() == 2)
			{
				for (auto pair : spanningTripletCoverage[i])
				{
					coveredPredecessors.insert(pair.first.first);
					coveredSuccessors.insert(pair.first.second);
					solids.insert(pair.first);
				}
			}
		}
		if (coveredPredecessors.size() != unitigCopyCount) continue;
		if (coveredSuccessors.size() != unitigCopyCount) continue;
		result[i].insert(result[i].end(), solids.begin(), solids.end());
	}
	return result;
}

std::vector<bool> estimateUniqueUnitigs(const ChunkUnitigGraph& graph, const size_t longUnitigThreshold, const double estimatedSingleCopyCoverage)
{
	std::vector<bool> result;
	result.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < longUnitigThreshold) continue;
		if (graph.coverages[i] > 1.5 * estimatedSingleCopyCoverage) continue;
		result[i] = true;
	}
	return result;
}

uint64_t getNextFromTriplet(const phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>>& triplets, const uint64_t previous, const uint64_t current)
{
	if (triplets.count(current & maskUint64_t) == 0) return std::numeric_limits<size_t>::max();
	if (current & firstBitUint64_t)
	{
		for (auto pair : triplets.at(current & maskUint64_t))
		{
			if (pair.first == previous)
			{
				return pair.second;
			}
		}
	}
	else
	{
		for (auto pair : triplets.at(current & maskUint64_t))
		{
			if (pair.second == (previous ^ firstBitUint64_t))
			{
				return pair.first ^ firstBitUint64_t;
			}
		}
	}
	return std::numeric_limits<size_t>::max();
}

std::vector<uint64_t> followTripletPath(const ChunkUnitigGraph& graph, const uint64_t startUnitig, const phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>>& triplets)
{
	std::vector<uint64_t> result;
	result.emplace_back(startUnitig);
	while (true)
	{
		if (graph.edges.getEdges(std::make_pair(result.back() & maskUint64_t, result.back() & firstBitUint64_t)).size() == 1)
		{
			auto nextNode = graph.edges.getEdges(std::make_pair(result.back() & maskUint64_t, result.back() & firstBitUint64_t))[0];
			uint64_t nextTripletNode = getNextFromTriplet(triplets, result.back(), nextNode.first + (nextNode.second ? firstBitUint64_t : 0));
			if (nextTripletNode != std::numeric_limits<size_t>::max())
			{
				result.emplace_back(nextNode.first + (nextNode.second ? firstBitUint64_t : 0));
				result.emplace_back(nextTripletNode);
				continue;
			}
		}
		if (triplets.count(result.back() & maskUint64_t) == 1)
		{
			if (result.size() >= 2)
			{
				uint64_t nextTripletNode = getNextFromTriplet(triplets, result[result.size()-2], result.back());
				if (nextTripletNode != std::numeric_limits<size_t>::max())
				{
					result.emplace_back(nextTripletNode);
					continue;
				}
			}
		}
		break;
	}
	return result;
}

size_t getPathLength(const ChunkUnitigGraph& graph, const std::vector<uint64_t>& path)
{
	size_t pathLength = 0;
	for (size_t i = 0; i < path.size(); i++)
	{
		pathLength += graph.unitigLengths[path[i] & maskUint64_t];
		if (i > 0)
		{
			auto key = canonNodePair(path[i-1], path[i]);
			size_t overlap = graph.edgeOverlaps.at(key);
			assert(pathLength >= overlap);
			pathLength -= overlap;
		}
	}
	return pathLength;
}

std::vector<std::vector<uint64_t>> getTripletExtendedPaths(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double estimatedSingleCopyCoverage)
{
	auto triplets = getNodeSpanningTriplets(graph, readPaths, estimatedSingleCopyCoverage);
	std::vector<std::vector<uint64_t>> result;
	std::vector<bool> alreadyInPath;
	alreadyInPath.resize(graph.unitigLengths.size(), false);
	for (const auto& pair : triplets)
	{
		for (const auto pair2 : pair.second)
		{
			if (alreadyInPath[pair2.first & maskUint64_t])
			{
				assert(alreadyInPath[pair2.second & maskUint64_t]);
				continue;
			}
			assert(!alreadyInPath[pair2.second & maskUint64_t]);
			auto bwPath = followTripletPath(graph, pair2.first ^ firstBitUint64_t, triplets);
			auto fwPath = followTripletPath(graph, pair2.second, triplets);
			assert(bwPath.size() >= 1);
			assert(fwPath.size() >= 1);
			std::reverse(bwPath.begin(), bwPath.end());
			for (size_t i = 0; i < bwPath.size(); i++)
			{
				bwPath[i] ^= firstBitUint64_t;
			}
			bwPath.emplace_back(pair.first + firstBitUint64_t);
			bwPath.insert(bwPath.end(), fwPath.begin(), fwPath.end());
			for (auto node : bwPath)
			{
				alreadyInPath[node & maskUint64_t] = true;
			}
			size_t pathLength = getPathLength(graph, bwPath);
			if (pathLength < 100000) continue;
			result.emplace_back(bwPath);
		}
	}
	return result;
}

std::vector<std::vector<uint64_t>> addUniqueNodePaths(const std::vector<std::vector<uint64_t>>& paths, const std::vector<bool>& isUniqueUnitig)
{
	std::vector<bool> alreadyInPath;
	alreadyInPath.resize(isUniqueUnitig.size(), false);
	std::vector<std::vector<uint64_t>> result = paths;
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 0; j < paths[i].size(); j++)
		{
			alreadyInPath[paths[i][j] & maskUint64_t] = true;
		}
	}
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (alreadyInPath[i]) continue;
		if (!isUniqueUnitig[i]) continue;
		result.emplace_back();
		result.back().emplace_back(i + firstBitUint64_t);
	}
	return result;
}

std::pair<std::vector<std::vector<uint64_t>>, std::vector<bool>> addRemainingSolidNodesAsPaths(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const std::vector<std::vector<uint64_t>>& contigPaths, const std::vector<bool>& removedNodes, const double estimatedSingleCopyCoverage)
{
	std::vector<bool> alreadyInPath;
	alreadyInPath.resize(graph.unitigLengths.size(), false);
	std::vector<std::vector<uint64_t>> result = contigPaths;
	for (size_t i = 0; i < contigPaths.size(); i++)
	{
		for (size_t j = 0; j < contigPaths[i].size(); j++)
		{
			alreadyInPath[contigPaths[i][j] & maskUint64_t] = true;
		}
	}
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (alreadyInPath[i]) continue;
		if (removedNodes[i]) continue;
		if (graph.unitigLengths[i] < 100000 && graph.coverages[i] < estimatedSingleCopyCoverage) continue;
		result.emplace_back();
		result.back().emplace_back(i + firstBitUint64_t);
	}
	return std::make_pair(result, removedNodes);
}

std::pair<std::vector<std::vector<uint64_t>>, std::vector<bool>> mergeContigPathsByReadPaths(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const std::vector<std::vector<uint64_t>>& contigPaths)
{
	phmap::flat_hash_map<size_t, uint64_t> unitigBelongsUniquelyToContigPath;
	phmap::flat_hash_set<size_t> unitigIsSharedInPaths;
	phmap::flat_hash_map<uint64_t, uint64_t> unitigIsContigPathTip;
	for (size_t i = 0; i < contigPaths.size(); i++)
	{
		for (size_t j = 0; j < contigPaths[i].size(); j++)
		{
			size_t node = contigPaths[i][j] & maskUint64_t;
			uint64_t path = i + (contigPaths[i][j] & firstBitUint64_t);
			if (unitigBelongsUniquelyToContigPath.count(node) == 1 && unitigBelongsUniquelyToContigPath.at(node) != path)
			{
				unitigIsSharedInPaths.insert(node);
				unitigBelongsUniquelyToContigPath[node] = std::numeric_limits<size_t>::max();
				continue;
			}
			if (unitigBelongsUniquelyToContigPath.count(node) == 0)
			{
				unitigBelongsUniquelyToContigPath[node] = path;
			}
		}
	}
	for (size_t i = 0; i < contigPaths.size(); i++)
	{
		assert(unitigIsSharedInPaths.count(contigPaths[i][0] & maskUint64_t) == 0);
		assert(unitigBelongsUniquelyToContigPath.at(contigPaths[i][0] & maskUint64_t) == (i + (contigPaths[i][0] & firstBitUint64_t)));
		assert(unitigIsSharedInPaths.count(contigPaths[i].back() & maskUint64_t) == 0);
		assert(unitigBelongsUniquelyToContigPath.at(contigPaths[i].back() & maskUint64_t) == (i + (contigPaths[i].back() & firstBitUint64_t)));
		assert(unitigIsContigPathTip.count(contigPaths[i][0] ^ firstBitUint64_t) == 0);
		assert(unitigIsContigPathTip.count(contigPaths[i].back()) == 0);
		assert((contigPaths[i][0] ^ firstBitUint64_t) != (contigPaths[i].back()));
		unitigIsContigPathTip[contigPaths[i][0] ^ firstBitUint64_t] = i;
		unitigIsContigPathTip[contigPaths[i].back()] = i + firstBitUint64_t;
	}
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<std::tuple<size_t, size_t, size_t, size_t>>> bridges; // (contig, contig) -> (read, readpathindex, pathstart, pathend INCLUSIVE)
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			uint64_t lastContigTip = std::numeric_limits<size_t>::max();
			uint64_t lastContigTipIndex = 0;
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (unitigIsContigPathTip.count(readPaths[i][j].path[k] ^ firstBitUint64_t) == 1)
				{
					if (lastContigTip != std::numeric_limits<size_t>::max())
					{
						uint64_t thisContigTip = unitigIsContigPathTip.at(readPaths[i][j].path[k] ^ firstBitUint64_t);
						bridges[std::make_pair(lastContigTip, thisContigTip ^ firstBitUint64_t)].emplace_back(i, j, lastContigTipIndex, k);
					}
					lastContigTip = std::numeric_limits<size_t>::max();
				}
				if (unitigIsContigPathTip.count(readPaths[i][j].path[k]) == 1)
				{
					lastContigTip = unitigIsContigPathTip.at(readPaths[i][j].path[k]);
					lastContigTipIndex = k;
				}
			}
		}
	}
	phmap::flat_hash_map<uint64_t, uint64_t> nodeTipBelongsToUnresolvedTangle;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		nodeTipBelongsToUnresolvedTangle[i + firstBitUint64_t] = i + firstBitUint64_t;
		nodeTipBelongsToUnresolvedTangle[i] = i;
	}
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		for (auto edge : graph.edges.getEdges(std::make_pair(i, true)))
		{
			merge(nodeTipBelongsToUnresolvedTangle, i + firstBitUint64_t, edge.first + (edge.second ? 0 : firstBitUint64_t)); // reverse edge!
		}
		for (auto edge : graph.edges.getEdges(std::make_pair(i, false)))
		{
			merge(nodeTipBelongsToUnresolvedTangle, i, edge.first + (edge.second ? 0 : firstBitUint64_t)); // reverse edge!
		}
	}
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (unitigIsSharedInPaths.count(i) == 1) continue;
		if (unitigBelongsUniquelyToContigPath.count(i) == 1) continue;
		merge(nodeTipBelongsToUnresolvedTangle, i, i + firstBitUint64_t);
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint64_t>> contigTipsTouchingTangle;
	phmap::flat_hash_map<uint64_t, uint64_t> contigTipInTangle;
	for (auto pair : unitigIsContigPathTip)
	{
		uint64_t tangle = find(nodeTipBelongsToUnresolvedTangle, pair.first);
		contigTipsTouchingTangle[tangle].insert(pair.second);
		contigTipInTangle[pair.second] = tangle;
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> bridgeCoveragesPerTangle;
	for (const auto& pair : bridges)
	{
		assert(contigTipInTangle.count(pair.first.first) == 1);
		assert(contigTipInTangle.count(pair.first.second ^ firstBitUint64_t) == 1);
		assert(contigTipInTangle.at(pair.first.first) == contigTipInTangle.at(pair.first.second ^ firstBitUint64_t));
		uint64_t tangle = contigTipInTangle.at(pair.first.first);
		auto key = canonNodePair(pair.first.first, pair.first.second);
		bridgeCoveragesPerTangle[tangle][key] += pair.second.size();
	}
	std::vector<bool> removedNodes;
	removedNodes.resize(graph.unitigLengths.size(), false);
	phmap::flat_hash_set<uint64_t> solvedTangles;
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> chosenMerges;
	for (const auto& tangle : bridgeCoveragesPerTangle)
	{
		if (contigTipsTouchingTangle.at(tangle.first).size() % 2 != 0)
		{
			continue;
		}
		phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, size_t>> bridgeCoverages; // from -> (to -> coverage)
		for (auto pair : tangle.second)
		{
			bridgeCoverages[pair.first.first][pair.first.second] = pair.second;
			bridgeCoverages[pair.first.second ^ firstBitUint64_t][pair.first.first ^ firstBitUint64_t] = pair.second;
		}
		phmap::flat_hash_map<uint64_t, uint64_t> uniqueBestMatch;
		for (const auto& pair : bridgeCoverages)
		{
			size_t highestCoverageMatch = std::numeric_limits<size_t>::max();
			size_t highestMatchCoverage = 0;
			for (auto pair2 : pair.second)
			{
				if (pair2.second == highestMatchCoverage)
				{
					highestCoverageMatch = std::numeric_limits<size_t>::max();
				}
				if (pair2.second > highestMatchCoverage)
				{
					highestMatchCoverage = pair2.second;
					highestCoverageMatch = pair2.first;
				}
			}
			if (highestCoverageMatch != std::numeric_limits<size_t>::max()) uniqueBestMatch[pair.first] = highestCoverageMatch;
		}
		if (uniqueBestMatch.size() != contigTipsTouchingTangle.at(tangle.first).size()) continue;
		bool valid = true;
		for (auto pair : uniqueBestMatch)
		{
			assert(uniqueBestMatch.count(pair.second ^ firstBitUint64_t) == 1);
			if (uniqueBestMatch.at(pair.second ^ firstBitUint64_t) != (pair.first ^ firstBitUint64_t))
			{
				valid = false;
				break;
			}
		}
		if (!valid) continue;
		solvedTangles.insert(tangle.first);
		for (auto pair : uniqueBestMatch)
		{
			if (std::make_pair(pair.first, pair.second) != canonNodePair(pair.first, pair.second)) continue;
			chosenMerges.insert(pair);
		}
	}
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<std::vector<uint64_t>>> mergeNodeSequences;
	for (const auto& pair : bridges)
	{
		bool fw = true;
		auto canoncheck = canonNodePair(pair.first.first, pair.first.second);
		if (canoncheck != pair.first)
		{
			fw = false;
		}
		if (chosenMerges.count(canoncheck) == 0) continue;
		for (auto t : pair.second)
		{
			std::vector<uint64_t> nodeSequence;
			assert(std::get<0>(t) < readPaths.size());
			assert(std::get<1>(t) < readPaths[std::get<0>(t)].size());
			assert(std::get<2>(t) <= std::get<3>(t));
			assert(std::get<3>(t) < readPaths[std::get<0>(t)][std::get<1>(t)].path.size());
			for (size_t i = std::get<2>(t); i <= std::get<3>(t); i++)
			{
				nodeSequence.emplace_back(readPaths[std::get<0>(t)][std::get<1>(t)].path[i]);
			}
			if (!fw)
			{
				std::reverse(nodeSequence.begin(), nodeSequence.end());
				for (size_t i = 0; i < nodeSequence.size(); i++)
				{
					nodeSequence[i] ^= firstBitUint64_t;
				}
			}
			mergeNodeSequences[canoncheck].emplace_back(nodeSequence);
		}
	}
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<uint64_t>> chosenMergeNodeSequence;
	for (const auto& pair : mergeNodeSequences)
	{
		std::map<std::vector<uint64_t>, size_t> pathCoverage;
		for (const auto& path : pair.second)
		{
			pathCoverage[path] += 1;
		}
		size_t highestPathCoverage = 0;
		for (const auto& pair2 : pathCoverage)
		{
			highestPathCoverage = std::max(highestPathCoverage, pair2.second);
		}
		for (const auto& pair2 : pathCoverage)
		{
			if (pair2.second == highestPathCoverage)
			{
				chosenMergeNodeSequence[pair.first] = pair2.first;
				break;
			}
		}
		assert(chosenMergeNodeSequence.count(pair.first) == 1);
	}
	phmap::flat_hash_map<uint64_t, uint64_t> pathOutEdge;
	for (const auto& pair : chosenMergeNodeSequence)
	{
		assert(pathOutEdge.count(pair.first.first) == 0);
		assert(pathOutEdge.count(pair.first.second ^ firstBitUint64_t) == 0);
		assert((pair.first.first & maskUint64_t) != (pair.first.second & maskUint64_t));
		pathOutEdge[pair.first.first] = pair.first.second;
		pathOutEdge[pair.first.second ^ firstBitUint64_t] = pair.first.first ^ firstBitUint64_t;
	}
	std::vector<std::vector<uint64_t>> newContigPathPaths;
	std::vector<bool> contigPathInPath;
	contigPathInPath.resize(contigPaths.size(), false);
	for (size_t i = 0; i < contigPaths.size(); i++)
	{
		if (contigPathInPath[i]) continue;
		std::vector<uint64_t> bwPath;
		bwPath.emplace_back(i);
		while (pathOutEdge.count(bwPath.back()) == 1)
		{
			if (pathOutEdge.at(bwPath.back()) == (bwPath.back() ^ firstBitUint64_t)) break;
			bwPath.emplace_back(pathOutEdge.at(bwPath.back()));
			if (bwPath.back() == i) break;
			assert(bwPath.size() <= contigPaths.size());
		}
		std::vector<uint64_t> fwPath;
		fwPath.emplace_back(i + firstBitUint64_t);
		while (pathOutEdge.count(fwPath.back()) == 1)
		{
			if (pathOutEdge.at(fwPath.back()) == (fwPath.back() ^ firstBitUint64_t)) break;
			fwPath.emplace_back(pathOutEdge.at(fwPath.back()));
			if (fwPath.back() == i + firstBitUint64_t) break;
			assert(fwPath.size() <= contigPaths.size());
		}
		if (bwPath.size() >= 2 && bwPath.back() == bwPath[0])
		{
			// circular path
			assert(fwPath.size() >= 2);
			assert(fwPath[0] == fwPath.back());
			assert(fwPath[0] == (bwPath[0] ^ firstBitUint64_t));
			fwPath.erase(fwPath.begin()+1, fwPath.end());
			bwPath.pop_back();
		}
		std::reverse(bwPath.begin(), bwPath.end());
		for (size_t j = 0; j < bwPath.size(); j++)
		{
			bwPath[j] ^= firstBitUint64_t;
		}
		assert(bwPath.size() >= 1);
		assert(fwPath.size() >= 1);
		assert(bwPath.back() == fwPath[0]);
		bwPath.insert(bwPath.end(), fwPath.begin()+1, fwPath.end());
		for (size_t j = 0; j < bwPath.size(); j++)
		{
			assert(!contigPathInPath[bwPath[j] & maskUint64_t]);
			contigPathInPath[bwPath[j] & maskUint64_t] = true;
		}
		assert(bwPath.size() >= 1);
		newContigPathPaths.emplace_back(bwPath);
		assert(newContigPathPaths.back().size() >= 1);
		assert(contigPathInPath[i]);
	}
	std::vector<std::vector<uint64_t>> result;
	for (size_t i = 0; i < newContigPathPaths.size(); i++)
	{
		result.emplace_back();
		for (size_t j = 0; j < newContigPathPaths[i].size(); j++)
		{
			std::vector<uint64_t> addedHere;
			addedHere = contigPaths[newContigPathPaths[i][j] & maskUint64_t];
			if (newContigPathPaths[i][j] & firstBitUint64_t)
			{
			}
			else
			{
				std::reverse(addedHere.begin(), addedHere.end());
				for (size_t k = 0; k < addedHere.size(); k++)
				{
					addedHere[k] ^= firstBitUint64_t;
				}
			}
			if (j >= 1)
			{
				auto key = canonNodePair(newContigPathPaths[i][j-1], newContigPathPaths[i][j]);
				assert(chosenMergeNodeSequence.count(key) == 1);
				std::vector<uint64_t> bridgeSequence = chosenMergeNodeSequence.at(key);
				assert(bridgeSequence.size() >= 2);
				if (key != std::make_pair(newContigPathPaths[i][j-1], newContigPathPaths[i][j]))
				{
					std::reverse(bridgeSequence.begin(), bridgeSequence.end());
					for (size_t k = 0; k < bridgeSequence.size(); k++)
					{
						bridgeSequence[k] ^= firstBitUint64_t;
					}
				}
				assert(bridgeSequence[0] == result[i].back());
				assert(addedHere[0] == bridgeSequence.back());
				result[i].insert(result[i].end(), bridgeSequence.begin()+1, bridgeSequence.end());
				result[i].insert(result[i].end(), addedHere.begin()+1, addedHere.end());
			}
			else
			{
				result[i].insert(result[i].end(), addedHere.begin(), addedHere.end());
			}
		}
		assert(result[i].size() >= 1);
	}
	removedNodes.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (unitigIsSharedInPaths.count(i) == 1) continue;
		if (unitigBelongsUniquelyToContigPath.count(i) == 1) continue;
		assert(find(nodeTipBelongsToUnresolvedTangle, i) == find(nodeTipBelongsToUnresolvedTangle, i + firstBitUint64_t));
		uint64_t tangle = find(nodeTipBelongsToUnresolvedTangle, i);
		if (solvedTangles.count(tangle) == 0) continue;
		removedNodes[i] = true;
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			removedNodes[result[i][j] & maskUint64_t] = false;
		}
	}
	return std::make_pair(result, removedNodes);
}

void getContigPathsAndConsensuses(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t kmerSize, const std::string& outPathsFile, const std::string& outContigsFasta)
{
	const size_t longUnitigThreshold = 100000;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	auto fixedChunksPerRead = getBidirectedChunks(chunksPerRead);
	std::tie(graph, readPaths) = getChunkUnitigGraph(fixedChunksPerRead, approxOneHapCoverage, kmerSize);
	double estimatedSingleCopyCoverage = getAverageLongNodeCoverage(graph, longUnitigThreshold);
	std::cerr << "estimated single copy coverage " << estimatedSingleCopyCoverage << std::endl;
	std::vector<bool> isUniqueUnitig = estimateUniqueUnitigs(graph, longUnitigThreshold, estimatedSingleCopyCoverage);
	std::cerr << "unique unitigs: ";
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (!isUniqueUnitig[i]) continue;
		std::cerr << i << ",";
	}
	std::cerr << std::endl;
	auto paths = getTripletExtendedPaths(graph, readPaths, estimatedSingleCopyCoverage);
	std::cerr << "step 1 paths:" << std::endl;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::cerr << "step 1 path id " << i << " length " << getPathLength(graph, paths[i]) << ": ";
		for (auto node : paths[i])
		{
			std::cerr << ((node & firstBitUint64_t) ? ">" : "<") << (node & maskUint64_t);
		}
		std::cerr << std::endl;
		for (auto node : paths[i])
		{
			std::cerr << (node & maskUint64_t) << ",";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	paths = addUniqueNodePaths(paths, isUniqueUnitig);
	std::cerr << "step 2 paths:" << std::endl;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::cerr << "step 2 path id " << i << " length " << getPathLength(graph, paths[i]) << ": ";
		for (auto node : paths[i])
		{
			std::cerr << ((node & firstBitUint64_t) ? ">" : "<") << (node & maskUint64_t);
		}
		std::cerr << std::endl;
		for (auto node : paths[i])
		{
			std::cerr << (node & maskUint64_t) << ",";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	std::vector<bool> removedNodes;
	std::tie(paths, removedNodes) = mergeContigPathsByReadPaths(graph, readPaths, paths);
	std::cerr << "step 3 paths:" << std::endl;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::cerr << "step 3 path id " << i << " length " << getPathLength(graph, paths[i]) << ": ";
		for (auto node : paths[i])
		{
			std::cerr << ((node & firstBitUint64_t) ? ">" : "<") << (node & maskUint64_t);
		}
		std::cerr << std::endl;
		for (auto node : paths[i])
		{
			std::cerr << (node & maskUint64_t) << ",";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	std::tie(paths, removedNodes) = addRemainingSolidNodesAsPaths(graph, readPaths, paths, removedNodes, estimatedSingleCopyCoverage);
	std::cerr << "step 4 paths:" << std::endl;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::cerr << "step 4 path id " << i << " length " << getPathLength(graph, paths[i]) << ": ";
		for (auto node : paths[i])
		{
			std::cerr << ((node & firstBitUint64_t) ? ">" : "<") << (node & maskUint64_t);
		}
		std::cerr << std::endl;
		for (auto node : paths[i])
		{
			std::cerr << (node & maskUint64_t) << ",";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	std::cerr << "removed nodes: ";
	for (size_t i = 0; i < removedNodes.size(); i++)
	{
		if (!removedNodes[i]) continue;
		std::cerr << i << ",";
	}
	std::cerr << std::endl;
}
