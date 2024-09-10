#include <mutex>
#include <fstream>
#include <map>
#include <limits>
#include <iostream>
#include "ChunkGraphWriter.h"
#include "PathWalker.h"
#include "Common.h"
#include "phmap.h"
#include "UnionFind.h"

std::vector<uint64_t> mergeBwFwPaths(const std::vector<uint64_t>& bwPath, const std::vector<uint64_t>& fwPath)
{
	assert(bwPath.size() >= 1);
	assert(fwPath.size() >= 1);
	assert(fwPath[0] == (bwPath[0] ^ firstBitUint64_t));
	if (bwPath.size() >= 2 && bwPath[0] == bwPath.back())
	{
		// circular path
		assert(fwPath[0] == (bwPath[0] ^ firstBitUint64_t));
		assert(fwPath[0] == fwPath.back());
		std::vector<uint64_t> result;
		result.insert(result.end(), fwPath.begin(), fwPath.end());
		result.pop_back();
		return result;
	}
	std::vector<uint64_t> result;
	result.insert(result.end(), bwPath.begin(), bwPath.end());
	std::reverse(result.begin(), result.end());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] ^= firstBitUint64_t;
	}
	assert(result.back() == fwPath[0]);
	result.insert(result.end(), fwPath.begin()+1, fwPath.end());
	return result;
}

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

std::vector<bool> estimateUniqueUnitigs(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>> readPaths, const size_t longUnitigThreshold, const double estimatedSingleCopyCoverage)
{
	auto triplets = getNodeSpanningTriplets(graph, readPaths, estimatedSingleCopyCoverage);
	std::vector<bool> result;
	result.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < longUnitigThreshold) continue;
		if (graph.coverages[i] > 1.5 * estimatedSingleCopyCoverage) continue;
		if (triplets.count(i) == 1) continue;
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

phmap::flat_hash_set<uint64_t> getTouchedUniqueTipsAfterCut(const ChunkUnitigGraph& graph, const std::vector<bool>& oldUniques, const uint64_t startNode)
{
	phmap::flat_hash_set<uint64_t> visited;
	phmap::flat_hash_set<uint64_t> result;
	std::vector<uint64_t> checkStack;
	checkStack.emplace_back(startNode);
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (visited.count(top) == 1) continue;
		visited.insert(top);
		if (oldUniques[top & maskUint64_t]) result.insert(top);
		if (!oldUniques[top & maskUint64_t] && (top & maskUint64_t) != (startNode & maskUint64_t)) checkStack.emplace_back(top ^ firstBitUint64_t);
		for (auto edge : graph.edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			checkStack.emplace_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	return result;
}

bool removingNodeCutsUniques(const ChunkUnitigGraph& graph, const std::vector<bool>& oldUniques, const size_t node)
{
	phmap::flat_hash_set<uint64_t> forwardUniquesAfterCut = getTouchedUniqueTipsAfterCut(graph, oldUniques, node + firstBitUint64_t);
	phmap::flat_hash_set<uint64_t> backwardUniquesAfterCut = getTouchedUniqueTipsAfterCut(graph, oldUniques, node);
	if (forwardUniquesAfterCut.size() == 0) return false;
	if (backwardUniquesAfterCut.size() == 0) return false;
	if (forwardUniquesAfterCut.size() != 1 && backwardUniquesAfterCut.size() != 1) return false;
	return true;
}

std::vector<bool> getTopologyExtendedUniques(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double estimatedSingleCopyCoverage, const std::vector<bool>& oldUniques)
{
	phmap::flat_hash_set<size_t> newUniques;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (oldUniques[i]) continue;
		if (!removingNodeCutsUniques(graph, oldUniques, i)) continue;
		newUniques.insert(i);
	}
	std::vector<bool> result = oldUniques;
	for (auto i : newUniques)
	{
		assert(!result[i]);
		result[i] = true;
	}
	if (newUniques.size() >= 1)
	{
		return getTopologyExtendedUniques(graph, readPaths, estimatedSingleCopyCoverage, result);
	}
	return result;
}

std::vector<std::vector<uint64_t>> connectTripletPaths(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double estimatedSingleCopyCoverage, const std::vector<std::vector<uint64_t>>& oldPaths)
{
	phmap::flat_hash_map<uint64_t, uint64_t> nodeIsPathTip;
	for (size_t i = 0; i < oldPaths.size(); i++)
	{
		assert(nodeIsPathTip.count(oldPaths[i].back()) == 0);
		assert(nodeIsPathTip.count(oldPaths[i][0] ^ firstBitUint64_t) == 0);
		nodeIsPathTip[oldPaths[i].back()] = i + firstBitUint64_t;
		nodeIsPathTip[oldPaths[i][0] ^ firstBitUint64_t] = i;
	}
	auto triplets = getNodeSpanningTriplets(graph, readPaths, estimatedSingleCopyCoverage);
	phmap::flat_hash_map<uint64_t, std::pair<uint64_t, uint64_t>> outEdges;
	for (const auto& pair : triplets)
	{
		bool allGood = true;
		for (const auto& pair2 : pair.second)
		{
			if (nodeIsPathTip.count(pair2.first) == 0)
			{
				allGood = false;
				break;
			}
			if (nodeIsPathTip.count(pair2.second ^ firstBitUint64_t) == 0)
			{
				allGood = false;
				break;
			}
			if ((nodeIsPathTip.at(pair2.first) & maskUint64_t) == (nodeIsPathTip.at(pair2.second ^ firstBitUint64_t) & maskUint64_t))
			{
				allGood = false;
				break;
			}
		}
		if (!allGood) continue;
		for (const auto& pair2 : pair.second)
		{
			assert(nodeIsPathTip.count(pair2.first) == 1);
			assert(nodeIsPathTip.count(pair2.second ^ firstBitUint64_t) == 1);
			uint64_t prevContig = nodeIsPathTip.at(pair2.first);
			uint64_t nextContig = nodeIsPathTip.at(pair2.second ^ firstBitUint64_t) ^ firstBitUint64_t;
			assert(outEdges.count(prevContig) == 0);
			assert(outEdges.count(nextContig ^ firstBitUint64_t) == 0);
			outEdges[prevContig] = std::make_pair(pair.first + firstBitUint64_t, nextContig);
			outEdges[nextContig ^ firstBitUint64_t] = std::make_pair(pair.first, prevContig ^ firstBitUint64_t);
		}
	}
	if (outEdges.size() == 0) return oldPaths;
	std::vector<bool> pathUsed;
	pathUsed.resize(oldPaths.size(), false);
	std::vector<std::vector<uint64_t>> result;
	for (size_t i = 0; i < oldPaths.size(); i++)
	{
		if (pathUsed[i]) continue;
		std::vector<uint64_t> bwPath;
		std::vector<uint64_t> fwPath;
		bwPath.emplace_back(i);
		while (outEdges.count(bwPath.back()) == 1)
		{
			bwPath.emplace_back(outEdges.at(bwPath.back()).second);
			assert(bwPath.size() <= oldPaths.size());
		}
		fwPath.emplace_back(i + firstBitUint64_t);
		while (outEdges.count(fwPath.back()) == 1)
		{
			fwPath.emplace_back(outEdges.at(fwPath.back()).second);
			assert(fwPath.size() <= oldPaths.size());
		}
		auto path = mergeBwFwPaths(bwPath, fwPath);
		assert(path.size() >= 1);
		for (size_t j = 0; j < path.size(); j++)
		{
			assert(!pathUsed[path[j] & maskUint64_t]);
			pathUsed[path[j] & maskUint64_t] = true;
		}
		result.emplace_back();
		for (size_t j = 0; j < path.size(); j++)
		{
			std::vector<uint64_t> addedHere = oldPaths[path[j] & maskUint64_t];
			if (path[j] & firstBitUint64_t)
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
			if (j > 0)
			{
				assert(outEdges.count(path[j-1]) == 1);
				result.back().emplace_back(outEdges.at(path[j-1]).first);
			}
			result.back().insert(result.back().end(), addedHere.begin(), addedHere.end());
		}
	}
	return result;
}

std::vector<bool> getTripletExtendedUniques(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double estimatedSingleCopyCoverage, const std::vector<bool>& oldUniques)
{
	auto triplets = getNodeSpanningTriplets(graph, readPaths, estimatedSingleCopyCoverage);
	std::vector<bool> result = oldUniques;
	for (const auto& pair : triplets)
	{
		if (graph.coverages[pair.first] < ((double)estimatedSingleCopyCoverage - 0.5) * pair.second.size()) continue;
		if (graph.coverages[pair.first] > ((double)estimatedSingleCopyCoverage + 0.5) * pair.second.size()) continue;
		bool allGood = true;
		for (const auto& pair2 : pair.second)
		{
			if (graph.coverages[pair2.first & maskUint64_t] > estimatedSingleCopyCoverage * 1.5) allGood = false;
			if (graph.coverages[pair2.second & maskUint64_t] > estimatedSingleCopyCoverage * 1.5) allGood = false;
		}
		if (!allGood) continue;
		double minFwCoverage = 10 * estimatedSingleCopyCoverage;
		double maxFwCoverage = 0;
		double minBwCoverage = 10 * estimatedSingleCopyCoverage;
		double maxBwCoverage = 0;
		for (const auto& pair2 : pair.second)
		{
			minBwCoverage = std::min(minBwCoverage, graph.coverages[pair2.first & maskUint64_t]);
			maxBwCoverage = std::max(maxBwCoverage, graph.coverages[pair2.first & maskUint64_t]);
			minFwCoverage = std::min(minFwCoverage, graph.coverages[pair2.second & maskUint64_t]);
			maxFwCoverage = std::max(maxFwCoverage, graph.coverages[pair2.second & maskUint64_t]);
		}
		bool bwValid = (maxBwCoverage < minBwCoverage * 2);
		bool fwValid = (maxFwCoverage < minFwCoverage * 2);
		for (const auto& pair2 : pair.second)
		{
			if (bwValid) result[pair2.first & maskUint64_t] = true;
			if (fwValid) result[pair2.second & maskUint64_t] = true;
		}
	}
	return result;
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
			auto path = mergeBwFwPaths(bwPath, fwPath);
			for (auto node : path)
			{
				alreadyInPath[node & maskUint64_t] = true;
			}
			size_t pathLength = getPathLength(graph, path);
			if (pathLength < 100000) continue;
			result.emplace_back(path);
		}
	}
	return result;
}

uint64_t getNextTripletCore(const ChunkUnitigGraph& graph, const phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>>& triplets, const uint64_t startNode)
{
	uint64_t result = std::numeric_limits<uint64_t>::max();
	for (auto pair : triplets.at(startNode & maskUint64_t))
	{
		uint64_t resultHere;
		if (startNode & firstBitUint64_t)
		{
			std::pair<size_t, bool> check { pair.second & maskUint64_t, pair.second & firstBitUint64_t };
			if (graph.edges.getEdges(check).size() != 1) return std::numeric_limits<uint64_t>::max();
			if (triplets.count(graph.edges.getEdges(check)[0].first) == 0) return std::numeric_limits<uint64_t>::max();
			resultHere = graph.edges.getEdges(check)[0].first + (graph.edges.getEdges(check)[0].second ? firstBitUint64_t : 0);
		}
		else
		{
			std::pair<size_t, bool> check { pair.first & maskUint64_t, (pair.first ^ firstBitUint64_t) & firstBitUint64_t };
			if (graph.edges.getEdges(check).size() != 1) return std::numeric_limits<uint64_t>::max();
			if (triplets.count(graph.edges.getEdges(check)[0].first) == 0) return std::numeric_limits<uint64_t>::max();
			resultHere = graph.edges.getEdges(check)[0].first + (graph.edges.getEdges(check)[0].second ? firstBitUint64_t : 0);
		}
		if (result == std::numeric_limits<uint64_t>::max()) result = resultHere;
		if (resultHere != result) return std::numeric_limits<uint64_t>::max();
	}
	return result;
}

std::vector<uint64_t> followTripletCorePath(const ChunkUnitigGraph& graph, const phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>>& triplets, const uint64_t startNode)
{
	assert(triplets.count(startNode & maskUint64_t) == 1);
	std::vector<uint64_t> result;
	result.emplace_back(startNode);
	while (true)
	{
		uint64_t next = getNextTripletCore(graph, triplets, result.back());
		if (next == std::numeric_limits<uint64_t>::max()) break;
		if (next == (result.back() ^ firstBitUint64_t)) break;
		result.emplace_back(next);
		if (next == startNode) break;
	}
	return result;
}

std::vector<bool> getTripletPathExtendedUniques(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double estimatedSingleCopyCoverage, const std::vector<bool>& oldUniques)
{
	auto triplets = getNodeSpanningTriplets(graph, readPaths, estimatedSingleCopyCoverage);
	std::vector<std::vector<uint64_t>> tripletChains;
	phmap::flat_hash_set<size_t> includedNodes;
	auto result = oldUniques;
	for (const auto& pair : triplets)
	{
		if (includedNodes.count(pair.first) == 1) continue;
		std::vector<uint64_t> fwCores = followTripletCorePath(graph, triplets, pair.first + firstBitUint64_t);
		std::vector<uint64_t> bwCores = followTripletCorePath(graph, triplets, pair.first);
		auto path = mergeBwFwPaths(bwCores, fwCores);
		if (path.size() == 1) continue;
		size_t length = 0;
		size_t ploidy = pair.second.size();
		bool hasUniqueCore = false;
		double coverageSum = 0;
		double coverageDivisor = 0;
		for (uint64_t node : path)
		{
			length += graph.unitigLengths[node & maskUint64_t];
			assert(triplets.count(node & maskUint64_t) == 1);
			assert(triplets.at(node & maskUint64_t).size() == ploidy);
			includedNodes.insert(node & maskUint64_t);
			if (oldUniques[node & maskUint64_t]) hasUniqueCore = true;
			coverageSum += graph.unitigLengths[node & maskUint64_t] * graph.coverages[node & maskUint64_t];
			coverageDivisor += graph.unitigLengths[node & maskUint64_t];
		}
		if (coverageSum / coverageDivisor < estimatedSingleCopyCoverage * 1.5 && hasUniqueCore) continue;
		bool alreadyHasUnique = false;
		for (uint64_t node : path)
		{
			bool allUniqueBw = true;
			bool allUniqueFw = true;
			for (auto pair2 : triplets.at(node & maskUint64_t))
			{
				if (!oldUniques[pair2.first & maskUint64_t]) allUniqueBw = false;
				if (!oldUniques[pair2.second & maskUint64_t]) allUniqueFw = false;
			}
			if (allUniqueBw || allUniqueFw) alreadyHasUnique = true;
		}
		if (length < 50000 && !alreadyHasUnique) continue;
		for (uint64_t node : path)
		{
			double minBwCoverage = estimatedSingleCopyCoverage * 10;
			double maxBwCoverage = 0;
			double minFwCoverage = estimatedSingleCopyCoverage * 10;
			double maxFwCoverage = 0;
			for (auto pair2 : triplets.at(node & maskUint64_t))
			{
				minBwCoverage = std::min(minBwCoverage, graph.coverages[pair2.first & maskUint64_t]);
				maxBwCoverage = std::max(maxBwCoverage, graph.coverages[pair2.first & maskUint64_t]);
				minFwCoverage = std::min(minFwCoverage, graph.coverages[pair2.second & maskUint64_t]);
				maxFwCoverage = std::max(maxFwCoverage, graph.coverages[pair2.second & maskUint64_t]);
			}
			bool bwValid = (maxBwCoverage < minBwCoverage * 2);
			bool fwValid = (maxFwCoverage < minFwCoverage * 2);
			for (auto pair2 : triplets.at(node & maskUint64_t))
			{
				if (bwValid) result[pair2.first & maskUint64_t] = true;
				if (fwValid) result[pair2.second & maskUint64_t] = true;
			}
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
		auto path = mergeBwFwPaths(bwPath, fwPath);
		for (size_t j = 0; j < path.size(); j++)
		{
			assert(!contigPathInPath[path[j] & maskUint64_t]);
			contigPathInPath[path[j] & maskUint64_t] = true;
		}
		assert(path.size() >= 1);
		newContigPathPaths.emplace_back(path);
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

std::pair<std::vector<std::vector<uint64_t>>, std::vector<std::vector<std::pair<size_t, size_t>>>> getPathChunkRepresentation(const ChunkUnitigGraph& graph, const std::vector<std::vector<uint64_t>>& paths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<std::vector<uint64_t>> contigPathChunks;
	std::vector<std::vector<std::pair<size_t, size_t>>> approxChunkPositionInContigPath;
	for (size_t i = 0; i < paths.size(); i++)
	{
		contigPathChunks.emplace_back();
		for (size_t j = 0; j < paths[i].size(); j++)
		{
			std::vector<uint64_t> addHere = graph.chunksInUnitig[paths[i][j] & maskUint64_t];
			if (paths[i][j] & firstBitUint64_t)
			{
			}
			else
			{
				std::reverse(addHere.begin(), addHere.end());
				for (size_t k = 0; k < addHere.size(); k++)
				{
					addHere[k] ^= firstBitUint64_t;
				}
			}
			assert(addHere.size() >= 1);
			contigPathChunks[i].insert(contigPathChunks[i].end(), addHere.begin(), addHere.end());
		}
	}
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
	for (size_t i = 0; i < contigPathChunks.size(); i++)
	{
		approxChunkPositionInContigPath.emplace_back();
		size_t endPos = 0;
		for (size_t j = 0; j < contigPathChunks[i].size(); j++)
		{
			size_t node = contigPathChunks[i][j] & maskUint64_t;
			size_t size = lengths[node][lengths[node].size()/2];
			endPos += size;
			if (j >= 1)
			{
				uint64_t from = contigPathChunks[i][j-1];
				uint64_t to = contigPathChunks[i][j];
				auto key = canonNodePair(from, to);
				assert(edgeOverlaps.count(key) == 1);
				size_t overlap = edgeOverlaps.at(key);
				if (overlap > size) overlap = size;
				size_t prevNode = contigPathChunks[i][j-1] & maskUint64_t;
				size_t prevSize = lengths[prevNode][lengths[prevNode].size()/2];
				if (overlap > prevSize) overlap = prevSize;
				assert(endPos >= overlap);
				endPos -= overlap;
			}
			assert(endPos >= size);
			approxChunkPositionInContigPath[i].emplace_back(endPos-size, endPos);
		}
	}
	return std::make_pair(contigPathChunks, approxChunkPositionInContigPath);
}

std::vector<std::vector<std::tuple<size_t, bool, size_t, size_t>>> getReadToContigPathAnchors(const std::vector<std::vector<uint64_t>>& contigPathChunks, const std::vector<std::vector<std::pair<size_t, size_t>>>& approxChunkPositionInContigPath, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> chunkPositionsInContigPaths; // contig, index, fw
	for (size_t i = 0; i < contigPathChunks.size(); i++)
	{
		for (size_t j = 0; j < contigPathChunks[i].size(); j++)
		{
			uint64_t node = contigPathChunks[i][j];
			while (chunkPositionsInContigPaths.size() <= (node & maskUint64_t)) chunkPositionsInContigPaths.emplace_back();
			chunkPositionsInContigPaths[node & maskUint64_t].emplace_back(i, j, node & firstBitUint64_t);
		}
	}
	std::vector<std::vector<std::tuple<size_t, bool, size_t, size_t>>> readToContigPathAnchors; // read, fw, readpos, contigpos
	readToContigPathAnchors.resize(contigPathChunks.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<std::tuple<size_t, bool, int, size_t, size_t, size_t, size_t>> readToChunkAnchors; // contigpath, fw, diagonal, readstart, readend, contigstart, contigend
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunkPositionsInContigPaths.size()) continue;
			bool readFw = std::get<2>(t) & firstBitUint64_t;
			size_t readPos;
			if (readFw)
			{
				readPos = std::get<0>(t);
			}
			else
			{
				readPos = std::get<1>(t);
			}
			for (auto t2 : chunkPositionsInContigPaths[std::get<2>(t) & maskUint64_t])
			{
				size_t contigPos;
				size_t contigChunkStart;
				size_t contigChunkEnd;
				if (std::get<2>(t2))
				{
					contigPos = approxChunkPositionInContigPath[std::get<0>(t2)][std::get<1>(t2)].first;
				}
				else
				{
					contigPos = approxChunkPositionInContigPath[std::get<0>(t2)][std::get<1>(t2)].second;
				}
				bool orientationMatch = !(readFw ^ std::get<2>(t2));
				if (orientationMatch)
				{
					contigChunkStart = approxChunkPositionInContigPath[std::get<0>(t2)][std::get<1>(t2)].first;
					contigChunkEnd = approxChunkPositionInContigPath[std::get<0>(t2)][std::get<1>(t2)].second;
				}
				else
				{
					contigChunkStart = approxChunkPositionInContigPath[std::get<0>(t2)][std::get<1>(t2)].second;
					contigChunkEnd = approxChunkPositionInContigPath[std::get<0>(t2)][std::get<1>(t2)].first;
				}
				readToChunkAnchors.emplace_back(std::get<0>(t2), orientationMatch, (int)readPos - (orientationMatch ? 1 : -1) * (int)contigPos, std::get<0>(t), std::get<1>(t), contigChunkStart, contigChunkEnd);
			}
		}
		std::sort(readToChunkAnchors.begin(), readToChunkAnchors.end());
		std::vector<std::tuple<size_t, bool, size_t, size_t, size_t, std::vector<std::pair<size_t, size_t>>>> unfilteredClusters; // path, fw, readclusterstart, readclusterend, matchbp, [readpos, pathpos]
		size_t clusterStart = 0;
		for (size_t j = 1; j <= readToChunkAnchors.size(); j++)
		{
			assert(j == readToChunkAnchors.size() || std::get<0>(readToChunkAnchors[j]) >= std::get<0>(readToChunkAnchors[j-1]));
			assert(j == readToChunkAnchors.size() || std::get<1>(readToChunkAnchors[j]) >= std::get<1>(readToChunkAnchors[j-1]) || std::get<0>(readToChunkAnchors[j]) > std::get<0>(readToChunkAnchors[j-1]));
			assert(j == readToChunkAnchors.size() || std::get<2>(readToChunkAnchors[j]) >= std::get<2>(readToChunkAnchors[j-1]) || std::get<1>(readToChunkAnchors[j]) > std::get<1>(readToChunkAnchors[j-1]) || std::get<0>(readToChunkAnchors[j]) > std::get<0>(readToChunkAnchors[j-1]));
			if (j == readToChunkAnchors.size() || std::get<0>(readToChunkAnchors[j]) != std::get<0>(readToChunkAnchors[j-1]) || std::get<1>(readToChunkAnchors[j]) != std::get<1>(readToChunkAnchors[j-1]) || std::get<2>(readToChunkAnchors[j]) > std::get<2>(readToChunkAnchors[j-1]) + 100)
			{
				std::vector<std::pair<size_t, size_t>> matchPositions;
				size_t clusterStartPosInRead = std::get<3>(readToChunkAnchors[clusterStart]);
				size_t clusterEndPosInRead = std::get<4>(readToChunkAnchors[clusterStart]);
				for (size_t k = clusterStart; k < j; k++)
				{
					clusterStartPosInRead = std::min(clusterStartPosInRead, std::get<3>(readToChunkAnchors[k]));
					clusterEndPosInRead = std::max(clusterEndPosInRead, std::get<4>(readToChunkAnchors[k]));
					matchPositions.emplace_back(std::get<3>(readToChunkAnchors[k]), std::get<4>(readToChunkAnchors[k]));
				}
				std::sort(matchPositions.begin(), matchPositions.end());
				size_t totalMatches = 0;
				size_t lastMatchEnd = 0;
				for (size_t k = 0; k < matchPositions.size(); k++)
				{
					if (matchPositions[k].first > lastMatchEnd)
					{
						totalMatches += matchPositions[k].second - matchPositions[k].first;
						lastMatchEnd = matchPositions[k].second;
					}
					else
					{
						if (matchPositions[k].second > lastMatchEnd)
						{
							totalMatches += matchPositions[k].second - lastMatchEnd;
							lastMatchEnd = matchPositions[k].second;
						}
					}
				}
				assert(totalMatches > 0);
				std::vector<std::pair<size_t, size_t>> matchPoses;
				for (size_t k = clusterStart; k < j; k++)
				{
					matchPoses.emplace_back(std::get<3>(readToChunkAnchors[k]), std::get<5>(readToChunkAnchors[k]));
					matchPoses.emplace_back(std::get<4>(readToChunkAnchors[k]), std::get<6>(readToChunkAnchors[k]));
				}
				assert(matchPoses.size() > 0);
				unfilteredClusters.emplace_back(std::get<0>(readToChunkAnchors[clusterStart]), std::get<1>(readToChunkAnchors[clusterStart]), clusterStartPosInRead, clusterEndPosInRead, totalMatches, matchPoses);
				clusterStart = j;
			}
		}
		std::vector<bool> removeCluster;
		removeCluster.resize(unfilteredClusters.size(), false);
		for (size_t j = 1; j < unfilteredClusters.size(); j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				if (std::get<2>(unfilteredClusters[j]) >= std::get<3>(unfilteredClusters[k])) continue;
				if (std::get<2>(unfilteredClusters[k]) >= std::get<3>(unfilteredClusters[j])) continue;
				if (std::get<4>(unfilteredClusters[j]) > std::get<4>(unfilteredClusters[k]))
				{
					removeCluster[k] = true;
				}
				if (std::get<4>(unfilteredClusters[k]) > std::get<4>(unfilteredClusters[j]))
				{
					removeCluster[j] = true;
				}
			}
		}
		for (size_t j = unfilteredClusters.size()-1; j < unfilteredClusters.size(); j--)
		{
			if (!removeCluster[j]) continue;
			std::swap(unfilteredClusters[j], unfilteredClusters.back());
			unfilteredClusters.pop_back();
		}
		for (size_t j = 0; j < unfilteredClusters.size(); j++)
		{
			for (size_t k = 0; k < std::get<5>(unfilteredClusters[j]).size(); k++)
			{
				readToContigPathAnchors[std::get<0>(unfilteredClusters[j])].emplace_back(i, std::get<1>(unfilteredClusters[j]), std::get<5>(unfilteredClusters[j])[k].first, std::get<5>(unfilteredClusters[j])[k].second);
			}
		}
	}
	for (size_t i = 0; i < readToContigPathAnchors.size(); i++)
	{
		std::sort(readToContigPathAnchors[i].begin(), readToContigPathAnchors[i].end());
	}
	return readToContigPathAnchors;
}

std::vector<std::vector<std::vector<std::tuple<size_t, size_t, bool>>>> getAnchorPaths(const std::vector<std::vector<std::tuple<size_t, bool, size_t, size_t>>>& readToContigPathAnchors)
{
	std::vector<std::vector<std::vector<std::tuple<size_t, size_t, bool>>>> result;
	for (size_t i = 0; i < readToContigPathAnchors.size(); i++)
	{
		if (readToContigPathAnchors[i].size() == 0)
		{
			result.emplace_back();
			continue;
		}
		phmap::flat_hash_map<std::tuple<size_t, size_t, bool>, std::tuple<size_t, size_t, bool>> anchorParent;
		phmap::flat_hash_map<size_t, std::vector<std::tuple<size_t, size_t, bool>>> anchorsPerPosition;
		for (auto t : readToContigPathAnchors[i])
		{
			auto key = std::make_tuple(std::get<0>(t), std::get<2>(t), std::get<1>(t));
			anchorParent[key] = key;
			anchorsPerPosition[std::get<3>(t)].emplace_back(key);
		}
		for (const auto& pair : anchorsPerPosition)
		{
			for (auto t : pair.second)
			{
				merge(anchorParent, t, *pair.second.begin());
			}
		}
		phmap::flat_hash_map<std::tuple<size_t, size_t, bool>, size_t> anchorToCluster;
		for (auto t : readToContigPathAnchors[i])
		{
			auto key = std::make_tuple(std::get<0>(t), std::get<2>(t), std::get<1>(t));
			auto p = find(anchorParent, key);
			if (anchorToCluster.count(p) == 1) continue;
			size_t nextNum = anchorToCluster.size();
			anchorToCluster[p] = nextNum;
		}
		std::vector<phmap::flat_hash_map<size_t, size_t>> outEdges;
		outEdges.resize(anchorToCluster.size());
		phmap::flat_hash_map<std::pair<size_t, bool>, phmap::flat_hash_set<size_t>> anchorsPerRead;
		for (auto t : readToContigPathAnchors[i])
		{
			anchorsPerRead[std::make_pair(std::get<0>(t), std::get<1>(t))].insert(std::get<2>(t));
		}
		phmap::flat_hash_set<size_t> forbiddenClusters;
		while (true)
		{
			bool removedAny = false;
			for (const auto& pair : anchorsPerRead)
			{
				std::vector<size_t> anchorsSorted { pair.second.begin(), pair.second.end() };
				for (size_t j = anchorsSorted.size()-1; j < anchorsSorted.size(); j--)
				{
					auto key = std::make_tuple(pair.first.first, anchorsSorted[j], pair.first.second);
					if (forbiddenClusters.count(anchorToCluster.at(find(anchorParent, key))) == 0) continue;
					std::swap(anchorsSorted[j], anchorsSorted.back());
					anchorsSorted.pop_back();
				}
				std::sort(anchorsSorted.begin(), anchorsSorted.end());
				if (pair.first.second)
				{
					for (size_t j = 1; j < anchorsSorted.size(); j++)
					{
						auto prevKey = std::make_tuple(pair.first.first, anchorsSorted[j-1], pair.first.second);
						auto thisKey = std::make_tuple(pair.first.first, anchorsSorted[j], pair.first.second);
						if (anchorToCluster.at(find(anchorParent, prevKey)) == anchorToCluster.at(find(anchorParent, thisKey)))
						{
							forbiddenClusters.insert(anchorToCluster.at(find(anchorParent, prevKey)));
							removedAny = true;
						}
					}
				}
				else
				{
					for (size_t j = 1; j < anchorsSorted.size(); j++)
					{
						auto prevKey = std::make_tuple(pair.first.first, anchorsSorted[j], pair.first.second);
						auto thisKey = std::make_tuple(pair.first.first, anchorsSorted[j-1], pair.first.second);
						if (anchorToCluster.at(find(anchorParent, prevKey)) == anchorToCluster.at(find(anchorParent, thisKey)))
						{
							forbiddenClusters.insert(anchorToCluster.at(find(anchorParent, prevKey)));
							removedAny = true;
						}
					}
				}
			}
			if (!removedAny) break;
		}
		for (const auto& pair : anchorsPerRead)
		{
			std::vector<size_t> anchorsSorted { pair.second.begin(), pair.second.end() };
			for (size_t j = anchorsSorted.size()-1; j < anchorsSorted.size(); j--)
			{
				auto key = std::make_tuple(pair.first.first, anchorsSorted[j], pair.first.second);
				if (forbiddenClusters.count(anchorToCluster.at(find(anchorParent, key))) == 0) continue;
				std::swap(anchorsSorted[j], anchorsSorted.back());
				anchorsSorted.pop_back();
			}
			std::sort(anchorsSorted.begin(), anchorsSorted.end());
			if (pair.first.second)
			{
				for (size_t j = 1; j < anchorsSorted.size(); j++)
				{
					auto prevKey = std::make_tuple(pair.first.first, anchorsSorted[j-1], pair.first.second);
					auto thisKey = std::make_tuple(pair.first.first, anchorsSorted[j], pair.first.second);
					outEdges[anchorToCluster.at(find(anchorParent, prevKey))][anchorToCluster.at(find(anchorParent, thisKey))] += 1;
				}
			}
			else
			{
				for (size_t j = 1; j < anchorsSorted.size(); j++)
				{
					auto prevKey = std::make_tuple(pair.first.first, anchorsSorted[j], pair.first.second);
					auto thisKey = std::make_tuple(pair.first.first, anchorsSorted[j-1], pair.first.second);
					outEdges[anchorToCluster.at(find(anchorParent, prevKey))][anchorToCluster.at(find(anchorParent, thisKey))] += 1;
				}
			}
		}
		std::vector<phmap::flat_hash_map<size_t, size_t>> inEdges;
		inEdges.resize(outEdges.size());
		for (size_t from = 0; from < outEdges.size(); from++)
		{
			for (auto pair : outEdges[from])
			{
				inEdges[pair.first][from] = pair.second;
			}
		}
		std::vector<size_t> anchorSequence;
		size_t startPos = 0;
		size_t iterations = 0;
		while (inEdges[startPos].size() >= 1)
		{
			size_t bestEdge = inEdges[startPos].begin()->first;
			for (auto pair : inEdges[startPos])
			{
				if (pair.second > inEdges[startPos].at(bestEdge)) bestEdge = pair.first;
			}
			startPos = bestEdge;
			assert(iterations < anchorToCluster.size()*2+10); // loops!
			iterations += 1;
		}
		anchorSequence.emplace_back(startPos);
		while (outEdges[startPos].size() >= 1)
		{
			size_t bestEdge = outEdges[startPos].begin()->first;
			for (auto pair : outEdges[startPos])
			{
				if (pair.second > outEdges[startPos].at(bestEdge)) bestEdge = pair.first;
			}
			startPos = bestEdge;
			anchorSequence.emplace_back(startPos);
			assert(iterations < anchorToCluster.size()*2+10); // loops!
			iterations += 1;
		}
		phmap::flat_hash_map<size_t, size_t> nodePositionInAnchorSequence;
		for (size_t j = 0; j < anchorSequence.size(); j++)
		{
			nodePositionInAnchorSequence[anchorSequence[j]] = j;
		}
		assert(anchorSequence.size() >= 2);
		result.emplace_back();
		result.back().resize(anchorSequence.size());
		for (auto t : readToContigPathAnchors[i])
		{
			auto key = std::make_tuple(std::get<0>(t), std::get<2>(t), std::get<1>(t));
			size_t cluster = anchorToCluster.at(find(anchorParent, key));
			if (nodePositionInAnchorSequence.count(cluster) == 0) continue;
			assert(nodePositionInAnchorSequence.at(cluster) < result.back().size());
			result.back()[nodePositionInAnchorSequence[cluster]].emplace_back(std::get<0>(t), std::get<2>(t), std::get<1>(t));
		}
		for (size_t j = 0; j < result.back().size(); j++)
		{
			assert(result.back()[j].size() >= 1);
		}
	}
	return result;
}

ConsensusString getAnchorPathConsensus(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& anchorPaths)
{
	assert(anchorPaths.size() >= 2);
	ConsensusString result;
	for (size_t i = 1; i < anchorPaths.size(); i++)
	{
		phmap::flat_hash_map<std::string, size_t> stringsHere;
		size_t totalCount = 0;
		for (auto t : anchorPaths[i-1])
		{
			for (auto t2 : anchorPaths[i])
			{
				if (std::get<0>(t) != std::get<0>(t2)) continue;
				if (std::get<2>(t) != std::get<2>(t2)) continue;
				if (std::get<2>(t))
				{
					if (std::get<1>(t) >= std::get<1>(t2)) continue;
					std::string seq = sequenceIndex.getSubstring(std::get<0>(t), std::get<1>(t), std::get<1>(t2) - std::get<1>(t) + 1);
					stringsHere[seq] += 1;
					totalCount += 1;
				}
				else
				{
					if (std::get<1>(t) <= std::get<1>(t2)) continue;
					std::string seq = sequenceIndex.getSubstring(std::get<0>(t), std::get<1>(t2), std::get<1>(t) - std::get<1>(t2) + 1);
					seq = revCompRaw(seq);
					stringsHere[seq] += 1;
					totalCount += 1;
				}
			}
		}
		if (stringsHere.size() == 0)
		{
			result.Nchunks.emplace_back(result.bases.size(), result.bases.size()+2);
			result.bases.emplace_back(0);
			result.bases.emplace_back(0);
			continue;
		}
		std::string consensusString;
		assert(stringsHere.size() <= totalCount);
		if (stringsHere.size() == 1 || totalCount == 2)
		{
			consensusString = stringsHere.begin()->first;
		}
		else
		{
			consensusString = getConsensus(stringsHere, totalCount);
		}
		size_t startPos = 0;
		if (i >= 2 && (result.Nchunks.size() == 0 || result.Nchunks.back().second < result.bases.size()))
		{
			startPos = 1;
		}
		for (size_t j = startPos; j < consensusString.size(); j++)
		{
			switch(consensusString[j])
			{
			case 'A':
				result.bases.emplace_back(0);
				break;
			case 'C':
				result.bases.emplace_back(1);
				break;
			case 'G':
				result.bases.emplace_back(2);
				break;
			case 'T':
				result.bases.emplace_back(3);
				break;
			default:
				assert(false);
				break;
			}
		}
	}
	return result;
}

std::vector<ConsensusString> getAnchorPathConsensuses(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<std::vector<std::tuple<size_t, size_t, bool>>>>& anchorPaths, const size_t numThreads)
{
	std::vector<ConsensusString> result;
	result.resize(anchorPaths.size());
	std::mutex printMutex;
	iterateMultithreaded(0, anchorPaths.size(), numThreads, [&sequenceIndex, &result, &anchorPaths, &printMutex](size_t i)
	{
		if (anchorPaths[i].size() == 0)
		{
			std::lock_guard<std::mutex> lock { printMutex };
			std::cerr << "skip consensus of path " << i << " due to no anchors" << std::endl;
			return;
		}
		{
			std::lock_guard<std::mutex> lock { printMutex };
			std::cerr << "start consensus of path " << i << std::endl;
		}
		result[i] = getAnchorPathConsensus(sequenceIndex, anchorPaths[i]);
		{
			std::lock_guard<std::mutex> lock { printMutex };
			std::cerr << "did get consensus of path " << i << ", length " << result[i].bases.size() << " bp" << std::endl;
		}
	});
	return result;
}

std::vector<ConsensusString> getPathSequences(const FastaCompressor::CompressedStringIndex& sequenceIndex, const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<uint64_t>>& paths, const size_t numThreads)
{
	std::vector<std::vector<uint64_t>> contigPathChunks;
	std::vector<std::vector<std::pair<size_t, size_t>>> approxChunkPositionInContigPath;
	std::tie(contigPathChunks, approxChunkPositionInContigPath) = getPathChunkRepresentation(graph, paths, chunksPerRead);
	std::vector<std::vector<std::tuple<size_t, bool, size_t, size_t>>> readToContigPathAnchors = getReadToContigPathAnchors(contigPathChunks, approxChunkPositionInContigPath, chunksPerRead);
	std::vector<std::vector<std::vector<std::tuple<size_t, size_t, bool>>>> anchorPaths = getAnchorPaths(readToContigPathAnchors);
	std::vector<ConsensusString> result = getAnchorPathConsensuses(sequenceIndex, anchorPaths, numThreads);
	return result;
}

void getContigPathsAndConsensuses(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t kmerSize, const std::string& outPathsFile, const std::string& outContigsFasta, const size_t numThreads)
{
	const size_t longUnitigThreshold = 100000;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	auto fixedChunksPerRead = getBidirectedChunks(chunksPerRead);
	std::tie(graph, readPaths) = getChunkUnitigGraph(fixedChunksPerRead, approxOneHapCoverage, kmerSize);
	double estimatedSingleCopyCoverage = getAverageLongNodeCoverage(graph, longUnitigThreshold);
	std::cerr << "estimated single copy coverage " << estimatedSingleCopyCoverage << std::endl;
	std::vector<bool> isUniqueUnitig = estimateUniqueUnitigs(graph, readPaths, longUnitigThreshold, estimatedSingleCopyCoverage);
	std::cerr << "step 1 unique unitigs: ";
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (!isUniqueUnitig[i]) continue;
		std::cerr << i << ",";
	}
	std::cerr << std::endl;
	isUniqueUnitig = getTripletExtendedUniques(graph, readPaths, estimatedSingleCopyCoverage, isUniqueUnitig);
	std::cerr << "step 1.5 unique unitigs: ";
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (!isUniqueUnitig[i]) continue;
		std::cerr << i << ",";
	}
	std::cerr << std::endl;
	isUniqueUnitig = getTripletPathExtendedUniques(graph, readPaths, estimatedSingleCopyCoverage, isUniqueUnitig);
	std::cerr << "step 2 unique unitigs: ";
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (!isUniqueUnitig[i]) continue;
		std::cerr << i << ",";
	}
	std::cerr << std::endl;
	isUniqueUnitig = getTopologyExtendedUniques(graph, readPaths, estimatedSingleCopyCoverage, isUniqueUnitig);
	std::cerr << "step 3 unique unitigs: ";
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (!isUniqueUnitig[i]) continue;
		std::cerr << i << ",";
	}
	std::cerr << std::endl;
	std::vector<std::vector<uint64_t>> paths;
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (!isUniqueUnitig[i]) continue;
		paths.emplace_back();
		paths.back().emplace_back(i + firstBitUint64_t);
	}
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
	paths = connectTripletPaths(graph, readPaths, estimatedSingleCopyCoverage, paths);
	std::cerr << "step 2.5 paths:" << std::endl;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::cerr << "step 2.5 path id " << i << " length " << getPathLength(graph, paths[i]) << ": ";
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
	{
		std::ofstream file { outPathsFile };
		for (size_t i = 0; i < paths.size(); i++)
		{
			file << i << "\t";
			for (auto node : paths[i])
			{
				file << ((node & firstBitUint64_t) ? ">" : "<") << (node & maskUint64_t);
			}
			file << std::endl;
		}
	}
	std::vector<ConsensusString> contigSequences = getPathSequences(sequenceIndex, graph, readPaths, fixedChunksPerRead, paths, numThreads);
	{
		std::ofstream file { outContigsFasta };
		for (size_t i = 0; i < paths.size(); i++)
		{
			file << ">contig_" << i << "_length_" << contigSequences[i].bases.size() << std::endl;
			std::string str = contigSequences[i].bases.toString();
			for (auto pair : contigSequences[i].Nchunks)
			{
				for (size_t j = pair.first; j < pair.first+pair.second; j++)
				{
					str[j] = 'N';
				}
			}
			file << str << std::endl;
		}
	}
}
