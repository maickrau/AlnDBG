#include <iostream>
#include <phmap.h>
#include "AnchorFinder.h"
#include "Common.h"
#include "MBGCommon.h"
#include "UnionFind.h"

std::tuple<uint64_t, uint64_t, uint64_t> canon(std::pair<size_t, bool> first, std::pair<size_t, bool> second, std::pair<size_t, bool> third)
{
	uint64_t one = first.first + (first.second ? firstBitUint64_t : 0);
	uint64_t two = second.first + (second.second ? firstBitUint64_t : 0);
	uint64_t three = third.first + (third.second ? firstBitUint64_t : 0);
	std::tuple<uint64_t, uint64_t, uint64_t> fw { one, two, three };
	std::tuple<uint64_t, uint64_t, uint64_t> bw { three ^ firstBitUint64_t, two ^ firstBitUint64_t, one ^ firstBitUint64_t };
	return std::min(fw, bw);
}

std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> canon(std::pair<size_t, bool> first, std::pair<size_t, bool> second, std::pair<size_t, bool> third, std::pair<size_t, bool> fourth)
{
	uint64_t one = first.first + (first.second ? firstBitUint64_t : 0);
	uint64_t two = second.first + (second.second ? firstBitUint64_t : 0);
	uint64_t three = third.first + (third.second ? firstBitUint64_t : 0);
	uint64_t four = fourth.first + (fourth.second ? firstBitUint64_t : 0);
	std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> fw { one, two, three, four };
	std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> bw { four ^ firstBitUint64_t, three ^ firstBitUint64_t, two ^ firstBitUint64_t, one ^ firstBitUint64_t };
	return std::min(fw, bw);
}

std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> getPossiblyTransitiveTriplets(const SparseEdgeContainer& rawEdgesWithTransitiveEdges)
{
	phmap::flat_hash_set<std::tuple<uint64_t, uint64_t, uint64_t>> result;
	for (size_t i = 0; i < rawEdgesWithTransitiveEdges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : rawEdgesWithTransitiveEdges.getEdges(fw))
		{
			for (auto edge2 : rawEdgesWithTransitiveEdges.getEdges(edge))
			{
				if (!rawEdgesWithTransitiveEdges.hasEdge(fw, edge2)) continue;
				result.emplace(canon(fw, edge, edge2));
			}
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : rawEdgesWithTransitiveEdges.getEdges(bw))
		{
			for (auto edge2 : rawEdgesWithTransitiveEdges.getEdges(edge))
			{
				if (!rawEdgesWithTransitiveEdges.hasEdge(bw, edge2)) continue;
				result.emplace(canon(bw, edge, edge2));
			}
		}
	}
	std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> resultvec { result.begin(), result.end() };
	return resultvec;
}

std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> getPossiblyTransitiveQuadruplets(const SparseEdgeContainer& rawEdgesWithTransitiveEdges)
{
	phmap::flat_hash_set<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> result;
	for (size_t i = 0; i < rawEdgesWithTransitiveEdges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : rawEdgesWithTransitiveEdges.getEdges(fw))
		{
			for (auto edge2 : rawEdgesWithTransitiveEdges.getEdges(edge))
			{
				for (auto edge3 : rawEdgesWithTransitiveEdges.getEdges(edge2))
				{
					if (!rawEdgesWithTransitiveEdges.hasEdge(fw, edge3)) continue;
					result.emplace(canon(fw, edge, edge2, edge3));
				}
			}
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : rawEdgesWithTransitiveEdges.getEdges(bw))
		{
			for (auto edge2 : rawEdgesWithTransitiveEdges.getEdges(edge))
			{
				for (auto edge3 : rawEdgesWithTransitiveEdges.getEdges(edge2))
				{
					if (!rawEdgesWithTransitiveEdges.hasEdge(bw, edge3)) continue;
					result.emplace(canon(bw, edge, edge2, edge3));
				}
			}
		}
	}
	std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> resultvec { result.begin(), result.end() };
	return resultvec;
}

bool isReachableWithoutNodeEvenBackwards(const uint64_t startNode, const uint64_t endNode, const uint64_t forbiddenNode, const SparseEdgeContainer& edges, const std::vector<bool>& anchor)
{
	phmap::flat_hash_set<uint64_t> checked;
	std::vector<uint64_t> stack;
	stack.push_back(startNode);
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (checked.count(top) == 1) continue;
		if ((top & maskUint64_t) == (endNode & maskUint64_t)) return true;
		checked.insert(top);
		if ((top & maskUint64_t) != startNode && (top & maskUint64_t) != (forbiddenNode & maskUint64_t) && !anchor[top & maskUint64_t]) stack.push_back(top ^ firstBitUint64_t);
		for (const auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	if (forbiddenNode != std::numeric_limits<size_t>::max())
	{
		assert(checked.count(forbiddenNode) == 1 || checked.count(forbiddenNode ^ firstBitUint64_t) == 1);
	}
	return false;
}

bool cutsGraph(const UnitigGraph& unitigGraph, const SparseEdgeContainer& edges, const uint64_t midNode, const std::vector<bool>& anchor)
{
	phmap::flat_hash_set<uint64_t> reachableFromBw;
	std::vector<uint64_t> checkStack;
	checkStack.emplace_back(midNode);
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (reachableFromBw.count(top) == 1) continue;
		if (top == midNode + firstBitUint64_t) return false;
		reachableFromBw.insert(top);
		if (!anchor[top & maskUint64_t] && (top & maskUint64_t) != midNode)
		{
			checkStack.emplace_back(top ^ firstBitUint64_t);
		}
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			checkStack.emplace_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	assert(reachableFromBw.count(midNode + firstBitUint64_t) == 0);
	return true;
}

bool blocksGraphPath(const UnitigGraph& unitigGraph, const SparseEdgeContainer& edges, const uint64_t startNode, const uint64_t midNode, const uint64_t endNode, const std::vector<bool>& anchor)
{
	if (!isReachableWithoutNodeEvenBackwards(startNode, endNode, std::numeric_limits<size_t>::max(), edges, anchor)) return false;
	if (isReachableWithoutNodeEvenBackwards(startNode, endNode, midNode, edges, anchor)) return false;
	return true;
}

SparseEdgeContainer getEdgesWithoutTransitiveEdges(const SparseEdgeContainer& rawEdgesWithTransitiveEdges, const std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>& maybeTransitiveTriplets, const std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& maybeTransitiveQuadruplets, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, std::vector<bool>& anchor, const size_t minCoverage)
{
	if (maybeTransitiveTriplets.size() == 0 && maybeTransitiveQuadruplets.size() == 0) return rawEdgesWithTransitiveEdges;
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> forbiddenTransitiveEdges;
	for (auto t : maybeTransitiveTriplets)
	{
		assert(anchor[std::get<1>(t) & maskUint64_t]);
		anchor[std::get<1>(t) & maskUint64_t] = false;
		if (!blocksGraphPath(unitigGraph, edges, std::get<0>(t), std::get<1>(t), std::get<2>(t), anchor))
		{
			anchor[std::get<1>(t) & maskUint64_t] = true;
			continue;
		}
		anchor[std::get<1>(t) & maskUint64_t] = true;
		forbiddenTransitiveEdges.emplace(std::get<0>(t), std::get<2>(t));
		std::cerr << "forbid transitive edge " << ((std::get<0>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<0>(t) & maskUint64_t) << " " << ((std::get<2>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(t) & maskUint64_t) << std::endl;
	}
	for (auto t : maybeTransitiveQuadruplets)
	{
		assert(anchor[std::get<1>(t) & maskUint64_t]);
		assert(anchor[std::get<2>(t) & maskUint64_t]);
		anchor[std::get<1>(t) & maskUint64_t] = false;
		anchor[std::get<2>(t) & maskUint64_t] = false;
		if (!blocksGraphPath(unitigGraph, edges, std::get<0>(t), std::get<1>(t), std::get<3>(t), anchor))
		{
			anchor[std::get<1>(t) & maskUint64_t] = true;
			anchor[std::get<2>(t) & maskUint64_t] = true;
			continue;
		}
		if (!blocksGraphPath(unitigGraph, edges, std::get<0>(t), std::get<2>(t), std::get<3>(t), anchor))
		{
			anchor[std::get<1>(t) & maskUint64_t] = true;
			anchor[std::get<2>(t) & maskUint64_t] = true;
			continue;
		}
		anchor[std::get<1>(t) & maskUint64_t] = true;
		anchor[std::get<2>(t) & maskUint64_t] = true;
		forbiddenTransitiveEdges.emplace(std::get<0>(t), std::get<3>(t));
		std::cerr << "forbid transitive edge " << ((std::get<0>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<0>(t) & maskUint64_t) << " " << ((std::get<3>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<3>(t) & maskUint64_t) << std::endl;
	}
	if (forbiddenTransitiveEdges.size() == 0) return rawEdgesWithTransitiveEdges;
	SparseEdgeContainer result;
	result.resize(rawEdgesWithTransitiveEdges.size());
	for (size_t i = 0; i < rawEdgesWithTransitiveEdges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (const auto edge : rawEdgesWithTransitiveEdges.getEdges(fw))
		{
			if (forbiddenTransitiveEdges.count(std::make_pair(fw.first + (fw.second ? firstBitUint64_t : 0), edge.first + (edge.second ? firstBitUint64_t : 0))) == 1) continue;
			if (forbiddenTransitiveEdges.count(std::make_pair(edge.first + (edge.second ? 0 : firstBitUint64_t), fw.first + (fw.second ? 0 : firstBitUint64_t))) == 1) continue;
			result.addEdge(fw, edge);
		}
		std::pair<size_t, bool> bw { i, false };
		for (const auto edge : rawEdgesWithTransitiveEdges.getEdges(bw))
		{
			if (forbiddenTransitiveEdges.count(std::make_pair(bw.first + (bw.second ? firstBitUint64_t : 0), edge.first + (edge.second ? firstBitUint64_t : 0))) == 1) continue;
			if (forbiddenTransitiveEdges.count(std::make_pair(edge.first + (edge.second ? 0 : firstBitUint64_t), bw.first + (bw.second ? 0 : firstBitUint64_t))) == 1) continue;
			result.addEdge(bw, edge);
		}
	}
	return result;
}

void addEdges(const std::pair<size_t, bool> start, SparseEdgeContainer& result, const SparseEdgeContainer& edges, const std::vector<bool>& anchor)
{
	assert(anchor[start.first]);
	phmap::flat_hash_set<std::pair<size_t, bool>> visited;
	std::vector<std::pair<size_t, bool>> stack;
	for (auto edge : edges.getEdges(start))
	{
		stack.push_back(edge);
	}
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		visited.insert(top);
		if (anchor[top.first])
		{
			result.addEdge(start, top);
			continue;
		}
		for (auto edge : edges.getEdges(top))
		{
			stack.push_back(edge);
		}
	}
}

SparseEdgeContainer getTopologicalEdges(const UnitigGraph& unitigGraph, const std::vector<bool>& anchor)
{
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	SparseEdgeContainer result;
	result.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < anchor.size(); i++)
	{
		if (!anchor[i]) continue;
		addEdges(std::make_pair(i, true), result, edges, anchor);
		addEdges(std::make_pair(i, false), result, edges, anchor);
	}
	return result;
}

SparseEdgeContainer getAnchorEdges(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<bool>& anchor, const size_t minCoverage)
{
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverage;
	edgeCoverage.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		uint64_t lastAnchor = std::numeric_limits<size_t>::max();
		size_t lastRightClip = 0;
		int lastAnchorImpliedStartPos = 0;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (!anchor[node & maskUint64_t]) continue;
				if (lastAnchor == std::numeric_limits<size_t>::max())
				{
					lastAnchor = node;
					lastAnchorImpliedStartPos = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t];
					if (k == readPaths[i].paths[j].path.size()-1)
					{
						lastRightClip = readPaths[i].paths[j].pathRightClipKmers;
					}
					else
					{
						lastRightClip = 0;
					}
					continue;
				}
				if (((int)readPos < lastAnchorImpliedStartPos + (int)unitigGraph.lengths[lastAnchor & maskUint64_t]) || ((int)readPos - (int)unitigGraph.lengths[node & maskUint64_t] < lastAnchorImpliedStartPos))
				{
					lastAnchor = node;
					lastAnchorImpliedStartPos = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t];
					if (k == readPaths[i].paths[j].path.size()-1)
					{
						lastRightClip = readPaths[i].paths[j].pathRightClipKmers;
					}
					else
					{
						lastRightClip = 0;
					}
					continue;
				}
				if (k == 0 && lastRightClip > 0 && lastAnchor == readPaths[i].paths[j].path[k])
				{
					if (unitigGraph.lengths[lastAnchor & maskUint64_t] - lastRightClip <= readPaths[i].paths[j].pathLeftClipKmers)
					{
						lastAnchor = node;
						lastAnchorImpliedStartPos = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t];
						if (k == readPaths[i].paths[j].path.size()-1)
						{
							lastRightClip = readPaths[i].paths[j].pathRightClipKmers;
						}
						else
						{
							lastRightClip = 0;
						}
						continue;
					}
				}
				std::pair<size_t, bool> fromEdge { lastAnchor & maskUint64_t, lastAnchor & firstBitUint64_t };
				std::pair<size_t, bool> toEdge { readPaths[i].paths[j].path[k] & maskUint64_t, readPaths[i].paths[j].path[k] & firstBitUint64_t };
				auto key = canon(fromEdge, toEdge);
				size_t coverage = 0;
				assert(key.first.first < edgeCoverage.size());
				assert(key.second.first < edgeCoverage.size());
				if (edgeCoverage.hasValue(key.first, key.second)) coverage = edgeCoverage.get(key.first, key.second);
				edgeCoverage.set(key.first, key.second, coverage+1);
				lastAnchor = node;
				lastAnchorImpliedStartPos = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t];
				if (k == readPaths[i].paths[j].path.size()-1)
				{
					lastRightClip = readPaths[i].paths[j].pathRightClipKmers;
				}
				else
				{
					lastRightClip = 0;
				}
				continue;
			}
		}
	}
	SparseEdgeContainer topologicalEdges = getTopologicalEdges(unitigGraph, anchor);
	SparseEdgeContainer result;
	result.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (!anchor[i]) continue;
		std::pair<size_t, bool> fw { i, true };
		if (topologicalEdges.getEdges(fw).size() == 1)
		{
			auto other = reverse(topologicalEdges.getEdges(fw)[0]);
			assert(anchor[other.first]);
			if (topologicalEdges.getEdges(other).size() == 1)
			{
				assert(reverse(topologicalEdges.getEdges(other)[0]) == fw);
				auto key = canon(fw, topologicalEdges.getEdges(fw)[0]);
				if (edgeCoverage.hasValue(key.first, key.second) && edgeCoverage.get(key.first, key.second) >= minCoverage)
				{
					result.addEdge(fw, topologicalEdges.getEdges(fw)[0]);
				}
			}
		}
		std::pair<size_t, bool> bw { i, false };
		if (topologicalEdges.getEdges(bw).size() == 1)
		{
			auto other = reverse(topologicalEdges.getEdges(bw)[0]);
			assert(anchor[other.first]);
			if (topologicalEdges.getEdges(other).size() == 1)
			{
				assert(reverse(topologicalEdges.getEdges(other)[0]) == bw);
				auto key = canon(bw, topologicalEdges.getEdges(bw)[0]);
				if (edgeCoverage.hasValue(key.first, key.second) && edgeCoverage.get(key.first, key.second) >= minCoverage)
				{
					result.addEdge(bw, topologicalEdges.getEdges(bw)[0]);
				}
			}
		}
	}
	return result;
}

std::vector<uint64_t> getReachableAnchors(const SparseEdgeContainer& edges, const uint64_t start, const std::vector<bool>& anchor)
{
	assert(!anchor[start & maskUint64_t]);
	std::vector<uint64_t> stack;
	stack.emplace_back(start);
	phmap::flat_hash_set<uint64_t> visited;
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		visited.insert(top);
		if (anchor[top & maskUint64_t]) continue;
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
		}
	}
	std::vector<uint64_t> result;
	for (uint64_t node : visited)
	{
		if (!anchor[node & maskUint64_t]) continue;
		result.emplace_back(node);
	}
	return result;
}

void setReachableAnchors(const SparseEdgeContainer& edges, std::vector<std::vector<uint64_t>>& forwardReachableAnchors, std::vector<std::vector<uint64_t>>& backwardReachableAnchors, const uint64_t startNode, const std::vector<bool>& anchor)
{
	std::vector<uint64_t> stack;
	phmap::flat_hash_set<uint64_t> checked;
	for (auto edge : edges.getEdges(std::make_pair(startNode & maskUint64_t, startNode & firstBitUint64_t)))
	{
		stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
	}
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (checked.count(top) == 1) continue;
		if (anchor[top & maskUint64_t]) continue;
		checked.insert(top);
		if (top & firstBitUint64_t)
		{
			backwardReachableAnchors[top & maskUint64_t].emplace_back(startNode);
		}
		else
		{
			forwardReachableAnchors[top & maskUint64_t].emplace_back(startNode);
		}
		for (const auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
		}
	}
}

void extendAnchors(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, std::vector<bool>& anchor)
{
	size_t maxTangleCheckSize = 100;
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<bool> tryExtendingThese;
	tryExtendingThese.resize(unitigGraph.nodeCount(), true);
	{
		std::vector<size_t> tipParent;
		tipParent.resize(unitigGraph.nodeCount()*2);
		for (size_t i = 0; i < unitigGraph.nodeCount()*2; i++)
		{
			tipParent[i] = i;
		}
		for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
		{
			if (!anchor[i])
			{
				merge(tipParent, i*2, i*2+1);
			}
			for (auto edge : edges.getEdges(std::make_pair(i, true)))
			{
				merge(tipParent, i*2+1, edge.first*2 + (edge.second ? 0 : 1));
			}
			for (auto edge : edges.getEdges(std::make_pair(i, false)))
			{
				merge(tipParent, i*2, edge.first*2 + (edge.second ? 0 : 1));
			}
		}
		phmap::flat_hash_map<size_t, size_t> clusterCoverage;
		for (size_t i = 0; i < unitigGraph.nodeCount()*2; i++)
		{
			clusterCoverage[find(tipParent, i)] += 1;
		}
		for (size_t i = 0; i < tryExtendingThese.size(); i++)
		{
			if (anchor[i]) continue;
			assert(find(tipParent, i*2) == find(tipParent, i*2+1));
			if (clusterCoverage[find(tipParent, i*2)] < maxTangleCheckSize) continue;
			tryExtendingThese[i] = false;
		}
	}
	std::vector<bool> anchorsBeforeExtension = anchor;
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (anchor[i]) continue;
		if (!tryExtendingThese[i]) continue;
		if (!cutsGraph(unitigGraph, edges, i, anchorsBeforeExtension)) continue;
		auto fwAnchors = getReachableAnchors(edges, i+firstBitUint64_t, anchorsBeforeExtension);
		if (fwAnchors.size() == 0) continue;
		auto bwAnchors = getReachableAnchors(edges, i, anchorsBeforeExtension);
		if (bwAnchors.size() == 0) continue;
		if (fwAnchors.size() != 1 && bwAnchors.size() != 1) continue;
		anchor[i] = true;
	}
}

size_t getPloidy(const std::vector<uint64_t>& nodes, const UnitigGraph& unitigGraph, const double approxOneHapCoverage)
{
	assert(nodes.size() >= 1);
	double coverageSum = 0;
	size_t coverageDivisor = 0;
	for (auto node : nodes)
	{
		assert((node & maskUint64_t) < unitigGraph.nodeCount());
		assert(unitigGraph.lengths[node & maskUint64_t] >= 1);
		coverageDivisor += unitigGraph.lengths[node & maskUint64_t];
		coverageSum += unitigGraph.lengths[node & maskUint64_t] * unitigGraph.coverages[node & maskUint64_t];
	}
	assert(coverageDivisor > 0);
	size_t estimatedPloidy = (coverageSum / coverageDivisor) / approxOneHapCoverage + 0.5;
	if (coverageDivisor > 1000 && estimatedPloidy == 0 && (coverageSum / coverageDivisor) >= approxOneHapCoverage*0.25) estimatedPloidy = 1;
	if (coverageDivisor > 10000 && estimatedPloidy == 0 && (coverageSum / coverageDivisor) >= 3) estimatedPloidy = 1;
	return estimatedPloidy;
}

std::vector<uint64_t> getChainOneWay(const uint64_t start, const SparseEdgeContainer& edges)
{
	uint64_t pos = start;
	std::vector<uint64_t> result;
	result.push_back(pos);
	while (edges.getEdges(std::make_pair(pos & maskUint64_t, pos & firstBitUint64_t)).size() == 1)
	{
		std::pair<size_t, bool> nextp = edges.getEdges(std::make_pair(pos & maskUint64_t, pos & firstBitUint64_t))[0];
		if (edges.getEdges(reverse(nextp)).size() != 1) break;
		uint64_t next = nextp.first + (nextp.second ? firstBitUint64_t : 0);
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next == start) break;
		pos = next;
		result.push_back(pos);
	}
	return result;
}

std::vector<uint64_t> getChain(std::vector<bool>& checked, const size_t startNode, const SparseEdgeContainer& edges)
{
	assert(edges.size() == checked.size());
	assert(startNode < checked.size());
	assert(!checked[startNode]);
	std::vector<uint64_t> fwChain = getChainOneWay(startNode + firstBitUint64_t, edges);
	std::vector<uint64_t> bwChain = getChainOneWay(startNode, edges);
	assert(fwChain.size() >= 1);
	assert(bwChain.size() >= 1);
	assert(fwChain[0] == startNode + firstBitUint64_t);
	assert(bwChain[0] == startNode);
	if (fwChain.size() == bwChain.size())
	{
		if (fwChain.size() >= 2)
		{
			if (fwChain.back() == (bwChain[1] ^ firstBitUint64_t))
			{
				for (size_t i = 0; i < fwChain.size(); i++)
				{
					assert(i == 0 || (fwChain[i] == (bwChain[bwChain.size()-i] ^ firstBitUint64_t)));
					assert(!checked[fwChain[i] & maskUint64_t]);
					checked[fwChain[i] & maskUint64_t] = true;
				}
				return fwChain;
			}
		}
	}
	std::reverse(bwChain.begin(), bwChain.end());
	for (size_t i = 0; i < bwChain.size(); i++)
	{
		bwChain[i] ^= firstBitUint64_t;
	}
	assert(bwChain.back() == fwChain[0]);
	bwChain.insert(bwChain.end(), fwChain.begin()+1, fwChain.end());
	for (size_t i = 0; i < bwChain.size(); i++)
	{
		assert((bwChain[i] & maskUint64_t) < checked.size());
		assert(!checked[bwChain[i] & maskUint64_t]);
		checked[bwChain[i] & maskUint64_t] = true;
	}
	return bwChain;
}

void getDistances(phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<int>>& distances, const std::vector<bool>& anchor, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths)
{
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		uint64_t lastAnchor = std::numeric_limits<uint64_t>::max();
		size_t lastAnchorOffset = 0;
		size_t lastAnchorPlus = 0;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readOffset = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				if (!anchor[node & maskUint64_t])
				{
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				if (lastAnchor == std::numeric_limits<uint64_t>::max())
				{
					lastAnchor = node;
					lastAnchorOffset = readOffset;
					lastAnchorPlus = 0;
					if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				if (distances.count(std::make_pair(lastAnchor, node)) == 1)
				{
					int distance = (int)readOffset - (int)lastAnchorOffset + (int)lastAnchorPlus;
					// assert(distance >= 0);
					if (k == 0) distance -= (int)readPaths[i].paths[j].pathLeftClipKmers;
					if (distance >= 0) distances.at(std::make_pair(lastAnchor, node)).emplace_back(distance);
					lastAnchor = node;
					lastAnchorOffset = readOffset;
					lastAnchorPlus = 0;
					if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				if (distances.count(std::make_pair(node ^ firstBitUint64_t, lastAnchor ^ firstBitUint64_t)) == 1)
				{
					int distance = (int)readOffset - (int)lastAnchorOffset + (int)lastAnchorPlus;
					// assert(distance >= 0);
					if (k == 0) distance -= (int)readPaths[i].paths[j].pathLeftClipKmers;
					distance += (int)unitigGraph.lengths[node & maskUint64_t] - (int)unitigGraph.lengths[lastAnchor & maskUint64_t];
					if (distance >= 0) distances.at(std::make_pair(node ^ firstBitUint64_t, lastAnchor ^ firstBitUint64_t)).emplace_back(distance);
					lastAnchor = node;
					lastAnchorOffset = readOffset;
					lastAnchorPlus = 0;
					if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				lastAnchor = node;
				lastAnchorOffset = readOffset;
				lastAnchorPlus = 0;
				if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
				readOffset += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
			}
		}
	}
}

double average(const std::vector<int>& values)
{
	assert(values.size() >= 1);
	double sum = 0;
	for (auto value : values)
	{
		sum += value;
	}
	return sum / values.size();
}

void addOffsets(std::vector<AnchorChain>& chains, const std::vector<bool>& anchor, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<int>> distances;
	for (size_t i = 0; i < chains.size(); i++)
	{
		for (size_t j = 1; j < chains[i].nodes.size(); j++)
		{
			assert(distances.count(std::make_pair(chains[i].nodes[j-1], chains[i].nodes[j])) == 0);
			distances[std::make_pair(chains[i].nodes[j-1], chains[i].nodes[j])];
		}
	}
	getDistances(distances, anchor, unitigGraph, readPaths);
	for (size_t i = 0; i < chains.size(); i++)
	{
		chains[i].nodeOffsets.push_back(0);
		for (size_t j = 1; j < chains[i].nodes.size(); j++)
		{
			auto key = std::make_pair(chains[i].nodes[j-1], chains[i].nodes[j]);
			assert(distances.at(key).size() >= 1);
			size_t distance = average(distances.at(key));
			// std::cerr << "chain " << i << " anchor " << j << " (" << (chains[i].nodes[j] & maskUint64_t) << ") distance " << distance << std::endl;
			assert(distance >= 0);
			assert(distance < 100000000);
			chains[i].nodeOffsets.push_back(chains[i].nodeOffsets.back() + distance);
		}
		assert(chains[i].nodeOffsets.size() == chains[i].nodes.size());
	}
}

std::vector<AnchorChain> getAnchorChains(const std::vector<bool>& anchor, const SparseEdgeContainer& edges, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	std::vector<AnchorChain> result;
	std::vector<bool> checked;
	checked.resize(anchor.size(), false);
	for (size_t i = 0; i < anchor.size(); i++)
	{
		if (!anchor[i]) continue;
		if (checked[i]) continue;
		result.emplace_back();
		result.back().nodes = getChain(checked, i, edges);
		result.back().ploidy = getPloidy(result.back().nodes, unitigGraph, approxOneHapCoverage);
	}
	addOffsets(result, anchor, unitigGraph, readPaths);
	return result;
}

SparseEdgeContainer getAnchorChainEdges(const std::vector<AnchorChain>& anchorChains, const SparseEdgeContainer& anchorEdges)
{
	SparseEdgeContainer result;
	result.resize(anchorChains.size());
	phmap::flat_hash_map<uint64_t, uint64_t> anchorStart;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		assert(anchorStart.count(anchorChains[i].nodes[0]) == 0);
		anchorStart[anchorChains[i].nodes[0]] = i + firstBitUint64_t;
		assert(anchorStart.count(anchorChains[i].nodes.back() ^ firstBitUint64_t) == 0);
		anchorStart[anchorChains[i].nodes.back() ^ firstBitUint64_t] = i;
	}
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::pair<size_t, bool> fw { anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t };
		for (auto edge : anchorEdges.getEdges(fw))
		{
			assert(anchorStart.count(edge.first + (edge.second ? firstBitUint64_t : 0)));
			uint64_t target = anchorStart.at(edge.first + (edge.second ? firstBitUint64_t : 0));
			result.addEdge(std::make_pair(i, true), std::make_pair(target & maskUint64_t, target & firstBitUint64_t));
		}
		std::pair<size_t, bool> bw { anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] & firstBitUint64_t) ^ firstBitUint64_t };
		for (auto edge : anchorEdges.getEdges(bw))
		{
			assert(anchorStart.count(edge.first + (edge.second ? firstBitUint64_t : 0)));
			uint64_t target = anchorStart.at(edge.first + (edge.second ? firstBitUint64_t : 0));
			result.addEdge(std::make_pair(i, false), std::make_pair(target & maskUint64_t, target & firstBitUint64_t));
		}
	}
	return result;
}

uint64_t getBubble(const uint64_t startPos, const SparseEdgeContainer& edges)
{
	std::vector<std::pair<size_t, bool>> S;
	phmap::flat_hash_set<std::pair<size_t, bool>> seen;
	phmap::flat_hash_set<std::pair<size_t, bool>> visited;
	S.emplace_back(startPos & maskUint64_t, startPos & firstBitUint64_t);
	seen.emplace(startPos & maskUint64_t, startPos & firstBitUint64_t);
	while (S.size() > 0)
	{
		auto v = S.back();
		S.pop_back();
		assert(seen.count(v) == 1);
		seen.erase(v);
		assert(visited.count(v) == 0);
		visited.insert(v);
		if (edges.getEdges(v).size() == 0) return std::numeric_limits<size_t>::max();
		for (auto u : edges.getEdges(v))
		{
			if (u.first == v.first) return std::numeric_limits<size_t>::max();
			if (u.first + (u.second ? firstBitUint64_t : 0) == startPos) return std::numeric_limits<size_t>::max();
			assert(visited.count(u) == 0);
			seen.insert(u);
			bool hasNonvisitedParent = false;
			assert(edges.getEdges(reverse(u)).size() >= 1);
			for (auto edge : edges.getEdges(reverse(u)))
			{
				if (visited.count(reverse(edge)) == 0) hasNonvisitedParent = true;
			}
			if (!hasNonvisitedParent) S.push_back(u);
		}
		if (S.size() == 1 && seen.size() == 1 && S[0] == *seen.begin())
		{
			auto t = S[0];
			for (auto edge : edges.getEdges(t))
			{
				if (edge.first + (edge.second ? firstBitUint64_t : 0) == startPos) return std::numeric_limits<size_t>::max();
			}
			return t.first + (t.second ? firstBitUint64_t : 0);
		}
	}
	return std::numeric_limits<size_t>::max();
}

std::vector<uint64_t> getChainOfChains(const uint64_t start, const SparseEdgeContainer& anchorEdges)
{
	uint64_t pos = start;
	while (true)
	{
		uint64_t next = getBubble(pos, anchorEdges);
		if (next == std::numeric_limits<size_t>::max()) break;
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next == start) break;
		pos = next;
	}
	pos = pos ^ firstBitUint64_t;
	const uint64_t resultStart = pos;
	std::vector<uint64_t> result;
	result.emplace_back(pos);
	while (true)
	{
		uint64_t next = getBubble(pos, anchorEdges);
		if (next == std::numeric_limits<size_t>::max()) break;
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next == resultStart) break;
		pos = next;
		result.emplace_back(pos);
	}
	return result;
}

std::vector<uint64_t> getSimpleBubble(const uint64_t bubbleStart, const uint64_t bubbleEnd, const SparseEdgeContainer& edges)
{
	std::vector<uint64_t> result;
	for (auto edge : edges.getEdges(std::make_pair(bubbleStart & maskUint64_t, bubbleStart & firstBitUint64_t)))
	{
		if (edges.getEdges(edge).size() != 1) return std::vector<uint64_t>{};
		if (edges.getEdges(edge)[0] != std::make_pair<size_t, bool>(bubbleEnd & maskUint64_t, bubbleEnd & firstBitUint64_t)) return std::vector<uint64_t>{};
		if (edges.getEdges(reverse(edge)).size() != 1) return std::vector<uint64_t>{};
		result.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
	}
	return result;
}

void refineAnchorChainPloidies(std::vector<AnchorChain>& anchorChains, const SparseEdgeContainer& anchorEdges, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	SparseEdgeContainer anchorChainEdges = getAnchorChainEdges(anchorChains, anchorEdges);
	std::vector<bool> checked;
	std::vector<bool> refined;
	refined.resize(anchorChains.size(), false);
	checked.resize(anchorChains.size(), false);
	std::vector<std::vector<uint64_t>> chainsOfChains;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (checked[i]) continue;
		auto chain = getChainOfChains(i, anchorChainEdges);
		for (auto node : chain)
		{
			assert((node & maskUint64_t) < checked.size());
			assert(!checked[node & maskUint64_t]);
			checked[node & maskUint64_t] = true;
		}
		chainsOfChains.emplace_back(chain);
	}
	std::sort(chainsOfChains.begin(), chainsOfChains.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
	for (size_t i = 0; i < chainsOfChains.size(); i++)
	{
		if (chainsOfChains[i].size() == 1 && refined[chainsOfChains[i][0] & maskUint64_t]) continue;
		double coverageSum = 0;
		size_t coverageDivisor = 0;
		for (const uint64_t chain : chainsOfChains[i])
		{
			for (const uint64_t node : anchorChains[chain & maskUint64_t].nodes)
			{
				coverageSum += unitigGraph.coverages[node & maskUint64_t] * unitigGraph.lengths[node & maskUint64_t];
				coverageDivisor += unitigGraph.lengths[node & maskUint64_t];
			}
		}
		size_t refinedPloidy = coverageSum / (double)coverageDivisor / approxOneHapCoverage + 0.5;
		if (refinedPloidy == 0 && coverageDivisor > 10000 && coverageSum / (double)coverageDivisor >= 3) refinedPloidy = 1;
		for (const uint64_t chain : chainsOfChains[i])
		{
			// if (anchorChains[chain & maskUint64_t].ploidy != refinedPloidy)
			// {
			// 	std::cerr << "refined chain of chains core chain " << (anchorChains[chain & maskUint64_t].nodes[0] & maskUint64_t) << "-" << (anchorChains[chain & maskUint64_t].nodes.back() & maskUint64_t) << " ploidy from " << anchorChains[chain & maskUint64_t].ploidy << " to " << refinedPloidy << std::endl;
			// }
			anchorChains[chain & maskUint64_t].ploidy = refinedPloidy;
			refined[chain & maskUint64_t] = true;
		}
		for (size_t j = 0; j+1 < chainsOfChains[i].size(); j++)
		{
			auto bubbleChains = getSimpleBubble(chainsOfChains[i][j], chainsOfChains[i][j+1], anchorChainEdges);
			if (bubbleChains.size() == 0)
			{
				// std::cerr << "could not make a refined coverage estimate for chain of anchor chains bubble between chains " << (anchorChains[chainsOfChains[i][j-1]].nodes[0] & maskUint64_t) << "-" << (anchorChains[chainsOfChains[i][j-1]].nodes.back() & maskUint64_t) << " and " << (anchorChains[chainsOfChains[i][j]].nodes[0] & maskUint64_t) << "-" << (anchorChains[chainsOfChains[i][j]].nodes.back() & maskUint64_t) << std::endl;
			}
			else
			{
				if (bubbleChains.size() == refinedPloidy)
				{
					for (const uint64_t chain : bubbleChains)
					{
						// if (anchorChains[chain & maskUint64_t].ploidy != 1)
						// {
						// 	std::cerr << "refined chain of chains bubble chain " << (anchorChains[chain & maskUint64_t].nodes[0] & maskUint64_t) << "-" << (anchorChains[chain & maskUint64_t].nodes.back() & maskUint64_t) << " ploidy from " << anchorChains[chain & maskUint64_t].ploidy << " to " << 1 << std::endl;
						// }
						anchorChains[chain & maskUint64_t].ploidy = 1;
						refined[chain & maskUint64_t] = true;
					}
				}
				else if (bubbleChains.size() > refinedPloidy)
				{
					// std::cerr << "could not distribute ploidy between too many chains (ploidy " << refinedPloidy << " vs chain count" << bubbleChains.size() << ")";
					// for (const uint64_t chain : bubbleChains)
					// {
					// 	std::cerr << " " << (anchorChains[chain & maskUint64_t].nodes[0] & maskUint64_t) << "-" << (anchorChains[chain & maskUint64_t].nodes.back() & maskUint64_t);
					// }
					// std::cerr << std::endl;
				}
				else
				{
					assert(bubbleChains.size() < refinedPloidy);
					std::cerr << "todo implement: could not distribute ploidy between too few chains (ploidy " << refinedPloidy << " vs chain count" << bubbleChains.size() << ")" << std::endl;
					// for (const uint64_t chain : bubbleChains)
					// {
					// 	std::cerr << " " << (anchorChains[chain & maskUint64_t].nodes[0] & maskUint64_t) << "-" << (anchorChains[chain & maskUint64_t].nodes.back() & maskUint64_t);
					// }
					// std::cerr << std::endl;
				}
			}
		}
	}
}

std::vector<AnchorChain> getAnchorChains(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t initialAnchorMinLength, const double approxOneHapCoverage)
{
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.lengths[i] >= initialAnchorMinLength) anchor[i] = true;
	}
	extendAnchors(unitigGraph, readPaths, anchor);
	SparseEdgeContainer edges = getAnchorEdges(unitigGraph, readPaths, anchor, approxOneHapCoverage*.25);
	for (size_t i = 0; i < edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : edges.getEdges(fw))
		{
			assert(edges.hasEdge(reverse(pair), reverse(fw)));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : edges.getEdges(bw))
		{
			assert(edges.hasEdge(reverse(pair), reverse(bw)));
		}
	}
	auto anchors = getAnchorChains(anchor, edges, unitigGraph, readPaths, approxOneHapCoverage);
	// for (size_t i = 0; i < anchors.size(); i++)
	// {
	// 	std::cerr << "raw chain " << i;
	// 	for (uint64_t node : anchors[i].nodes)
	// 	{
	// 		std::cerr << " " << (node & maskUint64_t);
	// 	}
	// 	std::cerr << std::endl;
	// }
	anchor.assign(unitigGraph.nodeCount(), false);
	for (const auto& chain : anchors)
	{
		if (chain.ploidy == 0) continue;
		for (auto node : chain.nodes)
		{
			anchor[node & maskUint64_t] = true;
		}
	}
	edges = getAnchorEdges(unitigGraph, readPaths, anchor, approxOneHapCoverage*.25);
	auto result = getAnchorChains(anchor, edges, unitigGraph, readPaths, approxOneHapCoverage);
	refineAnchorChainPloidies(result, edges, unitigGraph, readPaths, approxOneHapCoverage);
	return result;
}

std::vector<std::vector<ChainPosition>> getReadChainPositions(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains)
{
	std::vector<std::pair<size_t, size_t>> anchorLocator;
	anchorLocator.resize(unitigGraph.nodeCount(), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			anchorLocator[anchorChains[i].nodes[j] & maskUint64_t] = std::make_pair(i, j);
		}
	}
	std::vector<std::vector<ChainPosition>> result;
	result.resize(readPaths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, bool, int, size_t, int, int, size_t>> impliedPositions; // chain, fw, diagonal, readPos, readChainStartPos, readChainEndPos, kmerMatches
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				assert(readPos < readPaths[i].readLength);
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (anchorLocator[node & maskUint64_t].first == std::numeric_limits<size_t>::max()) continue;
				const size_t chain = anchorLocator[node & maskUint64_t].first;
				const size_t offset = anchorLocator[node & maskUint64_t].second;
				assert((node & maskUint64_t) == (anchorChains[chain].nodes[offset] & maskUint64_t));
				const bool fw = (node & firstBitUint64_t) == (anchorChains[chain].nodes[offset] & firstBitUint64_t);
				// std::cerr << "anchorchain read " << i << " matches chain " << chain << " node " << (node & maskUint64_t) << std::endl;
				size_t kmerMatches = unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) kmerMatches -= readPaths[i].paths[j].pathLeftClipKmers;
				if (k == readPaths[i].paths[j].path.size()-1) kmerMatches -= readPaths[i].paths[j].pathRightClipKmers;
				int diagonal = (int)readPos - (int)anchorChains[chain].nodeOffsets[offset] - (int)unitigGraph.lengths[anchorChains[chain].nodes[offset] & maskUint64_t];
				if (!fw) diagonal = (int)readPos + (int)anchorChains[chain].nodeOffsets[offset];
				int startPosInRead;
				int endPosInRead;
				if (fw)
				{
					// readPos is at the end of current node
					startPosInRead = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t] - (int)anchorChains[chain].nodeOffsets[offset];
					endPosInRead = (int)readPos + (int)anchorChains[chain].nodeOffsets.back() - (int)anchorChains[chain].nodeOffsets[offset] - (int)unitigGraph.lengths[node & maskUint64_t] + (int)unitigGraph.lengths[anchorChains[chain].nodes.back() & maskUint64_t];
					assert(startPosInRead < (int)readPaths[i].readLength);
				}
				else
				{
					// readPos is at the start of current node
					startPosInRead = (int)readPos - ((int)((int)unitigGraph.lengths[anchorChains[chain].nodes.back() & maskUint64_t] + (int)anchorChains[chain].nodeOffsets.back() - (int)anchorChains[chain].nodeOffsets[offset]));
					endPosInRead = (int)readPos + (int)anchorChains[chain].nodeOffsets[offset];
					assert(startPosInRead < (int)readPaths[i].readLength);
				}
				assert(startPosInRead < (int)readPaths[i].readLength);
				assert(endPosInRead > 0);
				impliedPositions.emplace_back(chain, fw, diagonal, readPos, startPosInRead, endPosInRead, kmerMatches);
			}
		}
		if (impliedPositions.size() == 0) continue;
		std::sort(impliedPositions.begin(), impliedPositions.end());
		size_t currentMinReadPos = std::numeric_limits<size_t>::max();
		size_t currentMaxReadPos = 0;
		for (size_t j = 0; j < impliedPositions.size(); j++)
		{
			if (j == 0 || std::get<0>(impliedPositions[j]) != std::get<0>(impliedPositions[j-1]) || std::get<1>(impliedPositions[j]) != std::get<1>(impliedPositions[j-1]) || std::get<2>(impliedPositions[j]) > std::get<2>(impliedPositions[j-1]) + (anchorChains[std::get<0>(impliedPositions[j])].nodeOffsets.back() + unitigGraph.lengths[anchorChains[std::get<0>(impliedPositions[j])].nodes.back() & maskUint64_t])*0.5)
			{
				result[i].emplace_back();
				result[i].back().chain = std::get<0>(impliedPositions[j]) + (std::get<1>(impliedPositions[j]) ? firstBitUint64_t : 0);
				currentMinReadPos = std::numeric_limits<size_t>::max();
				currentMaxReadPos = 0;
			}
			result[i].back().numKmerMatches += std::get<6>(impliedPositions[j]);
			if (std::get<3>(impliedPositions[j]) < currentMinReadPos)
			{
				result[i].back().chainStartPosInRead = std::get<4>(impliedPositions[j]);
				currentMinReadPos = std::get<3>(impliedPositions[j]);
			}
			if (std::get<3>(impliedPositions[j]) > currentMaxReadPos)
			{
				result[i].back().chainEndPosInRead = std::get<5>(impliedPositions[j]);
				currentMaxReadPos = std::get<3>(impliedPositions[j]);
			}
		}
		std::vector<bool> keep;
		keep.resize(result[i].size(), true);
		for (size_t j = 0; j < result[i].size(); j++)
		{
			for (size_t k = 0; k < result[i].size(); k++)
			{
				if (j == k) continue;
				if (result[i][j].chainStartPosInRead >= result[i][k].chainEndPosInRead) continue;
				if (result[i][j].chainEndPosInRead <= result[i][k].chainStartPosInRead) continue;
				assert(result[i][j].chainEndPosInRead > result[i][k].chainStartPosInRead);
				int overlap = std::min(result[i][j].chainEndPosInRead, result[i][k].chainEndPosInRead) - std::max(result[i][j].chainStartPosInRead, result[i][k].chainStartPosInRead);
				assert(overlap >= 1);
				if (overlap < (result[i][j].chainEndPosInRead - result[i][j].chainStartPosInRead) * 0.1) continue;
				if (overlap < (result[i][k].chainEndPosInRead - result[i][k].chainStartPosInRead) * 0.1) continue;
				if (result[i][j].numKmerMatches > result[i][k].numKmerMatches) keep[k] = false;
				if (result[i][j].numKmerMatches < result[i][k].numKmerMatches) keep[j] = false;
			}
		}
		for (size_t j = result[i].size()-1; j < result[i].size(); j--)
		{
			if (keep[j]) continue;
			std::swap(result[i][j], result[i].back());
			result[i].pop_back();
		}
		std::sort(result[i].begin(), result[i].end(), [](auto left, auto right) { return left.chainStartPosInRead < right.chainStartPosInRead; });
	}
	return result;
}
