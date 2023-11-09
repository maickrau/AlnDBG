#include <iostream>
#include "MBGCommon.h"
#include "GraphPhaser.h"
#include "Common.h"

std::pair<size_t, size_t> findSimpleBubble(const SparseEdgeContainer& activeEdges, const std::pair<size_t, bool> start, const UnitigGraph& graph, const double approxOneHapCoverage)
{
	if (activeEdges.getEdges(start).size() != 2) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.coverages[start.first] < approxOneHapCoverage * 1.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.coverages[start.first] > approxOneHapCoverage * 2.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	std::pair<size_t, bool> otherSideNode { std::numeric_limits<size_t>::max(), true };
	size_t otherSideNodeFound = 0;
	for (auto edge : activeEdges.getEdges(start))
	{
		if (graph.coverages[edge.first] < approxOneHapCoverage * 0.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (graph.coverages[edge.first] > approxOneHapCoverage * 1.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (activeEdges.getEdges(edge).size() > 1) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (activeEdges.getEdges(reverse(edge)).size() != 1) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (activeEdges.getEdges(edge).size() == 1)
		{
			otherSideNodeFound += 1;
			if (otherSideNode.first == std::numeric_limits<size_t>::max())
			{
				otherSideNode = activeEdges.getEdges(edge)[0];
			}
			else
			{
				if (activeEdges.getEdges(edge)[0] != otherSideNode) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
			}
		}
	}
	if (otherSideNodeFound == 1) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	auto edges = activeEdges.getEdges(start);
	if ((graph.lengths[edges[0].first] < 100 || graph.lengths[edges[1].first] < 100) && (otherSideNodeFound != 2)) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (otherSideNodeFound == 2)
	{
		if (graph.coverages[otherSideNode.first] < approxOneHapCoverage * 1.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (graph.coverages[otherSideNode.first] > approxOneHapCoverage * 2.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	}
	return std::make_pair(edges[0].first, edges[1].first);
}

std::vector<std::pair<size_t, size_t>> getSimpleBubbles(const SparseEdgeContainer& activeEdges, const UnitigGraph& graph, const double approxOneHapCoverage)
{
	phmap::flat_hash_set<std::pair<size_t, size_t>> result;
	for (size_t i = 0; i < activeEdges.size(); i++)
	{
		auto fwBubble = findSimpleBubble(activeEdges, std::make_pair(i, true), graph, approxOneHapCoverage);
		if (fwBubble.first < activeEdges.size())
		{
			assert(fwBubble.second < activeEdges.size());
			assert(fwBubble.first != fwBubble.second);
			result.emplace(std::min(fwBubble.first, fwBubble.second), std::max(fwBubble.first, fwBubble.second));
		}
		auto bwBubble = findSimpleBubble(activeEdges, std::make_pair(i, true), graph, approxOneHapCoverage);
		if (bwBubble.first < activeEdges.size())
		{
			assert(bwBubble.second < activeEdges.size());
			assert(bwBubble.first != bwBubble.second);
			result.emplace(std::min(bwBubble.first, bwBubble.second), std::max(bwBubble.first, bwBubble.second));
		}
	}
	return std::vector<std::pair<size_t, size_t>> { result.begin(), result.end() };
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> getReadsPerBubbleAllele(const std::vector<ReadPathBundle>& readPaths, const std::vector<std::pair<size_t, size_t>>& simpleBubbles)
{
	phmap::flat_hash_map<size_t, size_t> nodeToBubbleAllele;
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> result;
	result.resize(simpleBubbles.size());
	for (size_t i = 0; i < simpleBubbles.size(); i++)
	{
		nodeToBubbleAllele[simpleBubbles[i].first] = i;
		nodeToBubbleAllele[simpleBubbles[i].second] = i + firstBitUint64_t;
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				if (nodeToBubbleAllele.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 0) continue;
				size_t node = nodeToBubbleAllele.at(readPaths[i].paths[j].path[k] & maskUint64_t);
				if (node & firstBitUint64_t)
				{
					result[node & maskUint64_t].second.emplace_back(i);
				}
				else
				{
					result[node & maskUint64_t].first.emplace_back(i);
				}
			}
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		phmap::flat_hash_set<size_t> uniqs { result[i].first.begin(), result[i].first.end() };
		result[i].first.clear();
		result[i].first.insert(result[i].first.end(), uniqs.begin(), uniqs.end());
		std::sort(result[i].first.begin(), result[i].first.end());
		uniqs.clear();
		uniqs.insert(result[i].second.begin(), result[i].second.end());
		result[i].second.clear();
		result[i].second.insert(result[i].second.end(), uniqs.begin(), uniqs.end());
		std::sort(result[i].second.begin(), result[i].second.end());
	}
	return result;
}

size_t find(std::vector<size_t>& parent, const size_t i)
{
	while (parent[parent[i]] != parent[i]) parent[i] = parent[parent[i]];
	return parent[i];
}

void merge(std::vector<size_t>& parent, size_t i, size_t j)
{
	i = find(parent, i);
	j = find(parent, j);
	parent[j] = i;
}

bool bubblesMatch(const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& readsPerBubbleAllele, const size_t left, const size_t right)
{
	size_t normalOne = intersectSize(readsPerBubbleAllele[left].first, readsPerBubbleAllele[right].first);
	size_t normalTwo = intersectSize(readsPerBubbleAllele[left].second, readsPerBubbleAllele[right].second);
	size_t crossOne = intersectSize(readsPerBubbleAllele[left].first, readsPerBubbleAllele[right].second);
	size_t crossTwo = intersectSize(readsPerBubbleAllele[left].second, readsPerBubbleAllele[right].first);
	if (normalOne+normalTwo < crossOne+crossTwo)
	{
		std::swap(normalTwo, crossTwo);
		std::swap(normalOne, crossOne);
	}
	if (normalOne+normalTwo < 3) return false;
	if (normalOne == 0) return false;
	if (normalTwo == 0) return false;
	if (crossOne+crossTwo > 3) return false;
	if (normalOne+normalTwo < (crossOne+crossTwo)*3) return false;
	return true;
}

std::vector<std::vector<size_t>> mergeBubblesToClusters(const std::vector<std::pair<size_t, size_t>>& simpleBubbles, const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& readsPerBubbleAllele)
{
	std::vector<size_t> bubbleParent;
	for (size_t i = 0; i < simpleBubbles.size(); i++)
	{
		bubbleParent.emplace_back(i);
	}
	for (size_t i = 0; i < simpleBubbles.size(); i++)
	{
		for (size_t j = 1; j < simpleBubbles.size(); j++)
		{
			if (!bubblesMatch(readsPerBubbleAllele, i, j)) continue;
			merge(bubbleParent, i, j);
		}
	}
	RankBitvector clusterToIndex;
	clusterToIndex.resize(simpleBubbles.size());
	for (size_t i = 0; i < simpleBubbles.size(); i++)
	{
		if (find(bubbleParent, i) != i) continue;
		clusterToIndex.set(i, true);
	}
	clusterToIndex.buildRanks();
	std::vector<std::vector<size_t>> bubbleIndicesPerCluster;
	bubbleIndicesPerCluster.resize(clusterToIndex.getRank(clusterToIndex.size()-1) + (clusterToIndex.get(clusterToIndex.size()-1) ? 1 : 0));
	for (size_t i = 0; i < simpleBubbles.size(); i++)
	{
		bubbleIndicesPerCluster[clusterToIndex.getRank(find(bubbleParent, i))].push_back(i);
	}
	for (size_t i = bubbleIndicesPerCluster.size()-1; i < bubbleIndicesPerCluster.size(); i--)
	{
		if (bubbleIndicesPerCluster[i].size() != 1) continue;
		std::swap(bubbleIndicesPerCluster[i], bubbleIndicesPerCluster.back());
		bubbleIndicesPerCluster.pop_back();
	}
	return bubbleIndicesPerCluster;
}

// positive: good in same phase, negative: good in cross phase, zero: ambiguous
int getPhaseScore(const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& readsPerBubbleAllele, const size_t left, const size_t right)
{
	int result = 0;
	result += intersectSize(readsPerBubbleAllele[left].first, readsPerBubbleAllele[right].first) + intersectSize(readsPerBubbleAllele[left].second, readsPerBubbleAllele[right].second);
	result -= intersectSize(readsPerBubbleAllele[left].second, readsPerBubbleAllele[right].first) + intersectSize(readsPerBubbleAllele[left].first, readsPerBubbleAllele[right].second);
	return result;
}

std::vector<bool> getBubbleOrientations(const std::vector<size_t>& bubbleIndicesPerCluster, const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& readsPerBubbleAllele)
{
	std::vector<bool> result;
	result.resize(bubbleIndicesPerCluster.size(), true);
	std::vector<bool> handled;
	handled.resize(bubbleIndicesPerCluster.size(), false);
	size_t maxNode = 0;
	size_t maxTotal = 0;
	for (size_t i = 0; i < bubbleIndicesPerCluster.size(); i++)
	{
		size_t total = 0;
		for (size_t j = 0; j < bubbleIndicesPerCluster.size(); j++)
		{
			if (j == i) continue;
			total += abs(getPhaseScore(readsPerBubbleAllele, bubbleIndicesPerCluster[i], bubbleIndicesPerCluster[j]));
		}
		if (total > maxTotal)
		{
			maxTotal = total;
			maxNode = i;
		}
	}
	result[maxNode] = true;
	handled[maxNode] = true;
	for (size_t iter = 1; iter < bubbleIndicesPerCluster.size(); iter++)
	{
		int maxTotal = 0;
		size_t maxNode = 0;
		bool hasMaxNode = false;
		for (size_t i = 0; i < bubbleIndicesPerCluster.size(); i++)
		{
			if (handled[i]) continue;
			int total = 0;
			for (size_t j = 0; j < bubbleIndicesPerCluster.size(); j++)
			{
				if (!handled[i]) continue;
				total += getPhaseScore(readsPerBubbleAllele, i, j);
			}
			if (!hasMaxNode || abs(total) > abs(maxTotal))
			{
				maxTotal = total;
				maxNode = i;
				hasMaxNode = true;
			}
		}
		assert(hasMaxNode);
		assert(!handled[maxNode]);
		if (maxTotal > 0)
		{
			result[maxNode] = true;
		}
		else
		{
			result[maxNode] = false;
		}
		handled[maxNode] = true;
	}
	return result;
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> orientBubblesAndGetPhaseBlocks(const std::vector<std::vector<size_t>>& bubbleIndicesPerCluster, const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& readsPerBubbleAllele, const std::vector<std::pair<size_t, size_t>>& simpleBubbles)
{
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> result;
	result.resize(bubbleIndicesPerCluster.size());
	for (size_t i = 0; i < bubbleIndicesPerCluster.size(); i++)
	{
		std::vector<bool> bubbleOrientations = getBubbleOrientations(bubbleIndicesPerCluster[i], readsPerBubbleAllele);
		for (size_t j = 0; j < bubbleOrientations.size(); j++)
		{
			if (bubbleOrientations[j])
			{
				result[i].first.push_back(simpleBubbles[bubbleIndicesPerCluster[i][j]].first);
				result[i].second.push_back(simpleBubbles[bubbleIndicesPerCluster[i][j]].second);
			}
			else
			{
				result[i].first.push_back(simpleBubbles[bubbleIndicesPerCluster[i][j]].second);
				result[i].second.push_back(simpleBubbles[bubbleIndicesPerCluster[i][j]].first);
			}
		}
	}
	return result;
}

void clearInvalidBubbles(std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& readsPerBubbleAllele, std::vector<std::pair<size_t, size_t>>& simpleBubbles)
{
	for (size_t i = readsPerBubbleAllele.size()-1; i < readsPerBubbleAllele.size(); i++)
	{
		if (intersectSize(readsPerBubbleAllele[i].first, readsPerBubbleAllele[i].second) < 2) continue;
		std::swap(readsPerBubbleAllele[i], readsPerBubbleAllele.back());
		readsPerBubbleAllele.pop_back();
		std::swap(simpleBubbles[i], simpleBubbles.back());
		simpleBubbles.pop_back();
	}
}

size_t findbidirected(std::vector<size_t>& parent, const size_t i)
{
	assert((i & firstBitUint64_t) == 0);
	if (parent[i] == i) return i;
	auto result = findbidirected(parent, parent[i] & maskUint64_t) ^ (parent[i] & firstBitUint64_t);
	parent[i] = result;
	return result;
}

void merge(std::vector<size_t>& parent, const size_t i, const size_t j, const bool fw)
{
	auto leftParent = findbidirected(parent, i);
	auto rightParent = findbidirected(parent, j);
	if ((leftParent & maskUint64_t) == (rightParent & maskUint64_t))
	{
		assert((leftParent & firstBitUint64_t) ^ (rightParent & firstBitUint64_t) ^ fw);
		return;
	}
	parent[rightParent & maskUint64_t] = (leftParent & maskUint64_t) + (((leftParent & firstBitUint64_t) ^ (rightParent & firstBitUint64_t) ^ fw) ? 0 : firstBitUint64_t);
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> mergePhaseBlocks(const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& blocks, const std::vector<ReadPathBundle>& readPaths)
{
	phmap::flat_hash_map<size_t, size_t> nodeToBlockAllele;
	for (size_t i = 0; i < blocks.size(); i++)
	{
		for (auto node : blocks[i].first)
		{
			assert(nodeToBlockAllele.count(node) == 0);
			nodeToBlockAllele[node] = i + firstBitUint64_t;
		}
		for (auto node : blocks[i].second)
		{
			assert(nodeToBlockAllele.count(node) == 0);
			nodeToBlockAllele[node] = i;
		}
	}
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> connectionCounts;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		phmap::flat_hash_set<size_t> blocksHere;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				if (nodeToBlockAllele.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 0) continue;
				blocksHere.insert(nodeToBlockAllele.at(readPaths[i].paths[j].path[k] & maskUint64_t));
			}
		}
		if (blocksHere.size() >= 2)
		{
			for (auto b1 : blocksHere)
			{
				for (auto b2 : blocksHere)
				{
					if (b1 < b2) connectionCounts[std::make_pair(b1, b2)] += 1;
					if (b2 < b1) connectionCounts[std::make_pair(b2, b1)] += 1;
				}
			}
		}
	}
	std::vector<size_t> parent;
	for (size_t i = 0; i < blocks.size(); i++)
	{
		parent.push_back(i);
	}
	for (size_t i = 0; i < blocks.size(); i++)
	{
		for (size_t j = i+1; j < blocks.size(); j++)
		{
			bool hasBothNormal = false;
			bool hasBothCross = false;
			bool hasAnyNormal = false;
			bool hasAnyCross = false;
			if (connectionCounts.count(std::make_pair(i, j)) == 1 && connectionCounts.count(std::make_pair(i + firstBitUint64_t, j + firstBitUint64_t)) == 1)
			{
				if (connectionCounts.at(std::make_pair(i, j)) >= 2 && connectionCounts.at(std::make_pair(i + firstBitUint64_t, j + firstBitUint64_t)) >= 2)
				{
					hasBothNormal = true;
				}
			}
			if (connectionCounts.count(std::make_pair(i + firstBitUint64_t, j)) == 1 && connectionCounts.count(std::make_pair(i, j + firstBitUint64_t)) == 1)
			{
				if (connectionCounts.at(std::make_pair(i + firstBitUint64_t, j)) >= 2 && connectionCounts.at(std::make_pair(i, j + firstBitUint64_t)) >= 2)
				{
					hasBothCross = true;
				}
			}
			if (connectionCounts.count(std::make_pair(i, j)) == 1) hasAnyNormal = true;
			if (connectionCounts.count(std::make_pair(i + firstBitUint64_t, j + firstBitUint64_t)) == 1) hasAnyNormal = true;
			if (connectionCounts.count(std::make_pair(i, j + firstBitUint64_t)) == 1) hasAnyCross = true;
			if (connectionCounts.count(std::make_pair(i + firstBitUint64_t, j)) == 1) hasAnyCross = true;
			if (hasBothNormal && !hasAnyCross)
			{
				merge(parent, i, j, true);
			}
			if (hasBothCross && !hasAnyNormal)
			{
				merge(parent, i, j, false);
			}
		}
	}
	RankBitvector kept;
	kept.resize(parent.size());
	for (size_t i = 0; i < parent.size(); i++)
	{
		if (findbidirected(parent, i) == i) kept.set(i, true);
	}
	kept.buildRanks();
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> result;
	result.resize(kept.getRank(kept.size()-1) + (kept.get(kept.size()-1) ? 1 : 0));
	for (size_t i = 0; i < parent.size(); i++)
	{
		size_t target = findbidirected(parent, i);
		if (target & firstBitUint64_t)
		{
			result[kept.getRank(target & maskUint64_t)].first.insert(result[kept.getRank(target & maskUint64_t)].first.end(), blocks[i].first.begin(), blocks[i].first.end());
			result[kept.getRank(target & maskUint64_t)].second.insert(result[kept.getRank(target & maskUint64_t)].second.end(), blocks[i].second.begin(), blocks[i].second.end());
		}
		else
		{
			result[kept.getRank(target & maskUint64_t)].second.insert(result[kept.getRank(target & maskUint64_t)].second.end(), blocks[i].first.begin(), blocks[i].first.end());
			result[kept.getRank(target & maskUint64_t)].first.insert(result[kept.getRank(target & maskUint64_t)].first.end(), blocks[i].second.begin(), blocks[i].second.end());
		}
	}
	return result;
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> getGraphPhaseBlockNodes(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	SparseEdgeContainer activeEdges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<std::pair<size_t, size_t>> simpleBubbles = getSimpleBubbles(activeEdges, unitigGraph, approxOneHapCoverage);
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> readsPerBubbleAllele = getReadsPerBubbleAllele(readPaths, simpleBubbles);
	clearInvalidBubbles(readsPerBubbleAllele, simpleBubbles);
	std::vector<std::vector<size_t>> bubbleIndicesPerCluster = mergeBubblesToClusters(simpleBubbles, readsPerBubbleAllele);
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> result = orientBubblesAndGetPhaseBlocks(bubbleIndicesPerCluster, readsPerBubbleAllele, simpleBubbles);
	result = mergePhaseBlocks(result, readPaths);
	return result;
}

std::pair<std::vector<size_t>, std::vector<size_t>> getReadsPerPhaseBlock(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::pair<std::vector<size_t>, std::vector<size_t>>& phaseBlock)
{
	phmap::flat_hash_set<size_t> leftNodes { phaseBlock.first.begin(), phaseBlock.first.end() };
	phmap::flat_hash_set<size_t> rightNodes { phaseBlock.second.begin(), phaseBlock.second.end() };
	std::pair<std::vector<size_t>, std::vector<size_t>> result;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		size_t leftScore = 0;
		size_t rightScore = 0;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				if (leftNodes.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 1)
				{
					size_t matchSize = unitigGraph.lengths[readPaths[i].paths[j].path[k] & maskUint64_t];
					if (k == 0) matchSize -= readPaths[i].paths[j].pathLeftClipKmers;
					if (k == readPaths[i].paths[j].path.size()-1) matchSize -= readPaths[i].paths[j].pathRightClipKmers;
					assert(matchSize >= 1);
					assert(matchSize <= unitigGraph.lengths[readPaths[i].paths[j].path[k] & maskUint64_t]);
					leftScore += matchSize;
				}
				else if (rightNodes.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 1)
				{
					size_t matchSize = unitigGraph.lengths[readPaths[i].paths[j].path[k] & maskUint64_t];
					if (k == 0) matchSize -= readPaths[i].paths[j].pathLeftClipKmers;
					if (k == readPaths[i].paths[j].path.size()-1) matchSize -= readPaths[i].paths[j].pathRightClipKmers;
					assert(matchSize >= 1);
					assert(matchSize <= unitigGraph.lengths[readPaths[i].paths[j].path[k] & maskUint64_t]);
					rightScore += matchSize;
				}
			}
		}
		if (leftScore > rightScore)
		{
			result.first.emplace_back(i);
		}
		else if (rightScore > leftScore)
		{
			result.second.emplace_back(i);
		}
	}
	return result;
}
