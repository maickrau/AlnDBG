#include <iostream>
#include <tuple>
#include "MBGCommon.h"
#include "GraphPhaser.h"
#include "Common.h"

std::pair<size_t, size_t> findSimpleBubble(const SparseEdgeContainer& activeEdges, const std::pair<size_t, bool> start, const UnitigGraph& graph, const double approxOneHapCoverage)
{
	if (activeEdges.getEdges(start).size() != 2) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.coverages[start.first] < approxOneHapCoverage * 1.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.coverages[start.first] > approxOneHapCoverage * 2.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	std::pair<size_t, bool> otherSideNode { std::numeric_limits<size_t>::max(), true };
	for (auto edge : activeEdges.getEdges(start))
	{
		if (graph.coverages[edge.first] < approxOneHapCoverage * 0.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (graph.coverages[edge.first] > approxOneHapCoverage * 1.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (activeEdges.getEdges(edge).size() != 1) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (activeEdges.getEdges(reverse(edge)).size() != 1) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (otherSideNode.first == std::numeric_limits<size_t>::max())
		{
			otherSideNode = activeEdges.getEdges(edge)[0];
		}
		else
		{
			if (activeEdges.getEdges(edge)[0] != otherSideNode) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		}
	}
	if (activeEdges.getEdges(reverse(otherSideNode)).size() != 2) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.coverages[otherSideNode.first] < approxOneHapCoverage * 1.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.coverages[otherSideNode.first] > approxOneHapCoverage * 2.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.lengths[start.first] < 100 && graph.lengths[otherSideNode.first] < 100) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	auto edges = activeEdges.getEdges(start);
	return std::make_pair(edges[0].first, edges[1].first);
}

std::pair<size_t, size_t> findHighCoverageFork(const SparseEdgeContainer& activeEdges, const std::pair<size_t, bool> start, const UnitigGraph& graph, const double approxOneHapCoverage)
{
	if (activeEdges.getEdges(start).size() != 2) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.coverages[start.first] < approxOneHapCoverage * 1.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.coverages[start.first] > approxOneHapCoverage * 2.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	if (graph.lengths[start.first] < 100) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	for (auto edge : activeEdges.getEdges(start))
	{
		if (graph.coverages[edge.first] < approxOneHapCoverage * 0.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (graph.coverages[edge.first] > approxOneHapCoverage * 1.5) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (activeEdges.getEdges(reverse(edge)).size() != 1) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		if (graph.lengths[edge.first] < 100) return std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	}
	auto edges = activeEdges.getEdges(start);
	return std::make_pair(edges[0].first, edges[1].first);
}

std::pair<size_t, size_t> findPhaseableStructure(const SparseEdgeContainer& activeEdges, const std::pair<size_t, bool> start, const UnitigGraph& graph, const double approxOneHapCoverage)
{
	auto found = findSimpleBubble(activeEdges, start, graph, approxOneHapCoverage);
	if (found.first != std::numeric_limits<size_t>::max()) return found;
	found = findHighCoverageFork(activeEdges, start, graph, approxOneHapCoverage);
	return found;
}

std::vector<std::pair<std::pair<size_t, size_t>, size_t>> getTripletCounts(const SparseEdgeContainer& activeEdges, const std::vector<ReadPathBundle>& readPaths, const size_t startNode)
{
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> counts;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 1; k+1 < readPaths[i].paths[j].path.size(); k++)
			{
				if ((readPaths[i].paths[j].path[k] & maskUint64_t) != startNode) continue;
				size_t prev = readPaths[i].paths[j].path[k-1];
				size_t next = readPaths[i].paths[j].path[k+1];
				if ((readPaths[i].paths[j].path[k] & firstBitUint64_t) == 0)
				{
					std::swap(next, prev);
					prev ^= firstBitUint64_t;
					next ^= firstBitUint64_t;
				}
				counts[std::make_pair(prev, next)] += 1;
				if (counts.size() > 2) return std::vector<std::pair<std::pair<size_t, size_t>, size_t>>{};
			}
		}
	}
	std::vector<std::pair<std::pair<size_t, size_t>, size_t>> result { counts.begin(), counts.end() };
	return result;
}

std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>> findPhaseableXShape(const SparseEdgeContainer& activeEdges, const UnitigGraph& graph, const std::vector<ReadPathBundle>& readPaths, const size_t startNode, const double approxOneHapCoverage)
{
	if (graph.coverages[startNode] < approxOneHapCoverage * 1.5) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
//	if (graph.coverages[startNode] > approxOneHapCoverage * 2.5) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	if (activeEdges.getEdges(std::make_pair(startNode, false)).size() != 2) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	if (activeEdges.getEdges(std::make_pair(startNode, true)).size() != 2) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	for (auto edge : activeEdges.getEdges(std::make_pair(startNode, false)))
	{
//		if (graph.coverages[edge.first] < approxOneHapCoverage * 0.5) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
//		if (graph.coverages[edge.first] > approxOneHapCoverage * 1.5) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
		if (activeEdges.getEdges(reverse(edge)).size() != 1) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	}
	for (auto edge : activeEdges.getEdges(std::make_pair(startNode, true)))
	{
//		if (graph.coverages[edge.first] < approxOneHapCoverage * 0.5) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
//		if (graph.coverages[edge.first] > approxOneHapCoverage * 1.5) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
		if (activeEdges.getEdges(reverse(edge)).size() != 1) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	}
	auto tripletCounts = getTripletCounts(activeEdges, readPaths, startNode);
	if (tripletCounts.size() != 2) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	if (tripletCounts[0].second < 5) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	if (tripletCounts[1].second < 5) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	if (tripletCounts[0].first.first == tripletCounts[1].first.first) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	if (tripletCounts[0].first.second == tripletCounts[1].first.first) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	if (tripletCounts[0].first.first == tripletCounts[1].first.second) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	if (tripletCounts[0].first.second == tripletCounts[1].first.second) return std::make_pair(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	return std::make_pair(std::make_pair(tripletCounts[0].first.first & maskUint64_t, tripletCounts[1].first.first & maskUint64_t), std::make_pair(tripletCounts[0].first.second & maskUint64_t, tripletCounts[1].first.second & maskUint64_t));
}

std::vector<std::pair<size_t, size_t>> getSimpleBubbles(const SparseEdgeContainer& activeEdges, const UnitigGraph& graph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	phmap::flat_hash_set<std::pair<size_t, size_t>> result;
	for (size_t i = 0; i < activeEdges.size(); i++)
	{
		auto fwBubble = findPhaseableStructure(activeEdges, std::make_pair(i, true), graph, approxOneHapCoverage);
		if (fwBubble.first < activeEdges.size())
		{
			assert(fwBubble.second < activeEdges.size());
			assert(fwBubble.first != fwBubble.second);
			result.emplace(std::min(fwBubble.first, fwBubble.second), std::max(fwBubble.first, fwBubble.second));
		}
		auto bwBubble = findPhaseableStructure(activeEdges, std::make_pair(i, false), graph, approxOneHapCoverage);
		if (bwBubble.first < activeEdges.size())
		{
			assert(bwBubble.second < activeEdges.size());
			assert(bwBubble.first != bwBubble.second);
			result.emplace(std::min(bwBubble.first, bwBubble.second), std::max(bwBubble.first, bwBubble.second));
		}
	}
	for (size_t i = 0; i < graph.nodeCount(); i++)
	{
		auto bubbles = findPhaseableXShape(activeEdges, graph, readPaths, i, approxOneHapCoverage);
		if (bubbles.first.first < activeEdges.size())
		{
			result.emplace(std::min(bubbles.first.first, bubbles.first.second), std::max(bubbles.first.first, bubbles.first.second));
			result.emplace(std::min(bubbles.second.first, bubbles.second.second), std::max(bubbles.second.first, bubbles.second.second));
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

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> mergePhaseBlocks(const UnitigGraph& unitigGraph, const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& blocks, const std::vector<ReadPathBundle>& readPaths)
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
		phmap::flat_hash_map<size_t, size_t> blockCountsHere;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				if (nodeToBlockAllele.count(readPaths[i].paths[j].path[k] & maskUint64_t) == 0) continue;
				size_t lengthHere = unitigGraph.lengths[readPaths[i].paths[j].path[k] & maskUint64_t];
				if (k == 0) lengthHere -= readPaths[i].paths[j].pathLeftClipKmers;
				if (k == readPaths[i].paths[j].path.size()-1) lengthHere -= readPaths[i].paths[j].pathRightClipKmers;
				assert(lengthHere <= unitigGraph.lengths[readPaths[i].paths[j].path[k] & maskUint64_t]);
				blockCountsHere[nodeToBlockAllele.at(readPaths[i].paths[j].path[k] & maskUint64_t)] += lengthHere;
			}
		}
		phmap::flat_hash_set<size_t> blocksHere;
		for (auto pair : blockCountsHere)
		{
			if (blockCountsHere.count(pair.first ^ firstBitUint64_t) == 1 && blockCountsHere.at(pair.first ^ firstBitUint64_t)+1 >= pair.second) continue;
			blocksHere.insert(pair.first);
		}
		if (blocksHere.size() >= 2)
		{
			for (auto b1 : blocksHere)
			{
				for (auto b2 : blocksHere)
				{
					if ((b1 & maskUint64_t) < (b2 & maskUint64_t))
					{
						connectionCounts[std::make_pair(b1, b2)] += 1;
					}
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
			size_t normalCounts = 0;
			size_t crossCounts = 0;
			if (connectionCounts.count(std::make_pair(i, j)) == 1) normalCounts += connectionCounts.at(std::make_pair(i, j));
			if (connectionCounts.count(std::make_pair(i + firstBitUint64_t, j + firstBitUint64_t)) == 1) normalCounts += connectionCounts.at(std::make_pair(i + firstBitUint64_t, j + firstBitUint64_t));
			if (connectionCounts.count(std::make_pair(i + firstBitUint64_t, j)) == 1) crossCounts += connectionCounts.at(std::make_pair(i + firstBitUint64_t, j));
			if (connectionCounts.count(std::make_pair(i, j + firstBitUint64_t)) == 1) crossCounts += connectionCounts.at(std::make_pair(i, j + firstBitUint64_t));
			if (normalCounts >= 3 && crossCounts == 0)
			{
				merge(parent, i, j, true);
			}
			if (crossCounts >= 3 && normalCounts == 0)
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

void collectMiddleNodes(phmap::flat_hash_set<size_t>& result, const std::vector<ReadPathBundle>& readPaths, const size_t read, const phmap::flat_hash_set<size_t>& hetNodes)
{
	size_t lastMatchJ = std::numeric_limits<size_t>::max();
	size_t lastMatchK = std::numeric_limits<size_t>::max();
	for (size_t j = 0; j < readPaths[read].paths.size(); j++)
	{
		for (size_t k = 0; k < readPaths[read].paths[j].path.size(); k++)
		{
			if (hetNodes.count(readPaths[read].paths[j].path[k] & maskUint64_t) == 0) continue;
			if (lastMatchJ != std::numeric_limits<size_t>::max())
			{
				assert(lastMatchK != std::numeric_limits<size_t>::max());
				for (size_t otherj = lastMatchJ; otherj <= j; otherj++)
				{
					for (size_t otherk = (otherj == lastMatchJ ? lastMatchK : 0); otherk < (otherj == j ? k : readPaths[read].paths[otherj].path.size()); otherk++)
					{
						result.insert(readPaths[read].paths[otherj].path[otherk] & maskUint64_t);
					}
				}
			}
			lastMatchJ = j;
			lastMatchK = k;
		}
	}
}

phmap::flat_hash_set<size_t> getNodesInPhaseBlock(const std::vector<ReadPathBundle>& readPaths, const std::pair<std::vector<size_t>, std::vector<size_t>>& readsPerPhaseBlock, const std::pair<std::vector<size_t>, std::vector<size_t>>& phaseBlock)
{
	phmap::flat_hash_set<size_t> result;
	phmap::flat_hash_set<size_t> hetNodes { phaseBlock.first.begin(), phaseBlock.first.end() };
	hetNodes.insert(phaseBlock.second.begin(), phaseBlock.second.end());
	for (auto read : readsPerPhaseBlock.first)
	{
		collectMiddleNodes(result, readPaths, read, hetNodes);
	}
	for (auto read : readsPerPhaseBlock.second)
	{
		collectMiddleNodes(result, readPaths, read, hetNodes);
	}
	result.insert(hetNodes.begin(), hetNodes.end());
	return result;
}

phmap::flat_hash_map<size_t, size_t> getMidNodeCounts(const std::vector<size_t>& reads, const std::vector<ReadPathBundle>& readPaths, const phmap::flat_hash_set<size_t>& hetNodes)
{
	phmap::flat_hash_map<size_t, size_t> result;
	for (size_t read : reads)
	{
		size_t lastMatchJ = std::numeric_limits<size_t>::max();
		size_t lastMatchK = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < readPaths[read].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[read].paths[j].path.size(); k++)
			{
				if (hetNodes.count(readPaths[read].paths[j].path[k] & maskUint64_t) == 0) continue;
				if (lastMatchJ == std::numeric_limits<size_t>::max())
				{
					lastMatchJ = j;
					lastMatchK = k;
					continue;
				}
				for (size_t otherj = lastMatchJ; otherj <= j; otherj++)
				{
					for (size_t otherk = (otherj == lastMatchJ ? lastMatchK : 0); otherk < (otherj == j ? k : readPaths[read].paths[otherj].path.size()); otherk++)
					{
						result[readPaths[read].paths[otherj].path[otherk] & maskUint64_t] += 1;
					}
				}
				lastMatchJ = j;
				lastMatchK = k;
			}
		}
	}
	return result;
}

void insertInternalNodes(std::pair<std::vector<size_t>, std::vector<size_t>>& result, const std::vector<ReadPathBundle>& readPaths, const std::pair<std::vector<size_t>, std::vector<size_t>>& readsPerHaplotype, const std::pair<std::vector<size_t>, std::vector<size_t>>& rawNodes)
{
	phmap::flat_hash_set<size_t> hetNodes;
	hetNodes.insert(rawNodes.first.begin(), rawNodes.first.end());
	hetNodes.insert(rawNodes.second.begin(), rawNodes.second.end());
	phmap::flat_hash_map<size_t, size_t> countInHap1 = getMidNodeCounts(readsPerHaplotype.first, readPaths, hetNodes);
	phmap::flat_hash_map<size_t, size_t> countInHap2 = getMidNodeCounts(readsPerHaplotype.second, readPaths, hetNodes);
	for (auto pair : countInHap1)
	{
		if (pair.second < 3) continue;
		if (countInHap2.count(pair.first) == 1) continue;
		if (hetNodes.count(pair.first) == 1) continue;
		result.first.push_back(pair.first);
	}
	for (auto pair : countInHap2)
	{
		if (pair.second < 3) continue;
		if (countInHap1.count(pair.first) == 1) continue;
		if (hetNodes.count(pair.first) == 1) continue;
		result.first.push_back(pair.first);
	}
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> fillInternalPhaseNodes(const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& rawNodes, const std::vector<ReadPathBundle>& readPaths, const UnitigGraph& unitigGraph)
{
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> result = rawNodes;
	for (size_t i = 0; i < rawNodes.size(); i++)
	{
		auto readsPerHaplotype = getReadsPerPhaseBlock(unitigGraph, readPaths, rawNodes[i]);
		std::cerr << "check internals " << i << std::endl;
		insertInternalNodes(result[i], readPaths, readsPerHaplotype, rawNodes[i]);
	}
	return result;
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> getGraphPhaseBlockNodes(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	SparseEdgeContainer activeEdges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<std::pair<size_t, size_t>> simpleBubbles = getSimpleBubbles(activeEdges, unitigGraph, readPaths, approxOneHapCoverage);
	if (simpleBubbles.size() == 0) return std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> {};
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> readsPerBubbleAllele = getReadsPerBubbleAllele(readPaths, simpleBubbles);
	clearInvalidBubbles(readsPerBubbleAllele, simpleBubbles);
	std::vector<std::vector<size_t>> bubbleIndicesPerCluster = mergeBubblesToClusters(simpleBubbles, readsPerBubbleAllele);
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> result = orientBubblesAndGetPhaseBlocks(bubbleIndicesPerCluster, readsPerBubbleAllele, simpleBubbles);
	// result = fillInternalPhaseNodes(result, readPaths, unitigGraph);
	std::ofstream info { "colors.csv" };
	info << "node,color,info" << std::endl;
	std::vector<std::string> groupcolors { "#00AA00", "#AA0000", "#AA00AA", "#0000AA", "#AAAA00", "#00AAAA", "#FF0000", "#00FF00", "#0000FF", "#AAFF00", "#AA00FF", "#00AAFF", "#888888" };
	for (size_t i = 0; i < result.size(); i++)
	{
		for (auto node : result[i].first) info << "node_" << node << "," << groupcolors[i % groupcolors.size()] << ",block_" << i << "_allele1" << std::endl;
		for (auto node : result[i].second) info << "node_" << node << "," << groupcolors[i % groupcolors.size()] << ",block_" << i << "_allele2" << std::endl;
	}
	// result = mergePhaseBlocks(unitigGraph, result, readPaths);
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> remapNodesToAlleleNodes(const UnitigGraph& oldGraph, const std::vector<ReadPathBundle>& readPaths, const phmap::flat_hash_map<size_t, size_t>& phaseBlockPerNode, const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& readsPerPhaseBlock)
{
	std::vector<ReadPathBundle> result = readPaths;
	UnitigGraph resultGraph = oldGraph;
	result.resize(readPaths.size());
	size_t nextNode = resultGraph.nodeCount();
	std::vector<phmap::flat_hash_map<size_t, size_t>> nodeToAlleleNodeHap1;
	std::vector<phmap::flat_hash_map<size_t, size_t>> nodeToAlleleNodeHap2;
	nodeToAlleleNodeHap1.resize(readsPerPhaseBlock.size());
	nodeToAlleleNodeHap2.resize(readsPerPhaseBlock.size());
	for (const auto& pair : phaseBlockPerNode)
	{
		resultGraph.coverages.push_back(0);
		resultGraph.coverages.push_back(0);
		resultGraph.lengths.push_back(oldGraph.lengths[pair.first]);
		nodeToAlleleNodeHap1[pair.second][pair.first] = nextNode;
		nextNode += 1;
		resultGraph.lengths.push_back(oldGraph.lengths[pair.first]);
		nodeToAlleleNodeHap2[pair.second][pair.first] = nextNode;
		nextNode += 1;
	}
	for (size_t blocki = 0; blocki < readsPerPhaseBlock.size(); blocki++)
	{
		for (size_t read : readsPerPhaseBlock[blocki].first)
		{
			for (size_t j = 0; j < result[read].paths.size(); j++)
			{
				for (size_t k = 0; k < result[read].paths[j].path.size(); k++)
				{
					if (nodeToAlleleNodeHap1[blocki].count(result[read].paths[j].path[k] & maskUint64_t) == 0) continue;
					result[read].paths[j].path[k] = nodeToAlleleNodeHap1[blocki].at(result[read].paths[j].path[k] & maskUint64_t) + (result[read].paths[j].path[k] & firstBitUint64_t);
				}
			}
		}
		for (size_t read : readsPerPhaseBlock[blocki].second)
		{
			for (size_t j = 0; j < result[read].paths.size(); j++)
			{
				for (size_t k = 0; k < result[read].paths[j].path.size(); k++)
				{
					if (nodeToAlleleNodeHap2[blocki].count(result[read].paths[j].path[k] & maskUint64_t) == 0) continue;
					result[read].paths[j].path[k] = nodeToAlleleNodeHap2[blocki].at(result[read].paths[j].path[k] & maskUint64_t) + (result[read].paths[j].path[k] & firstBitUint64_t);
				}
			}
		}
	}
	return std::make_pair(std::move(resultGraph), std::move(result));
}

std::vector<size_t> getContainedSubgraph(std::vector<bool>& checked, const SparseEdgeContainer& edges, const size_t start, const phmap::flat_hash_map<size_t, size_t>& nodeToPhaseBlock)
{
	assert(!checked[start]);
	assert(nodeToPhaseBlock.count(start) == 0);
	std::vector<size_t> checkStack;
	checkStack.push_back(start);
	phmap::flat_hash_set<size_t> boundingPhaseBlocks;
	std::vector<size_t> subgraph;
	while (checkStack.size() >= 1)
	{
		size_t top = checkStack.back();
		checkStack.pop_back();
		if (checked[top]) continue;
		if (nodeToPhaseBlock.count(top) == 1)
		{
			boundingPhaseBlocks.insert(nodeToPhaseBlock.at(top));
			continue;
		}
		checked[top] = true;
		subgraph.push_back(top);
		for (auto edge : edges.getEdges(std::make_pair(start, true)))
		{
			checkStack.push_back(edge.first);
		}
		for (auto edge : edges.getEdges(std::make_pair(start, false)))
		{
			checkStack.push_back(edge.first);
		}
	}
	if (boundingPhaseBlocks.size() != 1) subgraph.clear();
	return subgraph;
}

void removeContainedNonphasedErrorNodes(const UnitigGraph& unphasedGraph, const phmap::flat_hash_map<size_t, size_t>& nodeToPhaseBlock, RankBitvector& kept)
{
	auto edges = getActiveEdges(unphasedGraph.edgeCoverages, unphasedGraph.nodeCount());
	std::vector<bool> checked;
	checked.resize(unphasedGraph.nodeCount(), false);
	for (size_t i = 0; i < unphasedGraph.nodeCount(); i++)
	{
		if (!kept.get(i)) continue;
		if (checked[i]) continue;
		std::vector<size_t> subgraphNodes = getContainedSubgraph(checked, edges, i, nodeToPhaseBlock);
		size_t totalSize = 0;
		if (subgraphNodes.size() >= 10) continue;
		for (size_t node : subgraphNodes) totalSize += unphasedGraph.lengths[node];
		if (totalSize > 200) continue;
		for (size_t node : subgraphNodes) kept.set(node, false);
	}
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphBubbles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& phaseBlocks)
{
	std::cerr << phaseBlocks.size() << " phase blocks" << std::endl;
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> readsPerPhaseBlock;
	for (size_t i = 0; i < phaseBlocks.size(); i++)
	{
		readsPerPhaseBlock.emplace_back(getReadsPerPhaseBlock(unitigGraph, readPaths, phaseBlocks[i]));
	}
	phmap::flat_hash_map<size_t, std::vector<size_t>> phaseBlocksPerRead;
	for (size_t i = 0; i < readsPerPhaseBlock.size(); i++)
	{
		for (auto read : readsPerPhaseBlock[i].first) phaseBlocksPerRead[read].emplace_back(i);
		for (auto read : readsPerPhaseBlock[i].second) phaseBlocksPerRead[read].emplace_back(i);
	}
	size_t countReadsInMultiplePhaseBlocks = 0;
	for (auto pair : phaseBlocksPerRead)
	{
		if (pair.second.size() != 1) countReadsInMultiplePhaseBlocks += 1;
	}
	std::cerr << phaseBlocksPerRead.size() << " reads in phase blocks" << std::endl;
	std::cerr << countReadsInMultiplePhaseBlocks << " reads in multiple phase blocks" << std::endl;
	std::vector<phmap::flat_hash_set<size_t>> nodesInPhaseBlock;
	std::vector<std::vector<size_t>> phaseBlocksPerNode;
	phaseBlocksPerNode.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < phaseBlocks.size(); i++)
	{
		nodesInPhaseBlock.emplace_back(getNodesInPhaseBlock(readPaths, readsPerPhaseBlock[i], phaseBlocks[i]));
		for (auto node : nodesInPhaseBlock.back())
		{
			phaseBlocksPerNode[node].emplace_back(i);
		}
	}
	size_t countNodesInPhaseBlocks = 0;
	size_t countNodesInMultiplePhaseBlocks = 0;
	phmap::flat_hash_map<size_t, size_t> uniquePhaseBlockPerNode;
	for (size_t i = 0; i < phaseBlocksPerNode.size(); i++)
	{
		if (phaseBlocksPerNode[i].size() == 0) continue;
		if (phaseBlocksPerNode[i].size() == 1)
		{
			uniquePhaseBlockPerNode[i] = phaseBlocksPerNode[i][0];
		}
		countNodesInPhaseBlocks += 1;
		if (phaseBlocksPerNode[i].size() > 1)
		{
			countNodesInMultiplePhaseBlocks += 1;
			std::cerr << "node: " << i << std::endl;
		}
	}
	std::cerr << countNodesInPhaseBlocks << " nodes in phase blocks" << std::endl;
	std::cerr << countNodesInMultiplePhaseBlocks << " nodes in multiple phase blocks" << std::endl;
	if (countNodesInMultiplePhaseBlocks > 0) std::make_pair(unitigGraph, readPaths);
	std::vector<ReadPathBundle> remappedReads;
	UnitigGraph newGraph;
	std::tie(newGraph, remappedReads) = remapNodesToAlleleNodes(unitigGraph, readPaths, uniquePhaseBlockPerNode, readsPerPhaseBlock);
	RankBitvector kept;
	kept.resize(newGraph.nodeCount());
	for (size_t i = 0; i < kept.size(); i++)
	{
		kept.set(i, true);
	}
	for (const auto& node : uniquePhaseBlockPerNode)
	{
		kept.set(node.first, false);
	}
	// removeContainedNonphasedErrorNodes(unitigGraph, uniquePhaseBlockPerNode, kept);
	return filterUnitigGraph(newGraph, remappedReads, kept);
}

std::vector<uint64_t> getReachableEvenBackwards(const SparseEdgeContainer& edges, const phmap::flat_hash_set<size_t>& allowedNodes, const uint64_t start, const phmap::flat_hash_set<uint64_t>& boundingCoreNodes, const size_t forbiddenNode)
{
	std::vector<uint64_t> stack;
	phmap::flat_hash_set<uint64_t> checked;
	stack.push_back(start);
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (checked.count(top) == 1) continue;
		checked.insert(top);
		if (allowedNodes.count(top & maskUint64_t) == 0 && boundingCoreNodes.count(top) == 0) continue;
		if ((top & maskUint64_t) != forbiddenNode && boundingCoreNodes.count(top) == 0) stack.push_back(top ^ firstBitUint64_t);
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back((edge.first & maskUint64_t) + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	std::vector<uint64_t> result;
	for (auto node : boundingCoreNodes)
	{
		if (checked.count(node) == 0) continue;
		result.push_back(node);
	}
	return result;
}

std::pair<uint64_t, uint64_t> canonPointingInwards(uint64_t left, uint64_t right)
{
	return std::make_pair(std::min(left, right), std::max(left, right));
}

phmap::flat_hash_set<std::pair<size_t, size_t>> getSubgraphBoundaryConnections(const SparseEdgeContainer& edges, const phmap::flat_hash_set<size_t>& insideNodes, const phmap::flat_hash_set<uint64_t>& boundingCoreNodes)
{
	phmap::flat_hash_set<std::pair<size_t, size_t>> result;
	for (auto node : boundingCoreNodes)
	{
		auto reachable = getReachableEvenBackwards(edges, insideNodes, node, boundingCoreNodes, std::numeric_limits<size_t>::max() & maskUint64_t);
		for (auto node2 : reachable)
		{
			result.insert(canonPointingInwards(node, node2));
		}
	}
	return result;
}

phmap::flat_hash_set<std::pair<size_t, size_t>> getSubgraphBoundaryConnectionsWithoutOneNode(const SparseEdgeContainer& edges, const phmap::flat_hash_set<size_t>& insideNodes, const phmap::flat_hash_set<uint64_t>& boundingCoreNodes, const size_t forbiddenNode)
{
	phmap::flat_hash_set<std::pair<size_t, size_t>> result;
	for (auto node : boundingCoreNodes)
	{
		auto reachable = getReachableEvenBackwards(edges, insideNodes, node, boundingCoreNodes, forbiddenNode);
		for (auto node2 : reachable)
		{
			result.insert(canonPointingInwards(node, node2));
		}
	}
	return result;
}

void addTangleCoreNodes(std::vector<bool>& fwChecked, std::vector<bool>& bwChecked, std::vector<bool>& coreNode, const SparseEdgeContainer& edges, const size_t startNode)
{
	assert(!fwChecked[startNode]);
	assert(!bwChecked[startNode]);
	assert(!coreNode[startNode]);
	std::vector<size_t> stack;
	stack.push_back(startNode);
	phmap::flat_hash_set<size_t> insideNodes;
	phmap::flat_hash_set<uint64_t> boundingCoreNodes;
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (!coreNode[top & maskUint64_t]) insideNodes.insert(top & maskUint64_t);
		if (coreNode[top & maskUint64_t])
		{
			boundingCoreNodes.insert(top);
		}
		if ((top & firstBitUint64_t) && fwChecked[top & maskUint64_t]) continue;
		if (((top & firstBitUint64_t) == 0) && bwChecked[top & maskUint64_t]) continue;
		if (top & firstBitUint64_t)
		{
			fwChecked[top & maskUint64_t] = true;
		}
		else
		{
			bwChecked[top & maskUint64_t] = true;
		}
		if (!coreNode[top & maskUint64_t]) stack.push_back(top ^ firstBitUint64_t);
		auto edgesHere = edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t));
		for (auto edge : edgesHere)
		{
			stack.push_back((edge.first & maskUint64_t) + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	if (insideNodes.size() < 2) return;
	if (boundingCoreNodes.size() != 2) return;
	phmap::flat_hash_set<std::pair<size_t, size_t>> normalEdges = getSubgraphBoundaryConnections(edges, insideNodes, boundingCoreNodes);
	// if (normalEdges.size() != 1) return;
	for (auto node : insideNodes)
	{
		phmap::flat_hash_set<std::pair<size_t, size_t>> edgesHere = getSubgraphBoundaryConnectionsWithoutOneNode(edges, insideNodes, boundingCoreNodes, node);
		assert(edgesHere.size() <= normalEdges.size());
		if (edgesHere.size() < normalEdges.size())
		{
			coreNode[node] = true;;
		}
	}
}

void addCoreNodes(std::vector<bool>& coreNode, const SparseEdgeContainer& edges)
{
	std::vector<bool> fwChecked;
	std::vector<bool> bwChecked;
	fwChecked.resize(coreNode.size(), false);
	bwChecked.resize(coreNode.size(), false);
	for (size_t i = 0; i < coreNode.size(); i++)
	{
		if (coreNode[i]) continue;
		if (fwChecked[i])
		{
			assert(bwChecked[i]);
			continue;
		}
		assert(!bwChecked[i]);
		addTangleCoreNodes(fwChecked, bwChecked, coreNode, edges, i);
	}
}

phmap::flat_hash_set<uint64_t> getReachableCoreNodes(const SparseEdgeContainer& edges, const std::vector<bool>& coreNode, const uint64_t start)
{
	assert(coreNode[start & maskUint64_t]);
	phmap::flat_hash_set<uint64_t> result;
	phmap::flat_hash_set<uint64_t> checked;
	std::vector<uint64_t> stack;
	for (auto edge : edges.getEdges(std::make_pair(start & maskUint64_t, start & firstBitUint64_t)))
	{
		stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
	}
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (checked.count(top) == 1) continue;
		checked.insert(top);
		if (coreNode[top & maskUint64_t])
		{
			result.insert(top);
			continue;
		}
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
		}
	}
	return result;
}

SparseEdgeContainer getCoreEdges(const SparseEdgeContainer& edges, const std::vector<bool>& coreNode)
{
	SparseEdgeContainer coreEdges;
	coreEdges.resize(edges.size());
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> uniqueEdges;
	for (size_t i = 0; i < coreNode.size(); i++)
	{
		if (!coreNode[i]) continue;
		auto fwreachable = getReachableCoreNodes(edges, coreNode, i + firstBitUint64_t);
		auto bwreachable = getReachableCoreNodes(edges, coreNode, i);
		for (uint64_t node : fwreachable)
		{
			assert(coreNode[node & maskUint64_t]);
			uniqueEdges.emplace(canonPointingInwards(i + firstBitUint64_t, node ^ firstBitUint64_t));
		}
		for (uint64_t node : bwreachable)
		{
			assert(coreNode[node & maskUint64_t]);
			uniqueEdges.emplace(canonPointingInwards(i, node ^ firstBitUint64_t));
		}
	}
	for (auto pair : uniqueEdges)
	{
		std::pair<size_t, bool> from { pair.first & maskUint64_t, pair.first & firstBitUint64_t };
		std::pair<size_t, bool> to { pair.second & maskUint64_t, (pair.second & firstBitUint64_t) ^ firstBitUint64_t };
		coreEdges.addEdge(from, to);
		coreEdges.addEdge(reverse(to), reverse(from));
	}
	return coreEdges;
}

std::vector<uint64_t> getChainOneWay(const SparseEdgeContainer& coreEdges, std::vector<bool>& checked, const std::vector<bool>& coreNode, const uint64_t start)
{
	assert(coreNode[start & maskUint64_t]);
	assert(!checked[start & maskUint64_t]);
	std::vector<uint64_t> result;
	uint64_t pos = start;
	result.emplace_back(pos);
	while (coreEdges.getEdges(std::make_pair(pos & maskUint64_t, pos & firstBitUint64_t)).size() == 1)
	{
		std::pair<size_t, bool> nextpair = coreEdges.getEdges(std::make_pair(pos & maskUint64_t, pos & firstBitUint64_t))[0];
		if (coreEdges.getEdges(reverse(nextpair)).size() != 1) break;
		uint64_t next = nextpair.first + (nextpair.second ? firstBitUint64_t : 0);
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next == start) break;
		assert(!checked[next & maskUint64_t]);
		assert(coreNode[next & maskUint64_t]);
		pos = next;
		result.emplace_back(pos);
	}
	return result;
}

std::vector<uint64_t> getCoreChain(const SparseEdgeContainer& coreEdges, std::vector<bool>& checked, const std::vector<bool>& coreNode, const uint64_t start)
{
	assert((start & firstBitUint64_t) == 0);
	assert(!checked[start]);
	assert(coreNode[start]);
	std::vector<uint64_t> chainFw = getChainOneWay(coreEdges, checked, coreNode, start + firstBitUint64_t);
	std::vector<uint64_t> chainBw = getChainOneWay(coreEdges, checked, coreNode, start);
	assert(chainFw.size() >= 1);
	assert(chainBw.size() >= 1);
	assert(chainFw[0] == start + firstBitUint64_t);
	assert(chainBw[0] == start);
	if (chainFw.size() >= 2 && chainBw.size() >= 2)
	{
		if (chainFw.size() == chainBw.size())
		{
			if ((chainFw.back() ^ firstBitUint64_t) == chainBw[1])
			{
				// circular chain
				assert((chainBw.back() ^ firstBitUint64_t) == chainFw[1]);
				return chainFw;
			}
		}
	}
	std::reverse(chainBw.begin(), chainBw.end());
	for (size_t i = 0; i < chainBw.size(); i++)
	{
		chainBw[i] ^= firstBitUint64_t;
	}
	chainBw.insert(chainBw.end(), chainFw.begin()+1, chainFw.end());
	for (size_t i = 0; i < chainBw.size(); i++)
	{
		checked[chainBw[i] & maskUint64_t] = true;
	}
	assert(chainBw.size() >= 1);
	return chainBw;
}

std::vector<std::vector<uint64_t>> getCoreNodeChains(const SparseEdgeContainer& edges, const std::vector<bool>& coreNode)
{
	auto coreEdges = getCoreEdges(edges, coreNode);
	std::vector<bool> checked;
	checked.resize(coreNode.size(), false);
	std::vector<std::vector<uint64_t>> result;
	for (size_t i = 0; i < coreNode.size(); i++)
	{
		if (!coreNode[i]) continue;
		if (checked[i]) continue;
		result.emplace_back(getCoreChain(coreEdges, checked, coreNode, i));
	}
	std::vector<bool> found;
	found.resize(coreNode.size(), false);
	for (size_t i = 0 ; i < result.size(); i++)
	{
		for (auto node : result[i])
		{
			assert(!found[node & maskUint64_t]);
			found[node & maskUint64_t] = true;
		}
	}
	for (size_t i = 0; i < coreNode.size(); i++)
	{
		if (!coreNode[i])
		{
			assert(!found[i]);
		}
		else
		{
			assert(found[i]);
		}
	}
	return result;
}

size_t estimateCoreNodeChainPloidy(const std::vector<uint64_t> chain, const UnitigGraph& unitigGraph, const double approxOneHapCoverage)
{
	double coverageSum = 0;
	size_t lengthSum = 0;
	for (auto node : chain)
	{
		coverageSum += unitigGraph.coverages[node & maskUint64_t] * unitigGraph.lengths[node & maskUint64_t];
		lengthSum += unitigGraph.lengths[node & maskUint64_t];
	}
	assert(lengthSum > 0);
	size_t result = ((coverageSum / (double)lengthSum) / approxOneHapCoverage) + 0.5;
	return result;
}

size_t getAlleleIndex(std::vector<std::vector<std::vector<std::vector<uint64_t>>>>& allelesPerChain, const size_t i, const size_t j, const std::vector<uint64_t>& allele)
{
	for (size_t k = 0; k < allelesPerChain[i][j].size(); k++)
	{
		if (allele == allelesPerChain[i][j][k]) return k;
	}
	allelesPerChain[i][j].emplace_back(allele);
	return allelesPerChain[i][j].size()-1;
}

std::pair<std::vector<std::vector<std::vector<std::vector<uint64_t>>>>, std::vector<std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>>> getRawReadInfoPerChain(const std::vector<std::vector<uint64_t>>& coreNodeChains, const std::vector<std::vector<size_t>>& coreNodeOffsetInChain, const std::vector<ReadPathBundle>& readPaths, const UnitigGraph& unitigGraph)
{
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>> allelesPerReadPerChain;
	phmap::flat_hash_map<uint64_t, std::pair<size_t, size_t>> coreNodeLocator;
	allelesPerChain.resize(coreNodeChains.size());
	allelesPerReadPerChain.resize(coreNodeChains.size());
	for (size_t i = 0; i < coreNodeChains.size(); i++)
	{
		allelesPerChain[i].resize(coreNodeChains[i].size()-1);
		for (size_t j = 0; j < coreNodeChains[i].size(); j++)
		{
			coreNodeLocator[coreNodeChains[i][j] & maskUint64_t] = std::make_pair(i, j);
		}
	}
	for (size_t readi = 0; readi < readPaths.size(); readi++)
	{
		std::vector<std::tuple<size_t, int, size_t, size_t>> allelesInThisRead;
		for (size_t j = 0; j < readPaths[readi].paths.size(); j++)
		{
			size_t lastCore = std::numeric_limits<size_t>::max();
			size_t readPos = readPaths[readi].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[readi].paths[j].path.size(); k++)
			{
				readPos += unitigGraph.lengths[readPaths[readi].paths[j].path[k] & maskUint64_t];
				if (k == 0) readPos -= readPaths[readi].paths[j].pathLeftClipKmers;
				if (coreNodeLocator.count(readPaths[readi].paths[j].path[k] & maskUint64_t) == 0) continue;
				std::cerr << "allele between " << (readPaths[readi].paths[j].path[lastCore] & maskUint64_t) << " and " << (readPaths[readi].paths[j].path[k] & maskUint64_t) << std::endl;
				if (lastCore == std::numeric_limits<size_t>::max())
				{
					lastCore = k;
					continue;
				}
				std::pair<size_t, size_t> prev = coreNodeLocator.at(readPaths[readi].paths[j].path[lastCore] & maskUint64_t);
				std::pair<size_t, size_t> curr = coreNodeLocator.at(readPaths[readi].paths[j].path[k] & maskUint64_t);
				if (prev.first != curr.first)
				{
					lastCore = k;
					continue;
				}
				if (prev.second != curr.second+1 && curr.second != prev.second+1)
				{
					lastCore = k;
					continue;
				}
				bool fw = true;
				if (curr.second == prev.second + 1)
				{
					if (coreNodeChains[prev.first][prev.second] != readPaths[readi].paths[j].path[lastCore])
					{
						lastCore = k;
						continue;
					}
					if (coreNodeChains[curr.first][curr.second] != readPaths[readi].paths[j].path[k])
					{
						lastCore = k;
						continue;
					}
				}
				if (prev.second == curr.second + 1)
				{
					fw = false;
					if (coreNodeChains[prev.first][prev.second] != (readPaths[readi].paths[j].path[lastCore] ^ firstBitUint64_t))
					{
						lastCore = k;
						continue;
					}
					if (coreNodeChains[curr.first][curr.second] != (readPaths[readi].paths[j].path[k] ^ firstBitUint64_t))
					{
						lastCore = k;
						continue;
					}
				}
				std::vector<uint64_t> allele { readPaths[readi].paths[j].path.begin() + lastCore + 1, readPaths[readi].paths[j].path.begin() + k };
				if (!fw)
				{
					std::reverse(allele.begin(), allele.end());
					for (size_t l = 0; l < allele.size(); l++)
					{
						allele[l] ^= firstBitUint64_t;
					}
				}
				size_t index = getAlleleIndex(allelesPerChain, curr.first, std::min(curr.second, prev.second), allele);
				int diagonal = (int)readPos - (int)coreNodeOffsetInChain[curr.first][std::min(curr.second, prev.second)];
				if (!fw) diagonal = (int)readPos + (int)coreNodeOffsetInChain[curr.first][std::min(curr.second, prev.second)];
				allelesInThisRead.emplace_back(curr.first + (fw ? firstBitUint64_t : 0), diagonal, std::min(curr.second, prev.second), index);
				lastCore = k;
			}
		}
		std::sort(allelesInThisRead.begin(), allelesInThisRead.end());
		for (size_t i = 0; i < allelesInThisRead.size(); i++)
		{
			if (i == 0 || std::get<0>(allelesInThisRead[i]) != std::get<0>(allelesInThisRead[i-1]) || std::get<1>(allelesInThisRead[i]) > std::get<1>(allelesInThisRead[i-1]) + (int)coreNodeOffsetInChain[std::get<0>(allelesInThisRead[i])].back()/2)
			{
				allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].emplace_back();
				std::get<0>(allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back()) = i;
				std::get<1>(allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back()) = std::get<1>(allelesInThisRead[i]);
			}
			std::get<2>(allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back()).emplace_back(std::get<2>(allelesInThisRead[i]), std::get<3>(allelesInThisRead[i]));
		}
	}
	return std::make_pair(std::move(allelesPerChain), std::move(allelesPerReadPerChain));
}

void heapify(std::vector<std::pair<size_t, size_t>>& heap)
{
	// good enough, todo implement better
	std::sort(heap.begin(), heap.end());
}

size_t getShortestPathLength(const uint64_t start, const uint64_t end, const SparseEdgeContainer& edges, const UnitigGraph& graph)
{
	std::vector<std::pair<size_t, size_t>> minHeap;
	phmap::flat_hash_set<uint64_t> checked;
	minHeap.emplace_back(start, 0);
	while (minHeap.size() >= 1)
	{
		auto top = minHeap[0];
		std::swap(minHeap[0], minHeap.back());
		minHeap.pop_back();
		heapify(minHeap);
		if (top.first == end) return top.second;
		if (checked.count(top.first) == 1) continue;
		checked.insert(top.first);
		for (auto edge : edges.getEdges(std::make_pair(top.first & maskUint64_t, top.first & firstBitUint64_t)))
		{
			if (checked.count((edge.first & maskUint64_t) + (edge.second ? firstBitUint64_t : 0)) == 1) continue;
			minHeap.emplace_back((edge.first & maskUint64_t) + (edge.second ? firstBitUint64_t : 0), top.second + graph.lengths[top.first & maskUint64_t]);
		}
		std::sort(minHeap.begin(), minHeap.end());
	}
	assert(false);
}

std::vector<std::vector<size_t>> getCoreNodeOffsetsPerChain(const std::vector<std::vector<uint64_t>>& coreNodeChains, const SparseEdgeContainer& edges, const UnitigGraph& unitigGraph)
{
	std::vector<std::vector<size_t>> result;
	result.resize(coreNodeChains.size());
	for (size_t i = 0; i < coreNodeChains.size(); i++)
	{
		result[i].push_back(0);
		for (size_t j = 1; j < coreNodeChains[i].size(); j++)
		{
			result[i].push_back(result[i].back() + getShortestPathLength(coreNodeChains[i][j-1], coreNodeChains[i][j], edges, unitigGraph));
		}
	}
	return result;
}

std::vector<size_t> getConsensusAlleles(const std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>& allelesPerRead, const size_t alleleCount)
{
	std::vector<phmap::flat_hash_map<size_t, size_t>> alleleCoverage;
	alleleCoverage.resize(alleleCount);
	for (size_t readi = 0; readi < allelesPerRead.size(); readi++)
	{
		for (auto pair : std::get<2>(allelesPerRead[readi]))
		{
			assert(pair.first < alleleCoverage.size());
			alleleCoverage[pair.first][pair.second] += 1;
		}
	}
	std::vector<size_t> result;
	for (size_t i = 0; i < alleleCount; i++)
	{
		std::pair<size_t, size_t> maxResult { std::numeric_limits<size_t>::max(), 0 };
		for (auto pair : alleleCoverage[i])
		{
			if (pair.second > maxResult.second) maxResult = pair;
		}
		result.push_back(maxResult.first);
	}
	return result;
}

std::pair<std::vector<std::vector<size_t>>, size_t> getMostCoveredAlleles(const std::vector<std::vector<size_t>>& readsPerAllele, const size_t ploidy, const double approxOneHapCoverage)
{
	std::vector<std::pair<size_t, size_t>> coveragePerAllele;
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		coveragePerAllele.emplace_back(i, readsPerAllele[i].size());
	}
	std::sort(coveragePerAllele.begin(), coveragePerAllele.end(), [](auto left, auto right) { return left.second > right.second; });
	size_t estimatedTotalPloidy = 0;
	for (size_t i = 0; i < coveragePerAllele.size(); i++)
	{
		estimatedTotalPloidy += (double)coveragePerAllele[i].second / approxOneHapCoverage + 0.5;
	}
	std::vector<std::vector<size_t>> result;
	if (estimatedTotalPloidy != ploidy)
	{
		return std::make_pair(result, 0);
	}
	size_t gotPloidy = 0;
	size_t i = 0;
	for (i = 0; i < coveragePerAllele.size(); i++)
	{
		gotPloidy += (double)coveragePerAllele[i].second / approxOneHapCoverage + 0.5;
		result.emplace_back(readsPerAllele[coveragePerAllele[i].first]);
		std::sort(result.back().begin(), result.back().end());
		if (gotPloidy == estimatedTotalPloidy)
		{
			i += 1;
			break;
		}
	}
	size_t minorCount = 0;
	if (i < coveragePerAllele.size())
	{
		for (; i < coveragePerAllele.size(); i++)
		{
			minorCount += readsPerAllele[coveragePerAllele[i].first].size();
		}
	}
	return std::make_pair(result, minorCount);
}

bool allelesCouldMatch(const std::vector<std::vector<size_t>>& bestAllelesLeft, const std::vector<std::vector<size_t>>& bestAllelesRight, const size_t ploidy, const double approxOneHapCoverage)
{
	assert(ploidy == 2);
	assert(bestAllelesLeft.size() >= 2);
	assert(bestAllelesRight.size() >= 2);
	size_t requiredCoverage = (double)approxOneHapCoverage * ((double)ploidy - 0.5);
	if (intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[1]) >= requiredCoverage) return true;
	if (intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[0]) >= requiredCoverage) return true;
	return false;
}

std::vector<bool> getPossiblyHaplotypeInformativeSites(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& coreNodeChain, const size_t ploidy, const double approxOneHapCoverage)
{
	std::vector<std::vector<std::vector<size_t>>> bestAlleles;
	bestAlleles.resize(readsPerAllele.size());
	std::vector<bool> notgood;
	notgood.resize(readsPerAllele.size(), false);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		size_t minorCount = 0;
		std::tie(bestAlleles[i], minorCount) = getMostCoveredAlleles(readsPerAllele[i], ploidy, approxOneHapCoverage);
		if (bestAlleles[i].size() == 0)
		{
			notgood[i] = true;
			continue;
		}
		if (bestAlleles[i].size() < ploidy)
		{
			notgood[i] = true;
			continue;
		}
		size_t totalCoverage = 0;
		for (size_t j = 0; j < bestAlleles[i].size(); j++)
		{
			totalCoverage += bestAlleles[i][j].size();
		}
		if (minorCount*10 > totalCoverage)
		{
			notgood[i] = true;
			continue;
		}
	}
	std::vector<bool> result;
	result.resize(readsPerAllele.size(), false);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		if (notgood[i]) continue;
		for (size_t j = i+1; j < readsPerAllele.size(); j++)
		{
			if (notgood[j]) continue;
			if (allelesCouldMatch(bestAlleles[i], bestAlleles[j], ploidy, approxOneHapCoverage))
			{
				result[i] = true;
				result[j] = true;
			}
		}
	}
	return result;
}

std::vector<std::vector<size_t>> getAllelesPerHaplotype(const std::vector<size_t>& readAssignment, const size_t ploidy, const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& readOrder)
{
	assert(ploidy >= 2);
	std::vector<std::vector<phmap::flat_hash_map<size_t, size_t>>> alleleCoveragePerHaplotype;
	alleleCoveragePerHaplotype.resize(ploidy);
	for (size_t i = 0; i < ploidy; i++)
	{
		alleleCoveragePerHaplotype[i].resize(readsPerAllele.size());
	}
	for (size_t site = 0; site < readsPerAllele.size(); site++)
	{
		for (size_t allele = 0; allele < readsPerAllele[site].size(); allele++)
		{
			for (size_t read : readsPerAllele[site][allele])
			{
				if (readOrder[read] >= readAssignment.size()) continue;
				if (readAssignment[readOrder[read]] == std::numeric_limits<size_t>::max()) continue;
				alleleCoveragePerHaplotype[readAssignment[readOrder[read]]][site][allele] += 1;
			}
		}
	}
	std::vector<std::vector<size_t>> result;
	result.resize(ploidy);
	for (size_t hap = 0; hap < ploidy; hap++)
	{
		for (size_t site = 0; site < readsPerAllele.size(); site++)
		{
			std::pair<size_t, size_t> bestHere { std::numeric_limits<size_t>::max(), 0 };
			size_t totalCount = 0;
			for (auto pair : alleleCoveragePerHaplotype[hap][site])
			{
				if (pair.second > bestHere.second) bestHere = pair;
				totalCount += pair.second;
			}
			if (bestHere.second*10 < totalCount)
			{
				// very ambiguous, assume site is sequencing error
				bestHere.first = std::numeric_limits<size_t>::max();
			}
			result[hap].push_back(bestHere.first);
		}
	}
	return result;
}

size_t getHaplotypeScore(const std::vector<size_t>& readAssignment, const size_t ploidy, const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& readOrder)
{
	assert(readOrder.size() >= 1);
	assert(ploidy >= 2);
	size_t score = 0;
	std::vector<std::vector<size_t>> allelesPerHaplotype = getAllelesPerHaplotype(readAssignment, ploidy, readsPerAllele, readOrder);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		for (size_t allele = 0; allele < readsPerAllele[i].size(); allele++)
		{
			for (size_t read : readsPerAllele[i][allele])
			{
				if (readOrder[read] >= readAssignment.size()) continue;
				if (readAssignment[readOrder[read]] == std::numeric_limits<size_t>::max())
				{
					score += 1;
					continue;
				}
				size_t haplotypeAllele = allelesPerHaplotype[readAssignment[readOrder[read]]][i];
				if (haplotypeAllele == std::numeric_limits<size_t>::max())
				{
					score += 1;
				}
				else if (haplotypeAllele != allele)
				{
					score += 10;
				}
			}
		}
	}
	return score;
}

class ReadPartition
{
public:
	std::vector<size_t> readAssignment;
	size_t maxAssignmentPlusOne;
	size_t score;
};

std::vector<size_t> getUnweightedHeuristicMEC(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const size_t ploidy, const std::vector<size_t>& readOrder, const size_t realReadCount)
{
	const size_t beamWidth = 10000;
	size_t addEverythingUntilHere = log(beamWidth)/log(ploidy+1);
	assert(ploidy >= 2);
	std::vector<ReadPartition> activeAssignments;
	activeAssignments.emplace_back();
	activeAssignments.back().maxAssignmentPlusOne = 0;
	activeAssignments.back().score = 0;
	for (size_t i = 0; i < addEverythingUntilHere && i < realReadCount; i++)
	{
		std::vector<ReadPartition> nextActiveAssignments;
		for (const auto& haplotype : activeAssignments)
		{
			for (size_t j = 0; j < std::min(haplotype.maxAssignmentPlusOne+1, ploidy); j++)
			{
				nextActiveAssignments.emplace_back();
				nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
				nextActiveAssignments.back().maxAssignmentPlusOne = std::max(haplotype.maxAssignmentPlusOne, j+1);
				nextActiveAssignments.back().readAssignment.push_back(j);
			}
			nextActiveAssignments.emplace_back();
			nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
			nextActiveAssignments.back().maxAssignmentPlusOne = haplotype.maxAssignmentPlusOne;
			nextActiveAssignments.back().readAssignment.push_back(std::numeric_limits<size_t>::max());
		}
		activeAssignments = nextActiveAssignments;
	}
	for (size_t i = 0; i < activeAssignments.size(); i++)
	{
		activeAssignments[i].score = getHaplotypeScore(activeAssignments[i].readAssignment, ploidy, readsPerAllele, readOrder);
	}
	std::sort(activeAssignments.begin(), activeAssignments.end(), [](const auto& left, const auto& right) { return left.score < right.score; });
	for (size_t i = addEverythingUntilHere; i < realReadCount; i++)
	{
		std::vector<ReadPartition> nextActiveAssignments;
		for (const auto& haplotype : activeAssignments)
		{
			for (size_t j = 0; j < std::min(haplotype.maxAssignmentPlusOne+1, ploidy); j++)
			{
				nextActiveAssignments.emplace_back();
				nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
				nextActiveAssignments.back().readAssignment.push_back(j);
				nextActiveAssignments.back().maxAssignmentPlusOne = std::max(haplotype.maxAssignmentPlusOne, j+1);
				nextActiveAssignments.back().score = getHaplotypeScore(nextActiveAssignments.back().readAssignment, ploidy, readsPerAllele, readOrder);
			}
			nextActiveAssignments.emplace_back();
			nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
			nextActiveAssignments.back().maxAssignmentPlusOne = haplotype.maxAssignmentPlusOne;
			nextActiveAssignments.back().readAssignment.push_back(std::numeric_limits<size_t>::max());
			nextActiveAssignments.back().score = getHaplotypeScore(nextActiveAssignments.back().readAssignment, ploidy, readsPerAllele, readOrder);
		}
		std::sort(nextActiveAssignments.begin(), nextActiveAssignments.end(), [](const auto& left, const auto& right) { return left.score < right.score; });
		if (nextActiveAssignments.size() > beamWidth)
		{
			nextActiveAssignments.erase(nextActiveAssignments.begin()+beamWidth, nextActiveAssignments.end());
		}
		activeAssignments = nextActiveAssignments;
	}
	assert(activeAssignments[0].readAssignment.size() == realReadCount);
	std::cerr << "MEC best scores:";
	for (size_t i = 0; i < 10 && i < activeAssignments.size(); i++)
	{
		std::cerr << " " << activeAssignments[i].score;
	}
	std::cerr << std::endl;
	std::cerr << "best two assignments:" << std::endl;
	for (size_t i = 0; i < activeAssignments[0].readAssignment.size(); i++)
	{
		if (activeAssignments[0].readAssignment[i] == std::numeric_limits<size_t>::max())
		{
			std::cerr << "_";
		}
		else
		{
			std::cerr << activeAssignments[0].readAssignment[i];
		}
	}
	std::cerr << std::endl;
	for (size_t i = 0; i < activeAssignments[1].readAssignment.size(); i++)
	{
		if (activeAssignments[1].readAssignment[i] == std::numeric_limits<size_t>::max())
		{
			std::cerr << "_";
		}
		else
		{
			std::cerr << activeAssignments[1].readAssignment[i];
		}
	}
	std::cerr << std::endl;
	std::cerr << "best alleles:" << std::endl;
	std::vector<std::vector<size_t>> allelesPerHaplotype = getAllelesPerHaplotype(activeAssignments[0].readAssignment, ploidy, readsPerAllele, readOrder);
	for (size_t i = 0; i < allelesPerHaplotype.size(); i++)
	{
		for (size_t j = 0; j < allelesPerHaplotype[i].size(); j++)
		{
			if (allelesPerHaplotype[i][j] == std::numeric_limits<size_t>::max())
			{
				std::cerr << "_,";
			}
			else
			{
				std::cerr << allelesPerHaplotype[i][j] << ",";
			}
		}
		std::cerr << std::endl;
	}
	std::cerr << "last read alleles:" << std::endl;
	for (size_t site = 0; site < readsPerAllele.size(); site++)
	{
		size_t lastReadAlleleHere = std::numeric_limits<size_t>::max();
		for (size_t allele = 0; allele < readsPerAllele[site].size(); allele++)
		{
			for (size_t read : readsPerAllele[site][allele])
			{
				if (readOrder[read] == activeAssignments[0].readAssignment.size()-1)
				{
					lastReadAlleleHere = allele;
				}
			}
		}
		if (lastReadAlleleHere == std::numeric_limits<size_t>::max())
		{
			std::cerr << "_,";
		}
		else
		{
			std::cerr << lastReadAlleleHere << ",";
		}
	}
	std::cerr << std::endl;
	return activeAssignments[0].readAssignment;
}

std::vector<size_t> getReadOrder(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const size_t readCount)
{
	std::vector<size_t> firstPerRead;
	std::vector<size_t> lastPerRead;
	firstPerRead.resize(readCount, std::numeric_limits<size_t>::max());
	lastPerRead.resize(readCount, std::numeric_limits<size_t>::max());
	for (size_t site = 0; site < readsPerAllele.size(); site++)
	{
		for (size_t allele = 0; allele < readsPerAllele[site].size(); allele++)
		{
			for (size_t read : readsPerAllele[site][allele])
			{
				if (firstPerRead[read] == std::numeric_limits<size_t>::max())
				{
					firstPerRead[read] = site;
				}
				lastPerRead[read] = site;
			}
		}
	}
	std::vector<std::tuple<size_t, size_t, size_t>> ordering;
	for (size_t i = 0; i < readCount; i++)
	{
		if (firstPerRead[i] == std::numeric_limits<size_t>::max()) continue;
		assert(lastPerRead[i] != std::numeric_limits<size_t>::max());
		ordering.emplace_back(i, firstPerRead[i], lastPerRead[i]);
	}
	std::sort(ordering.begin(), ordering.end(), [](auto left, auto right)
	{
		if (std::get<1>(left) < std::get<1>(right)) return true;
		if (std::get<1>(left) > std::get<1>(right)) return false;
		if (std::get<2>(left) < std::get<2>(right)) return true;
		if (std::get<2>(left) > std::get<2>(right)) return false;
		return false;
	});
	std::vector<size_t> result;
	result.resize(readCount, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < ordering.size(); i++)
	{
		assert(result[std::get<0>(ordering[i])] == std::numeric_limits<size_t>::max());
		result[std::get<0>(ordering[i])] = i;
	}
	std::vector<bool> found;
	found.resize(ordering.size(), false);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		for (size_t j = 0; j < readsPerAllele[i].size(); j++)
		{
			for (size_t read : readsPerAllele[i][j])
			{
				assert(result[read] != std::numeric_limits<size_t>::max());
				found[result[read]] = true;
			}
		}
	}
	for (size_t i = 0; i < found.size(); i++)
	{
		assert(found[i]);
	}
	return result;
}

std::vector<std::vector<std::vector<size_t>>> getInformativeReads(const std::vector<std::vector<std::vector<size_t>>>& unfiltered, const std::vector<bool>& maybeHaplotypeInformative)
{
	assert(unfiltered.size() == maybeHaplotypeInformative.size());
	phmap::flat_hash_map<size_t, size_t> readAlleleCount;
	for (size_t site = 0; site < maybeHaplotypeInformative.size(); site++)
	{
		if (!maybeHaplotypeInformative[site]) continue;
		for (size_t allele = 0; allele < unfiltered[site].size(); allele++)
		{
			for (auto read : unfiltered[site][allele])
			{
				readAlleleCount[read] += 1;
			}
		}
	}
	std::vector<std::vector<std::vector<size_t>>> result;
	for (size_t site = 0; site < maybeHaplotypeInformative.size(); site++)
	{
		if (!maybeHaplotypeInformative[site]) continue;
		result.emplace_back();
		for (size_t allele = 0; allele < unfiltered[site].size(); allele++)
		{
			result.back().emplace_back();
			for (auto read : unfiltered[site][allele])
			{
				if (readAlleleCount[read] < 2) continue;
				result.back().back().emplace_back(read);
			}
		}
	}
	return result;
}

std::vector<std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>> phaseReadsToHaplotypes(const std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>& allelesPerRead, const std::vector<uint64_t>& coreNodeChain, const size_t ploidy, const size_t alleleCount, const double approxOneHapCoverage)
{
	assert(ploidy >= 2);
	assert(alleleCount >= 1);
	size_t readCount = allelesPerRead.size();
	std::vector<std::vector<std::vector<size_t>>> unfilteredReadsPerAllele;
	unfilteredReadsPerAllele.resize(alleleCount);
	for (size_t readi = 0; readi < allelesPerRead.size(); readi++)
	{
		for (auto pair : std::get<2>(allelesPerRead[readi]))
		{
			assert(pair.first < unfilteredReadsPerAllele.size());
			while (pair.second >= unfilteredReadsPerAllele[pair.first].size()) unfilteredReadsPerAllele[pair.first].emplace_back();
			unfilteredReadsPerAllele[pair.first][pair.second].emplace_back(readi);
		}
	}
	std::vector<bool> maybeHaplotypeInformative = getPossiblyHaplotypeInformativeSites(unfilteredReadsPerAllele, coreNodeChain, ploidy, approxOneHapCoverage);
	size_t countMaybeInformativeSites = 0;
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		if (maybeHaplotypeInformative[i]) countMaybeInformativeSites += 1;
	}
	if (countMaybeInformativeSites < 2)
	{
		std::cerr << "num informative sites " << countMaybeInformativeSites << ", skipped" << std::endl;
		return std::vector<std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>>{};
	}
	std::vector<std::vector<std::vector<size_t>>> readsPerAllele = getInformativeReads(unfilteredReadsPerAllele, maybeHaplotypeInformative);
	std::vector<size_t> readOrder = getReadOrder(readsPerAllele, readCount);
	size_t realReadCount = 0;
	for (size_t i = 0; i < readOrder.size(); i++)
	{
		if (readOrder[i] != std::numeric_limits<size_t>::max()) realReadCount += 1;
	}
	std::cerr << realReadCount << " real reads" << std::endl;
	std::vector<size_t> readAssignments = getUnweightedHeuristicMEC(readsPerAllele, ploidy, readOrder, realReadCount);
	std::vector<std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>> splittedReads;
	splittedReads.resize(ploidy);
	for (size_t i = 0; i < readAssignments.size(); i++)
	{
		if (readAssignments[i] == std::numeric_limits<size_t>::max()) continue;
		splittedReads[readAssignments[i]].emplace_back(allelesPerRead[i]);
	}
	return splittedReads;
}

std::vector<std::vector<std::vector<size_t>>> phaseCoreChains(const std::vector<std::vector<uint64_t>>& coreNodeChains, const std::vector<size_t>& chainPloidies, const std::vector<std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>>& allelesPerReadPerChain, const double approxOneHapCoverage)
{
	assert(coreNodeChains.size() == chainPloidies.size());
	assert(coreNodeChains.size() == allelesPerReadPerChain.size());
	std::vector<std::vector<std::vector<size_t>>> result;
	result.resize(coreNodeChains.size());
	for (size_t i = 0; i < coreNodeChains.size(); i++)
	{
		if (coreNodeChains[i].size() < 2) continue;
		std::cerr << "ploidy " << chainPloidies[i] << std::endl;
		if (chainPloidies[i] == 1)
		{
			result[i].push_back(getConsensusAlleles(allelesPerReadPerChain[i], coreNodeChains[i].size()-1));
			continue;
		}
		else if (chainPloidies[i] == 2)
		{
			auto phasedReads = phaseReadsToHaplotypes(allelesPerReadPerChain[i], coreNodeChains[i], chainPloidies[i], coreNodeChains[i].size()-1, approxOneHapCoverage);
			if (phasedReads.size() == 0) continue;
			std::cerr << "phased" << std::endl;
			for (size_t j = 0; j < phasedReads.size(); j++)
			{
				result[i].push_back(getConsensusAlleles(phasedReads[j], coreNodeChains[i].size()-1));
			}
		}
		else
		{
			std::cerr << "skipped chain with ploidy " << chainPloidies[i] << std::endl;
		}
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLinearizable(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	auto edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<bool> coreNode;
	coreNode.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.lengths[i] < 100) continue;
		coreNode[i] = true;
	}
	addCoreNodes(coreNode, edges);
	std::vector<std::vector<uint64_t>> coreNodeChains = getCoreNodeChains(edges, coreNode);
	std::vector<size_t> chainPloidies;
	std::cerr << coreNodeChains.size() << " core node chains" << std::endl;
	for (size_t i = 0; i < coreNodeChains.size(); i++)
	{
		size_t ploidy = estimateCoreNodeChainPloidy(coreNodeChains[i], unitigGraph, approxOneHapCoverage);
		chainPloidies.emplace_back(ploidy);
	}
	std::vector<std::vector<size_t>> coreNodeOffsets = getCoreNodeOffsetsPerChain(coreNodeChains, edges, unitigGraph);
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>> allelesPerReadPerChain;
	std::tie(allelesPerChain, allelesPerReadPerChain) = getRawReadInfoPerChain(coreNodeChains, coreNodeOffsets, readPaths, unitigGraph);
	for (size_t i = 0; i < coreNodeChains.size(); i++)
	{
		for (size_t j = 0; j+1 < coreNodeChains[i].size(); j++)
		{
			std::cerr << "chain " << i << " ploidy " << chainPloidies[i] << " node_" << (coreNodeChains[i][j] & maskUint64_t) << " node_" << (coreNodeChains[i][j+1] & maskUint64_t) << " allele count " << allelesPerChain[i][j].size() << std::endl;
		}
	}
	std::vector<std::vector<std::vector<size_t>>> chainHaplotypes = phaseCoreChains(coreNodeChains, chainPloidies, allelesPerReadPerChain, approxOneHapCoverage);
	return std::make_pair(unitigGraph, readPaths);
}
