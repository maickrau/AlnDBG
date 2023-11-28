#include <iostream>
#include <tuple>
#include "MBGCommon.h"
#include "GraphPhaser.h"
#include "Common.h"
#include "AnchorFinder.h"

class PhaseBlock
{
public:
	size_t chainNumber;
	std::vector<size_t> bubbleIndices;
	std::vector<std::vector<size_t>> allelesPerHaplotype;
	bool chainStartPhased;
	bool chainEndPhased;
};

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

size_t getAlleleIndex(std::vector<std::vector<std::vector<std::vector<uint64_t>>>>& allelesPerChain, const size_t i, const size_t j, const std::vector<uint64_t>& allele)
{
	for (size_t k = 0; k < allelesPerChain[i][j].size(); k++)
	{
		if (allele == allelesPerChain[i][j][k]) return k;
	}
	allelesPerChain[i][j].emplace_back(allele);
	return allelesPerChain[i][j].size()-1;
}

std::pair<std::vector<std::vector<std::vector<std::vector<uint64_t>>>>, std::vector<std::vector<std::tuple<size_t, bool, int, std::vector<std::pair<size_t, size_t>>>>>> getRawReadInfoPerChain(const std::vector<AnchorChain>& anchorChains, const std::vector<ReadPathBundle>& readPaths, const UnitigGraph& unitigGraph, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads)
{
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<std::tuple<size_t, bool, int, std::vector<std::pair<size_t, size_t>>>>> allelesPerReadPerChain;
	phmap::flat_hash_map<uint64_t, std::pair<size_t, size_t>> coreNodeLocator;
	allelesPerChain.resize(anchorChains.size());
	allelesPerReadPerChain.resize(anchorChains.size());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		allelesPerChain[i].resize(anchorChains[i].nodes.size()+1);
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			coreNodeLocator[anchorChains[i].nodes[j] & maskUint64_t] = std::make_pair(i, j);
		}
	}
	phmap::flat_hash_map<size_t, std::tuple<size_t, size_t, size_t>> simpleBubbleAllele;
	{
		SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
		for (size_t i = 0; i < anchorChains.size(); i++)
		{
			for (size_t j = 1; j < anchorChains[i].nodes.size(); j++)
			{
				const uint64_t prevNode = anchorChains[i].nodes[j-1];
				const uint64_t thisNode = anchorChains[i].nodes[j];
				if (edges.getEdges(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)).size() != 2) continue;
				if (edges.getEdges(std::make_pair(thisNode & maskUint64_t, (thisNode & firstBitUint64_t) ^ firstBitUint64_t)).size() != 2) continue;
				bool valid = true;
				for (auto edge : edges.getEdges(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)))
				{
					if (edges.getEdges(edge).size() != 1)
					{
						valid = false;
						break;
					}
					if (edges.getEdges(reverse(edge)).size() != 1)
					{
						valid = false;
						break;
					}
					if (edges.getEdges(edge)[0] != std::pair<size_t, bool>{ thisNode & maskUint64_t, thisNode & firstBitUint64_t })
					{
						valid = false;
						break;
					}
				}
				if (!valid) continue;
				std::cerr << "found simple bubble between " << (prevNode & maskUint64_t) << " " << (thisNode & maskUint64_t) << std::endl;
				for (auto edge : edges.getEdges(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)))
				{
					std::vector<uint64_t> allele;
					allele.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
					size_t index = getAlleleIndex(allelesPerChain, i, j, allele);
					assert(simpleBubbleAllele.count(edge.first) == 0);
					simpleBubbleAllele[edge.first] = std::make_tuple(i, j, index);
				}
			}
		}
	}
	for (size_t readi = 0; readi < readPaths.size(); readi++)
	{
		std::vector<std::tuple<size_t, int, size_t, size_t>> allelesInThisRead;
		// bubble gapless alleles
		for (size_t j = 0; j < readPaths[readi].paths.size(); j++)
		{
			size_t lastCore = std::numeric_limits<size_t>::max();
			size_t readPos = readPaths[readi].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[readi].paths[j].path.size(); k++)
			{
				readPos += unitigGraph.lengths[readPaths[readi].paths[j].path[k] & maskUint64_t];
				if (k == 0)
				{
					assert(readPos > readPaths[readi].paths[j].pathLeftClipKmers);
					readPos -= readPaths[readi].paths[j].pathLeftClipKmers;
				}
				if (coreNodeLocator.count(readPaths[readi].paths[j].path[k] & maskUint64_t) == 0) continue;
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
					if (anchorChains[prev.first].nodes[prev.second] != readPaths[readi].paths[j].path[lastCore])
					{
						lastCore = k;
						continue;
					}
					if (anchorChains[curr.first].nodes[curr.second] != readPaths[readi].paths[j].path[k])
					{
						lastCore = k;
						continue;
					}
				}
				if (prev.second == curr.second + 1)
				{
					fw = false;
					if (anchorChains[prev.first].nodes[prev.second] != (readPaths[readi].paths[j].path[lastCore] ^ firstBitUint64_t))
					{
						lastCore = k;
						continue;
					}
					if (anchorChains[curr.first].nodes[curr.second] != (readPaths[readi].paths[j].path[k] ^ firstBitUint64_t))
					{
						lastCore = k;
						continue;
					}
				}
				std::vector<uint64_t> allele { readPaths[readi].paths[j].path.begin() + lastCore + 1, readPaths[readi].paths[j].path.begin() + k };
				if (!fw)
				{
					assert((readPaths[readi].paths[j].path[k] & firstBitUint64_t) != (anchorChains[curr.first].nodes[curr.second] & firstBitUint64_t));
					std::reverse(allele.begin(), allele.end());
					for (size_t l = 0; l < allele.size(); l++)
					{
						allele[l] ^= firstBitUint64_t;
					}
				}
				else
				{
					assert((readPaths[readi].paths[j].path[k] & firstBitUint64_t) == (anchorChains[curr.first].nodes[curr.second] & firstBitUint64_t));
				}
				size_t index = getAlleleIndex(allelesPerChain, curr.first, std::min(curr.second, prev.second)+1, allele);
				assert(curr.first < anchorChains.size());
				assert(curr.second < anchorChains[curr.first].nodeOffsets.size());
				int diagonal = (int)readPos - (int)anchorChains[curr.first].nodeOffsets[curr.second];
				if (!fw) diagonal = (int)readPos + (int)anchorChains[curr.first].nodeOffsets[curr.second];
				std::cerr << "insert allele read " << readi << " chain " << (curr.first & maskUint64_t) << (fw ? "+" : "-") << " bubble " << (std::min(curr.second, prev.second)+1) << " allele " << index << " diagonal " << diagonal << std::endl;
				assert(diagonal > -(int)(anchorChains[curr.first].nodeOffsets.back()*2+readPos*2));
				assert(diagonal < (int)(anchorChains[curr.first].nodeOffsets.back()*2+readPos*2));
				allelesInThisRead.emplace_back(curr.first + (fw ? firstBitUint64_t : 0), diagonal, std::min(curr.second, prev.second)+1, index);
				lastCore = k;
			}
		}
		// chain connection alleles
		for (size_t j = 1; j < chainPositionsInReads[readi].size(); j++)
		{
			std::pair<size_t, bool> lastChain = std::make_pair(chainPositionsInReads[readi][j-1].chain & maskUint64_t, chainPositionsInReads[readi][j-1].chain & firstBitUint64_t);
			std::pair<size_t, bool> currChain = std::make_pair(chainPositionsInReads[readi][j].chain & maskUint64_t, chainPositionsInReads[readi][j].chain & firstBitUint64_t);
			if (!validChainEdges.hasEdge(lastChain, currChain)) continue;
			bool prevFw = chainPositionsInReads[readi][j-1].chain & firstBitUint64_t;
			bool currFw = chainPositionsInReads[readi][j].chain & firstBitUint64_t;
			std::vector<uint64_t> prevAllele;
			prevAllele.emplace_back(chainPositionsInReads[readi][j].chain);
			size_t prevIndex = getAlleleIndex(allelesPerChain, chainPositionsInReads[readi][j-1].chain & maskUint64_t, prevFw ? anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodes.size() : 0, prevAllele);
			std::vector<uint64_t> currAllele;
			currAllele.emplace_back(chainPositionsInReads[readi][j-1].chain ^ firstBitUint64_t);
			size_t currIndex = getAlleleIndex(allelesPerChain, chainPositionsInReads[readi][j].chain & maskUint64_t, currFw ? 0 : anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.size(), currAllele);
			int prevDiagonal, currDiagonal;
			if (prevFw)
			{
				prevDiagonal = chainPositionsInReads[readi][j-1].chainEndPosInRead - anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodeOffsets.back();
			}
			else
			{
				prevDiagonal = chainPositionsInReads[readi][j-1].chainEndPosInRead;
			}
			if (currFw)
			{
				currDiagonal = chainPositionsInReads[readi][j].chainStartPosInRead + unitigGraph.lengths[anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes[0]];
			}
			else
			{
				currDiagonal = chainPositionsInReads[readi][j].chainStartPosInRead + anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodeOffsets.back() + unitigGraph.lengths[anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.back()];
			}
			std::cerr << "tangleconnection insert allele read " << readi << " chain " << (chainPositionsInReads[readi][j-1].chain & maskUint64_t) << (prevFw ? "+" : "-") << " bubble " << (prevFw ? anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodes.size() : 0) << " allele " << prevIndex << " diagonal " << prevDiagonal << std::endl;
			std::cerr << "tangleconnection insert allele read " << readi << " chain " << (chainPositionsInReads[readi][j].chain & maskUint64_t) << (currFw ? "+" : "-") << " bubble " << (currFw ? 0 : anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.size()) << " allele " << currIndex << " diagonal " << currDiagonal << std::endl;
			allelesInThisRead.emplace_back(chainPositionsInReads[readi][j-1].chain, prevDiagonal, prevFw ? anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodes.size() : 0, prevIndex);
			allelesInThisRead.emplace_back(chainPositionsInReads[readi][j].chain, currDiagonal, currFw ? 0 : anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.size(), currIndex);
		}
		// simple bubble alleles with gaps
		uint64_t lastAnchor = std::numeric_limits<size_t>::max();
		int lastAnchorDiagonal = 0;
		size_t lastAnchorJ = 0;
		size_t lastSimpleBubble = std::numeric_limits<size_t>::max();
		std::tuple<size_t, size_t, size_t> uniqueLastSimpleBubbleAllele { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
		for (size_t j = 0; j < readPaths[readi].paths.size(); j++)
		{
			size_t readPos = readPaths[readi].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[readi].paths[j].path.size(); k++)
			{
				const uint64_t node = readPaths[readi].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[readi].paths[j].pathLeftClipKmers;
				if (coreNodeLocator.count(node & maskUint64_t) == 0)
				{
					if (simpleBubbleAllele.count(node & maskUint64_t) == 1 && node != lastSimpleBubble)
					{
						if (std::get<0>(uniqueLastSimpleBubbleAllele) != std::numeric_limits<size_t>::max())
						{
							uniqueLastSimpleBubbleAllele = std::make_tuple(std::numeric_limits<size_t>::max()-1, std::numeric_limits<size_t>::max()-1, std::numeric_limits<size_t>::max()-1);
						}
						else
						{
							uniqueLastSimpleBubbleAllele = simpleBubbleAllele.at(node & maskUint64_t);
						}
						lastSimpleBubble = node;
					}
					continue;
				}
				std::pair<size_t, size_t> pos = coreNodeLocator.at(node & maskUint64_t);
				assert((node & maskUint64_t) == (anchorChains[pos.first].nodes[pos.second] & maskUint64_t));
				bool fw = (node & firstBitUint64_t) == (anchorChains[pos.first].nodes[pos.second] & firstBitUint64_t);
				int diagonal = (int)readPos - (int)anchorChains[pos.first].nodeOffsets[pos.second];
				if (!fw) diagonal = (int)readPos + (int)anchorChains[pos.first].nodeOffsets[pos.second];
				if (lastAnchor == std::numeric_limits<size_t>::max() || lastAnchorJ == j || std::get<0>(uniqueLastSimpleBubbleAllele) != pos.first || (fw && std::get<1>(uniqueLastSimpleBubbleAllele)+1 != pos.second) || (!fw && std::get<1>(uniqueLastSimpleBubbleAllele) != pos.second))
				{
					lastAnchor = node;
					lastAnchorJ = j;
					uniqueLastSimpleBubbleAllele = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
					lastAnchorDiagonal = diagonal;
					lastSimpleBubble = std::numeric_limits<size_t>::max();
					continue;
				}
				std::pair<size_t, size_t> prevpos = coreNodeLocator.at(lastAnchor & maskUint64_t);
				bool prevfw = (lastAnchor & firstBitUint64_t) == (anchorChains[prevpos.first].nodes[prevpos.second] & firstBitUint64_t);
				int maxDiagonalDifference = (int)(anchorChains[pos.first].nodeOffsets.back() + unitigGraph.lengths[anchorChains[pos.first].nodes.back() & maskUint64_t])*0.5;
				if (fw != prevfw || prevpos.first != pos.first || (fw && prevpos.second+1 != pos.second) || (!fw && prevpos.second != pos.second+1) || diagonal > lastAnchorDiagonal + maxDiagonalDifference || diagonal < lastAnchorDiagonal - maxDiagonalDifference)
				{
					lastAnchor = node;
					lastAnchorJ = j;
					uniqueLastSimpleBubbleAllele = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
					lastAnchorDiagonal = diagonal;
					lastSimpleBubble = std::numeric_limits<size_t>::max();
					continue;
				}
				assert(std::get<2>(uniqueLastSimpleBubbleAllele) < std::numeric_limits<size_t>::max()-1);
				std::cerr << "simplebubble insert allele read " << readi << " chain " << (pos.first) << (fw ? "+" : "-") << " bubble " << (std::min(pos.second, prevpos.second)+1) << " allele " << std::get<2>(uniqueLastSimpleBubbleAllele) << " diagonal " << diagonal << std::endl;
				allelesInThisRead.emplace_back(pos.first + (fw ? firstBitUint64_t : 0), diagonal, std::min(pos.second, prevpos.second)+1, std::get<2>(uniqueLastSimpleBubbleAllele));
				lastAnchor = node;
				lastAnchorJ = j;
				uniqueLastSimpleBubbleAllele = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
				lastAnchorDiagonal = diagonal;
				lastSimpleBubble = std::numeric_limits<size_t>::max();
			}
		}
		std::sort(allelesInThisRead.begin(), allelesInThisRead.end());
		for (size_t i = 0; i < allelesInThisRead.size(); i++)
		{
			if (i == 0 || std::get<0>(allelesInThisRead[i]) != std::get<0>(allelesInThisRead[i-1]) || std::get<1>(allelesInThisRead[i]) > std::get<1>(allelesInThisRead[i-1]) + (int)(anchorChains[std::get<0>(allelesInThisRead[i])].nodeOffsets.back() + unitigGraph.lengths[anchorChains[std::get<0>(allelesInThisRead[i])].nodes.back() & maskUint64_t])/2)
			{
				allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].emplace_back();
				std::get<0>(allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back()) = readi;
				std::get<1>(allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back()) = std::get<0>(allelesInThisRead[i]) & firstBitUint64_t;
				std::get<2>(allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back()) = std::get<1>(allelesInThisRead[i]);
			}
			std::cerr << "has allele chain " << (std::get<0>(allelesInThisRead[i]) & maskUint64_t) << ((std::get<0>(allelesInThisRead[i]) & firstBitUint64_t) ? "+" : "-") << " read " << readi << " diagonal " << std::get<1>(allelesInThisRead[i]) << " bubble " << std::get<2>(allelesInThisRead[i]) << " allele " << std::get<3>(allelesInThisRead[i]) << std::endl;
			std::get<3>(allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back()).emplace_back(std::get<2>(allelesInThisRead[i]), std::get<3>(allelesInThisRead[i]));
		}
	}
	return std::make_pair(std::move(allelesPerChain), std::move(allelesPerReadPerChain));
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

std::pair<std::vector<std::vector<size_t>>, size_t> getTopNCoveredAlleles(const std::vector<std::vector<size_t>>& readsPerAllele, const size_t numAlleles)
{
	std::vector<std::pair<size_t, size_t>> coveragePerAllele;
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		coveragePerAllele.emplace_back(i, readsPerAllele[i].size());
	}
	std::sort(coveragePerAllele.begin(), coveragePerAllele.end(), [](auto left, auto right) { return left.second > right.second; });
	size_t numMinors = 0;
	for (size_t i = numAlleles; i < coveragePerAllele.size(); i++)
	{
		numMinors += coveragePerAllele[i].second;
	}
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < coveragePerAllele.size() && i < numAlleles; i++)
	{
		result.emplace_back(readsPerAllele[coveragePerAllele[i].first]);
	}
	return std::make_pair(result, numMinors);
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
		size_t ploidyHere = (double)coveragePerAllele[i].second / approxOneHapCoverage + 0.5;
		std::cerr << i << " " << coveragePerAllele[i].second << " " << ploidyHere << std::endl;
		estimatedTotalPloidy += ploidyHere;
	}
	std::vector<std::vector<size_t>> result;
	std::cerr << "estimated total ploidy " << estimatedTotalPloidy << std::endl;
	if (estimatedTotalPloidy < ploidy || estimatedTotalPloidy > ploidy+1) //loosely allow one more found ploidy just in case haplotype coverage hovers around 1.5x average
	{
		return std::make_pair(result, 0);
	}
	size_t minorCount = 0;
	while (coveragePerAllele.size() >= 1)
	{
		size_t lastPloidy = (double)coveragePerAllele.back().second / approxOneHapCoverage + 0.5;
		if (lastPloidy > 0) break;
		minorCount += readsPerAllele[coveragePerAllele.back().first].size();
		coveragePerAllele.pop_back();
	}
	size_t gotPloidy = 0;
	for (size_t i = 0; i < coveragePerAllele.size(); i++)
	{
		gotPloidy += (double)coveragePerAllele[i].second / approxOneHapCoverage + 0.5;
		result.emplace_back(readsPerAllele[coveragePerAllele[i].first]);
		std::sort(result.back().begin(), result.back().end());
	}
	assert(gotPloidy == estimatedTotalPloidy);
	return std::make_pair(result, minorCount);
}

bool allelesCouldSomewhatMatch(const std::vector<std::vector<size_t>>& bestAllelesLeft, const std::vector<std::vector<size_t>>& bestAllelesRight, const size_t ploidy, const size_t extraNonMatchers)
{
	assert(ploidy == 2); // needs to handle polyploidy
	if (bestAllelesLeft.size() != ploidy) return false;
	if (bestAllelesRight.size() != ploidy) return false;
	size_t requiredCoveragePerHap = 3;
	size_t normalCounts = intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[1]);
	size_t crossCounts = intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[0]);
	if (intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) >= requiredCoveragePerHap && intersectSize(bestAllelesLeft[1], bestAllelesRight[1]) >= requiredCoveragePerHap && normalCounts > (crossCounts+extraNonMatchers)*9) return true;
	if (intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) >= requiredCoveragePerHap && intersectSize(bestAllelesLeft[1], bestAllelesRight[0]) >= requiredCoveragePerHap && crossCounts > (normalCounts+extraNonMatchers)*9) return true;
	return false;
}

bool allelesCouldMatch(const std::vector<std::vector<size_t>>& bestAllelesLeft, const std::vector<std::vector<size_t>>& bestAllelesRight, const size_t ploidy, const double approxOneHapCoverage)
{
	assert(ploidy == 2); // needs to handle polyploidy
	assert(bestAllelesLeft.size() >= 2);
	assert(bestAllelesRight.size() >= 2);
	size_t requiredCoverage = 4;
	// size_t requiredCoverage = (double)approxOneHapCoverage * ((double)ploidy - 0.5);
	size_t normalCounts = intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[1]);
	size_t crossCounts = intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[0]);
	if (normalCounts >= requiredCoverage && normalCounts > crossCounts*9) return true;
	if (crossCounts >= requiredCoverage && crossCounts > normalCounts*9) return true;
	return false;
}

std::vector<bool> getVeryLikelyHaplotypeInformativeSites(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& coreNodeChain, const size_t ploidy, const double approxOneHapCoverage)
{
	assert(ploidy == 2); // needs to handle polyploidy
	std::vector<std::vector<std::vector<size_t>>> bestAlleles;
	bestAlleles.resize(readsPerAllele.size());
	std::vector<bool> notgood;
	notgood.resize(readsPerAllele.size(), false);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << std::endl;
		size_t minorCount = 0;
		std::tie(bestAlleles[i], minorCount) = getMostCoveredAlleles(readsPerAllele[i], ploidy, approxOneHapCoverage);
		if (bestAlleles[i].size() == 0)
		{
			std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << " best alleles zero, not good" << std::endl;
			notgood[i] = true;
			continue;
		}
		if (bestAlleles[i].size() < ploidy)
		{
			std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << " best alleles less than ploidy, not good" << std::endl;
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
			std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << " minor allele coverage too high, not good" << std::endl;
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

std::vector<bool> extendMaybeHaplotypeInformativeSites(const std::vector<bool>& veryLikelyHaplotypeInformative, const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& coreNodeChain, const size_t ploidy, const double approxOneHapCoverage)
{
	std::vector<bool> result = veryLikelyHaplotypeInformative;
	std::vector<std::vector<std::vector<size_t>>> bestAlleles;
	std::vector<size_t> minorCount;
	bestAlleles.resize(readsPerAllele.size());
	minorCount.resize(readsPerAllele.size(), 0);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		std::tie(bestAlleles[i], minorCount[i]) = getTopNCoveredAlleles(readsPerAllele[i], ploidy);
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		if (result[i]) continue;
		for (size_t j = 0; j < veryLikelyHaplotypeInformative.size(); j++)
		{
			if (!veryLikelyHaplotypeInformative[j]) continue;
			if (allelesCouldSomewhatMatch(bestAlleles[i], bestAlleles[j], ploidy, minorCount[i]+minorCount[j]))
			{
				std::cerr << "extended haplotype informative to bubble index " << i << "/" << readsPerAllele.size() << std::endl;
				result[i] = true;
				break;
			}
		}
	}
	return result;
}

std::vector<bool> getPossiblyHaplotypeInformativeSites(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& coreNodeChain, const size_t ploidy, const double approxOneHapCoverage)
{
	std::vector<bool> veryLikelyHaplotypeInformative = getVeryLikelyHaplotypeInformativeSites(readsPerAllele, coreNodeChain, ploidy, approxOneHapCoverage);
	std::vector<bool> maybeHaplotypeInformative = extendMaybeHaplotypeInformativeSites(veryLikelyHaplotypeInformative, readsPerAllele, coreNodeChain, ploidy, approxOneHapCoverage);
	return maybeHaplotypeInformative;
}

std::vector<std::vector<size_t>> getAllelesPerHaplotype(const std::vector<size_t>& readAssignment, const size_t ploidy, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& allelesPerRead, const std::vector<size_t>& numAllelesPerSite)
{
	assert(ploidy >= 2);
	std::vector<std::vector<std::vector<size_t>>> alleleCoveragePerHaplotype;
	alleleCoveragePerHaplotype.resize(ploidy);
	for (size_t hap = 0; hap < ploidy; hap++)
	{
		alleleCoveragePerHaplotype[hap].resize(numAllelesPerSite.size());
		for (size_t site = 0; site < numAllelesPerSite.size(); site++)
		{
			alleleCoveragePerHaplotype[hap][site].resize(numAllelesPerSite[site], 0);
		}
	}
	for (size_t i = 0; i < readAssignment.size(); i++)
	{
		if (readAssignment[i] == std::numeric_limits<size_t>::max()) continue;
		for (std::tuple<size_t, size_t, size_t> t : allelesPerRead[i])
		{
			size_t site = std::get<0>(t);
			size_t allele = std::get<1>(t);
			size_t weight = std::get<2>(t);
			alleleCoveragePerHaplotype[readAssignment[i]][site][allele] += weight;
		}
	}
	std::vector<std::vector<size_t>> result;
	result.resize(ploidy);
	for (size_t hap = 0; hap < ploidy; hap++)
	{
		result[hap].resize(numAllelesPerSite.size());
	}
	for (size_t site = 0; site < numAllelesPerSite.size(); site++)
	{
		size_t totalCoverageOfSite = 0;
		size_t totalMismatchesInSite = 0;
		for (size_t hap = 0; hap < ploidy; hap++)
		{
			std::pair<size_t, size_t> bestHere { std::numeric_limits<size_t>::max(), 0 };
			size_t totalCount = 0;
			for (size_t allele = 0; allele < alleleCoveragePerHaplotype[hap][site].size(); allele++)
			{
				if (alleleCoveragePerHaplotype[hap][site][allele] > bestHere.second)
				{
					bestHere.first = allele;
					bestHere.second = alleleCoveragePerHaplotype[hap][site][allele];
				}
				totalCount += alleleCoveragePerHaplotype[hap][site][allele];
			}
			totalMismatchesInSite += totalCount-bestHere.second;
			totalCoverageOfSite += totalCount;
			result[hap][site] = bestHere.first;
		}
		if (totalMismatchesInSite*10 > totalCoverageOfSite)
		{
			for (size_t hap = 0; hap < ploidy; hap++)
			{
				result[hap][site] = std::numeric_limits<size_t>::max();
			}
		}
	}
	return result;
}

size_t getHaplotypeScore(const std::vector<size_t>& readAssignment, const size_t ploidy, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& allelesPerRead, const std::vector<size_t>& numAllelesPerSite)
{
	assert(ploidy >= 2);
	size_t score = 0;
	std::vector<std::vector<size_t>> allelesPerHaplotype = getAllelesPerHaplotype(readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
	for (size_t i = 0; i < readAssignment.size(); i++)
	{
		for (std::tuple<size_t, size_t, size_t> t : allelesPerRead[i])
		{
			size_t site = std::get<0>(t);
			size_t allele = std::get<1>(t);
			size_t weight = std::get<2>(t);
			if (readAssignment[i] == std::numeric_limits<size_t>::max())
			{
				score += weight;
				continue;
			}
			size_t haplotypeAllele = allelesPerHaplotype[readAssignment[i]][site];
			if (haplotypeAllele == std::numeric_limits<size_t>::max())
			{
				score += weight;
				continue;
			}
			if (allele != haplotypeAllele)
			{
				score += weight*10;
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

std::vector<size_t> getUnweightedHeuristicMEC(const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& allelesPerRead, const size_t ploidy, const std::vector<size_t>& numAllelesPerSite)
{
	const size_t beamWidth = 1000;
	size_t addEverythingUntilHere = log(beamWidth)/log(ploidy+1);
	assert(ploidy >= 2);
	std::vector<ReadPartition> activeAssignments;
	activeAssignments.emplace_back();
	activeAssignments.back().maxAssignmentPlusOne = 0;
	activeAssignments.back().score = 0;
	for (size_t i = 0; i < addEverythingUntilHere && i < allelesPerRead.size(); i++)
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
		activeAssignments[i].score = getHaplotypeScore(activeAssignments[i].readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
	}
	std::sort(activeAssignments.begin(), activeAssignments.end(), [](const auto& left, const auto& right) { return left.score < right.score; });
	for (size_t i = addEverythingUntilHere; i < allelesPerRead.size(); i++)
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
				nextActiveAssignments.back().score = getHaplotypeScore(nextActiveAssignments.back().readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
			}
			nextActiveAssignments.emplace_back();
			nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
			nextActiveAssignments.back().maxAssignmentPlusOne = haplotype.maxAssignmentPlusOne;
			nextActiveAssignments.back().readAssignment.push_back(std::numeric_limits<size_t>::max());
			nextActiveAssignments.back().score = getHaplotypeScore(nextActiveAssignments.back().readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
		}
		std::sort(nextActiveAssignments.begin(), nextActiveAssignments.end(), [](const auto& left, const auto& right) { return left.score < right.score; });
		if (nextActiveAssignments.size() > beamWidth)
		{
			nextActiveAssignments.erase(nextActiveAssignments.begin()+beamWidth, nextActiveAssignments.end());
		}
		activeAssignments = nextActiveAssignments;
	}
	assert(activeAssignments[0].readAssignment.size() == allelesPerRead.size());
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
	std::vector<std::vector<size_t>> allelesPerHaplotype = getAllelesPerHaplotype(activeAssignments[0].readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
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
	return activeAssignments[0].readAssignment;
}

std::vector<size_t> getReadOrder(const std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>>& readsPerAllele, const size_t readCount)
{
	std::vector<size_t> firstPerRead;
	std::vector<size_t> lastPerRead;
	firstPerRead.resize(readCount, std::numeric_limits<size_t>::max());
	lastPerRead.resize(readCount, std::numeric_limits<size_t>::max());
	for (size_t site = 0; site < readsPerAllele.size(); site++)
	{
		for (size_t allele = 0; allele < readsPerAllele[site].size(); allele++)
		{
			for (std::pair<size_t, size_t> pair : readsPerAllele[site][allele])
			{
				size_t read = pair.first;
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
			for (std::pair<size_t, size_t> pair : readsPerAllele[i][j])
			{
				size_t read = pair.first;
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

std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> mergeReadsPerAllele(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const size_t readCount)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> allelesPerRead;
	allelesPerRead.resize(readCount);
	for (size_t site = 0; site < readsPerAllele.size(); site++)
	{
		for (size_t allele = 0; allele < readsPerAllele[site].size(); allele++)
		{
			for (size_t read : readsPerAllele[site][allele])
			{
				allelesPerRead[read].emplace_back(site, allele);
			}
		}
	}
	std::sort(allelesPerRead.begin(), allelesPerRead.end());
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> result;
	size_t currentReadCount = 1;
	for (size_t i = 1; i < allelesPerRead.size(); i++)
	{
		if (allelesPerRead[i] == allelesPerRead[i-1])
		{
			currentReadCount += 1;
			continue;
		}
		if (allelesPerRead[i-1].size() < 2)
		{
			currentReadCount = 1;
			continue;
		}
		result.emplace_back();
		for (auto pair : allelesPerRead[i-1])
		{
			result.back().emplace_back(pair.first, pair.second, currentReadCount);
		}
		currentReadCount = 1;
	}
	if (allelesPerRead.back().size() >= 2)
	{
		result.emplace_back();
		for (auto pair : allelesPerRead.back())
		{
			result.back().emplace_back(pair.first, pair.second, currentReadCount);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		std::sort(result[i].begin(), result[i].end());
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right)
	{
		if (std::get<0>(left[0]) < std::get<0>(right[0])) return true;
		if (std::get<0>(left[0]) > std::get<0>(right[0])) return false;
		if (std::get<0>(left.back()) < std::get<0>(right.back())) return true;
		if (std::get<0>(left.back()) > std::get<0>(right.back())) return false;
		if (left.size() > right.size()) return true;
		if (left.size() < right.size()) return false;
		for (size_t i = 0; i < left.size(); i++)
		{
			if (std::get<0>(left[i]) < std::get<0>(right[i])) return true;
			if (std::get<0>(left[i]) > std::get<0>(right[i])) return false;
		}
		if (std::get<2>(left[0]) > std::get<2>(right[0])) return true;
		if (std::get<2>(left[0]) < std::get<2>(right[0])) return false;
		return false;
	});
	return result;
}

std::vector<PhaseBlock> splitPhaseBlocks(const PhaseBlock& raw, const std::vector<size_t> readAssignment, const std::vector<std::vector<size_t>>& allelesPerHaplotype, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& mergedAllelesPerRead, const size_t ploidy)
{
	assert(allelesPerHaplotype.size() == ploidy);
	assert(mergedAllelesPerRead.size() == readAssignment.size());
	std::vector<size_t> previousNonGarbageSite;
	previousNonGarbageSite.resize(allelesPerHaplotype[0].size(), std::numeric_limits<size_t>::max());
	size_t lastNonGarbage = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < allelesPerHaplotype[0].size(); i++)
	{
		if (allelesPerHaplotype[0][i] == std::numeric_limits<size_t>::max()) continue;
		previousNonGarbageSite[i] = lastNonGarbage;
		lastNonGarbage = i;
	}
	std::cerr << "pre-split phased at start: " << raw.chainStartPhased << " end: " << raw.chainEndPhased << std::endl;
	std::cerr << "bubble indices:";
	for (size_t i = 0; i < raw.bubbleIndices.size(); i++) std::cerr << " " << raw.bubbleIndices[i];
	std::cerr << std::endl;
	std::cerr << "number of sites: " << raw.bubbleIndices.size() << std::endl;
	{
		std::vector<bool> notGarbageSite;
		notGarbageSite.resize(raw.bubbleIndices.size(), true);
		for (size_t i = 0; i < allelesPerHaplotype[0].size(); i++)
		{
			if (allelesPerHaplotype[0][i] == std::numeric_limits<size_t>::max())
			{
				notGarbageSite[i] = false;
			}
		}
		for (size_t i = 0; i < notGarbageSite.size(); i++)
		{
			std::cerr << (notGarbageSite[i] ? "X" : "_");
		}
		std::cerr << " valid sites" << std::endl;
	}
	for (size_t i = 0; i < mergedAllelesPerRead.size(); i++)
	{
		if (readAssignment[i] == std::numeric_limits<size_t>::max()) continue;
		std::vector<bool> hasAllele;
		hasAllele.resize(raw.bubbleIndices.size(), false);
		for (auto t : mergedAllelesPerRead[i])
		{
			hasAllele[std::get<0>(t)] = true;
		}
		for (size_t i = 0; i < hasAllele.size(); i++)
		{
			std::cerr << (hasAllele[i] ? "X" : "_");
		}
		std::cerr << " read " << i << " hap " << readAssignment[i] << std::endl;
	}
	std::vector<std::vector<bool>> haplotypeHasCrossingReads;
	haplotypeHasCrossingReads.resize(ploidy);
	for (size_t i = 0; i < ploidy; i++)
	{
		haplotypeHasCrossingReads[i].resize(raw.bubbleIndices.size(), false);
	}
	for (size_t i = 0; i < mergedAllelesPerRead.size(); i++)
	{
		if (readAssignment[i] == std::numeric_limits<size_t>::max()) continue;
		for (size_t j = 1; j < mergedAllelesPerRead[i].size(); j++)
		{
			if (allelesPerHaplotype[0][std::get<0>(mergedAllelesPerRead[i][j])] == std::numeric_limits<size_t>::max()) continue;
			for (size_t k = 0; k < j; k++)
			{
				if (allelesPerHaplotype[0][std::get<0>(mergedAllelesPerRead[i][k])] == std::numeric_limits<size_t>::max()) continue;
				if (std::get<0>(mergedAllelesPerRead[i][k]) == previousNonGarbageSite[std::get<0>(mergedAllelesPerRead[i][j])])
				{
					haplotypeHasCrossingReads[readAssignment[i]][std::get<0>(mergedAllelesPerRead[i][j])] = true;
				}
			}
		}
	}
	std::vector<PhaseBlock> result;
	result.emplace_back();
	result.back().chainNumber = raw.chainNumber;
	result.back().chainStartPhased = raw.chainStartPhased;
	result.back().allelesPerHaplotype.resize(ploidy);
	bool lastPhased = raw.chainEndPhased;
	size_t bubblei = 0;
	for (size_t i = 0; i < allelesPerHaplotype[0].size(); i++)
	{
		if (allelesPerHaplotype[0][i] == std::numeric_limits<size_t>::max())
		{
			if (i == 0) result.back().chainStartPhased = false;
			if (i == allelesPerHaplotype[0].size()-1) lastPhased = false;
			continue;
		}
		if (previousNonGarbageSite[i] == std::numeric_limits<size_t>::max())
		{
			result.back().bubbleIndices.emplace_back(raw.bubbleIndices[bubblei]);
			for (size_t k = 0; k < ploidy; k++)
			{
				assert(raw.allelesPerHaplotype[k][bubblei] != std::numeric_limits<size_t>::max());
				result.back().allelesPerHaplotype[k].push_back(raw.allelesPerHaplotype[k][bubblei]);
			}
			// can't split before first site
			bubblei += 1;
			continue;
		}
		bool validSite = true;
		for (size_t k = 0; k < ploidy; k++)
		{
			if (!haplotypeHasCrossingReads[k][i])
			{
				validSite = false;
				break;
			}
		}
		if (!validSite)
		{
			std::cerr << "split block at chain " << raw.chainNumber << " index " << i << std::endl;
			result.emplace_back();
			result.back().chainNumber = raw.chainNumber;
			result.back().allelesPerHaplotype.resize(ploidy);
		}
		result.back().bubbleIndices.emplace_back(raw.bubbleIndices[bubblei]);
		for (size_t k = 0; k < ploidy; k++)
		{
			assert(raw.allelesPerHaplotype[k][bubblei] != std::numeric_limits<size_t>::max());
			result.back().allelesPerHaplotype[k].push_back(raw.allelesPerHaplotype[k][bubblei]);
		}
		bubblei += 1;
	}
	assert(bubblei == raw.bubbleIndices.size());
	result.back().chainEndPhased = lastPhased;
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		assert(result[i].bubbleIndices.size() == result[i].allelesPerHaplotype[0].size());
		if (result[i].bubbleIndices.size() < 2)
		{
			std::cerr << "remove too small block " << std::endl;
			std::swap(result[i], result.back());
			result.pop_back();
		}
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.bubbleIndices[0] < right.bubbleIndices[0]; });
	return result;
}

std::vector<PhaseBlock> getChainPhaseBlocks(const size_t coreNodeChainIndex, const std::vector<std::tuple<size_t, bool, int, std::vector<std::pair<size_t, size_t>>>>& allelesPerRead, const std::vector<uint64_t>& coreNodeChain, const size_t ploidy, const size_t alleleCount, const double approxOneHapCoverage)
{
	assert(ploidy >= 2);
	assert(alleleCount >= 1);
	size_t readCount = allelesPerRead.size();
	std::vector<std::vector<std::vector<size_t>>> unfilteredReadsPerAllele;
	unfilteredReadsPerAllele.resize(alleleCount);
	for (size_t readi = 0; readi < allelesPerRead.size(); readi++)
	{
		for (auto pair : std::get<3>(allelesPerRead[readi]))
		{
			assert(pair.first < unfilteredReadsPerAllele.size());
			while (pair.second >= unfilteredReadsPerAllele[pair.first].size()) unfilteredReadsPerAllele[pair.first].emplace_back();
			assert(std::find(unfilteredReadsPerAllele[pair.first][pair.second].begin(), unfilteredReadsPerAllele[pair.first][pair.second].end(), readi) == unfilteredReadsPerAllele[pair.first][pair.second].end());
			unfilteredReadsPerAllele[pair.first][pair.second].emplace_back(readi);
		}
	}
	std::cerr << "check chain " << coreNodeChainIndex << std::endl;
	std::vector<bool> maybeHaplotypeInformative = getPossiblyHaplotypeInformativeSites(unfilteredReadsPerAllele, coreNodeChain, ploidy, approxOneHapCoverage);
	std::cerr << "first informative: " << (maybeHaplotypeInformative[0] ? "yes" : "no") << " last: " << (maybeHaplotypeInformative.back() ? "yes" : "no") << std::endl;
	size_t countMaybeInformativeSites = 0;
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		if (maybeHaplotypeInformative[i]) countMaybeInformativeSites += 1;
	}
	if (countMaybeInformativeSites < 2)
	{
		std::cerr << "num informative sites " << countMaybeInformativeSites << ", skipped" << std::endl;
		return std::vector<PhaseBlock>{};
	}
	std::vector<std::vector<std::vector<size_t>>> readsPerAllele = getInformativeReads(unfilteredReadsPerAllele, maybeHaplotypeInformative);
	size_t informativeReadCount = 0;
	{
		phmap::flat_hash_set<size_t> foundRead;
		for (size_t i = 0; i < readsPerAllele.size(); i++)
		{
			for (size_t j = 0; j < readsPerAllele[i].size(); j++)
			{
				for (size_t read : readsPerAllele[i][j])
				{
					foundRead.insert(read);
				}
			}
		}
		informativeReadCount = foundRead.size();
	}
	std::cerr << informativeReadCount << " informative real reads" << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> mergedAllelesPerRead = mergeReadsPerAllele(readsPerAllele, readCount);
	std::cerr << mergedAllelesPerRead.size() << " merged reads" << std::endl;
	std::vector<size_t> numAllelesPerSite;
	numAllelesPerSite.resize(readsPerAllele.size());
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		numAllelesPerSite[i] = readsPerAllele[i].size();
	}
	std::vector<size_t> readAssignments = getUnweightedHeuristicMEC(mergedAllelesPerRead, ploidy, numAllelesPerSite);
	std::vector<std::vector<size_t>> allelesPerHaplotype = getAllelesPerHaplotype(readAssignments, ploidy, mergedAllelesPerRead, numAllelesPerSite);
	PhaseBlock rawResult;
	rawResult.chainStartPhased = false;
	rawResult.chainEndPhased = false;
	if (maybeHaplotypeInformative[0]) rawResult.chainStartPhased = true;
	if (maybeHaplotypeInformative.back()) rawResult.chainEndPhased = true;
	rawResult.chainNumber = coreNodeChainIndex;
	for (size_t i = 0; i < ploidy; i++)
	{
		rawResult.allelesPerHaplotype.emplace_back();
	}
	// todo split across haplo block boundaries
	size_t bubblei = 0;
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		if (!maybeHaplotypeInformative[i]) continue;
		if (allelesPerHaplotype[0][bubblei] == std::numeric_limits<size_t>::max())
		{
			bubblei += 1;
			continue;
		}
		rawResult.bubbleIndices.push_back(i);
		for (size_t j = 0; j < ploidy; j++)
		{
			rawResult.allelesPerHaplotype[j].push_back(allelesPerHaplotype[j][bubblei]);
		}
		bubblei += 1;
	}
	assert(bubblei == allelesPerHaplotype[0].size());
	return splitPhaseBlocks(rawResult, readAssignments, allelesPerHaplotype, mergedAllelesPerRead, ploidy);
}

std::vector<PhaseBlock> phaseCoreChains(const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<std::tuple<size_t, bool, int, std::vector<std::pair<size_t, size_t>>>>>& allelesPerReadPerChain, const double approxOneHapCoverage)
{
	assert(anchorChains.size() == allelesPerReadPerChain.size());
	std::vector<PhaseBlock> result;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::cerr << "check chain " << i << " start " << ((anchorChains[i].nodes[0] & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes[0] & maskUint64_t) << " end " << ((anchorChains[i].nodes.back() & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes.back() & maskUint64_t) << " ploidy " << anchorChains[i].ploidy << std::endl;
		if (anchorChains[i].ploidy == 1)
		{
			continue;
		}
		else if (anchorChains[i].ploidy == 2)
		{
			auto phaseBlocks = getChainPhaseBlocks(i, allelesPerReadPerChain[i], anchorChains[i].nodes, anchorChains[i].ploidy, anchorChains[i].nodes.size()+1, approxOneHapCoverage);
			if (phaseBlocks.size() == 0) continue;
			std::cerr << "phased with " << phaseBlocks.size() << " blocks. start: " << (phaseBlocks[0].chainStartPhased ? "yes" : "no") << ", end: " << (phaseBlocks.back().chainEndPhased ? "yes" : "no") << std::endl;
			result.insert(result.end(), phaseBlocks.begin(), phaseBlocks.end());
			for (size_t j = 0; j < phaseBlocks.size(); j++)
			{
				for (size_t sitei = 0; sitei < phaseBlocks[j].allelesPerHaplotype[0].size(); sitei++)
				{
					phmap::flat_hash_set<size_t> foundAlleles;
					for (size_t k = 0; k < anchorChains[i].ploidy; k++)
					{
						foundAlleles.insert(phaseBlocks[j].allelesPerHaplotype[k][sitei]);
					}
					if (foundAlleles.size() != anchorChains[i].ploidy)
					{
						std::cerr << "red flag: allele shared between haps, ploidy " << anchorChains[i].ploidy << ", num distinct alleles: " << foundAlleles.size() << std::endl;
					}
				}
			}
		}
		else
		{
			std::cerr << "skipped chain with ploidy " << anchorChains[i].ploidy << std::endl;
		}
	}
	return result;
}

void markNodePositions(std::vector<std::pair<size_t, size_t>>& result, const uint64_t startNode, const uint64_t endNode, const SparseEdgeContainer& edges, const size_t chain, const size_t offset)
{
	std::vector<uint64_t> stack;
	stack.push_back(startNode);
	phmap::flat_hash_set<uint64_t> visited;
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (visited.count(top & maskUint64_t) == 1) continue;
		visited.insert(top & maskUint64_t);
		assert(result[top & maskUint64_t].first == std::numeric_limits<size_t>::max() || (top & maskUint64_t) == (startNode & maskUint64_t));
		result[top & maskUint64_t] = std::make_pair(chain, offset);
		if (top != (startNode ^ firstBitUint64_t) && top != endNode)
		{
			for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
			{
				stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
			}
		}
		if (top != startNode && top != (endNode ^ firstBitUint64_t))
		{
			for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, (top & firstBitUint64_t) ^ firstBitUint64_t)))
			{
				stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
			}
		}
	}
}

bool canUnambiguouslyMark(const SparseEdgeContainer& edges, const uint64_t startNode, const uint64_t endNode, const std::vector<bool>& anchor)
{
	std::vector<uint64_t> stack;
	phmap::flat_hash_set<uint64_t> visited;
	size_t anchorCount = 0;
	stack.push_back(startNode);
	while (stack.size() >= 1)
	{
		uint64_t top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		visited.insert(top);
		if (anchor[top & maskUint64_t])
		{
			anchorCount += 1;
			if (anchorCount > 2) return false;
		}
		if (!anchor[top & maskUint64_t]) stack.push_back(top ^ firstBitUint64_t);
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	if (visited.count(endNode) == 0 && visited.count(endNode ^ firstBitUint64_t) == 0) return false;
	if (anchorCount != 2) return false;
	return true;
}

std::vector<std::pair<size_t, size_t>> getNodeLocationsWithinChain(const UnitigGraph& unitigGraph, const std::vector<AnchorChain>& anchorChains)
{
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<std::pair<size_t, size_t>> result;
	result.resize(unitigGraph.nodeCount(), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			anchor[anchorChains[i].nodes[j] & maskUint64_t] = true;
		}
	}
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 1; j < anchorChains[i].nodes.size(); j++)
		{
			if (canUnambiguouslyMark(edges, anchorChains[i].nodes[j-1], anchorChains[i].nodes[j], anchor))
			{
				markNodePositions(result, anchorChains[i].nodes[j-1], anchorChains[i].nodes[j], edges, i, j-1);
			}
			result[anchorChains[i].nodes[j-1] & maskUint64_t] = std::make_pair(i, j-1);
			result[anchorChains[i].nodes[j] & maskUint64_t] = std::make_pair(i, j);
		}
		result[anchorChains[i].nodes.back() & maskUint64_t] = std::make_pair(i, anchorChains[i].nodes.size()-1);
	}
	return result;
}

void unzipPhaseBlocks(UnitigGraph& resultGraph, std::vector<ReadPathBundle>& resultPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<PhaseBlock>& chainHaplotypes, const std::vector<std::vector<std::vector<std::tuple<size_t, bool, int>>>>& readsPerHaplotypePerChain, const std::vector<std::pair<size_t, size_t>>& nodeLocationInChain, const phmap::flat_hash_set<std::pair<size_t, size_t>>& solvedLocations, const phmap::flat_hash_set<std::pair<size_t, size_t>>& solvedAnchors, const phmap::flat_hash_set<size_t>& solvedInterchainTangles, const VectorWithDirection<size_t>& nodeLocationInInterchainTangles, const size_t tangleCount, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>>& haplotypeDiagonalsPerRead, const phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>>& validHaplotypeConnectionsPerChainEdge)
{
	std::vector<int> maxDiagonalDifferences;
	maxDiagonalDifferences.resize(anchorChains.size());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		maxDiagonalDifferences[i] = ((int)(anchorChains[i].nodeOffsets.back() + resultGraph.lengths[anchorChains[i].nodes.back() & maskUint64_t]) * .5);
	}
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> nodeReplacement;
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> interchainTangleConnectionHap;
	interchainTangleConnectionHap.resize(tangleCount);
	for (size_t i = 0; i < resultPaths.size(); i++)
	{
		// std::cerr << "new read index " << i << std::endl;
		// std::cerr << "diagonals:";
		// for (auto t : haplotypeDiagonalsPerRead[i]) std::cerr << " " << std::get<0>(t) << "(" << chainHaplotypes[std::get<0>(t)].chainNumber << ")" << (std::get<1>(t) ? "+" : "-") << " " << std::get<2>(t) << " " << std::get<3>(t);
		std::cerr << std::endl;
		std::vector<std::tuple<size_t, size_t, size_t, size_t, int, uint64_t>> readToAnchorMatches; // chain, bubble, readPos (extrapolated node end), haplotype, diagonal, node
		for (size_t j = 0; j < resultPaths[i].paths.size(); j++)
		{
			size_t readPos = resultPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < resultPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = resultPaths[i].paths[j].path[k];
				readPos += resultGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= resultPaths[i].paths[j].pathLeftClipKmers;
				assert((node & maskUint64_t) < nodeLocationInChain.size());
				size_t chain = nodeLocationInChain[node & maskUint64_t].first;
				size_t offset = nodeLocationInChain[node & maskUint64_t].second;
				if (chain == std::numeric_limits<size_t>::max()) continue;
				if ((node & maskUint64_t) != (anchorChains[chain].nodes[offset] & maskUint64_t))
				{
					continue;
				}
				bool fw = (node & firstBitUint64_t) == (anchorChains[chain].nodes[offset] & firstBitUint64_t);
				int diagonal = (int)readPos - (int)anchorChains[chain].nodeOffsets[offset];
				if (!fw) diagonal = (int)readPos + (int)anchorChains[chain].nodeOffsets[offset];
				size_t uniqueHap = std::numeric_limits<size_t>::max();
				if (anchorChains[chain].ploidy == 1)
				{
					uniqueHap = 0;
				}
				else
				{
					for (auto t : haplotypeDiagonalsPerRead[i])
					{
						if (chainHaplotypes[std::get<0>(t)].chainNumber != chain) continue;
						if (std::get<1>(t) != fw) continue;
						if (std::get<2>(t) > diagonal + (int)maxDiagonalDifferences[chain]) continue;
						if (std::get<2>(t) < diagonal - (int)maxDiagonalDifferences[chain]) continue;
						if (uniqueHap == std::numeric_limits<size_t>::max())
						{
							uniqueHap = std::get<3>(t);
						}
						else if (uniqueHap != std::get<3>(t))
						{
							uniqueHap = std::numeric_limits<size_t>::max()-1;
						}
					}
				}
				if (uniqueHap == std::numeric_limits<size_t>::max()-1) uniqueHap = std::numeric_limits<size_t>::max();
				readToAnchorMatches.emplace_back(j, k, readPos, uniqueHap, diagonal, node);
			}
		}
		bool previousAnchorSolved = false;
		std::vector<std::pair<size_t, size_t>> breakBeforeHere;
		for (size_t m = 1; m < readToAnchorMatches.size(); m++)
		{
			uint64_t lastNode = std::get<5>(readToAnchorMatches[m-1]);
			size_t lastChain = nodeLocationInChain[lastNode & maskUint64_t].first;
			size_t lastOffset = nodeLocationInChain[lastNode & maskUint64_t].second;
			assert((anchorChains[lastChain].nodes[lastOffset] & maskUint64_t) == (lastNode & maskUint64_t));
			bool lastFw = (lastNode & firstBitUint64_t) == (anchorChains[lastChain].nodes[lastOffset] & firstBitUint64_t);
			uint64_t currNode = std::get<5>(readToAnchorMatches[m]);
			size_t currChain = nodeLocationInChain[currNode & maskUint64_t].first;
			size_t currOffset = nodeLocationInChain[currNode & maskUint64_t].second;
			assert((anchorChains[currChain].nodes[currOffset] & maskUint64_t) == (currNode & maskUint64_t));
			bool currFw = (currNode & firstBitUint64_t) == (anchorChains[currChain].nodes[currOffset] & firstBitUint64_t);
			bool solveThisBlock = false;
			bool solveLastAnchor = false;
			bool solveCurrAnchor = false;
			size_t solutionHap = std::numeric_limits<size_t>::max();
			// std::cerr << "check " << m << " " << lastChain << " " << lastOffset << " " << (lastFw ? "fw" : "bw") << " " << (lastNode & maskUint64_t) << " " << std::get<3>(readToAnchorMatches[m-1]) << " " << currChain << " " << currOffset << " " << (currFw ? "fw" : "bw") << " " << (currNode & maskUint64_t) << " " << std::get<3>(readToAnchorMatches[m]) << std::endl;
			if (currChain == lastChain && currFw == lastFw && ((currFw && currOffset == lastOffset+1) || (!currFw && currOffset+1 == lastOffset) || (currOffset == lastOffset)))
			{
				// std::cerr << "has chain match" << std::endl;
				int lastDiagonal = std::get<4>(readToAnchorMatches[m-1]);
				int currDiagonal = std::get<4>(readToAnchorMatches[m]);
				// std::cerr << "diagonals " << lastDiagonal << " " << currDiagonal << std::endl;
				if (lastDiagonal > currDiagonal - maxDiagonalDifferences[currChain] && lastDiagonal < currDiagonal + maxDiagonalDifferences[currChain])
				{
					// std::cerr << "has curr last diagonal" << std::endl;
					if (std::get<3>(readToAnchorMatches[m-1]) == std::get<3>(readToAnchorMatches[m]) && std::get<3>(readToAnchorMatches[m]) != std::numeric_limits<size_t>::max())
					{
						// std::cerr << "has unique hap " << std::get<3>(readToAnchorMatches[m]) << " check solved location " << currChain << " " << std::min(lastOffset, currOffset) << std::endl;
						if (solvedLocations.count(std::make_pair(currChain, std::min(lastOffset, currOffset))) == 1)
						{
							// std::cerr << "has solved location" << std::endl;
							solveThisBlock = true;
							solveLastAnchor = (solvedAnchors.count(std::make_pair(lastChain, lastOffset)) == 1);
							if (previousAnchorSolved) solveLastAnchor = false;
							solveCurrAnchor = (solvedAnchors.count(std::make_pair(currChain, currOffset)) == 1);
							solutionHap = std::get<3>(readToAnchorMatches[m]);
						}
					}
				}
			}
			if (!solveThisBlock)
			{
				if (nodeLocationInInterchainTangles[std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t)] == nodeLocationInInterchainTangles[std::make_pair(currNode & maskUint64_t, (currNode ^ firstBitUint64_t) & firstBitUint64_t)])
				{
					if (nodeLocationInInterchainTangles[std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t)] != std::numeric_limits<size_t>::max())
					{
						size_t lastChain = nodeLocationInChain[lastNode & maskUint64_t].first;
						size_t lastOffset = nodeLocationInChain[lastNode & maskUint64_t].second;
						assert((anchorChains[lastChain].nodes[lastOffset] & maskUint64_t) == (lastNode & maskUint64_t));
						bool lastFw = (lastNode & firstBitUint64_t) == (anchorChains[lastChain].nodes[lastOffset] & firstBitUint64_t);
						size_t lastHap = std::get<3>(readToAnchorMatches[m-1]);
						uint64_t currNode = std::get<5>(readToAnchorMatches[m]);
						size_t currChain = nodeLocationInChain[currNode & maskUint64_t].first;
						size_t currOffset = nodeLocationInChain[currNode & maskUint64_t].second;
						assert((anchorChains[currChain].nodes[currOffset] & maskUint64_t) == (currNode & maskUint64_t));
						bool currFw = (currNode & firstBitUint64_t) == (anchorChains[currChain].nodes[currOffset] & firstBitUint64_t);
						size_t currHap = std::get<3>(readToAnchorMatches[m]);
						// std::cerr << "has interchain nodes" << std::endl;
						bool shouldBreakThis = false;
						if (validChainEdges.hasEdge(std::make_pair(lastChain, lastFw), std::make_pair(currChain, currFw)))
						{
							// std::cerr << "has interchain edge" << std::endl;
							if (solvedInterchainTangles.count(nodeLocationInInterchainTangles[std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t)]) == 1)
							{
								shouldBreakThis = true;
								// std::cerr << "is in solved tangle" << std::endl;
								if (validHaplotypeConnectionsPerChainEdge.count(lastChain + (lastFw ? firstBitUint64_t : 0)) == 1 && validHaplotypeConnectionsPerChainEdge.at(lastChain + (lastFw ? firstBitUint64_t : 0)).count(currChain + (currFw ? firstBitUint64_t : 0)) == 1)
								{
									// std::cerr << "has some valid haplotypes" << std::endl;
									if (validHaplotypeConnectionsPerChainEdge.at(lastChain + (lastFw ? firstBitUint64_t : 0)).at(currChain + (currFw ? firstBitUint64_t : 0)).count(std::make_pair(lastHap, currHap)) == 1)
									{
										// std::cerr << "has valid haplotype" << std::endl;
										shouldBreakThis = false;
										solveThisBlock = true;
										solveLastAnchor = false;
										solveCurrAnchor = false;
										if (anchorChains[lastChain].nodes.size() == 1)
										{
											solveLastAnchor = (solvedAnchors.count(std::make_pair(lastChain, lastOffset)) == 1);
											if (previousAnchorSolved) solveLastAnchor = false;
										}
										if (anchorChains[currChain].nodes.size() == 1)
										{
											solveCurrAnchor = (solvedAnchors.count(std::make_pair(currChain, currOffset)) == 1);
										}
										auto keypair = canon(std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t), std::make_pair(currNode & maskUint64_t, currNode & firstBitUint64_t));
										std::pair<uint64_t, uint64_t> key { keypair.first.first + (keypair.first.second ? firstBitUint64_t : 0), keypair.second.first + (keypair.second.second ? firstBitUint64_t : 0) };
										size_t tangle = nodeLocationInInterchainTangles[std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t)];
										assert(tangle < interchainTangleConnectionHap.size());
										if (interchainTangleConnectionHap[tangle].count(key) == 0)
										{
											solutionHap = interchainTangleConnectionHap[tangle].size();
											interchainTangleConnectionHap[tangle][key] = solutionHap;
										}
										else
										{
											solutionHap = interchainTangleConnectionHap[tangle].at(key);
										}
									}
								}
							}
						}
						if (shouldBreakThis && std::get<0>(readToAnchorMatches[m-1]) == std::get<0>(readToAnchorMatches[m]) && std::get<1>(readToAnchorMatches[m-1])+1 == std::get<1>(readToAnchorMatches[m]))
						{
							breakBeforeHere.emplace_back(std::get<0>(readToAnchorMatches[m]), std::get<1>(readToAnchorMatches[m]));
						}
					}
				}
			}
			if (solveThisBlock)
			{
				// std::cerr << "read " << i << " solve from " << std::get<0>(readToAnchorMatches[m-1]) << "," << std::get<1>(readToAnchorMatches[m-1]) << " to " << std::get<0>(readToAnchorMatches[m]) << "," << std::get<1>(readToAnchorMatches[m]) << std::endl;
				// if (solveLastAnchor) std::cerr << "read " << i << " solve anchor " << std::get<0>(readToAnchorMatches[m-1]) << " " << std::get<1>(readToAnchorMatches[m-1]) << " " << (std::get<5>(readToAnchorMatches[m-1]) & maskUint64_t) << std::endl;
				// if (solveCurrAnchor) std::cerr << "read " << i << " solve anchor " << std::get<0>(readToAnchorMatches[m]) << " " << std::get<1>(readToAnchorMatches[m]) << " " << (std::get<5>(readToAnchorMatches[m]) & maskUint64_t) << std::endl;
				for (size_t j = std::get<0>(readToAnchorMatches[m-1]); j <= std::get<0>(readToAnchorMatches[m]); j++)
				{
					for (size_t k = (j == std::get<0>(readToAnchorMatches[m-1]) ? std::get<1>(readToAnchorMatches[m-1]) : 0); k < (j == std::get<0>(readToAnchorMatches[m]) ? std::get<1>(readToAnchorMatches[m])+1 : resultPaths[i].paths[j].path.size()); k++)
					{
						if (!solveLastAnchor && j == std::get<0>(readToAnchorMatches[m-1]) && k == std::get<1>(readToAnchorMatches[m-1])) continue;
						if (!solveCurrAnchor && j == std::get<0>(readToAnchorMatches[m]) && k == std::get<1>(readToAnchorMatches[m])) continue;
						uint64_t node = resultPaths[i].paths[j].path[k];
						assert((node & maskUint64_t) < nodeLocationInChain.size());
						if (nodeReplacement.count(std::make_pair(node & maskUint64_t, solutionHap)) == 0)
						{
							nodeReplacement[std::make_pair(node & maskUint64_t, solutionHap)] = resultGraph.lengths.size();
							resultGraph.lengths.push_back(resultGraph.lengths[node & maskUint64_t]);
							resultGraph.coverages.emplace_back();
						}
						assert(nodeReplacement.count(std::make_pair(node & maskUint64_t, solutionHap)) == 1);
						// std::cerr << "read " << i << " " << j << " " << k << " replace " << (node & maskUint64_t) << " with " << nodeReplacement.at(std::make_pair(node & maskUint64_t, solutionHap)) << std::endl;
						resultPaths[i].paths[j].path[k] = nodeReplacement.at(std::make_pair(node & maskUint64_t, solutionHap)) + (node & firstBitUint64_t);
					}
				}
			}
			previousAnchorSolved = solveCurrAnchor;
		}
		if (breakBeforeHere.size() > 0)
		{
			for (size_t j = breakBeforeHere.size()-1; j < breakBeforeHere.size(); j--)
			{
				// std::cerr << "break read " << i << " index " << breakBeforeHere[j].first << " " << breakBeforeHere[j].second << " (" << (resultPaths[i].paths[breakBeforeHere[j].first].path[breakBeforeHere[j].second-1] & maskUint64_t) << " " << (resultPaths[i].paths[breakBeforeHere[j].first].path[breakBeforeHere[j].second] & maskUint64_t) << ")" << std::endl;
				assert(breakBeforeHere[j].first < resultPaths[i].paths.size());
				assert(breakBeforeHere[j].second < resultPaths[i].paths[breakBeforeHere[j].first].path.size());
				assert(breakBeforeHere[j].second > 0);
				resultPaths[i].paths.emplace_back();
				resultPaths[i].paths.back().pathLeftClipKmers = 0;
				resultPaths[i].paths.back().pathRightClipKmers = resultPaths[i].paths[breakBeforeHere[j].first].pathRightClipKmers;
				size_t startPos = resultPaths[i].paths[breakBeforeHere[j].first].readStartPos;
				for (size_t k = 0; k < breakBeforeHere[j].second; k++)
				{
					startPos += resultGraph.lengths[resultPaths[i].paths[breakBeforeHere[j].first].path[k] & maskUint64_t];
					if (k == 0) startPos -= resultPaths[i].paths[breakBeforeHere[j].first].pathLeftClipKmers;
				}
				resultPaths[i].paths.back().readStartPos = startPos;
				resultPaths[i].paths.back().path.insert(resultPaths[i].paths.back().path.end(), resultPaths[i].paths[breakBeforeHere[j].first].path.begin() + breakBeforeHere[j].second, resultPaths[i].paths[breakBeforeHere[j].first].path.end());
				resultPaths[i].paths[breakBeforeHere[j].first].path.erase(resultPaths[i].paths[breakBeforeHere[j].first].path.begin() + breakBeforeHere[j].second, resultPaths[i].paths[breakBeforeHere[j].first].path.end());
				resultPaths[i].paths[breakBeforeHere[j].first].pathRightClipKmers = 0;
			}
			std::sort(resultPaths[i].paths.begin(), resultPaths[i].paths.end(), [](const auto& left, const auto& right) { return left.readStartPos < right.readStartPos; });
		}
	}
}

phmap::flat_hash_set<uint64_t> getInterchainTangleNodes(const uint64_t start, const std::vector<bool>& anchor, const SparseEdgeContainer& edges)
{
	assert(anchor[start & maskUint64_t]);
	phmap::flat_hash_set<uint64_t> result;
	std::vector<uint64_t> stack;
	stack.push_back(start);
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (result.count(top) == 1) continue;
		result.insert(top);
		if (!anchor[top & maskUint64_t]) stack.push_back(top ^ firstBitUint64_t);
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	return result;
}

VectorWithDirection<size_t> getNodeLocationsWithinInterchainTangles(const UnitigGraph& unitigGraph, const std::vector<AnchorChain>& anchorChains)
{
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			anchor[anchorChains[i].nodes[j] & maskUint64_t] = true;
		}
	}
	VectorWithDirection<size_t> result;
	result.resize(unitigGraph.nodeCount(), std::numeric_limits<size_t>::max());
	phmap::flat_hash_set<uint64_t> checked;
	size_t tangleNumber = 0;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (checked.count(anchorChains[i].nodes[0] ^ firstBitUint64_t) == 0)
		{
			phmap::flat_hash_set<uint64_t> nodesHere = getInterchainTangleNodes(anchorChains[i].nodes[0] ^ firstBitUint64_t, anchor, edges);
			assert(nodesHere.size() >= 1);
			assert(nodesHere.count(anchorChains[i].nodes[0] ^ firstBitUint64_t) == 1);
			for (auto node : nodesHere)
			{
				checked.insert(node);
				assert(result[std::make_pair(node & maskUint64_t, node & firstBitUint64_t)] == std::numeric_limits<size_t>::max());
				result[std::make_pair(node & maskUint64_t, node & firstBitUint64_t)] = tangleNumber;
			}
			tangleNumber += 1;
		}
		if (checked.count(anchorChains[i].nodes.back()) == 0)
		{
			phmap::flat_hash_set<uint64_t> nodesHere = getInterchainTangleNodes(anchorChains[i].nodes.back(), anchor, edges);
			assert(nodesHere.size() >= 1);
			assert(nodesHere.count(anchorChains[i].nodes.back()) == 1);
			for (auto node : nodesHere)
			{
				checked.insert(node);
				assert(result[std::make_pair(node & maskUint64_t, node & firstBitUint64_t)] == std::numeric_limits<size_t>::max());
				result[std::make_pair(node & maskUint64_t, node & firstBitUint64_t)] = tangleNumber;
			}
			tangleNumber += 1;
		}
	}
	return result;
}

std::pair<std::vector<bool>, std::vector<bool>> getUnzippableChainEnds(const std::vector<AnchorChain>& anchorChains, const std::vector<PhaseBlock>& chainHaplotypes, const SparseEdgeContainer& validChainEdges, const phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>>& validHaplotypeConnectionsPerChainEdge)
{
	std::vector<bool> canUnzipStart;
	std::vector<bool> canUnzipEnd;
	canUnzipStart.resize(anchorChains.size(), false);
	canUnzipEnd.resize(anchorChains.size(), false);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy == 1)
		{
			canUnzipStart[i] = true;
			canUnzipEnd[i] = true;
		}
	}
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		if (chainHaplotypes[i].chainStartPhased) canUnzipStart[chainHaplotypes[i].chainNumber] = true;
		if (chainHaplotypes[i].chainEndPhased) canUnzipEnd[chainHaplotypes[i].chainNumber] = true;
	}
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		for (auto edge : validChainEdges.getEdges(std::make_pair(i, true)))
		{
			if (validHaplotypeConnectionsPerChainEdge.count(i + firstBitUint64_t) == 0)
			{
				canUnzipEnd[i] = false;
				break;
			}
			if (validHaplotypeConnectionsPerChainEdge.at(i + firstBitUint64_t).count(edge.first + (edge.second ? firstBitUint64_t : 0)) == 0)
			{
				canUnzipEnd[i] = false;
				break;
			}
			if (validHaplotypeConnectionsPerChainEdge.at(i + firstBitUint64_t).at(edge.first + (edge.second ? firstBitUint64_t : 0)).size() == 0)
			{
				canUnzipEnd[i] = false;
				break;
			}
		}
		for (auto edge : validChainEdges.getEdges(std::make_pair(i, false)))
		{
			if (validHaplotypeConnectionsPerChainEdge.count(i) == 0)
			{
				canUnzipStart[i] = false;
				break;
			}
			if (validHaplotypeConnectionsPerChainEdge.at(i).count(edge.first + (edge.second ? firstBitUint64_t : 0)) == 0)
			{
				canUnzipStart[i] = false;
				break;
			}
			if (validHaplotypeConnectionsPerChainEdge.at(i).at(edge.first + (edge.second ? firstBitUint64_t : 0)).size() == 0)
			{
				canUnzipStart[i] = false;
				break;
			}
		}
	}
	// could do smartly, but naive method is good enough
	while (true)
	{
		bool repeat = false;
		for (size_t i = 0; i < anchorChains.size(); i++)
		{
			if (!canUnzipStart[i])
			{
				for (auto edge : validChainEdges.getEdges(std::make_pair(i, false)))
				{
					if (edge.second)
					{
						if (canUnzipStart[edge.first]) repeat = true;
						canUnzipStart[edge.first] = false;
					}
					else
					{
						if (canUnzipEnd[edge.first]) repeat = true;
						canUnzipEnd[edge.first] = false;
					}
				}
			}
			if (!canUnzipEnd[i])
			{
				for (auto edge : validChainEdges.getEdges(std::make_pair(i, true)))
				{
					if (edge.second)
					{
						if (canUnzipStart[edge.first]) repeat = true;
						canUnzipStart[edge.first] = false;
					}
					else
					{
						if (canUnzipEnd[edge.first]) repeat = true;
						canUnzipEnd[edge.first] = false;
					}
				}
			}
		}
		if (!repeat) break;
	}
	return std::make_pair(canUnzipStart, canUnzipEnd);
}

phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> getValidHaplotypeConnectionsPerChainEdge(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>>& haplotypeDiagonalsPerRead, const std::vector<PhaseBlock>& chainHaplotypes, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
	std::vector<int> maxDiagonalDifferences;
	maxDiagonalDifferences.resize(anchorChains.size());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		maxDiagonalDifferences[i] = ((int)(anchorChains[i].nodeOffsets.back() + unitigGraph.lengths[anchorChains[i].nodes.back() & maskUint64_t]) * .5);
	}
	assert(chainPositionsInReads.size() == readPaths.size());
	assert(readPaths.size() == haplotypeDiagonalsPerRead.size());
	std::vector<std::vector<size_t>> chainStartHaplotypeInRead;
	std::vector<std::vector<size_t>> chainEndHaplotypeInRead;
	chainStartHaplotypeInRead.resize(chainPositionsInReads.size());
	chainEndHaplotypeInRead.resize(chainPositionsInReads.size());
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		chainStartHaplotypeInRead[i].resize(chainPositionsInReads[i].size(), std::numeric_limits<size_t>::max());
		chainEndHaplotypeInRead[i].resize(chainPositionsInReads[i].size(), std::numeric_limits<size_t>::max());
		for (size_t j = 0; j < chainPositionsInReads[i].size(); j++)
		{
			size_t uniqueStartHap = std::numeric_limits<size_t>::max();
			size_t uniqueEndHap = std::numeric_limits<size_t>::max();
			size_t chain = chainPositionsInReads[i][j].chain & maskUint64_t;
			bool fw = (chainPositionsInReads[i][j].chain & firstBitUint64_t);
			int startDiagonal = chainPositionsInReads[i][j].chainStartPosInRead + (int)unitigGraph.lengths[anchorChains[chain].nodes[0]];
			int endDiagonal = chainPositionsInReads[i][j].chainEndPosInRead - (int)anchorChains[chain].nodeOffsets.back();
			if (!fw)
			{
				startDiagonal = chainPositionsInReads[i][j].chainStartPosInRead + (int)unitigGraph.lengths[anchorChains[chain].nodes.back()] + (int)anchorChains[chain].nodeOffsets.back();
				endDiagonal = chainPositionsInReads[i][j].chainEndPosInRead;
			}
			for (size_t k = 0; k < haplotypeDiagonalsPerRead[i].size(); k++)
			{
				if (chainHaplotypes[std::get<0>(haplotypeDiagonalsPerRead[i][k])].chainNumber != chain) continue;
				if (std::get<1>(haplotypeDiagonalsPerRead[i][k]) != fw) continue;
				if (std::get<2>(haplotypeDiagonalsPerRead[i][k]) > startDiagonal - maxDiagonalDifferences[chain] && std::get<2>(haplotypeDiagonalsPerRead[i][k]) < startDiagonal + maxDiagonalDifferences[chain])
				{
					if (uniqueStartHap == std::numeric_limits<size_t>::max())
					{
						uniqueStartHap = std::get<3>(haplotypeDiagonalsPerRead[i][k]);
					}
					else if (uniqueStartHap != std::get<3>(haplotypeDiagonalsPerRead[i][k]))
					{
						uniqueStartHap = std::numeric_limits<size_t>::max()-1;
					}
				}
				if (std::get<2>(haplotypeDiagonalsPerRead[i][k]) > endDiagonal - maxDiagonalDifferences[chain] && std::get<2>(haplotypeDiagonalsPerRead[i][k]) < endDiagonal + maxDiagonalDifferences[chain])
				{
					if (uniqueEndHap == std::numeric_limits<size_t>::max())
					{
						uniqueEndHap = std::get<3>(haplotypeDiagonalsPerRead[i][k]);
					}
					else if (uniqueEndHap != std::get<3>(haplotypeDiagonalsPerRead[i][k]))
					{
						uniqueEndHap = std::numeric_limits<size_t>::max()-1;
					}
				}
			}
			if (anchorChains[chain].ploidy == 1)
			{
				assert(uniqueStartHap == std::numeric_limits<size_t>::max());
				assert(uniqueEndHap == std::numeric_limits<size_t>::max());
				uniqueStartHap = 0;
				uniqueEndHap = 0;
			}
			if (uniqueStartHap < std::numeric_limits<size_t>::max()-1) chainStartHaplotypeInRead[i][j] = uniqueStartHap;
			if (uniqueEndHap < std::numeric_limits<size_t>::max()-1) chainEndHaplotypeInRead[i][j] = uniqueEndHap;
		}
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>>> edgeHapCoverage;
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		for (size_t j = 1; j < chainPositionsInReads[i].size(); j++)
		{
			uint64_t prevChain = chainPositionsInReads[i][j-1].chain;
			uint64_t currChain = chainPositionsInReads[i][j].chain;
			size_t prevHap = chainEndHaplotypeInRead[i][j-1];
			size_t currHap = chainStartHaplotypeInRead[i][j];
			if (prevHap == std::numeric_limits<size_t>::max()) continue;
			if (currHap == std::numeric_limits<size_t>::max()) continue;
			if (!validChainEdges.hasEdge(std::make_pair(prevChain & maskUint64_t, prevChain & firstBitUint64_t), std::make_pair(currChain & maskUint64_t, currChain & firstBitUint64_t))) continue;
			edgeHapCoverage[prevChain][currChain][std::make_pair(prevHap, currHap)] += 1;
			edgeHapCoverage[currChain ^ firstBitUint64_t][prevChain ^firstBitUint64_t][std::make_pair(currHap, prevHap)] += 1;
		}
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> result;
	for (const auto& pair : edgeHapCoverage)
	{
		for (const auto& pair2 : pair.second)
		{
			for (const auto& pair3 : pair2.second)
			{
				std::cerr << "hap check chain " << (pair.first & maskUint64_t) << ((pair.first & firstBitUint64_t) ? "+" : "-") << " to chain " << (pair2.first & maskUint64_t) << ((pair2.first & firstBitUint64_t) ? "+" : "-") << " haps " << pair3.first.first << " " << pair3.first.second << " coverage " << pair3.second << std::endl;
				if (pair3.second < approxOneHapCoverage * 0.5) continue;
				result[pair.first][pair2.first].insert(pair3.first);
			}
		}
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipPhaseBlocks(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<std::tuple<size_t, bool, int, std::vector<std::pair<size_t, size_t>>>>>& allelesPerReadPerChain, const std::vector<PhaseBlock>& chainHaplotypes, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	std::vector<std::pair<size_t, size_t>> nodeLocationInChain = getNodeLocationsWithinChain(unitigGraph, anchorChains);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			assert(nodeLocationInChain[anchorChains[i].nodes[j] & maskUint64_t] == std::make_pair(i, j));
		}
	}
	VectorWithDirection<size_t> nodeLocationInInterchainTangles = getNodeLocationsWithinInterchainTangles(unitigGraph, anchorChains);
	size_t tangleCount = 0;
	for (size_t i = 0; i < nodeLocationInInterchainTangles.size(); i++)
	{
		if (nodeLocationInInterchainTangles[std::make_pair(i, true)] != std::numeric_limits<size_t>::max())
		{
			tangleCount = std::max(tangleCount, nodeLocationInInterchainTangles[std::make_pair(i, true)]);
		}
		if (nodeLocationInInterchainTangles[std::make_pair(i, false)] != std::numeric_limits<size_t>::max())
		{
			tangleCount = std::max(tangleCount, nodeLocationInInterchainTangles[std::make_pair(i, false)]);
		}
	}
	tangleCount += 1;
	std::vector<std::vector<std::vector<std::tuple<size_t, bool, int>>>> readsPerHaplotypePerChain;
	readsPerHaplotypePerChain.resize(chainHaplotypes.size());
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		const size_t chain = chainHaplotypes[i].chainNumber;
		const size_t ploidy = anchorChains[chain].ploidy;
		readsPerHaplotypePerChain[i].resize(ploidy);
		phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> alleleBelongsToUniqueHaplotype;
		for (size_t j = 0; j < chainHaplotypes[i].bubbleIndices.size(); j++)
		{
			for (size_t k = 0; k < ploidy; k++)
			{
				auto key = std::make_pair(chainHaplotypes[i].bubbleIndices[j], chainHaplotypes[i].allelesPerHaplotype[k][j]);
				if (alleleBelongsToUniqueHaplotype.count(key) == 1)
				{
					alleleBelongsToUniqueHaplotype[key] = std::numeric_limits<size_t>::max();
				}
				else
				{
					alleleBelongsToUniqueHaplotype[key] = k;
				}
				std::cerr << "block " << i << " chain " << chain << " bubble " << chainHaplotypes[i].bubbleIndices[j] << " allele " << chainHaplotypes[i].allelesPerHaplotype[k][j] << " belongs to hap " << k << std::endl;
			}
		}
		for (size_t j = 0; j < allelesPerReadPerChain[chain].size(); j++)
		{
			std::vector<size_t> countMatches;
			countMatches.resize(ploidy);
			size_t totalMatches = 0;
			for (auto pair : std::get<3>(allelesPerReadPerChain[chain][j]))
			{
				if (alleleBelongsToUniqueHaplotype.count(pair) == 0) continue;
				if (alleleBelongsToUniqueHaplotype.at(pair) == std::numeric_limits<size_t>::max()) continue;
				countMatches[alleleBelongsToUniqueHaplotype.at(pair)] += 1;
				totalMatches += 1;
			}
			if (totalMatches == 0)
			{
				std::cerr << "read " << std::get<0>(allelesPerReadPerChain[chain][j]) << " diagonal " << std::get<2>(allelesPerReadPerChain[chain][j]) << " has no het matches in block " << i << " chain " << chain << std::endl;
				continue;
			}
			size_t bestHap = 0;
			for (size_t k = 0; k < ploidy; k++)
			{
				if (countMatches[k] > countMatches[bestHap]) bestHap = k;
			}
			if (countMatches[bestHap] <= totalMatches * 0.9)
			{
				std::cerr << "read " << std::get<0>(allelesPerReadPerChain[chain][j]) << " diagonal " << std::get<2>(allelesPerReadPerChain[chain][j]) << " is ambiguous (" << countMatches[bestHap] << " vs " << totalMatches << ") in block " << i << " chain " << chain << std::endl;
				continue;
			}
			readsPerHaplotypePerChain[i][bestHap].emplace_back(std::get<0>(allelesPerReadPerChain[chain][j]), std::get<1>(allelesPerReadPerChain[chain][j]), std::get<2>(allelesPerReadPerChain[chain][j]));
			std::cerr << "read " << std::get<0>(allelesPerReadPerChain[chain][j]) << " diagonal " << std::get<2>(allelesPerReadPerChain[chain][j]) << " has hap " << bestHap << " in block " << i << " chain " << chain << std::endl;
		}
		for (size_t j = 0; j < ploidy; j++)
		{
			std::sort(readsPerHaplotypePerChain[i][j].begin(), readsPerHaplotypePerChain[i][j].end());
			assert(readsPerHaplotypePerChain[i][j].size() >= 1);
		}
	}
	std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>> haplotypeDiagonalsPerRead;
	haplotypeDiagonalsPerRead.resize(resultPaths.size());
	for (size_t block = 0; block < readsPerHaplotypePerChain.size(); block++)
	{
		for (size_t hap = 0; hap < readsPerHaplotypePerChain[block].size(); hap++)
		{
			for (auto t : readsPerHaplotypePerChain[block][hap])
			{
				size_t read = std::get<0>(t);
				bool fw = std::get<1>(t);
				int diagonal = std::get<2>(t);
				haplotypeDiagonalsPerRead[read].emplace_back(block, fw, diagonal, hap);
				std::cerr << "readhap " << read << " " << block << " (" << (chainHaplotypes[block].chainNumber & maskUint64_t) << ")" << (fw ? "+" : "-") << " " << diagonal << " " << hap << std::endl;
			}
		}
	}
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedLocations;
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedAnchors;
	std::vector<bool> canUnzipStart;
	std::vector<bool> canUnzipEnd;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> validHaplotypeConnectionsPerChainEdge = getValidHaplotypeConnectionsPerChainEdge(unitigGraph, readPaths, anchorChains, haplotypeDiagonalsPerRead, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	std::tie(canUnzipStart, canUnzipEnd) = getUnzippableChainEnds(anchorChains, chainHaplotypes, validChainEdges, validHaplotypeConnectionsPerChainEdge);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::cerr << "anchor chain " << i << " (" << (anchorChains[i].nodes[0] & maskUint64_t) << " to " << (anchorChains[i].nodes.back() & maskUint64_t) << ") unzippable at start: " << (canUnzipStart[i] ? "yes" : "no") << " at end: " << (canUnzipEnd[i] ? "yes" : "no") << std::endl;
	}
	phmap::flat_hash_set<size_t> solvedInterchainTangles;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		// canUnzipStart[i] = false;
		// canUnzipEnd[i] = false;
		if (canUnzipStart[i])
		{
			std::pair<size_t, bool> key { anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] ^ firstBitUint64_t) & firstBitUint64_t };
			assert(nodeLocationInInterchainTangles[key] != std::numeric_limits<size_t>::max());
			solvedInterchainTangles.insert(nodeLocationInInterchainTangles[key]);
		}
		if (canUnzipEnd[i])
		{
			std::pair<size_t, bool> key { anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t };
			assert(nodeLocationInInterchainTangles[key] != std::numeric_limits<size_t>::max());
			solvedInterchainTangles.insert(nodeLocationInInterchainTangles[key]);
		}
	}
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		const size_t chain = chainHaplotypes[i].chainNumber;
		size_t minSolvedIndex = chainHaplotypes[i].bubbleIndices[0];
		size_t maxSolvedIndex = chainHaplotypes[i].bubbleIndices.back();
		assert(chainHaplotypes[i].bubbleIndices.size() >= 2);
		std::cerr << "block " << i << " chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << " (" << anchorChains[chain].nodes.size() << ")" << std::endl;
		if (!canUnzipStart[chain] && minSolvedIndex == 0) minSolvedIndex = chainHaplotypes[i].bubbleIndices[1];
		if (!canUnzipEnd[chain] && maxSolvedIndex == anchorChains[chain].nodes.size()) maxSolvedIndex = chainHaplotypes[i].bubbleIndices[chainHaplotypes[i].bubbleIndices.size()-2];
		std::cerr << "block " << i << " chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << std::endl;
		for (size_t j = minSolvedIndex; j <= maxSolvedIndex; j++)
		{
			assert(solvedLocations.count(std::make_pair(chain, j)) == 0);
			std::cerr << "solve " << chain << " " << j << " (" << (j > 0 ? (anchorChains[chain].nodes[j-1] & maskUint64_t) : -1) << " to " << (j < anchorChains[chain].nodes.size() ? (anchorChains[chain].nodes[j] & maskUint64_t) : std::numeric_limits<size_t>::max()) << ")" << std::endl;
			solvedLocations.emplace(chain, j);
			if (j != maxSolvedIndex) solvedAnchors.emplace(chain, j);
		}
		if (minSolvedIndex > 0) solvedLocations.emplace(chain, minSolvedIndex-1);
	}
	unzipPhaseBlocks(resultGraph, resultPaths, anchorChains, chainHaplotypes, readsPerHaplotypePerChain, nodeLocationInChain, solvedLocations, solvedAnchors, solvedInterchainTangles, nodeLocationInInterchainTangles, tangleCount, validChainEdges, haplotypeDiagonalsPerRead, validHaplotypeConnectionsPerChainEdge);
	RankBitvector kept;
	kept.resize(resultGraph.nodeCount());
	for (size_t i = 0; i < kept.size(); i++)
	{
		kept.set(i, true);
		if (i >= nodeLocationInChain.size()) continue;
		if (nodeLocationInChain[i].first != std::numeric_limits<size_t>::max())
		{
			assert(nodeLocationInChain[i].first < anchorChains.size());
			assert(nodeLocationInChain[i].second < anchorChains[nodeLocationInChain[i].first].nodes.size());
			if ((anchorChains[nodeLocationInChain[i].first].nodes[nodeLocationInChain[i].second] & maskUint64_t) == i)
			{
				if (solvedAnchors.count(nodeLocationInChain[i]) == 0) continue;
			}
			else
			{
				if (solvedLocations.count(nodeLocationInChain[i]) == 0) continue;
			}
			std::cerr << "don't keep " << i << ", inside phased chain" << std::endl;
			kept.set(i, false);
		}
		if (nodeLocationInInterchainTangles[std::make_pair(i, true)] != std::numeric_limits<size_t>::max() && nodeLocationInInterchainTangles[std::make_pair(i, true)] == nodeLocationInInterchainTangles[std::make_pair(i, false)] && solvedInterchainTangles.count(nodeLocationInInterchainTangles[std::make_pair(i, true)]) == 1)
		{
			std::cerr << "don't keep " << i << ", inside phased tangle " << nodeLocationInInterchainTangles[std::make_pair(i, true)] << std::endl;
			kept.set(i, false);
		}
	}
	return filterUnitigGraph(resultGraph, resultPaths, kept);
}

SparseEdgeContainer getValidChainEdges(const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> edgeCoverage;
	size_t maxChain = 0;
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		if (chainPositionsInReads[i].size() == 0) continue;
		std::cerr << i << std::endl;
		maxChain = std::max(maxChain, chainPositionsInReads[i][0].chain & maskUint64_t);
		for (size_t j = 1; j < chainPositionsInReads[i].size(); j++)
		{
			maxChain = std::max(maxChain, chainPositionsInReads[i][j].chain & maskUint64_t);
			const std::pair<size_t, bool> prevChain = std::make_pair(chainPositionsInReads[i][j-1].chain & maskUint64_t, chainPositionsInReads[i][j-1].chain & firstBitUint64_t);
			const std::pair<size_t, bool> currChain = std::make_pair(chainPositionsInReads[i][j].chain & maskUint64_t, chainPositionsInReads[i][j].chain & firstBitUint64_t);
			auto key = canon(prevChain, currChain);
			edgeCoverage[key] += 1;
		}
	}
	SparseEdgeContainer result;
	result.resize(maxChain+1);
	for (auto pair : edgeCoverage)
	{
		if (pair.second < approxOneHapCoverage * 0.5) continue;
		std::cerr << "chain edge between chains " << (pair.first.first.first) << " and " << (pair.first.second.first) << std::endl;
		result.addEdge(pair.first.first, pair.first.second);
		result.addEdge(reverse(pair.first.second), reverse(pair.first.first));
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLinearizable(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		for (auto chain : chainPositionsInReads[i])
		{
			std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
		}
	}
	SparseEdgeContainer validChainEdges = getValidChainEdges(chainPositionsInReads, approxOneHapCoverage);
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<std::tuple<size_t, bool, int, std::vector<std::pair<size_t, size_t>>>>> allelesPerReadPerChain;
	std::tie(allelesPerChain, allelesPerReadPerChain) = getRawReadInfoPerChain(anchorChains, readPaths, unitigGraph, validChainEdges, chainPositionsInReads);
	std::vector<PhaseBlock> chainHaplotypes = phaseCoreChains(anchorChains, allelesPerReadPerChain, approxOneHapCoverage);
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipPhaseBlocks(unitigGraph, readPaths, anchorChains, allelesPerReadPerChain, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	return std::make_pair(std::move(unzippedGraph), std::move(unzippedReadPaths));
}
