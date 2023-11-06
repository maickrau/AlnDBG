#include <iostream>
#include "GraphCleaner.h"
#include "Common.h"
#include "MBGCommon.h"
#include "SparseEdgeContainer.h"

bool canReach(const std::pair<size_t, bool> start, const std::pair<size_t, bool> end, const SparseEdgeContainer& edges, const std::vector<bool>& keptNodes)
{
	phmap::flat_hash_set<uint64_t> reached;
	std::vector<uint64_t> checkStack;
	checkStack.push_back(start.first + (start.second ? firstBitUint64_t : 0));
	for (size_t i = 0; i < 5; i++)
	{
		std::vector<uint64_t> nextStack;
		for (auto pos : checkStack)
		{
			if (reached.count(pos) == 1) continue;
			reached.insert(pos);
			std::pair<size_t, bool> thispos { pos & maskUint64_t, pos & firstBitUint64_t };
			if (thispos == end) return true;
			for (auto edge : edges.getEdges(thispos))
			{
				if (!keptNodes[edge.first]) continue;
				nextStack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
			}
		}
		checkStack = nextStack;
	}
	return false;
}

bool canTopologicallyRemove(const SparseEdgeContainer& edges, const size_t node, const std::vector<bool>& keptNodes)
{
	if (edges.getEdges(std::make_pair(node, true)).size() == 0) return true;
	if (edges.getEdges(std::make_pair(node, false)).size() == 0) return true;
	if (edges.getEdges(std::make_pair(node, true)).size() > 1) return false;
	if (edges.getEdges(std::make_pair(node, false)).size() > 1) return false;
	std::pair<size_t, bool> prev = reverse(edges.getEdges(std::make_pair(node, false))[0]);
	std::pair<size_t, bool> next = edges.getEdges(std::make_pair(node, true))[0];
	if (!keptNodes[prev.first]) return false;
	if (!keptNodes[next.first]) return false;
	if (!canReach(prev, next, edges, keptNodes)) return false;
	return true;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> cleanUnitigGraph(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double averageHaplotypeCoverage)
{
	std::vector<bool> hasFwCoverage;
	std::vector<bool> hasBwCoverage;
	std::vector<bool> certainlyKept;
	hasFwCoverage.resize(unitigGraph.nodeCount(), false);
	hasBwCoverage.resize(unitigGraph.nodeCount(), false);
	certainlyKept.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				if (readPaths[i].paths[j].path[k] & firstBitUint64_t)
				{
					hasFwCoverage[readPaths[i].paths[j].path[k] & maskUint64_t] = true;
				}
				else
				{
					hasBwCoverage[readPaths[i].paths[j].path[k] & maskUint64_t] = true;
				}
			}
		}
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (!hasBwCoverage[i]) continue;
		if (!hasFwCoverage[i]) continue;
		if (unitigGraph.coverages[i] < averageHaplotypeCoverage * 0.75) continue;
		certainlyKept[i] = true;
	}
	RankBitvector keep;
	keep.resize(unitigGraph.nodeCount());
	bool removeAny = false;
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		keep.set(i, true);
		if (!hasBwCoverage[i] || !hasFwCoverage[i])
		{
			if ((double)unitigGraph.coverages[i] < averageHaplotypeCoverage * 0.5)
			{
				if (canTopologicallyRemove(edges, i, certainlyKept))
				{
					keep.set(i, false);
					removeAny = true;
				}
			}
		}
		else if ((double)unitigGraph.coverages[i] < averageHaplotypeCoverage * 0.25)
		{
			if (canTopologicallyRemove(edges, i, certainlyKept))
			{
				keep.set(i, false);
				removeAny = true;
			}
		}
	}
	if (!removeAny)
	{
		return std::make_pair(unitigGraph, readPaths);
	}
	keep.buildRanks();
	size_t countNewNodes = keep.getRank(keep.size()-1) + (keep.get(keep.size()-1) ? 1 : 0);
	return filterUnitigGraph(unitigGraph, readPaths, keep);
}
