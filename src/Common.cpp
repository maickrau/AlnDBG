#include "MBGCommon.h"
#include "Common.h"

SparseEdgeContainer getActiveEdges(const MostlySparse2DHashmap<uint8_t, size_t>& edges, const size_t numNodes)
{
	SparseEdgeContainer result;
	result.resize(numNodes);
	for (size_t i = 0; i < numNodes; i++)
	{
		for (auto pair : edges.getValues(std::make_pair(i, true)))
		{
			result.addEdge(std::make_pair(i, true), pair.first);
			result.addEdge(reverse(pair.first), std::make_pair(i, false));
		}
		for (auto pair : edges.getValues(std::make_pair(i, false)))
		{
			result.addEdge(std::make_pair(i, false), pair.first);
			result.addEdge(reverse(pair.first), std::make_pair(i, true));
		}
	}
	return result;
}

size_t intersectSize(const std::vector<size_t>& left, const std::vector<size_t>& right)
{
	size_t lefti = 0;
	size_t righti = 0;
	size_t count = 0;
	while (lefti < left.size() && righti < right.size())
	{
		if (left[lefti] < right[righti])
		{
			lefti += 1;
		}
		else if (right[righti] < left[lefti])
		{
			righti += 1;
		}
		else
		{
			assert(left[lefti] == right[righti]);
			lefti += 1;
			righti += 1;
			count += 1;
		}
	}
	return count;
}
