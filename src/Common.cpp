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

std::pair<size_t, bool> reverse(std::pair<size_t, bool> pos)
{
	return std::make_pair(pos.first, !pos.second);
}

std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> canon(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	if (to.first < from.first)
	{
		return std::make_pair(reverse(to), reverse(from));
	}
	if (to.first == from.first && !to.second && !from.second)
	{
		return std::make_pair(reverse(to), reverse(from));
	}
	return std::make_pair(from, to);
}

std::string revCompRaw(const std::string& raw)
{
	std::string result { raw.rbegin(), raw.rend() };
	for (size_t i = 0; i < result.size(); i++)
	{
		switch (result[i])
		{
			case 'a':
			case 'A':
				result[i] = 'T';
				break;
			case 'c':
			case 'C':
				result[i] = 'G';
				break;
			case 'g':
			case 'G':
				result[i] = 'C';
				break;
			case 't':
			case 'T':
				result[i] = 'A';
				break;
			default:
				assert(false);
				break;
		}
	}
	return result;
}
