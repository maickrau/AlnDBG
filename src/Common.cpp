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

size_t intersectSize(const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right)
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

std::pair<uint64_t, uint64_t> canonNodePair(const uint64_t from, const uint64_t to)
{
	if ((to & maskUint64_t) < (from & maskUint64_t))
	{
		return std::make_pair(to ^ firstBitUint64_t, from ^ firstBitUint64_t);
	}
	if ((to & maskUint64_t) == (from & maskUint64_t) && (to & firstBitUint64_t) == 0 && (from & firstBitUint64_t) == 0)
	{
		return std::make_pair(to ^ firstBitUint64_t, from ^ firstBitUint64_t);
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

bool NonexistantChunk(const uint64_t chunk)
{
	if (chunk == std::numeric_limits<uint64_t>::max()) return true;
	if (chunk == std::numeric_limits<uint64_t>::max()-1) return true;
	if (chunk == (std::numeric_limits<uint64_t>::max() ^ firstBitUint64_t)) return true;
	if (chunk == (std::numeric_limits<uint64_t>::max() ^ firstBitUint64_t)-1) return true;
	return false;
}

std::chrono::time_point<std::chrono::steady_clock> getTime()
{
	return std::chrono::steady_clock::now();
}

std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
	size_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	return std::to_string(milliseconds / 1000) + "," + std::to_string(milliseconds % 1000) + " s";
}

std::chrono::time_point<std::chrono::steady_clock> programStartTime = getTime();
