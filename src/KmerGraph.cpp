#include <iostream>
#include <cassert>
#include "MBGCommon.h"
#include "KmerGraph.h"
#include "Common.h"

size_t popcount(uint64_t x);

std::pair<bool, bool> extendBreakpointsFwFw(const std::vector<size_t>& readLengths, std::vector<RankBitvector>& breakpoints, size_t leftRead, size_t leftStart, size_t leftEnd, size_t rightRead, size_t rightStart, size_t rightEnd)
{
	assert(leftRead != rightRead);
	assert(leftEnd - leftStart == rightEnd - rightStart);
	std::vector<uint64_t>& leftBits = breakpoints[leftRead].getBits();
	std::vector<uint64_t>& rightBits = breakpoints[rightRead].getBits();
	size_t pos = 0;
	bool addedAnyLeft = false;
	bool addedAnyRight = false;
	while (pos <= leftEnd - leftStart)
	{
		size_t leftIndex = (leftStart + pos) / 64;
		size_t leftOffset = (leftStart + pos) % 64;
		size_t rightIndex = (rightStart + pos) / 64;
		size_t rightOffset = (rightStart + pos) % 64;
		size_t availableLeft = 64 - leftOffset;
		size_t availableRight = 64 - rightOffset;
		size_t take = std::min(availableLeft, availableRight);
		take = std::min(take, (leftEnd-leftStart) - pos + 1);
		assert(take <= 64);
		assert(take >= 1);
		uint64_t leftGot = leftBits[leftIndex] >> leftOffset;
		uint64_t rightGot = rightBits[rightIndex] >> rightOffset;
		uint64_t mask = 0xFFFFFFFFFFFFFFFFull;
		if (take < 64) mask = (1ull << (take)) - 1;
		assert(popcount(mask) == take);
		assert((mask & 1) == (1));
		leftGot &= mask;
		rightGot &= mask;
		if ((leftGot ^ rightGot) == 0)
		{
			pos += take;
			continue;
		}
		if ((leftGot & rightGot) != rightGot)
		{
			addedAnyLeft = true;
		}
		if ((leftGot & rightGot) != leftGot)
		{
			addedAnyRight = true;
		}
		uint64_t addThese = leftGot | rightGot;
		leftBits[leftIndex] |= (addThese << leftOffset);
		rightBits[rightIndex] |= (addThese << rightOffset);
		pos += take;
	}
	assert(pos == leftEnd - leftStart + 1);
	return std::make_pair(addedAnyLeft, addedAnyRight);
}

std::pair<bool, bool> extendBreakpointsFwBw(const std::vector<size_t>& readLengths, std::vector<RankBitvector>& breakpoints, size_t leftRead, size_t leftStart, size_t leftEnd, size_t rightRead, size_t rightStart, size_t rightEnd)
{
	assert(leftRead != rightRead);
	assert(leftEnd - leftStart == rightEnd - rightStart);
	bool addedAnyLeft = false;
	bool addedAnyRight = false;
	size_t leftIndex = leftStart;
	size_t rightIndex = readLengths[rightRead] - rightStart;
	for (size_t i = 0; i <= leftEnd-leftStart; i++)
	{
		bool leftbit = breakpoints[leftRead].get(leftIndex);
		bool rightbit = breakpoints[rightRead].get(rightIndex);
		if (leftbit && !rightbit)
		{
			breakpoints[rightRead].set(rightIndex, true);
			addedAnyRight = true;
		}
		else if (!leftbit && rightbit)
		{
			breakpoints[leftRead].set(leftIndex, true);
			addedAnyLeft = true;
		}
		leftIndex += 1;
		rightIndex -= 1;
	}
	return std::make_pair(addedAnyLeft, addedAnyRight);
}

std::vector<RankBitvector> extendBreakpoints(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches)
{
	std::vector<RankBitvector> breakpoints;
	breakpoints.resize(readLengths.size());
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		breakpoints[i].resize(readLengths[i]+1);
		breakpoints[i].set(0, true);
		breakpoints[i].set(readLengths[i], true);
	}
	for (size_t groupi = 0; groupi < matches.size(); groupi++)
	{
		for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
		{
			breakpoints[matches[groupi].leftRead].set(matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, true);
			breakpoints[matches[groupi].leftRead].set(matches[groupi].leftStart + matches[groupi].matches[posi].leftStart + matches[groupi].matches[posi].length, true);
			if (matches[groupi].rightFw)
			{
				breakpoints[matches[groupi].rightRead].set(matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, true);
				breakpoints[matches[groupi].rightRead].set(matches[groupi].rightStart + matches[groupi].matches[posi].rightStart + matches[groupi].matches[posi].length, true);
			}
			else
			{
				breakpoints[matches[groupi].rightRead].set(readLengths[matches[groupi].rightRead] - (matches[groupi].rightStart + matches[groupi].matches[posi].rightStart), true);
				breakpoints[matches[groupi].rightRead].set(readLengths[matches[groupi].rightRead] - (matches[groupi].rightStart + matches[groupi].matches[posi].rightStart + matches[groupi].matches[posi].length), true);
			}
		}
	}
	while (true)
	{
		bool changed = false;
		for (size_t groupi = 0; groupi < matches.size(); groupi++)
		{
			for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
			{
				std::pair<bool, bool> addedAny;
				if (matches[groupi].rightFw)
				{
					addedAny = extendBreakpointsFwFw(readLengths, breakpoints, matches[groupi].leftRead, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
				}
				else
				{
					addedAny = extendBreakpointsFwBw(readLengths, breakpoints, matches[groupi].leftRead, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
				}
				if (addedAny.first || addedAny.second) changed = true;
			}
		}
		if (!changed) break;
	}
	for (size_t i = 0; i < breakpoints.size(); i++) breakpoints[i].buildRanks();
	return breakpoints;
}

uint64_t find(std::vector<uint64_t>& result, const size_t pos)
{
	assert((pos & firstBitUint64_t) == 0);
	if (result[pos] == (pos & maskUint64_t))
	{
		assert(result[pos] == pos);
		return pos;
	}
	uint64_t foundPos = result[pos];
	while (result[foundPos & maskUint64_t] != (foundPos & maskUint64_t))
	{
		uint64_t nextFoundPos = result[foundPos & maskUint64_t];
		if (foundPos & firstBitUint64_t) nextFoundPos ^= firstBitUint64_t;
		foundPos = nextFoundPos;
	}
	assert(result[foundPos & maskUint64_t] == (foundPos & maskUint64_t));
	uint64_t finalPos = foundPos;
	foundPos = result[pos];
	result[pos] = finalPos;
	while (result[foundPos & maskUint64_t] != (foundPos & maskUint64_t))
	{
		uint64_t nextFoundPos = result[foundPos & maskUint64_t];
		if (foundPos & firstBitUint64_t) finalPos ^= firstBitUint64_t;
		result[foundPos & maskUint64_t] = finalPos;
		foundPos = nextFoundPos;
	}
	return result[pos];
}

void mergeSegments(std::vector<uint64_t>& result, const size_t leftPos, const size_t rightPos, const bool fw)
{
	auto leftParent = find(result, leftPos);
	auto rightParent = find(result, rightPos);
	assert((result[leftParent & maskUint64_t] & maskUint64_t) == (leftParent & maskUint64_t));
	assert((result[rightParent & maskUint64_t] & maskUint64_t) == (rightParent & maskUint64_t));
	if ((leftParent & maskUint64_t) == (rightParent & maskUint64_t))
	{
		assert(((leftParent & firstBitUint64_t) == 0) ^ ((rightParent & firstBitUint64_t) == 0) ^ fw);
		return;
	}
	result[rightParent & maskUint64_t] = (leftParent & maskUint64_t) + ((((leftParent & firstBitUint64_t) == 0) ^ ((rightParent & firstBitUint64_t) == 0) ^ fw) ? 0 : firstBitUint64_t);
}

void mergeSegments(const std::vector<size_t>& readLengths, const std::vector<size_t>& countBeforeRead, std::vector<uint64_t>& result, const std::vector<RankBitvector>& breakpoints, const size_t leftRead, const bool leftFw, const size_t leftStart, const size_t leftEnd, const size_t rightRead, const bool rightFw, const size_t rightStart, const size_t rightEnd)
{
	assert(leftFw);
	assert(breakpoints[leftRead].get(leftStart));
	assert(breakpoints[leftRead].get(leftEnd));
	size_t leftFirst = breakpoints[leftRead].getRank(leftStart);
	size_t leftLast = breakpoints[leftRead].getRank(leftEnd);
	assert(leftLast > leftFirst);
	size_t rightFirst;
	size_t rightLast;
	if (rightFw)
	{
		assert(breakpoints[rightRead].get(rightStart));
		assert(breakpoints[rightRead].get(rightEnd));
		rightFirst = breakpoints[rightRead].getRank(rightStart);
		rightLast = breakpoints[rightRead].getRank(rightEnd);
	}
	else
	{
		assert(breakpoints[rightRead].get(readLengths[rightRead]-rightStart));
		assert(breakpoints[rightRead].get(readLengths[rightRead]-rightEnd));
		rightLast = breakpoints[rightRead].getRank(readLengths[rightRead]-rightStart);
		rightFirst = breakpoints[rightRead].getRank(readLengths[rightRead]-rightEnd);
	}
	assert(rightLast > rightFirst);
	assert(rightLast-rightFirst == leftLast-leftFirst);
	for (size_t i = 0; i < leftLast-leftFirst; i++)
	{
		if (rightFw)
		{
			mergeSegments(result, countBeforeRead[leftRead] + leftFirst + i, countBeforeRead[rightRead] + rightFirst+i, true);
		}
		else
		{
			mergeSegments(result, countBeforeRead[leftRead] + leftFirst + i, countBeforeRead[rightRead] + rightLast-i-1, false);
		}
	}
}

std::vector<uint64_t> mergeSegments(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const std::vector<RankBitvector>& breakpoints, const size_t countBreakpoints)
{
	assert(breakpoints.size() == readLengths.size());
	std::vector<uint64_t> result;
	assert(countBreakpoints >= 2*readLengths.size());
	result.resize(countBreakpoints - readLengths.size());
	std::vector<size_t> countBeforeRead;
	countBeforeRead.resize(readLengths.size(), 0);
	size_t countBeforeNow = 0;
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		countBeforeRead[i] = countBeforeNow;
		size_t index = 0;
		for (size_t j = 1; j < breakpoints[i].size(); j++)
		{
			if (!breakpoints[i].get(j)) continue;
			result[countBeforeNow + index] = countBeforeNow + index;
			index += 1;
		}
		countBeforeNow += index;
	}
	assert(countBeforeNow == result.size());
	for (size_t groupi = 0; groupi < matches.size(); groupi++)
	{
		for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
		{
			mergeSegments(readLengths, countBeforeRead,result, breakpoints, matches[groupi].leftRead, true, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightFw, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		find(result, i);
	}
	return result;
}

RankBitvector getSegmentToNode(const std::vector<uint64_t>& segments, const size_t minCoverage)
{
	RankBitvector result;
	result.resize(segments.size());
	for (size_t i = 0; i < segments.size(); i++)
	{
		assert((segments[i] & maskUint64_t) < segments.size());
		if ((segments[i] & maskUint64_t) != i) continue;
		assert(segments[i] == i);
		result.set(i, true);
	}
	result.buildRanks();
	return result;
}

MostlySparse2DHashmap<uint8_t, size_t> getEdgeCoverages(const std::vector<size_t>& readLengths, const RankBitvector& segmentToNode, const std::vector<uint64_t>& segments, const std::vector<RankBitvector>& breakpoints, const size_t minCoverage, const std::vector<size_t>& nodeCoverage, const size_t countNodes)
{
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverage;
	edgeCoverage.resize(countNodes);
	size_t i = 0;
	for (size_t readi = 0; readi < breakpoints.size(); readi++)
	{
		i += 1;
		for (size_t readpos = 1; readpos+1 < breakpoints[readi].size(); readpos++)
		{
			if (!breakpoints[readi].get(readpos)) continue;
			std::pair<size_t, bool> edgeFrom;
			edgeFrom.first = segmentToNode.getRank(segments[i-1] & maskUint64_t);
			edgeFrom.second = (segments[i-1] & firstBitUint64_t) == 0;
			if (nodeCoverage[edgeFrom.first] < minCoverage)
			{
				i += 1;
				continue;
			}
			std::pair<size_t, bool> edgeTo;
			edgeTo.first = segmentToNode.getRank(segments[i] & maskUint64_t);
			if (nodeCoverage[edgeTo.first] < minCoverage)
			{
				i += 1;
				continue;
			}
			edgeTo.second = (segments[i] & firstBitUint64_t) == 0;
			auto key = canon(edgeFrom, edgeTo);
			assert(key.first.first < edgeCoverage.size());
			assert(key.second.first < edgeCoverage.size());
			size_t coverage = 0;
			if (edgeCoverage.hasValue(key.first, key.second)) coverage = edgeCoverage.get(key.first, key.second);
			edgeCoverage.set(key.first, key.second, coverage+1);
			i += 1;
		}
	}
	assert(i == segments.size());
	return edgeCoverage;
}

std::vector<size_t> getNodeCoverage(const std::vector<uint64_t>& segments, const RankBitvector& segmentToNode, const size_t countNodes)
{
	std::vector<size_t> result;
	result.resize(countNodes, 0);
	for (size_t i = 0; i < segments.size(); i++)
	{
		uint64_t node = segmentToNode.getRank(segments[i] & maskUint64_t);
		assert(node < result.size());
		result[node] += 1;
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i] != 0);
	}
	return result;
}

std::vector<size_t> getNodeLengths(const std::vector<uint64_t>& segments, const RankBitvector& segmentToNode, const std::vector<RankBitvector>& breakpoints, const size_t countNodes)
{
	std::vector<size_t> result;
	result.resize(countNodes, std::numeric_limits<size_t>::max());
	size_t i = 0;
	for (size_t readi = 0; readi < breakpoints.size(); readi++)
	{
		size_t lastBreakpoint = 0;
		assert(breakpoints[readi].get(0));
		for (size_t pos = 1; pos < breakpoints[readi].size(); pos++)
		{
			if (!breakpoints[readi].get(pos)) continue;
			const size_t node = segmentToNode.getRank(segments[i] & maskUint64_t);
			const size_t length = pos - lastBreakpoint;
			assert(node < result.size());
			assert(length >= 1);
			if (!(result[node] == std::numeric_limits<size_t>::max() || result[node] == length))
			{
				std::cerr << length << " " << result[node] << std::endl;
				std::cerr << i << " " << segments[i] << " " << node << std::endl;
			}
			assert(result[node] == std::numeric_limits<size_t>::max() || result[node] == length);
			result[node] = length;
			i += 1;
			lastBreakpoint = pos;
		}
	}
	assert(i == segments.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i] != std::numeric_limits<size_t>::max());
	}
	return result;
}
