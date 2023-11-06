#include <cassert>
#include <phmap.h>
#include "AlnHaploFilter.h"
#include "Common.h"

class AlnSpan
{
public:
	size_t alignmentIndex;
	std::vector<std::pair<uint32_t, uint32_t>> matchSpans;
};

std::vector<std::pair<uint32_t, uint32_t>> getMergedSpans(const std::vector<std::pair<uint32_t, uint32_t>>& spans)
{
	assert(spans.size() >= 1);
	std::vector<std::pair<uint32_t, uint32_t>> result;
	result.emplace_back(spans[0].first, spans[0].second);
	assert(spans[0].second > spans[0].first);
	for (size_t i = 1; i < spans.size(); i++)
	{
		assert(spans[i].second > spans[i].first);
		assert(spans[i].first >= spans[i-1].first);
		if (spans[i].first > result.back().second)
		{
			result.emplace_back(spans[i].first, spans[i].second);
		}
		else
		{
			result.back().second = std::max(result.back().second, spans[i].second);
		}
	}
	return result;
}

void forbidPerRead(const size_t readi, const std::vector<AlnSpan>& spansHere, std::vector<bool>& kept, const size_t minHapCoverage, const size_t maxCrossCoverage)
{
	assert(spansHere.size() >= 1);
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> gapCount;
	for (size_t i = 0; i < spansHere.size(); i++)
	{
		assert(spansHere[i].matchSpans.size() >= 1);
		assert(spansHere[i].matchSpans[0].first < spansHere[i].matchSpans[0].second);
		for (size_t j = 1; j < spansHere[i].matchSpans.size(); j++)
		{
			assert(spansHere[i].matchSpans[j-1].second < spansHere[i].matchSpans[j].first);
			assert(spansHere[i].matchSpans[j].first < spansHere[i].matchSpans[j].second);
			gapCount[std::make_pair(spansHere[i].matchSpans[j-1].second, spansHere[i].matchSpans[j].first)] += 1;
		}
	}
	std::vector<std::pair<size_t, size_t>> potentialGaps;
	std::vector<std::vector<size_t>> gapRefs;
	std::vector<std::vector<size_t>> gapAlts;
	for (auto pair : gapCount)
	{
		if (pair.second < minHapCoverage) continue;
		potentialGaps.emplace_back(pair.first);
	}
	std::sort(potentialGaps.begin(), potentialGaps.end(), [](auto left, auto right){ return left.first < right.first; });
	gapRefs.resize(potentialGaps.size());
	gapAlts.resize(potentialGaps.size());
	for (size_t i = 0; i < spansHere.size(); i++)
	{
		size_t potentialGapIndex = 0;
		for (size_t j = 0; j < spansHere[i].matchSpans.size(); j++)
		{
			while (potentialGapIndex < potentialGaps.size() && spansHere[i].matchSpans[j].first > potentialGaps[potentialGapIndex].second) potentialGapIndex += 1;
			if (j >= 1)
			{
				if (potentialGapIndex < potentialGaps.size() && spansHere[i].matchSpans[j-1].second == potentialGaps[potentialGapIndex].first && spansHere[i].matchSpans[j].first == potentialGaps[potentialGapIndex].second)
				{
					gapAlts[potentialGapIndex].emplace_back(i);
					potentialGapIndex += 1;
				}
			}
			while (potentialGapIndex < potentialGaps.size() && spansHere[i].matchSpans[j].first > potentialGaps[potentialGapIndex].first) potentialGapIndex += 1;
			while (potentialGapIndex < potentialGaps.size() && spansHere[i].matchSpans[j].first <= potentialGaps[potentialGapIndex].first && spansHere[i].matchSpans[j].second >= potentialGaps[potentialGapIndex].second)
			{
				gapRefs[potentialGapIndex].emplace_back(i);
				potentialGapIndex += 1;
			}
			while (potentialGapIndex < potentialGaps.size() && spansHere[i].matchSpans[j].second > potentialGaps[potentialGapIndex].first) potentialGapIndex += 1;
		}
		assert(potentialGapIndex == potentialGaps.size() || potentialGaps[potentialGapIndex].first >= spansHere[i].matchSpans.back().second);
	}
	for (size_t i = 0; i < potentialGaps.size(); i++)
	{
		std::sort(gapRefs[i].begin(), gapRefs[i].end());
		std::sort(gapAlts[i].begin(), gapAlts[i].end());
	}
	for (size_t i = 0; i < potentialGaps.size(); i++)
	{
		if (gapRefs[i].size() < minHapCoverage) continue;
		if (gapAlts[i].size() < minHapCoverage) continue;
		for (size_t j = i+1; j < potentialGaps.size(); j++)
		{
			if (gapRefs[j].size() < minHapCoverage) continue;
			if (gapAlts[j].size() < minHapCoverage) continue;
			if (intersectSize(gapRefs[i], gapRefs[j]) < minHapCoverage) continue;
			if (intersectSize(gapAlts[i], gapAlts[j]) < minHapCoverage) continue;
			if (intersectSize(gapRefs[i], gapAlts[j]) + intersectSize(gapAlts[i], gapRefs[j]) >= maxCrossCoverage) continue;
			size_t ji = 0;
			for (auto n : gapAlts[i])
			{
				while (ji < gapAlts[j].size() && gapAlts[j][ji] < n) ji++;
				if (gapAlts[j][ji] == n) kept[n] = false;
			}
		}
	}
}

std::vector<bool> getValidAlignments(const std::vector<MatchGroup>& matches, const std::vector<size_t>& readLengths, const size_t minHapCoverage, const size_t maxCrossCoverage)
{
	std::vector<std::vector<AlnSpan>> spansPerRead;
	spansPerRead.resize(readLengths.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		spansPerRead[matches[i].leftRead].emplace_back();
		spansPerRead[matches[i].leftRead].back().alignmentIndex = i;
		std::vector<std::pair<uint32_t, uint32_t>> lefts;
		for (auto match : matches[i].matches)
		{
			lefts.emplace_back(matches[i].leftStart + (size_t)match.leftStart, matches[i].leftStart + (size_t)match.leftStart + (size_t)match.length);
		}
		std::sort(lefts.begin(), lefts.end(), [](auto left, auto right){ return left.first < right.first; });
		spansPerRead[matches[i].leftRead].back().matchSpans = getMergedSpans(lefts);
		spansPerRead[matches[i].rightRead].emplace_back();
		spansPerRead[matches[i].rightRead].back().alignmentIndex = i;
		std::vector<std::pair<uint32_t, uint32_t>> rights;
		for (auto match : matches[i].matches)
		{
			if (matches[i].rightFw)
			{
				rights.emplace_back(matches[i].rightStart + (size_t)match.rightStart, matches[i].rightStart + (size_t)match.rightStart + (size_t)match.length);
			}
			else
			{
				rights.emplace_back(readLengths[matches[i].rightRead] - (matches[i].rightStart + (size_t)match.rightStart + (size_t)match.length), readLengths[matches[i].rightRead] - (matches[i].rightStart + (size_t)match.rightStart));
			}
		}
		std::sort(rights.begin(), rights.end(), [](auto left, auto right){ return left.first < right.first; });
		spansPerRead[matches[i].rightRead].back().matchSpans = getMergedSpans(rights);
	}
	std::vector<bool> kept;
	kept.resize(matches.size(), true);
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		if (spansPerRead[i].size() == 0) continue;
		forbidPerRead(i, spansPerRead[i], kept, minHapCoverage, maxCrossCoverage);
	}
	return kept;
}


