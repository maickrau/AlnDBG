#include <limits>
#include "OverlapMatcher.h"

std::vector<std::pair<size_t, size_t>> getLocallyUniqueKmers(const std::string& readSeq, const size_t kmerSize)
{
	std::vector<std::pair<size_t, size_t>> result;
	iterateLocallyUniqueKmers(readSeq, kmerSize, 100, [&result](const size_t kmer, const size_t pos)
	{
		result.emplace_back(pos, kmer);
	});
	std::sort(result.begin(), result.end());
	return result;
}

phmap::flat_hash_map<size_t, std::vector<size_t>> getLowCountKmers(const std::vector<std::pair<size_t, size_t>>& readKmers)
{
	const size_t maxSize = 5;
	phmap::flat_hash_map<size_t, std::vector<size_t>> result;
	for (auto pair : readKmers)
	{
		if (result.count(pair.second) == 0 || result.at(pair.second).size() < maxSize)
		{
			result[pair.second].emplace_back(pair.first);
		}
	}
	phmap::flat_hash_set<size_t> removeThese;
	for (auto pair : result)
	{
		assert(pair.second.size() <= maxSize);
		if (pair.second.size() == maxSize)
		{
			removeThese.insert(pair.first);
		}
	}
	for (auto kmer : removeThese)
	{
		result.erase(kmer);
	}
	return result;
}

int getBestDiagonal(const phmap::flat_hash_map<size_t, std::vector<size_t>>& refKmers, const phmap::flat_hash_map<size_t, std::vector<size_t>>& queryKmers)
{
	std::vector<int> diagonalMatches;
	for (const auto& pair : refKmers)
	{
		if (queryKmers.count(pair.first) == 0) continue;
		for (auto pos : pair.second)
		{
			for (auto pos2 : queryKmers.at(pair.first))
			{
				int diagonal = (int)pos - (int)pos2;
				diagonalMatches.emplace_back(diagonal);
			}
		}
	}
	if (diagonalMatches.size() < 10) return std::numeric_limits<int>::max();
	std::sort(diagonalMatches.begin(), diagonalMatches.end());
	size_t currentBandStart = 0;
	size_t bestBandStart = 0;
	size_t bestBandEnd = 0;
	for (size_t bandEnd = 1; bandEnd < diagonalMatches.size(); bandEnd++)
	{
		assert(currentBandStart < bandEnd);
		while (diagonalMatches[currentBandStart] + 100 < diagonalMatches[bandEnd]) currentBandStart += 1;
		assert(currentBandStart <= bandEnd);
		if (bandEnd - currentBandStart > bestBandEnd - bestBandStart)
		{
			bestBandEnd = bandEnd;
			bestBandStart = currentBandStart;
		}
	}
	size_t index = (bestBandEnd - bestBandStart) / 2 + bestBandStart;
	assert(index < diagonalMatches.size());
	return diagonalMatches[index];
}

void filterOutDoubleMatchedPositions(std::vector<std::pair<size_t, size_t>>& matches)
{
	phmap::flat_hash_set<size_t> matchedAtFirst;
	phmap::flat_hash_set<size_t> matchedAtSecond;
	phmap::flat_hash_set<size_t> needsToRemoveFirst;
	phmap::flat_hash_set<size_t> needsToRemoveSecond;
	for (auto pair : matches)
	{
		if (matchedAtFirst.count(pair.first) == 1)
		{
			needsToRemoveFirst.insert(pair.first);
		}
		if (matchedAtSecond.count(pair.second) == 1)
		{
			needsToRemoveSecond.insert(pair.second);
		}
		matchedAtFirst.insert(pair.first);
		matchedAtSecond.insert(pair.second);
	}
	for (size_t i = matches.size()-1; i < matches.size(); i--)
	{
		if (needsToRemoveFirst.count(matches[i].first) == 1 || needsToRemoveSecond.count(matches[i].second) == 1)
		{
			std::swap(matches[i], matches.back());
			matches.pop_back();
		}
	}
}
