#include <limits>
#include "OverlapMatcher.h"

size_t arraySize(const std::array<size_t, 5>& array)
{
	if (array[0] == std::numeric_limits<size_t>::max()) return 0;
	if (array[1] == std::numeric_limits<size_t>::max()) return 1;
	if (array[2] == std::numeric_limits<size_t>::max()) return 2;
	if (array[3] == std::numeric_limits<size_t>::max()) return 3;
	if (array[4] == std::numeric_limits<size_t>::max()) return 4;
	return 5;
}

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

phmap::flat_hash_map<size_t, std::array<size_t, 5>> getLowCountKmers(const std::vector<std::pair<size_t, size_t>>& readKmers)
{
	const size_t maxSize = 5;
	phmap::flat_hash_map<size_t, std::array<size_t, 5>> result;
	for (auto pair : readKmers)
	{
		if (result.count(pair.second) == 0)
		{
			result[pair.second][0] = pair.first;
			result[pair.second][1] = std::numeric_limits<size_t>::max();
			result[pair.second][2] = std::numeric_limits<size_t>::max();
			result[pair.second][3] = std::numeric_limits<size_t>::max();
			result[pair.second][4] = std::numeric_limits<size_t>::max();
		}
		else
		{
			std::array<size_t, 5>& arr = result[pair.second];
			if (arraySize(arr) < 5) arr[arraySize(arr)] = pair.first;
		}
	}
	phmap::flat_hash_set<size_t> removeThese;
	for (auto pair : result)
	{
		size_t s = arraySize(pair.second);
		assert(s <= maxSize);
		if (s == maxSize)
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

int getBestDiagonal(const phmap::flat_hash_map<size_t, std::array<size_t, 5>>& refKmers, const phmap::flat_hash_map<size_t, std::array<size_t, 5>>& queryKmers)
{
	std::vector<int> diagonalMatches;
	for (const auto& pair : refKmers)
	{
		if (queryKmers.count(pair.first) == 0) continue;
		for (auto pos : pair.second)
		{
			if (pos == std::numeric_limits<size_t>::max()) break;
			for (auto pos2 : queryKmers.at(pair.first))
			{
				if (pos2 == std::numeric_limits<size_t>::max()) break;
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

ReadOverlapInformation getReadOverlapInformation(const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t readIndex, const size_t kmerSize)
{
	ReadOverlapInformation result;
	std::string readSeq;
	if (readIndex < sequenceIndex.size())
	{
		readSeq = sequenceIndex.getSequence(readIndex);
	}
	else
	{
		readSeq = revCompRaw(sequenceIndex.getSequence(readIndex - sequenceIndex.size()));
	}
	result.readKmers = getLocallyUniqueKmers(readSeq, kmerSize);
	result.lowCountKmers = getLowCountKmers(result.readKmers);
	result.readIndex = readIndex;
	result.readLength = readSeq.size();
	return result;
}
