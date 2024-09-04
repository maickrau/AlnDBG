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
	result.readIndex = readIndex;
	result.readLength = readSeq.size();
	return result;
}
