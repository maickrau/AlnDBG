#ifndef OverlapMatcher_h
#define OverlapMatcher_h

#include <iostream>
#include <cstddef>
#include <vector>
#include <string>
#include "Common.h"
#include "CompressedStringIndex.h"
#include "phmap.h"

class ReadOverlapInformation
{
public:
	std::vector<std::pair<size_t, size_t>> readKmers;
	size_t readIndex;
	size_t readLength;
};

std::vector<std::pair<size_t, size_t>> getLocallyUniqueKmers(const std::string& readSeq, const size_t kmerSize);
void filterOutDoubleMatchedPositions(std::vector<std::pair<size_t, size_t>>& matches);
ReadOverlapInformation getReadOverlapInformation(const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t readIndex, const size_t kmerSize);

// kmers with same sequence are iterated ordered by their position, but kmers with different sequence will not be
template <typename F>
void iterateLocallyUniqueKmers(const std::string& sequence, const size_t kmerSize, const size_t localitySize, F callback)
{
	static_assert(sizeof(size_t) == 8, "size_t should be 8 bytes");
	//if size_t is not 8 bytes, this constant needs to be changed
	const size_t kmerIsNotLocallyUniqueBit = 0x8000000000000000ull;
	phmap::flat_hash_map<size_t, size_t> lastKmerPosition;
	size_t kmer = 0;
	const size_t mask = (1ull << (2ull*kmerSize)) - 1ull;
	for (size_t i = 0; i < kmerSize; i++)
	{
		kmer <<= 2;
		switch(sequence[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
		}
	}
	lastKmerPosition[kmer] = 0;
	for (size_t i = kmerSize; i < sequence.size(); i++)
	{
		kmer <<= 2;
		switch(sequence[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
		}
		kmer &= mask;
		size_t posHere = i-kmerSize+1;
		if (lastKmerPosition.count(kmer) == 0)
		{
			lastKmerPosition[kmer] = posHere;
		}
		else if ((lastKmerPosition.at(kmer) & ~kmerIsNotLocallyUniqueBit) + localitySize > posHere)
		{
			lastKmerPosition[kmer] = posHere | kmerIsNotLocallyUniqueBit;
		}
		else
		{
			if ((lastKmerPosition.at(kmer) & kmerIsNotLocallyUniqueBit) == 0) callback(kmer, lastKmerPosition.at(kmer));
			lastKmerPosition[kmer] = posHere;
		}
	}
	for (auto pair : lastKmerPosition)
	{
		if (pair.second & kmerIsNotLocallyUniqueBit) continue;
		callback(pair.first, pair.second);
	}
}

template <typename F>
void iterateMatchingKmersInDiagonal(const std::vector<std::pair<size_t, size_t>>& refReadKmers, const std::vector<std::pair<size_t, size_t>>& queryReadKmers, const int diagonal, F callback)
{
	phmap::flat_hash_map<size_t, size_t> activeRefKmers;
	size_t refIndex = 0;
	size_t queryIndex = 0;
	while (refIndex < refReadKmers.size() && queryIndex < queryReadKmers.size())
	{
		int diagonalHere = (int)refReadKmers[refIndex].first - (int)queryReadKmers[queryIndex].first;
		activeRefKmers[refReadKmers[refIndex].second] = refReadKmers[refIndex].first;
		size_t queryKmer = queryReadKmers[queryIndex].second;
		if (activeRefKmers.count(queryKmer) == 1)
		{
			int matchDiagonal = (int)activeRefKmers.at(queryKmer) - (int)queryReadKmers[queryIndex].first;
			if (matchDiagonal >= diagonal-100 && matchDiagonal <= diagonal+100)
			{
				callback(activeRefKmers.at(queryKmer), queryReadKmers[queryIndex].first);
			}
		}
		if (diagonalHere < diagonal+100)
		{
			refIndex += 1;
			continue;
		}
		if (diagonalHere > diagonal+100)
		{
			queryIndex += 1;
			continue;
		}
		refIndex += 1;
		queryIndex += 1;
		continue;
	}
}

template <typename F>
void iterateUniqueKmerMatches(const ReadOverlapInformation& refRead, const ReadOverlapInformation& queryRead, const int diagonal, F callback)
{
	if (diagonal == std::numeric_limits<int>::max()) return;
	std::vector<std::pair<size_t, size_t>> matches;
	iterateMatchingKmersInDiagonal(refRead.readKmers, queryRead.readKmers, diagonal, [&matches](const size_t refPos, const size_t queryPos)
	{
		matches.emplace_back(refPos, queryPos);
	});
	filterOutDoubleMatchedPositions(matches);
	if (matches.size() < 10) return;
	bool hasBigGap = false;
	bool matchesStart = false;
	bool matchesEnd = false;
	std::sort(matches.begin(), matches.end(), [](auto left, auto right) { return left.second < right.second; });
	size_t minMaxMatchLength = matches.back().second - matches[0].second;
	if (matches[0].second < 1000) matchesStart = true;
	if (matches.back().second + 1000 >= queryRead.readLength) matchesEnd = true;
	for (size_t i = 1; i < matches.size(); i++)
	{
		assert(matches[i].second >= matches[i-1].second);
		if (matches[i].second - matches[i-1].second > 1000)
		{
			hasBigGap = true;
			break;
		}
	}
	if (hasBigGap) return;
	if (minMaxMatchLength < 5000) return;
	std::sort(matches.begin(), matches.end());
	minMaxMatchLength = std::min(minMaxMatchLength, matches.back().first - matches[0].first);
	if (matches[0].first < 1000) matchesStart = true;
	if (matches.back().first + 1000 >= refRead.readLength) matchesEnd = true;
	for (size_t i = 1; i < matches.size(); i++)
	{
		assert(matches[i].first >= matches[i-1].first);
		if (matches[i].first - matches[i-1].first > 1000)
		{
			hasBigGap = true;
			break;
		}
	}
	if (hasBigGap) return;
	if (!matchesStart) return;
	if (!matchesEnd) return;
	if (minMaxMatchLength < 5000) return;
	for (size_t i = 0; i < matches.size(); i++)
	{
		callback(queryRead.readIndex, matches[i].first, matches[i].second);
	}
}

template <typename F>
void iterateUniqueKmerMatches(const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t refreadIndex, const std::vector<size_t>& matchingReads, const std::vector<int>& matchingDiagonals, const size_t kmerSize, F callback)
{
	ReadOverlapInformation refRead = getReadOverlapInformation(sequenceIndex, refreadIndex, kmerSize);
	assert(matchingReads.size() == matchingDiagonals.size());
	for (size_t i = 0; i < matchingReads.size(); i++)
	{
		const size_t read = matchingReads[i];
		const int diagonal = matchingDiagonals[i];
		if (read == refreadIndex) continue;
		if (refreadIndex < sequenceIndex.size() && read == refreadIndex + sequenceIndex.size()) continue;
		if (refreadIndex >= sequenceIndex.size() && read == refreadIndex - sequenceIndex.size()) continue;
		ReadOverlapInformation queryRead = getReadOverlapInformation(sequenceIndex, read, kmerSize);
		iterateUniqueKmerMatches(refRead, queryRead, diagonal, callback);
	}
}

#endif
