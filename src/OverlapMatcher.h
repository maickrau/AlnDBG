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
		if (refReadKmers[refIndex].first < queryReadKmers[queryIndex].first)
		{
			refIndex += 1;
			continue;
		}
		if (refReadKmers[refIndex].first > queryReadKmers[queryIndex].first)
		{
			queryIndex += 1;
			continue;
		}
		assert(refReadKmers[refIndex].first == queryReadKmers[queryIndex].first);
		if ((refIndex+1 == refReadKmers.size() || refReadKmers[refIndex+1].first != refReadKmers[refIndex].first) && (queryIndex+1 == queryReadKmers.size() || queryReadKmers[queryIndex+1].first != queryReadKmers[queryIndex].first))
		{
			int diagonalHere = (int)refReadKmers[refIndex].second - (int)queryReadKmers[queryIndex].second;
			if (diagonalHere >= diagonal-100 && diagonalHere <= diagonal+100)
			{
				callback(refReadKmers[refIndex].second, queryReadKmers[queryIndex].second);
			}
			refIndex += 1;
			queryIndex += 1;
			continue;
		}
		while (refIndex < refReadKmers.size() && queryIndex < queryReadKmers.size() && refReadKmers[refIndex].first == queryReadKmers[queryIndex].first)
		{
			int diagonalHere = (int)refReadKmers[refIndex].second - (int)queryReadKmers[queryIndex].second;
			if (diagonalHere >= diagonal-100 && diagonalHere <= diagonal+100)
			{
				bool matchedByAnother = false;
				if (refIndex > 0 && refReadKmers[refIndex-1].first == refReadKmers[refIndex].first)
				{
					int otherDiagonal = (int)refReadKmers[refIndex-1].second - (int)queryReadKmers[queryIndex].second;
					if (otherDiagonal >= diagonal-100 && otherDiagonal <= diagonal+100) matchedByAnother = true;
				}
				if (refIndex+1 < refReadKmers.size() && refReadKmers[refIndex+1].first == refReadKmers[refIndex].first)
				{
					int otherDiagonal = (int)refReadKmers[refIndex+1].second - (int)queryReadKmers[queryIndex].second;
					if (otherDiagonal >= diagonal-100 && otherDiagonal <= diagonal+100) matchedByAnother = true;
				}
				if (queryIndex > 0 && queryReadKmers[queryIndex-1].first == queryReadKmers[queryIndex].first)
				{
					int otherDiagonal = (int)refReadKmers[refIndex].second - (int)queryReadKmers[queryIndex-1].second;
					if (otherDiagonal >= diagonal-100 && otherDiagonal <= diagonal+100) matchedByAnother = true;
				}
				if (queryIndex+1 < queryReadKmers.size() && queryReadKmers[queryIndex+1].first == queryReadKmers[queryIndex].first)
				{
					int otherDiagonal = (int)refReadKmers[refIndex].second - (int)queryReadKmers[queryIndex+1].second;
					if (otherDiagonal >= diagonal-100 && otherDiagonal <= diagonal+100) matchedByAnother = true;
				}
				if (!matchedByAnother) callback(refReadKmers[refIndex].second, queryReadKmers[queryIndex].second);
			}
			if (diagonalHere < diagonal)
			{
				refIndex += 1;
				continue;
			}
			if (diagonalHere > diagonal)
			{
				queryIndex += 1;
				continue;
			}
			assert(diagonalHere == diagonal);
			refIndex += 1;
			queryIndex += 1;
			if (refIndex == refReadKmers.size() || refReadKmers[refIndex].first != refReadKmers[refIndex-1].first)
			{
				break;
			}
		}
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
	if (matches.size() < 10) return;
	std::sort(matches.begin(), matches.end(), [](auto left, auto right) { return left.second < right.second; });
	bool matchesStart = false;
	bool matchesEnd = false;
	size_t lastValidMatchPos = std::numeric_limits<size_t>::max();
	size_t lastAnyMatchPos = std::numeric_limits<size_t>::max();
	size_t maxMatchPos = std::numeric_limits<size_t>::max();
	size_t minMatchPos = std::numeric_limits<size_t>::max();
	for (size_t i = matches.size()-1; i > 0; i--)
	{
		assert(matches[i].second >= matches[i-1].second);
		if (matches[i].second == matches[i-1].second || matches[i].second == lastAnyMatchPos)
		{
			lastAnyMatchPos = matches[i].second;
			std::swap(matches[i], matches.back());
			matches.pop_back();
			if (matches.size() < 10) return;
			continue;
		}
		assert(matches[i].second > matches[i-1].second);
		assert(matches[i].second < lastAnyMatchPos);
		if (maxMatchPos == std::numeric_limits<size_t>::max())
		{
			maxMatchPos = matches[i].second;
		}
		assert(minMatchPos > matches[i].second);
		minMatchPos = matches[i].second;
		if (lastValidMatchPos != std::numeric_limits<size_t>::max())
		{
			assert(lastValidMatchPos > matches[i].second);
			if (lastValidMatchPos - matches[i].second > 1000)
			{
				return; // has big gap
			}
		}
		lastValidMatchPos = matches[i].second;
		lastAnyMatchPos = matches[i].second;
	}
	if (lastAnyMatchPos > matches[0].second)
	{
		assert(minMatchPos > matches[0].second);
		minMatchPos = matches[0].second;
		assert(lastValidMatchPos != std::numeric_limits<size_t>::max());
		assert(lastValidMatchPos > matches[0].second);
		if (lastValidMatchPos - matches[0].second > 1000)
		{
			return; // has big gap
		}
	}
	if (minMatchPos < 1000) matchesStart = true;
	if (maxMatchPos+1000 > queryRead.readLength) matchesEnd = true;
	size_t minMaxMatchLength = maxMatchPos - minMatchPos;
	if (minMaxMatchLength < 5000) return;
	std::sort(matches.begin(), matches.end());
	std::vector<size_t> skipMatchSites; // so no need to sort matches again
	lastValidMatchPos = std::numeric_limits<size_t>::max();
	lastAnyMatchPos = std::numeric_limits<size_t>::max();
	maxMatchPos = std::numeric_limits<size_t>::max();
	minMatchPos = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i+1 < matches.size(); i++)
	{
		assert(matches[i].first <= matches[i+1].first);
		if (matches[i].first == matches[i+1].first || matches[i].first == lastAnyMatchPos)
		{
			lastAnyMatchPos = matches[i].first;
			skipMatchSites.emplace_back(i);
			assert(matches.size() > skipMatchSites.size());
			if (matches.size() - skipMatchSites.size() < 10) return;
			continue;
		}
		assert(matches[i].first < matches[i+1].first);
		assert(lastAnyMatchPos == std::numeric_limits<size_t>::max() || matches[i].first > lastAnyMatchPos);
		assert(lastValidMatchPos == std::numeric_limits<size_t>::max() || matches[i].first > lastValidMatchPos);
		if (lastValidMatchPos != std::numeric_limits<size_t>::max())
		{
			assert(matches[i].first > lastValidMatchPos);
			if (matches[i].first - lastValidMatchPos > 1000)
			{
				return; // has big gap
			}
		}
		if (minMatchPos == std::numeric_limits<size_t>::max())
		{
			minMatchPos = matches[i].first;
			if (minMatchPos > 1000 && !matchesStart) return;
		}
		assert(maxMatchPos == std::numeric_limits<size_t>::max() || matches[i].first > maxMatchPos);
		maxMatchPos = matches[i].first;
		lastValidMatchPos = matches[i].first;
		lastAnyMatchPos = matches[i].first;
	}
	if (matches.back().first > lastAnyMatchPos)
	{
		assert(maxMatchPos == std::numeric_limits<size_t>::max() || matches.back().first > maxMatchPos);
		maxMatchPos = matches.back().first;
		assert(lastValidMatchPos != std::numeric_limits<size_t>::max());
		assert(matches.back().first > lastValidMatchPos);
		if (matches.back().first - lastValidMatchPos > 1000)
		{
			return; // has big gap
		}
	}
	else
	{
		skipMatchSites.push_back(matches.size()-1);
	}
	if (minMatchPos < 1000) matchesStart = true;
	if (maxMatchPos + 1000 > refRead.readLength) matchesEnd = true;
	minMaxMatchLength = std::min(minMaxMatchLength, maxMatchPos - minMatchPos);
	if (!matchesStart) return;
	if (!matchesEnd) return;
	if (minMaxMatchLength < 5000) return;
	assert(matches.size() > skipMatchSites.size());
	if (matches.size() - skipMatchSites.size() < 10) return;
	size_t skipIndex = 0;
	for (size_t i = 0; i < matches.size(); i++)
	{
		if (skipIndex < skipMatchSites.size())
		{
			while (skipMatchSites[skipIndex] < i)
			{
				assert(skipIndex+1 < skipMatchSites.size());
				skipIndex += 1;
			}
			if (skipMatchSites[skipIndex] == i)
			{
				skipIndex += 1;
				continue;
			}
		}
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
