#ifndef OverlapMatcher_h
#define OverlapMatcher_h

#include <iostream>
#include <cstddef>
#include <vector>
#include <string>
#include "Common.h"
#include "CompressedStringIndex.h"
#include "phmap.h"

std::vector<std::pair<size_t, size_t>> getLocallyUniqueKmers(const std::string& readSeq, const size_t kmerSize);
phmap::flat_hash_map<size_t, std::vector<size_t>> getLowCountKmers(const std::vector<std::pair<size_t, size_t>>& readKmers);
// plus: starts at positive location in ref. minus: starts at negative location in ref. starts at 0 in query
int getBestDiagonal(const phmap::flat_hash_map<size_t, std::vector<size_t>>& refKmers, const phmap::flat_hash_map<size_t, std::vector<size_t>>& queryKmers);
void filterOutDoubleMatchedPositions(std::vector<std::pair<size_t, size_t>>& matches);

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
	phmap::flat_hash_map<size_t, size_t> kmers;
	kmers[kmer] += 1;
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
		kmers[kmer] += 1;
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
	size_t uniqueKmers = 0;
	for (auto pair : kmers)
	{
		if (pair.second == 1) uniqueKmers += 1;
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
void iterateUniqueKmerMatches(const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t refreadIndex, const phmap::flat_hash_set<size_t>& matchingReads, const size_t kmerSize, F callback)
{
	std::string refSeq;
	if (refreadIndex < sequenceIndex.size())
	{
		refSeq = sequenceIndex.getSequence(refreadIndex);
	}
	else
	{
		refSeq = revCompRaw(sequenceIndex.getSequence(refreadIndex - sequenceIndex.size()));
	}
	std::vector<std::pair<size_t, size_t>> refReadKmers = getLocallyUniqueKmers(refSeq, kmerSize);
	if (refReadKmers.size() < 10) return;
	phmap::flat_hash_map<size_t, std::vector<size_t>> lowCountRefKmers = getLowCountKmers(refReadKmers);
	for (size_t read : matchingReads)
	{
		if (read == refreadIndex) continue;
		if (refreadIndex < sequenceIndex.size() && read == refreadIndex + sequenceIndex.size()) continue;
		if (refreadIndex >= sequenceIndex.size() && read == refreadIndex - sequenceIndex.size()) continue;
		std::string querySeq;
		if (read < sequenceIndex.size())
		{
			querySeq = sequenceIndex.getSequence(read);
		}
		else
		{
			querySeq = revCompRaw(sequenceIndex.getSequence(read - sequenceIndex.size()));
		}
		std::vector<std::pair<size_t, size_t>> queryReadKmers = getLocallyUniqueKmers(querySeq, kmerSize);
		if (queryReadKmers.size() < 10) continue;
		phmap::flat_hash_map<size_t, std::vector<size_t>> lowCountQueryKmers = getLowCountKmers(queryReadKmers);
		int bestDiagonal = getBestDiagonal(lowCountRefKmers, lowCountQueryKmers);
		if (bestDiagonal == std::numeric_limits<int>::max()) continue;
		std::vector<std::pair<size_t, size_t>> matches;
		iterateMatchingKmersInDiagonal(refReadKmers, queryReadKmers, bestDiagonal, [&matches](const size_t refPos, const size_t queryPos)
		{
			matches.emplace_back(refPos, queryPos);
		});
		filterOutDoubleMatchedPositions(matches);
		if (matches.size() < 10) continue;
		bool hasBigGap = false;
		bool matchesStart = false;
		bool matchesEnd = false;
		std::sort(matches.begin(), matches.end(), [](auto left, auto right) { return left.second < right.second; });
		size_t minMaxMatchLength = matches.back().second - matches[0].second;
		if (matches[0].second < 1000) matchesStart = true;
		if (matches.back().second + 1000 >= querySeq.size()) matchesEnd = true;
		for (size_t i = 1; i < matches.size(); i++)
		{
			assert(matches[i].second >= matches[i-1].second);
			if (matches[i].second - matches[i-1].second > 1000)
			{
				hasBigGap = true;
				break;
			}
		}
		std::sort(matches.begin(), matches.end());
		minMaxMatchLength = std::min(minMaxMatchLength, matches.back().first - matches[0].first);
		if (matches[0].first < 1000) matchesStart = true;
		if (matches.back().first + 1000 >= refSeq.size()) matchesEnd = true;
		for (size_t i = 1; i < matches.size(); i++)
		{
			assert(matches[i].first >= matches[i-1].first);
			if (matches[i].first - matches[i-1].first > 1000)
			{
				hasBigGap = true;
				break;
			}
		}
		if (hasBigGap) continue;
		if (!matchesStart) continue;
		if (!matchesEnd) continue;
		if (minMaxMatchLength < 5000) continue;
		for (size_t i = 0; i < matches.size(); i++)
		{
			callback(read, matches[i].first, matches[i].second);
		}
	}
}

#endif
