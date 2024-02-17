#include <cassert>
#include <mutex>
#include <thread>
#include <phmap.h>
#include <iostream>
#include "KmerMatcher.h"
#include "FastHasher.h"

template <typename F>
void iterateKmerMatchPositions(const uint64_t kmer, const phmap::flat_hash_map<uint64_t, uint64_t>& firstPositions, const phmap::flat_hash_map<uint64_t, std::vector<uint16_t>>& extraPositions, F callback)
{
	auto found = firstPositions.find(kmer);
	if (found == firstPositions.end()) return;
	callback(found->second & 0x000000000000FFFFull);
	if ((found->second & 0x00000000FFFF0000ull) == 0x00000000FFFF0000ull) return;
	callback((found->second >> 16ull) & 0x000000000000FFFFull);
	if ((found->second & 0x0000FFFF00000000ull) == 0x0000FFFF00000000ull) return;
	callback((found->second >> 32ull) & 0x000000000000FFFFull);
	if ((found->second & 0xFFFF000000000000ull) == 0xFFFF000000000000ull) return;
	callback((found->second >> 48ull) & 0x000000000000FFFFull);
	auto found2 = extraPositions.find(kmer);
	if (found2 == extraPositions.end()) return;
	for (auto pos : found2->second) callback(pos);
}

template <typename F1, typename F2>
void iterateSyncmersBigK(const std::vector<TwobitString>& readSequences, const size_t k, const size_t w, const size_t maxLen, F1 callback, F2 getchar)
{
	thread_local std::vector<std::tuple<size_t, uint64_t>> smerOrder;
	assert(w < k);
	assert(w >= 3);
	const size_t s = k-w+1;
	MBG::FastHasher kmer { k };
	MBG::FastHasher smer { s };
	assert(kmer.hash() == 0);
	assert(kmer.getFwHash() == 0);
	assert(smer.hash() == 0);
	assert(smer.getFwHash() == 0);
	for (size_t i = 0; i < s; i++)
	{
		uint64_t c = getchar(i);
		smer.addChar(c);
		kmer.addChar(c);
	}
	assert(smerOrder.size() == 0);
	smerOrder.emplace_back(0, smer.getFwHash());
	for (size_t i = s; i < k; i++)
	{
		uint64_t c = getchar(i);
		kmer.addChar(c);
		smer.addChar(c);
		smer.removeChar(getchar(i-s));
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > smer.getFwHash()) smerOrder.pop_back();
		smerOrder.emplace_back(i-s+1, smer.getFwHash());
	}
	if ((std::get<0>(smerOrder.front()) == 0) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == w-1))
	{
		callback(kmer.getFwHash(), 0);
	}
	for (size_t i = k; i < maxLen; i++)
	{
		uint64_t c = getchar(i);
		kmer.addChar(c);
		smer.addChar(c);
		smer.removeChar(getchar(i-s));
		kmer.removeChar(getchar(i-k));
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > smer.getFwHash()) smerOrder.pop_back();
		while (smerOrder.size() > 0 && std::get<0>(smerOrder.front()) <= i-s+1-w) smerOrder.erase(smerOrder.begin());
		smerOrder.emplace_back(i-s+1, smer.getFwHash());
		if ((std::get<0>(smerOrder.front()) == i-s+2-w) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == i-s+1))
		{
			callback(kmer.getFwHash(), i-k+1);
		}
	}
	smerOrder.clear();
}

template <typename F1, typename F2>
void iterateSyncmers(const std::vector<TwobitString>& readSequences, const size_t k, const size_t w, const size_t maxLen, F1 callback, F2 getchar)
{
	thread_local std::vector<std::tuple<size_t, uint64_t>> smerOrder;
	assert(k <= 31);
	assert(w < k);
	assert(w >= 3);
	const uint64_t mask = (1ull << (2ull*k)) - 1;
	const size_t s = k-w+1;
	const uint64_t smask = (1ull << (2ull*s)) - 1;
	uint64_t kmer = 0;
	uint64_t smer = 0;
	for (size_t i = 0; i < s; i++)
	{
		uint64_t c = getchar(i);
		smer <<= 2;
		smer += c;
	}
	kmer = smer;
	smerOrder.emplace_back(0, smer);
	for (size_t i = s; i < k; i++)
	{
		uint64_t c = getchar(i);
		kmer <<= 2;
		kmer += c;
		smer <<= 2;
		smer += c;
		smer &= smask;
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > smer) smerOrder.pop_back();
		smerOrder.emplace_back(i-s+1, smer);
	}
	if ((std::get<0>(smerOrder.front()) == 0) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == w-1))
	{
		callback(kmer, 0);
	}
	for (size_t i = k; i < maxLen; i++)
	{
		uint64_t c = getchar(i);
		kmer <<= 2;
		kmer += c;
		kmer &= mask;
		smer <<= 2;
		smer += c;
		smer &= smask;
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > smer) smerOrder.pop_back();
		while (smerOrder.size() > 0 && std::get<0>(smerOrder.front()) <= i-s+1-w) smerOrder.erase(smerOrder.begin());
		smerOrder.emplace_back(i-s+1, smer);
		if ((std::get<0>(smerOrder.front()) == i-s+2-w) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == i-s+1))
		{
			callback(kmer, i-k+1);
		}
	}
	smerOrder.clear();
}

template <typename F>
void iterateSyncmers(const std::vector<TwobitString>& readSequences, const size_t k, const size_t w, const size_t read, const size_t readStart, const size_t readEnd, const bool fw, F callback)
{
	assert(readStart < readEnd);
	assert(readEnd <= readSequences[read].size());
	if (k <= 31)
	{
		if (fw)
		{
			iterateSyncmers(readSequences, k, w, readEnd-readStart, callback, [&readSequences, read, readStart, readEnd](size_t index){ return readSequences[read].get(readStart+index); });
		}
		else
		{
			iterateSyncmers(readSequences, k, w, readEnd-readStart, callback, [&readSequences, read, readStart, readEnd](size_t index){ return 3-readSequences[read].get(readSequences[read].size() - 1 - (readStart+index)); });
		}
	}
	else
	{
		if (fw)
		{
			iterateSyncmersBigK(readSequences, k, w, readEnd-readStart, callback, [&readSequences, read, readStart, readEnd](size_t index){ return readSequences[read].get(readStart+index); });
		}
		else
		{
			iterateSyncmersBigK(readSequences, k, w, readEnd-readStart, callback, [&readSequences, read, readStart, readEnd](size_t index){ return 3-readSequences[read].get(readSequences[read].size() - 1 - (readStart+index)); });
		}
	}
}

void addKmer(phmap::flat_hash_map<uint64_t, uint64_t>& firstKmerPositionInLeft, phmap::flat_hash_map<uint64_t, std::vector<uint16_t>>& extraKmerPositionsInLeft, uint64_t kmer, size_t pos)
{
	if (firstKmerPositionInLeft.count(kmer) == 0)
	{
		firstKmerPositionInLeft[kmer] = 0xFFFFFFFFFFFF0000ull + pos;
	}
	else
	{
		auto& val = firstKmerPositionInLeft[kmer];
		if ((val & 0x00000000FFFF0000ull) == 0x00000000FFFF0000ull)
		{
			val &= 0xFFFFFFFF0000FFFFull;
			val += (uint64_t)pos << 16ull;
		}
		else if ((val & 0x0000FFFF00000000ull) == 0x0000FFFF00000000ull)
		{
			val &= 0xFFFF0000FFFFFFFFull;
			val += (uint64_t)pos << 32ull;
		}
		else if ((val & 0xFFFF000000000000ull) == 0xFFFF000000000000ull)
		{
			val &= 0x0000FFFFFFFFFFFFull;
			val += (uint64_t)pos << 48ull;
		}
		else
		{
			extraKmerPositionsInLeft[kmer].push_back(pos);
		}
	}
}

void heapify(std::vector<uint16_t>& vec)
{
	size_t pos = 0;
	while (pos*2+2 < vec.size())
	{
		if (vec[pos] < vec[pos*2+2] && vec[pos] < vec[pos*2+1]) return;
		if (vec[pos*2+2] < vec[pos*2+1])
		{
			std::swap(vec[pos], vec[pos*2+2]);
			pos = pos*2+2;
		}
		else
		{
			assert(vec[pos*2+1] < vec[pos*2+2]);
			std::swap(vec[pos], vec[pos*2+1]);
			pos = pos*2+1;
		}
	}
	if (pos*2+1 < vec.size())
	{
		if (vec[pos*2+1] < vec[pos])
		{
			std::swap(vec[pos], vec[pos*2+1]);
		}
	}
}

void removeKmer(phmap::flat_hash_map<uint64_t, uint64_t>& firstKmerPositionInLeft, phmap::flat_hash_map<uint64_t, std::vector<uint16_t>>& extraKmerPositionsInLeft, uint64_t kmer, size_t pos)
{
	auto found = firstKmerPositionInLeft.find(kmer);
	assert(found != firstKmerPositionInLeft.end());
	bool removed = false;
	if ((found->second & 0x000000000000FFFFull) == pos)
	{
		found->second >>= 16;
		found->second |= 0xFFFF000000000000ull;
		if (found->second == 0xFFFFFFFFFFFFFFFFull)
		{
			assert(extraKmerPositionsInLeft.count(kmer) == 0);
			firstKmerPositionInLeft.erase(found);
			return;
		}
		removed = true;
	}
	auto found2 = extraKmerPositionsInLeft.find(kmer);
	if (found2 == extraKmerPositionsInLeft.end()) return;
	assert(found2->second.size() >= 1);
	if (removed) found->second = (found->second & 0x0000FFFFFFFFFFFFull) + ((uint64_t)found2->second[0] << 48ull);
	std::swap(found2->second[0], found2->second.back());
	found2->second.pop_back();
	if (found2->second.size() == 0)
	{
		extraKmerPositionsInLeft.erase(found2);
	}
	else
	{
		heapify(found2->second);
	}
}

void removeHashCollisionMatches(const std::vector<TwobitString>& readSequences, MatchGroup& mappingMatch, const size_t k)
{
	for (size_t matchi = mappingMatch.matches.size()-1; matchi < mappingMatch.matches.size(); matchi--)
	{
		for (size_t i = mappingMatch.matches[matchi].length-1+k-1; i < mappingMatch.matches[matchi].length+k; i--)
		{
			uint8_t leftChar = readSequences[mappingMatch.leftRead].get(mappingMatch.leftStart+(size_t)mappingMatch.matches[matchi].leftStart+i);
			uint8_t rightChar;
			if (mappingMatch.rightFw)
			{
				rightChar = readSequences[mappingMatch.rightRead].get(mappingMatch.rightStart+(size_t)mappingMatch.matches[matchi].rightStart+i);
			}
			else
			{
				rightChar = 3-readSequences[mappingMatch.rightRead].get(readSequences[mappingMatch.rightRead].size()-1-(mappingMatch.rightStart+(size_t)mappingMatch.matches[matchi].rightStart+i));
			}
			if (leftChar == rightChar) continue;
			// std::cerr << "chop match leftread " << mappingMatch.leftRead << " rightread " << mappingMatch.rightRead << " " << (mappingMatch.rightFw ? "fw" : "bw") << " left " << mappingMatch.leftStart+(size_t)mappingMatch.matches[matchi].leftStart+i << " (" << (int)leftChar << ") right " << mappingMatch.rightStart+(size_t)mappingMatch.matches[matchi].rightStart+i << " (" << (int)rightChar << ")" << std::endl;
			if (i+1 < mappingMatch.matches[matchi].length)
			{
				mappingMatch.matches.emplace_back();
				mappingMatch.matches.back().leftStart = mappingMatch.matches[matchi].leftStart + i + 1;
				mappingMatch.matches.back().rightStart = mappingMatch.matches[matchi].rightStart + i + 1;
				mappingMatch.matches.back().length = mappingMatch.matches[matchi].length - (i + 1);
			}
			if (i >= k)
			{
				mappingMatch.matches[matchi].length = i-(k-1);
			}
			else
			{
				std::swap(mappingMatch.matches[matchi], mappingMatch.matches.back());
				mappingMatch.matches.pop_back();
				break;
			}
		}
	}
}

void getKmerMatches(const std::vector<TwobitString>& readSequences, MatchGroup& mappingMatch, const size_t k, const size_t d)
{
	thread_local size_t lastLeftRead = std::numeric_limits<size_t>::max();
	thread_local size_t lastLeftStart = std::numeric_limits<size_t>::max();
	thread_local size_t lastLeftEnd = std::numeric_limits<size_t>::max();
	thread_local std::vector<std::pair<uint64_t, size_t>> leftSyncmers;
	assert(mappingMatch.leftEnd - mappingMatch.leftStart < std::numeric_limits<uint16_t>::max());
	assert(mappingMatch.rightEnd - mappingMatch.rightStart < std::numeric_limits<uint16_t>::max());
	const size_t syncmerw = k > 11 ? 7 : (k/2);
	assert(k % 2 == 1);
	size_t leftStart = mappingMatch.leftStart;
	size_t leftEnd = mappingMatch.leftEnd;
	size_t rightStart = mappingMatch.rightStart;
	size_t rightEnd = mappingMatch.rightEnd;
	bool rightFw = mappingMatch.rightFw;
	size_t left = mappingMatch.leftRead;
	size_t right = mappingMatch.rightRead;
	phmap::flat_hash_map<uint64_t, uint64_t> firstKmerPositionInLeft;
	phmap::flat_hash_map<uint64_t, std::vector<uint16_t>> extraKmerPositionsInLeft;
	if (left != lastLeftRead || leftStart != lastLeftStart || leftEnd > lastLeftEnd)
	{
		leftSyncmers.clear();
		iterateSyncmers(readSequences, k, syncmerw, left, leftStart, leftEnd, true, [&leftSyncmers](const size_t kmer, const size_t pos)
		{
			leftSyncmers.emplace_back(kmer, pos);
		});
		lastLeftRead = left;
		lastLeftStart = leftStart;
		lastLeftEnd = leftEnd;
	}
	std::vector<std::pair<size_t, size_t>> currentMatchesPerDiagonal;
	size_t diagonalCount = 2*d+1;
	size_t zeroDiagonal = d;
	if (leftEnd-leftStart > rightEnd-rightStart)
	{
		diagonalCount += (leftEnd-leftStart)-(rightEnd-rightStart);
		zeroDiagonal = d + (leftEnd-leftStart)-(rightEnd-rightStart);
	}
	if (rightEnd-rightStart > leftEnd-leftStart) diagonalCount += (rightEnd-rightStart)-(leftEnd-leftStart);
	currentMatchesPerDiagonal.resize(diagonalCount, std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	size_t leftLeadingIndex = 0;
	size_t leftTrailingIndex = 0;
	iterateSyncmers(readSequences, k, syncmerw, right, rightStart, rightEnd, rightFw, [&firstKmerPositionInLeft, &extraKmerPositionsInLeft, &currentMatchesPerDiagonal, diagonalCount, zeroDiagonal, rightFw, left, right, leftStart, rightStart, leftEnd, rightEnd, d, k, &mappingMatch, &leftTrailingIndex, &leftLeadingIndex, &leftSyncmers](const size_t kmer, const size_t rightPos)
	{
		size_t interpolatedLeftPos = (double)(rightPos) / (double)(rightEnd-rightStart) * (double)(leftEnd-leftStart);
		while (leftLeadingIndex < leftSyncmers.size() && leftSyncmers[leftLeadingIndex].second <= interpolatedLeftPos + d && leftSyncmers[leftLeadingIndex].second < leftEnd)
		{
			addKmer(firstKmerPositionInLeft, extraKmerPositionsInLeft, leftSyncmers[leftLeadingIndex].first, leftSyncmers[leftLeadingIndex].second);
			leftLeadingIndex += 1;
		}
		while (leftTrailingIndex < leftSyncmers.size() && leftSyncmers[leftTrailingIndex].second + d < interpolatedLeftPos && leftSyncmers[leftTrailingIndex].second < leftEnd)
		{
			removeKmer(firstKmerPositionInLeft, extraKmerPositionsInLeft, leftSyncmers[leftTrailingIndex].first, leftSyncmers[leftTrailingIndex].second);
			leftTrailingIndex += 1;
		}
		assert(rightPos + zeroDiagonal >= interpolatedLeftPos + d);
		assert(rightPos + zeroDiagonal + d >= interpolatedLeftPos);
		size_t minDiagonal = rightPos + zeroDiagonal - interpolatedLeftPos - d;
		size_t maxDiagonal = rightPos + zeroDiagonal - interpolatedLeftPos + d;
		assert(maxDiagonal <= diagonalCount);
		iterateKmerMatchPositions(kmer, firstKmerPositionInLeft, extraKmerPositionsInLeft, [zeroDiagonal, rightPos, minDiagonal, maxDiagonal, &currentMatchesPerDiagonal, &mappingMatch, rightFw, left, right, leftStart, rightStart, diagonalCount, k](const size_t leftPos)
		{
			if (leftPos > zeroDiagonal + rightPos) return;
			if (zeroDiagonal + rightPos - leftPos >= diagonalCount) return;
			size_t diagonal = zeroDiagonal + rightPos - leftPos;
			if (diagonal < minDiagonal || diagonal > maxDiagonal) return;
			if (currentMatchesPerDiagonal[diagonal].first != std::numeric_limits<size_t>::max() && currentMatchesPerDiagonal[diagonal].second + k > rightPos)
			{
				currentMatchesPerDiagonal[diagonal].second = rightPos+1;
				return;
			}
			if (currentMatchesPerDiagonal[diagonal].first != std::numeric_limits<size_t>::max())
			{
				assert(currentMatchesPerDiagonal[diagonal].second != std::numeric_limits<size_t>::max());
				assert(currentMatchesPerDiagonal[diagonal].second > currentMatchesPerDiagonal[diagonal].first);
				assert(currentMatchesPerDiagonal[diagonal].first + zeroDiagonal >= diagonal);
				size_t leftMatchStart = currentMatchesPerDiagonal[diagonal].first + zeroDiagonal - diagonal;
				size_t rightMatchStart = currentMatchesPerDiagonal[diagonal].first;
				size_t length = currentMatchesPerDiagonal[diagonal].second - currentMatchesPerDiagonal[diagonal].first;
				assert(leftMatchStart < std::numeric_limits<uint16_t>::max());
				assert(rightMatchStart < std::numeric_limits<uint16_t>::max());
				assert(length < std::numeric_limits<uint16_t>::max());
				mappingMatch.matches.emplace_back();
				mappingMatch.matches.back().leftStart = leftMatchStart;
				mappingMatch.matches.back().rightStart = rightMatchStart;
				mappingMatch.matches.back().length = length;
			}
			currentMatchesPerDiagonal[diagonal].first = rightPos;
			currentMatchesPerDiagonal[diagonal].second = rightPos+1;
		});
	});
	for (size_t diagonal = 0; diagonal < diagonalCount; diagonal++)
	{
		if (currentMatchesPerDiagonal[diagonal].first == std::numeric_limits<size_t>::max()) continue;
		assert(currentMatchesPerDiagonal[diagonal].second != std::numeric_limits<size_t>::max());
		assert(currentMatchesPerDiagonal[diagonal].second > currentMatchesPerDiagonal[diagonal].first);
		assert(currentMatchesPerDiagonal[diagonal].first + zeroDiagonal >= diagonal);
		size_t leftMatchStart = currentMatchesPerDiagonal[diagonal].first + zeroDiagonal - diagonal;
		size_t rightMatchStart = currentMatchesPerDiagonal[diagonal].first;
		size_t length = currentMatchesPerDiagonal[diagonal].second - currentMatchesPerDiagonal[diagonal].first;
		assert(leftMatchStart < std::numeric_limits<uint16_t>::max());
		assert(rightMatchStart < std::numeric_limits<uint16_t>::max());
		assert(length < std::numeric_limits<uint16_t>::max());
		mappingMatch.matches.emplace_back();
		mappingMatch.matches.back().leftStart = leftMatchStart;
		mappingMatch.matches.back().rightStart = rightMatchStart;
		mappingMatch.matches.back().length = length;
	}
}

bool matchContained(MatchGroup::Match smaller, MatchGroup::Match bigger)
{
	if (smaller.leftStart <= bigger.leftStart) return false;
	if (smaller.rightStart <= bigger.rightStart) return false;
	if (smaller.leftStart + smaller.length >= bigger.leftStart + bigger.length) return false;
	if (smaller.rightStart + smaller.length >= bigger.rightStart + bigger.length) return false;
	return true;
}

void removeContainedKmerMatches(MatchGroup& matches)
{
	std::sort(matches.matches.begin(), matches.matches.end(), [](auto left, auto right) { return left.length < right.length; });
	std::vector<bool> remove;
	remove.resize(matches.matches.size(), false);
	for (size_t i = 0; i < matches.matches.size(); i++)
	{
		for (size_t j = matches.matches.size()-1; j > i; j--)
		{
			if (matches.matches[j].length < matches.matches[i].length+2) break;
			if (!matchContained(matches.matches[i], matches.matches[j])) continue;
			remove[i] = true;
			break;
		}
	}
	for (size_t i = matches.matches.size()-1; i < matches.matches.size(); i--)
	{
		if (!remove[i]) continue;
		std::swap(matches.matches[i], matches.matches.back());
		matches.matches.pop_back();
	}
}

std::vector<MatchGroup> addKmerMatches(const size_t numThreads, const std::vector<TwobitString>& readSequences, const std::vector<MatchGroup>& matches, const size_t graphk, const size_t graphd)
{
	const size_t SVLengthThreshold = 10000 + graphk;
	std::atomic<size_t> kmerMatchCount;
	kmerMatchCount = 0;
	size_t nextIndex = 0;
	std::mutex indexMutex;
	std::mutex resultMutex;
	std::vector<std::thread> threads;
	std::vector<MatchGroup> result;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&matches, &readSequences, graphk, graphd, &nextIndex, &indexMutex, &resultMutex, &kmerMatchCount, &result, SVLengthThreshold](){
			while (true)
			{
				size_t startIndex = 0;
				size_t endIndex = 0;
				{
					std::lock_guard<std::mutex> lock { indexMutex };
					startIndex = nextIndex;
					nextIndex += 1;
					while (nextIndex < matches.size() && matches[nextIndex].leftRead == matches[nextIndex-1].leftRead && matches[nextIndex].leftStart == matches[nextIndex-1].leftStart) nextIndex += 1;
					endIndex = nextIndex;
				}
				if (startIndex >= matches.size()) break;
				for (size_t i = startIndex; i < endIndex; i++)
				{
					std::vector<MatchGroup> tmp;
					size_t leftstart = matches[i].leftStart;
					size_t leftend = matches[i].leftEnd;
					size_t rightstart = matches[i].rightStart;
					size_t rightend = matches[i].rightEnd;
					size_t left = matches[i].leftRead;
					size_t right = matches[i].rightRead;
					bool rightFw = matches[i].rightFw;
					size_t numAlnChunks = std::max(leftend+1-leftstart, rightend+1-rightstart)/(std::numeric_limits<uint16_t>::max()-graphk-graphd-100)+1;
					assert(numAlnChunks >= 1);
					if (numAlnChunks == 1)
					{
						tmp.emplace_back();
						tmp.back().leftRead = left;
						tmp.back().rightRead = right;
						tmp.back().rightFw = rightFw;
						tmp.back().leftStart = leftstart;
						tmp.back().leftEnd = leftend+1;
						tmp.back().rightStart = rightstart;
						tmp.back().rightEnd = rightend+1;
					}
					else
					{
						assert(numAlnChunks >= 2);
						size_t leftPerChunk = (double)(leftend+1-leftstart)/(double)numAlnChunks;
						size_t rightPerChunk = (double)(rightend+1-rightstart)/(double)numAlnChunks;
						for (size_t l = 0; l < numAlnChunks; l++)
						{
							tmp.emplace_back();
							tmp.back().leftRead = left;
							tmp.back().rightRead = right;
							tmp.back().rightFw = rightFw;
							tmp.back().leftStart = leftstart + l * leftPerChunk;
							tmp.back().leftEnd = leftstart + (l+1) * leftPerChunk + graphk + graphd;
							tmp.back().rightStart = rightstart + l * rightPerChunk;
							tmp.back().rightEnd = rightstart + (l+1) * rightPerChunk + graphk + graphd;
						}
						tmp.back().leftEnd = leftend+1;
						tmp.back().rightEnd = rightend+1;
					}
					bool valid = true;
					size_t lastLeftEnd = leftstart;
					size_t lastRightEnd = rightstart;
					for (size_t j = 0; j < tmp.size(); j++)
					{
						getKmerMatches(readSequences, tmp[j], graphk, graphd);
						if (graphk > 31) removeHashCollisionMatches(readSequences, tmp[j], graphk);
						removeContainedKmerMatches(tmp[j]);
						std::sort(tmp[j].matches.begin(), tmp[j].matches.end(), [](auto left, auto right) { return left.leftStart < right.leftStart; });
						for (auto match : tmp[j].matches)
						{
							if (tmp[j].leftStart + (size_t)match.leftStart > lastLeftEnd + SVLengthThreshold)
							{
								valid = false;
								break;
							}
							if (tmp[j].rightStart + (size_t)match.rightStart > lastRightEnd + SVLengthThreshold)
							{
								valid = false;
								break;
							}
							lastLeftEnd = tmp[j].leftStart + (size_t)match.leftStart + (size_t)match.length;
							lastRightEnd = tmp[j].rightStart + (size_t)match.rightStart + (size_t)match.length;
						}
						if (!valid) break;
					}
					if (lastLeftEnd+SVLengthThreshold < leftend) valid = false;
					if (lastRightEnd+SVLengthThreshold < rightend) valid = false;
					if (valid)
					{
						std::lock_guard<std::mutex> lock { resultMutex };
						kmerMatchCount += matches[i].matches.size();
						result.insert(result.end(), tmp.begin(), tmp.end());
						for (size_t j = 0; j < tmp.size(); j++)
						{
							kmerMatchCount += tmp[j].matches.size();
						}
					}
				}
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	std::cerr << kmerMatchCount << " kmer matches" << std::endl;
	return result;
}
