#include <iostream>
#include <cassert>
#include <thread>
#include <chrono>
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

void extendBreakpoints(std::vector<RankBitvector>& breakpoints, const std::vector<size_t>& readOrder, const std::vector<std::vector<size_t>>& fwMatchChunks, const std::vector<std::vector<size_t>>& bwMatchChunks, const std::vector<size_t>& leftoverMatchesFw, const std::vector<size_t>& leftoverMatchesBw, const size_t chunkSize, const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t numThreads)
{
	std::vector<bool> shouldDoChunk;
	shouldDoChunk.resize(fwMatchChunks.size(), true);
	std::vector<std::thread> threads;
	std::atomic<bool> allDone;
	std::atomic<size_t> threadsComputingNow;
	threadsComputingNow = 0;
	allDone = false;
	std::mutex chunkPickerMutex;
	std::mutex leftoverMutex;
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([thread, &threads, &allDone, &fwMatchChunks, &bwMatchChunks, &shouldDoChunk, &chunkPickerMutex, &leftoverMutex, &threadsComputingNow, &matches, &breakpoints, &readLengths, numThreads, chunkSize, &readOrder, &leftoverMatchesBw, &leftoverMatchesFw]()
		{
			while (true)
			{
				size_t i = std::numeric_limits<size_t>::max();
				{
					std::lock_guard<std::mutex> lock { chunkPickerMutex };
					for (size_t chunk = 0; chunk < fwMatchChunks.size(); chunk++)
					{
						if (!shouldDoChunk[chunk]) continue;
						shouldDoChunk[chunk] = false;
						i = chunk;
						break;
					}
				}
				if (i == std::numeric_limits<size_t>::max())
				{
					if (allDone) return;
					if (leftoverMutex.try_lock())
					{
						while (threadsComputingNow > 0)
						{
							std::this_thread::sleep_for(std::chrono::milliseconds(10));
						}
						std::lock_guard<std::mutex> lock { chunkPickerMutex };
						assert(threadsComputingNow == 0);
						bool changed = false;
						for (size_t groupi : leftoverMatchesFw)
						{
							for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
							{
								std::pair<bool, bool> addedAny;
								addedAny = extendBreakpointsFwFw(readLengths, breakpoints, matches[groupi].leftRead, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
								if (addedAny.first || addedAny.second)
								{
									changed = true;
									shouldDoChunk[readOrder[matches[groupi].leftRead] / chunkSize] = true;
									shouldDoChunk[readOrder[matches[groupi].rightRead] / chunkSize] = true;
								}
							}
						}
						for (size_t groupi : leftoverMatchesBw)
						{
							for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
							{
								std::pair<bool, bool> addedAny;
								addedAny = extendBreakpointsFwBw(readLengths, breakpoints, matches[groupi].leftRead, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
								if (addedAny.first || addedAny.second)
								{
									changed = true;
									shouldDoChunk[readOrder[matches[groupi].leftRead] / chunkSize] = true;
									shouldDoChunk[readOrder[matches[groupi].rightRead] / chunkSize] = true;
								}
							}
						}
						if (!changed) allDone = true;
						leftoverMutex.unlock();
					}
					if (allDone) return;
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
					continue;
				}
				assert(threadsComputingNow < numThreads);
				threadsComputingNow += 1;
				while (true)
				{
					while (true)
					{
						bool fwChanged = false;
						for (size_t groupi : fwMatchChunks[i])
						{
							for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
							{
								std::pair<bool, bool> addedAny;
								addedAny = extendBreakpointsFwFw(readLengths, breakpoints, matches[groupi].leftRead, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
								if (addedAny.first || addedAny.second) fwChanged = true;
							}
						}
						if (!fwChanged) break;
					}
					bool bwChanged = false;
					for (size_t groupi : bwMatchChunks[i])
					{
						for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
						{
							std::pair<bool, bool> addedAny;
							addedAny = extendBreakpointsFwBw(readLengths, breakpoints, matches[groupi].leftRead, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
							if (addedAny.first || addedAny.second) bwChanged = true;
						}
					}
					if (!bwChanged) break;
				}
				assert(threadsComputingNow >= 1);
				threadsComputingNow -= 1;
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
}

bool isPalindrome(const TwobitString& string, const size_t start, const size_t end, const size_t k)
{
	assert((end - start) % 2 == 0);
	for (size_t i = 0; i < (end - start + k)/2; i++)
	{
		uint8_t fwchar = string.get(start+i);
		uint8_t bwchar = string.get(end+k-2-i);
		if (bwchar != 3-fwchar) return false;
	}
	return true;
}

bool fixPalindromeBreakpoints(std::vector<RankBitvector>& breakpoints, const std::vector<TwobitString>& readSequences, const size_t k)
{
	bool fixedAny = false;
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		size_t chunkStart = 0;
		for (size_t j = 0; j < breakpoints[i].size(); j++)
		{
			if (!breakpoints[i].get(j)) continue;
			if (j == chunkStart) continue;
			if ((j - chunkStart) % 2 == 0)
			{
				if (isPalindrome(readSequences[i], chunkStart, j, k))
				{
					assert((chunkStart - j)/2 >= 1);
					assert(!breakpoints[i].get(chunkStart + (j - chunkStart)/2));
					breakpoints[i].set(chunkStart + (j - chunkStart)/2, true);
					assert(!breakpoints[i].get(chunkStart + (j - chunkStart)/2 - 1) || j-chunkStart == 2);
					breakpoints[i].set(chunkStart + (j - chunkStart)/2 - 1, true);
					fixedAny = true;
				}
			}
			chunkStart = j;
		}
	}
	return fixedAny;
}

std::vector<RankBitvector> extendBreakpoints(const std::vector<TwobitString>& readSequences, const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t numThreads, const size_t k)
{
	size_t numChunks = 25;
	std::vector<RankBitvector> breakpoints;
	breakpoints.resize(readLengths.size());
	for (size_t i = 0; i < readLengths.size(); i++)
	{
		breakpoints[i].resize(readLengths[i]+1);
		breakpoints[i].set(0, true);
		breakpoints[i].set(readLengths[i], true);
	}
	std::vector<phmap::flat_hash_set<size_t>> readHasMatch;
	readHasMatch.resize(readLengths.size());
	for (size_t groupi = 0; groupi < matches.size(); groupi++)
	{
		readHasMatch[matches[groupi].leftRead].emplace(matches[groupi].rightRead);
		readHasMatch[matches[groupi].rightRead].emplace(matches[groupi].leftRead);
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
	std::vector<size_t> readOrder;
	readOrder.resize(readLengths.size(), std::numeric_limits<size_t>::max());
	std::vector<size_t> stackone;
	std::vector<size_t> stacktwo;
	size_t nextOrder = 0;
	for (size_t i = 0; i < readOrder.size(); i++)
	{
		if (readOrder[i] != std::numeric_limits<size_t>::max()) continue;
		if (readHasMatch[i].size() == 0) continue;
		stackone.emplace_back(i);
		while (stackone.size() > 0 || stacktwo.size() > 0)
		{
			if (stackone.size() == 0)
			{
				while (stacktwo.size() > 0)
				{
					stackone.emplace_back(stacktwo.back());
					stacktwo.pop_back();
				}
			}
			assert(stackone.size() >= 1);
			auto top = stackone.back();
			stackone.pop_back();
			if (readOrder[top] != std::numeric_limits<size_t>::max()) continue;
			readOrder[top] = nextOrder;
			nextOrder += 1;
			for (auto read : readHasMatch[top])
			{
				if (readOrder[read] != std::numeric_limits<size_t>::max()) continue;
				stacktwo.emplace_back(read);
			}
		}
	}
	assert(nextOrder <= readOrder.size());
	size_t chunkSize = nextOrder / numChunks + 1;
	std::vector<std::vector<size_t>> fwMatchChunks;
	fwMatchChunks.resize(numChunks);
	std::vector<std::vector<size_t>> bwMatchChunks;
	bwMatchChunks.resize(numChunks);
	std::vector<size_t> leftoverMatchesFw;
	std::vector<size_t> leftoverMatchesBw;
	for (size_t i = 0; i < matches.size(); i++)
	{
		if (readOrder[matches[i].leftRead] / chunkSize == readOrder[matches[i].rightRead] / chunkSize)
		{
			size_t chunk = readOrder[matches[i].leftRead] / chunkSize;
			if (matches[i].rightFw)
			{
				fwMatchChunks[chunk].emplace_back(i);
			}
			else
			{
				bwMatchChunks[chunk].emplace_back(i);
			}
		}
		else
		{
			if (matches[i].rightFw)
			{
				leftoverMatchesFw.emplace_back(i);
			}
			else
			{
				leftoverMatchesBw.emplace_back(i);
			}
		}
	}
	extendBreakpoints(breakpoints, readOrder, fwMatchChunks, bwMatchChunks, leftoverMatchesFw, leftoverMatchesBw, chunkSize, readLengths, matches, numThreads);
	bool fixedAny = fixPalindromeBreakpoints(breakpoints, readSequences, k);
	if (fixedAny)
	{
		extendBreakpoints(breakpoints, readOrder, fwMatchChunks, bwMatchChunks, leftoverMatchesFw, leftoverMatchesBw, chunkSize, readLengths, matches, numThreads);
	}
	for (size_t i = 0; i < breakpoints.size(); i++) breakpoints[i].buildRanks();
	return breakpoints;
}

uint64_t find(std::vector<std::vector<uint64_t>>& result, const std::vector<size_t>& countBeforeRead, const RankBitvector& segmentToRead, const size_t pos)
{
	assert((pos & firstBitUint64_t) == 0);
	size_t readi = segmentToRead.getRank(pos);
	size_t offset = pos - countBeforeRead[readi];
	if (result[readi][offset] == (pos & maskUint64_t))
	{
		assert(result[readi][offset] == pos);
		return pos;
	}
	uint64_t foundPos = result[readi][offset];
	while (true)
	{
		size_t nextreadi = segmentToRead.getRank(foundPos & maskUint64_t);
		size_t nextoffset = (foundPos & maskUint64_t) - countBeforeRead[nextreadi];
		if (result[nextreadi][nextoffset] == (foundPos & maskUint64_t)) break;
		uint64_t nextFoundPos = result[nextreadi][nextoffset];
		if (foundPos & firstBitUint64_t) nextFoundPos ^= firstBitUint64_t;
		foundPos = nextFoundPos;
	}
	uint64_t finalPos = foundPos;
	foundPos = result[readi][offset];
	result[readi][offset] = finalPos;
	while (true)
	{
		size_t nextreadi = segmentToRead.getRank(foundPos & maskUint64_t);
		size_t nextoffset = (foundPos & maskUint64_t) - countBeforeRead[nextreadi];
		if (result[nextreadi][nextoffset] == (foundPos & maskUint64_t)) break;
		uint64_t nextFoundPos = result[nextreadi][nextoffset];
		if (foundPos & firstBitUint64_t) finalPos ^= firstBitUint64_t;
		result[nextreadi][nextoffset] = finalPos;
		foundPos = nextFoundPos;
	}
	return result[readi][offset];
}

void mergeSegments(std::vector<std::vector<uint64_t>>& result, const std::vector<size_t>& countBeforeRead, const RankBitvector& segmentToRead, const size_t leftPos, const size_t rightPos, const bool fw)
{
	auto leftParent = find(result, countBeforeRead, segmentToRead, leftPos);
	auto rightParent = find(result, countBeforeRead, segmentToRead, rightPos);
	size_t leftreadi = segmentToRead.getRank(leftParent & maskUint64_t);
	size_t leftoffset = (leftParent & maskUint64_t) - countBeforeRead[leftreadi];
	size_t rightreadi = segmentToRead.getRank(rightParent & maskUint64_t);
	size_t rightoffset = (rightParent & maskUint64_t) - countBeforeRead[rightreadi];
	assert(result[leftreadi][leftoffset] == (leftParent & maskUint64_t));
	assert(result[rightreadi][rightoffset] == (rightParent & maskUint64_t));
	if ((leftParent & maskUint64_t) == (rightParent & maskUint64_t))
	{
		assert(((leftParent & firstBitUint64_t) == 0) ^ ((rightParent & firstBitUint64_t) == 0) ^ fw);
		return;
	}
	result[rightreadi][rightoffset] = (leftParent & maskUint64_t) + ((((leftParent & firstBitUint64_t) == 0) ^ ((rightParent & firstBitUint64_t) == 0) ^ fw) ? 0 : firstBitUint64_t);
}

void mergeSegments(const std::vector<size_t>& readLengths, const std::vector<size_t>& countBeforeRead, const RankBitvector& segmentToRead, std::vector<std::vector<uint64_t>>& result, const std::vector<RankBitvector>& breakpoints, const size_t leftRead, const bool leftFw, const size_t leftStart, const size_t leftEnd, const size_t rightRead, const bool rightFw, const size_t rightStart, const size_t rightEnd)
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
			mergeSegments(result, countBeforeRead, segmentToRead, countBeforeRead[leftRead] + leftFirst + i, countBeforeRead[rightRead] + rightFirst+i, true);
		}
		else
		{
			mergeSegments(result, countBeforeRead, segmentToRead, countBeforeRead[leftRead] + leftFirst + i, countBeforeRead[rightRead] + rightLast-i-1, false);
		}
	}
}

std::pair<std::vector<std::vector<uint64_t>>, std::vector<size_t>> mergeSegments(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const std::vector<RankBitvector>& breakpoints, const size_t countBreakpoints)
{
	assert(breakpoints.size() == readLengths.size());
	std::vector<std::vector<uint64_t>> result;
	assert(countBreakpoints >= 2*readLengths.size());
	result.resize(breakpoints.size());
	std::vector<size_t> countBeforeRead;
	countBeforeRead.resize(readLengths.size(), 0);
	RankBitvector segmentToRead;
	segmentToRead.resize(countBreakpoints - readLengths.size());
	size_t countBeforeNow = 0;
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		countBeforeRead[i] = countBeforeNow;
		size_t count = 0;
		for (size_t j = 1; j < breakpoints[i].size(); j++)
		{
			if (!breakpoints[i].get(j)) continue;
			count += 1;
		}
		result[i].resize(count);
		for (size_t j = 0; j < count; j++)
		{
			result[i][j] = countBeforeNow + j;
		}
		assert(count >= 1);
		segmentToRead.set(countBeforeNow + count - 1, true);
		countBeforeNow += count;
	}
	segmentToRead.buildRanks();
	assert(countBeforeNow == segmentToRead.size());
	for (size_t groupi = 0; groupi < matches.size(); groupi++)
	{
		for (size_t posi = 0; posi < matches[groupi].matches.size(); posi++)
		{
			mergeSegments(readLengths, countBeforeRead, segmentToRead, result, breakpoints, matches[groupi].leftRead, true, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart, matches[groupi].leftStart + matches[groupi].matches[posi].leftStart+matches[groupi].matches[posi].length, matches[groupi].rightRead, matches[groupi].rightFw, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart, matches[groupi].rightStart + matches[groupi].matches[posi].rightStart+matches[groupi].matches[posi].length);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			find(result, countBeforeRead, segmentToRead, countBeforeRead[i]+j);
		}
	}
	return std::make_pair(std::move(result), std::move(countBeforeRead));
}

RankBitvector getSegmentToNode(const std::vector<std::vector<uint64_t>>& segments, const std::vector<uint64_t>& segmentCountBeforeRead)
{
	RankBitvector result;
	result.resize(segmentCountBeforeRead.back() + segments.back().size());
	for (size_t i = 0; i < segments.size(); i++)
	{
		for (size_t j = 0; j < segments[i].size(); j++)
		{
			assert((segments[i][j] & maskUint64_t) < result.size());
			if ((segments[i][j] & maskUint64_t) != segmentCountBeforeRead[i] + j) continue;
			assert(segments[i][j] == segmentCountBeforeRead[i] + j);
			result.set(segments[i][j], true);
		}
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

RankBitvector keepCoveredNodes(const std::vector<std::vector<uint64_t>>& segments, const RankBitvector& segmentToNode, const size_t countNodes, const size_t minCoverage)
{
	std::vector<size_t> coverages;
	coverages.resize(countNodes, 0);
	for (size_t i = 0; i < segments.size(); i++)
	{
		for (size_t j = 0; j < segments[i].size(); j++)
		{
			uint64_t node = segmentToNode.getRank(segments[i][j] & maskUint64_t);
			assert(node < coverages.size());
			coverages[node] += 1;
		}
	}
	RankBitvector result;
	result.resize(countNodes);
	for (size_t i = 0; i < coverages.size(); i++)
	{
		assert(coverages[i] != 0);
		if (coverages[i] >= minCoverage) result.set(i, true);
	}
	result.buildRanks();
	return result;
}

std::vector<size_t> getNodeLengths(const std::vector<std::vector<uint64_t>>& segments, const RankBitvector& segmentToNode, const RankBitvector& keptNodes, const std::vector<RankBitvector>& breakpoints, const size_t countKeptNodes)
{
	std::vector<size_t> result;
	result.resize(countKeptNodes, std::numeric_limits<size_t>::max());
	for (size_t readi = 0; readi < breakpoints.size(); readi++)
	{
		size_t lastBreakpoint = 0;
		assert(breakpoints[readi].get(0));
		size_t i = 0;
		for (size_t pos = 1; pos < breakpoints[readi].size(); pos++)
		{
			if (!breakpoints[readi].get(pos)) continue;
			size_t node = segmentToNode.getRank(segments[readi][i] & maskUint64_t);
			if (!keptNodes.get(node))
			{
				i += 1;
				lastBreakpoint = pos;
				continue;
			}
			const size_t length = pos - lastBreakpoint;
			node = keptNodes.getRank(node);
			assert(node < result.size());
			assert(length >= 1);
			assert(result[node] == std::numeric_limits<size_t>::max() || result[node] == length);
			result[node] = length;
			i += 1;
			lastBreakpoint = pos;
		}
		assert(i == segments[readi].size());
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i] != std::numeric_limits<size_t>::max());
	}
	return result;
}

// destroy segments to save memory, otherwise segments and result both use huge memory to represent the same thing
std::vector<ReadPathBundle> getReadPathsAndDestroySegments(std::vector<std::vector<uint64_t>>& segments, const std::vector<RankBitvector>& breakpoints, const RankBitvector& segmentToNode, const RankBitvector& keptNodes, const std::vector<size_t>& readLengths)
{
	std::vector<ReadPathBundle> result;
	result.resize(breakpoints.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i].readName = i;
		result[i].readLength = readLengths[i];
		uint64_t lastNode = std::numeric_limits<size_t>::max();
		size_t segmenti = 0;
		for (size_t j = 0; j < breakpoints[i].size()-1; j++)
		{
			if (!breakpoints[i].get(j)) continue;
			size_t node = segmentToNode.getRank(segments[i][segmenti] & maskUint64_t);
			bool fw = (segments[i][segmenti] & firstBitUint64_t) == 0;
			segmenti += 1;
			if (!keptNodes.get(node))
			{
				lastNode = std::numeric_limits<size_t>::max();
				continue;
			}
			node = keptNodes.getRank(node);
			if (lastNode == std::numeric_limits<size_t>::max())
			{
				result[i].paths.emplace_back();
				result[i].paths.back().pathLeftClipKmers = 0;
				result[i].paths.back().pathRightClipKmers = 0;
				result[i].paths.back().readStartPos = j;
			}
			result[i].paths.back().path.emplace_back(node + (fw ? firstBitUint64_t : 0));
			lastNode = node + (fw ? firstBitUint64_t : 0);
		}
		assert(segmenti == segments[i].size());
		std::vector<uint64_t> tmp;
		std::swap(segments[i], tmp);
	}
	return result;
}

std::pair<KmerGraph, std::vector<ReadPathBundle>> makeKmerGraph(const std::vector<TwobitString>& readSequences, const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t minCoverage, const size_t numThreads, const size_t graphk)
{
	KmerGraph result;
	std::vector<RankBitvector> breakpoints = extendBreakpoints(readSequences, readLengths, matches, numThreads, graphk);
	size_t countBreakpoints = 0;
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		for (size_t j = 0; j < breakpoints[i].size(); j++)
		{
			if (breakpoints[i].get(j)) countBreakpoints += 1;
		}
	}
//	std::cerr << countBreakpoints << " breakpoints" << std::endl;
	std::vector<size_t> segmentCountBeforeRead;
	std::vector<std::vector<uint64_t>> segments;
	std::tie(segments, segmentCountBeforeRead) = mergeSegments(readLengths, matches, breakpoints, countBreakpoints);
//	std::cerr << segmentCountBeforeRead.back() + segments.back().size() << " segments" << std::endl;
	RankBitvector segmentToNode = getSegmentToNode(segments, segmentCountBeforeRead);
	{
		std::vector<size_t> tmp;
		std::swap(tmp, segmentCountBeforeRead);
	}
	size_t countNodes = (segmentToNode.getRank(segmentToNode.size()-1) + (segmentToNode.get(segmentToNode.size()-1) ? 1 : 0));
//	std::cerr << countNodes << " kmer-nodes pre coverage filter" << std::endl;
	RankBitvector keptNodes = keepCoveredNodes(segments, segmentToNode, countNodes, minCoverage);
	size_t countKeptNodes = (keptNodes.getRank(keptNodes.size()-1) + (keptNodes.get(keptNodes.size()-1) ? 1 : 0));
//	std::cerr << countKeptNodes << " kmer-nodes post coverage filter" << std::endl;
	result.lengths = getNodeLengths(segments, segmentToNode, keptNodes, breakpoints, countKeptNodes);
	std::vector<ReadPathBundle> readPaths = getReadPathsAndDestroySegments(segments, breakpoints, segmentToNode, keptNodes, readLengths);
	size_t countNodeMatches = 0;
	size_t countReadPaths = 0;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		countReadPaths += readPaths[i].paths.size();
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			countNodeMatches += readPaths[i].paths[j].path.size();
		}
	}
//	std::cerr << countReadPaths << " read paths" << std::endl;
//	std::cerr << countNodeMatches << " read-node matches" << std::endl;
	return std::make_pair(std::move(result), std::move(readPaths));
}

size_t KmerGraph::nodeCount() const
{
	return lengths.size();
}
