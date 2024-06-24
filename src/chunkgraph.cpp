#include <map>
#include <queue>
#include <mutex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "Common.h"
#include "RankBitvector.h"
#include "UnionFind.h"
#include "EdlibWrapper.h"
#include "CompressedStringIndex.h"
#include "TwobitString.h"
#include "ChunkUnitigGraph.h"
#include "ChunkGraphWriter.h"
#include "KmerIterator.h"
#include "ConsensusMaker.h"
#include "SequenceHelper.h"
#include "GraphCleaner.h"
#include "ChunkExtractor.h"

double mismatchFraction;

size_t popcount(uint64_t x);

size_t findBiggestSmallerIndex(const std::vector<size_t>& nums, const size_t value)
{
	if (nums.size() < 5)
	{
		assert(nums[0] < value);
		for (size_t i = 1; i < nums.size(); i++)
		{
			if (nums[i] > value) return i-1;
		}
		return nums.size()-1;
	}
	size_t low = 0;
	size_t high = nums.size();
	while (high > low)
	{
		size_t mid = (low + high) / 2;
		assert(nums[mid] != value);
		if (nums[mid] < value)
		{
			low = mid+1;
		}
		else if (nums[mid] > value)
		{
			high = mid;
		}
	}
	assert(low == high);
	assert(low <= nums.size());
	assert(low == nums.size() || nums[low] > value);
	assert(low > 0);
	assert(nums[low-1] < value);
	return low-1;
}

void splitPerLength(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by length" << std::endl;
	const double differenceFraction = 0.02;
	const size_t differenceConstant = 50;
	std::vector<std::vector<size_t>> lengthsPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			while ((std::get<2>(t) & maskUint64_t) >= lengthsPerChunk.size()) lengthsPerChunk.emplace_back();
			lengthsPerChunk[std::get<2>(t) & maskUint64_t].emplace_back(std::get<1>(t) - std::get<0>(t));
		}
	}
	std::vector<std::vector<size_t>> splitters;
	splitters.resize(lengthsPerChunk.size());
	iterateMultithreaded(0, lengthsPerChunk.size(), numThreads, [&lengthsPerChunk, &splitters, differenceFraction, differenceConstant](const size_t i)
	{
		std::sort(lengthsPerChunk[i].begin(), lengthsPerChunk[i].end());
		splitters[i].emplace_back(0);
		for (size_t j = 1; j < lengthsPerChunk[i].size(); j++)
		{
			size_t distance = lengthsPerChunk[i][j] - lengthsPerChunk[i][j-1];
			if (distance < std::max((size_t)(lengthsPerChunk[i][j] * differenceFraction), (size_t)differenceConstant)) continue;
			splitters[i].emplace_back((lengthsPerChunk[i][j] + lengthsPerChunk[i][j-1])/2);
		}
	});
	std::vector<std::vector<size_t>> splitterToNode;
	splitterToNode.resize(lengthsPerChunk.size());
	size_t nextNum = 0;
	for (size_t i = 0; i < splitters.size(); i++)
	{
		for (size_t j = 0; j < splitters[i].size(); j++)
		{
			splitterToNode[i].emplace_back(nextNum);
			nextNum += 1;
		}
	}
	iterateMultithreaded(0, chunksPerRead.size(), numThreads, [&splitters, &splitterToNode, &chunksPerRead](const size_t i)
	{
		for (auto& t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			assert(std::get<1>(t) > std::get<0>(t));
			size_t distance = std::get<1>(t) - std::get<0>(t);
			assert((std::get<2>(t) & maskUint64_t) < splitters.size());
			size_t index = findBiggestSmallerIndex(splitters[std::get<2>(t) & maskUint64_t], distance);
			assert(index < splitterToNode[std::get<2>(t) & maskUint64_t].size());
			std::get<2>(t) = (std::get<2>(t) & firstBitUint64_t) + splitterToNode[std::get<2>(t) & maskUint64_t][index];
		}
	});
}

void checkPhasablePair(const phmap::flat_hash_map<size_t, size_t>& bwForks, const phmap::flat_hash_map<size_t, size_t>& fwForks, std::vector<std::vector<size_t>>& phaseIdentities)
{
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> pairCoverage;
	phmap::flat_hash_set<size_t> foundOccurrences;
	phmap::flat_hash_set<size_t> foundBw;
	phmap::flat_hash_set<size_t> foundFw;
	for (auto bwfork : bwForks)
	{
		foundBw.insert(bwfork.second);
		foundOccurrences.insert(bwfork.first);
		if (fwForks.count(bwfork.first) == 0) continue;
		size_t bw = bwfork.second;
		size_t fw = fwForks.at(bwfork.first);
		pairCoverage[std::make_pair(bw, fw)] += 1;
	}
	for (auto fwfork : fwForks)
	{
		foundFw.insert(fwfork.second);
		foundOccurrences.insert(fwfork.first);
	}
	if (foundFw.size() != foundBw.size()) return;
	if (foundOccurrences.size() < phaseIdentities.size()) return;
	size_t coverageInPairs = 0;
	for (auto pair : pairCoverage)
	{
		coverageInPairs += pair.second;
		if (pair.second < 5) return;
	}
	if (coverageInPairs < phaseIdentities.size() * 0.95) return;
	phmap::flat_hash_map<size_t, size_t> parent;
	for (auto pair : pairCoverage)
	{
		assert((pair.first.first & firstBitUint64_t) == 0);
		assert((pair.first.second & firstBitUint64_t) == 0);
		if (parent.count(pair.first.first) == 0) parent[pair.first.first] = pair.first.first;
		if (parent.count(pair.first.second + firstBitUint64_t) == 0) parent[pair.first.second + firstBitUint64_t] = pair.first.second + firstBitUint64_t;
		merge(parent, pair.first.first, pair.first.second + firstBitUint64_t);
	}
	phmap::flat_hash_map<size_t, size_t> clusterCount;
	for (auto pair : pairCoverage)
	{
		auto cluster = find(parent, pair.first.first);
		clusterCount[cluster] += pair.second;
	}
	if (clusterCount.size() < 2) return;
	for (auto kmer : foundBw)
	{
		if (parent.count(kmer) == 0) return;
	}
	for (auto kmer : foundFw)
	{
		if (parent.count(kmer + firstBitUint64_t) == 0) return;
	}
	size_t inserted = 0;
	for (auto bwfork : bwForks)
	{
		if (fwForks.count(bwfork.first) == 1)
		{
			assert(find(parent, bwfork.second) == find(parent, fwForks.at(bwfork.first) + firstBitUint64_t));
		}
		phaseIdentities[bwfork.first].emplace_back(find(parent, bwfork.second));
		inserted += 1;
	}
	for (auto fwfork : fwForks)
	{
		if (bwForks.count(fwfork.first) == 1) continue;
		phaseIdentities[fwfork.first].emplace_back(find(parent, fwfork.second + firstBitUint64_t));
		inserted += 1;
	}
	assert(inserted == phaseIdentities.size());
}

std::vector<std::vector<std::vector<size_t>>> getAlleleOccurrences(const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>>& alleles, const size_t numOccurrences)
{
	std::vector<std::vector<std::vector<size_t>>> result;
	size_t clusterStart = 0;
	for (size_t i = 1; i <= alleles.size(); i++)
	{
		if (i < alleles.size() && std::get<0>(alleles[i]) == std::get<0>(alleles[i-1]) && std::get<1>(alleles[i]) == std::get<1>(alleles[i-1]) && std::get<2>(alleles[i]) == std::get<2>(alleles[i-1]) && std::get<3>(alleles[i]) == std::get<3>(alleles[i-1])) continue;
		if (i - clusterStart != numOccurrences)
		{
			clusterStart = i;
			continue;
		}
		std::vector<bool> foundOccurrences;
		foundOccurrences.resize(numOccurrences, false);
		bool valid = true;
		for (size_t j = clusterStart; j < i; j++)
		{
			assert(std::get<4>(alleles[j]) < foundOccurrences.size());
			if (foundOccurrences[std::get<4>(alleles[j])])
			{
				valid = false;
				break;
			}
			foundOccurrences[std::get<4>(alleles[j])] = true;
		}
		for (size_t j = 0; j < foundOccurrences.size(); j++)
		{
			if (!foundOccurrences[j])
			{
				valid = false;
				break;
			}
		}
		if (!valid)
		{
			clusterStart = i;
			continue;
		}
		phmap::flat_hash_map<size_t, size_t> alleleCounts;
		for (size_t j = clusterStart; j < i; j++)
		{
			alleleCounts[std::get<5>(alleles[j])] += 1;
		}
		if (alleleCounts.size() < 2)
		{
			clusterStart = i;
			continue;
		}
		for (auto pair : alleleCounts)
		{
			if (pair.second < 5)
			{
				valid = false;
				break;
			}
		}
		if (!valid)
		{
			clusterStart = i;
			continue;
		}
		result.emplace_back();
		phmap::flat_hash_map<size_t, size_t> alleleMap;
		for (auto pair : alleleCounts)
		{
			alleleMap[pair.first] = result.back().size();
//			std::cerr << "allele occurrence " << result.size()-1 << " allele " << result.back().size() << " is " << pair.first << std::endl;
			result.back().emplace_back();
		}
		for (size_t j = clusterStart; j < i; j++)
		{
			assert(alleleMap.count(std::get<5>(alleles[j])) == 1);
			assert(alleleMap.at(std::get<5>(alleles[j])) < result.back().size());
			result.back()[alleleMap.at(std::get<5>(alleles[j]))].emplace_back(std::get<4>(alleles[j]));
		}
		clusterStart = i;
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			std::sort(result[i][j].begin(), result[i][j].end());
		}
	}
	return result;
}

bool allelesMatchTwoVariantsThreeHaps(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right)
{
	if (left.size() != 2) return false;
	if (right.size() != 2) return false;
	if (left[0].size() < 10) return false;
	if (left[1].size() < 10) return false;
	if (right[0].size() < 10) return false;
	if (right[1].size() < 10) return false;
	size_t zerozero = intersectSize(left[0], right[0]);
	size_t zeroone = intersectSize(left[0], right[1]);
	size_t onezero = intersectSize(left[1], right[0]);
	size_t oneone = intersectSize(left[1], right[1]);
	if (zerozero > 0 && zerozero < 10) return false;
	if (zeroone > 0 && zeroone < 10) return false;
	if (onezero > 0 && onezero < 10) return false;
	if (oneone > 0 && oneone < 10) return false;
	if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	return true;
}

bool allelesMatchPerfectly(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right)
{
	if (left.size() != right.size()) return false;
	phmap::flat_hash_set<size_t> rightMatched;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i].size() < 5) return false;
		size_t uniqueMatch = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < right.size(); j++)
		{
			if (right[j].size() < 5) return false;
			if (left[i].size() != right[j].size()) continue;
			bool match = true;
			for (size_t k = 0; k < left[i].size(); k++)
			{
				if (left[i][k] != right[j][k])
				{
					match = false;
					break;
				}
			}
			if (!match) continue;
			if (uniqueMatch != std::numeric_limits<size_t>::max()) return false;
			uniqueMatch = j;
		}
		if (uniqueMatch == std::numeric_limits<size_t>::max()) return false;
		if (rightMatched.count(uniqueMatch) == 1) return false;
		rightMatched.insert(uniqueMatch);
	}
	assert(rightMatched.size() == right.size());
	return true;
}

std::string getAlleleStr(const TwobitString& readSequence, const size_t readStart, const size_t readEnd, const bool fw)
{
	assert(readEnd > readStart);
	std::string str = readSequence.substr(readStart, readEnd - readStart);
	if (!fw) str = revCompRaw(str);
	return str;
}

size_t getAllele(const TwobitString& readSequence, const size_t readStart, const size_t readEnd, phmap::flat_hash_map<std::string, size_t>& alleleNumbers, const bool fw)
{
	std::string allele = getAlleleStr(readSequence, readStart, readEnd, fw);
	if (alleleNumbers.count(allele) == 0)
	{
		size_t result = alleleNumbers.size();
		alleleNumbers[allele] = result;
	}
	return alleleNumbers.at(allele);
}

size_t getHammingdistance(const RankBitvector& left, const RankBitvector& right, const size_t maxEdits)
{
	assert(left.size() == right.size());
	if (left.size() == 0) return 0;
	const std::vector<uint64_t>& leftbits = left.getBits();
	const std::vector<uint64_t>& rightbits = right.getBits();
	assert(leftbits.size() == rightbits.size());
	assert(left.size() <= leftbits.size()*64);
	size_t result = 0;
	for (size_t i = 0; i*64+64 < left.size(); i++)
	{
		uint64_t notEqual = leftbits[i] ^ rightbits[i];
		result += popcount(notEqual);
		if (result > maxEdits) return result;
	}
	for (size_t i = ((size_t)((left.size()-1)/64))*64; i < left.size(); i++)
	{
		if (left.get(i) != right.get(i))
		{
			result += 1;
		}
	}
	return result;
}

bool siteIsPerfectlyPhased(const std::vector<RankBitvector>& columns, const size_t left, const size_t right)
{
	size_t zerozero = 0;
	size_t zeroone = 0;
	size_t onezero = 0;
	size_t oneone = 0;
	assert(columns[left].size() == columns[right].size());
	const std::vector<uint64_t>& leftbits = columns[left].getBits();
	const std::vector<uint64_t>& rightbits = columns[right].getBits();
	assert(leftbits.size() == rightbits.size());
	for (size_t i = 0; i*64+64 < columns[left].size(); i++)
	{
		oneone += popcount(leftbits[i] & rightbits[i]);
		zeroone += popcount(~leftbits[i] & rightbits[i]);
		zerozero += popcount(~leftbits[i] & ~rightbits[i]);
		onezero += popcount(leftbits[i] & ~rightbits[i]);
		if (oneone > 0 || zerozero > 0) return false;
	}
	size_t remainingBits = columns[left].size() % 64;
	uint64_t mask = 0xFFFFFFFFFFFFFFFFull;
	if (remainingBits > 0) // remainingbits 0 means the whole block is valid
	{
		mask >>= (64-remainingBits);
	}
	oneone += popcount((leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zeroone += popcount((~leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zerozero += popcount((~leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	onezero += popcount((leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	if (oneone == 0 && zerozero == 0 && onezero >= 5 && zeroone >= 5) return true;
	return false;
}

bool siteIsPhasedTwoVariantsThreeHaps(const std::vector<RankBitvector>& columns, const size_t left, const size_t right)
{
	size_t zerozero = 0;
	size_t zeroone = 0;
	size_t onezero = 0;
	size_t oneone = 0;
	assert(columns[left].size() == columns[right].size());
	const std::vector<uint64_t>& leftbits = columns[left].getBits();
	const std::vector<uint64_t>& rightbits = columns[right].getBits();
	assert(leftbits.size() == rightbits.size());
	for (size_t i = 0; i*64+64 < columns[left].size(); i++)
	{
		oneone += popcount(leftbits[i] & rightbits[i]);
		zeroone += popcount(~leftbits[i] & rightbits[i]);
		zerozero += popcount(~leftbits[i] & ~rightbits[i]);
		onezero += popcount(leftbits[i] & ~rightbits[i]);
		if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	}
	size_t remainingBits = columns[left].size() % 64;
	uint64_t mask = 0xFFFFFFFFFFFFFFFFull;
	if (remainingBits > 0) // remainingbits 0 means the whole block is valid
	{
		mask >>= (64-remainingBits);
	}
	oneone += popcount((leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zeroone += popcount((~leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zerozero += popcount((~leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	onezero += popcount((leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	if (zerozero > 0 && zerozero < 10) return false;
	if (onezero > 0 && onezero < 10) return false;
	if (zeroone > 0 && zeroone < 10) return false;
	if (oneone > 0 && oneone < 10) return false;
	if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	return true;
}

bool siteIsInformative(const std::vector<RankBitvector>& columns, const size_t left, const size_t right)
{
	size_t zerozero = 0;
	size_t zeroone = 0;
	size_t onezero = 0;
	size_t oneone = 0;
	assert(columns[left].size() == columns[right].size());
	const std::vector<uint64_t>& leftbits = columns[left].getBits();
	const std::vector<uint64_t>& rightbits = columns[right].getBits();
	assert(leftbits.size() == rightbits.size());
	for (size_t i = 0; i*64+64 < columns[left].size(); i++)
	{
		oneone += popcount(leftbits[i] & rightbits[i]);
		zeroone += popcount(~leftbits[i] & rightbits[i]);
		zerozero += popcount(~leftbits[i] & ~rightbits[i]);
		onezero += popcount(leftbits[i] & ~rightbits[i]);
		if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	}
	size_t remainingBits = columns[left].size() % 64;
	uint64_t mask = 0xFFFFFFFFFFFFFFFFull;
	if (remainingBits > 0) // remainingbits 0 means the whole block is valid
	{
		mask >>= (64-remainingBits);
	}
	oneone += popcount((leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zeroone += popcount((~leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zerozero += popcount((~leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	onezero += popcount((leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	if (zerozero + zeroone < 5) return false;
	if (onezero + oneone < 5) return false;
	if (zerozero + onezero < 5) return false;
	if (zeroone + oneone < 5) return false;
	if (zerozero > 0 && zerozero < 5) return false;
	if (zeroone > 0 && zeroone < 5) return false;
	if (onezero > 0 && onezero < 5) return false;
	if (oneone > 0 && oneone < 5) return false;
	if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	return true;
}

template <typename T>
void filterVector(std::vector<T>& vec, const std::vector<bool>& keep)
{
	std::vector<size_t> getIndexFrom;
	for (size_t i = 0; i < keep.size(); i++)
	{
		if (!keep[i]) continue;
		getIndexFrom.emplace_back(i);
	}
	assert(getIndexFrom.size() <= keep.size());
	assert(getIndexFrom.size() == 0 || getIndexFrom.back() < keep.size());
	for (size_t i = 0; i < getIndexFrom.size(); i++)
	{
		vec[i] = vec[getIndexFrom[i]];
	}
	vec.resize(getIndexFrom.size());
}

void filterMatrix(std::vector<RankBitvector>& matrix, const std::vector<bool>& keep)
{
	std::vector<size_t> getIndexFrom;
	for (size_t i = 0; i < keep.size(); i++)
	{
		if (!keep[i]) continue;
		getIndexFrom.emplace_back(i);
	}
	assert(getIndexFrom.size() <= keep.size());
	assert(getIndexFrom.size() == 0 || getIndexFrom.back() < keep.size());
	for (size_t i = 0; i < matrix.size(); i++)
	{
		assert(matrix[i].size() == keep.size());
		for (size_t k = 0; k < getIndexFrom.size(); k++)
		{
			matrix[i].set(k, matrix[i].get(getIndexFrom[k]));
		}
		matrix[i].resize(getIndexFrom.size());
	}
}

template <typename F>
void iterateChunksByCoverage(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, F callback)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	std::vector<size_t> iterationOrder;
	iterationOrder.reserve(occurrencesPerChunk.size());
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		iterationOrder.emplace_back(i);
	}
	std::sort(iterationOrder.begin(), iterationOrder.end(), [&occurrencesPerChunk](size_t left, size_t right) { return occurrencesPerChunk[left].size() > occurrencesPerChunk[right].size(); });
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&occurrencesPerChunk, &iterationOrder, callback](const size_t iterationIndex)
	{
		callback(iterationOrder[iterationIndex], occurrencesPerChunk);
	});
}

void checkAddMismatch(std::vector<std::vector<std::pair<size_t, size_t>>>& closestNeighborAndMismatchesPerLine, std::vector<size_t>& maxEditsHerePerLine, std::vector<size_t>& countEqualToMaxEditsPerLine, const size_t addHere, const size_t addThis, const size_t mismatches, const size_t countNeighbors)
{
	if (closestNeighborAndMismatchesPerLine[addHere].size() < countNeighbors)
	{
		closestNeighborAndMismatchesPerLine[addHere].emplace_back(mismatches, addThis);
		if (closestNeighborAndMismatchesPerLine[addHere].size() == countNeighbors)
		{
			maxEditsHerePerLine[addHere] = 0;
			for (size_t l = 0; l < closestNeighborAndMismatchesPerLine[addHere].size(); l++)
			{
				if (closestNeighborAndMismatchesPerLine[addHere][l].first < maxEditsHerePerLine[addHere]) continue;
				if (closestNeighborAndMismatchesPerLine[addHere][l].first > maxEditsHerePerLine[addHere]) countEqualToMaxEditsPerLine[addHere] = 0;
				maxEditsHerePerLine[addHere] = closestNeighborAndMismatchesPerLine[addHere][l].first;
				countEqualToMaxEditsPerLine[addHere] += 1;
			}
		}
	}
	else
	{
		if (mismatches <= maxEditsHerePerLine[addHere])
		{
			closestNeighborAndMismatchesPerLine[addHere].emplace_back(mismatches, addThis);
			if (mismatches == maxEditsHerePerLine[addHere])
			{
				countEqualToMaxEditsPerLine[addHere] += 1;
			}
			else
			{
				if (closestNeighborAndMismatchesPerLine[addHere].size() > countNeighbors+countEqualToMaxEditsPerLine[addHere])
				{
					size_t newMaxEdits = 0;
					countEqualToMaxEditsPerLine[addHere] = 0;
					for (size_t l = closestNeighborAndMismatchesPerLine[addHere].size()-1; l < closestNeighborAndMismatchesPerLine[addHere].size(); l--)
					{
						if (closestNeighborAndMismatchesPerLine[addHere][l].first == maxEditsHerePerLine[addHere])
						{
							std::swap(closestNeighborAndMismatchesPerLine[addHere][l], closestNeighborAndMismatchesPerLine[addHere].back());
							closestNeighborAndMismatchesPerLine[addHere].pop_back();
						}
						else
						{
							assert(closestNeighborAndMismatchesPerLine[addHere][l].first < maxEditsHerePerLine[addHere]);
							if (closestNeighborAndMismatchesPerLine[addHere][l].first > newMaxEdits) countEqualToMaxEditsPerLine[addHere] = 0;
							newMaxEdits = std::max(newMaxEdits, closestNeighborAndMismatchesPerLine[addHere][l].first);
							if (closestNeighborAndMismatchesPerLine[addHere][l].first == newMaxEdits) countEqualToMaxEditsPerLine[addHere] += 1;
						}
					}
					assert(newMaxEdits < maxEditsHerePerLine[addHere]);
					maxEditsHerePerLine[addHere] = newMaxEdits;
				}
			}
		}
	}
}

template <typename F>
std::vector<size_t> getFastTransitiveClosure(const size_t itemCount, const size_t maxDistance, F distanceFunction)
{
	std::vector<size_t> parent;
	parent.resize(itemCount);
	for (size_t i = 0; i < itemCount; i++)
	{
		parent[i] = i;
	}
	std::vector<size_t> clusterExample;
	std::vector<std::vector<size_t>> clusterAdditionals;
	std::vector<size_t> clusterMaxDistance;
	for (size_t i = 0; i < itemCount; i++)
	{
		bool found = false;
		for (size_t j = 0; j < clusterExample.size(); j++)
		{
			size_t distance = distanceFunction(i, clusterExample[j], maxDistance);
			if (distance <= maxDistance)
			{
				clusterAdditionals[j].emplace_back(i);
				found = true;
				clusterMaxDistance[j] = std::max(distance, clusterMaxDistance[j]);
				merge(parent, i, clusterExample[j]);
				break;
			}
		}
		if (!found)
		{
			clusterExample.emplace_back(i);
			clusterAdditionals.emplace_back();
			clusterMaxDistance.emplace_back(0);
		}
	}
	for (size_t i = 1; i < clusterExample.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (find(parent, clusterExample[i]) == find(parent, clusterExample[j])) continue;
			size_t distance = distanceFunction(clusterExample[i], clusterExample[j], maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j]);
			if (distance <= maxDistance)
			{
				merge(parent, clusterExample[i], clusterExample[j]);
				continue;
			}
			else if (distance <= maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j])
			{
				bool found = false;
				for (size_t k : clusterAdditionals[i])
				{
					if (distanceFunction(k, clusterExample[j], maxDistance) <= maxDistance)
					{
						found = true;
						merge(parent, k, clusterExample[j]);
						break;
					}
				}
				if (!found)
				{
					for (size_t l : clusterAdditionals[j])
					{
						if (distanceFunction(l, clusterExample[i], maxDistance) <= maxDistance)
						{
							found = true;
							merge(parent, l, clusterExample[i]);
							break;
						}
					}
				}
				if (!found)
				{
					for (size_t k : clusterAdditionals[i])
					{
						for (size_t l : clusterAdditionals[j])
						{
							if (distanceFunction(k, l, maxDistance) <= maxDistance)
							{
								found = true;
								merge(parent, k, l);
								break;
							}
						}
						if (found) break;
					}
				}
			}
		}
	}
	return parent;
}

std::vector<RankBitvector> getCorrectedMatrix(const std::vector<RankBitvector>& matrix, const size_t countNeighbors)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> closestNeighborAndMismatchesPerLine;
	std::vector<size_t> maxEditsHerePerLine;
	std::vector<size_t> countEqualToMaxEditsPerLine;
	closestNeighborAndMismatchesPerLine.resize(matrix.size());
	maxEditsHerePerLine.resize(matrix.size(), std::numeric_limits<size_t>::max());
	countEqualToMaxEditsPerLine.resize(matrix.size(), 0);
	for (size_t distance = 1; distance < matrix.size(); distance++)
	{
		for (size_t high = distance; high < matrix.size(); high++)
		{
			size_t low = high - distance;
			size_t mismatches = getHammingdistance(matrix[low], matrix[high], std::max(maxEditsHerePerLine[low], maxEditsHerePerLine[high]));
			checkAddMismatch(closestNeighborAndMismatchesPerLine, maxEditsHerePerLine, countEqualToMaxEditsPerLine, low, high, mismatches, countNeighbors);
			checkAddMismatch(closestNeighborAndMismatchesPerLine, maxEditsHerePerLine, countEqualToMaxEditsPerLine, high, low, mismatches, countNeighbors);
		}
	}
	std::vector<RankBitvector> correctedMatrix;
	correctedMatrix.resize(matrix.size());
	for (size_t j = 0; j < matrix.size(); j++)
	{
		correctedMatrix[j].resize(matrix[j].size());
		closestNeighborAndMismatchesPerLine[j].emplace_back(0, j);
		std::sort(closestNeighborAndMismatchesPerLine[j].begin(), closestNeighborAndMismatchesPerLine[j].end());
		size_t endIndex = countNeighbors;
		while (endIndex+1 < closestNeighborAndMismatchesPerLine[j].size() && closestNeighborAndMismatchesPerLine[j][endIndex+1].first == closestNeighborAndMismatchesPerLine[j][countNeighbors].first)
		{
			endIndex += 1;
		}
		for (size_t k = 0; k < matrix[j].size(); k++)
		{
			size_t ones = 0;
			for (size_t m = 0; m <= endIndex; m++)
			{
				if (matrix[closestNeighborAndMismatchesPerLine[j][m].second].get(k)) ones += 1;
			}
			correctedMatrix[j].set(k, ones >= (endIndex+1)/2);
		}
	}
	return correctedMatrix;
}

bool rowEqual(const RankBitvector& left, const RankBitvector& right)
{
	assert(left.size() == right.size());
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left.get(i) != right.get(i)) return false;
	}
	return true;
}

void sortAndInterleave(std::vector<RankBitvector>& correctedMatrix)
{
	std::sort(correctedMatrix.begin(), correctedMatrix.end(), [](const RankBitvector& left, const RankBitvector& right)
	{
		for(size_t i = 0; i < left.size(); i++)
		{
			if (!left.get(i) && right.get(i)) return true;
			if (left.get(i) && !right.get(i)) return false;
		}
		return false;
	});
	size_t countDistinctRows = 1;
	for (size_t i = 1; i < correctedMatrix.size(); i++)
	{
		if (!rowEqual(correctedMatrix[i], correctedMatrix[i-1])) countDistinctRows += 1;
	}
	std::swap(correctedMatrix[1], correctedMatrix.back());
	size_t pos = 2;
	for (size_t div = 2; div <= 64; div *= 2)
	{
		for (size_t off = 1; off < div; off += 2)
		{
			std::swap(correctedMatrix[pos], correctedMatrix[correctedMatrix.size()/div*off]);
			pos += 1;
		}
	}
	assert(pos == 65);
}

void splitPerCorrectedKmerPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by corrected kmer phasing" << std::endl;
	const size_t countNeighbors = 5;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, countNeighbors, kmerSize](size_t dummy)
	{
		while (true)
		{
			std::vector<std::pair<size_t, size_t>> chunkBeingDone;
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunksNeedProcessing.size() == 0) return;
				std::swap(chunkBeingDone, chunksNeedProcessing.back());
				chunksNeedProcessing.pop_back();
			}
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				chunksDoneProcessing.emplace_back();
				std::cerr << "corrected kmer skipped chunk with coverage " << chunkBeingDone.size() << std::endl;
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerClusterToColumn;
			std::vector<RankBitvector> matrix; // doesn't actually use rank in any way, but exposes uint64 for bitparallel hamming distance calculation
			std::vector<std::pair<size_t, size_t>> clusters;
			matrix.resize(chunkBeingDone.size());
			std::vector<TwobitString> chunkSequences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
			}
			auto validClusters = iterateSolidKmers(chunkSequences, kmerSize, 5, true, true, [&kmerClusterToColumn, &clusters, &matrix](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				assert(occurrenceID < matrix.size());
				size_t column = std::numeric_limits<size_t>::max();
				if (kmerClusterToColumn.count(std::make_pair(kmer, clusterIndex)) == 1)
				{
					column = kmerClusterToColumn.at(std::make_pair(kmer, clusterIndex));
				}
				else
				{
					column = kmerClusterToColumn.size();
					kmerClusterToColumn[std::make_pair(kmer, clusterIndex)] = column;
					clusters.emplace_back(std::make_pair(kmer, clusterIndex));
				}
				if (matrix[occurrenceID].size() <= column)
				{
					assert(column < kmerClusterToColumn.size());
					matrix[occurrenceID].resize(kmerClusterToColumn.size());
				}
				assert(!matrix[occurrenceID].get(column));
				matrix[occurrenceID].set(column, true);
			});
			std::vector<std::pair<double, double>> kmerLocation;
			kmerLocation.resize(clusters.size(), std::make_pair(-1, -1));
			for (auto pair : kmerClusterToColumn)
			{
				size_t index = pair.second;
				assert(kmerLocation[index].first == -1);
				kmerLocation[index] = validClusters.at(pair.first.first)[pair.first.second];
			}
			for (size_t i = 0; i < kmerLocation.size(); i++)
			{
				assert(kmerLocation[i].first != -1);
			}
			for (size_t j = 0; j < matrix.size(); j++)
			{
				if (matrix[j].size() == kmerClusterToColumn.size()) continue;
				assert(matrix[j].size() < kmerClusterToColumn.size());
				matrix[j].resize(kmerClusterToColumn.size());
			}
			size_t columnsInUnfiltered = matrix[0].size();
			std::vector<bool> covered;
			covered.resize(matrix[0].size(), false);
			size_t columnsInCovered = 0;
			std::vector<RankBitvector> columns;
			columns.resize(matrix[0].size());
			for (size_t j = 0; j < matrix[0].size(); j++)
			{
				columns[j].resize(matrix.size());
				size_t zeros = 0;
				size_t ones = 0;
				for (size_t k = 0; k < matrix.size(); k++)
				{
					assert(matrix[k].size() == matrix[0].size());
					columns[j].set(k, matrix[k].get(j));
					if (matrix[k].get(j))
					{
						ones += 1;
					}
					else
					{
						zeros += 1;
					}
				}
				if (zeros >= 5 && ones >= 5)
				{
					covered[j] = true;
					columnsInCovered += 1;
				}
			}
			if (columnsInCovered < 2)
			{
				auto endTime = getTime();
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "corrected kmer skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " due to covered " << columnsInCovered << " time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<bool> informativeSite;
			informativeSite.resize(matrix[0].size(), false);
			assert(clusters.size() == informativeSite.size());
			assert(clusters.size() == kmerClusterToColumn.size());
			for (size_t j = 0; j < informativeSite.size(); j++)
			{
				if (!covered[j]) continue;
				for (size_t k = 0; k < j; k++)
				{
					if (!covered[k]) continue;
					if (informativeSite[j] && informativeSite[k]) continue;
					assert(validClusters.at(clusters[j].first).size() > clusters[j].second);
					assert(validClusters.at(clusters[k].first).size() > clusters[k].second);
					if (validClusters.at(clusters[j].first)[clusters[j].second].first < validClusters.at(clusters[k].first)[clusters[k].second].second + 1)
					{
						if (validClusters.at(clusters[j].first)[clusters[j].second].second + 1 > validClusters.at(clusters[k].first)[clusters[k].second].first)
						{
							continue;
						}
					}
					if (siteIsInformative(columns, j, k) || siteIsInformative(columns, k, j))
					{
						informativeSite[j] = true;
						informativeSite[k] = true;
					}
				}
			}
			size_t columnsInInformative = 0;
			for (size_t j = 0; j < informativeSite.size(); j++)
			{
				if (informativeSite[j]) columnsInInformative += 1;
			}
			if (columnsInInformative < 2)
			{
				auto endTime = getTime();
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "corrected kmer skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " due to informative " << columnsInInformative << " time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			filterVector(columns, informativeSite);
			filterVector(kmerLocation, informativeSite);
			filterMatrix(matrix, informativeSite);
			std::vector<RankBitvector> correctedMatrix = getCorrectedMatrix(matrix, countNeighbors);
			assert(correctedMatrix.size() == matrix.size());
			assert(correctedMatrix.size() == chunkBeingDone.size());
			if (correctedMatrix.size() >= 1024)
			{
				std::vector<RankBitvector> columnMatrix = correctedMatrix;
				sortAndInterleave(columnMatrix);
				for (size_t j = 0; j < columnMatrix[0].size(); j++)
				{
					for (size_t k = 0; k < columnMatrix.size(); k++)
					{
						assert(columnMatrix[k].size() == columnMatrix[0].size());
						columns[j].set(k, columnMatrix[k].get(j));
					}
				}
			}
			std::vector<bool> phasingSite;
			phasingSite.resize(matrix[0].size(), false);
			std::vector<size_t> phasedSoFar;
			std::vector<size_t> notPhasedSoFar;
			notPhasedSoFar.emplace_back(0);
			for (size_t j = 1; j < matrix[0].size(); j++)
			{
				for (auto k : phasedSoFar)
				{
					if (kmerLocation[j].first < kmerLocation[k].second+1 && kmerLocation[j].second + 1 > kmerLocation[k].first) continue;
					if (siteIsPerfectlyPhased(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						break;
					}
					else if (siteIsPhasedTwoVariantsThreeHaps(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						break;
					}
				}
				for (size_t kindex = notPhasedSoFar.size()-1; kindex < notPhasedSoFar.size(); kindex--)
				{
					size_t k = notPhasedSoFar[kindex];
					if (kmerLocation[j].first < kmerLocation[k].second+1 && kmerLocation[j].second + 1 > kmerLocation[k].first) continue;
					if (siteIsPerfectlyPhased(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						phasedSoFar.emplace_back(k);
						std::swap(notPhasedSoFar[kindex], notPhasedSoFar.back());
						notPhasedSoFar.pop_back();
					}
					else if (siteIsPhasedTwoVariantsThreeHaps(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						phasedSoFar.emplace_back(k);
						std::swap(notPhasedSoFar[kindex], notPhasedSoFar.back());
						notPhasedSoFar.pop_back();
					}
				}
				if (phasingSite[j])
				{
					phasedSoFar.emplace_back(j);
				}
				else
				{
					notPhasedSoFar.emplace_back(j);
				}
			}
			filterMatrix(correctedMatrix, phasingSite);
			size_t numPhasingSites = correctedMatrix[0].size();
			std::vector<size_t> parent = getFastTransitiveClosure(correctedMatrix.size(), 0, [&correctedMatrix](const size_t i, const size_t j, const size_t maxDist) { return getHammingdistance(correctedMatrix[i], correctedMatrix[j], maxDist); });
			auto endTime = getTime();
			size_t nextNum = 0;
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
			if (keyToNode.size() == 1)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "corrected kmer splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " phasing sites " << numPhasingSites << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.resize(keyToNode.size());
			for (size_t j = 0; j < parent.size(); j++)
			{
				chunkResult[keyToNode.at(find(parent, j))].emplace_back(chunkBeingDone[j]);
			}
			// sort smallest last, so emplace-pop-swap puts biggest on top of chunksNeedProcessing
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (keyToNode.size() >= 2) countSplitted += 1;
				std::cerr << "corrected kmer splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " phasing sites " << numPhasingSites << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				while (chunkResult.size() > 0)
				{
					chunksNeedProcessing.emplace_back();
					std::swap(chunksNeedProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
			}
		}
	});
	for (size_t i = 0; i < chunksDoneProcessing.size(); i++)
	{
		for (auto pair : chunksDoneProcessing[i])
		{
			std::get<2>(chunksPerRead[pair.first][pair.second]) = i + (std::get<2>(chunksPerRead[pair.first][pair.second]) & firstBitUint64_t);
		}
	}
	std::cerr << "corrected kmer phasing splitted " << countSplitted << " chunks" << std::endl;
}

class ScopedCounterIncrementer
{
public:
	ScopedCounterIncrementer(std::atomic<size_t>& counter) :
	counter(counter)
	{
		counter += 1;
	}
	~ScopedCounterIncrementer()
	{
		counter -= 1;
	}
private:
	std::atomic<size_t>& counter;
};

void splitPerNearestNeighborPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by nearest neighbor phasing" << std::endl;
	const size_t countNeighbors = 5;
	const size_t countDifferences = 100;
	size_t countSplitted = 0;
	std::atomic<size_t> threadsRunning;
	threadsRunning = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, &threadsRunning, countNeighbors, countDifferences, kmerSize](size_t dummy)
	{
		while (true)
		{
			std::vector<std::pair<size_t, size_t>> chunkBeingDone;
			bool gotOne = false;
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunksNeedProcessing.size() >= 1)
				{
					std::swap(chunkBeingDone, chunksNeedProcessing.back());
					chunksNeedProcessing.pop_back();
					gotOne = true;
				}
				else
				{
					if (threadsRunning == 0) return;
				}
			}
			if (!gotOne)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				continue;
			}
			ScopedCounterIncrementer threadCounter { threadsRunning };
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "skip chunk with size " << chunkBeingDone.size() << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerClusterToColumn;
			std::vector<RankBitvector> matrix; // doesn't actually use rank in any way, but exposes uint64 for bitparallel hamming distance calculation
			std::vector<std::pair<size_t, size_t>> clusters;
			matrix.resize(chunkBeingDone.size());
			std::vector<TwobitString> chunkSequences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
			}
			auto validClusters = iterateSolidKmers(chunkSequences, kmerSize, 5, true, true, [&kmerClusterToColumn, &clusters, &matrix](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				assert(occurrenceID < matrix.size());
				size_t column = std::numeric_limits<size_t>::max();
				if (kmerClusterToColumn.count(std::make_pair(kmer, clusterIndex)) == 1)
				{
					column = kmerClusterToColumn.at(std::make_pair(kmer, clusterIndex));
				}
				else
				{
					column = kmerClusterToColumn.size();
					kmerClusterToColumn[std::make_pair(kmer, clusterIndex)] = column;
					clusters.emplace_back(std::make_pair(kmer, clusterIndex));
				}
				if (matrix[occurrenceID].size() <= column)
				{
					assert(column < kmerClusterToColumn.size());
					matrix[occurrenceID].resize(kmerClusterToColumn.size());
				}
				assert(!matrix[occurrenceID].get(column));
				matrix[occurrenceID].set(column, true);
			});
			for (size_t j = 0; j < matrix.size(); j++)
			{
				if (matrix[j].size() == kmerClusterToColumn.size()) continue;
				assert(matrix[j].size() < kmerClusterToColumn.size());
				matrix[j].resize(kmerClusterToColumn.size());
			}
			size_t columnsInUnfiltered = matrix[0].size();
			std::vector<bool> covered;
			covered.resize(matrix[0].size(), false);
			size_t columnsInCovered = 0;
			std::vector<RankBitvector> columns;
			columns.resize(matrix[0].size());
			for (size_t j = 0; j < matrix[0].size(); j++)
			{
				columns[j].resize(matrix.size());
				size_t zeros = 0;
				size_t ones = 0;
				for (size_t k = 0; k < matrix.size(); k++)
				{
					assert(matrix[k].size() == matrix[0].size());
					columns[j].set(k, matrix[k].get(j));
				}
				const std::vector<uint64_t>& bits = columns[j].getBits();
				for (size_t k = 0; k*64+64 <= columns[j].size(); k++)
				{
					ones += popcount(bits[k]);
					zeros += 64-popcount(bits[k]);
					if (zeros >= 5 && ones >= 5) break;
				}
				if (zeros < 5 || ones < 5)
				{
					for (size_t k = ((int)(columns[j].size()/64))*64; k < columns[j].size(); k++)
					{
						if (columns[j].get(k))
						{
							ones += 1;
						}
						else
						{
							zeros += 1;
						}
					}
				}
				if (zeros >= 5 && ones >= 5)
				{
					covered[j] = true;
					columnsInCovered += 1;
				}
			}
			if (columnsInCovered < countDifferences)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				auto endTime = getTime();
				std::cerr << "nearest neighbor skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " due to covered " << columnsInCovered << ", time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<bool> informativeSite;
			informativeSite.resize(matrix[0].size(), false);
			assert(clusters.size() == informativeSite.size());
			assert(clusters.size() == kmerClusterToColumn.size());
			for (size_t j = 0; j < informativeSite.size(); j++)
			{
				if (!covered[j]) continue;
				for (size_t k = 0; k < j; k++)
				{
					if (!covered[k]) continue;
					if (informativeSite[j] && informativeSite[k]) continue;
					assert(validClusters.at(clusters[j].first).size() > clusters[j].second);
					assert(validClusters.at(clusters[k].first).size() > clusters[k].second);
					if (validClusters.at(clusters[j].first)[clusters[j].second].first < validClusters.at(clusters[k].first)[clusters[k].second].second + 1)
					{
						if (validClusters.at(clusters[j].first)[clusters[j].second].second + 1 > validClusters.at(clusters[k].first)[clusters[k].second].first)
						{
							continue;
						}
					}
					if (siteIsInformative(columns, j, k) || siteIsInformative(columns, k, j))
					{
						informativeSite[j] = true;
						informativeSite[k] = true;
					}
				}
			}
			size_t columnsInInformative = 0;
			for (size_t i = 0; i < informativeSite.size(); i++)
			{
				if (informativeSite[i]) columnsInInformative += 1;
			}
			if (columnsInInformative < countDifferences)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				auto endTime = getTime();
				std::cerr << "nearest neighbor skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " due to informative " << columnsInInformative << ", time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			filterMatrix(matrix, informativeSite);
			std::vector<RankBitvector> correctedMatrix = getCorrectedMatrix(matrix, countNeighbors);
			assert(correctedMatrix.size() == matrix.size());
			assert(correctedMatrix.size() == chunkBeingDone.size());
			std::vector<size_t> parent = getFastTransitiveClosure(correctedMatrix.size(), countDifferences, [&correctedMatrix](const size_t i, const size_t j, const size_t maxDist) { return getHammingdistance(correctedMatrix[i], correctedMatrix[j], maxDist); });
			auto endTime = getTime();
			size_t nextNum = 0;
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
			if (keyToNode.size() == 1)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "nearest neighbor splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.resize(keyToNode.size());
			for (size_t j = 0; j < parent.size(); j++)
			{
				chunkResult[keyToNode.at(find(parent, j))].emplace_back(chunkBeingDone[j]);
			}
			// sort smallest last, so emplace-pop-swap puts biggest on top of chunksNeedProcessing
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (keyToNode.size() >= 2) countSplitted += 1;
				std::cerr << "nearest neighbor splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				while (chunkResult.size() > 0)
				{
					chunksNeedProcessing.emplace_back();
					std::swap(chunksNeedProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
			}
		}
	});
	for (size_t i = 0; i < chunksDoneProcessing.size(); i++)
	{
		for (auto pair : chunksDoneProcessing[i])
		{
			std::get<2>(chunksPerRead[pair.first][pair.second]) = i + (std::get<2>(chunksPerRead[pair.first][pair.second]) & firstBitUint64_t);
		}
	}
	std::cerr << "nearest neighbor phasing splitted " << countSplitted << " chunks" << std::endl;
}

void splitPerAllelePhasingWithinChunk(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by allele phasing" << std::endl;
	size_t countSplitted = 0;
	size_t countSplittedTo = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &countSplittedTo, &rawReadLengths, kmerSize](size_t dummy)
	{
		while (true)
		{
			std::vector<std::pair<size_t, size_t>> chunkBeingDone;
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunksNeedProcessing.size() == 0) return;
				std::swap(chunkBeingDone, chunksNeedProcessing.back());
				chunksNeedProcessing.pop_back();
			}
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
//				std::cerr << "skip chunk with coverage " << chunkBeingDone.size() << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>> alleles; // startkmer, startcluster, endkmer, endcluster, occurrenceID, allele
			size_t lastOccurrence = std::numeric_limits<size_t>::max();
			size_t lastKmer = std::numeric_limits<size_t>::max();
			size_t lastCluster = std::numeric_limits<size_t>::max();
			size_t lastKmerPos = std::numeric_limits<size_t>::max();
			std::vector<TwobitString> chunkSequences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
			}
			phmap::flat_hash_map<std::string, size_t> alleleNumbers;
			iterateSolidKmers(chunkSequences, kmerSize, chunkBeingDone.size(), false, false, [&chunkBeingDone, &chunkSequences, &alleles, &lastOccurrence, &lastKmer, &lastCluster, &lastKmerPos, &alleleNumbers, kmerSize](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				if (occurrenceID != lastOccurrence || lastKmer == std::numeric_limits<size_t>::max())
				{
					lastOccurrence = occurrenceID;
					lastKmer = kmer;
					lastCluster = clusterIndex;
					lastKmerPos = pos;
					return;
				}
				size_t allele;
				allele = getAllele(chunkSequences[occurrenceID], lastKmerPos, pos + kmerSize, alleleNumbers, true);
				alleles.emplace_back(lastKmer, lastCluster, kmer, clusterIndex, occurrenceID, allele);
				lastOccurrence = occurrenceID;
				lastKmer = kmer;
				lastCluster = clusterIndex;
				lastKmerPos = pos;
			});
			std::sort(alleles.begin(), alleles.end());
			std::vector<std::vector<std::vector<size_t>>> occurrencesPerAlleleSite;
			occurrencesPerAlleleSite = getAlleleOccurrences(alleles, chunkBeingDone.size());
			std::vector<bool> informativeOccurrence;
			informativeOccurrence.resize(occurrencesPerAlleleSite.size(), false);
			for (size_t j = 1; j < occurrencesPerAlleleSite.size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					if (informativeOccurrence[j] && informativeOccurrence[k]) continue;
					if (allelesMatchPerfectly(occurrencesPerAlleleSite[j], occurrencesPerAlleleSite[k]))
					{
						informativeOccurrence[j] = true;
						informativeOccurrence[k] = true;
					}
					else if (allelesMatchTwoVariantsThreeHaps(occurrencesPerAlleleSite[j], occurrencesPerAlleleSite[k]))
					{
						informativeOccurrence[j] = true;
						informativeOccurrence[k] = true;
					}
				}
			}
			std::vector<std::vector<size_t>> occurrenceAlleles;
			occurrenceAlleles.resize(chunkBeingDone.size());
			for (size_t site = 0; site < occurrencesPerAlleleSite.size(); site++)
			{
				if (!informativeOccurrence[site]) continue;
				for (size_t allele = 0; allele < occurrencesPerAlleleSite[site].size(); allele++)
				{
					for (auto occurrence : occurrencesPerAlleleSite[site][allele])
					{
						occurrenceAlleles[occurrence].emplace_back(allele);
					}
				}
			}
			std::vector<size_t> orderedOccurrences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				orderedOccurrences.emplace_back(j);
			}
			std::sort(orderedOccurrences.begin(), orderedOccurrences.end(), [&occurrenceAlleles](size_t left, size_t right)
			{
				assert(occurrenceAlleles[left].size() == occurrenceAlleles[right].size());
				for (size_t i = 0; i < occurrenceAlleles[left].size(); i++)
				{
					if (occurrenceAlleles[left][i] < occurrenceAlleles[right][i]) return true;
					if (occurrenceAlleles[left][i] > occurrenceAlleles[right][i]) return false;
				}
				return false;
			});
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.emplace_back();
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				if (j > 0 && occurrenceAlleles[orderedOccurrences[j]] != occurrenceAlleles[orderedOccurrences[j-1]])
				{
					chunkResult.emplace_back();
				}
				chunkResult.back().emplace_back(chunkBeingDone[orderedOccurrences[j]]);
			}
			if (chunkResult.size() == 1)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
//				std::cerr << "chunk with coverage " << chunkBeingDone.size() << " not split" << std::endl;
//				std::cerr << "alleles " << alleles.size() << " occurrences " << occurrencesPerAlleleSite.size() << " occurrenceAlleles " << occurrenceAlleles.size() << " informative sites " << occurrenceAlleles[0].size() << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkResult[0]);
				continue;
			}
			assert(chunkResult.size() >= 2);
			// sort smallest last, so emplace-swap-pop puts biggest on top of chunksNeedProcessing
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
//				std::cerr << "chunk with coverage " << chunkBeingDone.size() << " split to " << chunkResult.size() << std::endl;
//				std::cerr << "alleles " << alleles.size() << " occurrences " << occurrencesPerAlleleSite.size() << " occurrenceAlleles " << occurrenceAlleles.size() << " informative sites " << occurrenceAlleles[0].size() << std::endl;
				countSplittedTo += chunkResult.size();
				countSplitted += 1;
				while (chunkResult.size() > 0)
				{
					chunksNeedProcessing.emplace_back();
					std::swap(chunksNeedProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
			}
		}
	});
	for (size_t i = 0; i < chunksDoneProcessing.size(); i++)
	{
		for (auto pair : chunksDoneProcessing[i])
		{
			std::get<2>(chunksPerRead[pair.first][pair.second]) = i + (std::get<2>(chunksPerRead[pair.first][pair.second]) & firstBitUint64_t);
		}
	}
	std::cerr << "allele phasing splitted " << countSplitted << " chunks to " << countSplittedTo << std::endl;
}

void splitPerPhasingKmersWithinChunk(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by phasing kmers" << std::endl;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, kmerSize](size_t dummy)
	{
		while (true)
		{
			std::vector<std::pair<size_t, size_t>> chunkBeingDone;
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunksNeedProcessing.size() == 0) return;
				std::swap(chunkBeingDone, chunksNeedProcessing.back());
				chunksNeedProcessing.pop_back();
			}
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "skip chunk with coverage " << chunkBeingDone.size() << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			phmap::flat_hash_map<std::pair<size_t, size_t>, phmap::flat_hash_map<size_t, size_t>> fwForks; // (hom kmer, hom cluster) -> (occurrence ID -> seq char)
			phmap::flat_hash_map<std::pair<size_t, size_t>, phmap::flat_hash_map<size_t, size_t>> bwForks; // (hom kmer, hom cluster) -> (occurrence ID -> seq char)
			size_t lastOccurrenceID = std::numeric_limits<size_t>::max();
			size_t lastKmer = std::numeric_limits<size_t>::max();
			size_t lastCluster = std::numeric_limits<size_t>::max();
			size_t lastKmerPos = std::numeric_limits<size_t>::max();
			std::vector<TwobitString> chunkSequences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
			}
			auto validClusters = iterateSolidKmers(chunkSequences, kmerSize, chunkBeingDone.size() * 0.95, false, true, [&fwForks, &bwForks, &lastOccurrenceID, &lastKmer, &lastKmerPos, &lastCluster, &chunkBeingDone, &chunkSequences, &chunksPerRead, kmerSize](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				if (occurrenceID != lastOccurrenceID || pos != lastKmerPos+1)
				{
					if (lastOccurrenceID != std::numeric_limits<size_t>::max())
					{
						if (lastKmerPos + kmerSize < chunkSequences[lastOccurrenceID].size())
						{
							assert(fwForks.count(std::make_pair(lastKmer, lastCluster)) == 0 || fwForks.at(std::make_pair(lastKmer, lastCluster)).count(lastOccurrenceID) == 0);
							fwForks[std::make_pair(lastKmer, lastCluster)][lastOccurrenceID] = chunkSequences[lastOccurrenceID].get(lastKmerPos + kmerSize);
						}
					}
					if (pos > 0)
					{
						assert(bwForks.count(std::make_pair(kmer, clusterIndex)) == 0 || bwForks.at(std::make_pair(kmer, clusterIndex)).count(occurrenceID) == 0);
						bwForks[std::make_pair(kmer, clusterIndex)][occurrenceID] = chunkSequences[occurrenceID].get(pos - 1);
					}
				}
				lastOccurrenceID = occurrenceID;
				lastCluster = clusterIndex;
				lastKmerPos = pos;
				lastKmer = kmer;
			});
			if (lastOccurrenceID != std::numeric_limits<size_t>::max())
			{
				if (lastKmerPos + kmerSize < chunkSequences[lastOccurrenceID].size())
				{
					assert(fwForks.count(std::make_pair(lastKmer, lastCluster)) == 0 || fwForks.at(std::make_pair(lastKmer, lastCluster)).count(lastOccurrenceID) == 0);
					fwForks[std::make_pair(lastKmer, lastCluster)][lastOccurrenceID] = chunkSequences[lastOccurrenceID].get(lastKmerPos + kmerSize);
				}
			}
			size_t checkedPairs = 0;
			std::vector<std::vector<size_t>> phaseIdentities;
			phaseIdentities.resize(chunkBeingDone.size());
			size_t checkedFwBw = 0;
			size_t checkedBwBw = 0;
			size_t checkedFwFw = 0;
			for (const auto& bwFork : bwForks)
			{
				for (const auto& fwFork : fwForks)
				{
					if (validClusters.at(bwFork.first.first)[bwFork.first.second].second > validClusters.at(fwFork.first.first)[fwFork.first.second].first) continue;
					checkedPairs += 1;
					checkedFwBw += 1;
					checkPhasablePair(bwFork.second, fwFork.second, phaseIdentities);
				}
			}
			for (const auto& bwFork : bwForks)
			{
				for (const auto& bwFork2 : bwForks)
				{
					if (validClusters.at(bwFork.first.first)[bwFork.first.second].second+1 > validClusters.at(bwFork2.first.first)[bwFork2.first.second].first && validClusters.at(bwFork2.first.first)[bwFork2.first.second].second+1 > validClusters.at(bwFork.first.first)[bwFork.first.second].first) continue;
					checkedPairs += 1;
					checkedBwBw += 1;
					checkPhasablePair(bwFork.second, bwFork2.second, phaseIdentities);
				}
			}
			for (const auto& fwFork : fwForks)
			{
				for (const auto& fwFork2 : fwForks)
				{
					if (validClusters.at(fwFork.first.first)[fwFork.first.second].second+1 > validClusters.at(fwFork2.first.first)[fwFork2.first.second].first && validClusters.at(fwFork2.first.first)[fwFork2.first.second].second+1 > validClusters.at(fwFork.first.first)[fwFork.first.second].first) continue;
					checkedPairs += 1;
					checkedFwFw += 1;
					checkPhasablePair(fwFork.second, fwFork2.second, phaseIdentities);
				}
			}
			std::vector<size_t> parent;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				parent.emplace_back(j);
			}
			for (size_t j = 1; j < phaseIdentities.size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					if (phaseIdentities[k] != phaseIdentities[j]) continue;
					merge(parent, k, j);
				}
			}
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			size_t nextNum = 0;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
			if (keyToNode.size() == 1)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				auto endTime = getTime();
				std::cerr << "phasing kmer splitted chunk with coverage " << chunkBeingDone.size() << " forkpairs " << phaseIdentities[0].size() << " checked " << checkedPairs << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				// std::cerr << "bwforks " << bwForks.size() << " fwforks " << fwForks.size() << std::endl;
				// std::cerr << "bwfw " << checkedFwBw << " bwbw " << checkedBwBw << " fwfw " << checkedFwFw << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.resize(keyToNode.size());
			for (size_t j = 0; j < parent.size(); j++)
			{
				chunkResult[keyToNode.at(find(parent, j))].emplace_back(chunkBeingDone[j]);
			}
			// sort smallest last, so emplace-pop-swap puts biggest on top of chunksNeedProcessing
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				auto endTime = getTime();
				if (keyToNode.size() >= 2) countSplitted += 1;
				std::cerr << "phasing kmer splitted chunk with coverage " << chunkBeingDone.size() << " forkpairs " << phaseIdentities[0].size() << " checked " << checkedPairs << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				// std::cerr << "bwforks " << bwForks.size() << " fwforks " << fwForks.size() << std::endl;
				// std::cerr << "bwfw " << checkedFwBw << " bwbw " << checkedBwBw << " fwfw " << checkedFwFw << std::endl;
				while (chunkResult.size() > 0)
				{
					chunksNeedProcessing.emplace_back();
					std::swap(chunksNeedProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
			}
		}
	});
	for (size_t i = 0; i < chunksDoneProcessing.size(); i++)
	{
		for (auto pair : chunksDoneProcessing[i])
		{
			std::get<2>(chunksPerRead[pair.first][pair.second]) = i + (std::get<2>(chunksPerRead[pair.first][pair.second]) & firstBitUint64_t);
		}
	}
	std::cerr << "phasing kmers splitted " << countSplitted << " chunks" << std::endl;
}

void splitPerSequenceIdentityRoughly(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting roughly by sequence identity" << std::endl;
	const size_t mismatchFloor = 10;
	size_t nextNum = 0;
	size_t totalSplitted = 0;
	size_t totalSplittedTo = 0;
	std::mutex resultMutex;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, &totalSplitted, &totalSplittedTo, mismatchFloor](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		if (occurrencesPerChunk[i].size() < 1000)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			// std::cerr << "skip chunk with coverage " << occurrencesPerChunk[i].size() << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		auto startTime = getTime();
		std::vector<std::string> sequences;
		sequences.reserve(occurrencesPerChunk[i].size());
		size_t longestLength = 0;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
			longestLength = std::max(longestLength, sequences.back().size());
		}
		std::vector<size_t> parent = getFastTransitiveClosure(sequences.size(), std::max(longestLength * mismatchFraction, mismatchFraction), [&sequences](const size_t i, const size_t j, const size_t maxDist) { return getNumMismatches(sequences[i], sequences[j], maxDist); });
		phmap::flat_hash_map<size_t, size_t> parentToCluster;
		for (size_t j = 0; j < parent.size(); j++)
		{
			find(parent, j);
		}
		size_t nextCluster = 0;
		for (size_t j = 0; j < parent.size(); j++)
		{
			if (parent[j] != j) continue;
			parentToCluster[j] = nextCluster;
			nextCluster += 1;
		}
		auto endTime = getTime();
		assert(nextCluster >= 1);
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			if (nextCluster > 1)
			{
				totalSplitted += 1;
				totalSplittedTo += nextCluster;
			}
			// std::cerr << "split chunk with coverage " << occurrencesPerChunk[i].size() << " into " << nextCluster << " chunks time " << formatTime(startTime, endTime) << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + parentToCluster.at(parent[j]) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += nextCluster;
		}
	});
	std::cerr << "rough sequence identity splitting splitted " << totalSplitted << " chunks to " << totalSplittedTo << " chunks" << std::endl;
}

void splitPerSequenceIdentity(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by sequence identity" << std::endl;
	const size_t mismatchFloor = 10;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	std::vector<size_t> iterationOrder;
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		iterationOrder.emplace_back(i);
	}
	std::sort(iterationOrder.begin(), iterationOrder.end(), [&occurrencesPerChunk](size_t left, size_t right) { return occurrencesPerChunk[left].size() > occurrencesPerChunk[right].size(); });
	size_t firstSmall = iterationOrder.size();
	for (size_t i = 0; i < iterationOrder.size(); i++)
	{
		if (occurrencesPerChunk[iterationOrder[i]].size() > 200 || occurrencesPerChunk[iterationOrder[i]].size() > 200*numThreads) continue;
		firstSmall = i;
		break;
	}
	size_t nextNum = 0;
	for (size_t iterationIndex = 0; iterationIndex < firstSmall; iterationIndex++)
	{
		const size_t i = iterationOrder[iterationIndex];
		// std::cerr << "try split chunk " << i << std::endl;
		// auto startTime = getTime();
		std::vector<std::pair<std::string, size_t>> sequencesPerOccurrence;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sequencesPerOccurrence.emplace_back();
			sequencesPerOccurrence.back().second = j;
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequencesPerOccurrence.back().first = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
		}
		std::sort(sequencesPerOccurrence.begin(), sequencesPerOccurrence.end(), [](const auto& left, const auto& right) {
			if (left.first.size() < right.first.size()) return true;
			if (left.first.size() > right.first.size()) return false;
			return left.first < right.first;
		});
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		const size_t blockSize = 50;
		const size_t numBlocks = (sequencesPerOccurrence.size()+blockSize-1)/blockSize;
		// std::cerr << "split part 1" << std::endl;
		iterateMultithreaded(0, numBlocks, numThreads, [&sequencesPerOccurrence, &parent, blockSize, i, mismatchFloor](size_t index)
		{
			for (size_t j = index * blockSize+1; j < (index+1) * blockSize && j < sequencesPerOccurrence.size(); j++)
			{
				for (size_t k = j-1; k < j && k >= index*blockSize; k--)
				{
					if (find(parent, k) == find(parent, j)) continue;
					assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
					size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
					if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches) break;
					size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
					if (mismatches > maxMismatches) continue;
					merge(parent, j, k);
				}
			}
		});
		for (size_t index = 1; index < numBlocks; index++)
		{
			size_t j = index * blockSize;
			size_t k = j-1;
			assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
			size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
			size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
			if (mismatches > maxMismatches) continue;
			merge(parent, j, k);
		}
		for (size_t i = 0; i < parent.size(); i++)
		{
			find(parent, i);
		}
		std::vector<size_t> activeBlocks;
		for (size_t i = 1; i < numBlocks; i++)
		{
			activeBlocks.emplace_back(i);
		}
		for (size_t blockDistance = 1; blockDistance < numBlocks; blockDistance++)
		{
			// std::cerr << "split part 2 dist " << blockDistance << std::endl;
			std::vector<std::pair<size_t, size_t>> merges;
			std::mutex mergeMutex;
			phmap::flat_hash_set<size_t> removeBlocks;
			iterateMultithreaded(0, activeBlocks.size(), numThreads, [&sequencesPerOccurrence, &parent, &merges, &mergeMutex, &activeBlocks, &removeBlocks, blockSize, blockDistance, i, mismatchFloor](size_t blockindex)
			{
				size_t index = activeBlocks[blockindex];
				phmap::flat_hash_map<size_t, size_t> extraParent;
				std::vector<std::pair<size_t, size_t>> mergesHere;
				for (size_t j = index * blockSize; j < (index+1) * blockSize && j < sequencesPerOccurrence.size(); j++)
				{
					for (size_t k = (index-blockDistance+1) * blockSize - 1; k < j && k >= (index-blockDistance) * blockSize; k--)
					{
						size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
						if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches)
						{
							if (j == index * blockSize && k == (index-blockDistance+1) * blockSize - 1)
							{
								assert(mergesHere.size() == 0);
								std::lock_guard<std::mutex> lock { mergeMutex };
								removeBlocks.emplace(index);
								return;
							}
							break;
						}
						if (parent[k] == parent[j]) continue;
						if (extraParent.count(parent[k]) == 0) extraParent[parent[k]] = parent[k];
						if (extraParent.count(parent[j]) == 0) extraParent[parent[j]] = parent[j];
						if (find(extraParent, parent[k]) == find(extraParent, parent[j])) continue;
						assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
						size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
						if (mismatches > maxMismatches) continue;
						merge(extraParent, parent[j], parent[k]);
						mergesHere.emplace_back(parent[j], parent[k]);
					}
				}
				{
					std::lock_guard<std::mutex> lock { mergeMutex };
					merges.insert(merges.end(), mergesHere.begin(), mergesHere.end());
				}
			});
			// std::cerr << merges.size() << " merges" << std::endl;
			for (auto pair : merges)
			{
				merge(parent, pair.first, pair.second);
			}
			if (merges.size() >= 1)
			{
				for (size_t i = 0; i < parent.size(); i++)
				{
					find(parent, i);
				}
			}
			removeBlocks.emplace(blockDistance);
			for (size_t i = activeBlocks.size()-1; i < activeBlocks.size(); i--)
			{
				if (removeBlocks.count(activeBlocks[i]) == 0) continue;
				std::swap(activeBlocks[i], activeBlocks.back());
				activeBlocks.pop_back();
			}
			// std::cerr << "remove " << removeBlocks.size() << ", " << activeBlocks.size() << " remaining" << std::endl;
			if (activeBlocks.size() == 0) break;
		}
		// auto endTime = getTime();
		phmap::flat_hash_map<size_t, size_t> keyToNode;
		for (size_t j = 0; j < parent.size(); j++)
		{
			size_t key = find(parent, j);
			if (keyToNode.count(key) == 1) continue;
			keyToNode[key] = nextNum;
			nextNum += 1;
		}
		// std::cerr << "sequence identity big chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " time " << formatTime(startTime, endTime) << " splitted to " << keyToNode.size() << std::endl;
		// std::cerr << "sequence identity splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			size_t index = sequencesPerOccurrence[j].second;
			std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
		}
	}
	std::mutex resultMutex;
	iterateMultithreaded(firstSmall, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &sequenceIndex, &iterationOrder, &rawReadLengths, mismatchFloor](const size_t iterationIndex)
	{
		const size_t i = iterationOrder[iterationIndex];
		// {
		// 	std::lock_guard<std::mutex> lock { resultMutex };
		// 	std::cerr << "sequence identity split chunk " << i << " coverage " << occurrencesPerChunk[i].size() << std::endl;
		// }
		// size_t coverage = occurrencesPerChunk[i].size();
		// size_t totalCheckablePairs = 0;
		// size_t totalMergedPairs = 0;
		// size_t totalNotMergedPairs = 0;
		// size_t totalDistinct = 0;
		// auto startTime = getTime();
		std::vector<std::pair<std::string, size_t>> sequencesPerOccurrence;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sequencesPerOccurrence.emplace_back();
			sequencesPerOccurrence.back().second = j;
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequencesPerOccurrence.back().first = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
		}
		std::sort(sequencesPerOccurrence.begin(), sequencesPerOccurrence.end(), [](const auto& left, const auto& right) {
			if (left.first.size() < right.first.size()) return true;
			if (left.first.size() > right.first.size()) return false;
			return left.first < right.first;
		});
		std::vector<bool> differentFromPredecessor;
		differentFromPredecessor.resize(occurrencesPerChunk[i].size(), false);
		differentFromPredecessor[0] = true;
		// totalDistinct += 1;
		for (size_t j = 1; j < differentFromPredecessor.size(); j++)
		{
			if (sequencesPerOccurrence[j].first == sequencesPerOccurrence[j-1].first) continue;
			differentFromPredecessor[j] = true;
			// totalDistinct += 1;
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
		{
			if (!differentFromPredecessor[j])
			{
				merge(parent, j, j-1);
				continue;
			}
			for (size_t k = j-1; k < occurrencesPerChunk[i].size(); k--)
			{
				if (!differentFromPredecessor[k]) continue;
				// totalCheckablePairs += 1;
				if (find(parent, k) == find(parent, j)) continue;
				assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
				size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
				if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches) break;
				size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
				if (mismatches > maxMismatches)
				{
					// totalNotMergedPairs += 1;
					continue;
				}
				// totalMergedPairs += 1;
				merge(parent, j, k);
			}
		}
		// auto endTime = getTime();
		phmap::flat_hash_map<size_t, size_t> keyToNode;
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
			// std::cerr << "sequence identity chunk " << i << " coverage " << coverage << " distinct " << totalDistinct << " checkable pairs " << totalCheckablePairs << " merged " << totalMergedPairs << " not merged " << totalNotMergedPairs << " time " << formatTime(startTime, endTime) << " splitted to " << keyToNode.size() << std::endl;
			// std::cerr << "sequence identity splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t index = sequencesPerOccurrence[j].second;
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
}

std::vector<size_t> getMinHashes(const std::string& sequence, const size_t k, const size_t count)
{
	assert(sequence.size() > k + count);
	uint64_t kmer = 0;
	uint64_t mask = (1ull << 2ull*k) - 1;
	for (size_t i = 0; i < k; i++)
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
	std::priority_queue<size_t> queue;
	queue.emplace(hashFwAndBw(kmer, k));
	for (size_t i = k; i < sequence.size(); i++)
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
		size_t h = hashFwAndBw(kmer, k);
		if (queue.size() < count)
		{
			queue.emplace(h);
		}
		else if (h < queue.top())
		{
			queue.pop();
			queue.emplace(h);
		}
	}
	std::vector<size_t> result;
	while (queue.size() > 0)
	{
		result.emplace_back(queue.top());
		queue.pop();
	}
	return result;
}

void splitPerFirstLastKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	assert(kmerSize <= 31);
	size_t nextNum = 0;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> endKmersToNumber;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		std::string readseq = getReadSequence(sequenceIndex, i);
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<0>(chunksPerRead[i][j]) < std::get<1>(chunksPerRead[i][j]));
			assert(std::get<1>(chunksPerRead[i][j]) < readseq.size());
			uint64_t firstKmer = 0;
			uint64_t lastKmer = 0;
			for (size_t k = 0; k < kmerSize; k++)
			{
				firstKmer <<= 2;
				lastKmer <<= 2;
				switch(readseq[std::get<0>(chunksPerRead[i][j])+k])
				{
					case 'A':
						firstKmer += 0;
						break;
					case 'C':
						firstKmer += 1;
						break;
					case 'G':
						firstKmer += 2;
						break;
					case 'T':
						firstKmer += 3;
						break;
					default:
						assert(false);
				}
				switch(readseq[std::get<1>(chunksPerRead[i][j])-k])
				{
					case 'A':
						lastKmer += 3;
						break;
					case 'C':
						lastKmer += 2;
						break;
					case 'G':
						lastKmer += 1;
						break;
					case 'T':
						lastKmer += 0;
						break;
					default:
						assert(false);
				}
			}
			auto key = std::make_pair(firstKmer, lastKmer);
			if (endKmersToNumber.count(key) == 0)
			{
				endKmersToNumber[key] = nextNum;
				nextNum += 1;
			}
			std::get<2>(chunksPerRead[i][j]) = endKmersToNumber.at(key) + firstBitUint64_t;
		}
	}
}

void splitPerBaseCounts(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t mismatchFloor = 10;
	std::cerr << "splitting by base counts" << std::endl;
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, mismatchFloor](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		std::vector<std::vector<size_t>> countsPerOccurrence;
		countsPerOccurrence.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			std::string chunkSequence = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
			assert(chunkSequence.size() == std::get<1>(t)-std::get<0>(t)+1);
			countsPerOccurrence[j].resize(4);
			for (size_t k = 0; k < chunkSequence.size(); k++)
			{
				switch(chunkSequence[k])
				{
				case 'A':
					countsPerOccurrence[j][0] += 1;
					break;
				case 'C':
					countsPerOccurrence[j][1] += 1;
					break;
				case 'G':
					countsPerOccurrence[j][2] += 1;
					break;
				case 'T':
					countsPerOccurrence[j][3] += 1;
					break;
				default:
					assert(false);
				}
			}
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				size_t maxEdits = std::min(std::get<1>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) - std::get<0>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]), std::get<1>(chunksPerRead[occurrencesPerChunk[i][k].first][occurrencesPerChunk[i][k].second]) - std::get<0>(chunksPerRead[occurrencesPerChunk[i][k].first][occurrencesPerChunk[i][k].second]));
				maxEdits *= mismatchFraction;
				maxEdits = std::max(maxEdits, mismatchFloor);
				size_t edits = 0;
				for (size_t c = 0; c < 4; c++)
				{
					if (countsPerOccurrence[j][c] > countsPerOccurrence[k][c])
					{
						edits += countsPerOccurrence[j][c] - countsPerOccurrence[k][c];
					}
					else
					{
						edits += countsPerOccurrence[k][c] - countsPerOccurrence[j][c];
					}
				}
				if (edits < maxEdits) merge(parent, j, k);
			}
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "base count splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
}

void splitPerInterchunkPhasedKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize)
{
	std::cerr << "splitting by interchunk phasing kmers" << std::endl;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t) + 1);
			occurrencesPerChunk[std::get<2>(t)].emplace_back(i, j);
		}
	}
	std::vector<bool> repetitive;
	repetitive.resize(occurrencesPerChunk.size(), false);
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		phmap::flat_hash_set<size_t> readsHere;
		size_t repeats = 0;
		for (auto t : occurrencesPerChunk[i])
		{
			if (readsHere.count(t.first) == 1) repeats += 1;
			readsHere.insert(t.first);
		}
		if (repeats >= 2) repetitive[i] = true;
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint64_t>> edges;
	{
		phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, size_t>> edgeCoverage;
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 1; j < chunksPerRead[i].size(); j++)
			{
				uint64_t prevNode = std::get<2>(chunksPerRead[i][j-1]);
				uint64_t thisNode = std::get<2>(chunksPerRead[i][j]);
				if (NonexistantChunk(prevNode)) continue;
				if (NonexistantChunk(thisNode)) continue;
				if (repetitive[prevNode & maskUint64_t]) continue;
				if (repetitive[thisNode & maskUint64_t]) continue;
				if (occurrencesPerChunk[prevNode & maskUint64_t].size() == 1) continue;
				if (occurrencesPerChunk[thisNode & maskUint64_t].size() == 1) continue;
				edgeCoverage[prevNode][thisNode] += 1;
				edgeCoverage[thisNode ^ firstBitUint64_t][prevNode ^ firstBitUint64_t] += 1;
			}
		}
		for (const auto& pair : edgeCoverage)
		{
			for (const auto pair2 : pair.second)
			{
				if (pair2.second == 1) continue;
				edges[pair.first].emplace(pair2.first);
			}
		}
	}
	std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> maybePhaseGroups;
	{
		phmap::flat_hash_map<uint64_t, std::tuple<size_t, uint64_t, uint64_t>> edgeToGroupMapping;
		size_t nextPhaseGroup = 0;
		for (const auto& pair : edges)
		{
			if (pair.second.size() != 2) continue;
			assert(!NonexistantChunk(pair.first));
			assert(!repetitive[pair.first & maskUint64_t]);
			uint64_t first = *pair.second.begin();
			uint64_t second = *(++pair.second.begin());
			assert(!NonexistantChunk(first));
			assert(!NonexistantChunk(second));
			assert(!repetitive[first & maskUint64_t]);
			assert(!repetitive[second & maskUint64_t]);
			assert(edges.count(first ^ firstBitUint64_t) == 1);
			assert(edges.count(second ^ firstBitUint64_t) == 1);
			if (edges.at(first ^ firstBitUint64_t).size() != 1) continue;
			if (edges.at(second ^ firstBitUint64_t).size() != 1) continue;
			assert(first != second);
			if ((first & maskUint64_t) == (second & maskUint64_t)) continue;
			assert(*edges.at(first ^ firstBitUint64_t).begin() == (pair.first ^ firstBitUint64_t));
			assert(*edges.at(second ^ firstBitUint64_t).begin() == (pair.first ^ firstBitUint64_t));
			edgeToGroupMapping[pair.first] = std::make_tuple(nextPhaseGroup, first, second);
			nextPhaseGroup += 1;
			maybePhaseGroups.emplace_back();
		}
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 1; j < chunksPerRead[i].size(); j++)
			{
				uint64_t prevNode = std::get<2>(chunksPerRead[i][j-1]);
				uint64_t thisNode = std::get<2>(chunksPerRead[i][j]);
				if (NonexistantChunk(prevNode)) continue;
				if (NonexistantChunk(thisNode)) continue;
				uint64_t allele = std::numeric_limits<size_t>::max();
				uint64_t fork = std::numeric_limits<size_t>::max();
				if (edgeToGroupMapping.count(prevNode) == 1 && edgeToGroupMapping.count(thisNode ^ firstBitUint64_t) == 1) continue;
				if (edgeToGroupMapping.count(prevNode) == 1)
				{
					assert(edgeToGroupMapping.count(thisNode ^ firstBitUint64_t) == 0);
					fork = prevNode;
					allele = thisNode;
				}
				if (edgeToGroupMapping.count(thisNode ^ firstBitUint64_t) == 1)
				{
					assert(edgeToGroupMapping.count(prevNode) == 0);
					fork = (thisNode ^ firstBitUint64_t);
					allele = (prevNode ^ firstBitUint64_t);
				}
				if (fork == std::numeric_limits<size_t>::max()) continue;
				assert(allele != std::numeric_limits<size_t>::max());
				assert(allele != std::get<1>(edgeToGroupMapping.at(fork)) || allele != std::get<2>(edgeToGroupMapping.at(fork)));
				if (allele == std::get<1>(edgeToGroupMapping.at(fork)))
				{
					maybePhaseGroups[std::get<0>(edgeToGroupMapping.at(fork))].first.emplace(i);
				}
				else if (allele == std::get<2>(edgeToGroupMapping.at(fork)))
				{
					maybePhaseGroups[std::get<0>(edgeToGroupMapping.at(fork))].second.emplace(i);
				}
			}
		}
		for (size_t i = maybePhaseGroups.size()-1; i < maybePhaseGroups.size(); i--)
		{
			phmap::flat_hash_set<size_t> inconsistents;
			for (size_t read : maybePhaseGroups[i].first)
			{
				if (maybePhaseGroups[i].second.count(read) == 1) inconsistents.insert(read);
			}
			if (inconsistents.size() >= 3)
			{
				std::swap(maybePhaseGroups[i], maybePhaseGroups.back());
				maybePhaseGroups.pop_back();
				continue;
			}
			for (size_t read : inconsistents)
			{
				maybePhaseGroups[i].first.erase(read);
				maybePhaseGroups[i].second.erase(read);
			}
		}
		for (size_t i = maybePhaseGroups.size()-1; i < maybePhaseGroups.size(); i--)
		{
			if (maybePhaseGroups[i].first.size() >= 4 && maybePhaseGroups[i].second.size() >= 4) continue;
			std::swap(maybePhaseGroups.back(), maybePhaseGroups[i]);
			maybePhaseGroups.pop_back();
		}
	}
	size_t nextNum = 0;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &sequenceIndex, &repetitive, &maybePhaseGroups, &countSplitted, &rawReadLengths, kmerSize](const size_t i)
	{
		std::vector<size_t> applicablePhasingGroups;
		std::vector<size_t> readsHere;
		{
			phmap::flat_hash_set<size_t> readsHereUnique;
			for (auto t : occurrencesPerChunk[i])
			{
				readsHereUnique.emplace(t.first);
			}
			readsHere.insert(readsHere.end(), readsHereUnique.begin(), readsHereUnique.end());
		}
		for (size_t j = 0; j < maybePhaseGroups.size(); j++)
		{
			size_t firsts = 0;
			size_t seconds = 0;
			for (size_t read : readsHere)
			{
				if (maybePhaseGroups[j].first.count(read) == 1) firsts += 1;
				if (maybePhaseGroups[j].second.count(read) == 1) seconds += 1;
			}
			if (firsts < 3) continue;
			if (seconds < 3) continue;
			if (firsts+seconds < readsHere.size() * 0.75 && (firsts < 5 || seconds < 5)) continue;
			if (firsts+seconds < readsHere.size() * 0.25) continue;
			applicablePhasingGroups.push_back(j);
		}
		if (applicablePhasingGroups.size() == 0)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		phmap::flat_hash_map<size_t, std::vector<size_t>> phaseGroupsPerRead;
		size_t nextGroupNum = 0;
		std::vector<phmap::flat_hash_set<size_t>> solidKmersPerOccurrence;
		solidKmersPerOccurrence.resize(occurrencesPerChunk[i].size());
		phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerClusterToNumber;
		std::vector<TwobitString> chunkSequences;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
		}
		iterateSolidKmers(chunkSequences, kmerSize, occurrencesPerChunk[i].size(), true, false, [&solidKmersPerOccurrence, &kmerClusterToNumber](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
		{
			size_t kmerkey = 0;
			std::pair<size_t, size_t> key { kmer, clusterIndex };
			if (kmerClusterToNumber.count(key) == 0)
			{
				kmerkey = kmerClusterToNumber.size();
				kmerClusterToNumber[key] = kmerkey;
			}
			else
			{
				kmerkey = kmerClusterToNumber.at(key);
			}
			solidKmersPerOccurrence[occurrenceID].emplace(kmerkey);
		});
		phmap::flat_hash_map<size_t, size_t> kmerToNumber;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			phmap::flat_hash_set<size_t> kmersHere;
			iterateKmers(chunkSequences[j], 0, chunkSequences[j].size()-1, true, kmerSize, [&kmersHere](const size_t kmer, const size_t pos)
			{
				kmersHere.insert(kmer);
			});
			for (auto kmer : kmersHere)
			{
				size_t kmerkey = 0;
				if (kmerToNumber.count(kmer) == 0)
				{
					kmerkey = kmerToNumber.size() + kmerClusterToNumber.size();
					kmerToNumber[kmer] = kmerkey;
				}
				else
				{
					kmerkey = kmerToNumber.at(kmer);
				}
				solidKmersPerOccurrence[j].emplace(kmerkey);
			}
		}
		{
			std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>> alleles; // startkmer, startcluster, endkmer, endcluster, occurrenceID, allele
			size_t lastOccurrence = std::numeric_limits<size_t>::max();
			size_t lastKmer = std::numeric_limits<size_t>::max();
			size_t lastCluster = std::numeric_limits<size_t>::max();
			size_t lastKmerPos = std::numeric_limits<size_t>::max();
			phmap::flat_hash_map<std::string, size_t> alleleNumbers;
			iterateSolidKmers(chunkSequences, kmerSize, occurrencesPerChunk[i].size(), false, false, [&occurrencesPerChunk, &chunkSequences, &alleles, &lastOccurrence, &lastKmer, &lastCluster, &lastKmerPos, &alleleNumbers, kmerSize, i](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				if (occurrenceID != lastOccurrence || lastKmer == std::numeric_limits<size_t>::max())
				{
					lastOccurrence = occurrenceID;
					lastKmer = kmer;
					lastCluster = clusterIndex;
					lastKmerPos = pos;
					return;
				}
				size_t allele;
				allele = getAllele(chunkSequences[occurrenceID], lastKmerPos, pos + kmerSize, alleleNumbers, true);
				alleles.emplace_back(lastKmer, lastCluster, kmer, clusterIndex, occurrenceID, allele);
				lastOccurrence = occurrenceID;
				lastKmer = kmer;
				lastCluster = clusterIndex;
				lastKmerPos = pos;
			});
			std::sort(alleles.begin(), alleles.end());
			std::vector<std::vector<std::vector<size_t>>> occurrencesPerAlleleSite;
			occurrencesPerAlleleSite = getAlleleOccurrences(alleles, occurrencesPerChunk[i].size());
			size_t nextNum = kmerToNumber.size() + kmerClusterToNumber.size();
			for (size_t j = 0; j < occurrencesPerAlleleSite.size(); j++)
			{
				if (occurrencesPerAlleleSite[j].size() < 2) continue;
				for (size_t allele = 0; allele < occurrencesPerAlleleSite[j].size(); allele++)
				{
					for (size_t k : occurrencesPerAlleleSite[j][allele])
					{
						solidKmersPerOccurrence[k].emplace(nextNum);
					}
					nextNum += 1;
				}
			}
		}
		{
			phmap::flat_hash_map<size_t, size_t> kmerCoverage;
			for (size_t j = 0; j < solidKmersPerOccurrence.size(); j++)
			{
				for (auto kmer : solidKmersPerOccurrence[j])
				{
					kmerCoverage[kmer] += 1;
				}
			}
			phmap::flat_hash_set<size_t> removeKmers;
			for (auto pair : kmerCoverage)
			{
				if (pair.second < 3 || pair.second + 3 > occurrencesPerChunk[i].size())
				{
					removeKmers.insert(pair.first);
				}
			}
			if (removeKmers.size() >= 1)
			{
				for (size_t j = 0; j < solidKmersPerOccurrence.size(); j++)
				{
					for (auto kmer : removeKmers)
					{
						if (solidKmersPerOccurrence[j].count(kmer) == 1)
						{
							solidKmersPerOccurrence[j].erase(kmer);
						}
					}
				}
			}
		}
		for (const size_t phaseGroup : applicablePhasingGroups)
		{
			phmap::flat_hash_set<size_t> kmersInFirst;
			phmap::flat_hash_set<size_t> kmersInSecond;
			phmap::flat_hash_set<size_t> kmersInEvenOneFirst;
			phmap::flat_hash_set<size_t> kmersInEvenOneSecond;
			bool hasFirst = false;
			bool hasSecond = false;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 0 && maybePhaseGroups[phaseGroup].second.count(read) == 0) continue;
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				const phmap::flat_hash_set<size_t> kmersHere = solidKmersPerOccurrence[j];
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1)
				{
					kmersInEvenOneFirst.insert(kmersHere.begin(), kmersHere.end());
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
					if (!hasFirst)
					{
						hasFirst = true;
						kmersInFirst = kmersHere;
					}
					else
					{
						phmap::flat_hash_set<size_t> removeKmers;
						for (size_t kmer : kmersInFirst)
						{
							if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
						}
						for (size_t kmer : removeKmers)
						{
							kmersInFirst.erase(kmer);
						}
					}
				}
				else
				{
					kmersInEvenOneSecond.insert(kmersHere.begin(), kmersHere.end());
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 1);
					if (!hasSecond)
					{
						hasSecond = true;
						kmersInSecond = kmersHere;
					}
					else
					{
						phmap::flat_hash_set<size_t> removeKmers;
						for (size_t kmer : kmersInSecond)
						{
							if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
						}
						for (size_t kmer : removeKmers)
						{
							kmersInSecond.erase(kmer);
						}
					}
				}
			}
			assert(hasFirst);
			assert(hasSecond);
			{
				for (size_t kmer : kmersInEvenOneFirst)
				{
					if (kmersInSecond.count(kmer) == 0) continue;
					kmersInSecond.erase(kmer);
				}
				for (size_t kmer : kmersInEvenOneSecond)
				{
					if (kmersInFirst.count(kmer) == 0) continue;
					kmersInFirst.erase(kmer);
				}
			}
			{
				phmap::flat_hash_set<size_t> sharedKmers;
				for (auto kmer : kmersInFirst)
				{
					if (kmersInSecond.count(kmer) == 0) continue;
					sharedKmers.insert(kmer);
				}
				for (auto kmer : sharedKmers)
				{
					kmersInFirst.erase(kmer);
					kmersInSecond.erase(kmer);
				}
			}
			if (kmersInFirst.size() == 0 || kmersInSecond.size() == 0) continue;
			bool possiblyValid = true;
			phmap::flat_hash_set<size_t> removeKmers;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1 || maybePhaseGroups[phaseGroup].second.count(read) == 1) continue;
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				const phmap::flat_hash_set<size_t>& kmersHere = solidKmersPerOccurrence[j];
				size_t firstMatches = 0;
				size_t secondMatches = 0;
				for (size_t kmer : kmersHere)
				{
					if (kmersInFirst.count(kmer) == 1)
					{
						assert(kmersInSecond.count(kmer) == 0);
						firstMatches += 1;
					}
					if (kmersInSecond.count(kmer) == 1)
					{
						assert(kmersInFirst.count(kmer) == 0);
						secondMatches += 1;
					}
				}
				if (firstMatches == secondMatches)
				{
					possiblyValid = false;
					break;
				}
				if (firstMatches > secondMatches)
				{
					for (size_t kmer : kmersInFirst)
					{
						if (kmersHere.count(kmer) == 1) continue;
						removeKmers.insert(kmer);
					}
					for (size_t kmer : kmersInSecond)
					{
						if (kmersHere.count(kmer) == 0) continue;
						removeKmers.insert(kmer);
					}
				}
				else
				{
					assert(firstMatches < secondMatches);
					for (size_t kmer : kmersInSecond)
					{
						if (kmersHere.count(kmer) == 1) continue;
						removeKmers.insert(kmer);
					}
					for (size_t kmer : kmersInFirst)
					{
						if (kmersHere.count(kmer) == 0) continue;
						removeKmers.insert(kmer);
					}
				}
			}
			if (!possiblyValid) continue;
			for (auto kmer : removeKmers)
			{
				assert(kmersInFirst.count(kmer) == 1 || kmersInSecond.count(kmer) == 1);
				if (kmersInFirst.count(kmer) == 1)
				{
					assert(kmersInSecond.count(kmer) == 0);
					kmersInFirst.erase(kmer);
				}
				else
				{
					assert(kmersInSecond.count(kmer) == 1);
					kmersInSecond.erase(kmer);
				}
			}
			phmap::flat_hash_map<size_t, size_t> assignments;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1)
				{
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
					assignments[read] = 0;
					continue;
				}
				if (maybePhaseGroups[phaseGroup].second.count(read) == 1)
				{
					assignments[read] = 1;
					continue;
				}
				assert(maybePhaseGroups[phaseGroup].first.count(read) == 0);
				assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				const phmap::flat_hash_set<size_t>& kmersHere = solidKmersPerOccurrence[j];
				size_t firstMatches = 0;
				size_t secondMatches = 0;
				for (size_t kmer : kmersHere)
				{
					if (kmersInFirst.count(kmer) == 1)
					{
						assert(kmersInSecond.count(kmer) == 0);
						firstMatches += 1;
					}
					if (kmersInSecond.count(kmer) == 1)
					{
						assert(kmersInFirst.count(kmer) == 0);
						secondMatches += 1;
					}
				}
				if (firstMatches == secondMatches)
				{
					possiblyValid = false;
					break;
				}
				if (firstMatches > secondMatches)
				{
					if (assignments.count(read) == 1 && assignments.at(read) != 0)
					{
						possiblyValid = false;
						break;
					}
					assert(assignments.count(read) == 0 || assignments.at(read) == 0);
					assignments[read] = 0;
				}
				else
				{
					if (assignments.count(read) == 1 && assignments.at(read) != 1)
					{
						possiblyValid = false;
						break;
					}
					assert(assignments.count(read) == 0 || assignments.at(read) == 1);
					assignments[read] = 1;
				}
			}
			if (!possiblyValid) continue;
			for (auto pair : assignments)
			{
				phaseGroupsPerRead[pair.first].push_back(nextGroupNum + pair.second);
			}
			nextGroupNum += 2;
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		if (phaseGroupsPerRead.size() == 0)
		{
			for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
			{
				merge(parent, 0, j);
			}
		}
		else
		{
			for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					size_t readj = occurrencesPerChunk[i][j].first;
					size_t readk = occurrencesPerChunk[i][k].first;
					if (phaseGroupsPerRead.at(readj) != phaseGroupsPerRead.at(readk)) continue;
					merge(parent, j, k);
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> clusterToNode;
			// size_t first = nextNum;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				if (clusterToNode.count(find(parent, j)) == 1) continue;
				clusterToNode[find(parent, j)] = nextNum;
				nextNum += 1;
			}
			// size_t last = nextNum-1;
			if (clusterToNode.size() >= 2) countSplitted += 1;
//			std::cerr << "interchunk phasing kmers splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks range " << first << " - " << last << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, j)) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
	std::cerr << "interchunk phasing kmers splitted " << countSplitted << " chunks" << std::endl;
}

void splitPerMinHashes(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t kmerSize = 11;
	const size_t hashCount = 10;
	std::cerr << "splitting by minhash" << std::endl;
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, kmerSize, hashCount](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		phmap::flat_hash_map<size_t, size_t> parent;
		std::vector<size_t> oneHashPerLocation;
		size_t smallestSequence = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			smallestSequence = std::min(smallestSequence, std::get<1>(t)-std::get<0>(t)+1);
		}
		if (smallestSequence <= kmerSize + hashCount + 1)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			std::vector<uint64_t> minHashes;
			auto chunkSequence = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
			minHashes = getMinHashes(chunkSequence, kmerSize, hashCount);
			assert(minHashes.size() >= 1);
			for (auto hash : minHashes)
			{
				if (parent.count(hash) == 0) parent[hash] = hash;
			}
			for (auto hash : minHashes)
			{
				merge(parent, hash, minHashes[0]);
			}
			oneHashPerLocation.emplace_back(minHashes[0]);
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> clusterToNode;
			for (auto hash : oneHashPerLocation)
			{
				if (clusterToNode.count(find(parent, hash)) == 1) continue;
				clusterToNode[find(parent, hash)] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "minhash splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, oneHashPerLocation[j])) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
}

void splitPerDiploidChunkWithNeighbors(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "splitting by diploid chunk with neighbors" << std::endl;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			maxChunk = std::max(maxChunk, std::get<2>(t) & maskUint64_t);
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t) + 1);
			occurrencesPerChunk[std::get<2>(t)].emplace_back(i, j);
		}
	}
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<size_t> chunkBelongsToUnitig;
	chunkBelongsToUnitig.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			assert((node & maskUint64_t) < chunkBelongsToUnitig.size());
			assert(chunkBelongsToUnitig[node & maskUint64_t] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[node & maskUint64_t] = i;
		}
	}
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> fwForksPerUnitig;
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> bwForksPerUnitig;
	fwForksPerUnitig.resize(graph.unitigLengths.size());
	bwForksPerUnitig.resize(graph.unitigLengths.size());
	std::vector<bool> hasFwFork;
	std::vector<bool> hasBwFork;
	hasFwFork.resize(graph.unitigLengths.size(), false);
	hasBwFork.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		if (graph.edges.getEdges(fw).size() == 2)
		{
			hasFwFork[i] = true;
		}
		std::pair<size_t, bool> bw { i, false };
		if (graph.edges.getEdges(bw).size() == 2)
		{
			hasBwFork[i] = true;
		}
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 1; k < readPaths[i][j].path.size(); k++)
			{
				uint64_t prev = readPaths[i][j].path[k-1];
				uint64_t curr = readPaths[i][j].path[k];
				if ((prev & firstBitUint64_t) && hasFwFork[prev & maskUint64_t])
				{
					std::pair<size_t, bool> currPair { curr & maskUint64_t, curr & firstBitUint64_t };
					if (currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, true))[0])
					{
						fwForksPerUnitig[prev & maskUint64_t].first.emplace_back(i);
					}
					else
					{
						assert(currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, true))[1]);
						fwForksPerUnitig[prev & maskUint64_t].second.emplace_back(i);
					}
				}
				if (((prev ^ firstBitUint64_t) & firstBitUint64_t) && hasBwFork[prev & maskUint64_t])
				{
					std::pair<size_t, bool> currPair { curr & maskUint64_t, curr & firstBitUint64_t };
					if (currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, false))[0])
					{
						bwForksPerUnitig[prev & maskUint64_t].first.emplace_back(i);
					}
					else
					{
						assert(currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, false))[1]);
						bwForksPerUnitig[prev & maskUint64_t].second.emplace_back(i);
					}
				}
				if ((curr & firstBitUint64_t) && hasBwFork[curr & maskUint64_t])
				{
					std::pair<size_t, bool> prevPair { prev & maskUint64_t, prev & firstBitUint64_t };
					prevPair.second = !prevPair.second;
					if (prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, false))[0])
					{
						bwForksPerUnitig[curr & maskUint64_t].first.emplace_back(i);
					}
					else
					{
						assert(prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, false))[1]);
						bwForksPerUnitig[curr & maskUint64_t].second.emplace_back(i);
					}
				}
				if (((curr ^ firstBitUint64_t) & firstBitUint64_t) && hasFwFork[curr & maskUint64_t])
				{
					std::pair<size_t, bool> prevPair { prev & maskUint64_t, prev & firstBitUint64_t };
					prevPair.second = !prevPair.second;
					if (prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, true))[0])
					{
						fwForksPerUnitig[curr & maskUint64_t].first.emplace_back(i);
					}
					else
					{
						assert(prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, true))[1]);
						fwForksPerUnitig[curr & maskUint64_t].second.emplace_back(i);
					}
				}
			}
		}
	}
	for (size_t i = 0; i < fwForksPerUnitig.size(); i++)
	{
		std::sort(fwForksPerUnitig[i].first.begin(), fwForksPerUnitig[i].first.end());
		std::sort(fwForksPerUnitig[i].second.begin(), fwForksPerUnitig[i].second.end());
		std::sort(bwForksPerUnitig[i].first.begin(), bwForksPerUnitig[i].first.end());
		std::sort(bwForksPerUnitig[i].second.begin(), bwForksPerUnitig[i].second.end());
		if (hasFwFork[i])
		{
			if (intersectSize(fwForksPerUnitig[i].first, fwForksPerUnitig[i].second) >= 1)
			{
				hasFwFork[i] = false;
			}
			if (phmap::flat_hash_set<size_t>(fwForksPerUnitig[i].first.begin(), fwForksPerUnitig[i].first.end()).size() != fwForksPerUnitig[i].first.size()) hasFwFork[i] = false;
			if (phmap::flat_hash_set<size_t>(fwForksPerUnitig[i].second.begin(), fwForksPerUnitig[i].second.end()).size() != fwForksPerUnitig[i].second.size()) hasFwFork[i] = false;
			if (fwForksPerUnitig[i].first.size() < 2) hasFwFork[i] = false;
			if (fwForksPerUnitig[i].second.size() < 2) hasFwFork[i] = false;
		}
		if (hasBwFork[i])
		{
			if (intersectSize(bwForksPerUnitig[i].first, bwForksPerUnitig[i].second) >= 1)
			{
				hasBwFork[i] = false;
			}
			if (phmap::flat_hash_set<size_t>(bwForksPerUnitig[i].first.begin(), bwForksPerUnitig[i].first.end()).size() != bwForksPerUnitig[i].first.size()) hasBwFork[i] = false;
			if (phmap::flat_hash_set<size_t>(bwForksPerUnitig[i].second.begin(), bwForksPerUnitig[i].second.end()).size() != bwForksPerUnitig[i].second.size()) hasBwFork[i] = false;
			if (bwForksPerUnitig[i].first.size() < 2) hasBwFork[i] = false;
			if (bwForksPerUnitig[i].second.size() < 2) hasBwFork[i] = false;
		}
	}
	std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> fwForksPerUnitigSet;
	std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> bwForksPerUnitigSet;
	fwForksPerUnitigSet.resize(graph.unitigLengths.size());
	bwForksPerUnitigSet.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < fwForksPerUnitig.size(); i++)
	{
		if (hasFwFork[i])
		{
			fwForksPerUnitigSet[i].first.insert(fwForksPerUnitig[i].first.begin(), fwForksPerUnitig[i].first.end());
			fwForksPerUnitigSet[i].second.insert(fwForksPerUnitig[i].second.begin(), fwForksPerUnitig[i].second.end());
		}
		if (hasBwFork[i])
		{
			bwForksPerUnitigSet[i].first.insert(bwForksPerUnitig[i].first.begin(), bwForksPerUnitig[i].first.end());
			bwForksPerUnitigSet[i].second.insert(bwForksPerUnitig[i].second.begin(), bwForksPerUnitig[i].second.end());
		}
	}
	size_t nextNum = 0;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &sequenceIndex, &hasFwFork, &hasBwFork, &fwForksPerUnitigSet, &bwForksPerUnitigSet,  &countSplitted, &chunkBelongsToUnitig, &rawReadLengths, kmerSize](const size_t i)
	{
		if (chunkBelongsToUnitig[i] == std::numeric_limits<size_t>::max())
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		size_t unitig = chunkBelongsToUnitig[i];
		bool anythingToPhase = false;
		std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> maybePhaseGroups;
		if (hasFwFork[unitig])
		{
			if (fwForksPerUnitigSet[unitig].first.size() >= 2 && fwForksPerUnitigSet[unitig].second.size() >= 2)
			{
				size_t firstsHere = 0;
				size_t secondsHere = 0;
				bool valid = true;
				for (auto pair : occurrencesPerChunk[i])
				{
					if (fwForksPerUnitigSet[unitig].first.count(pair.first) == 1)
					{
						firstsHere += 1;
						if (fwForksPerUnitigSet[unitig].second.count(pair.first) == 1)
						{
							valid = false;
							break;
						}
					}
					if (fwForksPerUnitigSet[unitig].second.count(pair.first) == 1)
					{
						secondsHere += 1;
					}
				}
				if (firstsHere >= 2 && secondsHere >= 2 && valid)
				{
					maybePhaseGroups.emplace_back(fwForksPerUnitigSet[unitig]);
				}
			}
		}
		if (hasBwFork[unitig])
		{
			if (bwForksPerUnitigSet[unitig].first.size() >= 2 && bwForksPerUnitigSet[unitig].second.size() >= 2)
			{
				size_t firstsHere = 0;
				size_t secondsHere = 0;
				bool valid = true;
				for (auto pair : occurrencesPerChunk[i])
				{
					if (bwForksPerUnitigSet[unitig].first.count(pair.first) == 1)
					{
						firstsHere += 1;
						if (bwForksPerUnitigSet[unitig].second.count(pair.first) == 1)
						{
							valid = false;
							break;
						}
					}
					if (bwForksPerUnitigSet[unitig].second.count(pair.first) == 1)
					{
						secondsHere += 1;
					}
				}
				if (firstsHere >= 2 && secondsHere >= 2 && valid)
				{
					maybePhaseGroups.emplace_back(bwForksPerUnitigSet[unitig]);
				}
			}
		}
		if (maybePhaseGroups.size() == 0)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		size_t nextGroupNum = 0;
		std::vector<phmap::flat_hash_set<size_t>> solidKmersPerOccurrence;
		solidKmersPerOccurrence.resize(occurrencesPerChunk[i].size());
		phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerClusterToNumber;
		std::vector<TwobitString> chunkSequences;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
		}
		iterateSolidKmers(chunkSequences, kmerSize, occurrencesPerChunk[i].size(), true, false, [&solidKmersPerOccurrence, &kmerClusterToNumber](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
		{
			size_t kmerkey = 0;
			std::pair<size_t, size_t> key { kmer, clusterIndex };
			if (kmerClusterToNumber.count(key) == 0)
			{
				kmerkey = kmerClusterToNumber.size();
				kmerClusterToNumber[key] = kmerkey;
			}
			else
			{
				kmerkey = kmerClusterToNumber.at(key);
			}
			solidKmersPerOccurrence[occurrenceID].emplace(kmerkey);
		});
		phmap::flat_hash_map<size_t, size_t> kmerToNumber;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			phmap::flat_hash_set<size_t> kmersHere;
			iterateKmers(chunkSequences[j], 0, chunkSequences[j].size()-1, true, kmerSize, [&kmersHere](const size_t kmer, const size_t pos)
			{
				kmersHere.insert(kmer);
			});
			for (auto kmer : kmersHere)
			{
				size_t kmerkey = 0;
				if (kmerToNumber.count(kmer) == 0)
				{
					kmerkey = kmerToNumber.size() + kmerClusterToNumber.size();
					kmerToNumber[kmer] = kmerkey;
				}
				else
				{
					kmerkey = kmerToNumber.at(kmer);
				}
				solidKmersPerOccurrence[j].emplace(kmerkey);
			}
		}
		{
			std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>> alleles; // startkmer, startcluster, endkmer, endcluster, occurrenceID, allele
			size_t lastOccurrence = std::numeric_limits<size_t>::max();
			size_t lastKmer = std::numeric_limits<size_t>::max();
			size_t lastCluster = std::numeric_limits<size_t>::max();
			size_t lastKmerPos = std::numeric_limits<size_t>::max();
			phmap::flat_hash_map<std::string, size_t> alleleNumbers;
			iterateSolidKmers(chunkSequences, kmerSize, occurrencesPerChunk[i].size(), false, false, [&occurrencesPerChunk, &chunkSequences, &alleles, &lastOccurrence, &lastKmer, &lastCluster, &lastKmerPos, &alleleNumbers, kmerSize, i](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				if (occurrenceID != lastOccurrence || lastKmer == std::numeric_limits<size_t>::max())
				{
					lastOccurrence = occurrenceID;
					lastKmer = kmer;
					lastCluster = clusterIndex;
					lastKmerPos = pos;
					return;
				}
				size_t allele;
				allele = getAllele(chunkSequences[occurrenceID], lastKmerPos, pos + kmerSize, alleleNumbers, true);
				alleles.emplace_back(lastKmer, lastCluster, kmer, clusterIndex, occurrenceID, allele);
				lastOccurrence = occurrenceID;
				lastKmer = kmer;
				lastCluster = clusterIndex;
				lastKmerPos = pos;
			});
			std::sort(alleles.begin(), alleles.end());
			std::vector<std::vector<std::vector<size_t>>> occurrencesPerAlleleSite;
			occurrencesPerAlleleSite = getAlleleOccurrences(alleles, occurrencesPerChunk[i].size());
			size_t nextNum = kmerToNumber.size() + kmerClusterToNumber.size();
			for (size_t j = 0; j < occurrencesPerAlleleSite.size(); j++)
			{
				if (occurrencesPerAlleleSite[j].size() < 2) continue;
				for (size_t allele = 0; allele < occurrencesPerAlleleSite[j].size(); allele++)
				{
					for (size_t k : occurrencesPerAlleleSite[j][allele])
					{
						solidKmersPerOccurrence[k].emplace(nextNum);
					}
					nextNum += 1;
				}
			}
		}
		{
			phmap::flat_hash_map<size_t, size_t> kmerCoverage;
			for (size_t j = 0; j < solidKmersPerOccurrence.size(); j++)
			{
				for (auto kmer : solidKmersPerOccurrence[j])
				{
					kmerCoverage[kmer] += 1;
				}
			}
			phmap::flat_hash_set<size_t> removeKmers;
			for (auto pair : kmerCoverage)
			{
				if (pair.second < 2 || pair.second + 2 > occurrencesPerChunk[i].size())
				{
					removeKmers.insert(pair.first);
				}
			}
			if (removeKmers.size() >= 1)
			{
				for (size_t j = 0; j < solidKmersPerOccurrence.size(); j++)
				{
					for (auto kmer : removeKmers)
					{
						if (solidKmersPerOccurrence[j].count(kmer) == 1)
						{
							solidKmersPerOccurrence[j].erase(kmer);
						}
					}
				}
			}
		}
		phmap::flat_hash_map<size_t, std::vector<size_t>> phaseGroupsPerRead;
		for (size_t phaseGroup = 0; phaseGroup < maybePhaseGroups.size(); phaseGroup++)
		{
			phmap::flat_hash_set<size_t> kmersInFirst;
			phmap::flat_hash_set<size_t> kmersInSecond;
			phmap::flat_hash_set<size_t> kmersInEvenOneFirst;
			phmap::flat_hash_set<size_t> kmersInEvenOneSecond;
			bool hasFirst = false;
			bool hasSecond = false;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 0 && maybePhaseGroups[phaseGroup].second.count(read) == 0) continue;
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				const phmap::flat_hash_set<size_t> kmersHere = solidKmersPerOccurrence[j];
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1)
				{
					kmersInEvenOneFirst.insert(kmersHere.begin(), kmersHere.end());
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
					if (!hasFirst)
					{
						hasFirst = true;
						kmersInFirst = kmersHere;
					}
					else
					{
						phmap::flat_hash_set<size_t> removeKmers;
						for (size_t kmer : kmersInFirst)
						{
							if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
						}
						for (size_t kmer : removeKmers)
						{
							kmersInFirst.erase(kmer);
						}
					}
				}
				else
				{
					kmersInEvenOneSecond.insert(kmersHere.begin(), kmersHere.end());
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 1);
					if (!hasSecond)
					{
						hasSecond = true;
						kmersInSecond = kmersHere;
					}
					else
					{
						phmap::flat_hash_set<size_t> removeKmers;
						for (size_t kmer : kmersInSecond)
						{
							if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
						}
						for (size_t kmer : removeKmers)
						{
							kmersInSecond.erase(kmer);
						}
					}
				}
			}
			assert(hasFirst);
			assert(hasSecond);
			{
				for (size_t kmer : kmersInEvenOneFirst)
				{
					if (kmersInSecond.count(kmer) == 0) continue;
					kmersInSecond.erase(kmer);
				}
				for (size_t kmer : kmersInEvenOneSecond)
				{
					if (kmersInFirst.count(kmer) == 0) continue;
					kmersInFirst.erase(kmer);
				}
			}
			{
				phmap::flat_hash_set<size_t> sharedKmers;
				for (auto kmer : kmersInFirst)
				{
					if (kmersInSecond.count(kmer) == 0) continue;
					sharedKmers.insert(kmer);
				}
				for (auto kmer : sharedKmers)
				{
					kmersInFirst.erase(kmer);
					kmersInSecond.erase(kmer);
				}
			}
			if (kmersInFirst.size() == 0 || kmersInSecond.size() == 0) continue;
			bool possiblyValid = true;
			phmap::flat_hash_set<size_t> removeKmers;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1 || maybePhaseGroups[phaseGroup].second.count(read) == 1) continue;
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				const phmap::flat_hash_set<size_t>& kmersHere = solidKmersPerOccurrence[j];
				size_t firstMatches = 0;
				size_t secondMatches = 0;
				for (size_t kmer : kmersHere)
				{
					if (kmersInFirst.count(kmer) == 1)
					{
						assert(kmersInSecond.count(kmer) == 0);
						firstMatches += 1;
					}
					if (kmersInSecond.count(kmer) == 1)
					{
						assert(kmersInFirst.count(kmer) == 0);
						secondMatches += 1;
					}
				}
				if (firstMatches == secondMatches)
				{
					possiblyValid = false;
					break;
				}
				if (firstMatches > secondMatches)
				{
					for (size_t kmer : kmersInFirst)
					{
						if (kmersHere.count(kmer) == 1) continue;
						removeKmers.insert(kmer);
					}
					for (size_t kmer : kmersInSecond)
					{
						if (kmersHere.count(kmer) == 0) continue;
						removeKmers.insert(kmer);
					}
				}
				else
				{
					assert(firstMatches < secondMatches);
					for (size_t kmer : kmersInSecond)
					{
						if (kmersHere.count(kmer) == 1) continue;
						removeKmers.insert(kmer);
					}
					for (size_t kmer : kmersInFirst)
					{
						if (kmersHere.count(kmer) == 0) continue;
						removeKmers.insert(kmer);
					}
				}
			}
			if (!possiblyValid) continue;
			for (auto kmer : removeKmers)
			{
				assert(kmersInFirst.count(kmer) == 1 || kmersInSecond.count(kmer) == 1);
				if (kmersInFirst.count(kmer) == 1)
				{
					assert(kmersInSecond.count(kmer) == 0);
					kmersInFirst.erase(kmer);
				}
				else
				{
					assert(kmersInSecond.count(kmer) == 1);
					kmersInSecond.erase(kmer);
				}
			}
			phmap::flat_hash_map<size_t, size_t> assignments;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1)
				{
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
					assignments[read] = 0;
					continue;
				}
				if (maybePhaseGroups[phaseGroup].second.count(read) == 1)
				{
					assignments[read] = 1;
					continue;
				}
				assert(maybePhaseGroups[phaseGroup].first.count(read) == 0);
				assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				const phmap::flat_hash_set<size_t>& kmersHere = solidKmersPerOccurrence[j];
				size_t firstMatches = 0;
				size_t secondMatches = 0;
				for (size_t kmer : kmersHere)
				{
					if (kmersInFirst.count(kmer) == 1)
					{
						assert(kmersInSecond.count(kmer) == 0);
						firstMatches += 1;
					}
					if (kmersInSecond.count(kmer) == 1)
					{
						assert(kmersInFirst.count(kmer) == 0);
						secondMatches += 1;
					}
				}
				if (firstMatches == secondMatches)
				{
					possiblyValid = false;
					break;
				}
				if (firstMatches > secondMatches)
				{
					if (assignments.count(read) == 1 && assignments.at(read) != 0)
					{
						possiblyValid = false;
						break;
					}
					assert(assignments.count(read) == 0 || assignments.at(read) == 0);
					assignments[read] = 0;
				}
				else
				{
					if (assignments.count(read) == 1 && assignments.at(read) != 1)
					{
						possiblyValid = false;
						break;
					}
					assert(assignments.count(read) == 0 || assignments.at(read) == 1);
					assignments[read] = 1;
				}
			}
			if (!possiblyValid) continue;
			for (auto pair : assignments)
			{
				phaseGroupsPerRead[pair.first].push_back(nextGroupNum + pair.second);
			}
			nextGroupNum += 2;
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		if (phaseGroupsPerRead.size() == 0)
		{
			for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
			{
				merge(parent, 0, j);
			}
		}
		else
		{
			for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					size_t readj = occurrencesPerChunk[i][j].first;
					size_t readk = occurrencesPerChunk[i][k].first;
					if (phaseGroupsPerRead.at(readj) != phaseGroupsPerRead.at(readk)) continue;
					merge(parent, j, k);
				}
			}
		}
		phmap::flat_hash_map<size_t, size_t> groupCoverage;
		for (size_t i = 0; i < parent.size(); i++)
		{
			groupCoverage[find(parent, i)] += 1;
		}
		bool valid = true;
		if (groupCoverage.size() != 2) valid = false;
		for (auto pair : groupCoverage)
		{
			if (pair.second < 3) valid = false;
		}
		if (!valid)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> clusterToNode;
			size_t first = nextNum;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				if (clusterToNode.count(find(parent, j)) == 1) continue;
				clusterToNode[find(parent, j)] = nextNum;
				nextNum += 1;
			}
			size_t last = nextNum-1;
			if (clusterToNode.size() >= 2) countSplitted += 1;
			std::cerr << "diploid chunk with neighbors splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks range " << first << " - " << last << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, j)) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
	std::cerr << "diploid chunk with neighbors splitted " << countSplitted << " chunks" << std::endl;
}

void resolveVerySmallNodes(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const size_t maxResolveSize, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving very small nodes" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	phmap::flat_hash_set<size_t> resolvableTinies;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.chunksInUnitig[i].size() != 1) continue;
		if (graph.unitigLengths[i] >= maxResolveSize) continue;
		if (graph.edges.getEdges(std::make_pair(i, true)).size() < 2 && graph.edges.getEdges(std::make_pair(i, false)).size() < 2) continue;
		resolvableTinies.insert(graph.chunksInUnitig[i][0] & maskUint64_t);
	}
	size_t nextNum = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (!NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) nextNum = std::max(nextNum, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			if (resolvableTinies.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			if (std::get<0>(chunksPerRead[i][j]) == 0)
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (std::get<1>(chunksPerRead[i][j])+1 == rawReadLengths[i])
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1])+1)
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (j+1 < chunksPerRead[i].size() && std::get<1>(chunksPerRead[i][j])+1 == std::get<1>(chunksPerRead[i][j-1]))
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
		}
	}
	for (size_t chunk : resolvableTinies)
	{
		std::cerr << "resolve tiny chunk " << chunk << std::endl;
	}
	nextNum += 1;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeToNumber;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<size_t> eraseIndices;
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunks;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (resolvableTinies.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			eraseIndices.emplace_back(j);
			if (j > 0)
			{
				uint64_t here = std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t;
				uint64_t neighbor = std::get<2>(chunksPerRead[i][j-1]) ^ firstBitUint64_t;
				if (!NonexistantChunk(neighbor))
				{
					std::pair<uint64_t, uint64_t> key { here, neighbor };
					if (edgeToNumber.count(key) == 0)
					{
						edgeToNumber[key] = nextNum;
						nextNum += 1;
					}
					uint64_t node = edgeToNumber.at(key);
					newChunks.emplace_back(std::get<0>(chunksPerRead[i][j])-1, std::get<1>(chunksPerRead[i][j]), node);
				}
			}
			if (j+1 < chunksPerRead[i].size())
			{
				uint64_t here = std::get<2>(chunksPerRead[i][j]);
				uint64_t neighbor = std::get<2>(chunksPerRead[i][j+1]);
				if (!NonexistantChunk(neighbor))
				{
					std::pair<uint64_t, uint64_t> key { here, neighbor };
					if (edgeToNumber.count(key) == 0)
					{
						edgeToNumber[key] = nextNum;
						nextNum += 1;
					}
					uint64_t node = edgeToNumber.at(key) + firstBitUint64_t;
					newChunks.emplace_back(std::get<0>(chunksPerRead[i][j]), std::get<1>(chunksPerRead[i][j])+1, node);
				}
			}
		}
		if (eraseIndices.size() == 0 && newChunks.size() == 0) continue;
		for (size_t j = eraseIndices.size()-1; j < eraseIndices.size(); j--)
		{
			size_t index = eraseIndices[j];
			chunksPerRead[i].erase(chunksPerRead[i].begin()+index);
		}
		chunksPerRead[i].insert(chunksPerRead[i].end(), newChunks.begin(), newChunks.end());
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
	}
}

void resolveSemiAmbiguousUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving semiambiguous resolvable unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> tripletsPerUnitig;
	tripletsPerUnitig.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				uint64_t prev, after;
				if (mid & firstBitUint64_t)
				{
					prev = readPaths[i][j].path[k-1];
					after = readPaths[i][j].path[k+1];
				}
				else
				{
					prev = readPaths[i][j].path[k+1] ^ firstBitUint64_t;
					after = readPaths[i][j].path[k-1] ^ firstBitUint64_t;
				}
				tripletsPerUnitig[mid & maskUint64_t][std::make_pair(prev, after)] += 1;
			}
		}
	}
	std::vector<bool> canResolveUnitig;
	canResolveUnitig.resize(graph.unitigLengths.size(), false);
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> prevToAllele;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> afterToAllele;
	prevToAllele.resize(graph.unitigLengths.size());
	afterToAllele.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (tripletsPerUnitig[i].size() < 2) continue;
		phmap::flat_hash_set<uint64_t> coveredPrev;
		phmap::flat_hash_set<uint64_t> coveredAfter;
		phmap::flat_hash_map<uint64_t, uint64_t> parent;
		bool valid = true;
		for (auto triplet : tripletsPerUnitig[i])
		{
			if (triplet.second < 2)
			{
				valid = false;
				break;
			}
			coveredPrev.insert(triplet.first.first);
			coveredAfter.insert(triplet.first.second);
			assert((triplet.first.first & (firstBitUint64_t >> 1)) == 0);
			uint64_t prev = triplet.first.first + (firstBitUint64_t >> 1);
			uint64_t after = triplet.first.second;
			if (parent.count(prev) == 0) parent[prev] = prev;
			if (parent.count(after) == 0) parent[after] = after;
			merge(parent, prev, after);
		}
		if (!valid) continue;
		assert(coveredPrev.size() <= graph.edges.getEdges(std::make_pair(i, false)).size());
		assert(coveredAfter.size() <= graph.edges.getEdges(std::make_pair(i, true)).size());
		if (coveredPrev.size() < graph.edges.getEdges(std::make_pair(i, false)).size()) continue;
		if (coveredAfter.size() < graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		phmap::flat_hash_map<uint64_t, size_t> clusterToAllele;
		size_t alleleNum = 0;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t after = triplet.first.second;
			uint64_t cluster = find(parent, after);
			if (clusterToAllele.count(cluster) == 0)
			{
				clusterToAllele[cluster] = alleleNum;
				alleleNum += 1;
			}
		}
		assert(clusterToAllele.size() == alleleNum);
		assert(clusterToAllele.size() >= 1);
		if (clusterToAllele.size() == 1) continue;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t prev = triplet.first.first;
			uint64_t after = triplet.first.second;
			prevToAllele[i][prev] = clusterToAllele.at(find(parent, after));
			afterToAllele[i][after] = clusterToAllele.at(find(parent, after));
		}
		canResolveUnitig[i] = true;
	}
	phmap::flat_hash_set<size_t> readsAllowedToDiscard;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				if (!canResolveUnitig[mid & maskUint64_t]) continue;
				uint64_t prev = std::numeric_limits<size_t>::max();
				uint64_t after = std::numeric_limits<size_t>::max();
				if (k > 0)
				{
					prev = readPaths[i][j].path[k-1];
				}
				else if (j > 0)
				{
					prev = readPaths[i][j-1].path.back();
				}
				if (k+1 < readPaths[i][j].path.size())
				{
					after = readPaths[i][j].path[k+1];
				}
				else if (j+1 < readPaths[i].size())
				{
					after = readPaths[i][j+1].path[0];
				}
				if ((mid ^ firstBitUint64_t) & firstBitUint64_t)
				{
					std::swap(prev, after);
					if (prev != std::numeric_limits<size_t>::max()) prev ^= firstBitUint64_t;
					if (after != std::numeric_limits<size_t>::max()) after ^= firstBitUint64_t;
				}
				if (prev == std::numeric_limits<size_t>::max() && after == std::numeric_limits<size_t>::max())
				{
					if (readPaths[i].size() == 1 && readPaths[i][0].path.size() == 1)
					{
						readsAllowedToDiscard.insert(i);
						continue;
					}
					else
					{
						canResolveUnitig[mid & maskUint64_t] = false;
						continue;
					}
				}
				size_t prevAllele = std::numeric_limits<size_t>::max();
				size_t afterAllele = std::numeric_limits<size_t>::max();
				if (prev != std::numeric_limits<size_t>::max() && prevToAllele[mid & maskUint64_t].count(prev) == 1)
				{
					prevAllele = prevToAllele[mid & maskUint64_t].at(prev);
				}
				if (after != std::numeric_limits<size_t>::max() && afterToAllele[mid & maskUint64_t].count(after) == 1)
				{
					afterAllele = afterToAllele[mid & maskUint64_t].at(after);
				}
				if (prevAllele == std::numeric_limits<size_t>::max() && afterAllele == std::numeric_limits<size_t>::max())
				{
					canResolveUnitig[mid & maskUint64_t] = false;
					continue;
				}
				if (prevAllele != std::numeric_limits<size_t>::max() && afterAllele != std::numeric_limits<size_t>::max() && prevAllele != afterAllele)
				{
					canResolveUnitig[mid & maskUint64_t] = false;
					continue;
				}
			}
		}
	}
	assert(readPaths.size() == chunksPerRead.size());
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			assert(chunkBelongsToUnitig[chunk] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, size_t>> chunkAlleleReplacement;
	size_t nextNum = 0;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> unitigAllelesInThisRead;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (!canResolveUnitig[readPaths[i][j].path[k] & maskUint64_t]) continue;
				uint64_t prev = std::numeric_limits<size_t>::max();
				uint64_t after = std::numeric_limits<size_t>::max();
				if (k > 0)
				{
					prev = readPaths[i][j].path[k-1];
				}
				else if (j > 0)
				{
					prev = readPaths[i][j-1].path.back();
				}
				if (k+1 < readPaths[i][j].path.size())
				{
					after = readPaths[i][j].path[k+1];
				}
				else if (j+1 < readPaths[i].size())
				{
					after = readPaths[i][j+1].path[0];
				}
				if ((readPaths[i][j].path[k] ^ firstBitUint64_t) & firstBitUint64_t)
				{
					std::swap(prev, after);
					if (prev != std::numeric_limits<size_t>::max()) prev ^= firstBitUint64_t;
					if (after != std::numeric_limits<size_t>::max()) after ^= firstBitUint64_t;
				}
				if (prev == std::numeric_limits<size_t>::max() && after == std::numeric_limits<size_t>::max())
				{
					assert(readsAllowedToDiscard.count(i) == 1);
					assert(readPaths[i].size() == 1);
					assert(readPaths[i][0].path.size() == 1);
					continue;
				}
				assert(prev != std::numeric_limits<size_t>::max() || after != std::numeric_limits<size_t>::max());
				size_t alleleBefore = std::numeric_limits<size_t>::max();
				size_t alleleAfter = std::numeric_limits<size_t>::max();
				if (prev != std::numeric_limits<size_t>::max() && prevToAllele[readPaths[i][j].path[k] & maskUint64_t].count(prev) == 1) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(prev);
				if (after != std::numeric_limits<size_t>::max() && afterToAllele[readPaths[i][j].path[k] & maskUint64_t].count(after) == 1) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(after);
				assert(alleleBefore != std::numeric_limits<size_t>::max() || alleleAfter != std::numeric_limits<size_t>::max());
				assert(alleleBefore == alleleAfter || alleleBefore == std::numeric_limits<size_t>::max() || alleleAfter == std::numeric_limits<size_t>::max());
				size_t allele = alleleBefore;
				if (alleleBefore == std::numeric_limits<size_t>::max()) allele = alleleAfter;
				unitigAllelesInThisRead.emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k] & maskUint64_t, allele);
			}
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			uint64_t node = std::get<2>(chunksPerRead[i][j]);
			if (NonexistantChunk(node)) continue;
			size_t allele = std::numeric_limits<size_t>::max();
			if ((node & maskUint64_t) < chunkBelongsToUnitig.size() && chunkBelongsToUnitig[node & maskUint64_t] != std::numeric_limits<size_t>::max())
			{
				size_t unitig = chunkBelongsToUnitig[node & maskUint64_t];
				if (canResolveUnitig[unitig])
				{
					if (readsAllowedToDiscard.count(i) == 1)
					{
						assert(readPaths[i].size() == 1);
						assert(readPaths[i][0].path.size() == 1);
						std::get<2>(chunksPerRead[i][j]) = std::numeric_limits<size_t>::max();
						continue;
					}
					for (auto t : unitigAllelesInThisRead)
					{
						if (std::get<2>(t) != unitig) continue;
						if (std::get<0>(t) > std::get<0>(chunksPerRead[i][j])) continue;
						if (std::get<1>(t) < std::get<1>(chunksPerRead[i][j])) continue;
						assert(allele == std::numeric_limits<size_t>::max());
						allele = std::get<3>(t);
					}
					assert(allele != std::numeric_limits<size_t>::max());
				}
			}
			if (chunkAlleleReplacement.count(node & maskUint64_t) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t];
			}
			if (chunkAlleleReplacement.at(node & maskUint64_t).count(allele) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t][allele] = nextNum;
				nextNum += 1;
			}
			uint64_t replacement = chunkAlleleReplacement.at(node & maskUint64_t).at(allele);
			assert((replacement & firstBitUint64_t) == 0);
			replacement += (node & firstBitUint64_t);
			std::get<2>(chunksPerRead[i][j]) = replacement;
		}
	}
	size_t countUnitigsReplaced = 0;
	size_t countChunksReplaced = 0;
	for (size_t i = 0; i < canResolveUnitig.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		countUnitigsReplaced += 1;
		countChunksReplaced += graph.chunksInUnitig[i].size();
	}
	std::cerr << "semiambiguously resolvable unitig resolution resolved " << countUnitigsReplaced << " unitigs, " << countChunksReplaced << " chunks" << std::endl;
}

void resolveBetweenLongUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t longNodeThreshold, const size_t maxDistance, const size_t kmerSize)
{
	const uint64_t noNodeWithinDistance = (std::numeric_limits<size_t>::max() ^ (firstBitUint64_t>>1)) - 1;
	std::cerr << "resolving between long unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	phmap::flat_hash_set<size_t> chunkInsideLongUnitig;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < longNodeThreshold) continue;
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			chunkInsideLongUnitig.emplace(node & maskUint64_t);
		}
	}
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> longNodesInReads;
	longNodesInReads.resize(readPaths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (graph.unitigLengths[readPaths[i][j].path[k] & maskUint64_t] < longNodeThreshold) continue;
				longNodesInReads[i].emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k]);
			}
		}
		for (size_t j = 1; j < longNodesInReads[i].size(); j++)
		{
			assert(std::get<0>(longNodesInReads[i][j]) > std::get<0>(longNodesInReads[i][j-1]));
			assert(std::get<1>(longNodesInReads[i][j]) > std::get<1>(longNodesInReads[i][j-1]));
		}
	}
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	size_t countSplitted = 0;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&chunksPerRead, &occurrencesPerChunk, &longNodesInReads, &rawReadLengths, &nextNum, &resultMutex, &chunkInsideLongUnitig, &countSplitted, maxDistance, noNodeWithinDistance](const size_t chunki)
	{
		if (chunkInsideLongUnitig.count(chunki) == 1 || occurrencesPerChunk[chunki].size() < 4)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		std::vector<std::tuple<uint64_t, uint64_t, size_t, size_t>> predecessorsAndSuccessors;
		for (auto t : occurrencesPerChunk[chunki])
		{
			size_t read = t.first;
			size_t offset = t.second;
			bool fw = std::get<2>(chunksPerRead[read][offset]) & firstBitUint64_t;
			size_t chunkstartpos = std::get<0>(chunksPerRead[read][offset]);
			size_t chunkendpos = std::get<1>(chunksPerRead[read][offset]);
			uint64_t predecessor = std::numeric_limits<size_t>::max();
			uint64_t successor = std::numeric_limits<size_t>::max();
			for (size_t j = 0; j < longNodesInReads[read].size(); j++)
			{
				if (std::get<1>(longNodesInReads[read][j]) < chunkendpos && std::get<1>(longNodesInReads[read][j])+maxDistance > chunkstartpos)
				{
					predecessor = std::get<2>(longNodesInReads[read][j]);
				}
				if (std::get<0>(longNodesInReads[read][j]) > chunkstartpos && std::get<0>(longNodesInReads[read][j]) < chunkendpos+maxDistance && successor == std::numeric_limits<size_t>::max())
				{
					successor = std::get<2>(longNodesInReads[read][j]);
				}
			}
			if (chunkstartpos > maxDistance && predecessor == std::numeric_limits<size_t>::max()) predecessor = noNodeWithinDistance;
			if (chunkendpos + maxDistance < rawReadLengths[read] && successor == std::numeric_limits<size_t>::max()) successor = noNodeWithinDistance;
			if (!fw)
			{
				std::swap(predecessor, successor);
				if (predecessor != std::numeric_limits<size_t>::max() && predecessor != noNodeWithinDistance) predecessor ^= firstBitUint64_t;
				if (successor != std::numeric_limits<size_t>::max() && successor != noNodeWithinDistance) successor ^= firstBitUint64_t;
			}
			predecessorsAndSuccessors.emplace_back(predecessor, successor, read, offset);
		}
		assert(predecessorsAndSuccessors.size() == occurrencesPerChunk[chunki].size());
		bool hasAny = false;
		phmap::flat_hash_set<uint64_t> mergedSomewhere;
		phmap::flat_hash_map<uint64_t, uint64_t> parent;
		for (auto t : predecessorsAndSuccessors)
		{
			assert(std::get<1>(t) == std::numeric_limits<size_t>::max() || ((std::get<1>(t) & (firstBitUint64_t>>1)) == 0));
			if (std::get<0>(t) != std::numeric_limits<size_t>::max() && parent.count(std::get<0>(t)) == 0)
			{
				hasAny = true;
				parent[std::get<0>(t)] = std::get<0>(t);
			}
			if (std::get<1>(t) != std::numeric_limits<size_t>::max() && parent.count(std::get<1>(t) + (firstBitUint64_t>>1)) == 0)
			{
				hasAny = true;
				parent[std::get<1>(t) + (firstBitUint64_t>>1)] = std::get<1>(t) + (firstBitUint64_t>>1);
			}
			if (std::get<0>(t) != std::numeric_limits<size_t>::max() && std::get<1>(t) != std::numeric_limits<size_t>::max())
			{
				mergedSomewhere.insert(std::get<0>(t));
				mergedSomewhere.insert(std::get<1>(t) + (firstBitUint64_t>>1));
				merge(parent, std::get<0>(t), std::get<1>(t) + (firstBitUint64_t>>1));
			}
		}
		if (!hasAny)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "chunk " << chunki << " doesn't have any" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		bool canSplit = true;
		for (auto pair : parent)
		{
			if (mergedSomewhere.count(pair.first) == 0)
			{
				canSplit = false;
				break;
			}
		}
		if (!canSplit)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "can't split chunk " << chunki << " because of missing merge" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		phmap::flat_hash_map<uint64_t, size_t> clusterToNode;
		size_t nextId = 0;
		for (uint64_t node : mergedSomewhere)
		{
			uint64_t found = find(parent, node);
			if (clusterToNode.count(found) == 1)
			{
				continue;
			}
			clusterToNode[found] = nextId;
			nextId += 1;
		}
		phmap::flat_hash_map<uint64_t, size_t> clusterCoverage;
		for (auto t : predecessorsAndSuccessors)
		{
			if (std::get<0>(t) != std::numeric_limits<size_t>::max())
			{
				clusterCoverage[clusterToNode.at(find(parent, std::get<0>(t)))] += 1;
			}
			else if (std::get<1>(t) != std::numeric_limits<size_t>::max())
			{
				clusterCoverage[clusterToNode.at(find(parent, std::get<1>(t) + (firstBitUint64_t>>1)))] += 1;
			}
		}
		for (auto pair : clusterCoverage)
		{
			if (pair.second < 2)
			{
				canSplit = false;
			}
		}
		if (!canSplit)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "can't split chunk " << chunki << " because of low coverage" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		assert(clusterToNode.size() >= 1);
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			if (clusterToNode.size() >= 2) countSplitted += 1;
//			std::cerr << "chunk " << chunki << " split to " << clusterToNode.size() << " chunks" << std::endl;
			for (auto t : predecessorsAndSuccessors)
			{
				if (std::get<0>(t) == std::numeric_limits<size_t>::max() && std::get<1>(t) == std::numeric_limits<size_t>::max())
				{
					std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) = std::numeric_limits<size_t>::max();
					continue;
				}
				size_t cluster;
				if (std::get<0>(t) != std::numeric_limits<size_t>::max())
				{
					cluster = clusterToNode.at(find(parent, std::get<0>(t)));
				}
				else
				{
					assert(std::get<1>(t) != std::numeric_limits<size_t>::max());
					cluster = clusterToNode.at(find(parent, std::get<1>(t) + (firstBitUint64_t>>1)));
				}
				std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) = nextNum + cluster + (std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) & firstBitUint64_t);
			}
			nextNum += clusterToNode.size();
			return;
		}
	});
	std::cerr << "between long unitigs resolved " << countSplitted << " chunks" << std::endl;
}

void resolveUnambiguouslyResolvableUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving unambiguously resolvable unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<bool> edgesBalanced;
	std::vector<bool> unitigHasContainedPath;
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> tripletsPerUnitig;
	edgesBalanced.resize(graph.unitigLengths.size(), false);
	unitigHasContainedPath.resize(graph.unitigLengths.size(), false);
	tripletsPerUnitig.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.edges.getEdges(std::make_pair(i, true)).size() == graph.edges.getEdges(std::make_pair(i, false)).size())
		{
			if (graph.edges.getEdges(std::make_pair(i, true)).size() >= 2)
			{
				edgesBalanced[i] = true;
			}
		}
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			if (readPaths[i][j].path.size() == 1)
			{
				unitigHasContainedPath[readPaths[i][j].path[0] & maskUint64_t] = true;
				continue;
			}
			assert(readPaths[i][j].path.size() >= 2);
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				if (!edgesBalanced[mid & maskUint64_t]) continue;
				if (unitigHasContainedPath[mid & maskUint64_t]) continue;
				uint64_t prev, after;
				if (mid & firstBitUint64_t)
				{
					prev = readPaths[i][j].path[k-1];
					after = readPaths[i][j].path[k+1];
				}
				else
				{
					prev = readPaths[i][j].path[k+1] ^ firstBitUint64_t;
					after = readPaths[i][j].path[k-1] ^ firstBitUint64_t;
				}
				tripletsPerUnitig[mid & maskUint64_t][std::make_pair(prev, after)] += 1;
			}
		}
	}
	std::vector<bool> canResolveUnitig;
	canResolveUnitig.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (!edgesBalanced[i]) continue;
		if (unitigHasContainedPath[i]) continue;
		if (tripletsPerUnitig[i].size() != graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		phmap::flat_hash_set<uint64_t> nodesPrev;
		phmap::flat_hash_set<uint64_t> nodesAfter;
		bool valid = true;
		for (auto triplet : tripletsPerUnitig[i])
		{
			if (triplet.second < 2)
			{
				valid = false;
				break;
			}
			nodesPrev.insert(triplet.first.first);
			nodesAfter.insert(triplet.first.second);
		}
		if (!valid) continue;
		if (nodesPrev.size() != nodesAfter.size()) continue;
		if (nodesPrev.size() != graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		canResolveUnitig[i] = true;
	}
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> prevToAllele;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> afterToAllele;
	prevToAllele.resize(graph.unitigLengths.size());
	afterToAllele.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		size_t alleleNum = 0;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t prev = triplet.first.first;
			uint64_t after = triplet.first.second;
			assert(prevToAllele[i].count(prev) == 0);
			assert(afterToAllele[i].count(after) == 0);
			prevToAllele[i][prev] = alleleNum;
			afterToAllele[i][after] = alleleNum;
			alleleNum += 1;
		}
	}
	assert(readPaths.size() == chunksPerRead.size());
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			assert(chunkBelongsToUnitig[chunk] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, size_t>> chunkAlleleReplacement;
	size_t nextNum = 0;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> unitigAllelesInThisRead;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (!canResolveUnitig[readPaths[i][j].path[k] & maskUint64_t]) continue;
				size_t alleleBefore = std::numeric_limits<size_t>::max();
				size_t alleleAfter = std::numeric_limits<size_t>::max();
				if (readPaths[i][j].path[k] & firstBitUint64_t)
				{
					if (k > 0) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k-1]);
					if (k+1 < readPaths[i][j].path.size()) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k+1]);
				}
				else
				{
					if (k > 0) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k-1] ^ firstBitUint64_t);
					if (k+1 < readPaths[i][j].path.size()) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k+1] ^ firstBitUint64_t);
				}
				assert(alleleBefore != std::numeric_limits<size_t>::max() || alleleAfter != std::numeric_limits<size_t>::max());
				assert(alleleBefore == alleleAfter || alleleBefore == std::numeric_limits<size_t>::max() || alleleAfter == std::numeric_limits<size_t>::max());
				size_t allele = alleleBefore;
				if (alleleBefore == std::numeric_limits<size_t>::max()) allele = alleleAfter;
				unitigAllelesInThisRead.emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k] & maskUint64_t, allele);
			}
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			uint64_t node = std::get<2>(chunksPerRead[i][j]);
			if (NonexistantChunk(node)) continue;
			size_t allele = std::numeric_limits<size_t>::max();
			if ((node & maskUint64_t) < chunkBelongsToUnitig.size() && chunkBelongsToUnitig[node & maskUint64_t] != std::numeric_limits<size_t>::max())
			{
				size_t unitig = chunkBelongsToUnitig[node & maskUint64_t];
				if (canResolveUnitig[unitig])
				{
					for (auto t : unitigAllelesInThisRead)
					{
						if (std::get<2>(t) != unitig) continue;
						if (std::get<0>(t) > std::get<0>(chunksPerRead[i][j])) continue;
						if (std::get<1>(t) < std::get<1>(chunksPerRead[i][j])) continue;
						assert(allele == std::numeric_limits<size_t>::max());
						allele = std::get<3>(t);
					}
					assert(allele != std::numeric_limits<size_t>::max());
				}
			}
			if (chunkAlleleReplacement.count(node & maskUint64_t) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t];
			}
			if (chunkAlleleReplacement.at(node & maskUint64_t).count(allele) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t][allele] = nextNum;
				nextNum += 1;
			}
			uint64_t replacement = chunkAlleleReplacement.at(node & maskUint64_t).at(allele);
			assert((replacement & firstBitUint64_t) == 0);
			replacement += (node & firstBitUint64_t);
			std::get<2>(chunksPerRead[i][j]) = replacement;
		}
	}
	size_t countUnitigsReplaced = 0;
	size_t countChunksReplaced = 0;
	for (size_t i = 0; i < canResolveUnitig.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		countUnitigsReplaced += 1;
		countChunksReplaced += graph.chunksInUnitig[i].size();
	}
	std::cerr << "unambiguously resolvable unitig resolution resolved " << countUnitigsReplaced << " unitigs, " << countChunksReplaced << " chunks" << std::endl;
}

void countReadRepetitiveUnitigs(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "counting read-repetitive unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	std::vector<uint8_t> chunkRepeatCount;
	chunkRepeatCount.resize(chunkBelongsToUnitig.size(), 0);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		phmap::flat_hash_set<size_t> foundHere;
		phmap::flat_hash_set<size_t> repetitiveHere;
		for (auto t : chunksPerRead[i])
		{
			size_t chunk = std::get<2>(t) & maskUint64_t;
			if (foundHere.count(chunk) == 1) repetitiveHere.insert(chunk);
			foundHere.insert(chunk);
		}
		for (size_t chunk : repetitiveHere)
		{
			if (chunk >= chunkRepeatCount.size()) continue;
			if (chunkRepeatCount[chunk] >= 2) continue;
			chunkRepeatCount[chunk] += 1;
		}
	}
	std::vector<bool> unitigRepetitive;
	unitigRepetitive.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < chunkRepeatCount.size(); i++)
	{
		if (chunkRepeatCount[i] < 2) continue;
		if (chunkBelongsToUnitig[i] == std::numeric_limits<size_t>::max()) continue;
		unitigRepetitive[chunkBelongsToUnitig[i]] = true;
	}
	size_t countRepetitive = 0;
	std::ofstream file { "repetitive.txt" };
	for (size_t i = 0; i < unitigRepetitive.size(); i++)
	{
		if (!unitigRepetitive[i]) continue;
		countRepetitive += 1;
		file << i << std::endl;
	}
	std::cerr << countRepetitive << " read-repetitive unitigs" << std::endl;
}

std::vector<std::pair<size_t, size_t>> getNoncontainedUniqueRanges(std::vector<std::pair<size_t, size_t>> raw)
{
	assert(raw.size() >= 1);
	std::sort(raw.begin(), raw.end());
	for (size_t j = raw.size()-1; j < raw.size() && j > 0; j--)
	{
		if (raw[j] != raw[j-1]) continue;
		std::swap(raw[j], raw.back());
		raw.pop_back();
	}
	assert(raw.size() >= 1);
	std::sort(raw.begin(), raw.end());
	size_t longestExtent = raw[0].second;
	std::vector<bool> contained;
	contained.resize(raw.size(), false);
	for (size_t j = 1; j < raw.size(); j++)
	{
		assert(raw[j].first >= raw[j-1].first);
		if (raw[j].second <= longestExtent) contained[j] = true;
		if (raw[j].first == raw[j-1].first)
		{
			assert(raw[j].second > raw[j-1].second);
			contained[j-1] = true;
		}
		longestExtent = std::max(longestExtent, raw[j].second);
	}
	for (size_t j = raw.size()-1; j < raw.size(); j--)
	{
		if (!contained[j]) continue;
		std::swap(raw[j], raw.back());
		raw.pop_back();
	}
	std::sort(raw.begin(), raw.end());
	return raw;
}

void expandChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t expandedSize)
{
	std::vector<size_t> chunkLengths;
	{
		std::vector<std::vector<size_t>> rawChunkLengths;
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
				size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
				size_t length = std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) + 1;
				while (chunk >= rawChunkLengths.size()) rawChunkLengths.emplace_back();
				rawChunkLengths[chunk].emplace_back(length);
			}
		}
		chunkLengths.resize(rawChunkLengths.size());
		for (size_t i = 0; i < rawChunkLengths.size(); i++)
		{
			assert(rawChunkLengths[i].size() >= 1);
			std::sort(rawChunkLengths[i].begin(), rawChunkLengths[i].end());
			chunkLengths[i] = rawChunkLengths[i][rawChunkLengths[i].size()/2];
		}
	}
	std::map<std::vector<uint64_t>, size_t> chunkmerToNewChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<std::pair<size_t, size_t>> newChunkmerRanges;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) newChunkmerRanges.emplace_back(j, j);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t end = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			while (end < chunksPerRead[i].size() && length < expandedSize)
			{
				end += 1;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end]))) break;
				if (end == chunksPerRead[i].size()) break;
				length += chunkLengths[std::get<2>(chunksPerRead[i][end]) & maskUint64_t];
				length -= kmerSize;
			}
			if (length < expandedSize) continue;
			if (end == chunksPerRead[i].size()) continue;
			newChunkmerRanges.emplace_back(j, end);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t start = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			while (start > 0 && length < expandedSize)
			{
				start -= 1;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][start]))) break;
				length += chunkLengths[std::get<2>(chunksPerRead[i][start]) & maskUint64_t];
				length -= kmerSize;
			}
			if (length < expandedSize) continue;
			newChunkmerRanges.emplace_back(start, j);
		}
		if (newChunkmerRanges.size() == 0)
		{
			chunksPerRead[i].clear();
			continue;
		}
		newChunkmerRanges = getNoncontainedUniqueRanges(newChunkmerRanges);
		assert(newChunkmerRanges.size() >= 1);
		for (size_t j = 1; j < newChunkmerRanges.size(); j++)
		{
			assert(newChunkmerRanges[j].first > newChunkmerRanges[j-1].first);
			assert(newChunkmerRanges[j].second > newChunkmerRanges[j-1].second);
		}
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunks;
		for (size_t j = 0; j < newChunkmerRanges.size(); j++)
		{
			size_t newchunk = std::numeric_limits<size_t>::max();
			std::vector<uint64_t> chunkmer;
			bool valid = true;
			for (size_t k = newChunkmerRanges[j].first; k <= newChunkmerRanges[j].second; k++)
			{
				chunkmer.emplace_back(std::get<2>(chunksPerRead[i][k]));
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][k]))) valid = false;
			}
			if (valid)
			{
				if (chunkmerToNewChunk.count(chunkmer) == 1)
				{
					newchunk = chunkmerToNewChunk.at(chunkmer) + firstBitUint64_t;
				}
				else
				{
					newchunk = chunkmerToNewChunk.size();
					chunkmerToNewChunk[chunkmer] = newchunk;
					newchunk += firstBitUint64_t;
				}
			}
			else
			{
				newchunk = std::numeric_limits<size_t>::max();
			}
			newChunks.emplace_back(std::get<0>(chunksPerRead[i][newChunkmerRanges[j].first]), std::get<1>(chunksPerRead[i][newChunkmerRanges[j].second]), newchunk);
		}
		std::swap(chunksPerRead[i], newChunks);
	}
}

void resolveTinyNodesRecklessly(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize)
{
	size_t nextNum = 0;
	std::mutex resultMutex;
	size_t countResolved = 0;
	size_t countResolvedTo = 0;
	std::vector<std::tuple<size_t, size_t, size_t>> replacements;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &countResolved, &countResolvedTo, &replacements, kmerSize](const size_t chunki, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		if (occurrencesPerChunk[chunki].size() <= 2)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				replacements.emplace_back(occurrencesPerChunk[chunki][j].first, occurrencesPerChunk[chunki][j].second, nextNum);
			}
			nextNum += 1;
			return;
		}
		bool isSmall = true;
		for (auto pair : occurrencesPerChunk[chunki])
		{
			auto t = chunksPerRead[pair.first][pair.second];
			if (std::get<1>(t) - std::get<0>(t) >= kmerSize*2)
			{
				isSmall = false;
				break;
			}
		}
		if (!isSmall)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				replacements.emplace_back(occurrencesPerChunk[chunki][j].first, occurrencesPerChunk[chunki][j].second, nextNum);
			}
			nextNum += 1;
			return;
		}
		phmap::flat_hash_map<size_t, size_t> parent;
		for (auto pair : occurrencesPerChunk[chunki])
		{
			size_t prev = std::numeric_limits<size_t>::max();
			size_t next = std::numeric_limits<size_t>::max();
			if (pair.second > 0 && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second-1])))
			{
				prev = std::get<2>(chunksPerRead[pair.first][pair.second-1]);
				if (parent.count(prev) == 0) parent[prev] = prev;
			}
			if (pair.second+1 < chunksPerRead[pair.first].size() && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second+1])))
			{
				next = std::get<2>(chunksPerRead[pair.first][pair.second+1]) + (firstBitUint64_t >> 1);
				if (parent.count(next) == 0) parent[next] = next;
			}
			if (prev != std::numeric_limits<size_t>::max() && next != std::numeric_limits<size_t>::max())
			{
				merge(parent, prev, next);
			}
		}
		phmap::flat_hash_map<size_t, size_t> clusterToNode;
		size_t nextClusterNum = 1;
		for (auto pair : parent)
		{
			if (pair.second != pair.first) continue;
			clusterToNode[pair.second] = nextClusterNum;
			nextClusterNum += 1;
		}
		size_t zeroOffset = 0;
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			zeroOffset = nextNum;
			nextNum += nextClusterNum;
			if (nextClusterNum > 2)
			{
				countResolved += 1;
				countResolvedTo += nextClusterNum-1;
			}
			for (auto pair : occurrencesPerChunk[chunki])
			{
				size_t clusterNum = 0;
				if (pair.second > 0 && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second-1])))
				{
					clusterNum = clusterToNode.at(find(parent, std::get<2>(chunksPerRead[pair.first][pair.second-1])));
				}
				if (pair.second+1 < chunksPerRead[pair.first].size() && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second+1])))
				{
					size_t foundClusterNum = clusterToNode.at(find(parent, std::get<2>(chunksPerRead[pair.first][pair.second+1]) + (firstBitUint64_t >> 1)));
					assert(clusterNum == 0 || foundClusterNum == clusterNum);
					clusterNum = foundClusterNum;
				}
				replacements.emplace_back(pair.first, pair.second, zeroOffset + clusterNum);
			}
		}
	});
	for (auto t : replacements)
	{
		std::get<2>(chunksPerRead[std::get<0>(t)][std::get<1>(t)]) = std::get<2>(t) + (std::get<2>(chunksPerRead[std::get<0>(t)][std::get<1>(t)]) & firstBitUint64_t);
	}
	std::cerr << "tiny node resolution resolved " << countResolved << " chunks to " << countResolvedTo << " chunks" << std::endl;
}

void removeTinyProblemNodes(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			assert(std::get<2>(t) & firstBitUint64_t);
			maxChunk = std::max(maxChunk, std::get<2>(t) & maskUint64_t);
		}
	}
	std::vector<size_t> uniquePredecessor;
	std::vector<size_t> uniqueSuccessor;
	uniquePredecessor.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	uniqueSuccessor.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1]))) continue;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t prev = std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t;
			size_t curr = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			assert(prev < uniqueSuccessor.size());
			assert(curr < uniquePredecessor.size());
			if (uniqueSuccessor[prev] == std::numeric_limits<size_t>::max())
			{
				uniqueSuccessor[prev] = curr;
			}
			else if (uniqueSuccessor[prev] != curr)
			{
				uniqueSuccessor[prev] = std::numeric_limits<size_t>::max()-1;
			}
			if (uniquePredecessor[curr] == std::numeric_limits<size_t>::max())
			{
				uniquePredecessor[curr] = prev;
			}
			else if (uniquePredecessor[curr] != prev)
			{
				uniquePredecessor[curr] = std::numeric_limits<size_t>::max()-1;
			}
		}
	}
	std::vector<bool> canBeRemoved;
	canBeRemoved.resize(maxChunk+1, true);
	for (size_t i = 0; i < uniquePredecessor.size(); i++)
	{
		if (uniquePredecessor[i] != std::numeric_limits<size_t>::max()-1 && uniqueSuccessor[i] != std::numeric_limits<size_t>::max()-1)
		{
			canBeRemoved[i] = false;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i][0])) && std::get<1>(chunksPerRead[i][0]) - std::get<0>(chunksPerRead[i][0]) >= 2*kmerSize) canBeRemoved[std::get<2>(chunksPerRead[i][0]) & maskUint64_t] = false;
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i].back())) && std::get<1>(chunksPerRead[i].back()) - std::get<0>(chunksPerRead[i].back()) >= 2*kmerSize) canBeRemoved[std::get<2>(chunksPerRead[i].back()) & maskUint64_t] = false;
		for (size_t j = 1; j+1 < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			if (std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) >= 2*kmerSize)
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j])] = false;
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j+1])))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j])] = false;
				continue;
			}
			if (canBeRemoved[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (canBeRemoved[std::get<2>(chunksPerRead[i][j+1]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (std::get<1>(chunksPerRead[i][j-1])+1 <= std::get<0>(chunksPerRead[i][j+1]))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
		}
	}
	for (size_t i = 0; i < chunksPerRead.size()/2; i++)
	{
		size_t other = i + chunksPerRead.size()/2;
		assert(chunksPerRead[i].size() == chunksPerRead[other].size());
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) - std::get<0>(chunksPerRead[other][chunksPerRead[other].size()-1-j]));
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				assert(NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])));
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])))
			{
				std::cerr << i << " " << other << std::endl;
				std::cerr << j << " " << chunksPerRead[i].size() << std::endl;
				std::cerr << std::get<2>(chunksPerRead[i][j]) << " " << ((std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(chunksPerRead[i][j]) & maskUint64_t) << std::endl;
				std::cerr << std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) << " " << ((std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t) << std::endl;
			}
			assert(!NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])));
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t] = false;
			}
			if (!canBeRemoved[std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
			}
		}
	}
	size_t countRemovableChunks = 0;
	for (size_t i = 0; i < canBeRemoved.size(); i++)
	{
		if (!canBeRemoved[i]) continue;
		countRemovableChunks += 1;
	}
	size_t countRemoved = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = chunksPerRead[i].size()-1; j < chunksPerRead[i].size(); j--)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			std::swap(chunksPerRead[i][j], chunksPerRead[i].back());
			chunksPerRead[i].pop_back();
			countRemoved += 1;
		}
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<0>(chunksPerRead[i][j]) > std::get<0>(chunksPerRead[i][j-1]));
			assert(std::get<1>(chunksPerRead[i][j]) > std::get<1>(chunksPerRead[i][j-1]));
		}
	}
	std::cerr << countRemovableChunks << " distinct removable chunks, removed " << countRemoved << " total chunks" << std::endl;
}

void makeGraph(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize, const size_t windowSize, const size_t middleSkip, const size_t startStage)
{
	std::cerr << "start at stage " << startStage << std::endl;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead;
	if (startStage > 0)
	{
		chunksPerRead = readChunksFromFakePathsFile("fakepaths" + std::to_string(startStage) + ".txt");
		while (chunksPerRead.size() < sequenceIndex.size()*2) chunksPerRead.emplace_back();
		std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
	switch(startStage)
	{
		case 0:
			std::cerr << "getting chunks from reads" << std::endl;
			chunksPerRead = getMinimizerBoundedChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize);
			//addMissingPiecesBetweenChunks(chunksPerRead, k);
			{
				size_t numChunks = 0;
				for (size_t i = 0; i < chunksPerRead.size(); i++)
				{
					numChunks += chunksPerRead[i].size();
				}
				std::cerr << numChunks << " chunks" << std::endl;
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
			splitPerFirstLastKmers(sequenceIndex, chunksPerRead, kmerSize);
			splitPerLength(chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(1, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 1:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerBaseCounts(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(2, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 2:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerMinHashes(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(3, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 3:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(4, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 4:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentityRoughly(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentity(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(5, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 5:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			expandChunks(chunksPerRead, kmerSize, 5000);
			writeStage(51, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 51:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(6, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 6:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(7, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 7:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerNearestNeighborPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(8, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 8:
			splitPerCorrectedKmerPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(9, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 9:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(10, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 10:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(11, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 11:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(12, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 12:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(13, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 13:
//			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(14, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 14:
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, 2, kmerSize);
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.25, kmerSize);
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.5, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(15, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 15:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(16, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 16:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(17, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 17:
//			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			writeStage(18, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 18:
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			countReadRepetitiveUnitigs(chunksPerRead, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadChunkSequences("sequences-chunk19.txt", rawReadLengths, chunksPerRead, sequenceIndex);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(19, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			writeBidirectedUnitigGraph("graph-dbg-final.gfa", "paths-dbg-final.gaf", "unitigs-dbg-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadUnitigSequences("sequences-dbg-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 19:
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 50000, 100000);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(20, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 20:
			writeBidirectedUnitigGraph("graph-final.gfa", "paths-final.gaf", "unitigs-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadUnitigSequences("sequences-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
}

int main(int argc, char** argv)
{
	const size_t numThreads = std::stoi(argv[1]);
	const size_t k = std::stoull(argv[2]);
	const size_t windowSize = std::stoull(argv[3]);
	const size_t middleSkip = std::stoull(argv[4]);
	const double approxOneHapCoverage = std::stod(argv[5]);
	mismatchFraction = std::stod(argv[6]); // try 2-3x average error rate
	const size_t startStage = std::stoull(argv[7]);
	std::vector<std::string> readFiles;
	for (size_t i = 8; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	FastaCompressor::CompressedStringIndex sequenceIndex { 5, 100 };
	std::vector<size_t> readBasepairLengths;
	for (auto file : readFiles)
	{
		std::cerr << "reading from file " << file << std::endl;
		readFilesAndAddToSequenceIndex(file, sequenceIndex, readBasepairLengths, numThreads);
	}
	sequenceIndex.removeConstructionVariables();
	size_t lastReal = readBasepairLengths.size();
	for (size_t i = 0; i < lastReal; i++)
	{
		readBasepairLengths.emplace_back(readBasepairLengths[i]);
	}
	std::cerr << sequenceIndex.size() << " reads" << std::endl;
	makeGraph(sequenceIndex, readBasepairLengths, numThreads, approxOneHapCoverage, k, windowSize, middleSkip, startStage);
}
