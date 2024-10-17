#include <map>
#include <cassert>
#include <iostream>
#include <fstream>
#include "ChunkPhasing.h"
#include "Common.h"
#include "SequenceHelper.h"
#include "UnionFind.h"
#include "KmerIterator.h"
#include "ChunkUnitigGraph.h"
#include "TransitiveClosure.h"
#include "CanonHelper.h"
#include "edlib.h"
#include "ConsensusMaker.h"

class TripletSplittingQueue
{
public:
	TripletSplittingQueue() = default;
	TripletSplittingQueue(const phmap::flat_hash_map<std::string, size_t>* filteredCoverages, const std::vector<size_t>* coveredIndices, const std::vector<std::vector<uint8_t>>* readFakeMSABases, std::vector<std::vector<uint8_t>>* result, std::vector<bool>* indexUsed, std::atomic<size_t>* countRunning) :
		queueMutex(),
		resultMutex(),
		filteredCoverages(filteredCoverages),
		coveredIndices(coveredIndices),
		readFakeMSABases(readFakeMSABases),
		result(result),
		indexUsed(indexUsed),
		countRunning(countRunning),
		nextChunki(0),
		nextChunkj(0),
		nextChunkk(0),
		maxChunkIndex((coveredIndices->size()+9)/10)
	{
	}
	struct Item
	{
	public:
		Item() :
			filteredCoverages(nullptr),
			coveredIndices(nullptr),
			readFakeMSABases(nullptr),
			result(nullptr),
			indexUsed(nullptr),
			countRunning(nullptr),
			resultMutex(nullptr),
			chunki(std::numeric_limits<size_t>::max()),
			chunkj(std::numeric_limits<size_t>::max()),
			chunkk(std::numeric_limits<size_t>::max())
		{
		}
		const phmap::flat_hash_map<std::string, size_t>* filteredCoverages;
		const std::vector<size_t>* coveredIndices;
		const std::vector<std::vector<uint8_t>>* readFakeMSABases;
		std::vector<std::vector<uint8_t>>* result;
		std::vector<bool>* indexUsed;
		std::atomic<size_t>* countRunning;
		std::mutex* resultMutex;
		size_t chunki;
		size_t chunkj;
		size_t chunkk;
	};
	bool empty() const
	{
		std::lock_guard<std::mutex> lock { queueMutex };
		assert(nextChunki <= maxChunkIndex);
		return nextChunki == maxChunkIndex;
	}
	bool popItem(Item* item)
	{
		std::lock_guard<std::mutex> lock { queueMutex };
		if (nextChunki == maxChunkIndex) return false;
		assert(nextChunki < maxChunkIndex);
		assert(nextChunkj < maxChunkIndex);
		assert(nextChunkk < maxChunkIndex);
		item->resultMutex = &resultMutex;
		item->filteredCoverages = filteredCoverages;
		item->coveredIndices = coveredIndices;
		item->readFakeMSABases = readFakeMSABases;
		item->result = result;
		item->indexUsed = indexUsed;
		item->countRunning = countRunning;
		item->chunki = nextChunki;
		item->chunkj = nextChunkj;
		item->chunkk = nextChunkk;
		(*countRunning) += 1;
		nextChunkk += 1;
		if (nextChunkk == maxChunkIndex)
		{
			nextChunkj += 1;
			if (nextChunkj == maxChunkIndex)
			{
				nextChunki += 1;
				nextChunkj = nextChunki;
			}
			nextChunkk = nextChunkj;
		}
		return true;
	}
private:
	mutable std::mutex queueMutex;
	std::mutex resultMutex;
	const phmap::flat_hash_map<std::string, size_t>* filteredCoverages;
	const std::vector<size_t>* coveredIndices;
	const std::vector<std::vector<uint8_t>>* readFakeMSABases;
	std::vector<std::vector<uint8_t>>* result;
	std::vector<bool>* indexUsed;
	std::atomic<size_t>* countRunning;
	size_t nextChunki;
	size_t nextChunkj;
	size_t nextChunkk;
	size_t maxChunkIndex;
};

class SNPAlignmentQueue
{
public:
	struct Item
	{
	public:
		Item() = default;
		Item(const std::vector<std::string>* strings, std::string* consensus, size_t minIndex, size_t maxIndex, std::vector<std::vector<uint8_t>>* result, std::atomic<size_t>* totalMismatches, std::atomic<size_t>* totalSequenceLength, std::atomic<bool>* anyAlignmentError, std::atomic<size_t>* countRunning) :
			strings(strings),
			consensus(consensus),
			minIndex(minIndex),
			maxIndex(maxIndex),
			result(result),
			totalMismatches(totalMismatches),
			totalSequenceLength(totalSequenceLength),
			anyAlignmentError(anyAlignmentError),
			countRunning(countRunning)
		{
		}
		const std::vector<std::string>* strings;
		std::string* consensus;
		size_t minIndex;
		size_t maxIndex;
		std::vector<std::vector<uint8_t>>* result;
		std::atomic<size_t>* totalMismatches;
		std::atomic<size_t>* totalSequenceLength;
		std::atomic<bool>* anyAlignmentError;
		std::atomic<size_t>* countRunning;
	};
	void insertItems(const std::vector<std::string>* strings, std::string* consensus, size_t minIndex, size_t maxIndex, size_t splitCount, std::vector<std::vector<uint8_t>>* result, std::atomic<size_t>* totalMismatches, std::atomic<size_t>* totalSequenceLength, std::atomic<bool>* anyAlignmentError, std::atomic<size_t>* countRunning)
	{
		std::lock_guard<std::mutex> lock { mutex };
		size_t lastEnd = minIndex;
		for (size_t i = 0; i < splitCount; i++)
		{
			size_t nextIndex = (maxIndex-minIndex)*(i+1)/splitCount+minIndex;
			assert(nextIndex > lastEnd);
			assert(nextIndex <= maxIndex);
			items.emplace_back(strings, consensus, lastEnd, nextIndex, result, totalMismatches, totalSequenceLength, anyAlignmentError, countRunning);
			lastEnd = nextIndex;
		}
	}
	bool popItem(Item* result)
	{
		std::lock_guard<std::mutex> lock { mutex };
		if (items.size() == 0) return false;
		std::swap(items.back(), *result);
		items.pop_back();
		(*result->countRunning) += 1;
		return true;
	}
	size_t size() const
	{
		return items.size();
	}
private:
	std::mutex mutex;
	std::vector<Item> items;
};

size_t getHammingdistance(const std::vector<uint8_t>& left, const std::vector<uint8_t>& right, const size_t maxDistance)
{
	assert(left.size() == right.size());
	size_t result = 0;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] != right[i])
		{
			result += 1;
			if (result > maxDistance) return result;
		}
	}
	return result;
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
	assert(vec.size() == keep.size());
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

bool allelesMatchApproximately(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right)
{
	assert(left.size() >= 1);
	assert(right.size() >= 1);
	if (left.size() < 2) return false;
	if (right.size() < 2) return false;
	if (left[0].size() < 10) return false;
	if (right[0].size() < 10) return false;
	if (left[1].size() < 10) return false;
	if (right[1].size() < 10) return false;
	if (left.size() > right.size() && left[right.size()].size() >= 10) return false;
	if (right.size() > left.size() && right[left.size()].size() >= 10) return false;
	for (size_t i = 0; i < left.size() && i < right.size(); i++)
	{
		if (left[i].size() < right[i].size() * 0.9) return false;
		if (right[i].size() < left[i].size() * 0.9) return false;
		if (left[i].size() < 10) break;
		if (right[i].size() < 10) break;
	}
	for (size_t i = 0; i < left.size() && i < right.size(); i++)
	{
		size_t intersect = intersectSize(left[i], right[i]);
		if (intersect < left[i].size() * 0.9) return false;
		if (intersect < right[i].size() * 0.9) return false;
		if (left[i].size() < 10) break;
		if (right[i].size() < 10) break;
	}
	return true;
}

std::vector<std::vector<uint8_t>> filterByOccurrenceLinkage(const std::vector<std::vector<uint8_t>>& unfiltered)
{
	assert(unfiltered.size() >= 1);
	assert(unfiltered[0].size() >= 1);
	std::vector<std::vector<std::vector<size_t>>> readsPerAllelePerSite;
	readsPerAllelePerSite.resize(unfiltered[0].size());
	for (size_t i = 0; i < unfiltered[0].size(); i++)
	{
		phmap::flat_hash_map<uint8_t, uint8_t> alleleToIndex;
		for (size_t j = 0; j < unfiltered.size(); j++)
		{
			if (unfiltered[j][i] == '-') continue;
			assert(unfiltered[j].size() == unfiltered[0].size());
			if (alleleToIndex.count(unfiltered[j][i]) == 0)
			{
				alleleToIndex[unfiltered[j][i]] = readsPerAllelePerSite[i].size();
				readsPerAllelePerSite[i].emplace_back();
			}
		}
		for (size_t j = 0; j < unfiltered.size(); j++)
		{
			if (unfiltered[j][i] == '-') continue;
			readsPerAllelePerSite[i][alleleToIndex.at(unfiltered[j][i])].emplace_back(j);
		}
		std::sort(readsPerAllelePerSite[i].begin(), readsPerAllelePerSite[i].end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
		for (size_t j = 0; j < readsPerAllelePerSite[i].size(); j++)
		{
			std::sort(readsPerAllelePerSite[i][j].begin(), readsPerAllelePerSite[i][j].end());
		}
	}
	std::vector<size_t> possibleSites;
	for (size_t i = 0; i < unfiltered[0].size(); i++)
	{
		if (readsPerAllelePerSite[i].size() < 2) continue;
		if (readsPerAllelePerSite[i][1].size() < 5) continue;
		possibleSites.emplace_back(i);
	}
	std::sort(possibleSites.begin(), possibleSites.end(), [&readsPerAllelePerSite](const size_t left, const size_t right) { return readsPerAllelePerSite[left][0].size() > readsPerAllelePerSite[right][0].size(); });
	std::vector<bool> covered;
	covered.resize(unfiltered[0].size(), false);
	for (size_t i = 0; i < possibleSites.size(); i++)
	{
		for (size_t j = i+1; j < possibleSites.size(); j++)
		{
			assert(readsPerAllelePerSite[possibleSites[j]][0].size() <= readsPerAllelePerSite[possibleSites[i]][0].size());
			if (readsPerAllelePerSite[possibleSites[j]][0].size() < readsPerAllelePerSite[possibleSites[i]][0].size() * 0.9) break;
			if (possibleSites[j] + 10 > possibleSites[i] && possibleSites[j] < possibleSites[i] + 10) continue;
			if (covered[possibleSites[i]] && covered[possibleSites[j]]) continue;
			if (allelesMatchApproximately(readsPerAllelePerSite[possibleSites[i]], readsPerAllelePerSite[possibleSites[j]]) || allelesMatchPerfectly(readsPerAllelePerSite[possibleSites[i]], readsPerAllelePerSite[possibleSites[j]]))
			{
				covered[possibleSites[i]] = true;
				covered[possibleSites[j]] = true;
			}
		}
	}
	std::vector<std::vector<uint8_t>> result = unfiltered;
	for (size_t i = 0; i < result.size(); i++)
	{
		filterVector(result[i], covered);
	}
	return result;
}

void checkPairPhasingGroupsDiscardSingletons(const std::vector<std::vector<uint8_t>>& readFakeMSABases, const size_t firstIndex, const size_t secondIndex, std::vector<std::vector<size_t>>& pairPhasingGroups)
{
	phmap::flat_hash_set<std::pair<uint8_t, uint8_t>> alleles;
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, size_t> alleleCoverage;
	for (size_t i = 0; i < readFakeMSABases.size(); i++)
	{
		alleles.emplace(readFakeMSABases[i][firstIndex], readFakeMSABases[i][secondIndex]);
		alleleCoverage[std::make_pair(readFakeMSABases[i][firstIndex], readFakeMSABases[i][secondIndex])] += 1;
	}
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, std::pair<uint8_t, uint8_t>> parent;
	for (auto pair : alleles)
	{
		if (alleleCoverage.at(pair) == 1) continue;
		parent[pair] = pair;
	}
	for (auto pair : alleles)
	{
		if (alleleCoverage.at(pair) == 1) continue;
		for (auto pair2 : alleles)
		{
			if (alleleCoverage.at(pair2) == 1) continue;
			if (pair == pair2) continue;
			if (pair.first != pair2.first && pair.second != pair2.second) continue;
			merge(parent, pair, pair2);
		}
	}
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, size_t> keyToNode;
	size_t nextNum = 0;
	for (auto pair : alleles)
	{
		if (alleleCoverage.at(pair) == 1) continue;
		if (keyToNode.count(find(parent, pair)) == 0)
		{
			keyToNode[find(parent, pair)] = nextNum;
			nextNum += 1;
		}
	}
	assert(nextNum >= 1);
	if (nextNum < 2) return;
	for (size_t i = 0; i < readFakeMSABases.size(); i++)
	{
		std::pair<uint8_t, uint8_t> key { readFakeMSABases[i][firstIndex], readFakeMSABases[i][secondIndex] };
		if (alleleCoverage.at(key) == 1)
		{
			pairPhasingGroups[i].emplace_back(nextNum);
			nextNum += 1;
		}
		else
		{
			pairPhasingGroups[i].emplace_back(keyToNode.at(find(parent, key)));
		}
	}
}

void checkPairPhasingGroups(const std::vector<const std::vector<uint8_t>*>& readFakeMSABases, const phmap::flat_hash_map<std::string, size_t>& blockSequenceCoverage, const size_t firstIndex, const size_t secondIndex, std::vector<std::vector<size_t>>& pairPhasingGroups, const size_t originalFirstIndx, const size_t originalSecondIndex)
{
	phmap::flat_hash_set<std::pair<uint8_t, uint8_t>> alleles;
	for (const auto& pair : blockSequenceCoverage)
	{
		alleles.emplace(pair.first[firstIndex], pair.first[secondIndex]);
	}
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, std::pair<uint8_t, uint8_t>> parent;
	for (auto pair : alleles)
	{
		parent[pair] = pair;
	}
	for (auto pair : alleles)
	{
		for (auto pair2 : alleles)
		{
			if (pair == pair2) continue;
			if (pair.first != pair2.first && pair.second != pair2.second) continue;
			merge(parent, pair, pair2);
		}
	}
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, size_t> keyToNode;
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, size_t> clusterCoverage;
	size_t nextNum = 0;
	for (auto pair : alleles)
	{
		if (keyToNode.count(find(parent, pair)) == 0)
		{
			keyToNode[find(parent, pair)] = nextNum;
			nextNum += 1;
		}
	}
	for (const auto& pair : blockSequenceCoverage)
	{
		clusterCoverage[find(parent, std::make_pair<uint8_t, uint8_t>(pair.first[firstIndex], pair.first[secondIndex]))] += pair.second;
	}
	for (auto pair : clusterCoverage)
	{
		if (pair.second < 3) return;
	}
	assert(nextNum >= 1);
	if (nextNum < 2) return;
	for (size_t i = 0; i < readFakeMSABases.size(); i++)
	{
		std::pair<uint8_t, uint8_t> key { (*readFakeMSABases[i])[originalFirstIndx], (*readFakeMSABases[i])[originalSecondIndex] };
		assert(parent.count(key) == 1);
		auto p = find(parent, key);
		assert(keyToNode.count(p) == 1);
		pairPhasingGroups[i].emplace_back(keyToNode.at(p));
	}
}

void checkPairPhasingGroups(const std::vector<std::vector<uint8_t>>& readFakeMSABases, const size_t firstIndex, const size_t secondIndex, std::vector<std::vector<size_t>>& pairPhasingGroups)
{
	phmap::flat_hash_set<std::pair<uint8_t, uint8_t>> alleles;
	for (size_t i = 0; i < readFakeMSABases.size(); i++)
	{
		alleles.emplace(readFakeMSABases[i][firstIndex], readFakeMSABases[i][secondIndex]);
	}
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, std::pair<uint8_t, uint8_t>> parent;
	for (auto pair : alleles)
	{
		parent[pair] = pair;
	}
	for (auto pair : alleles)
	{
		for (auto pair2 : alleles)
		{
			if (pair == pair2) continue;
			if (pair.first != pair2.first && pair.second != pair2.second) continue;
			merge(parent, pair, pair2);
		}
	}
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, size_t> keyToNode;
	size_t nextNum = 0;
	for (auto pair : alleles)
	{
		if (keyToNode.count(find(parent, pair)) == 0)
		{
			keyToNode[find(parent, pair)] = nextNum;
			nextNum += 1;
		}
	}
	assert(nextNum >= 1);
	if (nextNum < 2) return;
	for (size_t i = 0; i < readFakeMSABases.size(); i++)
	{
		std::pair<uint8_t, uint8_t> key { readFakeMSABases[i][firstIndex], readFakeMSABases[i][secondIndex] };
		pairPhasingGroups[i].emplace_back(keyToNode.at(find(parent, key)));
	}
}

std::vector<std::vector<size_t>> getPairPhasingGroupsDiscardSingletons(const std::vector<std::vector<uint8_t>>& readFakeMSABases)
{
	std::vector<std::vector<size_t>> result;
	result.resize(readFakeMSABases.size());
	for (size_t j = 0; j < readFakeMSABases[0].size(); j++)
	{
		for (size_t k = j+1; k < readFakeMSABases[0].size(); k++)
		{
			checkPairPhasingGroupsDiscardSingletons(readFakeMSABases, j, k, result);
			if (result[0].size() > 0) return result;
		}
	}
	return result;
}

std::vector<std::vector<size_t>> getPairPhasingGroupsRec(const std::vector<const std::vector<uint8_t>*>& readFakeMSABases)
{
	assert(readFakeMSABases.size() >= 1);
	std::vector<std::vector<size_t>> result;
	result.resize(readFakeMSABases.size());
	for (size_t blockj = 0; blockj < (readFakeMSABases[0]->size()+9)/10; blockj++)
	{
		for (size_t blockk = blockj; blockk < (readFakeMSABases[0]->size()+9)/10; blockk++)
		{
			phmap::flat_hash_map<std::string, size_t> blockSequenceCoverage;
			for (size_t k = 0; k < readFakeMSABases.size(); k++)
			{
				std::string sequenceHere;
				for (size_t l = 0; l < 10; l++)
				{
					if (blockj*10+l < readFakeMSABases[0]->size())
					{
						sequenceHere += (*readFakeMSABases[k])[blockj*10+l];
					}
					else
					{
						sequenceHere += 'N';
					}
				}
				for (size_t l = 0; l < 10; l++)
				{
					if (blockk*10+l < readFakeMSABases[0]->size())
					{
						sequenceHere += (*readFakeMSABases[k])[blockk*10+l];
					}
					else
					{
						sequenceHere += 'N';
					}
				}
				blockSequenceCoverage[sequenceHere] += 1;
			}
			for (size_t j = 0; j < 10; j++)
			{
				if (blockj*10 + j >= readFakeMSABases[0]->size()) break;
				for (size_t k = (blockj == blockk ? j+1 : 0); k < 10; k++)
				{
					if (blockk*10 + k >= readFakeMSABases[0]->size()) break;
					checkPairPhasingGroups(readFakeMSABases, blockSequenceCoverage, j, 10+k, result, blockj*10+j, blockk*10+k);
				}
			}
		}
	}
	if (result[0].size() == 0) return result;
	std::vector<std::vector<size_t>> linesPerCluster;
	std::vector<std::vector<const std::vector<uint8_t>*>> sitesPerLinePerCluster;
	std::map<std::vector<size_t>, size_t> valuesToCluster;
	for (size_t i = 0; i < result.size(); i++)
	{
		std::vector<size_t>* key = &result[i];
		if (valuesToCluster.count(*key) == 0)
		{
			valuesToCluster[*key] = linesPerCluster.size();
			linesPerCluster.emplace_back();
			sitesPerLinePerCluster.emplace_back();
		}
		size_t cluster = valuesToCluster.at(*key);
		linesPerCluster[cluster].emplace_back(i);
		sitesPerLinePerCluster[cluster].emplace_back(readFakeMSABases[i]);
	}
	assert(linesPerCluster.size() >= 2);
	assert(linesPerCluster.size() == sitesPerLinePerCluster.size());
	assert(linesPerCluster.size() == valuesToCluster.size());
	for (size_t i = 0; i < linesPerCluster.size(); i++)
	{
		std::vector<std::vector<size_t>> recursiveResult = getPairPhasingGroupsRec(sitesPerLinePerCluster[i]);
		assert(recursiveResult.size() == linesPerCluster[i].size());
		for (size_t j = 0; j < recursiveResult.size(); j++)
		{
			result[linesPerCluster[i][j]].insert(result[linesPerCluster[i][j]].end(), recursiveResult[j].begin(), recursiveResult[j].end());
		}
	}
	return result;
}

std::vector<std::vector<size_t>> getPairPhasingGroups(const std::vector<std::vector<uint8_t>>& readFakeMSABases)
{
	std::vector<const std::vector<uint8_t>*> basesPointers;
	for (size_t i = 0; i < readFakeMSABases.size(); i++)
	{
		basesPointers.emplace_back(&(readFakeMSABases[i]));
	}
	return getPairPhasingGroupsRec(basesPointers);
}

void insertTripletClusters(std::vector<std::vector<uint8_t>>& result, const std::vector<std::vector<uint8_t>>& readFakeMSABases, const size_t i, const size_t j, const size_t k)
{
	phmap::flat_hash_set<std::tuple<uint8_t, uint8_t, uint8_t>> triplets;
	for (size_t l = 0; l < readFakeMSABases.size(); l++)
	{
		triplets.insert(std::make_tuple(readFakeMSABases[l][i], readFakeMSABases[l][j], readFakeMSABases[l][k]));
	}
	phmap::flat_hash_map<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>> parent;
	for (auto pair : triplets)
	{
		parent[pair] = pair;
	}
	for (auto pair : triplets)
	{
		for (auto pair2 : triplets)
		{
			size_t distance = 0;
			if (std::get<0>(pair) != std::get<0>(pair2)) distance += 1;
			if (std::get<1>(pair) != std::get<1>(pair2)) distance += 1;
			if (std::get<2>(pair) != std::get<2>(pair2)) distance += 1;
			if (distance <= 1)
			{
				merge(parent, pair, pair2);
			}
		}
	}
	phmap::flat_hash_map<std::tuple<size_t, size_t, size_t>, size_t> clusterToIndex;
	size_t nextNum = 0;
	for (auto t : triplets)
	{
		if (clusterToIndex.count(find(parent, t)) == 1) continue;
		clusterToIndex[find(parent, t)] = nextNum;
		nextNum += 1;
	}
	assert(nextNum >= 2);
	assert(nextNum <= 5);
	for (size_t l = 0; l < readFakeMSABases.size(); l++)
	{
		std::tuple<size_t, size_t, size_t> key { readFakeMSABases[l][i], readFakeMSABases[l][j], readFakeMSABases[l][k] };
		assert(parent.count(key) == 1);
		assert(clusterToIndex.count(parent.at(key)) == 1);
		result[l].emplace_back(clusterToIndex.at(parent.at(key)));
	}
}

bool tripletsArePhasable(const phmap::flat_hash_map<std::string, size_t>& alleleCoverages, const size_t i, const size_t j, const size_t k)
{
	phmap::flat_hash_map<std::tuple<uint8_t, uint8_t, uint8_t>, size_t> tripletCoverage;
	for (const auto& pair : alleleCoverages)
	{
		tripletCoverage[std::make_tuple(pair.first[i], pair.first[j], pair.first[k])] += pair.second;
	}
	phmap::flat_hash_map<std::tuple<uint8_t, uint8_t, uint8_t>, std::tuple<uint8_t, uint8_t, uint8_t>> parent;
	for (auto pair : tripletCoverage)
	{
		parent[pair.first] = pair.first;
	}
	for (auto pair : tripletCoverage)
	{
		for (auto pair2 : tripletCoverage)
		{
			size_t distance = 0;
			if (std::get<0>(pair.first) != std::get<0>(pair2.first)) distance += 1;
			if (std::get<1>(pair.first) != std::get<1>(pair2.first)) distance += 1;
			if (std::get<2>(pair.first) != std::get<2>(pair2.first)) distance += 1;
			if (distance <= 1)
			{
				merge(parent, pair.first, pair2.first);
			}
		}
	}
	phmap::flat_hash_map<std::tuple<uint8_t, uint8_t, uint8_t>, size_t> clusterCoverage;
	for (auto pair : tripletCoverage)
	{
		clusterCoverage[find(parent, pair.first)] += pair.second;
	}
	if (clusterCoverage.size() < 2) return false;
	for (auto pair : clusterCoverage)
	{
		if (pair.second < 5) return false;
	}
	return true;
}

void runTripletItem(const TripletSplittingQueue::Item& item)
{
	assert(item.chunki != std::numeric_limits<size_t>::max());
	assert(item.filteredCoverages != nullptr);
	phmap::flat_hash_map<std::string, size_t> alleleCoverages;
	for (const auto& pair : (*item.filteredCoverages))
	{
		std::string allele;
		for (size_t i = 0; i < 10; i++)
		{
			if (item.chunki*10+i < item.coveredIndices->size())
			{
				allele.push_back(pair.first[item.chunki*10+i]);
			}
			else
			{
				allele.push_back('N');
			}
		}
		for (size_t i = 0; i < 10; i++)
		{
			if (item.chunkj*10+i < item.coveredIndices->size())
			{
				allele.push_back(pair.first[item.chunkj*10+i]);
			}
			else
			{
				allele.push_back('N');
			}
		}
		for (size_t i = 0; i < 10; i++)
		{
			if (item.chunkk*10+i < item.coveredIndices->size())
			{
				allele.push_back(pair.first[item.chunkk*10+i]);
			}
			else
			{
				allele.push_back('N');
			}
		}
		alleleCoverages[allele] += pair.second;
	}
	assert(alleleCoverages.size() >= 2);
	std::vector<bool> indexUsedHere;
	indexUsedHere.resize(30, false);
	{
		std::lock_guard<std::mutex> lock { *item.resultMutex };
		for (size_t i = 0; i < 10; i++)
		{
			if (item.chunki*10+i < item.coveredIndices->size()) indexUsedHere[i] = (*item.indexUsed)[item.chunki*10+i];
			if (item.chunkj*10+i < item.coveredIndices->size()) indexUsedHere[10+i] = (*item.indexUsed)[item.chunkj*10+i];
			if (item.chunkk*10+i < item.coveredIndices->size()) indexUsedHere[20+i] = (*item.indexUsed)[item.chunkk*10+i];
		}
	}
	bool didSomething = false;
	std::vector<std::vector<uint8_t>> resultHere;
	for (size_t offi = 0; offi < 10; offi++)
	{
		size_t i = item.chunki*10+offi;
		if (i >= item.coveredIndices->size()) break;
		if (indexUsedHere[offi]) continue;
		for (size_t offj = 0; offj < 10; offj++)
		{
			if (indexUsedHere[offi]) break;
			size_t j = item.chunkj*10+offj;
			if (j >= item.coveredIndices->size()) break;
			if (indexUsedHere[10+offj]) continue;
			if ((*item.coveredIndices)[i]+10 > (*item.coveredIndices)[j] && (*item.coveredIndices)[j]+10 > (*item.coveredIndices)[i]) continue;
			for (size_t offk = 0; offk < 10; offk++)
			{
				if (indexUsedHere[offi]) break;
				if (indexUsedHere[10+offj]) break;
				size_t k = item.chunkk*10+offk;
				if (k >= item.coveredIndices->size()) break;
				if (indexUsedHere[20+offk]) continue;
				if ((*item.coveredIndices)[i]+10 > (*item.coveredIndices)[k] && (*item.coveredIndices)[k]+10 > (*item.coveredIndices)[i]) continue;
				if ((*item.coveredIndices)[k]+10 > (*item.coveredIndices)[j] && (*item.coveredIndices)[j]+10 > (*item.coveredIndices)[k]) continue;
				if (!tripletsArePhasable(alleleCoverages, offi, 10+offj, 20+offk)) continue;
				didSomething = true;
				if (resultHere.size() == 0) resultHere.resize(item.readFakeMSABases->size());
				insertTripletClusters(resultHere, *item.readFakeMSABases, (*item.coveredIndices)[i], (*item.coveredIndices)[j], (*item.coveredIndices)[k]);
				indexUsedHere[offi] = true;
				indexUsedHere[10+offj] = true;
				indexUsedHere[20+offk] = true;
			}
		}
	}
	if (didSomething)
	{
		std::lock_guard<std::mutex> lock { *item.resultMutex };
		for (size_t i = 0; i < 10; i++)
		{
			if (indexUsedHere[i]) (*item.indexUsed)[item.chunki*10+i] = true;
			if (indexUsedHere[10+i]) (*item.indexUsed)[item.chunkj*10+i] = true;
			if (indexUsedHere[20+i]) (*item.indexUsed)[item.chunkk*10+i] = true;
		}
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			(*item.result)[i].insert((*item.result)[i].end(), resultHere[i].begin(), resultHere[i].end());
		}
	}
	assert((*item.countRunning) >= 1);
	(*item.countRunning) -= 1;
}

std::vector<std::vector<uint8_t>> filterByTriplets(const std::vector<std::vector<uint8_t>>& readFakeMSABases, const double estimatedAverageErrorRate, std::vector<std::shared_ptr<TripletSplittingQueue>>& tripletQueues, std::mutex& tripletQueuesMutex)
{
	const size_t maxItems = 1000000;
	assert(readFakeMSABases.size() >= 1);
	assert(readFakeMSABases[0].size() >= 1);
	std::vector<size_t> coveredIndices;
	for (size_t i = 0; i < readFakeMSABases[0].size(); i++)
	{
		phmap::flat_hash_map<uint8_t, size_t> alleles;
		for (size_t j = 0; j < readFakeMSABases.size(); j++)
		{
			alleles[readFakeMSABases[j][i]] += 1;
		}
		size_t countCovered = 0;
		for (auto pair : alleles)
		{
			if (pair.second < 3) continue;
			if (pair.second < readFakeMSABases.size() * estimatedAverageErrorRate) continue;
			if (pair.first == '-') continue;
			countCovered += 1;
		}
		if (countCovered < 2) continue;
		coveredIndices.emplace_back(i);
	}
	std::vector<std::vector<uint8_t>> result;
	result.resize(readFakeMSABases.size());
	if (coveredIndices.size() < 3) return result;
	phmap::flat_hash_map<std::string, size_t> filteredCoverages;
	for (size_t i = 0; i < readFakeMSABases.size(); i++)
	{
		std::string allele;
		for (size_t j = 0; j < coveredIndices.size(); j++)
		{
			allele.push_back(readFakeMSABases[i][coveredIndices[j]]);
		}
		filteredCoverages[allele] += 1;
	}
	std::vector<bool> indexUsed;
	indexUsed.resize(coveredIndices.size(), false);
	std::atomic<size_t> countRunning;
	countRunning = 0;
	size_t countItems = coveredIndices.size() * coveredIndices.size() * coveredIndices.size() * readFakeMSABases.size();
	std::shared_ptr<TripletSplittingQueue> queue = std::make_shared<TripletSplittingQueue>(&filteredCoverages, &coveredIndices, &readFakeMSABases, &result, &indexUsed, &countRunning);
	if (countItems > maxItems*2)
	{
		std::lock_guard<std::mutex> lock { tripletQueuesMutex };
		tripletQueues.emplace_back(queue);
	}
	TripletSplittingQueue::Item item;
	while (queue->popItem(&item))
	{
		runTripletItem(item);
	}
	while (countRunning >= 1)
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
	return result;
}

void runAlignmentItem(SNPAlignmentQueue::Item& item)
{
	assert(item.consensus != nullptr);
	assert(item.result != nullptr);
	assert(item.strings != nullptr);
	std::string& consensusSeq = *item.consensus;
	std::vector<std::vector<uint8_t>>& readSNPMSA = *item.result;
	const std::vector<std::string>& strings = *item.strings;
	size_t totalMismatches = 0;
	size_t sequenceLengthSum = 0;
	bool anyAlignmentError = false;
	for (size_t j = item.minIndex; j < item.maxIndex; j++)
	{
		sequenceLengthSum += strings[j].size();
		EdlibAlignResult result = edlibAlign(strings[j].data(), strings[j].size(), consensusSeq.data(), consensusSeq.size(), edlibNewAlignConfig(std::max(consensusSeq.size(), strings[j].size()), EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
		if (result.status != EDLIB_STATUS_OK)
		{
			anyAlignmentError = true;
			edlibFreeAlignResult(result);
			break;
		}
		totalMismatches += result.editDistance;
		size_t consensusPos = 0;
		size_t seqPos = 0;
		readSNPMSA[j].resize(consensusSeq.size(), 'N');
		assert(result.alignmentLength >= std::max(consensusSeq.size(), strings[j].size()));
		for (size_t k = 0; k < result.alignmentLength; k++)
		{
			switch(result.alignment[k])
			{
			case 0:
				assert(seqPos < strings[j].size());
				assert(consensusPos < consensusSeq.size());
				assert(strings[j][seqPos] == consensusSeq[consensusPos]);
				readSNPMSA[j][consensusPos] = strings[j][seqPos];
				seqPos += 1;
				consensusPos += 1;
				continue;
			case 1:
				assert(seqPos < strings[j].size());
				seqPos += 1;
				continue;
			case 2:
				assert(consensusPos < consensusSeq.size());
				readSNPMSA[j][consensusPos] = '-';
				consensusPos += 1;
				continue;
			case 3:
				assert(seqPos < strings[j].size());
				assert(consensusPos < consensusSeq.size());
				assert(strings[j][seqPos] != consensusSeq[consensusPos]);
				readSNPMSA[j][consensusPos] = strings[j][seqPos];
				seqPos += 1;
				consensusPos += 1;
				continue;
			default:
				assert(false);
			}
		}
		assert(seqPos == strings[j].size());
		assert(consensusPos == consensusSeq.size());
		edlibFreeAlignResult(result);
		for (size_t k = 0; k < readSNPMSA[j].size(); k++)
		{
			assert(readSNPMSA[j][k] != 'N');
		}
	}
	if (anyAlignmentError) (*item.anyAlignmentError) = true;
	(*item.totalSequenceLength) += sequenceLengthSum;
	(*item.totalMismatches) += totalMismatches;
	assert((*item.countRunning) >= 1);
	(*item.countRunning) -= 1;
}

std::tuple<std::vector<std::vector<uint8_t>>, size_t, size_t> getSNPMSA(const std::vector<std::string>& strings, std::vector<std::shared_ptr<SNPAlignmentQueue>>& alignmentQueues, std::mutex& alignmentQueuesMutex)
{
	const size_t bpsPerBlock = 2000000;
	std::vector<std::vector<uint8_t>> readSNPMSA;
	std::string consensusSeq;
	{
		phmap::flat_hash_map<std::string, size_t> sequenceCount;
		for (size_t j = 0; j < strings.size(); j++)
		{
			sequenceCount[strings[j]] += 1;
		}
		if (sequenceCount.size() >= 2)
		{
			consensusSeq = getConsensus(sequenceCount, strings.size());
		}
		else
		{
			return std::make_tuple(readSNPMSA, 0, 0);
		}
	}
	bool hasNonATCG = false;
	for (size_t j = 0; j < consensusSeq.size(); j++)
	{
		if (consensusSeq[j] != 'A' && consensusSeq[j] != 'C' && consensusSeq[j] != 'T' && consensusSeq[j] != 'G')
		{
			hasNonATCG = true;
			break;
		}
	}
	if (hasNonATCG)
	{
		std::cerr << "ERROR SNP transitive closure clustering skipped chunk with coverage " << strings.size() << " NON-ATCG IN CONSENSUS: " << consensusSeq << std::endl;
		return std::make_tuple(readSNPMSA, 0, 0);
	}
	readSNPMSA.resize(strings.size());
	std::atomic<bool> anyAlignmentError = false;
	std::atomic<size_t> totalMismatches = 0;
	std::atomic<size_t> totalSequenceLength = 0;
	std::atomic<size_t> countRunning = 0;
	size_t totalBps = strings.size() * consensusSeq.size();
	if (totalBps < bpsPerBlock*2)
	{
		SNPAlignmentQueue::Item item;
		item.strings = &strings;
		item.consensus = &consensusSeq;
		item.minIndex = 0;
		item.maxIndex = strings.size();
		item.result = &readSNPMSA;
		item.totalMismatches = &totalMismatches;
		item.totalSequenceLength = &totalSequenceLength;
		item.anyAlignmentError = &anyAlignmentError;
		item.countRunning = &countRunning;
		(*item.countRunning) += 1;
		runAlignmentItem(item);
	}
	else
	{
		std::shared_ptr<SNPAlignmentQueue> queue = std::make_shared<SNPAlignmentQueue>();
		queue->insertItems(&strings, &consensusSeq, 0, strings.size(), totalBps/bpsPerBlock, &readSNPMSA, &totalMismatches, &totalSequenceLength, &anyAlignmentError, &countRunning);
		{
			std::lock_guard<std::mutex> lock { alignmentQueuesMutex };
			alignmentQueues.emplace_back(queue);
		}
		while (true)
		{
			SNPAlignmentQueue::Item item;
			bool got = queue->popItem(&item);
			if (!got) break;
			runAlignmentItem(item);
		}
		while (countRunning >= 1)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
	}
	if (anyAlignmentError)
	{
		std::cerr << "ERROR SNP transitive closure clustering skipped chunk with coverage " << strings.size() << " EDLIB ERROR" << std::endl;
		readSNPMSA.clear();
	}
	return std::make_tuple(std::move(readSNPMSA), (size_t)totalMismatches, (size_t)totalSequenceLength);
}

std::tuple<std::vector<std::vector<uint8_t>>, size_t, size_t> getSNPMSA(const std::vector<std::string>& strings)
{
	std::vector<std::shared_ptr<SNPAlignmentQueue>> fakeQueue;
	std::mutex fakeMutex;
	return getSNPMSA(strings, fakeQueue, fakeMutex);
}

std::vector<std::vector<size_t>> getSNPTransitiveClosure(const std::vector<std::vector<uint8_t>>& readSNPMSA, const size_t edits)
{
	std::vector<size_t> parent = getFastTransitiveClosure(readSNPMSA.size(), edits, [&readSNPMSA](const size_t i, const size_t j, const size_t maxDist) { return getHammingdistance(readSNPMSA[i], readSNPMSA[j], maxDist); });
	size_t nextNum = 0;
	phmap::flat_hash_map<size_t, size_t> keyToNode;
	for (size_t j = 0; j < parent.size(); j++)
	{
		size_t key = find(parent, j);
		if (keyToNode.count(key) == 1) continue;
		keyToNode[key] = nextNum;
		nextNum += 1;
	}
	std::vector<std::vector<size_t>> result;
	result.resize(nextNum);
	for (size_t j = 0; j < parent.size(); j++)
	{
		result[keyToNode.at(find(parent, j))].push_back(j);
	}
	return result;
}

std::vector<std::vector<uint8_t>> filterMSAByCoverage(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const double estimatedAverageErrorRate)
{
	const size_t minSolidBaseCoverage = 10;
	std::vector<bool> coveredSite;
	assert(unfilteredReadFakeMSABases.size() >= 1);
	size_t consensusLength = unfilteredReadFakeMSABases[0].size();
	coveredSite.resize(consensusLength, false);
	size_t countCovered = 0;
	// {
	// 	std::lock_guard<std::mutex> lock { resultMutex };
	// 	std::cerr << "chunk with coverage " << chunkBeingDone.size() << " get covered bases time " << formatTime(startTime, getTime()) << std::endl;
	// }
	for (size_t j = 0; j < consensusLength; j++)
	{
		phmap::flat_hash_map<uint8_t, uint32_t> counts;
		for (size_t k = 0; k < unfilteredReadFakeMSABases.size(); k++)
		{
			assert(unfilteredReadFakeMSABases[k].size() == consensusLength);
			counts[unfilteredReadFakeMSABases[k][j]] += 1;
		}
		if (counts.size() < 2) continue;
		size_t countSolid = 0;
		for (auto pair : counts)
		{
			if (pair.second < minSolidBaseCoverage) continue;
			if (pair.second < unfilteredReadFakeMSABases.size() * estimatedAverageErrorRate) continue;
			if (pair.first == '-') continue;
			countSolid += 1;
		}
		if (countSolid < 2) continue;
		coveredSite[j] = true;
		countCovered += 1;
	}
	std::vector<std::vector<uint8_t>> readFakeMSABases = unfilteredReadFakeMSABases;
	for (size_t j = 0; j < unfilteredReadFakeMSABases.size(); j++)
	{
		filterVector(readFakeMSABases[j], coveredSite);
	}
	return readFakeMSABases;
}

bool sitePairIsMultinomiallySignificant(const std::vector<std::vector<uint8_t>>& unfiltered, const size_t left, const size_t right)
{
	phmap::flat_hash_map<uint8_t, size_t> leftCoverage;
	phmap::flat_hash_map<uint8_t, size_t> rightCoverage;
	phmap::flat_hash_map<uint8_t, phmap::flat_hash_map<uint8_t, size_t>> leftCoverageConditionalOnRight;
	phmap::flat_hash_map<uint8_t, phmap::flat_hash_map<uint8_t, size_t>> rightCoverageConditionalOnLeft;
	for (size_t i = 0; i < unfiltered.size(); i++)
	{
		leftCoverage[unfiltered[i][left]] += 1;
		rightCoverage[unfiltered[i][right]] += 1;
		leftCoverageConditionalOnRight[unfiltered[i][right]][unfiltered[i][left]] += 1;
		rightCoverageConditionalOnLeft[unfiltered[i][left]][unfiltered[i][right]] += 1;
	}
	// test if P(left) != P(left|right) and P(right) != P(right|left)
	// chi squared test p=0.001, null is that P(left) == P(left|right) and P(right) == P(right|left)
	// distributions are multinomial with five categories: A C G T -
	const double chiSquaredCutoff = 18.467; // 4 degrees of freedom, critical value=0.999
	for (const auto& pair : leftCoverageConditionalOnRight)
	{
		size_t totalCoverageHere = 0;
		for (const auto& pair2 : pair.second)
		{
			totalCoverageHere += pair2.second;
		}
		double chiSquaredValue = 0;
		for (const auto& pair2 : pair.second)
		{
			assert(leftCoverage.count(pair2.first) == 1);
			assert(leftCoverage.at(pair2.first) > 0);
			assert(pair2.second <= leftCoverage.at(pair2.first));
			double expectedCount = (double)leftCoverage.at(pair2.first) * (double)totalCoverageHere / (double)unfiltered.size();
			assert(expectedCount > 0);
			chiSquaredValue += ((double)pair2.second - expectedCount) * ((double)pair2.second - expectedCount) / expectedCount;
		}
		if (chiSquaredValue >= chiSquaredCutoff) return true;
	}
	for (const auto& pair : rightCoverageConditionalOnLeft)
	{
		size_t totalCoverageHere = 0;
		for (const auto& pair2 : pair.second)
		{
			totalCoverageHere += pair2.second;
		}
		double chiSquaredValue = 0;
		for (const auto& pair2 : pair.second)
		{
			assert(rightCoverage.count(pair2.first) == 1);
			assert(rightCoverage.at(pair2.first) > 0);
			assert(pair2.second <= rightCoverage.at(pair2.first));
			double expectedCount = (double)rightCoverage.at(pair2.first) * (double)totalCoverageHere / (double)unfiltered.size();
			assert(expectedCount > 0);
			chiSquaredValue += ((double)pair2.second - expectedCount) * ((double)pair2.second - expectedCount) / expectedCount;
		}
		if (chiSquaredValue >= chiSquaredCutoff) return true;
	}
	return false;
}

std::vector<std::vector<uint8_t>> filterByMultinomiallySignificantSNPs(const std::vector<std::vector<uint8_t>>& unfiltered, const double estimatedAverageErrorRate)
{
	std::vector<std::vector<uint8_t>> filtered = filterMSAByCoverage(unfiltered, estimatedAverageErrorRate);
	size_t oldSize = filtered[0].size();
	std::vector<bool> multinomialSignificant;
	multinomialSignificant.resize(filtered[0].size(), false);
	for (size_t i = 1; i < filtered[0].size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (multinomialSignificant[i] && multinomialSignificant[j]) continue;
			if (!sitePairIsMultinomiallySignificant(filtered, i, j)) continue;
			multinomialSignificant[i] = true;
			multinomialSignificant[j] = true;
		}
	}
	for (size_t i = 0; i < filtered.size(); i++)
	{
		filterVector(filtered[i], multinomialSignificant);
	}
	return filtered;
}

std::vector<std::vector<size_t>> trySNPSplittingLowerMismatchRate(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const size_t consensusLength, const double estimatedAverageErrorRate)
{
	std::vector<std::vector<uint8_t>> readFakeMSABases = filterMSAByCoverage(unfilteredReadFakeMSABases, estimatedAverageErrorRate);
	size_t biggestEdits = readFakeMSABases[0].size() * mismatchFraction + 1;
	size_t smallestEdits = readFakeMSABases[0].size() * estimatedAverageErrorRate + 1;
	assert(biggestEdits >= smallestEdits);
	assert(smallestEdits >= 1);
	for (size_t edits = biggestEdits; edits >= smallestEdits; edits--)
	{
		auto result = getSNPTransitiveClosure(readFakeMSABases, edits);
		if (result.size() >= 2) return result;
		if (edits == smallestEdits) return result;
	}
	assert(false);
}

std::vector<std::vector<size_t>> trySNPSplitting(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const size_t consensusLength, const double estimatedAverageErrorRate)
{
	std::vector<std::vector<uint8_t>> readFakeMSABases = filterMSAByCoverage(unfilteredReadFakeMSABases, estimatedAverageErrorRate);
	size_t maxEdits = readFakeMSABases[0].size() * mismatchFraction + 1;
	return getSNPTransitiveClosure(readFakeMSABases, maxEdits);
}

std::vector<std::vector<size_t>> tryMultinomiallySignificantSNPSplittingLowerMismatchRate(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const size_t consensusLength, const double estimatedAverageErrorRate)
{
	std::vector<std::vector<uint8_t>> readFakeMSABases = filterByMultinomiallySignificantSNPs(unfilteredReadFakeMSABases, estimatedAverageErrorRate);
	size_t biggestEdits = readFakeMSABases[0].size() * mismatchFraction + 1;
	size_t smallestEdits = readFakeMSABases[0].size() * estimatedAverageErrorRate + 1;
	assert(biggestEdits >= smallestEdits);
	assert(smallestEdits >= 1);
	for (size_t edits = biggestEdits; edits >= smallestEdits; edits--)
	{
		auto result = getSNPTransitiveClosure(readFakeMSABases, edits);
		if (result.size() >= 2) return result;
		if (edits == smallestEdits) return result;
	}
	assert(false);
}

std::vector<std::vector<size_t>> tryMultinomiallySignificantSNPSplitting(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const size_t consensusLength, const double estimatedAverageErrorRate)
{
	std::vector<std::vector<uint8_t>> readFakeMSABases = filterByMultinomiallySignificantSNPs(unfilteredReadFakeMSABases, estimatedAverageErrorRate);
	size_t maxEdits = readFakeMSABases[0].size() * mismatchFraction + 1;
	return getSNPTransitiveClosure(readFakeMSABases, maxEdits);
}

std::vector<std::vector<size_t>> tryOccurrenceLinkageSNPSplitting(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const size_t consensusLength, const double estimatedAverageErrorRate)
{
	std::vector<std::vector<uint8_t>> readFakeMSABases = filterByOccurrenceLinkage(unfilteredReadFakeMSABases);
	size_t maxEdits = readFakeMSABases[0].size() * mismatchFraction + 1;
	return getSNPTransitiveClosure(readFakeMSABases, maxEdits);
}

std::vector<std::vector<size_t>> tryTripletSplitting(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const size_t consensusLength, const double estimatedAverageErrorRate, std::vector<std::shared_ptr<TripletSplittingQueue>>& tripletQueues, std::mutex& tripletQueuesMutex)
{
	std::vector<std::vector<uint8_t>> readFakeMSABases = filterByTriplets(unfilteredReadFakeMSABases, estimatedAverageErrorRate, tripletQueues, tripletQueuesMutex);
	return getSNPTransitiveClosure(readFakeMSABases, 0);
}

std::vector<std::vector<size_t>> tryPairPhasingGroupSplitting(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const size_t consensusLength, const double estimatedAverageErrorRate)
{
	std::vector<std::vector<uint8_t>> readFakeMSABases = filterByOccurrenceLinkage(unfilteredReadFakeMSABases);
	std::vector<std::vector<size_t>> pairPhasingGroups = getPairPhasingGroups(readFakeMSABases);
	std::vector<size_t> parent = getFastTransitiveClosure(unfilteredReadFakeMSABases.size(), 0, [&pairPhasingGroups](const size_t i, const size_t j, const size_t maxEdits) { return (pairPhasingGroups[i] == pairPhasingGroups[j]) ? 0 : 1; });
	size_t nextNum = 0;
	phmap::flat_hash_map<size_t, size_t> keyToNode;
	for (size_t j = 0; j < parent.size(); j++)
	{
		size_t key = find(parent, j);
		if (keyToNode.count(key) == 1) continue;
		keyToNode[key] = nextNum;
		nextNum += 1;
	}
	std::vector<std::vector<size_t>> result;
	result.resize(nextNum);
	for (size_t j = 0; j < parent.size(); j++)
	{
		result[keyToNode.at(find(parent, j))].push_back(j);
	}
	return result;
}

std::vector<std::vector<size_t>> tryPairPhasingGroupSplittingDiscardSingletons(const std::vector<std::vector<uint8_t>>& unfilteredReadFakeMSABases, const size_t consensusLength, const double estimatedAverageErrorRate)
{
	std::vector<std::vector<uint8_t>> readFakeMSABases = filterByOccurrenceLinkage(unfilteredReadFakeMSABases);
	std::vector<std::vector<size_t>> pairPhasingGroups = getPairPhasingGroupsDiscardSingletons(readFakeMSABases);
	std::vector<size_t> parent = getFastTransitiveClosure(unfilteredReadFakeMSABases.size(), 0, [&pairPhasingGroups](const size_t i, const size_t j, const size_t maxEdits) { return (pairPhasingGroups[i] == pairPhasingGroups[j]) ? 0 : 1; });
	size_t nextNum = 0;
	phmap::flat_hash_map<size_t, size_t> keyToNode;
	for (size_t j = 0; j < parent.size(); j++)
	{
		size_t key = find(parent, j);
		if (keyToNode.count(key) == 1) continue;
		keyToNode[key] = nextNum;
		nextNum += 1;
	}
	std::vector<std::vector<size_t>> result;
	result.resize(nextNum);
	for (size_t j = 0; j < parent.size(); j++)
	{
		result[keyToNode.at(find(parent, j))].push_back(j);
	}
	return result;
}

void checkAndFixPalindromicChunk(std::vector<std::vector<size_t>>& clusters, const std::vector<std::pair<size_t, size_t>>& chunkBeingDone, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	if (clusters.size() == 1) return;
	phmap::flat_hash_map<size_t, size_t> clusterReverse;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> chunkInCluster;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		for (size_t j : clusters[i])
		{
			assert(j < chunkBeingDone.size());
			std::pair<size_t, size_t> chunk = chunkBeingDone[j];
			chunkInCluster[chunk] = i;
		}
	}
	bool chunksConsistent = true;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		for (size_t j : clusters[i])
		{
			assert(j < chunkBeingDone.size());
			std::pair<size_t, size_t> chunk = chunkBeingDone[j];
			std::pair<size_t, size_t> reverseChunk = chunk;
			if (chunk.first >= chunksPerRead.size() / 2)
			{
				reverseChunk.first -= chunksPerRead.size() / 2;
			}
			else
			{
				reverseChunk.first += chunksPerRead.size() / 2;
			}
			reverseChunk.second = chunksPerRead[chunk.first].size()-1-chunk.second;
			assert(chunkInCluster.count(reverseChunk) == 1);
			if (clusterReverse.count(i) == 0)
			{
				clusterReverse[i] = chunkInCluster.at(reverseChunk);
			}
			else if (clusterReverse.at(i) != chunkInCluster.at(reverseChunk))
			{
				chunksConsistent = false;
				break;
			}
		}
		if (!chunksConsistent) break;
	}
	if (!chunksConsistent)
	{
		std::cerr << "clusters not palindrome-consistent\n";
		for (size_t i = 1; i < clusters.size(); i++)
		{
			clusters[0].insert(clusters[0].end(), clusters[i].begin(), clusters[i].end());
		}
		clusters.resize(1);
	}
	else
	{
		std::cerr << "clusters are palindrome-consistent\n";
	}
}

void splitPerSNPTransitiveClosureClustering(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t minSolidBaseCoverage = 10;
	std::cerr << "splitting by SNP transitive closure clustering" << std::endl;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	std::atomic<size_t> threadsRunning;
	threadsRunning = 0;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	for (size_t i = chunksNeedProcessing.size()-1; i < chunksNeedProcessing.size(); i--)
	{
		if (chunksNeedProcessing[i].size() == 0)
		{
			std::swap(chunksNeedProcessing[i], chunksNeedProcessing.back());
			chunksNeedProcessing.pop_back();
			continue;
		}
		if (chunksNeedProcessing[i].size() < 10 || chunksNeedProcessing[i].size() < 2*minSolidBaseCoverage)
		{
			chunksDoneProcessing.emplace_back();
			std::swap(chunksDoneProcessing.back(), chunksNeedProcessing[i]);
			std::swap(chunksNeedProcessing[i], chunksNeedProcessing.back());
			chunksNeedProcessing.pop_back();
			continue;
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	auto oldChunks = chunksPerRead;
	std::vector<std::shared_ptr<SNPAlignmentQueue>> alignmentQueues;
	std::vector<std::shared_ptr<TripletSplittingQueue>> tripletQueues;
	std::mutex alignmentQueuesMutex;
	std::mutex tripletQueuesMutex;
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, &threadsRunning, &alignmentQueues, &alignmentQueuesMutex, &tripletQueues, &tripletQueuesMutex, minSolidBaseCoverage](size_t dummy)
	{
		while (true)
		{
			{
				SNPAlignmentQueue::Item item;
				bool got = false;
				{
					std::lock_guard<std::mutex> lock { alignmentQueuesMutex };
					if (alignmentQueues.size() >= 1)
					{
						got = alignmentQueues.back()->popItem(&item);
						if (alignmentQueues.back()->size() == 0) alignmentQueues.pop_back();
					}
				}
				if (got)
				{
					runAlignmentItem(item);
					continue;
				}
			}
			{
				TripletSplittingQueue::Item item;
				assert(item.chunki == std::numeric_limits<size_t>::max());
				bool got = false;
				{
					std::lock_guard<std::mutex> lock { tripletQueuesMutex };
					if (tripletQueues.size() >= 1)
					{
						got = tripletQueues.back()->popItem(&item);
						if (tripletQueues.back()->empty()) tripletQueues.pop_back();
					}
				}
				if (got)
				{
					assert(item.chunki != std::numeric_limits<size_t>::max());
					runTripletItem(item);
					continue;
				}
			}
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
			if (chunkBeingDone.size() < 10 || chunkBeingDone.size() < 2*minSolidBaseCoverage)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				chunksDoneProcessing.emplace_back();
				std::cerr << "SNP transitive closure clustering skipped chunk with coverage " << chunkBeingDone.size() << std::endl;
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "SNP transitive closure clustering begin chunk with coverage " << chunkBeingDone.size() << std::endl;
			}
			size_t lengthSum = 0;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				lengthSum += std::get<1>(chunksPerRead[chunkBeingDone[j].first][chunkBeingDone[j].second]) - std::get<0>(chunksPerRead[chunkBeingDone[j].first][chunkBeingDone[j].second]);
			}
			size_t chunkLength = lengthSum / chunkBeingDone.size();
			if (chunkLength < 100)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				chunksDoneProcessing.emplace_back();
				std::cerr << "SNP transitive closure clustering skipped short chunk with coverage " << chunkBeingDone.size() << " length " << chunkLength << std::endl;
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			bool chunkIsPalindrome = false;
			{
				phmap::flat_hash_set<std::pair<size_t, size_t>> extantChunks { chunkBeingDone.begin(), chunkBeingDone.end() };
				for (auto pair : extantChunks)
				{
					if (pair.first >= chunksPerRead.size()/2) continue;
					size_t otheri = pair.first + chunksPerRead.size()/2;
					size_t otherj = chunksPerRead[pair.first].size()-1-pair.second;
					if (extantChunks.count(std::make_pair(otheri, otherj)) == 1)
					{
						chunkIsPalindrome = true;
						break;
					}
				}
			}
			if (chunkIsPalindrome)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "SNP transitive closure clustering chunk with coverage " << chunkBeingDone.size() << " is palindrome" << std::endl;
			}
			std::vector<std::vector<uint8_t>> unfilteredReadFakeMSABases;
			size_t totalMismatches = 0;
			size_t sequenceLengthSum = 0;
			{
				std::vector<std::string> strings;
				for (size_t j = 0; j < chunkBeingDone.size(); j++)
				{
					strings.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
				}
				std::tie(unfilteredReadFakeMSABases, totalMismatches, sequenceLengthSum) = getSNPMSA(strings, alignmentQueues, alignmentQueuesMutex);
			}
			if (unfilteredReadFakeMSABases.size() == 0)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				chunksDoneProcessing.emplace_back();
				std::cerr << "ERROR SNP transitive closure clustering skipped chunk with coverage " << chunkBeingDone.size() << std::endl;
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			size_t consensusLength = unfilteredReadFakeMSABases[0].size();
			assert(unfilteredReadFakeMSABases.size() == chunkBeingDone.size());
			double estimatedAverageErrorRate = (double)totalMismatches/(double)sequenceLengthSum; // not really error rate, actually error+divergence rate
			assert(unfilteredReadFakeMSABases.size() == chunkBeingDone.size());
			std::vector<std::vector<size_t>> clusters = trySNPSplitting(unfilteredReadFakeMSABases, consensusLength, estimatedAverageErrorRate);
			if (chunkIsPalindrome) checkAndFixPalindromicChunk(clusters, chunkBeingDone, chunksPerRead);
			assert(clusters.size() >= 1);
			if (clusters.size() == 1)
			{
				clusters = tryPairPhasingGroupSplitting(unfilteredReadFakeMSABases, consensusLength, estimatedAverageErrorRate);
				if (chunkIsPalindrome) checkAndFixPalindromicChunk(clusters, chunkBeingDone, chunksPerRead);
				assert(clusters.size() >= 1);
			}
			if (clusters.size() == 1)
			{
				clusters = tryTripletSplitting(unfilteredReadFakeMSABases, consensusLength, estimatedAverageErrorRate, tripletQueues, tripletQueuesMutex);
				if (chunkIsPalindrome) checkAndFixPalindromicChunk(clusters, chunkBeingDone, chunksPerRead);
				assert(clusters.size() >= 1);
			}
			if (clusters.size() == 1 && chunkBeingDone.size() > 1000)
			{
				clusters = tryPairPhasingGroupSplittingDiscardSingletons(unfilteredReadFakeMSABases, consensusLength, estimatedAverageErrorRate);
				if (chunkIsPalindrome) checkAndFixPalindromicChunk(clusters, chunkBeingDone, chunksPerRead);
				assert(clusters.size() >= 1);
			}
			if (clusters.size() == 1)
			{
				clusters = tryOccurrenceLinkageSNPSplitting(unfilteredReadFakeMSABases, consensusLength, estimatedAverageErrorRate);
				if (chunkIsPalindrome) checkAndFixPalindromicChunk(clusters, chunkBeingDone, chunksPerRead);
				assert(clusters.size() >= 1);
			}
			if (clusters.size() == 1)
			{
				clusters = tryMultinomiallySignificantSNPSplitting(unfilteredReadFakeMSABases, consensusLength, estimatedAverageErrorRate);
				if (chunkIsPalindrome) checkAndFixPalindromicChunk(clusters, chunkBeingDone, chunksPerRead);
				assert(clusters.size() >= 1);
			}
			if (clusters.size() == 1 && chunkBeingDone.size() > 1000)
			{
				clusters = tryMultinomiallySignificantSNPSplittingLowerMismatchRate(unfilteredReadFakeMSABases, consensusLength, estimatedAverageErrorRate);
				if (chunkIsPalindrome) checkAndFixPalindromicChunk(clusters, chunkBeingDone, chunksPerRead);
				assert(clusters.size() >= 1);
			}
			if (clusters.size() == 1 && chunkBeingDone.size() > 1000)
			{
				clusters = trySNPSplittingLowerMismatchRate(unfilteredReadFakeMSABases, consensusLength, estimatedAverageErrorRate);
				if (chunkIsPalindrome) checkAndFixPalindromicChunk(clusters, chunkBeingDone, chunksPerRead);
				assert(clusters.size() >= 1);
			}
			auto endTime = getTime();
			// sort smallest last, so emplace-pop-swap puts biggest on top of chunksNeedProcessing
			assert(clusters.size() >= 1);
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.resize(clusters.size());
			for (size_t i = 0; i < clusters.size(); i++)
			{
				for (size_t j = 0; j < clusters[i].size(); j++)
				{
					assert(clusters[i][j] < chunkBeingDone.size());
					chunkResult[i].emplace_back(chunkBeingDone[clusters[i][j]]);
				}
			}
			assert(chunkResult.size() >= 1);
			bool removedNoncanonicals = false;
			if (chunkIsPalindrome)
			{
				std::vector<size_t> clusterReverse;
				clusterReverse.resize(clusters.size(), std::numeric_limits<size_t>::max());
				phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> chunkInCluster;
				for (size_t cluster = 0; cluster < clusters.size(); cluster++)
				{
					for (size_t j : clusters[cluster])
					{
						std::pair<size_t, size_t> chunk = chunkBeingDone[j];
						chunkInCluster[chunk] = cluster;
					}
				}
				for (size_t cluster = 0; cluster < clusters.size(); cluster++)
				{
					for (size_t j : clusters[cluster])
					{
						std::pair<size_t, size_t> chunk = chunkBeingDone[j];
						std::pair<size_t, size_t> reverseChunk = chunk;
						if (chunk.first >= chunksPerRead.size()/2)
						{
							reverseChunk.first -= chunksPerRead.size()/2;
						}
						else
						{
							reverseChunk.first += chunksPerRead.size()/2;
						}
						reverseChunk.second = chunksPerRead[chunk.first].size()-1-chunk.second;
						assert(chunkInCluster.count(reverseChunk) == 1);
						assert(clusterReverse[cluster] == std::numeric_limits<size_t>::max() || clusterReverse[cluster] == chunkInCluster.at(reverseChunk));
						clusterReverse[cluster] = chunkInCluster.at(reverseChunk);
					}
				}
				std::vector<size_t> noncanonicalClusters;
				for (size_t j = 0; j < clusters.size(); j++)
				{
					assert(clusterReverse[j] != std::numeric_limits<size_t>::max());
					if (clusterReverse[j] > j) noncanonicalClusters.emplace_back(j);
				}
				if (noncanonicalClusters.size() >= 1)
				{
					for (size_t j = 1; j < noncanonicalClusters.size(); j++)
					{
						assert(noncanonicalClusters[j] < noncanonicalClusters[j-1]);
					}
					removedNoncanonicals = true;
					std::lock_guard<std::mutex> lock { resultMutex };
					std::cerr << "SNP transitive closure clustering discard noncanonical clusters from chunk with coverage " << chunkBeingDone.size() << ", count noncanonicals " << noncanonicalClusters.size() << std::endl;
					for (size_t j : noncanonicalClusters)
					{
						chunksDoneProcessing.emplace_back();
						std::swap(chunksDoneProcessing.back(), chunkResult[j]);
						std::swap(chunkResult[j], chunkResult.back());
						chunkResult.pop_back();
					}
				}
				assert(chunkResult.size() >= 1);
			}
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunkResult.size() >= 2) countSplitted += 1;
				std::cerr << "SNP transitive closure clustering splitted chunk with coverage " << chunkBeingDone.size() << " to " << chunkResult.size() << " chunks, biggest chunk size: " << chunkResult[0].size() << ", sites: " << consensusLength << " estimated average error rate: " << estimatedAverageErrorRate << " time " << formatTime(startTime, endTime) << std::endl;
				if (chunkResult.size() == 1 && !removedNoncanonicals)
				{
					chunksDoneProcessing.emplace_back();
					std::swap(chunksDoneProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
				// std::cerr << "times " << formatTime(startTime, getSequencesTime) << " " << formatTime(getSequencesTime, getConsensusTime) << " " << formatTime(getConsensusTime, alignTime) << " " << formatTime(alignTime, getCoveredTime) << " " << formatTime(getCoveredTime, filterTime) << " " << formatTime(filterTime, transitiveClosureTime) << " " << formatTime(transitiveClosureTime, endTime) << std::endl;
/*				std::cerr << "sites:";
				for (size_t k = 0; k < coveredSite.size(); k++)
				{
					if (coveredSite[k]) std::cerr << " " << k;
				}
				std::cerr << std::endl;
				std::cerr << "new chunk sizes:";
				for (size_t k = 0; k < chunkResult.size(); k++)
				{
					std::cerr << " " << chunkResult[k].size();
				}
				std::cerr << std::endl;*/
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
	std::cerr << "SNP transitive closure clustering splitted " << countSplitted << " chunks" << std::endl;
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerCorrectedKmerClustering(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by corrected kmer clustering" << std::endl;
	const size_t countNeighbors = 5;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	std::atomic<size_t> threadsRunning;
	threadsRunning = 0;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	auto oldChunks = chunksPerRead;
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, &threadsRunning, countNeighbors, kmerSize](size_t dummy)
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
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				chunksDoneProcessing.emplace_back();
				std::cerr << "corrected kmer clustering skipped chunk with coverage " << chunkBeingDone.size() << std::endl;
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			ScopedCounterIncrementer threadCounter { threadsRunning };
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
				std::cerr << "corrected kmer clustering skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " due to covered " << columnsInCovered << " time " << formatTime(startTime, endTime) << std::endl;
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
				std::cerr << "corrected kmer clustering skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " due to informative " << columnsInInformative << " time " << formatTime(startTime, endTime) << std::endl;
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
			std::vector<size_t> pairwiseHammingDistances;
			pairwiseHammingDistances.resize(columns.size()+1);
			for (size_t j = 1; j < matrix.size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					size_t distance = getHammingdistance(correctedMatrix[j], correctedMatrix[k], pairwiseHammingDistances.size());
					assert(distance < pairwiseHammingDistances.size());
					pairwiseHammingDistances[distance] += 1;
				}
			}
			size_t clusteringDistance = 0;
			for (size_t j = 0; j < pairwiseHammingDistances.size(); j++)
			{
				if (pairwiseHammingDistances[j] > pairwiseHammingDistances[clusteringDistance]) clusteringDistance = j;
			}
			clusteringDistance *= 2;
			clusteringDistance += 1;
			std::vector<size_t> parent;
			for (size_t i = 0; i < matrix.size(); i++)
			{
				parent.emplace_back(i);
			}
			for (size_t j = 1; j < matrix.size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					if (find(parent, j) == find(parent, k)) continue;
					size_t distance = getHammingdistance(correctedMatrix[j], correctedMatrix[k], clusteringDistance+1);
					if (distance <= clusteringDistance) merge(parent, j, k);
				}
			}
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
				std::cerr << "corrected kmer clustering splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " clustering distance " << clusteringDistance << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
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
				std::cerr << "corrected kmer clustering splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " clustering distance " << clusteringDistance << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
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
	std::cerr << "corrected kmer clustering splitted " << countSplitted << " chunks" << std::endl;
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerCorrectedKmerPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by corrected kmer phasing" << std::endl;
	const size_t countNeighbors = 5;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	std::atomic<size_t> threadsRunning;
	threadsRunning = 0;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	auto oldChunks = chunksPerRead;
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, &threadsRunning, countNeighbors, kmerSize](size_t dummy)
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
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				chunksDoneProcessing.emplace_back();
				std::cerr << "corrected kmer skipped chunk with coverage " << chunkBeingDone.size() << std::endl;
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			ScopedCounterIncrementer threadCounter { threadsRunning };
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
			std::vector<bool> coveredCorrectedColumn;
			size_t columnsInCorrectedCovered = 0;
			coveredCorrectedColumn.resize(columns.size(), false);
			for (size_t j = 0; j < columns.size(); j++)
			{
				size_t zeros = 0;
				size_t ones = 0;
				for (size_t k = 0; k < correctedMatrix.size(); k++)
				{
					assert(correctedMatrix[k].size() == correctedMatrix[0].size());
					if (correctedMatrix[k].get(j))
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
					coveredCorrectedColumn[j] = true;
					columnsInCorrectedCovered += 1;
				}
			}
			for (size_t j = 0; j < matrix[0].size(); j++)
			{
				if (!coveredCorrectedColumn[j]) continue;
				for (auto k : phasedSoFar)
				{
					if (kmerLocation[j].first < kmerLocation[k].second+1 && kmerLocation[j].second + 1 > kmerLocation[k].first) continue;
					if (siteIsPerfectlyPhased(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						break;
					}
/*					else if (siteIsPhasedTwoVariantsThreeHaps(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						break;
					}*/
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
/*					else if (siteIsPhasedTwoVariantsThreeHaps(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						phasedSoFar.emplace_back(k);
						std::swap(notPhasedSoFar[kindex], notPhasedSoFar.back());
						notPhasedSoFar.pop_back();
					}*/
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
				std::cerr << "corrected kmer splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " corrected covered " << columnsInCorrectedCovered << " phasing sites " << numPhasingSites << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
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
				std::cerr << "corrected kmer splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " corrected covered " << columnsInCorrectedCovered << " phasing sites " << numPhasingSites << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
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
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerNearestNeighborPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by nearest neighbor phasing" << std::endl;
	const size_t countNeighbors = 5;
	const size_t countDifferences = 100;
	size_t countSplitted = 0;
	std::atomic<size_t> threadsRunning;
	threadsRunning = 0;
	std::mutex resultMutex;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	auto oldChunks = chunksPerRead;
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
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerAllelePhasingWithinChunk(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by allele phasing" << std::endl;
	size_t countSplitted = 0;
	size_t countSplittedTo = 0;
	std::mutex resultMutex;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	auto oldChunks = chunksPerRead;
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
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerPhasingKmersWithinChunk(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by phasing kmers" << std::endl;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	auto oldChunks = chunksPerRead;
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
			{
				phmap::flat_hash_set<std::pair<size_t, size_t>> removeFwForks;
				phmap::flat_hash_set<std::pair<size_t, size_t>> removeBwForks;
				for (const auto& pair : fwForks)
				{
					phmap::flat_hash_map<size_t, size_t> alleleCoverage;
					for (auto pair2 : pair.second)
					{
						alleleCoverage[pair2.second] += 1;
					}
					for (auto pair2 : alleleCoverage)
					{
						if (pair2.second < 5)
						{
							removeFwForks.insert(pair.first);
							break;
						}
						if (pair2.second < chunkBeingDone.size() * mismatchFraction)
						{
							removeFwForks.insert(pair.first);
							break;
						}
					}
				}
				for (const auto& pair : bwForks)
				{
					phmap::flat_hash_map<size_t, size_t> alleleCoverage;
					for (auto pair2 : pair.second)
					{
						alleleCoverage[pair2.second] += 1;
					}
					for (auto pair2 : alleleCoverage)
					{
						if (pair2.second < 5)
						{
							removeBwForks.insert(pair.first);
							break;
						}
						if (pair2.second < chunkBeingDone.size() * mismatchFraction)
						{
							removeBwForks.insert(pair.first);
							break;
						}
					}
				}
				for (auto key : removeFwForks)
				{
					fwForks.erase(key);
				}
				for (auto key : removeBwForks)
				{
					bwForks.erase(key);
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
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerInterchunkPhasedKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize)
{
	std::cerr << "splitting by interchunk phasing kmers" << std::endl;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
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
	auto oldChunks = chunksPerRead;
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
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		if (canonical[i]) continue;
		occurrencesPerChunk[i].clear();
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
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerDiploidChunkWithNeighbors(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "splitting by diploid chunk with neighbors" << std::endl;
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			maxChunk = std::max(maxChunk, std::get<2>(t) & maskUint64_t);
			if (!canonical[std::get<2>(t) & maskUint64_t]) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t) + 1);
			occurrencesPerChunk[std::get<2>(t)].emplace_back(i, j);
		}
	}
	auto oldChunks = chunksPerRead;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<size_t> chunkBelongsToUnitig;
	std::vector<size_t> chunkStartPosInUnitig;
	std::vector<size_t> chunkEndPosInUnitig;
	chunkBelongsToUnitig.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	chunkStartPosInUnitig.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	chunkEndPosInUnitig.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	assert(graph.unitigChunkBreakpointPositions.size() == graph.chunksInUnitig.size());
	assert(graph.unitigChunkBreakpointPositions.size() == graph.unitigLengths.size());
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		assert(graph.unitigChunkBreakpointPositions[i].size() == graph.chunksInUnitig[i].size());
		for (size_t j = 0; j < graph.chunksInUnitig[i].size(); j++)
		{
			uint64_t node = graph.chunksInUnitig[i][j];
			assert(node & firstBitUint64_t);
			size_t startPos = graph.unitigChunkBreakpointPositions[i][j].first;
			size_t endPos = graph.unitigChunkBreakpointPositions[i][j].second;
			assert(endPos > startPos);
			assert(endPos <= graph.unitigLengths[i]);
			assert((node & maskUint64_t) < chunkBelongsToUnitig.size());
			assert(chunkBelongsToUnitig[node & maskUint64_t] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[node & maskUint64_t] = i;
			chunkStartPosInUnitig[node & maskUint64_t] = startPos;
			assert(startPos < graph.unitigLengths[i]);
			chunkEndPosInUnitig[node & maskUint64_t] = endPos;
			assert(chunkEndPosInUnitig[node & maskUint64_t] <= graph.unitigLengths[i]);
		}
	}
	std::vector<std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<std::pair<size_t, size_t>>>> fwForksPerUnitig;
	std::vector<std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<std::pair<size_t, size_t>>>> bwForksPerUnitig;
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
				assert(prev & firstBitUint64_t);
				assert(curr & firstBitUint64_t);
				if (hasFwFork[prev & maskUint64_t])
				{
					std::pair<size_t, bool> currPair { curr & maskUint64_t, curr & firstBitUint64_t };
					if (currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, true))[0])
					{
						fwForksPerUnitig[prev & maskUint64_t].first.emplace_back(i, readPaths[i][j].readPartInPathnode[k-1].second);
					}
					else
					{
						assert(currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, true))[1]);
						fwForksPerUnitig[prev & maskUint64_t].second.emplace_back(i, readPaths[i][j].readPartInPathnode[k-1].second);
					}
				}
				if (hasBwFork[curr & maskUint64_t])
				{
					std::pair<size_t, bool> prevPair { prev & maskUint64_t, prev & firstBitUint64_t };
					prevPair.second = !prevPair.second;
					if (prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, false))[0])
					{
						bwForksPerUnitig[curr & maskUint64_t].first.emplace_back(i, readPaths[i][j].readPartInPathnode[k].first);
					}
					else
					{
						assert(prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, false))[1]);
						bwForksPerUnitig[curr & maskUint64_t].second.emplace_back(i, readPaths[i][j].readPartInPathnode[k].first);
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
			if (phmap::flat_hash_set<std::pair<size_t, size_t>>(fwForksPerUnitig[i].first.begin(), fwForksPerUnitig[i].first.end()).size() != fwForksPerUnitig[i].first.size()) hasFwFork[i] = false;
			if (phmap::flat_hash_set<std::pair<size_t, size_t>>(fwForksPerUnitig[i].second.begin(), fwForksPerUnitig[i].second.end()).size() != fwForksPerUnitig[i].second.size()) hasFwFork[i] = false;
			if (fwForksPerUnitig[i].first.size() < 2) hasFwFork[i] = false;
			if (fwForksPerUnitig[i].second.size() < 2) hasFwFork[i] = false;
		}
		if (hasBwFork[i])
		{
			if (intersectSize(bwForksPerUnitig[i].first, bwForksPerUnitig[i].second) >= 1)
			{
				hasBwFork[i] = false;
			}
			if (phmap::flat_hash_set<std::pair<size_t, size_t>>(bwForksPerUnitig[i].first.begin(), bwForksPerUnitig[i].first.end()).size() != bwForksPerUnitig[i].first.size()) hasBwFork[i] = false;
			if (phmap::flat_hash_set<std::pair<size_t, size_t>>(bwForksPerUnitig[i].second.begin(), bwForksPerUnitig[i].second.end()).size() != bwForksPerUnitig[i].second.size()) hasBwFork[i] = false;
			if (bwForksPerUnitig[i].first.size() < 2) hasBwFork[i] = false;
			if (bwForksPerUnitig[i].second.size() < 2) hasBwFork[i] = false;
		}
	}
	std::vector<std::pair<phmap::flat_hash_set<std::pair<size_t, size_t>>, phmap::flat_hash_set<std::pair<size_t, size_t>>>> fwForksPerUnitigSet;
	std::vector<std::pair<phmap::flat_hash_set<std::pair<size_t, size_t>>, phmap::flat_hash_set<std::pair<size_t, size_t>>>> bwForksPerUnitigSet;
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
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &sequenceIndex, &hasFwFork, &hasBwFork, &fwForksPerUnitigSet, &bwForksPerUnitigSet,  &countSplitted, &chunkBelongsToUnitig, &rawReadLengths, &graph, &chunkStartPosInUnitig, &chunkEndPosInUnitig, &canonical, kmerSize](const size_t i)
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
		if (occurrencesPerChunk[i].size() == 0)
		{
			assert(!canonical[i]);
			return;
		}
		size_t unitig = chunkBelongsToUnitig[i];
		assert(unitig < graph.unitigLengths.size());
		bool anythingToPhase = false;
		std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> maybePhaseGroups;
		if (hasFwFork[unitig])
		{
			if (fwForksPerUnitigSet[unitig].first.size() >= 2 && fwForksPerUnitigSet[unitig].second.size() >= 2)
			{
				std::vector<size_t> firstsHere;
				std::vector<size_t> secondsHere;
				bool valid = true;
				assert(chunkEndPosInUnitig[i] <= graph.unitigLengths[unitig]);
				size_t extraUntilEnd = graph.unitigLengths[unitig] - chunkEndPosInUnitig[i];
				assert(extraUntilEnd < graph.unitigLengths[unitig]);
				for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
				{
					auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
					size_t read = occurrencesPerChunk[i][j].first;
					size_t endPos = std::get<1>(t) + extraUntilEnd;
					bool matchesFirst = false;
					bool matchesSecond = false;
					for (auto pair : fwForksPerUnitigSet[unitig].first)
					{
						if (pair.first != read) continue;
						if (endPos > pair.second+100) continue;
						if (pair.second > endPos+100) continue;
						matchesFirst = true;
					}
					for (auto pair : fwForksPerUnitigSet[unitig].second)
					{
						if (pair.first != read) continue;
						if (endPos > pair.second+100) continue;
						if (pair.second > endPos+100) continue;
						matchesSecond = true;
					}
					if (matchesFirst && matchesSecond)
					{
						valid = false;
						break;
					}
					if (matchesFirst) firstsHere.emplace_back(j);
					if (matchesSecond) secondsHere.emplace_back(j);
				}
				if (firstsHere.size() >= 2 && secondsHere.size() >= 2 && valid)
				{
					maybePhaseGroups.emplace_back();
					maybePhaseGroups.back().first.insert(firstsHere.begin(), firstsHere.end());
					maybePhaseGroups.back().second.insert(secondsHere.begin(), secondsHere.end());
				}
			}
		}
		if (hasBwFork[unitig])
		{
			if (bwForksPerUnitigSet[unitig].first.size() >= 2 && bwForksPerUnitigSet[unitig].second.size() >= 2)
			{
				std::vector<size_t> firstsHere;
				std::vector<size_t> secondsHere;
				bool valid = true;
				size_t extraUntilStart = chunkStartPosInUnitig[i];
				assert(extraUntilStart < graph.unitigLengths[unitig]);
				for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
				{
					auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
					size_t read = occurrencesPerChunk[i][j].first;
					if (std::get<0>(t) < extraUntilStart) continue;
					size_t startPos = std::get<0>(t) - extraUntilStart;
					bool matchesFirst = false;
					bool matchesSecond = false;
					for (auto pair : bwForksPerUnitigSet[unitig].first)
					{
						if (pair.first != read) continue;
						if (startPos > pair.second+100) continue;
						if (pair.second > startPos+100) continue;
						matchesFirst = true;
					}
					for (auto pair : bwForksPerUnitigSet[unitig].second)
					{
						if (pair.first != read) continue;
						if (startPos > pair.second+100) continue;
						if (pair.second > startPos+100) continue;
						matchesSecond = true;
					}
					if (matchesFirst && matchesSecond)
					{
						valid = false;
						break;
					}
					if (matchesFirst) firstsHere.emplace_back(j);
					if (matchesSecond) secondsHere.emplace_back(j);
				}
				if (firstsHere.size() >= 2 && secondsHere.size() >= 2 && valid)
				{
					maybePhaseGroups.emplace_back();
					maybePhaseGroups.back().first.insert(firstsHere.begin(), firstsHere.end());
					maybePhaseGroups.back().second.insert(secondsHere.begin(), secondsHere.end());
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
		size_t nextKmerNum = kmerToNumber.size() + kmerClusterToNumber.size();
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
			for (size_t j = 0; j < occurrencesPerAlleleSite.size(); j++)
			{
				if (occurrencesPerAlleleSite[j].size() < 2) continue;
				for (size_t allele = 0; allele < occurrencesPerAlleleSite[j].size(); allele++)
				{
					for (size_t k : occurrencesPerAlleleSite[j][allele])
					{
						solidKmersPerOccurrence[k].emplace(nextKmerNum);
					}
					nextKmerNum += 1;
				}
			}
		}
		{
			std::vector<std::vector<uint8_t>> unfilteredReadFakeMSABases;
			size_t totalMismatches = 0;
			size_t sequenceLengthSum = 0;
			{
				std::vector<std::string> strings;
				for (size_t j = 0; j < chunkSequences.size(); j++)
				{
					strings.emplace_back(chunkSequences[j].toString());
				}
				std::tie(unfilteredReadFakeMSABases, totalMismatches, sequenceLengthSum) = getSNPMSA(strings);
			}
			if (unfilteredReadFakeMSABases.size() > 0)
			{
				assert(unfilteredReadFakeMSABases.size() == chunkSequences.size());
				std::vector<bool> covered;
				covered.resize(unfilteredReadFakeMSABases[0].size(), false);
				for (size_t j = 0; j < unfilteredReadFakeMSABases[0].size(); j++)
				{
					phmap::flat_hash_map<uint8_t, size_t> alleleCoverages;
					for (size_t k = 0; k < unfilteredReadFakeMSABases.size(); k++)
					{
						alleleCoverages[unfilteredReadFakeMSABases[k][j]] += 1;
					}
					size_t countCovered = 0;
					for (auto pair : alleleCoverages)
					{
						if (pair.first == '-') continue;
						if (pair.second < 3) continue;
						if (pair.second < chunkSequences.size() * mismatchFraction) continue;
						countCovered += 1;
					}
					if (countCovered >= 2)
					{
						covered[j] = true;
					}
				}
				for (size_t j = 0; j < covered.size(); j++)
				{
					if (!covered[j]) continue;
					phmap::flat_hash_map<uint8_t, size_t> alleleToIndex;
					for (size_t k = 0; k < unfilteredReadFakeMSABases.size(); k++)
					{
						if (alleleToIndex.count(unfilteredReadFakeMSABases[k][j]) == 0)
						{
							alleleToIndex[unfilteredReadFakeMSABases[k][j]] = nextKmerNum;
							nextKmerNum += 1;
						}
						solidKmersPerOccurrence[k].emplace(alleleToIndex.at(unfilteredReadFakeMSABases[k][j]));
					}
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
		phmap::flat_hash_map<size_t, std::vector<size_t>> phaseGroupsPerOccurrence;
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
				if (maybePhaseGroups[phaseGroup].first.count(j) == 0 && maybePhaseGroups[phaseGroup].second.count(j) == 0) continue;
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				const phmap::flat_hash_set<size_t> kmersHere = solidKmersPerOccurrence[j];
				if (maybePhaseGroups[phaseGroup].first.count(j) == 1)
				{
					kmersInEvenOneFirst.insert(kmersHere.begin(), kmersHere.end());
					assert(maybePhaseGroups[phaseGroup].second.count(j) == 0);
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
					assert(maybePhaseGroups[phaseGroup].second.count(j) == 1);
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
				if (maybePhaseGroups[phaseGroup].first.count(j) == 1 || maybePhaseGroups[phaseGroup].second.count(j) == 1) continue;
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
				if (maybePhaseGroups[phaseGroup].first.count(j) == 1)
				{
					assert(maybePhaseGroups[phaseGroup].second.count(j) == 0);
					assignments[j] = 0;
					continue;
				}
				if (maybePhaseGroups[phaseGroup].second.count(j) == 1)
				{
					assignments[j] = 1;
					continue;
				}
				assert(maybePhaseGroups[phaseGroup].first.count(j) == 0);
				assert(maybePhaseGroups[phaseGroup].second.count(j) == 0);
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
					if (assignments.count(j) == 1 && assignments.at(j) != 0)
					{
						possiblyValid = false;
						break;
					}
					assert(assignments.count(j) == 0 || assignments.at(j) == 0);
					assignments[j] = 0;
				}
				else
				{
					if (assignments.count(j) == 1 && assignments.at(j) != 1)
					{
						possiblyValid = false;
						break;
					}
					assert(assignments.count(j) == 0 || assignments.at(j) == 1);
					assignments[j] = 1;
				}
			}
			if (!possiblyValid) continue;
			for (auto pair : assignments)
			{
				phaseGroupsPerOccurrence[pair.first].push_back(nextGroupNum + pair.second);
			}
			nextGroupNum += 2;
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		if (phaseGroupsPerOccurrence.size() == 0)
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
					if (phaseGroupsPerOccurrence.at(j) != phaseGroupsPerOccurrence.at(k)) continue;
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
		if (groupCoverage.size() < 2) valid = false;
		if (groupCoverage.size() > 5) valid = false; // sanity check
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
	chunksPerRead = extrapolateCanonInformation(oldChunks, chunksPerRead);
}

void splitPerTrioPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const TrioKmerCounter& counter, const size_t numThreads)
{
	size_t nextNum = 0;
	std::mutex resultMutex;
	size_t countSplitted = 0;
	size_t countSplittedTo = 0;
	iterateCanonicalChunksByCoverage(chunksPerRead, numThreads, [&chunksPerRead, &rawReadLengths, &sequenceIndex, &counter, &nextNum, &resultMutex, &countSplitted, &countSplittedTo](const size_t chunkindex, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		bool allPhasable = true;
		if (occurrencesPerChunk[chunkindex].size() == 0) return;
		if (occurrencesPerChunk[chunkindex].size() < 2)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunkindex].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunkindex][j].first][occurrencesPerChunk[chunkindex][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunkindex][j].first][occurrencesPerChunk[chunkindex][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		size_t countHap1 = 0;
		size_t countHap2 = 0;
		std::vector<bool> isHapOne;
		isHapOne.resize(occurrencesPerChunk[chunkindex].size(), false);
		for (size_t j = 0; j < occurrencesPerChunk[chunkindex].size(); j++)
		{
			std::string chunkSequence = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[chunkindex][j].first, occurrencesPerChunk[chunkindex][j].second);
			std::pair<size_t, size_t> hapmers = counter.getHaplotypeMatchCounts(chunkSequence);
			if (hapmers.first > 0 && hapmers.second > 0)
			{
				allPhasable = false;
				break;
			}
			if (hapmers.first == 0 && hapmers.second == 0)
			{
				allPhasable = false;
				break;
			}
			if (hapmers.first > 0)
			{
				assert(hapmers.second == 0);
				countHap1 += 1;
				isHapOne[j] = true;
			}
			else
			{
				assert(hapmers.second > 0);
				countHap2 += 1;
			}
		}
		if (!allPhasable || countHap1 < 2 || countHap2 < 2)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunkindex].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunkindex][j].first][occurrencesPerChunk[chunkindex][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunkindex][j].first][occurrencesPerChunk[chunkindex][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			countSplitted += 1;
			countSplittedTo += 1;
			for (size_t j = 0; j < occurrencesPerChunk[chunkindex].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunkindex][j].first][occurrencesPerChunk[chunkindex][j].second]) = nextNum + (isHapOne[j] ? 1 : 0) + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunkindex][j].first][occurrencesPerChunk[chunkindex][j].second]) & firstBitUint64_t);
			}
			nextNum += 2;
			return;
		}
	});
	std::cerr << "trio phasing splitted " << countSplitted << " chunks to " << countSplittedTo << " chunks" << std::endl;
}
