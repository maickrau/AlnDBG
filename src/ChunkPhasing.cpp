#include <cassert>
#include <iostream>
#include "ChunkPhasing.h"
#include "Common.h"
#include "SequenceHelper.h"
#include "UnionFind.h"
#include "KmerIterator.h"
#include "ChunkUnitigGraph.h"
#include "TransitiveClosure.h"

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