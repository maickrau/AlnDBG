#include <cmath>
#include "VariableWidthIntVector.h"
#include "TrioKmerCounter.h"
#include "fastqloader.h"
#include "KmerIterator.h"

class EliasFanoEncodedVector
{
public:
	void set(const std::vector<size_t>& sortedValues, const size_t maxValue)
	{
		if (sortedValues.size() == 0)
		{
			lowerBits.resize(0);
			return;
		}
		size_t lowerBitCount = ceil(log2((double)maxValue / (double)sortedValues.size()));
		if (lowerBitCount == log2(maxValue)) lowerBitCount -= 1;
		size_t lowerBitsMask = (1ull << (lowerBitCount)) - 1ull;
		assert(lowerBitCount < log2(maxValue));
		upperBitCount = log2(maxValue) - lowerBitCount;
		assert(upperBitCount >= 1);
		assert(upperBitCount < log2(maxValue));
		lowerBits.resize(0);
		lowerBits.setWidth(lowerBitCount);
		lowerBits.resize(sortedValues.size());
		upperBits.clear();
		size_t currentUpperBits = 0;
		for (size_t i = 0; i < sortedValues.size(); i++)
		{
			size_t upperBitsHere = sortedValues[i] >> lowerBitCount;
			assert(upperBitsHere >= currentUpperBits);
			for (size_t j = currentUpperBits; j < upperBitsHere; j++)
			{
				upperBits.emplace_back(1);
			}
			currentUpperBits = upperBitsHere;
			upperBits.emplace_back(0);
			lowerBits.set(i, sortedValues[i] & lowerBitsMask);
		}
		assert(currentUpperBits < pow(2, upperBitCount));
		for (size_t j = currentUpperBits; j < pow(2, upperBitCount); j++)
		{
			upperBits.emplace_back(1);
		}
		assert(upperBits.size() == pow(2, upperBitCount) + sortedValues.size());
	}
	template <typename F>
	void iterateIndicesInSortedOrder(F callback) const
	{
		if (lowerBits.size() == 0) return;
		size_t index = 0;
		size_t valueUpperBits = 0;
		for (size_t i = 0; i < upperBits.size(); i++)
		{
			if (upperBits[i])
			{
				valueUpperBits += 1;
				continue;
			}
			callback(index, valueUpperBits * pow(2, lowerBits.width()) + lowerBits.get(index));
			index += 1;
		}
		assert(index == lowerBits.size());
		assert(valueUpperBits == pow(2, upperBitCount));
	}
	size_t size() const
	{
		return lowerBits.size();
	}
private:
	size_t upperBitCount;
	std::vector<bool> upperBits;
	FastaCompressor::VariableWidthIntVector lowerBits;
};

class TemporaryKmerCounterGroup
{
	const size_t prefixLength = 5;
	const size_t maxKmersBeforePacking = 50000000;
public:
	TemporaryKmerCounterGroup(const size_t kmerSize) :
		kmerSize(kmerSize)
	{
		existingKmers.resize(pow(4, prefixLength));
		existingCoverages.resize(pow(4, prefixLength));
		hap1KmersNeedAdding.reserve(maxKmersBeforePacking);
		hap2KmersNeedAdding.reserve(maxKmersBeforePacking);
	}
	void addHap1Kmer(const uint64_t kmer)
	{
		hap1KmersNeedAdding.emplace_back(kmer);
		if (hap1KmersNeedAdding.size() >= maxKmersBeforePacking)
		{
			pack();
		}
	}
	void addHap2Kmer(const uint64_t kmer)
	{
		hap2KmersNeedAdding.emplace_back(kmer);
		if (hap2KmersNeedAdding.size() >= maxKmersBeforePacking)
		{
			pack();
		}
	}
	std::pair<phmap::flat_hash_set<uint64_t>, phmap::flat_hash_set<uint64_t>> getHaplotypeSpecificKmers()
	{
		pack();
		assert(hap1KmersNeedAdding.size() == 0);
		assert(hap2KmersNeedAdding.size() == 0);
		std::pair<phmap::flat_hash_set<uint64_t>, phmap::flat_hash_set<uint64_t>> result;
		for (size_t block = 0; block < existingCoverages.size(); block++)
		{
			assert(existingKmers[block].size() == existingCoverages[block].size());
			existingKmers[block].iterateIndicesInSortedOrder([this, &result, block](const size_t i, const size_t value)
			{
				if ((existingCoverages[block][i] & 0x0F) >= 4 && (existingCoverages[block][i] >> 4) == 0)
				{
					result.first.insert((value << (prefixLength*2)) + block);
				}
				if ((existingCoverages[block][i] & 0x0F) == 0 && (existingCoverages[block][i] >> 4) >= 4)
				{
					result.second.insert((value << (prefixLength*2)) + block);
				}
			});
		}
		return result;
	}
private:
	void mergeCoveragesToBucket(phmap::flat_hash_map<size_t, uint8_t>& addCoverages, const size_t bucketIndex, std::vector<uint64_t>& newKeysTmpArray)
	{
		existingKmers[bucketIndex].iterateIndicesInSortedOrder([this, &bucketIndex, &addCoverages](const size_t index, const size_t value)
		{
			if (addCoverages.count(value) == 0)
			{
				addCoverages[value] = existingCoverages[bucketIndex][index];
			}
			else
			{
				uint8_t& newCoverage = addCoverages[value];
				uint8_t oldCoverage = existingCoverages[bucketIndex][index];
				uint8_t hap1Coverage = (newCoverage & 0x0F) + (oldCoverage & 0x0F);
				if (hap1Coverage > 15) hap1Coverage = 15;
				uint8_t hap2Coverage = (newCoverage >> 4) + (oldCoverage >> 4);
				if (hap2Coverage > 15) hap2Coverage = 15;
				newCoverage = (hap2Coverage << 4) + hap1Coverage;
			}
		});
		assert(addCoverages.size() >= existingKmers[bucketIndex].size());
		if (addCoverages.size() > existingKmers[bucketIndex].size())
		{
			newKeysTmpArray.reserve(addCoverages.size());
			for (auto pair : addCoverages)
			{
				newKeysTmpArray.emplace_back(pair.first);
			}
			assert(newKeysTmpArray.size() == addCoverages.size());
			std::sort(newKeysTmpArray.begin(), newKeysTmpArray.end());
			existingKmers[bucketIndex].set(newKeysTmpArray, pow(4, (kmerSize-prefixLength)));
			existingCoverages[bucketIndex].resize(newKeysTmpArray.size(), 0);
			newKeysTmpArray.clear();
		}
		existingKmers[bucketIndex].iterateIndicesInSortedOrder([this, &bucketIndex, &addCoverages](const size_t index, const size_t value)
		{
			assert(addCoverages.count(value) == 1);
			existingCoverages[bucketIndex][index] = addCoverages.at(value);
		});
	}
	void pack()
	{
		std::cerr << "pack" << std::endl;
		const size_t prefixMask = ((1ull << (2ull*prefixLength)) - 1ull);
		if (hap1KmersNeedAdding.size() == 0 && hap2KmersNeedAdding.size() == 0) return;
		std::sort(hap1KmersNeedAdding.begin(), hap1KmersNeedAdding.end(), [prefixMask](const size_t left, const size_t right) { return (left & prefixMask) < (right & prefixMask); });
		std::sort(hap2KmersNeedAdding.begin(), hap2KmersNeedAdding.end(), [prefixMask](const size_t left, const size_t right) { return (left & prefixMask) < (right & prefixMask); });
		size_t lastPrefix = 0;
		size_t hap1Index = 0;
		size_t hap2Index = 0;
		phmap::flat_hash_map<size_t, uint8_t> addCoverages; // keep these outside loop to reduce malloc/free, less memory fragmentation
		std::vector<uint64_t> newKeysTmpArray; // keep these outside loop to reduce malloc/free, less memory fragmentation
		while (hap1Index < hap1KmersNeedAdding.size() || hap2Index < hap2KmersNeedAdding.size())
		{
			size_t prefixHere;
			if (hap2Index == hap2KmersNeedAdding.size())
			{
				assert(hap1Index < hap1KmersNeedAdding.size());
				prefixHere = hap1KmersNeedAdding[hap1Index] & prefixMask;
			}
			else if (hap1Index == hap1KmersNeedAdding.size())
			{
				assert(hap2Index < hap2KmersNeedAdding.size());
				prefixHere = hap2KmersNeedAdding[hap2Index] & prefixMask;
			}
			else
			{
				assert(hap1Index < hap1KmersNeedAdding.size());
				assert(hap2Index < hap2KmersNeedAdding.size());
				prefixHere = std::min(hap1KmersNeedAdding[hap1Index] & prefixMask, hap2KmersNeedAdding[hap2Index] & prefixMask);
			}
			assert(prefixHere > lastPrefix || lastPrefix == 0);
			lastPrefix = prefixHere;
			while (hap1Index < hap1KmersNeedAdding.size() && (hap1KmersNeedAdding[hap1Index] & prefixMask) == prefixHere)
			{
				uint64_t suffix = hap1KmersNeedAdding[hap1Index] >> (prefixLength*2ull);
				uint8_t& coverage = addCoverages[suffix];
				if ((coverage & 0x0F) < 15)
				{
					coverage += 1;
				}
				hap1Index += 1;
			}
			while (hap2Index < hap2KmersNeedAdding.size() && (hap2KmersNeedAdding[hap2Index] & prefixMask) == prefixHere)
			{
				uint64_t suffix = hap2KmersNeedAdding[hap2Index] >> (prefixLength*2ull);
				uint8_t& coverage = addCoverages[suffix];
				if ((coverage >> 4) < 15)
				{
					coverage += 0x10;
				}
				hap2Index += 1;
			}
			mergeCoveragesToBucket(addCoverages, prefixHere, newKeysTmpArray);
			addCoverages.clear();
		}
		assert(addCoverages.size() == 0);
		hap1KmersNeedAdding.clear();
		hap2KmersNeedAdding.clear();
	}
	size_t kmerSize;
	std::vector<EliasFanoEncodedVector> existingKmers;
	std::vector<std::vector<uint8_t>> existingCoverages;
	std::vector<uint64_t> hap1KmersNeedAdding;
	std::vector<uint64_t> hap2KmersNeedAdding;
};

template <typename F>
void iterateNSplittedStrings(const std::string& str, F callback)
{
	size_t lastValid = 0;
	for (size_t i = 0; i < str.size(); i++)
	{
		if (str[i] == 'N')
		{
			if (i > lastValid)
			{
				callback(str.substr(lastValid, i-lastValid));
			}
			lastValid = i+1;
		}
	}
	if (lastValid > 0 && lastValid < str.size())
	{
		callback(str.substr(lastValid));
	}
	else if (lastValid == 0)
	{
		callback(str);
	}
}

void TrioKmerCounter::initialize(const std::string& hap1File, const std::string& hap2File, const size_t k, const size_t w)
{
	this->k = k;
	this->w = w;
	assert(hap1Kmers.size() == 0);
	assert(hap2Kmers.size() == 0);
	TemporaryKmerCounterGroup temporaryCounters { k };
	std::cerr << "start reading file " << hap1File << std::endl;
	size_t countReads = 0;
	size_t countSyncmers = 0;
	FastQ::streamFastqFromFile(hap1File, false, [&temporaryCounters, k, w, &countReads, &countSyncmers](const FastQ& read)
	{
		iterateNSplittedStrings(read.sequence, [&temporaryCounters, k, w, &countReads, &countSyncmers](const std::string& substring)
		{
			countReads += 1;
			if (countReads % 100000 == 0)
			{
				std::cerr << "count reads " << countReads << std::endl;
			}
			iterateSyncmers(substring, k, w, [&substring, &temporaryCounters, k, &countSyncmers](const size_t pos)
			{
				uint64_t kmer = 0;
				for (size_t i = 0; i < k; i++)
				{
					kmer <<= 2;
					kmer += charToInt(substring[pos + i]);
				}
				temporaryCounters.addHap1Kmer(kmer);
			});
		});
	});
	std::cerr << "start reading file " << hap2File << std::endl;
	countReads = 0;
	FastQ::streamFastqFromFile(hap2File, false, [&temporaryCounters, k, w, &countReads, &countSyncmers](const FastQ& read)
	{
		iterateNSplittedStrings(read.sequence, [&temporaryCounters, k, w, &countReads, &countSyncmers](const std::string& substring)
		{
			countReads += 1;
			if (countReads % 100000 == 0)
			{
				std::cerr << "count reads " << countReads << std::endl;
			}
			iterateSyncmers(substring, k, w, [&substring, &temporaryCounters, k, &countSyncmers](const size_t pos)
			{
				uint64_t kmer = 0;
				for (size_t i = 0; i < k; i++)
				{
					kmer <<= 2;
					kmer += charToInt(substring[pos + i]);
				}
				temporaryCounters.addHap2Kmer(kmer);
			});
		});
	});
	std::cerr << "merge buckets" << std::endl;
	std::tie(hap1Kmers, hap2Kmers) = temporaryCounters.getHaplotypeSpecificKmers();
}

size_t TrioKmerCounter::hap1KmerCount() const
{
	return hap1Kmers.size();
}

size_t TrioKmerCounter::hap2KmerCount() const
{
	return hap2Kmers.size();
}

std::pair<size_t, size_t> TrioKmerCounter::getHaplotypeMatchCounts(const std::string& str) const
{
	size_t hap1Matches = 0;
	size_t hap2Matches = 0;
	iterateSyncmers(str, k, w, [this, &str, &hap1Matches, &hap2Matches](const size_t pos)
	{
		uint64_t kmer = 0;
		for (size_t i = 0; i < k; i++)
		{
			kmer <<= 2;
			kmer += charToInt(str[pos + i]);
		}
		if (hap1Kmers.count(kmer) == 1) hap1Matches += 1;
		if (hap2Kmers.count(kmer) == 1) hap2Matches += 1;
	});
	return std::make_pair(hap1Matches, hap2Matches);
}

bool TrioKmerCounter::notEmpty() const
{
	return hap1KmerCount() > 0 || hap2KmerCount() > 0;
}