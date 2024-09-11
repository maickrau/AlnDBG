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
		size_t lowerBitsMask = (1ull << (lowerBitCount)) - 1;
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
	const size_t maxKmersBeforePacking = 50000000;
public:
	TemporaryKmerCounterGroup(const size_t kmerSize) :
		kmerSize(kmerSize)
	{
		existingKmers.resize(pow(4, 11));
		existingCoverages.resize(pow(4, 11));
	}
	void addHap1Kmer(const uint64_t kmer)
	{
		hap1KmersNeedAdding.emplace_back(kmer);
		if (hap1KmersNeedAdding.size() + hap2KmersNeedAdding.size() > maxKmersBeforePacking)
		{
			pack();
		}
	}
	void addHap2Kmer(const uint64_t kmer)
	{
		hap1KmersNeedAdding.emplace_back(kmer);
		if (hap1KmersNeedAdding.size() + hap2KmersNeedAdding.size() > maxKmersBeforePacking)
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
					result.first.insert((value << 22) + block);
				}
				if ((existingCoverages[block][i] & 0x0F) == 0 && (existingCoverages[block][i] >> 4) >= 4)
				{
					result.second.insert((value << 22) + block);
				}
			});
		}
		return result;
	}
private:
	void pack()
	{
		std::cerr << "pack" << std::endl;
		if (hap1KmersNeedAdding.size() == 0 && hap2KmersNeedAdding.size() == 0) return;
		phmap::flat_hash_map<size_t, phmap::flat_hash_map<uint64_t, uint8_t>> addCoverages;
		for (uint64_t kmer : hap1KmersNeedAdding)
		{
			uint64_t prefix = kmer & 0b00000000001111111111111111111111;
			uint64_t suffix = kmer >> 22;
			uint8_t& coverage = addCoverages[prefix][suffix];
			if ((coverage & 0x0F) < 15)
			{
				coverage += 1;
			}
		}
		hap1KmersNeedAdding.clear();
		for (uint64_t kmer : hap2KmersNeedAdding)
		{
			uint64_t prefix = kmer & 0b00000000001111111111111111111111;
			uint64_t suffix = kmer >> 22;
			uint8_t& coverage = addCoverages[prefix][suffix];
			if ((coverage >> 4) < 15)
			{
				coverage += 0x10;
			}
		}
		hap2KmersNeedAdding.clear();
		size_t countAdded = 0;
		for (auto& pair : addCoverages)
		{
			existingKmers[pair.first].iterateIndicesInSortedOrder([this, &pair](const size_t index, const size_t value)
			{
				if (pair.second.count(value) == 0)
				{
					pair.second[value] = existingCoverages[pair.first][index];
				}
				else
				{
					uint8_t& newCoverage = pair.second[value];
					uint8_t oldCoverage = existingCoverages[pair.first][index];
					uint8_t hap1Coverage = (newCoverage & 0x0F) + (oldCoverage & 0x0F);
					if (hap1Coverage > 15) hap1Coverage = 15;
					uint8_t hap2Coverage = (newCoverage >> 4) + (oldCoverage >> 4);
					if (hap2Coverage > 15) hap2Coverage = 15;
					newCoverage = (hap2Coverage << 4) + hap1Coverage;
				}
			});
			std::vector<std::pair<uint64_t, uint8_t>> newValues;
			newValues.insert(newValues.end(), pair.second.begin(), pair.second.end());
			assert(newValues.size() >= existingKmers[pair.first].size());
			countAdded += newValues.size() - existingKmers[pair.first].size();
			std::sort(newValues.begin(), newValues.end());
			std::vector<uint64_t> newKeys;
			existingCoverages[pair.first].clear();
			for (size_t i = 0; i < newValues.size(); i++)
			{
				newKeys.emplace_back(newValues[i].first);
				existingCoverages[pair.first].emplace_back(newValues[i].second);
			}
			existingKmers[pair.first].set(newKeys, pow(4, (kmerSize-11)));
			assert(existingCoverages[pair.first].size() == newValues.size());
			assert(existingKmers[pair.first].size() == newValues.size());
		}
		std::cerr << "added " << countAdded << " new kmers" << std::endl;
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
	FastQ::streamFastqFromFile(hap1File, false, [&temporaryCounters, k, w, &countReads](const FastQ& read)
	{
		iterateNSplittedStrings(read.sequence, [&temporaryCounters, k, w, &countReads](const std::string& substring)
		{
			countReads += 1;
			if (countReads % 100000 == 0)
			{
				std::cerr << "count reads " << countReads << std::endl;
			}
			iterateSyncmers(substring, k, w, [&substring, &temporaryCounters, k](const size_t pos)
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
	FastQ::streamFastqFromFile(hap2File, false, [&temporaryCounters, k, w, &countReads](const FastQ& read)
	{
		iterateNSplittedStrings(read.sequence, [&temporaryCounters, k, w, &countReads](const std::string& substring)
		{
			countReads += 1;
			if (countReads % 100000 == 0)
			{
				std::cerr << "count reads " << countReads << std::endl;
			}
			iterateSyncmers(substring, k, w, [&substring, &temporaryCounters, k](const size_t pos)
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