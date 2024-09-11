#include <cmath>
#include "TrioKmerCounter.h"
#include "fastqloader.h"
#include "KmerIterator.h"

class TemporaryKmerCounter
{
public:
	void addHap1Kmer(const uint64_t kmer)
	{
		hap1KmersNeedAdding.emplace_back(kmer);
		if (hap1KmersNeedAdding.size() + hap2KmersNeedAdding.size() > 100)
		{
			pack();
		}
	}
	void addHap2Kmer(const uint64_t kmer)
	{
		hap2KmersNeedAdding.emplace_back(kmer);
		if (hap1KmersNeedAdding.size() + hap2KmersNeedAdding.size() > 100)
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
		assert(existingKmers.size() == existingCoverages.size());
		for (size_t i = 0; i < existingKmers.size(); i++)
		{
			if ((existingCoverages[i] & 0x0F) >= 4 && (existingCoverages[i] >> 4) == 0)
			{
				result.first.insert(i);
			}
			if ((existingCoverages[i] & 0x0F) == 0 && (existingCoverages[i] >> 4) >= 4)
			{
				result.second.insert(i);
			}
		}
		return result;
	}
private:
	void pack()
	{
		if (hap1KmersNeedAdding.size() == 0 && hap2KmersNeedAdding.size() == 0) return;
		phmap::flat_hash_map<uint64_t, uint8_t> coverages;
		assert(existingKmers.size() == existingCoverages.size());
		for (size_t i = 0; i < existingKmers.size(); i++)
		{
			coverages[i] = existingCoverages[i];
		}
		for (uint64_t kmer : hap1KmersNeedAdding)
		{
			if (coverages.count(kmer) == 0)
			{
				coverages[kmer] = 0;
			}
			size_t oldCoverage = coverages.at(kmer) & 0x0F;
			if (oldCoverage < 15)
			{
				coverages[kmer] += 1;
			}
		}
		for (uint64_t kmer : hap2KmersNeedAdding)
		{
			if (coverages.count(kmer) == 0)
			{
				coverages[kmer] = 0;
			}
			size_t oldCoverage = coverages.at(kmer) >> 4;
			if (oldCoverage < 15)
			{
				coverages[kmer] += 0x10;
			}
		}
		existingKmers.clear();
		existingCoverages.clear();
		for (auto pair : coverages)
		{
			existingKmers.emplace_back(pair.first);
			existingCoverages.emplace_back(pair.second);
		}
		hap1KmersNeedAdding.clear();
		hap2KmersNeedAdding.clear();
	}
	std::vector<uint64_t> existingKmers;
	std::vector<uint8_t> existingCoverages;
	std::vector<uint64_t> hap1KmersNeedAdding;
	std::vector<uint64_t> hap2KmersNeedAdding;
};

void TrioKmerCounter::initialize(const std::string& hap1File, const std::string& hap2File, const size_t k, const size_t w)
{
	this->k = k;
	this->w = w;
	assert(hap1Kmers.size() == 0);
	assert(hap2Kmers.size() == 0);
	std::vector<TemporaryKmerCounter> temporaryCounters;
	temporaryCounters.resize(pow(4, 11));
	FastQ::streamFastqFromFile(hap1File, false, [&temporaryCounters, k, w](const FastQ& read)
	{
		iterateSyncmers(read.sequence, k, w, [&read, &temporaryCounters, k](const size_t pos)
		{
			uint64_t prefix = 0;
			uint64_t suffix = 0;
			for (size_t i = 0; i < 11; i++)
			{
				prefix <<= 2;
				prefix += charToInt(read.sequence[pos + i]);
			}
			for (size_t i = 11; i < k; i++)
			{
				suffix <<= 2;
				suffix += charToInt(read.sequence[pos + i]);
			}
			temporaryCounters[prefix].addHap1Kmer(suffix);
		});
	});
	FastQ::streamFastqFromFile(hap2File, false, [&temporaryCounters, k, w](const FastQ& read)
	{
		iterateSyncmers(read.sequence, k, w, [&read, &temporaryCounters, k](const size_t pos)
		{
			uint64_t prefix = 0;
			uint64_t suffix = 0;
			for (size_t i = 0; i < 11; i++)
			{
				prefix <<= 2;
				prefix += charToInt(read.sequence[pos + i]);
			}
			for (size_t i = 11; i < k; i++)
			{
				suffix <<= 2;
				suffix += charToInt(read.sequence[pos + i]);
			}
			temporaryCounters[prefix].addHap2Kmer(suffix);
		});
	});
	for (size_t i = 0; i < temporaryCounters.size(); i++)
	{
		size_t prefix = pow(4, 11);
		auto partialResult = temporaryCounters[i].getHaplotypeSpecificKmers();
		for (auto suffix : partialResult.first)
		{
			uint64_t kmer = (prefix << ((k-11ull)*2ull)) + suffix;
			assert(hap1Kmers.count(kmer) == 0);
			assert(hap2Kmers.count(kmer) == 0);
			hap1Kmers.insert(kmer);
		}
		for (auto suffix : partialResult.second)
		{
			uint64_t kmer = (prefix << ((k-11ull)*2ull)) + suffix;
			assert(hap1Kmers.count(kmer) == 0);
			assert(hap2Kmers.count(kmer) == 0);
			hap2Kmers.insert(kmer);
		}
	}
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