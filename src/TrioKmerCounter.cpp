#include <cmath>
#include "TrioKmerCounter.h"
#include "fastqloader.h"
#include "KmerIterator.h"

class TemporaryKmerCounterGroup
{
public:
	TemporaryKmerCounterGroup()
	{
		existingKmers.resize(pow(4, 11));
		existingCoverages.resize(pow(4, 11));
	}
	void addHap1Kmer(const uint64_t kmer)
	{
		hap1KmersNeedAdding.emplace_back(kmer);
		if (hap1KmersNeedAdding.size() + hap2KmersNeedAdding.size() > 1000000)
		{
			pack();
		}
	}
	void addHap2Kmer(const uint64_t kmer)
	{
		hap1KmersNeedAdding.emplace_back(kmer);
		if (hap1KmersNeedAdding.size() + hap2KmersNeedAdding.size() > 1000000)
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
			for (size_t i = 0; i < existingKmers[block].size(); i++)
			{
				if ((existingCoverages[block][i] & 0x0F) >= 4 && (existingCoverages[block][i] >> 4) == 0)
				{
					result.first.insert((existingKmers[block][i] << 22) + block);
				}
				if ((existingCoverages[block][i] & 0x0F) == 0 && (existingCoverages[block][i] >> 4) >= 4)
				{
					result.second.insert((existingKmers[block][i] << 22) + block);
				}
			}
		}
		return result;
	}
private:
	void pack()
	{
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
		for (const auto& pair : addCoverages)
		{
			phmap::flat_hash_set<uint64_t> found;
			for (size_t i = 0; i < existingKmers[pair.first].size(); i++)
			{
				if (pair.second.count(existingKmers[pair.first][i]) == 0) continue;
				found.insert(existingKmers[pair.first][i]);
				uint8_t oldCoverage = existingCoverages[pair.first][i];
				uint8_t addCoverage = pair.second.at(existingKmers[pair.first][i]);
				uint64_t hap1Coverage = (oldCoverage & 0x0F) + (addCoverage & 0x0F);
				if (hap1Coverage > 15) hap1Coverage = 15;
				uint64_t hap2Coverage = (oldCoverage >> 4) + (addCoverage >> 4);
				if (hap2Coverage > 15) hap2Coverage = 15;
				assert(hap1Coverage <= 15);
				assert(hap2Coverage <= 15);
				assert(hap1Coverage > 0 || hap2Coverage > 0);
				existingKmers[pair.first][i] = (hap2Coverage << 4) + hap1Coverage;
			}
			assert(found.size() <= pair.second.size());
			if (found.size() < pair.second.size())
			{
				for (auto pair2 : pair.second)
				{
					if (found.count(pair2.first) == 1) continue;
					existingKmers[pair.first].emplace_back(pair2.first);
					existingCoverages[pair.first].emplace_back(pair2.second);
				}
			}
		}
	}
	std::vector<std::vector<uint64_t>> existingKmers;
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
	TemporaryKmerCounterGroup temporaryCounters;
	std::cerr << "start reading file " << hap1File << std::endl;
	FastQ::streamFastqFromFile(hap1File, false, [&temporaryCounters, k, w](const FastQ& read)
	{
		iterateNSplittedStrings(read.sequence, [&temporaryCounters, k, w](const std::string& substring)
		{
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
	FastQ::streamFastqFromFile(hap2File, false, [&temporaryCounters, k, w](const FastQ& read)
	{
		iterateNSplittedStrings(read.sequence, [&temporaryCounters, k, w](const std::string& substring)
		{
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