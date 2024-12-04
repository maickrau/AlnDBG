#ifndef TrioKmerCounter_h
#define TrioKmerCounter_h

#include <string>
#include <phmap.h>
#include <tuple>

class TrioKmerCounter
{
public:
	void initialize(const std::string& hap1File, const std::string& hap2File, const size_t k, const size_t w);
	std::pair<size_t, size_t> getHaplotypeMatchCounts(const std::string& str) const;
	size_t hap1KmerCount() const;
	size_t hap2KmerCount() const;
	bool notEmpty() const;
private:
	size_t k;
	size_t w;
	phmap::flat_hash_set<uint64_t> hap1Kmers;
	phmap::flat_hash_set<uint64_t> hap2Kmers;
};

#endif
