#ifndef KmerGraph_h
#define KmerGraph_h

#include <vector>
#include <cstdint>
#include <cstddef>
#include "RankBitvector.h"
#include "MostlySparse2DHashmap.h"
#include "MatchGroup.h"
#include "Common.h"

class KmerGraph
{
public:
	size_t nodeCount() const;
	std::vector<size_t> lengths;
};

std::pair<KmerGraph, std::vector<ReadPathBundle>> makeKmerGraph(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t minCoverage, const size_t numThreads);

#endif
