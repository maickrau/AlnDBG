#ifndef KmerGraph_h
#define KmerGraph_h

#include <vector>
#include <cstdint>
#include <cstddef>
#include "RankBitvector.h"
#include "MostlySparse2DHashmap.h"
#include "MatchGroup.h"

class KmerGraph
{
public:
	size_t nodeCount() const;
	std::vector<uint64_t> readSegmentPaths;
	std::vector<size_t> coverages;
	std::vector<size_t> lengths;
	RankBitvector segmentToNode;
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverages;
};

KmerGraph makeKmerGraph(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t minCoverage);

#endif
