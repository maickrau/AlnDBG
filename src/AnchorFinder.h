#ifndef AnchorFinder_h
#define AnchorFinder_h

#include <vector>
#include "UnitigGraph.h"
#include "Common.h"

class AnchorChain
{
public:
	std::vector<uint64_t> nodes;
	std::vector<size_t> nodeOffsets;
	size_t ploidy;
};

std::vector<AnchorChain> getAnchorChains(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage);

#endif
