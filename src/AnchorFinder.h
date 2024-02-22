#ifndef AnchorFinder_h
#define AnchorFinder_h

#include <vector>
#include "UnitigGraph.h"
#include "Common.h"

class AnchorChain
{
public:
	std::vector<uint64_t> nodes;
	std::vector<size_t> nodeOffsets; // distance from start of first node to start of this node
	size_t ploidy;
};

class ChainPosition
{
public:
	uint64_t chain;
	int chainStartPosInRead;
	int chainEndPosInRead;
	size_t numKmerMatches;
};

std::vector<AnchorChain> getAnchorChains(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t initialAnchorMinLength, const double approxOneHapCoverage);
std::vector<std::vector<ChainPosition>> getReadChainPositions(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains);

#endif
