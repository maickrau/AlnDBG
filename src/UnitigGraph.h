#ifndef UnitigGraph_h
#define UnitigGraph_h

#include <vector>
#include <string>
#include <cstdint>
#include <cstddef>
#include "RankBitvector.h"
#include "MostlySparse2DHashmap.h"
#include "KmerGraph.h"
#include "Common.h"

class UnitigGraph
{
public:
	size_t nodeCount() const;
	std::vector<double> coverages;
	std::vector<size_t> lengths;
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverages;
};

void writeGraph(std::string outputFileName, const UnitigGraph& unitigGraph, const size_t minCoverage, const size_t k);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> makeUnitigGraph(const KmerGraph& kmerGraph, const std::vector<ReadPathBundle>& kmerGraphReadPaths, const size_t minCoverage);

#endif
