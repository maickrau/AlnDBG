#ifndef UnitigGraph_h
#define UnitigGraph_h

#include <vector>
#include <string>
#include <cstdint>
#include <cstddef>
#include "RankBitvector.h"
#include "MostlySparse2DHashmap.h"

std::tuple<std::vector<std::vector<uint64_t>>, std::vector<uint64_t>, std::vector<uint64_t>> getUnitigs(const size_t countNodes, const size_t minCoverage, const MostlySparse2DHashmap<uint8_t, size_t>& edgeCoverages);
std::pair<std::vector<size_t>, std::vector<double>> getUnitigLengthAndCoverage(const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<size_t>& nodeCoverage, const std::vector<size_t>& nodeLength);
MostlySparse2DHashmap<uint8_t, size_t> getUnitigEdgeCoverages(const std::vector<uint64_t>& unitigLeftmostNode, const std::vector<uint64_t>& unitigRightmostNode, const size_t minCoverage, const MostlySparse2DHashmap<uint8_t, size_t>& edgeCoverages);
void writeGraph(std::string outputFileName, const std::vector<double>& nodeCoverage, const std::vector<size_t>& nodeLength, const MostlySparse2DHashmap<uint8_t, size_t>& edgeCoverage, const size_t minCoverage, const size_t k);

#endif
