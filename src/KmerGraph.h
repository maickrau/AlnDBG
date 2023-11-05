#ifndef KmerGraph_h
#define KmerGraph_h

#include <vector>
#include <cstdint>
#include <cstddef>
#include "RankBitvector.h"
#include "MostlySparse2DHashmap.h"
#include "MatchGroup.h"

std::vector<RankBitvector> extendBreakpoints(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches);
std::vector<uint64_t> mergeSegments(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const std::vector<RankBitvector>& breakpoints, const size_t countBreakpoints);
RankBitvector getSegmentToNode(const std::vector<uint64_t>& segments, const size_t minCoverage);
std::vector<size_t> getNodeCoverage(const std::vector<uint64_t>& segments, const RankBitvector& segmentToNode, const size_t countNodes);
std::vector<size_t> getNodeLengths(const std::vector<uint64_t>& segments, const RankBitvector& segmentToNode, const std::vector<RankBitvector>& breakpoints, const size_t countNodes);
MostlySparse2DHashmap<uint8_t, size_t> getEdgeCoverages(const std::vector<size_t>& readLengths, const RankBitvector& segmentToNode, const std::vector<uint64_t>& segments, const std::vector<RankBitvector>& breakpoints, const size_t minCoverage, const std::vector<size_t>& nodeCoverage, const size_t countNodes);

#endif
