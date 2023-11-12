#ifndef GraphPhaser_h
#define GraphPhaser_h

#include <vector>
#include <tuple>
#include "UnitigGraph.h"
#include "MatchGroup.h"

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> getGraphPhaseBlockNodes(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphBubbles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& phaseBlocks);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLinearizable(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage);

#endif
