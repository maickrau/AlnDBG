#ifndef GraphPhaser_h
#define GraphPhaser_h

#include <vector>
#include <tuple>
#include "UnitigGraph.h"
#include "MatchGroup.h"

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLinearizable(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphHapmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t resolveLength);

#endif
