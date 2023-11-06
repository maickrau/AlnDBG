#ifndef GraphCleaner_h
#define GraphCleaner_h

#include <vector>
#include <tuple>
#include "Common.h"
#include "UnitigGraph.h"

std::pair<UnitigGraph, std::vector<ReadPathBundle>> cleanUnitigGraph(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double averageHaplotypeCoverage);

#endif
