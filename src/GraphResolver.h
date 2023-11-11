#ifndef GraphResolver_h
#define GraphResolver_h

#include <tuple>
#include <vector>
#include "Common.h"
#include "UnitigGraph.h"

std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolveSimpleStructures(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double averageOneHaplotypeCoverage);

#endif
