#ifndef MultiplexResolverCaller_h
#define MultiplexResolverCaller_h

#include "Common.h"
#include "UnitigGraph.h"

std::pair<UnitigGraph, std::vector<ReadPathBundle>> runMBGMultiplexResolution(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const size_t maxk);

#endif
