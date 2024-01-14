#ifndef ChunkmerFilter_h
#define ChunkmerFilter_h

#include <vector>
#include "MatchIndex.h"

std::vector<bool> getFilteredValidChunks(const MatchIndex& index, const std::vector<size_t>& readLengths, const size_t maxCoverage, const size_t repetitiveLength, const size_t repetitiveCoverage);

#endif
