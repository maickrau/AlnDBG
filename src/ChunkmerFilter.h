#ifndef ChunkmerFilter_h
#define ChunkmerFilter_h

#include <vector>
#include "MatchIndex.h"

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getChunksPerRead(const MatchIndex& matchIndex, const std::vector<size_t>& rawReadLengths, const std::vector<bool>& useTheseChunks);
std::vector<bool> getFilteredValidChunks(const MatchIndex& index, const std::vector<size_t>& readLengths, const size_t maxCoverage, const size_t repetitiveLength, const size_t repetitiveCoverage);
std::vector<std::pair<size_t, bool>> getParent(const MatchIndex& matchIndex, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead);

#endif
