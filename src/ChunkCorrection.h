#ifndef ChunkCorrection_h
#define ChunkCorrection_h

#include <vector>
#include <tuple>
#include "CompressedStringIndex.h"

void correctChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const std::vector<size_t>& rawReadLengths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t kmerSize, const size_t numThreads);

#endif
