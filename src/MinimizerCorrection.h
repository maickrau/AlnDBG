#ifndef MinimizerCorrection_h
#define MinimizerCorrection_h

#include <vector>
#include <tuple>
#include "CompressedStringIndex.h"

std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>> getCorrectedMinimizers(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t kmerSize, const size_t windowSize, const size_t numThreads);
std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>> getChunksPerReadFromCorrectedMinimizers(const std::vector<std::vector<size_t>>& correctedMinimizers, const std::vector<std::vector<size_t>>& minimizerKmers, const std::vector<size_t>& rawReadLengths, const size_t kmerSize);

#endif
