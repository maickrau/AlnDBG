#ifndef PathWalker_h
#define PathWalker_h

#include <vector>
#include <tuple>
#include <string>
#include "CompressedStringIndex.h"

void getContigPathsAndConsensuses(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t kmerSize, const std::string& outPathsFile, const std::string& outContigsFasta, const size_t numThreads);

#endif
