#ifndef SequenceIdentitySplitter_h
#define SequenceIdentitySplitter_h

#include <vector>
#include <tuple>
#include "CompressedStringIndex.h"

void splitPerSequenceIdentityRoughly(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads);
void splitPerSequenceIdentity(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads);

#endif
