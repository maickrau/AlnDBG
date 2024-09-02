#ifndef ChunkExtractor_h
#define ChunkExtractor_h

#include <vector>
#include <tuple>
#include "CompressedStringIndex.h"

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getGapEnclosingChunksPerRead(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const size_t k, const size_t windowSize, const size_t gapSize);
std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getMinimizerBoundedChunksPerRead(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const size_t k, const size_t windowSize);
void addMissingPiecesBetweenChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize);
void readFilesAndAddToSequenceIndex(const std::string& filename, FastaCompressor::CompressedStringIndex& sequenceIndex, std::vector<size_t>& readBasepairLengths, const size_t numThreads);

#endif
