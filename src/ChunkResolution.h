#ifndef ChunkResolution_h
#define ChunkResolution_h

#include <vector>
#include <tuple>
#include <string>
#include "CompressedStringIndex.h"

size_t findBiggestSmallerIndex(const std::vector<size_t>& nums, const size_t value);
void splitPerLength(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double differenceFraction, const size_t differenceConstant, const size_t kmerSize, const size_t numThreads);
std::vector<size_t> getMinHashes(const std::string& sequence, const size_t k, const size_t count);
void splitPerFirstLastKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads);
void splitPerBaseCounts(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads);
void splitPerMinHashes(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads);
void resolveVerySmallNodes(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const size_t maxResolveSize, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize);
void resolveSemiAmbiguousUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize);
void resolveBetweenLongUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t longNodeThreshold, const size_t maxDistance, const size_t kmerSize);
void resolveUnambiguouslyResolvableUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize);
void resolveTinyNodesRecklessly(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize);
void removeTinyProblemNodes(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize);
void resolveBetweenTangles(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t longNodeThreshold, const size_t kmerSize);
void resolveBetweenTanglesAllowGaps(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t longNodeThreshold, const size_t kmerSize);
void expandResolvableChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize, const size_t expandedSize);
void contextResolve(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t expandedSize);
void resegmentChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize);
void expandChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t expandedSize);
void addGapChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize);
void expandChunksUntilSolids(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize, const size_t expandedSize);
void fragmentChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const size_t kmerSize);
void mergeNonexistentChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead);
void removeBadShortHighCoverageChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize);
void resplitFalselyMergedChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& oldChunks);
void fixYForks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize);

#endif
