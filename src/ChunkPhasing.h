#ifndef ChunkPhasing_h
#define ChunkPhasing_h

#include <vector>
#include <tuple>
#include <phmap.h>
#include "RankBitvector.h"
#include "CompressedStringIndex.h"
#include "TwobitString.h"

void checkPhasablePair(const phmap::flat_hash_map<size_t, size_t>& bwForks, const phmap::flat_hash_map<size_t, size_t>& fwForks, std::vector<std::vector<size_t>>& phaseIdentities);
std::vector<std::vector<std::vector<size_t>>> getAlleleOccurrences(const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>>& alleles, const size_t numOccurrences);
bool allelesMatchTwoVariantsThreeHaps(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right);
bool allelesMatchPerfectly(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right);
std::string getAlleleStr(const TwobitString& readSequence, const size_t readStart, const size_t readEnd, const bool fw);
size_t getAllele(const TwobitString& readSequence, const size_t readStart, const size_t readEnd, phmap::flat_hash_map<std::string, size_t>& alleleNumbers, const bool fw);
size_t getHammingdistance(const RankBitvector& left, const RankBitvector& right, const size_t maxEdits);
bool siteIsPerfectlyPhased(const std::vector<RankBitvector>& columns, const size_t left, const size_t right);
bool siteIsPhasedTwoVariantsThreeHaps(const std::vector<RankBitvector>& columns, const size_t left, const size_t right);
bool siteIsInformative(const std::vector<RankBitvector>& columns, const size_t left, const size_t right);
void filterMatrix(std::vector<RankBitvector>& matrix, const std::vector<bool>& keep);
void checkAddMismatch(std::vector<std::vector<std::pair<size_t, size_t>>>& closestNeighborAndMismatchesPerLine, std::vector<size_t>& maxEditsHerePerLine, std::vector<size_t>& countEqualToMaxEditsPerLine, const size_t addHere, const size_t addThis, const size_t mismatches, const size_t countNeighbors);
std::vector<RankBitvector> getCorrectedMatrix(const std::vector<RankBitvector>& matrix, const size_t countNeighbors);
bool rowEqual(const RankBitvector& left, const RankBitvector& right);
void sortAndInterleave(std::vector<RankBitvector>& correctedMatrix);
void splitPerCorrectedKmerPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads);
void splitPerNearestNeighborPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads);
void splitPerAllelePhasingWithinChunk(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads);
void splitPerPhasingKmersWithinChunk(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads);
void splitPerInterchunkPhasedKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize);
void splitPerDiploidChunkWithNeighbors(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize);

#endif
