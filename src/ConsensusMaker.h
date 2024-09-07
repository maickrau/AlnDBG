#ifndef ConsensusMaker_h
#define ConsensusMaker_h

#include <vector>
#include <string>
#include <tuple>
#include <phmap.h>
#include "CompressedStringIndex.h"

void removeKmer(std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& kmersPerString, std::pair<size_t, size_t> removeThis);
std::vector<std::pair<size_t, size_t>> getConsensusSolidKmerPath(std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& kmersPerString);
std::string pickMostCentralString(const std::vector<std::string>& options, const phmap::flat_hash_map<std::string, size_t>& pieceCounts);
std::string getConsensusFromSolidKmers(const phmap::flat_hash_map<std::string, size_t>& sequenceCount, const size_t totalCount);
std::string getConsensusPickArbitrary(const phmap::flat_hash_map<std::string, size_t>& sequenceCount, const size_t totalCount);
std::string getConsensus(const phmap::flat_hash_map<std::string, size_t>& sequenceCount, const size_t totalCount);
std::string getConsensusSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const phmap::flat_hash_map<size_t, std::vector<std::tuple<size_t, bool, size_t>>>& readAnchorPoints, const size_t start, const size_t end, const size_t kmerSize);

#endif
