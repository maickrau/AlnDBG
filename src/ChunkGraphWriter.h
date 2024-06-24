#ifndef ChunkGraphWriter_h
#define ChunkGraphWriter_h

#include <vector>
#include <string>
#include <phmap.h>
#include "SparseEdgeContainer.h"
#include "ChunkUnitigGraph.h"
#include "CompressedStringIndex.h"

void writeStage(const size_t stage, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize);
std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> readChunksFromFakePathsFile(const std::string& filename);
void writeGraph(const std::string& filename, const std::vector<bool>& allowedNode, const SparseEdgeContainer& allowedEdges, const std::vector<std::vector<size_t>>& lengths, const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeCoverage, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps);
void writeUnitigGraph(const std::string& filename, const ChunkUnitigGraph& unitigGraph);
void writeUnitigPaths(const std::string& filename, const ChunkUnitigGraph& unitigGraph, const std::vector<std::vector<UnitigPath>>& readPaths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths);
void writeBidirectedUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::string& unitigSequencesFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t numThreads, const size_t kmerSize);
void writeBidirectedUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize);
void writeUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize);
std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getBidirectedChunks(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead);
void writeUnitigSequences(const std::string& filename, const std::vector<ConsensusString>& sequences);
void writeGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead);
void writeReadChunkSequences(const std::string& filename, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex);
void writeReadUnitigSequences(const std::string& filename, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t kmerSize);

#endif
