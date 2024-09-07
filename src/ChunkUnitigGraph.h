#ifndef ChunkUnitigGraph_h
#define ChunkUnitigGraph_h

#include <vector>
#include <tuple>
#include <phmap.h>
#include "TwobitString.h"
#include "SparseEdgeContainer.h"
#include "CompressedStringIndex.h"
#include "ConsensusMaker.h"

class ChunkUnitigGraph
{
public:
	std::vector<size_t> unitigLengths;
	std::vector<std::vector<std::pair<size_t, size_t>>> unitigChunkBreakpointPositions;
	std::vector<std::vector<size_t>> chunksInUnitig;
	std::vector<double> coverages;
	SparseEdgeContainer edges;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverages;
private:
};

class UnitigPath
{
public:
	std::vector<uint64_t> path;
	std::vector<std::pair<size_t, size_t>> readPartInPathnode;
	std::vector<std::pair<size_t, size_t>> chunksInPathnode;
	size_t pathLeftClip;
	size_t pathRightClip;
};

class ConsensusString
{
public:
	TwobitString bases;
	std::vector<std::pair<size_t, size_t>> Nchunks;
};

std::pair<std::vector<bool>, SparseEdgeContainer> getAllowedNodesAndEdges(const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeCoverage, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead);
std::pair<std::vector<std::vector<size_t>>, std::vector<size_t>> getLengthsAndCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead);
phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getEdgeOverlaps(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead);
phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getEdgeCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead);
std::vector<uint64_t> getUnitig(const uint64_t startNode, const std::vector<bool>& allowedNode, const SparseEdgeContainer& allowedEdges);
std::tuple<std::vector<std::vector<uint64_t>>, std::vector<size_t>, std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>, std::vector<std::vector<std::pair<size_t, size_t>>>> getUnitigs(const std::vector<bool>& allowedNode, const SparseEdgeContainer& allowedEdges, const std::vector<std::vector<size_t>>& lengths, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps);
std::vector<std::vector<UnitigPath>> getUnitigPaths(const ChunkUnitigGraph& graph, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<bool>& allowedNode, const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>& chunkLocationInUnitig);
std::vector<ConsensusString> getUnitigConsensuses(const ChunkUnitigGraph& unitigGraph, const std::vector<std::vector<UnitigPath>>& readPaths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<size_t>>& chunkLengths, const std::vector<size_t>& chunkCoverages, const std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>& chunkLocationInUnitig, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps, const size_t numThreads, const size_t kmerSize);
std::vector<double> getUnitigCoverages(const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<std::vector<size_t>>& chunkLengths, const std::vector<size_t>& chunkCoverages);
SparseEdgeContainer filterOutZEdges(const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<size_t>& unitigLengths, const SparseEdgeContainer& allowedEdges, const std::vector<std::vector<size_t>>& lengths, const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverages, const std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>& chunkLocationInUnitig, const double approxOneHapCoverage);
std::tuple<ChunkUnitigGraph, std::vector<std::vector<UnitigPath>>, std::vector<ConsensusString>> getChunkUnitigGraphInner(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const bool alsoConsensuses, const double approxOneHapCoverage, const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t numThreads, const size_t kmerSize);
std::tuple<ChunkUnitigGraph, std::vector<std::vector<UnitigPath>>, std::vector<ConsensusString>> getChunkUnitigGraphWithConsensuses(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t numThreads, const size_t kmerSize);
std::tuple<ChunkUnitigGraph, std::vector<std::vector<UnitigPath>>> getChunkUnitigGraph(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize);

#endif
