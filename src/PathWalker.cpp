#include <limits>
#include <iostream>
#include "ChunkGraphWriter.h"
#include "PathWalker.h"
#include "Common.h"
#include "phmap.h"

double getAverageLongNodeCoverage(const ChunkUnitigGraph& graph, const size_t minLength)
{
	double totalLength = 0;
	double totalCoverage = 0;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < minLength) continue;
		totalLength += graph.unitigLengths[i];
		totalCoverage += graph.unitigLengths[i] * graph.coverages[i];
	}
	return totalCoverage / totalLength;
}

phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>> getNodeSpanningTriplets(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double estimatedSingleCopyCoverage)
{
	std::vector<bool> canBeSpanned;
	canBeSpanned.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.edges.getEdges(std::make_pair(i, true)).size() < 2) continue;
		if (graph.edges.getEdges(std::make_pair(i, false)).size() != graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		canBeSpanned[i] = true;
	}
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> spanningTripletCoverage;
	spanningTripletCoverage.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				if (!canBeSpanned[readPaths[i][j].path[k] & maskUint64_t]) continue;
				uint64_t prev = readPaths[i][j].path[k-1];
				uint64_t after = readPaths[i][j].path[k+1];
				if (readPaths[i][j].path[k] & firstBitUint64_t)
				{
				}
				else
				{
					std::swap(prev, after);
					prev ^= firstBitUint64_t;
					after ^= firstBitUint64_t;
				}
				spanningTripletCoverage[readPaths[i][j].path[k] & maskUint64_t][std::make_pair(prev, after)] += 1;
			}
		}
	}
	phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>> result;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (!canBeSpanned[i]) continue;
		size_t unitigCopyCount = graph.edges.getEdges(std::make_pair(i, true)).size();
		assert(unitigCopyCount >= 2);
		phmap::flat_hash_set<uint64_t> coveredPredecessors;
		phmap::flat_hash_set<uint64_t> coveredSuccessors;
		phmap::flat_hash_set<std::pair<size_t, size_t>> solids;
		for (auto pair : spanningTripletCoverage[i])
		{
			if (pair.second < 2) continue;
			coveredPredecessors.insert(pair.first.first);
			coveredSuccessors.insert(pair.first.second);
			solids.insert(pair.first);
		}
		if (solids.size() > unitigCopyCount) continue;
		if (spanningTripletCoverage[i].size() > unitigCopyCount+1) continue;
		if (solids.size() == unitigCopyCount-1 && spanningTripletCoverage[i].size() == unitigCopyCount)
		{
			for (auto pair : spanningTripletCoverage[i])
			{
				coveredPredecessors.insert(pair.first.first);
				coveredSuccessors.insert(pair.first.second);
				solids.insert(pair.first);
			}
		}
		if (solids.size() < 2 && unitigCopyCount == 2)
		{
			if (solids.size() == 1 && spanningTripletCoverage[i].size() == 1)
			{
				phmap::flat_hash_set<uint64_t> uncoveredPredecessor;
				phmap::flat_hash_set<uint64_t> uncoveredSuccessor;
				for (auto edge : graph.edges.getEdges(std::make_pair(i, true)))
				{
					uncoveredSuccessor.insert(edge.first + (edge.second ? firstBitUint64_t : 0));
				}
				for (auto edge : graph.edges.getEdges(std::make_pair(i, false)))
				{
					uncoveredPredecessor.insert(edge.first + (edge.second ? 0 : firstBitUint64_t));
				}
				for (auto pair : solids)
				{
					assert(uncoveredPredecessor.count(pair.first) == 1);
					assert(uncoveredSuccessor.count(pair.second) == 1);
					uncoveredPredecessor.erase(pair.first);
					uncoveredSuccessor.erase(pair.second);
				}
				if (uncoveredSuccessor.size() == 1 && uncoveredPredecessor.size() == 1)
				{
					std::pair<uint64_t, uint64_t> newUncoveredTriplet { *uncoveredPredecessor.begin(), *uncoveredSuccessor.begin() };
					solids.insert(newUncoveredTriplet);
					coveredPredecessors.insert(newUncoveredTriplet.first);
					coveredSuccessors.insert(newUncoveredTriplet.second);
				}
			}
			if (solids.size() == 0 && spanningTripletCoverage[i].size() == 2)
			{
				for (auto pair : spanningTripletCoverage[i])
				{
					coveredPredecessors.insert(pair.first.first);
					coveredSuccessors.insert(pair.first.second);
					solids.insert(pair.first);
				}
			}
		}
		if (coveredPredecessors.size() != unitigCopyCount) continue;
		if (coveredSuccessors.size() != unitigCopyCount) continue;
		result[i].insert(result[i].end(), solids.begin(), solids.end());
	}
	return result;
}

std::vector<bool> estimateUniqueUnitigs(const ChunkUnitigGraph& graph, const size_t longUnitigThreshold, const double estimatedSingleCopyCoverage)
{
	std::vector<bool> result;
	result.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < longUnitigThreshold) continue;
		if (graph.coverages[i] > 1.5 * estimatedSingleCopyCoverage) continue;
		result[i] = true;
	}
	return result;
}

uint64_t getNextFromTriplet(const phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>>& triplets, const uint64_t previous, const uint64_t current)
{
	if (triplets.count(current & maskUint64_t) == 0) return std::numeric_limits<size_t>::max();
	if (current & firstBitUint64_t)
	{
		for (auto pair : triplets.at(current & maskUint64_t))
		{
			if (pair.first == previous)
			{
				return pair.second;
			}
		}
	}
	else
	{
		for (auto pair : triplets.at(current & maskUint64_t))
		{
			if (pair.second == (previous ^ firstBitUint64_t))
			{
				return pair.first ^ firstBitUint64_t;
			}
		}
	}
	return std::numeric_limits<size_t>::max();
}

std::vector<uint64_t> followTripletPath(const ChunkUnitigGraph& graph, const uint64_t startUnitig, const phmap::flat_hash_map<size_t, std::vector<std::pair<uint64_t, uint64_t>>>& triplets)
{
	std::vector<uint64_t> result;
	result.emplace_back(startUnitig);
	while (true)
	{
		if (graph.edges.getEdges(std::make_pair(result.back() & maskUint64_t, result.back() & firstBitUint64_t)).size() == 1)
		{
			auto nextNode = graph.edges.getEdges(std::make_pair(result.back() & maskUint64_t, result.back() & firstBitUint64_t))[0];
			uint64_t nextTripletNode = getNextFromTriplet(triplets, result.back(), nextNode.first + (nextNode.second ? firstBitUint64_t : 0));
			if (nextTripletNode != std::numeric_limits<size_t>::max())
			{
				result.emplace_back(nextNode.first + (nextNode.second ? firstBitUint64_t : 0));
				result.emplace_back(nextTripletNode);
				continue;
			}
		}
		if (triplets.count(result.back() & maskUint64_t) == 1)
		{
			if (result.size() >= 2)
			{
				uint64_t nextTripletNode = getNextFromTriplet(triplets, result[result.size()-2], result.back());
				if (nextTripletNode != std::numeric_limits<size_t>::max())
				{
					result.emplace_back(nextTripletNode);
					continue;
				}
			}
		}
		break;
	}
	return result;
}

std::vector<std::vector<uint64_t>> getTripletExtendedPaths(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const double estimatedSingleCopyCoverage)
{
	auto triplets = getNodeSpanningTriplets(graph, readPaths, estimatedSingleCopyCoverage);
	std::vector<std::vector<uint64_t>> result;
	std::vector<bool> alreadyInPath;
	alreadyInPath.resize(graph.unitigLengths.size(), false);
	for (const auto& pair : triplets)
	{
		for (const auto pair2 : pair.second)
		{
			if (alreadyInPath[pair2.first & maskUint64_t])
			{
				assert(alreadyInPath[pair2.second & maskUint64_t]);
				continue;
			}
			assert(!alreadyInPath[pair2.second & maskUint64_t]);
			auto bwPath = followTripletPath(graph, pair2.first ^ firstBitUint64_t, triplets);
			auto fwPath = followTripletPath(graph, pair2.second, triplets);
			assert(bwPath.size() >= 1);
			assert(fwPath.size() >= 1);
			std::reverse(bwPath.begin(), bwPath.end());
			for (size_t i = 0; i < bwPath.size(); i++)
			{
				bwPath[i] ^= firstBitUint64_t;
			}
			bwPath.emplace_back(pair.first + firstBitUint64_t);
			bwPath.insert(bwPath.end(), fwPath.begin(), fwPath.end());
			for (auto node : bwPath)
			{
				alreadyInPath[node & maskUint64_t] = true;
			}
			size_t pathLength = 0;
			for (size_t i = 0; i < bwPath.size(); i++)
			{
				pathLength += graph.unitigLengths[bwPath[i] & maskUint64_t];
				if (i > 0)
				{
					auto unitigpairkey = canon(std::make_pair(bwPath[i-1] & maskUint64_t, bwPath[i-1] & firstBitUint64_t), std::make_pair(bwPath[i] & maskUint64_t, bwPath[i] & firstBitUint64_t));
					std::pair<uint64_t, uint64_t> unitigkey { unitigpairkey.first.first + (unitigpairkey.first.second ? firstBitUint64_t : 0), unitigpairkey.second.first + (unitigpairkey.second.second ? firstBitUint64_t : 0) };
					size_t overlap = graph.edgeOverlaps.at(unitigkey);
					assert(pathLength >= overlap);
					pathLength -= overlap;
				}
			}
			if (pathLength < 100000) continue;
			result.emplace_back(bwPath);
		}
	}
	return result;
}

void addUniqueNodePaths(std::vector<std::vector<uint64_t>>& paths, const std::vector<bool>& isUniqueUnitig)
{
	std::vector<bool> alreadyInPath;
	alreadyInPath.resize(isUniqueUnitig.size(), false);
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 0; j < paths[i].size(); j++)
		{
			alreadyInPath[paths[i][j] & maskUint64_t] = true;
		}
	}
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (alreadyInPath[i]) continue;
		if (!isUniqueUnitig[i]) continue;
		paths.emplace_back();
		paths.back().emplace_back(i + firstBitUint64_t);
	}
}

void getContigPathsAndConsensuses(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t kmerSize, const std::string& outPathsFile, const std::string& outContigsFasta)
{
	const size_t longUnitigThreshold = 100000;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	auto fixedChunksPerRead = getBidirectedChunks(chunksPerRead);
	std::tie(graph, readPaths) = getChunkUnitigGraph(fixedChunksPerRead, approxOneHapCoverage, kmerSize);
	double estimatedSingleCopyCoverage = getAverageLongNodeCoverage(graph, longUnitigThreshold);
	std::cerr << "estimated single copy coverage " << estimatedSingleCopyCoverage << std::endl;
	std::vector<bool> isUniqueUnitig = estimateUniqueUnitigs(graph, longUnitigThreshold, estimatedSingleCopyCoverage);
	std::cerr << "unique unitigs: ";
	for (size_t i = 0; i < isUniqueUnitig.size(); i++)
	{
		if (!isUniqueUnitig[i]) continue;
		std::cerr << i << ",";
	}
	std::cerr << std::endl;
	auto paths = getTripletExtendedPaths(graph, readPaths, estimatedSingleCopyCoverage);
	std::cerr << "step 1 paths:" << std::endl;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::cerr << "step 1 path " << i << ": ";
		for (auto node : paths[i])
		{
			std::cerr << ((node & firstBitUint64_t) ? ">" : "<") << (node & maskUint64_t);
		}
		std::cerr << std::endl;
		for (auto node : paths[i])
		{
			std::cerr << (node & maskUint64_t) << ",";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	addUniqueNodePaths(paths, isUniqueUnitig);
	std::cerr << "step 2 paths:" << std::endl;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::cerr << "step 2 path " << i << ": ";
		for (auto node : paths[i])
		{
			std::cerr << ((node & firstBitUint64_t) ? ">" : "<") << (node & maskUint64_t);
		}
		std::cerr << std::endl;
		for (auto node : paths[i])
		{
			std::cerr << (node & maskUint64_t) << ",";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
}
