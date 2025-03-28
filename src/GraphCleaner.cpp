#include <iostream>
#include "GraphCleaner.h"
#include "Common.h"

std::vector<std::pair<uint64_t, std::vector<std::vector<size_t>>>> getUnitigForkReads(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& unitigPaths)
{
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> forkAlleleToIndex;
	phmap::flat_hash_map<uint64_t, size_t> forkToIndex;
	std::vector<std::pair<uint64_t, std::vector<std::vector<size_t>>>> result;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.edges.getEdges(std::make_pair(i, true)).size() >= 2)
		{
			size_t forkNumber = forkToIndex.size();
			size_t alleleNumber = 0;
			result.emplace_back();
			result.back().first = i + firstBitUint64_t;
			forkToIndex[i + firstBitUint64_t] = forkNumber;
			forkAlleleToIndex.emplace_back();
			assert(forkAlleleToIndex.size() == forkNumber+1);
			for (auto edge : graph.edges.getEdges(std::make_pair(i, true)))
			{
				forkAlleleToIndex[forkNumber][edge.first + (edge.second ? firstBitUint64_t : 0)] = alleleNumber;
				result.back().second.emplace_back();
				alleleNumber += 1;
			}
		}
		if (graph.edges.getEdges(std::make_pair(i, false)).size() >= 2)
		{
			size_t forkNumber = forkToIndex.size();
			size_t alleleNumber = 0;
			result.emplace_back();
			result.back().first = i;
			forkToIndex[i] = forkNumber;
			forkAlleleToIndex.emplace_back();
			assert(forkAlleleToIndex.size() == forkNumber+1);
			for (auto edge : graph.edges.getEdges(std::make_pair(i, false)))
			{
				forkAlleleToIndex[forkNumber][edge.first + (edge.second ? firstBitUint64_t : 0)] = alleleNumber;
				result.back().second.emplace_back();
				alleleNumber += 1;
			}
		}
	}
	for (size_t i = 0; i < unitigPaths.size(); i++)
	{
		for (size_t j = 0; j < unitigPaths[i].size(); j++)
		{
			for (size_t k = 1; k < unitigPaths[i][j].path.size(); k++)
			{
				uint64_t prev = unitigPaths[i][j].path[k-1];
				if (forkToIndex.count(prev) == 1)
				{
					size_t forkIndex = forkToIndex.at(prev);
					uint64_t curr = unitigPaths[i][j].path[k];
					assert(forkAlleleToIndex[forkIndex].count(curr) == 1);
					result[forkIndex].second[forkAlleleToIndex[forkIndex].at(curr)].emplace_back(i);
				}
				uint64_t currRev = unitigPaths[i][j].path[k] ^ firstBitUint64_t;
				if (forkToIndex.count(currRev) == 1)
				{
					size_t forkIndex = forkToIndex.at(currRev);
					uint64_t prevRev = unitigPaths[i][j].path[k-1] ^ firstBitUint64_t;
					assert(forkAlleleToIndex[forkIndex].count(prevRev) == 1);
					result[forkIndex].second[forkAlleleToIndex[forkIndex].at(prevRev)].emplace_back(i);
				}
			}
		}
	}
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		bool valid = true;
		for (size_t j = 0; j < result[i].second.size(); j++)
		{
			std::sort(result[i].second[j].begin(), result[i].second[j].end());
			for (size_t k = 1; k < result[i].second[j].size(); k++)
			{
				if (result[i].second[j][k-1] == result[i].second[j][k])
				{
					valid = false;
					break;
				}
			}
			if (!valid) break;
		}
		if (valid)
		{
			for (size_t j = 1; j < result[i].second.size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					if (intersectSize(result[i].second[j], result[i].second[k]) > 0)
					{
						valid = false;
						break;
					}
				}
				if (!valid) break;
			}
		}
		if (!valid)
		{
			std::swap(result[i], result.back());
			result.pop_back();
		}
	}
	return result;
}

bool forkAllelesMatch(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right, const size_t minMatch)
{
	assert(left.size() == right.size());
	phmap::flat_hash_set<size_t> uniqueMatchesFound;
	for (size_t i = 0; i < left.size(); i++)
	{
		size_t uniqueMatch = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < right.size(); j++)
		{
			if (intersectSize(left[i], right[j]) == 0) continue;
			if (intersectSize(left[i], right[j]) < minMatch) return false;
			if (uniqueMatch != std::numeric_limits<size_t>::max()) return false;
			uniqueMatch = j;
		}
		if (uniqueMatchesFound.count(uniqueMatch) == 1) return false;
		uniqueMatchesFound.insert(uniqueMatch);
	}
	return true;
}

phmap::flat_hash_set<uint64_t> getSolidForks(const std::vector<std::pair<uint64_t, std::vector<std::vector<size_t>>>>& forkReads, const size_t minCoverage)
{
	std::vector<size_t> checkableIndices;
	std::vector<std::vector<size_t>> forksPerRead;
	for (size_t i = 0; i < forkReads.size(); i++)
	{
		bool canCheckThis = true;
		for (size_t j = 0; j < forkReads[i].second.size(); j++)
		{
			if (forkReads[i].second[j].size() < minCoverage)
			{
				canCheckThis = false;
				break;
			}
		}
		if (canCheckThis)
		{
			for (size_t j = 0; j < forkReads[i].second.size(); j++)
			{
				for (size_t k : forkReads[i].second[j])
				{
					while (forksPerRead.size() <= k) forksPerRead.emplace_back();
					forksPerRead[k].emplace_back(checkableIndices.size());
				}
			}
			checkableIndices.emplace_back(i);
		}
	}
	std::vector<bool> solid;
	solid.resize(checkableIndices.size(), false);
	for (size_t i = 0; i < checkableIndices.size(); i++)
	{
		phmap::flat_hash_set<size_t> possiblyMatchingForks;
		{
			phmap::flat_hash_set<size_t> readsHere;
			for (size_t j = 0; j < forkReads[checkableIndices[i]].second.size(); j++)
			{
				for (auto k : forkReads[checkableIndices[i]].second[j])
				{
					readsHere.emplace(k);
				}
			}
			for (auto k : readsHere)
			{
				possiblyMatchingForks.insert(forksPerRead[k].begin(), forksPerRead[k].end());
			}
		}
		for (auto j : possiblyMatchingForks)
		{
			if (j >= i) continue;
			if (solid[i] && solid[j]) continue;
			if (forkReads[checkableIndices[i]].second.size() != forkReads[checkableIndices[j]].second.size()) continue;
			if (!forkAllelesMatch(forkReads[checkableIndices[i]].second, forkReads[checkableIndices[j]].second, minCoverage)) continue;
			solid[i] = true;
			solid[j] = true;
		}
	}
	phmap::flat_hash_set<uint64_t> result;
	for (size_t i = 0; i < checkableIndices.size(); i++)
	{
		if (!solid[i]) continue;
		result.emplace(forkReads[checkableIndices[i]].first);
	}
	return result;
}

phmap::flat_hash_set<uint64_t> getAcceptableForks(const std::vector<std::pair<uint64_t, std::vector<std::vector<size_t>>>>& forkReads, const phmap::flat_hash_set<uint64_t>& solids, const size_t minCoverage)
{
	std::vector<size_t> solidIndices;
	std::vector<size_t> checkableIndices;
	for (size_t i = 0; i < forkReads.size(); i++)
	{
		if (solids.count(forkReads[i].first) == 1)
		{
			solidIndices.emplace_back(i);
			continue;
		}
		bool canCheckThis = true;
		for (size_t j = 0; j < forkReads[i].second.size(); j++)
		{
			if (forkReads[i].second[j].size() < minCoverage)
			{
				canCheckThis = false;
				break;
			}
		}
		if (canCheckThis) checkableIndices.emplace_back(i);
	}
	std::vector<bool> acceptable;
	acceptable.resize(checkableIndices.size(), false);
	for (size_t i = 0; i < checkableIndices.size(); i++)
	{
		for (size_t j : solidIndices)
		{
			if (forkReads[checkableIndices[i]].second.size() != forkReads[j].second.size()) continue;
			if (!forkAllelesMatch(forkReads[checkableIndices[i]].second, forkReads[j].second, minCoverage)) continue;
			acceptable[i] = true;
			break;
		}
	}
	phmap::flat_hash_set<uint64_t> result;
	for (size_t i = 0; i < checkableIndices.size(); i++)
	{
		if (!acceptable[i]) continue;
		result.emplace(forkReads[checkableIndices[i]].first);
	}
	return result;
}

bool hasOtherHigherCoverageEdge(const uint64_t fork, const uint64_t edge, const ChunkUnitigGraph& graph)
{
	std::pair<size_t, bool> forkpair { fork & maskUint64_t, fork & firstBitUint64_t };
	std::pair<size_t, bool> edgepair { edge & maskUint64_t, edge & firstBitUint64_t };
	assert(graph.edges.hasEdge(forkpair, edgepair));
	auto keypair = canon(forkpair, edgepair);
	std::pair<uint64_t, uint64_t> key { keypair.first.first + (keypair.first.second ? firstBitUint64_t : 0), keypair.second.first + (keypair.second.second ? firstBitUint64_t : 0) };
	size_t compareCoverage = graph.edgeCoverages.at(key);
	for (auto edge : graph.edges.getEdges(forkpair))
	{
		if (edge == edgepair) continue;
		keypair = canon(forkpair, edge);
		std::pair<uint64_t, uint64_t> key { keypair.first.first + (keypair.first.second ? firstBitUint64_t : 0), keypair.second.first + (keypair.second.second ? firstBitUint64_t : 0) };
		if (graph.edgeCoverages.at(key) >= compareCoverage * 2) return true;
	}
	return false;
}

double estimateCoverage(const ChunkUnitigGraph& graph)
{
	double sum = 0;
	double div = 0;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < 100000) continue;
		sum += graph.unitigLengths[i] * graph.coverages[i];
		div += graph.unitigLengths[i];
	}
	return div/sum;
}

phmap::flat_hash_set<uint64_t> getValidBubbleForks(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths)
{
	phmap::flat_hash_set<size_t> relevantNode;
	phmap::flat_hash_map<size_t, std::pair<uint64_t, uint64_t>> bubbleContainingNodes;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.edges.getEdges(std::make_pair(i, true)).size() != 1) continue;
		if (graph.edges.getEdges(std::make_pair(i, false)).size() != 1) continue;
		std::pair<size_t, bool> fwEdge;
		std::pair<size_t, bool> bwEdge;
		fwEdge = graph.edges.getEdges(std::make_pair(i, true))[0];
		bwEdge = graph.edges.getEdges(std::make_pair(i, false))[0];
		if (graph.edges.getEdges(fwEdge).size() != 2) continue;
		if (graph.edges.getEdges(bwEdge).size() != 2) continue;
		if (graph.edges.getEdges(reverse(fwEdge)).size() != 2) continue;
		if (graph.edges.getEdges(reverse(bwEdge)).size() != 2) continue;
		std::pair<size_t, bool> otherNode { std::numeric_limits<size_t>::max(), true };
		for (auto edge : graph.edges.getEdges(reverse(bwEdge)))
		{
			if (edge.first == i) continue;
			assert(otherNode.first == std::numeric_limits<size_t>::max());
			otherNode = edge;
		}
		if (otherNode.first == std::numeric_limits<size_t>::max()) continue;
		if (graph.edges.getEdges(otherNode).size() != 1) continue;
		if (graph.edges.getEdges(reverse(otherNode)).size() != 1) continue;
		if (graph.edges.getEdges(otherNode)[0] != fwEdge) continue;
		relevantNode.insert(fwEdge.first);
		relevantNode.insert(bwEdge.first);
		assert(bubbleContainingNodes.count(i) == 0);
		bubbleContainingNodes[i] = std::make_pair(bwEdge.first + (bwEdge.second ? 0 : firstBitUint64_t), fwEdge.first + (fwEdge.second ? 0 : firstBitUint64_t));
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> tripletsPerUnitig;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				if (relevantNode.count(readPaths[i][j].path[k] & maskUint64_t) == 0) continue;
				uint64_t prevNode = readPaths[i][j].path[k-1];
				uint64_t nextNode = readPaths[i][j].path[k+1];
				if (readPaths[i][j].path[k] & firstBitUint64_t)
				{
				}
				else
				{
					std::swap(prevNode, nextNode);
					prevNode ^= firstBitUint64_t;
					nextNode ^= firstBitUint64_t;
				}
				tripletsPerUnitig[readPaths[i][j].path[k] & maskUint64_t][std::make_pair(prevNode, nextNode)] += 1;
			}
		}
	}
	phmap::flat_hash_set<size_t> unitigsWithValidTriplets;
	for (const auto& pair : tripletsPerUnitig)
	{
		if (pair.second.size() != 2) continue;
		bool valid = true;
		size_t totalCoverage = 0;
		for (auto pair2 : pair.second)
		{
			if (pair2.second < 2) valid = false;
			totalCoverage += pair2.second;
		}
		if (!valid && totalCoverage < 6) continue;
		phmap::flat_hash_set<uint64_t> prevs;
		phmap::flat_hash_set<uint64_t> nexts;
		for (auto pair2 : pair.second)
		{
			prevs.insert(pair2.first.first);
			nexts.insert(pair2.first.second);
		}
		assert(prevs.size() >= 1);
		assert(nexts.size() >= 1);
		assert(prevs.size() <= 2);
		assert(nexts.size() <= 2);
		if (prevs.size() < 2) continue;
		if (nexts.size() < 2) continue;
		unitigsWithValidTriplets.insert(pair.first);
	}
	phmap::flat_hash_set<uint64_t> result;
	for (auto pair : bubbleContainingNodes)
	{
		uint64_t bwNodePointingInwards = pair.second.first;
		uint64_t fwNodePointingInwards = pair.second.second;
		if (unitigsWithValidTriplets.count(bwNodePointingInwards & maskUint64_t) == 0) continue;
		if (unitigsWithValidTriplets.count(fwNodePointingInwards & maskUint64_t) == 0) continue;
		result.insert(bwNodePointingInwards);
		result.insert(fwNodePointingInwards);
	}
	return result;
}

void cleanTips(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const double maxRemoveCoverage, const size_t kmerSize)
{
	std::cerr << "cleaning tips" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	assert(chunksPerRead.size()%2 == 0);
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize, chunksPerRead.size()/2);
	double reestimatedCoverage = estimateCoverage(graph);
	if (reestimatedCoverage < 1) reestimatedCoverage = approxOneHapCoverage;
	phmap::flat_hash_set<uint64_t> solidFork;
	phmap::flat_hash_set<uint64_t> acceptableFork;
	{
		std::vector<std::pair<uint64_t, std::vector<std::vector<size_t>>>> forkReads = getUnitigForkReads(graph, readPaths);
		solidFork = getSolidForks(forkReads, std::min(approxOneHapCoverage, reestimatedCoverage) * 0.5);
		acceptableFork = getAcceptableForks(forkReads, solidFork, 3);
		phmap::flat_hash_set<uint64_t> bubbleForks = getValidBubbleForks(graph, readPaths);
		acceptableFork.insert(bubbleForks.begin(), bubbleForks.end());
	}
	phmap::flat_hash_set<size_t> removeChunks;
	size_t removeUnitigCount = 0;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.coverages[i] >= maxRemoveCoverage) continue;
		if (graph.unitigLengths[i] >= 50000) continue;
		if (graph.edges.getEdges(std::make_pair(i, true)).size() >= 1)
		{
			bool allGood = true;
			for (auto otherpair : graph.edges.getEdges(std::make_pair(i, true)))
			{
				if (otherpair.first == i) continue;
				uint64_t other = otherpair.first + (otherpair.second ? firstBitUint64_t : 0);
				other ^= firstBitUint64_t;
				if (solidFork.count(other) == 0)
				{
					if (acceptableFork.count(other) == 0)
					{
						if (hasOtherHigherCoverageEdge(other, i, graph))
						{
							continue;
						}
					}
				}
				allGood = false;
				break;
			}
			if (!allGood) continue;
		}
		if (graph.edges.getEdges(std::make_pair(i, false)).size() >= 1)
		{
			bool allGood = true;
			for (auto otherpair : graph.edges.getEdges(std::make_pair(i, false)))
			{
				if (otherpair.first == i) continue;
				uint64_t other = otherpair.first + (otherpair.second ? firstBitUint64_t : 0);
				other ^= firstBitUint64_t;
				if (solidFork.count(other) == 0)
				{
					if (acceptableFork.count(other) == 0)
					{
						if (hasOtherHigherCoverageEdge(other, i + firstBitUint64_t, graph))
						{
							continue;
						}
					}
				}
				allGood = false;
				break;
			}
			if (!allGood) continue;
		}
		removeUnitigCount += 1;
		for (uint64_t chunk : graph.chunksInUnitig[i])
		{
			removeChunks.insert(chunk & maskUint64_t);
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (removeChunks.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			std::get<2>(chunksPerRead[i][j]) = std::numeric_limits<size_t>::max();
		}
	}
	std::cerr << "removed " << removeUnitigCount << " unitigs, " << removeChunks.size() << " chunks" << std::endl;
}
