#include <queue>
#include <mutex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "UnionFind.h"
#include "CommonUtils.h"
#include "ReadStorage.h"
#include "MatchIndex.h"
#include "KmerGraph.h"
#include "KmerMatcher.h"
#include "MatchGroup.h"
#include "UnitigGraph.h"
#include "GraphCleaner.h"
#include "AlnHaploFilter.h"
#include "GraphPhaser.h"
#include "GraphResolver.h"
#include "AnchorFinder.h"
#include "ChunkmerFilter.h"
#include "edlib.h"

const double mismatchFraction = 0.03; // try 2-3x avg error rate

bool NonexistantChunk(const uint64_t chunk)
{
	if (chunk == std::numeric_limits<uint64_t>::max()) return true;
	if (chunk == std::numeric_limits<uint64_t>::max()-1) return true;
	if (chunk == (std::numeric_limits<uint64_t>::max() ^ firstBitUint64_t)) return true;
	if (chunk == (std::numeric_limits<uint64_t>::max() ^ firstBitUint64_t)-1) return true;
	return false;
}

class ChunkUnitigGraph
{
public:
	std::vector<size_t> unitigLengths;
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
	size_t pathLeftClip;
	size_t pathRightClip;
	size_t readStartPos;
	size_t readEndPos;
};

void writeUnitigPaths(const std::string& filename, const ChunkUnitigGraph& unitigGraph, const std::vector<std::vector<UnitigPath>>& readPaths, const std::vector<std::string>& readNames, const std::vector<size_t>& rawReadLengths)
{
	std::ofstream file { filename };
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			std::string pathstr;
			size_t pathLength = 0;
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				uint64_t node = readPaths[i][j].path[k];
				pathstr += ((node & firstBitUint64_t) ? ">" : "<") + std::to_string(node & maskUint64_t);
				pathLength += unitigGraph.unitigLengths[node & maskUint64_t];
				if (k >= 1)
				{
					uint64_t prevNode = readPaths[i][j].path[k-1];
					auto pairkey = MBG::canon(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t), std::make_pair(node & maskUint64_t, node & firstBitUint64_t));
					std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
					pathLength -= unitigGraph.edgeOverlaps.at(key);
				}
			}
			file << readNames[i] << "\t" << rawReadLengths[i] << "\t" << readPaths[i][j].readStartPos << "\t" << readPaths[i][j].readEndPos << "\t+\t" << pathstr << "\t" << pathLength << "\t" << readPaths[i][j].pathLeftClip << "\t" << (pathLength - readPaths[i][j].pathRightClip) << "\t" << (pathLength - readPaths[i][j].pathLeftClip- readPaths[i][j].pathRightClip) << "\t" << (pathLength - readPaths[i][j].pathLeftClip- readPaths[i][j].pathRightClip) << "\t60" << std::endl;
		}
	}
}

void writeUnitigGraph(const std::string& filename, const ChunkUnitigGraph& unitigGraph)
{
	assert(unitigGraph.unitigLengths.size() == unitigGraph.coverages.size());
	std::ofstream file { filename };
	for (size_t i = 0; i < unitigGraph.unitigLengths.size(); i++)
	{
		file << "S\t" << i << "\t*\tLN:i:" << unitigGraph.unitigLengths[i] << "\tll:f:" << unitigGraph.coverages[i] << "\tFC:f:" << (unitigGraph.unitigLengths[i] * unitigGraph.coverages[i]) << std::endl;
	}
	for (auto pair : unitigGraph.edgeOverlaps)
	{
		assert(unitigGraph.edges.hasEdge(std::make_pair(pair.first.first & maskUint64_t, pair.first.first & firstBitUint64_t), std::make_pair(pair.first.second & maskUint64_t, pair.first.second & firstBitUint64_t)) == 1);
		assert(unitigGraph.edges.hasEdge(std::make_pair(pair.first.second & maskUint64_t, (pair.first.second ^ firstBitUint64_t) & firstBitUint64_t), std::make_pair(pair.first.first & maskUint64_t, (pair.first.first ^ firstBitUint64_t) & firstBitUint64_t)) == 1);
		file << "L\t" << (pair.first.first & maskUint64_t) << "\t" << ((pair.first.first & firstBitUint64_t) ? "+" : "-") << "\t" << (pair.first.second & maskUint64_t) << "\t" << ((pair.first.second & firstBitUint64_t) ? "+" : "-") << "\t" << pair.second << "M\tec:i:" << unitigGraph.edgeCoverages.at(pair.first) << std::endl; 
	}
}

void writeGraph(const std::string& filename, const std::vector<bool>& allowedNode, const SparseEdgeContainer& allowedEdges, const std::vector<std::vector<size_t>>& lengths, const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeCoverage, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps)
{
	std::ofstream graph { filename };
	for (size_t i = 0; i < coverages.size(); i++)
	{
		assert(lengths[i].size() >= coverages[i]);
		if (!allowedNode[i]) continue;
		size_t length = lengths[i][lengths[i].size()/2];
		size_t coverage = coverages[i];
		graph << "S\t" << i << "\t*\tLN:i:" << length << "\tll:i:" << coverage << "\tFC:i:" << length*coverage << std::endl;
	}
	for (auto pair : edgeCoverage)
	{
		if (!allowedEdges.hasEdge(std::make_pair(pair.first.first & maskUint64_t, pair.first.first & firstBitUint64_t), std::make_pair(pair.first.second & maskUint64_t, pair.first.second & firstBitUint64_t))) continue;
		assert(edgeOverlaps.count(pair.first) == 1);
		graph << "L\t" << (pair.first.first & maskUint64_t) << "\t" << ((pair.first.first & firstBitUint64_t) ? "+" : "-") << "\t" << (pair.first.second & maskUint64_t) << "\t" << ((pair.first.second & firstBitUint64_t) ? "+" : "-") << "\t" << edgeOverlaps.at(pair.first) << "M\tec:i:" << pair.second << std::endl;
	}
}

std::pair<std::vector<bool>, SparseEdgeContainer> getAllowedNodesAndEdges(const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeCoverage, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<bool> allowedNode;
	allowedNode.resize(coverages.size(), false);
	SparseEdgeContainer allowedEdges;
	allowedEdges.resize(coverages.size());
	for (size_t i = 0; i < coverages.size(); i++)
	{
		allowedNode[i] = (coverages[i] >= 2);
	}
	for (auto pair : edgeCoverage)
	{
		if (pair.second < 2) continue;
		allowedEdges.addEdge(std::make_pair(pair.first.first & maskUint64_t, pair.first.first & firstBitUint64_t), std::make_pair(pair.first.second & maskUint64_t, pair.first.second & firstBitUint64_t));
		allowedEdges.addEdge(std::make_pair(pair.first.second & maskUint64_t, (pair.first.second ^ firstBitUint64_t) & firstBitUint64_t), std::make_pair(pair.first.first & maskUint64_t, (pair.first.first ^ firstBitUint64_t) & firstBitUint64_t));
	}
	phmap::flat_hash_set<uint64_t> tip;
	for (size_t i = 0; i < coverages.size(); i++)
	{
		if (!allowedNode[i]) continue;
		if (allowedEdges.getEdges(std::make_pair(i, true)).size() == 0) tip.insert(i + firstBitUint64_t);
		if (allowedEdges.getEdges(std::make_pair(i, false)).size() == 0) tip.insert(i);
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, size_t>> tipConnections;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		size_t lastTip = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!allowedNode[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			if (lastTip != std::numeric_limits<size_t>::max() && tip.count(std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t) == 1)
			{
				uint64_t tipHere = std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t;
				tipConnections[lastTip][tipHere] += 1;
				tipConnections[tipHere][lastTip] += 1;
			}
			lastTip = std::numeric_limits<size_t>::max();
			if (tip.count(std::get<2>(chunksPerRead[i][j])) == 1)
			{
				lastTip = std::get<2>(chunksPerRead[i][j]);
			}
		}
	}
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> allowedTips;
	for (uint64_t t : tip)
	{
		assert(!NonexistantChunk(t));
		if (tipConnections.count(t) == 0) continue;
		if (tipConnections.at(t).size() != 1) continue;
		uint64_t otherEnd = tipConnections.at(t).begin()->first;
		assert(!NonexistantChunk(otherEnd));
		assert(tipConnections.count(otherEnd) == 1);
		if (tipConnections.at(otherEnd).size() != 1) continue;
		assert((tipConnections.at(otherEnd).begin()->first) == t);
		allowedTips.emplace(t, otherEnd);
		allowedTips.emplace(otherEnd, t);
//		std::cerr << "allowed tip from " << ((t & firstBitUint64_t) ? ">" : "<") << (t & maskUint64_t) << " to " << ((otherEnd & firstBitUint64_t) ? ">" : "<") << (otherEnd & maskUint64_t) << std::endl;
	}
	phmap::flat_hash_set<size_t> newlyAllowedNodes;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		size_t lastTip = std::numeric_limits<size_t>::max();
		size_t lastTipIndex = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!allowedNode[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			if (lastTip != std::numeric_limits<size_t>::max() && tip.count(std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t) == 1)
			{
				assert(lastTipIndex < j);
				uint64_t tipHere = std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t;
				if (allowedTips.count(std::make_pair(lastTip, tipHere)) == 1)
				{
					for (size_t k = lastTipIndex; k < j; k++)
					{
						if (NonexistantChunk(std::get<2>(chunksPerRead[i][k]))) continue;
						newlyAllowedNodes.insert(std::get<2>(chunksPerRead[i][k]) & maskUint64_t);
						std::pair<size_t, bool> fromnode { std::get<2>(chunksPerRead[i][k]) & maskUint64_t, std::get<2>(chunksPerRead[i][k]) & firstBitUint64_t };
						std::pair<size_t, bool> tonode { std::get<2>(chunksPerRead[i][k+1]) & maskUint64_t, std::get<2>(chunksPerRead[i][k+1]) & firstBitUint64_t };
						allowedEdges.addEdge(fromnode, tonode);
						allowedEdges.addEdge(reverse(tonode), reverse(fromnode));
					}
				}
			}
			lastTip = std::numeric_limits<size_t>::max();
			lastTipIndex = std::numeric_limits<size_t>::max();
			if (tip.count(std::get<2>(chunksPerRead[i][j])) == 1)
			{
				lastTip = std::get<2>(chunksPerRead[i][j]);
				lastTipIndex = j;
			}
		}
	}
	for (auto node : newlyAllowedNodes)
	{
		assert(!NonexistantChunk(node));
		allowedNode[node] = true;
	}
	return std::make_pair(allowedNode, allowedEdges);
}

std::pair<std::vector<std::vector<size_t>>, std::vector<size_t>> getLengthsAndCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if ((std::get<2>(chunksPerRead[i][j]) & maskUint64_t) >= coverages.size())
			{
				coverages.resize((std::get<2>(chunksPerRead[i][j]) & maskUint64_t)+1, 0);
				lengths.resize((std::get<2>(chunksPerRead[i][j]) & maskUint64_t)+1);
			}
			assert(j == 0 || std::get<0>(chunksPerRead[i][j]) > std::get<0>(chunksPerRead[i][j-1]));
			assert(j == 0 || std::get<1>(chunksPerRead[i][j]) > std::get<1>(chunksPerRead[i][j-1]));
			lengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t].emplace_back(std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]));
		}
	}
	for (size_t i = 0; i < lengths.size(); i++)
	{
		std::sort(lengths[i].begin(), lengths[i].end());
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			assert(j == 0 || std::get<0>(chunksPerRead[i][j]) > std::get<0>(chunksPerRead[i][j-1]));
			assert(j == 0 || std::get<1>(chunksPerRead[i][j]) > std::get<1>(chunksPerRead[i][j-1]));
			coverages[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] += 1;
		}
	}
	return std::make_pair(lengths, coverages);
}

phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getEdgeOverlaps(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<size_t>> overlaps;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1]))) continue;
			assert(std::get<0>(chunksPerRead[i][j]) > std::get<0>(chunksPerRead[i][j-1]));
			assert(std::get<1>(chunksPerRead[i][j]) > std::get<1>(chunksPerRead[i][j-1]));
			auto prev = std::get<2>(chunksPerRead[i][j-1]);
			auto curr = std::get<2>(chunksPerRead[i][j]);
			auto pairkey = MBG::canon(std::make_pair(prev & maskUint64_t, prev & firstBitUint64_t), std::make_pair(curr & maskUint64_t, curr & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
			size_t overlap = 0;
			if (std::get<1>(chunksPerRead[i][j-1]) > std::get<0>(chunksPerRead[i][j])) overlap = std::get<1>(chunksPerRead[i][j-1]) - std::get<0>(chunksPerRead[i][j]);
			overlaps[key].push_back(overlap);
		}
	}
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> result;
	for (auto& pair : overlaps)
	{
		std::sort(pair.second.begin(), pair.second.end());
		result[pair.first] = pair.second[pair.second.size()/2];
	}
	return result;
}

phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getEdgeCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1]))) continue;
			assert(std::get<0>(chunksPerRead[i][j]) > std::get<0>(chunksPerRead[i][j-1]));
			assert(std::get<1>(chunksPerRead[i][j]) > std::get<1>(chunksPerRead[i][j-1]));
			auto prev = std::get<2>(chunksPerRead[i][j-1]);
			auto curr = std::get<2>(chunksPerRead[i][j]);
			auto pairkey = MBG::canon(std::make_pair(prev & maskUint64_t, prev & firstBitUint64_t), std::make_pair(curr & maskUint64_t, curr & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
			edgeCoverage[key] += 1;
		}
	}
	return edgeCoverage;
}

void removeContainedChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = chunksPerRead[i].size()-1; j < chunksPerRead[i].size(); j--)
		{
			if (j > 0)
			{
				assert(std::get<1>(chunksPerRead[i][j-1]) <= std::get<1>(chunksPerRead[i][j]));
				assert(std::get<0>(chunksPerRead[i][j-1]) <= std::get<0>(chunksPerRead[i][j]));
				if (std::get<1>(chunksPerRead[i][j-1]) == std::get<1>(chunksPerRead[i][j]))
				{
					assert(std::get<0>(chunksPerRead[i][j-1]) < std::get<0>(chunksPerRead[i][j]));
					chunksPerRead[i].erase(chunksPerRead[i].begin()+j);
					continue;
				}
			}
			if (j+1 < chunksPerRead[i].size())
			{
				assert(std::get<1>(chunksPerRead[i][j]) <= std::get<1>(chunksPerRead[i][j+1]));
				assert(std::get<0>(chunksPerRead[i][j]) <= std::get<0>(chunksPerRead[i][j+1]));
				if (std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j+1]))
				{
					assert(std::get<1>(chunksPerRead[i][j]) < std::get<1>(chunksPerRead[i][j+1]));
					chunksPerRead[i].erase(chunksPerRead[i].begin()+j);
					continue;
				}
			}
		}
	}
}

void splitPerLength(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::cerr << "splitting by length" << std::endl;
	const double differenceFraction = 0.02;
	const size_t differenceConstant = 50;
	phmap::flat_hash_map<size_t, std::vector<size_t>> lengthsPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			lengthsPerChunk[std::get<2>(t) & maskUint64_t].emplace_back(std::get<1>(t) - std::get<0>(t));
		}
	}
	phmap::flat_hash_map<size_t, std::vector<size_t>> splitters;
	for (auto& pair : lengthsPerChunk)
	{
		std::sort(pair.second.begin(), pair.second.end());
		splitters[pair.first].emplace_back(0);
		for (size_t i = 1; i < pair.second.size(); i++)
		{
			size_t distance = pair.second[i] - pair.second[i-1];
			if (distance < std::max((size_t)(pair.second[i] * differenceFraction), (size_t)differenceConstant)) continue;
			splitters[pair.first].emplace_back((pair.second[i] + pair.second[i-1])/2);
		}
	}
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> splitterToNode;
	size_t nextNum = 0;
	for (const auto& pair : splitters)
	{
		for (const size_t len : pair.second)
		{
			auto key = std::make_pair(pair.first, len);
			assert(splitterToNode.count(key) == 0);
			assert(nextNum == splitterToNode.size());
			splitterToNode[key] = nextNum;
			nextNum += 1;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto& t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			assert(std::get<1>(t) > std::get<0>(t));
			size_t distance = std::get<1>(t) - std::get<0>(t);
			size_t dist = std::numeric_limits<size_t>::max();
			for (size_t len : splitters[std::get<2>(t) & maskUint64_t])
			{
				if (len >= distance) break;
				dist = len;
			}
			assert(dist != std::numeric_limits<size_t>::max());
			assert(dist < distance);
			auto key = std::make_pair(std::get<2>(t) & maskUint64_t, dist);
			assert(splitterToNode.count(key) == 1);
			std::get<2>(t) = (std::get<2>(t) & firstBitUint64_t) + splitterToNode.at(key);
		}
	}
}

size_t getNumMismatches(const std::string& leftSequence, const std::string& rightSequence, const size_t maxMismatches)
{
	EdlibAlignResult result = edlibAlign(leftSequence.data(), leftSequence.size(), rightSequence.data(), rightSequence.size(), edlibNewAlignConfig(maxMismatches+1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	size_t mismatches = maxMismatches+1;
	if (result.status == EDLIB_STATUS_OK)
	{
		mismatches = result.editDistance;
	}
	edlibFreeAlignResult(result);
	return mismatches;
}

template <typename F>
void iterateKmers(const TwobitString& baseSequence, const size_t start, const size_t end, const bool fw, const size_t k, F callback)
{
	const size_t mask = (1ull << (2ull*k))-1;
	size_t kmer = 0;
	if (fw)
	{
		for (size_t m = start; m <= end; m++)
		{
			kmer <<= 2;
			kmer += baseSequence.get(m);
			kmer &= mask;
			if (m < start+k) continue;
			callback(kmer, m - start);
		}
	}
	else
	{
		for (size_t m = end; m >= start && m < baseSequence.size(); m--)
		{
			kmer <<= 2;
			kmer += 3 - baseSequence.get(m);
			kmer &= mask;
			if (m+k > end) continue;
			callback(kmer, end - m);
		}
	}
}

void checkPhasablePair(const std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t>>>& bwForks, const std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t>>>& fwForks, const size_t bwHom, const size_t fwHom, std::vector<std::vector<size_t>>& phaseIdentities)
{
	assert(bwForks.size() == fwForks.size());
	std::vector<std::pair<size_t, size_t>> pairs;
	for (size_t i = 0; i < bwForks.size(); i++)
	{
		size_t bw = std::numeric_limits<size_t>::max();
		size_t fw = std::numeric_limits<size_t>::max();
		size_t bwPos = std::numeric_limits<size_t>::max();
		for (auto t : bwForks[i])
		{
			if (std::get<1>(t) != bwHom) continue;
			assert(bw == std::numeric_limits<size_t>::max());
			bw = std::get<2>(t);
			bwPos = std::get<0>(t);
		}
		if (bw == std::numeric_limits<size_t>::max()) return;
		for (auto t : fwForks[i])
		{
			if (std::get<1>(t) != fwHom) continue;
			assert(fw == std::numeric_limits<size_t>::max());
			if (std::get<0>(t) > bwPos) return;
			fw = std::get<2>(t);
		}
		if (fw == std::numeric_limits<size_t>::max()) return;
		pairs.emplace_back(bw, fw);
	}
	assert(pairs.size() == bwForks.size());
	std::sort(pairs.begin(), pairs.end());
	phmap::flat_hash_map<size_t, size_t> parent;
	for (auto pair : pairs)
	{
		assert((pair.first & firstBitUint64_t) == 0);
		assert((pair.second & firstBitUint64_t) == 0);
		if (parent.count(pair.first) == 0) parent[pair.first] = pair.first;
		if (parent.count(pair.second + firstBitUint64_t) == 0) parent[pair.second + firstBitUint64_t] = pair.second + firstBitUint64_t;
		merge(parent, pair.first, pair.second + firstBitUint64_t);
	}
	phmap::flat_hash_map<size_t, size_t> clusterCount;
	for (auto pair : pairs)
	{
		auto cluster = find(parent, pair.first);
		clusterCount[cluster] += 1;
	}
	if (clusterCount.size() < 2) return;
	for (auto pair : clusterCount)
	{
		if (pair.second < 4) return;
	}
	for (size_t i = 0; i < bwForks.size(); i++)
	{
		bool found = false;
		for (auto t : bwForks[i])
		{
			if (std::get<1>(t) != bwHom) continue;
			assert(!found);
			phaseIdentities[std::get<3>(t)].emplace_back(find(parent, std::get<2>(t)));
			found = true;
		}
		assert(found);
	}
}

void splitPerPhasingKmersWithinChunk(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by phasing kmers" << std::endl;
	const size_t k = 11;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, k](const size_t i)
	{
		phmap::flat_hash_set<size_t> kmersEverywhere;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			phmap::flat_hash_set<size_t> kmersHere;
			phmap::flat_hash_set<size_t> kmersRepeatingHere;
			iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&occurrencesPerChunk, &kmersHere, &kmersRepeatingHere, j, i](const size_t kmer, const size_t pos)
			{
				if (kmersHere.count(kmer) == 0)
				{
					kmersHere.insert(kmer);
					return;
				}
				else
				{
					kmersRepeatingHere.insert(kmer);
				}
			});
			if (j == 0)
			{
				for (size_t kmer : kmersHere)
				{
					if (kmersRepeatingHere.count(kmer) == 1) continue;
					kmersEverywhere.insert(kmer);
				}
			}
			else
			{
				phmap::flat_hash_set<size_t> removeKmers;
				for (size_t kmer : kmersEverywhere)
				{
					if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
					if (kmersRepeatingHere.count(kmer) == 1) removeKmers.insert(kmer);
				}
				for (size_t kmer : removeKmers) kmersEverywhere.erase(kmer);
			}
		}
		std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t>>> fwForks; // hom pos, hom kmer, het kmer, occurrence ID
		std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t>>> bwForks; // hom pos, hom kmer, het kmer, occurrence ID
		fwForks.resize(occurrencesPerChunk[i].size());
		bwForks.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			size_t lastKmer = std::numeric_limits<size_t>::max();
			iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [i, j, &fwForks, &bwForks, &lastKmer, &kmersEverywhere](const size_t kmer, const size_t pos)
			{
				if (kmersEverywhere.count(lastKmer) == 1 && kmersEverywhere.count(kmer) == 0)
				{
					fwForks[j].emplace_back(pos-1, lastKmer, kmer, j);
				}
				else if (lastKmer != std::numeric_limits<size_t>::max() && kmersEverywhere.count(lastKmer) == 0 && kmersEverywhere.count(kmer) == 1)
				{
					bwForks[j].emplace_back(pos, kmer, lastKmer, j);
				}
				lastKmer = kmer;
			});
		}
		phmap::flat_hash_set<std::pair<size_t, size_t>> checkedPairs;
		std::vector<std::vector<size_t>> phaseIdentities;
		phaseIdentities.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			for (size_t m = 0; m < fwForks[j].size(); m++)
			{
				for (size_t n = 0; n < bwForks[j].size(); n++)
				{
					if (std::get<0>(bwForks[j][n]) > std::get<0>(fwForks[j][m])) break;
					std::pair<size_t, size_t> key { std::get<1>(bwForks[j][n]), std::get<1>(fwForks[j][m]) };
					if (checkedPairs.count(key) == 1) continue;
					checkPhasablePair(bwForks, fwForks, key.first, key.second, phaseIdentities);
					checkedPairs.insert(key);
				}
			}
		}
		std::vector<size_t> sortedPhaseIdentityIndices;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sortedPhaseIdentityIndices.emplace_back(j);
		}
		std::sort(sortedPhaseIdentityIndices.begin(), sortedPhaseIdentityIndices.end(), [&phaseIdentities](size_t left, size_t right) { return phaseIdentities[left] < phaseIdentities[right]; });
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < sortedPhaseIdentityIndices.size(); j++)
		{
			if (phaseIdentities[sortedPhaseIdentityIndices[j-1]] != phaseIdentities[sortedPhaseIdentityIndices[j]]) continue;
			merge(parent, sortedPhaseIdentityIndices[j-1], sortedPhaseIdentityIndices[j]);
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "phasable kmer splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
}

bool sequenceAndDistanceApproxMatch(const std::vector<size_t>& leftKmers, const std::vector<size_t>& rightKmers, const std::vector<size_t>& leftDistances, const std::vector<size_t>& rightDistances, const size_t maxDistanceDifference)
{
	assert(leftKmers.size() == rightKmers.size());
	assert(leftDistances.size() == rightDistances.size());
	assert(leftDistances.size() == leftKmers.size()+1);
	for (size_t i = 0; i < leftKmers.size(); i++)
	{
		if (leftKmers[i] != rightKmers[i]) return false;
	}
	for (size_t i = 0; i < leftDistances.size(); i++)
	{
		if (leftDistances[i] > rightDistances[i]+maxDistanceDifference) return false;
		if (rightDistances[i] > leftDistances[i]+maxDistanceDifference) return false;
	}
	return true;
}

void splitPerAllUniqueKmerSVs(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by all-unique kmers" << std::endl;
	const size_t maxDistanceDifference = 50;
	const size_t k = 11;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, k](const size_t i)
	{
		phmap::flat_hash_set<size_t> kmersEverywhere;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			phmap::flat_hash_set<size_t> kmersHere;
			phmap::flat_hash_set<size_t> kmersRepeatingHere;
			iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&occurrencesPerChunk, &kmersHere, &kmersRepeatingHere, j, i](const size_t kmer, const size_t pos)
			{
				if (kmersHere.count(kmer) == 0)
				{
					kmersHere.insert(kmer);
					return;
				}
				else
				{
					kmersRepeatingHere.insert(kmer);
				}
			});
			if (j == 0)
			{
				for (size_t kmer : kmersHere)
				{
					if (kmersRepeatingHere.count(kmer) == 1) continue;
					kmersEverywhere.insert(kmer);
				}
			}
			else
			{
				phmap::flat_hash_set<size_t> removeKmers;
				for (size_t kmer : kmersEverywhere)
				{
					if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
					if (kmersRepeatingHere.count(kmer) == 1) removeKmers.insert(kmer);
				}
				for (size_t kmer : removeKmers) kmersEverywhere.erase(kmer);
			}
		}
		std::vector<std::vector<size_t>> kmerSequence;
		std::vector<std::vector<size_t>> kmerDistances;
		kmerSequence.resize(occurrencesPerChunk[i].size());
		kmerDistances.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			size_t lastPos = 0;
			iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&kmersEverywhere, &kmerSequence, &kmerDistances, &lastPos, j, i](const size_t kmer, const size_t pos)
			{
				if (kmersEverywhere.count(kmer) == 0) return;
				kmerSequence[j].emplace_back(kmer);
				kmerDistances[j].emplace_back(pos - lastPos);
				lastPos = pos;
			});
			kmerDistances[j].emplace_back(std::get<1>(t)+1 - lastPos);
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				if (find(parent, j) == find(parent, k)) continue;
				if (!sequenceAndDistanceApproxMatch(kmerSequence[j], kmerSequence[k], kmerDistances[j], kmerDistances[j], maxDistanceDifference)) continue;
				merge(parent, j, k);
			}
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "all-unique kmer splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
}

void countGoodishKmersInChunks(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "counting goodish kmers" << std::endl;
	const size_t k = 11;
	const size_t distance = 100;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, k, distance](const size_t i)
	{
		phmap::flat_hash_set<size_t> kmersEverywhere;
		size_t averageSize = 0;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(std::get<1>(t) > std::get<0>(t));
			averageSize += std::get<1>(t) - std::get<0>(t);
		}
		averageSize /= occurrencesPerChunk[i].size();
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			phmap::flat_hash_map<size_t, size_t> kmerLastOccurrence;
			phmap::flat_hash_set<size_t> kmersRepeatingHere;
			iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&occurrencesPerChunk, &kmerLastOccurrence, &kmersRepeatingHere, j, i, distance](const size_t kmer, const size_t pos)
			{
				if (kmerLastOccurrence.count(kmer) == 0)
				{
					kmerLastOccurrence[kmer] = pos;
					return;
				}
				else if (kmerLastOccurrence.at(kmer) + distance >= pos)
				{
					kmersRepeatingHere.insert(kmer);
				}
				kmerLastOccurrence[kmer] = pos;
			});
			if (j == 0)
			{
				for (auto pair : kmerLastOccurrence)
				{
					if (kmersRepeatingHere.count(pair.first) == 1) continue;
					kmersEverywhere.insert(pair.first);
				}
			}
			else
			{
				phmap::flat_hash_set<size_t> removeKmers;
				for (size_t kmer : kmersEverywhere)
				{
					if (kmerLastOccurrence.count(kmer) == 0) removeKmers.insert(kmer);
					if (kmersRepeatingHere.count(kmer) == 1) removeKmers.insert(kmer);
				}
				for (size_t kmer : removeKmers) kmersEverywhere.erase(kmer);
			}
		}
		phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, size_t>>> maybeGoodishKmerPositions;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&kmersEverywhere, &maybeGoodishKmerPositions, averageSize, j, i, distance, t](const size_t kmer, const size_t pos)
			{
				if (kmersEverywhere.count(kmer) == 0) return;
				size_t extrapolatedPosition = (double)pos / (double)(std::get<1>(t) - std::get<0>(t)) * (double)averageSize;
				maybeGoodishKmerPositions[kmer].emplace_back(extrapolatedPosition, j);
			});
		}
		size_t countGoodishClusters = 0;
		size_t countGoodishKmers = 0;
		for (auto& pair : maybeGoodishKmerPositions)
		{
			std::sort(pair.second.begin(), pair.second.end());
			phmap::flat_hash_set<size_t> readsSinceLastBreak;
			assert(pair.second.size() >= 1);
			readsSinceLastBreak.insert(pair.second[0].second);
			bool hasAny = false;
			bool thisValid = true;
			for (size_t j = 1; j < pair.second.size(); j++)
			{
				if (pair.second[j].first < pair.second[j-1].first + distance)
				{
					if (readsSinceLastBreak.count(pair.second[j].second) == 1)
					{
						thisValid = false;
					}
					readsSinceLastBreak.insert(pair.second[j].second);
					continue;
				}
				if (thisValid && readsSinceLastBreak.size() == occurrencesPerChunk[i].size())
				{
					hasAny = true;
					countGoodishClusters += 1;
				}
				thisValid = true;
				readsSinceLastBreak.clear();
				readsSinceLastBreak.insert(pair.second[j].second);
			}
			if (thisValid && readsSinceLastBreak.size() == occurrencesPerChunk[i].size())
			{
				hasAny = true;
				countGoodishClusters += 1;
			}
			if (hasAny) countGoodishKmers += 1;
		}
		std::cerr << "chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " has " << countGoodishClusters << " goodish clusters, " << countGoodishKmers << " goodish kmers (" << kmersEverywhere.size() << " pre-filter kmers)" << std::endl;
	});
}

void countGoodKmersInChunks(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "counting good kmers" << std::endl;
	const size_t k = 11;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, k](const size_t i)
	{
		phmap::flat_hash_set<size_t> kmersEverywhere;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			phmap::flat_hash_set<size_t> kmersHere;
			phmap::flat_hash_set<size_t> kmersRepeatingHere;
			iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&occurrencesPerChunk, &kmersHere, &kmersRepeatingHere, j, i](const size_t kmer, const size_t pos)
			{
				if (kmersHere.count(kmer) == 0)
				{
					kmersHere.insert(kmer);
					return;
				}
				else
				{
					kmersRepeatingHere.insert(kmer);
				}
			});
			if (j == 0)
			{
				for (size_t kmer : kmersHere)
				{
					if (kmersRepeatingHere.count(kmer) == 1) continue;
					kmersEverywhere.insert(kmer);
				}
			}
			else
			{
				phmap::flat_hash_set<size_t> removeKmers;
				for (size_t kmer : kmersEverywhere)
				{
					if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
					if (kmersRepeatingHere.count(kmer) == 1) removeKmers.insert(kmer);
				}
				for (size_t kmer : removeKmers) kmersEverywhere.erase(kmer);
			}
		}
		std::cerr << "chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " has " << kmersEverywhere.size() << " good kmers" << std::endl;

	});
}

void splitPerSequenceIdentity(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by sequence identity" << std::endl;
	const size_t mismatchFloor = 10;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	std::vector<size_t> iterationOrder;
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		iterationOrder.emplace_back(i);
	}
	std::sort(iterationOrder.begin(), iterationOrder.end(), [&occurrencesPerChunk](size_t left, size_t right) { return occurrencesPerChunk[left].size() > occurrencesPerChunk[right].size(); });
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, &iterationOrder, mismatchFloor](const size_t iterationIndex)
	{
		const size_t i = iterationOrder[iterationIndex];
//		{
//			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "sequence identity split chunk " << i << " coverage " << occurrencesPerChunk[i].size() << std::endl;
//		}
		std::vector<std::pair<std::string, size_t>> sequencesPerOccurrence;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sequencesPerOccurrence.emplace_back();
			sequencesPerOccurrence.back().second = j;
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			if (std::get<2>(t) & firstBitUint64_t)
			{
				for (size_t k = std::get<0>(t); k <= std::get<1>(t); k++)
				{
					sequencesPerOccurrence.back().first.push_back("ACGT"[readSequences[occurrencesPerChunk[i][j].first].get(k)]);
				}
			}
			else
			{
				for (size_t k = std::get<1>(t); k >= std::get<0>(t) && k < readSequences[occurrencesPerChunk[i][j].first].size(); k--)
				{
					sequencesPerOccurrence.back().first.push_back("ACGT"[3-readSequences[occurrencesPerChunk[i][j].first].get(k)]);
				}
			}
		}
		std::sort(sequencesPerOccurrence.begin(), sequencesPerOccurrence.end(), [](const auto& left, const auto& right) {
			if (left.first.size() < right.first.size()) return true;
			if (left.first.size() > right.first.size()) return false;
			return left.first < right.first;
		});
		std::vector<bool> differentFromPredecessor;
		differentFromPredecessor.resize(occurrencesPerChunk[i].size(), false);
		differentFromPredecessor[0] = true;
		for (size_t j = 1; j < differentFromPredecessor.size(); j++)
		{
			if (sequencesPerOccurrence[j].first == sequencesPerOccurrence[j-1].first) continue;
			differentFromPredecessor[j] = true;
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
		{
			if (!differentFromPredecessor[j])
			{
				merge(parent, j, j-1);
				continue;
			}
			for (size_t k = j-1; k < occurrencesPerChunk[i].size(); k--)
			{
				if (!differentFromPredecessor[k]) continue;
				if (find(parent, k) == find(parent, j)) continue;
				assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
				size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
				if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches) break;
				size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
				if (mismatches > maxMismatches) continue;
				merge(parent, j, k);
			}
		}
		phmap::flat_hash_map<size_t, size_t> keyToNode;
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "sequence identity splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t index = sequencesPerOccurrence[j].second;
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
}

// https://naml.us/post/inverse-of-a-hash-function/
uint64_t hash(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

std::vector<size_t> getMinHashes(const std::string& sequence, const size_t k, const size_t count)
{
	assert(sequence.size() > k + count);
	uint64_t kmer = 0;
	uint64_t mask = (1ull << 2ull*k) - 1;
	for (size_t i = 0; i < k; i++)
	{
		kmer <<= 2;
		switch(sequence[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
		}
	}
	std::priority_queue<size_t> queue;
	queue.emplace(hash(kmer));
	for (size_t i = k; i < sequence.size(); i++)
	{
		kmer <<= 2;
		switch(sequence[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
		}
		kmer &= mask;
		size_t h = hash(kmer);
		if (queue.size() < count)
		{
			queue.emplace(h);
		}
		else if (h < queue.top())
		{
			queue.pop();
			queue.emplace(h);
		}
	}
	std::vector<size_t> result;
	while (queue.size() > 0)
	{
		result.emplace_back(queue.top());
		queue.pop();
	}
	return result;
}

void splitPerBaseCounts(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t mismatchFloor = 10;
	std::cerr << "splitting by base counts" << std::endl;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t) + 1);
			occurrencesPerChunk[std::get<2>(t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, mismatchFloor](const size_t i)
	{
		std::vector<std::vector<size_t>> countsPerOccurrence;
		countsPerOccurrence.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			countsPerOccurrence[j].resize(4);
			for (size_t k = std::get<0>(t); k < std::get<1>(t); k++)
			{
				countsPerOccurrence[j][readSequences[occurrencesPerChunk[i][j].first].get(k)] += 1;
			}
			if ((std::get<2>(t) & firstBitUint64_t) ^ firstBitUint64_t)
			{
				std::swap(countsPerOccurrence[j][0], countsPerOccurrence[j][3]);
				std::swap(countsPerOccurrence[j][1], countsPerOccurrence[j][2]);
			}
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				size_t maxEdits = std::min(std::get<1>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) - std::get<0>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]), std::get<1>(chunksPerRead[occurrencesPerChunk[i][k].first][occurrencesPerChunk[i][k].second]) - std::get<0>(chunksPerRead[occurrencesPerChunk[i][k].first][occurrencesPerChunk[i][k].second]));
				maxEdits *= mismatchFraction;
				maxEdits = std::max(maxEdits, mismatchFloor);
				size_t edits = 0;
				for (size_t c = 0; c < 4; c++)
				{
					if (countsPerOccurrence[j][c] > countsPerOccurrence[k][c])
					{
						edits += countsPerOccurrence[j][c] - countsPerOccurrence[k][c];
					}
					else
					{
						edits += countsPerOccurrence[k][c] - countsPerOccurrence[j][c];
					}
				}
				if (edits < maxEdits) merge(parent, j, k);
			}
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "base count splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
}

void splitPerInterchunkPhasedKmers(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by interchunk phasing kmers" << std::endl;
	const size_t k = 11;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t) + 1);
			occurrencesPerChunk[std::get<2>(t)].emplace_back(i, j);
		}
	}
	std::vector<bool> repetitive;
	repetitive.resize(occurrencesPerChunk.size(), false);
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		phmap::flat_hash_set<size_t> readsHere;
		size_t repeats = 0;
		for (auto t : occurrencesPerChunk[i])
		{
			if (readsHere.count(t.first) == 1) repeats += 1;
			readsHere.insert(t.first);
		}
		if (repeats >= 2) repetitive[i] = true;
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint64_t>> edges;
	{
		phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, size_t>> edgeCoverage;
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 1; j < chunksPerRead[i].size(); j++)
			{
				uint64_t prevNode = std::get<2>(chunksPerRead[i][j-1]);
				uint64_t thisNode = std::get<2>(chunksPerRead[i][j]);
				if (NonexistantChunk(prevNode)) continue;
				if (NonexistantChunk(thisNode)) continue;
				if (repetitive[prevNode & maskUint64_t]) continue;
				if (repetitive[thisNode & maskUint64_t]) continue;
				if (occurrencesPerChunk[prevNode & maskUint64_t].size() == 1) continue;
				if (occurrencesPerChunk[thisNode & maskUint64_t].size() == 1) continue;
				edgeCoverage[prevNode][thisNode] += 1;
				edgeCoverage[thisNode ^ firstBitUint64_t][prevNode ^ firstBitUint64_t] += 1;
			}
		}
		for (const auto& pair : edgeCoverage)
		{
			for (const auto pair2 : pair.second)
			{
				if (pair2.second == 1) continue;
				edges[pair.first].emplace(pair2.first);
			}
		}
	}
	std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> maybePhaseGroups;
	{
		phmap::flat_hash_map<uint64_t, std::tuple<size_t, uint64_t, uint64_t>> edgeToGroupMapping;
		size_t nextPhaseGroup = 0;
		for (const auto& pair : edges)
		{
			if (pair.second.size() != 2) continue;
			assert(!NonexistantChunk(pair.first));
			assert(!repetitive[pair.first & maskUint64_t]);
			uint64_t first = *pair.second.begin();
			uint64_t second = *(++pair.second.begin());
			assert(!NonexistantChunk(first));
			assert(!NonexistantChunk(second));
			assert(!repetitive[first & maskUint64_t]);
			assert(!repetitive[second & maskUint64_t]);
			assert(edges.count(first ^ firstBitUint64_t) == 1);
			assert(edges.count(second ^ firstBitUint64_t) == 1);
			if (edges.at(first ^ firstBitUint64_t).size() != 1) continue;
			if (edges.at(second ^ firstBitUint64_t).size() != 1) continue;
			assert(first != second);
			if ((first & maskUint64_t) == (second & maskUint64_t)) continue;
			assert(*edges.at(first ^ firstBitUint64_t).begin() == (pair.first ^ firstBitUint64_t));
			assert(*edges.at(second ^ firstBitUint64_t).begin() == (pair.first ^ firstBitUint64_t));
			edgeToGroupMapping[pair.first] = std::make_tuple(nextPhaseGroup, first, second);
			nextPhaseGroup += 1;
			maybePhaseGroups.emplace_back();
		}
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 1; j < chunksPerRead[i].size(); j++)
			{
				uint64_t prevNode = std::get<2>(chunksPerRead[i][j-1]);
				uint64_t thisNode = std::get<2>(chunksPerRead[i][j]);
				if (NonexistantChunk(prevNode)) continue;
				if (NonexistantChunk(thisNode)) continue;
				uint64_t allele = std::numeric_limits<size_t>::max();
				uint64_t fork = std::numeric_limits<size_t>::max();
				if (edgeToGroupMapping.count(prevNode) == 1 && edgeToGroupMapping.count(thisNode ^ firstBitUint64_t) == 1) continue;
				if (edgeToGroupMapping.count(prevNode) == 1)
				{
					assert(edgeToGroupMapping.count(thisNode ^ firstBitUint64_t) == 0);
					fork = prevNode;
					allele = thisNode;
				}
				if (edgeToGroupMapping.count(thisNode ^ firstBitUint64_t) == 1)
				{
					assert(edgeToGroupMapping.count(prevNode) == 0);
					fork = (thisNode ^ firstBitUint64_t);
					allele = (prevNode ^ firstBitUint64_t);
				}
				if (fork == std::numeric_limits<size_t>::max()) continue;
				assert(allele != std::numeric_limits<size_t>::max());
				assert(allele != std::get<1>(edgeToGroupMapping.at(fork)) || allele != std::get<2>(edgeToGroupMapping.at(fork)));
				if (allele == std::get<1>(edgeToGroupMapping.at(fork)))
				{
					maybePhaseGroups[std::get<0>(edgeToGroupMapping.at(fork))].first.emplace(i);
				}
				else if (allele == std::get<2>(edgeToGroupMapping.at(fork)))
				{
					maybePhaseGroups[std::get<0>(edgeToGroupMapping.at(fork))].second.emplace(i);
				}
			}
		}
		for (size_t i = maybePhaseGroups.size()-1; i < maybePhaseGroups.size(); i--)
		{
			phmap::flat_hash_set<size_t> inconsistents;
			for (size_t read : maybePhaseGroups[i].first)
			{
				if (maybePhaseGroups[i].second.count(read) == 1) inconsistents.insert(read);
			}
			if (inconsistents.size() >= 3)
			{
				std::swap(maybePhaseGroups[i], maybePhaseGroups.back());
				maybePhaseGroups.pop_back();
				continue;
			}
			for (size_t read : inconsistents)
			{
				maybePhaseGroups[i].first.erase(read);
				maybePhaseGroups[i].second.erase(read);
			}
		}
		for (size_t i = maybePhaseGroups.size()-1; i < maybePhaseGroups.size(); i--)
		{
			if (maybePhaseGroups[i].first.size() >= 4 && maybePhaseGroups[i].second.size() >= 4) continue;
			std::swap(maybePhaseGroups.back(), maybePhaseGroups[i]);
			maybePhaseGroups.pop_back();
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, &repetitive, &maybePhaseGroups, k](const size_t i)
	{
		std::vector<size_t> applicablePhasingGroups;
		std::vector<size_t> readsHere;
		{
			phmap::flat_hash_set<size_t> readsHereUnique;
			for (auto t : occurrencesPerChunk[i])
			{
				readsHereUnique.emplace(t.first);
			}
			readsHere.insert(readsHere.end(), readsHereUnique.begin(), readsHereUnique.end());
		}
		for (size_t j = 0; j < maybePhaseGroups.size(); j++)
		{
			size_t firsts = 0;
			size_t seconds = 0;
			for (size_t read : readsHere)
			{
				if (maybePhaseGroups[j].first.count(read) == 1) firsts += 1;
				if (maybePhaseGroups[j].second.count(read) == 1) seconds += 1;
			}
			if (firsts < 3) continue;
			if (seconds < 3) continue;
			if (firsts+seconds < readsHere.size() * 0.75) continue;
			applicablePhasingGroups.push_back(j);
		}
		phmap::flat_hash_map<size_t, std::vector<size_t>> phaseGroupsPerRead;
		size_t nextGroupNum = 0;
		for (const size_t phaseGroup : applicablePhasingGroups)
		{
			phmap::flat_hash_set<size_t> kmersInFirst;
			phmap::flat_hash_set<size_t> kmersInSecond;
			phmap::flat_hash_set<size_t> kmersInEvenOneFirst;
			phmap::flat_hash_set<size_t> kmersInEvenOneSecond;
			bool hasFirst = false;
			bool hasSecond = false;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 0 && maybePhaseGroups[phaseGroup].second.count(read) == 0) continue;
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				phmap::flat_hash_set<size_t> kmersHere;
				iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&kmersHere](const size_t kmer, const size_t pos)
				{
					kmersHere.insert(kmer);
				});
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1)
				{
					kmersInEvenOneFirst.insert(kmersHere.begin(), kmersHere.end());
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
					if (!hasFirst)
					{
						hasFirst = true;
						kmersInFirst = kmersHere;
					}
					else
					{
						phmap::flat_hash_set<size_t> removeKmers;
						for (size_t kmer : kmersInFirst)
						{
							if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
						}
						for (size_t kmer : removeKmers)
						{
							kmersInFirst.erase(kmer);
						}
					}
				}
				else
				{
					kmersInEvenOneSecond.insert(kmersHere.begin(), kmersHere.end());
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 1);
					if (!hasSecond)
					{
						hasSecond = true;
						kmersInSecond = kmersHere;
					}
					else
					{
						phmap::flat_hash_set<size_t> removeKmers;
						for (size_t kmer : kmersInSecond)
						{
							if (kmersHere.count(kmer) == 0) removeKmers.insert(kmer);
						}
						for (size_t kmer : removeKmers)
						{
							kmersInSecond.erase(kmer);
						}
					}
				}
			}
			{
				for (size_t kmer : kmersInEvenOneFirst)
				{
					if (kmersInSecond.count(kmer) == 0) continue;
					kmersInSecond.erase(kmer);
				}
				for (size_t kmer : kmersInEvenOneSecond)
				{
					if (kmersInFirst.count(kmer) == 0) continue;
					kmersInFirst.erase(kmer);
				}
			}
			{
				phmap::flat_hash_set<size_t> sharedKmers;
				for (auto kmer : kmersInFirst)
				{
					if (kmersInSecond.count(kmer) == 0) continue;
					sharedKmers.insert(kmer);
				}
				for (auto kmer : sharedKmers)
				{
					kmersInFirst.erase(kmer);
					kmersInSecond.erase(kmer);
				}
			}
			if (kmersInFirst.size() == 0 || kmersInSecond.size() == 0) continue;
			bool possiblyValid = true;
			phmap::flat_hash_set<size_t> removeKmers;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1 || maybePhaseGroups[phaseGroup].second.count(read) == 1) continue;
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				phmap::flat_hash_set<size_t> kmersHere;
				iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&kmersHere](const size_t kmer, const size_t pos)
				{
					kmersHere.insert(kmer);
				});
				size_t firstMatches = 0;
				size_t secondMatches = 0;
				for (size_t kmer : kmersHere)
				{
					if (kmersInFirst.count(kmer) == 1)
					{
						assert(kmersInSecond.count(kmer) == 0);
						firstMatches += 1;
					}
					if (kmersInSecond.count(kmer) == 1)
					{
						assert(kmersInFirst.count(kmer) == 0);
						secondMatches += 1;
					}
				}
				if (firstMatches == secondMatches)
				{
					possiblyValid = false;
					break;
				}
				if (firstMatches > secondMatches)
				{
					for (size_t kmer : kmersInFirst)
					{
						if (kmersHere.count(kmer) == 1) continue;
						removeKmers.insert(kmer);
					}
					for (size_t kmer : kmersInSecond)
					{
						if (kmersHere.count(kmer) == 0) continue;
						removeKmers.insert(kmer);
					}
				}
				else
				{
					assert(firstMatches < secondMatches);
					for (size_t kmer : kmersInSecond)
					{
						if (kmersHere.count(kmer) == 1) continue;
						removeKmers.insert(kmer);
					}
					for (size_t kmer : kmersInFirst)
					{
						if (kmersHere.count(kmer) == 0) continue;
						removeKmers.insert(kmer);
					}
				}
			}
			if (!possiblyValid) continue;
			for (auto kmer : removeKmers)
			{
				assert(kmersInFirst.count(kmer) == 1 || kmersInSecond.count(kmer) == 1);
				if (kmersInFirst.count(kmer) == 1)
				{
					assert(kmersInSecond.count(kmer) == 0);
					kmersInFirst.erase(kmer);
				}
				else
				{
					assert(kmersInSecond.count(kmer) == 1);
					kmersInSecond.erase(kmer);
				}
			}
			phmap::flat_hash_map<size_t, size_t> assignments;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				size_t read = occurrencesPerChunk[i][j].first;
				if (maybePhaseGroups[phaseGroup].first.count(read) == 1)
				{
					assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
					assignments[read] = 0;
					continue;
				}
				if (maybePhaseGroups[phaseGroup].second.count(read) == 1)
				{
					assignments[read] = 1;
					continue;
				}
				assert(maybePhaseGroups[phaseGroup].first.count(read) == 0);
				assert(maybePhaseGroups[phaseGroup].second.count(read) == 0);
				auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
				assert(!NonexistantChunk(std::get<2>(t)));
				phmap::flat_hash_set<size_t> kmersHere;
				iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & firstBitUint64_t, k, [&kmersHere](const size_t kmer, const size_t pos)
				{
					kmersHere.insert(kmer);
				});
				size_t firstMatches = 0;
				size_t secondMatches = 0;
				for (size_t kmer : kmersHere)
				{
					if (kmersInFirst.count(kmer) == 1)
					{
						assert(kmersInSecond.count(kmer) == 0);
						firstMatches += 1;
					}
					if (kmersInSecond.count(kmer) == 1)
					{
						assert(kmersInFirst.count(kmer) == 0);
						secondMatches += 1;
					}
				}
				if (firstMatches == secondMatches)
				{
					possiblyValid = false;
					break;
				}
				if (firstMatches > secondMatches)
				{
					if (assignments.count(read) == 1 && assignments.at(read) != 0)
					{
						possiblyValid = false;
						break;
					}
					assert(assignments.count(read) == 0 || assignments.at(read) == 0);
					assignments[read] = 0;
				}
				else
				{
					if (assignments.count(read) == 1 && assignments.at(read) != 1)
					{
						possiblyValid = false;
						break;
					}
					assert(assignments.count(read) == 0 || assignments.at(read) == 1);
					assignments[read] = 1;
				}
			}
			if (!possiblyValid) continue;
			for (auto pair : assignments)
			{
				phaseGroupsPerRead[pair.first].push_back(nextGroupNum + pair.second);
			}
			nextGroupNum += 2;
		}
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		if (phaseGroupsPerRead.size() == 0)
		{
			for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
			{
				merge(parent, 0, j);
			}
		}
		else
		{
			for (size_t j = 1; j < occurrencesPerChunk[i].size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					size_t readj = occurrencesPerChunk[i][j].first;
					size_t readk = occurrencesPerChunk[i][k].first;
					if (phaseGroupsPerRead.at(readj) != phaseGroupsPerRead.at(readk)) continue;
					merge(parent, j, k);
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> clusterToNode;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				if (clusterToNode.count(find(parent, j)) == 1) continue;
				clusterToNode[find(parent, j)] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "interchunk phasing kmers splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, j)) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
}

void splitPerMinHashes(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by minhash" << std::endl;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t) + 1);
			occurrencesPerChunk[std::get<2>(t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences](const size_t i)
	{
		phmap::flat_hash_map<size_t, size_t> parent;
		std::vector<size_t> oneHashPerLocation;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			std::vector<uint64_t> minHashes;
			assert(std::get<0>(t) < std::get<1>(t));
			assert(std::get<1>(t) < readSequences[occurrencesPerChunk[i][j].first].size());
			if (std::get<2>(t) & firstBitUint64_t)
			{
				minHashes = getMinHashes(readSequences[occurrencesPerChunk[i][j].first].substr(std::get<0>(t), std::get<1>(t)-std::get<0>(t)+1), 11, 10);
			}
			else
			{
				minHashes = getMinHashes(MBG::revCompRaw(readSequences[occurrencesPerChunk[i][j].first].substr(std::get<0>(t), std::get<1>(t)-std::get<0>(t)+1)), 11, 10);
			}
			assert(minHashes.size() >= 1);
			for (auto hash : minHashes)
			{
				if (parent.count(hash) == 0) parent[hash] = hash;
			}
			for (auto hash : minHashes)
			{
				merge(parent, hash, minHashes[0]);
			}
			oneHashPerLocation.emplace_back(minHashes[0]);
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> clusterToNode;
			for (auto hash : oneHashPerLocation)
			{
				if (clusterToNode.count(find(parent, hash)) == 1) continue;
				clusterToNode[find(parent, hash)] = nextNum;
				nextNum += 1;
			}
//			std::cerr << "minhash splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, oneHashPerLocation[j])) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
}

void removeSingleCopyChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	phmap::flat_hash_map<size_t, size_t> coverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			coverage[std::get<2>(t) & maskUint64_t] += 1;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = chunksPerRead[i].size()-1; j < chunksPerRead[i].size(); j--)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (coverage.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) > 1) continue;
			chunksPerRead[i].erase(chunksPerRead[i].begin()+j);
		}
	}
}

void removeHighCoverageChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t maxCoverage)
{
	phmap::flat_hash_map<size_t, size_t> coverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			coverage[std::get<2>(t) & maskUint64_t] += 1;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = chunksPerRead[i].size()-1; j < chunksPerRead[i].size(); j--)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (coverage.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) < maxCoverage) continue;
			chunksPerRead[i].erase(chunksPerRead[i].begin()+j);
		}
	}
}

void writeGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::cerr << "writing graph " << graphFile << std::endl;
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
	std::vector<bool> allowedNode;
	SparseEdgeContainer allowedEdges;
	std::tie(allowedNode, allowedEdges) = getAllowedNodesAndEdges(coverages, edgeCoverage, chunksPerRead);
	writeGraph(graphFile, allowedNode, allowedEdges, lengths, coverages, edgeCoverage, edgeOverlaps);
	std::cerr << "writing paths " << pathsFile << std::endl;
	std::ofstream pathfile { pathsFile };
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			uint64_t rawnode = std::get<2>(t);
			pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
		}
	}
}

std::vector<uint64_t> getUnitig(const uint64_t startNode, const std::vector<bool>& allowedNode, const SparseEdgeContainer& allowedEdges)
{
	assert(!NonexistantChunk(startNode));
	assert(allowedNode[startNode & maskUint64_t]);
	uint64_t pos = startNode;
	while (true)
	{
		std::pair<size_t, bool> pairpos { pos & maskUint64_t, pos & firstBitUint64_t };
		if (allowedEdges.getEdges(pairpos).size() != 1) break;
		std::pair<size_t, bool> next = allowedEdges.getEdges(pairpos)[0];
		assert(allowedEdges.hasEdge(MBG::reverse(next), MBG::reverse(pairpos)));
		if (allowedEdges.getEdges(MBG::reverse(next)).size() != 1) break;
		if (next.first == pairpos.first) break;
		if (next.first == (startNode & maskUint64_t)) break;
		pos = next.first + (next.second ? firstBitUint64_t : 0);
		assert(!NonexistantChunk(pos));
		assert(allowedNode[pos & maskUint64_t]);
	}
	pos = pos ^ firstBitUint64_t;
	const uint64_t unitigStart = pos;
	std::vector<uint64_t> result;
	while (true)
	{
		assert(!NonexistantChunk(pos));
		result.emplace_back(pos);
		std::pair<size_t, bool> pairpos { pos & maskUint64_t, pos & firstBitUint64_t };
		if (allowedEdges.getEdges(pairpos).size() != 1) break;
		std::pair<size_t, bool> next = allowedEdges.getEdges(pairpos)[0];
		assert(allowedEdges.hasEdge(MBG::reverse(next), MBG::reverse(pairpos)));
		if (allowedEdges.getEdges(MBG::reverse(next)).size() != 1) break;
		if (next.first == pairpos.first) break;
		if (next.first == (unitigStart & maskUint64_t)) break;
		pos = next.first + (next.second ? firstBitUint64_t : 0);
		assert(allowedNode[pos & maskUint64_t]);
	}
	return result;
}

std::tuple<std::vector<std::vector<uint64_t>>, std::vector<size_t>, std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>> getUnitigs(const std::vector<bool>& allowedNode, const SparseEdgeContainer& allowedEdges, const std::vector<std::vector<size_t>>& lengths, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps)
{
	std::vector<std::vector<uint64_t>> unitigs;
	std::vector<bool> checked;
	checked.resize(allowedNode.size(), false);
	std::vector<std::tuple<uint64_t, size_t, size_t, size_t>> chunkLocationInUnitig;
	chunkLocationInUnitig.resize(allowedNode.size(), std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	for (size_t i = 0; i < allowedNode.size(); i++)
	{
		if (!allowedNode[i]) continue;
		if (checked[i])
		{
			assert(std::get<0>(chunkLocationInUnitig[i]) != std::numeric_limits<size_t>::max());
			continue;
		}
		assert(std::get<0>(chunkLocationInUnitig[i]) == std::numeric_limits<size_t>::max());
		auto unitig = getUnitig(i, allowedNode, allowedEdges);
		assert(unitig.size() >= 1);
		unitigs.emplace_back(unitig);
		for (size_t j = 0; j < unitig.size(); j++)
		{
			uint64_t node = unitig[j];
			assert(!NonexistantChunk(node));
			assert(!checked[node & maskUint64_t]);
			assert(std::get<0>(chunkLocationInUnitig[node & maskUint64_t]) == std::numeric_limits<size_t>::max());
			chunkLocationInUnitig[node & maskUint64_t] = std::make_tuple(unitigs.size()-1 + (node & firstBitUint64_t), j, 0, 0);
			checked[node & maskUint64_t] = true;
		}
	}
	std::vector<size_t> unitigLengths;
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		unitigLengths.emplace_back(0);
		for (size_t j = 0; j < unitigs[i].size(); j++)
		{
			assert(lengths[unitigs[i][j] & maskUint64_t].size() >= 1);
			size_t length = lengths[unitigs[i][j] & maskUint64_t][lengths[unitigs[i][j] & maskUint64_t].size()/2];
			assert(length >= 0);
			assert(length <= 1000000);
			assert(std::get<2>(chunkLocationInUnitig[unitigs[i][j] & maskUint64_t]) == 0);
			assert(std::get<3>(chunkLocationInUnitig[unitigs[i][j] & maskUint64_t]) == 0);
			if (j > 0)
			{
				std::pair<size_t, bool> fromnode { unitigs[i][j-1] & maskUint64_t, unitigs[i][j-1] & firstBitUint64_t };
				std::pair<size_t, bool> tonode { unitigs[i][j] & maskUint64_t, unitigs[i][j] & firstBitUint64_t };
				auto pairkey = MBG::canon(fromnode, tonode);
				std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
				unitigLengths.back() -= edgeOverlaps.at(key);
			}
			std::get<2>(chunkLocationInUnitig[unitigs[i][j] & maskUint64_t]) = unitigLengths.back();
			std::get<3>(chunkLocationInUnitig[unitigs[i][j] & maskUint64_t]) = unitigLengths.back() + length;
			unitigLengths.back() += length;
		}
	}
	for (size_t i = 0; i < chunkLocationInUnitig.size(); i++)
	{
		if (!allowedNode[i])
		{
			assert(std::get<0>(chunkLocationInUnitig[i]) == std::numeric_limits<uint64_t>::max());
			continue;
		}
		size_t unitig = std::get<0>(chunkLocationInUnitig[i]) & maskUint64_t;
		assert(std::get<0>(chunkLocationInUnitig[i]) != std::numeric_limits<uint64_t>::max());
		assert(std::get<3>(chunkLocationInUnitig[i]) != 0);
		assert(std::get<2>(chunkLocationInUnitig[i]) < std::get<3>(chunkLocationInUnitig[i]));
		assert(std::get<3>(chunkLocationInUnitig[i]) <= unitigLengths[unitig]);
		if (std::get<0>(chunkLocationInUnitig[i]) & firstBitUint64_t) continue;
		std::get<1>(chunkLocationInUnitig[i]) = unitigs[unitig].size() - 1 - std::get<1>(chunkLocationInUnitig[i]);
		std::swap(std::get<2>(chunkLocationInUnitig[i]), std::get<3>(chunkLocationInUnitig[i]));
		std::get<2>(chunkLocationInUnitig[i]) = unitigLengths[unitig] - std::get<2>(chunkLocationInUnitig[i]);
		std::get<3>(chunkLocationInUnitig[i]) = unitigLengths[unitig] - std::get<3>(chunkLocationInUnitig[i]);
		assert(std::get<2>(chunkLocationInUnitig[i]) < std::get<3>(chunkLocationInUnitig[i]));
		assert(std::get<3>(chunkLocationInUnitig[i]) <= unitigLengths[unitig]);
	}
	return std::make_tuple(unitigs, unitigLengths, chunkLocationInUnitig);
}

std::vector<std::vector<UnitigPath>> getUnitigPaths(const ChunkUnitigGraph& graph, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<bool>& allowedNode, const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>& chunkLocationInUnitig)
{
	std::vector<std::vector<UnitigPath>> result;
	result.resize(chunksPerRead.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<uint64_t> path;
		uint64_t currentUnitig;
		size_t currentUnitigIndex;
		size_t currentPathLeftClip;
		size_t currentPathRightClip;
		size_t currentReadStart;
		size_t currentReadEnd;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || !allowedNode[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
			{
				if (path.size() > 0)
				{
					result[i].emplace_back();
					result[i].back().path = path;
					result[i].back().pathLeftClip = currentPathLeftClip;
					result[i].back().pathRightClip = currentPathRightClip;
					result[i].back().readStartPos = currentReadStart;
					result[i].back().readEndPos = currentReadEnd;
				}
				path.clear();
				currentUnitig = std::numeric_limits<size_t>::max();
				currentUnitigIndex = 0;
				currentPathLeftClip= 0;
				currentPathRightClip = 0;
				currentReadStart = 0;
				currentReadEnd = 0;
				continue;
			}
			uint64_t chunk = std::get<2>(chunksPerRead[i][j]);
			uint64_t unitig = std::get<0>(chunkLocationInUnitig[chunk]);
			size_t unitigIndex = std::get<1>(chunkLocationInUnitig[chunk]);
			size_t unitigStartPos = std::get<2>(chunkLocationInUnitig[chunk]);
			size_t unitigEndPos = std::get<3>(chunkLocationInUnitig[chunk]);
			size_t readStart = std::get<0>(chunksPerRead[i][j]);
			size_t readEnd = std::get<1>(chunksPerRead[i][j]);
			assert(unitigEndPos > unitigStartPos);
			assert(unitigEndPos <= graph.unitigLengths[unitig & maskUint64_t]);
			if ((chunk ^ firstBitUint64_t) & firstBitUint64_t)
			{
				unitig ^= firstBitUint64_t;
				unitigIndex = unitigs[unitig & maskUint64_t].size()-1-unitigIndex;
				std::swap(unitigStartPos, unitigEndPos);
				unitigStartPos = graph.unitigLengths[unitig & maskUint64_t] - unitigStartPos;
				unitigEndPos = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
			}
			assert(unitigEndPos > unitigStartPos);
			assert(unitigEndPos <= graph.unitigLengths[unitig & maskUint64_t]);
			assert(readEnd > readStart);
			if (path.size() == 0)
			{
				path.emplace_back(unitig);
				currentUnitig = unitig;
				currentUnitigIndex = unitigIndex;
				currentPathLeftClip = unitigStartPos;
				currentPathRightClip = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
				currentReadStart = readStart;
				currentReadEnd = readEnd;
				continue;
			}
			if (currentUnitig == unitig && currentUnitigIndex+1 == unitigIndex)
			{
				assert(readEnd > currentReadEnd);
				assert(graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos < currentPathRightClip);
				currentUnitigIndex = unitigIndex;
				currentReadEnd = readEnd;
				currentPathRightClip = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
				continue;
			}
			if (currentUnitigIndex+1 == unitigs[currentUnitig & maskUint64_t].size() && unitigIndex == 0)
			{
				std::pair<size_t, bool> fromUnitig { currentUnitig & maskUint64_t, currentUnitig & firstBitUint64_t };
				std::pair<size_t, bool> thisUnitig { unitig & maskUint64_t, unitig & firstBitUint64_t };
				if (graph.edges.hasEdge(fromUnitig, thisUnitig))
				{
					assert(readEnd > currentReadEnd);
					path.emplace_back(unitig);
					currentUnitig = unitig;
					currentUnitigIndex = unitigIndex;
					currentReadEnd = readEnd;
					currentPathRightClip = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
					continue;
				}
			}
			result[i].emplace_back();
			result[i].back().path = path;
			result[i].back().pathLeftClip = currentPathLeftClip;
			result[i].back().pathRightClip = currentPathRightClip;
			result[i].back().readStartPos = currentReadStart;
			result[i].back().readEndPos = currentReadEnd;
			path.clear();
			path.emplace_back(unitig);
			currentUnitig = unitig;
			currentUnitigIndex = unitigIndex;
			currentPathLeftClip = unitigStartPos;
			currentPathRightClip = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
			currentReadStart = readStart;
			currentReadEnd = readEnd;
		}
		if (path.size() >= 1)
		{
			result[i].emplace_back();
			result[i].back().path = path;
			result[i].back().pathLeftClip = currentPathLeftClip;
			result[i].back().pathRightClip = currentPathRightClip;
			result[i].back().readStartPos = currentReadStart;
			result[i].back().readEndPos = currentReadEnd;
		}
	}
	return result;
}

std::tuple<ChunkUnitigGraph, std::vector<std::vector<UnitigPath>>> getChunkUnitigGraph(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	ChunkUnitigGraph result;
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
	std::vector<bool> allowedNode;
	SparseEdgeContainer allowedEdges;
	std::tie(allowedNode, allowedEdges) = getAllowedNodesAndEdges(coverages, edgeCoverage, chunksPerRead);
	std::vector<std::vector<uint64_t>> unitigs;
	std::vector<size_t> unitigLengths;
	std::vector<std::tuple<uint64_t, size_t, size_t, size_t>> chunkLocationInUnitig;
	std::tie(unitigs, unitigLengths, chunkLocationInUnitig) = getUnitigs(allowedNode, allowedEdges, lengths, edgeOverlaps);
	result.unitigLengths = unitigLengths;
	result.edges.resize(result.unitigLengths.size());
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		std::pair<size_t, bool> lastNode { unitigs[i].back() & maskUint64_t, unitigs[i].back() & firstBitUint64_t };
		std::pair<size_t, bool> firstNode { unitigs[i][0] & maskUint64_t, (unitigs[i][0] ^ firstBitUint64_t) & firstBitUint64_t };
		assert(!NonexistantChunk(unitigs[i].back()));
		assert(!NonexistantChunk(unitigs[i][0]));
		for (auto edge : allowedEdges.getEdges(lastNode))
		{
			uint64_t targetUnitig = std::get<0>(chunkLocationInUnitig[edge.first]);
			if (!edge.second) targetUnitig ^= firstBitUint64_t;
			auto nodepairkey = MBG::canon(lastNode, edge);
			std::pair<uint64_t, uint64_t> nodekey { nodepairkey.first.first + (nodepairkey.first.second ? firstBitUint64_t : 0), nodepairkey.second.first + (nodepairkey.second.second ? firstBitUint64_t : 0) };
			auto unitigpairkey = MBG::canon(std::make_pair(i, true), std::make_pair(targetUnitig & maskUint64_t, targetUnitig & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> unitigkey { unitigpairkey.first.first + (unitigpairkey.first.second ? firstBitUint64_t : 0), unitigpairkey.second.first + (unitigpairkey.second.second ? firstBitUint64_t : 0) };
			result.edges.addEdge(unitigpairkey.first, unitigpairkey.second);
			result.edges.addEdge(MBG::reverse(unitigpairkey.second), MBG::reverse(unitigpairkey.first));
			result.edgeOverlaps[unitigkey] = edgeOverlaps.at(nodekey);
			result.edgeCoverages[unitigkey] = edgeCoverage.at(nodekey);
		}
		for (auto edge : allowedEdges.getEdges(firstNode))
		{
			uint64_t targetUnitig = std::get<0>(chunkLocationInUnitig[edge.first]);
			if (!edge.second) targetUnitig ^= firstBitUint64_t;
			auto nodepairkey = MBG::canon(firstNode, edge);
			std::pair<uint64_t, uint64_t> nodekey { nodepairkey.first.first + (nodepairkey.first.second ? firstBitUint64_t : 0), nodepairkey.second.first + (nodepairkey.second.second ? firstBitUint64_t : 0) };
			auto unitigpairkey = MBG::canon(std::make_pair(i, false), std::make_pair(targetUnitig & maskUint64_t, targetUnitig & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> unitigkey { unitigpairkey.first.first + (unitigpairkey.first.second ? firstBitUint64_t : 0), unitigpairkey.second.first + (unitigpairkey.second.second ? firstBitUint64_t : 0) };
			result.edges.addEdge(unitigpairkey.first, unitigpairkey.second);
			result.edges.addEdge(MBG::reverse(unitigpairkey.second), MBG::reverse(unitigpairkey.first));
			result.edgeOverlaps[unitigkey] = edgeOverlaps.at(nodekey);
			result.edgeCoverages[unitigkey] = edgeCoverage.at(nodekey);
		}
	}
	result.coverages.resize(result.unitigLengths.size(), 0);
	for (size_t i = 0; i < result.coverages.size(); i++)
	{
		double sum = 0;
		double divisor = 0;
		for (uint64_t node : unitigs[i])
		{
			size_t length = lengths[node & maskUint64_t][lengths[node & maskUint64_t].size() / 2];
			sum += length * coverages[node & maskUint64_t];
			divisor += length;
		}
		result.coverages[i] = sum / divisor;
	}
	result.chunksInUnitig.resize(result.unitigLengths.size());
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		result.chunksInUnitig[i] = unitigs[i];
	}
	auto paths = getUnitigPaths(result, chunksPerRead, allowedNode, unitigs, chunkLocationInUnitig);
	return std::make_tuple(result, paths);
}

void writeUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::string>& readNames, const std::vector<size_t>& rawReadLengths)
{
	std::cerr << "writing unitig graph " << graphFile << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead);
	writeUnitigGraph(graphFile, graph);
	std::cerr << "writing unitig paths " << pathsFile << std::endl;
	writeUnitigPaths(pathsFile, graph, readPaths, readNames, rawReadLengths);
}

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
		if (canCheckThis) checkableIndices.emplace_back(i);
	}
	std::vector<bool> solid;
	solid.resize(checkableIndices.size(), false);
	for (size_t i = 0; i < checkableIndices.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
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
	auto keypair = MBG::canon(forkpair, edgepair);
	std::pair<uint64_t, uint64_t> key { keypair.first.first + (keypair.first.second ? firstBitUint64_t : 0), keypair.second.first + (keypair.second.second ? firstBitUint64_t : 0) };
	size_t compareCoverage = graph.edgeCoverages.at(key);
	for (auto edge : graph.edges.getEdges(forkpair))
	{
		if (edge == edgepair) continue;
		keypair = MBG::canon(forkpair, edge);
		std::pair<uint64_t, uint64_t> key { keypair.first.first + (keypair.first.second ? firstBitUint64_t : 0), keypair.second.first + (keypair.second.second ? firstBitUint64_t : 0) };
		if (graph.edgeCoverages.at(key) >= compareCoverage * 2) return true;
	}
	return false;
}

void cleanTips(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage)
{
	std::cerr << "cleaning tips" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead);
	phmap::flat_hash_set<uint64_t> solidFork;
	phmap::flat_hash_set<uint64_t> acceptableFork;
	{
		std::vector<std::pair<uint64_t, std::vector<std::vector<size_t>>>> forkReads = getUnitigForkReads(graph, readPaths);
		solidFork = getSolidForks(forkReads, approxOneHapCoverage * 0.5);
		acceptableFork = getAcceptableForks(forkReads, solidFork, 3);
	}
	phmap::flat_hash_set<size_t> removeChunks;
	size_t removeUnitigCount = 0;
	size_t removeChunkCount = 0;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.coverages[i] >= approxOneHapCoverage*0.5) continue;
		if (graph.edges.getEdges(std::make_pair(i, true)).size() >= 2) continue;
		if (graph.edges.getEdges(std::make_pair(i, false)).size() >= 2) continue;
		if (graph.edges.getEdges(std::make_pair(i, true)).size() == 1)
		{
			std::pair<size_t, bool> otherpair = graph.edges.getEdges(std::make_pair(i, true))[0];
			uint64_t other = otherpair.first + (otherpair.second ? firstBitUint64_t : 0);
			other ^= firstBitUint64_t;
			if (solidFork.count(other) == 1) continue;
			if (acceptableFork.count(other) == 1) continue;
			if (!hasOtherHigherCoverageEdge(other, i, graph)) continue;
		}
		if (graph.edges.getEdges(std::make_pair(i, false)).size() == 1)
		{
			std::pair<size_t, bool> otherpair = graph.edges.getEdges(std::make_pair(i, false))[0];
			uint64_t other = otherpair.first + (otherpair.second ? firstBitUint64_t : 0);
			other ^= firstBitUint64_t;
			if (solidFork.count(other) == 1) continue;
			if (acceptableFork.count(other) == 1) continue;
			if (!hasOtherHigherCoverageEdge(other, i + firstBitUint64_t, graph)) continue;
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

void makeGraph(const MatchIndex& matchIndex, const std::vector<std::string>& readNames, const std::vector<size_t>& rawReadLengths, const std::vector<TwobitString>& readSequences, const size_t numThreads)
{
	std::vector<bool> useTheseChunks;
	useTheseChunks.resize(matchIndex.numWindowChunks() - matchIndex.numUniqueChunks(), true);
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead = getChunksPerRead(matchIndex, rawReadLengths, useTheseChunks);
	removeContainedChunks(chunksPerRead);
	splitPerLength(chunksPerRead);
	writeUnitigGraph("graph-round1.gfa", "paths1.gaf", chunksPerRead, readNames, rawReadLengths);
	splitPerBaseCounts(readSequences, chunksPerRead, numThreads);
	writeUnitigGraph("graph-round2.gfa", "paths2.gaf", chunksPerRead, readNames, rawReadLengths);
	splitPerMinHashes(readSequences, chunksPerRead, numThreads);
	writeUnitigGraph("graph-round3.gfa", "paths3.gaf", chunksPerRead, readNames, rawReadLengths);
	splitPerPhasingKmersWithinChunk(readSequences, chunksPerRead, numThreads);
	writeUnitigGraph("graph-round4.gfa", "paths4.gaf", chunksPerRead, readNames, rawReadLengths);
	splitPerSequenceIdentity(readSequences, chunksPerRead, numThreads);
	writeUnitigGraph("graph-round5.gfa", "paths5.gaf", chunksPerRead, readNames, rawReadLengths);
	splitPerPhasingKmersWithinChunk(readSequences, chunksPerRead, numThreads);
	writeUnitigGraph("graph-round6.gfa", "paths6.gaf", chunksPerRead, readNames, rawReadLengths);
//	splitPerAllUniqueKmerSVs(readSequences, chunksPerRead, numThreads);
	splitPerPhasingKmersWithinChunk(readSequences, chunksPerRead, numThreads);
	writeUnitigGraph("graph-round7.gfa", "paths7.gaf", chunksPerRead, readNames, rawReadLengths);
	splitPerInterchunkPhasedKmers(readSequences, chunksPerRead, numThreads);
	writeUnitigGraph("graph-round8.gfa", "paths8.gaf", chunksPerRead, readNames, rawReadLengths);
	splitPerInterchunkPhasedKmers(readSequences, chunksPerRead, numThreads);
//	countGoodKmersInChunks(readSequences, chunksPerRead, 1);
//	countGoodishKmersInChunks(readSequences, chunksPerRead, 1);
	writeUnitigGraph("graph-round9.gfa", "paths9.gaf", chunksPerRead, readNames, rawReadLengths);
	cleanTips(readSequences, chunksPerRead, numThreads, 20);
	writeUnitigGraph("graph-round10.gfa", "paths10.gaf", chunksPerRead, readNames, rawReadLengths);
	cleanTips(readSequences, chunksPerRead, numThreads, 20);
	writeUnitigGraph("graph-round11.gfa", "paths11.gaf", chunksPerRead, readNames, rawReadLengths);
	removeSingleCopyChunks(chunksPerRead);
	writeUnitigGraph("graph-round12.gfa", "paths12.gaf", chunksPerRead, readNames, rawReadLengths);
	removeHighCoverageChunks(chunksPerRead, 60);
	writeUnitigGraph("graph-round13.gfa", "paths13.gaf", chunksPerRead, readNames, rawReadLengths);
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	std::vector<std::string> readFiles;
	for (size_t i = 5; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	MatchIndex matchIndex { k, numWindows, windowSize };
	ReadStorage storage;
	std::vector<TwobitString> readSequences;
	for (auto file : readFiles)
	{
		std::mutex indexMutex;
		std::mutex sequenceMutex;
		storage.iterateReadsFromFile(file, numThreads, false, [&matchIndex, &indexMutex, &sequenceMutex, &readSequences](size_t readName, const std::string& sequence)
		{
			matchIndex.addMatchesFromRead(readName, indexMutex, sequence);
			{
				std::lock_guard<std::mutex> lock { sequenceMutex };
				while (readSequences.size() <= readName) readSequences.emplace_back();
				readSequences[readName] = sequence;
			}
		});
	}
	std::cerr << readSequences.size() << " reads" << std::endl;
	std::cerr << matchIndex.numWindowChunks() << " distinct windowchunks" << std::endl;
	std::cerr << matchIndex.numUniqueChunks() << " windowchunks have only one read" << std::endl;
	matchIndex.clearConstructionVariablesAndCompact();
	const std::vector<size_t>& readBasepairLengths = storage.getRawReadLengths();
	const std::vector<std::string>& readNames = storage.getNames();
	makeGraph(matchIndex, readNames, readBasepairLengths, readSequences, numThreads);
}
