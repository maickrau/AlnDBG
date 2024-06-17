#include <queue>
#include <mutex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "RankBitvector.h"
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
#include "CompressedStringIndex.h"

double mismatchFraction;

size_t popcount(uint64_t x);

auto getTime()
{
	return std::chrono::steady_clock::now();
}

std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
	size_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	return std::to_string(milliseconds / 1000) + "," + std::to_string(milliseconds % 1000) + " s";
}

auto programStartTime = getTime();

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
	std::vector<std::vector<size_t>> unitigChunkBreakpointPositions;
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
	size_t pathLeftClip;
	size_t pathRightClip;
};

class ConsensusString
{
public:
	TwobitString bases;
	std::vector<std::pair<size_t, size_t>> Nchunks;
};

std::string getReadSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t index)
{
	if (index < sequenceIndex.size())
	{
		return sequenceIndex.getSequence(index);
	}
	assert(index < sequenceIndex.size()*2);
	return MBG::revCompRaw(sequenceIndex.getSequence(index - sequenceIndex.size()));
}

std::string getChunkSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t read, const size_t index)
{
	if (read < sequenceIndex.size())
	{
		auto t = chunksPerRead[read][index];
		std::string result = sequenceIndex.getSubstring(read, std::get<0>(t), std::get<1>(t)-std::get<0>(t)+1);
		if (std::get<2>(t) & firstBitUint64_t)
		{
		}
		else
		{
			result = MBG::revCompRaw(result);
		}
		return result;
	}
	else
	{
		assert(read < sequenceIndex.size()*2);
		auto t = chunksPerRead[read][index];
		std::string result = sequenceIndex.getSubstring(read - sequenceIndex.size(), rawReadLengths[read-sequenceIndex.size()] - 1 - std::get<1>(t), std::get<1>(t) - std::get<0>(t) + 1);
		if (std::get<2>(t) & firstBitUint64_t)
		{
			result = MBG::revCompRaw(result);
		}
		return result;
	}
}

void writeUnitigPaths(const std::string& filename, const ChunkUnitigGraph& unitigGraph, const std::vector<std::vector<UnitigPath>>& readPaths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths)
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
			std::string name;
			size_t length;
			if (i < sequenceIndex.size())
			{
				name = sequenceIndex.getName(i);
				length = rawReadLengths[i];
			}
			else
			{
				name = sequenceIndex.getName(i - sequenceIndex.size()) + "_bw";
				length = rawReadLengths[i - sequenceIndex.size()];
			}
			file << name << "\t" << length << "\t" << readPaths[i][j].readPartInPathnode[0].first << "\t" << readPaths[i][j].readPartInPathnode.back().second << "\t+\t" << pathstr << "\t" << pathLength << "\t" << readPaths[i][j].pathLeftClip << "\t" << (pathLength - readPaths[i][j].pathRightClip) << "\t" << (pathLength - readPaths[i][j].pathLeftClip- readPaths[i][j].pathRightClip) << "\t" << (pathLength - readPaths[i][j].pathLeftClip- readPaths[i][j].pathRightClip) << "\t60" << std::endl;
		}
	}
}

void writeUnitigGraph(const std::string& filename, const ChunkUnitigGraph& unitigGraph)
{
	assert(unitigGraph.unitigLengths.size() == unitigGraph.coverages.size());
	std::ofstream file { filename };
	for (size_t i = 0; i < unitigGraph.unitigLengths.size(); i++)
	{
		file << "S\t" << i << "\t*\tLN:i:" << unitigGraph.unitigLengths[i] << "\tll:f:" << unitigGraph.coverages[i] << "\tFC:f:" << (unitigGraph.unitigLengths[i] * unitigGraph.coverages[i]);
		file << "\tsg:Z:";
		for (size_t chunk : unitigGraph.chunksInUnitig[i])
		{
			file << ((chunk & firstBitUint64_t) ? ">" : "<") << (chunk & maskUint64_t);
		}
		file << "\tbr:Z:";
		for (size_t j = 0; j < unitigGraph.unitigChunkBreakpointPositions[i].size(); j++)
		{
			if (j > 0) file << ",";
			file << unitigGraph.unitigChunkBreakpointPositions[i][j];
		}
		file << std::endl;
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
	phmap::flat_hash_map<uint64_t, size_t> maxCoverageEdge;
	for (auto pair : edgeCoverage)
	{
		maxCoverageEdge[pair.first.first] = std::max(maxCoverageEdge[pair.first.first], pair.second);
		maxCoverageEdge[pair.first.second ^ firstBitUint64_t] = std::max(maxCoverageEdge[pair.first.second ^ firstBitUint64_t], pair.second);
	}
	for (auto pair : edgeCoverage)
	{
		if (pair.second < 2) continue;
		if (pair.second <= 3)
		{
			if (maxCoverageEdge.at(pair.first.first) > pair.second*3)
			{
				if (maxCoverageEdge.at(pair.first.second ^ firstBitUint64_t) > pair.second*3)
				{
					continue;
				}
			}
		}
		assert(!NonexistantChunk(pair.first.first));
		assert(!NonexistantChunk(pair.first.second));
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
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				lastTip = std::numeric_limits<size_t>::max();
				continue;
			}
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
		if (otherEnd == t)
		{
			if (coverages[t & maskUint64_t] == 2) continue; // bunch of fake inverted should-be-duplex reads in ONT, skip them. todo this should also check the coverage is all from just one read
		}
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
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				lastTip = std::numeric_limits<size_t>::max();
				lastTipIndex = std::numeric_limits<size_t>::max();
				continue;
			}
			if (!allowedNode[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			if (lastTip != std::numeric_limits<size_t>::max() && tip.count(std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t) == 1)
			{
				assert(lastTipIndex < j);
				uint64_t tipHere = std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t;
				if (allowedTips.count(std::make_pair(lastTip, tipHere)) == 1)
				{
					for (size_t k = lastTipIndex; k < j; k++)
					{
						assert(!NonexistantChunk(std::get<2>(chunksPerRead[i][k])));
						assert(!NonexistantChunk(std::get<2>(chunksPerRead[i][k+1])));
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
		assert(node < allowedNode.size());
		allowedNode[node] = true;
	}
	for (size_t i = 0; i < allowedEdges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : allowedEdges.getEdges(fw))
		{
			assert(edge.first < allowedNode.size());
			assert(allowedNode[i]);
			assert(allowedNode[edge.first]);
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : allowedEdges.getEdges(bw))
		{
			assert(edge.first < allowedNode.size());
			assert(allowedNode[i]);
			assert(allowedNode[edge.first]);
		}
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
			size_t overlap = 0;
			if (std::get<1>(chunksPerRead[i][j-1]) > std::get<0>(chunksPerRead[i][j])) overlap = std::get<1>(chunksPerRead[i][j-1]) - std::get<0>(chunksPerRead[i][j]);
			auto pairkey = MBG::canon(std::make_pair(prev & maskUint64_t, prev & firstBitUint64_t), std::make_pair(curr & maskUint64_t, curr & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
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

size_t findBiggestSmallerIndex(const std::vector<size_t>& nums, const size_t value)
{
	if (nums.size() < 5)
	{
		assert(nums[0] < value);
		for (size_t i = 1; i < nums.size(); i++)
		{
			if (nums[i] > value) return i-1;
		}
		return nums.size()-1;
	}
	size_t low = 0;
	size_t high = nums.size();
	while (high > low)
	{
		size_t mid = (low + high) / 2;
		assert(nums[mid] != value);
		if (nums[mid] < value)
		{
			low = mid+1;
		}
		else if (nums[mid] > value)
		{
			high = mid;
		}
	}
	assert(low == high);
	assert(low <= nums.size());
	assert(low == nums.size() || nums[low] > value);
	assert(low > 0);
	assert(nums[low-1] < value);
	return low-1;
}

void splitPerLength(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting by length" << std::endl;
	const double differenceFraction = 0.02;
	const size_t differenceConstant = 50;
	std::vector<std::vector<size_t>> lengthsPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			while ((std::get<2>(t) & maskUint64_t) >= lengthsPerChunk.size()) lengthsPerChunk.emplace_back();
			lengthsPerChunk[std::get<2>(t) & maskUint64_t].emplace_back(std::get<1>(t) - std::get<0>(t));
		}
	}
	std::vector<std::vector<size_t>> splitters;
	splitters.resize(lengthsPerChunk.size());
	iterateMultithreaded(0, lengthsPerChunk.size(), numThreads, [&lengthsPerChunk, &splitters, differenceFraction, differenceConstant](const size_t i)
	{
		std::sort(lengthsPerChunk[i].begin(), lengthsPerChunk[i].end());
		splitters[i].emplace_back(0);
		for (size_t j = 1; j < lengthsPerChunk[i].size(); j++)
		{
			size_t distance = lengthsPerChunk[i][j] - lengthsPerChunk[i][j-1];
			if (distance < std::max((size_t)(lengthsPerChunk[i][j] * differenceFraction), (size_t)differenceConstant)) continue;
			splitters[i].emplace_back((lengthsPerChunk[i][j] + lengthsPerChunk[i][j-1])/2);
		}
	});
	std::vector<std::vector<size_t>> splitterToNode;
	splitterToNode.resize(lengthsPerChunk.size());
	size_t nextNum = 0;
	for (size_t i = 0; i < splitters.size(); i++)
	{
		for (size_t j = 0; j < splitters[i].size(); j++)
		{
			splitterToNode[i].emplace_back(nextNum);
			nextNum += 1;
		}
	}
	iterateMultithreaded(0, chunksPerRead.size(), numThreads, [&splitters, &splitterToNode, &chunksPerRead](const size_t i)
	{
		for (auto& t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			assert(std::get<1>(t) > std::get<0>(t));
			size_t distance = std::get<1>(t) - std::get<0>(t);
			assert((std::get<2>(t) & maskUint64_t) < splitters.size());
			size_t index = findBiggestSmallerIndex(splitters[std::get<2>(t) & maskUint64_t], distance);
			assert(index < splitterToNode[std::get<2>(t) & maskUint64_t].size());
			std::get<2>(t) = (std::get<2>(t) & firstBitUint64_t) + splitterToNode[std::get<2>(t) & maskUint64_t][index];
		}
	});
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
			if (m < start+k-1) continue;
			callback(kmer, m - (start+k-1));
		}
	}
	else
	{
		for (size_t m = end; m >= start && m < baseSequence.size(); m--)
		{
			kmer <<= 2;
			kmer += 3 - baseSequence.get(m);
			kmer &= mask;
			if (m > end-k+1) continue;
			callback(kmer, end - m + 1 - k);
		}
	}
}

void checkPhasablePair(const phmap::flat_hash_map<size_t, size_t>& bwForks, const phmap::flat_hash_map<size_t, size_t>& fwForks, std::vector<std::vector<size_t>>& phaseIdentities)
{
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> pairCoverage;
	phmap::flat_hash_set<size_t> foundOccurrences;
	phmap::flat_hash_set<size_t> foundBw;
	phmap::flat_hash_set<size_t> foundFw;
	for (auto bwfork : bwForks)
	{
		foundBw.insert(bwfork.second);
		foundOccurrences.insert(bwfork.first);
		if (fwForks.count(bwfork.first) == 0) continue;
		size_t bw = bwfork.second;
		size_t fw = fwForks.at(bwfork.first);
		pairCoverage[std::make_pair(bw, fw)] += 1;
	}
	for (auto fwfork : fwForks)
	{
		foundFw.insert(fwfork.second);
		foundOccurrences.insert(fwfork.first);
	}
	if (foundFw.size() != foundBw.size()) return;
	if (foundOccurrences.size() < phaseIdentities.size()) return;
	size_t coverageInPairs = 0;
	for (auto pair : pairCoverage)
	{
		coverageInPairs += pair.second;
		if (pair.second < 5) return;
	}
	if (coverageInPairs < phaseIdentities.size() * 0.95) return;
	phmap::flat_hash_map<size_t, size_t> parent;
	for (auto pair : pairCoverage)
	{
		assert((pair.first.first & firstBitUint64_t) == 0);
		assert((pair.first.second & firstBitUint64_t) == 0);
		if (parent.count(pair.first.first) == 0) parent[pair.first.first] = pair.first.first;
		if (parent.count(pair.first.second + firstBitUint64_t) == 0) parent[pair.first.second + firstBitUint64_t] = pair.first.second + firstBitUint64_t;
		merge(parent, pair.first.first, pair.first.second + firstBitUint64_t);
	}
	phmap::flat_hash_map<size_t, size_t> clusterCount;
	for (auto pair : pairCoverage)
	{
		auto cluster = find(parent, pair.first.first);
		clusterCount[cluster] += pair.second;
	}
	if (clusterCount.size() < 2) return;
	for (auto kmer : foundBw)
	{
		if (parent.count(kmer) == 0) return;
	}
	for (auto kmer : foundFw)
	{
		if (parent.count(kmer + firstBitUint64_t) == 0) return;
	}
	size_t inserted = 0;
	for (auto bwfork : bwForks)
	{
		if (fwForks.count(bwfork.first) == 1)
		{
			assert(find(parent, bwfork.second) == find(parent, fwForks.at(bwfork.first) + firstBitUint64_t));
		}
		phaseIdentities[bwfork.first].emplace_back(find(parent, bwfork.second));
		inserted += 1;
	}
	for (auto fwfork : fwForks)
	{
		if (bwForks.count(fwfork.first) == 1) continue;
		phaseIdentities[fwfork.first].emplace_back(find(parent, fwfork.second + firstBitUint64_t));
		inserted += 1;
	}
	assert(inserted == phaseIdentities.size());
}

std::vector<std::vector<std::vector<size_t>>> getAlleleOccurrences(const std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>>& alleles, const size_t numOccurrences)
{
	std::vector<std::vector<std::vector<size_t>>> result;
	size_t clusterStart = 0;
	for (size_t i = 1; i <= alleles.size(); i++)
	{
		if (i < alleles.size() && std::get<0>(alleles[i]) == std::get<0>(alleles[i-1]) && std::get<1>(alleles[i]) == std::get<1>(alleles[i-1]) && std::get<2>(alleles[i]) == std::get<2>(alleles[i-1]) && std::get<3>(alleles[i]) == std::get<3>(alleles[i-1])) continue;
		if (i - clusterStart != numOccurrences)
		{
			clusterStart = i;
			continue;
		}
		std::vector<bool> foundOccurrences;
		foundOccurrences.resize(numOccurrences, false);
		bool valid = true;
		for (size_t j = clusterStart; j < i; j++)
		{
			assert(std::get<4>(alleles[j]) < foundOccurrences.size());
			if (foundOccurrences[std::get<4>(alleles[j])])
			{
				valid = false;
				break;
			}
			foundOccurrences[std::get<4>(alleles[j])] = true;
		}
		for (size_t j = 0; j < foundOccurrences.size(); j++)
		{
			if (!foundOccurrences[j])
			{
				valid = false;
				break;
			}
		}
		if (!valid)
		{
			clusterStart = i;
			continue;
		}
		phmap::flat_hash_map<size_t, size_t> alleleCounts;
		for (size_t j = clusterStart; j < i; j++)
		{
			alleleCounts[std::get<5>(alleles[j])] += 1;
		}
		if (alleleCounts.size() < 2)
		{
			clusterStart = i;
			continue;
		}
		for (auto pair : alleleCounts)
		{
			if (pair.second < 5)
			{
				valid = false;
				break;
			}
		}
		if (!valid)
		{
			clusterStart = i;
			continue;
		}
		result.emplace_back();
		phmap::flat_hash_map<size_t, size_t> alleleMap;
		for (auto pair : alleleCounts)
		{
			alleleMap[pair.first] = result.back().size();
//			std::cerr << "allele occurrence " << result.size()-1 << " allele " << result.back().size() << " is " << pair.first << std::endl;
			result.back().emplace_back();
		}
		for (size_t j = clusterStart; j < i; j++)
		{
			assert(alleleMap.count(std::get<5>(alleles[j])) == 1);
			assert(alleleMap.at(std::get<5>(alleles[j])) < result.back().size());
			result.back()[alleleMap.at(std::get<5>(alleles[j]))].emplace_back(std::get<4>(alleles[j]));
		}
		clusterStart = i;
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			std::sort(result[i][j].begin(), result[i][j].end());
		}
	}
	return result;
}

bool allelesMatchTwoVariantsThreeHaps(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right)
{
	if (left.size() != 2) return false;
	if (right.size() != 2) return false;
	if (left[0].size() < 10) return false;
	if (left[1].size() < 10) return false;
	if (right[0].size() < 10) return false;
	if (right[1].size() < 10) return false;
	size_t zerozero = intersectSize(left[0], right[0]);
	size_t zeroone = intersectSize(left[0], right[1]);
	size_t onezero = intersectSize(left[1], right[0]);
	size_t oneone = intersectSize(left[1], right[1]);
	if (zerozero > 0 && zerozero < 10) return false;
	if (zeroone > 0 && zeroone < 10) return false;
	if (onezero > 0 && onezero < 10) return false;
	if (oneone > 0 && oneone < 10) return false;
	if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	return true;
}

bool allelesMatchPerfectly(const std::vector<std::vector<size_t>>& left, const std::vector<std::vector<size_t>>& right)
{
	if (left.size() != right.size()) return false;
	phmap::flat_hash_set<size_t> rightMatched;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i].size() < 5) return false;
		size_t uniqueMatch = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < right.size(); j++)
		{
			if (right[j].size() < 5) return false;
			if (left[i].size() != right[j].size()) continue;
			bool match = true;
			for (size_t k = 0; k < left[i].size(); k++)
			{
				if (left[i][k] != right[j][k])
				{
					match = false;
					break;
				}
			}
			if (!match) continue;
			if (uniqueMatch != std::numeric_limits<size_t>::max()) return false;
			uniqueMatch = j;
		}
		if (uniqueMatch == std::numeric_limits<size_t>::max()) return false;
		if (rightMatched.count(uniqueMatch) == 1) return false;
		rightMatched.insert(uniqueMatch);
	}
	assert(rightMatched.size() == right.size());
	return true;
}

std::string getAlleleStr(const TwobitString& readSequence, const size_t readStart, const size_t readEnd, const bool fw)
{
	assert(readEnd > readStart);
	std::string str = readSequence.substr(readStart, readEnd - readStart);
	if (!fw) str = MBG::revCompRaw(str);
	return str;
}

size_t getAllele(const TwobitString& readSequence, const size_t readStart, const size_t readEnd, phmap::flat_hash_map<std::string, size_t>& alleleNumbers, const bool fw)
{
	std::string allele = getAlleleStr(readSequence, readStart, readEnd, fw);
	if (alleleNumbers.count(allele) == 0)
	{
		size_t result = alleleNumbers.size();
		alleleNumbers[allele] = result;
	}
	return alleleNumbers.at(allele);
}

template <typename KmerType, typename F>
phmap::flat_hash_map<size_t, std::vector<std::pair<double, double>>> iterateSolidKmersWithKmerType(const std::vector<TwobitString>& chunkSequences, const size_t kmerSize, const size_t minSolidThreshold, const bool allowTwoAlleleRepeats, const bool needClusters, F callback)
{
	assert(chunkSequences.size() < (size_t)std::numeric_limits<uint32_t>::max());
	assert(chunkSequences.size() >= 1);
	double repetitiveClusteringDistance = 1.0;
	// if one sequence is too small then don't do anything even if theoretically other sequences could have solid kmers
	size_t smallestSequence = chunkSequences[0].size();
	for (size_t i = 0; i < chunkSequences.size(); i++)
	{
		if (chunkSequences[i].size() < kmerSize+2) return phmap::flat_hash_map<size_t, std::vector<std::pair<double, double>>> {};
		smallestSequence = std::min(smallestSequence, chunkSequences[i].size());
	}
	std::vector<size_t> currentPosPerOccurrence;
	// 40 blocks so it advances 2.5% of sequence every block so kmers seen before block cannot cluster with kmers seen after block if no kmers in block
	size_t numBlocks = 40;
	// unless the strings are small, then no blocks
	if (smallestSequence < 2000)
	{
		repetitiveClusteringDistance = 100 * 20.0/(double)smallestSequence;
	}
	if (smallestSequence < 2000 || smallestSequence < 40*kmerSize*2)
	{
		numBlocks = 1;
	}
	for (size_t j = 0; j < chunkSequences.size(); j++)
	{
		assert(chunkSequences[j].size() < (size_t)std::numeric_limits<uint16_t>::max());
		currentPosPerOccurrence.emplace_back(0);
	}
	std::vector<std::vector<std::pair<double, double>>> validClusters;
	std::vector<std::vector<std::tuple<uint16_t, uint16_t, KmerType>>> solidKmers; // pos, cluster, kmer
	solidKmers.resize(chunkSequences.size());
	phmap::flat_hash_map<KmerType, size_t> kmerToIndex;
	std::vector<std::vector<std::tuple<float, uint32_t, uint16_t>>> kmers; // kmer -> (extrapolated-pos, occurrence, occurrence-pos)
	std::vector<size_t> nextCluster;
	std::vector<size_t> kmersDroppingOutOfContext;
	std::vector<KmerType> actualKmer;
	for (size_t offsetBreakpoint = 0; offsetBreakpoint < numBlocks; offsetBreakpoint++)
	{
		phmap::flat_hash_set<size_t> inContext;
		for (size_t j = 0; j < chunkSequences.size(); j++)
		{
			size_t nextOffset = (double)(offsetBreakpoint+1)/(double)numBlocks * (chunkSequences[j].size()-1);
			assert(offsetBreakpoint != (numBlocks-1) || nextOffset == chunkSequences[j].size()-1);
			assert(nextOffset > currentPosPerOccurrence[j] + kmerSize);
			size_t startPos, endPos;
			startPos = currentPosPerOccurrence[j];
			endPos = nextOffset;
			iterateKmers(chunkSequences[j], startPos, endPos, true, kmerSize, [&kmers, &nextCluster, &kmerToIndex, &inContext, &currentPosPerOccurrence, &validClusters, &actualKmer, &chunkSequences, kmerSize, j, offsetBreakpoint, startPos, endPos](const size_t kmer, const size_t pos)
			{
				assert(pos < endPos-startPos);
				assert(pos + currentPosPerOccurrence[j] < chunkSequences[j].size());
				assert((size_t)(pos+currentPosPerOccurrence[j]) < (size_t)std::numeric_limits<uint16_t>::max());
				double extrapolatedPos = 100.0 * (double)(pos+currentPosPerOccurrence[j]+((double)kmerSize/2.0)) / (double)chunkSequences[j].size();
				size_t index = kmerToIndex.size();
				auto found = kmerToIndex.find(kmer);
				if (found != kmerToIndex.end())
				{
					index = found->second;
				}
				else
				{
					kmerToIndex[kmer] = index;
					kmers.emplace_back();
					nextCluster.emplace_back(0);
					validClusters.emplace_back();
					actualKmer.emplace_back(kmer);
				}
				kmers[index].emplace_back(extrapolatedPos, j, (pos+currentPosPerOccurrence[j]));
				inContext.insert(index);
			});
			currentPosPerOccurrence[j] = nextOffset - kmerSize + 2;
		}
		if (offsetBreakpoint == (numBlocks-1))
		{
			kmersDroppingOutOfContext.insert(kmersDroppingOutOfContext.end(), inContext.begin(), inContext.end());
			inContext.clear();
		}
		for (size_t checkindex = kmersDroppingOutOfContext.size()-1; checkindex < kmersDroppingOutOfContext.size(); checkindex--)
		{
			size_t index = kmersDroppingOutOfContext[checkindex];
			assert(index < kmers.size());
			if (inContext.count(index) == 1) continue;
			if (kmers[index].size() >= minSolidThreshold)
			{
				std::sort(kmers[index].begin(), kmers[index].end());
				size_t clusterStart = 0;
				std::vector<bool> occurrencesHere;
				occurrencesHere.resize(chunkSequences.size(), false);
				bool currentlyNonrepetitive = true;
				occurrencesHere[std::get<1>(kmers[index][0])] = true;
				size_t clusterNum = nextCluster[index];
				for (size_t j = 1; j <= kmers[index].size(); j++)
				{
					if (j < kmers[index].size() && std::get<0>(kmers[index][j]) < std::get<0>(kmers[index][j-1])+repetitiveClusteringDistance)
					{
						if (occurrencesHere[std::get<1>(kmers[index][j])])
						{
							currentlyNonrepetitive = false;
						}
						occurrencesHere[std::get<1>(kmers[index][j])] = true;
						continue;
					}
					if (j-clusterStart >= minSolidThreshold && (currentlyNonrepetitive || allowTwoAlleleRepeats))
					{
						double minPos = std::get<0>(kmers[index][clusterStart]);
						double maxPos = std::get<0>(kmers[index][j-1]);
						if (minPos + repetitiveClusteringDistance > maxPos)
						{
							if (!currentlyNonrepetitive && allowTwoAlleleRepeats)
							{
								size_t countOccurrences = 0;
								for (size_t k = 0; k < occurrencesHere.size(); k++)
								{
									countOccurrences += occurrencesHere[k] ? 1 : 0;
								}
								if (countOccurrences >= minSolidThreshold)
								{
									phmap::flat_hash_map<size_t, size_t> kmerCountPerOccurrence;
									for (size_t k = clusterStart; k < j; k++)
									{
										kmerCountPerOccurrence[std::get<1>(kmers[index][k])] += 1;
									}
									phmap::flat_hash_set<size_t> foundCounts;
									for (auto pair2 : kmerCountPerOccurrence)
									{
										foundCounts.insert(pair2.second);
									}
									if (foundCounts.size() == 2)
									{
										size_t firstCount = *foundCounts.begin();
										size_t secondCount = *(++foundCounts.begin());
										assert(firstCount != secondCount);
										for (auto pair2 : kmerCountPerOccurrence)
										{
											if (pair2.second == firstCount)
											{
												// todo pos should be something reasonable not this. but nothing uses poses of repetitive kmers as of comment time
												assert(clusterNum < (size_t)std::numeric_limits<uint16_t>::max());
												solidKmers[pair2.first].emplace_back(std::numeric_limits<uint16_t>::max(), clusterNum, index);
											}
											else
											{
												assert(pair2.second == secondCount);
												// todo pos should be something reasonable not this. but nothing uses poses of repetitive kmers as of comment time
												assert(clusterNum+1 < (size_t)std::numeric_limits<uint16_t>::max());
												solidKmers[pair2.first].emplace_back(std::numeric_limits<uint16_t>::max(), clusterNum+1, index);
											}
										}
										clusterNum += 2;
										validClusters[index].emplace_back(minPos-0.01, maxPos+0.01);
										validClusters[index].emplace_back(minPos-0.01, maxPos+0.01);
									}
								}
							}
							else
							{
								for (size_t k = clusterStart; k < j; k++)
								{
									assert((size_t)std::get<2>(kmers[index][k]) < (size_t)std::numeric_limits<uint16_t>::max());
									assert(std::get<0>(kmers[index][k]) >= minPos);
									assert(std::get<0>(kmers[index][k]) <= maxPos);
									assert(clusterNum < (size_t)std::numeric_limits<uint16_t>::max());
									solidKmers[std::get<1>(kmers[index][k])].emplace_back(std::get<2>(kmers[index][k]), clusterNum, index);
								}
								clusterNum += 1;
								validClusters[index].emplace_back(minPos-0.01, maxPos+0.01);
							}
						}
					}
					currentlyNonrepetitive = true;
					occurrencesHere.assign(occurrencesHere.size(), false);
					if (j < kmers[index].size())
					{
						occurrencesHere[std::get<1>(kmers[index][j])] = true;
					}
					clusterStart = j;
				}
				nextCluster[index] = clusterNum;
			}
			std::swap(kmersDroppingOutOfContext[checkindex], kmersDroppingOutOfContext.back());
			kmersDroppingOutOfContext.pop_back();
			{
				std::remove_reference<decltype(kmers[index])>::type tmp;
				std::swap(tmp, kmers[index]);
			}
		}
		kmersDroppingOutOfContext.insert(kmersDroppingOutOfContext.end(), inContext.begin(), inContext.end());
	}
	for (size_t i = 0; i < kmers.size(); i++)
	{
		assert(kmers[i].size() == 0);
	}
	phmap::flat_hash_map<size_t, std::vector<std::pair<double, double>>> fixedClusters;
	if (needClusters)
	{
		for (size_t i = 0; i < validClusters.size(); i++)
		{
			fixedClusters[actualKmer[i]];
			std::swap(fixedClusters[actualKmer[i]], validClusters[i]);
		}
	}
	for (size_t i = 0; i < solidKmers.size(); i++)
	{
		std::sort(solidKmers[i].begin(), solidKmers[i].end());
		for (auto kmer : solidKmers[i])
		{
			// occurrence ID, chunk startpos in read, chunk endpos in read, chunk, kmer, clusternum, kmerpos in chunk
			callback(i, 0, chunkSequences[i].size(), std::numeric_limits<size_t>::max(), (size_t)actualKmer[std::get<2>(kmer)], (size_t)std::get<1>(kmer), (size_t)std::get<0>(kmer));
		}
	}
	return fixedClusters;
}

template <typename F>
phmap::flat_hash_map<size_t, std::vector<std::pair<double, double>>> iterateSolidKmers(const std::vector<TwobitString>& chunkSequences, const size_t kmerSize, const size_t minSolidThreshold, const bool allowTwoAlleleRepeats, const bool needClusters, F callback)
{
	if (kmerSize < 16)
	{
		return iterateSolidKmersWithKmerType<uint32_t>(chunkSequences, kmerSize, minSolidThreshold, allowTwoAlleleRepeats, needClusters, callback);
	}
	else
	{
		return iterateSolidKmersWithKmerType<size_t>(chunkSequences, kmerSize, minSolidThreshold, allowTwoAlleleRepeats, needClusters, callback);
	}
}

size_t getHammingdistance(const RankBitvector& left, const RankBitvector& right, const size_t maxEdits)
{
	assert(left.size() == right.size());
	if (left.size() == 0) return 0;
	const std::vector<uint64_t>& leftbits = left.getBits();
	const std::vector<uint64_t>& rightbits = right.getBits();
	assert(leftbits.size() == rightbits.size());
	assert(left.size() <= leftbits.size()*64);
	size_t result = 0;
	for (size_t i = 0; i*64+64 < left.size(); i++)
	{
		uint64_t notEqual = leftbits[i] ^ rightbits[i];
		result += popcount(notEqual);
		if (result > maxEdits) return result;
	}
	for (size_t i = ((size_t)((left.size()-1)/64))*64; i < left.size(); i++)
	{
		if (left.get(i) != right.get(i))
		{
			result += 1;
		}
	}
	return result;
}

bool siteIsPerfectlyPhased(const std::vector<RankBitvector>& columns, const size_t left, const size_t right)
{
	size_t zerozero = 0;
	size_t zeroone = 0;
	size_t onezero = 0;
	size_t oneone = 0;
	assert(columns[left].size() == columns[right].size());
	const std::vector<uint64_t>& leftbits = columns[left].getBits();
	const std::vector<uint64_t>& rightbits = columns[right].getBits();
	assert(leftbits.size() == rightbits.size());
	for (size_t i = 0; i*64+64 < columns[left].size(); i++)
	{
		oneone += popcount(leftbits[i] & rightbits[i]);
		zeroone += popcount(~leftbits[i] & rightbits[i]);
		zerozero += popcount(~leftbits[i] & ~rightbits[i]);
		onezero += popcount(leftbits[i] & ~rightbits[i]);
		if (oneone > 0 || zerozero > 0) return false;
	}
	size_t remainingBits = columns[left].size() % 64;
	uint64_t mask = 0xFFFFFFFFFFFFFFFFull;
	if (remainingBits > 0) // remainingbits 0 means the whole block is valid
	{
		mask >>= (64-remainingBits);
	}
	oneone += popcount((leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zeroone += popcount((~leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zerozero += popcount((~leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	onezero += popcount((leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	if (oneone == 0 && zerozero == 0 && onezero >= 5 && zeroone >= 5) return true;
	return false;
}

bool siteIsPhasedTwoVariantsThreeHaps(const std::vector<RankBitvector>& columns, const size_t left, const size_t right)
{
	size_t zerozero = 0;
	size_t zeroone = 0;
	size_t onezero = 0;
	size_t oneone = 0;
	assert(columns[left].size() == columns[right].size());
	const std::vector<uint64_t>& leftbits = columns[left].getBits();
	const std::vector<uint64_t>& rightbits = columns[right].getBits();
	assert(leftbits.size() == rightbits.size());
	for (size_t i = 0; i*64+64 < columns[left].size(); i++)
	{
		oneone += popcount(leftbits[i] & rightbits[i]);
		zeroone += popcount(~leftbits[i] & rightbits[i]);
		zerozero += popcount(~leftbits[i] & ~rightbits[i]);
		onezero += popcount(leftbits[i] & ~rightbits[i]);
		if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	}
	size_t remainingBits = columns[left].size() % 64;
	uint64_t mask = 0xFFFFFFFFFFFFFFFFull;
	if (remainingBits > 0) // remainingbits 0 means the whole block is valid
	{
		mask >>= (64-remainingBits);
	}
	oneone += popcount((leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zeroone += popcount((~leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zerozero += popcount((~leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	onezero += popcount((leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	if (zerozero > 0 && zerozero < 10) return false;
	if (onezero > 0 && onezero < 10) return false;
	if (zeroone > 0 && zeroone < 10) return false;
	if (oneone > 0 && oneone < 10) return false;
	if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	return true;
}

bool siteIsInformative(const std::vector<RankBitvector>& columns, const size_t left, const size_t right)
{
	size_t zerozero = 0;
	size_t zeroone = 0;
	size_t onezero = 0;
	size_t oneone = 0;
	assert(columns[left].size() == columns[right].size());
	const std::vector<uint64_t>& leftbits = columns[left].getBits();
	const std::vector<uint64_t>& rightbits = columns[right].getBits();
	assert(leftbits.size() == rightbits.size());
	for (size_t i = 0; i*64+64 < columns[left].size(); i++)
	{
		oneone += popcount(leftbits[i] & rightbits[i]);
		zeroone += popcount(~leftbits[i] & rightbits[i]);
		zerozero += popcount(~leftbits[i] & ~rightbits[i]);
		onezero += popcount(leftbits[i] & ~rightbits[i]);
		if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	}
	size_t remainingBits = columns[left].size() % 64;
	uint64_t mask = 0xFFFFFFFFFFFFFFFFull;
	if (remainingBits > 0) // remainingbits 0 means the whole block is valid
	{
		mask >>= (64-remainingBits);
	}
	oneone += popcount((leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zeroone += popcount((~leftbits[(columns[left].size()-1)/64] & rightbits[(columns[left].size()-1)/64]) & mask);
	zerozero += popcount((~leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	onezero += popcount((leftbits[(columns[left].size()-1)/64] & ~rightbits[(columns[left].size()-1)/64]) & mask);
	if (zerozero + zeroone < 5) return false;
	if (onezero + oneone < 5) return false;
	if (zerozero + onezero < 5) return false;
	if (zeroone + oneone < 5) return false;
	if (zerozero > 0 && zerozero < 5) return false;
	if (zeroone > 0 && zeroone < 5) return false;
	if (onezero > 0 && onezero < 5) return false;
	if (oneone > 0 && oneone < 5) return false;
	if (zerozero > 0 && zeroone > 0 && onezero > 0 && oneone > 0) return false;
	return true;
}

template <typename T>
void filterVector(std::vector<T>& vec, const std::vector<bool>& keep)
{
	std::vector<size_t> getIndexFrom;
	for (size_t i = 0; i < keep.size(); i++)
	{
		if (!keep[i]) continue;
		getIndexFrom.emplace_back(i);
	}
	assert(getIndexFrom.size() <= keep.size());
	assert(getIndexFrom.size() == 0 || getIndexFrom.back() < keep.size());
	for (size_t i = 0; i < getIndexFrom.size(); i++)
	{
		vec[i] = vec[getIndexFrom[i]];
	}
	vec.resize(getIndexFrom.size());
}

void filterMatrix(std::vector<RankBitvector>& matrix, const std::vector<bool>& keep)
{
	std::vector<size_t> getIndexFrom;
	for (size_t i = 0; i < keep.size(); i++)
	{
		if (!keep[i]) continue;
		getIndexFrom.emplace_back(i);
	}
	assert(getIndexFrom.size() <= keep.size());
	assert(getIndexFrom.size() == 0 || getIndexFrom.back() < keep.size());
	for (size_t i = 0; i < matrix.size(); i++)
	{
		assert(matrix[i].size() == keep.size());
		for (size_t k = 0; k < getIndexFrom.size(); k++)
		{
			matrix[i].set(k, matrix[i].get(getIndexFrom[k]));
		}
		matrix[i].resize(getIndexFrom.size());
	}
}

template <typename F>
void iterateChunksByCoverage(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, F callback)
{
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
	iterationOrder.reserve(occurrencesPerChunk.size());
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		iterationOrder.emplace_back(i);
	}
	std::sort(iterationOrder.begin(), iterationOrder.end(), [&occurrencesPerChunk](size_t left, size_t right) { return occurrencesPerChunk[left].size() > occurrencesPerChunk[right].size(); });
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&occurrencesPerChunk, &iterationOrder, callback](const size_t iterationIndex)
	{
		callback(iterationOrder[iterationIndex], occurrencesPerChunk);
	});
}

void checkAddMismatch(std::vector<std::vector<std::pair<size_t, size_t>>>& closestNeighborAndMismatchesPerLine, std::vector<size_t>& maxEditsHerePerLine, std::vector<size_t>& countEqualToMaxEditsPerLine, const size_t addHere, const size_t addThis, const size_t mismatches, const size_t countNeighbors)
{
	if (closestNeighborAndMismatchesPerLine[addHere].size() < countNeighbors)
	{
		closestNeighborAndMismatchesPerLine[addHere].emplace_back(mismatches, addThis);
		if (closestNeighborAndMismatchesPerLine[addHere].size() == countNeighbors)
		{
			maxEditsHerePerLine[addHere] = 0;
			for (size_t l = 0; l < closestNeighborAndMismatchesPerLine[addHere].size(); l++)
			{
				if (closestNeighborAndMismatchesPerLine[addHere][l].first < maxEditsHerePerLine[addHere]) continue;
				if (closestNeighborAndMismatchesPerLine[addHere][l].first > maxEditsHerePerLine[addHere]) countEqualToMaxEditsPerLine[addHere] = 0;
				maxEditsHerePerLine[addHere] = closestNeighborAndMismatchesPerLine[addHere][l].first;
				countEqualToMaxEditsPerLine[addHere] += 1;
			}
		}
	}
	else
	{
		if (mismatches <= maxEditsHerePerLine[addHere])
		{
			closestNeighborAndMismatchesPerLine[addHere].emplace_back(mismatches, addThis);
			if (mismatches == maxEditsHerePerLine[addHere])
			{
				countEqualToMaxEditsPerLine[addHere] += 1;
			}
			else
			{
				if (closestNeighborAndMismatchesPerLine[addHere].size() > countNeighbors+countEqualToMaxEditsPerLine[addHere])
				{
					size_t newMaxEdits = 0;
					countEqualToMaxEditsPerLine[addHere] = 0;
					for (size_t l = closestNeighborAndMismatchesPerLine[addHere].size()-1; l < closestNeighborAndMismatchesPerLine[addHere].size(); l--)
					{
						if (closestNeighborAndMismatchesPerLine[addHere][l].first == maxEditsHerePerLine[addHere])
						{
							std::swap(closestNeighborAndMismatchesPerLine[addHere][l], closestNeighborAndMismatchesPerLine[addHere].back());
							closestNeighborAndMismatchesPerLine[addHere].pop_back();
						}
						else
						{
							assert(closestNeighborAndMismatchesPerLine[addHere][l].first < maxEditsHerePerLine[addHere]);
							if (closestNeighborAndMismatchesPerLine[addHere][l].first > newMaxEdits) countEqualToMaxEditsPerLine[addHere] = 0;
							newMaxEdits = std::max(newMaxEdits, closestNeighborAndMismatchesPerLine[addHere][l].first);
							if (closestNeighborAndMismatchesPerLine[addHere][l].first == newMaxEdits) countEqualToMaxEditsPerLine[addHere] += 1;
						}
					}
					assert(newMaxEdits < maxEditsHerePerLine[addHere]);
					maxEditsHerePerLine[addHere] = newMaxEdits;
				}
			}
		}
	}
}

template <typename F>
std::vector<size_t> getFastTransitiveClosure(const size_t itemCount, const size_t maxDistance, F distanceFunction)
{
	std::vector<size_t> parent;
	parent.resize(itemCount);
	for (size_t i = 0; i < itemCount; i++)
	{
		parent[i] = i;
	}
	std::vector<size_t> clusterExample;
	std::vector<std::vector<size_t>> clusterAdditionals;
	std::vector<size_t> clusterMaxDistance;
	for (size_t i = 0; i < itemCount; i++)
	{
		bool found = false;
		for (size_t j = 0; j < clusterExample.size(); j++)
		{
			size_t distance = distanceFunction(i, clusterExample[j], maxDistance);
			if (distance <= maxDistance)
			{
				clusterAdditionals[j].emplace_back(i);
				found = true;
				clusterMaxDistance[j] = std::max(distance, clusterMaxDistance[j]);
				merge(parent, i, clusterExample[j]);
				break;
			}
		}
		if (!found)
		{
			clusterExample.emplace_back(i);
			clusterAdditionals.emplace_back();
			clusterMaxDistance.emplace_back(0);
		}
	}
	for (size_t i = 1; i < clusterExample.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (find(parent, clusterExample[i]) == find(parent, clusterExample[j])) continue;
			size_t distance = distanceFunction(clusterExample[i], clusterExample[j], maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j]);
			if (distance <= maxDistance)
			{
				merge(parent, clusterExample[i], clusterExample[j]);
				continue;
			}
			else if (distance <= maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j])
			{
				bool found = false;
				for (size_t k : clusterAdditionals[i])
				{
					if (distanceFunction(k, clusterExample[j], maxDistance) <= maxDistance)
					{
						found = true;
						merge(parent, k, clusterExample[j]);
						break;
					}
				}
				if (!found)
				{
					for (size_t l : clusterAdditionals[j])
					{
						if (distanceFunction(l, clusterExample[i], maxDistance) <= maxDistance)
						{
							found = true;
							merge(parent, l, clusterExample[i]);
							break;
						}
					}
				}
				if (!found)
				{
					for (size_t k : clusterAdditionals[i])
					{
						for (size_t l : clusterAdditionals[j])
						{
							if (distanceFunction(k, l, maxDistance) <= maxDistance)
							{
								found = true;
								merge(parent, k, l);
								break;
							}
						}
						if (found) break;
					}
				}
			}
		}
	}
	return parent;
}

std::vector<RankBitvector> getCorrectedMatrix(const std::vector<RankBitvector>& matrix, const size_t countNeighbors)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> closestNeighborAndMismatchesPerLine;
	std::vector<size_t> maxEditsHerePerLine;
	std::vector<size_t> countEqualToMaxEditsPerLine;
	closestNeighborAndMismatchesPerLine.resize(matrix.size());
	maxEditsHerePerLine.resize(matrix.size(), std::numeric_limits<size_t>::max());
	countEqualToMaxEditsPerLine.resize(matrix.size(), 0);
	for (size_t distance = 1; distance < matrix.size(); distance++)
	{
		for (size_t high = distance; high < matrix.size(); high++)
		{
			size_t low = high - distance;
			size_t mismatches = getHammingdistance(matrix[low], matrix[high], std::max(maxEditsHerePerLine[low], maxEditsHerePerLine[high]));
			checkAddMismatch(closestNeighborAndMismatchesPerLine, maxEditsHerePerLine, countEqualToMaxEditsPerLine, low, high, mismatches, countNeighbors);
			checkAddMismatch(closestNeighborAndMismatchesPerLine, maxEditsHerePerLine, countEqualToMaxEditsPerLine, high, low, mismatches, countNeighbors);
		}
	}
	std::vector<RankBitvector> correctedMatrix;
	correctedMatrix.resize(matrix.size());
	for (size_t j = 0; j < matrix.size(); j++)
	{
		correctedMatrix[j].resize(matrix[j].size());
		closestNeighborAndMismatchesPerLine[j].emplace_back(0, j);
		std::sort(closestNeighborAndMismatchesPerLine[j].begin(), closestNeighborAndMismatchesPerLine[j].end());
		size_t endIndex = countNeighbors;
		while (endIndex+1 < closestNeighborAndMismatchesPerLine[j].size() && closestNeighborAndMismatchesPerLine[j][endIndex+1].first == closestNeighborAndMismatchesPerLine[j][countNeighbors].first)
		{
			endIndex += 1;
		}
		for (size_t k = 0; k < matrix[j].size(); k++)
		{
			size_t ones = 0;
			for (size_t m = 0; m <= endIndex; m++)
			{
				if (matrix[closestNeighborAndMismatchesPerLine[j][m].second].get(k)) ones += 1;
			}
			correctedMatrix[j].set(k, ones >= (endIndex+1)/2);
		}
	}
	return correctedMatrix;
}

bool rowEqual(const RankBitvector& left, const RankBitvector& right)
{
	assert(left.size() == right.size());
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left.get(i) != right.get(i)) return false;
	}
	return true;
}

void sortAndInterleave(std::vector<RankBitvector>& correctedMatrix)
{
	std::sort(correctedMatrix.begin(), correctedMatrix.end(), [](const RankBitvector& left, const RankBitvector& right)
	{
		for(size_t i = 0; i < left.size(); i++)
		{
			if (!left.get(i) && right.get(i)) return true;
			if (left.get(i) && !right.get(i)) return false;
		}
		return false;
	});
	size_t countDistinctRows = 1;
	for (size_t i = 1; i < correctedMatrix.size(); i++)
	{
		if (!rowEqual(correctedMatrix[i], correctedMatrix[i-1])) countDistinctRows += 1;
	}
	std::swap(correctedMatrix[1], correctedMatrix.back());
	size_t pos = 2;
	for (size_t div = 2; div <= 64; div *= 2)
	{
		for (size_t off = 1; off < div; off += 2)
		{
			std::swap(correctedMatrix[pos], correctedMatrix[correctedMatrix.size()/div*off]);
			pos += 1;
		}
	}
	assert(pos == 65);
}

void splitPerCorrectedKmerPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by corrected kmer phasing" << std::endl;
	const size_t countNeighbors = 5;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, countNeighbors, kmerSize](size_t dummy)
	{
		while (true)
		{
			std::vector<std::pair<size_t, size_t>> chunkBeingDone;
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunksNeedProcessing.size() == 0) return;
				std::swap(chunkBeingDone, chunksNeedProcessing.back());
				chunksNeedProcessing.pop_back();
			}
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				chunksDoneProcessing.emplace_back();
				std::cerr << "corrected kmer skipped chunk with coverage " << chunkBeingDone.size() << std::endl;
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerClusterToColumn;
			std::vector<RankBitvector> matrix; // doesn't actually use rank in any way, but exposes uint64 for bitparallel hamming distance calculation
			std::vector<std::pair<size_t, size_t>> clusters;
			matrix.resize(chunkBeingDone.size());
			std::vector<TwobitString> chunkSequences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
			}
			auto validClusters = iterateSolidKmers(chunkSequences, kmerSize, 5, true, true, [&kmerClusterToColumn, &clusters, &matrix](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				assert(occurrenceID < matrix.size());
				size_t column = std::numeric_limits<size_t>::max();
				if (kmerClusterToColumn.count(std::make_pair(kmer, clusterIndex)) == 1)
				{
					column = kmerClusterToColumn.at(std::make_pair(kmer, clusterIndex));
				}
				else
				{
					column = kmerClusterToColumn.size();
					kmerClusterToColumn[std::make_pair(kmer, clusterIndex)] = column;
					clusters.emplace_back(std::make_pair(kmer, clusterIndex));
				}
				if (matrix[occurrenceID].size() <= column)
				{
					assert(column < kmerClusterToColumn.size());
					matrix[occurrenceID].resize(kmerClusterToColumn.size());
				}
				assert(!matrix[occurrenceID].get(column));
				matrix[occurrenceID].set(column, true);
			});
			std::vector<std::pair<double, double>> kmerLocation;
			kmerLocation.resize(clusters.size(), std::make_pair(-1, -1));
			for (auto pair : kmerClusterToColumn)
			{
				size_t index = pair.second;
				assert(kmerLocation[index].first == -1);
				kmerLocation[index] = validClusters.at(pair.first.first)[pair.first.second];
			}
			for (size_t i = 0; i < kmerLocation.size(); i++)
			{
				assert(kmerLocation[i].first != -1);
			}
			for (size_t j = 0; j < matrix.size(); j++)
			{
				if (matrix[j].size() == kmerClusterToColumn.size()) continue;
				assert(matrix[j].size() < kmerClusterToColumn.size());
				matrix[j].resize(kmerClusterToColumn.size());
			}
			size_t columnsInUnfiltered = matrix[0].size();
			std::vector<bool> covered;
			covered.resize(matrix[0].size(), false);
			size_t columnsInCovered = 0;
			std::vector<RankBitvector> columns;
			columns.resize(matrix[0].size());
			for (size_t j = 0; j < matrix[0].size(); j++)
			{
				columns[j].resize(matrix.size());
				size_t zeros = 0;
				size_t ones = 0;
				for (size_t k = 0; k < matrix.size(); k++)
				{
					assert(matrix[k].size() == matrix[0].size());
					columns[j].set(k, matrix[k].get(j));
					if (matrix[k].get(j))
					{
						ones += 1;
					}
					else
					{
						zeros += 1;
					}
				}
				if (zeros >= 5 && ones >= 5)
				{
					covered[j] = true;
					columnsInCovered += 1;
				}
			}
			if (columnsInCovered < 2)
			{
				auto endTime = getTime();
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "corrected kmer skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " due to covered " << columnsInCovered << " time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<bool> informativeSite;
			informativeSite.resize(matrix[0].size(), false);
			assert(clusters.size() == informativeSite.size());
			assert(clusters.size() == kmerClusterToColumn.size());
			for (size_t j = 0; j < informativeSite.size(); j++)
			{
				if (!covered[j]) continue;
				for (size_t k = 0; k < j; k++)
				{
					if (!covered[k]) continue;
					if (informativeSite[j] && informativeSite[k]) continue;
					assert(validClusters.at(clusters[j].first).size() > clusters[j].second);
					assert(validClusters.at(clusters[k].first).size() > clusters[k].second);
					if (validClusters.at(clusters[j].first)[clusters[j].second].first < validClusters.at(clusters[k].first)[clusters[k].second].second + 1)
					{
						if (validClusters.at(clusters[j].first)[clusters[j].second].second + 1 > validClusters.at(clusters[k].first)[clusters[k].second].first)
						{
							continue;
						}
					}
					if (siteIsInformative(columns, j, k) || siteIsInformative(columns, k, j))
					{
						informativeSite[j] = true;
						informativeSite[k] = true;
					}
				}
			}
			size_t columnsInInformative = 0;
			for (size_t j = 0; j < informativeSite.size(); j++)
			{
				if (informativeSite[j]) columnsInInformative += 1;
			}
			if (columnsInInformative < 2)
			{
				auto endTime = getTime();
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "corrected kmer skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " due to informative " << columnsInInformative << " time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			filterVector(columns, informativeSite);
			filterVector(kmerLocation, informativeSite);
			filterMatrix(matrix, informativeSite);
			std::vector<RankBitvector> correctedMatrix = getCorrectedMatrix(matrix, countNeighbors);
			assert(correctedMatrix.size() == matrix.size());
			assert(correctedMatrix.size() == chunkBeingDone.size());
			if (correctedMatrix.size() >= 1024)
			{
				std::vector<RankBitvector> columnMatrix = correctedMatrix;
				sortAndInterleave(columnMatrix);
				for (size_t j = 0; j < columnMatrix[0].size(); j++)
				{
					for (size_t k = 0; k < columnMatrix.size(); k++)
					{
						assert(columnMatrix[k].size() == columnMatrix[0].size());
						columns[j].set(k, columnMatrix[k].get(j));
					}
				}
			}
			std::vector<bool> phasingSite;
			phasingSite.resize(matrix[0].size(), false);
			std::vector<size_t> phasedSoFar;
			std::vector<size_t> notPhasedSoFar;
			notPhasedSoFar.emplace_back(0);
			for (size_t j = 1; j < matrix[0].size(); j++)
			{
				for (auto k : phasedSoFar)
				{
					if (kmerLocation[j].first < kmerLocation[k].second+1 && kmerLocation[j].second + 1 > kmerLocation[k].first) continue;
					if (siteIsPerfectlyPhased(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						break;
					}
					else if (siteIsPhasedTwoVariantsThreeHaps(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						break;
					}
				}
				for (size_t kindex = notPhasedSoFar.size()-1; kindex < notPhasedSoFar.size(); kindex--)
				{
					size_t k = notPhasedSoFar[kindex];
					if (kmerLocation[j].first < kmerLocation[k].second+1 && kmerLocation[j].second + 1 > kmerLocation[k].first) continue;
					if (siteIsPerfectlyPhased(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						phasedSoFar.emplace_back(k);
						std::swap(notPhasedSoFar[kindex], notPhasedSoFar.back());
						notPhasedSoFar.pop_back();
					}
					else if (siteIsPhasedTwoVariantsThreeHaps(columns, j, k))
					{
						phasingSite[j] = true;
						phasingSite[k] = true;
						phasedSoFar.emplace_back(k);
						std::swap(notPhasedSoFar[kindex], notPhasedSoFar.back());
						notPhasedSoFar.pop_back();
					}
				}
				if (phasingSite[j])
				{
					phasedSoFar.emplace_back(j);
				}
				else
				{
					notPhasedSoFar.emplace_back(j);
				}
			}
			filterMatrix(correctedMatrix, phasingSite);
			size_t numPhasingSites = correctedMatrix[0].size();
			std::vector<size_t> parent = getFastTransitiveClosure(correctedMatrix.size(), 0, [&correctedMatrix](const size_t i, const size_t j, const size_t maxDist) { return getHammingdistance(correctedMatrix[i], correctedMatrix[j], maxDist); });
			auto endTime = getTime();
			size_t nextNum = 0;
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
			if (keyToNode.size() == 1)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "corrected kmer splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " phasing sites " << numPhasingSites << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.resize(keyToNode.size());
			for (size_t j = 0; j < parent.size(); j++)
			{
				chunkResult[keyToNode.at(find(parent, j))].emplace_back(chunkBeingDone[j]);
			}
			// sort smallest last, so emplace-pop-swap puts biggest on top of chunksNeedProcessing
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (keyToNode.size() >= 2) countSplitted += 1;
				std::cerr << "corrected kmer splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " phasing sites " << numPhasingSites << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				while (chunkResult.size() > 0)
				{
					chunksNeedProcessing.emplace_back();
					std::swap(chunksNeedProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
			}
		}
	});
	for (size_t i = 0; i < chunksDoneProcessing.size(); i++)
	{
		for (auto pair : chunksDoneProcessing[i])
		{
			std::get<2>(chunksPerRead[pair.first][pair.second]) = i + (std::get<2>(chunksPerRead[pair.first][pair.second]) & firstBitUint64_t);
		}
	}
	std::cerr << "corrected kmer phasing splitted " << countSplitted << " chunks" << std::endl;
}

class ScopedCounterIncrementer
{
public:
	ScopedCounterIncrementer(std::atomic<size_t>& counter) :
	counter(counter)
	{
		counter += 1;
	}
	~ScopedCounterIncrementer()
	{
		counter -= 1;
	}
private:
	std::atomic<size_t>& counter;
};

void splitPerNearestNeighborPhasing(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by nearest neighbor phasing" << std::endl;
	const size_t countNeighbors = 5;
	const size_t countDifferences = 100;
	size_t countSplitted = 0;
	std::atomic<size_t> threadsRunning;
	threadsRunning = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, &threadsRunning, countNeighbors, countDifferences, kmerSize](size_t dummy)
	{
		while (true)
		{
			std::vector<std::pair<size_t, size_t>> chunkBeingDone;
			bool gotOne = false;
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunksNeedProcessing.size() >= 1)
				{
					std::swap(chunkBeingDone, chunksNeedProcessing.back());
					chunksNeedProcessing.pop_back();
					gotOne = true;
				}
				else
				{
					if (threadsRunning == 0) return;
				}
			}
			if (!gotOne)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				continue;
			}
			ScopedCounterIncrementer threadCounter { threadsRunning };
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "skip chunk with size " << chunkBeingDone.size() << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerClusterToColumn;
			std::vector<RankBitvector> matrix; // doesn't actually use rank in any way, but exposes uint64 for bitparallel hamming distance calculation
			std::vector<std::pair<size_t, size_t>> clusters;
			matrix.resize(chunkBeingDone.size());
			std::vector<TwobitString> chunkSequences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
			}
			auto validClusters = iterateSolidKmers(chunkSequences, kmerSize, 5, true, true, [&kmerClusterToColumn, &clusters, &matrix](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				assert(occurrenceID < matrix.size());
				size_t column = std::numeric_limits<size_t>::max();
				if (kmerClusterToColumn.count(std::make_pair(kmer, clusterIndex)) == 1)
				{
					column = kmerClusterToColumn.at(std::make_pair(kmer, clusterIndex));
				}
				else
				{
					column = kmerClusterToColumn.size();
					kmerClusterToColumn[std::make_pair(kmer, clusterIndex)] = column;
					clusters.emplace_back(std::make_pair(kmer, clusterIndex));
				}
				if (matrix[occurrenceID].size() <= column)
				{
					assert(column < kmerClusterToColumn.size());
					matrix[occurrenceID].resize(kmerClusterToColumn.size());
				}
				assert(!matrix[occurrenceID].get(column));
				matrix[occurrenceID].set(column, true);
			});
			for (size_t j = 0; j < matrix.size(); j++)
			{
				if (matrix[j].size() == kmerClusterToColumn.size()) continue;
				assert(matrix[j].size() < kmerClusterToColumn.size());
				matrix[j].resize(kmerClusterToColumn.size());
			}
			size_t columnsInUnfiltered = matrix[0].size();
			std::vector<bool> covered;
			covered.resize(matrix[0].size(), false);
			size_t columnsInCovered = 0;
			std::vector<RankBitvector> columns;
			columns.resize(matrix[0].size());
			for (size_t j = 0; j < matrix[0].size(); j++)
			{
				columns[j].resize(matrix.size());
				size_t zeros = 0;
				size_t ones = 0;
				for (size_t k = 0; k < matrix.size(); k++)
				{
					assert(matrix[k].size() == matrix[0].size());
					columns[j].set(k, matrix[k].get(j));
				}
				const std::vector<uint64_t>& bits = columns[j].getBits();
				for (size_t k = 0; k*64+64 <= columns[j].size(); k++)
				{
					ones += popcount(bits[k]);
					zeros += 64-popcount(bits[k]);
					if (zeros >= 5 && ones >= 5) break;
				}
				if (zeros < 5 || ones < 5)
				{
					for (size_t k = ((int)(columns[j].size()/64))*64; k < columns[j].size(); k++)
					{
						if (columns[j].get(k))
						{
							ones += 1;
						}
						else
						{
							zeros += 1;
						}
					}
				}
				if (zeros >= 5 && ones >= 5)
				{
					covered[j] = true;
					columnsInCovered += 1;
				}
			}
			if (columnsInCovered < countDifferences)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				auto endTime = getTime();
				std::cerr << "nearest neighbor skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " due to covered " << columnsInCovered << ", time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<bool> informativeSite;
			informativeSite.resize(matrix[0].size(), false);
			assert(clusters.size() == informativeSite.size());
			assert(clusters.size() == kmerClusterToColumn.size());
			for (size_t j = 0; j < informativeSite.size(); j++)
			{
				if (!covered[j]) continue;
				for (size_t k = 0; k < j; k++)
				{
					if (!covered[k]) continue;
					if (informativeSite[j] && informativeSite[k]) continue;
					assert(validClusters.at(clusters[j].first).size() > clusters[j].second);
					assert(validClusters.at(clusters[k].first).size() > clusters[k].second);
					if (validClusters.at(clusters[j].first)[clusters[j].second].first < validClusters.at(clusters[k].first)[clusters[k].second].second + 1)
					{
						if (validClusters.at(clusters[j].first)[clusters[j].second].second + 1 > validClusters.at(clusters[k].first)[clusters[k].second].first)
						{
							continue;
						}
					}
					if (siteIsInformative(columns, j, k) || siteIsInformative(columns, k, j))
					{
						informativeSite[j] = true;
						informativeSite[k] = true;
					}
				}
			}
			size_t columnsInInformative = 0;
			for (size_t i = 0; i < informativeSite.size(); i++)
			{
				if (informativeSite[i]) columnsInInformative += 1;
			}
			if (columnsInInformative < countDifferences)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				auto endTime = getTime();
				std::cerr << "nearest neighbor skipped chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " due to informative " << columnsInInformative << ", time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			filterMatrix(matrix, informativeSite);
			std::vector<RankBitvector> correctedMatrix = getCorrectedMatrix(matrix, countNeighbors);
			assert(correctedMatrix.size() == matrix.size());
			assert(correctedMatrix.size() == chunkBeingDone.size());
			std::vector<size_t> parent = getFastTransitiveClosure(correctedMatrix.size(), countDifferences, [&correctedMatrix](const size_t i, const size_t j, const size_t maxDist) { return getHammingdistance(correctedMatrix[i], correctedMatrix[j], maxDist); });
			auto endTime = getTime();
			size_t nextNum = 0;
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
			if (keyToNode.size() == 1)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "nearest neighbor splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.resize(keyToNode.size());
			for (size_t j = 0; j < parent.size(); j++)
			{
				chunkResult[keyToNode.at(find(parent, j))].emplace_back(chunkBeingDone[j]);
			}
			// sort smallest last, so emplace-pop-swap puts biggest on top of chunksNeedProcessing
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (keyToNode.size() >= 2) countSplitted += 1;
				std::cerr << "nearest neighbor splitted chunk with coverage " << chunkBeingDone.size() << " columns " << columnsInUnfiltered << " covered " << columnsInCovered << " informative " << columnsInInformative << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				while (chunkResult.size() > 0)
				{
					chunksNeedProcessing.emplace_back();
					std::swap(chunksNeedProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
			}
		}
	});
	for (size_t i = 0; i < chunksDoneProcessing.size(); i++)
	{
		for (auto pair : chunksDoneProcessing[i])
		{
			std::get<2>(chunksPerRead[pair.first][pair.second]) = i + (std::get<2>(chunksPerRead[pair.first][pair.second]) & firstBitUint64_t);
		}
	}
	std::cerr << "nearest neighbor phasing splitted " << countSplitted << " chunks" << std::endl;
}

void splitPerAllelePhasingWithinChunk(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by allele phasing" << std::endl;
	size_t countSplitted = 0;
	size_t countSplittedTo = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &countSplittedTo, &rawReadLengths, kmerSize](size_t dummy)
	{
		while (true)
		{
			std::vector<std::pair<size_t, size_t>> chunkBeingDone;
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunksNeedProcessing.size() == 0) return;
				std::swap(chunkBeingDone, chunksNeedProcessing.back());
				chunksNeedProcessing.pop_back();
			}
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
//				std::cerr << "skip chunk with coverage " << chunkBeingDone.size() << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>> alleles; // startkmer, startcluster, endkmer, endcluster, occurrenceID, allele
			size_t lastOccurrence = std::numeric_limits<size_t>::max();
			size_t lastKmer = std::numeric_limits<size_t>::max();
			size_t lastCluster = std::numeric_limits<size_t>::max();
			size_t lastKmerPos = std::numeric_limits<size_t>::max();
			std::vector<TwobitString> chunkSequences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
			}
			phmap::flat_hash_map<std::string, size_t> alleleNumbers;
			iterateSolidKmers(chunkSequences, kmerSize, chunkBeingDone.size(), false, false, [&chunkBeingDone, &chunkSequences, &alleles, &lastOccurrence, &lastKmer, &lastCluster, &lastKmerPos, &alleleNumbers, kmerSize](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				if (occurrenceID != lastOccurrence || lastKmer == std::numeric_limits<size_t>::max())
				{
					lastOccurrence = occurrenceID;
					lastKmer = kmer;
					lastCluster = clusterIndex;
					lastKmerPos = pos;
					return;
				}
				size_t allele;
				allele = getAllele(chunkSequences[occurrenceID], lastKmerPos, pos + kmerSize, alleleNumbers, true);
				alleles.emplace_back(lastKmer, lastCluster, kmer, clusterIndex, occurrenceID, allele);
				lastOccurrence = occurrenceID;
				lastKmer = kmer;
				lastCluster = clusterIndex;
				lastKmerPos = pos;
			});
			std::sort(alleles.begin(), alleles.end());
			std::vector<std::vector<std::vector<size_t>>> occurrencesPerAlleleSite;
			occurrencesPerAlleleSite = getAlleleOccurrences(alleles, chunkBeingDone.size());
			std::vector<bool> informativeOccurrence;
			informativeOccurrence.resize(occurrencesPerAlleleSite.size(), false);
			for (size_t j = 1; j < occurrencesPerAlleleSite.size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					if (informativeOccurrence[j] && informativeOccurrence[k]) continue;
					if (allelesMatchPerfectly(occurrencesPerAlleleSite[j], occurrencesPerAlleleSite[k]))
					{
						informativeOccurrence[j] = true;
						informativeOccurrence[k] = true;
					}
					else if (allelesMatchTwoVariantsThreeHaps(occurrencesPerAlleleSite[j], occurrencesPerAlleleSite[k]))
					{
						informativeOccurrence[j] = true;
						informativeOccurrence[k] = true;
					}
				}
			}
			std::vector<std::vector<size_t>> occurrenceAlleles;
			occurrenceAlleles.resize(chunkBeingDone.size());
			for (size_t site = 0; site < occurrencesPerAlleleSite.size(); site++)
			{
				if (!informativeOccurrence[site]) continue;
				for (size_t allele = 0; allele < occurrencesPerAlleleSite[site].size(); allele++)
				{
					for (auto occurrence : occurrencesPerAlleleSite[site][allele])
					{
						occurrenceAlleles[occurrence].emplace_back(allele);
					}
				}
			}
			std::vector<size_t> orderedOccurrences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				orderedOccurrences.emplace_back(j);
			}
			std::sort(orderedOccurrences.begin(), orderedOccurrences.end(), [&occurrenceAlleles](size_t left, size_t right)
			{
				assert(occurrenceAlleles[left].size() == occurrenceAlleles[right].size());
				for (size_t i = 0; i < occurrenceAlleles[left].size(); i++)
				{
					if (occurrenceAlleles[left][i] < occurrenceAlleles[right][i]) return true;
					if (occurrenceAlleles[left][i] > occurrenceAlleles[right][i]) return false;
				}
				return false;
			});
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.emplace_back();
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				if (j > 0 && occurrenceAlleles[orderedOccurrences[j]] != occurrenceAlleles[orderedOccurrences[j-1]])
				{
					chunkResult.emplace_back();
				}
				chunkResult.back().emplace_back(chunkBeingDone[orderedOccurrences[j]]);
			}
			if (chunkResult.size() == 1)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
//				std::cerr << "chunk with coverage " << chunkBeingDone.size() << " not split" << std::endl;
//				std::cerr << "alleles " << alleles.size() << " occurrences " << occurrencesPerAlleleSite.size() << " occurrenceAlleles " << occurrenceAlleles.size() << " informative sites " << occurrenceAlleles[0].size() << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkResult[0]);
				continue;
			}
			assert(chunkResult.size() >= 2);
			// sort smallest last, so emplace-swap-pop puts biggest on top of chunksNeedProcessing
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
//				std::cerr << "chunk with coverage " << chunkBeingDone.size() << " split to " << chunkResult.size() << std::endl;
//				std::cerr << "alleles " << alleles.size() << " occurrences " << occurrencesPerAlleleSite.size() << " occurrenceAlleles " << occurrenceAlleles.size() << " informative sites " << occurrenceAlleles[0].size() << std::endl;
				countSplittedTo += chunkResult.size();
				countSplitted += 1;
				while (chunkResult.size() > 0)
				{
					chunksNeedProcessing.emplace_back();
					std::swap(chunksNeedProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
			}
		}
	});
	for (size_t i = 0; i < chunksDoneProcessing.size(); i++)
	{
		for (auto pair : chunksDoneProcessing[i])
		{
			std::get<2>(chunksPerRead[pair.first][pair.second]) = i + (std::get<2>(chunksPerRead[pair.first][pair.second]) & firstBitUint64_t);
		}
	}
	std::cerr << "allele phasing splitted " << countSplitted << " chunks to " << countSplittedTo << std::endl;
}

void splitPerPhasingKmersWithinChunk(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "splitting by phasing kmers" << std::endl;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksNeedProcessing;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunksDoneProcessing;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= chunksNeedProcessing.size()) chunksNeedProcessing.resize((std::get<2>(t) & maskUint64_t)+1);
			chunksNeedProcessing[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	// biggest on top so starts processing first
	std::sort(chunksNeedProcessing.begin(), chunksNeedProcessing.end(), [](const std::vector<std::pair<size_t, size_t>>& left, const std::vector<std::pair<size_t, size_t>>& right) { return left.size() < right.size(); });
	iterateMultithreaded(0, numThreads, numThreads, [&sequenceIndex, &chunksPerRead, &chunksNeedProcessing, &chunksDoneProcessing, &resultMutex, &countSplitted, &rawReadLengths, kmerSize](size_t dummy)
	{
		while (true)
		{
			std::vector<std::pair<size_t, size_t>> chunkBeingDone;
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				if (chunksNeedProcessing.size() == 0) return;
				std::swap(chunkBeingDone, chunksNeedProcessing.back());
				chunksNeedProcessing.pop_back();
			}
			auto startTime = getTime();
			if (chunkBeingDone.size() < 10)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				std::cerr << "skip chunk with coverage " << chunkBeingDone.size() << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			phmap::flat_hash_map<std::pair<size_t, size_t>, phmap::flat_hash_map<size_t, size_t>> fwForks; // (hom kmer, hom cluster) -> (occurrence ID -> seq char)
			phmap::flat_hash_map<std::pair<size_t, size_t>, phmap::flat_hash_map<size_t, size_t>> bwForks; // (hom kmer, hom cluster) -> (occurrence ID -> seq char)
			size_t lastOccurrenceID = std::numeric_limits<size_t>::max();
			size_t lastKmer = std::numeric_limits<size_t>::max();
			size_t lastCluster = std::numeric_limits<size_t>::max();
			size_t lastKmerPos = std::numeric_limits<size_t>::max();
			std::vector<TwobitString> chunkSequences;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, chunkBeingDone[j].first, chunkBeingDone[j].second));
			}
			auto validClusters = iterateSolidKmers(chunkSequences, kmerSize, chunkBeingDone.size() * 0.95, false, true, [&fwForks, &bwForks, &lastOccurrenceID, &lastKmer, &lastKmerPos, &lastCluster, &chunkBeingDone, &chunkSequences, &chunksPerRead, kmerSize](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				if (occurrenceID != lastOccurrenceID || pos != lastKmerPos+1)
				{
					if (lastOccurrenceID != std::numeric_limits<size_t>::max())
					{
						if (lastKmerPos + kmerSize < chunkSequences[lastOccurrenceID].size())
						{
							assert(fwForks.count(std::make_pair(lastKmer, lastCluster)) == 0 || fwForks.at(std::make_pair(lastKmer, lastCluster)).count(lastOccurrenceID) == 0);
							fwForks[std::make_pair(lastKmer, lastCluster)][lastOccurrenceID] = chunkSequences[lastOccurrenceID].get(lastKmerPos + kmerSize);
						}
					}
					if (pos > 0)
					{
						assert(bwForks.count(std::make_pair(kmer, clusterIndex)) == 0 || bwForks.at(std::make_pair(kmer, clusterIndex)).count(occurrenceID) == 0);
						bwForks[std::make_pair(kmer, clusterIndex)][occurrenceID] = chunkSequences[occurrenceID].get(pos - 1);
					}
				}
				lastOccurrenceID = occurrenceID;
				lastCluster = clusterIndex;
				lastKmerPos = pos;
				lastKmer = kmer;
			});
			if (lastOccurrenceID != std::numeric_limits<size_t>::max())
			{
				if (lastKmerPos + kmerSize < chunkSequences[lastOccurrenceID].size())
				{
					assert(fwForks.count(std::make_pair(lastKmer, lastCluster)) == 0 || fwForks.at(std::make_pair(lastKmer, lastCluster)).count(lastOccurrenceID) == 0);
					fwForks[std::make_pair(lastKmer, lastCluster)][lastOccurrenceID] = chunkSequences[lastOccurrenceID].get(lastKmerPos + kmerSize);
				}
			}
			size_t checkedPairs = 0;
			std::vector<std::vector<size_t>> phaseIdentities;
			phaseIdentities.resize(chunkBeingDone.size());
			size_t checkedFwBw = 0;
			size_t checkedBwBw = 0;
			size_t checkedFwFw = 0;
			for (const auto& bwFork : bwForks)
			{
				for (const auto& fwFork : fwForks)
				{
					if (validClusters.at(bwFork.first.first)[bwFork.first.second].second > validClusters.at(fwFork.first.first)[fwFork.first.second].first) continue;
					checkedPairs += 1;
					checkedFwBw += 1;
					checkPhasablePair(bwFork.second, fwFork.second, phaseIdentities);
				}
			}
			for (const auto& bwFork : bwForks)
			{
				for (const auto& bwFork2 : bwForks)
				{
					if (validClusters.at(bwFork.first.first)[bwFork.first.second].second+1 > validClusters.at(bwFork2.first.first)[bwFork2.first.second].first && validClusters.at(bwFork2.first.first)[bwFork2.first.second].second+1 > validClusters.at(bwFork.first.first)[bwFork.first.second].first) continue;
					checkedPairs += 1;
					checkedBwBw += 1;
					checkPhasablePair(bwFork.second, bwFork2.second, phaseIdentities);
				}
			}
			for (const auto& fwFork : fwForks)
			{
				for (const auto& fwFork2 : fwForks)
				{
					if (validClusters.at(fwFork.first.first)[fwFork.first.second].second+1 > validClusters.at(fwFork2.first.first)[fwFork2.first.second].first && validClusters.at(fwFork2.first.first)[fwFork2.first.second].second+1 > validClusters.at(fwFork.first.first)[fwFork.first.second].first) continue;
					checkedPairs += 1;
					checkedFwFw += 1;
					checkPhasablePair(fwFork.second, fwFork2.second, phaseIdentities);
				}
			}
			std::vector<size_t> parent;
			for (size_t j = 0; j < chunkBeingDone.size(); j++)
			{
				parent.emplace_back(j);
			}
			for (size_t j = 1; j < phaseIdentities.size(); j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					if (phaseIdentities[k] != phaseIdentities[j]) continue;
					merge(parent, k, j);
				}
			}
			phmap::flat_hash_map<size_t, size_t> keyToNode;
			size_t nextNum = 0;
			for (size_t j = 0; j < parent.size(); j++)
			{
				size_t key = find(parent, j);
				if (keyToNode.count(key) == 1) continue;
				keyToNode[key] = nextNum;
				nextNum += 1;
			}
			if (keyToNode.size() == 1)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				auto endTime = getTime();
				std::cerr << "phasing kmer splitted chunk with coverage " << chunkBeingDone.size() << " forkpairs " << phaseIdentities[0].size() << " checked " << checkedPairs << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				// std::cerr << "bwforks " << bwForks.size() << " fwforks " << fwForks.size() << std::endl;
				// std::cerr << "bwfw " << checkedFwBw << " bwbw " << checkedBwBw << " fwfw " << checkedFwFw << std::endl;
				chunksDoneProcessing.emplace_back();
				std::swap(chunksDoneProcessing.back(), chunkBeingDone);
				continue;
			}
			std::vector<std::vector<std::pair<size_t, size_t>>> chunkResult;
			chunkResult.resize(keyToNode.size());
			for (size_t j = 0; j < parent.size(); j++)
			{
				chunkResult[keyToNode.at(find(parent, j))].emplace_back(chunkBeingDone[j]);
			}
			// sort smallest last, so emplace-pop-swap puts biggest on top of chunksNeedProcessing
			std::sort(chunkResult.begin(), chunkResult.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				auto endTime = getTime();
				if (keyToNode.size() >= 2) countSplitted += 1;
				std::cerr << "phasing kmer splitted chunk with coverage " << chunkBeingDone.size() << " forkpairs " << phaseIdentities[0].size() << " checked " << checkedPairs << " to " << keyToNode.size() << " chunks time " << formatTime(startTime, endTime) << std::endl;
				// std::cerr << "bwforks " << bwForks.size() << " fwforks " << fwForks.size() << std::endl;
				// std::cerr << "bwfw " << checkedFwBw << " bwbw " << checkedBwBw << " fwfw " << checkedFwFw << std::endl;
				while (chunkResult.size() > 0)
				{
					chunksNeedProcessing.emplace_back();
					std::swap(chunksNeedProcessing.back(), chunkResult.back());
					chunkResult.pop_back();
				}
			}
		}
	});
	for (size_t i = 0; i < chunksDoneProcessing.size(); i++)
	{
		for (auto pair : chunksDoneProcessing[i])
		{
			std::get<2>(chunksPerRead[pair.first][pair.second]) = i + (std::get<2>(chunksPerRead[pair.first][pair.second]) & firstBitUint64_t);
		}
	}
	std::cerr << "phasing kmers splitted " << countSplitted << " chunks" << std::endl;
}

void splitPerSequenceIdentityRoughly(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting roughly by sequence identity" << std::endl;
	const size_t mismatchFloor = 10;
	size_t nextNum = 0;
	size_t totalSplitted = 0;
	size_t totalSplittedTo = 0;
	std::mutex resultMutex;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, &totalSplitted, &totalSplittedTo, mismatchFloor](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		if (occurrencesPerChunk[i].size() < 1000)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			// std::cerr << "skip chunk with coverage " << occurrencesPerChunk[i].size() << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		auto startTime = getTime();
		std::vector<std::string> sequences;
		sequences.reserve(occurrencesPerChunk[i].size());
		size_t longestLength = 0;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
			longestLength = std::max(longestLength, sequences.back().size());
		}
		std::vector<size_t> parent = getFastTransitiveClosure(sequences.size(), std::max(longestLength * mismatchFraction, mismatchFraction), [&sequences](const size_t i, const size_t j, const size_t maxDist) { return getNumMismatches(sequences[i], sequences[j], maxDist); });
		phmap::flat_hash_map<size_t, size_t> parentToCluster;
		for (size_t j = 0; j < parent.size(); j++)
		{
			find(parent, j);
		}
		size_t nextCluster = 0;
		for (size_t j = 0; j < parent.size(); j++)
		{
			if (parent[j] != j) continue;
			parentToCluster[j] = nextCluster;
			nextCluster += 1;
		}
		auto endTime = getTime();
		assert(nextCluster >= 1);
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			if (nextCluster > 1)
			{
				totalSplitted += 1;
				totalSplittedTo += nextCluster;
			}
			// std::cerr << "split chunk with coverage " << occurrencesPerChunk[i].size() << " into " << nextCluster << " chunks time " << formatTime(startTime, endTime) << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + parentToCluster.at(parent[j]) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += nextCluster;
		}
	});
	std::cerr << "rough sequence identity splitting splitted " << totalSplitted << " chunks to " << totalSplittedTo << " chunks" << std::endl;
}

void splitPerSequenceIdentity(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
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
	size_t firstSmall = iterationOrder.size();
	for (size_t i = 0; i < iterationOrder.size(); i++)
	{
		if (occurrencesPerChunk[iterationOrder[i]].size() > 200 || occurrencesPerChunk[iterationOrder[i]].size() > 200*numThreads) continue;
		firstSmall = i;
		break;
	}
	size_t nextNum = 0;
	for (size_t iterationIndex = 0; iterationIndex < firstSmall; iterationIndex++)
	{
		const size_t i = iterationOrder[iterationIndex];
		// std::cerr << "try split chunk " << i << std::endl;
		// auto startTime = getTime();
		std::vector<std::pair<std::string, size_t>> sequencesPerOccurrence;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sequencesPerOccurrence.emplace_back();
			sequencesPerOccurrence.back().second = j;
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequencesPerOccurrence.back().first = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
		}
		std::sort(sequencesPerOccurrence.begin(), sequencesPerOccurrence.end(), [](const auto& left, const auto& right) {
			if (left.first.size() < right.first.size()) return true;
			if (left.first.size() > right.first.size()) return false;
			return left.first < right.first;
		});
		std::vector<size_t> parent;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			parent.emplace_back(j);
		}
		const size_t blockSize = 50;
		const size_t numBlocks = (sequencesPerOccurrence.size()+blockSize-1)/blockSize;
		// std::cerr << "split part 1" << std::endl;
		iterateMultithreaded(0, numBlocks, numThreads, [&sequencesPerOccurrence, &parent, blockSize, i, mismatchFloor](size_t index)
		{
			for (size_t j = index * blockSize+1; j < (index+1) * blockSize && j < sequencesPerOccurrence.size(); j++)
			{
				for (size_t k = j-1; k < j && k >= index*blockSize; k--)
				{
					if (find(parent, k) == find(parent, j)) continue;
					assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
					size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
					if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches) break;
					size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
					if (mismatches > maxMismatches) continue;
					merge(parent, j, k);
				}
			}
		});
		for (size_t index = 1; index < numBlocks; index++)
		{
			size_t j = index * blockSize;
			size_t k = j-1;
			assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
			size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
			size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
			if (mismatches > maxMismatches) continue;
			merge(parent, j, k);
		}
		for (size_t i = 0; i < parent.size(); i++)
		{
			find(parent, i);
		}
		std::vector<size_t> activeBlocks;
		for (size_t i = 1; i < numBlocks; i++)
		{
			activeBlocks.emplace_back(i);
		}
		for (size_t blockDistance = 1; blockDistance < numBlocks; blockDistance++)
		{
			// std::cerr << "split part 2 dist " << blockDistance << std::endl;
			std::vector<std::pair<size_t, size_t>> merges;
			std::mutex mergeMutex;
			phmap::flat_hash_set<size_t> removeBlocks;
			iterateMultithreaded(0, activeBlocks.size(), numThreads, [&sequencesPerOccurrence, &parent, &merges, &mergeMutex, &activeBlocks, &removeBlocks, blockSize, blockDistance, i, mismatchFloor](size_t blockindex)
			{
				size_t index = activeBlocks[blockindex];
				phmap::flat_hash_map<size_t, size_t> extraParent;
				std::vector<std::pair<size_t, size_t>> mergesHere;
				for (size_t j = index * blockSize; j < (index+1) * blockSize && j < sequencesPerOccurrence.size(); j++)
				{
					for (size_t k = (index-blockDistance+1) * blockSize - 1; k < j && k >= (index-blockDistance) * blockSize; k--)
					{
						size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
						if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches)
						{
							if (j == index * blockSize && k == (index-blockDistance+1) * blockSize - 1)
							{
								assert(mergesHere.size() == 0);
								std::lock_guard<std::mutex> lock { mergeMutex };
								removeBlocks.emplace(index);
								return;
							}
							break;
						}
						if (parent[k] == parent[j]) continue;
						if (extraParent.count(parent[k]) == 0) extraParent[parent[k]] = parent[k];
						if (extraParent.count(parent[j]) == 0) extraParent[parent[j]] = parent[j];
						if (find(extraParent, parent[k]) == find(extraParent, parent[j])) continue;
						assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
						size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
						if (mismatches > maxMismatches) continue;
						merge(extraParent, parent[j], parent[k]);
						mergesHere.emplace_back(parent[j], parent[k]);
					}
				}
				{
					std::lock_guard<std::mutex> lock { mergeMutex };
					merges.insert(merges.end(), mergesHere.begin(), mergesHere.end());
				}
			});
			// std::cerr << merges.size() << " merges" << std::endl;
			for (auto pair : merges)
			{
				merge(parent, pair.first, pair.second);
			}
			if (merges.size() >= 1)
			{
				for (size_t i = 0; i < parent.size(); i++)
				{
					find(parent, i);
				}
			}
			removeBlocks.emplace(blockDistance);
			for (size_t i = activeBlocks.size()-1; i < activeBlocks.size(); i--)
			{
				if (removeBlocks.count(activeBlocks[i]) == 0) continue;
				std::swap(activeBlocks[i], activeBlocks.back());
				activeBlocks.pop_back();
			}
			// std::cerr << "remove " << removeBlocks.size() << ", " << activeBlocks.size() << " remaining" << std::endl;
			if (activeBlocks.size() == 0) break;
		}
		// auto endTime = getTime();
		phmap::flat_hash_map<size_t, size_t> keyToNode;
		for (size_t j = 0; j < parent.size(); j++)
		{
			size_t key = find(parent, j);
			if (keyToNode.count(key) == 1) continue;
			keyToNode[key] = nextNum;
			nextNum += 1;
		}
		// std::cerr << "sequence identity big chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " time " << formatTime(startTime, endTime) << " splitted to " << keyToNode.size() << std::endl;
		// std::cerr << "sequence identity splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			size_t index = sequencesPerOccurrence[j].second;
			std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][index].first][occurrencesPerChunk[i][index].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
		}
	}
	std::mutex resultMutex;
	iterateMultithreaded(firstSmall, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &sequenceIndex, &iterationOrder, &rawReadLengths, mismatchFloor](const size_t iterationIndex)
	{
		const size_t i = iterationOrder[iterationIndex];
		// {
		// 	std::lock_guard<std::mutex> lock { resultMutex };
		// 	std::cerr << "sequence identity split chunk " << i << " coverage " << occurrencesPerChunk[i].size() << std::endl;
		// }
		// size_t coverage = occurrencesPerChunk[i].size();
		// size_t totalCheckablePairs = 0;
		// size_t totalMergedPairs = 0;
		// size_t totalNotMergedPairs = 0;
		// size_t totalDistinct = 0;
		// auto startTime = getTime();
		std::vector<std::pair<std::string, size_t>> sequencesPerOccurrence;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sequencesPerOccurrence.emplace_back();
			sequencesPerOccurrence.back().second = j;
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			sequencesPerOccurrence.back().first = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
		}
		std::sort(sequencesPerOccurrence.begin(), sequencesPerOccurrence.end(), [](const auto& left, const auto& right) {
			if (left.first.size() < right.first.size()) return true;
			if (left.first.size() > right.first.size()) return false;
			return left.first < right.first;
		});
		std::vector<bool> differentFromPredecessor;
		differentFromPredecessor.resize(occurrencesPerChunk[i].size(), false);
		differentFromPredecessor[0] = true;
		// totalDistinct += 1;
		for (size_t j = 1; j < differentFromPredecessor.size(); j++)
		{
			if (sequencesPerOccurrence[j].first == sequencesPerOccurrence[j-1].first) continue;
			differentFromPredecessor[j] = true;
			// totalDistinct += 1;
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
				// totalCheckablePairs += 1;
				if (find(parent, k) == find(parent, j)) continue;
				assert(sequencesPerOccurrence[j].first.size() >= sequencesPerOccurrence[k].first.size());
				size_t maxMismatches = std::max(mismatchFloor, (size_t)(sequencesPerOccurrence[k].first.size()*mismatchFraction));
				if (sequencesPerOccurrence[j].first.size() - sequencesPerOccurrence[k].first.size() >= maxMismatches) break;
				size_t mismatches = getNumMismatches(sequencesPerOccurrence[j].first, sequencesPerOccurrence[k].first, maxMismatches);
				if (mismatches > maxMismatches)
				{
					// totalNotMergedPairs += 1;
					continue;
				}
				// totalMergedPairs += 1;
				merge(parent, j, k);
			}
		}
		// auto endTime = getTime();
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
			// std::cerr << "sequence identity chunk " << i << " coverage " << coverage << " distinct " << totalDistinct << " checkable pairs " << totalCheckablePairs << " merged " << totalMergedPairs << " not merged " << totalNotMergedPairs << " time " << formatTime(startTime, endTime) << " splitted to " << keyToNode.size() << std::endl;
			// std::cerr << "sequence identity splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
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

bool hasHomopolymerWithLengthThreeOrMore(uint64_t kmer, const size_t kmerSize)
{
	kmer ^= kmer >> 2;
	kmer = (~(kmer | (kmer >> 1))) & 0x5555555555555555ull;
	kmer &= kmer >> 2;
	kmer &= (1ull << (2ull * (kmerSize-2)))-1;
	return kmer != 0;
}

bool hasMicrosatelliteMotifLengthTwo(uint64_t kmer, const size_t kmerSize)
{
	kmer ^= kmer >> 4;
	kmer = (~(kmer | (kmer >> 1))) & 0x5555555555555555ull;
	kmer &= kmer >> 2;
	kmer &= (1ull << (2ull * (kmerSize-3)))-1;
	return kmer != 0;
}

uint64_t weightedHash(uint64_t kmer, const size_t kmerSize)
{
	const size_t hashmask = 0x1FFFFFFFFFFFFFFFull;
	uint64_t result = hash(kmer) & hashmask;
	if (hasHomopolymerWithLengthThreeOrMore(kmer, kmerSize)) result |= 0x8000000000000000ull;
	if (hasMicrosatelliteMotifLengthTwo(kmer, kmerSize)) result |= 0x4000000000000000ull;
	return result;
}

uint64_t reverseKmer(uint64_t kmer, const size_t kmerSize)
{
	kmer = ~kmer;
	kmer = ((kmer >> 2) & 0x3333333333333333ull) + ((kmer << 2) & 0xCCCCCCCCCCCCCCCCull);
	kmer = ((kmer >> 4) & 0x0F0F0F0F0F0F0F0Full) + ((kmer << 4) & 0xF0F0F0F0F0F0F0F0ull);
	kmer = ((kmer >> 8) & 0x00FF00FF00FF00FFull) + ((kmer << 8) & 0xFF00FF00FF00FF00ull);
	kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFull) + ((kmer << 16) & 0xFFFF0000FFFF0000ull);
	kmer = ((kmer >> 32) & 0x00000000FFFFFFFFull) + ((kmer << 32) & 0xFFFFFFFF00000000ull);
	kmer >>= (32ull-kmerSize)*2ull;
	return kmer;
}

uint64_t hashFwAndBw(const uint64_t kmer, const size_t kmerSize)
{
	uint64_t fwHash = weightedHash(kmer, kmerSize);
	uint64_t bwHash = weightedHash(reverseKmer(kmer, kmerSize), kmerSize);
	return std::min(fwHash, bwHash);
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
	queue.emplace(hashFwAndBw(kmer, k));
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
		size_t h = hashFwAndBw(kmer, k);
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

void splitPerFirstLastKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	assert(kmerSize <= 31);
	size_t nextNum = 0;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> endKmersToNumber;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		std::string readseq = getReadSequence(sequenceIndex, i);
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<0>(chunksPerRead[i][j]) < std::get<1>(chunksPerRead[i][j]));
			assert(std::get<1>(chunksPerRead[i][j]) < readseq.size());
			uint64_t firstKmer = 0;
			uint64_t lastKmer = 0;
			for (size_t k = 0; k < kmerSize; k++)
			{
				firstKmer <<= 2;
				lastKmer <<= 2;
				switch(readseq[std::get<0>(chunksPerRead[i][j])+k])
				{
					case 'A':
						firstKmer += 0;
						break;
					case 'C':
						firstKmer += 1;
						break;
					case 'G':
						firstKmer += 2;
						break;
					case 'T':
						firstKmer += 3;
						break;
					default:
						assert(false);
				}
				switch(readseq[std::get<1>(chunksPerRead[i][j])-k])
				{
					case 'A':
						lastKmer += 3;
						break;
					case 'C':
						lastKmer += 2;
						break;
					case 'G':
						lastKmer += 1;
						break;
					case 'T':
						lastKmer += 0;
						break;
					default:
						assert(false);
				}
			}
			auto key = std::make_pair(firstKmer, lastKmer);
			if (endKmersToNumber.count(key) == 0)
			{
				endKmersToNumber[key] = nextNum;
				nextNum += 1;
			}
			std::get<2>(chunksPerRead[i][j]) = endKmersToNumber.at(key) + firstBitUint64_t;
		}
	}
}

void splitPerBaseCounts(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t mismatchFloor = 10;
	std::cerr << "splitting by base counts" << std::endl;
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, mismatchFloor](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		std::vector<std::vector<size_t>> countsPerOccurrence;
		countsPerOccurrence.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			std::string chunkSequence = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
			assert(chunkSequence.size() == std::get<1>(t)-std::get<0>(t)+1);
			countsPerOccurrence[j].resize(4);
			for (size_t k = 0; k < chunkSequence.size(); k++)
			{
				switch(chunkSequence[k])
				{
				case 'A':
					countsPerOccurrence[j][0] += 1;
					break;
				case 'C':
					countsPerOccurrence[j][1] += 1;
					break;
				case 'G':
					countsPerOccurrence[j][2] += 1;
					break;
				case 'T':
					countsPerOccurrence[j][3] += 1;
					break;
				default:
					assert(false);
				}
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

void splitPerInterchunkPhasedKmers(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize)
{
	std::cerr << "splitting by interchunk phasing kmers" << std::endl;
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
	size_t countSplitted = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &sequenceIndex, &repetitive, &maybePhaseGroups, &countSplitted, &rawReadLengths, kmerSize](const size_t i)
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
			if (firsts+seconds < readsHere.size() * 0.75 && (firsts < 5 || seconds < 5)) continue;
			if (firsts+seconds < readsHere.size() * 0.25) continue;
			applicablePhasingGroups.push_back(j);
		}
		if (applicablePhasingGroups.size() == 0)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		phmap::flat_hash_map<size_t, std::vector<size_t>> phaseGroupsPerRead;
		size_t nextGroupNum = 0;
		std::vector<phmap::flat_hash_set<size_t>> solidKmersPerOccurrence;
		solidKmersPerOccurrence.resize(occurrencesPerChunk[i].size());
		phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerClusterToNumber;
		std::vector<TwobitString> chunkSequences;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
		}
		iterateSolidKmers(chunkSequences, kmerSize, occurrencesPerChunk[i].size(), true, false, [&solidKmersPerOccurrence, &kmerClusterToNumber](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
		{
			size_t kmerkey = 0;
			std::pair<size_t, size_t> key { kmer, clusterIndex };
			if (kmerClusterToNumber.count(key) == 0)
			{
				kmerkey = kmerClusterToNumber.size();
				kmerClusterToNumber[key] = kmerkey;
			}
			else
			{
				kmerkey = kmerClusterToNumber.at(key);
			}
			solidKmersPerOccurrence[occurrenceID].emplace(kmerkey);
		});
		phmap::flat_hash_map<size_t, size_t> kmerToNumber;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			phmap::flat_hash_set<size_t> kmersHere;
			iterateKmers(chunkSequences[j], 0, chunkSequences[j].size()-1, true, kmerSize, [&kmersHere](const size_t kmer, const size_t pos)
			{
				kmersHere.insert(kmer);
			});
			for (auto kmer : kmersHere)
			{
				size_t kmerkey = 0;
				if (kmerToNumber.count(kmer) == 0)
				{
					kmerkey = kmerToNumber.size() + kmerClusterToNumber.size();
					kmerToNumber[kmer] = kmerkey;
				}
				else
				{
					kmerkey = kmerToNumber.at(kmer);
				}
				solidKmersPerOccurrence[j].emplace(kmerkey);
			}
		}
		{
			std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>> alleles; // startkmer, startcluster, endkmer, endcluster, occurrenceID, allele
			size_t lastOccurrence = std::numeric_limits<size_t>::max();
			size_t lastKmer = std::numeric_limits<size_t>::max();
			size_t lastCluster = std::numeric_limits<size_t>::max();
			size_t lastKmerPos = std::numeric_limits<size_t>::max();
			phmap::flat_hash_map<std::string, size_t> alleleNumbers;
			iterateSolidKmers(chunkSequences, kmerSize, occurrencesPerChunk[i].size(), false, false, [&occurrencesPerChunk, &chunkSequences, &alleles, &lastOccurrence, &lastKmer, &lastCluster, &lastKmerPos, &alleleNumbers, kmerSize, i](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				if (occurrenceID != lastOccurrence || lastKmer == std::numeric_limits<size_t>::max())
				{
					lastOccurrence = occurrenceID;
					lastKmer = kmer;
					lastCluster = clusterIndex;
					lastKmerPos = pos;
					return;
				}
				size_t allele;
				allele = getAllele(chunkSequences[occurrenceID], lastKmerPos, pos + kmerSize, alleleNumbers, true);
				alleles.emplace_back(lastKmer, lastCluster, kmer, clusterIndex, occurrenceID, allele);
				lastOccurrence = occurrenceID;
				lastKmer = kmer;
				lastCluster = clusterIndex;
				lastKmerPos = pos;
			});
			std::sort(alleles.begin(), alleles.end());
			std::vector<std::vector<std::vector<size_t>>> occurrencesPerAlleleSite;
			occurrencesPerAlleleSite = getAlleleOccurrences(alleles, occurrencesPerChunk[i].size());
			size_t nextNum = kmerToNumber.size() + kmerClusterToNumber.size();
			for (size_t j = 0; j < occurrencesPerAlleleSite.size(); j++)
			{
				if (occurrencesPerAlleleSite[j].size() < 2) continue;
				for (size_t allele = 0; allele < occurrencesPerAlleleSite[j].size(); allele++)
				{
					for (size_t k : occurrencesPerAlleleSite[j][allele])
					{
						solidKmersPerOccurrence[k].emplace(nextNum);
					}
					nextNum += 1;
				}
			}
		}
		{
			phmap::flat_hash_map<size_t, size_t> kmerCoverage;
			for (size_t j = 0; j < solidKmersPerOccurrence.size(); j++)
			{
				for (auto kmer : solidKmersPerOccurrence[j])
				{
					kmerCoverage[kmer] += 1;
				}
			}
			phmap::flat_hash_set<size_t> removeKmers;
			for (auto pair : kmerCoverage)
			{
				if (pair.second < 3 || pair.second + 3 > occurrencesPerChunk[i].size())
				{
					removeKmers.insert(pair.first);
				}
			}
			if (removeKmers.size() >= 1)
			{
				for (size_t j = 0; j < solidKmersPerOccurrence.size(); j++)
				{
					for (auto kmer : removeKmers)
					{
						if (solidKmersPerOccurrence[j].count(kmer) == 1)
						{
							solidKmersPerOccurrence[j].erase(kmer);
						}
					}
				}
			}
		}
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
				const phmap::flat_hash_set<size_t> kmersHere = solidKmersPerOccurrence[j];
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
			assert(hasFirst);
			assert(hasSecond);
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
				const phmap::flat_hash_set<size_t>& kmersHere = solidKmersPerOccurrence[j];
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
				const phmap::flat_hash_set<size_t>& kmersHere = solidKmersPerOccurrence[j];
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
			// size_t first = nextNum;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				if (clusterToNode.count(find(parent, j)) == 1) continue;
				clusterToNode[find(parent, j)] = nextNum;
				nextNum += 1;
			}
			// size_t last = nextNum-1;
			if (clusterToNode.size() >= 2) countSplitted += 1;
//			std::cerr << "interchunk phasing kmers splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks range " << first << " - " << last << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, j)) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
	std::cerr << "interchunk phasing kmers splitted " << countSplitted << " chunks" << std::endl;
}

void splitPerMinHashes(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t kmerSize = 11;
	const size_t hashCount = 10;
	std::cerr << "splitting by minhash" << std::endl;
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &sequenceIndex, &rawReadLengths, kmerSize, hashCount](const size_t i, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		phmap::flat_hash_map<size_t, size_t> parent;
		std::vector<size_t> oneHashPerLocation;
		size_t smallestSequence = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			smallestSequence = std::min(smallestSequence, std::get<1>(t)-std::get<0>(t)+1);
		}
		if (smallestSequence <= kmerSize + hashCount + 1)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			std::vector<uint64_t> minHashes;
			auto chunkSequence = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second);
			minHashes = getMinHashes(chunkSequence, kmerSize, hashCount);
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
			if (NonexistantChunk(std::get<2>(t)))
			{
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << "-" << std::endl;
				continue;
			}
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
			assert(length >= 1);
			assert(length <= 1000000);
			assert(std::get<0>(chunkLocationInUnitig[unitigs[i][j] & maskUint64_t]) != std::numeric_limits<size_t>::max());
			assert((std::get<0>(chunkLocationInUnitig[unitigs[i][j] & maskUint64_t]) & maskUint64_t) == i);
			assert(std::get<2>(chunkLocationInUnitig[unitigs[i][j] & maskUint64_t]) == 0);
			assert(std::get<3>(chunkLocationInUnitig[unitigs[i][j] & maskUint64_t]) == 0);
			if (j > 0)
			{
				std::pair<size_t, bool> fromnode { unitigs[i][j-1] & maskUint64_t, unitigs[i][j-1] & firstBitUint64_t };
				std::pair<size_t, bool> tonode { unitigs[i][j] & maskUint64_t, unitigs[i][j] & firstBitUint64_t };
				auto pairkey = MBG::canon(fromnode, tonode);
				std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
				size_t overlap = edgeOverlaps.at(key);
				size_t prevLength = lengths[unitigs[i][j-1] & maskUint64_t][lengths[unitigs[i][j-1] & maskUint64_t].size()/2];
				if (overlap >= prevLength) overlap = prevLength-1; // wrong but do it anyway to prevent crash later
				if (overlap >= length) overlap = length-1; // wrong but do it anyway to prevent crash later
				assert(unitigLengths.back() > overlap);
				unitigLengths.back() -= overlap;
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
		assert(std::get<3>(chunkLocationInUnitig[i]) - std::get<2>(chunkLocationInUnitig[i]) == lengths[i][lengths[i].size()/2]);
		assert(std::get<3>(chunkLocationInUnitig[i]) <= unitigLengths[unitig]);
		if (std::get<0>(chunkLocationInUnitig[i]) & firstBitUint64_t) continue;
		std::get<1>(chunkLocationInUnitig[i]) = unitigs[unitig].size() - 1 - std::get<1>(chunkLocationInUnitig[i]);
		std::swap(std::get<2>(chunkLocationInUnitig[i]), std::get<3>(chunkLocationInUnitig[i]));
		std::get<2>(chunkLocationInUnitig[i]) = unitigLengths[unitig] - std::get<2>(chunkLocationInUnitig[i]);
		std::get<3>(chunkLocationInUnitig[i]) = unitigLengths[unitig] - std::get<3>(chunkLocationInUnitig[i]);
		assert(std::get<2>(chunkLocationInUnitig[i]) < std::get<3>(chunkLocationInUnitig[i]));
		assert(std::get<3>(chunkLocationInUnitig[i]) - std::get<2>(chunkLocationInUnitig[i]) == lengths[i][lengths[i].size()/2]);
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
		std::vector<std::pair<size_t, size_t>> readPartInPathnode;
		uint64_t currentUnitig = std::numeric_limits<uint64_t>::max();
		size_t currentUnitigIndex = std::numeric_limits<size_t>::max();
		size_t currentPathLeftClip = std::numeric_limits<size_t>::max();
		size_t currentPathRightClip = std::numeric_limits<size_t>::max();
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
					result[i].back().readPartInPathnode = readPartInPathnode;
				}
				path.clear();
				readPartInPathnode.clear();
				currentUnitig = std::numeric_limits<size_t>::max();
				currentUnitigIndex = 0;
				currentPathLeftClip= 0;
				currentPathRightClip = 0;
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
				readPartInPathnode.emplace_back(readStart, readEnd);
				continue;
			}
			if (currentUnitig == unitig && currentUnitigIndex+1 == unitigIndex)
			{
				assert(readEnd > readPartInPathnode.back().second);
				assert(graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos < currentPathRightClip);
				currentUnitigIndex = unitigIndex;
				readPartInPathnode.back().second = readEnd;
				currentPathRightClip = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
				continue;
			}
			if (currentUnitigIndex+1 == unitigs[currentUnitig & maskUint64_t].size() && unitigIndex == 0)
			{
				std::pair<size_t, bool> fromUnitig { currentUnitig & maskUint64_t, currentUnitig & firstBitUint64_t };
				std::pair<size_t, bool> thisUnitig { unitig & maskUint64_t, unitig & firstBitUint64_t };
				if (graph.edges.hasEdge(fromUnitig, thisUnitig))
				{
					assert(readEnd > readPartInPathnode.back().second);
					path.emplace_back(unitig);
					readPartInPathnode.emplace_back(readStart, readEnd);
					currentUnitig = unitig;
					currentUnitigIndex = unitigIndex;
					currentPathRightClip = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
					continue;
				}
			}
			result[i].emplace_back();
			result[i].back().path = path;
			result[i].back().pathLeftClip = currentPathLeftClip;
			result[i].back().pathRightClip = currentPathRightClip;
			result[i].back().readPartInPathnode = readPartInPathnode;
			path.clear();
			readPartInPathnode.clear();
			path.emplace_back(unitig);
			currentUnitig = unitig;
			currentUnitigIndex = unitigIndex;
			currentPathLeftClip = unitigStartPos;
			currentPathRightClip = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
			readPartInPathnode.emplace_back(readStart, readEnd);
		}
		if (path.size() >= 1)
		{
			result[i].emplace_back();
			result[i].back().path = path;
			result[i].back().pathLeftClip = currentPathLeftClip;
			result[i].back().pathRightClip = currentPathRightClip;
			result[i].back().readPartInPathnode = readPartInPathnode;
		}
	}
	return result;
}

void removeKmer(std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& kmersPerString, std::pair<size_t, size_t> removeThis)
{
	assert(removeThis.first != std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < kmersPerString.size(); i++)
	{
		assert(kmersPerString[i].size() >= 2);
		assert(std::get<0>(kmersPerString[i][0]) != removeThis.first || std::get<1>(kmersPerString[i][0]) != removeThis.second);
		assert(std::get<0>(kmersPerString[i].back()) != removeThis.first || std::get<1>(kmersPerString[i].back()) != removeThis.second);
		for (size_t j = kmersPerString[i].size()-2; j > 0; j--)
		{
			if (std::get<0>(kmersPerString[i][j]) != removeThis.first || std::get<1>(kmersPerString[i][j]) != removeThis.second) continue;
			kmersPerString[i].erase(kmersPerString[i].begin()+j);
		}
	}
}

std::vector<std::pair<size_t, size_t>> getConsensusSolidKmerPath(std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& kmersPerString)
{
start:
	phmap::flat_hash_map<std::pair<size_t, size_t>, phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>> outEdges;
	for (size_t i = 0; i < kmersPerString.size(); i++)
	{
		for (size_t j = 1; j < kmersPerString[i].size(); j++)
		{
			assert(std::get<2>(kmersPerString[i][j]) > std::get<2>(kmersPerString[i][j-1]));
			std::pair<size_t, size_t> from { std::get<0>(kmersPerString[i][j-1]), std::get<1>(kmersPerString[i][j-1]) };
			std::pair<size_t, size_t> to { std::get<0>(kmersPerString[i][j]), std::get<1>(kmersPerString[i][j]) };
			outEdges[from][to] += 1;
		}
	}
	std::vector<std::pair<size_t, size_t>> greedyChosenPath;
	greedyChosenPath.emplace_back(std::numeric_limits<size_t>::max(), 0);
	assert(outEdges.count(greedyChosenPath.back()) == 1);
	assert(outEdges.at(greedyChosenPath.back()).size() >= 1);
	phmap::flat_hash_set<std::pair<size_t, size_t>> visited;
	visited.emplace(std::numeric_limits<size_t>::max(), 0);
	while (outEdges.count(greedyChosenPath.back()) == 1)
	{
		assert(greedyChosenPath.size() <= outEdges.size()+5);
		size_t maxcov = 0;
		for (auto edge : outEdges.at(greedyChosenPath.back()))
		{
			maxcov = std::max(maxcov, edge.second);
		}
		assert(maxcov >= 1);
		for (auto edge : outEdges.at(greedyChosenPath.back()))
		{
			if (edge.second == maxcov)
			{
				if (visited.count(edge.first) == 1)
				{
					removeKmer(kmersPerString, edge.first);
					goto start; // manual tail call recursion optimization because apparently gcc doesn't
				}
				visited.emplace(edge.first);
				greedyChosenPath.emplace_back(edge.first);
				break;
			}
		}
	}
	assert(greedyChosenPath.size() >= 2);
	assert(greedyChosenPath.back().first == std::numeric_limits<size_t>::max());
	assert(greedyChosenPath.back().second == 1);
	return greedyChosenPath;
}

std::string pickMostCentralString(const std::vector<std::string>& options, const phmap::flat_hash_map<std::string, size_t>& pieceCounts)
{
	std::vector<size_t> editSum;
	editSum.resize(options.size(), 0);
	for (size_t i = 0; i < options.size(); i++)
	{
		for (const auto& pair : pieceCounts)
		{
			size_t edits = getNumMismatches(options[i], pair.first, std::max(options[i].size(), pair.first.size()));
			editSum[i] += edits * pair.second;
		}
	}
	size_t minIndex = 0;
	for (size_t i = 1; i < editSum.size(); i++)
	{
		if (editSum[i] < editSum[minIndex]) minIndex = i;
	}
	return options[minIndex];
}

std::string getConsensusFromSolidKmers(const phmap::flat_hash_map<std::string, size_t>& sequenceCount, const size_t totalCount)
{
	const size_t kmerSize = 11;
	assert(sequenceCount.size() >= 2);
	std::vector<TwobitString> strings;
	for (auto pair : sequenceCount)
	{
		for (size_t i = 0; i < pair.second; i++)
		{
			strings.emplace_back(pair.first);
		}
	}
	assert(totalCount >= 3);
	assert(strings.size() == totalCount);
	assert(strings.size() >= sequenceCount.size());
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> kmersPerString;
	kmersPerString.resize(strings.size());
	for (size_t i = 0; i < kmersPerString.size(); i++)
	{
		kmersPerString[i].emplace_back(std::numeric_limits<size_t>::max(), 0, 0);
	}
	iterateSolidKmers(strings, kmerSize, (strings.size()+1)*0.5, false, true, [&kmersPerString](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
	{
		if (pos == 0) return;
		assert(pos > std::get<2>(kmersPerString[occurrenceID].back()));
		kmersPerString[occurrenceID].emplace_back(kmer, clusterIndex, pos);
	});
	for (size_t i = 0; i < kmersPerString.size(); i++)
	{
		assert(std::get<2>(kmersPerString[i].back()) < strings[i].size());
		kmersPerString[i].emplace_back(std::numeric_limits<size_t>::max(), 1, strings[i].size());
	}
	std::vector<std::pair<size_t, size_t>> chosenPath = getConsensusSolidKmerPath(kmersPerString);
	assert(chosenPath.size() >= 2);
	assert(chosenPath[0].first == std::numeric_limits<size_t>::max());
	assert(chosenPath[0].second == 0);
	assert(chosenPath.back().first == std::numeric_limits<size_t>::max());
	assert(chosenPath.back().second == 1);
	std::vector<phmap::flat_hash_map<std::string, size_t>> pieceCounts;
	pieceCounts.resize(chosenPath.size()-1);
	phmap::flat_hash_map<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>, size_t> pieceLocation;
	for (size_t i = 1; i < chosenPath.size(); i++)
	{
		pieceLocation[std::make_pair(chosenPath[i-1], chosenPath[i])] = i-1;
	}
	for (size_t i = 0; i < kmersPerString.size(); i++)
	{
		for (size_t j = 1; j < kmersPerString[i].size(); j++)
		{
			std::pair<size_t, size_t> from { std::get<0>(kmersPerString[i][j-1]), std::get<1>(kmersPerString[i][j-1]) };
			std::pair<size_t, size_t> to { std::get<0>(kmersPerString[i][j]), std::get<1>(kmersPerString[i][j]) };
			if (pieceLocation.count(std::make_pair(from, to)) == 1)
			{
				std::string seq = strings[i].substr(std::get<2>(kmersPerString[i][j-1]), std::get<2>(kmersPerString[i][j])-std::get<2>(kmersPerString[i][j-1]));
				pieceCounts[pieceLocation.at(std::make_pair(from, to))][seq] += 1;
			}
		}
	}
	std::string result;
	for (size_t i = 0; i < pieceCounts.size(); i++)
	{
		size_t maxCount = 0;
		for (const auto& pair : pieceCounts[i])
		{
			maxCount = std::max(maxCount, pair.second);
		}
		assert(maxCount >= 1);
		std::vector<std::string> stringsWithMaxCount;
		for (const auto& pair : pieceCounts[i])
		{
			if (pair.second == maxCount) stringsWithMaxCount.emplace_back(pair.first);
		}
		assert(stringsWithMaxCount.size() >= 1);
		if (stringsWithMaxCount.size() == 1)
		{
			result += stringsWithMaxCount[0];
		}
		else
		{
			result += pickMostCentralString(stringsWithMaxCount, pieceCounts[i]);
		}
	}
	return result;
}

std::string getConsensusPickArbitrary(const phmap::flat_hash_map<std::string, size_t>& sequenceCount, const size_t totalCount)
{
	size_t maxCount = 0;
	for (const auto& pair : sequenceCount)
	{
		maxCount = std::max(pair.second, maxCount);
	}
	for (const auto& pair : sequenceCount)
	{
		if (pair.second == maxCount) return pair.first;
	}
	assert(false);
}

std::string getConsensus(const phmap::flat_hash_map<std::string, size_t>& sequenceCount, const size_t totalCount)
{
	std::vector<size_t> lengths;
	for (const auto& pair : sequenceCount)
	{
		for (size_t i = 0; i < pair.second; i++)
		{
			lengths.emplace_back(pair.first.size());
		}
	}
	std::sort(lengths.begin(), lengths.end());
	if (lengths[lengths.size()/2] < 20)
	{
		return getConsensusPickArbitrary(sequenceCount, totalCount);
	}
	return getConsensusFromSolidKmers(sequenceCount, totalCount);
}

std::string getConsensusSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const phmap::flat_hash_map<size_t, std::vector<std::tuple<size_t, bool, size_t>>>& readAnchorPoints, const size_t start, const size_t end)
{
	phmap::flat_hash_map<std::string, size_t> sequenceCount;
	size_t totalCount = 0;
	for (auto t1 : readAnchorPoints.at(start))
	{
		for (auto t2 : readAnchorPoints.at(end))
		{
			if (std::get<0>(t1) != std::get<0>(t2)) continue;
			if (std::get<1>(t1) != std::get<1>(t2)) continue;
			std::string sequenceHere;
			if (std::get<1>(t1))
			{
				if (std::get<2>(t1) >= std::get<2>(t2)) continue;
				if (std::get<2>(t2)-std::get<2>(t1) > end - start + 50) continue;
				if (std::get<2>(t2)-std::get<2>(t1) + 50 < end - start) continue;
				sequenceHere = sequenceIndex.getSubstring(std::get<0>(t1), std::get<2>(t1), std::get<2>(t2)-std::get<2>(t1)+1);
			}
			else
			{
				if (std::get<2>(t1) <= std::get<2>(t2)) continue;
				if (std::get<2>(t1)-std::get<2>(t2) > end - start + 50) continue;
				if (std::get<2>(t1)-std::get<2>(t2) + 50 < end - start) continue;
				sequenceHere = sequenceIndex.getSubstring(std::get<0>(t1), std::get<2>(t2), std::get<2>(t1)-std::get<2>(t2)+1);
				sequenceHere = MBG::revCompRaw(sequenceHere);
			}
			sequenceCount[sequenceHere] += 1;
			totalCount += 1;
		}
	}
	if (sequenceCount.size() == 0)
	{
		return "";
	}
	size_t maxCount = 0;
	for (const auto& pair : sequenceCount)
	{
		maxCount = std::max(maxCount, pair.second);
	}
	if (maxCount*2 > totalCount || totalCount == 2)
	{
		for (const auto& pair : sequenceCount)
		{
			if (pair.second == maxCount) return pair.first;
		}
	}
	return getConsensus(sequenceCount, totalCount);
}

std::vector<ConsensusString> getUnitigConsensuses(const ChunkUnitigGraph& unitigGraph, const std::vector<std::vector<UnitigPath>>& readPaths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<size_t>>& chunkLengths, const std::vector<size_t>& chunkCoverages, const std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>& chunkLocationInUnitig, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps, const size_t numThreads)
{
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "getting unitig consensuses" << std::endl;
	std::vector<phmap::flat_hash_map<size_t, std::vector<std::tuple<size_t, bool, size_t>>>> readAnchorPoints; // unitigpos -> read, fw, readpos
	readAnchorPoints.resize(unitigGraph.unitigLengths.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			size_t readStart = std::get<0>(t);
			size_t readEnd = std::get<1>(t);
			uint64_t node = std::get<2>(t);
			assert((node & maskUint64_t) < chunkLocationInUnitig.size());
			uint64_t unitig = std::get<0>(chunkLocationInUnitig[node & maskUint64_t]);
			if (unitig == std::numeric_limits<size_t>::max()) continue;
			size_t unitigStartPos = std::get<2>(chunkLocationInUnitig[node & maskUint64_t]);
			size_t unitigEndPos = std::get<3>(chunkLocationInUnitig[node & maskUint64_t]);
			bool fw = true;
			if (unitig & firstBitUint64_t)
			{
			}
			else
			{
				std::swap(unitigStartPos, unitigEndPos);
				unitigStartPos = unitigGraph.unitigLengths[unitig & maskUint64_t] - unitigStartPos;
				unitigEndPos = unitigGraph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
				node ^= firstBitUint64_t;
			}
			unitig = unitig & maskUint64_t;
			if (node & firstBitUint64_t)
			{
			}
			else
			{
				std::swap(unitigStartPos, unitigEndPos);
				fw = !fw;
			}
			assert(unitig < readAnchorPoints.size());
			readAnchorPoints[unitig][unitigStartPos].emplace_back(i, fw, readStart);
			readAnchorPoints[unitig][unitigEndPos].emplace_back(i, fw, readEnd);
		}
	}
	std::vector<ConsensusString> result;
	result.resize(unitigGraph.unitigLengths.size());
	iterateMultithreaded(0, unitigGraph.unitigLengths.size(), numThreads, [&readAnchorPoints, &sequenceIndex, &result, &unitigGraph](const size_t unitigi)
	{
		std::vector<size_t> anchorPositions;
		for (const auto& pair : readAnchorPoints[unitigi])
		{
			anchorPositions.emplace_back(pair.first);
		}
		assert(anchorPositions.size() >= 2);
		std::sort(anchorPositions.begin(), anchorPositions.end());
		assert(anchorPositions[0] == 0);
		assert(anchorPositions.back() == unitigGraph.unitigLengths[unitigi]);
		for (size_t i = 1; i < anchorPositions.size(); i++)
		{
			std::string seq = getConsensusSequence(sequenceIndex, readAnchorPoints[unitigi], anchorPositions[i-1], anchorPositions[i]);
			if (seq.size() >= 1)
			{
				if (i+1 < anchorPositions.size()) seq.pop_back();
				for (size_t j = 0; j < seq.size(); j++)
				{
					switch(seq[j])
					{
					case 'A':
						result[unitigi].bases.emplace_back(0);
						break;
					case 'C':
						result[unitigi].bases.emplace_back(1);
						break;
					case 'G':
						result[unitigi].bases.emplace_back(2);
						break;
					case 'T':
						result[unitigi].bases.emplace_back(3);
						break;
					default:
						assert(false);
						break;
					}
				}
			}
			else
			{
				result[unitigi].Nchunks.emplace_back(result[unitigi].bases.size(), anchorPositions[i]-anchorPositions[i-1]);
				for(size_t j = anchorPositions[i-1]; j < anchorPositions[i]; j++)
				{
					result[unitigi].bases.emplace_back(0);
				}
			}
		}
	});
	std::cerr << "got unitig consensuses" << std::endl;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	return result;
}

std::vector<double> getUnitigCoverages(const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<std::vector<size_t>>& chunkLengths, const std::vector<size_t>& chunkCoverages)
{
	std::vector<double> result;
	result.resize(unitigs.size(), 0);
	for (size_t i = 0; i < result.size(); i++)
	{
		double sum = 0;
		double divisor = 0;
		for (uint64_t node : unitigs[i])
		{
			size_t length = chunkLengths[node & maskUint64_t][chunkLengths[node & maskUint64_t].size() / 2];
			sum += length * chunkCoverages[node & maskUint64_t];
			divisor += length;
		}
		result[i] = sum / divisor;
	}
	return result;
}

SparseEdgeContainer filterOutZEdges(const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<size_t>& unitigLengths, const SparseEdgeContainer& allowedEdges, const std::vector<std::vector<size_t>>& lengths, const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverages, const std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>& chunkLocationInUnitig, const double approxOneHapCoverage)
{
	phmap::flat_hash_set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>> oneSideWantsToRemoveEdge;
	auto unitigCoverages = getUnitigCoverages(unitigs, lengths, coverages);
	std::vector<bool> singleCopyUnitigs;
	singleCopyUnitigs.resize(unitigs.size(), false);
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (unitigLengths[i] < 50000) continue;
		if (unitigCoverages[i] >= approxOneHapCoverage * 1.5) continue;
		singleCopyUnitigs[i] = true;
	}
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		if (!singleCopyUnitigs[i]) continue;
		std::pair<size_t, bool> lastNode { unitigs[i].back() & maskUint64_t, unitigs[i].back() & firstBitUint64_t };
		std::pair<size_t, bool> firstNode { unitigs[i][0] & maskUint64_t, (unitigs[i][0] ^ firstBitUint64_t) & firstBitUint64_t };
		if (allowedEdges.getEdges(lastNode).size() == 2)
		{
			std::pair<size_t, bool> first = allowedEdges.getEdges(lastNode)[0];
			std::pair<size_t, bool> second = allowedEdges.getEdges(lastNode)[1];
			auto canonkey = canon(lastNode, first);
			size_t firstEdgeCoverage = edgeCoverages.at(std::make_pair(canonkey.first.first + (canonkey.first.second ? firstBitUint64_t : 0), canonkey.second.first + (canonkey.second.second ? firstBitUint64_t : 0)));
			canonkey = canon(lastNode, second);
			size_t secondEdgeCoverage = edgeCoverages.at(std::make_pair(canonkey.first.first + (canonkey.first.second ? firstBitUint64_t : 0), canonkey.second.first + (canonkey.second.second ? firstBitUint64_t : 0)));
			if (singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[first.first]) & maskUint64_t] && singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[second.first]) & maskUint64_t])
			{
				if (firstEdgeCoverage < secondEdgeCoverage && allowedEdges.getEdges(reverse(first)).size() == 2 && allowedEdges.getEdges(reverse(second)).size() == 1)
				{
					std::pair<size_t, bool> third { std::numeric_limits<size_t>::max(), false };
					assert(allowedEdges.hasEdge(reverse(first), reverse(lastNode)));
					for (auto edge : allowedEdges.getEdges(reverse(first)))
					{
						if (edge.first != lastNode.first)
						{
							assert(third.first == std::numeric_limits<size_t>::max());
							third = edge;
						}
					}
					if (third.first != std::numeric_limits<size_t>::max())
					{
						if (singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[third.first]) & maskUint64_t])
						{
							canonkey = canon(reverse(first), third);
							size_t thirdEdgeCoverage = edgeCoverages.at(std::make_pair(canonkey.first.first + (canonkey.first.second ? firstBitUint64_t : 0), canonkey.second.first + (canonkey.second.second ? firstBitUint64_t : 0)));
							if (thirdEdgeCoverage > firstEdgeCoverage)
							{
								oneSideWantsToRemoveEdge.emplace(lastNode, first);
							}
						}
					}
				}
			}
			std::swap(first, second);
			std::swap(firstEdgeCoverage, secondEdgeCoverage);
			if (singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[first.first]) & maskUint64_t] && singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[second.first]) & maskUint64_t])
			{
				if (firstEdgeCoverage < secondEdgeCoverage && allowedEdges.getEdges(reverse(first)).size() == 2 && allowedEdges.getEdges(reverse(second)).size() == 1)
				{
					std::pair<size_t, bool> third { std::numeric_limits<size_t>::max(), false };
					assert(allowedEdges.hasEdge(reverse(first), reverse(lastNode)));
					for (auto edge : allowedEdges.getEdges(reverse(first)))
					{
						if (edge.first != lastNode.first)
						{
							assert(third.first == std::numeric_limits<size_t>::max());
							third = edge;
						}
					}
					if (third.first != std::numeric_limits<size_t>::max())
					{
						if (singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[third.first]) & maskUint64_t])
						{
							canonkey = canon(reverse(first), third);
							size_t thirdEdgeCoverage = edgeCoverages.at(std::make_pair(canonkey.first.first + (canonkey.first.second ? firstBitUint64_t : 0), canonkey.second.first + (canonkey.second.second ? firstBitUint64_t : 0)));
							if (thirdEdgeCoverage > firstEdgeCoverage)
							{
								oneSideWantsToRemoveEdge.emplace(lastNode, first);
							}
						}
					}
				}
			}
		}
		if (allowedEdges.getEdges(firstNode).size() == 2)
		{
			std::pair<size_t, bool> first = allowedEdges.getEdges(firstNode)[0];
			std::pair<size_t, bool> second = allowedEdges.getEdges(firstNode)[1];
			auto canonkey = canon(firstNode, first);
			size_t firstEdgeCoverage = edgeCoverages.at(std::make_pair(canonkey.first.first + (canonkey.first.second ? firstBitUint64_t : 0), canonkey.second.first + (canonkey.second.second ? firstBitUint64_t : 0)));
			canonkey = canon(firstNode, second);
			size_t secondEdgeCoverage = edgeCoverages.at(std::make_pair(canonkey.first.first + (canonkey.first.second ? firstBitUint64_t : 0), canonkey.second.first + (canonkey.second.second ? firstBitUint64_t : 0)));
			if (singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[first.first]) & maskUint64_t] && singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[second.first]) & maskUint64_t])
			{
				if (firstEdgeCoverage < secondEdgeCoverage && allowedEdges.getEdges(reverse(first)).size() == 2 && allowedEdges.getEdges(reverse(second)).size() == 1)
				{
					std::pair<size_t, bool> third { std::numeric_limits<size_t>::max(), false };
					assert(allowedEdges.hasEdge(reverse(first), reverse(firstNode)));
					for (auto edge : allowedEdges.getEdges(reverse(first)))
					{
						if (edge.first != firstNode.first)
						{
							assert(third.first == std::numeric_limits<size_t>::max());
							third = edge;
						}
					}
					if (third.first != std::numeric_limits<size_t>::max())
					{
						if (singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[third.first]) & maskUint64_t])
						{
							canonkey = canon(reverse(first), third);
							size_t thirdEdgeCoverage = edgeCoverages.at(std::make_pair(canonkey.first.first + (canonkey.first.second ? firstBitUint64_t : 0), canonkey.second.first + (canonkey.second.second ? firstBitUint64_t : 0)));
							if (thirdEdgeCoverage > firstEdgeCoverage)
							{
								oneSideWantsToRemoveEdge.emplace(firstNode, first);
							}
						}
					}
				}
			}
			std::swap(first, second);
			std::swap(firstEdgeCoverage, secondEdgeCoverage);
			if (singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[first.first]) & maskUint64_t] && singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[second.first]) & maskUint64_t])
			{
				if (firstEdgeCoverage < secondEdgeCoverage && allowedEdges.getEdges(reverse(first)).size() == 2 && allowedEdges.getEdges(reverse(second)).size() == 1)
				{
					std::pair<size_t, bool> third { std::numeric_limits<size_t>::max(), false };
					assert(allowedEdges.hasEdge(reverse(first), reverse(firstNode)));
					for (auto edge : allowedEdges.getEdges(reverse(first)))
					{
						if (edge.first != firstNode.first)
						{
							assert(third.first == std::numeric_limits<size_t>::max());
							third = edge;
						}
					}
					if (third.first != std::numeric_limits<size_t>::max())
					{
						if (singleCopyUnitigs[std::get<0>(chunkLocationInUnitig[third.first]) & maskUint64_t])
						{
							canonkey = canon(reverse(first), third);
							size_t thirdEdgeCoverage = edgeCoverages.at(std::make_pair(canonkey.first.first + (canonkey.first.second ? firstBitUint64_t : 0), canonkey.second.first + (canonkey.second.second ? firstBitUint64_t : 0)));
							if (thirdEdgeCoverage > firstEdgeCoverage)
							{
								oneSideWantsToRemoveEdge.emplace(firstNode, first);
							}
						}
					}
				}
			}
		}
	}
	SparseEdgeContainer result;
	result.resize(allowedEdges.size());
	for (size_t i = 0; i < allowedEdges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : allowedEdges.getEdges(fw))
		{
			if (oneSideWantsToRemoveEdge.count(std::make_pair(fw, edge)) == 1)
			{
				if (oneSideWantsToRemoveEdge.count(std::make_pair(reverse(edge), reverse(fw))) == 1)
				{
					continue;
				}
			}
			result.addEdge(fw, edge);
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : allowedEdges.getEdges(bw))
		{
			if (oneSideWantsToRemoveEdge.count(std::make_pair(bw, edge)) == 1)
			{
				if (oneSideWantsToRemoveEdge.count(std::make_pair(reverse(edge), reverse(bw))) == 1)
				{
					continue;
				}
			}
			result.addEdge(bw, edge);
		}
	}
	return result;
}

std::tuple<ChunkUnitigGraph, std::vector<std::vector<UnitigPath>>, std::vector<ConsensusString>> getChunkUnitigGraphInner(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const bool alsoConsensuses, const double approxOneHapCoverage, const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t numThreads, const size_t kmerSize)
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
	std::vector<std::tuple<uint64_t, size_t, size_t, size_t>> chunkLocationInUnitig; // unitig, index, startpos, endpos. reverse also reverses index and poses.
	std::tie(unitigs, unitigLengths, chunkLocationInUnitig) = getUnitigs(allowedNode, allowedEdges, lengths, edgeOverlaps);
	{
		allowedEdges = filterOutZEdges(unitigs, unitigLengths, allowedEdges, lengths, coverages, edgeCoverage, chunkLocationInUnitig, approxOneHapCoverage);
		std::tie(unitigs, unitigLengths, chunkLocationInUnitig) = getUnitigs(allowedNode, allowedEdges, lengths, edgeOverlaps);
	}
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
			assert(!NonexistantChunk(edge.first + (edge.second ? firstBitUint64_t : 0)));
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
			assert(!NonexistantChunk(edge.first + (edge.second ? firstBitUint64_t : 0)));
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
	result.coverages = getUnitigCoverages(unitigs, lengths, coverages);
	result.chunksInUnitig.resize(result.unitigLengths.size());
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		result.chunksInUnitig[i] = unitigs[i];
	}
	auto paths = getUnitigPaths(result, chunksPerRead, allowedNode, unitigs, chunkLocationInUnitig);
	std::vector<ConsensusString> consensuses;
	if (alsoConsensuses)
	{
		consensuses = getUnitigConsensuses(result, paths, sequenceIndex, chunksPerRead, lengths, coverages, chunkLocationInUnitig, edgeOverlaps, numThreads);
	}
	result.unitigChunkBreakpointPositions.resize(result.chunksInUnitig.size());
	for (size_t i = 0; i < result.chunksInUnitig.size(); i++)
	{
		result.unitigChunkBreakpointPositions[i].emplace_back(0);
		for (size_t j = 1; j < result.chunksInUnitig[i].size(); j++)
		{
			size_t chunkLength = lengths[result.chunksInUnitig[i][j-1] & maskUint64_t][lengths[result.chunksInUnitig[i][j-1] & maskUint64_t].size()/2];
			result.unitigChunkBreakpointPositions[i].emplace_back(result.unitigChunkBreakpointPositions[i].back() + chunkLength - kmerSize);
		}
	}
	return std::make_tuple(result, paths, consensuses);
}

std::tuple<ChunkUnitigGraph, std::vector<std::vector<UnitigPath>>, std::vector<ConsensusString>> getChunkUnitigGraphWithConsensuses(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t numThreads, const size_t kmerSize)
{
	return getChunkUnitigGraphInner(chunksPerRead, true, approxOneHapCoverage, sequenceIndex, numThreads, kmerSize);
}

std::tuple<ChunkUnitigGraph, std::vector<std::vector<UnitigPath>>> getChunkUnitigGraph(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize)
{
	FastaCompressor::CompressedStringIndex fakeSequences { 0, 0 };
	auto result = getChunkUnitigGraphInner(chunksPerRead, false, approxOneHapCoverage, fakeSequences, 1, kmerSize);
	return std::make_tuple(std::get<0>(result), std::get<1>(result));
}

void splitPerDiploidChunkWithNeighbors(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "splitting by diploid chunk with neighbors" << std::endl;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			maxChunk = std::max(maxChunk, std::get<2>(t) & maskUint64_t);
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t) + 1);
			occurrencesPerChunk[std::get<2>(t)].emplace_back(i, j);
		}
	}
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<size_t> chunkBelongsToUnitig;
	chunkBelongsToUnitig.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			assert((node & maskUint64_t) < chunkBelongsToUnitig.size());
			assert(chunkBelongsToUnitig[node & maskUint64_t] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[node & maskUint64_t] = i;
		}
	}
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> fwForksPerUnitig;
	std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> bwForksPerUnitig;
	fwForksPerUnitig.resize(graph.unitigLengths.size());
	bwForksPerUnitig.resize(graph.unitigLengths.size());
	std::vector<bool> hasFwFork;
	std::vector<bool> hasBwFork;
	hasFwFork.resize(graph.unitigLengths.size(), false);
	hasBwFork.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		if (graph.edges.getEdges(fw).size() == 2)
		{
			hasFwFork[i] = true;
		}
		std::pair<size_t, bool> bw { i, false };
		if (graph.edges.getEdges(bw).size() == 2)
		{
			hasBwFork[i] = true;
		}
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 1; k < readPaths[i][j].path.size(); k++)
			{
				uint64_t prev = readPaths[i][j].path[k-1];
				uint64_t curr = readPaths[i][j].path[k];
				if ((prev & firstBitUint64_t) && hasFwFork[prev & maskUint64_t])
				{
					std::pair<size_t, bool> currPair { curr & maskUint64_t, curr & firstBitUint64_t };
					if (currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, true))[0])
					{
						fwForksPerUnitig[prev & maskUint64_t].first.emplace_back(i);
					}
					else
					{
						assert(currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, true))[1]);
						fwForksPerUnitig[prev & maskUint64_t].second.emplace_back(i);
					}
				}
				if (((prev ^ firstBitUint64_t) & firstBitUint64_t) && hasBwFork[prev & maskUint64_t])
				{
					std::pair<size_t, bool> currPair { curr & maskUint64_t, curr & firstBitUint64_t };
					if (currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, false))[0])
					{
						bwForksPerUnitig[prev & maskUint64_t].first.emplace_back(i);
					}
					else
					{
						assert(currPair == graph.edges.getEdges(std::make_pair(prev & maskUint64_t, false))[1]);
						bwForksPerUnitig[prev & maskUint64_t].second.emplace_back(i);
					}
				}
				if ((curr & firstBitUint64_t) && hasBwFork[curr & maskUint64_t])
				{
					std::pair<size_t, bool> prevPair { prev & maskUint64_t, prev & firstBitUint64_t };
					prevPair.second = !prevPair.second;
					if (prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, false))[0])
					{
						bwForksPerUnitig[curr & maskUint64_t].first.emplace_back(i);
					}
					else
					{
						assert(prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, false))[1]);
						bwForksPerUnitig[curr & maskUint64_t].second.emplace_back(i);
					}
				}
				if (((curr ^ firstBitUint64_t) & firstBitUint64_t) && hasFwFork[curr & maskUint64_t])
				{
					std::pair<size_t, bool> prevPair { prev & maskUint64_t, prev & firstBitUint64_t };
					prevPair.second = !prevPair.second;
					if (prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, true))[0])
					{
						fwForksPerUnitig[curr & maskUint64_t].first.emplace_back(i);
					}
					else
					{
						assert(prevPair == graph.edges.getEdges(std::make_pair(curr & maskUint64_t, true))[1]);
						fwForksPerUnitig[curr & maskUint64_t].second.emplace_back(i);
					}
				}
			}
		}
	}
	for (size_t i = 0; i < fwForksPerUnitig.size(); i++)
	{
		std::sort(fwForksPerUnitig[i].first.begin(), fwForksPerUnitig[i].first.end());
		std::sort(fwForksPerUnitig[i].second.begin(), fwForksPerUnitig[i].second.end());
		std::sort(bwForksPerUnitig[i].first.begin(), bwForksPerUnitig[i].first.end());
		std::sort(bwForksPerUnitig[i].second.begin(), bwForksPerUnitig[i].second.end());
		if (hasFwFork[i])
		{
			if (intersectSize(fwForksPerUnitig[i].first, fwForksPerUnitig[i].second) >= 1)
			{
				hasFwFork[i] = false;
			}
			if (phmap::flat_hash_set<size_t>(fwForksPerUnitig[i].first.begin(), fwForksPerUnitig[i].first.end()).size() != fwForksPerUnitig[i].first.size()) hasFwFork[i] = false;
			if (phmap::flat_hash_set<size_t>(fwForksPerUnitig[i].second.begin(), fwForksPerUnitig[i].second.end()).size() != fwForksPerUnitig[i].second.size()) hasFwFork[i] = false;
			if (fwForksPerUnitig[i].first.size() < 2) hasFwFork[i] = false;
			if (fwForksPerUnitig[i].second.size() < 2) hasFwFork[i] = false;
		}
		if (hasBwFork[i])
		{
			if (intersectSize(bwForksPerUnitig[i].first, bwForksPerUnitig[i].second) >= 1)
			{
				hasBwFork[i] = false;
			}
			if (phmap::flat_hash_set<size_t>(bwForksPerUnitig[i].first.begin(), bwForksPerUnitig[i].first.end()).size() != bwForksPerUnitig[i].first.size()) hasBwFork[i] = false;
			if (phmap::flat_hash_set<size_t>(bwForksPerUnitig[i].second.begin(), bwForksPerUnitig[i].second.end()).size() != bwForksPerUnitig[i].second.size()) hasBwFork[i] = false;
			if (bwForksPerUnitig[i].first.size() < 2) hasBwFork[i] = false;
			if (bwForksPerUnitig[i].second.size() < 2) hasBwFork[i] = false;
		}
	}
	std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> fwForksPerUnitigSet;
	std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> bwForksPerUnitigSet;
	fwForksPerUnitigSet.resize(graph.unitigLengths.size());
	bwForksPerUnitigSet.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < fwForksPerUnitig.size(); i++)
	{
		if (hasFwFork[i])
		{
			fwForksPerUnitigSet[i].first.insert(fwForksPerUnitig[i].first.begin(), fwForksPerUnitig[i].first.end());
			fwForksPerUnitigSet[i].second.insert(fwForksPerUnitig[i].second.begin(), fwForksPerUnitig[i].second.end());
		}
		if (hasBwFork[i])
		{
			bwForksPerUnitigSet[i].first.insert(bwForksPerUnitig[i].first.begin(), bwForksPerUnitig[i].first.end());
			bwForksPerUnitigSet[i].second.insert(bwForksPerUnitig[i].second.begin(), bwForksPerUnitig[i].second.end());
		}
	}
	size_t nextNum = 0;
	size_t countSplitted = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &sequenceIndex, &hasFwFork, &hasBwFork, &fwForksPerUnitigSet, &bwForksPerUnitigSet,  &countSplitted, &chunkBelongsToUnitig, &rawReadLengths, kmerSize](const size_t i)
	{
		if (chunkBelongsToUnitig[i] == std::numeric_limits<size_t>::max())
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		size_t unitig = chunkBelongsToUnitig[i];
		bool anythingToPhase = false;
		std::vector<std::pair<phmap::flat_hash_set<size_t>, phmap::flat_hash_set<size_t>>> maybePhaseGroups;
		if (hasFwFork[unitig])
		{
			if (fwForksPerUnitigSet[unitig].first.size() >= 2 && fwForksPerUnitigSet[unitig].second.size() >= 2)
			{
				size_t firstsHere = 0;
				size_t secondsHere = 0;
				bool valid = true;
				for (auto pair : occurrencesPerChunk[i])
				{
					if (fwForksPerUnitigSet[unitig].first.count(pair.first) == 1)
					{
						firstsHere += 1;
						if (fwForksPerUnitigSet[unitig].second.count(pair.first) == 1)
						{
							valid = false;
							break;
						}
					}
					if (fwForksPerUnitigSet[unitig].second.count(pair.first) == 1)
					{
						secondsHere += 1;
					}
				}
				if (firstsHere >= 2 && secondsHere >= 2 && valid)
				{
					maybePhaseGroups.emplace_back(fwForksPerUnitigSet[unitig]);
				}
			}
		}
		if (hasBwFork[unitig])
		{
			if (bwForksPerUnitigSet[unitig].first.size() >= 2 && bwForksPerUnitigSet[unitig].second.size() >= 2)
			{
				size_t firstsHere = 0;
				size_t secondsHere = 0;
				bool valid = true;
				for (auto pair : occurrencesPerChunk[i])
				{
					if (bwForksPerUnitigSet[unitig].first.count(pair.first) == 1)
					{
						firstsHere += 1;
						if (bwForksPerUnitigSet[unitig].second.count(pair.first) == 1)
						{
							valid = false;
							break;
						}
					}
					if (bwForksPerUnitigSet[unitig].second.count(pair.first) == 1)
					{
						secondsHere += 1;
					}
				}
				if (firstsHere >= 2 && secondsHere >= 2 && valid)
				{
					maybePhaseGroups.emplace_back(bwForksPerUnitigSet[unitig]);
				}
			}
		}
		if (maybePhaseGroups.size() == 0)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		size_t nextGroupNum = 0;
		std::vector<phmap::flat_hash_set<size_t>> solidKmersPerOccurrence;
		solidKmersPerOccurrence.resize(occurrencesPerChunk[i].size());
		phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> kmerClusterToNumber;
		std::vector<TwobitString> chunkSequences;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			chunkSequences.emplace_back(getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, occurrencesPerChunk[i][j].first, occurrencesPerChunk[i][j].second));
		}
		iterateSolidKmers(chunkSequences, kmerSize, occurrencesPerChunk[i].size(), true, false, [&solidKmersPerOccurrence, &kmerClusterToNumber](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
		{
			size_t kmerkey = 0;
			std::pair<size_t, size_t> key { kmer, clusterIndex };
			if (kmerClusterToNumber.count(key) == 0)
			{
				kmerkey = kmerClusterToNumber.size();
				kmerClusterToNumber[key] = kmerkey;
			}
			else
			{
				kmerkey = kmerClusterToNumber.at(key);
			}
			solidKmersPerOccurrence[occurrenceID].emplace(kmerkey);
		});
		phmap::flat_hash_map<size_t, size_t> kmerToNumber;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			assert(!NonexistantChunk(std::get<2>(t)));
			phmap::flat_hash_set<size_t> kmersHere;
			iterateKmers(chunkSequences[j], 0, chunkSequences[j].size()-1, true, kmerSize, [&kmersHere](const size_t kmer, const size_t pos)
			{
				kmersHere.insert(kmer);
			});
			for (auto kmer : kmersHere)
			{
				size_t kmerkey = 0;
				if (kmerToNumber.count(kmer) == 0)
				{
					kmerkey = kmerToNumber.size() + kmerClusterToNumber.size();
					kmerToNumber[kmer] = kmerkey;
				}
				else
				{
					kmerkey = kmerToNumber.at(kmer);
				}
				solidKmersPerOccurrence[j].emplace(kmerkey);
			}
		}
		{
			std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, size_t>> alleles; // startkmer, startcluster, endkmer, endcluster, occurrenceID, allele
			size_t lastOccurrence = std::numeric_limits<size_t>::max();
			size_t lastKmer = std::numeric_limits<size_t>::max();
			size_t lastCluster = std::numeric_limits<size_t>::max();
			size_t lastKmerPos = std::numeric_limits<size_t>::max();
			phmap::flat_hash_map<std::string, size_t> alleleNumbers;
			iterateSolidKmers(chunkSequences, kmerSize, occurrencesPerChunk[i].size(), false, false, [&occurrencesPerChunk, &chunkSequences, &alleles, &lastOccurrence, &lastKmer, &lastCluster, &lastKmerPos, &alleleNumbers, kmerSize, i](size_t occurrenceID, size_t chunkStartPos, size_t chunkEndPos, uint64_t node, size_t kmer, size_t clusterIndex, size_t pos)
			{
				if (occurrenceID != lastOccurrence || lastKmer == std::numeric_limits<size_t>::max())
				{
					lastOccurrence = occurrenceID;
					lastKmer = kmer;
					lastCluster = clusterIndex;
					lastKmerPos = pos;
					return;
				}
				size_t allele;
				allele = getAllele(chunkSequences[occurrenceID], lastKmerPos, pos + kmerSize, alleleNumbers, true);
				alleles.emplace_back(lastKmer, lastCluster, kmer, clusterIndex, occurrenceID, allele);
				lastOccurrence = occurrenceID;
				lastKmer = kmer;
				lastCluster = clusterIndex;
				lastKmerPos = pos;
			});
			std::sort(alleles.begin(), alleles.end());
			std::vector<std::vector<std::vector<size_t>>> occurrencesPerAlleleSite;
			occurrencesPerAlleleSite = getAlleleOccurrences(alleles, occurrencesPerChunk[i].size());
			size_t nextNum = kmerToNumber.size() + kmerClusterToNumber.size();
			for (size_t j = 0; j < occurrencesPerAlleleSite.size(); j++)
			{
				if (occurrencesPerAlleleSite[j].size() < 2) continue;
				for (size_t allele = 0; allele < occurrencesPerAlleleSite[j].size(); allele++)
				{
					for (size_t k : occurrencesPerAlleleSite[j][allele])
					{
						solidKmersPerOccurrence[k].emplace(nextNum);
					}
					nextNum += 1;
				}
			}
		}
		{
			phmap::flat_hash_map<size_t, size_t> kmerCoverage;
			for (size_t j = 0; j < solidKmersPerOccurrence.size(); j++)
			{
				for (auto kmer : solidKmersPerOccurrence[j])
				{
					kmerCoverage[kmer] += 1;
				}
			}
			phmap::flat_hash_set<size_t> removeKmers;
			for (auto pair : kmerCoverage)
			{
				if (pair.second < 2 || pair.second + 2 > occurrencesPerChunk[i].size())
				{
					removeKmers.insert(pair.first);
				}
			}
			if (removeKmers.size() >= 1)
			{
				for (size_t j = 0; j < solidKmersPerOccurrence.size(); j++)
				{
					for (auto kmer : removeKmers)
					{
						if (solidKmersPerOccurrence[j].count(kmer) == 1)
						{
							solidKmersPerOccurrence[j].erase(kmer);
						}
					}
				}
			}
		}
		phmap::flat_hash_map<size_t, std::vector<size_t>> phaseGroupsPerRead;
		for (size_t phaseGroup = 0; phaseGroup < maybePhaseGroups.size(); phaseGroup++)
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
				const phmap::flat_hash_set<size_t> kmersHere = solidKmersPerOccurrence[j];
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
			assert(hasFirst);
			assert(hasSecond);
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
				const phmap::flat_hash_set<size_t>& kmersHere = solidKmersPerOccurrence[j];
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
				const phmap::flat_hash_set<size_t>& kmersHere = solidKmersPerOccurrence[j];
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
		phmap::flat_hash_map<size_t, size_t> groupCoverage;
		for (size_t i = 0; i < parent.size(); i++)
		{
			groupCoverage[find(parent, i)] += 1;
		}
		bool valid = true;
		if (groupCoverage.size() != 2) valid = false;
		for (auto pair : groupCoverage)
		{
			if (pair.second < 3) valid = false;
		}
		if (!valid)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			phmap::flat_hash_map<size_t, size_t> clusterToNode;
			size_t first = nextNum;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				if (clusterToNode.count(find(parent, j)) == 1) continue;
				clusterToNode[find(parent, j)] = nextNum;
				nextNum += 1;
			}
			size_t last = nextNum-1;
			if (clusterToNode.size() >= 2) countSplitted += 1;
			std::cerr << "diploid chunk with neighbors splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks range " << first << " - " << last << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, j)) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
	std::cerr << "diploid chunk with neighbors splitted " << countSplitted << " chunks" << std::endl;
}

void writeUnitigSequences(const std::string& filename, const std::vector<ConsensusString>& sequences)
{
	std::ofstream file { filename };
	for (size_t i = 0; i < sequences.size(); i++)
	{
		file << ">unitig_" << i << std::endl;
		std::string str = sequences[i].bases.toString();
		for (auto pair : sequences[i].Nchunks)
		{
			for (size_t j = pair.first; j < pair.first+pair.second; j++)
			{
				str[j] = 'N';
			}
		}
		file << str << std::endl;
	}
}

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getBidirectedChunks(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, uint64_t> pairToNode;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> fixedChunksPerRead;
	fixedChunksPerRead.resize(rawChunksPerRead.size()/2);
	for (size_t i = 0; i < rawChunksPerRead.size()/2; i++)
	{
		size_t other = i + rawChunksPerRead.size()/2;
		fixedChunksPerRead[i] = rawChunksPerRead[i];
		if (!(rawChunksPerRead[i].size() == rawChunksPerRead[other].size()))
		{
			std::cerr << i << " " << other << std::endl;
		}
		assert(rawChunksPerRead[i].size() == rawChunksPerRead[other].size());
		for (size_t j = 0; j < rawChunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(rawChunksPerRead[i][j])))
			{
				if (!NonexistantChunk(std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j])))
				{
					std::cerr << i << " " << other << std::endl;
					std::cerr << j << " " << rawChunksPerRead[other].size()-1-j << std::endl;
					std::cerr << ((std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & maskUint64_t) << std::endl;
				}
				assert(NonexistantChunk(std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j])));
				continue;
			}
			if (NonexistantChunk(std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j])))
			{
				std::cerr << i << " " << other << std::endl;
				std::cerr << j << " " << rawChunksPerRead[other].size()-1-j << std::endl;
				std::cerr << ((std::get<2>(rawChunksPerRead[i][j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(rawChunksPerRead[i][j]) & maskUint64_t) << std::endl;
			}
			assert(!NonexistantChunk(std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j])));
			assert((std::get<2>(rawChunksPerRead[i][j]) & firstBitUint64_t) == firstBitUint64_t);
			assert((std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & firstBitUint64_t) == firstBitUint64_t);
			if ((std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & maskUint64_t) < (std::get<2>(rawChunksPerRead[i][j]) & maskUint64_t))
			{
				std::get<2>(fixedChunksPerRead[i][j]) = (std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & maskUint64_t);
			}
			/*
			uint64_t fwnode = std::get<2>(rawChunksPerRead[i][j]);
			uint64_t bwnode = std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) ^ firstBitUint64_t;
			std::pair<uint64_t, uint64_t> fwpair { fwnode, bwnode };
			std::pair<uint64_t, uint64_t> bwpair { bwnode ^ firstBitUint64_t, fwnode ^ firstBitUint64_t };
			bool fw = true;
			if (bwpair < fwpair)
			{
				std::swap(bwpair, fwpair);
				fw = false;
			}
			if (pairToNode.count(fwpair) == 0)
			{
				uint64_t next = pairToNode.size();
				pairToNode[fwpair] = next;
			}
			std::get<2>(fixedChunksPerRead[i][j]) = pairToNode.at(fwpair) + (fw ? firstBitUint64_t : 0);
			*/
		}
	}
	return fixedChunksPerRead;
}

void writeBidirectedUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::string& unitigSequencesFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t numThreads, const size_t kmerSize)
{
	std::cerr << "writing unitig graph " << graphFile << std::endl;
	auto fixedChunksPerRead = getBidirectedChunks(rawChunksPerRead);
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::vector<ConsensusString> unitigSequences;
	std::tie(graph, readPaths, unitigSequences) = getChunkUnitigGraphWithConsensuses(fixedChunksPerRead, sequenceIndex, approxOneHapCoverage, numThreads, kmerSize);
	writeUnitigGraph(graphFile, graph);
	std::cerr << "writing unitig paths " << pathsFile << std::endl;
	writeUnitigPaths(pathsFile, graph, readPaths, sequenceIndex, rawReadLengths);
	std::cerr << "writing unitig sequences " << unitigSequencesFile << std::endl;
	writeUnitigSequences(unitigSequencesFile, unitigSequences);
}

void writeBidirectedUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "writing unitig graph " << graphFile << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	auto fixedChunksPerRead = getBidirectedChunks(chunksPerRead);
	std::tie(graph, readPaths) = getChunkUnitigGraph(fixedChunksPerRead, approxOneHapCoverage, kmerSize);
	writeUnitigGraph(graphFile, graph);
	std::cerr << "writing unitig paths " << pathsFile << std::endl;
	writeUnitigPaths(pathsFile, graph, readPaths, sequenceIndex, rawReadLengths);
}

void writeUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "writing unitig digraph " << graphFile << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	writeUnitigGraph(graphFile, graph);
	std::cerr << "writing unitig paths " << pathsFile << std::endl;
	writeUnitigPaths(pathsFile, graph, readPaths, sequenceIndex, rawReadLengths);
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

void cleanTips(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const double maxRemoveCoverage, const size_t kmerSize)
{
	std::cerr << "cleaning tips" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	double reestimatedCoverage = estimateCoverage(graph);
	if (reestimatedCoverage < 1) reestimatedCoverage = approxOneHapCoverage;
	phmap::flat_hash_set<uint64_t> solidFork;
	phmap::flat_hash_set<uint64_t> acceptableFork;
	{
		std::vector<std::pair<uint64_t, std::vector<std::vector<size_t>>>> forkReads = getUnitigForkReads(graph, readPaths);
		solidFork = getSolidForks(forkReads, std::min(approxOneHapCoverage, reestimatedCoverage) * 0.5);
		acceptableFork = getAcceptableForks(forkReads, solidFork, 3);
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

void resolveVerySmallNodes(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const size_t maxResolveSize, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving very small nodes" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	phmap::flat_hash_set<size_t> resolvableTinies;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.chunksInUnitig[i].size() != 1) continue;
		if (graph.unitigLengths[i] >= maxResolveSize) continue;
		if (graph.edges.getEdges(std::make_pair(i, true)).size() < 2 && graph.edges.getEdges(std::make_pair(i, false)).size() < 2) continue;
		resolvableTinies.insert(graph.chunksInUnitig[i][0] & maskUint64_t);
	}
	size_t nextNum = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (!NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) nextNum = std::max(nextNum, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			if (resolvableTinies.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			if (std::get<0>(chunksPerRead[i][j]) == 0)
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (std::get<1>(chunksPerRead[i][j])+1 == rawReadLengths[i])
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1])+1)
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
			if (j+1 < chunksPerRead[i].size() && std::get<1>(chunksPerRead[i][j])+1 == std::get<1>(chunksPerRead[i][j-1]))
			{
				resolvableTinies.erase(std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
		}
	}
	for (size_t chunk : resolvableTinies)
	{
		std::cerr << "resolve tiny chunk " << chunk << std::endl;
	}
	nextNum += 1;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeToNumber;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<size_t> eraseIndices;
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunks;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (resolvableTinies.count(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) == 0) continue;
			eraseIndices.emplace_back(j);
			if (j > 0)
			{
				uint64_t here = std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t;
				uint64_t neighbor = std::get<2>(chunksPerRead[i][j-1]) ^ firstBitUint64_t;
				if (!NonexistantChunk(neighbor))
				{
					std::pair<uint64_t, uint64_t> key { here, neighbor };
					if (edgeToNumber.count(key) == 0)
					{
						edgeToNumber[key] = nextNum;
						nextNum += 1;
					}
					uint64_t node = edgeToNumber.at(key);
					newChunks.emplace_back(std::get<0>(chunksPerRead[i][j])-1, std::get<1>(chunksPerRead[i][j]), node);
				}
			}
			if (j+1 < chunksPerRead[i].size())
			{
				uint64_t here = std::get<2>(chunksPerRead[i][j]);
				uint64_t neighbor = std::get<2>(chunksPerRead[i][j+1]);
				if (!NonexistantChunk(neighbor))
				{
					std::pair<uint64_t, uint64_t> key { here, neighbor };
					if (edgeToNumber.count(key) == 0)
					{
						edgeToNumber[key] = nextNum;
						nextNum += 1;
					}
					uint64_t node = edgeToNumber.at(key) + firstBitUint64_t;
					newChunks.emplace_back(std::get<0>(chunksPerRead[i][j]), std::get<1>(chunksPerRead[i][j])+1, node);
				}
			}
		}
		if (eraseIndices.size() == 0 && newChunks.size() == 0) continue;
		for (size_t j = eraseIndices.size()-1; j < eraseIndices.size(); j--)
		{
			size_t index = eraseIndices[j];
			chunksPerRead[i].erase(chunksPerRead[i].begin()+index);
		}
		chunksPerRead[i].insert(chunksPerRead[i].end(), newChunks.begin(), newChunks.end());
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
	}
}

void resolveSemiAmbiguousUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving semiambiguous resolvable unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> tripletsPerUnitig;
	tripletsPerUnitig.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				uint64_t prev, after;
				if (mid & firstBitUint64_t)
				{
					prev = readPaths[i][j].path[k-1];
					after = readPaths[i][j].path[k+1];
				}
				else
				{
					prev = readPaths[i][j].path[k+1] ^ firstBitUint64_t;
					after = readPaths[i][j].path[k-1] ^ firstBitUint64_t;
				}
				tripletsPerUnitig[mid & maskUint64_t][std::make_pair(prev, after)] += 1;
			}
		}
	}
	std::vector<bool> canResolveUnitig;
	canResolveUnitig.resize(graph.unitigLengths.size(), false);
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> prevToAllele;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> afterToAllele;
	prevToAllele.resize(graph.unitigLengths.size());
	afterToAllele.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (tripletsPerUnitig[i].size() < 2) continue;
		phmap::flat_hash_set<uint64_t> coveredPrev;
		phmap::flat_hash_set<uint64_t> coveredAfter;
		phmap::flat_hash_map<uint64_t, uint64_t> parent;
		bool valid = true;
		for (auto triplet : tripletsPerUnitig[i])
		{
			if (triplet.second < 2)
			{
				valid = false;
				break;
			}
			coveredPrev.insert(triplet.first.first);
			coveredAfter.insert(triplet.first.second);
			assert((triplet.first.first & (firstBitUint64_t >> 1)) == 0);
			uint64_t prev = triplet.first.first + (firstBitUint64_t >> 1);
			uint64_t after = triplet.first.second;
			if (parent.count(prev) == 0) parent[prev] = prev;
			if (parent.count(after) == 0) parent[after] = after;
			merge(parent, prev, after);
		}
		if (!valid) continue;
		assert(coveredPrev.size() <= graph.edges.getEdges(std::make_pair(i, false)).size());
		assert(coveredAfter.size() <= graph.edges.getEdges(std::make_pair(i, true)).size());
		if (coveredPrev.size() < graph.edges.getEdges(std::make_pair(i, false)).size()) continue;
		if (coveredAfter.size() < graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		phmap::flat_hash_map<uint64_t, size_t> clusterToAllele;
		size_t alleleNum = 0;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t after = triplet.first.second;
			uint64_t cluster = find(parent, after);
			if (clusterToAllele.count(cluster) == 0)
			{
				clusterToAllele[cluster] = alleleNum;
				alleleNum += 1;
			}
		}
		assert(clusterToAllele.size() == alleleNum);
		assert(clusterToAllele.size() >= 1);
		if (clusterToAllele.size() == 1) continue;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t prev = triplet.first.first;
			uint64_t after = triplet.first.second;
			prevToAllele[i][prev] = clusterToAllele.at(find(parent, after));
			afterToAllele[i][after] = clusterToAllele.at(find(parent, after));
		}
		canResolveUnitig[i] = true;
	}
	phmap::flat_hash_set<size_t> readsAllowedToDiscard;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				if (!canResolveUnitig[mid & maskUint64_t]) continue;
				uint64_t prev = std::numeric_limits<size_t>::max();
				uint64_t after = std::numeric_limits<size_t>::max();
				if (k > 0)
				{
					prev = readPaths[i][j].path[k-1];
				}
				else if (j > 0)
				{
					prev = readPaths[i][j-1].path.back();
				}
				if (k+1 < readPaths[i][j].path.size())
				{
					after = readPaths[i][j].path[k+1];
				}
				else if (j+1 < readPaths[i].size())
				{
					after = readPaths[i][j+1].path[0];
				}
				if ((mid ^ firstBitUint64_t) & firstBitUint64_t)
				{
					std::swap(prev, after);
					if (prev != std::numeric_limits<size_t>::max()) prev ^= firstBitUint64_t;
					if (after != std::numeric_limits<size_t>::max()) after ^= firstBitUint64_t;
				}
				if (prev == std::numeric_limits<size_t>::max() && after == std::numeric_limits<size_t>::max())
				{
					if (readPaths[i].size() == 1 && readPaths[i][0].path.size() == 1)
					{
						readsAllowedToDiscard.insert(i);
						continue;
					}
					else
					{
						canResolveUnitig[mid & maskUint64_t] = false;
						continue;
					}
				}
				size_t prevAllele = std::numeric_limits<size_t>::max();
				size_t afterAllele = std::numeric_limits<size_t>::max();
				if (prev != std::numeric_limits<size_t>::max() && prevToAllele[mid & maskUint64_t].count(prev) == 1)
				{
					prevAllele = prevToAllele[mid & maskUint64_t].at(prev);
				}
				if (after != std::numeric_limits<size_t>::max() && afterToAllele[mid & maskUint64_t].count(after) == 1)
				{
					afterAllele = afterToAllele[mid & maskUint64_t].at(after);
				}
				if (prevAllele == std::numeric_limits<size_t>::max() && afterAllele == std::numeric_limits<size_t>::max())
				{
					canResolveUnitig[mid & maskUint64_t] = false;
					continue;
				}
				if (prevAllele != std::numeric_limits<size_t>::max() && afterAllele != std::numeric_limits<size_t>::max() && prevAllele != afterAllele)
				{
					canResolveUnitig[mid & maskUint64_t] = false;
					continue;
				}
			}
		}
	}
	assert(readPaths.size() == chunksPerRead.size());
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			assert(chunkBelongsToUnitig[chunk] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, size_t>> chunkAlleleReplacement;
	size_t nextNum = 0;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> unitigAllelesInThisRead;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (!canResolveUnitig[readPaths[i][j].path[k] & maskUint64_t]) continue;
				uint64_t prev = std::numeric_limits<size_t>::max();
				uint64_t after = std::numeric_limits<size_t>::max();
				if (k > 0)
				{
					prev = readPaths[i][j].path[k-1];
				}
				else if (j > 0)
				{
					prev = readPaths[i][j-1].path.back();
				}
				if (k+1 < readPaths[i][j].path.size())
				{
					after = readPaths[i][j].path[k+1];
				}
				else if (j+1 < readPaths[i].size())
				{
					after = readPaths[i][j+1].path[0];
				}
				if ((readPaths[i][j].path[k] ^ firstBitUint64_t) & firstBitUint64_t)
				{
					std::swap(prev, after);
					if (prev != std::numeric_limits<size_t>::max()) prev ^= firstBitUint64_t;
					if (after != std::numeric_limits<size_t>::max()) after ^= firstBitUint64_t;
				}
				if (prev == std::numeric_limits<size_t>::max() && after == std::numeric_limits<size_t>::max())
				{
					assert(readsAllowedToDiscard.count(i) == 1);
					assert(readPaths[i].size() == 1);
					assert(readPaths[i][0].path.size() == 1);
					continue;
				}
				assert(prev != std::numeric_limits<size_t>::max() || after != std::numeric_limits<size_t>::max());
				size_t alleleBefore = std::numeric_limits<size_t>::max();
				size_t alleleAfter = std::numeric_limits<size_t>::max();
				if (prev != std::numeric_limits<size_t>::max() && prevToAllele[readPaths[i][j].path[k] & maskUint64_t].count(prev) == 1) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(prev);
				if (after != std::numeric_limits<size_t>::max() && afterToAllele[readPaths[i][j].path[k] & maskUint64_t].count(after) == 1) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(after);
				assert(alleleBefore != std::numeric_limits<size_t>::max() || alleleAfter != std::numeric_limits<size_t>::max());
				assert(alleleBefore == alleleAfter || alleleBefore == std::numeric_limits<size_t>::max() || alleleAfter == std::numeric_limits<size_t>::max());
				size_t allele = alleleBefore;
				if (alleleBefore == std::numeric_limits<size_t>::max()) allele = alleleAfter;
				unitigAllelesInThisRead.emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k] & maskUint64_t, allele);
			}
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			uint64_t node = std::get<2>(chunksPerRead[i][j]);
			if (NonexistantChunk(node)) continue;
			size_t allele = std::numeric_limits<size_t>::max();
			if ((node & maskUint64_t) < chunkBelongsToUnitig.size() && chunkBelongsToUnitig[node & maskUint64_t] != std::numeric_limits<size_t>::max())
			{
				size_t unitig = chunkBelongsToUnitig[node & maskUint64_t];
				if (canResolveUnitig[unitig])
				{
					if (readsAllowedToDiscard.count(i) == 1)
					{
						assert(readPaths[i].size() == 1);
						assert(readPaths[i][0].path.size() == 1);
						std::get<2>(chunksPerRead[i][j]) = std::numeric_limits<size_t>::max();
						continue;
					}
					for (auto t : unitigAllelesInThisRead)
					{
						if (std::get<2>(t) != unitig) continue;
						if (std::get<0>(t) > std::get<0>(chunksPerRead[i][j])) continue;
						if (std::get<1>(t) < std::get<1>(chunksPerRead[i][j])) continue;
						assert(allele == std::numeric_limits<size_t>::max());
						allele = std::get<3>(t);
					}
					assert(allele != std::numeric_limits<size_t>::max());
				}
			}
			if (chunkAlleleReplacement.count(node & maskUint64_t) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t];
			}
			if (chunkAlleleReplacement.at(node & maskUint64_t).count(allele) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t][allele] = nextNum;
				nextNum += 1;
			}
			uint64_t replacement = chunkAlleleReplacement.at(node & maskUint64_t).at(allele);
			assert((replacement & firstBitUint64_t) == 0);
			replacement += (node & firstBitUint64_t);
			std::get<2>(chunksPerRead[i][j]) = replacement;
		}
	}
	size_t countUnitigsReplaced = 0;
	size_t countChunksReplaced = 0;
	for (size_t i = 0; i < canResolveUnitig.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		countUnitigsReplaced += 1;
		countChunksReplaced += graph.chunksInUnitig[i].size();
	}
	std::cerr << "semiambiguously resolvable unitig resolution resolved " << countUnitigsReplaced << " unitigs, " << countChunksReplaced << " chunks" << std::endl;
}

void resolveBetweenLongUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t longNodeThreshold, const size_t maxDistance, const size_t kmerSize)
{
	const uint64_t noNodeWithinDistance = (std::numeric_limits<size_t>::max() ^ (firstBitUint64_t>>1)) - 1;
	std::cerr << "resolving between long unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	phmap::flat_hash_set<size_t> chunkInsideLongUnitig;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.unitigLengths[i] < longNodeThreshold) continue;
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			chunkInsideLongUnitig.emplace(node & maskUint64_t);
		}
	}
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> longNodesInReads;
	longNodesInReads.resize(readPaths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (graph.unitigLengths[readPaths[i][j].path[k] & maskUint64_t] < longNodeThreshold) continue;
				longNodesInReads[i].emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k]);
			}
		}
		for (size_t j = 1; j < longNodesInReads[i].size(); j++)
		{
			assert(std::get<0>(longNodesInReads[i][j]) > std::get<0>(longNodesInReads[i][j-1]));
			assert(std::get<1>(longNodesInReads[i][j]) > std::get<1>(longNodesInReads[i][j-1]));
		}
	}
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
	size_t countSplitted = 0;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&chunksPerRead, &occurrencesPerChunk, &longNodesInReads, &rawReadLengths, &nextNum, &resultMutex, &chunkInsideLongUnitig, &countSplitted, maxDistance, noNodeWithinDistance](const size_t chunki)
	{
		if (chunkInsideLongUnitig.count(chunki) == 1 || occurrencesPerChunk[chunki].size() < 4)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		std::vector<std::tuple<uint64_t, uint64_t, size_t, size_t>> predecessorsAndSuccessors;
		for (auto t : occurrencesPerChunk[chunki])
		{
			size_t read = t.first;
			size_t offset = t.second;
			bool fw = std::get<2>(chunksPerRead[read][offset]) & firstBitUint64_t;
			size_t chunkstartpos = std::get<0>(chunksPerRead[read][offset]);
			size_t chunkendpos = std::get<1>(chunksPerRead[read][offset]);
			uint64_t predecessor = std::numeric_limits<size_t>::max();
			uint64_t successor = std::numeric_limits<size_t>::max();
			for (size_t j = 0; j < longNodesInReads[read].size(); j++)
			{
				if (std::get<1>(longNodesInReads[read][j]) < chunkendpos && std::get<1>(longNodesInReads[read][j])+maxDistance > chunkstartpos)
				{
					predecessor = std::get<2>(longNodesInReads[read][j]);
				}
				if (std::get<0>(longNodesInReads[read][j]) > chunkstartpos && std::get<0>(longNodesInReads[read][j]) < chunkendpos+maxDistance && successor == std::numeric_limits<size_t>::max())
				{
					successor = std::get<2>(longNodesInReads[read][j]);
				}
			}
			if (chunkstartpos > maxDistance && predecessor == std::numeric_limits<size_t>::max()) predecessor = noNodeWithinDistance;
			if (chunkendpos + maxDistance < rawReadLengths[read] && successor == std::numeric_limits<size_t>::max()) successor = noNodeWithinDistance;
			if (!fw)
			{
				std::swap(predecessor, successor);
				if (predecessor != std::numeric_limits<size_t>::max() && predecessor != noNodeWithinDistance) predecessor ^= firstBitUint64_t;
				if (successor != std::numeric_limits<size_t>::max() && successor != noNodeWithinDistance) successor ^= firstBitUint64_t;
			}
			predecessorsAndSuccessors.emplace_back(predecessor, successor, read, offset);
		}
		assert(predecessorsAndSuccessors.size() == occurrencesPerChunk[chunki].size());
		bool hasAny = false;
		phmap::flat_hash_set<uint64_t> mergedSomewhere;
		phmap::flat_hash_map<uint64_t, uint64_t> parent;
		for (auto t : predecessorsAndSuccessors)
		{
			assert(std::get<1>(t) == std::numeric_limits<size_t>::max() || ((std::get<1>(t) & (firstBitUint64_t>>1)) == 0));
			if (std::get<0>(t) != std::numeric_limits<size_t>::max() && parent.count(std::get<0>(t)) == 0)
			{
				hasAny = true;
				parent[std::get<0>(t)] = std::get<0>(t);
			}
			if (std::get<1>(t) != std::numeric_limits<size_t>::max() && parent.count(std::get<1>(t) + (firstBitUint64_t>>1)) == 0)
			{
				hasAny = true;
				parent[std::get<1>(t) + (firstBitUint64_t>>1)] = std::get<1>(t) + (firstBitUint64_t>>1);
			}
			if (std::get<0>(t) != std::numeric_limits<size_t>::max() && std::get<1>(t) != std::numeric_limits<size_t>::max())
			{
				mergedSomewhere.insert(std::get<0>(t));
				mergedSomewhere.insert(std::get<1>(t) + (firstBitUint64_t>>1));
				merge(parent, std::get<0>(t), std::get<1>(t) + (firstBitUint64_t>>1));
			}
		}
		if (!hasAny)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "chunk " << chunki << " doesn't have any" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		bool canSplit = true;
		for (auto pair : parent)
		{
			if (mergedSomewhere.count(pair.first) == 0)
			{
				canSplit = false;
				break;
			}
		}
		if (!canSplit)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "can't split chunk " << chunki << " because of missing merge" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		phmap::flat_hash_map<uint64_t, size_t> clusterToNode;
		size_t nextId = 0;
		for (uint64_t node : mergedSomewhere)
		{
			uint64_t found = find(parent, node);
			if (clusterToNode.count(found) == 1)
			{
				continue;
			}
			clusterToNode[found] = nextId;
			nextId += 1;
		}
		phmap::flat_hash_map<uint64_t, size_t> clusterCoverage;
		for (auto t : predecessorsAndSuccessors)
		{
			if (std::get<0>(t) != std::numeric_limits<size_t>::max())
			{
				clusterCoverage[clusterToNode.at(find(parent, std::get<0>(t)))] += 1;
			}
			else if (std::get<1>(t) != std::numeric_limits<size_t>::max())
			{
				clusterCoverage[clusterToNode.at(find(parent, std::get<1>(t) + (firstBitUint64_t>>1)))] += 1;
			}
		}
		for (auto pair : clusterCoverage)
		{
			if (pair.second < 2)
			{
				canSplit = false;
			}
		}
		if (!canSplit)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
//			std::cerr << "can't split chunk " << chunki << " because of low coverage" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) = nextNum + (std::get<2>(chunksPerRead[occurrencesPerChunk[chunki][j].first][occurrencesPerChunk[chunki][j].second]) & firstBitUint64_t);
			}
			nextNum += 1;
			return;
		}
		assert(clusterToNode.size() >= 1);
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			if (clusterToNode.size() >= 2) countSplitted += 1;
//			std::cerr << "chunk " << chunki << " split to " << clusterToNode.size() << " chunks" << std::endl;
			for (auto t : predecessorsAndSuccessors)
			{
				if (std::get<0>(t) == std::numeric_limits<size_t>::max() && std::get<1>(t) == std::numeric_limits<size_t>::max())
				{
					std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) = std::numeric_limits<size_t>::max();
					continue;
				}
				size_t cluster;
				if (std::get<0>(t) != std::numeric_limits<size_t>::max())
				{
					cluster = clusterToNode.at(find(parent, std::get<0>(t)));
				}
				else
				{
					assert(std::get<1>(t) != std::numeric_limits<size_t>::max());
					cluster = clusterToNode.at(find(parent, std::get<1>(t) + (firstBitUint64_t>>1)));
				}
				std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) = nextNum + cluster + (std::get<2>(chunksPerRead[std::get<2>(t)][std::get<3>(t)]) & firstBitUint64_t);
			}
			nextNum += clusterToNode.size();
			return;
		}
	});
	std::cerr << "between long unitigs resolved " << countSplitted << " chunks" << std::endl;
}

void resolveUnambiguouslyResolvableUnitigs(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "resolving unambiguously resolvable unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<bool> edgesBalanced;
	std::vector<bool> unitigHasContainedPath;
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> tripletsPerUnitig;
	edgesBalanced.resize(graph.unitigLengths.size(), false);
	unitigHasContainedPath.resize(graph.unitigLengths.size(), false);
	tripletsPerUnitig.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (graph.edges.getEdges(std::make_pair(i, true)).size() == graph.edges.getEdges(std::make_pair(i, false)).size())
		{
			if (graph.edges.getEdges(std::make_pair(i, true)).size() >= 2)
			{
				edgesBalanced[i] = true;
			}
		}
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			if (readPaths[i][j].path.size() == 1)
			{
				unitigHasContainedPath[readPaths[i][j].path[0] & maskUint64_t] = true;
				continue;
			}
			assert(readPaths[i][j].path.size() >= 2);
			for (size_t k = 1; k+1 < readPaths[i][j].path.size(); k++)
			{
				uint64_t mid = readPaths[i][j].path[k];
				if (!edgesBalanced[mid & maskUint64_t]) continue;
				if (unitigHasContainedPath[mid & maskUint64_t]) continue;
				uint64_t prev, after;
				if (mid & firstBitUint64_t)
				{
					prev = readPaths[i][j].path[k-1];
					after = readPaths[i][j].path[k+1];
				}
				else
				{
					prev = readPaths[i][j].path[k+1] ^ firstBitUint64_t;
					after = readPaths[i][j].path[k-1] ^ firstBitUint64_t;
				}
				tripletsPerUnitig[mid & maskUint64_t][std::make_pair(prev, after)] += 1;
			}
		}
	}
	std::vector<bool> canResolveUnitig;
	canResolveUnitig.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (!edgesBalanced[i]) continue;
		if (unitigHasContainedPath[i]) continue;
		if (tripletsPerUnitig[i].size() != graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		phmap::flat_hash_set<uint64_t> nodesPrev;
		phmap::flat_hash_set<uint64_t> nodesAfter;
		bool valid = true;
		for (auto triplet : tripletsPerUnitig[i])
		{
			if (triplet.second < 2)
			{
				valid = false;
				break;
			}
			nodesPrev.insert(triplet.first.first);
			nodesAfter.insert(triplet.first.second);
		}
		if (!valid) continue;
		if (nodesPrev.size() != nodesAfter.size()) continue;
		if (nodesPrev.size() != graph.edges.getEdges(std::make_pair(i, true)).size()) continue;
		canResolveUnitig[i] = true;
	}
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> prevToAllele;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> afterToAllele;
	prevToAllele.resize(graph.unitigLengths.size());
	afterToAllele.resize(graph.unitigLengths.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		size_t alleleNum = 0;
		for (auto triplet : tripletsPerUnitig[i])
		{
			uint64_t prev = triplet.first.first;
			uint64_t after = triplet.first.second;
			assert(prevToAllele[i].count(prev) == 0);
			assert(afterToAllele[i].count(after) == 0);
			prevToAllele[i][prev] = alleleNum;
			afterToAllele[i][after] = alleleNum;
			alleleNum += 1;
		}
	}
	assert(readPaths.size() == chunksPerRead.size());
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			assert(chunkBelongsToUnitig[chunk] == std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	phmap::flat_hash_map<size_t, phmap::flat_hash_map<size_t, size_t>> chunkAlleleReplacement;
	size_t nextNum = 0;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> unitigAllelesInThisRead;
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				if (!canResolveUnitig[readPaths[i][j].path[k] & maskUint64_t]) continue;
				size_t alleleBefore = std::numeric_limits<size_t>::max();
				size_t alleleAfter = std::numeric_limits<size_t>::max();
				if (readPaths[i][j].path[k] & firstBitUint64_t)
				{
					if (k > 0) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k-1]);
					if (k+1 < readPaths[i][j].path.size()) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k+1]);
				}
				else
				{
					if (k > 0) alleleAfter = afterToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k-1] ^ firstBitUint64_t);
					if (k+1 < readPaths[i][j].path.size()) alleleBefore = prevToAllele[readPaths[i][j].path[k] & maskUint64_t].at(readPaths[i][j].path[k+1] ^ firstBitUint64_t);
				}
				assert(alleleBefore != std::numeric_limits<size_t>::max() || alleleAfter != std::numeric_limits<size_t>::max());
				assert(alleleBefore == alleleAfter || alleleBefore == std::numeric_limits<size_t>::max() || alleleAfter == std::numeric_limits<size_t>::max());
				size_t allele = alleleBefore;
				if (alleleBefore == std::numeric_limits<size_t>::max()) allele = alleleAfter;
				unitigAllelesInThisRead.emplace_back(readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second, readPaths[i][j].path[k] & maskUint64_t, allele);
			}
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			uint64_t node = std::get<2>(chunksPerRead[i][j]);
			if (NonexistantChunk(node)) continue;
			size_t allele = std::numeric_limits<size_t>::max();
			if ((node & maskUint64_t) < chunkBelongsToUnitig.size() && chunkBelongsToUnitig[node & maskUint64_t] != std::numeric_limits<size_t>::max())
			{
				size_t unitig = chunkBelongsToUnitig[node & maskUint64_t];
				if (canResolveUnitig[unitig])
				{
					for (auto t : unitigAllelesInThisRead)
					{
						if (std::get<2>(t) != unitig) continue;
						if (std::get<0>(t) > std::get<0>(chunksPerRead[i][j])) continue;
						if (std::get<1>(t) < std::get<1>(chunksPerRead[i][j])) continue;
						assert(allele == std::numeric_limits<size_t>::max());
						allele = std::get<3>(t);
					}
					assert(allele != std::numeric_limits<size_t>::max());
				}
			}
			if (chunkAlleleReplacement.count(node & maskUint64_t) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t];
			}
			if (chunkAlleleReplacement.at(node & maskUint64_t).count(allele) == 0)
			{
				chunkAlleleReplacement[node & maskUint64_t][allele] = nextNum;
				nextNum += 1;
			}
			uint64_t replacement = chunkAlleleReplacement.at(node & maskUint64_t).at(allele);
			assert((replacement & firstBitUint64_t) == 0);
			replacement += (node & firstBitUint64_t);
			std::get<2>(chunksPerRead[i][j]) = replacement;
		}
	}
	size_t countUnitigsReplaced = 0;
	size_t countChunksReplaced = 0;
	for (size_t i = 0; i < canResolveUnitig.size(); i++)
	{
		if (!canResolveUnitig[i]) continue;
		countUnitigsReplaced += 1;
		countChunksReplaced += graph.chunksInUnitig[i].size();
	}
	std::cerr << "unambiguously resolvable unitig resolution resolved " << countUnitigsReplaced << " unitigs, " << countChunksReplaced << " chunks" << std::endl;
}

void countReadRepetitiveUnitigs(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "counting read-repetitive unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	std::vector<uint8_t> chunkRepeatCount;
	chunkRepeatCount.resize(chunkBelongsToUnitig.size(), 0);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		phmap::flat_hash_set<size_t> foundHere;
		phmap::flat_hash_set<size_t> repetitiveHere;
		for (auto t : chunksPerRead[i])
		{
			size_t chunk = std::get<2>(t) & maskUint64_t;
			if (foundHere.count(chunk) == 1) repetitiveHere.insert(chunk);
			foundHere.insert(chunk);
		}
		for (size_t chunk : repetitiveHere)
		{
			if (chunk >= chunkRepeatCount.size()) continue;
			if (chunkRepeatCount[chunk] >= 2) continue;
			chunkRepeatCount[chunk] += 1;
		}
	}
	std::vector<bool> unitigRepetitive;
	unitigRepetitive.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < chunkRepeatCount.size(); i++)
	{
		if (chunkRepeatCount[i] < 2) continue;
		if (chunkBelongsToUnitig[i] == std::numeric_limits<size_t>::max()) continue;
		unitigRepetitive[chunkBelongsToUnitig[i]] = true;
	}
	size_t countRepetitive = 0;
	std::ofstream file { "repetitive.txt" };
	for (size_t i = 0; i < unitigRepetitive.size(); i++)
	{
		if (!unitigRepetitive[i]) continue;
		countRepetitive += 1;
		file << i << std::endl;
	}
	std::cerr << countRepetitive << " read-repetitive unitigs" << std::endl;
}

void writeReadChunkSequences(const std::string& filename, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex)
{
	std::cerr << "writing read chunk sequences" << std::endl;
	std::ofstream file { filename };
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			size_t chunk = std::get<2>(t) & maskUint64_t;
			bool fw = std::get<2>(t) & firstBitUint64_t;
			std::string seq = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, i, j);
			std::string name;
			if (i < sequenceIndex.size())
			{
				name = sequenceIndex.getName(i);
			}
			else
			{
				name = sequenceIndex.getName(i - sequenceIndex.size()) + "_bw";
			}
			file << chunk << " " << name << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << seq << std::endl;
		}
	}
}

void writeReadUnitigSequences(const std::string& filename, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "writing read unitig sequences" << std::endl;
	auto chunksPerRead = getBidirectedChunks(rawChunksPerRead);
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::ofstream file { filename };
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				size_t unitig = readPaths[i][j].path[k] & maskUint64_t;
				bool fw = readPaths[i][j].path[k] & firstBitUint64_t;
				std::string seq = sequenceIndex.getSubstring(i, readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second - readPaths[i][j].readPartInPathnode[k].first+1);
				if (!fw) seq = MBG::revCompRaw(seq);
				size_t unitigStart = 0;
				size_t unitigEnd = graph.unitigLengths[unitig];
				if (k == 0) unitigStart += readPaths[i][j].pathLeftClip;
				if (k+1 == readPaths[i][j].path.size()) unitigEnd -= readPaths[i][j].pathRightClip;
				file << unitig << " " << unitigStart << " " << unitigEnd << " " << sequenceIndex.getName(i) << " " << readPaths[i][j].readPartInPathnode[k].first << " " << (readPaths[i][j].readPartInPathnode[k].second - readPaths[i][j].readPartInPathnode[k].first) << " " << seq << std::endl;
			}
		}
	}
}

size_t charToInt(const unsigned char c)
{
	switch(c)
	{
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	}
	assert(false);
	return 0;
}

template <typename CallbackF>
void iterateMinimizers(const std::string& str, const size_t kmerSize, const size_t windowSize, CallbackF callback)
{
	assert(kmerSize * 2 <= sizeof(size_t) * 8);
	assert(kmerSize <= windowSize);
	if (str.size() < kmerSize + windowSize) return;
	const size_t mask = ~(0xFFFFFFFFFFFFFFFF << (kmerSize * 2));
	assert(mask == pow(4, kmerSize)-1);
	thread_local std::vector<std::tuple<size_t, size_t>> window;
	window.clear();
	size_t kmer = 0;
	for (size_t i = 0; i < kmerSize; i++)
	{
		kmer <<= 2;
		kmer |= charToInt(str[i]);
	}
	window.clear();
	window.emplace_back(0, hashFwAndBw(kmer, kmerSize));
	for (size_t i = kmerSize; i < kmerSize + windowSize; i++)
	{
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[i]);
		size_t hashHere = hashFwAndBw(kmer, kmerSize);
		while (!window.empty() && std::get<1>(window.back()) > hashHere) window.pop_back();
		window.emplace_back(i-kmerSize+1, hashHere);
	}
	callback(std::get<0>(window[0]));
	size_t lastCallbackPos = std::get<0>(window[0]);
	for (size_t i = kmerSize+windowSize; i < str.size(); i++)
	{
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[i]);
		auto hashHere = hashFwAndBw(kmer, kmerSize);
		assert(window.size() >= 1);
		bool frontPopped = false;
		if (std::get<0>(window[0]) <= i - windowSize - kmerSize)
		{
			frontPopped = true;
			window.erase(window.begin());
		}
		assert(window.size() == 0 || std::get<0>(window[0]) > i - windowSize - kmerSize);
		if (hashHere < std::get<1>(window[0]))
		{
			for (size_t j = 1; j < window.size(); j++)
			{
				if (std::get<1>(window[j]) == std::get<1>(window[0]))
				{
					assert(std::get<0>(window[j]) > lastCallbackPos);
					callback(std::get<0>(window[j]));
					lastCallbackPos = std::get<0>(window[j]);
				}
			}
		}
		while (window.size() >= 1 && std::get<1>(window.back()) > hashHere) window.pop_back();
		if (window.size() == 0) frontPopped = true;
		window.emplace_back(i-kmerSize+1, hashHere);
		assert(window.size() >= 1);
		if (frontPopped)
		{
			assert(std::get<0>(window[0]) > lastCallbackPos);
			callback(std::get<0>(window[0]));
			lastCallbackPos = std::get<0>(window[0]);
		}
		else
		{
			assert(std::get<0>(window[0]) == lastCallbackPos);
		}
	}
	assert(window.size() >= 1);
	assert(lastCallbackPos == std::get<0>(window[0]));
	for (size_t j = 1; j < window.size(); j++)
	{
		if (std::get<1>(window[j]) == std::get<0>(window[j]))
		{
			assert(std::get<0>(window[j]) > lastCallbackPos);
			callback(std::get<0>(window[j]));
			lastCallbackPos = std::get<0>(window[j]);
		}
	}
}

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getMinimizerBoundedChunksPerRead(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const size_t k, const size_t windowSize)
{
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> result;
	result.resize(sequenceIndex.size()*2);
	const size_t backwardOffset = sequenceIndex.size();
	iterateMultithreaded(0, sequenceIndex.size(), numThreads, [&result, &sequenceIndex, backwardOffset, k, windowSize](const size_t readIndex)
	{
		std::string readSequence = sequenceIndex.getSequence(readIndex);
		std::vector<size_t> minimizerPositions;
		iterateMinimizers(readSequence, k, windowSize, [&minimizerPositions](const size_t pos)
		{
			assert(minimizerPositions.size() == 0 || pos > minimizerPositions.back());
			minimizerPositions.emplace_back(pos);
		});
		for (size_t i = 1; i < minimizerPositions.size(); i++)
		{
			assert(minimizerPositions[i] > minimizerPositions[i-1]);
			result[readIndex].emplace_back(minimizerPositions[i-1], minimizerPositions[i]+k-1, 0);
		};
		for (size_t i = 0; i < result[readIndex].size(); i++)
		{
			result[readIndex+backwardOffset].emplace_back(result[readIndex][result[readIndex].size()-1-i]);
			std::swap(std::get<0>(result[readIndex+backwardOffset].back()), std::get<1>(result[readIndex+backwardOffset].back()));
			std::get<0>(result[readIndex+backwardOffset].back()) = readSequence.size()-1-std::get<0>(result[readIndex+backwardOffset].back());
			std::get<1>(result[readIndex+backwardOffset].back()) = readSequence.size()-1-std::get<1>(result[readIndex+backwardOffset].back());
		}
	});
	return result;
}

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getBetterChunksPerRead(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const size_t k, const size_t windowSize, const size_t middleSkip)
{
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> result;
	result.resize(sequenceIndex.size());
	std::vector<std::string> readFiles { };
	MBG::ReadpartIterator partIterator { 31, 1, MBG::ErrorMasking::No, 1, readFiles, false, "" };
	phmap::flat_hash_map<uint64_t, size_t> hashToNode;
	phmap::flat_hash_set<MBG::HashType> removedHashes;
	std::mutex resultMutex;
	iterateMultithreaded(0, sequenceIndex.size(), numThreads, [&result, &sequenceIndex, &partIterator, &hashToNode, &resultMutex, &removedHashes, k, windowSize, middleSkip](const size_t readIndex)
	{
		std::string readSequence = sequenceIndex.getSequence(readIndex);
		std::vector<std::tuple<size_t, size_t, MBG::HashType>> hashes;
		partIterator.iteratePartsOfRead("", readSequence, [&hashToNode, &hashes, &resultMutex, readIndex, k, windowSize, middleSkip](const MBG::ReadInfo& read, const MBG::SequenceCharType& seq, const MBG::SequenceLengthType& poses, const std::string& raw)
		{
			iterateWindowchunks(seq, k, windowSize, middleSkip, [&hashToNode, &hashes, &resultMutex, &poses, &seq, readIndex, k](const std::vector<uint64_t>& hashPositions)
			{
				assert(hashPositions.size() == 2);
				size_t startPos = hashPositions[0];
				size_t endPos = hashPositions.back() + k - 1;
				size_t realStartPos = poses[startPos];
				size_t realEndPos = poses[endPos+1]-1;
				if (readIndex == 58059)
				{
					std::cerr << "!!! " << readIndex << " " << startPos << " " << endPos << " " << realStartPos << " " << realEndPos << " !!!" << std::endl;
				}
				MBG::SequenceCharType hashableSequence;
				for (size_t i = 0; i < hashPositions.size()/2; i++)
				{
					size_t pos = hashPositions[i];
					hashableSequence.insert(hashableSequence.end(), seq.begin() + pos, seq.begin() + pos + k);
				}
				hashableSequence.emplace_back(10);
				for (size_t i = hashPositions.size()/2; i < hashPositions.size(); i++)
				{
					size_t pos = hashPositions[i];
					MBG::SequenceCharType insertSubstr { seq.begin() + pos, seq.begin() + pos + k };
					//insertSubstr = MBG::revCompRLE(insertSubstr);
					hashableSequence.insert(hashableSequence.end(), insertSubstr.begin(), insertSubstr.end());
				}
				MBG::HashType totalhashfw = MBG::hash(hashableSequence);
				hashes.emplace_back(realStartPos, realEndPos, totalhashfw);
			});
		});
		std::vector<bool> removeHash;
		removeHash.resize(hashes.size(), false);
		phmap::flat_hash_set<MBG::HashType> removeHashesHere;
		for (size_t i = 1; i < hashes.size(); i++)
		{
			assert(std::get<0>(hashes[i]) > std::get<0>(hashes[i-1]) || std::get<1>(hashes[i]) > std::get<1>(hashes[i-1]));
			assert(std::get<0>(hashes[i]) >= std::get<0>(hashes[i-1]));
			assert(std::get<1>(hashes[i]) >= std::get<1>(hashes[i-1]));
			if (std::get<0>(hashes[i]) != std::get<0>(hashes[i-1]) && std::get<1>(hashes[i]) != std::get<1>(hashes[i-1])) continue;
			size_t prevLen = std::get<1>(hashes[i-1]) - std::get<0>(hashes[i-1]);
			size_t currLen = std::get<1>(hashes[i]) - std::get<0>(hashes[i]);
			size_t prevDistFromMid;
			size_t currDistFromMid;
			if (prevLen < middleSkip + windowSize)
			{
				prevDistFromMid = middleSkip + windowSize - prevLen;
			}
			else
			{
				prevDistFromMid = prevLen - (middleSkip + windowSize);
			}
			if (currLen < middleSkip + windowSize)
			{
				currDistFromMid = middleSkip + windowSize - currLen;
			}
			else
			{
				currDistFromMid = currLen - (middleSkip + windowSize);
			}
			if (prevDistFromMid < currDistFromMid)
			{
				if (readIndex == 58059)
				{
					std::cerr << "remove " << std::get<0>(hashes[i]) << " " << std::get<1>(hashes[i]) << " due to prev" << std::endl;
				}
				removeHash[i] = true;
				removeHashesHere.insert(std::get<2>(hashes[i]));
			}
			if (currDistFromMid < prevDistFromMid)
			{
				if (readIndex == 58059)
				{
					std::cerr << "remove " << std::get<0>(hashes[i-1]) << " " << std::get<1>(hashes[i-1]) << " due to next" << std::endl;
				}
				removeHash[i-1] = true;
				removeHashesHere.insert(std::get<2>(hashes[i-1]));
			}
			if (currDistFromMid == prevDistFromMid)
			{
				assert(prevLen != currLen);
				if (prevLen > currLen)
				{
					removeHash[i] = true;
					removeHashesHere.insert(std::get<2>(hashes[i]));
				}
				else
				{
					removeHash[i-1] = true;
					removeHashesHere.insert(std::get<2>(hashes[i-1]));
				}
			}
		}
		for (size_t i = hashes.size()-1; i < hashes.size(); i--)
		{
			if (!removeHash[i]) continue;
			std::swap(hashes[i], hashes.back());
			hashes.pop_back();
		}
		std::sort(hashes.begin(), hashes.end(), [](auto left, auto right) { return std::get<0>(left) < std::get<0>(right); });
		for (size_t i = 1; i < hashes.size(); i++)
		{
			assert(std::get<0>(hashes[i]) > std::get<0>(hashes[i-1]));
			assert(std::get<1>(hashes[i]) > std::get<1>(hashes[i-1]));
		}
		std::lock_guard<std::mutex> lock { resultMutex };
		for (size_t j = 0; j < hashes.size(); j++)
		{
			removedHashes.insert(removeHashesHere.begin(), removeHashesHere.end());
			MBG::HashType totalhashfw = std::get<2>(hashes[j]);
			MBG::HashType totalhashbw = (totalhashfw << (MBG::HashType)64) + (totalhashfw >> (MBG::HashType)64);
//			if (totalhashfw == totalhashbw) continue;
			uint64_t hash;
			bool fw;
			if (totalhashfw <= totalhashbw)
			{
				hash = totalhashfw + 3*(totalhashfw >> 64);
				fw = true;
			}
			else
			{
				hash = totalhashbw + 3*(totalhashbw >> 64);
				fw = false;
			}
			size_t node;
			if (hashToNode.count(hash) == 0)
			{
				node = hashToNode.size();
				hashToNode[hash] = node;
			}
			else
			{
				node = hashToNode.at(hash);
			}
			result[readIndex].emplace_back(std::get<0>(hashes[j]), std::get<1>(hashes[j]), node + (fw ? firstBitUint64_t : 0));
		}
	});
	phmap::flat_hash_set<uint64_t> removedNodes;
	for (MBG::HashType totalhashfw : removedHashes)
	{
		MBG::HashType totalhashbw = (totalhashfw << (MBG::HashType)64) + (totalhashfw >> (MBG::HashType)64);
		uint64_t hash;
		if (totalhashfw <= totalhashbw)
		{
			hash = totalhashfw + 3*(totalhashfw >> 64);
		}
		else
		{
			hash = totalhashbw + 3*(totalhashbw >> 64);
		}
		if (hashToNode.count(hash) == 0) continue;
		uint64_t node = hashToNode.at(hash);
		removedNodes.insert(node);
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		if (result[i].size() == 0) continue;
		if (std::get<1>(result[i].back()) + windowSize >= rawReadLengths[i])
		{
			if (i == 58059)
			{
				std::cerr << "end remove " << std::get<0>(result[i].back()) << " " << std::get<1>(result[i].back()) << " end" << std::endl;
			}
			if (removedNodes.count(std::get<2>(result[i].back()) & maskUint64_t) == 1)
			{
				result[i].pop_back();
			}
		}
		if (result[i].size() == 0) continue;
		if (std::get<0>(result[i][0]) < windowSize)
		{
			if (i == 58059)
			{
				std::cerr << "end remove " << std::get<0>(result[i][0]) << " " << std::get<1>(result[i][0]) << " start" << std::endl;
			}
			if (removedNodes.count(std::get<2>(result[i][0]) & maskUint64_t) == 1)
			{
				result[i].erase(result[i].begin());
			}
		}
	}
	return result;
}

void addMissingPiecesBetweenChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	std::vector<std::tuple<size_t, size_t, size_t>> additionalPieces;
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		maxChunk = std::max(std::get<2>(chunksPerRead[i][0]) & maskUint64_t, maxChunk);
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			maxChunk = std::max(std::get<2>(chunksPerRead[i][j]) & maskUint64_t, maxChunk);
			if (std::get<1>(chunksPerRead[i][j-1]) < std::get<0>(chunksPerRead[i][j]))
			{
				size_t start = std::get<1>(chunksPerRead[i][j-1])-kmerSize+1;
				size_t end = std::get<0>(chunksPerRead[i][j])+kmerSize-1;
				assert(start > std::get<0>(chunksPerRead[i][j-1]));
				assert(end > std::get<1>(chunksPerRead[i][j-1]));
				assert(start < std::get<0>(chunksPerRead[i][j]));
				assert(end < std::get<1>(chunksPerRead[i][j]));
				additionalPieces.emplace_back(start, end, i);
			}
		}
	}
	for (auto t : additionalPieces)
	{
		chunksPerRead[std::get<2>(t)].emplace_back(std::get<0>(t), std::get<1>(t), maxChunk+1);
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
	}
}

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> readChunksFromFakePathsFile(const std::string& filename)
{
	std::cerr << "reading chunks from file " << filename << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> result;
	std::ifstream file { filename };
	if (!file.good())
	{
		std::cerr << "file " << filename << " could not be opened!" << std::endl;
		std::abort();
	}
	size_t countChunks = 0;
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (!file.good()) break;
		std::stringstream sstr { line };
		size_t readindex = 0;
		size_t chunkindex = 0;
		size_t readstart = 0;
		size_t readend = 0;
		sstr >> readindex >> chunkindex >> readstart >> readend;
		assert((readindex >= result.size() && chunkindex == 0) || (readindex == result.size()-1 && chunkindex == result.back().size()));
		assert(readend > readstart);
		while (readindex >= result.size()) result.emplace_back();
		assert(result.back().size() == chunkindex);
		std::string nodestr;
		sstr >> nodestr;
		uint64_t node;
		assert(nodestr.size() >= 1);
		assert(nodestr.size() >= 2 || nodestr == "-");
		if (nodestr == "-")
		{
			node = std::numeric_limits<size_t>::max();
		}
		else if (nodestr[0] == '>')
		{
			node = std::stoull(nodestr.substr(1)) + firstBitUint64_t;
		}
		else
		{
			assert(nodestr[0] == '<');
			node = std::stoull(nodestr.substr(1));
		}
		result.back().emplace_back(readstart, readend, node);
		countChunks += 1;
	}
	std::cerr << "read " << countChunks << " chunks" << std::endl;
	return result;
}

void writeStage(const size_t stage, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize)
{
	writeUnitigGraph("digraph-round" + std::to_string(stage) + ".gfa", "dipaths" + std::to_string(stage) + ".gaf", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
	writeGraph("fakegraph" + std::to_string(stage) + ".gfa", "fakepaths" + std::to_string(stage) + ".txt", chunksPerRead);
	writeBidirectedUnitigGraph("graph-round" + std::to_string(stage) + ".gfa", "paths" + std::to_string(stage) + ".gaf", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
}

void resolveTinyNodesRecklessly(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, const size_t kmerSize)
{
	size_t nextNum = 0;
	std::mutex resultMutex;
	size_t countResolved = 0;
	size_t countResolvedTo = 0;
	std::vector<std::tuple<size_t, size_t, size_t>> replacements;
	iterateChunksByCoverage(chunksPerRead, numThreads, [&nextNum, &resultMutex, &chunksPerRead, &countResolved, &countResolvedTo, &replacements, kmerSize](const size_t chunki, const std::vector<std::vector<std::pair<size_t, size_t>>>& occurrencesPerChunk)
	{
		if (occurrencesPerChunk[chunki].size() <= 2)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				replacements.emplace_back(occurrencesPerChunk[chunki][j].first, occurrencesPerChunk[chunki][j].second, nextNum);
			}
			nextNum += 1;
			return;
		}
		bool isSmall = true;
		for (auto pair : occurrencesPerChunk[chunki])
		{
			auto t = chunksPerRead[pair.first][pair.second];
			if (std::get<1>(t) - std::get<0>(t) >= kmerSize*2)
			{
				isSmall = false;
				break;
			}
		}
		if (!isSmall)
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			for (size_t j = 0; j < occurrencesPerChunk[chunki].size(); j++)
			{
				replacements.emplace_back(occurrencesPerChunk[chunki][j].first, occurrencesPerChunk[chunki][j].second, nextNum);
			}
			nextNum += 1;
			return;
		}
		phmap::flat_hash_map<size_t, size_t> parent;
		for (auto pair : occurrencesPerChunk[chunki])
		{
			size_t prev = std::numeric_limits<size_t>::max();
			size_t next = std::numeric_limits<size_t>::max();
			if (pair.second > 0 && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second-1])))
			{
				prev = std::get<2>(chunksPerRead[pair.first][pair.second-1]);
				if (parent.count(prev) == 0) parent[prev] = prev;
			}
			if (pair.second+1 < chunksPerRead[pair.first].size() && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second+1])))
			{
				next = std::get<2>(chunksPerRead[pair.first][pair.second+1]) + (firstBitUint64_t >> 1);
				if (parent.count(next) == 0) parent[next] = next;
			}
			if (prev != std::numeric_limits<size_t>::max() && next != std::numeric_limits<size_t>::max())
			{
				merge(parent, prev, next);
			}
		}
		phmap::flat_hash_map<size_t, size_t> clusterToNode;
		size_t nextClusterNum = 1;
		for (auto pair : parent)
		{
			if (pair.second != pair.first) continue;
			clusterToNode[pair.second] = nextClusterNum;
			nextClusterNum += 1;
		}
		size_t zeroOffset = 0;
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			zeroOffset = nextNum;
			nextNum += nextClusterNum;
			if (nextClusterNum > 2)
			{
				countResolved += 1;
				countResolvedTo += nextClusterNum-1;
			}
			for (auto pair : occurrencesPerChunk[chunki])
			{
				size_t clusterNum = 0;
				if (pair.second > 0 && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second-1])))
				{
					clusterNum = clusterToNode.at(find(parent, std::get<2>(chunksPerRead[pair.first][pair.second-1])));
				}
				if (pair.second+1 < chunksPerRead[pair.first].size() && !NonexistantChunk(std::get<2>(chunksPerRead[pair.first][pair.second+1])))
				{
					size_t foundClusterNum = clusterToNode.at(find(parent, std::get<2>(chunksPerRead[pair.first][pair.second+1]) + (firstBitUint64_t >> 1)));
					assert(clusterNum == 0 || foundClusterNum == clusterNum);
					clusterNum = foundClusterNum;
				}
				replacements.emplace_back(pair.first, pair.second, zeroOffset + clusterNum);
			}
		}
	});
	for (auto t : replacements)
	{
		std::get<2>(chunksPerRead[std::get<0>(t)][std::get<1>(t)]) = std::get<2>(t) + (std::get<2>(chunksPerRead[std::get<0>(t)][std::get<1>(t)]) & firstBitUint64_t);
	}
	std::cerr << "tiny node resolution resolved " << countResolved << " chunks to " << countResolvedTo << " chunks" << std::endl;
}

void removeTinyProblemNodes(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			assert(std::get<2>(t) & firstBitUint64_t);
			maxChunk = std::max(maxChunk, std::get<2>(t) & maskUint64_t);
		}
	}
	std::vector<size_t> uniquePredecessor;
	std::vector<size_t> uniqueSuccessor;
	uniquePredecessor.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	uniqueSuccessor.resize(maxChunk+1, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1]))) continue;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t prev = std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t;
			size_t curr = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			assert(prev < uniqueSuccessor.size());
			assert(curr < uniquePredecessor.size());
			if (uniqueSuccessor[prev] == std::numeric_limits<size_t>::max())
			{
				uniqueSuccessor[prev] = curr;
			}
			else if (uniqueSuccessor[prev] != curr)
			{
				uniqueSuccessor[prev] = std::numeric_limits<size_t>::max()-1;
			}
			if (uniquePredecessor[curr] == std::numeric_limits<size_t>::max())
			{
				uniquePredecessor[curr] = prev;
			}
			else if (uniquePredecessor[curr] != prev)
			{
				uniquePredecessor[curr] = std::numeric_limits<size_t>::max()-1;
			}
		}
	}
	std::vector<bool> canBeRemoved;
	canBeRemoved.resize(maxChunk+1, true);
	for (size_t i = 0; i < uniquePredecessor.size(); i++)
	{
		if (uniquePredecessor[i] != std::numeric_limits<size_t>::max()-1 && uniqueSuccessor[i] != std::numeric_limits<size_t>::max()-1)
		{
			canBeRemoved[i] = false;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i][0])) && std::get<1>(chunksPerRead[i][0]) - std::get<0>(chunksPerRead[i][0]) >= 2*kmerSize) canBeRemoved[std::get<2>(chunksPerRead[i][0]) & maskUint64_t] = false;
		if (!NonexistantChunk(std::get<2>(chunksPerRead[i].back())) && std::get<1>(chunksPerRead[i].back()) - std::get<0>(chunksPerRead[i].back()) >= 2*kmerSize) canBeRemoved[std::get<2>(chunksPerRead[i].back()) & maskUint64_t] = false;
		for (size_t j = 1; j+1 < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			if (std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) >= 2*kmerSize)
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1])))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j])] = false;
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j+1])))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j])] = false;
				continue;
			}
			if (canBeRemoved[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (canBeRemoved[std::get<2>(chunksPerRead[i][j+1]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
			if (std::get<1>(chunksPerRead[i][j-1])+1 <= std::get<0>(chunksPerRead[i][j+1]))
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
				continue;
			}
		}
	}
	for (size_t i = 0; i < chunksPerRead.size()/2; i++)
	{
		size_t other = i + chunksPerRead.size()/2;
		assert(chunksPerRead[i].size() == chunksPerRead[other].size());
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) - std::get<0>(chunksPerRead[other][chunksPerRead[other].size()-1-j]));
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])))
			{
				assert(NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])));
				continue;
			}
			if (NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])))
			{
				std::cerr << i << " " << other << std::endl;
				std::cerr << j << " " << chunksPerRead[i].size() << std::endl;
				std::cerr << std::get<2>(chunksPerRead[i][j]) << " " << ((std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(chunksPerRead[i][j]) & maskUint64_t) << std::endl;
				std::cerr << std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) << " " << ((std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t) << std::endl;
			}
			assert(!NonexistantChunk(std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j])));
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t] = false;
			}
			if (!canBeRemoved[std::get<2>(chunksPerRead[other][chunksPerRead[other].size()-1-j]) & maskUint64_t])
			{
				canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
			}
		}
	}
	size_t countRemovableChunks = 0;
	for (size_t i = 0; i < canBeRemoved.size(); i++)
	{
		if (!canBeRemoved[i]) continue;
		countRemovableChunks += 1;
	}
	size_t countRemoved = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = chunksPerRead[i].size()-1; j < chunksPerRead[i].size(); j--)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			if (!canBeRemoved[std::get<2>(chunksPerRead[i][j]) & maskUint64_t]) continue;
			std::swap(chunksPerRead[i][j], chunksPerRead[i].back());
			chunksPerRead[i].pop_back();
			countRemoved += 1;
		}
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			assert(std::get<0>(chunksPerRead[i][j]) > std::get<0>(chunksPerRead[i][j-1]));
			assert(std::get<1>(chunksPerRead[i][j]) > std::get<1>(chunksPerRead[i][j-1]));
		}
	}
	std::cerr << countRemovableChunks << " distinct removable chunks, removed " << countRemoved << " total chunks" << std::endl;
}

void makeGraph(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize, const size_t windowSize, const size_t middleSkip, const size_t startStage)
{
	std::cerr << "start at stage " << startStage << std::endl;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead;
	if (startStage > 0)
	{
		chunksPerRead = readChunksFromFakePathsFile("fakepaths" + std::to_string(startStage) + ".txt");
		while (chunksPerRead.size() < sequenceIndex.size()*2) chunksPerRead.emplace_back();
		std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
	switch(startStage)
	{
		case 0:
			std::cerr << "getting chunks from reads" << std::endl;
			chunksPerRead = getMinimizerBoundedChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize);
			//addMissingPiecesBetweenChunks(chunksPerRead, k);
			{
				size_t numChunks = 0;
				for (size_t i = 0; i < chunksPerRead.size(); i++)
				{
					numChunks += chunksPerRead[i].size();
				}
				std::cerr << numChunks << " chunks" << std::endl;
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
			splitPerFirstLastKmers(sequenceIndex, chunksPerRead, kmerSize);
			splitPerLength(chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(1, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 1:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerBaseCounts(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(2, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 2:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerMinHashes(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(3, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 3:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(4, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 4:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentityRoughly(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentity(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(5, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 5:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(6, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 6:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(7, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 7:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerNearestNeighborPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(8, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 8:
			splitPerCorrectedKmerPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(9, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 9:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(10, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 10:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(11, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 11:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(12, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 12:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(13, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 13:
//			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(14, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 14:
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, 2, kmerSize);
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.25, kmerSize);
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.5, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(15, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 15:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(16, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 16:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(17, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 17:
//			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			writeStage(18, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 18:
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			countReadRepetitiveUnitigs(chunksPerRead, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadChunkSequences("sequences-chunk19.txt", rawReadLengths, chunksPerRead, sequenceIndex);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(19, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			writeBidirectedUnitigGraph("graph-dbg-final.gfa", "paths-dbg-final.gaf", "unitigs-dbg-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadUnitigSequences("sequences-dbg-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 19:
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 50000, 100000);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(20, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 20:
			writeBidirectedUnitigGraph("graph-final.gfa", "paths-final.gaf", "unitigs-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadUnitigSequences("sequences-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
}

void readFilesAndAddToSequenceIndex(const std::string& filename, FastaCompressor::CompressedStringIndex& sequenceIndex, std::vector<size_t>& readBasepairLengths, const size_t numThreads)
{
	std::vector<std::tuple<size_t, std::string*, std::string*>> sequenceStack;
	std::vector<std::thread> threads;
	std::atomic<bool> allDone;
	std::mutex stackMutex;
	size_t nextNum = 0;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&allDone, &stackMutex, &sequenceStack, &sequenceIndex]()
		{
			while (true)
			{
				std::tuple<size_t, std::string*, std::string*> item { std::numeric_limits<size_t>::max(), nullptr, nullptr };
				{
					std::lock_guard<std::mutex> lock { stackMutex };
					if (sequenceStack.size() == 0)
					{
						if (allDone) break;
					}
					else
					{
						item = sequenceStack.back();
						sequenceStack.pop_back();
					}
				}
				if (std::get<0>(item) == std::numeric_limits<size_t>::max())
				{
					assert(std::get<1>(item) == nullptr);
					assert(std::get<2>(item) == nullptr);
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
					continue;
				}
				sequenceIndex.addString(std::get<0>(item), *std::get<1>(item), *std::get<2>(item));
				delete std::get<1>(item);
				delete std::get<2>(item);
			}
		});
	}
	FastQ::streamFastqFromFile(filename, false, [&sequenceStack, &stackMutex, &nextNum, &readBasepairLengths](FastQ& read)
	{
		std::tuple<size_t, std::string*, std::string*> item;
		readBasepairLengths.emplace_back(read.sequence.size());
		std::get<0>(item) = nextNum;
		std::get<1>(item) = new std::string;
		std::get<2>(item) = new std::string;
		std::swap(*std::get<1>(item), read.seq_id);
		std::swap(*std::get<2>(item), read.sequence);
		nextNum += 1;
		while (true)
		{
			{
				std::lock_guard<std::mutex> lock { stackMutex };
				if (sequenceStack.size() < 100)
				{
					sequenceStack.emplace_back(item);
					break;
				}
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			continue;
		}
	});
	allDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

int main(int argc, char** argv)
{
	const size_t numThreads = std::stoi(argv[1]);
	const size_t k = std::stoull(argv[2]);
	const size_t windowSize = std::stoull(argv[3]);
	const size_t middleSkip = std::stoull(argv[4]);
	const double approxOneHapCoverage = std::stod(argv[5]);
	mismatchFraction = std::stod(argv[6]); // try 2-3x average error rate
	const size_t startStage = std::stoull(argv[7]);
	std::vector<std::string> readFiles;
	for (size_t i = 8; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	FastaCompressor::CompressedStringIndex sequenceIndex { 5, 100 };
	std::vector<size_t> readBasepairLengths;
	for (auto file : readFiles)
	{
		std::cerr << "reading from file " << file << std::endl;
		readFilesAndAddToSequenceIndex(file, sequenceIndex, readBasepairLengths, numThreads);
	}
	sequenceIndex.removeConstructionVariables();
	size_t lastReal = readBasepairLengths.size();
	for (size_t i = 0; i < lastReal; i++)
	{
		readBasepairLengths.emplace_back(readBasepairLengths[i]);
	}
	std::cerr << sequenceIndex.size() << " reads" << std::endl;
	makeGraph(sequenceIndex, readBasepairLengths, numThreads, approxOneHapCoverage, k, windowSize, middleSkip, startStage);
}
