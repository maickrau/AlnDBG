#include <iostream>
#include "ChunkUnitigGraph.h"
#include "Common.h"
#include "UnionFind.h"

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
			assert(j == 0 || std::get<0>(chunksPerRead[i][j]) >= std::get<0>(chunksPerRead[i][j-1]));
			assert(j == 0 || std::get<1>(chunksPerRead[i][j]) >= std::get<1>(chunksPerRead[i][j-1]));
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
			assert(j == 0 || std::get<0>(chunksPerRead[i][j]) >= std::get<0>(chunksPerRead[i][j-1]));
			assert(j == 0 || std::get<1>(chunksPerRead[i][j]) >= std::get<1>(chunksPerRead[i][j-1]));
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
			assert(std::get<0>(chunksPerRead[i][j]) >= std::get<0>(chunksPerRead[i][j-1]));
			assert(std::get<1>(chunksPerRead[i][j]) >= std::get<1>(chunksPerRead[i][j-1]));
			auto prev = std::get<2>(chunksPerRead[i][j-1]);
			auto curr = std::get<2>(chunksPerRead[i][j]);
			size_t overlap = 0;
			if (std::get<1>(chunksPerRead[i][j-1]) > std::get<0>(chunksPerRead[i][j])) overlap = std::get<1>(chunksPerRead[i][j-1]) - std::get<0>(chunksPerRead[i][j]);
			auto pairkey = canon(std::make_pair(prev & maskUint64_t, prev & firstBitUint64_t), std::make_pair(curr & maskUint64_t, curr & firstBitUint64_t));
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
			assert(std::get<0>(chunksPerRead[i][j]) >= std::get<0>(chunksPerRead[i][j-1]));
			assert(std::get<1>(chunksPerRead[i][j]) >= std::get<1>(chunksPerRead[i][j-1]));
			auto prev = std::get<2>(chunksPerRead[i][j-1]);
			auto curr = std::get<2>(chunksPerRead[i][j]);
			auto pairkey = canon(std::make_pair(prev & maskUint64_t, prev & firstBitUint64_t), std::make_pair(curr & maskUint64_t, curr & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
			edgeCoverage[key] += 1;
		}
	}
	return edgeCoverage;
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
		assert(allowedEdges.hasEdge(reverse(next), reverse(pairpos)));
		if (allowedEdges.getEdges(reverse(next)).size() != 1) break;
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
		assert(allowedEdges.hasEdge(reverse(next), reverse(pairpos)));
		if (allowedEdges.getEdges(reverse(next)).size() != 1) break;
		if (next.first == pairpos.first) break;
		if (next.first == (unitigStart & maskUint64_t)) break;
		pos = next.first + (next.second ? firstBitUint64_t : 0);
		assert(allowedNode[pos & maskUint64_t]);
	}
	return result;
}

std::tuple<std::vector<std::vector<uint64_t>>, std::vector<size_t>, std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>, std::vector<std::vector<std::pair<size_t, size_t>>>> getUnitigs(const std::vector<bool>& allowedNode, const SparseEdgeContainer& allowedEdges, const std::vector<std::vector<size_t>>& lengths, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps)
{
	std::vector<std::vector<uint64_t>> unitigs;
	std::vector<bool> checked;
	checked.resize(allowedNode.size(), false);
	std::vector<std::tuple<uint64_t, size_t, size_t, size_t>> chunkLocationInUnitig;
	std::vector<std::vector<std::pair<size_t, size_t>>> unitigChunkBreakpointPositions;
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
	unitigChunkBreakpointPositions.resize(unitigs.size());
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		unitigLengths.emplace_back(0);
		unitigChunkBreakpointPositions[i].reserve(unitigs[i].size());
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
				auto pairkey = canon(fromnode, tonode);
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
			unitigChunkBreakpointPositions[i].emplace_back(unitigLengths.back(), unitigLengths.back() + length);
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
	return std::make_tuple(unitigs, unitigLengths, chunkLocationInUnitig, unitigChunkBreakpointPositions);
}

std::vector<std::vector<UnitigPath>> getUnitigPaths(const ChunkUnitigGraph& graph, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<bool>& allowedNode, const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>& chunkLocationInUnitig)
{
	std::vector<std::vector<UnitigPath>> result;
	result.resize(chunksPerRead.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<uint64_t> path;
		std::vector<std::pair<size_t, size_t>> readPartInPathnode;
		std::vector<std::pair<size_t, size_t>> chunksInPathnode;
		uint64_t currentUnitig = std::numeric_limits<uint64_t>::max();
		size_t currentUnitigIndex = std::numeric_limits<size_t>::max();
		size_t currentPathLeftClipBases = std::numeric_limits<size_t>::max();
		size_t currentPathRightClipBases = std::numeric_limits<size_t>::max();
		size_t currentPathLeftClipChunks = std::numeric_limits<size_t>::max();
		size_t currentPathRightClipChunks = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j])) || !allowedNode[std::get<2>(chunksPerRead[i][j]) & maskUint64_t])
			{
				if (path.size() > 0)
				{
					result[i].emplace_back();
					result[i].back().path = path;
					result[i].back().pathLeftClipBases = currentPathLeftClipBases;
					result[i].back().pathRightClipBases = currentPathRightClipBases;
					result[i].back().readPartInPathnode = readPartInPathnode;
					result[i].back().chunksInPathnode = chunksInPathnode;
					result[i].back().pathLeftClipChunks = currentPathLeftClipChunks;
					result[i].back().pathRightClipChunks = currentPathRightClipChunks;
				}
				path.clear();
				readPartInPathnode.clear();
				chunksInPathnode.clear();
				currentUnitig = std::numeric_limits<size_t>::max();
				currentUnitigIndex = 0;
				currentPathLeftClipBases = 0;
				currentPathRightClipBases = 0;
				currentPathLeftClipChunks = 0;
				currentPathRightClipChunks = 0;
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
				currentPathLeftClipChunks = unitigIndex;
				currentPathRightClipChunks = graph.chunksInUnitig[unitig & maskUint64_t].size() - 1 - unitigIndex;
				currentPathLeftClipBases = unitigStartPos;
				currentPathRightClipBases = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
				readPartInPathnode.emplace_back(readStart, readEnd);
				chunksInPathnode.emplace_back(j, j);
				continue;
			}
			if (currentUnitig == unitig && currentUnitigIndex+1 == unitigIndex)
			{
				assert(readEnd >= readPartInPathnode.back().second);
				assert(graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos < currentPathRightClipBases);
				currentUnitigIndex = unitigIndex;
				readPartInPathnode.back().second = readEnd;
				chunksInPathnode.back().second = j;
				currentPathRightClipChunks = graph.chunksInUnitig[unitig & maskUint64_t].size() - 1 - unitigIndex;
				currentPathRightClipBases = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
				continue;
			}
			if (currentUnitigIndex+1 == unitigs[currentUnitig & maskUint64_t].size() && unitigIndex == 0)
			{
				std::pair<size_t, bool> fromUnitig { currentUnitig & maskUint64_t, currentUnitig & firstBitUint64_t };
				std::pair<size_t, bool> thisUnitig { unitig & maskUint64_t, unitig & firstBitUint64_t };
				if (graph.edges.hasEdge(fromUnitig, thisUnitig))
				{
					assert(readEnd >= readPartInPathnode.back().second);
					path.emplace_back(unitig);
					readPartInPathnode.emplace_back(readStart, readEnd);
					chunksInPathnode.emplace_back(j, j);
					currentUnitig = unitig;
					currentUnitigIndex = unitigIndex;
					currentPathRightClipChunks = graph.chunksInUnitig[unitig & maskUint64_t].size() - 1 - unitigIndex;
					currentPathRightClipBases = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
					continue;
				}
			}
			result[i].emplace_back();
			result[i].back().path = path;
			result[i].back().pathLeftClipBases = currentPathLeftClipBases;
			result[i].back().pathRightClipBases = currentPathRightClipBases;
			result[i].back().pathLeftClipChunks = currentPathLeftClipChunks;
			result[i].back().pathRightClipChunks = currentPathRightClipChunks;
			result[i].back().readPartInPathnode = readPartInPathnode;
			result[i].back().chunksInPathnode = chunksInPathnode;
			path.clear();
			readPartInPathnode.clear();
			chunksInPathnode.clear();
			path.emplace_back(unitig);
			currentUnitig = unitig;
			currentUnitigIndex = unitigIndex;
			currentPathLeftClipChunks = unitigIndex;
			currentPathRightClipChunks = graph.chunksInUnitig[unitig & maskUint64_t].size() - 1 - unitigIndex;
			currentPathLeftClipBases = unitigStartPos;
			currentPathRightClipBases = graph.unitigLengths[unitig & maskUint64_t] - unitigEndPos;
			readPartInPathnode.emplace_back(readStart, readEnd);
			chunksInPathnode.emplace_back(j, j);
		}
		if (path.size() >= 1)
		{
			result[i].emplace_back();
			result[i].back().path = path;
			result[i].back().pathLeftClipBases = currentPathLeftClipBases;
			result[i].back().pathRightClipBases = currentPathRightClipBases;
			result[i].back().pathLeftClipChunks = currentPathLeftClipChunks;
			result[i].back().pathRightClipChunks = currentPathRightClipChunks;
			result[i].back().readPartInPathnode = readPartInPathnode;
			result[i].back().chunksInPathnode = chunksInPathnode;
		}
	}
	return result;
}

std::vector<ConsensusString> getUnitigConsensuses(const ChunkUnitigGraph& unitigGraph, const std::vector<std::vector<UnitigPath>>& readPaths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<size_t>>& chunkLengths, const std::vector<size_t>& chunkCoverages, const std::vector<std::tuple<uint64_t, size_t, size_t, size_t>>& chunkLocationInUnitig, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps, const size_t numThreads, const size_t kmerSize)
{
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "getting unitig consensuses" << std::endl;
	std::vector<phmap::flat_hash_map<size_t, std::vector<std::tuple<size_t, bool, size_t>>>> readAnchorPoints; // unitigpos -> read, fw, readpos
	readAnchorPoints.resize(unitigGraph.chunksInUnitig.size());
	std::vector<std::vector<size_t>> chunkEndToKmerIndex;
	chunkEndToKmerIndex.resize(unitigGraph.chunksInUnitig.size());
	{
		assert(chunksPerRead.size() == sequenceIndex.size()); // bidirected chunks! not directed as usual
		size_t maxChunk = 0;
		for (size_t i = 0; i < readPaths.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
				maxChunk = std::max(maxChunk, std::get<2>(chunksPerRead[i][j]) & maskUint64_t);
			}
		}
		maxChunk += 1;
		std::vector<size_t> chunkForwardNeighbor;
		std::vector<size_t> chunkBackwardNeighbor;
		std::vector<bool> chunkStartPosEqualToBack;
		std::vector<bool> chunkEndPosEqualToFront;
		chunkForwardNeighbor.resize(maxChunk, std::numeric_limits<size_t>::max());
		chunkBackwardNeighbor.resize(maxChunk, std::numeric_limits<size_t>::max());
		chunkStartPosEqualToBack.resize(maxChunk, true);
		chunkEndPosEqualToFront.resize(maxChunk, true);
		for (size_t i = 0; i < unitigGraph.chunksInUnitig.size(); i++)
		{
			for (size_t j = 1; j < unitigGraph.chunksInUnitig[i].size(); j++)
			{
				if (unitigGraph.chunksInUnitig[i][j-1] & firstBitUint64_t)
				{
					chunkForwardNeighbor[unitigGraph.chunksInUnitig[i][j-1] & maskUint64_t] = unitigGraph.chunksInUnitig[i][j];
				}
				else
				{
					chunkBackwardNeighbor[unitigGraph.chunksInUnitig[i][j-1] & maskUint64_t] = unitigGraph.chunksInUnitig[i][j] ^ firstBitUint64_t;
				}
				if (unitigGraph.chunksInUnitig[i][j] & firstBitUint64_t)
				{
					chunkBackwardNeighbor[unitigGraph.chunksInUnitig[i][j] & maskUint64_t] = unitigGraph.chunksInUnitig[i][j-1];
				}
				else
				{
					chunkForwardNeighbor[unitigGraph.chunksInUnitig[i][j] & maskUint64_t] = unitigGraph.chunksInUnitig[i][j-1] ^ firstBitUint64_t;
				}
			}
		}
		std::vector<bool> chunkCheckedForward;
		std::vector<bool> chunkCheckedBackward;
		chunkCheckedForward.resize(maxChunk, false);
		chunkCheckedBackward.resize(maxChunk, false);
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 1; j < chunksPerRead[i].size(); j++)
			{
				if (std::get<2>(chunksPerRead[i][j-1]) & firstBitUint64_t)
				{
					if (chunkForwardNeighbor[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] != std::numeric_limits<size_t>::max())
					{
						if (chunkForwardNeighbor[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] == std::get<2>(chunksPerRead[i][j]))
						{
							chunkCheckedForward[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = true;
							if (std::get<1>(chunksPerRead[i][j-1]) != std::get<1>(chunksPerRead[i][j]))
							{
								chunkEndPosEqualToFront[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;
							}
						}
					}
				}
				else
				{
					if (chunkBackwardNeighbor[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] != std::numeric_limits<size_t>::max())
					{
						if (chunkBackwardNeighbor[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] == (std::get<2>(chunksPerRead[i][j]) ^ firstBitUint64_t))
						{
							chunkCheckedBackward[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = true;
							if (std::get<1>(chunksPerRead[i][j-1]) != std::get<1>(chunksPerRead[i][j]))
							{
								chunkStartPosEqualToBack[std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t] = false;
							}
						}
					}
				}
				if (std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t)
				{
					if (chunkBackwardNeighbor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] != std::numeric_limits<size_t>::max())
					{
						if (chunkBackwardNeighbor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] == std::get<2>(chunksPerRead[i][j-1]))
						{
							chunkCheckedBackward[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = true;
							if (std::get<0>(chunksPerRead[i][j-1]) != std::get<0>(chunksPerRead[i][j]))
							{
								chunkStartPosEqualToBack[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
							}
						}
					}
				}
				else
				{
					if (chunkForwardNeighbor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] != std::numeric_limits<size_t>::max())
					{
						if (chunkForwardNeighbor[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] == (std::get<2>(chunksPerRead[i][j-1]) ^ firstBitUint64_t))
						{
							chunkCheckedForward[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = true;
							if (std::get<0>(chunksPerRead[i][j-1]) != std::get<0>(chunksPerRead[i][j]))
							{
								chunkEndPosEqualToFront[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
							}
						}
					}
				}
			}
		}
		for (size_t i = 0; i < maxChunk; i++)
		{
			assert(chunkForwardNeighbor[i] == std::numeric_limits<size_t>::max() || chunkCheckedForward[i]);
			assert(chunkBackwardNeighbor[i] == std::numeric_limits<size_t>::max() || chunkCheckedBackward[i]);
			assert(chunkForwardNeighbor[i] != std::numeric_limits<size_t>::max() || !chunkCheckedForward[i]);
			assert(chunkBackwardNeighbor[i] != std::numeric_limits<size_t>::max() || !chunkCheckedBackward[i]);
		}
		for (size_t i = 0; i < unitigGraph.chunksInUnitig.size(); i++)
		{
			assert(unitigGraph.unitigChunkBreakpointPositions[i].back().second == unitigGraph.unitigLengths[i]);
			chunkEndToKmerIndex[i].emplace_back(unitigGraph.unitigChunkBreakpointPositions[i][0].second - kmerSize + 1);
			assert(unitigGraph.unitigChunkBreakpointPositions[i][0].second <= unitigGraph.unitigLengths[i]);
			size_t index = unitigGraph.unitigChunkBreakpointPositions[i][0].second - kmerSize + 1;
			for (size_t j = 1; j < unitigGraph.chunksInUnitig[i].size(); j++)
			{
				if (unitigGraph.chunksInUnitig[i][j-1] & firstBitUint64_t)
				{
					if (chunkEndPosEqualToFront[unitigGraph.chunksInUnitig[i][j-1] & maskUint64_t])
					{
						chunkEndToKmerIndex[i].emplace_back(index);
						continue;
					}
				}
				else
				{
					if (chunkStartPosEqualToBack[unitigGraph.chunksInUnitig[i][j-1] & maskUint64_t])
					{
						chunkEndToKmerIndex[i].emplace_back(index);
						continue;
					}
				}
				assert(unitigGraph.unitigChunkBreakpointPositions[i][j].second > kmerSize - 1);
				assert(unitigGraph.unitigChunkBreakpointPositions[i][j].second - kmerSize + 1 > index);
				assert(unitigGraph.unitigChunkBreakpointPositions[i][j].second <= unitigGraph.unitigLengths[i]);
				index = unitigGraph.unitigChunkBreakpointPositions[i][j].second - kmerSize + 1;
				chunkEndToKmerIndex[i].emplace_back(index);
			}
			assert(chunkEndToKmerIndex[i].size() >= 1);
			assert(chunkEndToKmerIndex[i].back() <= unitigGraph.unitigLengths[i] - kmerSize + 1);
			chunkEndToKmerIndex[i].back() = unitigGraph.unitigLengths[i] - kmerSize + 1;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			if (NonexistantChunk(std::get<2>(t))) continue;
			size_t readStart = std::get<0>(t);
			size_t readEnd = std::get<1>(t) - kmerSize + 1;
			uint64_t node = std::get<2>(t);
			assert((node & maskUint64_t) < chunkLocationInUnitig.size());
			uint64_t unitig = std::get<0>(chunkLocationInUnitig[node & maskUint64_t]);
			if (unitig == std::numeric_limits<size_t>::max()) continue;
			size_t unitigPos = std::get<1>(chunkLocationInUnitig[node & maskUint64_t]);
			assert(unitigPos < unitigGraph.chunksInUnitig[unitig & maskUint64_t].size());
			assert((unitig & maskUint64_t) < readAnchorPoints.size());
			bool fw = true;
			if (node & firstBitUint64_t)
			{
			}
			else
			{
				std::swap(readEnd, readStart);
				fw = false;
			}
			if (unitig & firstBitUint64_t)
			{
				if (unitigPos == 0)
				{
					readAnchorPoints[unitig & maskUint64_t][0].emplace_back(i, fw, readStart);
				}
				readAnchorPoints[unitig & maskUint64_t][chunkEndToKmerIndex[unitig & maskUint64_t][unitigPos]].emplace_back(i, fw, readEnd);
			}
			else
			{
				fw = !fw;
				unitigPos = unitigGraph.chunksInUnitig[unitig & maskUint64_t].size() - 1 - unitigPos;
				readAnchorPoints[unitig & maskUint64_t][chunkEndToKmerIndex[unitig & maskUint64_t][unitigPos]].emplace_back(i, fw, readStart);
				if (unitigPos == 0)
				{
					readAnchorPoints[unitig & maskUint64_t][0].emplace_back(i, fw, readEnd);
				}
			}
		}
	}
	std::vector<ConsensusString> result;
	result.resize(unitigGraph.unitigLengths.size());
	iterateMultithreaded(0, unitigGraph.unitigLengths.size(), numThreads, [&readAnchorPoints, &sequenceIndex, &result, &unitigGraph, kmerSize](const size_t unitigi)
	{
		std::vector<size_t> anchorPositions;
		for (const auto& pair : readAnchorPoints[unitigi])
		{
			anchorPositions.emplace_back(pair.first);
		}
		assert(anchorPositions.size() >= 2);
		std::sort(anchorPositions.begin(), anchorPositions.end());
		assert(anchorPositions[0] == 0);
		assert(anchorPositions.back() == unitigGraph.unitigLengths[unitigi] - kmerSize+1);
		for (size_t i = 1; i < anchorPositions.size(); i++)
		{
			std::string seq = getConsensusSequence(sequenceIndex, readAnchorPoints[unitigi], anchorPositions[i-1], anchorPositions[i], kmerSize);
			if (seq.size() == 0)
			{
				result[unitigi].Nchunks.emplace_back(result[unitigi].bases.size(), anchorPositions[i]-anchorPositions[i-1]);
				for(size_t j = anchorPositions[i-1]; j < anchorPositions[i]; j++)
				{
					result[unitigi].bases.emplace_back(0);
				}
				continue;
			}
			assert(seq.size() >= kmerSize+1);
			size_t startPos = 0;
			if (i != 1 && (result[unitigi].Nchunks.size() == 0 || result[unitigi].Nchunks.back().first + result[unitigi].Nchunks.back().second + kmerSize < result[unitigi].bases.size()))
			{
				assert(result[unitigi].bases.size() >= kmerSize);
				std::string previousEnd = result[unitigi].bases.substr(result[unitigi].bases.size()-kmerSize, kmerSize);
				std::string thisStart = seq.substr(0, kmerSize);
				assert(previousEnd == thisStart);
				startPos = kmerSize;
			}
			for (size_t j = startPos; j < seq.size(); j++)
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
	std::vector<std::vector<std::pair<size_t, size_t>>> unitigChunkBreakpointPositions;
	std::vector<std::tuple<uint64_t, size_t, size_t, size_t>> chunkLocationInUnitig; // unitig, index, startpos, endpos. reverse also reverses index and poses.
	std::tie(unitigs, unitigLengths, chunkLocationInUnitig, unitigChunkBreakpointPositions) = getUnitigs(allowedNode, allowedEdges, lengths, edgeOverlaps);
	{
		allowedEdges = filterOutZEdges(unitigs, unitigLengths, allowedEdges, lengths, coverages, edgeCoverage, chunkLocationInUnitig, approxOneHapCoverage);
		std::tie(unitigs, unitigLengths, chunkLocationInUnitig, unitigChunkBreakpointPositions) = getUnitigs(allowedNode, allowedEdges, lengths, edgeOverlaps);
	}
	std::swap(result.unitigChunkBreakpointPositions, unitigChunkBreakpointPositions);
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
			auto nodepairkey = canon(lastNode, edge);
			std::pair<uint64_t, uint64_t> nodekey { nodepairkey.first.first + (nodepairkey.first.second ? firstBitUint64_t : 0), nodepairkey.second.first + (nodepairkey.second.second ? firstBitUint64_t : 0) };
			auto unitigpairkey = canon(std::make_pair(i, true), std::make_pair(targetUnitig & maskUint64_t, targetUnitig & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> unitigkey { unitigpairkey.first.first + (unitigpairkey.first.second ? firstBitUint64_t : 0), unitigpairkey.second.first + (unitigpairkey.second.second ? firstBitUint64_t : 0) };
			result.edges.addEdge(unitigpairkey.first, unitigpairkey.second);
			result.edges.addEdge(reverse(unitigpairkey.second), reverse(unitigpairkey.first));
			result.edgeOverlaps[unitigkey] = edgeOverlaps.at(nodekey);
			result.edgeCoverages[unitigkey] = edgeCoverage.at(nodekey);
		}
		for (auto edge : allowedEdges.getEdges(firstNode))
		{
			assert(!NonexistantChunk(edge.first + (edge.second ? firstBitUint64_t : 0)));
			uint64_t targetUnitig = std::get<0>(chunkLocationInUnitig[edge.first]);
			if (!edge.second) targetUnitig ^= firstBitUint64_t;
			auto nodepairkey = canon(firstNode, edge);
			std::pair<uint64_t, uint64_t> nodekey { nodepairkey.first.first + (nodepairkey.first.second ? firstBitUint64_t : 0), nodepairkey.second.first + (nodepairkey.second.second ? firstBitUint64_t : 0) };
			auto unitigpairkey = canon(std::make_pair(i, false), std::make_pair(targetUnitig & maskUint64_t, targetUnitig & firstBitUint64_t));
			std::pair<uint64_t, uint64_t> unitigkey { unitigpairkey.first.first + (unitigpairkey.first.second ? firstBitUint64_t : 0), unitigpairkey.second.first + (unitigpairkey.second.second ? firstBitUint64_t : 0) };
			result.edges.addEdge(unitigpairkey.first, unitigpairkey.second);
			result.edges.addEdge(reverse(unitigpairkey.second), reverse(unitigpairkey.first));
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
		consensuses = getUnitigConsensuses(result, paths, sequenceIndex, chunksPerRead, lengths, coverages, chunkLocationInUnitig, edgeOverlaps, numThreads, kmerSize);
	}
	assert(result.unitigChunkBreakpointPositions.size() == result.unitigLengths.size());
	for (size_t i = 0; i < result.unitigLengths.size(); i++)
	{
		assert(result.unitigChunkBreakpointPositions[i].size() == result.chunksInUnitig[i].size());
		for (size_t j = 0; j < result.unitigChunkBreakpointPositions[i].size(); j++)
		{
			assert(result.unitigChunkBreakpointPositions[i][j].first < result.unitigChunkBreakpointPositions[i][j].second);
			assert(result.unitigChunkBreakpointPositions[i][j].second <= result.unitigLengths[i]);
			if (j > 0)
			{
				assert(result.unitigChunkBreakpointPositions[i][j].first >= result.unitigChunkBreakpointPositions[i][j-1].first);
				assert(result.unitigChunkBreakpointPositions[i][j].second >= result.unitigChunkBreakpointPositions[i][j-1].second);
			}
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

std::pair<size_t, size_t> firstAndLastMinimizerWithinSpan(const std::vector<size_t>& minimizerPositions, const size_t start, const size_t end, const size_t kmerSize)
{
	assert(end >= start + kmerSize-1);
	size_t first = std::numeric_limits<size_t>::max();
	size_t last = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < minimizerPositions.size(); i++)
	{
		if (minimizerPositions[i] == start)
		{
			assert(first == std::numeric_limits<size_t>::max());
			first = i;
		}
		if (minimizerPositions[i]+kmerSize-1 == end)
		{
			assert(last == std::numeric_limits<size_t>::max());
			last = i;
			break;
		}
	}
	assert(first != std::numeric_limits<size_t>::max());
	assert(last != std::numeric_limits<size_t>::max());
	assert(last >= first);
	assert(last < minimizerPositions.size());
	return std::make_pair(first, last);
}

std::pair<std::pair<size_t, size_t>, bool> findBidirected(phmap::flat_hash_map<std::pair<size_t, size_t>, std::pair<std::pair<size_t, size_t>, bool>>& bidirectedParent, std::pair<std::pair<size_t, size_t>, bool> key)
{
	while (true)
	{
		auto parent = bidirectedParent.at(key.first);
		auto parent2 = bidirectedParent.at(parent.first);
		if (parent2.first == parent.first) break;
		bidirectedParent[key.first] = parent2;
		if (!parent.second) bidirectedParent[key.first].second = !bidirectedParent[key.first].second;
	}
	auto result = bidirectedParent.at(key.first);
	if (!key.second) result.second = !result.second;
	return result;
}

void mergeBidirected(phmap::flat_hash_map<std::pair<size_t, size_t>, std::pair<std::pair<size_t, size_t>, bool>>& bidirectedParent, std::pair<std::pair<size_t, size_t>, bool> left, std::pair<std::pair<size_t, size_t>, bool> right)
{
	auto leftp = findBidirected(bidirectedParent, left);
	auto rightp = findBidirected(bidirectedParent, right);
	assert(bidirectedParent.at(leftp.first).first == leftp.first);
	assert(bidirectedParent.at(rightp.first).first == rightp.first);
	if (leftp.first == rightp.first) return;
	if (leftp.second == rightp.second)
	{
		bidirectedParent[rightp.first] = std::make_pair(leftp.first, true);
	}
	else
	{
		bidirectedParent[rightp.first] = std::make_pair(leftp.first, false);
	}
}

std::pair<std::vector<size_t>, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> getSmallChunkSizesAndOverlaps(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& fixedChunksPerRead, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const size_t kmerSize)
{
	std::vector<size_t> countSmallchunksPerChunk;
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> smallChunkOverlap;
	for (size_t read = 0; read < fixedChunksPerRead.size(); read++)
	{
		for (size_t j = 0; j < fixedChunksPerRead[read].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(fixedChunksPerRead[read][j]))) continue;
			size_t chunk = std::get<2>(fixedChunksPerRead[read][j]) & maskUint64_t;
			bool fw = std::get<2>(fixedChunksPerRead[read][j]) & firstBitUint64_t;
			std::pair<size_t, size_t> firstAndLast = firstAndLastMinimizerWithinSpan(minimizerPositionsPerRead[read], std::get<0>(fixedChunksPerRead[read][j]), std::get<1>(fixedChunksPerRead[read][j]), kmerSize);
			size_t countMinimizers = firstAndLast.second - firstAndLast.first + 1;
			assert(countMinimizers >= 2);
			assert(countMinimizers <= minimizerPositionsPerRead[read].size());
			size_t countSmallchunks = countMinimizers - 1;
			while (countSmallchunksPerChunk.size() <= chunk) countSmallchunksPerChunk.emplace_back(std::numeric_limits<size_t>::max());
			assert(countSmallchunksPerChunk[chunk] == countSmallchunks || countSmallchunksPerChunk[chunk] == std::numeric_limits<size_t>::max());
			countSmallchunksPerChunk[chunk] = countSmallchunks;
		}
		for (size_t j = 1; j < fixedChunksPerRead[read].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(fixedChunksPerRead[read][j-1]))) continue;
			if (NonexistantChunk(std::get<2>(fixedChunksPerRead[read][j]))) continue;
			size_t currchunk = std::get<2>(fixedChunksPerRead[read][j]) & maskUint64_t;
			size_t prevchunk = std::get<2>(fixedChunksPerRead[read][j-1]) & maskUint64_t;
			bool currfw = std::get<2>(fixedChunksPerRead[read][j]) & firstBitUint64_t;
			bool prevfw = std::get<2>(fixedChunksPerRead[read][j-1]) & firstBitUint64_t;
			assert(std::get<1>(fixedChunksPerRead[read][j-1]) >= std::get<0>(fixedChunksPerRead[read][j]));
			std::pair<size_t, size_t> firstAndLast = firstAndLastMinimizerWithinSpan(minimizerPositionsPerRead[read], std::get<0>(fixedChunksPerRead[read][j]), std::get<1>(fixedChunksPerRead[read][j-1]), kmerSize);
			assert(firstAndLast.second >= firstAndLast.first);
			size_t overlap = firstAndLast.second - firstAndLast.first;
			assert(prevchunk < countSmallchunksPerChunk.size());
			assert(currchunk < countSmallchunksPerChunk.size());
			assert(overlap <= countSmallchunksPerChunk[prevchunk]);
			assert(overlap <= countSmallchunksPerChunk[currchunk]);
			auto key = canonNodePair(std::get<2>(fixedChunksPerRead[read][j-1]), std::get<2>(fixedChunksPerRead[read][j]));
			assert(smallChunkOverlap.count(key) == 0 || smallChunkOverlap.at(key) == overlap);
			smallChunkOverlap[key] = overlap;
		}
	}
	return std::make_pair(countSmallchunksPerChunk, smallChunkOverlap);
}

std::tuple<std::vector<std::vector<std::pair<size_t, bool>>>, std::vector<std::vector<std::tuple<size_t, size_t, bool>>>, phmap::flat_hash_set<size_t>> getSmallchunkInfo(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const std::vector<size_t>& countSmallchunksPerChunk, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& smallChunkOverlap, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& fixedChunksPerRead, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const size_t kmerSize)
{
	assert(fixedChunksPerRead.size()*2 == minimizerPositionsPerRead.size());
	std::vector<std::pair<size_t, bool>> smallChunkParent;
	std::vector<std::vector<std::pair<size_t, bool>>> smallChunksPerUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		size_t countSmallChunksHere = 0;
		for (size_t j = 0; j < graph.chunksInUnitig[i].size(); j++)
		{
			assert((graph.chunksInUnitig[i][j] & maskUint64_t) < countSmallchunksPerChunk.size());
			countSmallChunksHere += countSmallchunksPerChunk[graph.chunksInUnitig[i][j] & maskUint64_t];
			if (j >= 1)
			{
				auto key = canonNodePair(graph.chunksInUnitig[i][j-1], graph.chunksInUnitig[i][j]);
				assert(smallChunkOverlap.count(key) == 1);
				assert(smallChunkOverlap.at(key) <= countSmallChunksHere);
				countSmallChunksHere -= smallChunkOverlap.at(key);
			}
		}
		smallChunksPerUnitig.emplace_back();
		for (size_t j = 0; j < countSmallChunksHere; j++)
		{
			smallChunksPerUnitig.back().emplace_back(smallChunkParent.size(), true);
			smallChunkParent.emplace_back(smallChunksPerUnitig.back().back().first, true);
		}
	}
	assert(smallChunkParent.back().first+1 == smallChunkParent.size());
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : graph.edges.getEdges(fw))
		{
			uint64_t fromChunk = graph.chunksInUnitig[i].back();
			uint64_t toChunk;
			if (edge.second)
			{
				toChunk = graph.chunksInUnitig[edge.first][0];
			}
			else
			{
				toChunk = graph.chunksInUnitig[edge.first].back() ^ firstBitUint64_t;
			}
			auto key = canonNodePair(fromChunk, toChunk);
			assert(smallChunkOverlap.count(key) == 1);
			size_t overlap = smallChunkOverlap.at(key);
			if (overlap >= 1)
			{
				assert(i < smallChunksPerUnitig.size());
				assert(edge.first < smallChunksPerUnitig.size());
				assert(overlap <= smallChunksPerUnitig[i].size());
				assert(overlap <= smallChunksPerUnitig[edge.first].size());
				if (edge.second)
				{
					for (size_t k = 0; k < overlap; k++)
					{
						mergeAllowPalindrome(smallChunkParent, smallChunksPerUnitig[i][smallChunksPerUnitig[i].size()-overlap+k], smallChunksPerUnitig[edge.first][k], true);
					}
				}
				else
				{
					for (size_t k = 0; k < overlap; k++)
					{
						mergeAllowPalindrome(smallChunkParent, smallChunksPerUnitig[i][smallChunksPerUnitig[i].size()-overlap+k], smallChunksPerUnitig[edge.first][smallChunksPerUnitig[edge.first].size()-1-k], false);
					}
				}
			}
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : graph.edges.getEdges(bw))
		{
			uint64_t fromChunk = graph.chunksInUnitig[i][0] ^ firstBitUint64_t;
			uint64_t toChunk;
			if (edge.second)
			{
				toChunk = graph.chunksInUnitig[edge.first][0];
			}
			else
			{
				toChunk = graph.chunksInUnitig[edge.first].back() ^ firstBitUint64_t;
			}
			auto key = canonNodePair(fromChunk, toChunk);
			assert(smallChunkOverlap.count(key) == 1);
			size_t overlap = smallChunkOverlap.at(key);
			if (overlap >= 1)
			{
				assert(i < smallChunksPerUnitig.size());
				assert(edge.first < smallChunksPerUnitig.size());
				assert(overlap <= smallChunksPerUnitig[i].size());
				assert(overlap <= smallChunksPerUnitig[edge.first].size());
				if (edge.second)
				{
					for (size_t k = 0; k < overlap; k++)
					{
						mergeAllowPalindrome(smallChunkParent, smallChunksPerUnitig[i][overlap-1-k], smallChunksPerUnitig[edge.first][k], false);
					}
				}
				else
				{
					for (size_t k = 0; k < overlap; k++)
					{
						mergeAllowPalindrome(smallChunkParent, smallChunksPerUnitig[i][overlap-1-k], smallChunksPerUnitig[edge.first][smallChunksPerUnitig[edge.first].size()-1-k], true);
					}
				}
			}
		}
	}
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (size_t j = 0; j < smallChunksPerUnitig[i].size(); j++)
		{
			smallChunksPerUnitig[i][j] = find(smallChunkParent, smallChunksPerUnitig[i][j]);
		}
	}
	phmap::flat_hash_set<size_t> palindromicSmallChunks;
	for (size_t i = 0; i < smallChunkParent.size(); i++)
	{
		if (find(smallChunkParent, i).first != i) continue;
		if (find(smallChunkParent, i).second == false)
		{
			palindromicSmallChunks.emplace(i);
		}
	}
	std::vector<std::tuple<size_t, size_t, size_t, bool>> chunkSmallChunkPositionInUnitig; // unitig, startpos, endpos, fw. startpos inclusive endpos exclusive even when not fw.
	chunkSmallChunkPositionInUnitig.resize(countSmallchunksPerChunk.size(), std::make_tuple(std::numeric_limits<size_t>::max(), 0, 0, true));
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		size_t startpos = 0;
		for (size_t j = 0; j < graph.chunksInUnitig[i].size(); j++)
		{
			if (j >= 1)
			{
				uint64_t prev = graph.chunksInUnitig[i][j-1];
				uint64_t curr = graph.chunksInUnitig[i][j];
				auto key = canonNodePair(prev, curr);
				assert(smallChunkOverlap.count(key) == 1);
				size_t overlap = smallChunkOverlap.at(key);
				assert(startpos >= overlap);
				startpos -= overlap;
			}
			size_t chunk = graph.chunksInUnitig[i][j] & maskUint64_t;
			bool fw = graph.chunksInUnitig[i][j] & firstBitUint64_t;
			assert(chunk < countSmallchunksPerChunk.size());
			size_t endpos = startpos + countSmallchunksPerChunk[chunk];
			assert(std::get<0>(chunkSmallChunkPositionInUnitig[chunk]) == std::numeric_limits<size_t>::max());
			chunkSmallChunkPositionInUnitig[chunk] = std::make_tuple(i, startpos, endpos, fw);
			startpos = endpos;
		}
		assert(startpos == smallChunksPerUnitig[i].size());
	}
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> smallchunksPerRead; // readminimizerindex, smallchunkindex, smallchunkfw
	smallchunksPerRead.resize(fixedChunksPerRead.size());
	for (size_t i = 0; i < fixedChunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < fixedChunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(fixedChunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(fixedChunksPerRead[i][j]) & maskUint64_t;
			bool fw = std::get<2>(fixedChunksPerRead[i][j]) & firstBitUint64_t;
			if (std::get<0>(chunkSmallChunkPositionInUnitig[chunk]) == std::numeric_limits<size_t>::max()) continue;
			std::pair<size_t, size_t> firstAndLast = firstAndLastMinimizerWithinSpan(minimizerPositionsPerRead[i], std::get<0>(fixedChunksPerRead[i][j]), std::get<1>(fixedChunksPerRead[i][j]), kmerSize);
			assert(firstAndLast.second > firstAndLast.first);
			size_t unitig = std::get<0>(chunkSmallChunkPositionInUnitig[chunk]);
			size_t start = std::get<1>(chunkSmallChunkPositionInUnitig[chunk]);
			size_t end = std::get<2>(chunkSmallChunkPositionInUnitig[chunk]);
			assert(end-start == countSmallchunksPerChunk[chunk]);
			assert(end-start == firstAndLast.second-firstAndLast.first);
			size_t chunkfw = std::get<3>(chunkSmallChunkPositionInUnitig[chunk]);
			bool readAndUnitigOrientationMatch = (fw == chunkfw);
			if (readAndUnitigOrientationMatch)
			{
				for (size_t k = 0; k < end-start; k++)
				{
					std::pair<size_t, bool> key = smallChunksPerUnitig[unitig][start+k];
					smallchunksPerRead[i].emplace_back(firstAndLast.first+k, key.first, key.second);
				}
			}
			else
			{
				for (size_t k = 0; k < end-start; k++)
				{
					std::pair<size_t, bool> key = reverse(smallChunksPerUnitig[unitig][end-1-k]);
					smallchunksPerRead[i].emplace_back(firstAndLast.first+k, key.first, key.second);
				}
			}
		}
	}
	return std::make_tuple(smallChunksPerUnitig, smallchunksPerRead, palindromicSmallChunks);
}

phmap::flat_hash_map<std::pair<size_t, size_t>, std::vector<std::tuple<size_t, size_t, size_t, bool>>> getReadPiecesPerSmallchunk(const std::vector<std::vector<std::pair<std::pair<size_t, size_t>, bool>>>& smallchunksPerRead, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const size_t kmerSize)
{
	phmap::flat_hash_map<std::pair<size_t, size_t>, std::vector<std::tuple<size_t, size_t, size_t, bool>>> result;
	for (size_t i = 0; i < smallchunksPerRead.size(); i++)
	{
		if (minimizerPositionsPerRead[i].size() < 2) continue;
		assert(minimizerPositionsPerRead[i].size() == smallchunksPerRead[i].size()+1);
		for (size_t j = 0; j < smallchunksPerRead[i].size(); j++)
		{
			if (smallchunksPerRead[i][j].first.first == std::numeric_limits<size_t>::max()) continue;
			std::pair<size_t, size_t> key = smallchunksPerRead[i][j].first;
			result[key].emplace_back(i, minimizerPositionsPerRead[i][j], minimizerPositionsPerRead[i][j+1]+kmerSize-1, smallchunksPerRead[i][j].second);
		}
	}
	return result;
}

std::vector<TwobitString> getSmallchunkSequences(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<std::tuple<size_t, size_t, bool>>>& smallchunksPerRead, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const phmap::flat_hash_set<size_t>& palindromicSmallChunks, const size_t numThreads, const size_t kmerSize)
{
	assert(smallchunksPerRead.size() == sequenceIndex.size());
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> readPiecesPerSmallChunk; // read, readminimizerindex, fw
	for (size_t i = 0; i < smallchunksPerRead.size(); i++)
	{
		for (auto t : smallchunksPerRead[i])
		{
			while (readPiecesPerSmallChunk.size() <= std::get<1>(t)) readPiecesPerSmallChunk.emplace_back();
			readPiecesPerSmallChunk[std::get<1>(t)].emplace_back(i, std::get<0>(t), std::get<2>(t));
		}
	}
	std::mutex indexMutex;
	std::mutex resultMutex;
	size_t nextnum = 0;
	std::vector<TwobitString> result;
	result.resize(readPiecesPerSmallChunk.size());
	std::vector<std::thread> threads;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&indexMutex, &sequenceIndex, &resultMutex, &nextnum, &result, &minimizerPositionsPerRead, &readPiecesPerSmallChunk, kmerSize]()
		{
			while (true)
			{
				size_t indexHere = 0;
				{
					std::lock_guard<std::mutex> lock { indexMutex };
					indexHere = nextnum;
					nextnum += 1;
				}
				if (indexHere >= readPiecesPerSmallChunk.size()) return;
				if (readPiecesPerSmallChunk[indexHere].size() == 0) continue;
				phmap::flat_hash_set<std::tuple<size_t, size_t, bool>> uniquePieces { readPiecesPerSmallChunk[indexHere].begin(), readPiecesPerSmallChunk[indexHere].end() };
				phmap::flat_hash_map<std::string, size_t> sequenceCounts;
				for (auto t : uniquePieces)
				{
					assert(std::get<0>(t) < minimizerPositionsPerRead.size());
					assert(std::get<1>(t)+1 < minimizerPositionsPerRead[std::get<0>(t)].size());
					auto sequence = sequenceIndex.getSubstring(std::get<0>(t), minimizerPositionsPerRead[std::get<0>(t)][std::get<1>(t)], minimizerPositionsPerRead[std::get<0>(t)][std::get<1>(t)+1] - minimizerPositionsPerRead[std::get<0>(t)][std::get<1>(t)] + kmerSize - 1 + 1);
					assert(sequence.size() >= kmerSize+1);
					if (!std::get<2>(t))
					{
						sequence = revCompRaw(sequence);
					}
					sequenceCounts[sequence] += 1;
				}
				auto first = sequenceCounts.begin()->first;
				for (const auto& pair : sequenceCounts)
				{
					assert(pair.first.size() >= kmerSize);
					assert(pair.first.substr(0, kmerSize) == first.substr(0, kmerSize));
					assert(pair.first.substr(pair.first.size()-kmerSize) == first.substr(first.size()-kmerSize));
				}
				std::string consensus;
				if (sequenceCounts.size() >= 2 && uniquePieces.size() >= 3)
				{
					consensus = getConsensus(sequenceCounts, uniquePieces.size());
				}
				else
				{
					consensus = first;
				}
				assert(consensus.size() >= kmerSize+1);
				assert(consensus.substr(0, kmerSize) == first.substr(0, kmerSize));
				assert(consensus.substr(consensus.size()-kmerSize) == first.substr(first.size()-kmerSize));
				{
					std::lock_guard<std::mutex> lock { resultMutex };
					result[indexHere] = consensus;
				}
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	for (auto i : palindromicSmallChunks)
	{
		std::string seq = result[i].toString();
		std::string revcomp = revCompRaw(seq);
		for (size_t j = seq.size()/2; j < seq.size(); j++)
		{
			if (seq[j] != revcomp[j]) seq[j] = revcomp[j];
		}
		result[i] = seq;
	}
	return result;
}

phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getSmallchunkOverlapPerchunk(const std::vector<std::vector<std::pair<std::pair<size_t, size_t>, bool>>>& smallchunksPerRead, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& fixedChunksPerRead, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const size_t kmerSize)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> result;
	for (size_t i = 0; i < fixedChunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < fixedChunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(fixedChunksPerRead[i][j-1]))) continue;
			if (NonexistantChunk(std::get<2>(fixedChunksPerRead[i][j]))) continue;
			std::pair<size_t, size_t> prev = firstAndLastMinimizerWithinSpan(minimizerPositionsPerRead[i], std::get<0>(fixedChunksPerRead[i][j-1]), std::get<1>(fixedChunksPerRead[i][j-1]), kmerSize);
			std::pair<size_t, size_t> curr = firstAndLastMinimizerWithinSpan(minimizerPositionsPerRead[i], std::get<0>(fixedChunksPerRead[i][j]), std::get<1>(fixedChunksPerRead[i][j]), kmerSize);
			assert(prev.second >= curr.first);
			size_t overlap = prev.second - curr.first;
			auto key = canonNodePair(std::get<2>(fixedChunksPerRead[i][j-1]), std::get<2>(fixedChunksPerRead[i][j]));
			assert(result.count(key) == 0 || result.at(key) == overlap);
			result[key] = overlap;
		}
	}
	return result;
}

std::tuple<std::vector<TwobitString>, std::vector<std::vector<std::pair<size_t, size_t>>>, std::vector<std::vector<size_t>>> getUnitigSequences(const ChunkUnitigGraph& graph, const std::vector<size_t>& countSmallchunksPerChunk, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& smallChunkOverlap, const std::vector<std::vector<std::pair<size_t, bool>>>& smallChunksPerUnitig, const std::vector<TwobitString>& smallchunkSequences, const size_t kmerSize)
{
	std::vector<TwobitString> result;
	std::vector<std::vector<size_t>> minimizerPositionsInUnitigs;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunkPositions;
	minimizerPositionsInUnitigs.resize(graph.chunksInUnitig.size());
	result.resize(graph.chunksInUnitig.size());
	chunkPositions.resize(graph.chunksInUnitig.size());
	for (size_t unitig = 0; unitig < graph.chunksInUnitig.size(); unitig++)
	{
		std::string resultHere;
		minimizerPositionsInUnitigs[unitig].emplace_back(0);
		for (size_t j = 0; j < smallChunksPerUnitig[unitig].size(); j++)
		{
			std::string add = smallchunkSequences[smallChunksPerUnitig[unitig][j].first].toString();
			assert(add.size() >= kmerSize+1);
			if (!smallChunksPerUnitig[unitig][j].second) add = revCompRaw(add);
			if (j > 0)
			{
				assert(resultHere.size() >= kmerSize+1);
				assert(add.substr(0, kmerSize) == resultHere.substr(resultHere.size()-kmerSize));
				resultHere.erase(resultHere.begin()+resultHere.size()-kmerSize, resultHere.end());
			}
			resultHere += add;
			minimizerPositionsInUnitigs[unitig].emplace_back(resultHere.size() - kmerSize);
		}
		assert(resultHere.size() >= kmerSize+1);
		result[unitig] = resultHere;
		assert(result[unitig].size() >= kmerSize+1);
		size_t startMinimizerIndex = 0;
		for (size_t j = 0; j < graph.chunksInUnitig[unitig].size(); j++)
		{
			if (j >= 1)
			{
				uint64_t prev = graph.chunksInUnitig[unitig][j-1];
				uint64_t curr = graph.chunksInUnitig[unitig][j];
				auto key = canonNodePair(prev, curr);
				assert(smallChunkOverlap.count(key) == 1);
				size_t overlap = smallChunkOverlap.at(key);
				assert(startMinimizerIndex >= overlap);
				startMinimizerIndex -= overlap;
			}
			size_t chunk = graph.chunksInUnitig[unitig][j] & maskUint64_t;
			size_t endMinimizerIndex = startMinimizerIndex+countSmallchunksPerChunk[chunk];
			chunkPositions[unitig].emplace_back(minimizerPositionsInUnitigs[unitig][startMinimizerIndex], minimizerPositionsInUnitigs[unitig][endMinimizerIndex]+kmerSize-1);
			startMinimizerIndex = endMinimizerIndex;
		}
		assert(startMinimizerIndex+1 == minimizerPositionsInUnitigs[unitig].size());
	}
	return std::make_tuple(std::move(result), std::move(chunkPositions), std::move(minimizerPositionsInUnitigs));
}

phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> getFixedOverlaps(const ChunkUnitigGraph& graph, const std::vector<std::vector<size_t>>& minimizerPositionsInUnitigs, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& smallChunkOverlap, const size_t kmerSize)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> result;
	for (const auto& unitigpair : graph.edgeOverlaps)
	{
		uint64_t chunkFrom;
		uint64_t chunkTo;
		if (unitigpair.first.first & firstBitUint64_t)
		{
			chunkFrom = graph.chunksInUnitig[unitigpair.first.first & maskUint64_t].back();
		}
		else
		{
			chunkFrom = graph.chunksInUnitig[unitigpair.first.first & maskUint64_t][0] ^ firstBitUint64_t;
		}
		if (unitigpair.first.second & firstBitUint64_t)
		{
			chunkTo = graph.chunksInUnitig[unitigpair.first.second & maskUint64_t][0];
		}
		else
		{
			chunkTo = graph.chunksInUnitig[unitigpair.first.second & maskUint64_t].back() ^ firstBitUint64_t;
		}
		auto chunkpair = canonNodePair(chunkFrom, chunkTo);
		assert(smallChunkOverlap.count(chunkpair) == 1);
		size_t overlap = smallChunkOverlap.at(chunkpair);
		size_t overlapBases = 0;
		assert((unitigpair.first.first & maskUint64_t) < minimizerPositionsInUnitigs.size());
		assert((unitigpair.first.second & maskUint64_t) < minimizerPositionsInUnitigs.size());
		assert(overlap < minimizerPositionsInUnitigs[unitigpair.first.second & maskUint64_t].size());
		if (unitigpair.first.second & firstBitUint64_t)
		{
			overlapBases = minimizerPositionsInUnitigs[unitigpair.first.second & maskUint64_t][overlap] + kmerSize;
		}
		else
		{
			overlapBases = minimizerPositionsInUnitigs[unitigpair.first.second & maskUint64_t][minimizerPositionsInUnitigs[unitigpair.first.second & maskUint64_t].size()-1] - minimizerPositionsInUnitigs[unitigpair.first.second & maskUint64_t][minimizerPositionsInUnitigs[unitigpair.first.second & maskUint64_t].size()-1-overlap] + kmerSize;
		}
		if (unitigpair.first.first & firstBitUint64_t)
		{
			assert(overlapBases == (minimizerPositionsInUnitigs[unitigpair.first.first & maskUint64_t][minimizerPositionsInUnitigs[unitigpair.first.first & maskUint64_t].size()-1] - minimizerPositionsInUnitigs[unitigpair.first.first & maskUint64_t][minimizerPositionsInUnitigs[unitigpair.first.first & maskUint64_t].size()-1-overlap] + kmerSize));
		}
		else
		{
			assert(overlapBases == (minimizerPositionsInUnitigs[unitigpair.first.first & maskUint64_t][overlap] + kmerSize));
		}
		result[unitigpair.first] = overlapBases;
	}
	return result;
}

std::tuple<std::vector<TwobitString>, std::vector<std::vector<std::pair<size_t, size_t>>>, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> getUnitigDBGSequences(const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const size_t kmerSize, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& fixedChunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const size_t numThreads)
{
	std::vector<size_t> countSmallchunksPerChunk;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> smallChunkOverlap;
	std::tie(countSmallchunksPerChunk, smallChunkOverlap) = getSmallChunkSizesAndOverlaps(fixedChunksPerRead, minimizerPositionsPerRead, kmerSize);
	std::vector<std::vector<std::pair<size_t, bool>>> smallChunksPerUnitig;
	std::vector<std::vector<std::tuple<size_t, size_t, bool>>> smallchunksPerRead; // readminimizerindex, smallchunkindex, smallchunkfw
	phmap::flat_hash_set<size_t> palindromicSmallChunks;
	std::tie(smallChunksPerUnitig, smallchunksPerRead, palindromicSmallChunks) = getSmallchunkInfo(graph, readPaths, countSmallchunksPerChunk, smallChunkOverlap, fixedChunksPerRead, minimizerPositionsPerRead, kmerSize);
	std::vector<TwobitString> smallchunkSequences = getSmallchunkSequences(sequenceIndex, smallchunksPerRead, minimizerPositionsPerRead, palindromicSmallChunks, numThreads, kmerSize);
	for (size_t i = 0; i < smallChunksPerUnitig.size(); i++)
	{
		for (size_t j = 0; j < smallChunksPerUnitig[i].size(); j++)
		{
			assert(smallChunksPerUnitig[i][j].first < smallchunkSequences.size());
			assert(smallchunkSequences[smallChunksPerUnitig[i][j].first].size() >= kmerSize+1);
		}
	}
	std::vector<TwobitString> unitigSequences;
	std::vector<std::vector<size_t>> minimizerPositionsInUnitigs;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunkSequencePositionsWithinUnitigs; // end inclusive!
	std::tie(unitigSequences, chunkSequencePositionsWithinUnitigs, minimizerPositionsInUnitigs) = getUnitigSequences(graph, countSmallchunksPerChunk, smallChunkOverlap, smallChunksPerUnitig, smallchunkSequences, kmerSize);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> fixedOverlaps = getFixedOverlaps(graph, minimizerPositionsInUnitigs, smallChunkOverlap, kmerSize);
	return std::make_tuple(std::move(unitigSequences), std::move(chunkSequencePositionsWithinUnitigs), std::move(fixedOverlaps));
}
