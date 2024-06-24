#include <iostream>
#include "ChunkUnitigGraph.h"
#include "Common.h"

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
			assert(std::get<0>(chunksPerRead[i][j]) > std::get<0>(chunksPerRead[i][j-1]));
			assert(std::get<1>(chunksPerRead[i][j]) > std::get<1>(chunksPerRead[i][j-1]));
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
