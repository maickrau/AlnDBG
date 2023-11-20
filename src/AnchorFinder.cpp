#include <iostream>
#include <phmap.h>
#include "AnchorFinder.h"
#include "Common.h"
#include "MBGCommon.h"

std::tuple<uint64_t, uint64_t, uint64_t> canon(std::pair<size_t, bool> first, std::pair<size_t, bool> second, std::pair<size_t, bool> third)
{
	uint64_t one = first.first + (first.second ? firstBitUint64_t : 0);
	uint64_t two = second.first + (second.second ? firstBitUint64_t : 0);
	uint64_t three = third.first + (third.second ? firstBitUint64_t : 0);
	std::tuple<uint64_t, uint64_t, uint64_t> fw { one, two, three };
	std::tuple<uint64_t, uint64_t, uint64_t> bw { three ^ firstBitUint64_t, two ^ firstBitUint64_t, one ^ firstBitUint64_t };
	return std::min(fw, bw);
}

std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> getPossiblyTransitiveTriplets(const SparseEdgeContainer& rawEdgesWithTransitiveEdges)
{
	phmap::flat_hash_set<std::tuple<uint64_t, uint64_t, uint64_t>> result;
	for (size_t i = 0; i < rawEdgesWithTransitiveEdges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : rawEdgesWithTransitiveEdges.getEdges(fw))
		{
			for (auto edge2 : rawEdgesWithTransitiveEdges.getEdges(edge))
			{
				if (!rawEdgesWithTransitiveEdges.hasEdge(fw, edge2)) continue;
				result.emplace(canon(fw, edge, edge2));
			}
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : rawEdgesWithTransitiveEdges.getEdges(bw))
		{
			for (auto edge2 : rawEdgesWithTransitiveEdges.getEdges(edge))
			{
				if (!rawEdgesWithTransitiveEdges.hasEdge(bw, edge2)) continue;
				result.emplace(canon(bw, edge, edge2));
			}
		}
	}
	std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> resultvec { result.begin(), result.end() };
	return resultvec;
}

SparseEdgeContainer getEdgesWithoutTransitiveEdges(const SparseEdgeContainer& rawEdgesWithTransitiveEdges, const std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>& maybeTransitiveEdges, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<bool>& anchor, const double approxOneHapCoverage)
{
	SparseEdgeContainer result { rawEdgesWithTransitiveEdges };
	// todo implement
	return result;
}

SparseEdgeContainer getAnchorEdges(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<bool>& anchor, const double approxOneHapCoverage)
{
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverage;
	edgeCoverage.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		uint64_t lastAnchor = std::numeric_limits<size_t>::max();
		size_t lastRightClip = 0;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				if (!anchor[readPaths[i].paths[j].path[k] & maskUint64_t]) continue;
				if (lastAnchor == std::numeric_limits<size_t>::max())
				{
					lastAnchor = readPaths[i].paths[j].path[k];
					if (k == readPaths[i].paths[j].path.size()-1)
					{
						lastRightClip = readPaths[i].paths[j].pathRightClipKmers;
					}
					else
					{
						lastRightClip = 0;
					}
					continue;
				}
				if (k == 0 && lastRightClip > 0 && lastAnchor == readPaths[i].paths[j].path[k])
				{
					if (unitigGraph.lengths[lastAnchor & maskUint64_t] - lastRightClip <= readPaths[i].paths[j].pathLeftClipKmers)
					{
						lastAnchor = readPaths[i].paths[j].path[k];
						if (k == readPaths[i].paths[j].path.size()-1)
						{
							lastRightClip = readPaths[i].paths[j].pathRightClipKmers;
						}
						else
						{
							lastRightClip = 0;
						}
						continue;
					}
				}
				std::pair<size_t, bool> fromEdge { lastAnchor & maskUint64_t, lastAnchor & firstBitUint64_t };
				std::pair<size_t, bool> toEdge { readPaths[i].paths[j].path[k] & maskUint64_t, readPaths[i].paths[j].path[k] & firstBitUint64_t };
				auto key = canon(fromEdge, toEdge);
				size_t coverage = 0;
				assert(key.first.first < edgeCoverage.size());
				assert(key.second.first < edgeCoverage.size());
				if (edgeCoverage.hasValue(key.first, key.second)) coverage = edgeCoverage.get(key.first, key.second);
				edgeCoverage.set(key.first, key.second, coverage+1);
				lastAnchor = readPaths[i].paths[j].path[k];
				if (k == readPaths[i].paths[j].path.size()-1)
				{
					lastRightClip = readPaths[i].paths[j].pathRightClipKmers;
				}
				else
				{
					lastRightClip = 0;
				}
				continue;
			}
		}
	}
	SparseEdgeContainer rawEdgesWithTransitiveEdges;
	rawEdgesWithTransitiveEdges.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (!anchor[i]) continue;
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : edgeCoverage.getValues(fw))
		{
			if ((double)pair.second < approxOneHapCoverage * 0.5) continue;
			assert(anchor[pair.first.first]);
			assert(fw.first < rawEdgesWithTransitiveEdges.size());
			assert(pair.first.first < rawEdgesWithTransitiveEdges.size());
			rawEdgesWithTransitiveEdges.addEdge(fw, pair.first);
			rawEdgesWithTransitiveEdges.addEdge(reverse(pair.first), reverse(fw));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : edgeCoverage.getValues(bw))
		{
			if ((double)pair.second < approxOneHapCoverage * 0.5) continue;
			assert(anchor[pair.first.first]);
			assert(bw.first < rawEdgesWithTransitiveEdges.size());
			assert(pair.first.first < rawEdgesWithTransitiveEdges.size());
			rawEdgesWithTransitiveEdges.addEdge(bw, pair.first);
			rawEdgesWithTransitiveEdges.addEdge(reverse(pair.first), reverse(bw));
		}
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : rawEdgesWithTransitiveEdges.getEdges(fw))
		{
			assert(rawEdgesWithTransitiveEdges.hasEdge(reverse(pair), reverse(fw)));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : rawEdgesWithTransitiveEdges.getEdges(bw))
		{
			assert(rawEdgesWithTransitiveEdges.hasEdge(reverse(pair), reverse(bw)));
		}
	}
	auto maybeTransitiveEdges = getPossiblyTransitiveTriplets(rawEdgesWithTransitiveEdges);
	std::cerr << maybeTransitiveEdges.size() << " maybe transitive edges" << std::endl;
	for (auto t : maybeTransitiveEdges)
	{
		std::cerr << "transitive triplet " << ((std::get<0>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<0>(t) & maskUint64_t) << " "  << ((std::get<1>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<1>(t) & maskUint64_t) << " "  << ((std::get<2>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(t) & maskUint64_t) << std::endl;
	}
	// return rawEdgesWithTransitiveEdges;
	return getEdgesWithoutTransitiveEdges(rawEdgesWithTransitiveEdges, maybeTransitiveEdges, unitigGraph, readPaths, anchor, approxOneHapCoverage);
	// SparseEdgeContainer edges = getEdgesWithoutTransitiveEdges(rawEdgesWithTransitiveEdges, maybeTransitiveEdges, unitigGraph, readPaths, anchor, approxOneHapCoverage);
	// return edges;
}

void extendAnchors(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, std::vector<bool>& anchor)
{
	// todo implement
	return;
}

size_t getPloidy(const std::vector<uint64_t>& nodes, const UnitigGraph& unitigGraph, const double approxOneHapCoverage)
{
	assert(nodes.size() >= 1);
	double coverageSum = 0;
	size_t coverageDivisor = 0;
	for (auto node : nodes)
	{
		assert((node & maskUint64_t) < unitigGraph.nodeCount());
		assert(unitigGraph.lengths[node & maskUint64_t] >= 1);
		coverageDivisor += unitigGraph.lengths[node & maskUint64_t];
		coverageSum += unitigGraph.lengths[node & maskUint64_t] * unitigGraph.coverages[node & maskUint64_t];
	}
	assert(coverageDivisor > 0);
	return (coverageSum / coverageDivisor) / approxOneHapCoverage + 0.5;
}

std::vector<uint64_t> getChainOneWay(const uint64_t start, const SparseEdgeContainer& edges)
{
	uint64_t pos = start;
	std::vector<uint64_t> result;
	result.push_back(pos);
	while (edges.getEdges(std::make_pair(pos & maskUint64_t, pos & firstBitUint64_t)).size() == 1)
	{
		std::pair<size_t, bool> nextp = edges.getEdges(std::make_pair(pos & maskUint64_t, pos & firstBitUint64_t))[0];
		if (edges.getEdges(reverse(nextp)).size() != 1) break;
		uint64_t next = nextp.first + (nextp.second ? firstBitUint64_t : 0);
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next == start) break;
		pos = next;
		result.push_back(pos);
	}
	return result;
}

std::vector<uint64_t> getChain(std::vector<bool>& checked, const size_t startNode, const SparseEdgeContainer& edges)
{
	assert(edges.size() == checked.size());
	assert(startNode < checked.size());
	for (size_t i = 0; i < edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : edges.getEdges(fw))
		{
			assert(edges.hasEdge(reverse(pair), reverse(fw)));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : edges.getEdges(bw))
		{
			assert(edges.hasEdge(reverse(pair), reverse(bw)));
		}
	}
	assert(!checked[startNode]);
	std::vector<uint64_t> fwChain = getChainOneWay(startNode + firstBitUint64_t, edges);
	std::vector<uint64_t> bwChain = getChainOneWay(startNode, edges);
	assert(fwChain.size() >= 1);
	assert(bwChain.size() >= 1);
	assert(fwChain[0] == startNode + firstBitUint64_t);
	assert(bwChain[0] == startNode);
	if (fwChain.size() == bwChain.size())
	{
		if (fwChain.size() >= 2)
		{
			if (fwChain.back() == (bwChain[1] ^ firstBitUint64_t))
			{
				for (size_t i = 0; i < fwChain.size(); i++)
				{
					assert(i == 0 || (fwChain[i] == (bwChain[bwChain.size()-i] ^ firstBitUint64_t)));
					assert(!checked[fwChain[i] & maskUint64_t]);
					checked[fwChain[i] & maskUint64_t] = true;
				}
				return fwChain;
			}
		}
	}
	std::reverse(bwChain.begin(), bwChain.end());
	for (size_t i = 0; i < bwChain.size(); i++)
	{
		bwChain[i] ^= firstBitUint64_t;
	}
	assert(bwChain.back() == fwChain[0]);
	bwChain.insert(bwChain.end(), fwChain.begin()+1, fwChain.end());
	for (size_t i = 0; i < bwChain.size(); i++)
	{
		assert((bwChain[i] & maskUint64_t) < checked.size());
		assert(!checked[bwChain[i] & maskUint64_t]);
		checked[bwChain[i] & maskUint64_t] = true;
	}
	return bwChain;
}

void getDistances(phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<size_t>>& distances, const std::vector<bool> anchor, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths)
{
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		uint64_t lastAnchor = std::numeric_limits<uint64_t>::max();
		size_t lastAnchorOffset = 0;
		size_t lastAnchorPlus = 0;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readOffset = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				if (!anchor[node & maskUint64_t])
				{
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				if (lastAnchor == std::numeric_limits<uint64_t>::max())
				{
					lastAnchor = node;
					lastAnchorOffset = readOffset;
					lastAnchorPlus = 0;
					if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				if (distances.count(std::make_pair(lastAnchor, node)) == 1)
				{
					distances.at(std::make_pair(lastAnchor, node)).emplace_back(readOffset - lastAnchorOffset + lastAnchorPlus);
					if (k == 0) distances.at(std::make_pair(lastAnchor, node)).back() -= readPaths[i].paths[j].pathLeftClipKmers;
					lastAnchor = node;
					lastAnchorOffset = readOffset;
					lastAnchorPlus = 0;
					if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				if (distances.count(std::make_pair(node ^ firstBitUint64_t, lastAnchor ^ firstBitUint64_t)) == 1)
				{
					distances.at(std::make_pair(node ^ firstBitUint64_t, lastAnchor ^ firstBitUint64_t)).emplace_back(readOffset - lastAnchorOffset + lastAnchorPlus);
					if (k == 0) distances.at(std::make_pair(node ^ firstBitUint64_t, lastAnchor ^ firstBitUint64_t)).back() -= readPaths[i].paths[j].pathLeftClipKmers;
					distances.at(std::make_pair(node ^ firstBitUint64_t, lastAnchor ^ firstBitUint64_t)).back() += unitigGraph.lengths[node & maskUint64_t] - unitigGraph.lengths[lastAnchor & maskUint64_t];
					lastAnchor = node;
					lastAnchorOffset = readOffset;
					lastAnchorPlus = 0;
					if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				lastAnchor = node;
				readOffset += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
			}
		}
	}
}

double average(const std::vector<size_t>& values)
{
	assert(values.size() >= 1);
	double sum = 0;
	for (auto value : values)
	{
		sum += value;
	}
	return sum / values.size();
}

void addOffsets(std::vector<AnchorChain>& chains, const std::vector<bool> anchor, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<size_t>> distances;
	for (size_t i = 0; i < chains.size(); i++)
	{
		for (size_t j = 1; j < chains[i].nodes.size(); j++)
		{
			assert(distances.count(std::make_pair(chains[i].nodes[j-1], chains[i].nodes[j])) == 0);
			distances[std::make_pair(chains[i].nodes[j-1], chains[i].nodes[j])];
		}
	}
	getDistances(distances, anchor, unitigGraph, readPaths);
	for (size_t i = 0; i < chains.size(); i++)
	{
		chains[i].nodeOffsets.push_back(0);
		for (size_t j = 1; j < chains[i].nodes.size(); j++)
		{
			auto key = std::make_pair(chains[i].nodes[j-1], chains[i].nodes[j]);
			assert(distances.at(key).size() >= 1);
			chains[i].nodeOffsets.push_back(chains[i].nodeOffsets.back() + average(distances.at(key)));
		}
		assert(chains[i].nodeOffsets.size() == chains[i].nodes.size());
	}
}

std::vector<AnchorChain> getAnchorChains(const std::vector<bool>& anchor, const SparseEdgeContainer& edges, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	std::vector<AnchorChain> result;
	std::vector<bool> checked;
	checked.resize(anchor.size(), false);
	for (size_t i = 0; i < anchor.size(); i++)
	{
		if (!anchor[i]) continue;
		if (checked[i]) continue;
		result.emplace_back();
		result.back().nodes = getChain(checked, i, edges);
		result.back().ploidy = getPloidy(result.back().nodes, unitigGraph, approxOneHapCoverage);
	}
	addOffsets(result, anchor, unitigGraph, readPaths);
	return result;
}

std::vector<AnchorChain> getAnchorChains(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.lengths[i] >= 100) anchor[i] = true;
	}
	// extendAnchors(unitigGraph, readPaths, anchor);
	SparseEdgeContainer edges = getAnchorEdges(unitigGraph, readPaths, anchor, approxOneHapCoverage);
	for (size_t i = 0; i < edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : edges.getEdges(fw))
		{
			assert(edges.hasEdge(reverse(pair), reverse(fw)));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : edges.getEdges(bw))
		{
			assert(edges.hasEdge(reverse(pair), reverse(bw)));
		}
	}
	auto anchors = getAnchorChains(anchor, edges, unitigGraph, readPaths, approxOneHapCoverage);
	anchor.assign(unitigGraph.nodeCount(), false);
	for (const auto& chain : anchors)
	{
		if (chain.ploidy == 0) continue;
		for (auto node : chain.nodes)
		{
			anchor[node & maskUint64_t] = true;
		}
	}
	edges = getAnchorEdges(unitigGraph, readPaths, anchor, approxOneHapCoverage);
	auto result = getAnchorChains(anchor, edges, unitigGraph, readPaths, approxOneHapCoverage);
	return result;
}
