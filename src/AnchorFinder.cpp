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

bool isReachableWithoutNodeEvenBackwards(const uint64_t startNode, const uint64_t endNode, const uint64_t forbiddenNode, const SparseEdgeContainer& edges, const std::vector<bool>& anchor)
{
	phmap::flat_hash_set<uint64_t> checked;
	std::vector<uint64_t> stack;
	stack.push_back(startNode);
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (checked.count(top) == 1) continue;
		if ((top & maskUint64_t) == (endNode & maskUint64_t)) return true;
		checked.insert(top);
		if ((top & maskUint64_t) != startNode && (top & maskUint64_t) != (forbiddenNode & maskUint64_t) && !anchor[top & maskUint64_t]) stack.push_back(top ^ firstBitUint64_t);
		for (const auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	if (forbiddenNode != std::numeric_limits<size_t>::max())
	{
		assert(checked.count(forbiddenNode) == 1 || checked.count(forbiddenNode ^ firstBitUint64_t) == 1);
	}
	return false;
}

bool blocksGraphPath(const UnitigGraph& unitigGraph, const SparseEdgeContainer& edges, const uint64_t startNode, const uint64_t midNode, const uint64_t endNode, const std::vector<bool>& anchor)
{
	if (!isReachableWithoutNodeEvenBackwards(startNode, endNode, std::numeric_limits<size_t>::max(), edges, anchor)) return false;
	if (isReachableWithoutNodeEvenBackwards(startNode, endNode, midNode, edges, anchor)) return false;
	return true;
}

SparseEdgeContainer getEdgesWithoutTransitiveEdges(const SparseEdgeContainer& rawEdgesWithTransitiveEdges, const std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>& maybeTransitiveEdges, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, std::vector<bool>& anchor, const double approxOneHapCoverage)
{
	if (maybeTransitiveEdges.size() == 0) return rawEdgesWithTransitiveEdges;
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> forbiddenTransitiveEdges;
	for (auto t : maybeTransitiveEdges)
	{
		assert(anchor[std::get<1>(t) & maskUint64_t]);
		anchor[std::get<1>(t) & maskUint64_t] = false;
		if (!blocksGraphPath(unitigGraph, edges, std::get<0>(t), std::get<1>(t), std::get<2>(t), anchor))
		{
			anchor[std::get<1>(t) & maskUint64_t] = true;
			continue;
		}
		anchor[std::get<1>(t) & maskUint64_t] = true;
		forbiddenTransitiveEdges.emplace(std::get<0>(t), std::get<2>(t));
		std::cerr << "forbid transitive edge " << ((std::get<0>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<0>(t) & maskUint64_t) << " " << ((std::get<2>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(t) & maskUint64_t) << std::endl;
	}
	if (forbiddenTransitiveEdges.size() == 0) return rawEdgesWithTransitiveEdges;
	SparseEdgeContainer result;
	result.resize(rawEdgesWithTransitiveEdges.size());
	for (size_t i = 0; i < rawEdgesWithTransitiveEdges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (const auto edge : rawEdgesWithTransitiveEdges.getEdges(fw))
		{
			if (forbiddenTransitiveEdges.count(std::make_pair(fw.first + (fw.second ? firstBitUint64_t : 0), edge.first + (edge.second ? firstBitUint64_t : 0))) == 1) continue;
			if (forbiddenTransitiveEdges.count(std::make_pair(edge.first + (edge.second ? 0 : firstBitUint64_t), fw.first + (fw.second ? 0 : firstBitUint64_t))) == 1) continue;
			result.addEdge(fw, edge);
		}
		std::pair<size_t, bool> bw { i, false };
		for (const auto edge : rawEdgesWithTransitiveEdges.getEdges(bw))
		{
			if (forbiddenTransitiveEdges.count(std::make_pair(bw.first + (bw.second ? firstBitUint64_t : 0), edge.first + (edge.second ? firstBitUint64_t : 0))) == 1) continue;
			if (forbiddenTransitiveEdges.count(std::make_pair(edge.first + (edge.second ? 0 : firstBitUint64_t), bw.first + (bw.second ? 0 : firstBitUint64_t))) == 1) continue;
			result.addEdge(bw, edge);
		}
	}
	return result;
}

SparseEdgeContainer getAnchorEdges(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, std::vector<bool>& anchor, const double approxOneHapCoverage)
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
			if ((double)pair.second < approxOneHapCoverage * 0.5 && (pair.second < unitigGraph.coverages[i] * 0.5 || pair.second < unitigGraph.coverages[pair.first.first] * 0.5)) continue;
			assert(anchor[pair.first.first]);
			assert(fw.first < rawEdgesWithTransitiveEdges.size());
			assert(pair.first.first < rawEdgesWithTransitiveEdges.size());
			rawEdgesWithTransitiveEdges.addEdge(fw, pair.first);
			rawEdgesWithTransitiveEdges.addEdge(reverse(pair.first), reverse(fw));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : edgeCoverage.getValues(bw))
		{
			if ((double)pair.second < approxOneHapCoverage * 0.5 && (pair.second < unitigGraph.coverages[i] * 0.5 || pair.second < unitigGraph.coverages[pair.first.first] * 0.5)) continue;
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

void setReachableAnchors(const SparseEdgeContainer& edges, std::vector<std::vector<uint64_t>>& forwardReachableAnchors, std::vector<std::vector<uint64_t>>& backwardReachableAnchors, const uint64_t startNode, const std::vector<bool>& anchor)
{
	std::vector<uint64_t> stack;
	phmap::flat_hash_set<uint64_t> checked;
	for (auto edge : edges.getEdges(std::make_pair(startNode & maskUint64_t, startNode & firstBitUint64_t)))
	{
		stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
	}
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (checked.count(top) == 1) continue;
		if (anchor[top & maskUint64_t]) continue;
		checked.insert(top);
		if (top & firstBitUint64_t)
		{
			backwardReachableAnchors[top & maskUint64_t].emplace_back(startNode);
		}
		else
		{
			forwardReachableAnchors[top & maskUint64_t].emplace_back(startNode);
		}
		for (const auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
		}
	}
}

void extendAnchors(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, std::vector<bool>& anchor)
{
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<std::vector<uint64_t>> forwardReachableAnchors;
	std::vector<std::vector<uint64_t>> backwardReachableAnchors;
	forwardReachableAnchors.resize(unitigGraph.nodeCount());
	backwardReachableAnchors.resize(unitigGraph.nodeCount());
	std::vector<bool> anchorsBeforeExtension = anchor;
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (!anchor[i]) continue;
		setReachableAnchors(edges, forwardReachableAnchors, backwardReachableAnchors, i + firstBitUint64_t, anchor);
		setReachableAnchors(edges, forwardReachableAnchors, backwardReachableAnchors, i, anchor);
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (anchor[i]) continue;
		if (forwardReachableAnchors[i].size() != 1 && backwardReachableAnchors[i].size() != 1) continue;
		if (forwardReachableAnchors[i].size() == 0 || backwardReachableAnchors[i].size() == 0) continue;
		bool anyNotBlocked = false;
		for (size_t j = 0; j < forwardReachableAnchors[i].size(); j++)
		{
			for (size_t k = 0; k < backwardReachableAnchors[i].size(); k++)
			{
				if (!blocksGraphPath(unitigGraph, edges, backwardReachableAnchors[i][k], i, forwardReachableAnchors[i][j] ^ firstBitUint64_t, anchorsBeforeExtension))
				{
					anyNotBlocked = true;
					break;
				}
			}
			if (anyNotBlocked) break;
		}
		if (anyNotBlocked) continue;
		anchor[i] = true;
	}
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

void getDistances(phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<int>>& distances, const std::vector<bool> anchor, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths)
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
					int distance = readOffset - lastAnchorOffset + lastAnchorPlus;
					assert(distance >= 0);
					if (k == 0) distance -= (int)readPaths[i].paths[j].pathLeftClipKmers;
					distances.at(std::make_pair(lastAnchor, node)).emplace_back(distance);
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
					int distance = readOffset - lastAnchorOffset + lastAnchorPlus;
					assert(distance >= 0);
					if (k == 0) distance -= (int)readPaths[i].paths[j].pathLeftClipKmers;
					distance += (int)unitigGraph.lengths[node & maskUint64_t] - (int)unitigGraph.lengths[lastAnchor & maskUint64_t];
					distances.at(std::make_pair(node ^ firstBitUint64_t, lastAnchor ^ firstBitUint64_t)).emplace_back(distance);
					lastAnchor = node;
					lastAnchorOffset = readOffset;
					lastAnchorPlus = 0;
					if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
					readOffset += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
					continue;
				}
				lastAnchor = node;
				lastAnchorOffset = readOffset;
				lastAnchorPlus = 0;
				if (k == 0) lastAnchorPlus = readPaths[i].paths[j].pathLeftClipKmers;
				readOffset += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readOffset -= readPaths[i].paths[j].pathLeftClipKmers;
			}
		}
	}
}

double average(const std::vector<int>& values)
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
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<int>> distances;
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
			size_t distance = average(distances.at(key));
			std::cerr << "chain " << i << " anchor " << j << " (" << (chains[i].nodes[j] & maskUint64_t) << ") distance " << distance << std::endl;
			assert(distance >= 0);
			assert(distance < 100000000);
			chains[i].nodeOffsets.push_back(chains[i].nodeOffsets.back() + distance);
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

SparseEdgeContainer getAnchorChainEdges(const std::vector<AnchorChain>& anchorChains, const SparseEdgeContainer& anchorEdges)
{
	SparseEdgeContainer result;
	result.resize(anchorChains.size());
	phmap::flat_hash_map<uint64_t, uint64_t> anchorStart;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		assert(anchorStart.count(anchorChains[i].nodes[0]) == 0);
		anchorStart[anchorChains[i].nodes[0]] = i + firstBitUint64_t;
		assert(anchorStart.count(anchorChains[i].nodes.back() ^ firstBitUint64_t) == 0);
		anchorStart[anchorChains[i].nodes.back() ^ firstBitUint64_t] = i;
	}
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::pair<size_t, bool> fw { anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t };
		for (auto edge : anchorEdges.getEdges(fw))
		{
			assert(anchorStart.count(edge.first + (edge.second ? firstBitUint64_t : 0)));
			uint64_t target = anchorStart.at(edge.first + (edge.second ? firstBitUint64_t : 0));
			result.addEdge(std::make_pair(i, true), std::make_pair(target & maskUint64_t, target & firstBitUint64_t));
		}
		std::pair<size_t, bool> bw { anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] & firstBitUint64_t) ^ firstBitUint64_t };
		for (auto edge : anchorEdges.getEdges(bw))
		{
			assert(anchorStart.count(edge.first + (edge.second ? firstBitUint64_t : 0)));
			uint64_t target = anchorStart.at(edge.first + (edge.second ? firstBitUint64_t : 0));
			result.addEdge(std::make_pair(i, false), std::make_pair(target & maskUint64_t, target & firstBitUint64_t));
		}
	}
	return result;
}

uint64_t getBubble(const uint64_t startPos, const SparseEdgeContainer& edges)
{
	std::vector<std::pair<size_t, bool>> S;
	phmap::flat_hash_set<std::pair<size_t, bool>> seen;
	phmap::flat_hash_set<std::pair<size_t, bool>> visited;
	S.emplace_back(startPos & maskUint64_t, startPos & firstBitUint64_t);
	seen.emplace(startPos & maskUint64_t, startPos & firstBitUint64_t);
	while (S.size() > 0)
	{
		auto v = S.back();
		S.pop_back();
		assert(seen.count(v) == 1);
		seen.erase(v);
		assert(visited.count(v) == 0);
		visited.insert(v);
		if (edges.getEdges(v).size() == 0) return std::numeric_limits<size_t>::max();
		for (auto u : edges.getEdges(v))
		{
			if (u.first == v.first) return std::numeric_limits<size_t>::max();
			if (u.first + (u.second ? firstBitUint64_t : 0) == startPos) return std::numeric_limits<size_t>::max();
			assert(visited.count(u) == 0);
			seen.insert(u);
			bool hasNonvisitedParent = false;
			assert(edges.getEdges(reverse(u)).size() >= 1);
			for (auto edge : edges.getEdges(reverse(u)))
			{
				if (visited.count(reverse(edge)) == 0) hasNonvisitedParent = true;
			}
			if (!hasNonvisitedParent) S.push_back(u);
		}
		if (S.size() == 1 && seen.size() == 1 && S[0] == *seen.begin())
		{
			auto t = S[0];
			for (auto edge : edges.getEdges(t))
			{
				if (edge.first + (edge.second ? firstBitUint64_t : 0) == startPos) return std::numeric_limits<size_t>::max();
			}
			return t.first + (t.second ? firstBitUint64_t : 0);
		}
	}
	return std::numeric_limits<size_t>::max();
}

std::vector<uint64_t> getChainOfChains(const uint64_t start, const SparseEdgeContainer& anchorEdges)
{
	uint64_t pos = start;
	while (true)
	{
		uint64_t next = getBubble(pos, anchorEdges);
		if (next == std::numeric_limits<size_t>::max()) break;
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next == start) break;
		pos = next;
	}
	pos = pos ^ firstBitUint64_t;
	const uint64_t resultStart = pos;
	std::vector<uint64_t> result;
	result.emplace_back(pos);
	while (true)
	{
		uint64_t next = getBubble(pos, anchorEdges);
		if (next == std::numeric_limits<size_t>::max()) break;
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next == resultStart) break;
		pos = next;
		result.emplace_back(pos);
	}
	return result;
}

std::vector<uint64_t> getSimpleBubble(const uint64_t bubbleStart, const uint64_t bubbleEnd, const SparseEdgeContainer& edges)
{
	std::vector<uint64_t> result;
	for (auto edge : edges.getEdges(std::make_pair(bubbleStart & maskUint64_t, bubbleStart & firstBitUint64_t)))
	{
		if (edges.getEdges(edge).size() != 1) return std::vector<uint64_t>{};
		if (edges.getEdges(edge)[0] != std::make_pair<size_t, bool>(bubbleEnd & maskUint64_t, bubbleEnd & firstBitUint64_t)) return std::vector<uint64_t>{};
		if (edges.getEdges(reverse(edge)).size() != 1) return std::vector<uint64_t>{};
		result.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
	}
	return result;
}

void refineAnchorChainPloidies(std::vector<AnchorChain>& anchorChains, const SparseEdgeContainer& anchorEdges, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	SparseEdgeContainer anchorChainEdges = getAnchorChainEdges(anchorChains, anchorEdges);
	std::vector<bool> checked;
	std::vector<bool> refined;
	refined.resize(anchorChains.size(), false);
	checked.resize(anchorChains.size(), false);
	std::vector<std::vector<uint64_t>> chainsOfChains;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (checked[i]) continue;
		auto chain = getChainOfChains(i, anchorChainEdges);
		for (auto node : chain)
		{
			assert((node & maskUint64_t) < checked.size());
			assert(!checked[node & maskUint64_t]);
			checked[node & maskUint64_t] = true;
		}
		chainsOfChains.emplace_back(chain);
	}
	std::sort(chainsOfChains.begin(), chainsOfChains.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
	for (size_t i = 0; i < chainsOfChains.size(); i++)
	{
		if (chainsOfChains[i].size() == 1 && refined[chainsOfChains[i][0] & maskUint64_t]) continue;
		double coverageSum = 0;
		size_t coverageDivisor = 0;
		for (const uint64_t chain : chainsOfChains[i])
		{
			for (const uint64_t node : anchorChains[chain & maskUint64_t].nodes)
			{
				coverageSum += unitigGraph.coverages[node & maskUint64_t] * unitigGraph.lengths[node & maskUint64_t];
				coverageDivisor += unitigGraph.lengths[node & maskUint64_t];
			}
		}
		size_t refinedPloidy = coverageSum / (double)coverageDivisor / approxOneHapCoverage + 0.5;
		for (const uint64_t chain : chainsOfChains[i])
		{
			if (anchorChains[chain & maskUint64_t].ploidy != refinedPloidy)
			{
				std::cerr << "refined chain of chains core chain " << (anchorChains[chain & maskUint64_t].nodes[0] & maskUint64_t) << "-" << (anchorChains[chain & maskUint64_t].nodes.back() & maskUint64_t) << " ploidy from " << anchorChains[chain & maskUint64_t].ploidy << " to " << refinedPloidy << std::endl;
			}
			anchorChains[chain & maskUint64_t].ploidy = refinedPloidy;
			refined[chain & maskUint64_t] = true;
		}
		for (size_t j = 0; j+1 < chainsOfChains[i].size(); j++)
		{
			auto bubbleChains = getSimpleBubble(chainsOfChains[i][j], chainsOfChains[i][j+1], anchorChainEdges);
			if (bubbleChains.size() == 0)
			{
				std::cerr << "could not make a refined coverage estimate for chain of anchor chains bubble between chains " << (anchorChains[chainsOfChains[i][j-1]].nodes[0] & maskUint64_t) << "-" << (anchorChains[chainsOfChains[i][j-1]].nodes.back() & maskUint64_t) << " and " << (anchorChains[chainsOfChains[i][j]].nodes[0] & maskUint64_t) << "-" << (anchorChains[chainsOfChains[i][j]].nodes.back() & maskUint64_t) << std::endl;
			}
			else
			{
				if (bubbleChains.size() == refinedPloidy)
				{
					for (const uint64_t chain : bubbleChains)
					{
						if (anchorChains[chain & maskUint64_t].ploidy != 1)
						{
							std::cerr << "refined chain of chains bubble chain " << (anchorChains[chain & maskUint64_t].nodes[0] & maskUint64_t) << "-" << (anchorChains[chain & maskUint64_t].nodes.back() & maskUint64_t) << " ploidy from " << anchorChains[chain & maskUint64_t].ploidy << " to " << 1 << std::endl;
						}
						anchorChains[chain & maskUint64_t].ploidy = 1;
						refined[chain & maskUint64_t] = true;
					}
				}
				else if (bubbleChains.size() > refinedPloidy)
				{
					std::cerr << "could not distribute ploidy between too many chains (ploidy " << refinedPloidy << " vs chain count" << bubbleChains.size() << ")";
					for (const uint64_t chain : bubbleChains)
					{
						std::cerr << " " << (anchorChains[chain & maskUint64_t].nodes[0] & maskUint64_t) << "-" << (anchorChains[chain & maskUint64_t].nodes.back() & maskUint64_t);
					}
					std::cerr << std::endl;
				}
				else
				{
					assert(bubbleChains.size() < refinedPloidy);
					std::cerr << "todo implement: could not distribute ploidy between too few chains (ploidy " << refinedPloidy << " vs chain count" << bubbleChains.size() << ")";
					for (const uint64_t chain : bubbleChains)
					{
						std::cerr << " " << (anchorChains[chain & maskUint64_t].nodes[0] & maskUint64_t) << "-" << (anchorChains[chain & maskUint64_t].nodes.back() & maskUint64_t);
					}
					std::cerr << std::endl;
				}
			}
		}
	}
}

std::vector<AnchorChain> getAnchorChains(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.lengths[i] >= 100) anchor[i] = true;
	}
	extendAnchors(unitigGraph, readPaths, anchor);
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
	for (size_t i = 0; i < anchors.size(); i++)
	{
		std::cerr << "raw chain " << i;
		for (uint64_t node : anchors[i].nodes)
		{
			std::cerr << " " << (node & maskUint64_t);
		}
		std::cerr << std::endl;
	}
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
	refineAnchorChainPloidies(result, edges, unitigGraph, readPaths, approxOneHapCoverage);
	return result;
}

std::vector<std::vector<ChainPosition>> getReadChainPositions(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains)
{
	std::vector<std::pair<size_t, size_t>> anchorLocator;
	anchorLocator.resize(unitigGraph.nodeCount(), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			anchorLocator[anchorChains[i].nodes[j] & maskUint64_t] = std::make_pair(i, j);
		}
	}
	std::vector<std::vector<ChainPosition>> result;
	result.resize(readPaths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, bool, int, size_t, int, int>> impliedPositions; // chain, fw, diagonal, readPos, readChainStartPos, readChainEndPos
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (anchorLocator[node & maskUint64_t].first == std::numeric_limits<size_t>::max()) continue;
				const size_t chain = anchorLocator[node & maskUint64_t].first;
				const size_t offset = anchorLocator[node & maskUint64_t].second;
				assert((node & maskUint64_t) == (anchorChains[chain].nodes[offset] & maskUint64_t));
				const bool fw = (node & firstBitUint64_t) == (anchorChains[chain].nodes[offset] & firstBitUint64_t);
				int diagonal = (int)readPos - (int)anchorChains[chain].nodeOffsets[offset] - (int)unitigGraph.lengths[anchorChains[chain].nodes[offset] & maskUint64_t];
				if (!fw) diagonal = (int)readPos + (int)anchorChains[chain].nodeOffsets[offset];
				if (fw)
				{
					// readPos is at the end of current node
					impliedPositions.emplace_back(chain, fw, diagonal, readPos, readPos - unitigGraph.lengths[node & maskUint64_t] - anchorChains[chain].nodeOffsets[offset], readPos + anchorChains[chain].nodeOffsets.back() - anchorChains[chain].nodeOffsets[offset] - unitigGraph.lengths[node & maskUint64_t] + unitigGraph.lengths[anchorChains[chain].nodes.back() & maskUint64_t]);
				}
				else
				{
					// readPos is at the start of current node
					impliedPositions.emplace_back(chain, fw, diagonal, readPos, readPos - (unitigGraph.lengths[anchorChains[chain].nodes.back() & maskUint64_t] + anchorChains[chain].nodeOffsets.back() - anchorChains[chain].nodeOffsets[offset]), readPos + anchorChains[chain].nodeOffsets[offset]);
				}
			}
		}
		if (impliedPositions.size() == 0) continue;
		std::sort(impliedPositions.begin(), impliedPositions.end());
		size_t currentMinReadPos = std::numeric_limits<size_t>::max();
		size_t currentMaxReadPos = 0;
		for (size_t j = 0; j < impliedPositions.size(); j++)
		{
			if (j == 0 || std::get<0>(impliedPositions[j]) != std::get<0>(impliedPositions[j-1]) || std::get<1>(impliedPositions[j]) != std::get<1>(impliedPositions[j-1]) || std::get<2>(impliedPositions[j]) > std::get<2>(impliedPositions[j-1]) + (anchorChains[std::get<0>(impliedPositions[j])].nodeOffsets.back() + unitigGraph.lengths[anchorChains[std::get<0>(impliedPositions[j])].nodes.back() & maskUint64_t])*0.5)
			{
				result[i].emplace_back();
				result[i].back().chain = std::get<0>(impliedPositions[j]) + (std::get<1>(impliedPositions[j]) ? firstBitUint64_t : 0);
				currentMinReadPos = std::numeric_limits<size_t>::max();
				currentMaxReadPos = 0;
			}
			if (std::get<3>(impliedPositions[j]) < currentMinReadPos)
			{
				result[i].back().chainStartPosInRead = std::get<4>(impliedPositions[j]);
				currentMinReadPos = std::get<3>(impliedPositions[j]);
			}
			if (std::get<3>(impliedPositions[j]) > currentMaxReadPos)
			{
				result[i].back().chainEndPosInRead = std::get<5>(impliedPositions[j]);
				currentMaxReadPos = std::get<3>(impliedPositions[j]);
			}
		}
		std::sort(result[i].begin(), result[i].end(), [](auto left, auto right) { return left.chainStartPosInRead < right.chainStartPosInRead; });
	}
	return result;
}
