#include <fstream>
#include <iostream>
#include <cassert>
#include "UnitigGraph.h"
#include "Common.h"

void startUnitig(const size_t startPos, const std::vector<uint64_t>& uniqueLeftEdge, const std::vector<uint64_t>& uniqueRightEdge, std::vector<bool>& belongsToUnitig, std::vector<std::vector<uint64_t>>& unitigs, std::vector<uint64_t>& unitigLeftmostNode, std::vector<uint64_t>& unitigRightmostNode)
{
	assert(!belongsToUnitig[startPos & maskUint64_t]);
	uint64_t unitigStartPos = startPos;
	while (true)
	{
		uint64_t next;
		if (unitigStartPos & firstBitUint64_t)
		{
			next = uniqueLeftEdge[unitigStartPos & maskUint64_t];
		}
		else
		{
			next = uniqueRightEdge[unitigStartPos & maskUint64_t];
		}
		if (next >= std::numeric_limits<size_t>::max()-1) break;
		next ^= firstBitUint64_t;
		if ((next & maskUint64_t) == (unitigStartPos & maskUint64_t)) break;
		if (next & firstBitUint64_t)
		{
			if (uniqueRightEdge[next & maskUint64_t] >= std::numeric_limits<size_t>::max()-1) break;
			assert(uniqueRightEdge[next & maskUint64_t] == unitigStartPos);
		}
		else
		{
			if (uniqueLeftEdge[next & maskUint64_t] >= std::numeric_limits<size_t>::max()-1) break;
			assert(uniqueLeftEdge[next & maskUint64_t] == unitigStartPos);
		}
		if (next == startPos) break;
		assert(!belongsToUnitig[next & maskUint64_t]);
		unitigStartPos = next;
	}
	unitigs.emplace_back();
	unitigs.back().push_back(unitigStartPos);
	unitigLeftmostNode.emplace_back(unitigStartPos ^ firstBitUint64_t);
	uint64_t pos = unitigStartPos;
	belongsToUnitig[pos & maskUint64_t] = true;
	while (true)
	{
		uint64_t next;
		if (pos & firstBitUint64_t)
		{
			next = uniqueRightEdge[pos & maskUint64_t];
		}
		else
		{
			next = uniqueLeftEdge[pos & maskUint64_t];
		}
		if (next >= std::numeric_limits<size_t>::max()-1) break;
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next & firstBitUint64_t)
		{
			if (uniqueLeftEdge[next & maskUint64_t] >= std::numeric_limits<size_t>::max()-1) break;
			assert(uniqueLeftEdge[next & maskUint64_t] == (pos ^ firstBitUint64_t));
		}
		else
		{
			if (uniqueRightEdge[next & maskUint64_t] >= std::numeric_limits<size_t>::max()-1) break;
			assert(uniqueRightEdge[next & maskUint64_t] == (pos ^ firstBitUint64_t));
		}
		if (next == unitigStartPos) break;
		pos = next;
		assert(!belongsToUnitig[pos & maskUint64_t]);
		belongsToUnitig[pos & maskUint64_t] = true;
		unitigs.back().push_back(pos);
	}
	unitigRightmostNode.emplace_back(pos);
}

std::tuple<std::vector<std::vector<uint64_t>>, std::vector<uint64_t>, std::vector<uint64_t>> getUnitigs(const size_t countNodes, const size_t minCoverage, const MostlySparse2DHashmap<uint8_t, size_t>& edgeCoverages)
{
	std::vector<uint64_t> uniqueRightEdge;
	std::vector<uint64_t> uniqueLeftEdge;
	uniqueRightEdge.resize(countNodes, std::numeric_limits<size_t>::max());
	uniqueLeftEdge.resize(countNodes, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < countNodes; i++)
	{
		for (auto pair : edgeCoverages.getValues(std::make_pair(i, true)))
		{
			if (pair.second < minCoverage) continue;
			if (uniqueRightEdge[i] == std::numeric_limits<size_t>::max())
			{
				uniqueRightEdge[i] = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			}
			else
			{
				uniqueRightEdge[i] = std::numeric_limits<size_t>::max()-1;
			}
			if (pair.first.second)
			{
				if (uniqueLeftEdge[pair.first.first] == std::numeric_limits<size_t>::max())
				{
					uniqueLeftEdge[pair.first.first] = i;
				}
				else
				{
					uniqueLeftEdge[pair.first.first] = std::numeric_limits<size_t>::max()-1;
				}
			}
			else
			{
				if (uniqueRightEdge[pair.first.first] == std::numeric_limits<size_t>::max())
				{
					uniqueRightEdge[pair.first.first] = i;
				}
				else
				{
					uniqueRightEdge[pair.first.first] = std::numeric_limits<size_t>::max()-1;
				}
			}
		}
		for (auto pair : edgeCoverages.getValues(std::make_pair(i, false)))
		{
			if (pair.second < minCoverage) continue;
			if (uniqueLeftEdge[i] == std::numeric_limits<size_t>::max())
			{
				uniqueLeftEdge[i] = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			}
			else
			{
				uniqueLeftEdge[i] = std::numeric_limits<size_t>::max()-1;
			}
			if (pair.first.second)
			{
				if (uniqueLeftEdge[pair.first.first] == std::numeric_limits<size_t>::max())
				{
					uniqueLeftEdge[pair.first.first] = i + firstBitUint64_t;
				}
				else
				{
					uniqueLeftEdge[pair.first.first] = std::numeric_limits<size_t>::max()-1;
				}
			}
			else
			{
				if (uniqueRightEdge[pair.first.first] == std::numeric_limits<size_t>::max())
				{
					uniqueRightEdge[pair.first.first] = i + firstBitUint64_t;
				}
				else
				{
					uniqueRightEdge[pair.first.first] = std::numeric_limits<size_t>::max()-1;
				}
			}
		}
	}
	std::vector<bool> belongsToUnitig;
	belongsToUnitig.resize(countNodes, false);
	std::vector<std::vector<uint64_t>> unitigs;
	std::vector<uint64_t> unitigLeftmostNode;
	std::vector<uint64_t> unitigRightmostNode;
	for (size_t i = 0; i < countNodes; i++)
	{
		if (belongsToUnitig[i]) continue;
		startUnitig(i, uniqueLeftEdge, uniqueRightEdge, belongsToUnitig, unitigs, unitigLeftmostNode, unitigRightmostNode);
	}
	return std::make_tuple(std::move(unitigs), std::move(unitigLeftmostNode), std::move(unitigRightmostNode));
}

std::pair<std::vector<size_t>, std::vector<double>> getUnitigLengthAndCoverage(const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<size_t>& nodeCoverage, const std::vector<size_t>& nodeLength)
{
	std::vector<double> coverage;
	std::vector<size_t> length;
	coverage.resize(unitigs.size(), 0);
	length.resize(unitigs.size(), 0);
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		for (uint64_t node : unitigs[i])
		{
			coverage[i] += nodeLength[node & maskUint64_t] * nodeCoverage[node & maskUint64_t];
			length[i] += nodeLength[node & maskUint64_t];
		}
		coverage[i] /= (double)length[i];
	}
	return std::make_pair(std::move(length), std::move(coverage));
}

MostlySparse2DHashmap<uint8_t, size_t> getUnitigEdgeCoverages(const std::vector<uint64_t>& unitigLeftmostNode, const std::vector<uint64_t>& unitigRightmostNode, const size_t minCoverage, const MostlySparse2DHashmap<uint8_t, size_t>& edgeCoverages)
{
	phmap::flat_hash_map<uint64_t, uint64_t> nodeToUnitig;
	assert(unitigLeftmostNode.size() == unitigRightmostNode.size());
	MostlySparse2DHashmap<uint8_t, size_t> result;
	result.resize(unitigLeftmostNode.size());
	for (size_t i = 0; i < unitigLeftmostNode.size(); i++)
	{
		assert(nodeToUnitig.count(unitigLeftmostNode[i] ^ firstBitUint64_t) == 0);
		nodeToUnitig[unitigLeftmostNode[i] ^ firstBitUint64_t] = i + firstBitUint64_t;
		assert(nodeToUnitig.count(unitigRightmostNode[i] ^ firstBitUint64_t) == 0);
		nodeToUnitig[unitigRightmostNode[i] ^ firstBitUint64_t] = i;
	}
	for (size_t i = 0; i < unitigLeftmostNode.size(); i++)
	{
		for (auto pair : edgeCoverages.getValues(std::make_pair(unitigLeftmostNode[i] & maskUint64_t, (unitigLeftmostNode[i] & firstBitUint64_t) == firstBitUint64_t)))
		{
			if (pair.second < minCoverage) continue;
			uint64_t kmerKey = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			assert(nodeToUnitig.count(kmerKey) == 1);
			uint64_t targetUnitig = nodeToUnitig.at(kmerKey);
			result.set(std::make_pair(i, false), std::make_pair(targetUnitig & maskUint64_t, (targetUnitig & firstBitUint64_t) == firstBitUint64_t), pair.second);
		}
		for (auto pair : edgeCoverages.getValues(std::make_pair(unitigRightmostNode[i] & maskUint64_t, (unitigRightmostNode[i] & firstBitUint64_t) == firstBitUint64_t)))
		{
			if (pair.second < minCoverage) continue;
			uint64_t kmerKey = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			assert(nodeToUnitig.count(kmerKey) == 1);
			uint64_t targetUnitig = nodeToUnitig.at(kmerKey);
			result.set(std::make_pair(i, true), std::make_pair(targetUnitig & maskUint64_t, (targetUnitig & firstBitUint64_t) == firstBitUint64_t), pair.second);
		}
	}
	return result;
}

void writeGraph(std::string outputFileName, const UnitigGraph& unitigGraph, const size_t minCoverage, const size_t k)
{
	std::ofstream graph { outputFileName };
	size_t countNodes = 0;
	size_t countEdges = 0;
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.coverages[i] < minCoverage) continue;
		graph << "S\tnode_" << i << "\t*\tLN:i:" << (unitigGraph.lengths[i]+k-1) << "\tll:f:" << unitigGraph.coverages[i] << "\tFC:i:" << unitigGraph.coverages[i] * (unitigGraph.lengths[i]+k-1) << std::endl;
		countNodes += 1;
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		auto fw = std::make_pair(i, true);
		for (auto pair : unitigGraph.edgeCoverages.getValues(fw))
		{
			if (pair.second < minCoverage) continue;
			graph << "L\tnode_" << i << "\t+\tnode_" << pair.first.first << "\t" << (pair.first.second ? "+" : "-") << "\t" << (k-1) << "M\tec:i:" << pair.second << std::endl;
			countEdges += 1;
		}
		auto bw = std::make_pair(i, false);
		for (auto pair : unitigGraph.edgeCoverages.getValues(bw))
		{
			if (pair.second < minCoverage) continue;
			graph << "L\tnode_" << i << "\t-\tnode_" << pair.first.first << "\t" << (pair.first.second ? "+" : "-") << "\t" << (k-1) << "M\tec:i:" << pair.second << std::endl;
			countEdges += 1;
		}
	}
	std::cerr << countNodes << " graph nodes" << std::endl;
	std::cerr << countEdges << " graph edges" << std::endl;
}

UnitigGraph makeUnitigGraph(const KmerGraph& kmerGraph, const size_t minCoverage)
{
	UnitigGraph result;
	std::vector<std::vector<uint64_t>> unitigs;
	std::vector<uint64_t> unitigLeftmostNode;
	std::vector<uint64_t> unitigRightmostNode;
	std::tie(unitigs, unitigLeftmostNode, unitigRightmostNode) = getUnitigs(kmerGraph.nodeCount(), minCoverage, kmerGraph.edgeCoverages);
	std::tie(result.lengths, result.coverages) = getUnitigLengthAndCoverage(unitigs, kmerGraph.coverages, kmerGraph.lengths);
	result.edgeCoverages = getUnitigEdgeCoverages(unitigLeftmostNode, unitigRightmostNode, minCoverage, kmerGraph.edgeCoverages);
	return result;
}

size_t UnitigGraph::nodeCount() const
{
	return coverages.size();
}
