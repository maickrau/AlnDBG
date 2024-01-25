#include <fstream>
#include <iostream>
#include <cassert>
#include "MBGCommon.h"
#include "UnitigGraph.h"
#include "Common.h"

void startUnitig(const size_t startPos, const VectorWithDirection<uint64_t>& uniqueEdges, std::vector<bool>& belongsToUnitig, std::vector<std::vector<uint64_t>>& unitigs)
{
	assert(!belongsToUnitig[startPos & maskUint64_t]);
	uint64_t unitigStartPos = startPos;
	while (true)
	{
		uint64_t next;
		if (unitigStartPos & firstBitUint64_t)
		{
			next = uniqueEdges[std::make_pair(unitigStartPos & maskUint64_t, false)];
		}
		else
		{
			next = uniqueEdges[std::make_pair(unitigStartPos & maskUint64_t, true)];
		}
		if (next >= std::numeric_limits<size_t>::max()-1) break;
		next ^= firstBitUint64_t;
		if ((next & maskUint64_t) == (unitigStartPos & maskUint64_t)) break;
		if (next & firstBitUint64_t)
		{
			if (uniqueEdges[std::make_pair(next & maskUint64_t, true)] >= std::numeric_limits<size_t>::max()-1) break;
			assert(uniqueEdges[std::make_pair(next & maskUint64_t, true)] == unitigStartPos);
		}
		else
		{
			if (uniqueEdges[std::make_pair(next & maskUint64_t, false)] >= std::numeric_limits<size_t>::max()-1) break;
			assert(uniqueEdges[std::make_pair(next & maskUint64_t, false)] == unitigStartPos);
		}
		if (next == startPos) break;
		assert(!belongsToUnitig[next & maskUint64_t]);
		unitigStartPos = next;
	}
	unitigs.emplace_back();
	unitigs.back().push_back(unitigStartPos);
	uint64_t pos = unitigStartPos;
	belongsToUnitig[pos & maskUint64_t] = true;
	while (true)
	{
		uint64_t next;
		if (pos & firstBitUint64_t)
		{
			next = uniqueEdges[std::make_pair(pos & maskUint64_t, true)];
		}
		else
		{
			next = uniqueEdges[std::make_pair(pos & maskUint64_t, false)];
		}
		if (next >= std::numeric_limits<size_t>::max()-1) break;
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next & firstBitUint64_t)
		{
			if (uniqueEdges[std::make_pair(next & maskUint64_t, false)] >= std::numeric_limits<size_t>::max()-1) break;
			assert(uniqueEdges[std::make_pair(next & maskUint64_t, false)] == (pos ^ firstBitUint64_t));
		}
		else
		{
			if (uniqueEdges[std::make_pair(next & maskUint64_t, true)] >= std::numeric_limits<size_t>::max()-1) break;
			assert(uniqueEdges[std::make_pair(next & maskUint64_t, true)] == (pos ^ firstBitUint64_t));
		}
		if (next == unitigStartPos) break;
		pos = next;
		assert(!belongsToUnitig[pos & maskUint64_t]);
		belongsToUnitig[pos & maskUint64_t] = true;
		unitigs.back().push_back(pos);
	}
}

std::vector<std::vector<uint64_t>> getUnitigs(const size_t countNodes, const VectorWithDirection<uint64_t>& uniqueEdges)
{
	std::vector<bool> belongsToUnitig;
	belongsToUnitig.resize(countNodes, false);
	std::vector<std::vector<uint64_t>> unitigs;
	for (size_t i = 0; i < countNodes; i++)
	{
		if (belongsToUnitig[i]) continue;
		startUnitig(i, uniqueEdges, belongsToUnitig, unitigs);
	}
	return unitigs;
}

std::vector<std::vector<uint64_t>> getUnitigs(const size_t countNodes, const RankBitvector& keptNodes, const VectorWithDirection<uint64_t>& uniqueEdges)
{
	std::vector<bool> belongsToUnitig;
	belongsToUnitig.resize(countNodes, false);
	std::vector<std::vector<uint64_t>> unitigs;
	for (size_t i = 0; i < countNodes; i++)
	{
		if (!keptNodes.get(i)) continue;
		if (belongsToUnitig[i]) continue;
		startUnitig(i, uniqueEdges, belongsToUnitig, unitigs);
	}
	return unitigs;
}

std::pair<std::vector<size_t>, std::vector<double>> getUnitigLengthAndCoverage(const std::vector<std::vector<uint64_t>>& unitigs, const std::vector<ReadPathBundle>& kmerGraphReadPaths, const std::vector<std::pair<uint64_t, size_t>>& kmerNodeToUnitig, const std::vector<size_t>& nodeLength)
{
	std::vector<double> coverage;
	std::vector<size_t> length;
	coverage.resize(unitigs.size(), 0);
	length.resize(unitigs.size(), 0);
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		for (uint64_t node : unitigs[i])
		{
			length[i] += nodeLength[node & maskUint64_t];
		}
	}
	for (size_t i = 0; i < kmerGraphReadPaths.size(); i++)
	{
		for (size_t j = 0; j < kmerGraphReadPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < kmerGraphReadPaths[i].paths[j].path.size(); k++)
			{
				size_t node = kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t;
				assert(node < kmerNodeToUnitig.size());
				size_t unitig = kmerNodeToUnitig[node].first;
				if (unitig == std::numeric_limits<uint64_t>::max()) continue;
				unitig = unitig & maskUint64_t;
				assert(unitig < coverage.size());
				assert(node < nodeLength.size());
				coverage[unitig] += nodeLength[node];
				if (k == 0 && kmerGraphReadPaths[i].paths[j].pathLeftClipKmers > 0)
				{
					assert(nodeLength[node] > kmerGraphReadPaths[i].paths[j].pathLeftClipKmers);
					coverage[unitig] -= kmerGraphReadPaths[i].paths[j].pathLeftClipKmers;
				}
				if (k == kmerGraphReadPaths[i].paths[j].path.size()-1 && kmerGraphReadPaths[i].paths[j].pathRightClipKmers > 0)
				{
					assert(nodeLength[node] > kmerGraphReadPaths[i].paths[j].pathRightClipKmers);
					coverage[unitig] -= kmerGraphReadPaths[i].paths[j].pathRightClipKmers;
				}
			}
		}
	}
	for (size_t i = 0; i < kmerGraphReadPaths.size(); i++)
	{
		for (size_t j = 1; j < kmerGraphReadPaths[i].paths.size(); j++)
		{
			if (kmerGraphReadPaths[i].paths[j-1].path.back() != kmerGraphReadPaths[i].paths[j].path[0]) continue;
			size_t thisNode = kmerGraphReadPaths[i].paths[j].path[0] & maskUint64_t;
			assert(thisNode < kmerNodeToUnitig.size());
			size_t unitig = kmerNodeToUnitig[thisNode].first;
			if (unitig == std::numeric_limits<uint64_t>::max()) continue;
			unitig = unitig & maskUint64_t;
			assert(thisNode < nodeLength.size());
			assert(unitig < coverage.size());
			size_t prevNodeEnd = nodeLength[thisNode] - kmerGraphReadPaths[i].paths[j-1].pathRightClipKmers;
			size_t currNodeStart = kmerGraphReadPaths[i].paths[j].pathLeftClipKmers;
			if (currNodeStart <= prevNodeEnd) continue;
			size_t currReadStart = kmerGraphReadPaths[i].paths[j].readStartPos;
			size_t prevReadEnd = kmerGraphReadPaths[i].paths[j-1].readStartPos;
			for (auto node : kmerGraphReadPaths[i].paths[j-1].path)
			{
				prevReadEnd += nodeLength[node & maskUint64_t];
			}
			prevReadEnd -= kmerGraphReadPaths[i].paths[j-1].pathLeftClipKmers;
			prevReadEnd -= kmerGraphReadPaths[i].paths[j-1].pathRightClipKmers;
			assert(prevReadEnd <= currReadStart);
			size_t readDistance = currReadStart - prevReadEnd;
			size_t nodeDistance = currNodeStart - prevNodeEnd;
			if (readDistance > nodeDistance + 50) continue;
			if (nodeDistance > readDistance + 50) continue;
			coverage[unitig] += nodeDistance;
		}
	}
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		assert(length[i] != 0);
		// assert(coverage[i] != 0);
		coverage[i] /= (double)length[i];
	}
	return std::make_pair(std::move(length), std::move(coverage));
}

MostlySparse2DHashmap<uint8_t, size_t> getUnitigEdgeCoverages(const std::vector<std::vector<uint64_t>>& unitigs, const size_t minCoverage, const std::vector<ReadPathBundle>& kmerGraphReadPaths, const std::vector<std::pair<uint64_t, size_t>>& kmerNodeToUnitig)
{
	MostlySparse2DHashmap<uint8_t, size_t> result;
	result.resize(unitigs.size());
	for (size_t i = 0; i < kmerGraphReadPaths.size(); i++)
	{
		for (size_t j = 0; j < kmerGraphReadPaths[i].paths.size(); j++)
		{
			for (size_t k = 1; k < kmerGraphReadPaths[i].paths[j].path.size(); k++)
			{
				std::pair<uint64_t, size_t> beforePos = kmerNodeToUnitig[kmerGraphReadPaths[i].paths[j].path[k-1] & maskUint64_t];
				if (beforePos.first == std::numeric_limits<uint64_t>::max()) continue;
				if ((kmerGraphReadPaths[i].paths[j].path[k-1] & firstBitUint64_t) == 0)
				{
					beforePos.first ^= firstBitUint64_t;
					beforePos.second = unitigs[beforePos.first & maskUint64_t].size() - 1 - beforePos.second;
				}
				std::pair<uint64_t, size_t> afterPos = kmerNodeToUnitig[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
				if (afterPos.first == std::numeric_limits<uint64_t>::max()) continue;
				if ((kmerGraphReadPaths[i].paths[j].path[k] & firstBitUint64_t) == 0)
				{
					afterPos.first ^= firstBitUint64_t;
					afterPos.second = unitigs[afterPos.first & maskUint64_t].size() - 1 - afterPos.second;
				}
				if (beforePos.second != unitigs[beforePos.first & maskUint64_t].size()-1) continue;
				if (afterPos.second != 0) continue;
				std::pair<size_t, bool> from { beforePos.first & maskUint64_t, beforePos.first & firstBitUint64_t };
				std::pair<size_t, bool> to { afterPos.first & maskUint64_t, afterPos.first & firstBitUint64_t };
				auto key = canon(from, to);
				size_t coverage = 0;
				if (result.hasValue(key.first, key.second)) coverage = result.get(key.first, key.second);
				coverage += 1;
				result.set(key.first, key.second, coverage);
			}
		}
	}
	if (minCoverage >= 1)
	{
		MostlySparse2DHashmap<uint8_t, size_t> filteredResult;
		filteredResult.resize(result.size());
		for (size_t i = 0; i < unitigs.size(); i++)
		{
			for (auto pair : result.getValues(std::make_pair(i, true)))
			{
				if (pair.second < minCoverage) continue;
				filteredResult.set(std::make_pair(i, true), pair.first, pair.second);
			}
			for (auto pair : result.getValues(std::make_pair(i, false)))
			{
				if (pair.second < minCoverage) continue;
				filteredResult.set(std::make_pair(i, false), pair.first, pair.second);
			}
		}
		std::swap(result,filteredResult);
	}
	return result;
}

void writeGraph(std::string outputFileName, const UnitigGraph& unitigGraph, const size_t k)
{
	std::ofstream graph { outputFileName };
	size_t countNodes = 0;
	size_t countEdges = 0;
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		graph << "S\tnode_" << i << "\t*\tLN:i:" << (unitigGraph.lengths[i]+k-1) << "\tll:f:" << unitigGraph.coverages[i] << "\tFC:i:" << unitigGraph.coverages[i] * (unitigGraph.lengths[i]+k-1) << std::endl;
		countNodes += 1;
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		auto fw = std::make_pair(i, true);
		for (auto pair : unitigGraph.edgeCoverages.getValues(fw))
		{
			graph << "L\tnode_" << i << "\t+\tnode_" << pair.first.first << "\t" << (pair.first.second ? "+" : "-") << "\t" << (k-1) << "M\tec:i:" << pair.second << std::endl;
			countEdges += 1;
		}
		auto bw = std::make_pair(i, false);
		for (auto pair : unitigGraph.edgeCoverages.getValues(bw))
		{
			graph << "L\tnode_" << i << "\t-\tnode_" << pair.first.first << "\t" << (pair.first.second ? "+" : "-") << "\t" << (k-1) << "M\tec:i:" << pair.second << std::endl;
			countEdges += 1;
		}
	}
	std::cerr << countNodes << " graph nodes" << std::endl;
	std::cerr << countEdges << " graph edges" << std::endl;
}

void writeGraph(std::string outputFileName, const UnitigGraph& unitigGraph, const std::vector<TwobitString>& nodeSequences, const size_t k)
{
	std::ofstream graph { outputFileName };
	size_t countNodes = 0;
	size_t countEdges = 0;
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		assert(nodeSequences[i].size() == (unitigGraph.lengths[i]+k-1));
		graph << "S\tnode_" << i << "\t" << nodeSequences[i].toString() << "\tLN:i:" << (unitigGraph.lengths[i]+k-1) << "\tll:f:" << unitigGraph.coverages[i] << "\tFC:i:" << unitigGraph.coverages[i] * (unitigGraph.lengths[i]+k-1) << std::endl;
		countNodes += 1;
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		auto fw = std::make_pair(i, true);
		for (auto pair : unitigGraph.edgeCoverages.getValues(fw))
		{
			graph << "L\tnode_" << i << "\t+\tnode_" << pair.first.first << "\t" << (pair.first.second ? "+" : "-") << "\t" << (k-1) << "M\tec:i:" << pair.second << std::endl;
			countEdges += 1;
		}
		auto bw = std::make_pair(i, false);
		for (auto pair : unitigGraph.edgeCoverages.getValues(bw))
		{
			graph << "L\tnode_" << i << "\t-\tnode_" << pair.first.first << "\t" << (pair.first.second ? "+" : "-") << "\t" << (k-1) << "M\tec:i:" << pair.second << std::endl;
			countEdges += 1;
		}
	}
	std::cerr << countNodes << " graph nodes" << std::endl;
	std::cerr << countEdges << " graph edges" << std::endl;
}

VectorWithDirection<uint64_t> getUniqueEdges(const MostlySparse2DHashmap<uint8_t, size_t>& edgeCoverage, const RankBitvector& keptNodes, const size_t numKmerNodes, const size_t minCoverage)
{
	VectorWithDirection<uint64_t> result;
	result.resize(numKmerNodes, std::numeric_limits<uint64_t>::max());
	for (size_t i = 0; i < numKmerNodes; i++)
	{
		if (!keptNodes.get(i)) continue;
		assert(i < numKmerNodes);
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : edgeCoverage.getValues(fw))
		{
			if (pair.second < minCoverage) continue;
			if (!keptNodes.get(pair.first.first)) continue;
			if (result[std::make_pair(i, true)] == std::numeric_limits<uint64_t>::max())
			{
				result[std::make_pair(i, true)] = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			}
			else
			{
				result[std::make_pair(i, true)] = std::numeric_limits<uint64_t>::max()-1;
			}
			if (result[reverse(pair.first)] == std::numeric_limits<uint64_t>::max())
			{
				result[reverse(pair.first)] = i;
			}
			else
			{
				result[reverse(pair.first)] = std::numeric_limits<uint64_t>::max()-1;
			}
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : edgeCoverage.getValues(bw))
		{
			if (pair.second < minCoverage) continue;
			if (!keptNodes.get(pair.first.first)) continue;
			if (result[std::make_pair(i, false)] == std::numeric_limits<uint64_t>::max())
			{
				result[std::make_pair(i, false)] = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			}
			else
			{
				result[std::make_pair(i, false)] = std::numeric_limits<uint64_t>::max()-1;
			}
			if (result[reverse(pair.first)] == std::numeric_limits<uint64_t>::max())
			{
				result[reverse(pair.first)] = i + firstBitUint64_t;
			}
			else
			{
				result[reverse(pair.first)] = std::numeric_limits<uint64_t>::max()-1;
			}
		}
	}
	return result;
}

VectorWithDirection<uint64_t> getUniqueEdges(const MostlySparse2DHashmap<uint8_t, size_t>& edgeCoverage, const size_t numKmerNodes, const size_t minCoverage)
{
	VectorWithDirection<uint64_t> result;
	result.resize(numKmerNodes, std::numeric_limits<uint64_t>::max());
	for (size_t i = 0; i < numKmerNodes; i++)
	{
		assert(i < numKmerNodes);
		std::pair<size_t, bool> fw { i, true };
		for (auto pair : edgeCoverage.getValues(fw))
		{
			if (pair.second < minCoverage) continue;
			if (result[std::make_pair(i, true)] == std::numeric_limits<uint64_t>::max())
			{
				result[std::make_pair(i, true)] = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			}
			else
			{
				result[std::make_pair(i, true)] = std::numeric_limits<uint64_t>::max()-1;
			}
			if (result[reverse(pair.first)] == std::numeric_limits<uint64_t>::max())
			{
				result[reverse(pair.first)] = i;
			}
			else
			{
				result[reverse(pair.first)] = std::numeric_limits<uint64_t>::max()-1;
			}
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : edgeCoverage.getValues(bw))
		{
			if (pair.second < minCoverage) continue;
			if (result[std::make_pair(i, false)] == std::numeric_limits<uint64_t>::max())
			{
				result[std::make_pair(i, false)] = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			}
			else
			{
				result[std::make_pair(i, false)] = std::numeric_limits<uint64_t>::max()-1;
			}
			if (result[reverse(pair.first)] == std::numeric_limits<uint64_t>::max())
			{
				result[reverse(pair.first)] = i + firstBitUint64_t;
			}
			else
			{
				result[reverse(pair.first)] = std::numeric_limits<uint64_t>::max()-1;
			}
		}
	}
	return result;
}

VectorWithDirection<uint64_t> getUniqueEdges(const std::vector<ReadPathBundle>& kmerGraphReadPaths, const size_t numKmerNodes, const size_t minCoverage)
{
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverage;
	edgeCoverage.resize(numKmerNodes);
	for (size_t i = 0; i < kmerGraphReadPaths.size(); i++)
	{
		for (size_t j = 0; j < kmerGraphReadPaths[i].paths.size(); j++)
		{
			for (size_t k = 1; k < kmerGraphReadPaths[i].paths[j].path.size(); k++)
			{
				std::pair<size_t, bool> fromNode { kmerGraphReadPaths[i].paths[j].path[k-1] & maskUint64_t, kmerGraphReadPaths[i].paths[j].path[k-1] & firstBitUint64_t };
				std::pair<size_t, bool> toNode { kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t, kmerGraphReadPaths[i].paths[j].path[k] & firstBitUint64_t };
				auto key = canon(fromNode, toNode);
				assert(key.first.first < numKmerNodes);
				assert(key.second.first < numKmerNodes);
				size_t coverage = 0;
				if (edgeCoverage.hasValue(key.first, key.second)) coverage = edgeCoverage.get(key.first, key.second);
				if (coverage < minCoverage) coverage += 1;
				edgeCoverage.set(key.first, key.second, coverage);
			}
		}
	}
	return getUniqueEdges(edgeCoverage, numKmerNodes, minCoverage);
}

VectorWithDirection<uint64_t> getUniqueEdges(const std::vector<ReadPathBundle>& kmerGraphReadPaths, const size_t numKmerNodes, const RankBitvector& keptNodes)
{
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverage;
	edgeCoverage.resize(numKmerNodes);
	for (size_t i = 0; i < kmerGraphReadPaths.size(); i++)
	{
		for (size_t j = 0; j < kmerGraphReadPaths[i].paths.size(); j++)
		{
			for (size_t k = 1; k < kmerGraphReadPaths[i].paths[j].path.size(); k++)
			{
				std::pair<size_t, bool> fromNode { kmerGraphReadPaths[i].paths[j].path[k-1] & maskUint64_t, kmerGraphReadPaths[i].paths[j].path[k-1] & firstBitUint64_t };
				std::pair<size_t, bool> toNode { kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t, kmerGraphReadPaths[i].paths[j].path[k] & firstBitUint64_t };
				if (!keptNodes.get(fromNode.first)) continue;
				if (!keptNodes.get(toNode.first)) continue;
				auto key = canon(fromNode, toNode);
				edgeCoverage.set(key.first, key.second, 1);
			}
		}
	}
	return getUniqueEdges(edgeCoverage, numKmerNodes, 1);
}

std::vector<std::pair<uint64_t, size_t>> getKmerPosToUnitigPos(const std::vector<std::vector<uint64_t>>& unitigs, const size_t kmerNodeCount)
{
	std::vector<std::pair<uint64_t, size_t>> result;
	result.resize(kmerNodeCount, std::make_pair(std::numeric_limits<uint64_t>::max(), std::numeric_limits<size_t>::max()));
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		for (size_t j = 0; j < unitigs[i].size(); j++)
		{
			size_t node = unitigs[i][j] & maskUint64_t;
			assert(node < result.size());
			assert(result[node].first == std::numeric_limits<uint64_t>::max());
			result[node].first = i;
			result[node].second = j;
			if (unitigs[i][j] & firstBitUint64_t)
			{
				result[node].first += firstBitUint64_t;
			}
			else
			{
				result[node].second = unitigs[i].size()-1-j;
			}
		}
	}
	return result;
}

std::vector<ReadPathBundle> getReadPaths(const std::vector<size_t>& kmerGraphNodeLengths, const std::vector<ReadPathBundle>& kmerGraphReadPaths, const std::vector<std::pair<uint64_t, size_t>>& kmerNodeToUnitig, const std::vector<std::vector<uint64_t>>& unitigs, const UnitigGraph& unitigGraph)
{
	std::vector<ReadPathBundle> result;
	result.resize(kmerGraphReadPaths.size());
	std::vector<std::vector<size_t>> unitigNodeStartPosKmers;
	unitigNodeStartPosKmers.resize(unitigs.size());
	for (size_t i = 0; i < unitigs.size(); i++)
	{
		unitigNodeStartPosKmers[i].resize(unitigs[i].size()+1);
		size_t startPos = 0;
		for (size_t j = 0; j < unitigs[i].size(); j++)
		{
			unitigNodeStartPosKmers[i][j] = startPos;
			startPos += kmerGraphNodeLengths[unitigs[i][j] & maskUint64_t];
		}
		unitigNodeStartPosKmers[i].back() = startPos;
		assert(startPos == unitigGraph.lengths[i]);
	}
	for (size_t i = 0; i < kmerGraphReadPaths.size(); i++)
	{
		result[i].readName = kmerGraphReadPaths[i].readName;
		result[i].readLength = kmerGraphReadPaths[i].readLength;
		for (size_t j = 0; j < kmerGraphReadPaths[i].paths.size(); j++)
		{
			std::pair<uint64_t, size_t> lastPos { std::numeric_limits<uint64_t>::max(), std::numeric_limits<size_t>::max() };
			size_t readPos = kmerGraphReadPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < kmerGraphReadPaths[i].paths[j].path.size(); k++)
			{
				std::pair<uint64_t, size_t> thisPos = kmerNodeToUnitig[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
				if (thisPos.first == std::numeric_limits<uint64_t>::max())
				{
					if (lastPos.first != std::numeric_limits<uint64_t>::max())
					{
						assert(result[i].paths.back().path.size() >= 1);
						if (lastPos.first & firstBitUint64_t)
						{
							result[i].paths.back().pathRightClipKmers = unitigGraph.lengths[lastPos.first & maskUint64_t] - unitigNodeStartPosKmers[lastPos.first & maskUint64_t][lastPos.second+1];
						}
						else
						{
							result[i].paths.back().pathRightClipKmers = unitigNodeStartPosKmers[lastPos.first & maskUint64_t][unitigs[lastPos.first & maskUint64_t].size() - lastPos.second - 1];
						}
					}
					readPos += kmerGraphNodeLengths[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
					if (k == 0) readPos -= kmerGraphReadPaths[i].paths[j].pathLeftClipKmers;
					lastPos = thisPos;
					continue;
				}
				if ((kmerGraphReadPaths[i].paths[j].path[k] & firstBitUint64_t) == 0)
				{
					thisPos.first ^= firstBitUint64_t;
					thisPos.second = unitigs[thisPos.first & maskUint64_t].size() - 1 - thisPos.second;
				}
				if (lastPos.first != std::numeric_limits<uint64_t>::max() && thisPos.first == lastPos.first && thisPos.second == lastPos.second+1)
				{
					readPos += kmerGraphNodeLengths[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
					if (k == 0) readPos -= kmerGraphReadPaths[i].paths[j].pathLeftClipKmers;
					lastPos = thisPos;
					continue;
				}
				if (lastPos.first != std::numeric_limits<uint64_t>::max() && thisPos.second == 0 && lastPos.second == unitigs[lastPos.first & maskUint64_t].size()-1)
				{
					auto from = std::make_pair(lastPos.first & maskUint64_t, lastPos.first & firstBitUint64_t);
					auto to = std::make_pair(thisPos.first & maskUint64_t, thisPos.first & firstBitUint64_t);
					auto key = canon(from, to);
					if (unitigGraph.edgeCoverages.hasValue(key.first, key.second))
					{
						result[i].paths.back().path.emplace_back(thisPos.first);
						readPos += kmerGraphNodeLengths[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
						if (k == 0) readPos -= kmerGraphReadPaths[i].paths[j].pathLeftClipKmers;
						lastPos = thisPos;
						continue;
					}
				}
				if (lastPos.first != std::numeric_limits<uint64_t>::max())
				{
					assert(result[i].paths.back().path.size() >= 1);
					if (lastPos.first & firstBitUint64_t)
					{
						result[i].paths.back().pathRightClipKmers = unitigGraph.lengths[lastPos.first & maskUint64_t] - unitigNodeStartPosKmers[lastPos.first & maskUint64_t][lastPos.second+1];
					}
					else
					{
						result[i].paths.back().pathRightClipKmers = unitigNodeStartPosKmers[lastPos.first & maskUint64_t][unitigs[lastPos.first & maskUint64_t].size() - lastPos.second - 1];
					}
				}
				result[i].paths.emplace_back();
				result[i].paths.back().readStartPos = readPos;
				if (thisPos.first & firstBitUint64_t)
				{
					result[i].paths.back().pathLeftClipKmers = unitigNodeStartPosKmers[thisPos.first & maskUint64_t][thisPos.second];
				}
				else
				{
					assert(unitigNodeStartPosKmers[thisPos.first & maskUint64_t][unitigs[thisPos.first].size()-1-thisPos.second] <= unitigGraph.lengths[thisPos.first & maskUint64_t] - 1);
					result[i].paths.back().pathLeftClipKmers = unitigGraph.lengths[thisPos.first & maskUint64_t] - unitigNodeStartPosKmers[thisPos.first & maskUint64_t][unitigs[thisPos.first].size()-thisPos.second];
				}
				if (k == 0) result[i].paths.back().pathLeftClipKmers += kmerGraphReadPaths[i].paths[j].pathLeftClipKmers;
				result[i].paths.back().path.emplace_back(thisPos.first);
				lastPos = thisPos;
				readPos += kmerGraphNodeLengths[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
				if (k == 0) readPos -= kmerGraphReadPaths[i].paths[j].pathLeftClipKmers;
			}
			if (lastPos.first != std::numeric_limits<uint64_t>::max())
			{
				assert(result[i].paths.back().path.size() >= 1);
				if (lastPos.first & firstBitUint64_t)
				{
					result[i].paths.back().pathRightClipKmers = unitigGraph.lengths[lastPos.first & maskUint64_t] - unitigNodeStartPosKmers[lastPos.first & maskUint64_t][lastPos.second+1];
				}
				else
				{
					result[i].paths.back().pathRightClipKmers = unitigNodeStartPosKmers[lastPos.first & maskUint64_t][unitigs[lastPos.first & maskUint64_t].size() - lastPos.second - 1];
				}
				result[i].paths.back().pathRightClipKmers += kmerGraphReadPaths[i].paths[j].pathRightClipKmers;
			}
		}
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> makeUnitigGraph(const KmerGraph& kmerGraph, const std::vector<ReadPathBundle>& kmerGraphReadPaths, const size_t minCoverage)
{
	UnitigGraph result;
	std::vector<std::vector<uint64_t>> unitigs = getUnitigs(kmerGraph.nodeCount(), getUniqueEdges(kmerGraphReadPaths, kmerGraph.nodeCount(), minCoverage));
	std::vector<std::pair<uint64_t, size_t>> kmerNodeToUnitig = getKmerPosToUnitigPos(unitigs, kmerGraph.nodeCount());
	for (size_t i = 0; i < kmerNodeToUnitig.size(); i++)
	{
		assert(kmerNodeToUnitig[i].first != std::numeric_limits<uint64_t>::max());
		if (kmerNodeToUnitig[i].first & firstBitUint64_t)
		{
			assert(unitigs[kmerNodeToUnitig[i].first & maskUint64_t][kmerNodeToUnitig[i].second] == (i + firstBitUint64_t));
		}
		else
		{
			assert(unitigs[kmerNodeToUnitig[i].first & maskUint64_t][unitigs[kmerNodeToUnitig[i].first & maskUint64_t].size() - 1 - kmerNodeToUnitig[i].second] == i);
		}
	}
	std::tie(result.lengths, result.coverages) = getUnitigLengthAndCoverage(unitigs, kmerGraphReadPaths, kmerNodeToUnitig, kmerGraph.lengths);
	result.edgeCoverages = getUnitigEdgeCoverages(unitigs, minCoverage, kmerGraphReadPaths, kmerNodeToUnitig);
	auto readPaths = getReadPaths(kmerGraph.lengths, kmerGraphReadPaths, kmerNodeToUnitig, unitigs, result);
	return std::make_pair(std::move(result), std::move(readPaths));
}

size_t UnitigGraph::nodeCount() const
{
	return coverages.size();
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> filterUnitigGraph(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const RankBitvector& keptNodes)
{
	assert(keptNodes.size() == unitigGraph.nodeCount());
	UnitigGraph result;
	std::vector<ReadPathBundle> resultPaths;
	std::vector<std::vector<uint64_t>> unitigs = getUnitigs(unitigGraph.nodeCount(), keptNodes, getUniqueEdges(readPaths, unitigGraph.nodeCount(), keptNodes));
	std::vector<std::pair<uint64_t, size_t>> kmerNodeToUnitig = getKmerPosToUnitigPos(unitigs, unitigGraph.nodeCount());
	std::tie(result.lengths, result.coverages) = getUnitigLengthAndCoverage(unitigs, readPaths, kmerNodeToUnitig, unitigGraph.lengths);
	result.edgeCoverages = getUnitigEdgeCoverages(unitigs, 1, readPaths, kmerNodeToUnitig);
	auto resultReadPaths = getReadPaths(unitigGraph.lengths, readPaths, kmerNodeToUnitig, unitigs, result);
	return std::make_pair(std::move(result), std::move(resultReadPaths));
}

void addSequence(std::vector<TwobitString>& sequences, std::vector<std::vector<bool>>& bpSequenceGotten, const uint64_t node, const size_t nodeKmerStartPos, const size_t nodeKmerEndPos, const size_t readI, const size_t readStartPos, const std::vector<TwobitString>& readSequences, const size_t k)
{
	for (size_t i = 0; i < nodeKmerEndPos - nodeKmerStartPos + k - 1; i++)
	{
		uint8_t c = readSequences[readI].get(readStartPos + i);
		size_t nodePos = nodeKmerStartPos + i;
		if ((node & firstBitUint64_t) == 0)
		{
			c = 3-c;
			nodePos = sequences[node & maskUint64_t].size() - 1 - nodePos;
		}
		if (!bpSequenceGotten[node & maskUint64_t][nodePos])
		{
			sequences[node & maskUint64_t].set(nodePos, c);
			bpSequenceGotten[node & maskUint64_t][nodePos] = true;
		}
		else
		{
			assert(sequences[node & maskUint64_t].get(nodePos) == c);
		}
	}
}

std::vector<TwobitString> getNodeSequences(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t kmerSize, const std::vector<TwobitString>& readSequences)
{
	std::vector<std::vector<bool>> bpSequenceGotten;
	std::vector<TwobitString> sequences;
	bpSequenceGotten.resize(unitigGraph.nodeCount());
	sequences.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < unitigGraph.nodeCount();i++)
	{
		bpSequenceGotten[i].resize(unitigGraph.lengths[i] + kmerSize - 1, false);
		sequences[i].resize(unitigGraph.lengths[i] + kmerSize - 1);
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readStartPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				size_t nodeKmerStartPos = 0;
				size_t nodeKmerEndPos = unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) nodeKmerStartPos = readPaths[i].paths[j].pathLeftClipKmers;
				if (k == readPaths[i].paths[j].path.size()-1) nodeKmerEndPos -= readPaths[i].paths[j].pathRightClipKmers;
				addSequence(sequences, bpSequenceGotten, node, nodeKmerStartPos, nodeKmerEndPos, i, readStartPos, readSequences, kmerSize);
				readStartPos += nodeKmerEndPos - nodeKmerStartPos;
			}
		}
	}
	for (size_t i = 0; i < bpSequenceGotten.size(); i++)
	{
		for (size_t j = 0; j < bpSequenceGotten[i].size(); j++)
		{
			// assert(bpSequenceGotten[i][j]);
		}
	}
	return sequences;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unitigify(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths)
{
	UnitigGraph result;
	std::vector<ReadPathBundle> resultPaths;
	std::vector<std::vector<uint64_t>> unitigs = getUnitigs(unitigGraph.nodeCount(), getUniqueEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount(), 1));
	std::vector<std::pair<uint64_t, size_t>> kmerNodeToUnitig = getKmerPosToUnitigPos(unitigs, unitigGraph.nodeCount());
	std::tie(result.lengths, result.coverages) = getUnitigLengthAndCoverage(unitigs, readPaths, kmerNodeToUnitig, unitigGraph.lengths);
	result.edgeCoverages = getUnitigEdgeCoverages(unitigs, 1, readPaths, kmerNodeToUnitig);
	auto resultReadPaths = getReadPaths(unitigGraph.lengths, readPaths, kmerNodeToUnitig, unitigs, result);
	return std::make_pair(std::move(result), std::move(resultReadPaths));
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unitigifyWithFilter(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const RankBitvector& keptNodes)
{
	UnitigGraph result;
	std::vector<ReadPathBundle> resultPaths;
	std::vector<std::vector<uint64_t>> unitigs = getUnitigs(unitigGraph.nodeCount(), keptNodes, getUniqueEdges(unitigGraph.edgeCoverages, keptNodes, unitigGraph.nodeCount(), 1));
	std::vector<std::pair<uint64_t, size_t>> kmerNodeToUnitig = getKmerPosToUnitigPos(unitigs, unitigGraph.nodeCount());
	std::tie(result.lengths, result.coverages) = getUnitigLengthAndCoverage(unitigs, readPaths, kmerNodeToUnitig, unitigGraph.lengths);
	result.edgeCoverages = getUnitigEdgeCoverages(unitigs, 1, readPaths, kmerNodeToUnitig);
	auto resultReadPaths = getReadPaths(unitigGraph.lengths, readPaths, kmerNodeToUnitig, unitigs, result);
	return std::make_pair(std::move(result), std::move(resultReadPaths));
}
