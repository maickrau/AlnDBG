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

std::vector<std::vector<uint64_t>> getUnitigs(const size_t countNodes, const size_t minCoverage, const VectorWithDirection<uint64_t>& uniqueEdges)
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
				size_t unitig = kmerNodeToUnitig[node].first & maskUint64_t;
				assert(unitig < coverage.size());
				assert(node < nodeLength.size());
				coverage[unitig] += nodeLength[node];
			}
		}
	}
	for (size_t i = 0; i < unitigs.size(); i++)
	{
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
				if ((kmerGraphReadPaths[i].paths[j].path[k-1] & firstBitUint64_t) == 0)
				{
					beforePos.first ^= firstBitUint64_t;
					beforePos.second = unitigs[beforePos.first & maskUint64_t].size() - 1 - beforePos.second;
				}
				std::pair<uint64_t, size_t> afterPos = kmerNodeToUnitig[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
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
				size_t coverage = 0;
				if (edgeCoverage.hasValue(key.first, key.second)) coverage = edgeCoverage.get(key.first, key.second);
				if (coverage < minCoverage) coverage += 1;
				edgeCoverage.set(key.first, key.second, coverage);
			}
		}
	}
	VectorWithDirection<uint64_t> result;
	result.resize(numKmerNodes, std::numeric_limits<uint64_t>::max());
	for (size_t i = 0; i < numKmerNodes; i++)
	{
		for (auto pair : edgeCoverage.getValues(std::make_pair(i, true)))
		{
			if (pair.second < minCoverage) continue;
			if (result[std::make_pair(i, true)] == std::numeric_limits<size_t>::max())
			{
				result[std::make_pair(i, true)] = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			}
			else
			{
				result[std::make_pair(i, true)] = std::numeric_limits<size_t>::max()-1;
			}
			if (result[reverse(pair.first)] == std::numeric_limits<size_t>::max())
			{
				result[reverse(pair.first)] = i;
			}
			else
			{
				result[reverse(pair.first)] = std::numeric_limits<size_t>::max()-1;
			}
		}
		for (auto pair : edgeCoverage.getValues(std::make_pair(i, false)))
		{
			if (pair.second < minCoverage) continue;
			if (result[std::make_pair(i, false)] == std::numeric_limits<size_t>::max())
			{
				result[std::make_pair(i, false)] = pair.first.first + (pair.first.second ? firstBitUint64_t : 0);
			}
			else
			{
				result[std::make_pair(i, false)] = std::numeric_limits<size_t>::max()-1;
			}
			if (result[reverse(pair.first)] == std::numeric_limits<size_t>::max())
			{
				result[reverse(pair.first)] = i + firstBitUint64_t;
			}
			else
			{
				result[reverse(pair.first)] = std::numeric_limits<size_t>::max()-1;
			}
		}
	}
	return result;
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
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i].first != std::numeric_limits<uint64_t>::max());
		if (result[i].first & firstBitUint64_t)
		{
			assert(unitigs[result[i].first & maskUint64_t][result[i].second] == (i + firstBitUint64_t));
		}
		else
		{
			assert(unitigs[result[i].first & maskUint64_t][unitigs[result[i].first & maskUint64_t].size() - 1 - result[i].second] == i);
		}
	}
	return result;
}

std::vector<ReadPathBundle> getReadPaths(const KmerGraph& kmerGraph, const std::vector<ReadPathBundle>& kmerGraphReadPaths, const std::vector<std::pair<uint64_t, size_t>>& kmerNodeToUnitig, const std::vector<std::vector<uint64_t>>& unitigs, const UnitigGraph& unitigGraph)
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
			startPos += kmerGraph.lengths[unitigs[i][j] & maskUint64_t];
		}
		unitigNodeStartPosKmers[i].back() = startPos;
	}
	for (size_t i = 0; i < kmerGraphReadPaths.size(); i++)
	{
		for (size_t j = 0; j < kmerGraphReadPaths[i].paths.size(); j++)
		{
			std::pair<uint64_t, size_t> lastPos { std::numeric_limits<uint64_t>::max(), std::numeric_limits<size_t>::max() };
			size_t readPos = kmerGraphReadPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < kmerGraphReadPaths[i].paths[j].path.size(); k++)
			{
				std::pair<uint64_t, size_t> thisPos = kmerNodeToUnitig[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
				if (!(kmerGraphReadPaths[i].paths[j].path[k] & firstBitUint64_t))
				{
					thisPos.first ^= firstBitUint64_t;
					thisPos.second = unitigs[thisPos.first & maskUint64_t].size() - 1 - thisPos.second;
				}
				if (lastPos.first != std::numeric_limits<uint64_t>::max() && thisPos.first == lastPos.first && thisPos.second == lastPos.second+1)
				{
					readPos += kmerGraph.lengths[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
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
						readPos += kmerGraph.lengths[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
						lastPos = thisPos;
						continue;
					}
				}
				if (lastPos.first != std::numeric_limits<uint64_t>::max())
				{
					assert(lastPos.first != std::numeric_limits<uint64_t>::max());
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
				result[i].paths.back().path.emplace_back(thisPos.first);
				lastPos = thisPos;
				readPos += kmerGraph.lengths[kmerGraphReadPaths[i].paths[j].path[k] & maskUint64_t];
			}
			if (lastPos.first != std::numeric_limits<uint64_t>::max())
			{
				assert(lastPos.first != std::numeric_limits<uint64_t>::max());
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
		}
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> makeUnitigGraph(const KmerGraph& kmerGraph, const std::vector<ReadPathBundle>& kmerGraphReadPaths, const size_t minCoverage)
{
	UnitigGraph result;
	std::vector<std::vector<uint64_t>> unitigs = getUnitigs(kmerGraph.nodeCount(), minCoverage, getUniqueEdges(kmerGraphReadPaths, kmerGraph.nodeCount(), minCoverage));
	std::vector<std::pair<uint64_t, size_t>> kmerNodeToUnitig = getKmerPosToUnitigPos(unitigs, kmerGraph.nodeCount());
	std::tie(result.lengths, result.coverages) = getUnitigLengthAndCoverage(unitigs, kmerGraphReadPaths, kmerNodeToUnitig, kmerGraph.lengths);
	result.edgeCoverages = getUnitigEdgeCoverages(unitigs, minCoverage, kmerGraphReadPaths, kmerNodeToUnitig);
	auto readPaths = getReadPaths(kmerGraph, kmerGraphReadPaths, kmerNodeToUnitig, unitigs, result);
	return std::make_pair(std::move(result), std::move(readPaths));
}

size_t UnitigGraph::nodeCount() const
{
	return coverages.size();
}
