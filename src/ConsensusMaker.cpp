#include "ConsensusMaker.h"
#include "Common.h"
#include "TwobitString.h"
#include "EdlibWrapper.h"
#include "KmerIterator.h"

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

std::string getConsensusSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const phmap::flat_hash_map<size_t, std::vector<std::tuple<size_t, bool, size_t>>>& readAnchorPoints, const size_t start, const size_t end, const size_t kmerSize)
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
				sequenceHere = sequenceIndex.getSubstring(std::get<0>(t1), std::get<2>(t1), std::get<2>(t2)-std::get<2>(t1)+kmerSize);
			}
			else
			{
				if (std::get<2>(t1) <= std::get<2>(t2)) continue;
				if (std::get<2>(t1)-std::get<2>(t2) > end - start + 50) continue;
				if (std::get<2>(t1)-std::get<2>(t2) + 50 < end - start) continue;
				sequenceHere = sequenceIndex.getSubstring(std::get<0>(t1), std::get<2>(t2), std::get<2>(t1)-std::get<2>(t2)+kmerSize);
				sequenceHere = revCompRaw(sequenceHere);
			}
			sequenceCount[sequenceHere] += 1;
			totalCount += 1;
		}
	}
	if (sequenceCount.size() == 0)
	{
		return "";
	}
	{
		const std::string& exampleRead = sequenceCount.begin()->first;
		assert(exampleRead.size() >= kmerSize);
		std::string firstKmer = exampleRead.substr(0, kmerSize);
		std::string lastKmer = exampleRead.substr(exampleRead.size() - kmerSize, kmerSize);
		for (const auto& pair : sequenceCount)
		{
			assert(pair.first.substr(0, kmerSize) == firstKmer);
			assert(pair.first.substr(pair.first.size()-kmerSize, kmerSize) == lastKmer);
		}
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
