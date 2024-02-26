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

void writeGraph(const std::string& filename, const std::vector<std::vector<size_t>>& lengths, const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeCoverage, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps)
{
	std::ofstream graph { filename };
	for (size_t i = 0; i < coverages.size(); i++)
	{
		assert(lengths[i].size() >= coverages[i]);
		if (coverages[i] < 2) continue;
		size_t length = lengths[i][lengths[i].size()/2];
		size_t coverage = coverages[i];
		graph << "S\t" << i << "\t*\tLN:i:" << length << "\tll:i:" << coverage << "\tFC:i:" << length*coverage << std::endl;
	}
	for (auto pair : edgeCoverage)
	{
		if (pair.second < 2) continue;
		assert(edgeOverlaps.count(pair.first) == 1);
		graph << "L\t" << (pair.first.first & maskUint64_t) << "\t" << ((pair.first.first & firstBitUint64_t) ? "+" : "-") << "\t" << (pair.first.second & maskUint64_t) << "\t" << ((pair.first.second & firstBitUint64_t) ? "+" : "-") << "\t" << edgeOverlaps.at(pair.first) << "M\tec:i:" << pair.second << std::endl;
	}
}

std::pair<std::vector<std::vector<size_t>>, std::vector<size_t>> getLengthsAndCoverages(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
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
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
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
	std::cerr << "splitting per length" << std::endl;
	const double differenceFraction = 0.02;
	const size_t differenceConstant = 50;
	phmap::flat_hash_map<size_t, std::vector<size_t>> lengthsPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
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

size_t getNumMismatches(const TwobitString& leftSequence, const TwobitString& rightSequence, const size_t maxMismatches)
{
	assert(leftSequence.size() >= rightSequence.size());
	assert(rightSequence.size() + maxMismatches >= leftSequence.size());
	// left sequence: left to right
	// right sequence: up to down
	// offset: diagonally to down-right
	// diagonal: up to down
	// zero diagonal (top left corner) at maxBefore
	std::vector<size_t> offset;
	size_t maxBefore = (leftSequence.size() - rightSequence.size()) + (maxMismatches - (leftSequence.size() - rightSequence.size()))/2;
	size_t maxAfter = (maxMismatches - (leftSequence.size() - rightSequence.size()))/2;
	assert(maxBefore < maxMismatches);
	assert(maxAfter < maxMismatches);
	offset.resize(maxBefore + maxAfter + 1, std::numeric_limits<size_t>::max());
	offset[maxBefore] = 0;
	std::vector<size_t> prevOffset;
	prevOffset.resize(offset.size());
	while (offset[maxBefore] < leftSequence.size() && offset[maxBefore] < rightSequence.size() && leftSequence.get(offset[maxBefore]) == rightSequence.get(offset[maxBefore])) offset[maxBefore] += 1;
	if (leftSequence.size() == rightSequence.size() && offset[maxBefore] == leftSequence.size()) return 0;
	for (size_t score = 1; score < maxMismatches; score++)
	{
		std::swap(offset, prevOffset);
		if (score <= maxBefore) prevOffset[maxBefore-score] = score;
		if (maxBefore+score < offset.size()) prevOffset[maxBefore+score] = score;
		for (size_t diagonal = (score < maxBefore ? (maxBefore - score) : 0); diagonal <= maxBefore + score && diagonal < offset.size(); diagonal++)
		{
			assert(prevOffset[diagonal] != std::numeric_limits<size_t>::max());
			offset[diagonal] = prevOffset[diagonal]+1;
			if (diagonal > 0 && prevOffset[diagonal-1] != std::numeric_limits<size_t>::max())
			{
				offset[diagonal] = std::max(offset[diagonal], prevOffset[diagonal-1]);
			}
			if (diagonal+1 < offset.size() && prevOffset[diagonal+1] != std::numeric_limits<size_t>::max())
			{
				offset[diagonal] = std::max(offset[diagonal], prevOffset[diagonal+1]+1);
			}
			offset[diagonal] = std::min(offset[diagonal], std::min(leftSequence.size(), rightSequence.size()+maxBefore-diagonal));
		}
		for (size_t diagonal = (score < maxBefore ? (maxBefore - score) : 0); diagonal <= maxBefore + score && diagonal < offset.size(); diagonal++)
		{
			assert(offset[diagonal] <= leftSequence.size());
			assert(offset[diagonal]+diagonal-maxBefore <= rightSequence.size());
			while (offset[diagonal] < leftSequence.size() && offset[diagonal]+diagonal-maxBefore < rightSequence.size() && leftSequence.get(offset[diagonal]) == rightSequence.get(offset[diagonal]+diagonal-maxBefore)) offset[diagonal] += 1;
			if (offset[diagonal] == leftSequence.size() && offset[diagonal]+diagonal-maxBefore == rightSequence.size()) return score;
		}
	}
	return maxMismatches+1;
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

bool checkPhasablePair(const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& bwForks, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& fwForks, const size_t bwHom, const size_t fwHom)
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
			bw = std::get<2>(t);
			bwPos = std::get<0>(t);
		}
		if (bw == std::numeric_limits<size_t>::max()) return false;
		for (auto t : fwForks[i])
		{
			if (std::get<1>(t) != fwHom) continue;
			if (std::get<0>(t) > bwPos) return false;
			fw = std::get<2>(t);
		}
		if (fw == std::numeric_limits<size_t>::max()) return false;
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
	if (clusterCount.size() < 2) return false;
	for (auto pair : clusterCount)
	{
		if (pair.second < 4) return false;
	}
	return true;
}

void splitPerPhasingKmersWithinChunk(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting per phasing kmers" << std::endl;
	const size_t k = 11;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
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
		std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> fwForks; // hom pos, hom kmer, het kmer
		std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> bwForks; // hom pos, hom kmer, het kmer
		fwForks.resize(occurrencesPerChunk[i].size());
		bwForks.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			size_t lastKmer = std::numeric_limits<size_t>::max();
			iterateKmers(readSequences[occurrencesPerChunk[i][j].first], std::get<0>(t), std::get<1>(t), std::get<2>(t) & maskUint64_t, k, [i, j, &fwForks, &bwForks, &lastKmer, &kmersEverywhere](const size_t kmer, const size_t pos)
			{
				if (kmersEverywhere.count(lastKmer) == 1 && kmersEverywhere.count(kmer) == 0)
				{
					fwForks[j].emplace_back(pos-1, lastKmer, kmer);
				}
				else if (lastKmer != std::numeric_limits<size_t>::max() && kmersEverywhere.count(lastKmer) == 0 && kmersEverywhere.count(kmer) == 1)
				{
					bwForks[j].emplace_back(pos, kmer, lastKmer);
				}
				lastKmer = kmer;
			});
		}
		phmap::flat_hash_set<std::pair<size_t, size_t>> validPairs;
		phmap::flat_hash_set<std::pair<size_t, size_t>> invalidPairs;
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
					if (invalidPairs.count(key) == 1) continue;
					if (validPairs.count(key) == 0)
					{
						bool phasable = checkPhasablePair(bwForks, fwForks, key.first, key.second);
						if (!phasable)
						{
							invalidPairs.insert(key);
							continue;
						}
						else
						{
							validPairs.insert(key);
						}
					}
					if (validPairs.count(key) == 1)
					{
						phaseIdentities[j].emplace_back(std::get<2>(bwForks[j][n]));
						phaseIdentities[j].emplace_back(std::get<2>(fwForks[j][m]));
					}
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
			std::cerr << "phasable kmer splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
			}
		}
	});
}

void splitPerSequenceIdentity(const std::vector<TwobitString>& readSequences, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	std::cerr << "splitting per sequence identity" << std::endl;
	const double mismatchFraction = 0.05;
	const size_t mismatchFloor = 10;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
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
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, &iterationOrder, mismatchFraction, mismatchFloor](const size_t iterationIndex)
	{
		const size_t i = iterationOrder[iterationIndex];
		{
			std::lock_guard<std::mutex> lock { resultMutex };
			std::cerr << "split chunk " << i << " coverage " << occurrencesPerChunk[i].size() << std::endl;
		}
		std::vector<std::pair<TwobitString, size_t>> sequencesPerOccurrence;
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			sequencesPerOccurrence.emplace_back();
			sequencesPerOccurrence.back().second = j;
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
			if (std::get<2>(t) & firstBitUint64_t)
			{
				for (size_t k = std::get<0>(t); k <= std::get<1>(t); k++)
				{
					sequencesPerOccurrence.back().first.emplace_back(readSequences[occurrencesPerChunk[i][j].first].get(k));
				}
			}
			else
			{
				for (size_t k = std::get<1>(t); k >= std::get<0>(t) && k < readSequences[occurrencesPerChunk[i][j].first].size(); k--)
				{
					sequencesPerOccurrence.back().first.emplace_back(3-readSequences[occurrencesPerChunk[i][j].first].get(k));
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
			std::cerr << "splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
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
	const double mismatchFraction = 0.05;
	const size_t mismatchFloor = 10;
	std::cerr << "splitting by base counts" << std::endl;
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t) + 1);
			occurrencesPerChunk[std::get<2>(t)].emplace_back(i, j);
		}
	}
	size_t nextNum = 0;
	std::mutex resultMutex;
	iterateMultithreaded(0, occurrencesPerChunk.size(), numThreads, [&nextNum, &resultMutex, &chunksPerRead, &occurrencesPerChunk, &readSequences, mismatchFraction, mismatchFloor](const size_t i)
	{
		std::vector<std::vector<size_t>> countsPerOccurrence;
		countsPerOccurrence.resize(occurrencesPerChunk[i].size());
		for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
		{
			auto t = chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second];
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
			std::cerr << "base count splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << keyToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t) + keyToNode.at(find(parent, j));
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
			std::cerr << "minhash splitted chunk " << i << " coverage " << occurrencesPerChunk[i].size() << " to " << clusterToNode.size() << " chunks" << std::endl;
			for (size_t j = 0; j < occurrencesPerChunk[i].size(); j++)
			{
				std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) = clusterToNode.at(find(parent, oneHashPerLocation[j])) + (std::get<2>(chunksPerRead[occurrencesPerChunk[i][j].first][occurrencesPerChunk[i][j].second]) & firstBitUint64_t);
			}
		}
	});
}

void removeHighCoverageChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t maxCoverage)
{
	phmap::flat_hash_map<size_t, size_t> coverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (auto t : chunksPerRead[i])
		{
			coverage[std::get<2>(t) & maskUint64_t] += 1;
		}
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = chunksPerRead[i].size()-1; j < chunksPerRead[i].size(); j--)
		{
			if (coverage.at(std::get<2>(chunksPerRead[i][j]) & maskUint64_t) < maxCoverage) continue;
			chunksPerRead[i].erase(chunksPerRead[i].begin()+j);
		}
	}
}

void makeGraph(const MatchIndex& matchIndex, const std::vector<size_t>& rawReadLengths, const std::vector<TwobitString>& readSequences, const size_t numThreads)
{
	std::vector<bool> useTheseChunks;
	useTheseChunks.resize(matchIndex.numWindowChunks() - matchIndex.numUniqueChunks(), true);
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead = getChunksPerRead(matchIndex, rawReadLengths, useTheseChunks);
	removeContainedChunks(chunksPerRead);
	splitPerLength(chunksPerRead);
	{
		std::cerr << "writing first graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round1.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
	splitPerBaseCounts(readSequences, chunksPerRead, numThreads);
	{
		std::cerr << "writing second graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round2.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths2.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
	splitPerMinHashes(readSequences, chunksPerRead, numThreads);
	{
		std::cerr << "writing third graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round3.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths3.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
	splitPerPhasingKmersWithinChunk(readSequences, chunksPerRead, numThreads);
	{
		std::cerr << "writing fourth graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round4.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths4.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
	splitPerSequenceIdentity(readSequences, chunksPerRead, numThreads);
	{
		std::cerr << "writing fifth graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round5.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths5.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
	splitPerPhasingKmersWithinChunk(readSequences, chunksPerRead, numThreads);
	{
		std::cerr << "writing sixth graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round6.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths6.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
	removeHighCoverageChunks(chunksPerRead, 60);
	{
		std::cerr << "writing seventh graph" << std::endl;
		std::vector<std::vector<size_t>> lengths;
		std::vector<size_t> coverages;
		std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
		writeGraph("graph-round7.gfa", lengths, coverages, edgeCoverage, edgeOverlaps);
		std::ofstream pathfile { "fakepaths7.txt" };
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
				auto t = chunksPerRead[i][j];
				uint64_t rawnode = std::get<2>(t);
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
			}
		}
	}
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
	for (size_t i = 0; i < readNames.size(); i++)
	{
		std::cerr << "readname " << i << " " << readNames[i] << std::endl;
	}
	makeGraph(matchIndex, readBasepairLengths, readSequences, numThreads);
}
