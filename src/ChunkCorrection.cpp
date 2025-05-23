#include <iostream>
#include <map>
#include <string>
#include <phmap.h>
#include "ChunkCorrection.h"
#include "Common.h"
#include "KmerIterator.h"

struct Overlap
{
	size_t referenceRead;
	size_t queryRead;
	std::vector<std::pair<size_t, size_t>> matches; // pair order ref, query
};

struct Correction
{
	size_t startIndex; // index of solid before replacement
	size_t endIndex; // index of solid after replacement
	std::vector<size_t> chunkReplacement; // does not include bounding solids
	std::vector<size_t> minimizerReplacement; // INCLUDES startIndex last minimizer and endIndex first minimizer
	bool operator<(const Correction& other) const
	{
		if (startIndex < other.startIndex) return true;
		if (startIndex > other.startIndex) return false;
		if (endIndex < other.endIndex) return true;
		if (endIndex > other.endIndex) return false;
		if (chunkReplacement < other.chunkReplacement) return true;
		if (chunkReplacement > other.chunkReplacement) return false;
		if (minimizerReplacement < other.minimizerReplacement) return true;
		if (minimizerReplacement > other.minimizerReplacement) return false;
		return false;
	}
};

std::vector<std::vector<size_t>> getReadsWhichHaveChunk(const std::vector<phmap::flat_hash_map<size_t, std::vector<size_t>>>& chunkIndicesWithinReads)
{
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < chunkIndicesWithinReads.size(); i++)
	{
		for (const auto& pair : chunkIndicesWithinReads[i])
		{
			while (result.size() <= pair.first) result.emplace_back();
			result[pair.first].emplace_back(i);
		}
	}
	return result;
}

std::vector<phmap::flat_hash_map<size_t, std::vector<size_t>>> getChunkIndicesWithinReads(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<phmap::flat_hash_map<size_t, std::vector<size_t>>> result;
	result.resize(chunksPerRead.size());
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			result[i][chunk].emplace_back(j);
		}
	}
	return result;
}

Overlap getUniqueBestOverlap(const std::vector<phmap::flat_hash_map<size_t, std::vector<size_t>>>& chunkIndicesWithinReads, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t referenceRead, const size_t queryRead)
{
	const size_t maxRepeatCoverage = 100;
	Overlap result;
	result.referenceRead = referenceRead;
	result.queryRead = queryRead;
	std::vector<std::pair<size_t, size_t>> matches;
	for (const auto& pair : chunkIndicesWithinReads[referenceRead])
	{
		if (pair.second.size() >= maxRepeatCoverage) continue;
		if (chunkIndicesWithinReads[queryRead].count(pair.first) == 0) continue;
		if (chunkIndicesWithinReads[queryRead].at(pair.first).size() >= maxRepeatCoverage) continue;
		for (const auto& pos2 : chunkIndicesWithinReads[queryRead].at(pair.first))
		{
			size_t queryChunkLength = std::get<1>(chunksPerRead[queryRead][pos2]) - std::get<0>(chunksPerRead[queryRead][pos2]);
			for (const size_t pos : pair.second)
			{
				assert(std::get<2>(chunksPerRead[referenceRead][pos]) == std::get<2>(chunksPerRead[queryRead][pos2]));
				size_t refChunkLength = std::get<1>(chunksPerRead[referenceRead][pos]) - std::get<0>(chunksPerRead[referenceRead][pos]);
				if (refChunkLength > queryChunkLength * 1.02 + 10) continue;
				if (queryChunkLength > refChunkLength * 1.02 + 10) continue;
				matches.emplace_back(pos, pos2);
			}
		}
	}
	if (matches.size() == 0)
	{
		return result;
	}
	std::sort(matches.begin(), matches.end());
	std::vector<size_t> bestPredecessor;
	bestPredecessor.resize(matches.size(), std::numeric_limits<size_t>::max());
	std::vector<size_t> bestMatchChainLength;
	bestMatchChainLength.resize(matches.size(), 0);
	for (size_t i = 1; i < matches.size(); i++)
	{
		for (size_t j = i-1; j < i; j--)
		{
			assert(matches[i].first >= matches[j].first);
			if (matches[i].first == matches[j].first) continue;
			if (matches[i].second <= matches[j].second) continue;
			int refDistance = (int)std::get<0>(chunksPerRead[referenceRead][matches[i].first]) - (int)std::get<1>(chunksPerRead[referenceRead][matches[j].first]);
			int queryDistance = (int)std::get<0>(chunksPerRead[queryRead][matches[i].second]) - (int)std::get<1>(chunksPerRead[queryRead][matches[j].second]);
			if (refDistance > queryDistance+50) continue;
			if (queryDistance > refDistance+50) continue;
			if (bestPredecessor[i] == std::numeric_limits<size_t>::max() || bestMatchChainLength[j] > bestMatchChainLength[bestPredecessor[i]])
			{
				bestPredecessor[i] = j;
			}
		}
		bestMatchChainLength[i] = std::get<1>(chunksPerRead[referenceRead][matches[i].first]) - std::get<0>(chunksPerRead[referenceRead][matches[i].first]);
		if (bestPredecessor[i] != std::numeric_limits<size_t>::max())
		{
			bestMatchChainLength[i] += bestMatchChainLength[bestPredecessor[i]];
		}
	}
	std::vector<bool> matchIsChainEnd;
	matchIsChainEnd.resize(matches.size(), true);
	for (size_t i = 0; i < matches.size(); i++)
	{
		if (bestPredecessor[i] == std::numeric_limits<size_t>::max()) continue;
		assert(bestPredecessor[i] < i);
		matchIsChainEnd[bestPredecessor[i]] = false;
	}
	bool bestChainIsUnique = true;
	size_t bestChainIndex = 0;
	for (size_t i = 0; i < matches.size(); i++)
	{
		if (!matchIsChainEnd[i]) continue;
		if (bestMatchChainLength[i] > bestMatchChainLength[bestChainIndex]) bestChainIndex = i;
	}
	for (size_t i = 0; i < matches.size(); i++)
	{
		if (i == bestChainIndex) continue;
		if (!matchIsChainEnd[i]) continue;
		if (bestMatchChainLength[i]*1.05+100 > bestMatchChainLength[bestChainIndex]) bestChainIsUnique = false;
	}
	if (!bestChainIsUnique)
	{
		return result;
	}
	std::vector<std::pair<size_t, size_t>> bestChain;
	size_t currentChainIndex = bestChainIndex;
	std::vector<bool> usedInBestChain;
	usedInBestChain.resize(matches.size(), false);
	while (true)
	{
		bestChain.emplace_back(matches[currentChainIndex]);
		usedInBestChain[currentChainIndex] = true;
		if (bestPredecessor[currentChainIndex] == std::numeric_limits<size_t>::max()) break;
		currentChainIndex = bestPredecessor[currentChainIndex];
	}
	currentChainIndex = bestChainIndex;
	while (true)
	{
		if (bestPredecessor[currentChainIndex] == std::numeric_limits<size_t>::max()) break;
		for (size_t j = 0; j < currentChainIndex; j++)
		{
			if (usedInBestChain[j]) continue;
			assert(matches[j].first <= matches[currentChainIndex].first);
			if (matches[j].first == matches[currentChainIndex].first) continue;
			if (matches[j].second >= matches[currentChainIndex].second) continue;
			if (bestMatchChainLength[j]+100 > bestMatchChainLength[bestPredecessor[currentChainIndex]])
			{
				bestChainIsUnique = false;
				return result;
			}
		}
		currentChainIndex = bestPredecessor[currentChainIndex];
	}
	assert(bestChain.size() >= 1);
	bool bestChainCoversEnough = true;
	size_t missingAtRefEnd = std::get<1>(chunksPerRead[referenceRead].back()) - std::get<1>(chunksPerRead[referenceRead][bestChain[0].first]);
	size_t missingAtQueryEnd = std::get<1>(chunksPerRead[queryRead].back()) - std::get<1>(chunksPerRead[queryRead][bestChain[0].second]);
	size_t missingAtRefStart = std::get<0>(chunksPerRead[referenceRead][bestChain.back().first]) - std::get<0>(chunksPerRead[referenceRead][0]);
	size_t missingAtQueryStart = std::get<0>(chunksPerRead[queryRead][bestChain.back().second]) - std::get<0>(chunksPerRead[queryRead][0]);
	if (missingAtRefEnd > 20000 && missingAtQueryEnd > 20000) bestChainCoversEnough = false;
	if (missingAtRefStart > 20000 && missingAtQueryStart > 20000) bestChainCoversEnough = false;
	size_t biggestMiddleGap = 0;
	for (size_t i = 1; i < bestChain.size(); i++)
	{
		assert(bestChain[i].first < bestChain[i-1].first);
		if (std::get<0>(chunksPerRead[referenceRead][bestChain[i-1].first]) > std::get<1>(chunksPerRead[referenceRead][bestChain[i].first]))
		{
			biggestMiddleGap = std::max(biggestMiddleGap, std::get<0>(chunksPerRead[referenceRead][bestChain[i-1].first]) - std::get<1>(chunksPerRead[referenceRead][bestChain[i].first]));
		}
		if (std::get<0>(chunksPerRead[queryRead][bestChain[i-1].second]) > std::get<1>(chunksPerRead[queryRead][bestChain[i].second]))
		{
			biggestMiddleGap = std::max(biggestMiddleGap, std::get<0>(chunksPerRead[queryRead][bestChain[i-1].second]) - std::get<1>(chunksPerRead[queryRead][bestChain[i].second]));
		}
	}
	if (biggestMiddleGap > 20000) bestChainCoversEnough = false;
	if (!bestChainCoversEnough)
	{
		return result;
	}
	std::swap(result.matches, bestChain);
	std::reverse(result.matches.begin(), result.matches.end());
	return result;
}

std::vector<std::vector<Overlap>> getOverlaps(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t maxRepeatChunkCoverage = 10000;
	assert(chunksPerRead.size() % 2 == 0);
	auto chunkIndicesWithinReads = getChunkIndicesWithinReads(chunksPerRead);
	auto readsWhichHaveChunk = getReadsWhichHaveChunk(chunkIndicesWithinReads);
	std::vector<std::vector<Overlap>> result;
	result.resize(chunksPerRead.size()/2);
	iterateMultithreaded(0, chunksPerRead.size()/2, numThreads, [&chunksPerRead, &chunkIndicesWithinReads, &readsWhichHaveChunk, &result, maxRepeatChunkCoverage](const size_t i)
	{
		phmap::flat_hash_set<size_t> possibleOverlapReads;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			assert(chunk < readsWhichHaveChunk.size());
			if (readsWhichHaveChunk[chunk].size() >= maxRepeatChunkCoverage) continue;
			for (auto read : readsWhichHaveChunk[chunk])
			{
				possibleOverlapReads.emplace(read);
			}
		}
		for (size_t read : possibleOverlapReads)
		{
			auto overlap = getUniqueBestOverlap(chunkIndicesWithinReads, chunksPerRead, i, read);
			if (overlap.matches.size() > 0)
			{
				result[i].emplace_back();
				std::swap(result[i].back(), overlap);
			}
		}
	});
	size_t countOverlaps = 0;
	for (size_t i = 0; i < result.size(); i++)
	{
		countOverlaps += result[i].size();
	}
	std::cerr << "count overlaps: " << countOverlaps << std::endl;
	return result;
}

std::vector<Correction> getSolidEdits(const std::vector<Overlap>& overlaps, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	const size_t solidThreshold = 5;
	std::map<Correction, size_t> editCoverage;
	for (const auto& overlap : overlaps)
	{
		assert(overlap.matches.size() >= 1);
		for (size_t i = 1; i < overlap.matches.size(); i++)
		{
			assert(overlap.matches[i].first > overlap.matches[i-1].first);
			assert(overlap.matches[i].second > overlap.matches[i-1].second);
			if (overlap.matches[i].first == overlap.matches[i-1].first+1 && overlap.matches[i].second == overlap.matches[i-1].second+1) continue;
			Correction replacement;
			replacement.startIndex = overlap.matches[i-1].first;
			replacement.endIndex = overlap.matches[i].first;
			assert(replacement.endIndex > replacement.startIndex);
			bool anyInvalid = false;
			for (size_t j = overlap.matches[i-1].second+1; j < overlap.matches[i].second; j++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[overlap.queryRead][j]))) anyInvalid = true;
				replacement.chunkReplacement.emplace_back(std::get<2>(chunksPerRead[overlap.queryRead][j]) & maskUint64_t);
			}
			if (anyInvalid) continue;
			for (size_t j = overlap.matches[i-1].first; j <= overlap.matches[i].second; j++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[overlap.referenceRead][j]))) anyInvalid = true;
			}
			if (anyInvalid) continue;
			editCoverage[replacement] += 1;
		}
	}
	std::vector<Correction> result;
	for (const auto& pair : editCoverage)
	{
		if (pair.second < solidThreshold) continue;
		result.emplace_back(pair.first);
	}
	return result;
}

std::vector<std::vector<Correction>> getPossibleCorrections(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const std::vector<std::vector<Overlap>>& overlaps, const size_t numThreads)
{
	assert(chunksPerRead.size() % 2 == 0);
	std::vector<std::vector<Correction>> result;
	result.resize(chunksPerRead.size()/2);
	iterateMultithreaded(0, overlaps.size(), numThreads, [&chunksPerRead, &overlaps, &result](const size_t i)
	{
		auto solidEdits = getSolidEdits(overlaps[i], chunksPerRead);
		result[i].insert(result[i].end(), solidEdits.begin(), solidEdits.end());
	});
	size_t countCorrections = 0;
	for (size_t i = 0; i < result.size(); i++)
	{
		countCorrections += result[i].size();
	}
	std::cerr << "count corrections: " << countCorrections << std::endl;
	return result;
}

std::pair<std::vector<uint64_t>, std::vector<uint64_t>> getChunkBoundingKmers(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numChunks, const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t kmerSize, const size_t numThreads)
{
	std::vector<uint64_t> chunkFirstKmer;
	std::vector<uint64_t> chunkLastKmer;
	chunkFirstKmer.resize(numChunks, std::numeric_limits<uint64_t>::max());
	chunkLastKmer.resize(numChunks, std::numeric_limits<uint64_t>::max());
	std::mutex chunkKmerMutex;
	iterateMultithreaded(0, chunksPerRead.size(), numThreads, [&chunksPerRead, &chunkFirstKmer, &chunkLastKmer, &chunkKmerMutex, &sequenceIndex, kmerSize](const size_t i)
	{
		std::string readSequence;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			{
				std::lock_guard<std::mutex> lock { chunkKmerMutex };
				if (chunkFirstKmer[chunk] != std::numeric_limits<uint64_t>::max() && chunkLastKmer[chunk] != std::numeric_limits<uint64_t>::max()) continue;
				if (chunkFirstKmer[chunk] == std::numeric_limits<uint64_t>::max()) chunkFirstKmer[chunk] = std::numeric_limits<uint64_t>::max()-1;
				if (chunkLastKmer[chunk] == std::numeric_limits<uint64_t>::max()) chunkLastKmer[chunk] = std::numeric_limits<uint64_t>::max()-1;
			}
			if (readSequence.size() == 0)
			{
				if (i >= sequenceIndex.size())
				{
					readSequence = sequenceIndex.getSequence(i-sequenceIndex.size());
					readSequence = revCompRaw(readSequence);
				}
				else
				{
					readSequence = sequenceIndex.getSequence(i);
				}
			}
			uint64_t firstKmer = 0;
			uint64_t lastKmer = 0;
			for (size_t l = 0; l < kmerSize; l++)
			{
				firstKmer <<= 2;
				lastKmer <<= 2;
				switch(readSequence[std::get<0>(chunksPerRead[i][j])+l])
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
				switch(readSequence[std::get<1>(chunksPerRead[i][j])-kmerSize+1+l])
				{
					case 'A':
						lastKmer += 0;
						break;
					case 'C':
						lastKmer += 1;
						break;
					case 'G':
						lastKmer += 2;
						break;
					case 'T':
						lastKmer += 3;
						break;
					default:
						assert(false);
				}
			}
			{
				std::lock_guard<std::mutex> lock { chunkKmerMutex };
				assert(chunkFirstKmer[chunk] != std::numeric_limits<uint64_t>::max());
				assert(chunkLastKmer[chunk] != std::numeric_limits<uint64_t>::max());
				if (chunkFirstKmer[chunk] == std::numeric_limits<uint64_t>::max()-1)
				{
					chunkFirstKmer[chunk] = firstKmer;
				}
				else
				{
					assert(chunkFirstKmer[chunk] == firstKmer);
				}
				if (chunkLastKmer[chunk] == std::numeric_limits<uint64_t>::max()-1)
				{
					chunkLastKmer[chunk] = lastKmer;
				}
				else
				{
					assert(chunkLastKmer[chunk] == lastKmer);
				}
			}
		}
	});
	return std::make_pair(chunkFirstKmer, chunkLastKmer);
}

void addEditKmers(Correction& correction, const std::string& readSequence, const std::vector<std::tuple<size_t, size_t, uint64_t>>& chunksInThisRead, const std::vector<size_t>& minimizersInThisRead, const size_t kmerSize, const std::vector<size_t>& chunkFirstKmer, const std::vector<size_t>& chunkLastKmer, const std::vector<size_t>& chunkMeanLength)
{
	assert(minimizersInThisRead.size() == chunksInThisRead.size()+1);
	size_t replacementLength = 0;
	for (size_t i = 0; i < correction.chunkReplacement.size(); i++)
	{
		replacementLength += chunkMeanLength[correction.chunkReplacement[i]];
	}
	if (correction.chunkReplacement.size() == 0)
	{
		if (correction.endIndex == correction.startIndex+1)
		{
			correction.minimizerReplacement.emplace_back(minimizersInThisRead[correction.startIndex+1]);
		}
		return;
	}
	if (!(replacementLength > (kmerSize-1)*(correction.chunkReplacement.size()-1)))
	{
		std::cerr << replacementLength << " " << ((kmerSize-1)*(correction.chunkReplacement.size()-1)) << " " << correction.chunkReplacement.size() << std::endl;
	}
	assert(replacementLength > (kmerSize-1)*(correction.chunkReplacement.size()-1));
	replacementLength -= (kmerSize-1)*(correction.chunkReplacement.size()-1);
	std::vector<uint64_t> wantedKmerSequence;
	phmap::flat_hash_map<uint64_t, std::vector<std::pair<double, size_t>>> wantedKmerFractionalPositionAndIndex;
	size_t currentPos = 0;
	for (size_t i = 0; i < correction.chunkReplacement.size(); i++)
	{
		uint64_t kmer = chunkFirstKmer[correction.chunkReplacement[i]];
		wantedKmerSequence.emplace_back(kmer);
		wantedKmerFractionalPositionAndIndex[kmer].emplace_back((double)currentPos/(double)replacementLength, i);
	}
	wantedKmerSequence.emplace_back(chunkLastKmer[correction.chunkReplacement.back()]);
	wantedKmerFractionalPositionAndIndex[correction.chunkReplacement.back()].emplace_back(1, correction.chunkReplacement.size());
	std::vector<size_t> replacementPositions;
	replacementPositions.resize(wantedKmerSequence.size(), std::numeric_limits<size_t>::max());
	replacementPositions[0] = minimizersInThisRead[correction.startIndex+1];
	replacementPositions.back() = minimizersInThisRead[correction.endIndex];
	if (replacementPositions.size() > 2)
	{
		iterateKmersChar(readSequence, replacementPositions[0]+1, replacementPositions.back()-1, kmerSize, [&wantedKmerFractionalPositionAndIndex, &replacementPositions](const size_t kmer, const size_t pos)
		{
			if (wantedKmerFractionalPositionAndIndex.count(kmer) == 0) return;
			for (const auto& pair : wantedKmerFractionalPositionAndIndex.at(kmer))
			{
				size_t wantedBaseLocation = pair.first * (replacementPositions.back() - replacementPositions[0]) + replacementPositions[0];
				if (pos+100 < wantedBaseLocation) continue;
				if (wantedBaseLocation+100 < pos) continue;
				if (replacementPositions[pair.second] != std::numeric_limits<size_t>::max())
				{
					replacementPositions[pair.second] = std::numeric_limits<size_t>::max()-1;
					return;
				}
				replacementPositions[pair.second] = pos;
			}
		});
	}
	bool anyInvalid = false;
	for (size_t i = 0; i < replacementPositions.size(); i++)
	{
		if (replacementPositions[i] == std::numeric_limits<size_t>::max()) anyInvalid = true;
		if (replacementPositions[i] == std::numeric_limits<size_t>::max()-1) anyInvalid = true;
	}
	for (size_t i = 1; i < replacementPositions.size(); i++)
	{
		if (replacementPositions[i] <= replacementPositions[i-1]) anyInvalid = true;
	}
	if (anyInvalid) return;
	correction.minimizerReplacement = replacementPositions;
}

void addCorrectionMinimizerSequences(std::vector<std::vector<Correction>>& edits, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const size_t kmerSize, const size_t numThreads)
{
	std::vector<size_t> chunkMeanLength;
	std::vector<size_t> chunkCoverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			while (chunkMeanLength.size() <= chunk) chunkMeanLength.emplace_back(0);
			while (chunkCoverage.size() <= chunk) chunkCoverage.emplace_back(0);
			chunkMeanLength[chunk] += std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]);
			chunkCoverage[chunk] += 1;
		}
	}
	for (size_t i = 0; i < chunkMeanLength.size(); i++)
	{
		assert(chunkCoverage[i] >= 1);
		chunkMeanLength[i] /= chunkCoverage[i];
	}
	std::vector<uint64_t> chunkFirstKmer;
	std::vector<uint64_t> chunkLastKmer;
	std::tie(chunkFirstKmer, chunkLastKmer) = getChunkBoundingKmers(chunksPerRead, chunkMeanLength.size(), sequenceIndex, kmerSize, numThreads);
	iterateMultithreaded(0, edits.size(), numThreads, [&edits, &chunksPerRead, &minimizerPositionsPerRead, &sequenceIndex, &chunkFirstKmer, &chunkLastKmer, &chunkMeanLength, kmerSize](const size_t i)
	{
		if (edits[i].size() == 0) return;
		std::string readSequence = sequenceIndex.getSequence(i);
		for (size_t j = 0; j < edits[i].size(); j++)
		{
			addEditKmers(edits[i][j], readSequence, chunksPerRead[i], minimizerPositionsPerRead[i], kmerSize, chunkFirstKmer, chunkLastKmer, chunkMeanLength);
		}
	});
}

void filterPossibleCorrections(std::vector<std::vector<Correction>>& edits, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads)
{
	const size_t solidChunkThreshold = 20;
	const size_t solidChunkEdgeThreshold = 5;
	std::vector<size_t> chunkCoverage;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			while (chunkCoverage.size() <= chunk) chunkCoverage.emplace_back(0);
			chunkCoverage[chunk] += 1;
		}
	}
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> adjacentChunkCoverages;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j-1]))) continue;
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			size_t prevchunk = std::get<2>(chunksPerRead[i][j-1]) & maskUint64_t;
			size_t currchunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
			std::pair<size_t, size_t> key { prevchunk, currchunk };
			adjacentChunkCoverages[key] += 1;
		}
	}
	iterateMultithreaded(0, edits.size(), numThreads, [&edits, &chunkCoverage, &chunksPerRead, &adjacentChunkCoverages](const size_t i)
	{
		if (edits[i].size() == 0) return;
		std::vector<bool> editCanBeDone;
		editCanBeDone.resize(edits[i].size(), true);
		for (size_t j = 0; j < edits[i].size(); j++)
		{
			if (edits[i][j].minimizerReplacement.size() == 0)
			{
				editCanBeDone[j] = false;
				continue;
			}
			bool allSolids = true;
			for (size_t k = edits[i][j].startIndex; k <= edits[i][j].endIndex; k++)
			{
				size_t chunk = std::get<2>(chunksPerRead[i][k]) & maskUint64_t;
				if (chunkCoverage[chunk] < solidChunkThreshold) allSolids = false;
			}
			if (allSolids) editCanBeDone[j] = false;
			if (!editCanBeDone[j]) continue;
			bool readPathSolid = true;
			for (size_t k = edits[i][j].startIndex+1; k < edits[i][j].endIndex; k++)
			{
				size_t prevchunk = std::get<2>(chunksPerRead[i][k-1]) & maskUint64_t;
				size_t currchunk = std::get<2>(chunksPerRead[i][k]) & maskUint64_t;
				std::pair<size_t, size_t> key { prevchunk, currchunk };
				if (adjacentChunkCoverages.at(key) < solidChunkEdgeThreshold) readPathSolid = false;
			}
			if (readPathSolid) editCanBeDone[j] = false;
			if (!editCanBeDone[j]) continue;
			for (size_t k = 0; k < edits[i].size(); k++)
			{
				if (k == j) continue;
				if (edits[i][k].startIndex < edits[i][j].endIndex && edits[i][k].endIndex > edits[i][j].startIndex)
				{
					editCanBeDone[j] = false;
					editCanBeDone[k] = false;
					break;
				}
			}
		}
		for (size_t j = edits[i].size()-1; j < edits[i].size(); j--)
		{
			if (editCanBeDone[j]) continue;
			std::swap(edits[i][j], edits[i].back());
			edits[i].pop_back();
		}
		std::sort(edits[i].begin(), edits[i].end(), [](const auto& left, const auto& right)
		{
			return left.startIndex < right.startIndex;
		});
	});
	size_t countCorrections = 0;
	for (size_t i = 0; i < edits.size(); i++)
	{
		countCorrections += edits[i].size();
	}
	std::cerr << "count solid corrections: " << countCorrections << std::endl;
}

void applyCorrections(const std::vector<std::vector<Correction>>& corrections, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const std::vector<size_t>& rawReadLengths, const size_t kmerSize, const size_t numThreads)
{
	iterateMultithreaded(0, corrections.size(), numThreads, [&corrections, &chunksPerRead, &minimizerPositionsPerRead, kmerSize](const size_t i)
	{
		if (corrections[i].size() == 0)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				std::get<2>(chunksPerRead[i][j]) = 0;
			}
			return;
		}
		chunksPerRead[i].clear();
		for (size_t j = 1; j < corrections[i].size(); j++)
		{
			assert(corrections[i][j].startIndex >= corrections[i][j-1].endIndex);
		}
		std::vector<size_t> newMinimizerPositions;
		size_t lastMatch = 0;
		for (size_t j = 0; j < corrections[i].size(); j++)
		{
			if (corrections[i][j].startIndex >= lastMatch)
			{
				newMinimizerPositions.insert(newMinimizerPositions.end(), minimizerPositionsPerRead[i].begin()+lastMatch, minimizerPositionsPerRead[i].begin()+corrections[i][j].startIndex);
			}
			newMinimizerPositions.insert(newMinimizerPositions.end(), corrections[i][j].minimizerReplacement.begin(), corrections[i][j].minimizerReplacement.end());
			lastMatch = corrections[i][j].endIndex+1;
		}
		newMinimizerPositions.insert(newMinimizerPositions.end(), minimizerPositionsPerRead[i].begin()+lastMatch, minimizerPositionsPerRead[i].end());
		std::swap(newMinimizerPositions, minimizerPositionsPerRead[i]);
		for (size_t j = 1; j < minimizerPositionsPerRead[i].size(); j++)
		{
			chunksPerRead[i].emplace_back(minimizerPositionsPerRead[i][j-1], minimizerPositionsPerRead[i][j-1]+kmerSize-1, 0);
		}
	});
	iterateMultithreaded(0, corrections.size(), numThreads, [&chunksPerRead, &corrections, &minimizerPositionsPerRead, &rawReadLengths, kmerSize](const size_t readIndex)
	{
		chunksPerRead[readIndex+corrections.size()].clear();
		minimizerPositionsPerRead[readIndex+corrections.size()].clear();
		for (size_t i = 0; i < minimizerPositionsPerRead[readIndex].size(); i++)
		{
			minimizerPositionsPerRead[readIndex+corrections.size()].emplace_back(rawReadLengths[readIndex]-kmerSize-minimizerPositionsPerRead[readIndex][minimizerPositionsPerRead[readIndex].size()-1-i]);
		}
		for (size_t i = 0; i < chunksPerRead[readIndex].size(); i++)
		{
			chunksPerRead[readIndex+corrections.size()].emplace_back(chunksPerRead[readIndex][chunksPerRead[readIndex].size()-1-i]);
			std::swap(std::get<0>(chunksPerRead[readIndex+corrections.size()].back()), std::get<1>(chunksPerRead[readIndex+corrections.size()].back()));
			std::get<0>(chunksPerRead[readIndex+corrections.size()].back()) = rawReadLengths[readIndex]-1-std::get<0>(chunksPerRead[readIndex+corrections.size()].back());
			std::get<1>(chunksPerRead[readIndex+corrections.size()].back()) = rawReadLengths[readIndex]-1-std::get<1>(chunksPerRead[readIndex+corrections.size()].back());
		}
	});
}

void correctChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const std::vector<size_t>& rawReadLengths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t kmerSize, const size_t numThreads)
{
	std::cerr << "correct chunks" << std::endl;
	assert(chunksPerRead.size() == 2*sequenceIndex.size());
	std::cerr << "get overlaps" << std::endl;
	auto overlaps = getOverlaps(chunksPerRead, numThreads);
	std::cerr << "get corrections" << std::endl;
	auto possibleCorrections = getPossibleCorrections(chunksPerRead, overlaps, numThreads);
	std::cerr << "get correction minimizers" << std::endl;
	addCorrectionMinimizerSequences(possibleCorrections, chunksPerRead, sequenceIndex, minimizerPositionsPerRead, kmerSize, numThreads);
	std::cerr << "filter corrections" << std::endl;
	filterPossibleCorrections(possibleCorrections, chunksPerRead, numThreads);
	std::cerr << "apply corrections" << std::endl;
	applyCorrections(possibleCorrections, chunksPerRead, minimizerPositionsPerRead, rawReadLengths, kmerSize, numThreads);
}
