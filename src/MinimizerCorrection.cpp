#include <iostream>
#include <string>
#include <phmap.h>
#include "Common.h"
#include "KmerIterator.h"
#include "MinimizerCorrection.h"
#include "UnionFind.h"

std::vector<std::pair<size_t, size_t>> getRawMinimizers(const std::string& sequence, const size_t kmerSize, const size_t windowSize)
{
	std::vector<std::pair<size_t, size_t>> result;
	iterateMinimizers(sequence, kmerSize, windowSize, [&result](const size_t pos)
	{
		result.emplace_back(pos, 0);
	});
	for (size_t j = 0; j < result.size(); j++)
	{
		uint64_t kmer = 0;
		for (size_t l = 0; l < kmerSize; l++)
		{
			kmer <<= 2;
			kmer += charToInt(sequence[result[j].first+l]);
		}
		result[j].first += kmerSize/2;
		result[j].second = kmer;
	}
	return result;
}

struct Overlap
{
	size_t referenceRead;
	size_t queryRead;
	std::vector<std::pair<size_t, size_t>> matches;
};

void removeConflictMatches(std::vector<std::pair<size_t, size_t>>& matches)
{
	std::vector<bool> conflict;
	conflict.resize(matches.size(), false);
	bool removedAny = false;
	for (size_t i = 1; i < matches.size(); i++)
	{
		if (matches[i].first > matches[i-1].first && matches[i].second > matches[i-1].second) continue;
		conflict[i] = true;
		conflict[i-1] = true;
		removedAny = true;
	}
	if (!removedAny) return;
	for (size_t i = matches.size()-1; i < matches.size(); i--)
	{
		if (!conflict[i]) continue;
		std::swap(matches[i], matches.back());
		matches.pop_back();
	}
	std::sort(matches.begin(), matches.end());
	removeConflictMatches(matches);
}

Overlap getBestOverlap(const std::vector<phmap::flat_hash_map<std::pair<size_t, size_t>, std::vector<size_t>>>& minimizerPairIndex, const std::vector<std::vector<std::pair<size_t, size_t>>>& rawMinimizers, const size_t referenceRead, const size_t queryRead)
{
	Overlap result;
	result.referenceRead = referenceRead;
	result.queryRead = queryRead;
	std::vector<std::pair<int, size_t>> diagonalMatches;
	for (const auto& pair : minimizerPairIndex[referenceRead])
	{
		if (minimizerPairIndex[queryRead].count(pair.first) == 0) continue;
		if (pair.second.size()*minimizerPairIndex[queryRead].at(pair.first).size() > 10000) continue;
		for (const auto pos2 : minimizerPairIndex[queryRead].at(pair.first))
		{
			for (const size_t pos : pair.second)
			{
				size_t refLength = rawMinimizers[referenceRead][pos].first - rawMinimizers[referenceRead][pos-1].first;
				size_t queryLength = rawMinimizers[queryRead][pos2].first - rawMinimizers[queryRead][pos2-1].first;
				if (refLength > queryLength*1.02+50) continue;
				if (queryLength > refLength*1.02+50) continue;
				int diagonal = (int)rawMinimizers[referenceRead][pos].first - (int)rawMinimizers[queryRead][pos2].first;
				diagonalMatches.emplace_back(diagonal, refLength);
			}
		}
	}
	if (diagonalMatches.size() == 0) return result;
	std::sort(diagonalMatches.begin(), diagonalMatches.end());
	int bestDiagonalStart = diagonalMatches[0].first;
	int bestDiagonalEnd = diagonalMatches[0].first;
	size_t bestDiagonalLengthSum = diagonalMatches[0].second;
	size_t clusterStart = 0;
	for (size_t i = 1; i <= diagonalMatches.size(); i++)
	{
		assert(i == diagonalMatches.size() || diagonalMatches[i].first >= diagonalMatches[i-1].first);
		if (i < diagonalMatches.size() && diagonalMatches[i].first < diagonalMatches[i-1].first+100) continue;
		size_t lengthSumHere = 0;
		for (size_t j = clusterStart; j < i; j++)
		{
			lengthSumHere += diagonalMatches[j].second;
		}
		if (lengthSumHere == bestDiagonalLengthSum)
		{
			bestDiagonalStart = 0;
			bestDiagonalEnd = 0;
		}
		if (lengthSumHere > bestDiagonalLengthSum)
		{
			bestDiagonalStart = diagonalMatches[clusterStart].first;
			bestDiagonalEnd = diagonalMatches[i-1].first;
			bestDiagonalLengthSum = lengthSumHere;
		}
		clusterStart = i;
	}
	if (bestDiagonalEnd == bestDiagonalStart) return result;
	if (bestDiagonalLengthSum < 5000) return result;
	std::vector<std::pair<size_t, size_t>> matches;
	for (const auto& pair : minimizerPairIndex[referenceRead])
	{
		if (minimizerPairIndex[queryRead].count(pair.first) == 0) continue;
		for (const auto pos2 : minimizerPairIndex[queryRead].at(pair.first))
		{
			for (const auto pos : pair.second)
			{
				size_t refLength = rawMinimizers[referenceRead][pos].first - rawMinimizers[referenceRead][pos-1].first;
				size_t queryLength = rawMinimizers[queryRead][pos2].first - rawMinimizers[queryRead][pos2-1].first;
				if (refLength > queryLength*1.02+50) continue;
				if (queryLength > refLength*1.02+50) continue;
				int diagonal = (int)rawMinimizers[referenceRead][pos].first - (int)rawMinimizers[queryRead][pos2].first;
				if (diagonal < bestDiagonalStart) continue;
				if (diagonal > bestDiagonalEnd) continue;
				matches.emplace_back(pos, pos2);
			}
		}
	}
	phmap::flat_hash_set<std::pair<size_t, size_t>> matchPoses;
	for (auto pair : matches)
	{
		matchPoses.emplace(pair.first, pair.second);
		matchPoses.emplace(pair.first-1, pair.second-1);
	}
	std::vector<std::pair<size_t, size_t>> matchPosesVec { matchPoses.begin(), matchPoses.end() };
	std::sort(matchPosesVec.begin(), matchPosesVec.end());
	removeConflictMatches(matchPosesVec);
	result.matches = matchPosesVec;
	return result;
}

std::vector<std::vector<Overlap>> getOverlaps(const std::vector<std::vector<std::pair<size_t, size_t>>> rawMinimizers, const size_t numThreads)
{
	assert(rawMinimizers.size() % 2 == 0);
	std::vector<phmap::flat_hash_map<std::pair<size_t, size_t>, std::vector<size_t>>> minimizerPairIndex;
	minimizerPairIndex.resize(rawMinimizers.size());
	iterateMultithreaded(0, rawMinimizers.size(), numThreads, [&minimizerPairIndex, &rawMinimizers](const size_t i)
	{
		for (size_t j = 1; j < rawMinimizers[i].size(); j++)
		{
			std::pair<size_t, size_t> key { rawMinimizers[i][j-1].second, rawMinimizers[i][j].second };
			minimizerPairIndex[i][key].emplace_back(j);
		}
	});
	phmap::flat_hash_map<std::pair<size_t, size_t>, std::vector<size_t>> readsWithMinimizerPair;
	for (size_t i = 0; i < rawMinimizers.size(); i++)
	{
		for (const auto& pair : minimizerPairIndex[i])
		{
			readsWithMinimizerPair[pair.first].emplace_back(i);
		}
	}
	std::vector<std::vector<Overlap>> result;
	result.resize(rawMinimizers.size()/2);
	iterateMultithreaded(0, rawMinimizers.size()/2, numThreads, [&result, &rawMinimizers, &minimizerPairIndex, &readsWithMinimizerPair](const size_t i)
	{
		phmap::flat_hash_set<size_t> readsPossiblyOverlapping;
		for (const auto& pair : minimizerPairIndex[i])
		{
			assert(readsWithMinimizerPair.count(pair.first) == 1);
			if (readsWithMinimizerPair.at(pair.first).size() > 1000) continue;
			readsPossiblyOverlapping.insert(readsWithMinimizerPair.at(pair.first).begin(), readsWithMinimizerPair.at(pair.first).end());
		}
		for (size_t read : readsPossiblyOverlapping)
		{
			auto overlap = getBestOverlap(minimizerPairIndex, rawMinimizers, i, read);
			if (overlap.matches.size() >= 2)
			{
				result[i].emplace_back(overlap);
			}
		}
		std::string info = "read " + std::to_string(i) + " has " + std::to_string(readsPossiblyOverlapping.size()) + " reads possibly overlapping, " + std::to_string(result[i].size()) + " overlaps\n";
		std::cerr << info;
	});
	size_t countOverlaps = 0;
	for (size_t i = 0; i < result.size(); i++)
	{
		countOverlaps += result[i].size();
	}
	std::cerr << "count overlaps: " << countOverlaps << std::endl;
	return result;
}

void addReverseMinimizers(std::vector<std::vector<size_t>>& rawMinimizers, const std::vector<size_t>& rawReadLengths)
{
	size_t oldSize = rawMinimizers.size();
	rawMinimizers.resize(oldSize*2);
	for (size_t i = 0; i < oldSize; i++)
	{
		for (const size_t pos : rawMinimizers[i])
		{
			rawMinimizers[i+oldSize].emplace_back(rawReadLengths[i]-1-pos);
		}
		std::sort(rawMinimizers[i+oldSize].begin(), rawMinimizers[i+oldSize].end());
	}
}

void makeReverseReads(std::vector<std::vector<std::pair<size_t, size_t>>>& rawMinimizers, const std::vector<size_t>& rawReadLengths, const size_t kmerSize)
{
	size_t oldSize = rawMinimizers.size();
	rawMinimizers.resize(oldSize*2);
	for (size_t i = 0; i < oldSize; i++)
	{
		for (auto pair : rawMinimizers[i])
		{
			rawMinimizers[i+oldSize].emplace_back(rawReadLengths[i]-1-pair.first, reverseKmer(pair.second, kmerSize));
		}
		std::sort(rawMinimizers[i+oldSize].begin(), rawMinimizers[i+oldSize].end());
	}
}

struct Edit
{
	size_t startIndex; // index of pre-anchor
	size_t endIndex; // inclusive, ie index of post-anchor
	std::vector<size_t> replacementKmers; // includes startIndex, endIndex
	std::vector<double> replacementKmerPositions;
	std::vector<size_t> actualKmerPositionsInReference;
	bool operator<(const Edit& other) const
	{
		if (startIndex < other.startIndex) return true;
		if (startIndex > other.startIndex) return false;
		if (endIndex < other.endIndex) return true;
		if (endIndex > other.endIndex) return false;
		if (replacementKmers < other.replacementKmers) return true;
		if (replacementKmers > other.replacementKmers) return false;
		if (replacementKmerPositions < other.replacementKmerPositions) return true;
		if (replacementKmerPositions > other.replacementKmerPositions) return false;
		return false;
	}
};

bool kmerPositionsApproxMatch(const std::vector<double>& left, const std::vector<double>& right)
{
	return true;
	assert(left.size() == right.size());
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] > right[i]+0.1) return false;
		if (left[i] < right[i]-0.1) return false;
	}
	return true;
}

std::vector<Edit> getCoveredEdits(const std::vector<std::vector<std::pair<size_t, size_t>>>& rawMinimizers, const std::vector<Overlap>& overlaps)
{
	if (overlaps.size() == 0) return std::vector<Edit> {};
	const size_t minCoverage = 5;
	std::vector<Edit> edits;
	for (const auto& overlap : overlaps)
	{
		for (size_t i = 1; i < overlap.matches.size(); i++)
		{
			Edit thisEdit;
			thisEdit.startIndex = overlap.matches[i-1].first;
			thisEdit.endIndex = overlap.matches[i].first;
			size_t lengthInQuery = rawMinimizers[overlap.queryRead][overlap.matches[i].second].first - rawMinimizers[overlap.queryRead][overlap.matches[i-1].second].first;
			for (size_t j = overlap.matches[i-1].second; j <= overlap.matches[i].second; j++)
			{
				thisEdit.replacementKmers.emplace_back(rawMinimizers[overlap.queryRead][j].second);
				thisEdit.replacementKmerPositions.emplace_back((double)(rawMinimizers[overlap.queryRead][j].first - rawMinimizers[overlap.queryRead][overlap.matches[i-1].second].first) / (double)lengthInQuery);
				assert(thisEdit.replacementKmerPositions.back() <= 1.0);
				assert(thisEdit.replacementKmerPositions.back() >= 0.0);
			}
			edits.emplace_back(thisEdit);
		}
	}
	std::sort(edits.begin(), edits.end());
	std::vector<size_t> parent;
	for (size_t i = 0; i < edits.size(); i++)
	{
		parent.emplace_back(i);
	}
	size_t clusterStart = 0;
	for (size_t i = 1; i <= edits.size(); i++)
	{
		if (i < edits.size() && edits[i].startIndex == edits[i-1].startIndex && edits[i].endIndex == edits[i-1].endIndex && edits[i].replacementKmers == edits[i-1].replacementKmers) continue;
		for (size_t j = clusterStart+1; j < i; j++)
		{
			for (size_t k = clusterStart; k < j; k++)
			{
				if (kmerPositionsApproxMatch(edits[j].replacementKmerPositions, edits[k].replacementKmerPositions))
				{
					merge(parent, j, k);
				}
			}
		}
		clusterStart = i;
	}
	phmap::flat_hash_map<size_t, size_t> clusterCoverage;
	for (size_t i = 0; i < parent.size(); i++)
	{
		clusterCoverage[find(parent, i)] += 1;
	}
	phmap::flat_hash_map<size_t, size_t> clusterRenaming;
	std::vector<std::vector<size_t>> indicesInCluster;
	for (auto pair : clusterCoverage)
	{
		if (pair.second < minCoverage) continue;
		clusterRenaming[pair.first] = indicesInCluster.size();
		indicesInCluster.emplace_back();
	}
	for (size_t i = 0; i < parent.size(); i++)
	{
		if (clusterRenaming.count(find(parent, i)) == 0) continue;
		indicesInCluster[clusterRenaming.at(find(parent, i))].emplace_back(i);
	}
	std::vector<Edit> result;
	for (size_t i = 0; i < indicesInCluster.size(); i++)
	{
		result.emplace_back(edits[indicesInCluster[i][0]]);
		for (size_t j = 0; j < result.back().replacementKmerPositions.size(); j++)
		{
			double positionSum = 0;
			for (size_t k : indicesInCluster[i])
			{
				positionSum += edits[k].replacementKmerPositions[j];
			}
			positionSum /= indicesInCluster[i].size();
			assert(positionSum >= 0.0);
			assert(positionSum <= 1.0);
			result.back().replacementKmerPositions[j] = positionSum;
		}
	}
	return result;
}

std::vector<std::vector<Edit>> getValidEdits(const std::vector<std::vector<std::pair<size_t, size_t>>>& rawMinimizers, const std::vector<std::vector<Overlap>>& overlaps, const size_t numThreads)
{
	assert(rawMinimizers.size() == 2*overlaps.size());
	std::vector<std::vector<Edit>> result;
	result.resize(overlaps.size());
	iterateMultithreaded(0, overlaps.size(), numThreads, [&result, &rawMinimizers, &overlaps](const size_t i)
	{
		auto edits = getCoveredEdits(rawMinimizers, overlaps[i]);
		std::vector<Edit> actualEdits;
		std::vector<bool> coveredEdges; // points at sink node of edge
		coveredEdges.resize(rawMinimizers[i].size(), false);
		for (const auto& edit : edits)
		{
			assert(edit.endIndex > edit.startIndex);
			if (edit.endIndex == edit.startIndex+1 && edit.replacementKmers.size() == 2)
			{
				coveredEdges[edit.endIndex] = true;
				continue;
			}
			actualEdits.emplace_back(edit);
		}
		std::vector<bool> actualEditCanBeMade;
		actualEditCanBeMade.resize(actualEdits.size(), true);
		for (size_t l = 0; l < actualEdits.size(); l++)
		{
			for (size_t j = actualEdits[l].startIndex+1; j <= actualEdits[l].endIndex; j++)
			{
				if (coveredEdges[j]) actualEditCanBeMade[l] = false;
			}
			if (!actualEditCanBeMade[l]) continue;
			for (size_t j = 0; j < actualEdits.size(); j++)
			{
				if (j == l) continue;
				if (actualEdits[l].startIndex < actualEdits[j].endIndex && actualEdits[j].startIndex < actualEdits[l].endIndex)
				{
					actualEditCanBeMade[l] = false;
					actualEditCanBeMade[j] = false;
				}
			}
		}
		for (size_t j = actualEdits.size()-1; j < actualEdits.size(); j--)
		{
			if (actualEditCanBeMade[j]) continue;
			std::swap(actualEdits[j], actualEdits.back());
			actualEdits.pop_back();
		}
		std::sort(actualEdits.begin(), actualEdits.end(), [](const auto& left, const auto& right)
		{
			assert(left.startIndex != right.startIndex);
			return left.startIndex < right.startIndex;
		});
		result[i].insert(result[i].end(), actualEdits.begin(), actualEdits.end());
	});
	size_t countEdits = 0;
	for (size_t i = 0; i < result.size(); i++)
	{
		countEdits += result[i].size();
	}
	std::cerr << "count edits " << countEdits << std::endl;
	return result;
}

// return position of MIDDLE base in k-mer
// startPosition, endPosition also refer to middle base
size_t findKmerPosition(const std::string& sequence, const size_t startPosition, const size_t endPosition, const size_t kmer, const size_t kmerSize)
{
	const uint64_t mask = (1ull << (2ull*kmerSize)) - 1;
	assert(startPosition >= kmerSize/2);
	assert(endPosition > startPosition);
	assert(sequence.size() >= endPosition+kmerSize/2);
	uint64_t kmerHere = 0;
	for (size_t i = 0; i < kmerSize; i++)
	{
		kmerHere <<= 2;
		assert(startPosition+1+i-kmerSize/2 < sequence.size());
		kmerHere += charToInt(sequence[startPosition+1+i-kmerSize/2]);
	}
	size_t result = std::numeric_limits<size_t>::max();
	if (kmerHere == kmer) result = startPosition+1;
	for (size_t i = startPosition+1; i < endPosition; i++)
	{
		kmerHere <<= 2;
		kmerHere &= mask;
		assert(i+kmerSize/2 < sequence.size());
		kmerHere += charToInt(sequence[i+kmerSize/2]);
		if (kmerHere == kmer)
		{
			if (result == std::numeric_limits<size_t>::max())
			{
				result = i;
			}
			else
			{
				return std::numeric_limits<size_t>::max()-1;
			}
		}
	}
	return result;
}

size_t getUintKmerChar(size_t kmer, size_t posFromRightmost)
{
	return (kmer >> (2ull*posFromRightmost)) & 3ull;
}

size_t getFwEdits(const std::string& sequence, const size_t startPosition, const size_t kmer, const size_t kmerSize, const size_t maxEdits)
{
	std::vector<size_t> editColumn;
	editColumn.resize(kmerSize+1);
	for (size_t i = 0; i < editColumn.size(); i++)
	{
		editColumn[i] = i;
	}
	size_t minEdits = editColumn.back();
	assert(minEdits > 0);
	for (size_t row = 0; row < kmerSize+maxEdits && row+startPosition < sequence.size(); row++)
	{
		std::vector<size_t> newColumn;
		newColumn.resize(editColumn.size(), kmerSize+1);
		newColumn[0] = editColumn[0]+1;
		for (size_t i = 1; i < editColumn.size(); i++)
		{
			newColumn[i] = std::min(newColumn[i-1]+1, editColumn[i]+1);
			if (charToInt(sequence[startPosition+row]) == getUintKmerChar(kmer, kmerSize-i))
			{
				newColumn[i] = std::min(newColumn[i], editColumn[i-1]);
			}
			else
			{
				assert(row != 0 || i != 1);
				newColumn[i] = std::min(newColumn[i], editColumn[i-1]+1);
			}
		}
		std::swap(editColumn, newColumn);
		minEdits = std::min(minEdits, editColumn.back());
		assert(row+1 >= kmerSize || minEdits > 0);
	}
	return minEdits;
}

size_t getBwEdits(const std::string& sequence, const size_t startPosition, const size_t kmer, const size_t kmerSize, const size_t maxEdits)
{
	std::vector<size_t> editColumn;
	editColumn.resize(kmerSize+1);
	for (size_t i = 0; i < editColumn.size(); i++)
	{
		editColumn[i] = i;
	}
	size_t minEdits = editColumn.back();
	assert(minEdits > 0);
	for (size_t row = 0; row < kmerSize+maxEdits && row < startPosition; row++)
	{
		std::vector<size_t> newColumn;
		newColumn.resize(editColumn.size(), kmerSize+1);
		newColumn[0] = editColumn[0]+1;
		for (size_t i = 1; i < editColumn.size(); i++)
		{
			newColumn[i] = std::min(newColumn[i-1]+1, editColumn[i]+1);
			if (charToInt(sequence[startPosition-row]) == getUintKmerChar(kmer, i-1))
			{
				newColumn[i] = std::min(newColumn[i], editColumn[i-1]);
			}
			else
			{
				assert(row != 0 || i != 1);
				newColumn[i] = std::min(newColumn[i], editColumn[i-1]+1);
			}
		}
		std::swap(editColumn, newColumn);
		minEdits = std::min(minEdits, editColumn.back());
		assert(row+1 >= kmerSize || minEdits > 0);
	}
	return minEdits;
}

// return position of MIDDLE base in k-mer
// startPosition, endPosition also refer to middle base
size_t findKmerPositionWithEdits(const std::string& sequence, const size_t startPosition, const size_t endPosition, const size_t kmer, const size_t kmerSize, const size_t maxEditsTotal, const size_t maxEditsPerSide)
{
	assert(kmerSize % 2 == 1);
	const size_t halfmerSize = (kmerSize+1)/2;
	const size_t halfmerMask = (1ull << (2ull*halfmerSize))-1ull;
	const size_t forwardHalfmer = kmer & halfmerMask;
	const size_t backwardHalfmer = (kmer >> 2ull*(kmerSize-halfmerSize));
	assert((backwardHalfmer & halfmerMask) == backwardHalfmer);
	assert(getUintKmerChar(forwardHalfmer, halfmerSize-1) == getUintKmerChar(backwardHalfmer, 0));
	size_t result = std::numeric_limits<size_t>::max();
	for (size_t midPos = startPosition+1; midPos < endPosition; midPos++)
	{
		if (charToInt(sequence[midPos]) != (backwardHalfmer & 3)) continue; // middle base must always match regardless of edits
		size_t fwEdits = getFwEdits(sequence, midPos, forwardHalfmer, halfmerSize, maxEditsPerSide);
		if (fwEdits > maxEditsPerSide) continue;
		size_t bwEdits = getBwEdits(sequence, midPos, backwardHalfmer, halfmerSize, maxEditsPerSide);
		if (bwEdits > maxEditsPerSide) continue;
		if (fwEdits+bwEdits > maxEditsTotal) continue;
		if (result == std::numeric_limits<size_t>::max())
		{
			result = midPos;
		}
		else
		{
			result = std::numeric_limits<size_t>::max()-1;
		}
	}
	return result;
}

std::vector<size_t> findKmerPositions(const std::string& sequence, const size_t startPosition, const size_t endPosition, const std::vector<size_t>& replacementKmers, const std::vector<double>& replacementKmerPositions, const size_t kmerSize)
{
	assert(endPosition > startPosition);
	std::vector<size_t> result;
	result.resize(replacementKmers.size(), std::numeric_limits<size_t>::max());
	result[0] = startPosition;
	result.back() = endPosition;
	size_t acceptableDistance = (endPosition-startPosition) * 0.1 + 10;
	for (size_t j = 1; j+1 < replacementKmers.size(); j++)
	{
		assert(replacementKmerPositions[j] > 0);
		assert(replacementKmerPositions[j] < 1.0);
		size_t checkStart = startPosition + replacementKmerPositions[j] * (endPosition-startPosition);
		if (checkStart > startPosition+acceptableDistance+1)
		{
			checkStart -= acceptableDistance;
		}
		else
		{
			checkStart = startPosition+1;
		}
		size_t checkEnd = startPosition + replacementKmerPositions[j] * (endPosition-startPosition);
		if (checkEnd + acceptableDistance < endPosition)
		{
			checkEnd += acceptableDistance;
		}
		else
		{
			checkEnd = endPosition-1;
		}
		size_t kmerPos = findKmerPosition(sequence, checkStart, checkEnd, replacementKmers[j], kmerSize);
		if (kmerPos == std::numeric_limits<size_t>::max())
		{
			kmerPos = findKmerPositionWithEdits(sequence, checkStart, checkEnd, replacementKmers[j], kmerSize, 1, 1);
		}
		if (kmerPos == std::numeric_limits<size_t>::max())
		{
			kmerPos = findKmerPositionWithEdits(sequence, checkStart, checkEnd, replacementKmers[j], kmerSize, 2, 2);
		}
		result[j] = kmerPos;
	}
	return result;
}

void filterOutEditsWhereKmerIsNotFound(std::vector<std::vector<Edit>>& validEdits, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<std::pair<size_t, size_t>>>& rawMinimizers, const size_t kmerSize, const size_t numThreads)
{
	assert(rawMinimizers.size() == sequenceIndex.size()*2);
	assert(validEdits.size() == sequenceIndex.size());
	iterateMultithreaded(0, validEdits.size(), numThreads, [&validEdits, &sequenceIndex, &rawMinimizers, kmerSize](const size_t i)
	{
		if (validEdits[i].size() == 0) return;
		std::string sequence = sequenceIndex.getSequence(i);
		for (size_t j = validEdits[i].size()-1; j < validEdits[i].size(); j--)
		{
			validEdits[i][j].actualKmerPositionsInReference = findKmerPositions(sequence, rawMinimizers[i][validEdits[i][j].startIndex].first, rawMinimizers[i][validEdits[i][j].endIndex].first, validEdits[i][j].replacementKmers, validEdits[i][j].replacementKmerPositions, kmerSize);
			assert(validEdits[i][j].actualKmerPositionsInReference.size() >= 2);
			assert(validEdits[i][j].actualKmerPositionsInReference.size() == validEdits[i][j].replacementKmers.size());
			assert(validEdits[i][j].actualKmerPositionsInReference[0] == rawMinimizers[i][validEdits[i][j].startIndex].first);
			assert(validEdits[i][j].actualKmerPositionsInReference.back() == rawMinimizers[i][validEdits[i][j].endIndex].first);
			bool allGood = true;
			for (size_t k = 0; k < validEdits[i][j].actualKmerPositionsInReference.size(); k++)
			{
				if (validEdits[i][j].actualKmerPositionsInReference[k] == std::numeric_limits<size_t>::max()) allGood = false;
			}
			for (size_t k = 1; k < validEdits[i][j].actualKmerPositionsInReference.size(); k++)
			{
				if (validEdits[i][j].actualKmerPositionsInReference[k] <= validEdits[i][j].actualKmerPositionsInReference[k-1]) allGood = false;
			}
			if (!allGood)
			{
				std::string info = "invalid edit read " + std::to_string(i) + " poses";
				for (size_t pos = validEdits[i][j].startIndex; pos <= validEdits[i][j].endIndex; pos++)
				{
					info += " " + std::to_string(rawMinimizers[i][pos].first);
				}
				info += " vs";
				for (size_t pos : validEdits[i][j].actualKmerPositionsInReference)
				{
					info += " " + std::to_string(pos);
				}
				info += "\n";
				std::cerr << info;
				std::swap(validEdits[i][j], validEdits[i].back());
				validEdits[i].pop_back();
			}
		}
		std::sort(validEdits[i].begin(), validEdits[i].end());
	});
}

std::vector<std::vector<std::pair<size_t, size_t>>> getCorrections(const std::vector<std::vector<std::pair<size_t, size_t>>>& rawMinimizers, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<std::vector<Overlap>>& overlaps, const size_t kmerSize, const size_t numThreads)
{
	assert(rawMinimizers.size() == 2*overlaps.size());
	std::vector<std::vector<Edit>> validEdits = getValidEdits(rawMinimizers, overlaps, numThreads);
	filterOutEditsWhereKmerIsNotFound(validEdits, sequenceIndex, rawMinimizers, kmerSize, numThreads);
	size_t countFilteredEdits = 0;
	for (size_t i = 0; i < validEdits.size(); i++)
	{
		countFilteredEdits += validEdits[i].size();
	}
	std::cerr << "count filtered edits " << countFilteredEdits << std::endl;
	std::vector<std::vector<std::pair<size_t, size_t>>> result;
	result.resize(sequenceIndex.size());
	iterateMultithreaded(0, result.size(), numThreads, [&result, &rawMinimizers, &validEdits](const size_t i)
	{
		for (size_t j = 1; j < validEdits[i].size(); j++)
		{
			assert(validEdits[i][j].startIndex > validEdits[i][j-1].startIndex);
			assert(validEdits[i][j].startIndex >= validEdits[i][j-1].endIndex);
		}
		size_t lastMatch = 0;
		for (size_t j = 0; j < validEdits[i].size(); j++)
		{
			result[i].insert(result[i].end(), rawMinimizers[i].begin()+lastMatch, rawMinimizers[i].begin()+validEdits[i][j].startIndex+1);
			assert(validEdits[i][j].actualKmerPositionsInReference.size() == validEdits[i][j].replacementKmers.size());
			assert(validEdits[i][j].actualKmerPositionsInReference[0] == rawMinimizers[i][validEdits[i][j].startIndex].first);
			assert(validEdits[i][j].actualKmerPositionsInReference.back() == rawMinimizers[i][validEdits[i][j].endIndex].first);
			for (size_t l = 1; l+1 < validEdits[i][j].replacementKmers.size(); l++)
			{
				result[i].emplace_back(validEdits[i][j].actualKmerPositionsInReference[l], validEdits[i][j].replacementKmers[l]);
			}
			std::string info = "correct read " + std::to_string(i) + " minimizer poses";
			for (size_t l = validEdits[i][j].startIndex; l <= validEdits[i][j].endIndex; l++)
			{
				info += " " + std::to_string(rawMinimizers[i][l].first);
			}
			info += " to";
			for (size_t l = 0; l < validEdits[i][j].replacementKmers.size(); l++)
			{
				info += " " + std::to_string(validEdits[i][j].actualKmerPositionsInReference[l]);
			}
			info += "\n";
			std::cerr << info;
			lastMatch = validEdits[i][j].endIndex;
		}
		result[i].insert(result[i].end(), rawMinimizers[i].begin()+lastMatch, rawMinimizers[i].end());
		for (size_t j = 1; j < result[i].size(); j++)
		{
			assert(result[i][j].first > result[i][j-1].first);
		}
	});
	return result;
}

std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>> getCorrectedMinimizers(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t kmerSize, const size_t windowSize, const size_t numThreads)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> rawMinimizers;
	rawMinimizers.resize(sequenceIndex.size());
	std::cerr << "get raw minimizers" << std::endl;
	iterateMultithreaded(0, sequenceIndex.size(), numThreads, [&rawMinimizers, &sequenceIndex, kmerSize, windowSize](const size_t i)
	{
		std::string seq = sequenceIndex.getSequence(i);
		rawMinimizers[i] = getRawMinimizers(seq, kmerSize, windowSize);
	});
	makeReverseReads(rawMinimizers, rawReadLengths, kmerSize);
	std::cerr << "get overlaps" << std::endl;
	auto overlaps = getOverlaps(rawMinimizers, numThreads);
	std::cerr << "get corrected minimizers" << std::endl;
	auto correctedMinimizers = getCorrections(rawMinimizers, sequenceIndex, overlaps, kmerSize, numThreads);
	assert(correctedMinimizers.size() == sequenceIndex.size());
	std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>> result;
	result.first.resize(correctedMinimizers.size());
	result.second.resize(correctedMinimizers.size());
	for (size_t i = 0; i < correctedMinimizers.size(); i++)
	{
		for (auto pair : correctedMinimizers[i])
		{
			result.first[i].emplace_back(pair.first);
			result.second[i].emplace_back(pair.second);
		}
	}
	return result;
}

std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>> getChunksPerReadFromCorrectedMinimizers(const std::vector<std::vector<size_t>>& correctedMinimizers, const std::vector<std::vector<size_t>>& minimizerKmers, const std::vector<size_t>& rawReadLengths, const size_t kmerSize)
{
	assert(correctedMinimizers.size()*2 == rawReadLengths.size());
	assert(minimizerKmers.size()*2 == rawReadLengths.size());
	std::vector<std::vector<size_t>> resultMinimizers;
	resultMinimizers = correctedMinimizers;
	addReverseMinimizers(resultMinimizers, rawReadLengths);
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> result;
	result.resize(correctedMinimizers.size()*2);
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> minimizerPairToChunk;
	for (size_t i = 0; i < correctedMinimizers.size(); i++)
	{
		for (size_t j = 1; j < correctedMinimizers[i].size(); j++)
		{
			std::pair<size_t, size_t> key { minimizerKmers[i][j-1], minimizerKmers[i][j] };
			if (minimizerPairToChunk.count(key) == 0)
			{
				size_t chunk = minimizerPairToChunk.size();
				minimizerPairToChunk[key] = chunk;
			}
			size_t chunk = minimizerPairToChunk.at(key);
			result[i].emplace_back(correctedMinimizers[i][j-1], correctedMinimizers[i][j], chunk + firstBitUint64_t);
			std::pair<size_t, size_t> reverseKey  { reverseKmer(key.second, kmerSize), reverseKmer(key.first, kmerSize) };
			if (minimizerPairToChunk.count(reverseKey) == 0)
			{
				size_t reverseChunk = minimizerPairToChunk.size();
				minimizerPairToChunk[reverseKey] = reverseChunk;
			}
			size_t reverseChunk = minimizerPairToChunk.at(reverseKey);
			result[i+correctedMinimizers.size()].emplace_back(rawReadLengths[i] - 1 - correctedMinimizers[i][j], rawReadLengths[i] - 1 - correctedMinimizers[i][j-1], reverseChunk + firstBitUint64_t);
		}
		std::reverse(result[i+correctedMinimizers.size()].begin(), result[i+correctedMinimizers.size()].end());
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			assert(std::get<1>(result[i][j]) > std::get<0>(result[i][j]));
			assert(j == 0 || std::get<0>(result[i][j]) > std::get<0>(result[i][j-1]));
			assert(j == 0 || std::get<1>(result[i][j]) > std::get<1>(result[i][j-1]));
		}
	}
	return std::make_pair(resultMinimizers, result);
}
