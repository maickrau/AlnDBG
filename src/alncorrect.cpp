#include <mutex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
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

std::pair<UnitigGraph, std::vector<ReadPathBundle>> makeGraph(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t minCoverage, const size_t numThreads)
{
	KmerGraph kmerGraph;
	std::vector<ReadPathBundle> readKmerGraphPaths;
	std::tie(kmerGraph, readKmerGraphPaths) = makeKmerGraph(readLengths, matches, minCoverage, numThreads);
	UnitigGraph unitigGraph;
	std::vector<ReadPathBundle> readUnitigGraphPaths;
	std::tie(unitigGraph, readUnitigGraphPaths) = makeUnitigGraph(kmerGraph, readKmerGraphPaths, minCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 10);
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 10);
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 10);
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 10);
	return std::make_pair(std::move(unitigGraph), std::move(readUnitigGraphPaths));
}

std::vector<MatchGroup> getMatches(const MatchIndex& matchIndex, const std::vector<TwobitString>& readSequences, size_t numThreads, const std::vector<size_t>& rawReadLengths, const size_t minAlignmentLength, const size_t k, const size_t graphk, const size_t graphd, const std::vector<bool>& usableChunkmers, const size_t maxIndexCoverage)
{
	std::mutex printMutex;
	std::vector<MatchGroup> matches;
	auto result = matchIndex.iterateMatchChains(numThreads, 2, maxIndexCoverage, 50, rawReadLengths, usableChunkmers, [&printMutex, &matches, minAlignmentLength, &rawReadLengths, k, graphd](const size_t left, const size_t leftstart, const size_t leftend, const bool leftFw, const size_t right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		if (leftend+1-leftstart < minAlignmentLength) return;
		if (rightend+1-rightstart < minAlignmentLength) return;
		std::lock_guard<std::mutex> lock { printMutex };
		assert(leftFw);
		matches.emplace_back();
		matches.back().leftRead = left;
		matches.back().rightRead = right;
		matches.back().rightFw = rightFw;
		matches.back().leftStart = leftstart;
		matches.back().leftEnd = leftend;
		matches.back().rightStart = rightstart;
		matches.back().rightEnd = rightend;
		double slope = (double)(rightend-rightstart) / (double)(leftend-leftstart);
		if (rightstart*slope < leftstart)
		{
			matches.back().leftStart -= rightstart*slope;
			matches.back().rightStart = 0;
		}
		else
		{
			assert(rightstart >= (size_t)(leftstart/slope));
			matches.back().rightStart -= leftstart/slope;
			matches.back().leftStart = 0;
		}
		if ((rawReadLengths[right]-rightend-1)*slope < (rawReadLengths[left]-leftend-1))
		{
			matches.back().leftEnd += (rawReadLengths[right]-rightend-1)*slope;
			assert(matches.back().leftEnd < rawReadLengths[left]);
			matches.back().rightEnd = rawReadLengths[right]-1;
		}
		else
		{
			matches.back().rightEnd += (rawReadLengths[left]-leftend-1)/slope;
			assert(matches.back().rightEnd < rawReadLengths[right]);
			matches.back().leftEnd = rawReadLengths[left]-1;
		}
	});
	std::stable_sort(matches.begin(), matches.end(), [](const auto& left, const auto& right){
		if (left.leftRead < right.leftRead) return true;
		if (left.leftRead > right.leftRead) return false;
		if (left.leftStart < right.leftStart) return true;
		if (left.leftStart > right.leftStart) return false;
		if (left.leftEnd > right.leftEnd) return true;
		return false;
	});
	std::cerr << result.totalReadChunkMatches << " read-windowchunk matches (except unique)" << std::endl;
	std::cerr << result.readsWithMatch << " reads with a match" << std::endl;
	std::cerr << result.readPairMatches << " read-read matches" << std::endl;
	std::cerr << result.readChainMatches << " chain matches" << std::endl;
	std::cerr << result.totalMatches << " window matches" << std::endl;
	std::cerr << result.maxPerChunk << " max windowchunk size" << std::endl;
	std::cerr << matches.size() << " mapping matches" << std::endl;
	matches = addKmerMatches(numThreads, readSequences, matches, graphk, graphd);
	std::cerr << matches.size() << " mapping matches" << std::endl;
	return matches;
}

phmap::flat_hash_set<size_t> getNonrepetitiveNodes(const UnitigGraph& graph, const std::vector<ReadPathBundle>& paths)
{
	std::vector<bool> repetitive;
	repetitive.resize(graph.nodeCount(), false);
	for (size_t i = 0; i < paths.size(); i++)
	{
		phmap::flat_hash_set<size_t> seenHere;
		for (size_t j = 0; j < paths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < paths[i].paths[j].path.size(); k++)
			{
				uint64_t node = paths[i].paths[j].path[k];
				if (j > 0 && k == 0 && seenHere.count(node & maskUint64_t) == 1)
				{
					if (paths[i].paths[j-1].path.back() == paths[i].paths[j].path[k])
					{
						bool couldBeSame = true;
						size_t thisNode = paths[i].paths[j].path[k] & maskUint64_t;
						size_t prevNodeEnd = graph.lengths[thisNode] - paths[i].paths[j-1].pathRightClipKmers;
						size_t currNodeStart = paths[i].paths[j].pathLeftClipKmers;
						if (currNodeStart <= prevNodeEnd) couldBeSame = false;
						if (couldBeSame)
						{
							size_t currReadStart = paths[i].paths[j].readStartPos;
							size_t prevReadEnd = paths[i].paths[j-1].readStartPos;
							for (auto node : paths[i].paths[j-1].path)
							{
								prevReadEnd += graph.lengths[node & maskUint64_t];
							}
							assert(prevReadEnd >= paths[i].paths[j-1].pathLeftClipKmers + paths[i].paths[j-1].pathRightClipKmers);
							prevReadEnd -= paths[i].paths[j-1].pathLeftClipKmers;
							prevReadEnd -= paths[i].paths[j-1].pathRightClipKmers;
							assert(prevReadEnd <= currReadStart);
							size_t readDistance = currReadStart - prevReadEnd;
							size_t nodeDistance = currNodeStart - prevNodeEnd;
							if (readDistance > nodeDistance + 50) couldBeSame = false;
							if (nodeDistance > readDistance + 50) couldBeSame =	false;
							if (couldBeSame) continue;
						}
					}
				}
				if (seenHere.count(node & maskUint64_t) == 1) repetitive[node & maskUint64_t] = true;
				seenHere.insert(node & maskUint64_t);
			}
		}
	}
	phmap::flat_hash_set<size_t> result;
	for (size_t i = 0; i < repetitive.size(); i++)
	{
		if (!repetitive[i]) result.emplace(i);
	}
	return result;
}

std::vector<std::pair<uint64_t, uint64_t>> getPhasedAlleles(const std::vector<std::pair<uint64_t, size_t>>& bw, const std::vector<std::pair<uint64_t, size_t>>& fw)
{
	phmap::flat_hash_map<uint64_t, std::vector<size_t>> bwReadsPerAllele;
	phmap::flat_hash_map<uint64_t, std::vector<size_t>> fwReadsPerAllele;
	for (auto pair : bw)
	{
		bwReadsPerAllele[pair.first].emplace_back(pair.second);
	}
	for (auto pair : fw)
	{
		fwReadsPerAllele[pair.first].emplace_back(pair.second);
	}
	for (auto& pair : bwReadsPerAllele)
	{
		std::sort(pair.second.begin(), pair.second.end());
	}
	for (auto& pair : fwReadsPerAllele)
	{
		std::sort(pair.second.begin(), pair.second.end());
	}
	if (bwReadsPerAllele.size() != fwReadsPerAllele.size()) return std::vector<std::pair<uint64_t, uint64_t>>{};
	if (bwReadsPerAllele.size() < 2) return std::vector<std::pair<uint64_t, uint64_t>>{};
	std::vector<std::pair<uint64_t, uint64_t>> result;
	phmap::flat_hash_set<uint64_t> foundBw;
	phmap::flat_hash_set<uint64_t> foundFw;
	for (const auto& pair : bwReadsPerAllele)
	{
		for (const auto& pair2 : fwReadsPerAllele)
		{
			size_t matching = intersectSize(pair.second, pair2.second);
			if (matching == 0) continue;
			if (matching < 4) return std::vector<std::pair<uint64_t, uint64_t>>{};
			if (matching < pair.second.size() * 0.5) return std::vector<std::pair<uint64_t, uint64_t>>{};
			if (matching < pair2.second.size() * 0.5) return std::vector<std::pair<uint64_t, uint64_t>>{};
			if (foundBw.count(pair.first) == 1) return std::vector<std::pair<uint64_t, uint64_t>>{};
			if (foundFw.count(pair2.first) == 1) return std::vector<std::pair<uint64_t, uint64_t>>{};
			result.emplace_back(pair.first, pair2.first);
			foundBw.emplace(pair.first);
			foundFw.emplace(pair2.first);
		}
	}
	if (result.size() != bwReadsPerAllele.size()) return std::vector<std::pair<uint64_t, uint64_t>>{};
	return result;
}

std::pair<TwobitString, size_t> getCorrectedSequence(const std::vector<size_t>& readKmerLengths, const std::vector<TwobitString>& readSequences, const std::vector<MatchGroup>& matches, const size_t graphk)
{
	UnitigGraph graph;
	std::vector<ReadPathBundle> paths;
	std::tie(graph, paths) = makeGraph(readKmerLengths, matches, 2, 1);
	phmap::flat_hash_set<size_t> uniqueNodes = getNonrepetitiveNodes(graph, paths);
	phmap::flat_hash_map<uint64_t, size_t> nodeIndexInRead;
	for (size_t j = 0; j < paths[0].paths.size(); j++)
	{
		for (size_t k = 0; k < paths[0].paths[j].path.size(); k++)
		{
			uint64_t node = paths[0].paths[j].path[k];
			if (uniqueNodes.count(node & maskUint64_t) == 0) continue;
			if (nodeIndexInRead.count(node) == 1) continue;
			if (nodeIndexInRead.count(node ^ firstBitUint64_t) == 1) continue;
			size_t index = nodeIndexInRead.size();
			nodeIndexInRead[node] = index;
		}
	}
	std::vector<std::vector<std::pair<uint64_t, size_t>>> bwForkAlleles;
	std::vector<std::vector<std::pair<uint64_t, size_t>>> fwForkAlleles;
	bwForkAlleles.resize(nodeIndexInRead.size());
	fwForkAlleles.resize(nodeIndexInRead.size());
	for (size_t i = 0; i < paths.size(); i++)
	{
		for (size_t j = 0; j < paths[i].paths.size(); j++)
		{
			for (size_t k = 1; k < paths[i].paths[j].path.size(); k++)
			{
				if (nodeIndexInRead.count(paths[i].paths[j].path[k-1]) == 1)
				{
					size_t index = nodeIndexInRead.at(paths[i].paths[j].path[k-1]);
					fwForkAlleles[index].emplace_back(paths[i].paths[j].path[k], i);
				}
				if (nodeIndexInRead.count(paths[i].paths[j].path[k-1] ^ firstBitUint64_t) == 1)
				{
					size_t index = nodeIndexInRead.at(paths[i].paths[j].path[k-1] ^ firstBitUint64_t);
					bwForkAlleles[index].emplace_back(paths[i].paths[j].path[k] ^ firstBitUint64_t, i);
				}
				if (nodeIndexInRead.count(paths[i].paths[j].path[k]) == 1)
				{
					size_t index = nodeIndexInRead.at(paths[i].paths[j].path[k]);
					bwForkAlleles[index].emplace_back(paths[i].paths[j].path[k-1], i);
				}
				if (nodeIndexInRead.count(paths[i].paths[j].path[k] ^ firstBitUint64_t) == 1)
				{
					size_t index = nodeIndexInRead.at(paths[i].paths[j].path[k] ^ firstBitUint64_t);
					fwForkAlleles[index].emplace_back(paths[i].paths[j].path[k-1] ^ firstBitUint64_t, i);
				}
			}
		}
	}
	phmap::flat_hash_set<size_t> forbiddenReads;
	for (size_t i = 0; i < bwForkAlleles.size(); i++)
	{
		for (size_t j = i; j < fwForkAlleles.size(); j++)
		{
			auto phasedAlleles = getPhasedAlleles(bwForkAlleles[i], fwForkAlleles[j]);
			if (phasedAlleles.size() < 2) continue;
			uint64_t bwNode = std::numeric_limits<size_t>::max();
			uint64_t fwNode = std::numeric_limits<size_t>::max();
			bool valid = true;
			for (auto pair : bwForkAlleles[i])
			{
				if (pair.second != 0) continue;
				if (bwNode != std::numeric_limits<size_t>::max()) valid = false;
				bwNode = pair.first;
			}
			for (auto pair : fwForkAlleles[j])
			{
				if (pair.second != 0) continue;
				if (fwNode != std::numeric_limits<size_t>::max()) valid = false;
				fwNode = pair.first;
			}
			if (!valid) continue;
			size_t bwAllele = std::numeric_limits<size_t>::max();
			size_t fwAllele = std::numeric_limits<size_t>::max();
			for (size_t i = 0; i < phasedAlleles.size(); i++)
			{
				if (phasedAlleles[i].first == bwNode)
				{
					assert(bwAllele == std::numeric_limits<size_t>::max());
					bwAllele = i;
				}
				if (phasedAlleles[i].second == fwNode)
				{
					assert(fwAllele == std::numeric_limits<size_t>::max());
					fwAllele = i;
				}
			}
			if (bwAllele == std::numeric_limits<size_t>::max() && fwAllele == std::numeric_limits<size_t>::max()) continue;
			if (bwAllele != std::numeric_limits<size_t>::max() && fwAllele != std::numeric_limits<size_t>::max() && bwAllele != fwAllele) continue;
			size_t allele = bwAllele;
			if (bwAllele == std::numeric_limits<size_t>::max())
			{
				assert(fwAllele != std::numeric_limits<size_t>::max());
				allele = fwAllele;
			}
			assert(bwNode == phasedAlleles[allele].first || bwNode == std::numeric_limits<size_t>::max());
			assert(fwNode == phasedAlleles[allele].second || fwNode == std::numeric_limits<size_t>::max());
			bwNode = phasedAlleles[allele].first;
			fwNode = phasedAlleles[allele].second;
			for (auto pair : bwForkAlleles[i])
			{
				if (pair.first == bwNode) continue;
				assert(pair.second != 0);
				forbiddenReads.emplace(pair.second);
			}
			for (auto pair : fwForkAlleles[j])
			{
				if (pair.first == fwNode) continue;
				assert(pair.second != 0);
				forbiddenReads.emplace(pair.second);
			}
		}
	}
	if (forbiddenReads.size() >= 1)
	{
		std::vector<size_t> subsetReadKmerLengths;
		std::vector<TwobitString> subsetReadSequences;
		std::vector<MatchGroup> subsetMatches;
		phmap::flat_hash_map<size_t, size_t> newIndex;
		for (size_t i = 0; i < readKmerLengths.size(); i++)
		{
			if (forbiddenReads.count(i) == 1) continue;
			size_t index = newIndex.size();
			newIndex[i] = index;
			subsetReadKmerLengths.emplace_back(readKmerLengths[i]);
			subsetReadSequences.emplace_back(readSequences[i]);
		}
		for (size_t i = 0; i < matches.size(); i++)
		{
			if (forbiddenReads.count(matches[i].leftRead) == 1) continue;
			if (forbiddenReads.count(matches[i].rightRead) == 1) continue;
			subsetMatches.emplace_back(matches[i]);
			subsetMatches.back().leftRead = newIndex.at(subsetMatches.back().leftRead);
			subsetMatches.back().rightRead = newIndex.at(subsetMatches.back().rightRead);
		}
		return getCorrectedSequence(subsetReadKmerLengths, subsetReadSequences, subsetMatches, graphk);
	}
	std::vector<std::tuple<size_t, size_t, uint64_t, size_t, size_t>> corrections;
	for (size_t j = 1; j < paths[0].paths.size(); j++)
	{
		if (paths[0].paths[j-1].path.back() != paths[0].paths[j].path[0]) continue;
		size_t thisNode = paths[0].paths[j].path[0] & maskUint64_t;
		size_t prevNodeEnd = graph.lengths[thisNode] - paths[0].paths[j-1].pathRightClipKmers;
		size_t currNodeStart = paths[0].paths[j].pathLeftClipKmers;
		if (currNodeStart <= prevNodeEnd) continue;
		size_t currReadStart = paths[0].paths[j].readStartPos;
		size_t prevReadEnd = paths[0].paths[j-1].readStartPos;
		for (auto node : paths[0].paths[j-1].path)
		{
			prevReadEnd += graph.lengths[node & maskUint64_t];
		}
		assert(prevReadEnd >= paths[0].paths[j-1].pathLeftClipKmers + paths[0].paths[j-1].pathRightClipKmers);
		prevReadEnd -= paths[0].paths[j-1].pathLeftClipKmers;
		prevReadEnd -= paths[0].paths[j-1].pathRightClipKmers;
		assert(prevReadEnd <= currReadStart);
		size_t readDistance = currReadStart - prevReadEnd;
		size_t nodeDistance = currNodeStart - prevNodeEnd;
		if (readDistance > nodeDistance + 50) continue;
		if (nodeDistance > readDistance + 50) continue;
		corrections.emplace_back(prevReadEnd, currReadStart, paths[0].paths[j].path[0], prevNodeEnd, currNodeStart);
	}
	if (corrections.size() == 0) return std::make_pair(readSequences[0], 0);
	auto nodeSequences = getNodeSequences(graph, paths, graphk, readSequences);
	std::vector<std::tuple<size_t, size_t, std::string>> chosenReplacements;
	for (auto t : corrections)
	{
		std::string replacementSequence;
		assert(std::get<0>(t) <= std::get<1>(t));
		assert(std::get<0>(t) < readSequences[0].size()-graphk+1);
		assert(std::get<1>(t) <= readSequences[0].size()-graphk+1);
		assert(std::get<3>(t) <= std::get<4>(t));
		assert(std::get<3>(t) < nodeSequences[std::get<2>(t) & maskUint64_t].size()-graphk+1);
		assert(std::get<4>(t) <= nodeSequences[std::get<2>(t) & maskUint64_t].size()-graphk+1);
		if (std::get<2>(t) & firstBitUint64_t)
		{
			replacementSequence = nodeSequences[std::get<2>(t) & maskUint64_t].substr(std::get<3>(t), std::get<4>(t)-std::get<3>(t)+graphk-1);
		}
		else
		{
			size_t size = nodeSequences[std::get<2>(t) & maskUint64_t].size();
			replacementSequence = nodeSequences[std::get<2>(t) & maskUint64_t].substr(size-(std::get<4>(t)+graphk-1), std::get<4>(t)-std::get<3>(t)+graphk-1);
			replacementSequence = CommonUtils::ReverseComplement(replacementSequence);
		}
		assert(replacementSequence.size() >= graphk-1);
		if (replacementSequence.substr(0, graphk-2) != readSequences[0].substr(std::get<0>(t), graphk-2)) continue;
		if (replacementSequence.substr(replacementSequence.size()-(graphk-1)) != readSequences[0].substr(std::get<1>(t), graphk-1)) continue;
		chosenReplacements.emplace_back(std::get<0>(t), std::get<1>(t), replacementSequence.substr(0, replacementSequence.size()-graphk+1));
		assert(std::get<1>(chosenReplacements.back()) >= std::get<0>(chosenReplacements.back()));
//		std::cerr << "correct " << readSequences[0].substr(std::get<0>(chosenReplacements.back()), std::get<1>(chosenReplacements.back()) - std::get<0>(chosenReplacements.back())) << " with " << std::get<2>(chosenReplacements.back()) << " (" << ((std::get<2>(t) & firstBitUint64_t) ? "fw" : "bw") << ")" << std::endl;
	}
	if (chosenReplacements.size() == 0) return std::make_pair(readSequences[0], 0);
	std::sort(chosenReplacements.begin(), chosenReplacements.end());
	TwobitString result;
	size_t lastIncluded = 1;
	assert(std::get<0>(chosenReplacements[0]) >= 1);
	result.resize(1);
	result.set(0, readSequences[0].get(0));
	for (size_t i = 0; i < chosenReplacements.size(); i++)
	{
		assert(i == 0 || std::get<0>(chosenReplacements[i]) >= std::get<1>(chosenReplacements[i-1]));
		assert(std::get<0>(chosenReplacements[i]) >= lastIncluded);
		for (size_t j = lastIncluded; j < std::get<0>(chosenReplacements[i]); j++)
		{
			result.emplace_back(readSequences[0].get(j));
		}
		for (size_t j = 0; j < std::get<2>(chosenReplacements[i]).size(); j++)
		{
			switch(std::get<2>(chosenReplacements[i])[j])
			{
			case 'A':
				result.emplace_back(0);
				break;
			case 'C':
				result.emplace_back(1);
				break;
			case 'G':
				result.emplace_back(2);
				break;
			case 'T':
				result.emplace_back(3);
				break;
			default:
				assert(false);
			}
		}
		lastIncluded = std::get<1>(chosenReplacements[i]);
	}
	for (size_t i = lastIncluded; i < readSequences[0].size(); i++)
	{
		result.emplace_back(readSequences[0].get(i));
	}
	return std::make_pair(result, chosenReplacements.size());
}

void correctReads(const std::vector<size_t>& readKmerLengths, const std::vector<std::string>& readNames, const std::vector<TwobitString>& readSequences, const std::vector<MatchGroup>& matches, const size_t graphk, const size_t numThreads)
{
	assert(readKmerLengths.size() == readNames.size());
	assert(readKmerLengths.size() == readSequences.size());
	std::vector<std::vector<size_t>> alignmentsPerRead;
	alignmentsPerRead.resize(readKmerLengths.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		assert(matches[i].leftRead != matches[i].rightRead);
		alignmentsPerRead[matches[i].leftRead].emplace_back(i);
		alignmentsPerRead[matches[i].rightRead].emplace_back(i);
	}
	std::mutex printMutex;
	std::mutex indexMutex;
	size_t nextIndex = 0;
	std::vector<std::thread> threads;
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([&alignmentsPerRead, &printMutex, &indexMutex, &nextIndex, &readKmerLengths, &readNames, &readSequences, &matches, graphk]()
		{
			while (true)
			{
				size_t readi = std::numeric_limits<size_t>::max();
				{
					std::lock_guard<std::mutex> lock { indexMutex };
					readi = nextIndex;
					nextIndex += 1;
				}
				if (readi >= readKmerLengths.size()) return;
				phmap::flat_hash_map<size_t, size_t> relevantReads;
				relevantReads[readi] = 0;
				if (alignmentsPerRead[readi].size() == 0)
				{
					std::string sequence = readSequences[readi].toString();
					std::lock_guard<std::mutex> lock { printMutex };
					std::cout << ">" << readNames[readi] << std::endl;
					std::cout << sequence << std::endl;
					std::cerr << "read " << readNames[readi] << " had 0 corrections (no overlaps)" << std::endl;
					continue;
				}
				std::vector<size_t> subsetReadKmerLengths;
				std::vector<TwobitString> subsetReadSequences;
				subsetReadKmerLengths.emplace_back(readKmerLengths[readi]);
				subsetReadSequences.emplace_back(readSequences[readi]);
				std::vector<MatchGroup> subsetMatches;
				phmap::flat_hash_set<size_t> includedAlignments;
				for (auto aln : alignmentsPerRead[readi])
				{
					if (relevantReads.count(matches[aln].leftRead) == 0)
					{
						size_t index = relevantReads.size();
						relevantReads[matches[aln].leftRead] = index;
						subsetReadKmerLengths.emplace_back(readKmerLengths[matches[aln].leftRead]);
						subsetReadSequences.emplace_back(readSequences[matches[aln].leftRead]);
					}
					if (relevantReads.count(matches[aln].rightRead) == 0)
					{
						size_t index = relevantReads.size();
						relevantReads[matches[aln].rightRead] = index;
						subsetReadKmerLengths.emplace_back(readKmerLengths[matches[aln].rightRead]);
						subsetReadSequences.emplace_back(readSequences[matches[aln].rightRead]);
					}
					assert(includedAlignments.count(aln) == 0);
					includedAlignments.emplace(aln);
					subsetMatches.emplace_back(matches[aln]);
					subsetMatches.back().leftRead = relevantReads.at(subsetMatches.back().leftRead);
					subsetMatches.back().rightRead = relevantReads.at(subsetMatches.back().rightRead);
				}
				for (auto pair : relevantReads)
				{
					for (size_t aln : alignmentsPerRead[pair.first])
					{
						if (matches[aln].leftRead == readi) continue;
						if (matches[aln].rightRead == readi) continue;
						if (relevantReads.count(matches[aln].leftRead) == 0) continue;
						if (relevantReads.count(matches[aln].rightRead) == 0) continue;
						if (includedAlignments.count(aln) == 1) continue;
						includedAlignments.emplace(aln);
						subsetMatches.emplace_back(matches[aln]);
						subsetMatches.back().leftRead = relevantReads.at(subsetMatches.back().leftRead);
						subsetMatches.back().rightRead = relevantReads.at(subsetMatches.back().rightRead);
					}
				}
				assert(subsetReadKmerLengths.size() == relevantReads.size());
				assert(subsetReadSequences.size() == relevantReads.size());
				TwobitString correctedSequence;
				size_t numCorrections = 0;
				std::tie(correctedSequence, numCorrections) = getCorrectedSequence(subsetReadKmerLengths, subsetReadSequences, subsetMatches, graphk);
				{
					std::string sequence = correctedSequence.toString();
					std::lock_guard<std::mutex> lock { printMutex };
					std::cout << ">" << readNames[readi] << std::endl;
					std::cout << sequence << std::endl;
					std::cerr << "read " << readNames[readi] << " had " << numCorrections << " corrections" << std::endl;
				}
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t searchk = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	size_t minAlignmentLength = std::stoull(argv[5]);
	const size_t graphk = std::stoull(argv[6]);
	const size_t minCoverage = 2;
	const size_t graphd = 50;
	std::vector<std::string> readFiles;
	for (size_t i = 7; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	MatchIndex matchIndex { searchk, numWindows, windowSize };
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
	std::vector<size_t> readKmerLengths;
	readKmerLengths.resize(readBasepairLengths.size());
	for (size_t i = 0; i < readBasepairLengths.size(); i++)
	{
		readKmerLengths[i] = readBasepairLengths[i] + 1 - graphk;
	}
	const std::vector<std::string>& readNames = storage.getNames();
	for (size_t i = 0; i < readNames.size(); i++)
	{
		std::cerr << "readname " << i << " " << readNames[i] << std::endl;
	}
	auto rawReadLengths = storage.getRawReadLengths();
	std::vector<bool> usableChunkmers = getFilteredValidChunks(matchIndex, rawReadLengths, std::numeric_limits<size_t>::max(), 10000, 10);
	std::vector<MatchGroup> matches = getMatches(matchIndex, readSequences, numThreads, rawReadLengths, minAlignmentLength, searchk, graphk, graphd, usableChunkmers, std::numeric_limits<size_t>::max());
	correctReads(readKmerLengths, readNames, readSequences, matches, graphk, numThreads);
}
