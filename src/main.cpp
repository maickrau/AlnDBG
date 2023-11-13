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

void writePaths(const std::string& filename, const std::vector<size_t>& readLengths, const std::vector<std::string>& readNames, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readUnitigGraphPaths, const size_t k)
{
	std::ofstream file { filename };
	for (size_t i = 0; i < readUnitigGraphPaths.size(); i++)
	{
		for (size_t j = 0; j < readUnitigGraphPaths[i].paths.size(); j++)
		{
			std::string pathstr;
			size_t pathLengthKmers = 0;
			for (size_t k = 0; k < readUnitigGraphPaths[i].paths[j].path.size(); k++)
			{
				pathstr += (readUnitigGraphPaths[i].paths[j].path[k] & firstBitUint64_t) ? ">" : "<";
				pathstr += std::to_string(readUnitigGraphPaths[i].paths[j].path[k] & maskUint64_t);
				pathLengthKmers += unitigGraph.lengths[readUnitigGraphPaths[i].paths[j].path[k] & maskUint64_t];
			}
			assert(pathLengthKmers > readUnitigGraphPaths[i].paths[j].pathLeftClipKmers + readUnitigGraphPaths[i].paths[j].pathRightClipKmers);
			size_t alnLengthKmers = pathLengthKmers - readUnitigGraphPaths[i].paths[j].pathLeftClipKmers - readUnitigGraphPaths[i].paths[j].pathRightClipKmers;
			file << readNames[i] << "\t" << readLengths[i] << "\t" << readUnitigGraphPaths[i].paths[j].readStartPos << "\t" << (readUnitigGraphPaths[i].paths[j].readStartPos + alnLengthKmers + k-1) << "\t+\t" << pathstr << "\t" << (pathLengthKmers + k-1) << "\t" << readUnitigGraphPaths[i].paths[j].pathLeftClipKmers << "\t" << (pathLengthKmers - readUnitigGraphPaths[i].paths[j].pathRightClipKmers + k-1) << "\t" << (alnLengthKmers + k-1) << "\t" << (alnLengthKmers + k-1) << "\t60" << std::endl;
		}
	}
}

bool haplotypeMismatch(const size_t leftRead, const size_t rightRead, const phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, bool>>>& readBelongsToPhaseBlock)
{
	if (readBelongsToPhaseBlock.count(leftRead) == 0) return false;
	if (readBelongsToPhaseBlock.count(rightRead) == 0) return false;
	size_t lefti = 0;
	size_t righti = 0;
	const std::vector<std::pair<size_t, bool>>& leftvec = readBelongsToPhaseBlock.at(leftRead);
	const std::vector<std::pair<size_t, bool>>& rightvec = readBelongsToPhaseBlock.at(rightRead);
	while (lefti < leftvec.size() && righti < rightvec.size())
	{
		if (leftvec[lefti].first < rightvec[righti].first)
		{
			lefti += 1;
			continue;
		}
		if (rightvec[righti].first < leftvec[lefti].first)
		{
			righti += 1;
			continue;
		}
		assert(leftvec[lefti].first == rightvec[righti].first);
		if (leftvec[lefti].second != rightvec[righti].second) return true;
		lefti += 1;
		righti += 1;
	}
	return false;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> makeGraph(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t minCoverage)
{
	KmerGraph kmerGraph;
	std::vector<ReadPathBundle> readKmerGraphPaths;
	std::tie(kmerGraph, readKmerGraphPaths) = makeKmerGraph(readLengths, matches, minCoverage);
	UnitigGraph unitigGraph;
	std::vector<ReadPathBundle> readUnitigGraphPaths;
	std::tie(unitigGraph, readUnitigGraphPaths) = makeUnitigGraph(kmerGraph, readKmerGraphPaths, minCoverage);
	std::cerr << unitigGraph.nodeCount() << " nodes before cleaning" << std::endl;
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 10);
	std::cerr << unitigGraph.nodeCount() << " nodes after cleaning" << std::endl;
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 10);
	std::cerr << unitigGraph.nodeCount() << " nodes after cleaning" << std::endl;
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 10);
	std::cerr << unitigGraph.nodeCount() << " nodes after cleaning" << std::endl;
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 10);
	std::cerr << unitigGraph.nodeCount() << " nodes after cleaning" << std::endl;
	return std::make_pair(std::move(unitigGraph), std::move(readUnitigGraphPaths));
}
/*
void forbidAlnsFromDifferentHaplotypes(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readUnitigGraphPaths, std::vector<MatchGroup>& matches, const std::vector<std::string>& readNames)
{
	auto phaseBlocks = getGraphPhaseBlockNodes(unitigGraph, readUnitigGraphPaths, 17);
	std::cerr << phaseBlocks.size() << " phase blocks" << std::endl;
	phmap::flat_hash_map<size_t, std::vector<std::pair<size_t, bool>>> readBelongsToPhaseBlock;
	for (size_t i = 0; i < phaseBlocks.size(); i++)
	{
		auto reads = getReadsPerPhaseBlock(unitigGraph, readUnitigGraphPaths, phaseBlocks[i]);
		for (auto read : reads.first)
		{
			readBelongsToPhaseBlock[read].emplace_back(i, true);
		}
		for (auto read : reads.second)
		{
			readBelongsToPhaseBlock[read].emplace_back(i, false);
		}
	}
	size_t countForbidden = 0;
	for (size_t i = matches.size()-1; i < matches.size(); i--)
	{
		if (haplotypeMismatch(matches[i].leftRead, matches[i].rightRead, readBelongsToPhaseBlock))
		{
			std::swap(matches[i], matches.back());
			matches.pop_back();
			countForbidden += 1;
		}
	}
	std::cerr << countForbidden << " alns forbidden by phasing" << std::endl;
}
*/
void makeGraph(const std::vector<size_t>& readLengths, const std::vector<std::string>& readNames, const std::vector<TwobitString>& readSequences, std::vector<MatchGroup>& matches, const size_t minCoverage, const std::string& outputFileName, const size_t k)
{
	UnitigGraph unitigGraph;
	std::vector<ReadPathBundle> readUnitigGraphPaths;
	std::tie(unitigGraph, readUnitigGraphPaths) = makeGraph(readLengths, matches, minCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, 20);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, 20);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, 20);
	unzipGraphLinearizable(unitigGraph, readUnitigGraphPaths, 20);
//	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraph(unitigGraph, readUnitigGraphPaths, getGraphPhaseBlockNodes(unitigGraph, readUnitigGraphPaths, 17));
//	forbidAlnsFromDifferentHaplotypes(unitigGraph, readUnitigGraphPaths, matches, readNames);
//	std::tie(unitigGraph, readUnitigGraphPaths) = makeGraph(readLengths, matches, minCoverage);
	auto nodeSequences  = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	writeGraph(outputFileName, unitigGraph, nodeSequences, k);
//	writeGraph(outputFileName, unitigGraph, k);
	writePaths("paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
}

void doHaplofilter(std::vector<MatchGroup>& matches, const std::vector<size_t>& readLengths)
{
	std::vector<bool> kept = getValidAlignments(matches, readLengths, 5, 2);
	size_t numRemoved = 0;
	assert(kept.size() == matches.size());
	for (size_t i = matches.size()-1; i < matches.size(); i--)
	{
		if (kept[i]) continue;
		std::swap(matches[i], matches.back());
		matches.pop_back();
		numRemoved += 1;
	}
	std::cerr << numRemoved << " mapping matches removed by haplofilter" << std::endl;
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	size_t minAlignmentLength = std::stoull(argv[5]);
	const size_t graphk = 11;
	const size_t minCoverage = 2;
	const size_t graphd = 50;
	std::vector<std::string> readFiles;
	for (size_t i = 6; i < argc; i++)
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
	std::vector<size_t> readKmerLengths;
	readKmerLengths.resize(readBasepairLengths.size());
	for (size_t i = 0; i < readBasepairLengths.size(); i++)
	{
		readKmerLengths[i] = readBasepairLengths[i] + 1 - graphk;
	}
	const std::vector<std::string>& readNames = storage.getNames();
	std::mutex printMutex;
	std::vector<MatchGroup> matches;
	auto result = matchIndex.iterateMatchChains(numThreads, 2, 1000, 50, storage.getRawReadLengths(), [&printMutex, &matches, minAlignmentLength, k, graphd](const size_t left, const size_t leftstart, const size_t leftend, const bool leftFw, const size_t right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		if (leftend+1-leftstart < minAlignmentLength) return;
		if (rightend+1-rightstart < minAlignmentLength) return;
		std::lock_guard<std::mutex> lock { printMutex };
		assert(leftFw);
		size_t numAlnChunks = std::max(leftend+1-leftstart, rightend+1-rightstart)/(std::numeric_limits<uint16_t>::max()-k-graphd-100)+1;
		assert(numAlnChunks >= 1);
		if (numAlnChunks == 1)
		{
			matches.emplace_back();
			matches.back().leftRead = left;
			matches.back().rightRead = right;
			matches.back().rightFw = rightFw;
			matches.back().leftStart = leftstart;
			matches.back().leftEnd = leftend+1;
			matches.back().rightStart = rightstart;
			matches.back().rightEnd = rightend+1;
		}
		else
		{
			assert(numAlnChunks >= 2);
			size_t leftPerChunk = (double)(leftend+1-leftstart)/(double)numAlnChunks;
			size_t rightPerChunk = (double)(rightend+1-rightstart)/(double)numAlnChunks;
			for (size_t i = 0; i < numAlnChunks; i++)
			{
				matches.emplace_back();
				matches.back().leftRead = left;
				matches.back().rightRead = right;
				matches.back().rightFw = rightFw;
				matches.back().leftStart = leftstart + i * leftPerChunk;
				matches.back().leftEnd = leftstart + (i+1) * leftPerChunk + k + graphd;
				matches.back().rightStart = rightstart + i * rightPerChunk;
				matches.back().rightEnd = rightstart + (i+1) * rightPerChunk + k + graphd;
			}
			matches.back().leftEnd = leftend+1;
			matches.back().rightEnd = rightend+1;
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
	addKmerMatches(numThreads, readSequences, matches, graphk, graphd);
//	doHaplofilter(matches, readKmerLengths);
	makeGraph(readKmerLengths, readNames, readSequences, matches, minCoverage, "graph.gfa", graphk);
}
