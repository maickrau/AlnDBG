#include <mutex>
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

void makeGraph(const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t minCoverage, const std::string& outputFileName, const size_t k)
{
	std::vector<RankBitvector> breakpoints = extendBreakpoints(readLengths, matches);
	size_t countBreakpoints = 0;
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		for (size_t j = 0; j < breakpoints[i].size(); j++)
		{
			if (breakpoints[i].get(j)) countBreakpoints += 1;
		}
	}
	std::cerr << countBreakpoints << " breakpoints" << std::endl;
	std::vector<uint64_t> segments = mergeSegments(readLengths, matches, breakpoints, countBreakpoints);
	std::cerr << segments.size() << " segments" << std::endl;
	RankBitvector segmentToNode = getSegmentToNode(segments, minCoverage);
	size_t countNodes = (segmentToNode.getRank(segmentToNode.size()-1) + (segmentToNode.get(segmentToNode.size()-1) ? 1 : 0));
	std::cerr << countNodes << " nodes pre coverage filter" << std::endl;
	std::vector<size_t> nodeCoverage = getNodeCoverage(segments, segmentToNode, countNodes);
	std::vector<size_t> nodeLength = getNodeLengths(segments, segmentToNode, breakpoints, countNodes);
	MostlySparse2DHashmap<uint8_t, size_t> edgeCoverages = getEdgeCoverages(readLengths, segmentToNode, segments, breakpoints, minCoverage, nodeCoverage, countNodes);
	std::cerr << edgeCoverages.size() << " edges pre coverage filter" << std::endl;
	std::vector<std::vector<uint64_t>> unitigs;
	std::vector<uint64_t> unitigLeftmostNode;
	std::vector<uint64_t> unitigRightmostNode;
	std::tie(unitigs, unitigLeftmostNode, unitigRightmostNode) = getUnitigs(countNodes, minCoverage, edgeCoverages);
	std::vector<size_t> unitigLength;
	std::vector<double> unitigCoverage;
	std::tie(unitigLength, unitigCoverage) = getUnitigLengthAndCoverage(unitigs, nodeCoverage, nodeLength);
	MostlySparse2DHashmap<uint8_t, size_t> unitigEdgeCoverages = getUnitigEdgeCoverages(unitigLeftmostNode, unitigRightmostNode, minCoverage, edgeCoverages);
	// writeGraph(outputFileName, nodeCoverage, nodeLength, edgeCoverages, minCoverage, k);
	writeGraph(outputFileName, unitigCoverage, unitigLength, unitigEdgeCoverages, minCoverage, k);
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	size_t minAlignmentLength = std::stoull(argv[5]);
	const size_t graphk = 31;
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
	makeGraph(readKmerLengths, matches, minCoverage, "graph.gfa", graphk);
}
