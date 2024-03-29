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
#include "MultiplexResolverCaller.h"

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
				if (k > 0)
				{
					std::pair<size_t, bool> from { readUnitigGraphPaths[i].paths[j].path[k-1] & maskUint64_t, readUnitigGraphPaths[i].paths[j].path[k-1] & firstBitUint64_t };
					std::pair<size_t, bool> to { readUnitigGraphPaths[i].paths[j].path[k] & maskUint64_t, readUnitigGraphPaths[i].paths[j].path[k] & firstBitUint64_t };
					assert(unitigGraph.edgeKmerOverlaps.hasValue(from, to));
					size_t overlap = unitigGraph.edgeKmerOverlaps.get(from, to);
					pathLengthKmers -= overlap;
				}
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

std::pair<UnitigGraph, std::vector<ReadPathBundle>> makeGraph(const std::vector<TwobitString>& readSequences, const std::vector<size_t>& readLengths, const std::vector<MatchGroup>& matches, const size_t minCoverage, const size_t numThreads, const size_t graphk)
{
	KmerGraph kmerGraph;
	std::vector<ReadPathBundle> readKmerGraphPaths;
	std::tie(kmerGraph, readKmerGraphPaths) = makeKmerGraph(readSequences, readLengths, matches, minCoverage, numThreads, graphk);
	UnitigGraph unitigGraph;
	std::vector<ReadPathBundle> readUnitigGraphPaths;
	std::tie(unitigGraph, readUnitigGraphPaths) = makeUnitigGraph(kmerGraph, readKmerGraphPaths, minCoverage);
	std::cerr << unitigGraph.nodeCount() << " nodes before cleaning" << std::endl;
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 7);
	std::cerr << unitigGraph.nodeCount() << " nodes after cleaning" << std::endl;
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 7);
	std::cerr << unitigGraph.nodeCount() << " nodes after cleaning" << std::endl;
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 7);
	std::cerr << unitigGraph.nodeCount() << " nodes after cleaning" << std::endl;
	std::tie(unitigGraph, readUnitigGraphPaths) = cleanUnitigGraph(unitigGraph, readUnitigGraphPaths, 7);
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

std::pair<UnitigGraph, std::vector<ReadPathBundle>> getInitialGraph(const std::vector<size_t>& readLengths, const std::vector<std::string>& readNames, const std::vector<TwobitString>& readSequences, std::vector<MatchGroup>& matches, const size_t minCoverage, const size_t k, const double approxOneHapCoverage, const std::string& prefix, const size_t numThreads)
{
	UnitigGraph unitigGraph;
	std::vector<ReadPathBundle> readUnitigGraphPaths;
	std::tie(unitigGraph, readUnitigGraphPaths) = makeGraph(readSequences, readLengths, matches, minCoverage, numThreads, k);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
//	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraph(unitigGraph, readUnitigGraphPaths, getGraphPhaseBlockNodes(unitigGraph, readUnitigGraphPaths, 17));
//	forbidAlnsFromDifferentHaplotypes(unitigGraph, readUnitigGraphPaths, matches, readNames);
//	std::tie(unitigGraph, readUnitigGraphPaths) = makeGraph(readLengths, matches, minCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);
	auto nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	writeGraph(prefix + "graph.gfa", unitigGraph, nodeSequences, k);
//	writeGraph(outputFileName, unitigGraph, k);
	writePaths(prefix + "paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	return std::make_pair(std::move(unitigGraph), std::move(readUnitigGraphPaths));
}

void makeGraph(const std::vector<size_t>& readLengths, const std::vector<std::string>& readNames, const std::vector<TwobitString>& readSequences, std::vector<MatchGroup>& matches, const size_t minCoverage, const std::string& outputFileName, const size_t k, const double approxOneHapCoverage, const size_t numThreads)
{
	UnitigGraph unitigGraph;
	std::vector<ReadPathBundle> readUnitigGraphPaths;
	std::tie(unitigGraph, readUnitigGraphPaths) = makeGraph(readSequences, readLengths, matches, minCoverage, numThreads, k);
	// std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	// std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	// std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
//	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraph(unitigGraph, readUnitigGraphPaths, getGraphPhaseBlockNodes(unitigGraph, readUnitigGraphPaths, 17));
//	forbidAlnsFromDifferentHaplotypes(unitigGraph, readUnitigGraphPaths, matches, readNames);
//	std::tie(unitigGraph, readUnitigGraphPaths) = makeGraph(readLengths, matches, minCoverage);
	auto nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	writeGraph(outputFileName, unitigGraph, nodeSequences, k);
//	writeGraph(outputFileName, unitigGraph, k);
	writePaths("paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	// for (size_t i = 0; i < 5; i++)
	// {
		std::tie(unitigGraph, readUnitigGraphPaths) = runMBGMultiplexResolution(unitigGraph, readUnitigGraphPaths, k, 1000);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);

	//nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	//writeGraph("hmm1-graph.gfa", unitigGraph, nodeSequences, k);
	writeGraph("hmm1-graph.gfa", unitigGraph, k);
	writePaths("hmm1-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
	for (size_t i = 0; i < readUnitigGraphPaths.size(); i++)
	{
		for (size_t j = 0; j < readUnitigGraphPaths[i].paths.size(); j++)
		{
			assert(readUnitigGraphPaths[i].paths[j].readStartPos < readUnitigGraphPaths[i].readLength);
		}
	}
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmersLocation(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 100, 100000);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphPolyploidTransitiveClosure(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphPolyploidTransitiveClosure(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);

	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("hmm2-graph.gfa", unitigGraph, nodeSequences, k);
	writeGraph("hmm2-graph.gfa", unitigGraph, k);
	writePaths("hmm2-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmersLocation(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 100, 100000);
		// std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);
	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("hmm3-graph.gfa", unitigGraph, nodeSequences, k);
	writeGraph("hmm3-graph.gfa", unitigGraph, k);
	writePaths("hmm3-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	std::exit(0);

		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
	writeGraph("hmm4-graph.gfa", unitigGraph, k);
	writePaths("hmm4-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmersLocation(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 10000, 100000);
	writeGraph("hmm5-graph.gfa", unitigGraph, k);
	writePaths("hmm5-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmersLocation(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 20000, 100000);
	writeGraph("hmm6-graph.gfa", unitigGraph, k);
	writePaths("hmm6-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
std::exit(0);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmersLocation(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 20000, 100000);
	writeGraph("hmm7-graph.gfa", unitigGraph, k);
	writePaths("hmm7-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	std::exit(0);
		std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);

	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("one1phase-graph.gfa", unitigGraph, nodeSequences, k);
	writeGraph("one1phase-graph.gfa", unitigGraph, k);
	writePaths("one1phase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmersLocation(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 2000, 100000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 2000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 2000, 2000, 100000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphPolyploidTransitiveClosure(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 3, 2);
		std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 3, 2);
	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("one2phase-graph.gfa", unitigGraph, nodeSequences, k);
	writeGraph("one2phase-graph.gfa", unitigGraph, k);
	writePaths("one2phase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmersLocation(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 10000, 100000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 10000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 10000, 10000, 100000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphPolyploidTransitiveClosure(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 3, 2);
		std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 3, 2);
	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("one3phase-graph.gfa", unitigGraph, nodeSequences, k);
	writeGraph("one3phase-graph.gfa", unitigGraph, k);
	writePaths("one3phase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmersLocation(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 20000, 100000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 20000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 20000, 20000, 100000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphPolyploidTransitiveClosure(unitigGraph, readUnitigGraphPaths, k, numThreads, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 3, 2);
	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("one4phase-graph.gfa", unitigGraph, nodeSequences, k);
	writeGraph("one4phase-graph.gfa", unitigGraph, k);
	writePaths("one4phase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
		std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, k, approxOneHapCoverage, 3, 2);
		/*
	nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	writeGraph("one4phase-graph.gfa", unitigGraph, nodeSequences, k);
	writePaths("one4phase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	std::exit(0);

		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000);
		std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000);
		std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000, 2000, 5);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000, 10000, 5);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000, 10000, 2);
		// std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHapmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000);
		// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHapmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000);
		// std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
		std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 4, 2);
		std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	// 	if (i < 4)
	// 	{
	// 		nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// 		writeGraph("iteration-" + std::to_string(i) + "-graph.gfa", unitigGraph, nodeSequences, k);
	// 		writePaths("iteration-" + std::to_string(i) + "-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	// 	}
	// }
	nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	writeGraph("onephase-graph.gfa", unitigGraph, nodeSequences, k);
	writePaths("onephase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000, 2000, 5);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHaploforks(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 200, true);
	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("one1phase-graph.gfa", unitigGraph, nodeSequences, k);
	// writePaths("one1phase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000, 10000, 5);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHaploforks(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 1000, true);
	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("one2phase-graph.gfa", unitigGraph, nodeSequences, k);
	// writePaths("one2phase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	// std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000, 10000, 2);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHaploforks(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000, false);
	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("one3phase-graph.gfa", unitigGraph, nodeSequences, k);
	// writePaths("one3phase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);

	// std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 4, 2);
	// std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	// std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);

	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	// // nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// // writeGraph("twophase-graph.gfa", unitigGraph, nodeSequences, k);
	// // writePaths("twophase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHapmers(unitigGraph, readUnitigGraphPaths, 20, 2000);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, 20);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHapmers(unitigGraph, readUnitigGraphPaths, 20, 10000);
	// // std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, 20);

	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphChainmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHapmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHaploforks(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000, false);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphHapmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000);
	std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 4, 2);
	nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	writeGraph("twophase-graph.gfa", unitigGraph, nodeSequences, k);
	writePaths("twophase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000, 2000, 5);

	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 2000, 10000, 20);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 4, 2);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000, 10000, 2);
	nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	writeGraph("threephase-graph.gfa", unitigGraph, nodeSequences, k);
	writePaths("threephase-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphDiploidMEC(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSimpleStructures(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = connectChainGaps(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 4, 2);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 200, 10000, 2);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = unzipGraphLocalUniqmers(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage, 10000, 10000, 2);
	std::tie(unitigGraph, readUnitigGraphPaths) = resolveSpannedTangles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	std::tie(unitigGraph, readUnitigGraphPaths) = popHaploidChainBubbles(unitigGraph, readUnitigGraphPaths, approxOneHapCoverage);
	*/
	// nodeSequences = getNodeSequences(unitigGraph, readUnitigGraphPaths, k, readSequences);
	// writeGraph("phased-graph.gfa", unitigGraph, nodeSequences, k);
	writeGraph("phased-graph.gfa", unitigGraph, k);
	writePaths("phased-paths.gaf", readLengths, readNames, unitigGraph, readUnitigGraphPaths, k);
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

bool annotationsMismatch(const std::vector<std::vector<ChainPosition>>& readAnnotations, const std::vector<size_t>& readLengths, const size_t leftRead, const size_t rightRead, const size_t leftStart, const size_t rightStart, const bool rightFw)
{
	if (readAnnotations[leftRead].size() == 0) return false;
	if (readAnnotations[rightRead].size() == 0) return false;
	std::vector<ChainPosition> annotationsLeft;
	std::vector<ChainPosition> annotationsRight;
	for (const auto chain : readAnnotations[leftRead])
	{
		annotationsLeft.emplace_back(chain);
		annotationsLeft.back().chainStartPosInRead -= (int)leftStart;
		annotationsLeft.back().chainEndPosInRead -= (int)leftStart;
	}
	for (const auto chain : readAnnotations[rightRead])
	{
		if (rightFw)
		{
			annotationsRight.emplace_back(chain);
			annotationsRight.back().chainStartPosInRead -= (int)rightStart;
			annotationsRight.back().chainEndPosInRead -= (int)rightStart;
		}
		else
		{
			annotationsRight.emplace_back(chain);
			annotationsRight.back().chain ^= firstBitUint64_t;
			std::swap(annotationsRight.back().chainStartPosInRead, annotationsRight.back().chainEndPosInRead);
			annotationsRight.back().chainStartPosInRead = readLengths[rightRead] - annotationsRight.back().chainStartPosInRead;
			annotationsRight.back().chainEndPosInRead = readLengths[rightRead] - annotationsRight.back().chainEndPosInRead;
			annotationsRight.back().chainStartPosInRead -= (int)rightStart;
			annotationsRight.back().chainEndPosInRead -= (int)rightStart;
		}
	}
	int leftReadStart = -(int)leftStart;
	int rightReadStart = -(int)rightStart;
	int leftReadEnd = leftReadStart + (int)readLengths[leftRead];
	int rightReadEnd = rightReadStart + (int)readLengths[rightRead];
	assert(leftReadEnd > 0);
	assert(rightReadEnd > 0);
	for (size_t i = 0; i < annotationsLeft.size(); i++)
	{
		if (annotationsLeft[i].chainEndPosInRead < rightReadStart+1000) continue;
		if (annotationsLeft[i].chainStartPosInRead > rightReadEnd-1000) continue;
		if (annotationsLeft[i].numKmerMatches < 100) continue;
		bool hasMatch = false;
		for (size_t j = 0; j < annotationsRight.size(); j++)
		{
			if (annotationsRight[j].chain != annotationsLeft[i].chain) continue;
			if (annotationsRight[j].chainEndPosInRead < annotationsLeft[i].chainEndPosInRead-1000) continue;
			if (annotationsRight[j].chainEndPosInRead > annotationsLeft[i].chainEndPosInRead+1000) continue;
			if (annotationsRight[j].chainStartPosInRead < annotationsLeft[i].chainStartPosInRead-1000) continue;
			if (annotationsRight[j].chainStartPosInRead > annotationsLeft[i].chainStartPosInRead+1000) continue;
			hasMatch = true;
		}
		if (!hasMatch) return true;
	}
	for (size_t i = 0; i < annotationsRight.size(); i++)
	{
		if (annotationsRight[i].chainEndPosInRead < rightReadStart+1000) continue;
		if (annotationsRight[i].chainStartPosInRead > rightReadEnd-1000) continue;
		if (annotationsRight[i].numKmerMatches < 100) continue;
		bool hasMatch = false;
		for (size_t j = 0; j < annotationsLeft.size(); j++)
		{
			if (annotationsLeft[j].chain != annotationsRight[i].chain) continue;
			if (annotationsLeft[j].chainEndPosInRead < annotationsRight[i].chainEndPosInRead-1000) continue;
			if (annotationsLeft[j].chainEndPosInRead > annotationsRight[i].chainEndPosInRead+1000) continue;
			if (annotationsLeft[j].chainStartPosInRead < annotationsRight[i].chainStartPosInRead-1000) continue;
			if (annotationsLeft[j].chainStartPosInRead > annotationsRight[i].chainStartPosInRead+1000) continue;
			hasMatch = true;
		}
		if (!hasMatch) return true;
	}
	return false;
}

bool annotationsMismatch(const std::vector<std::vector<HaplotypeInformativeSite>>& readAnnotations, const std::vector<size_t>& readLengths, const size_t leftRead, const size_t rightRead, const size_t leftStart, const size_t rightStart, const bool rightFw)
{
	if (readAnnotations[leftRead].size() == 0) return false;
	if (readAnnotations[rightRead].size() == 0) return false;
	std::vector<HaplotypeInformativeSite> annotationsLeft;
	std::vector<HaplotypeInformativeSite> annotationsRight;
	for (const auto chain : readAnnotations[leftRead])
	{
		annotationsLeft.emplace_back(chain);
		annotationsLeft.back().readPos -= (int)leftStart;
	}
	for (const auto chain : readAnnotations[rightRead])
	{
		if (rightFw)
		{
			annotationsRight.emplace_back(chain);
			annotationsRight.back().readPos -= (int)rightStart;
		}
		else
		{
			annotationsRight.emplace_back(chain);
			annotationsRight.back().bubble ^= firstBitUint64_t;
			annotationsRight.back().allele ^= firstBitUint64_t;
			annotationsRight.back().readPos = readLengths[rightRead] - annotationsRight.back().readPos;
		}
	}
	int leftReadStart = -(int)leftStart;
	int rightReadStart = -(int)rightStart;
	int leftReadEnd = leftReadStart + (int)readLengths[leftRead];
	int rightReadEnd = rightReadStart + (int)readLengths[rightRead];
	assert(leftReadEnd > 0);
	assert(rightReadEnd > 0);
	size_t missingCount = 0;
	size_t mismatchCount = 0;
	for (size_t i = 0; i < annotationsLeft.size(); i++)
	{
		if (annotationsLeft[i].readPos < rightReadStart+1000) continue;
		if (annotationsLeft[i].readPos > rightReadEnd-1000) continue;
		bool hasMatch = false;
		bool hasMismatch = false;
		for (size_t j = 0; j < annotationsRight.size(); j++)
		{
			if (annotationsRight[j].bubble != annotationsLeft[i].bubble) continue;
			if (annotationsRight[j].readPos < annotationsLeft[i].readPos-100) continue;
			if (annotationsRight[j].readPos > annotationsLeft[i].readPos+100) continue;
			if (annotationsRight[j].allele != annotationsLeft[i].allele)
			{
				hasMismatch = true;
			}
			else
			{
				hasMatch = true;
			}
		}
		if (!hasMatch && hasMismatch) mismatchCount += 1;
		if (!hasMatch) missingCount += 1;
	}
	for (size_t i = 0; i < annotationsRight.size(); i++)
	{
		if (annotationsRight[i].readPos < rightReadStart+1000) continue;
		if (annotationsRight[i].readPos > rightReadEnd-1000) continue;
		bool hasMatch = false;
		bool hasMismatch = false;
		for (size_t j = 0; j < annotationsLeft.size(); j++)
		{
			if (annotationsLeft[j].bubble != annotationsRight[i].bubble) continue;
			if (annotationsLeft[j].readPos < annotationsRight[i].readPos-100) continue;
			if (annotationsLeft[j].readPos > annotationsRight[i].readPos+100) continue;
			if (annotationsLeft[j].allele != annotationsRight[i].allele)
			{
				hasMismatch = true;
			}
			else
			{
				hasMatch = true;
			}
		}
		if (!hasMatch && hasMismatch) mismatchCount += 1;
		if (!hasMatch) mismatchCount += 1;
	}
	if (mismatchCount >= 1) return true;
	if (missingCount >= 3) return true;
	return false;
}

void tryExtendTransitiveMatch(std::vector<MatchGroup>& result, const size_t firstRead, const size_t middleRead, const size_t thirdRead, const MatchGroup& firstMatch, const MatchGroup& secondMatch, const std::vector<size_t>& rawReadLengths, const size_t graphk)
{
	assert(firstRead != thirdRead);
	assert(firstRead != middleRead);
	assert(middleRead != thirdRead);
	int firstStartPosInMiddle;
	int firstEndPosInMiddle;
	bool firstFwInMiddle;
	if (firstMatch.rightFw)
	{
		firstFwInMiddle = true;
		if (firstMatch.leftRead == firstRead)
		{
			assert(firstMatch.rightRead == middleRead);
			firstStartPosInMiddle = (int)firstMatch.rightStart - (int)firstMatch.leftStart;
			firstEndPosInMiddle = (int)firstMatch.rightEnd + ((int)rawReadLengths[firstRead]-1) - (int)firstMatch.leftEnd;
		}
		else
		{
			assert(firstMatch.leftRead == middleRead);
			assert(firstMatch.rightRead == firstRead);
			firstStartPosInMiddle = (int)firstMatch.leftStart - (int)firstMatch.rightStart;
			firstEndPosInMiddle = (int)firstMatch.leftEnd + ((int)rawReadLengths[firstRead]-1) - (int)firstMatch.rightEnd;
		}
	}
	else
	{
		firstFwInMiddle = false;
		if (firstMatch.leftRead == firstRead)
		{
			assert(firstMatch.rightRead == middleRead);
			firstStartPosInMiddle = (int)firstMatch.rightStart + (int)firstMatch.leftEnd - ((int)rawReadLengths[firstRead]-1);
			firstEndPosInMiddle = (int)firstMatch.rightEnd + (int)firstMatch.leftStart;
		}
		else
		{
			assert(firstMatch.leftRead == middleRead);
			assert(firstMatch.rightRead == firstRead);
			firstStartPosInMiddle = (int)firstMatch.leftStart + (int)firstMatch.rightEnd - ((int)rawReadLengths[firstRead]-1);
			firstEndPosInMiddle = (int)firstMatch.leftEnd + (int)firstMatch.rightStart;
		}
	}
	assert(firstEndPosInMiddle > firstStartPosInMiddle);
	int thirdStartPosInMiddle;
	int thirdEndPosInMiddle;
	bool thirdFwInMiddle;
	if (secondMatch.rightFw)
	{
		thirdFwInMiddle = true;
		if (secondMatch.leftRead == thirdRead)
		{
			assert(secondMatch.rightRead == middleRead);
			thirdStartPosInMiddle = (int)secondMatch.rightStart - (int)secondMatch.leftStart;
			thirdEndPosInMiddle = (int)secondMatch.rightEnd + ((int)rawReadLengths[thirdRead]-1) - (int)secondMatch.leftEnd;
		}
		else
		{
			assert(secondMatch.leftRead == middleRead);
			assert(secondMatch.rightRead == thirdRead);
			thirdStartPosInMiddle = (int)secondMatch.leftStart - (int)secondMatch.rightStart;
			thirdEndPosInMiddle = (int)secondMatch.leftEnd + ((int)rawReadLengths[thirdRead]-1) - (int)secondMatch.rightEnd;
		}
	}
	else
	{
		thirdFwInMiddle = false;
		if (secondMatch.leftRead == thirdRead)
		{
			assert(secondMatch.rightRead == middleRead);
			thirdStartPosInMiddle = (int)secondMatch.rightStart + (int)secondMatch.leftEnd - ((int)rawReadLengths[thirdRead]-1);
			thirdEndPosInMiddle = (int)secondMatch.rightEnd + (int)secondMatch.leftStart;
		}
		else
		{
			assert(secondMatch.leftRead == middleRead);
			assert(secondMatch.rightRead == thirdRead);
			thirdStartPosInMiddle = (int)secondMatch.leftStart + (int)secondMatch.rightEnd - ((int)rawReadLengths[thirdRead]-1);
			thirdEndPosInMiddle = (int)secondMatch.leftEnd + (int)secondMatch.rightStart;
		}
	}
	assert(thirdEndPosInMiddle > thirdStartPosInMiddle);
	std::cerr << "reads first " << firstRead << " middle " << middleRead << " third " << thirdRead << std::endl;
	std::cerr << "first reads " << firstMatch.leftRead << " (len " << rawReadLengths[firstMatch.leftRead] << ") " << firstMatch.rightRead << " (len " << rawReadLengths[firstMatch.rightRead] << ")" << std::endl;
	std::cerr << "first match " << firstMatch.leftStart << " " << firstMatch.leftEnd << " " << firstMatch.rightStart << " " << firstMatch.rightEnd << " " << (firstMatch.rightFw ? "fw" : "bw") << std::endl;
	std::cerr << "second reads " << secondMatch.leftRead << " (len " << rawReadLengths[secondMatch.leftRead] << ") " << secondMatch.rightRead << " (len " << rawReadLengths[secondMatch.rightRead] << ")" << std::endl;
	std::cerr << "second match " << secondMatch.leftStart << " " << secondMatch.leftEnd << " " << secondMatch.rightStart << " " << secondMatch.rightEnd << " " << (secondMatch.rightFw ? "fw" : "bw") << std::endl;
	std::cerr << "poses " << firstStartPosInMiddle << " " << firstEndPosInMiddle << " " << thirdStartPosInMiddle << " " << thirdEndPosInMiddle << std::endl;
	if (thirdStartPosInMiddle >= firstEndPosInMiddle) return;
	if (firstStartPosInMiddle >= thirdEndPosInMiddle) return;
	int overlapStart = std::max(firstStartPosInMiddle, thirdStartPosInMiddle);
	int overlapEnd = std::min(firstEndPosInMiddle, thirdEndPosInMiddle);
	if (overlapEnd - overlapStart < 2*graphk) return;
	assert(overlapStart >= firstStartPosInMiddle);
	assert(overlapEnd <= firstEndPosInMiddle);
	assert(overlapStart >= thirdStartPosInMiddle);
	assert(overlapEnd <= thirdEndPosInMiddle);
	std::cerr << "found overlap " << overlapStart << " " << overlapEnd << std::endl;
	assert(overlapEnd > overlapStart);
	result.emplace_back();
	result.back().leftRead = firstRead;
	result.back().rightRead = thirdRead;
	result.back().rightFw = (firstFwInMiddle == thirdFwInMiddle);
	if (firstFwInMiddle)
	{
		result.back().leftStart = overlapStart - firstStartPosInMiddle;
		result.back().leftEnd = overlapEnd - firstStartPosInMiddle;
	}
	else
	{
		result.back().leftStart = firstEndPosInMiddle - overlapEnd;
		result.back().leftEnd = firstEndPosInMiddle - overlapStart;
	}
	result.back().leftStart *= rawReadLengths[result.back().leftRead] / (double)(firstEndPosInMiddle - firstStartPosInMiddle);
	result.back().leftEnd *= rawReadLengths[result.back().leftRead] / (double)(firstEndPosInMiddle - firstStartPosInMiddle);
	assert(result.back().leftStart < result.back().leftEnd);
	assert(result.back().leftEnd <= rawReadLengths[result.back().leftRead]);
	if (result.back().leftEnd == rawReadLengths[result.back().leftRead]) result.back().leftEnd = rawReadLengths[result.back().leftRead]-1;
	if (thirdFwInMiddle)
	{
		result.back().rightStart = overlapStart - thirdStartPosInMiddle;
		result.back().rightEnd = overlapEnd - thirdStartPosInMiddle;
	}
	else
	{
		result.back().rightStart = thirdEndPosInMiddle - overlapEnd;
		result.back().rightEnd = thirdEndPosInMiddle - overlapStart;
	}
	result.back().rightStart *= rawReadLengths[result.back().rightRead] / (double)(thirdEndPosInMiddle - thirdStartPosInMiddle);
	result.back().rightEnd *= rawReadLengths[result.back().rightRead] / (double)(thirdEndPosInMiddle - thirdStartPosInMiddle);
	assert(result.back().rightStart < result.back().rightEnd);
	assert(result.back().rightEnd <= rawReadLengths[result.back().rightRead]);
	if (result.back().rightEnd == rawReadLengths[result.back().rightRead]) result.back().rightEnd = rawReadLengths[result.back().rightRead]-1;
}

std::vector<MatchGroup> extendTransitiveMatches(const std::vector<MatchGroup>& oldMatches, const std::vector<size_t>& rawReadLengths, const size_t graphk)
{
	std::cerr << "try transitively extend matches" << std::endl;
	std::vector<MatchGroup> result = oldMatches;
	std::vector<std::vector<size_t>> matchesPerRead;
	matchesPerRead.resize(rawReadLengths.size());
	for (size_t i = 0; i < oldMatches.size(); i++)
	{
		assert(oldMatches[i].leftRead < matchesPerRead.size());
		assert(oldMatches[i].rightRead < matchesPerRead.size());
		matchesPerRead[oldMatches[i].leftRead].push_back(i);
		matchesPerRead[oldMatches[i].rightRead].push_back(i);
	}
	for (size_t i = 0; i < rawReadLengths.size(); i++)
	{
		for (size_t match : matchesPerRead[i])
		{
			assert(oldMatches[match].leftRead == i || oldMatches[match].rightRead == i);
			assert(oldMatches[match].leftRead != i || oldMatches[match].rightRead != i);
			size_t otherRead = i;
			if (oldMatches[match].leftRead == i)
			{
				otherRead = oldMatches[match].rightRead;
			}
			else
			{
				otherRead = oldMatches[match].leftRead;
			}
			assert(otherRead != i);
			for (size_t match2 : matchesPerRead[otherRead])
			{
				assert(oldMatches[match2].leftRead == otherRead || oldMatches[match2].rightRead == otherRead);
				assert(oldMatches[match2].leftRead != otherRead || oldMatches[match2].rightRead != otherRead);
				if (oldMatches[match2].leftRead == i || oldMatches[match2].rightRead == i) continue;
				size_t thirdRead;
				if (oldMatches[match2].leftRead == otherRead)
				{
					thirdRead = oldMatches[match2].rightRead;
				}
				else
				{
					thirdRead = oldMatches[match2].leftRead;
				}
				tryExtendTransitiveMatch(result, i, otherRead, thirdRead, oldMatches[match], oldMatches[match2], rawReadLengths, graphk);
			}
		}
	}
	return result;
}

std::vector<MatchGroup> getMatchesFilteredByAnnotation(const MatchIndex& matchIndex, const std::vector<TwobitString>& readSequences, size_t numThreads, const std::vector<size_t>& rawReadLengths, const size_t minAlignmentLength, const size_t k, const size_t graphk, const size_t graphd, const size_t maxIndexCoverage, const std::vector<bool>& usableChunkmers, const std::vector<std::vector<ChainPosition>>& readAnnotations)
{
	std::mutex printMutex;
	std::vector<MatchGroup> matches;
	size_t countFiltered = 0;
	auto result = matchIndex.iterateMatchChains(numThreads, 2, maxIndexCoverage, 50, rawReadLengths, usableChunkmers, [&printMutex, &matches, minAlignmentLength, &rawReadLengths, k, graphd, &readAnnotations, &countFiltered](const size_t left, const size_t leftstart, const size_t leftend, const bool leftFw, const size_t right, const size_t rightstart, const size_t rightend, const bool rightFw)
	{
		if (leftend+1-leftstart < minAlignmentLength) return;
		if (rightend+1-rightstart < minAlignmentLength) return;
		if (annotationsMismatch(readAnnotations, rawReadLengths, left, right, leftstart, rightstart, rightFw))
		{
			countFiltered += 1;
			return;
		}
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
	// matches = extendTransitiveMatches(matches, rawReadLengths, graphk);
	// std::cerr << matches.size() << " mapping matches after transitive extension" << std::endl;
	matches = addKmerMatches(numThreads, readSequences, matches, graphk, graphd);
	size_t filteredFromKmerGaps = 0;
	size_t countRemainingKmerMatches = 0;
	// for (size_t i = matches.size()-1; i < matches.size(); i--)
	// {
	// 	size_t lastKmerLeftMatchPos = 0;
	// 	size_t lastKmerRightMatchPos = 0;
	// 	bool remove = false;
	// 	for (size_t j = 0; j < matches[i].matches.size(); j++)
	// 	{
	// 		assert(j == 0 || matches[i].matches[j].leftStart >= matches[i].matches[j-1].leftStart);
	// 		if (matches[i].matches[j].leftStart > lastKmerLeftMatchPos + 500)
	// 		{
	// 			remove = true;
	// 			break;
	// 		}
	// 		if (matches[i].matches[j].rightStart > lastKmerRightMatchPos + 500)
	// 		{
	// 			remove = true;
	// 			break;
	// 		}
	// 		lastKmerLeftMatchPos = std::max(lastKmerLeftMatchPos, (size_t)matches[i].matches[j].leftStart + (size_t)matches[i].matches[j].length);
	// 		lastKmerRightMatchPos = std::max(lastKmerRightMatchPos, (size_t)matches[i].matches[j].rightStart + (size_t)matches[i].matches[j].length);
	// 	}
	// 	if (matches[i].leftEnd - matches[i].leftStart > lastKmerLeftMatchPos + 500) remove = true;
	// 	if (matches[i].rightEnd - matches[i].rightStart > lastKmerRightMatchPos + 500) remove = true;
	// 	if (remove)
	// 	{
	// 		std::swap(matches[i], matches.back());
	// 		matches.pop_back();
	// 		filteredFromKmerGaps += 1;
	// 	}
	// 	else
	// 	{
	// 		countRemainingKmerMatches += matches[i].matches.size();
	// 	}
	// }
	std::cerr << filteredFromKmerGaps << " mapping matches filtered by k-mer gaps" << std::endl;
	std::cerr << matches.size() << " mapping matches" << std::endl;
	std::cerr << countRemainingKmerMatches << " kmer matches" << std::endl;
	return matches;
}

std::vector<MatchGroup> getMatches(const MatchIndex& matchIndex, const std::vector<TwobitString>& readSequences, size_t numThreads, const std::vector<size_t>& rawReadLengths, const size_t minAlignmentLength, const size_t k, const size_t graphk, const size_t graphd, const std::vector<bool>& usableChunkmers, const size_t maxIndexCoverage)
{
	std::vector<std::vector<ChainPosition>> tmp;
	tmp.resize(rawReadLengths.size());
	return getMatchesFilteredByAnnotation(matchIndex, readSequences, numThreads, rawReadLengths, minAlignmentLength, k, graphk, graphd, maxIndexCoverage, usableChunkmers, tmp);
}

std::vector<MatchGroup> filterMatchesByAnnotations(const std::vector<MatchGroup>& initialMatches, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<ChainPosition>>& chainAnnotations, const std::vector<std::vector<HaplotypeInformativeSite>>& haplotypeEdgeAnnotations)
{
	std::vector<MatchGroup> result;
	size_t filteredByChain = 0;
	size_t filteredByHaplotypeEdge = 0;
	for (size_t i = 0; i < initialMatches.size(); i++)
	{
		if (annotationsMismatch(chainAnnotations, rawReadLengths, initialMatches[i].leftRead, initialMatches[i].rightRead, initialMatches[i].leftStart, initialMatches[i].rightStart, initialMatches[i].rightFw))
		{
			filteredByChain += 1;
			continue;
		}
		if (annotationsMismatch(haplotypeEdgeAnnotations, rawReadLengths, initialMatches[i].leftRead, initialMatches[i].rightRead, initialMatches[i].leftStart, initialMatches[i].rightStart, initialMatches[i].rightFw))
		{
			filteredByHaplotypeEdge += 1;
			continue;
		}
		result.emplace_back(initialMatches[i]);
	}
	std::cerr << "filtered out " << filteredByChain << " matches by chain" << std::endl;
	std::cerr << "filtered out " << filteredByHaplotypeEdge << " matches by haplotype edges" << std::endl;
	return result;
}

int main(int argc, char** argv)
{
	size_t numThreads = std::stoi(argv[1]);
	size_t k = std::stoull(argv[2]);
	size_t numWindows = std::stoull(argv[3]);
	size_t windowSize = std::stoull(argv[4]);
	size_t minAlignmentLength = std::stoull(argv[5]);
	const size_t graphk = std::stoull(argv[6]);
	const size_t minCoverage = 2;
	const size_t graphd = 50;
	const double approxOneHapCoverage = 20;
	std::vector<std::string> readFiles;
	for (size_t i = 7; i < argc; i++)
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
		if (readBasepairLengths[i] + 1 > graphk)
		{
			readKmerLengths[i] = readBasepairLengths[i] + 1 - graphk;
		}
		else
		{
			readKmerLengths[i] = 1;
		}
	}
	const std::vector<std::string>& readNames = storage.getNames();
	for (size_t i = 0; i < readNames.size(); i++)
	{
		std::cerr << "readname " << i << " " << readNames[i] << std::endl;
	}
	auto rawReadLengths = storage.getRawReadLengths();
	std::vector<bool> usableChunkmers = getFilteredValidChunks(matchIndex, rawReadLengths, std::numeric_limits<size_t>::max(), 1000, 10);
	std::vector<MatchGroup> initialMatches = getMatches(matchIndex, readSequences, numThreads, rawReadLengths, minAlignmentLength, k, graphk, graphd, usableChunkmers, std::numeric_limits<size_t>::max());
	std::vector<MatchGroup> filteredMatches = initialMatches;
//	doHaplofilter(matches, readKmerLengths);
	// {
	// 	auto initial = getInitialGraph(readKmerLengths, readNames, readSequences, filteredMatches, minCoverage, graphk, approxOneHapCoverage, "initial1");
	// 	auto chains = getAnchorChains(initial.first, initial.second, approxOneHapCoverage);
	// 	auto readAnnotations = getReadChainPositions(initial.first, initial.second, chains);
	// 	auto haplotypeEdgeAnnotations = getHaplotypeInformativeForks(initial.first, initial.second, approxOneHapCoverage);
	// 	assert(haplotypeEdgeAnnotations.size() == readNames.size());
	// 	filteredMatches = filterMatchesByAnnotations(initialMatches, rawReadLengths, readAnnotations, haplotypeEdgeAnnotations);
	// 	std::cerr << filteredMatches.size() << " filtered matches out of " << initialMatches.size() << " initial matches" << std::endl;
	// }
	// {
	// 	auto initial = getInitialGraph(readKmerLengths, readNames, readSequences, filteredMatches, minCoverage, graphk, approxOneHapCoverage, "initial2");
	// 	auto chains = getAnchorChains(initial.first, initial.second, approxOneHapCoverage);
	// 	auto readAnnotations = getReadChainPositions(initial.first, initial.second, chains);
	// 	filteredMatches = filterMatchesByAnnotations(initialMatches, rawReadLengths, readAnnotations);
	// 	std::cerr << filteredMatches.size() << " filtered matches out of " << initialMatches.size() << " initial matches" << std::endl;
	// }
	// {
	// 	auto initial = getInitialGraph(readKmerLengths, readNames, readSequences, filteredMatches, minCoverage, graphk, approxOneHapCoverage, "initial3");
	// 	auto chains = getAnchorChains(initial.first, initial.second, approxOneHapCoverage);
	// 	auto readAnnotations = getReadChainPositions(initial.first, initial.second, chains);
	// 	filteredMatches = filterMatchesByAnnotations(initialMatches, rawReadLengths, readAnnotations);
	// 	std::cerr << filteredMatches.size() << " filtered matches out of " << initialMatches.size() << " initial matches" << std::endl;
	// }
	// {
	// 	auto initial = getInitialGraph(readKmerLengths, readNames, readSequences, filteredMatches, minCoverage, graphk, approxOneHapCoverage, "initial4");
	// 	auto chains = getAnchorChains(initial.first, initial.second, approxOneHapCoverage);
	// 	auto readAnnotations = getReadChainPositions(initial.first, initial.second, chains);
	// 	filteredMatches = filterMatchesByAnnotations(initialMatches, rawReadLengths, readAnnotations);
	// 	std::cerr << filteredMatches.size() << " filtered matches out of " << initialMatches.size() << " initial matches" << std::endl;
	// }
	makeGraph(readKmerLengths, readNames, readSequences, initialMatches, minCoverage, "graph.gfa", graphk, approxOneHapCoverage, numThreads);
}
