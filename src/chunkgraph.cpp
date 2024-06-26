#include <map>
#include <queue>
#include <mutex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "Common.h"
#include "RankBitvector.h"
#include "UnionFind.h"
#include "EdlibWrapper.h"
#include "CompressedStringIndex.h"
#include "TwobitString.h"
#include "ChunkUnitigGraph.h"
#include "ChunkGraphWriter.h"
#include "KmerIterator.h"
#include "ConsensusMaker.h"
#include "SequenceHelper.h"
#include "GraphCleaner.h"
#include "ChunkExtractor.h"
#include "SequenceIdentitySplitter.h"
#include "ChunkHelper.h"
#include "TransitiveClosure.h"
#include "ChunkPhasing.h"
#include "ChunkResolution.h"

double mismatchFraction;

void countReadRepetitiveUnitigs(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "counting read-repetitive unitigs" << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<size_t> chunkBelongsToUnitig;
	for (size_t i = 0; i < graph.chunksInUnitig.size(); i++)
	{
		for (uint64_t node : graph.chunksInUnitig[i])
		{
			size_t chunk = node & maskUint64_t;
			if (chunk >= chunkBelongsToUnitig.size()) chunkBelongsToUnitig.resize(chunk+1, std::numeric_limits<size_t>::max());
			chunkBelongsToUnitig[chunk] = i;
		}
	}
	std::vector<uint8_t> chunkRepeatCount;
	chunkRepeatCount.resize(chunkBelongsToUnitig.size(), 0);
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		phmap::flat_hash_set<size_t> foundHere;
		phmap::flat_hash_set<size_t> repetitiveHere;
		for (auto t : chunksPerRead[i])
		{
			size_t chunk = std::get<2>(t) & maskUint64_t;
			if (foundHere.count(chunk) == 1) repetitiveHere.insert(chunk);
			foundHere.insert(chunk);
		}
		for (size_t chunk : repetitiveHere)
		{
			if (chunk >= chunkRepeatCount.size()) continue;
			if (chunkRepeatCount[chunk] >= 2) continue;
			chunkRepeatCount[chunk] += 1;
		}
	}
	std::vector<bool> unitigRepetitive;
	unitigRepetitive.resize(graph.unitigLengths.size(), false);
	for (size_t i = 0; i < chunkRepeatCount.size(); i++)
	{
		if (chunkRepeatCount[i] < 2) continue;
		if (chunkBelongsToUnitig[i] == std::numeric_limits<size_t>::max()) continue;
		unitigRepetitive[chunkBelongsToUnitig[i]] = true;
	}
	size_t countRepetitive = 0;
	std::ofstream file { "repetitive.txt" };
	for (size_t i = 0; i < unitigRepetitive.size(); i++)
	{
		if (!unitigRepetitive[i]) continue;
		countRepetitive += 1;
		file << i << std::endl;
	}
	std::cerr << countRepetitive << " read-repetitive unitigs" << std::endl;
}

std::vector<std::pair<size_t, size_t>> getUniqueRanges(std::vector<std::pair<size_t, size_t>> raw)
{
	assert(raw.size() >= 1);
	std::sort(raw.begin(), raw.end());
	for (size_t j = raw.size()-1; j < raw.size() && j > 0; j--)
	{
		if (raw[j] != raw[j-1]) continue;
		std::swap(raw[j], raw.back());
		raw.pop_back();
	}
	assert(raw.size() >= 1);
	std::sort(raw.begin(), raw.end());
	return raw;
}

void expandChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize, const size_t expandedSize)
{
	std::vector<size_t> chunkLengths;
	{
		std::vector<std::vector<size_t>> rawChunkLengths;
		for (size_t i = 0; i < chunksPerRead.size(); i++)
		{
			for (size_t j = 0; j < chunksPerRead[i].size(); j++)
			{
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
				size_t chunk = std::get<2>(chunksPerRead[i][j]) & maskUint64_t;
				size_t length = std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) + 1;
				while (chunk >= rawChunkLengths.size()) rawChunkLengths.emplace_back();
				rawChunkLengths[chunk].emplace_back(length);
			}
		}
		chunkLengths.resize(rawChunkLengths.size());
		for (size_t i = 0; i < rawChunkLengths.size(); i++)
		{
			assert(rawChunkLengths[i].size() >= 1);
			std::sort(rawChunkLengths[i].begin(), rawChunkLengths[i].end());
			chunkLengths[i] = rawChunkLengths[i][rawChunkLengths[i].size()/2];
		}
	}
	std::map<std::vector<uint64_t>, size_t> chunkmerToNewChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::vector<std::pair<size_t, size_t>> newChunkmerRanges;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) newChunkmerRanges.emplace_back(j, j);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t end = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			while (end < chunksPerRead[i].size() && length < expandedSize)
			{
				end += 1;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][end]))) break;
				if (end == chunksPerRead[i].size()) break;
				length += chunkLengths[std::get<2>(chunksPerRead[i][end]) & maskUint64_t];
				length -= kmerSize;
			}
			if (length < expandedSize) continue;
			if (end == chunksPerRead[i].size()) continue;
			newChunkmerRanges.emplace_back(j, end);
		}
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t start = j;
			size_t length = chunkLengths[std::get<2>(chunksPerRead[i][j]) & maskUint64_t];
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			while (start > 0 && length < expandedSize)
			{
				start -= 1;
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][start]))) break;
				length += chunkLengths[std::get<2>(chunksPerRead[i][start]) & maskUint64_t];
				length -= kmerSize;
			}
			if (length < expandedSize) continue;
			newChunkmerRanges.emplace_back(start, j);
		}
		if (newChunkmerRanges.size() == 0)
		{
			chunksPerRead[i].clear();
			continue;
		}
		newChunkmerRanges = getUniqueRanges(newChunkmerRanges);
		assert(newChunkmerRanges.size() >= 1);
		for (size_t j = 1; j < newChunkmerRanges.size(); j++)
		{
			assert(newChunkmerRanges[j].first >= newChunkmerRanges[j-1].first);
			assert(newChunkmerRanges[j].second >= newChunkmerRanges[j-1].second);
		}
		std::vector<std::tuple<size_t, size_t, uint64_t>> newChunks;
		for (size_t j = 0; j < newChunkmerRanges.size(); j++)
		{
			size_t newchunk = std::numeric_limits<size_t>::max();
			std::vector<uint64_t> chunkmer;
			bool valid = true;
			for (size_t k = newChunkmerRanges[j].first; k <= newChunkmerRanges[j].second; k++)
			{
				chunkmer.emplace_back(std::get<2>(chunksPerRead[i][k]));
				if (NonexistantChunk(std::get<2>(chunksPerRead[i][k]))) valid = false;
			}
			if (valid)
			{
				if (chunkmerToNewChunk.count(chunkmer) == 1)
				{
					newchunk = chunkmerToNewChunk.at(chunkmer) + firstBitUint64_t;
				}
				else
				{
					newchunk = chunkmerToNewChunk.size();
					chunkmerToNewChunk[chunkmer] = newchunk;
					newchunk += firstBitUint64_t;
				}
			}
			else
			{
				newchunk = std::numeric_limits<size_t>::max();
			}
			newChunks.emplace_back(std::get<0>(chunksPerRead[i][newChunkmerRanges[j].first]), std::get<1>(chunksPerRead[i][newChunkmerRanges[j].second]), newchunk);
		}
		std::swap(chunksPerRead[i], newChunks);
	}
}

void makeGraph(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize, const size_t windowSize, const size_t middleSkip, const size_t startStage)
{
	std::cerr << "start at stage " << startStage << std::endl;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead;
	if (startStage > 0)
	{
		chunksPerRead = readChunksFromFakePathsFile("fakepaths" + std::to_string(startStage) + ".txt");
		while (chunksPerRead.size() < sequenceIndex.size()*2) chunksPerRead.emplace_back();
		std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
	switch(startStage)
	{
		case 0:
			std::cerr << "getting chunks from reads" << std::endl;
			chunksPerRead = getMinimizerBoundedChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize);
			//addMissingPiecesBetweenChunks(chunksPerRead, k);
			{
				size_t numChunks = 0;
				for (size_t i = 0; i < chunksPerRead.size(); i++)
				{
					numChunks += chunksPerRead[i].size();
				}
				std::cerr << numChunks << " chunks" << std::endl;
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
			splitPerFirstLastKmers(sequenceIndex, chunksPerRead, kmerSize);
			splitPerLength(chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(1, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 1:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerBaseCounts(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(2, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 2:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerMinHashes(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(3, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 3:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(4, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 4:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentityRoughly(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerSequenceIdentity(sequenceIndex, rawReadLengths, chunksPerRead, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 1000, 3000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
//			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			removeTinyProblemNodes(chunksPerRead, k);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(5, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 5:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			expandChunks(chunksPerRead, kmerSize, 5000);
			writeStage(51, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 51:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(6, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 6:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(7, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 7:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerNearestNeighborPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(8, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 8:
			splitPerCorrectedKmerPhasing(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(9, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 9:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(10, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 10:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(11, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 11:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(12, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 12:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(13, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 13:
//			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(14, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 14:
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, 2, kmerSize);
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.25, kmerSize);
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.5, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(15, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 15:
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(16, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 16:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
//			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(17, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 17:
//			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 31);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			writeStage(18, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 18:
			resolveTinyNodesRecklessly(chunksPerRead, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			countReadRepetitiveUnitigs(chunksPerRead, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadChunkSequences("sequences-chunk19.txt", rawReadLengths, chunksPerRead, sequenceIndex);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(19, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			writeBidirectedUnitigGraph("graph-dbg-final.gfa", "paths-dbg-final.gaf", "unitigs-dbg-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadUnitigSequences("sequences-dbg-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 19:
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 30000, 60000);
//			resolveBetweenLongUnitigs(chunksPerRead, rawReadLengths, numThreads, approxOneHapCoverage, 50000, 100000);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(20, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 20:
			writeBidirectedUnitigGraph("graph-final.gfa", "paths-final.gaf", "unitigs-final.fa", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeReadUnitigSequences("sequences-final.txt", chunksPerRead, sequenceIndex, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	}
}

int main(int argc, char** argv)
{
	const size_t numThreads = std::stoi(argv[1]);
	const size_t k = std::stoull(argv[2]);
	const size_t windowSize = std::stoull(argv[3]);
	const size_t middleSkip = std::stoull(argv[4]);
	const double approxOneHapCoverage = std::stod(argv[5]);
	mismatchFraction = std::stod(argv[6]); // try 2-3x average error rate
	const size_t startStage = std::stoull(argv[7]);
	std::vector<std::string> readFiles;
	for (size_t i = 8; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	FastaCompressor::CompressedStringIndex sequenceIndex { 5, 100 };
	std::vector<size_t> readBasepairLengths;
	for (auto file : readFiles)
	{
		std::cerr << "reading from file " << file << std::endl;
		readFilesAndAddToSequenceIndex(file, sequenceIndex, readBasepairLengths, numThreads);
	}
	sequenceIndex.removeConstructionVariables();
	size_t lastReal = readBasepairLengths.size();
	for (size_t i = 0; i < lastReal; i++)
	{
		readBasepairLengths.emplace_back(readBasepairLengths[i]);
	}
	std::cerr << sequenceIndex.size() << " reads" << std::endl;
	makeGraph(sequenceIndex, readBasepairLengths, numThreads, approxOneHapCoverage, k, windowSize, middleSkip, startStage);
}
