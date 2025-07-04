#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cxxopts.hpp>
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
#include "OverlapMatcher.h"
#include "PathWalker.h"

const size_t RestartStageLatest = std::numeric_limits<size_t>::max()-1;

double mismatchFraction;

void makeGraph(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const double approxOneHapCoverage, const size_t kmerSize, const size_t windowSize, size_t startStage, const size_t resolveSize)
{
	if (startStage != RestartStageLatest)
	{
		std::cerr << "start at stage " << startStage << std::endl;
	}
	else
	{
		std::cerr << "restart from latest stage" << std::endl;
	}
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> chunksPerRead;
	std::vector<std::vector<size_t>> minimizerPositionsPerRead;
	if (startStage == RestartStageLatest)
	{
		startStage = 0;
		bool good = true;
		std::string minimizerFileName = "tmppaths_minimizers.txt";
		try
		{
			minimizerPositionsPerRead = readMinimizersFromFile(minimizerFileName);
		}
		catch (FileCorruptedException& e)
		{
			std::cerr << "Minimizer file " << minimizerFileName << " is corrupted or not found, restart from scratch" << std::endl;
			good = false;
		}
		if (good)
		{
			for (int i = 28; i >= 0; i--)
			{
				std::string chunkFilename = "tmppaths" + std::to_string(i) + ".txt";
				std::cerr << "try to read stage " << i << std::endl;
				try
				{
					chunksPerRead = readChunksFromTempPathsFile(chunkFilename);
					startStage = i;
					break;
				}
				catch (FileCorruptedException& e)
				{
				}
			}
			std::cerr << "start from stage " << startStage << std::endl;
		}
	}
	else if (startStage > 0)
	{
		std::string minimizerFileName = "tmppaths_minimizers.txt";
		std::string chunkFilename = "tmppaths" + std::to_string(startStage) + ".txt";
		try
		{
			chunksPerRead = readChunksFromTempPathsFile(chunkFilename);
		}
		catch (FileCorruptedException& e)
		{
			std::cerr << "Chunk file " << chunkFilename << " is corrupted" << std::endl;
			std::abort();
		}
		try
		{
			minimizerPositionsPerRead = readMinimizersFromFile(minimizerFileName);
		}
		catch (FileCorruptedException& e)
		{
			std::cerr << "Minimizer file " << minimizerFileName << " is corrupted" << std::endl;
			std::abort();
		}
		assert(chunksPerRead.size() == sequenceIndex.size()*2);
		assert(minimizerPositionsPerRead.size() == sequenceIndex.size()*2);
	}
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	switch(startStage)
	{
		case 0:
			std::cerr << "getting chunks from reads" << std::endl;
			std::tie(chunksPerRead, minimizerPositionsPerRead) = getMinimizerBoundedChunksPerRead(sequenceIndex, rawReadLengths, numThreads, kmerSize, windowSize);
			for (size_t i = 0; i < chunksPerRead.size(); i++)
			{
				for (size_t j = 0; j < chunksPerRead[i].size(); j++)
				{
					assert((std::get<2>(chunksPerRead[i][j]) & firstBitUint64_t) == 0);
					std::get<2>(chunksPerRead[i][j]) |= firstBitUint64_t;
				}
			}
			{
				size_t numChunks = 0;
				for (size_t i = 0; i < chunksPerRead.size(); i++)
				{
					numChunks += chunksPerRead[i].size();
				}
				std::cerr << numChunks << " chunks" << std::endl;
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			}
			writeStage(1, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			writeMinimizers("tmppaths_minimizers.txt", minimizerPositionsPerRead);
			[[fallthrough]];
		case 1:
			splitPerFirstLastKmers(sequenceIndex, chunksPerRead, kmerSize, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerLength(chunksPerRead, 0.02, 50, kmerSize, numThreads);
			splitPerLength(chunksPerRead, mismatchFraction, 10, kmerSize, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			removeBadShortHighCoverageChunks(chunksPerRead, kmerSize);
			mergeNonexistentChunks(chunksPerRead);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(2, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 2:
			splitPerBaseCounts(sequenceIndex, rawReadLengths, chunksPerRead, kmerSize, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(3, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 3:
			splitPerSequenceIdentityRoughly(sequenceIndex, rawReadLengths, chunksPerRead, kmerSize, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(4, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 4:
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(5, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 5:
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(6, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 6:
			splitPerSNPTransitiveClosureClustering(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, false);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(7, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 7:
			expandChunks(chunksPerRead, kmerSize, resolveSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			addGapChunks(chunksPerRead, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resegmentChunks(chunksPerRead, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(8, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 8:
			splitPerSNPTransitiveClosureClustering(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, true);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(9, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 9:
			splitPerAllelePhasingWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(10, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 10:
			splitPerPhasingKmersWithinChunk(sequenceIndex, rawReadLengths, chunksPerRead, 11, numThreads);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(11, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 11:
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(12, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 12:
			splitPerNeighborForksPolyploid(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(13, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 13:
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(14, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 14:
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, 2, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.25, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.5, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(15, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 15:
			splitPerNeighborForksPolyploid(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerInterchunkPhasedKmers(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(16, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 16:
			splitPerNeighborForksPolyploid(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(17, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 17:
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveUnambiguouslyResolvableUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			resolveSemiAmbiguousUnitigs(chunksPerRead, numThreads, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			splitPerNeighborForksPolyploid(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(18, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 18:
			{
				auto oldChunks = chunksPerRead;
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				fragmentChunks(chunksPerRead, minimizerPositionsPerRead, kmerSize);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				contextResolve(chunksPerRead, kmerSize, resolveSize);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				expandChunksUntilSolids(chunksPerRead, approxOneHapCoverage, kmerSize, resolveSize);
				std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
				resplitFalselyMergedChunks(chunksPerRead, oldChunks);
			}
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(19, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 19:
			expandResolvableChunks(chunksPerRead, approxOneHapCoverage, kmerSize, 50000);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(20, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 20:
			fixYForks(chunksPerRead, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(21, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 21:
			resolveBetweenTangles(chunksPerRead, approxOneHapCoverage, 50000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(22, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 22:
			resolveBetweenTanglesAllowGaps(chunksPerRead, approxOneHapCoverage, 50000, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(23, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 23:
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, 2, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.25, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			cleanTips(chunksPerRead, numThreads, approxOneHapCoverage, approxOneHapCoverage*0.5, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(24, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			[[fallthrough]];
		case 24:
			splitPerDiploidChunkWithNeighbors(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(25, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 25:
			splitPerNeighborForksPolyploid(sequenceIndex, rawReadLengths, chunksPerRead, numThreads, approxOneHapCoverage, 11);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(26, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 26:
			expandResolvableChunks(chunksPerRead, approxOneHapCoverage, kmerSize, 50000);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(27, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 27:
			fixYForks(chunksPerRead, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			writeStage(28, chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			[[fallthrough]];
		case 28:
			writeBidirectedUnitigGraphWithSequences("graph-resolved-final.gfa", "paths-resolved-final.gaf", chunksPerRead, minimizerPositionsPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, numThreads, kmerSize);
			std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
			break;
		default:
			std::cerr << "invalid start stage " << startStage << "!" << std::endl;
			std::abort();

	}
}

int main(int argc, char** argv)
{
	std::cerr << "chunkgraph " << VERSION << std::endl;
	cxxopts::Options options { "chunkgraph" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input reads. Multiple files can be input with -i file1.fa -i file2.fa etc (required)", cxxopts::value<std::vector<std::string>>())
		("t,threads", "Number of threads", cxxopts::value<size_t>()->default_value("1"))
		("k", "K-mer size", cxxopts::value<size_t>()->default_value("11"))
		("w", "Window size", cxxopts::value<size_t>()->default_value("5000"))
		("r", "Resolve size", cxxopts::value<size_t>()->default_value("10000"))
		("avg-hap-coverage", "Average single haplotype coverage", cxxopts::value<double>())
		("max-error-rate", "Maximum error rate. Try 2-3x average read error rate", cxxopts::value<double>()->default_value("0.03"))
		("restart", "Restart from stage", cxxopts::value<size_t>())
		("restart-latest", "Restart from latest stage")
		;
	auto params = options.parse(argc, argv);
	if (params.count("v") == 1)
	{
		std::cerr << "Version: " << VERSION << std::endl;
		exit(0);
	}
	if (params.count("h") == 1)
	{
		std::cerr << options.help();
		exit(0);
	}
	bool paramError = false;
	if (params.count("avg-hap-coverage") == 0)
	{
		std::cerr << "Average single haplotype coverage (--avg-hap-coverage) is required" << std::endl;
		paramError = true;
	}
	if (params.count("i") == 0)
	{
		std::cerr << "Input reads are required" << std::endl;
		paramError = true;
	}
	if (params.count("restart") == 1 && params.count("restart-latest") == 1)
	{
		std::cerr << "Cannot use both --restart and --restart-latest, pick only one" << std::endl;
		paramError = true;
	}
	if (paramError) std::abort();
	const size_t numThreads = params["t"].as<size_t>();
	const size_t k = params["k"].as<size_t>();
	const size_t windowSize = params["w"].as<size_t>();
	const size_t resolveSize = params["r"].as<size_t>();
	const double approxOneHapCoverage = params["avg-hap-coverage"].as<double>();
	mismatchFraction = params["max-error-rate"].as<double>();
	size_t startStage = 0;
	if (params.count("restart") == 1) startStage = params["restart"].as<size_t>();
	if (params.count("restart-latest") == 1) startStage = RestartStageLatest;
	std::vector<std::string> readFiles = params["i"].as<std::vector<std::string>>();
	FastaCompressor::CompressedStringIndex sequenceIndex { 5, 100 };
	std::vector<size_t> readBasepairLengths;
	for (auto file : readFiles)
	{
		std::cerr << "reading from file " << file << std::endl;
		readFilesAndAddToSequenceIndex(file, sequenceIndex, readBasepairLengths, numThreads);
	}
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "postprocessing sequence index" << std::endl;
	sequenceIndex.removeConstructionVariables(numThreads);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	size_t lastReal = readBasepairLengths.size();
	for (size_t i = 0; i < lastReal; i++)
	{
		readBasepairLengths.emplace_back(readBasepairLengths[i]);
	}
	std::cerr << sequenceIndex.size() << " reads" << std::endl;
	makeGraph(sequenceIndex, readBasepairLengths, numThreads, approxOneHapCoverage, k, windowSize, startStage, resolveSize);
}
