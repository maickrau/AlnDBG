#include <fstream>
#include <iostream>
#include <sstream>
#include "Common.h"
#include "ChunkGraphWriter.h"
#include "EdlibWrapper.h"
#include "KmerIterator.h"
#include "SequenceHelper.h"

void writeStage(const size_t stage, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize)
{
	writeGraph("tmpgraph" + std::to_string(stage) + ".gfa", "tmppaths" + std::to_string(stage) + ".txt", chunksPerRead);
	writeUnitigGraph("digraph-round" + std::to_string(stage) + ".gfa", "dipaths" + std::to_string(stage) + ".gaf", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
	writeBidirectedUnitigGraph("graph-round" + std::to_string(stage) + ".gfa", "paths" + std::to_string(stage) + ".gaf", chunksPerRead, sequenceIndex, rawReadLengths, approxOneHapCoverage, kmerSize);
}

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> readChunksFromTempPathsFile(const std::string& filename)
{
	std::cerr << "reading chunks from file " << filename << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> result;
	std::ifstream file { filename };
	if (!file.good()) throw FileCorruptedException {};
	size_t countChunks = 0;
	size_t countReads = 0;
	{
		std::string line;
		std::getline(file, line);
		std::stringstream sstr { line };
		sstr >> countReads;
	}
	result.resize(countReads);
	for (size_t i = 0; i < countReads; i++)
	{
		if (!file.good()) throw FileCorruptedException {};
		size_t countChunksHere = 0;
		{
			std::string line;
			std::getline(file, line);
			std::stringstream sstr { line };
			sstr >> countChunksHere;
		}
		if (!file.good()) throw FileCorruptedException {};
		result[i].reserve(countChunksHere);
		for (size_t j = 0; j < countChunksHere; j++)
		{
			std::string line;
			std::getline(file, line);
			if (line.size() == 0) throw FileCorruptedException {};
			if (!file.good()) throw FileCorruptedException {};
			std::stringstream sstr { line };
			size_t readindex = 0;
			size_t chunkindex = 0;
			size_t readstart = 0;
			size_t readend = 0;
			sstr >> readindex >> chunkindex >> readstart >> readend;
			if (readindex != i) throw FileCorruptedException {};
			if (chunkindex != result[i].size()) throw FileCorruptedException {};
			assert(readend > readstart);
			std::string nodestr;
			sstr >> nodestr;
			uint64_t node;
			assert(nodestr.size() >= 1);
			assert(nodestr.size() >= 2 || nodestr == "-");
			if (nodestr == "-")
			{
				node = std::numeric_limits<size_t>::max();
			}
			else if (nodestr[0] == '>')
			{
				node = std::stoull(nodestr.substr(1)) + firstBitUint64_t;
			}
			else
			{
				assert(nodestr[0] == '<');
				node = std::stoull(nodestr.substr(1));
			}
			result[i].emplace_back(readstart, readend, node);
			countChunks += 1;
		}
	}
	if (!file.good()) throw FileCorruptedException {};
	int dummy = 1;
	file >> dummy;
	if (dummy != -1) throw FileCorruptedException {};
	std::cerr << "read " << countChunks << " chunks in " << countReads << " fw/bw reads" << std::endl;
	return result;
}

void writeGraph(const std::string& filename, const std::vector<bool>& allowedNode, const SparseEdgeContainer& allowedEdges, const std::vector<std::vector<size_t>>& lengths, const std::vector<size_t>& coverages, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeCoverage, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& edgeOverlaps)
{
	std::ofstream graph { filename };
	for (size_t i = 0; i < coverages.size(); i++)
	{
		assert(lengths[i].size() >= coverages[i]);
		if (!allowedNode[i]) continue;
		size_t length = lengths[i][lengths[i].size()/2];
		size_t coverage = coverages[i];
		graph << "S\t" << i << "\t*\tLN:i:" << length << "\tll:i:" << coverage << "\tFC:i:" << length*coverage << std::endl;
	}
	for (auto pair : edgeCoverage)
	{
		if (!allowedEdges.hasEdge(std::make_pair(pair.first.first & maskUint64_t, pair.first.first & firstBitUint64_t), std::make_pair(pair.first.second & maskUint64_t, pair.first.second & firstBitUint64_t))) continue;
		assert(edgeOverlaps.count(pair.first) == 1);
		graph << "L\t" << (pair.first.first & maskUint64_t) << "\t" << ((pair.first.first & firstBitUint64_t) ? "+" : "-") << "\t" << (pair.first.second & maskUint64_t) << "\t" << ((pair.first.second & firstBitUint64_t) ? "+" : "-") << "\t" << edgeOverlaps.at(pair.first) << "M\tec:i:" << pair.second << std::endl;
	}
}

void writeUnitigGraph(const std::string& filename, const ChunkUnitigGraph& unitigGraph)
{
	assert(unitigGraph.unitigLengths.size() == unitigGraph.coverages.size());
	std::ofstream file { filename };
	for (size_t i = 0; i < unitigGraph.unitigLengths.size(); i++)
	{
		file << "S\t" << i << "\t*\tLN:i:" << unitigGraph.unitigLengths[i] << "\tll:f:" << unitigGraph.coverages[i] << "\tFC:f:" << (unitigGraph.unitigLengths[i] * unitigGraph.coverages[i]);
		file << "\tsg:Z:";
		for (size_t chunk : unitigGraph.chunksInUnitig[i])
		{
			file << ((chunk & firstBitUint64_t) ? ">" : "<") << (chunk & maskUint64_t);
		}
		file << "\tbr:Z:";
		for (size_t j = 0; j < unitigGraph.unitigChunkBreakpointPositions[i].size(); j++)
		{
			if (j > 0) file << ",";
			file << unitigGraph.unitigChunkBreakpointPositions[i][j].first << "-" << unitigGraph.unitigChunkBreakpointPositions[i][j].second;
		}
		file << std::endl;
	}
	for (auto pair : unitigGraph.edgeOverlaps)
	{
		assert(unitigGraph.edges.hasEdge(std::make_pair(pair.first.first & maskUint64_t, pair.first.first & firstBitUint64_t), std::make_pair(pair.first.second & maskUint64_t, pair.first.second & firstBitUint64_t)) == 1);
		assert(unitigGraph.edges.hasEdge(std::make_pair(pair.first.second & maskUint64_t, (pair.first.second ^ firstBitUint64_t) & firstBitUint64_t), std::make_pair(pair.first.first & maskUint64_t, (pair.first.first ^ firstBitUint64_t) & firstBitUint64_t)) == 1);
		file << "L\t" << (pair.first.first & maskUint64_t) << "\t" << ((pair.first.first & firstBitUint64_t) ? "+" : "-") << "\t" << (pair.first.second & maskUint64_t) << "\t" << ((pair.first.second & firstBitUint64_t) ? "+" : "-") << "\t" << pair.second << "M\tec:i:" << unitigGraph.edgeCoverages.at(pair.first) << std::endl; 
	}
}

void writeUnitigPaths(const std::string& filename, const ChunkUnitigGraph& unitigGraph, const std::vector<std::vector<UnitigPath>>& readPaths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths)
{
	std::ofstream file { filename };
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			std::string pathstr;
			size_t pathLength = 0;
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				uint64_t node = readPaths[i][j].path[k];
				pathstr += ((node & firstBitUint64_t) ? ">" : "<") + std::to_string(node & maskUint64_t);
				pathLength += unitigGraph.unitigLengths[node & maskUint64_t];
				if (k >= 1)
				{
					uint64_t prevNode = readPaths[i][j].path[k-1];
					auto pairkey = canon(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t), std::make_pair(node & maskUint64_t, node & firstBitUint64_t));
					std::pair<uint64_t, uint64_t> key { pairkey.first.first + (pairkey.first.second ? firstBitUint64_t : 0), pairkey.second.first + (pairkey.second.second ? firstBitUint64_t : 0) };
					pathLength -= unitigGraph.edgeOverlaps.at(key);
				}
			}
			std::string name;
			size_t length;
			if (i < sequenceIndex.size())
			{
				name = sequenceIndex.getName(i);
				length = rawReadLengths[i];
			}
			else
			{
				name = sequenceIndex.getName(i - sequenceIndex.size()) + "_bw";
				length = rawReadLengths[i - sequenceIndex.size()];
			}
			file << name << "\t" << length << "\t" << readPaths[i][j].readPartInPathnode[0].first << "\t" << readPaths[i][j].readPartInPathnode.back().second << "\t+\t" << pathstr << "\t" << pathLength << "\t" << readPaths[i][j].pathLeftClipBases << "\t" << (pathLength - readPaths[i][j].pathRightClipBases) << "\t" << (pathLength - readPaths[i][j].pathLeftClipBases - readPaths[i][j].pathRightClipBases) << "\t" << (pathLength - readPaths[i][j].pathLeftClipBases - readPaths[i][j].pathRightClipBases) << "\t60" << std::endl;
		}
	}
}

void writeUnitigGraph(const std::string& graphFile, const ChunkUnitigGraph& graph, const std::vector<TwobitString>& unitigDBGSequences, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& fixedOverlaps)
{
	std::ofstream file { graphFile };
	file << "H\tVN:Z:1" << std::endl;
	assert(graph.unitigLengths.size() == unitigDBGSequences.size());
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		size_t length = unitigDBGSequences[i].size();
		double coverage = graph.coverages[i];
		file << "S\t" << i << "\t";
		file << unitigDBGSequences[i].toString();
		file << "\tll:i:" << coverage << "\tFC:i:" << length*coverage << std::endl;
	}
	for (size_t i = 0; i < graph.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : graph.edges.getEdges(fw))
		{
			uint64_t from = i + firstBitUint64_t;
			uint64_t to = (edge.first) + (edge.second ? firstBitUint64_t : 0);
			auto key = canonNodePair(from, to);
			assert(fixedOverlaps.count(key) == 1);
			size_t overlapBases = fixedOverlaps.at(key);
			assert(graph.edgeCoverages.count(key) == 1);
			size_t coverage = graph.edgeCoverages.at(key);
			file << "L\t" << i << "\t+\t" << edge.first << "\t" << (edge.second ? "+" : "-") << "\t" << overlapBases << "M\tec:i:" << coverage << std::endl;
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : graph.edges.getEdges(bw))
		{
			uint64_t from = i;
			uint64_t to = (edge.first) + (edge.second ? firstBitUint64_t : 0);
			auto key = canonNodePair(from, to);
			assert(fixedOverlaps.count(key) == 1);
			size_t overlapBases = fixedOverlaps.at(key);
			assert(graph.edgeCoverages.count(key) == 1);
			size_t coverage = graph.edgeCoverages.at(key);
			file << "L\t" << i << "\t-\t" << edge.first << "\t" << (edge.second ? "+" : "-") << "\t" << overlapBases << "M\tec:i:" << coverage << std::endl;
		}
	}
}

void writeUnitigPaths(const std::string& pathsFile, const ChunkUnitigGraph& graph, const std::vector<std::vector<UnitigPath>>& readPaths, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<TwobitString>& unitigDBGSequences, const std::vector<size_t>& rawReadLengths, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& fixedOverlaps)
{
	std::ofstream file { pathsFile };
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			std::string pathstr;
			size_t pathLength = 0;
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				uint64_t node = readPaths[i][j].path[k];
				pathstr += ((node & firstBitUint64_t) ? ">" : "<") + std::to_string(node & maskUint64_t);
				pathLength += unitigDBGSequences[node & maskUint64_t].size();
				if (k >= 1)
				{
					auto key = canonNodePair(readPaths[i][j].path[k-1], readPaths[i][j].path[k]);
					assert(fixedOverlaps.count(key) == 1);
					pathLength -= fixedOverlaps.at(key);
				}
			}
			std::string name;
			size_t length;
			if (i < sequenceIndex.size())
			{
				name = sequenceIndex.getName(i);
				length = rawReadLengths[i];
			}
			else
			{
				name = sequenceIndex.getName(i - sequenceIndex.size()) + "_bw";
				length = rawReadLengths[i - sequenceIndex.size()];
			}
			size_t pathLeftClip = readPaths[i][j].pathLeftClipBases;
			size_t pathRightClip = readPaths[i][j].pathRightClipBases;
			file << name << "\t" << length << "\t" << readPaths[i][j].readPartInPathnode[0].first << "\t" << readPaths[i][j].readPartInPathnode.back().second << "\t+\t" << pathstr << "\t" << pathLength << "\t" << pathLeftClip << "\t" << (pathLength - pathRightClip) << "\t" << (pathLength - pathLeftClip - pathRightClip) << "\t" << (pathLength - pathLeftClip - pathRightClip) << "\t60" << std::endl;
		}
	}
}

void fixPathClips(std::vector<std::vector<UnitigPath>>& readPaths, const ChunkUnitigGraph& graph, const std::vector<TwobitString>& unitigDBGSequences, const std::vector<std::vector<std::pair<size_t, size_t>>>& chunkSequencePositionsWithinUnitigs)
{
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			size_t chunkLeftClip = readPaths[i][j].pathLeftClipChunks;
			size_t chunkRightClip = readPaths[i][j].pathRightClipChunks;
			assert(chunkLeftClip < chunkSequencePositionsWithinUnitigs[readPaths[i][j].path[0] & maskUint64_t].size());
			assert(chunkRightClip < chunkSequencePositionsWithinUnitigs[readPaths[i][j].path.back() & maskUint64_t].size());
			assert(readPaths[i][j].path.size() >= 2 || chunkLeftClip + chunkRightClip < chunkSequencePositionsWithinUnitigs[readPaths[i][j].path[0] & maskUint64_t].size());
			if (readPaths[i][j].path[0] & firstBitUint64_t)
			{
				readPaths[i][j].pathLeftClipBases = chunkSequencePositionsWithinUnitigs[readPaths[i][j].path[0] & maskUint64_t][chunkLeftClip].first;
			}
			else
			{
				readPaths[i][j].pathLeftClipBases = unitigDBGSequences[readPaths[i][j].path[0] & maskUint64_t].size() - chunkSequencePositionsWithinUnitigs[readPaths[i][j].path[0] & maskUint64_t][chunkSequencePositionsWithinUnitigs[readPaths[i][j].path[0] & maskUint64_t].size()-1-chunkLeftClip].second;
				assert(readPaths[i][j].pathLeftClipBases >= 1);
				readPaths[i][j].pathLeftClipBases -= 1;
			}
			if (readPaths[i][j].path.back() & firstBitUint64_t)
			{
				readPaths[i][j].pathRightClipBases = unitigDBGSequences[readPaths[i][j].path.back() & maskUint64_t].size() - chunkSequencePositionsWithinUnitigs[readPaths[i][j].path.back() & maskUint64_t][chunkSequencePositionsWithinUnitigs[readPaths[i][j].path.back() & maskUint64_t].size()-1-chunkRightClip].second;
				assert(readPaths[i][j].pathRightClipBases >= 1);
				readPaths[i][j].pathRightClipBases -= 1;
			}
			else
			{
				readPaths[i][j].pathRightClipBases = chunkSequencePositionsWithinUnitigs[readPaths[i][j].path.back() & maskUint64_t][chunkRightClip].first;
			}
		}
	}
}

bool trimUnitigsByOneBasePair(std::vector<std::vector<UnitigPath>>& readPaths, const ChunkUnitigGraph& graph, std::vector<TwobitString>& unitigDBGSequences, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& fixedOverlaps)
{
	std::vector<bool> shouldTrimForward;
	std::vector<bool> shouldTrimBackward;
	shouldTrimForward.resize(graph.unitigLengths.size(), false);
	shouldTrimBackward.resize(graph.unitigLengths.size(), false);
	bool trimAnything = false;
	for (size_t i = 0; i < graph.unitigLengths.size(); i++)
	{
		for (auto edge : graph.edges.getEdges(std::make_pair(i, true)))
		{
			auto key = canonNodePair(i + firstBitUint64_t, edge.first + (edge.second ? firstBitUint64_t : 0));
			assert(fixedOverlaps.count(key) == 1);
			size_t overlap = fixedOverlaps.at(key);
			if (overlap == unitigDBGSequences[edge.first].size())
			{
				shouldTrimForward[i] = true;
				trimAnything = true;
				break;
			}
		}
		for (auto edge : graph.edges.getEdges(std::make_pair(i, false)))
		{
			auto key = canonNodePair(i, edge.first + (edge.second ? firstBitUint64_t : 0));
			assert(fixedOverlaps.count(key) == 1);
			size_t overlap = fixedOverlaps.at(key);
			if (overlap == unitigDBGSequences[edge.first].size())
			{
				shouldTrimBackward[i] = true;
				trimAnything = true;
				break;
			}
		}
	}
	if (!trimAnything) return false;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			if (((readPaths[i][j].path[0] & firstBitUint64_t) && shouldTrimBackward[readPaths[i][j].path[0] & maskUint64_t]) || (((readPaths[i][j].path[0] ^ firstBitUint64_t) & firstBitUint64_t) && shouldTrimForward[readPaths[i][j].path[0] & maskUint64_t]))
			{
				if (readPaths[i][j].pathLeftClipBases == 0)
				{
					readPaths[i][j].readPartInPathnode[0].first += 1;
				}
				else
				{
					readPaths[i][j].pathLeftClipBases -= 1;
				}
			}
			if (((readPaths[i][j].path.back() & firstBitUint64_t) && shouldTrimForward[readPaths[i][j].path.back() & maskUint64_t]) || (((readPaths[i][j].path.back() ^ firstBitUint64_t) & firstBitUint64_t) && shouldTrimBackward[readPaths[i][j].path.back() & maskUint64_t]))
			{
				if (readPaths[i][j].pathRightClipBases == 0)
				{
					readPaths[i][j].readPartInPathnode.back().second -= 1;
				}
				else
				{
					readPaths[i][j].pathRightClipBases -= 1;
				}
			}
		}
	}
	for (size_t i = 0; i < unitigDBGSequences.size(); i++)
	{
		size_t start = 0;
		size_t size = unitigDBGSequences[i].size();
		if (shouldTrimBackward[i])
		{
			start += 1;
			size -= 1;
		}
		if (shouldTrimForward[i])
		{
			size -= 1;
		}
		if (size < unitigDBGSequences[i].size()) unitigDBGSequences[i] = unitigDBGSequences[i].substr(start, size);
	}
	for (auto& pair : fixedOverlaps)
	{
		assert(pair.second >= 2);
		if (((pair.first.first & firstBitUint64_t) && shouldTrimForward[pair.first.first & maskUint64_t]) || (((pair.first.first ^ firstBitUint64_t) & firstBitUint64_t) && shouldTrimBackward[pair.first.first & maskUint64_t]))
		{
			pair.second -= 1;
		}
		if (((pair.first.second & firstBitUint64_t) && shouldTrimBackward[pair.first.second & maskUint64_t]) || (((pair.first.second ^ firstBitUint64_t) & firstBitUint64_t) && shouldTrimForward[pair.first.second & maskUint64_t]))
		{
			pair.second -= 1;
		}
	}
	return true;
}

void writeBidirectedUnitigGraphWithSequences(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead, const std::vector<std::vector<size_t>>& minimizerPositionsPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t numThreads, const size_t kmerSize)
{
	std::cerr << "writing unitig graph with sequences " << graphFile << std::endl;
	auto fixedChunksPerRead = getBidirectedChunks(rawChunksPerRead);
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(fixedChunksPerRead, approxOneHapCoverage, kmerSize);
	std::vector<TwobitString> unitigDBGSequences;
	std::vector<std::vector<std::pair<size_t, size_t>>> chunkSequencePositionsWithinUnitigs;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> fixedOverlaps;
	std::tie(unitigDBGSequences, chunkSequencePositionsWithinUnitigs, fixedOverlaps) = getUnitigDBGSequences(graph, readPaths, kmerSize, fixedChunksPerRead, sequenceIndex, minimizerPositionsPerRead, numThreads);
	fixPathClips(readPaths, graph, unitigDBGSequences, chunkSequencePositionsWithinUnitigs);
	// ensure that edge overlap does not contain any unitigs
	for (size_t i = 0; i+1 < kmerSize; i++)
	{
		bool trimmed = trimUnitigsByOneBasePair(readPaths, graph, unitigDBGSequences, fixedOverlaps);
		if (!trimmed) break;
	}
	writeUnitigGraph(graphFile, graph, unitigDBGSequences, fixedOverlaps);
	std::cerr << "writing unitig paths " << pathsFile << std::endl;
	writeUnitigPaths(pathsFile, graph, readPaths, sequenceIndex, unitigDBGSequences, rawReadLengths, fixedOverlaps);
}

void writeBidirectedUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::string& unitigSequencesFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t numThreads, const size_t kmerSize)
{
	std::cerr << "writing unitig graph " << graphFile << std::endl;
	auto fixedChunksPerRead = getBidirectedChunks(rawChunksPerRead);
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::vector<ConsensusString> unitigSequences;
	std::tie(graph, readPaths, unitigSequences) = getChunkUnitigGraphWithConsensuses(fixedChunksPerRead, sequenceIndex, approxOneHapCoverage, numThreads, kmerSize);
	writeUnitigGraph(graphFile, graph);
	std::cerr << "writing unitig paths " << pathsFile << std::endl;
	writeUnitigPaths(pathsFile, graph, readPaths, sequenceIndex, rawReadLengths);
	std::cerr << "writing unitig sequences " << unitigSequencesFile << std::endl;
	writeUnitigSequences(unitigSequencesFile, unitigSequences);
}

void writeBidirectedUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "writing unitig graph " << graphFile << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	auto fixedChunksPerRead = getBidirectedChunks(chunksPerRead);
	std::tie(graph, readPaths) = getChunkUnitigGraph(fixedChunksPerRead, approxOneHapCoverage, kmerSize);
	writeUnitigGraph(graphFile, graph);
	std::cerr << "writing unitig paths " << pathsFile << std::endl;
	writeUnitigPaths(pathsFile, graph, readPaths, sequenceIndex, rawReadLengths);
}

void writeUnitigGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "writing unitig digraph " << graphFile << std::endl;
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	writeUnitigGraph(graphFile, graph);
	std::cerr << "writing unitig paths " << pathsFile << std::endl;
	writeUnitigPaths(pathsFile, graph, readPaths, sequenceIndex, rawReadLengths);
}

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getBidirectedChunks(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, uint64_t> pairToNode;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> fixedChunksPerRead;
	fixedChunksPerRead.resize(rawChunksPerRead.size()/2);
	for (size_t i = 0; i < rawChunksPerRead.size()/2; i++)
	{
		size_t other = i + rawChunksPerRead.size()/2;
		fixedChunksPerRead[i] = rawChunksPerRead[i];
		if (!(rawChunksPerRead[i].size() == rawChunksPerRead[other].size()))
		{
			std::cerr << i << " " << other << std::endl;
		}
		assert(rawChunksPerRead[i].size() == rawChunksPerRead[other].size());
		for (size_t j = 0; j < rawChunksPerRead[i].size(); j++)
		{
			if (NonexistantChunk(std::get<2>(rawChunksPerRead[i][j])))
			{
				if (!NonexistantChunk(std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j])))
				{
					std::cerr << i << " " << other << std::endl;
					std::cerr << j << " " << rawChunksPerRead[other].size()-1-j << std::endl;
					std::cerr << ((std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & maskUint64_t) << std::endl;
				}
				assert(NonexistantChunk(std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j])));
				continue;
			}
			if (NonexistantChunk(std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j])))
			{
				std::cerr << i << " " << other << std::endl;
				std::cerr << j << " " << rawChunksPerRead[other].size()-1-j << std::endl;
				std::cerr << ((std::get<2>(rawChunksPerRead[i][j]) & firstBitUint64_t) ? ">" : "<") << (std::get<2>(rawChunksPerRead[i][j]) & maskUint64_t) << std::endl;
			}
			assert(!NonexistantChunk(std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j])));
			assert((std::get<2>(rawChunksPerRead[i][j]) & firstBitUint64_t) == firstBitUint64_t);
			assert((std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & firstBitUint64_t) == firstBitUint64_t);
			if ((std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & maskUint64_t) < (std::get<2>(rawChunksPerRead[i][j]) & maskUint64_t))
			{
				std::get<2>(fixedChunksPerRead[i][j]) = (std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) & maskUint64_t);
			}
			/*
			uint64_t fwnode = std::get<2>(rawChunksPerRead[i][j]);
			uint64_t bwnode = std::get<2>(rawChunksPerRead[other][rawChunksPerRead[other].size()-1-j]) ^ firstBitUint64_t;
			std::pair<uint64_t, uint64_t> fwpair { fwnode, bwnode };
			std::pair<uint64_t, uint64_t> bwpair { bwnode ^ firstBitUint64_t, fwnode ^ firstBitUint64_t };
			bool fw = true;
			if (bwpair < fwpair)
			{
				std::swap(bwpair, fwpair);
				fw = false;
			}
			if (pairToNode.count(fwpair) == 0)
			{
				uint64_t next = pairToNode.size();
				pairToNode[fwpair] = next;
			}
			std::get<2>(fixedChunksPerRead[i][j]) = pairToNode.at(fwpair) + (fw ? firstBitUint64_t : 0);
			*/
		}
	}
	return fixedChunksPerRead;
}

void writeUnitigSequences(const std::string& filename, const std::vector<ConsensusString>& sequences)
{
	std::ofstream file { filename };
	for (size_t i = 0; i < sequences.size(); i++)
	{
		file << ">unitig_" << i << std::endl;
		std::string str = sequences[i].bases.toString();
		for (auto pair : sequences[i].Nchunks)
		{
			for (size_t j = pair.first; j < pair.first+pair.second; j++)
			{
				str[j] = 'N';
			}
		}
		file << str << std::endl;
	}
}

void writeGraph(const std::string& graphFile, const std::string& pathsFile, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::cerr << "writing paths " << pathsFile << std::endl;
	std::ofstream pathfile { pathsFile };
	pathfile << chunksPerRead.size() << std::endl;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		pathfile << chunksPerRead[i].size() << std::endl;
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			// if (j > 0 && std::get<0>(chunksPerRead[i][j]) == std::get<0>(chunksPerRead[i][j-1]) && std::get<1>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[i][j-1])) continue;
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t)))
			{
				pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << "-" << std::endl;
				continue;
			}
			uint64_t rawnode = std::get<2>(t);
			pathfile << i << " " << j << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << ((rawnode & firstBitUint64_t) ? ">" : "<") << (rawnode & maskUint64_t) << std::endl;
		}
	}
	pathfile << "-1";
	std::cerr << "writing graph " << graphFile << std::endl;
	std::vector<std::vector<size_t>> lengths;
	std::vector<size_t> coverages;
	std::tie(lengths, coverages) = getLengthsAndCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeCoverage = getEdgeCoverages(chunksPerRead);
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> edgeOverlaps = getEdgeOverlaps(chunksPerRead);
	std::vector<bool> allowedNode;
	SparseEdgeContainer allowedEdges;
	std::tie(allowedNode, allowedEdges) = getAllowedNodesAndEdges(coverages, edgeCoverage, chunksPerRead);
	writeGraph(graphFile, allowedNode, allowedEdges, lengths, coverages, edgeCoverage, edgeOverlaps);
}

void writeReadChunkSequences(const std::string& filename, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex)
{
	std::cerr << "writing read chunk sequences" << std::endl;
	std::ofstream file { filename };
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			size_t chunk = std::get<2>(t) & maskUint64_t;
			bool fw = std::get<2>(t) & firstBitUint64_t;
			std::string seq = getChunkSequence(sequenceIndex, rawReadLengths, chunksPerRead, i, j);
			std::string name;
			if (i < sequenceIndex.size())
			{
				name = sequenceIndex.getName(i);
			}
			else
			{
				name = sequenceIndex.getName(i - sequenceIndex.size()) + "_bw";
			}
			file << chunk << " " << name << " " << std::get<0>(t) << " " << std::get<1>(t) << " " << seq << std::endl;
		}
	}
}

void writeReadUnitigSequences(const std::string& filename, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& rawChunksPerRead, const FastaCompressor::CompressedStringIndex& sequenceIndex, const double approxOneHapCoverage, const size_t kmerSize)
{
	std::cerr << "writing read unitig sequences" << std::endl;
	auto chunksPerRead = getBidirectedChunks(rawChunksPerRead);
	ChunkUnitigGraph graph;
	std::vector<std::vector<UnitigPath>> readPaths;
	std::tie(graph, readPaths) = getChunkUnitigGraph(chunksPerRead, approxOneHapCoverage, kmerSize);
	std::ofstream file { filename };
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].size(); j++)
		{
			for (size_t k = 0; k < readPaths[i][j].path.size(); k++)
			{
				size_t unitig = readPaths[i][j].path[k] & maskUint64_t;
				bool fw = readPaths[i][j].path[k] & firstBitUint64_t;
				std::string seq = sequenceIndex.getSubstring(i, readPaths[i][j].readPartInPathnode[k].first, readPaths[i][j].readPartInPathnode[k].second - readPaths[i][j].readPartInPathnode[k].first+1);
				if (!fw) seq = revCompRaw(seq);
				size_t unitigStart = 0;
				size_t unitigEnd = graph.unitigLengths[unitig];
				if (k == 0) unitigStart += readPaths[i][j].pathLeftClipBases;
				if (k+1 == readPaths[i][j].path.size()) unitigEnd -= readPaths[i][j].pathRightClipBases;
				file << unitig << " " << unitigStart << " " << unitigEnd << " " << sequenceIndex.getName(i) << " " << readPaths[i][j].readPartInPathnode[k].first << " " << (readPaths[i][j].readPartInPathnode[k].second - readPaths[i][j].readPartInPathnode[k].first) << " " << seq << std::endl;
			}
		}
	}
}

void writeMinimizers(const std::string& filename, std::vector<std::vector<size_t>>& minimizerPositionsPerRead)
{
	std::ofstream file { filename };
	file << minimizerPositionsPerRead.size() << "\n";
	for (size_t i = 0; i < minimizerPositionsPerRead.size(); i++)
	{
		file << minimizerPositionsPerRead[i].size();
		for (size_t pos : minimizerPositionsPerRead[i])
		{
			file << " " << pos;
		}
		file << "\n";
	}
	file << "-1";
	assert(file.good());
}

std::vector<std::vector<size_t>> readMinimizersFromFile(const std::string& filename)
{
	std::ifstream file { filename };
	size_t count = 0;
	file >> count;
	std::vector<std::vector<size_t>> result;
	result.resize(count);
	for (size_t i = 0; i < result.size(); i++)
	{
		if (!file.good()) throw FileCorruptedException {};
		assert(file.good());
		file >> count;
		result[i].resize(count);
		for (size_t j = 0; j < result[i].size(); j++)
		{
			if (!file.good()) throw FileCorruptedException {};
			assert(file.good());
			file >> result[i][j];
		}
	}
	if (!file.good()) throw FileCorruptedException {};
	int dummy = 1;
	file >> dummy;
	if (dummy != -1) throw FileCorruptedException {};
	return result;
}
