#include <iostream>
#include <tuple>
#include <vector>
#include <phmap.h>
#include "MBGCommon.h"
#include "GraphPhaser.h"
#include "Common.h"
#include "AnchorFinder.h"
#include "UnionFind.h"

class PhaseBlock
{
public:
	size_t chainNumber;
	std::vector<size_t> bubbleIndices;
	std::vector<std::vector<size_t>> allelesPerHaplotype;
	std::vector<std::vector<std::vector<size_t>>> extraAllelesPerHaplotype;
	bool chainStartPhased;
	bool chainEndPhased;
};

class ReadDiagonalAlleles
{
public:
	class AlleleMatch
	{
	public:
		size_t site;
		size_t allele;
		size_t readStartPos;
		size_t readEndPos;
	};
	size_t readId;
	bool fw;
	int diagonal;
	std::vector<AlleleMatch> alleles;
};

class Context
{
public:
	Context() :
		preHets(),
		overlapHets(),
		postHets(),
		preChains(),
		overlapChains(),
		postChains(),
		mergedhash(std::numeric_limits<size_t>::max())
	{
	}
	Context reverse() const
	{
		Context result = *this;
		result.mergedhash = std::numeric_limits<size_t>::max();
		std::swap(result.preHets, result.postHets);
		std::swap(result.preChains, result.postChains);
		std::reverse(result.preHets.begin(), result.preHets.end());
		std::reverse(result.overlapHets.begin(), result.overlapHets.end());
		std::reverse(result.postHets.begin(), result.postHets.end());
		std::reverse(result.preChains.begin(), result.preChains.end());
		std::reverse(result.overlapChains.begin(), result.overlapChains.end());
		std::reverse(result.postChains.begin(), result.postChains.end());
		for (size_t i = 0; i < result.preHets.size(); i++) result.preHets[i] ^= firstBitUint64_t;
		for (size_t i = 0; i < result.overlapHets.size(); i++) result.overlapHets[i] ^= firstBitUint64_t;
		for (size_t i = 0; i < result.postHets.size(); i++) result.postHets[i] ^= firstBitUint64_t;
		for (size_t i = 0; i < result.preChains.size(); i++) result.preChains[i] ^= firstBitUint64_t;
		for (size_t i = 0; i < result.overlapChains.size(); i++) result.overlapChains[i] ^= firstBitUint64_t;
		for (size_t i = 0; i < result.postChains.size(); i++) result.postChains[i] ^= firstBitUint64_t;
		return result;
	}
	std::vector<uint64_t> preHets;
	std::vector<uint64_t> overlapHets;
	std::vector<uint64_t> postHets;
	std::vector<uint64_t> preChains;
	std::vector<uint64_t> overlapChains;
	std::vector<uint64_t> postChains;
	bool operator!=(const Context& other) const
	{
		return !(*this == other);
	}
	bool operator==(const Context& other) const
	{
		return preHets == other.preHets && overlapHets == other.overlapHets && postHets == other.postHets && preChains == other.preChains && overlapChains == other.overlapChains && postChains == other.postChains;
	}
	bool operator<(const Context& other) const
	{
		if (preHets < other.preHets) return true;
		if (preHets > other.preHets) return false;
		if (overlapHets < other.overlapHets) return true;
		if (overlapHets > other.overlapHets) return false;
		if (postHets < other.postHets) return true;
		if (postHets > other.postHets) return false;
		if (preChains < other.preChains) return true;
		if (preChains > other.preChains) return false;
		if (overlapChains < other.overlapChains) return true;
		if (overlapChains > other.overlapChains) return false;
		if (postChains < other.postChains) return true;
		if (postChains > other.postChains) return false;
		return false;
	}
	size_t getHash() const
	{
		if (mergedhash == std::numeric_limits<size_t>::max())
		{
			mergedhash = 1;
			for (auto n : preHets) mergedhash = mergedhash * 3 + n;
			mergedhash *= 3;
			for (auto n : overlapHets) mergedhash = mergedhash * 3 + n;
			mergedhash *= 3;
			for (auto n : postHets) mergedhash = mergedhash * 3 + n;
			mergedhash *= 3;
			for (auto n : preChains) mergedhash = mergedhash * 3 + n;
			mergedhash *= 3;
			for (auto n : overlapChains) mergedhash = mergedhash * 3 + n;
			mergedhash *= 3;
			for (auto n : postChains) mergedhash = mergedhash * 3 + n;
			if (mergedhash == std::numeric_limits<size_t>::max()) mergedhash = std::numeric_limits<size_t>::max()-1;
		}
		return mergedhash;
	}
private:
	mutable size_t mergedhash;
};

namespace std
{
	template<>
	struct hash<Context>
	{
	public:
		size_t operator()(const Context& context) const
		{
			return context.getHash();
		}
	};
}

std::pair<uint64_t, uint64_t> canon(uint64_t from, uint64_t to)
{
	std::pair<uint64_t, uint64_t> fw { from, to };
	std::pair<uint64_t, uint64_t> bw { to ^ firstBitUint64_t, from ^ firstBitUint64_t };
	if (bw < fw) return bw;
	return fw;
}

size_t getAlleleIndex(std::vector<std::vector<std::vector<std::vector<uint64_t>>>>& allelesPerChain, const size_t i, const size_t j, const std::vector<uint64_t>& allele)
{
	for (size_t k = 0; k < allelesPerChain[i][j].size(); k++)
	{
		if (allele == allelesPerChain[i][j][k]) return k;
	}
	allelesPerChain[i][j].emplace_back(allele);
	return allelesPerChain[i][j].size()-1;
}

std::pair<std::vector<std::vector<std::vector<std::vector<uint64_t>>>>, std::vector<std::vector<ReadDiagonalAlleles>>> getRawReadInfoPerChain(const std::vector<AnchorChain>& anchorChains, const std::vector<ReadPathBundle>& readPaths, const UnitigGraph& unitigGraph, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads)
{
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<ReadDiagonalAlleles>> allelesPerReadPerChain;
	phmap::flat_hash_map<uint64_t, std::pair<size_t, size_t>> coreNodeLocator;
	allelesPerChain.resize(anchorChains.size());
	allelesPerReadPerChain.resize(anchorChains.size());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		allelesPerChain[i].resize(anchorChains[i].nodes.size()+1);
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			coreNodeLocator[anchorChains[i].nodes[j] & maskUint64_t] = std::make_pair(i, j);
		}
	}
	phmap::flat_hash_map<size_t, std::tuple<size_t, size_t, size_t>> simpleBubbleAllele;
	{
		SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
		for (size_t i = 0; i < anchorChains.size(); i++)
		{
			for (size_t j = 1; j < anchorChains[i].nodes.size(); j++)
			{
				const uint64_t prevNode = anchorChains[i].nodes[j-1];
				const uint64_t thisNode = anchorChains[i].nodes[j];
				if (edges.getEdges(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)).size() != 2) continue;
				if (edges.getEdges(std::make_pair(thisNode & maskUint64_t, (thisNode & firstBitUint64_t) ^ firstBitUint64_t)).size() != 2) continue;
				bool valid = true;
				for (auto edge : edges.getEdges(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)))
				{
					if (edges.getEdges(edge).size() != 1)
					{
						valid = false;
						break;
					}
					if (edges.getEdges(reverse(edge)).size() != 1)
					{
						valid = false;
						break;
					}
					if (edges.getEdges(edge)[0] != std::pair<size_t, bool>{ thisNode & maskUint64_t, thisNode & firstBitUint64_t })
					{
						valid = false;
						break;
					}
				}
				if (!valid) continue;
				// std::cerr << "found simple bubble between " << (prevNode & maskUint64_t) << " " << (thisNode & maskUint64_t) << std::endl;
				for (auto edge : edges.getEdges(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)))
				{
					std::vector<uint64_t> allele;
					allele.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
					size_t index = getAlleleIndex(allelesPerChain, i, j, allele);
					assert(simpleBubbleAllele.count(edge.first) == 0);
					simpleBubbleAllele[edge.first] = std::make_tuple(i, j, index);
				}
			}
		}
	}
	for (size_t readi = 0; readi < readPaths.size(); readi++)
	{
		std::vector<std::tuple<size_t, int, size_t, size_t, size_t, size_t>> allelesInThisRead; // chain, diagonal, bubble, allele, readStartPos, readEndPos
		// bubble gapless alleles
		for (size_t j = 0; j < readPaths[readi].paths.size(); j++)
		{
			size_t lastCore = std::numeric_limits<size_t>::max();
			size_t readPos = readPaths[readi].paths[j].readStartPos;
			size_t lastReadPos = std::numeric_limits<size_t>::max();
			for (size_t k = 0; k < readPaths[readi].paths[j].path.size(); k++)
			{
				readPos += unitigGraph.lengths[readPaths[readi].paths[j].path[k] & maskUint64_t];
				if (k == 0)
				{
					assert(readPos > readPaths[readi].paths[j].pathLeftClipKmers);
					readPos -= readPaths[readi].paths[j].pathLeftClipKmers;
				}
				if (coreNodeLocator.count(readPaths[readi].paths[j].path[k] & maskUint64_t) == 0) continue;
				if (lastCore == std::numeric_limits<size_t>::max())
				{
					lastCore = k;
					lastReadPos = readPos;
					continue;
				}
				std::pair<size_t, size_t> prev = coreNodeLocator.at(readPaths[readi].paths[j].path[lastCore] & maskUint64_t);
				std::pair<size_t, size_t> curr = coreNodeLocator.at(readPaths[readi].paths[j].path[k] & maskUint64_t);
				if (prev.first != curr.first)
				{
					lastCore = k;
					lastReadPos = readPos;
					continue;
				}
				if (prev.second != curr.second+1 && curr.second != prev.second+1)
				{
					lastCore = k;
					lastReadPos = readPos;
					continue;
				}
				bool fw = true;
				if (curr.second == prev.second + 1)
				{
					if (anchorChains[prev.first].nodes[prev.second] != readPaths[readi].paths[j].path[lastCore])
					{
						lastCore = k;
						lastReadPos = readPos;
						continue;
					}
					if (anchorChains[curr.first].nodes[curr.second] != readPaths[readi].paths[j].path[k])
					{
						lastCore = k;
						lastReadPos = readPos;
						continue;
					}
				}
				if (prev.second == curr.second + 1)
				{
					fw = false;
					if (anchorChains[prev.first].nodes[prev.second] != (readPaths[readi].paths[j].path[lastCore] ^ firstBitUint64_t))
					{
						lastCore = k;
						lastReadPos = readPos;
						continue;
					}
					if (anchorChains[curr.first].nodes[curr.second] != (readPaths[readi].paths[j].path[k] ^ firstBitUint64_t))
					{
						lastCore = k;
						lastReadPos = readPos;
						continue;
					}
				}
				std::vector<uint64_t> allele { readPaths[readi].paths[j].path.begin() + lastCore + 1, readPaths[readi].paths[j].path.begin() + k };
				if (!fw)
				{
					assert((readPaths[readi].paths[j].path[k] & firstBitUint64_t) != (anchorChains[curr.first].nodes[curr.second] & firstBitUint64_t));
					std::reverse(allele.begin(), allele.end());
					for (size_t l = 0; l < allele.size(); l++)
					{
						allele[l] ^= firstBitUint64_t;
					}
				}
				else
				{
					assert((readPaths[readi].paths[j].path[k] & firstBitUint64_t) == (anchorChains[curr.first].nodes[curr.second] & firstBitUint64_t));
				}
				size_t index = getAlleleIndex(allelesPerChain, curr.first, std::min(curr.second, prev.second)+1, allele);
				assert(curr.first < anchorChains.size());
				assert(curr.second < anchorChains[curr.first].nodeOffsets.size());
				int diagonal = (int)readPos - (int)anchorChains[curr.first].nodeOffsets[curr.second] - (int)unitigGraph.lengths[anchorChains[curr.first].nodes[curr.second] & maskUint64_t];
				if (!fw) diagonal = (int)readPos + (int)anchorChains[curr.first].nodeOffsets[curr.second];
				// std::cerr << "insert allele read " << readi << " chain " << (curr.first & maskUint64_t) << (fw ? "+" : "-") << " bubble " << (std::min(curr.second, prev.second)+1) << " allele " << index << " diagonal " << diagonal << std::endl;
				assert(diagonal > -(int)(anchorChains[curr.first].nodeOffsets.back()*2+readPos*2));
				assert(diagonal < (int)(anchorChains[curr.first].nodeOffsets.back()*2+readPos*2));
				assert(lastReadPos != std::numeric_limits<size_t>::max());
				allelesInThisRead.emplace_back(curr.first + (fw ? firstBitUint64_t : 0), diagonal, std::min(curr.second, prev.second)+1, index, lastReadPos, readPos-unitigGraph.lengths[anchorChains[curr.first].nodes[curr.second] & maskUint64_t]);
				lastCore = k;
				lastReadPos = readPos;
			}
		}
		// chain connection alleles
		for (size_t j = 1; j < chainPositionsInReads[readi].size(); j++)
		{
			std::pair<size_t, bool> lastChain = std::make_pair(chainPositionsInReads[readi][j-1].chain & maskUint64_t, chainPositionsInReads[readi][j-1].chain & firstBitUint64_t);
			std::pair<size_t, bool> currChain = std::make_pair(chainPositionsInReads[readi][j].chain & maskUint64_t, chainPositionsInReads[readi][j].chain & firstBitUint64_t);
			if (!validChainEdges.hasEdge(lastChain, currChain)) continue;
			bool prevFw = chainPositionsInReads[readi][j-1].chain & firstBitUint64_t;
			bool currFw = chainPositionsInReads[readi][j].chain & firstBitUint64_t;
			std::vector<uint64_t> prevAllele;
			prevAllele.emplace_back(chainPositionsInReads[readi][j].chain);
			size_t prevIndex = getAlleleIndex(allelesPerChain, chainPositionsInReads[readi][j-1].chain & maskUint64_t, prevFw ? anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodes.size() : 0, prevAllele);
			// std::cerr << "read " << readi << " chain-allele " << (lastChain.second ? ">" : "<") << lastChain.first << " to " << (currChain.second ? ">" : "<") << currChain.first << " allele " << prevIndex << std::endl;
			std::vector<uint64_t> currAllele;
			currAllele.emplace_back(chainPositionsInReads[readi][j-1].chain ^ firstBitUint64_t);
			size_t currIndex = getAlleleIndex(allelesPerChain, chainPositionsInReads[readi][j].chain & maskUint64_t, currFw ? 0 : anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.size(), currAllele);
			// std::cerr << "read " << readi << " chain-allele " << (currChain.second ? "<" : ">") << currChain.first << " to " << (lastChain.second ? "<" : ">") << lastChain.first << " allele " << currIndex << std::endl;
			int prevDiagonal, currDiagonal;
			if (prevFw)
			{
				prevDiagonal = chainPositionsInReads[readi][j-1].chainEndPosInRead - (int)anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodeOffsets.back() -  (int)unitigGraph.lengths[anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodes.back() & maskUint64_t];
			}
			else
			{
				prevDiagonal = chainPositionsInReads[readi][j-1].chainEndPosInRead;
			}
			if (currFw)
			{
				currDiagonal = chainPositionsInReads[readi][j].chainStartPosInRead;
			}
			else
			{
				currDiagonal = chainPositionsInReads[readi][j].chainStartPosInRead + (int)anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodeOffsets.back() + (int)unitigGraph.lengths[anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.back()];
			}
			// std::cerr << "tangleconnection insert allele read " << readi << " chain " << (chainPositionsInReads[readi][j-1].chain & maskUint64_t) << (prevFw ? "+" : "-") << " bubble " << (prevFw ? anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodes.size() : 0) << " allele " << prevIndex << " diagonal " << prevDiagonal << std::endl;
			// std::cerr << "tangleconnection insert allele read " << readi << " chain " << (chainPositionsInReads[readi][j].chain & maskUint64_t) << (currFw ? "+" : "-") << " bubble " << (currFw ? 0 : anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.size()) << " allele " << currIndex << " diagonal " << currDiagonal << std::endl;
			allelesInThisRead.emplace_back(chainPositionsInReads[readi][j-1].chain, prevDiagonal, prevFw ? anchorChains[chainPositionsInReads[readi][j-1].chain & maskUint64_t].nodes.size() : 0, prevIndex, chainPositionsInReads[readi][j-1].chainEndPosInRead, chainPositionsInReads[readi][j-1].chainEndPosInRead);
			allelesInThisRead.emplace_back(chainPositionsInReads[readi][j].chain, currDiagonal, currFw ? 0 : anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.size(), currIndex, chainPositionsInReads[readi][j].chainStartPosInRead, chainPositionsInReads[readi][j].chainStartPosInRead);
		}
		// simple bubble alleles with gaps
		uint64_t lastAnchor = std::numeric_limits<size_t>::max();
		int lastAnchorDiagonal = 0;
		size_t lastAnchorJ = 0;
		size_t lastSimpleBubble = std::numeric_limits<size_t>::max();
		size_t lastReadPos = std::numeric_limits<size_t>::max();
		std::tuple<size_t, size_t, size_t> uniqueLastSimpleBubbleAllele { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
		for (size_t j = 0; j < readPaths[readi].paths.size(); j++)
		{
			size_t readPos = readPaths[readi].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[readi].paths[j].path.size(); k++)
			{
				const uint64_t node = readPaths[readi].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[readi].paths[j].pathLeftClipKmers;
				if (coreNodeLocator.count(node & maskUint64_t) == 0)
				{
					if (simpleBubbleAllele.count(node & maskUint64_t) == 1 && node != lastSimpleBubble)
					{
						if (std::get<0>(uniqueLastSimpleBubbleAllele) != std::numeric_limits<size_t>::max() && uniqueLastSimpleBubbleAllele != simpleBubbleAllele.at(node & maskUint64_t))
						{
							uniqueLastSimpleBubbleAllele = std::make_tuple(std::numeric_limits<size_t>::max()-1, std::numeric_limits<size_t>::max()-1, std::numeric_limits<size_t>::max()-1);
						}
						else
						{
							uniqueLastSimpleBubbleAllele = simpleBubbleAllele.at(node & maskUint64_t);
						}
						lastSimpleBubble = node;
					}
					continue;
				}
				std::pair<size_t, size_t> pos = coreNodeLocator.at(node & maskUint64_t);
				assert((node & maskUint64_t) == (anchorChains[pos.first].nodes[pos.second] & maskUint64_t));
				bool fw = (node & firstBitUint64_t) == (anchorChains[pos.first].nodes[pos.second] & firstBitUint64_t);
				int diagonal = (int)readPos - (int)anchorChains[pos.first].nodeOffsets[pos.second] - (int)unitigGraph.lengths[anchorChains[pos.first].nodes[pos.second] & maskUint64_t];
				if (!fw) diagonal = (int)readPos + (int)anchorChains[pos.first].nodeOffsets[pos.second];
				if (lastAnchor == std::numeric_limits<size_t>::max() || lastAnchorJ == j || std::get<0>(uniqueLastSimpleBubbleAllele) != pos.first || (fw && std::get<1>(uniqueLastSimpleBubbleAllele) != pos.second) || (!fw && std::get<1>(uniqueLastSimpleBubbleAllele) != pos.second+1))
				{
					lastAnchor = node;
					lastAnchorJ = j;
					uniqueLastSimpleBubbleAllele = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
					lastAnchorDiagonal = diagonal;
					lastSimpleBubble = std::numeric_limits<size_t>::max();
					lastReadPos = readPos;
					continue;
				}
				std::pair<size_t, size_t> prevpos = coreNodeLocator.at(lastAnchor & maskUint64_t);
				bool prevfw = (lastAnchor & firstBitUint64_t) == (anchorChains[prevpos.first].nodes[prevpos.second] & firstBitUint64_t);
				int maxDiagonalDifference = (int)(anchorChains[pos.first].nodeOffsets.back() + unitigGraph.lengths[anchorChains[pos.first].nodes.back() & maskUint64_t])*0.5;
				if (fw != prevfw || prevpos.first != pos.first || (fw && prevpos.second+1 != pos.second) || (!fw && prevpos.second != pos.second+1) || diagonal > lastAnchorDiagonal + maxDiagonalDifference || diagonal < lastAnchorDiagonal - maxDiagonalDifference)
				{
					lastAnchor = node;
					lastAnchorJ = j;
					uniqueLastSimpleBubbleAllele = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
					lastAnchorDiagonal = diagonal;
					lastSimpleBubble = std::numeric_limits<size_t>::max();
					lastReadPos = readPos;
					continue;
				}
				assert(std::get<2>(uniqueLastSimpleBubbleAllele) < std::numeric_limits<size_t>::max()-1);
				// std::cerr << "simplebubble insert allele read " << readi << " chain " << (pos.first) << (fw ? "+" : "-") << " bubble " << (std::min(pos.second, prevpos.second)+1) << " allele " << std::get<2>(uniqueLastSimpleBubbleAllele) << " diagonal " << diagonal << std::endl;
				assert(lastReadPos != std::numeric_limits<size_t>::max());
				allelesInThisRead.emplace_back(pos.first + (fw ? firstBitUint64_t : 0), diagonal, std::min(pos.second, prevpos.second)+1, std::get<2>(uniqueLastSimpleBubbleAllele), lastReadPos, readPos-unitigGraph.lengths[anchorChains[pos.first].nodes[pos.second] & maskUint64_t]);
				lastAnchor = node;
				lastAnchorJ = j;
				uniqueLastSimpleBubbleAllele = std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
				lastAnchorDiagonal = diagonal;
				lastSimpleBubble = std::numeric_limits<size_t>::max();
				lastReadPos = readPos;
			}
		}
		std::sort(allelesInThisRead.begin(), allelesInThisRead.end());
		for (size_t i = 0; i < allelesInThisRead.size(); i++)
		{
			if (i == 0 || std::get<0>(allelesInThisRead[i]) != std::get<0>(allelesInThisRead[i-1]) || std::get<1>(allelesInThisRead[i]) > std::get<1>(allelesInThisRead[i-1]) + (int)(anchorChains[std::get<0>(allelesInThisRead[i])].nodeOffsets.back() + unitigGraph.lengths[anchorChains[std::get<0>(allelesInThisRead[i])].nodes.back() & maskUint64_t])/2)
			{
				allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].emplace_back();
				allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back().readId = readi;
				allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back().fw = std::get<0>(allelesInThisRead[i]) & firstBitUint64_t;
				allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back().diagonal = std::get<1>(allelesInThisRead[i]);
			}
			// std::cerr << "has allele chain " << (std::get<0>(allelesInThisRead[i]) & maskUint64_t) << ((std::get<0>(allelesInThisRead[i]) & firstBitUint64_t) ? "+" : "-") << " read " << readi << " diagonal " << std::get<1>(allelesInThisRead[i]) << " bubble " << std::get<2>(allelesInThisRead[i]) << " allele " << std::get<3>(allelesInThisRead[i]) << std::endl;
			allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back().alleles.emplace_back();
			allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back().alleles.back().site = std::get<2>(allelesInThisRead[i]);
			allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back().alleles.back().allele = std::get<3>(allelesInThisRead[i]);
			allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back().alleles.back().readStartPos = std::get<4>(allelesInThisRead[i]);
			allelesPerReadPerChain[std::get<0>(allelesInThisRead[i]) & maskUint64_t].back().alleles.back().readEndPos = std::get<5>(allelesInThisRead[i]);
		}
	}
	return std::make_pair(std::move(allelesPerChain), std::move(allelesPerReadPerChain));
}

std::vector<size_t> getConsensusAlleles(const std::vector<std::tuple<size_t, int, std::vector<std::pair<size_t, size_t>>>>& allelesPerRead, const size_t alleleCount)
{
	std::vector<phmap::flat_hash_map<size_t, size_t>> alleleCoverage;
	alleleCoverage.resize(alleleCount);
	for (size_t readi = 0; readi < allelesPerRead.size(); readi++)
	{
		for (auto pair : std::get<2>(allelesPerRead[readi]))
		{
			assert(pair.first < alleleCoverage.size());
			alleleCoverage[pair.first][pair.second] += 1;
		}
	}
	std::vector<size_t> result;
	for (size_t i = 0; i < alleleCount; i++)
	{
		std::pair<size_t, size_t> maxResult { std::numeric_limits<size_t>::max(), 0 };
		for (auto pair : alleleCoverage[i])
		{
			if (pair.second > maxResult.second) maxResult = pair;
		}
		result.push_back(maxResult.first);
	}
	return result;
}

std::pair<std::vector<std::vector<size_t>>, size_t> getTopNCoveredAlleles(const std::vector<std::vector<size_t>>& readsPerAllele, const size_t numAlleles)
{
	std::vector<std::pair<size_t, size_t>> coveragePerAllele;
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		coveragePerAllele.emplace_back(i, readsPerAllele[i].size());
	}
	std::sort(coveragePerAllele.begin(), coveragePerAllele.end(), [](auto left, auto right) { return left.second > right.second; });
	size_t numMinors = 0;
	for (size_t i = numAlleles; i < coveragePerAllele.size(); i++)
	{
		numMinors += coveragePerAllele[i].second;
	}
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < coveragePerAllele.size() && i < numAlleles; i++)
	{
		result.emplace_back(readsPerAllele[coveragePerAllele[i].first]);
		for (size_t j = 1; j < result.back().size(); j++)
		{
			assert(result.back()[j-1] < result.back()[j]);
		}
	}
	return std::make_pair(result, numMinors);
}

std::pair<std::vector<std::vector<size_t>>, size_t> getMostCoveredAlleles(const std::vector<std::vector<size_t>>& readsPerAllele, const size_t ploidy, const double approxOneHapCoverage)
{
	std::vector<std::pair<size_t, size_t>> coveragePerAllele;
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		coveragePerAllele.emplace_back(i, readsPerAllele[i].size());
	}
	std::sort(coveragePerAllele.begin(), coveragePerAllele.end(), [](auto left, auto right) { return left.second > right.second; });
	size_t estimatedTotalPloidy = 0;
	for (size_t i = 0; i < coveragePerAllele.size(); i++)
	{
		size_t ploidyHere = (double)coveragePerAllele[i].second / approxOneHapCoverage + 0.5;
		// std::cerr << i << " " << coveragePerAllele[i].second << " " << ploidyHere << std::endl;
		estimatedTotalPloidy += ploidyHere;
	}
	std::vector<std::vector<size_t>> result;
	// std::cerr << "estimated total ploidy " << estimatedTotalPloidy << std::endl;
	if (estimatedTotalPloidy < ploidy || estimatedTotalPloidy > ploidy+1) //loosely allow one more found ploidy just in case haplotype coverage hovers around 1.5x average
	{
		return std::make_pair(result, 0);
	}
	size_t minorCount = 0;
	while (coveragePerAllele.size() >= 1)
	{
		size_t lastPloidy = (double)coveragePerAllele.back().second / approxOneHapCoverage + 0.5;
		if (lastPloidy > 0) break;
		minorCount += readsPerAllele[coveragePerAllele.back().first].size();
		coveragePerAllele.pop_back();
	}
	size_t gotPloidy = 0;
	for (size_t i = 0; i < coveragePerAllele.size(); i++)
	{
		gotPloidy += (double)coveragePerAllele[i].second / approxOneHapCoverage + 0.5;
		result.emplace_back(readsPerAllele[coveragePerAllele[i].first]);
		std::sort(result.back().begin(), result.back().end());
	}
	assert(gotPloidy == estimatedTotalPloidy);
	return std::make_pair(result, minorCount);
}

bool allelesCouldSomewhatMatch(const std::vector<std::vector<size_t>>& bestAllelesLeft, const std::vector<std::vector<size_t>>& bestAllelesRight, const size_t ploidy, const size_t extraNonMatchers)
{
	if (bestAllelesLeft.size() != bestAllelesRight.size()) return false;
	if (bestAllelesLeft.size() != 2) return false;
	size_t requiredCoveragePerHap = 3;
	size_t normalCounts = intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[1]);
	size_t crossCounts = intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[0]);
	if (intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) >= requiredCoveragePerHap && intersectSize(bestAllelesLeft[1], bestAllelesRight[1]) >= requiredCoveragePerHap && normalCounts > (crossCounts+extraNonMatchers)*9 && crossCounts <= 2) return true;
	if (intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) >= requiredCoveragePerHap && intersectSize(bestAllelesLeft[1], bestAllelesRight[0]) >= requiredCoveragePerHap && crossCounts > (normalCounts+extraNonMatchers)*9 && normalCounts <= 2) return true;
	return false;
}

bool allelesCouldMatch(const std::vector<std::vector<size_t>>& bestAllelesLeft, const std::vector<std::vector<size_t>>& bestAllelesRight, const size_t ploidy, const double approxOneHapCoverage)
{
	if (bestAllelesLeft.size() != bestAllelesRight.size()) return false;
	if (bestAllelesLeft.size() != 2) return false;
	assert(bestAllelesLeft.size() >= 2);
	assert(bestAllelesRight.size() >= 2);
	size_t requiredCoverage = 4;
	// size_t requiredCoverage = (double)approxOneHapCoverage * ((double)ploidy - 0.5);
	size_t normalCounts = intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[1]);
	size_t crossCounts = intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) + intersectSize(bestAllelesLeft[1], bestAllelesRight[0]);
	if (normalCounts >= requiredCoverage && intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) >= 1 && intersectSize(bestAllelesLeft[1], bestAllelesRight[1]) >= 1 && normalCounts > crossCounts*9 && crossCounts <= 2) return true;
	if (crossCounts >= requiredCoverage && intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) >= 1 && intersectSize(bestAllelesLeft[1], bestAllelesRight[0]) >= 1 && crossCounts > normalCounts*9 && normalCounts <= 2) return true;
	return false;
}

std::vector<bool> getVeryLikelyHaplotypeInformativeSites(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& coreNodeChain, const size_t ploidy, const double approxOneHapCoverage)
{
	std::vector<std::vector<std::vector<size_t>>> bestAlleles;
	bestAlleles.resize(readsPerAllele.size());
	std::vector<bool> notgood;
	notgood.resize(readsPerAllele.size(), false);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		// std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << std::endl;
		size_t minorCount = 0;
		std::tie(bestAlleles[i], minorCount) = getMostCoveredAlleles(readsPerAllele[i], ploidy, approxOneHapCoverage);
		if (bestAlleles[i].size() == 0)
		{
			// std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << " best alleles zero, not good" << std::endl;
			notgood[i] = true;
			continue;
		}
		if (bestAlleles[i].size() < 2)
		{
			// std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << " best alleles less 2, not good" << std::endl;
			notgood[i] = true;
			continue;
		}
		size_t bestAllelePloidySum = 0;
		for (size_t j = 0; j < bestAlleles[i].size(); j++)
		{
			bestAllelePloidySum += (double)readsPerAllele[i][j].size() / approxOneHapCoverage + 0.5;
		}
		if (bestAllelePloidySum < ploidy)
		{
			// std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << " best alleles estimated ploidy less than wanted ploidy, not good" << std::endl;
			notgood[i] = true;
			continue;
		}
		size_t totalCoverage = 0;
		for (size_t j = 0; j < bestAlleles[i].size(); j++)
		{
			totalCoverage += bestAlleles[i][j].size();
		}
		if (minorCount*10 > totalCoverage)
		{
			// std::cerr << "bubble index " << i << "/" << readsPerAllele.size() << " minor allele coverage too high, not good" << std::endl;
			notgood[i] = true;
			continue;
		}
	}
	std::cerr << "first notgood: " << (notgood[0] ? "notgood" : "good") << std::endl;
	std::cerr << "last notgood: " << (notgood[0] ? "notgood" : "good") << std::endl;
	if (notgood[0])
	{
		size_t minorCount;
		std::tie(bestAlleles[0], minorCount) = getTopNCoveredAlleles(readsPerAllele[0], ploidy);
		if (minorCount == 0) notgood[0] = bestAlleles[0].size() != ploidy;
		std::cerr << "first minorcount: " << minorCount << std::endl;
		std::cerr << "first bestalleles size: " << bestAlleles[0].size() << " vs ploidy " << ploidy << std::endl;
	}
	if (notgood.back())
	{
		size_t minorCount;
		std::tie(bestAlleles.back(), minorCount) = getTopNCoveredAlleles(readsPerAllele.back(), ploidy);
		if (minorCount == 0) notgood.back() = bestAlleles.back().size() != ploidy;
	}
	std::cerr << "first notgood: " << (notgood[0] ? "notgood" : "good") << std::endl;
	std::cerr << "last notgood: " << (notgood[0] ? "notgood" : "good") << std::endl;
	std::vector<bool> result;
	result.resize(readsPerAllele.size(), false);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		if (notgood[i]) continue;
		for (size_t j = i+1; j < readsPerAllele.size(); j++)
		{
			if (notgood[j]) continue;
			if (allelesCouldMatch(bestAlleles[i], bestAlleles[j], ploidy, approxOneHapCoverage))
			{
				result[i] = true;
				result[j] = true;
			}
		}
	}
	return result;
}

std::vector<bool> extendMaybeHaplotypeInformativeSites(const std::vector<bool>& veryLikelyHaplotypeInformative, const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& coreNodeChain, const size_t ploidy, const double approxOneHapCoverage)
{
	std::vector<bool> result = veryLikelyHaplotypeInformative;
	std::vector<std::vector<std::vector<size_t>>> bestAlleles;
	std::vector<size_t> minorCount;
	bestAlleles.resize(readsPerAllele.size());
	minorCount.resize(readsPerAllele.size(), 0);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		std::tie(bestAlleles[i], minorCount[i]) = getTopNCoveredAlleles(readsPerAllele[i], ploidy);
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		if (result[i]) continue;
		for (size_t j = 0; j < veryLikelyHaplotypeInformative.size(); j++)
		{
			if (!veryLikelyHaplotypeInformative[j]) continue;
			if (allelesCouldSomewhatMatch(bestAlleles[i], bestAlleles[j], ploidy, minorCount[i]+minorCount[j]))
			{
				result[i] = true;
				break;
			}
		}
	}
	return result;
}

std::vector<bool> getPossiblyHaplotypeInformativeSites(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const std::vector<size_t>& coreNodeChain, const size_t ploidy, const double approxOneHapCoverage)
{
	std::vector<bool> veryLikelyHaplotypeInformative = getVeryLikelyHaplotypeInformativeSites(readsPerAllele, coreNodeChain, ploidy, approxOneHapCoverage);
	std::vector<bool> maybeHaplotypeInformative = extendMaybeHaplotypeInformativeSites(veryLikelyHaplotypeInformative, readsPerAllele, coreNodeChain, ploidy, approxOneHapCoverage);
	return maybeHaplotypeInformative;
}

std::vector<std::vector<size_t>> getAllelesPerHaplotype(const std::vector<size_t>& readAssignment, const size_t ploidy, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& allelesPerRead, const std::vector<size_t>& numAllelesPerSite)
{
	assert(readAssignment.size() <= allelesPerRead.size());
	assert(ploidy >= 2);
	std::vector<std::vector<std::vector<size_t>>> alleleCoveragePerHaplotype;
	alleleCoveragePerHaplotype.resize(ploidy);
	for (size_t hap = 0; hap < ploidy; hap++)
	{
		alleleCoveragePerHaplotype[hap].resize(numAllelesPerSite.size());
		for (size_t site = 0; site < numAllelesPerSite.size(); site++)
		{
			alleleCoveragePerHaplotype[hap][site].resize(numAllelesPerSite[site], 0);
		}
	}
	for (size_t i = 0; i < readAssignment.size(); i++)
	{
		if (readAssignment[i] == std::numeric_limits<size_t>::max()) continue;
		for (std::tuple<size_t, size_t, size_t> t : allelesPerRead[i])
		{
			size_t site = std::get<0>(t);
			size_t allele = std::get<1>(t);
			size_t weight = std::get<2>(t);
			assert(readAssignment[i] < alleleCoveragePerHaplotype.size());
			assert(site < alleleCoveragePerHaplotype[readAssignment[i]].size());
			assert(allele < alleleCoveragePerHaplotype[readAssignment[i]][site].size());
			alleleCoveragePerHaplotype[readAssignment[i]][site][allele] += weight;
		}
	}
	std::vector<std::vector<size_t>> result;
	result.resize(ploidy);
	for (size_t hap = 0; hap < ploidy; hap++)
	{
		result[hap].resize(numAllelesPerSite.size());
	}
	for (size_t site = 0; site < numAllelesPerSite.size(); site++)
	{
		size_t totalCoverageOfSite = 0;
		size_t totalMismatchesInSite = 0;
		for (size_t hap = 0; hap < ploidy; hap++)
		{
			std::pair<size_t, size_t> bestHere { std::numeric_limits<size_t>::max()-1, 0 };
			size_t totalCount = 0;
			for (size_t allele = 0; allele < alleleCoveragePerHaplotype[hap][site].size(); allele++)
			{
				if (alleleCoveragePerHaplotype[hap][site][allele] > bestHere.second)
				{
					bestHere.first = allele;
					bestHere.second = alleleCoveragePerHaplotype[hap][site][allele];
				}
				totalCount += alleleCoveragePerHaplotype[hap][site][allele];
			}
			totalMismatchesInSite += totalCount-bestHere.second;
			totalCoverageOfSite += totalCount;
			result[hap][site] = bestHere.first;
		}
		if (totalMismatchesInSite*10 > totalCoverageOfSite)
		{
			for (size_t hap = 0; hap < ploidy; hap++)
			{
				result[hap][site] = std::numeric_limits<size_t>::max();
			}
		}
	}
	return result;
}

size_t getHaplotypeScore(const std::vector<size_t>& readAssignment, const size_t ploidy, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& allelesPerRead, const std::vector<size_t>& numAllelesPerSite)
{
	assert(readAssignment.size() <= allelesPerRead.size());
	assert(ploidy >= 2);
	size_t score = 0;
	std::vector<std::vector<size_t>> allelesPerHaplotype = getAllelesPerHaplotype(readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
	for (size_t i = 0; i < readAssignment.size(); i++)
	{
		for (std::tuple<size_t, size_t, size_t> t : allelesPerRead[i])
		{
			size_t site = std::get<0>(t);
			size_t allele = std::get<1>(t);
			size_t weight = std::get<2>(t);
			if (readAssignment[i] == std::numeric_limits<size_t>::max())
			{
				score += weight;
				continue;
			}
			size_t haplotypeAllele = allelesPerHaplotype[readAssignment[i]][site];
			assert(haplotypeAllele != std::numeric_limits<size_t>::max()-1);
			if (haplotypeAllele == std::numeric_limits<size_t>::max())
			{
				score += weight;
				continue;
			}
			if (allele != haplotypeAllele)
			{
				score += weight*10;
			}
		}
	}
	return score;
}

class ReadPartition
{
public:
	std::vector<size_t> readAssignment;
	size_t maxAssignmentPlusOne;
	size_t score;
};

std::vector<size_t> getUnweightedHeuristicMEC(const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& allelesPerRead, const size_t ploidy, const std::vector<size_t>& numAllelesPerSite)
{
	const size_t beamWidth = 1000;
	size_t addEverythingUntilHere = log(beamWidth)/log(ploidy+1);
	assert(ploidy >= 2);
	std::vector<ReadPartition> activeAssignments;
	activeAssignments.emplace_back();
	activeAssignments.back().maxAssignmentPlusOne = 0;
	activeAssignments.back().score = 0;
	for (size_t i = 0; i < addEverythingUntilHere && i < allelesPerRead.size(); i++)
	{
		std::vector<ReadPartition> nextActiveAssignments;
		for (const auto& haplotype : activeAssignments)
		{
			for (size_t j = 0; j < std::min(haplotype.maxAssignmentPlusOne+1, ploidy); j++)
			{
				nextActiveAssignments.emplace_back();
				nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
				nextActiveAssignments.back().maxAssignmentPlusOne = std::max(haplotype.maxAssignmentPlusOne, j+1);
				nextActiveAssignments.back().readAssignment.push_back(j);
			}
			nextActiveAssignments.emplace_back();
			nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
			nextActiveAssignments.back().maxAssignmentPlusOne = haplotype.maxAssignmentPlusOne;
			nextActiveAssignments.back().readAssignment.push_back(std::numeric_limits<size_t>::max());
		}
		activeAssignments = nextActiveAssignments;
	}
	for (size_t i = 0; i < activeAssignments.size(); i++)
	{
		activeAssignments[i].score = getHaplotypeScore(activeAssignments[i].readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
	}
	std::sort(activeAssignments.begin(), activeAssignments.end(), [](const auto& left, const auto& right) { return left.score < right.score; });
	for (size_t i = addEverythingUntilHere; i < allelesPerRead.size(); i++)
	{
		std::vector<ReadPartition> nextActiveAssignments;
		for (const auto& haplotype : activeAssignments)
		{
			for (size_t j = 0; j < std::min(haplotype.maxAssignmentPlusOne+1, ploidy); j++)
			{
				nextActiveAssignments.emplace_back();
				nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
				nextActiveAssignments.back().readAssignment.push_back(j);
				nextActiveAssignments.back().maxAssignmentPlusOne = std::max(haplotype.maxAssignmentPlusOne, j+1);
				nextActiveAssignments.back().score = getHaplotypeScore(nextActiveAssignments.back().readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
			}
			nextActiveAssignments.emplace_back();
			nextActiveAssignments.back().readAssignment = haplotype.readAssignment;
			nextActiveAssignments.back().maxAssignmentPlusOne = haplotype.maxAssignmentPlusOne;
			nextActiveAssignments.back().readAssignment.push_back(std::numeric_limits<size_t>::max());
			nextActiveAssignments.back().score = getHaplotypeScore(nextActiveAssignments.back().readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
		}
		std::sort(nextActiveAssignments.begin(), nextActiveAssignments.end(), [](const auto& left, const auto& right) { return left.score < right.score; });
		if (nextActiveAssignments.size() > beamWidth)
		{
			nextActiveAssignments.erase(nextActiveAssignments.begin()+beamWidth, nextActiveAssignments.end());
		}
		activeAssignments = nextActiveAssignments;
	}
	assert(activeAssignments.size() >= 2);
	assert(activeAssignments[0].readAssignment.size() == allelesPerRead.size());
	std::cerr << "MEC best scores:";
	for (size_t i = 0; i < 10 && i < activeAssignments.size(); i++)
	{
		std::cerr << " " << activeAssignments[i].score;
	}
	std::cerr << std::endl;
	std::cerr << "best two assignments:" << std::endl;
	for (size_t i = 0; i < activeAssignments[0].readAssignment.size(); i++)
	{
		if (activeAssignments[0].readAssignment[i] == std::numeric_limits<size_t>::max())
		{
			std::cerr << "_";
		}
		else
		{
			std::cerr << activeAssignments[0].readAssignment[i];
		}
	}
	std::cerr << std::endl;
	for (size_t i = 0; i < activeAssignments[1].readAssignment.size(); i++)
	{
		if (activeAssignments[1].readAssignment[i] == std::numeric_limits<size_t>::max())
		{
			std::cerr << "_";
		}
		else
		{
			std::cerr << activeAssignments[1].readAssignment[i];
		}
	}
	std::cerr << std::endl;
	std::cerr << "best alleles:" << std::endl;
	std::vector<std::vector<size_t>> allelesPerHaplotype = getAllelesPerHaplotype(activeAssignments[0].readAssignment, ploidy, allelesPerRead, numAllelesPerSite);
	for (size_t i = 0; i < allelesPerHaplotype.size(); i++)
	{
		for (size_t j = 0; j < allelesPerHaplotype[i].size(); j++)
		{
			if (allelesPerHaplotype[i][j] == std::numeric_limits<size_t>::max())
			{
				std::cerr << "_,";
			}
			else if (allelesPerHaplotype[i][j] == std::numeric_limits<size_t>::max()-1)
			{
				std::cerr << ".,";
			}
			else
			{
				std::cerr << allelesPerHaplotype[i][j] << ",";
			}
		}
		std::cerr << std::endl;
	}
	return activeAssignments[0].readAssignment;
}

std::vector<size_t> getReadOrder(const std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>>& readsPerAllele, const size_t readCount)
{
	std::vector<size_t> firstPerRead;
	std::vector<size_t> lastPerRead;
	firstPerRead.resize(readCount, std::numeric_limits<size_t>::max());
	lastPerRead.resize(readCount, std::numeric_limits<size_t>::max());
	for (size_t site = 0; site < readsPerAllele.size(); site++)
	{
		for (size_t allele = 0; allele < readsPerAllele[site].size(); allele++)
		{
			for (std::pair<size_t, size_t> pair : readsPerAllele[site][allele])
			{
				size_t read = pair.first;
				if (firstPerRead[read] == std::numeric_limits<size_t>::max())
				{
					firstPerRead[read] = site;
				}
				lastPerRead[read] = site;
			}
		}
	}
	std::vector<std::tuple<size_t, size_t, size_t>> ordering;
	for (size_t i = 0; i < readCount; i++)
	{
		if (firstPerRead[i] == std::numeric_limits<size_t>::max()) continue;
		assert(lastPerRead[i] != std::numeric_limits<size_t>::max());
		ordering.emplace_back(i, firstPerRead[i], lastPerRead[i]);
	}
	std::sort(ordering.begin(), ordering.end(), [](auto left, auto right)
	{
		if (std::get<1>(left) < std::get<1>(right)) return true;
		if (std::get<1>(left) > std::get<1>(right)) return false;
		if (std::get<2>(left) < std::get<2>(right)) return true;
		if (std::get<2>(left) > std::get<2>(right)) return false;
		return false;
	});
	std::vector<size_t> result;
	result.resize(readCount, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < ordering.size(); i++)
	{
		assert(result[std::get<0>(ordering[i])] == std::numeric_limits<size_t>::max());
		result[std::get<0>(ordering[i])] = i;
	}
	std::vector<bool> found;
	found.resize(ordering.size(), false);
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		for (size_t j = 0; j < readsPerAllele[i].size(); j++)
		{
			for (std::pair<size_t, size_t> pair : readsPerAllele[i][j])
			{
				size_t read = pair.first;
				assert(result[read] != std::numeric_limits<size_t>::max());
				found[result[read]] = true;
			}
		}
	}
	for (size_t i = 0; i < found.size(); i++)
	{
		assert(found[i]);
	}
	return result;
}

std::vector<std::vector<std::vector<size_t>>> getInformativeReads(const std::vector<std::vector<std::vector<size_t>>>& unfiltered, const std::vector<bool>& maybeHaplotypeInformative)
{
	assert(unfiltered.size() == maybeHaplotypeInformative.size());
	phmap::flat_hash_map<size_t, size_t> readAlleleCount;
	for (size_t site = 0; site < maybeHaplotypeInformative.size(); site++)
	{
		if (!maybeHaplotypeInformative[site]) continue;
		for (size_t allele = 0; allele < unfiltered[site].size(); allele++)
		{
			for (auto read : unfiltered[site][allele])
			{
				readAlleleCount[read] += 1;
			}
		}
	}
	std::vector<std::vector<std::vector<size_t>>> result;
	for (size_t site = 0; site < maybeHaplotypeInformative.size(); site++)
	{
		if (!maybeHaplotypeInformative[site]) continue;
		result.emplace_back();
		for (size_t allele = 0; allele < unfiltered[site].size(); allele++)
		{
			result.back().emplace_back();
			for (auto read : unfiltered[site][allele])
			{
				if (readAlleleCount[read] < 2) continue;
				result.back().back().emplace_back(read);
			}
		}
	}
	return result;
}

std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> mergeReadsPerAllele(const std::vector<std::vector<std::vector<size_t>>>& readsPerAllele, const size_t readCount)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> allelesPerRead;
	allelesPerRead.resize(readCount);
	for (size_t site = 0; site < readsPerAllele.size(); site++)
	{
		for (size_t allele = 0; allele < readsPerAllele[site].size(); allele++)
		{
			for (size_t read : readsPerAllele[site][allele])
			{
				assert(read < allelesPerRead.size());
				allelesPerRead[read].emplace_back(site, allele);
			}
		}
	}
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		// todo is this necessary? they're probably already sorted due to readsPerAllele
		std::sort(allelesPerRead[i].begin(), allelesPerRead[i].end());
	}
	std::sort(allelesPerRead.begin(), allelesPerRead.end());
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> result;
	size_t currentReadCount = 1;
	for (size_t i = 1; i < allelesPerRead.size(); i++)
	{
		if (allelesPerRead[i] == allelesPerRead[i-1])
		{
			currentReadCount += 1;
			continue;
		}
		if (allelesPerRead[i-1].size() < 2)
		{
			currentReadCount = 1;
			continue;
		}
		result.emplace_back();
		for (auto pair : allelesPerRead[i-1])
		{
			result.back().emplace_back(pair.first, pair.second, currentReadCount);
		}
		currentReadCount = 1;
	}
	if (allelesPerRead.back().size() >= 2)
	{
		result.emplace_back();
		for (auto pair : allelesPerRead.back())
		{
			result.back().emplace_back(pair.first, pair.second, currentReadCount);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i].size() >= 1);
		std::sort(result[i].begin(), result[i].end());
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right)
	{
		if (std::get<0>(left[0]) < std::get<0>(right[0])) return true;
		if (std::get<0>(left[0]) > std::get<0>(right[0])) return false;
		if (std::get<0>(left.back()) < std::get<0>(right.back())) return true;
		if (std::get<0>(left.back()) > std::get<0>(right.back())) return false;
		if (left.size() > right.size()) return true;
		if (left.size() < right.size()) return false;
		for (size_t i = 0; i < left.size(); i++)
		{
			if (std::get<0>(left[i]) < std::get<0>(right[i])) return true;
			if (std::get<0>(left[i]) > std::get<0>(right[i])) return false;
		}
		if (std::get<2>(left[0]) > std::get<2>(right[0])) return true;
		if (std::get<2>(left[0]) < std::get<2>(right[0])) return false;
		return false;
	});
	return result;
}

std::vector<PhaseBlock> splitPhaseBlocks(const PhaseBlock& raw, const std::vector<size_t> readAssignment, const std::vector<std::vector<size_t>>& allelesPerHaplotype, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& mergedAllelesPerRead, const size_t ploidy)
{
	assert(ploidy >= 2);
	assert(allelesPerHaplotype.size() == ploidy);
	assert(mergedAllelesPerRead.size() == readAssignment.size());
	assert(raw.allelesPerHaplotype.size() == ploidy);
	assert(raw.bubbleIndices.size() == raw.allelesPerHaplotype[0].size());
	assert(raw.bubbleIndices.size() <= allelesPerHaplotype[0].size());
	assert(raw.bubbleIndices.size() >= 2);
	assert(allelesPerHaplotype[0].size() >= 2);
	std::vector<size_t> previousNonGarbageSite;
	previousNonGarbageSite.resize(allelesPerHaplotype[0].size(), std::numeric_limits<size_t>::max());
	size_t lastNonGarbage = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < allelesPerHaplotype[0].size(); i++)
	{
		if (allelesPerHaplotype[0][i] == std::numeric_limits<size_t>::max()) continue;
		previousNonGarbageSite[i] = lastNonGarbage;
		lastNonGarbage = i;
	}
	std::cerr << "pre-split phased at start: " << raw.chainStartPhased << " end: " << raw.chainEndPhased << std::endl;
	std::cerr << "bubble indices:";
	for (size_t i = 0; i < raw.bubbleIndices.size(); i++) std::cerr << " " << raw.bubbleIndices[i];
	std::cerr << std::endl;
	std::cerr << "number of sites: " << raw.bubbleIndices.size() << std::endl;
	{
		std::vector<bool> notGarbageSite;
		notGarbageSite.resize(allelesPerHaplotype[0].size(), true);
		for (size_t i = 0; i < allelesPerHaplotype[0].size(); i++)
		{
			if (allelesPerHaplotype[0][i] == std::numeric_limits<size_t>::max())
			{
				notGarbageSite[i] = false;
			}
		}
		for (size_t i = 0; i < notGarbageSite.size(); i++)
		{
			std::cerr << (notGarbageSite[i] ? "X" : "_");
		}
		std::cerr << " valid sites" << std::endl;
	}
	for (size_t i = 0; i < mergedAllelesPerRead.size(); i++)
	{
		std::vector<bool> hasAllele;
		hasAllele.resize(allelesPerHaplotype[0].size(), false);
		std::vector<size_t> allele;
		allele.resize(allelesPerHaplotype[0].size(), std::numeric_limits<size_t>::max());
		for (auto t : mergedAllelesPerRead[i])
		{
			assert(std::get<0>(t) < hasAllele.size());
			hasAllele[std::get<0>(t)] = true;
			if (allele[std::get<0>(t)] == std::numeric_limits<size_t>::max())
			{
				allele[std::get<0>(t)] = std::get<1>(t);
			}
			else if (allele[std::get<0>(t)] != std::get<1>(t))
			{
				allele[std::get<0>(t)] = std::numeric_limits<size_t>::max()-1;
			}
		}
		for (size_t i = 0; i < hasAllele.size(); i++)
		{
			if (!hasAllele[i])
			{
				std::cerr << "_";
			}
			else if (allele[i] <= 9)
			{
				std::cerr << allele[i];
			}
			else if (allele[i] == std::numeric_limits<size_t>::max()-1)
			{
				std::cerr << "X";
			}
			else
			{
				std::cerr << "A";
			}
		}
		std::cerr << " read " << i << " hap " << (readAssignment[i] < std::numeric_limits<size_t>::max() ? std::to_string(readAssignment[i]) : "-") << " weight " << std::get<2>(mergedAllelesPerRead[i][0]);
		if (readAssignment[i] == std::numeric_limits<size_t>::max()) std::cerr << " (GARBAGE READ)";
		std::cerr << std::endl;
	}
	std::vector<std::vector<bool>> haplotypeHasCrossingReads;
	haplotypeHasCrossingReads.resize(ploidy);
	for (size_t i = 0; i < ploidy; i++)
	{
		haplotypeHasCrossingReads[i].resize(allelesPerHaplotype[0].size(), false);
	}
	for (size_t i = 0; i < mergedAllelesPerRead.size(); i++)
	{
		if (readAssignment[i] == std::numeric_limits<size_t>::max()) continue;
		for (size_t j = 1; j < mergedAllelesPerRead[i].size(); j++)
		{
			if (allelesPerHaplotype[0][std::get<0>(mergedAllelesPerRead[i][j])] == std::numeric_limits<size_t>::max()) continue;
			for (size_t k = 0; k < j; k++)
			{
				if (allelesPerHaplotype[0][std::get<0>(mergedAllelesPerRead[i][k])] == std::numeric_limits<size_t>::max()) continue;
				if (std::get<0>(mergedAllelesPerRead[i][k]) == previousNonGarbageSite[std::get<0>(mergedAllelesPerRead[i][j])])
				{
					haplotypeHasCrossingReads[readAssignment[i]][std::get<0>(mergedAllelesPerRead[i][j])] = true;
				}
			}
		}
	}
	std::vector<PhaseBlock> result;
	result.emplace_back();
	result.back().chainNumber = raw.chainNumber;
	result.back().chainStartPhased = raw.chainStartPhased;
	result.back().allelesPerHaplotype.resize(ploidy);
	result.back().extraAllelesPerHaplotype.resize(ploidy);
	bool lastPhased = raw.chainEndPhased;
	size_t bubblei = 0;
	for (size_t i = 0; i < allelesPerHaplotype[0].size(); i++)
	{
		if (allelesPerHaplotype[0][i] == std::numeric_limits<size_t>::max())
		{
			if (i == 0) result.back().chainStartPhased = false;
			if (i == allelesPerHaplotype[0].size()-1) lastPhased = false;
			continue;
		}
		if (previousNonGarbageSite[i] == std::numeric_limits<size_t>::max())
		{
			result.back().bubbleIndices.emplace_back(raw.bubbleIndices[bubblei]);
			for (size_t k = 0; k < ploidy; k++)
			{
				assert(raw.allelesPerHaplotype[k][bubblei] != std::numeric_limits<size_t>::max());
				result.back().allelesPerHaplotype[k].push_back(raw.allelesPerHaplotype[k][bubblei]);
				result.back().extraAllelesPerHaplotype[k].push_back(raw.extraAllelesPerHaplotype[k][bubblei]);
			}
			// can't split before first site
			bubblei += 1;
			continue;
		}
		bool validSite = true;
		for (size_t k = 0; k < ploidy; k++)
		{
			if (!haplotypeHasCrossingReads[k][i])
			{
				validSite = false;
				break;
			}
		}
		if (!validSite)
		{
			std::cerr << "split block at chain " << raw.chainNumber << " index " << i << std::endl;
			result.back().chainEndPhased = false;
			result.emplace_back();
			result.back().chainStartPhased = false;
			result.back().chainNumber = raw.chainNumber;
			result.back().allelesPerHaplotype.resize(ploidy);
			result.back().extraAllelesPerHaplotype.resize(ploidy);
		}
		result.back().bubbleIndices.emplace_back(raw.bubbleIndices[bubblei]);
		for (size_t k = 0; k < ploidy; k++)
		{
			if (i == 0 && raw.allelesPerHaplotype[k][bubblei] == std::numeric_limits<size_t>::max()-1) result[0].chainStartPhased = false;
			if (i == allelesPerHaplotype[0].size()-1 && raw.allelesPerHaplotype[k][bubblei] == std::numeric_limits<size_t>::max()-1) lastPhased = false;
			assert(raw.allelesPerHaplotype[k][bubblei] != std::numeric_limits<size_t>::max());
			result.back().allelesPerHaplotype[k].push_back(raw.allelesPerHaplotype[k][bubblei]);
			result.back().extraAllelesPerHaplotype[k].push_back(raw.extraAllelesPerHaplotype[k][bubblei]);
		}
		bubblei += 1;
	}
	assert(bubblei == raw.bubbleIndices.size());
	result.back().chainEndPhased = lastPhased;
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		assert(result[i].bubbleIndices.size() == result[i].allelesPerHaplotype[0].size());
		if (result[i].bubbleIndices.size() < 2)
		{
			std::cerr << "remove too small block " << std::endl;
			std::swap(result[i], result.back());
			result.pop_back();
		}
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.bubbleIndices[0] < right.bubbleIndices[0]; });
	return result;
}

std::tuple<std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>,std::vector<bool>,std::vector<size_t>> getPhasingTableInformation(const size_t coreNodeChainIndex, const std::vector<ReadDiagonalAlleles>& allelesPerRead, const std::vector<uint64_t>& coreNodeChain, const size_t ploidy, const size_t alleleCount, const double approxOneHapCoverage)
{
	assert(ploidy >= 2);
	assert(alleleCount >= 1);
	size_t readCount = allelesPerRead.size();
	std::vector<std::vector<std::vector<size_t>>> unfilteredReadsPerAllele;
	unfilteredReadsPerAllele.resize(alleleCount);
	for (size_t readi = 0; readi < allelesPerRead.size(); readi++)
	{
		for (auto pair : allelesPerRead[readi].alleles)
		{
			assert(pair.site < unfilteredReadsPerAllele.size());
			while (pair.allele >= unfilteredReadsPerAllele[pair.site].size()) unfilteredReadsPerAllele[pair.site].emplace_back();
			assert(pair.allele < unfilteredReadsPerAllele[pair.site].size());
			unfilteredReadsPerAllele[pair.site][pair.allele].emplace_back(readi);
		}
	}
	for (size_t i = 0; i < unfilteredReadsPerAllele.size(); i++)
	{
		for (size_t j = 0; j < unfilteredReadsPerAllele[i].size(); j++)
		{
			phmap::flat_hash_set<size_t> unique { unfilteredReadsPerAllele[i][j].begin(), unfilteredReadsPerAllele[i][j].end() };
			unfilteredReadsPerAllele[i][j].clear();
			unfilteredReadsPerAllele[i][j].insert(unfilteredReadsPerAllele[i][j].end(), unique.begin(), unique.end());
			std::sort(unfilteredReadsPerAllele[i][j].begin(), unfilteredReadsPerAllele[i][j].end());
		}
	}
	std::cerr << "check chain " << coreNodeChainIndex << std::endl;
	std::vector<bool> maybeHaplotypeInformative = getPossiblyHaplotypeInformativeSites(unfilteredReadsPerAllele, coreNodeChain, ploidy, approxOneHapCoverage);
	std::cerr << "first informative: " << (maybeHaplotypeInformative[0] ? "yes" : "no") << " last: " << (maybeHaplotypeInformative.back() ? "yes" : "no") << std::endl;
	size_t countMaybeInformativeSites = 0;
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		if (maybeHaplotypeInformative[i]) countMaybeInformativeSites += 1;
	}
	if (countMaybeInformativeSites < 2)
	{
		std::cerr << "num informative sites " << countMaybeInformativeSites << ", skipped" << std::endl;
		return std::tuple<std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>,std::vector<bool>,std::vector<size_t>>{};
	}
	std::vector<std::vector<std::vector<size_t>>> readsPerAllele = getInformativeReads(unfilteredReadsPerAllele, maybeHaplotypeInformative);
	size_t informativeReadCount = 0;
	{
		phmap::flat_hash_set<size_t> foundRead;
		for (size_t i = 0; i < readsPerAllele.size(); i++)
		{
			for (size_t j = 0; j < readsPerAllele[i].size(); j++)
			{
				for (size_t read : readsPerAllele[i][j])
				{
					foundRead.insert(read);
				}
			}
		}
		informativeReadCount = foundRead.size();
	}
	std::cerr << informativeReadCount << " informative real reads" << std::endl;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> mergedAllelesPerRead = mergeReadsPerAllele(readsPerAllele, readCount);
	std::cerr << mergedAllelesPerRead.size() << " merged reads" << std::endl;
	std::vector<size_t> numAllelesPerSite;
	numAllelesPerSite.resize(readsPerAllele.size());
	for (size_t i = 0; i < readsPerAllele.size(); i++)
	{
		numAllelesPerSite[i] = readsPerAllele[i].size();
	}
	return std::make_tuple(std::move(mergedAllelesPerRead), std::move(maybeHaplotypeInformative), std::move(numAllelesPerSite));
}

std::vector<PhaseBlock> getChainPhaseBlocksMEC(const size_t coreNodeChainIndex, const std::vector<ReadDiagonalAlleles>& allelesPerRead, const std::vector<uint64_t>& coreNodeChain, const size_t ploidy, const size_t alleleCount, const double approxOneHapCoverage)
{
	assert(ploidy == 2);
	assert(alleleCount >= 1);
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> mergedAllelesPerRead;
	std::vector<bool> maybeHaplotypeInformative;
	std::vector<size_t> numAllelesPerSite;
	std::tie(mergedAllelesPerRead, maybeHaplotypeInformative, numAllelesPerSite) = getPhasingTableInformation(coreNodeChainIndex, allelesPerRead, coreNodeChain, ploidy, alleleCount, approxOneHapCoverage);
	if (maybeHaplotypeInformative.size() == 0) return std::vector<PhaseBlock> {};
	std::vector<size_t> readAssignments = getUnweightedHeuristicMEC(mergedAllelesPerRead, ploidy, numAllelesPerSite);
	std::vector<std::vector<size_t>> allelesPerHaplotype = getAllelesPerHaplotype(readAssignments, ploidy, mergedAllelesPerRead, numAllelesPerSite);
	PhaseBlock rawResult;
	rawResult.chainStartPhased = false;
	rawResult.chainEndPhased = false;
	if (maybeHaplotypeInformative[0]) rawResult.chainStartPhased = true;
	if (maybeHaplotypeInformative.back()) rawResult.chainEndPhased = true;
	rawResult.chainNumber = coreNodeChainIndex;
	for (size_t i = 0; i < ploidy; i++)
	{
		rawResult.allelesPerHaplotype.emplace_back();
		rawResult.extraAllelesPerHaplotype.emplace_back();
	}
	size_t bubblei = 0;
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		if (!maybeHaplotypeInformative[i]) continue;
		if (allelesPerHaplotype[0][bubblei] == std::numeric_limits<size_t>::max())
		{
			if (i == 0) rawResult.chainStartPhased = false;
			if (i == maybeHaplotypeInformative.size()-1) rawResult.chainEndPhased = false;
			bubblei += 1;
			continue;
		}
		rawResult.bubbleIndices.push_back(i);
		for (size_t j = 0; j < ploidy; j++)
		{
			assert(allelesPerHaplotype[j][bubblei] != std::numeric_limits<size_t>::max());
			rawResult.allelesPerHaplotype[j].push_back(allelesPerHaplotype[j][bubblei]);
			rawResult.extraAllelesPerHaplotype[j].emplace_back();
		}
		bubblei += 1;
	}
	assert(bubblei == allelesPerHaplotype[0].size());
	if (rawResult.bubbleIndices.size() < 2) return std::vector<PhaseBlock> {};
	return splitPhaseBlocks(rawResult, readAssignments, allelesPerHaplotype, mergedAllelesPerRead, ploidy);
}

std::vector<PhaseBlock> getChainPhaseBlocksSpanners(const size_t coreNodeChainIndex, const std::vector<ReadDiagonalAlleles>& allelesPerRead, const std::vector<uint64_t>& coreNodeChain, const size_t ploidy, const size_t alleleCount, const double approxOneHapCoverage)
{
	std::cerr << "try spanning read phase chain " << coreNodeChainIndex << std::endl;
	phmap::flat_hash_map<std::pair<size_t, size_t>, std::vector<size_t>> spannerReads;
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		size_t startAllele = std::numeric_limits<size_t>::max();
		size_t endAllele = std::numeric_limits<size_t>::max();
		for (auto pair : allelesPerRead[i].alleles)
		{
			if (pair.site == 0)
			{
				if (startAllele == std::numeric_limits<size_t>::max())
				{
					startAllele = pair.allele;
				}
				else if (startAllele != pair.allele)
				{
					startAllele = std::numeric_limits<size_t>::max()-1;
				}
			}
			if (pair.site == alleleCount-1)
			{
				if (endAllele == std::numeric_limits<size_t>::max())
				{
					endAllele = pair.allele;
				}
				else if (endAllele != pair.allele)
				{
					endAllele = std::numeric_limits<size_t>::max()-1;
				}
			}
		}
		if (startAllele >= std::numeric_limits<size_t>::max()-1) continue;
		if (endAllele >= std::numeric_limits<size_t>::max()-1) continue;
		spannerReads[std::make_pair(startAllele, endAllele)].push_back(i);
	}
	size_t estimatedPloidy = 0;
	size_t minorCoverage = 0;
	size_t effectivePloidy = 0;
	for (const auto& pair : spannerReads)
	{
		std::cerr << "spanner triple " << pair.first.first << " " << pair.first.second << " " << pair.second.size() << std::endl;
		size_t ploidy = pair.second.size() / approxOneHapCoverage + 0.5;
		estimatedPloidy += ploidy;
		if (ploidy == 0) minorCoverage += pair.second.size();
		if (ploidy >= 1) effectivePloidy += 1;
	}
	std::cerr << "spanner estimated ploidy sum " << estimatedPloidy << " minor count " << minorCoverage << " effective ploidy " << effectivePloidy << std::endl;
	if (estimatedPloidy != ploidy) return std::vector<PhaseBlock> {};
	if (minorCoverage >= 3) return std::vector<PhaseBlock> {};
	if (effectivePloidy == 1) return std::vector<PhaseBlock> {};
	PhaseBlock result;
	result.chainNumber = coreNodeChainIndex;
	result.bubbleIndices.push_back(0);
	result.chainStartPhased = true;
	result.chainEndPhased = true;
	result.bubbleIndices.push_back(alleleCount-1);
	std::cerr << "phased with spanners, alleles:" << std::endl;
	for (const auto& pair : spannerReads)
	{
		size_t ploidy = pair.second.size() / approxOneHapCoverage + 0.5;
		if (ploidy < 1) continue;
		std::cerr << "hap " << result.allelesPerHaplotype.size() << " " << pair.first.first << " " << pair.first.second << std::endl;
		result.allelesPerHaplotype.emplace_back();
		result.allelesPerHaplotype.back().push_back(pair.first.first);
		result.allelesPerHaplotype.back().push_back(pair.first.second);
		result.extraAllelesPerHaplotype.emplace_back();
		result.extraAllelesPerHaplotype.back().emplace_back();
		result.extraAllelesPerHaplotype.back().emplace_back();
	}
	return std::vector<PhaseBlock> { result };
}

std::vector<PhaseBlock> getChainPhaseBlocksTransitiveClosure(const size_t coreNodeChainIndex, const std::vector<ReadDiagonalAlleles>& allelesPerRead, const std::vector<uint64_t>& coreNodeChain, const size_t ploidy, const size_t alleleCount, const double approxOneHapCoverage)
{
	std::cerr << "try transitive closure phase chain " << coreNodeChainIndex << std::endl;
	assert(ploidy >= 2);
	assert(alleleCount >= 1);
	const size_t requiredMatchCoverage = approxOneHapCoverage * 0.5;
	const size_t requiredAmbiguousCoverage = std::min((size_t)(approxOneHapCoverage * 0.25), (size_t)3);
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> mergedAllelesPerRead;
	std::vector<bool> maybeHaplotypeInformative;
	std::vector<size_t> numAllelesPerSite;
	std::tie(mergedAllelesPerRead, maybeHaplotypeInformative, numAllelesPerSite) = getPhasingTableInformation(coreNodeChainIndex, allelesPerRead, coreNodeChain, ploidy, alleleCount, approxOneHapCoverage);
	if (maybeHaplotypeInformative.size() == 0) return std::vector<PhaseBlock> {};
	assert(maybeHaplotypeInformative.size() >= 2);
	std::cerr << "informative at start: " << (maybeHaplotypeInformative[0] ? "yes" : "no") << " end: " << (maybeHaplotypeInformative.back() ? "yes" : "no") << std::endl;
	std::cerr << "informative sites: ";
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		std::cerr << (maybeHaplotypeInformative[i] ? "1" : "_");
	}
	std::cerr << std::endl;
	if (!maybeHaplotypeInformative[0]) return std::vector<PhaseBlock> {};
	if (!maybeHaplotypeInformative.back()) return std::vector<PhaseBlock> {};
	std::vector<std::vector<std::pair<size_t, size_t>>> alleleParent;
	std::vector<std::vector<size_t>> alleleCoverage;
	alleleParent.resize(numAllelesPerSite.size());
	alleleCoverage.resize(numAllelesPerSite.size());
	for (size_t i = 0; i < numAllelesPerSite.size(); i++)
	{
		alleleCoverage[i].resize(numAllelesPerSite[i], 0);
		for (size_t j = 0; j < numAllelesPerSite[i]; j++)
		{
			alleleParent[i].emplace_back(i, j);
		}
	}
	for (size_t i = 0; i < mergedAllelesPerRead.size(); i++)
	{
		for (size_t j = 0; j < mergedAllelesPerRead[i].size(); j++)
		{
			assert(std::get<0>(mergedAllelesPerRead[i][j]) < alleleCoverage.size());
			assert(std::get<1>(mergedAllelesPerRead[i][j]) < alleleCoverage[std::get<0>(mergedAllelesPerRead[i][j])].size());
			alleleCoverage[std::get<0>(mergedAllelesPerRead[i][j])][std::get<1>(mergedAllelesPerRead[i][j])] += std::get<2>(mergedAllelesPerRead[i][j]);
		}
	}
	assert(alleleCoverage.size() == alleleParent.size());
	for (size_t i = 0; i < alleleCoverage.size(); i++)
	{
		assert(alleleCoverage[i].size() == alleleParent[i].size());
	}
	phmap::flat_hash_map<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>, size_t> alleleMatchCount;
	for (size_t i = 0; i < mergedAllelesPerRead.size(); i++)
	{
		for (size_t j = 0; j < mergedAllelesPerRead[i].size(); j++)
		{
			std::pair<size_t, size_t> from { std::get<0>(mergedAllelesPerRead[i][j]), std::get<1>(mergedAllelesPerRead[i][j]) };
			for (size_t k = j+1; k < mergedAllelesPerRead[i].size(); k++)
			{
				assert(std::get<0>(mergedAllelesPerRead[i][k]) >= std::get<0>(mergedAllelesPerRead[i][j]));
				assert(std::get<1>(mergedAllelesPerRead[i][k]) > std::get<1>(mergedAllelesPerRead[i][j]) || std::get<0>(mergedAllelesPerRead[i][k]) > std::get<0>(mergedAllelesPerRead[i][j]));
				std::pair<size_t, size_t> to { std::get<0>(mergedAllelesPerRead[i][k]), std::get<1>(mergedAllelesPerRead[i][k]) };
				assert(std::get<2>(mergedAllelesPerRead[i][j]) == std::get<2>(mergedAllelesPerRead[i][k]));
				alleleMatchCount[std::make_pair(from, to)] += std::get<2>(mergedAllelesPerRead[i][j]);
			}
		}
	}
	for (auto pair : alleleMatchCount)
	{
		if (pair.second < requiredMatchCoverage) continue;
		merge(alleleParent, pair.first.first, pair.first.second);
	}
	bool valid = true;
	for (auto pair : alleleMatchCount)
	{
		if (pair.second < requiredAmbiguousCoverage) continue;
		auto left = find(alleleParent, pair.first.first);
		auto right = find(alleleParent, pair.first.second);
		if (left != right)
		{
			valid = false;
			break;
		}
	}
	std::cerr << "valid: " << (valid ? "yes" : "no") << std::endl;
	if (!valid) return std::vector<PhaseBlock> {};
	phmap::flat_hash_set<std::pair<size_t, size_t>> validClusters;
	phmap::flat_hash_set<std::pair<size_t, size_t>> ambiguousClusters;
	for (size_t i = 0; i < alleleParent.size(); i++)
	{
		for (size_t j = 0; j < alleleParent[i].size(); j++)
		{
			if (alleleCoverage[i][j] >= requiredMatchCoverage) validClusters.insert(find(alleleParent, std::make_pair(i, j)));
			if (alleleCoverage[i][j] >= requiredAmbiguousCoverage) ambiguousClusters.insert(find(alleleParent, std::make_pair(i, j)));
		}
	}
	std::cerr << "valid clusters: " << validClusters.size() << std::endl;
	std::cerr << "ambiguous clusters: " << validClusters.size() << std::endl;
	if (validClusters.size() < 2) return std::vector<PhaseBlock> {};
	if (ambiguousClusters.size() != validClusters.size()) return std::vector<PhaseBlock> {};
	phmap::flat_hash_set<std::pair<size_t, size_t>> clustersAtStart;
	phmap::flat_hash_set<std::pair<size_t, size_t>> clustersAtEnd;
	for (size_t i = 0; i < alleleParent[0].size(); i++)
	{
		if (alleleCoverage[0][i] >= requiredAmbiguousCoverage) clustersAtStart.insert(find(alleleParent, std::make_pair(0, i)));
	}
	for (size_t i = 0; i < alleleParent.back().size(); i++)
	{
		if (alleleCoverage.back()[i] >= requiredAmbiguousCoverage) clustersAtEnd.insert(find(alleleParent, std::make_pair(alleleParent.size()-1, i)));
	}
	std::cerr << "clusters at start: " << clustersAtStart.size() << std::endl;
	std::cerr << "clusters at end: " << clustersAtEnd.size() << std::endl;
	if (clustersAtStart.size() != validClusters.size()) return std::vector<PhaseBlock> {};
	if (clustersAtEnd.size() != validClusters.size()) return std::vector<PhaseBlock> {};
	assert(clustersAtStart == clustersAtEnd);
	assert(clustersAtStart == validClusters);
	phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> clusterToHaplotype;
	for (auto pair : validClusters)
	{
		size_t hap = clusterToHaplotype.size();
		clusterToHaplotype[pair] = hap;
	}
	PhaseBlock result;
	result.chainNumber = coreNodeChainIndex;
	size_t effectivePloidy = validClusters.size();
	result.allelesPerHaplotype.resize(effectivePloidy);
	result.extraAllelesPerHaplotype.resize(effectivePloidy);
	assert(maybeHaplotypeInformative[0]);
	assert(maybeHaplotypeInformative.back());
	result.chainStartPhased = maybeHaplotypeInformative[0];
	result.chainEndPhased = maybeHaplotypeInformative.back();
	std::vector<size_t> bubbleToSite;
	bubbleToSite.resize(maybeHaplotypeInformative.size(), std::numeric_limits<size_t>::max());
	size_t bubbleIndex = 0;
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		if (!maybeHaplotypeInformative[i]) continue;
		result.bubbleIndices.emplace_back(i);
		for (size_t j = 0; j < effectivePloidy; j++)
		{
			result.allelesPerHaplotype[j].emplace_back(std::numeric_limits<size_t>::max());
			result.extraAllelesPerHaplotype[j].emplace_back();
		}
		bubbleToSite[i] = bubbleIndex;
		bubbleIndex += 1;
	}
	assert(bubbleIndex >= 2);
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		if (!maybeHaplotypeInformative[i]) continue;
		assert(bubbleToSite[i] < result.bubbleIndices.size());
		size_t site = bubbleToSite[i];
		assert(alleleCoverage[site].size() == numAllelesPerSite[site]);
		for (size_t j = 0; j < numAllelesPerSite[site]; j++)
		{
			if (alleleCoverage[site][j] < requiredAmbiguousCoverage) continue;
			size_t hap = clusterToHaplotype.at(find(alleleParent, std::make_pair(site, j)));
			result.extraAllelesPerHaplotype[hap][site].emplace_back(j);
		}
	}
	std::cerr << "yes phased" << std::endl;
	std::cerr << "effective ploidy: " << effectivePloidy << std::endl;
	std::cerr << "bubble indices:";
	for (size_t index : result.bubbleIndices) std::cerr << " " << index;
	std::cerr << std::endl;
	for (size_t hap = 0; hap < effectivePloidy; hap++)
	{
		assert(result.extraAllelesPerHaplotype[hap].size() == result.bubbleIndices.size());
		std::cerr << "hap " << hap << std::endl;
		for (size_t site = 0; site < result.extraAllelesPerHaplotype[hap].size(); site++)
		{
			for (size_t allele : result.extraAllelesPerHaplotype[hap][site])
			{
				std::cerr << "site " << site << " allele " << allele << std::endl;
			}
		}
	}
	return std::vector<PhaseBlock> { result };
}

bool isDiploidChain(const size_t chainNumber, const AnchorChain& chain, const std::vector<ReadDiagonalAlleles>& allelesPerRead, const double approxOneHapCoverage)
{
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> mergedAllelesPerRead;
	std::vector<bool> maybeHaplotypeInformative;
	std::vector<size_t> numAllelesPerSite;
	std::tie(mergedAllelesPerRead, maybeHaplotypeInformative, numAllelesPerSite) = getPhasingTableInformation(chainNumber, allelesPerRead, chain.nodes, 2, chain.nodes.size()+1, approxOneHapCoverage);
	if (maybeHaplotypeInformative.size() == 0) return false;
	size_t numHaplotypeInformative = 0;
	for (size_t i = 0; i < maybeHaplotypeInformative.size(); i++)
	{
		if (maybeHaplotypeInformative[i]) numHaplotypeInformative += 1;
	}
	return (maybeHaplotypeInformative[0] && maybeHaplotypeInformative.back()) && (numHaplotypeInformative >= 2);
}

double getNormalizedChainCoverage(const AnchorChain& anchorChain, const double approxOneHapCoverage, const UnitigGraph& unitigGraph)
{
	double sum = 0;
	double divisor = 0;
	for (auto node : anchorChain.nodes)
	{
		sum += unitigGraph.coverages[node & maskUint64_t];
		divisor += unitigGraph.lengths[node & maskUint64_t];
	}
	return sum / divisor;
}

std::vector<PhaseBlock> phaseDiploidChains(std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ReadDiagonalAlleles>>& allelesPerReadPerChain, const UnitigGraph& unitigGraph, const size_t numThreads, const double approxOneHapCoverage)
{
	assert(anchorChains.size() == allelesPerReadPerChain.size());
	std::vector<PhaseBlock> result;
	std::mutex resultMutex;
	iterateMultithreaded(0, anchorChains.size(), numThreads, [&resultMutex, &result, &anchorChains, &unitigGraph, &allelesPerReadPerChain, approxOneHapCoverage](size_t i)
	{
//		std::cerr << "check chain " << i << " start " << ((anchorChains[i].nodes[0] & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes[0] & maskUint64_t) << " end " << ((anchorChains[i].nodes.back() & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes.back() & maskUint64_t) << " ploidy " << anchorChains[i].ploidy << std::endl;
		if (anchorChains[i].ploidy == 1)
		{
//			std::cerr << "check if chain " << i << " is low coverage miscalled diploid" << std::endl;
			if (isDiploidChain(i, anchorChains[i], allelesPerReadPerChain[i], approxOneHapCoverage))
			{
//				std::cerr << "chain " << i << " is actually diploid" << std::endl;
				anchorChains[i].ploidy = 2;
			}
			else
			{
				return;
			}
		}
		bool checkingDiploidFromThree = false;
		if (anchorChains[i].ploidy == 3)
		{
			double fractionalCopyCount = getNormalizedChainCoverage(anchorChains[i], approxOneHapCoverage, unitigGraph);
			if (fractionalCopyCount < 3 && isDiploidChain(i, anchorChains[i], allelesPerReadPerChain[i], approxOneHapCoverage))
			{
//				std::cerr << "check if chain " << i << " is high coverage miscalled diploid" << std::endl;
				anchorChains[i].ploidy = 2;
				checkingDiploidFromThree = true;
			}
		}
		if (anchorChains[i].ploidy == 2)
		{
//			std::cerr << "try phase diploid chain " << i << " start " << ((anchorChains[i].nodes[0] & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes[0] & maskUint64_t) << " end " << ((anchorChains[i].nodes.back() & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes.back() & maskUint64_t) << std::endl;
			auto phaseBlocks = getChainPhaseBlocksMEC(i, allelesPerReadPerChain[i], anchorChains[i].nodes, anchorChains[i].ploidy, anchorChains[i].nodes.size()+1, approxOneHapCoverage);
			if (phaseBlocks.size() == 0) return;
//			std::cerr << "phased with " << phaseBlocks.size() << " blocks. start: " << (phaseBlocks[0].chainStartPhased ? "yes" : "no") << ", end: " << (phaseBlocks.back().chainEndPhased ? "yes" : "no") << std::endl;
			bool include = true;
			for (size_t j = 0; j < phaseBlocks.size(); j++)
			{
				for (size_t sitei = 0; sitei < phaseBlocks[j].allelesPerHaplotype[0].size(); sitei++)
				{
					phmap::flat_hash_set<size_t> foundAlleles;
					for (size_t k = 0; k < anchorChains[i].ploidy; k++)
					{
						foundAlleles.insert(phaseBlocks[j].allelesPerHaplotype[k][sitei]);
					}
					if (foundAlleles.size() != anchorChains[i].ploidy)
					{
//						std::cerr << "red flag: allele shared between haps, ploidy " << anchorChains[i].ploidy << ", num distinct alleles: " << foundAlleles.size() << ", skipping this phasing" << std::endl;
						include = false;
					}
				}
			}
			if (!include && checkingDiploidFromThree) anchorChains[i].ploidy = 3;
			if (include)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				result.insert(result.end(), phaseBlocks.begin(), phaseBlocks.end());
			}
		}
		else
		{
//			std::cerr << "skipped chain with ploidy " << anchorChains[i].ploidy << std::endl;
		}
	});
	return result;
}

std::vector<PhaseBlock> phasePolyploidsTransitiveClosure(std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ReadDiagonalAlleles>>& allelesPerReadPerChain, const double approxOneHapCoverage)
{
	assert(anchorChains.size() == allelesPerReadPerChain.size());
	std::vector<PhaseBlock> result;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::cerr << "check chain " << i << " start " << ((anchorChains[i].nodes[0] & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes[0] & maskUint64_t) << " end " << ((anchorChains[i].nodes.back() & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes.back() & maskUint64_t) << " ploidy " << anchorChains[i].ploidy << std::endl;
		if (anchorChains[i].ploidy == 0) continue;
		if (anchorChains[i].ploidy == 1)
		{
			std::cerr << "check if chain " << i << " is low coverage miscalled diploid" << std::endl;
			if (isDiploidChain(i, anchorChains[i], allelesPerReadPerChain[i], approxOneHapCoverage))
			{
				std::cerr << "chain " << i << " is actually diploid" << std::endl;
				anchorChains[i].ploidy = 2;
			}
			else
			{
				continue;
			}
		}
		std::cerr << "try phase ploidy " << anchorChains[i].ploidy << " chain " << i << " start " << ((anchorChains[i].nodes[0] & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes[0] & maskUint64_t) << " end " << ((anchorChains[i].nodes.back() & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes.back() & maskUint64_t) << std::endl;
		auto phaseBlocks = getChainPhaseBlocksTransitiveClosure(i, allelesPerReadPerChain[i], anchorChains[i].nodes, anchorChains[i].ploidy, anchorChains[i].nodes.size()+1, approxOneHapCoverage);
		if (phaseBlocks.size() == 0) continue;
		if (phaseBlocks.size() == 1 && phaseBlocks[0].allelesPerHaplotype.size() == 1) continue;
		std::cerr << "phased with " << phaseBlocks.size() << " blocks, effective ploidy " << phaseBlocks[0].allelesPerHaplotype.size() << ". start: " << (phaseBlocks[0].chainStartPhased ? "yes" : "no") << ", end: " << (phaseBlocks.back().chainEndPhased ? "yes" : "no") << std::endl;
		bool include = true;
		if (include) result.insert(result.end(), phaseBlocks.begin(), phaseBlocks.end());
	}
	return result;
}

void markNodePositions(std::vector<std::pair<size_t, size_t>>& result, const uint64_t startNode, const uint64_t endNode, const SparseEdgeContainer& edges, const size_t chain, const size_t offset)
{
	std::vector<uint64_t> stack;
	stack.push_back(startNode);
	phmap::flat_hash_set<uint64_t> visited;
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (visited.count(top & maskUint64_t) == 1) continue;
		visited.insert(top & maskUint64_t);
		assert(result[top & maskUint64_t].first == std::numeric_limits<size_t>::max() || (top & maskUint64_t) == (startNode & maskUint64_t));
		result[top & maskUint64_t] = std::make_pair(chain, offset);
		if (top != (startNode ^ firstBitUint64_t) && top != endNode)
		{
			for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
			{
				stack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
			}
		}
		if (top != startNode && top != (endNode ^ firstBitUint64_t))
		{
			for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, (top & firstBitUint64_t) ^ firstBitUint64_t)))
			{
				stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
			}
		}
	}
}

bool canUnambiguouslyMark(const SparseEdgeContainer& edges, const uint64_t startNode, const uint64_t endNode, const std::vector<bool>& anchor)
{
	std::vector<uint64_t> stack;
	phmap::flat_hash_set<uint64_t> visited;
	size_t anchorCount = 0;
	stack.push_back(startNode);
	while (stack.size() >= 1)
	{
		uint64_t top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		visited.insert(top);
		if (anchor[top & maskUint64_t])
		{
			anchorCount += 1;
			if (anchorCount > 2) return false;
		}
		if (!anchor[top & maskUint64_t]) stack.push_back(top ^ firstBitUint64_t);
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	if (visited.count(endNode) == 0 && visited.count(endNode ^ firstBitUint64_t) == 0) return false;
	if (anchorCount != 2) return false;
	return true;
}

std::vector<std::pair<size_t, size_t>> getNodeLocationsWithinChain(const UnitigGraph& unitigGraph, const std::vector<AnchorChain>& anchorChains)
{
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<std::pair<size_t, size_t>> result;
	result.resize(unitigGraph.nodeCount(), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			anchor[anchorChains[i].nodes[j] & maskUint64_t] = true;
		}
	}
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 1; j < anchorChains[i].nodes.size(); j++)
		{
			if (canUnambiguouslyMark(edges, anchorChains[i].nodes[j-1], anchorChains[i].nodes[j], anchor))
			{
				markNodePositions(result, anchorChains[i].nodes[j-1], anchorChains[i].nodes[j], edges, i, j-1);
			}
			result[anchorChains[i].nodes[j-1] & maskUint64_t] = std::make_pair(i, j-1);
			result[anchorChains[i].nodes[j] & maskUint64_t] = std::make_pair(i, j);
		}
		result[anchorChains[i].nodes.back() & maskUint64_t] = std::make_pair(i, anchorChains[i].nodes.size()-1);
	}
	return result;
}

// pathindex, nodeindex, readPos (extrapolated node end), haplotype, diagonal, node
std::vector<std::tuple<size_t, size_t, size_t, size_t, int, uint64_t>> getReadToAnchorMatches(const UnitigGraph& resultGraph, const std::vector<ReadPathBundle>& resultPaths, const std::vector<std::pair<size_t, size_t>>& nodeLocationInChain, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>>& haplotypeDiagonalsPerRead, const std::vector<PhaseBlock>& chainHaplotypes, const std::vector<int>& maxDiagonalDifferences, const size_t i)
{
	std::vector<std::tuple<size_t, size_t, size_t, size_t, int, uint64_t>> readToAnchorMatches;
	for (size_t j = 0; j < resultPaths[i].paths.size(); j++)
	{
		size_t readPos = resultPaths[i].paths[j].readStartPos;
		for (size_t k = 0; k < resultPaths[i].paths[j].path.size(); k++)
		{
			uint64_t node = resultPaths[i].paths[j].path[k];
			readPos += resultGraph.lengths[node & maskUint64_t];
			if (k == 0) readPos -= resultPaths[i].paths[j].pathLeftClipKmers;
			assert((node & maskUint64_t) < nodeLocationInChain.size());
			size_t chain = nodeLocationInChain[node & maskUint64_t].first;
			size_t offset = nodeLocationInChain[node & maskUint64_t].second;
			if (chain == std::numeric_limits<size_t>::max()) continue;
			if ((node & maskUint64_t) != (anchorChains[chain].nodes[offset] & maskUint64_t))
			{
				continue;
			}
			bool fw = (node & firstBitUint64_t) == (anchorChains[chain].nodes[offset] & firstBitUint64_t);
			int diagonal = (int)readPos - (int)anchorChains[chain].nodeOffsets[offset] - (int)resultGraph.lengths[anchorChains[chain].nodes[offset] & maskUint64_t];
			if (!fw) diagonal = (int)readPos + (int)anchorChains[chain].nodeOffsets[offset];
			size_t uniqueHap = std::numeric_limits<size_t>::max();
			for (auto t : haplotypeDiagonalsPerRead[i])
			{
				if (chainHaplotypes[std::get<0>(t)].chainNumber != chain) continue;
				if (std::get<1>(t) != fw) continue;
				if (std::get<2>(t) > diagonal + (int)maxDiagonalDifferences[chain]) continue;
				if (std::get<2>(t) < diagonal - (int)maxDiagonalDifferences[chain]) continue;
				if (uniqueHap == std::numeric_limits<size_t>::max())
				{
					uniqueHap = std::get<3>(t);
				}
				else if (uniqueHap != std::get<3>(t))
				{
					uniqueHap = std::numeric_limits<size_t>::max()-1;
				}
			}
			if (anchorChains[chain].ploidy == 1)
			{
				assert(uniqueHap == std::numeric_limits<size_t>::max());
				uniqueHap = 0;
			}
			readToAnchorMatches.emplace_back(j, k, readPos, uniqueHap, diagonal, node);
		}
	}
	return readToAnchorMatches;
}

void unzipDiploidPhaseBlocks(UnitigGraph& resultGraph, std::vector<ReadPathBundle>& resultPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<PhaseBlock>& chainHaplotypes, const std::vector<std::pair<size_t, size_t>>& nodeLocationInChain, const phmap::flat_hash_set<std::pair<size_t, size_t>>& solvedAnchors, const phmap::flat_hash_set<std::pair<size_t, size_t>>& solvedBetweenAnchorLocations, const phmap::flat_hash_set<size_t>& solvedInterchainTangles, const VectorWithDirection<size_t>& nodeLocationInInterchainTangles, const size_t tangleCount, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>>& haplotypeDiagonalsPerRead, const phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>>& validHaplotypeConnectionsPerChainEdge, const phmap::flat_hash_set<std::pair<size_t, size_t>>& phasedSites)
{
	std::vector<int> maxDiagonalDifferences;
	maxDiagonalDifferences.resize(anchorChains.size());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		maxDiagonalDifferences[i] = ((int)(anchorChains[i].nodeOffsets.back() + resultGraph.lengths[anchorChains[i].nodes.back() & maskUint64_t]) * .5);
	}
	phmap::flat_hash_map<std::tuple<uint64_t, uint64_t, size_t, size_t>, std::vector<size_t>> tangleGapLengths; // chain1, chain2, hap1, hap2
	phmap::flat_hash_map<std::tuple<size_t, size_t, size_t>, std::vector<size_t>> anchorGapLengths; // chain, node, hap
	for (size_t i = 0; i < resultPaths.size(); i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t, int, uint64_t>> readToAnchorMatches = getReadToAnchorMatches(resultGraph, resultPaths, nodeLocationInChain, anchorChains, haplotypeDiagonalsPerRead, chainHaplotypes, maxDiagonalDifferences, i); // pathindex, nodeindex, readPos (extrapolated node end), haplotype, diagonal, node
		for (size_t m = 1; m < readToAnchorMatches.size(); m++)
		{
			uint64_t prevNode = std::get<5>(readToAnchorMatches[m-1]);
			size_t prevChain = nodeLocationInChain[prevNode & maskUint64_t].first;
			size_t prevOffset = nodeLocationInChain[prevNode & maskUint64_t].second;
			uint64_t currNode = std::get<5>(readToAnchorMatches[m]);
			size_t currChain = nodeLocationInChain[currNode & maskUint64_t].first;
			size_t currOffset = nodeLocationInChain[currNode & maskUint64_t].second;
			if (prevChain == currChain)
			{
				assert((anchorChains[currChain].nodes[currOffset] & maskUint64_t) == (currNode & maskUint64_t));
				assert((anchorChains[prevChain].nodes[prevOffset] & maskUint64_t) == (prevNode & maskUint64_t));
				size_t prevhap = std::get<3>(readToAnchorMatches[m-1]);
				size_t currhap = std::get<3>(readToAnchorMatches[m]);
				bool fw = currNode == anchorChains[currChain].nodes[currOffset];
				if (fw && currOffset != prevOffset+1) continue;
				if (!fw && prevOffset != currOffset+1) continue;
				if (prevhap != currhap) continue;
				size_t prevEndPos = std::get<2>(readToAnchorMatches[m-1]);
				size_t currStartPos = std::get<2>(readToAnchorMatches[m]) - resultGraph.lengths[currNode];
				if (currStartPos < prevEndPos) continue;
				size_t distance = currStartPos - prevEndPos;
				anchorGapLengths[std::make_tuple(currChain, std::min(currOffset, prevOffset), currhap)].emplace_back(distance);
			}
			else
			{
				bool prevFw = prevNode == anchorChains[prevChain].nodes[prevOffset];
				bool currFw = currNode == anchorChains[currChain].nodes[currOffset];
				if (prevFw && prevOffset != anchorChains[prevChain].nodes.size()-1) continue;
				if (!prevFw && prevOffset != 0) continue;
				if (currFw && currOffset != anchorChains[currChain].nodes.size()-1) continue;
				if (!currFw && currOffset != 0) continue;
				uint64_t orientedPrevChain = prevChain + (prevFw ? firstBitUint64_t : 0);
				uint64_t orientedCurrChain = currChain + (currFw ? firstBitUint64_t : 0);
				if (validHaplotypeConnectionsPerChainEdge.count(orientedPrevChain) == 0) continue;
				if (validHaplotypeConnectionsPerChainEdge.at(orientedPrevChain).count(orientedCurrChain) == 0) continue;
				size_t prevhap = std::get<3>(readToAnchorMatches[m-1]);
				size_t currhap = std::get<3>(readToAnchorMatches[m]);
				if (validHaplotypeConnectionsPerChainEdge.at(orientedPrevChain).at(orientedCurrChain).count(std::make_pair(prevhap, currhap)) == 0) continue;
				if ((orientedCurrChain ^ firstBitUint64_t) < orientedPrevChain)
				{
					std::swap(orientedPrevChain, orientedCurrChain);
					std::swap(prevhap, currhap);
					orientedPrevChain ^= firstBitUint64_t;
					orientedCurrChain ^= firstBitUint64_t;
				}
				size_t prevEndPos = std::get<2>(readToAnchorMatches[m-1]);
				size_t currStartPos = std::get<2>(readToAnchorMatches[m]) - resultGraph.lengths[currNode];
				if (currStartPos < prevEndPos) continue;
				size_t distance = currStartPos - prevEndPos;
				tangleGapLengths[std::make_tuple(orientedPrevChain, orientedPrevChain, prevhap, currhap)].emplace_back(distance);
			}
		}
	}
	phmap::flat_hash_map<std::tuple<size_t, size_t, size_t, size_t>, size_t> nodeReplacement; // node, chain, offset, hap
	std::vector<std::pair<uint64_t, uint64_t>> newFakeEdges;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			if (solvedAnchors.count(std::make_pair(i, j)) == 0) continue;
			if (anchorChains[i].ploidy == 1) continue;
			assert(anchorChains[i].ploidy == 2);
			nodeReplacement[std::make_tuple(anchorChains[i].nodes[j] & maskUint64_t, i, j, 0)] = resultGraph.nodeCount();
			resultGraph.lengths.emplace_back(resultGraph.lengths[anchorChains[i].nodes[j] & maskUint64_t]);
			resultGraph.coverages.emplace_back(0);
			nodeReplacement[std::make_tuple(anchorChains[i].nodes[j] & maskUint64_t, i, j, 1)] = resultGraph.nodeCount();
			resultGraph.lengths.emplace_back(resultGraph.lengths[anchorChains[i].nodes[j] & maskUint64_t]);
			resultGraph.coverages.emplace_back(0);
		}
		for (size_t j = 1; j < anchorChains[i].nodes.size(); j++)
		{
			if (solvedBetweenAnchorLocations.count(std::make_pair(i, j-1)) == 0) continue;
			if (anchorChains[i].ploidy == 1) continue;
			if (phasedSites.count(std::make_pair(i, j-1)) == 1) continue;
			assert(anchorChains[i].ploidy == 2);
			size_t hap0length = 0;
			size_t hap1length = 0;
			if (anchorChains[i].nodeOffsets[j] >= (anchorChains[i].nodeOffsets[j-1] + resultGraph.lengths[anchorChains[i].nodes[j-1] & maskUint64_t]))
			{
				hap0length = anchorChains[i].nodeOffsets[j] >= (anchorChains[i].nodeOffsets[j-1] + resultGraph.lengths[anchorChains[i].nodes[j-1] & maskUint64_t]);
				hap1length = anchorChains[i].nodeOffsets[j] >= (anchorChains[i].nodeOffsets[j-1] + resultGraph.lengths[anchorChains[i].nodes[j-1] & maskUint64_t]);
			}
			if (anchorGapLengths.count(std::make_tuple(i, j-1, 0)) == 1)
			{
				std::vector<size_t>& vals = anchorGapLengths.at(std::make_tuple(i, j-1, 0));
				std::sort(vals.begin(), vals.end());
				hap0length = vals[vals.size()/2];
			}
			if (anchorGapLengths.count(std::make_tuple(i, j-1, 1)) == 1)
			{
				std::vector<size_t>& vals = anchorGapLengths.at(std::make_tuple(i, j-1, 1));
				std::sort(vals.begin(), vals.end());
				hap1length = vals[vals.size()/2];
			}
			size_t thisNode = resultGraph.lengths.size();
			resultGraph.lengths.emplace_back(hap0length);
			resultGraph.coverages.emplace_back(1);
			resultGraph.lengths.emplace_back(hap1length);
			resultGraph.coverages.emplace_back(1);
			size_t hap0Prev = anchorChains[i].nodes[j-1] & maskUint64_t;
			size_t hap1Prev = anchorChains[i].nodes[j-1] & maskUint64_t;
			size_t hap0Next = anchorChains[i].nodes[j] & maskUint64_t;
			size_t hap1Next = anchorChains[i].nodes[j] & maskUint64_t;
			assert(nodeReplacement.count(std::make_tuple(anchorChains[i].nodes[j-1] & maskUint64_t, i, j-1, 0)) == 1 || nodeReplacement.count(std::make_tuple(anchorChains[i].nodes[j] & maskUint64_t, i, j, 0)) == 1);
			if (nodeReplacement.count(std::make_tuple(anchorChains[i].nodes[j-1] & maskUint64_t, i, j-1, 0)) == 1)
			{
				hap0Prev = nodeReplacement.at(std::make_tuple(anchorChains[i].nodes[j-1] & maskUint64_t, i, j-1, 0));
				assert(nodeReplacement.count(std::make_tuple(anchorChains[i].nodes[j-1] & maskUint64_t, i, j-1, 1)) == 1);
				hap1Prev = nodeReplacement.at(std::make_tuple(anchorChains[i].nodes[j-1] & maskUint64_t, i, j-1, 1));
			}
			if (nodeReplacement.count(std::make_tuple(anchorChains[i].nodes[j] & maskUint64_t, i, j, 0)) == 1)
			{
				hap0Next = nodeReplacement.at(std::make_tuple(anchorChains[i].nodes[j] & maskUint64_t, i, j, 0));
				assert(nodeReplacement.count(std::make_tuple(anchorChains[i].nodes[j] & maskUint64_t, i, j, 1)) == 1);
				hap1Next = nodeReplacement.at(std::make_tuple(anchorChains[i].nodes[j] & maskUint64_t, i, j, 1));
			}
			newFakeEdges.emplace_back((anchorChains[i].nodes[j-1] & firstBitUint64_t) + hap0Prev, thisNode+0);
			newFakeEdges.emplace_back(thisNode+0, (anchorChains[i].nodes[j] & firstBitUint64_t) + hap0Next);
			newFakeEdges.emplace_back((anchorChains[i].nodes[j-1] & firstBitUint64_t) + hap1Prev, thisNode+1);
			newFakeEdges.emplace_back(thisNode+1, (anchorChains[i].nodes[j] & firstBitUint64_t) + hap1Next);
//			std::cerr << "insert chain hap gap of lengths " << hap0length << " " << hap1length << " between " << (anchorChains[i].nodes[j-1] & maskUint64_t) << " and " << (anchorChains[i].nodes[j] & maskUint64_t) << std::endl;
		}
	}
	for (const auto& pair : validHaplotypeConnectionsPerChainEdge)
	{
		uint64_t from = pair.first;
		uint64_t tip;
		if (from & firstBitUint64_t)
		{
			tip = anchorChains[from & maskUint64_t].nodes.back();
		}
		else
		{
			tip = anchorChains[from & maskUint64_t].nodes[0] ^ firstBitUint64_t;
		}
		size_t tangle = nodeLocationInInterchainTangles[std::make_pair(tip & maskUint64_t, tip & firstBitUint64_t)];
		if (solvedInterchainTangles.count(tangle) == 0) continue;
		for (const auto& pair2 : pair.second)
		{
			uint64_t to = pair2.first;
			uint64_t totip;
			if (to & firstBitUint64_t)
			{
				totip = anchorChains[to & maskUint64_t].nodes[0];
			}
			else
			{
				totip = anchorChains[to & maskUint64_t].nodes.back() ^ firstBitUint64_t;
			}
			if ((totip ^ firstBitUint64_t) < tip) continue;
			for (auto happair : pair2.second)
			{
				size_t distance = 1;
				if (tangleGapLengths.count(std::make_tuple(from, to, happair.first, happair.second)) == 1)
				{
					assert(tangleGapLengths.count(std::make_tuple(to ^ firstBitUint64_t, from ^ firstBitUint64_t, happair.second, happair.first)) == 0);
					std::vector<size_t>& distances = tangleGapLengths.at(std::make_tuple(from, to, happair.first, happair.second));
					std::sort(distances.begin(), distances.end());
					distance = distances[distances.size()/2];
				}
				else if (tangleGapLengths.count(std::make_tuple(to ^ firstBitUint64_t, from ^ firstBitUint64_t, happair.second, happair.first)) == 1)
				{
					std::vector<size_t>& distances = tangleGapLengths.at(std::make_tuple(from, to, happair.first, happair.second));
					std::sort(distances.begin(), distances.end());
					distance = distances[distances.size()/2];
				}
				if (distance == 0) distance = 1;
				uint64_t fromnode = tip;
				size_t fromchain = nodeLocationInChain[fromnode & maskUint64_t].first;
				size_t fromoffset = nodeLocationInChain[fromnode & maskUint64_t].second;
				if (nodeReplacement.count(std::make_tuple(fromnode & maskUint64_t, fromchain, fromoffset, happair.first)) == 1)
				{
					fromnode = (fromnode & firstBitUint64_t) + nodeReplacement.at(std::make_tuple(fromnode & maskUint64_t, fromchain, fromoffset, happair.first));
				}
				uint64_t tonode = totip;
				size_t tochain = nodeLocationInChain[tonode & maskUint64_t].first;
				size_t tooffset = nodeLocationInChain[tonode & maskUint64_t].second;
				if (nodeReplacement.count(std::make_tuple(tonode & maskUint64_t, tochain, tooffset, happair.second)) == 1)
				{
					tonode = (tonode & firstBitUint64_t) + nodeReplacement.at(std::make_tuple(tonode & maskUint64_t, tochain, tooffset, happair.second));
				}
				size_t newnode = resultGraph.lengths.size();
				resultGraph.lengths.emplace_back(distance);
				resultGraph.coverages.emplace_back(0);
				newFakeEdges.emplace_back(fromnode, newnode);
				newFakeEdges.emplace_back(newnode, tonode);
//				std::cerr << "insert tangle gap of length " << distance << " between " << (tip & maskUint64_t) << " (" << (fromnode & maskUint64_t) << ") and " << (totip & maskUint64_t) << " (" << (tonode & maskUint64_t) << ")" << std::endl;
			}
		}
	}
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> interchainTangleConnectionHap;
	interchainTangleConnectionHap.resize(tangleCount);
	for (size_t i = 0; i < resultPaths.size(); i++)
	{
		// std::cerr << "new read index " << i << std::endl;
		// std::cerr << "diagonals:";
		// for (auto t : haplotypeDiagonalsPerRead[i]) std::cerr << " " << std::get<0>(t) << "(" << chainHaplotypes[std::get<0>(t)].chainNumber << ")" << (std::get<1>(t) ? "+" : "-") << " " << std::get<2>(t) << " " << std::get<3>(t);
		// std::cerr << std::endl;
		std::vector<std::tuple<size_t, size_t, size_t, size_t, int, uint64_t>> readToAnchorMatches = getReadToAnchorMatches(resultGraph, resultPaths, nodeLocationInChain, anchorChains, haplotypeDiagonalsPerRead, chainHaplotypes, maxDiagonalDifferences, i); // pathindex, nodeindex, readPos (extrapolated node end), haplotype, diagonal, node
		std::vector<std::pair<size_t, size_t>> breakBeforeHere;
		for (size_t m = 0; m < readToAnchorMatches.size(); m++)
		{
			uint64_t currNode = std::get<5>(readToAnchorMatches[m]);
			size_t currChain = nodeLocationInChain[currNode & maskUint64_t].first;
			size_t currOffset = nodeLocationInChain[currNode & maskUint64_t].second;
			assert((anchorChains[currChain].nodes[currOffset] & maskUint64_t) == (currNode & maskUint64_t));
			size_t hap = std::get<3>(readToAnchorMatches[m]);
			size_t j = std::get<0>(readToAnchorMatches[m]);
			size_t k = std::get<1>(readToAnchorMatches[m]);
			if (anchorChains[currChain].ploidy == 2 && hap != std::numeric_limits<size_t>::max() && solvedAnchors.count(std::make_pair(currChain, currOffset)) == 1)
			{
				uint64_t node = resultPaths[i].paths[j].path[k];
				assert((node & maskUint64_t) < nodeLocationInChain.size());
				assert(nodeReplacement.count(std::make_tuple(node & maskUint64_t, currChain, currOffset, hap)) == 1);
//				std::cerr << "read " << i << " " << j << " " << k << " replace " << (node & maskUint64_t) << " with " << nodeReplacement.at(std::make_tuple(node & maskUint64_t, currChain, currOffset, hap)) << std::endl;
				resultPaths[i].paths[j].path[k] = nodeReplacement.at(std::make_tuple(node & maskUint64_t, currChain, currOffset, hap)) + (node & firstBitUint64_t);
			}
		}
		for (size_t m = 1; m < readToAnchorMatches.size(); m++)
		{
			uint64_t prevNode = std::get<5>(readToAnchorMatches[m-1]);
			size_t prevChain = nodeLocationInChain[prevNode & maskUint64_t].first;
			size_t prevOffset = nodeLocationInChain[prevNode & maskUint64_t].second;
			uint64_t currNode = std::get<5>(readToAnchorMatches[m]);
			size_t currChain = nodeLocationInChain[currNode & maskUint64_t].first;
			size_t currOffset = nodeLocationInChain[currNode & maskUint64_t].second;
			if (prevChain != currChain) continue;
			bool prevFw = prevNode == anchorChains[prevChain].nodes[prevOffset];
			bool currFw = currNode == anchorChains[currChain].nodes[currOffset];
			if (prevFw != currFw) continue;
			if (currFw && currOffset != prevOffset+1) continue;
			if (!currFw && currOffset+1 != prevOffset) continue;
			size_t prevhap = std::get<3>(readToAnchorMatches[m-1]);
			size_t currhap = std::get<3>(readToAnchorMatches[m]);
			if (prevhap != currhap) continue;
			if (phasedSites.count(std::make_pair(currChain, std::min(prevOffset, currOffset))) == 0) continue;
			size_t prevj = std::get<0>(readToAnchorMatches[m-1]);
			size_t prevk = std::get<1>(readToAnchorMatches[m-1]);
			size_t currj = std::get<0>(readToAnchorMatches[m]);
			size_t currk = std::get<1>(readToAnchorMatches[m]);
			for (size_t j = prevj; j <= currj; j++)
			{
				for (size_t k = (j == prevj ? prevk+1 : 0); k < (j == currj ? currk : resultPaths[i].paths[j].path.size()); k++)
				{
					uint64_t node = resultPaths[i].paths[j].path[k];
					if (nodeLocationInChain[node & maskUint64_t].first != currChain) continue;
					if (nodeLocationInChain[node & maskUint64_t].second != std::min(currOffset, prevOffset)) continue;
					assert(nodeLocationInChain[node & maskUint64_t].first == currChain);
					assert(nodeLocationInChain[node & maskUint64_t].second == std::min(currOffset, prevOffset));
					assert((anchorChains[nodeLocationInChain[node & maskUint64_t].first].nodes[nodeLocationInChain[node & maskUint64_t].second] & maskUint64_t) != (node & maskUint64_t));
					auto key = std::make_tuple(node & maskUint64_t, currChain, std::min(currOffset, prevOffset), currhap);
					if (nodeReplacement.count(key) == 0)
					{
						nodeReplacement[key] = resultGraph.lengths.size();
						resultGraph.lengths.emplace_back(resultGraph.lengths[node & maskUint64_t]);
						resultGraph.coverages.emplace_back(1);
					}
					resultPaths[i].paths[j].path[k] = (node & firstBitUint64_t) + nodeReplacement.at(key);
				}
			}
		}
		for (size_t m = 1; m < readToAnchorMatches.size(); m++)
		{
			uint64_t prevNode = std::get<5>(readToAnchorMatches[m-1]);
			size_t prevChain = nodeLocationInChain[prevNode & maskUint64_t].first;
			size_t prevOffset = nodeLocationInChain[prevNode & maskUint64_t].second;
			uint64_t currNode = std::get<5>(readToAnchorMatches[m]);
			size_t currChain = nodeLocationInChain[currNode & maskUint64_t].first;
			size_t currOffset = nodeLocationInChain[currNode & maskUint64_t].second;
			size_t prevj = std::get<0>(readToAnchorMatches[m-1]);
			size_t prevk = std::get<1>(readToAnchorMatches[m-1]);
			size_t currj = std::get<0>(readToAnchorMatches[m]);
			size_t currk = std::get<1>(readToAnchorMatches[m]);
			if (prevChain == currChain) continue;
			if (prevj != currj) continue;
			if (prevk+1 != currk) continue;
			size_t prevhap = std::get<3>(readToAnchorMatches[m-1]);
			size_t currhap = std::get<3>(readToAnchorMatches[m]);
			bool prevFw = prevNode == anchorChains[prevChain].nodes[prevOffset];
			bool currFw = currNode == anchorChains[currChain].nodes[currOffset];
			uint64_t orientedPrevChain = prevChain + (prevFw ? firstBitUint64_t : 0);
			uint64_t orientedCurrChain = currChain + (currFw ? firstBitUint64_t : 0);
			if (solvedInterchainTangles.count(nodeLocationInInterchainTangles[std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)]) == 1)
			{
				bool removeThis = true;
				if (validHaplotypeConnectionsPerChainEdge.count(orientedPrevChain) == 1)
				{
					if (validHaplotypeConnectionsPerChainEdge.at(orientedPrevChain).count(orientedCurrChain) == 1)
					{
						if (validHaplotypeConnectionsPerChainEdge.at(orientedPrevChain).at(orientedCurrChain).count(std::make_pair(prevhap, currhap)) == 1)
						{
							removeThis = false;
						}
					}
				}
				if (removeThis)
				{
					breakBeforeHere.emplace_back(currj, currk);
				}
			}
		}
		if (breakBeforeHere.size() > 0)
		{
			for (size_t j = breakBeforeHere.size()-1; j < breakBeforeHere.size(); j--)
			{
//				std::cerr << "break read " << i << " index " << breakBeforeHere[j].first << " " << breakBeforeHere[j].second << " (" << (resultPaths[i].paths[breakBeforeHere[j].first].path[breakBeforeHere[j].second-1] & maskUint64_t) << " " << (resultPaths[i].paths[breakBeforeHere[j].first].path[breakBeforeHere[j].second] & maskUint64_t) << ")" << std::endl;
				assert(breakBeforeHere[j].first < resultPaths[i].paths.size());
				assert(breakBeforeHere[j].second < resultPaths[i].paths[breakBeforeHere[j].first].path.size());
				assert(breakBeforeHere[j].second > 0);
				resultPaths[i].paths.emplace_back();
				resultPaths[i].paths.back().pathLeftClipKmers = 0;
				resultPaths[i].paths.back().pathRightClipKmers = resultPaths[i].paths[breakBeforeHere[j].first].pathRightClipKmers;
				size_t startPos = resultPaths[i].paths[breakBeforeHere[j].first].readStartPos;
				for (size_t k = 0; k < breakBeforeHere[j].second; k++)
				{
					startPos += resultGraph.lengths[resultPaths[i].paths[breakBeforeHere[j].first].path[k] & maskUint64_t];
					if (k == 0) startPos -= resultPaths[i].paths[breakBeforeHere[j].first].pathLeftClipKmers;
				}
				resultPaths[i].paths.back().readStartPos = startPos;
				resultPaths[i].paths.back().path.insert(resultPaths[i].paths.back().path.end(), resultPaths[i].paths[breakBeforeHere[j].first].path.begin() + breakBeforeHere[j].second, resultPaths[i].paths[breakBeforeHere[j].first].path.end());
				resultPaths[i].paths[breakBeforeHere[j].first].path.erase(resultPaths[i].paths[breakBeforeHere[j].first].path.begin() + breakBeforeHere[j].second, resultPaths[i].paths[breakBeforeHere[j].first].path.end());
				resultPaths[i].paths[breakBeforeHere[j].first].pathRightClipKmers = 0;
			}
			std::sort(resultPaths[i].paths.begin(), resultPaths[i].paths.end(), [](const auto& left, const auto& right) { return left.readStartPos < right.readStartPos; });
		}
	}
	assert(resultGraph.edgeCoverages.size() <= resultGraph.lengths.size());
	if (resultGraph.edgeCoverages.size() < resultGraph.lengths.size())
	{
		resultGraph.edgeCoverages.resize(resultGraph.lengths.size());
	}
	for (auto t : newFakeEdges)
	{
		std::pair<size_t, bool> from { t.first & maskUint64_t, t.first & firstBitUint64_t };
		std::pair<size_t, bool> to { t.second & maskUint64_t, t.second & firstBitUint64_t };
		resultGraph.edgeCoverages.set(from, to, 1);
	}
}

void unzipPhaseBlocks(UnitigGraph& resultGraph, std::vector<ReadPathBundle>& resultPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<PhaseBlock>& chainHaplotypes, const std::vector<std::pair<size_t, size_t>>& nodeLocationInChain, const phmap::flat_hash_set<std::pair<size_t, size_t>>& solvedAnchors, const phmap::flat_hash_set<std::pair<size_t, size_t>>& solvedBetweenAnchorLocations, const phmap::flat_hash_set<size_t>& solvedInterchainTangles, const VectorWithDirection<size_t>& nodeLocationInInterchainTangles, const size_t tangleCount, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>>& haplotypeDiagonalsPerRead, const phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>>& validHaplotypeConnectionsPerChainEdge)
{
	std::vector<int> maxDiagonalDifferences;
	maxDiagonalDifferences.resize(anchorChains.size());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		maxDiagonalDifferences[i] = ((int)(anchorChains[i].nodeOffsets.back() + resultGraph.lengths[anchorChains[i].nodes.back() & maskUint64_t]) * .5);
	}
	phmap::flat_hash_map<std::tuple<size_t, size_t, size_t, size_t>, size_t> nodeReplacement; // node, chain, offset, hap
	std::vector<phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>> interchainTangleConnectionHap;
	interchainTangleConnectionHap.resize(tangleCount);
	for (size_t i = 0; i < resultPaths.size(); i++)
	{
		// std::cerr << "new read index " << i << std::endl;
		// std::cerr << "diagonals:";
		// for (auto t : haplotypeDiagonalsPerRead[i]) std::cerr << " " << std::get<0>(t) << "(" << chainHaplotypes[std::get<0>(t)].chainNumber << ")" << (std::get<1>(t) ? "+" : "-") << " " << std::get<2>(t) << " " << std::get<3>(t);
		// std::cerr << std::endl;
		std::vector<std::tuple<size_t, size_t, size_t, size_t, int, uint64_t>> readToAnchorMatches; // chain, bubble, readPos (extrapolated node end), haplotype, diagonal, node
		for (size_t j = 0; j < resultPaths[i].paths.size(); j++)
		{
			size_t readPos = resultPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < resultPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = resultPaths[i].paths[j].path[k];
				readPos += resultGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= resultPaths[i].paths[j].pathLeftClipKmers;
				assert((node & maskUint64_t) < nodeLocationInChain.size());
				size_t chain = nodeLocationInChain[node & maskUint64_t].first;
				size_t offset = nodeLocationInChain[node & maskUint64_t].second;
				if (chain == std::numeric_limits<size_t>::max()) continue;
				if ((node & maskUint64_t) != (anchorChains[chain].nodes[offset] & maskUint64_t))
				{
					continue;
				}
				bool fw = (node & firstBitUint64_t) == (anchorChains[chain].nodes[offset] & firstBitUint64_t);
				int diagonal = (int)readPos - (int)anchorChains[chain].nodeOffsets[offset] - (int)resultGraph.lengths[anchorChains[chain].nodes[offset] & maskUint64_t];
				if (!fw) diagonal = (int)readPos + (int)anchorChains[chain].nodeOffsets[offset];
				size_t uniqueHap = std::numeric_limits<size_t>::max();
				for (auto t : haplotypeDiagonalsPerRead[i])
				{
					if (chainHaplotypes[std::get<0>(t)].chainNumber != chain) continue;
					if (std::get<1>(t) != fw) continue;
					if (std::get<2>(t) > diagonal + (int)maxDiagonalDifferences[chain]) continue;
					if (std::get<2>(t) < diagonal - (int)maxDiagonalDifferences[chain]) continue;
					if (uniqueHap == std::numeric_limits<size_t>::max())
					{
						uniqueHap = std::get<3>(t);
					}
					else if (uniqueHap != std::get<3>(t))
					{
						uniqueHap = std::numeric_limits<size_t>::max()-1;
					}
				}
				if (anchorChains[chain].ploidy == 1)
				{
					assert(uniqueHap == std::numeric_limits<size_t>::max());
					uniqueHap = 0;
				}
				readToAnchorMatches.emplace_back(j, k, readPos, uniqueHap, diagonal, node);
			}
		}
		std::vector<std::pair<size_t, size_t>> breakBeforeHere;
		for (size_t m = 1; m < readToAnchorMatches.size(); m++)
		{
			uint64_t lastNode = std::get<5>(readToAnchorMatches[m-1]);
			size_t lastChain = nodeLocationInChain[lastNode & maskUint64_t].first;
			size_t lastOffset = nodeLocationInChain[lastNode & maskUint64_t].second;
			assert((anchorChains[lastChain].nodes[lastOffset] & maskUint64_t) == (lastNode & maskUint64_t));
			bool lastFw = (lastNode & firstBitUint64_t) == (anchorChains[lastChain].nodes[lastOffset] & firstBitUint64_t);
			uint64_t currNode = std::get<5>(readToAnchorMatches[m]);
			size_t currChain = nodeLocationInChain[currNode & maskUint64_t].first;
			size_t currOffset = nodeLocationInChain[currNode & maskUint64_t].second;
			assert((anchorChains[currChain].nodes[currOffset] & maskUint64_t) == (currNode & maskUint64_t));
			bool currFw = (currNode & firstBitUint64_t) == (anchorChains[currChain].nodes[currOffset] & firstBitUint64_t);
			bool solveThisBlock = false;
			size_t solutionHap = std::numeric_limits<size_t>::max();
			size_t solutionChain = std::numeric_limits<size_t>::max();
			size_t solutionOffset = std::numeric_limits<size_t>::max();
			// std::cerr << "check " << m << " " << lastChain << " " << lastOffset << " " << (lastFw ? "fw" : "bw") << " " << (lastNode & maskUint64_t) << " " << std::get<3>(readToAnchorMatches[m-1]) << " " << currChain << " " << currOffset << " " << (currFw ? "fw" : "bw") << " " << (currNode & maskUint64_t) << " " << std::get<3>(readToAnchorMatches[m]) << std::endl;
			if (currChain == lastChain && currFw == lastFw && ((currFw && currOffset == lastOffset+1) || (!currFw && currOffset+1 == lastOffset) || (currOffset == lastOffset)))
			{
				// std::cerr << "has chain match" << std::endl;
				int lastDiagonal = std::get<4>(readToAnchorMatches[m-1]);
				int currDiagonal = std::get<4>(readToAnchorMatches[m]);
				// std::cerr << "diagonals " << lastDiagonal << " " << currDiagonal << std::endl;
				if (lastDiagonal > currDiagonal - maxDiagonalDifferences[currChain] && lastDiagonal < currDiagonal + maxDiagonalDifferences[currChain])
				{
					// std::cerr << "has curr last diagonal" << std::endl;
					if (solvedBetweenAnchorLocations.count(std::make_pair(currChain, std::min(lastOffset, currOffset))) == 1)
					{
						// std::cerr << "has unique hap " << std::get<3>(readToAnchorMatches[m]) << " check solved location " << currChain << " " << std::min(lastOffset, currOffset) << std::endl;
						if (std::get<3>(readToAnchorMatches[m-1]) == std::get<3>(readToAnchorMatches[m]) && std::get<3>(readToAnchorMatches[m]) != std::numeric_limits<size_t>::max()-1 && std::get<3>(readToAnchorMatches[m]) != std::numeric_limits<size_t>::max())
						{
							// std::cerr << "has solved location" << std::endl;
							solveThisBlock = true;
							solutionHap = std::get<3>(readToAnchorMatches[m]);
							solutionChain = currChain;
							solutionOffset = std::min(currOffset, lastOffset);
						}
					}
					else
					{
						solveThisBlock = true;
						solutionHap = 0;
						solutionChain = currChain;
						solutionOffset = std::min(currOffset, lastOffset);
					}
				}
			}
			if (!solveThisBlock)
			{
				if (nodeLocationInInterchainTangles[std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t)] == nodeLocationInInterchainTangles[std::make_pair(currNode & maskUint64_t, (currNode ^ firstBitUint64_t) & firstBitUint64_t)])
				{
					if (nodeLocationInInterchainTangles[std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t)] != std::numeric_limits<size_t>::max())
					{
						size_t lastHap = std::get<3>(readToAnchorMatches[m-1]);
						size_t currHap = std::get<3>(readToAnchorMatches[m]);
						// std::cerr << "has interchain nodes" << std::endl;
						bool shouldBreakThis = false;
						if (validChainEdges.hasEdge(std::make_pair(lastChain, lastFw), std::make_pair(currChain, currFw)))
						{
							// std::cerr << "has interchain edge" << std::endl;
							if (solvedInterchainTangles.count(nodeLocationInInterchainTangles[std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t)]) == 1)
							{
								shouldBreakThis = true;
								// std::cerr << "is in solved tangle" << std::endl;
								if (validHaplotypeConnectionsPerChainEdge.count(lastChain + (lastFw ? firstBitUint64_t : 0)) == 1 && validHaplotypeConnectionsPerChainEdge.at(lastChain + (lastFw ? firstBitUint64_t : 0)).count(currChain + (currFw ? firstBitUint64_t : 0)) == 1)
								{
									// std::cerr << "has some valid haplotypes" << std::endl;
									if (validHaplotypeConnectionsPerChainEdge.at(lastChain + (lastFw ? firstBitUint64_t : 0)).at(currChain + (currFw ? firstBitUint64_t : 0)).count(std::make_pair(lastHap, currHap)) == 1)
									{
										// std::cerr << "has valid haplotype" << std::endl;
										shouldBreakThis = false;
										solveThisBlock = true;
										auto keypair = canon(std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t), std::make_pair(currNode & maskUint64_t, currNode & firstBitUint64_t));
										std::pair<uint64_t, uint64_t> key { keypair.first.first + (keypair.first.second ? firstBitUint64_t : 0), keypair.second.first + (keypair.second.second ? firstBitUint64_t : 0) };
										size_t tangle = nodeLocationInInterchainTangles[std::make_pair(lastNode & maskUint64_t, lastNode & firstBitUint64_t)];
										assert(tangle < interchainTangleConnectionHap.size());
										if (interchainTangleConnectionHap[tangle].count(key) == 0)
										{
											solutionHap = interchainTangleConnectionHap[tangle].size();
											interchainTangleConnectionHap[tangle][key] = solutionHap;
										}
										else
										{
											solutionHap = interchainTangleConnectionHap[tangle].at(key);
										}
									}
								}
							}
						}
						if (shouldBreakThis && std::get<0>(readToAnchorMatches[m-1]) == std::get<0>(readToAnchorMatches[m]) && std::get<1>(readToAnchorMatches[m-1])+1 == std::get<1>(readToAnchorMatches[m]))
						{
							// std::cerr << "should break read " << i << " from node " << (lastNode & maskUint64_t) << " to " << (currNode & maskUint64_t) << std::endl;
							breakBeforeHere.emplace_back(std::get<0>(readToAnchorMatches[m]), std::get<1>(readToAnchorMatches[m]));
						}
					}
				}
			}
			if (solveThisBlock)
			{
				// std::cerr << "read " << i << " solve from " << std::get<0>(readToAnchorMatches[m-1]) << "," << std::get<1>(readToAnchorMatches[m-1]) << " to " << std::get<0>(readToAnchorMatches[m]) << "," << std::get<1>(readToAnchorMatches[m]) << std::endl;
				// if (solveLastAnchor) std::cerr << "read " << i << " solve anchor " << std::get<0>(readToAnchorMatches[m-1]) << " " << std::get<1>(readToAnchorMatches[m-1]) << " " << (std::get<5>(readToAnchorMatches[m-1]) & maskUint64_t) << std::endl;
				// if (solveCurrAnchor) std::cerr << "read " << i << " solve anchor " << std::get<0>(readToAnchorMatches[m]) << " " << std::get<1>(readToAnchorMatches[m]) << " " << (std::get<5>(readToAnchorMatches[m]) & maskUint64_t) << std::endl;
				for (size_t j = std::get<0>(readToAnchorMatches[m-1]); j <= std::get<0>(readToAnchorMatches[m]); j++)
				{
					for (size_t k = (j == std::get<0>(readToAnchorMatches[m-1]) ? std::get<1>(readToAnchorMatches[m-1])+1 : 0); k < (j == std::get<0>(readToAnchorMatches[m]) ? std::get<1>(readToAnchorMatches[m]) : resultPaths[i].paths[j].path.size()); k++)
					{
						uint64_t node = resultPaths[i].paths[j].path[k];
						assert((node & maskUint64_t) < nodeLocationInChain.size());
						if (nodeReplacement.count(std::make_tuple(node & maskUint64_t, solutionChain, solutionOffset, solutionHap)) == 0)
						{
							nodeReplacement[std::make_tuple(node & maskUint64_t, solutionChain, solutionOffset, solutionHap)] = resultGraph.lengths.size();
							resultGraph.lengths.push_back(resultGraph.lengths[node & maskUint64_t]);
							resultGraph.coverages.emplace_back();
						}
						assert(nodeReplacement.count(std::make_tuple(node & maskUint64_t, solutionChain, solutionOffset, solutionHap)) == 1);
						// std::cerr << "read " << i << " " << j << " " << k << " replace " << (node & maskUint64_t) << " with " << nodeReplacement.at(std::make_tuple(node & maskUint64_t, solutionChain, solutionOffset, solutionHap)) << std::endl;
						resultPaths[i].paths[j].path[k] = nodeReplacement.at(std::make_tuple(node & maskUint64_t, solutionChain, solutionOffset, solutionHap)) + (node & firstBitUint64_t);
					}
				}
			}
		}
		for (size_t m = 0; m < readToAnchorMatches.size(); m++)
		{
			uint64_t currNode = std::get<5>(readToAnchorMatches[m]);
			size_t currChain = nodeLocationInChain[currNode & maskUint64_t].first;
			size_t currOffset = nodeLocationInChain[currNode & maskUint64_t].second;
			assert((anchorChains[currChain].nodes[currOffset] & maskUint64_t) == (currNode & maskUint64_t));
			size_t hap = std::get<3>(readToAnchorMatches[m]);
			size_t j = std::get<0>(readToAnchorMatches[m]);
			size_t k = std::get<1>(readToAnchorMatches[m]);
			if (hap != std::numeric_limits<size_t>::max() && solvedAnchors.count(std::make_pair(currChain, currOffset)) == 1)
			{
				uint64_t node = resultPaths[i].paths[j].path[k];
				assert((node & maskUint64_t) < nodeLocationInChain.size());
				if (nodeReplacement.count(std::make_tuple(node & maskUint64_t, currChain, currOffset, hap)) == 0)
				{
					nodeReplacement[std::make_tuple(node & maskUint64_t, currChain, currOffset, hap)] = resultGraph.lengths.size();
					resultGraph.lengths.push_back(resultGraph.lengths[node & maskUint64_t]);
					resultGraph.coverages.emplace_back();
				}
				assert(nodeReplacement.count(std::make_tuple(node & maskUint64_t, currChain, currOffset, hap)) == 1);
				// std::cerr << "read " << i << " " << j << " " << k << " replace " << (node & maskUint64_t) << " with " << nodeReplacement.at(std::make_tuple(node & maskUint64_t, currChain, currOffset, hap)) << std::endl;
				resultPaths[i].paths[j].path[k] = nodeReplacement.at(std::make_tuple(node & maskUint64_t, currChain, currOffset, hap)) + (node & firstBitUint64_t);
			}
		}
		if (breakBeforeHere.size() > 0)
		{
			for (size_t j = breakBeforeHere.size()-1; j < breakBeforeHere.size(); j--)
			{
				// std::cerr << "break read " << i << " index " << breakBeforeHere[j].first << " " << breakBeforeHere[j].second << " (" << (resultPaths[i].paths[breakBeforeHere[j].first].path[breakBeforeHere[j].second-1] & maskUint64_t) << " " << (resultPaths[i].paths[breakBeforeHere[j].first].path[breakBeforeHere[j].second] & maskUint64_t) << ")" << std::endl;
				assert(breakBeforeHere[j].first < resultPaths[i].paths.size());
				assert(breakBeforeHere[j].second < resultPaths[i].paths[breakBeforeHere[j].first].path.size());
				assert(breakBeforeHere[j].second > 0);
				resultPaths[i].paths.emplace_back();
				resultPaths[i].paths.back().pathLeftClipKmers = 0;
				resultPaths[i].paths.back().pathRightClipKmers = resultPaths[i].paths[breakBeforeHere[j].first].pathRightClipKmers;
				size_t startPos = resultPaths[i].paths[breakBeforeHere[j].first].readStartPos;
				for (size_t k = 0; k < breakBeforeHere[j].second; k++)
				{
					startPos += resultGraph.lengths[resultPaths[i].paths[breakBeforeHere[j].first].path[k] & maskUint64_t];
					if (k == 0) startPos -= resultPaths[i].paths[breakBeforeHere[j].first].pathLeftClipKmers;
				}
				resultPaths[i].paths.back().readStartPos = startPos;
				resultPaths[i].paths.back().path.insert(resultPaths[i].paths.back().path.end(), resultPaths[i].paths[breakBeforeHere[j].first].path.begin() + breakBeforeHere[j].second, resultPaths[i].paths[breakBeforeHere[j].first].path.end());
				resultPaths[i].paths[breakBeforeHere[j].first].path.erase(resultPaths[i].paths[breakBeforeHere[j].first].path.begin() + breakBeforeHere[j].second, resultPaths[i].paths[breakBeforeHere[j].first].path.end());
				resultPaths[i].paths[breakBeforeHere[j].first].pathRightClipKmers = 0;
			}
			std::sort(resultPaths[i].paths.begin(), resultPaths[i].paths.end(), [](const auto& left, const auto& right) { return left.readStartPos < right.readStartPos; });
		}
	}
}

phmap::flat_hash_set<uint64_t> getInterchainTangleNodes(const uint64_t start, const std::vector<bool>& anchor, const SparseEdgeContainer& edges)
{
	assert(anchor[start & maskUint64_t]);
	phmap::flat_hash_set<uint64_t> result;
	std::vector<uint64_t> stack;
	stack.push_back(start);
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (result.count(top) == 1) continue;
		result.insert(top);
		if (!anchor[top & maskUint64_t]) stack.push_back(top ^ firstBitUint64_t);
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	return result;
}

VectorWithDirection<size_t> getNodeLocationsWithinInterchainTangles(const UnitigGraph& unitigGraph, const std::vector<AnchorChain>& anchorChains)
{
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			anchor[anchorChains[i].nodes[j] & maskUint64_t] = true;
		}
	}
	VectorWithDirection<size_t> result;
	result.resize(unitigGraph.nodeCount(), std::numeric_limits<size_t>::max());
	phmap::flat_hash_set<uint64_t> checked;
	size_t tangleNumber = 0;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (checked.count(anchorChains[i].nodes[0] ^ firstBitUint64_t) == 0)
		{
			phmap::flat_hash_set<uint64_t> nodesHere = getInterchainTangleNodes(anchorChains[i].nodes[0] ^ firstBitUint64_t, anchor, edges);
			assert(nodesHere.size() >= 1);
			assert(nodesHere.count(anchorChains[i].nodes[0] ^ firstBitUint64_t) == 1);
			for (auto node : nodesHere)
			{
				checked.insert(node);
				assert(result[std::make_pair(node & maskUint64_t, node & firstBitUint64_t)] == std::numeric_limits<size_t>::max());
				result[std::make_pair(node & maskUint64_t, node & firstBitUint64_t)] = tangleNumber;
			}
			tangleNumber += 1;
		}
		if (checked.count(anchorChains[i].nodes.back()) == 0)
		{
			phmap::flat_hash_set<uint64_t> nodesHere = getInterchainTangleNodes(anchorChains[i].nodes.back(), anchor, edges);
			assert(nodesHere.size() >= 1);
			assert(nodesHere.count(anchorChains[i].nodes.back()) == 1);
			for (auto node : nodesHere)
			{
				checked.insert(node);
				assert(result[std::make_pair(node & maskUint64_t, node & firstBitUint64_t)] == std::numeric_limits<size_t>::max());
				result[std::make_pair(node & maskUint64_t, node & firstBitUint64_t)] = tangleNumber;
			}
			tangleNumber += 1;
		}
	}
	return result;
}

std::vector<double> getChainCoverages(const UnitigGraph& unitigGraph, const std::vector<AnchorChain>& anchorChains)
{
	std::vector<double> chainCoverages;
	chainCoverages.resize(anchorChains.size(), 0);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		size_t length = 0;
		double coverageSum = 0;
		for (auto node : anchorChains[i].nodes)
		{
			length += unitigGraph.lengths[node & maskUint64_t];
			coverageSum += unitigGraph.lengths[node & maskUint64_t] * unitigGraph.coverages[node & maskUint64_t];
		}
		chainCoverages[i] = coverageSum / (double)length;
	}
	return chainCoverages;
}

std::tuple<phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>>, std::vector<bool>, std::vector<bool>> getValidHaplotypeConnectionsPerChainEdge(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>>& haplotypeDiagonalsPerRead, const std::vector<PhaseBlock>& chainHaplotypes, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
	std::vector<int> maxDiagonalDifferences;
	maxDiagonalDifferences.resize(anchorChains.size());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		maxDiagonalDifferences[i] = ((int)(anchorChains[i].nodeOffsets.back() + unitigGraph.lengths[anchorChains[i].nodes.back() & maskUint64_t]) * .5);
	}
	assert(chainPositionsInReads.size() == readPaths.size());
	assert(readPaths.size() == haplotypeDiagonalsPerRead.size());
	std::vector<bool> startPhased;
	std::vector<bool> endPhased;
	startPhased.resize(anchorChains.size(), false);
	endPhased.resize(anchorChains.size(), false);
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		size_t chain = chainHaplotypes[i].chainNumber;
		if (chainHaplotypes[i].chainStartPhased) startPhased[chain] = true;
		if (chainHaplotypes[i].chainEndPhased) endPhased[chain] = true;
	}
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::cerr << "chain " << i << " phased start: " << (startPhased[i] ? "yes" : "no") << " end: " << (endPhased[i] ? "yes" : "no") << std::endl;
	}
	std::vector<std::vector<size_t>> chainStartHaplotypeInRead;
	std::vector<std::vector<size_t>> chainEndHaplotypeInRead;
	chainStartHaplotypeInRead.resize(chainPositionsInReads.size());
	chainEndHaplotypeInRead.resize(chainPositionsInReads.size());
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		chainStartHaplotypeInRead[i].resize(chainPositionsInReads[i].size(), std::numeric_limits<size_t>::max());
		chainEndHaplotypeInRead[i].resize(chainPositionsInReads[i].size(), std::numeric_limits<size_t>::max());
		for (size_t j = 0; j < chainPositionsInReads[i].size(); j++)
		{
			size_t uniqueStartHap = std::numeric_limits<size_t>::max();
			size_t uniqueEndHap = std::numeric_limits<size_t>::max();
			size_t chain = chainPositionsInReads[i][j].chain & maskUint64_t;
			bool fw = (chainPositionsInReads[i][j].chain & firstBitUint64_t);
			int startDiagonal = chainPositionsInReads[i][j].chainStartPosInRead;
			int endDiagonal = chainPositionsInReads[i][j].chainEndPosInRead - (int)anchorChains[chain].nodeOffsets.back() - (int)unitigGraph.lengths[anchorChains[chain].nodes.back() & maskUint64_t];
			if (!fw)
			{
				startDiagonal = chainPositionsInReads[i][j].chainStartPosInRead + (int)unitigGraph.lengths[anchorChains[chain].nodes.back()] + (int)anchorChains[chain].nodeOffsets.back();
				endDiagonal = chainPositionsInReads[i][j].chainEndPosInRead;
			}
			for (size_t k = 0; k < haplotypeDiagonalsPerRead[i].size(); k++)
			{
				if (chainHaplotypes[std::get<0>(haplotypeDiagonalsPerRead[i][k])].chainNumber != chain) continue;
				if (std::get<1>(haplotypeDiagonalsPerRead[i][k]) != fw) continue;
				if (std::get<2>(haplotypeDiagonalsPerRead[i][k]) > startDiagonal - maxDiagonalDifferences[chain] && std::get<2>(haplotypeDiagonalsPerRead[i][k]) < startDiagonal + maxDiagonalDifferences[chain])
				{
					if (uniqueStartHap == std::numeric_limits<size_t>::max())
					{
						uniqueStartHap = std::get<3>(haplotypeDiagonalsPerRead[i][k]);
					}
					else if (uniqueStartHap != std::get<3>(haplotypeDiagonalsPerRead[i][k]))
					{
						uniqueStartHap = std::numeric_limits<size_t>::max()-1;
					}
				}
				if (std::get<2>(haplotypeDiagonalsPerRead[i][k]) > endDiagonal - maxDiagonalDifferences[chain] && std::get<2>(haplotypeDiagonalsPerRead[i][k]) < endDiagonal + maxDiagonalDifferences[chain])
				{
					if (uniqueEndHap == std::numeric_limits<size_t>::max())
					{
						uniqueEndHap = std::get<3>(haplotypeDiagonalsPerRead[i][k]);
					}
					else if (uniqueEndHap != std::get<3>(haplotypeDiagonalsPerRead[i][k]))
					{
						uniqueEndHap = std::numeric_limits<size_t>::max()-1;
					}
				}
			}
			if (anchorChains[chain].ploidy == 1)
			{
				assert(uniqueStartHap == std::numeric_limits<size_t>::max());
				assert(uniqueEndHap == std::numeric_limits<size_t>::max());
				uniqueStartHap = 0;
				uniqueEndHap = 0;
			}
			chainStartHaplotypeInRead[i][j] = uniqueStartHap;
			chainEndHaplotypeInRead[i][j] = uniqueEndHap;
		}
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<std::pair<size_t, size_t>, size_t>>> edgeHapCoverage;
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		for (size_t j = 1; j < chainPositionsInReads[i].size(); j++)
		{
			uint64_t prevChain = chainPositionsInReads[i][j-1].chain;
			uint64_t currChain = chainPositionsInReads[i][j].chain;
			size_t prevHap = chainEndHaplotypeInRead[i][j-1];
			size_t currHap = chainStartHaplotypeInRead[i][j];
			if (prevHap == std::numeric_limits<size_t>::max()-1) continue;
			if (currHap == std::numeric_limits<size_t>::max()-1) continue;
			if (prevHap == std::numeric_limits<size_t>::max() && (((prevChain & firstBitUint64_t) && endPhased[prevChain & maskUint64_t]) || (((prevChain ^ firstBitUint64_t) & firstBitUint64_t) && startPhased[prevChain & maskUint64_t]))) continue;
			if (currHap == std::numeric_limits<size_t>::max() && (((currChain & firstBitUint64_t) && startPhased[currChain & maskUint64_t]) || (((currChain ^ firstBitUint64_t) & firstBitUint64_t) && endPhased[currChain & maskUint64_t]))) continue;
			if (!validChainEdges.hasEdge(std::make_pair(prevChain & maskUint64_t, prevChain & firstBitUint64_t), std::make_pair(currChain & maskUint64_t, currChain & firstBitUint64_t))) continue;
			edgeHapCoverage[prevChain][currChain][std::make_pair(prevHap, currHap)] += 1;
			edgeHapCoverage[currChain ^ firstBitUint64_t][prevChain ^ firstBitUint64_t][std::make_pair(currHap, prevHap)] += 1;
		}
	}
	auto chainCoverages = getChainCoverages(unitigGraph, anchorChains);
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> unfilteredResult;
	for (const auto& pair : edgeHapCoverage)
	{
		for (const auto& pair2 : pair.second)
		{
			for (const auto& pair3 : pair2.second)
			{
				std::cerr << "hap check chain " << (pair.first & maskUint64_t) << ((pair.first & firstBitUint64_t) ? "+" : "-") << " to chain " << (pair2.first & maskUint64_t) << ((pair2.first & firstBitUint64_t) ? "+" : "-") << " haps " << pair3.first.first << " " << pair3.first.second << " coverage " << pair3.second << std::endl;
				if (pair3.second < approxOneHapCoverage * 0.25) continue;
				if (pair3.second < approxOneHapCoverage * 0.5 && pair3.second < chainCoverages[pair.first & maskUint64_t] * 0.5 && pair3.second < chainCoverages[pair2.first & maskUint64_t] * 0.5) continue;
				std::cerr << "add hap-edge chain " << (pair.first & maskUint64_t) << ((pair.first & firstBitUint64_t) ? "+" : "-") << " to chain " << (pair2.first & maskUint64_t) << ((pair2.first & firstBitUint64_t) ? "+" : "-") << " haps " << pair3.first.first << " " << pair3.first.second << " coverage " << pair3.second << std::endl;
				unfilteredResult[pair.first][pair2.first].insert(pair3.first);
			}
		}
	}
	std::vector<bool> canUnzipStart;
	std::vector<bool> canUnzipEnd;
	canUnzipStart.resize(anchorChains.size(), true);
	canUnzipEnd.resize(anchorChains.size(), true);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 1 && anchorChains[i].ploidy != 2)
		{
			canUnzipEnd[i] = false;
			canUnzipStart[i] = false;
			continue;
		}
		if (anchorChains[i].ploidy == 2 && !startPhased[i]) canUnzipStart[i] = false;
		if (anchorChains[i].ploidy == 2 && !endPhased[i]) canUnzipEnd[i] = false;
		size_t totalFoundPloidyFw = 0;
		for (auto edge : validChainEdges.getEdges(std::make_pair(i, true)))
		{
			uint64_t from = (i + firstBitUint64_t);
			uint64_t to = edge.first + (edge.second ? firstBitUint64_t : 0);
			if (unfilteredResult.count(from) == 0)
			{
				canUnzipEnd[i] = false;
				break;
			}
			if (unfilteredResult.at(from).count(to) == 0)
			{
				canUnzipEnd[i] = false;
				break;
			}
			assert(edgeHapCoverage.count(from) == 1);
			assert(edgeHapCoverage.at(from).count(to) == 1);
			for (auto happair : edgeHapCoverage.at(from).at(to))
			{
				size_t ploidyHere = happair.second / approxOneHapCoverage + 0.5;
				if (ploidyHere == 2 && (double)happair.second / approxOneHapCoverage < 1.75 && (anchorChains[to & maskUint64_t].ploidy == 1 || anchorChains[from & maskUint64_t].ploidy == 1))
				{
					ploidyHere = 1;
				}
				if (ploidyHere == 0 && (double)happair.second / approxOneHapCoverage > 0.25)
				{
					ploidyHere = 1;
				}
				totalFoundPloidyFw += ploidyHere;
			}
		}
		if (totalFoundPloidyFw != anchorChains[i].ploidy)
		{
			canUnzipEnd[i] = false;
		}
		size_t totalFoundPloidyBw = 0;
		for (auto edge : validChainEdges.getEdges(std::make_pair(i, false)))
		{
			uint64_t from = (i);
			uint64_t to = edge.first + (edge.second ? firstBitUint64_t : 0);
			if (unfilteredResult.count(from) == 0)
			{
				canUnzipStart[i] = false;
				break;
			}
			if (unfilteredResult.at(from).count(to) == 0)
			{
				canUnzipStart[i] = false;
				break;
			}
			assert(edgeHapCoverage.count(from) == 1);
			assert(edgeHapCoverage.at(from).count(to) == 1);
			for (auto happair : edgeHapCoverage.at(from).at(to))
			{
				size_t ploidyHere = happair.second / approxOneHapCoverage + 0.5;
				if (ploidyHere == 2 && (double)happair.second / approxOneHapCoverage < 1.75 && (anchorChains[to & maskUint64_t].ploidy == 1 || anchorChains[from & maskUint64_t].ploidy == 1))
				{
					ploidyHere = 1;
				}
				if (ploidyHere == 0 && (double)happair.second / approxOneHapCoverage > 0.25)
				{
					ploidyHere = 1;
				}
				totalFoundPloidyBw += ploidyHere;
			}
		}
		if (totalFoundPloidyBw != anchorChains[i].ploidy)
		{
			canUnzipStart[i] = false;
		}
	}
	// could do smartly, but naive method is good enough
	while (true)
	{
		bool repeat = false;
		for (size_t i = 0; i < anchorChains.size(); i++)
		{
			if (!canUnzipStart[i])
			{
				for (auto edge : validChainEdges.getEdges(std::make_pair(i, false)))
				{
					if (edge.second)
					{
						if (canUnzipStart[edge.first]) repeat = true;
						canUnzipStart[edge.first] = false;
					}
					else
					{
						if (canUnzipEnd[edge.first]) repeat = true;
						canUnzipEnd[edge.first] = false;
					}
				}
			}
			if (!canUnzipEnd[i])
			{
				for (auto edge : validChainEdges.getEdges(std::make_pair(i, true)))
				{
					if (edge.second)
					{
						if (canUnzipStart[edge.first]) repeat = true;
						canUnzipStart[edge.first] = false;
					}
					else
					{
						if (canUnzipEnd[edge.first]) repeat = true;
						canUnzipEnd[edge.first] = false;
					}
				}
			}
		}
		if (!repeat) break;
	}
	for (size_t i = 0; i < canUnzipEnd.size(); i++)
	{
		std::cerr << "chain " << i << " can unzip start: " << (canUnzipStart[i] ? "yes" : "no") << " end: " << (canUnzipEnd[i] ? "yes" : "no") << std::endl;
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> result;
	for (const auto& pair : unfilteredResult)
	{
		if ((pair.first & firstBitUint64_t) && !canUnzipEnd[pair.first & maskUint64_t]) continue;
		if (((pair.first ^ firstBitUint64_t) & firstBitUint64_t) && !canUnzipStart[pair.first & maskUint64_t]) continue;
		for (const auto& pair2 : pair.second)
		{
			if ((pair2.first & firstBitUint64_t) && !canUnzipStart[pair2.first & maskUint64_t]) continue;
			if (((pair2.first ^ firstBitUint64_t) & firstBitUint64_t) && !canUnzipEnd[pair2.first & maskUint64_t]) continue;
			for (const auto& pair3 : pair2.second)
			{
				result[pair.first][pair2.first].insert(pair3);
			}
		}
	}
	return std::make_tuple(std::move(result), std::move(canUnzipStart), std::move(canUnzipEnd));
}

std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>> getPhaseBlockReadHaplotypeDiagonals(const std::vector<PhaseBlock>& chainHaplotypes, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ReadDiagonalAlleles>>& allelesPerReadPerChain, const size_t pathCount)
{
	std::vector<std::vector<std::vector<std::tuple<size_t, bool, int>>>> readsPerHaplotypePerChain;
	readsPerHaplotypePerChain.resize(chainHaplotypes.size());
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		const size_t chain = chainHaplotypes[i].chainNumber;
		const size_t ploidy = chainHaplotypes[i].allelesPerHaplotype.size();
		assert(ploidy <= anchorChains[chain].ploidy);
		assert(ploidy == anchorChains[chain].ploidy || anchorChains[chain].ploidy != 2);
		readsPerHaplotypePerChain[i].resize(ploidy);
		phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> alleleBelongsToUniqueHaplotype;
		phmap::flat_hash_map<std::pair<size_t, size_t>, std::vector<size_t>> alleleBelongsToNonuniqueHaplotypes;
		for (size_t j = 0; j < chainHaplotypes[i].bubbleIndices.size(); j++)
		{
			for (size_t k = 0; k < ploidy; k++)
			{
				assert(chainHaplotypes[i].allelesPerHaplotype[k].size() == chainHaplotypes[i].bubbleIndices.size());
				assert(chainHaplotypes[i].extraAllelesPerHaplotype[k].size() == chainHaplotypes[i].bubbleIndices.size());
				if (chainHaplotypes[i].allelesPerHaplotype[k][j] != std::numeric_limits<size_t>::max())
				{
					auto key = std::make_pair(chainHaplotypes[i].bubbleIndices[j], chainHaplotypes[i].allelesPerHaplotype[k][j]);
					if (alleleBelongsToUniqueHaplotype.count(key) == 1)
					{
						if (alleleBelongsToUniqueHaplotype[key] != std::numeric_limits<size_t>::max())
						{
							alleleBelongsToNonuniqueHaplotypes[key].emplace_back(alleleBelongsToUniqueHaplotype.at(key));
							alleleBelongsToUniqueHaplotype[key] = std::numeric_limits<size_t>::max();
						}
						alleleBelongsToNonuniqueHaplotypes[key].emplace_back(k);
					}
					else
					{
						alleleBelongsToUniqueHaplotype[key] = k;
					}
					// std::cerr << "block " << i << " chain " << chain << " bubble " << chainHaplotypes[i].bubbleIndices[j] << " allele " << chainHaplotypes[i].allelesPerHaplotype[k][j] << " belongs to hap " << k << std::endl;
				}
				for (size_t l = 0; l < chainHaplotypes[i].extraAllelesPerHaplotype[k][j].size(); l++)
				{
					auto key = std::make_pair(chainHaplotypes[i].bubbleIndices[j], chainHaplotypes[i].extraAllelesPerHaplotype[k][j][l]);
					if (alleleBelongsToUniqueHaplotype.count(key) == 1)
					{
						if (alleleBelongsToUniqueHaplotype[key] != std::numeric_limits<size_t>::max())
						{
							alleleBelongsToNonuniqueHaplotypes[key].emplace_back(alleleBelongsToUniqueHaplotype.at(key));
							alleleBelongsToUniqueHaplotype[key] = std::numeric_limits<size_t>::max();
						}
						alleleBelongsToNonuniqueHaplotypes[key].emplace_back(k);
					}
					else
					{
						alleleBelongsToUniqueHaplotype[key] = k;
					}
					// std::cerr << "block " << i << " chain " << chain << " bubble " << chainHaplotypes[i].bubbleIndices[j] << " allele " << chainHaplotypes[i].extraAllelesPerHaplotype[k][j][l] << " belongs to hap " << k << std::endl;
				}
			}
		}
		for (size_t j = 0; j < allelesPerReadPerChain[chain].size(); j++)
		{
			std::vector<size_t> countMatches;
			countMatches.resize(ploidy);
			size_t totalMatches = 0;
			for (auto pair : allelesPerReadPerChain[chain][j].alleles)
			{
				bool hasMatch = false;
				if (alleleBelongsToUniqueHaplotype.count(std::make_pair(pair.site, pair.allele)) == 1 && alleleBelongsToUniqueHaplotype.at(std::make_pair(pair.site, pair.allele)) != std::numeric_limits<size_t>::max())
				{
					countMatches[alleleBelongsToUniqueHaplotype.at(std::make_pair(pair.site, pair.allele))] += 1;
					hasMatch = true;
				}
				if (alleleBelongsToNonuniqueHaplotypes.count(std::make_pair(pair.site, pair.allele)) == 1)
				{
					for (size_t k : alleleBelongsToNonuniqueHaplotypes.at(std::make_pair(pair.site, pair.allele)))
					{
						countMatches[k] += 1;
					}
					hasMatch = true;
				}
				if (hasMatch) totalMatches += 1;
			}
			if (totalMatches == 0)
			{
				// std::cerr << "read " << allelesPerReadPerChain[chain][j].readId << " diagonal " << allelesPerReadPerChain[chain][j].diagonal << " has no het matches in block " << i << " chain " << chain << std::endl;
				continue;
			}
			size_t bestHap = 0;
			bool ambiguous = false;
			for (size_t k = 1; k < ploidy; k++)
			{
				if (countMatches[k] == countMatches[bestHap]) ambiguous = true;
				if (countMatches[k] > countMatches[bestHap])
				{
					bestHap = k;
					ambiguous = false;
				}
			}
			if (countMatches[bestHap] <= totalMatches * 0.9)
			{
				// std::cerr << "read " << allelesPerReadPerChain[chain][j].readId << " diagonal " << allelesPerReadPerChain[chain][j].diagonal << " is ambiguous (" << countMatches[bestHap] << " vs " << totalMatches << ") in block " << i << " chain " << chain << std::endl;
				continue;
			}
			if (ambiguous)
			{
				// std::cerr << "read " << allelesPerReadPerChain[chain][j].readId << " diagonal " << allelesPerReadPerChain[chain][j].diagonal << " is ambiguous (same number) in block " << i << " chain " << chain << std::endl;
				continue;
			}
			readsPerHaplotypePerChain[i][bestHap].emplace_back(allelesPerReadPerChain[chain][j].readId, allelesPerReadPerChain[chain][j].fw, allelesPerReadPerChain[chain][j].diagonal);
			// std::cerr << "read " << allelesPerReadPerChain[chain][j].readId << " diagonal " << allelesPerReadPerChain[chain][j].diagonal << " has hap " << bestHap << " in block " << i << " chain " << chain << std::endl;
		}
		for (size_t j = 0; j < ploidy; j++)
		{
			std::sort(readsPerHaplotypePerChain[i][j].begin(), readsPerHaplotypePerChain[i][j].end());
			assert(readsPerHaplotypePerChain[i][j].size() >= 1);
		}
	}
	std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>> haplotypeDiagonalsPerRead;
	haplotypeDiagonalsPerRead.resize(pathCount);
	for (size_t block = 0; block < readsPerHaplotypePerChain.size(); block++)
	{
		for (size_t hap = 0; hap < readsPerHaplotypePerChain[block].size(); hap++)
		{
			for (auto t : readsPerHaplotypePerChain[block][hap])
			{
				size_t read = std::get<0>(t);
				bool fw = std::get<1>(t);
				int diagonal = std::get<2>(t);
				haplotypeDiagonalsPerRead[read].emplace_back(block, fw, diagonal, hap);
				// std::cerr << "readhap " << read << " " << block << " (" << (chainHaplotypes[block].chainNumber & maskUint64_t) << ")" << (fw ? "+" : "-") << " " << diagonal << " " << hap << std::endl;
			}
		}
	}
	return haplotypeDiagonalsPerRead;
}

phmap::flat_hash_set<size_t> getSolvedInterchainTangles(const std::vector<AnchorChain>& anchorChains, const std::vector<bool>& canUnzipStart, const std::vector<bool>& canUnzipEnd, const VectorWithDirection<size_t>& nodeLocationInInterchainTangles)
{
	phmap::flat_hash_set<size_t> solvedInterchainTangles;
	phmap::flat_hash_set<size_t> unsolvedInterchainTangles;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		// canUnzipStart[i] = false;
		// canUnzipEnd[i] = false;
		if (canUnzipStart[i])
		{
			std::pair<size_t, bool> key { anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] ^ firstBitUint64_t) & firstBitUint64_t };
			assert(nodeLocationInInterchainTangles[key] != std::numeric_limits<size_t>::max());
			solvedInterchainTangles.insert(nodeLocationInInterchainTangles[key]);
		}
		else
		{
			std::pair<size_t, bool> key { anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] ^ firstBitUint64_t) & firstBitUint64_t };
			assert(nodeLocationInInterchainTangles[key] != std::numeric_limits<size_t>::max());
			unsolvedInterchainTangles.insert(nodeLocationInInterchainTangles[key]);
		}
		if (canUnzipEnd[i])
		{
			std::pair<size_t, bool> key { anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t };
			assert(nodeLocationInInterchainTangles[key] != std::numeric_limits<size_t>::max());
			solvedInterchainTangles.insert(nodeLocationInInterchainTangles[key]);
		}
		else
		{
			std::pair<size_t, bool> key { anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t };
			assert(nodeLocationInInterchainTangles[key] != std::numeric_limits<size_t>::max());
			unsolvedInterchainTangles.insert(nodeLocationInInterchainTangles[key]);
		}
	}
	for (auto tangle : unsolvedInterchainTangles)
	{
		if (solvedInterchainTangles.count(tangle) == 1) solvedInterchainTangles.erase(tangle);
	}
	// no tangles with too many chains
	phmap::flat_hash_map<size_t, size_t> chainsPerTangle;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::pair<size_t, bool> startNode { anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] ^ firstBitUint64_t) & firstBitUint64_t };
		std::pair<size_t, bool> endNode { anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t };
		if (nodeLocationInInterchainTangles[startNode] != std::numeric_limits<size_t>::max())
		{
			chainsPerTangle[nodeLocationInInterchainTangles[startNode]] += 1;
		}
		if (nodeLocationInInterchainTangles[endNode] != std::numeric_limits<size_t>::max() && nodeLocationInInterchainTangles[endNode] != nodeLocationInInterchainTangles[startNode])
		{
			chainsPerTangle[nodeLocationInInterchainTangles[endNode]] += 1;
		}
	}
	for (auto pair : chainsPerTangle)
	{
		if (pair.second <= 4) continue;
		if (solvedInterchainTangles.count(pair.first) == 1) solvedInterchainTangles.erase(pair.first);
	}
	// arbitrarily say that too big tangles are not solved
	phmap::flat_hash_map<size_t, size_t> tangleSize;
	for (size_t i = 0; i < nodeLocationInInterchainTangles.size(); i++)
	{
		if (nodeLocationInInterchainTangles[std::make_pair(i, true)] == std::numeric_limits<size_t>::max()) continue;
		if (nodeLocationInInterchainTangles[std::make_pair(i, false)] == std::numeric_limits<size_t>::max()) continue;
		if (nodeLocationInInterchainTangles[std::make_pair(i, false)] != nodeLocationInInterchainTangles[std::make_pair(i, true)]) continue;
		tangleSize[nodeLocationInInterchainTangles[std::make_pair(i, true)]] += 1;
	}
	for (auto pair : tangleSize)
	{
		if (pair.second < 100) continue;
		if (solvedInterchainTangles.count(pair.first) == 1) solvedInterchainTangles.erase(pair.first);
	}
	for (auto tangle : solvedInterchainTangles)
	{
		std::cerr << "solved interchain tangle " << tangle << std::endl;
	}
	return solvedInterchainTangles;
}

std::pair<phmap::flat_hash_set<std::pair<size_t, size_t>>, phmap::flat_hash_set<std::pair<size_t, size_t>>> getSolvedLocationsInChain(const std::vector<AnchorChain>& anchorChains, const std::vector<PhaseBlock>& chainHaplotypes, const std::vector<bool>& canUnzipStart, const std::vector<bool>& canUnzipEnd)
{
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedAnchors;
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedBetweenAnchorLocations;
	for (size_t chain = 0; chain < anchorChains.size(); chain++)
	{
		if (anchorChains[chain].ploidy != 1) continue;
		size_t minSolvedIndex = 0;
		size_t maxSolvedIndex = anchorChains[chain].nodes.size();
		// std::cerr << "haploid chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << " (" << anchorChains[chain].nodes.size() << ")" << std::endl;
//		if (!canUnzipStart[chain] && minSolvedIndex == 0) minSolvedIndex = 1;
//		if (!canUnzipEnd[chain] && maxSolvedIndex == anchorChains[chain].nodes.size()) maxSolvedIndex = anchorChains[chain].nodes.size();
		// std::cerr << "haploid chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << std::endl;
		if (maxSolvedIndex < minSolvedIndex) continue;
		for (size_t j = minSolvedIndex; j <= maxSolvedIndex; j++)
		{
			assert(solvedBetweenAnchorLocations.count(std::make_pair(chain, j)) == 0);
			// std::cerr << "solve " << chain << " " << j << " (" << (j > 0 ? (anchorChains[chain].nodes[j-1] & maskUint64_t) : -1) << " to " << (j < anchorChains[chain].nodes.size() ? (anchorChains[chain].nodes[j] & maskUint64_t) : std::numeric_limits<size_t>::max()) << ")" << std::endl;
			solvedBetweenAnchorLocations.emplace(chain, j);
			if (j != maxSolvedIndex) solvedAnchors.emplace(chain, j);
		}
		if (minSolvedIndex > 0) solvedBetweenAnchorLocations.emplace(chain, minSolvedIndex-1);
	}
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		const size_t chain = chainHaplotypes[i].chainNumber;
		size_t minSolvedIndex = chainHaplotypes[i].bubbleIndices[0];
		size_t maxSolvedIndex = chainHaplotypes[i].bubbleIndices.back();
		assert(chainHaplotypes[i].bubbleIndices.size() >= 2);
		assert(anchorChains[chainHaplotypes[i].chainNumber].ploidy >= 2);
		// std::cerr << "block " << i << " chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << " (" << anchorChains[chain].nodes.size() << ")" << std::endl;
//		if (!canUnzipStart[chain] && minSolvedIndex == 0) minSolvedIndex = chainHaplotypes[i].bubbleIndices[1];
//		if (!canUnzipEnd[chain] && maxSolvedIndex == anchorChains[chain].nodes.size()) maxSolvedIndex = chainHaplotypes[i].bubbleIndices[chainHaplotypes[i].bubbleIndices.size()-2];
		// std::cerr << "block " << i << " chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << std::endl;
		for (size_t j = minSolvedIndex; j < maxSolvedIndex; j++)
		{
			assert(solvedBetweenAnchorLocations.count(std::make_pair(chain, j)) == 0);
			// std::cerr << "solve " << chain << " " << j << " (" << (j > 0 ? (anchorChains[chain].nodes[j-1] & maskUint64_t) : -1) << " to " << (j < anchorChains[chain].nodes.size() ? (anchorChains[chain].nodes[j] & maskUint64_t) : std::numeric_limits<size_t>::max()) << ")" << std::endl;
			solvedBetweenAnchorLocations.emplace(chain, j);
			if ((j != 0 && j+1 != anchorChains[chain].nodes.size()) || (j == 0 && canUnzipStart[chain]) || (j+1 == anchorChains[chain].nodes.size() && canUnzipEnd[chain])) solvedAnchors.emplace(chain, j);
		}
		if (maxSolvedIndex == anchorChains[chain].nodes.size())
		{
			assert(solvedBetweenAnchorLocations.count(std::make_pair(chain, maxSolvedIndex)) == 0);
			solvedBetweenAnchorLocations.emplace(chain, maxSolvedIndex);
		}
		if (minSolvedIndex > 0) solvedBetweenAnchorLocations.emplace(chain, minSolvedIndex-1);
	}
	return std::make_pair(std::move(solvedAnchors), std::move(solvedBetweenAnchorLocations));
}

phmap::flat_hash_set<std::pair<size_t, size_t>> getPhasedSimpleBubbleLocations(const UnitigGraph& unitigGraph, const std::vector<AnchorChain>& anchorChains, const std::vector<PhaseBlock>& chainHaplotypes)
{
	auto edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	phmap::flat_hash_set<std::pair<size_t, size_t>> topologicalSimpleBubbles;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 2) continue;
		for (size_t j = 0; j+1 < anchorChains[i].nodes.size(); j++)
		{
			std::pair<size_t, bool> start { anchorChains[i].nodes[j] & maskUint64_t, anchorChains[i].nodes[j] & firstBitUint64_t };
			std::pair<size_t, bool> end { anchorChains[i].nodes[j+1] & maskUint64_t, anchorChains[i].nodes[j+1] & firstBitUint64_t };
			if (edges.getEdges(start).size() != 2) continue;
			bool bubble = true;
			for (auto edge : edges.getEdges(start))
			{
				if (edges.getEdges(edge).size() != 1)
				{
					bubble = false;
					break;
				}
				if (edges.getEdges(reverse(edge)).size() != 1)
				{
					bubble = false;
					break;
				}
				if (edges.getEdges(edge)[0] != end)
				{
					bubble = false;
					break;
				}
			}
			if (!bubble) continue;
			topologicalSimpleBubbles.emplace(i, j);
		}
	}
	phmap::flat_hash_set<std::pair<size_t, size_t>> result;
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		size_t chain = chainHaplotypes[i].chainNumber;
		for (size_t bubble : chainHaplotypes[i].bubbleIndices)
		{
			if (bubble == 0) continue;
			if (topologicalSimpleBubbles.count(std::make_pair(chain, bubble-1)) == 0) continue;
			result.emplace(chain, bubble-1);
		}
	}
	return result;
}

void markUnsolvableUnzippable(std::vector<bool>& canUnzipStart, std::vector<bool>& canUnzipEnd, const VectorWithDirection<size_t> nodeLocationInInterchainTangles, const std::vector<AnchorChain>& anchorChains, const phmap::flat_hash_set<size_t>& solvedInterchainTangles)
{
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		uint64_t startNode = anchorChains[i].nodes[0] ^ firstBitUint64_t;
		std::pair<size_t, bool> start { startNode & maskUint64_t, startNode & firstBitUint64_t };
		if (nodeLocationInInterchainTangles[start] != std::numeric_limits<size_t>::max())
		{
			if (solvedInterchainTangles.count(nodeLocationInInterchainTangles[start]) == 0)
			{
				if (canUnzipStart[i]) std::cerr << "mark chain " << i << " start unresolvable" << std::endl;
				canUnzipStart[i] = false;
			}
		}
		uint64_t endNode = anchorChains[i].nodes.back();
		std::pair<size_t, bool> end { endNode & maskUint64_t, endNode & firstBitUint64_t };
		if (nodeLocationInInterchainTangles[end] != std::numeric_limits<size_t>::max())
		{
			if (solvedInterchainTangles.count(nodeLocationInInterchainTangles[end]) == 0)
			{
				if (canUnzipEnd[i]) std::cerr << "mark chain " << i << " end unresolvable" << std::endl;
				canUnzipEnd[i] = false;
			}
		}
	}
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipDiploidPhaseBlocks(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ReadDiagonalAlleles>>& allelesPerReadPerChain, std::vector<PhaseBlock>& chainHaplotypes, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	std::vector<std::pair<size_t, size_t>> nodeLocationInChain = getNodeLocationsWithinChain(unitigGraph, anchorChains);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			assert(nodeLocationInChain[anchorChains[i].nodes[j] & maskUint64_t] == std::make_pair(i, j));
		}
	}
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			anchor[anchorChains[i].nodes[j] & maskUint64_t] = true;
		}
	}
	VectorWithDirection<size_t> nodeLocationInInterchainTangles = getNodeLocationsWithinInterchainTangles(unitigGraph, anchorChains);
	size_t tangleCount = 0;
	for (size_t i = 0; i < nodeLocationInInterchainTangles.size(); i++)
	{
		if (nodeLocationInInterchainTangles[std::make_pair(i, true)] != std::numeric_limits<size_t>::max())
		{
			tangleCount = std::max(tangleCount, nodeLocationInInterchainTangles[std::make_pair(i, true)]);
		}
		if (nodeLocationInInterchainTangles[std::make_pair(i, false)] != std::numeric_limits<size_t>::max())
		{
			tangleCount = std::max(tangleCount, nodeLocationInInterchainTangles[std::make_pair(i, false)]);
		}
	}
	tangleCount += 1;
	std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>> haplotypeDiagonalsPerRead = getPhaseBlockReadHaplotypeDiagonals(chainHaplotypes, anchorChains, allelesPerReadPerChain, resultPaths.size());
	std::vector<bool> canUnzipStart;
	std::vector<bool> canUnzipEnd;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> validHaplotypeConnectionsPerChainEdge;
	std::tie(validHaplotypeConnectionsPerChainEdge, canUnzipStart, canUnzipEnd) = getValidHaplotypeConnectionsPerChainEdge(unitigGraph, readPaths, anchorChains, haplotypeDiagonalsPerRead, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	for (size_t i = 0; i < nodeLocationInInterchainTangles.size(); i++)
	{
		std::cerr << "node location in interchain tangle node " << i << " fw " << nodeLocationInInterchainTangles[std::make_pair(i, true)] << " bw " << nodeLocationInInterchainTangles[std::make_pair(i, false)] << std::endl;
	}
	phmap::flat_hash_set<size_t> solvedInterchainTangles = getSolvedInterchainTangles(anchorChains, canUnzipStart, canUnzipEnd, nodeLocationInInterchainTangles);
	markUnsolvableUnzippable(canUnzipStart, canUnzipEnd, nodeLocationInInterchainTangles, anchorChains, solvedInterchainTangles);
	// for (size_t i = 0; i < anchorChains.size(); i++)
	// {
	// 	std::cerr << "anchor chain " << i << " (" << (anchorChains[i].nodes[0] & maskUint64_t) << " to " << (anchorChains[i].nodes.back() & maskUint64_t) << ") unzippable at start: " << (canUnzipStart[i] ? "yes" : "no") << " at end: " << (canUnzipEnd[i] ? "yes" : "no") << std::endl;
	// }
	bool removedAny = false;
	for (size_t i = chainHaplotypes.size()-1; i < chainHaplotypes.size(); i--)
	{
		size_t chain = chainHaplotypes[i].chainNumber;
		if (chainHaplotypes[i].bubbleIndices.back() == anchorChains[chain].nodes.size() && !canUnzipEnd[chain])
		{
			chainHaplotypes[i].bubbleIndices.pop_back();
			assert(chainHaplotypes[i].allelesPerHaplotype.size() == anchorChains[chain].ploidy);
			for (size_t j = 0; j < chainHaplotypes[i].allelesPerHaplotype.size(); j++)
			{
				chainHaplotypes[i].allelesPerHaplotype[j].pop_back();
				chainHaplotypes[i].extraAllelesPerHaplotype[j].pop_back();
			}
		}
		if (chainHaplotypes[i].bubbleIndices[0] == 0 && !canUnzipStart[chain])
		{
			chainHaplotypes[i].bubbleIndices.erase(chainHaplotypes[i].bubbleIndices.begin());
			assert(chainHaplotypes[i].allelesPerHaplotype.size() == anchorChains[chain].ploidy);
			for (size_t j = 0; j < chainHaplotypes[i].allelesPerHaplotype.size(); j++)
			{
				chainHaplotypes[i].allelesPerHaplotype[j].erase(chainHaplotypes[i].allelesPerHaplotype[j].begin());
				chainHaplotypes[i].extraAllelesPerHaplotype[j].erase(chainHaplotypes[i].extraAllelesPerHaplotype[j].begin());
			}
		}
		if (chainHaplotypes[i].bubbleIndices.size() < 2)
		{
			std::swap(chainHaplotypes[i], chainHaplotypes.back());
			chainHaplotypes.pop_back();
			removedAny = true;
		}
	}
	if (removedAny) return unzipDiploidPhaseBlocks(unitigGraph, readPaths, anchorChains, allelesPerReadPerChain, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedAnchors;
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedBetweenAnchorLocations;
	std::tie(solvedAnchors, solvedBetweenAnchorLocations) = getSolvedLocationsInChain(anchorChains, chainHaplotypes, canUnzipStart, canUnzipEnd);
	phmap::flat_hash_set<std::pair<size_t, size_t>> phasedSites;
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		size_t chain = chainHaplotypes[i].chainNumber;
		for (size_t bubble : chainHaplotypes[i].bubbleIndices)
		{
			if (bubble != 0) phasedSites.emplace(chain, bubble-1);
		}
	}
	unzipDiploidPhaseBlocks(resultGraph, resultPaths, anchorChains, chainHaplotypes, nodeLocationInChain, solvedAnchors, solvedBetweenAnchorLocations, solvedInterchainTangles, nodeLocationInInterchainTangles, tangleCount, validChainEdges, haplotypeDiagonalsPerRead, validHaplotypeConnectionsPerChainEdge, phasedSites);
	RankBitvector kept;
	kept.resize(resultGraph.nodeCount());
	for (size_t i = 0; i < kept.size(); i++)
	{
		kept.set(i, true);
		if (i >= nodeLocationInChain.size()) continue;
		if (nodeLocationInChain[i].first != std::numeric_limits<size_t>::max())
		{
			assert(nodeLocationInChain[i].first < anchorChains.size());
			if (anchorChains[nodeLocationInChain[i].first].ploidy == 2)
			{
				assert(nodeLocationInChain[i].second < anchorChains[nodeLocationInChain[i].first].nodes.size());
				if ((anchorChains[nodeLocationInChain[i].first].nodes[nodeLocationInChain[i].second] & maskUint64_t) == i)
				{
					if (solvedAnchors.count(nodeLocationInChain[i]) == 0) continue;
				}
				else
				{
					if (solvedBetweenAnchorLocations.count(nodeLocationInChain[i]) == 0) continue;
					//if (phasedSites.count(std::make_pair(nodeLocationInChain[i].first, nodeLocationInChain[i].second)) == 1) continue;
				}
				std::cerr << "don't keep " << i << ", inside phased chain" << std::endl;
				kept.set(i, false);
			}
		}
		if (nodeLocationInInterchainTangles[std::make_pair(i, true)] != std::numeric_limits<size_t>::max() && nodeLocationInInterchainTangles[std::make_pair(i, true)] == nodeLocationInInterchainTangles[std::make_pair(i, false)] && solvedInterchainTangles.count(nodeLocationInInterchainTangles[std::make_pair(i, true)]) == 1)
		{
			if (!anchor[i])
			{
				std::cerr << "don't keep " << i << ", inside phased tangle " << nodeLocationInInterchainTangles[std::make_pair(i, true)] << std::endl;
				kept.set(i, false);
			}
		}
	}
	return unitigifyWithFilter(resultGraph, resultPaths, kept);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipPhaseBlocks(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ReadDiagonalAlleles>>& allelesPerReadPerChain, const std::vector<PhaseBlock>& chainHaplotypes, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	std::vector<std::pair<size_t, size_t>> nodeLocationInChain = getNodeLocationsWithinChain(unitigGraph, anchorChains);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			assert(nodeLocationInChain[anchorChains[i].nodes[j] & maskUint64_t] == std::make_pair(i, j));
		}
	}
	VectorWithDirection<size_t> nodeLocationInInterchainTangles = getNodeLocationsWithinInterchainTangles(unitigGraph, anchorChains);
	size_t tangleCount = 0;
	for (size_t i = 0; i < nodeLocationInInterchainTangles.size(); i++)
	{
		if (nodeLocationInInterchainTangles[std::make_pair(i, true)] != std::numeric_limits<size_t>::max())
		{
			tangleCount = std::max(tangleCount, nodeLocationInInterchainTangles[std::make_pair(i, true)]);
		}
		if (nodeLocationInInterchainTangles[std::make_pair(i, false)] != std::numeric_limits<size_t>::max())
		{
			tangleCount = std::max(tangleCount, nodeLocationInInterchainTangles[std::make_pair(i, false)]);
		}
	}
	tangleCount += 1;
	std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>> haplotypeDiagonalsPerRead = getPhaseBlockReadHaplotypeDiagonals(chainHaplotypes, anchorChains, allelesPerReadPerChain, resultPaths.size());
	std::vector<bool> canUnzipStart;
	std::vector<bool> canUnzipEnd;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> validHaplotypeConnectionsPerChainEdge;
	std::tie(validHaplotypeConnectionsPerChainEdge, canUnzipStart, canUnzipEnd) = getValidHaplotypeConnectionsPerChainEdge(unitigGraph, readPaths, anchorChains, haplotypeDiagonalsPerRead, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	// for (size_t i = 0; i < anchorChains.size(); i++)
	// {
	// 	std::cerr << "anchor chain " << i << " (" << (anchorChains[i].nodes[0] & maskUint64_t) << " to " << (anchorChains[i].nodes.back() & maskUint64_t) << ") unzippable at start: " << (canUnzipStart[i] ? "yes" : "no") << " at end: " << (canUnzipEnd[i] ? "yes" : "no") << std::endl;
	// }
	phmap::flat_hash_set<size_t> solvedInterchainTangles = getSolvedInterchainTangles(anchorChains, canUnzipStart, canUnzipEnd, nodeLocationInInterchainTangles);
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedAnchors;
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedBetweenAnchorLocations;
	std::tie(solvedAnchors, solvedBetweenAnchorLocations) = getSolvedLocationsInChain(anchorChains, chainHaplotypes, canUnzipStart, canUnzipEnd);
	unzipPhaseBlocks(resultGraph, resultPaths, anchorChains, chainHaplotypes, nodeLocationInChain, solvedAnchors, solvedBetweenAnchorLocations, solvedInterchainTangles, nodeLocationInInterchainTangles, tangleCount, validChainEdges, haplotypeDiagonalsPerRead, validHaplotypeConnectionsPerChainEdge);
	RankBitvector kept;
	kept.resize(resultGraph.nodeCount());
	for (size_t i = 0; i < kept.size(); i++)
	{
		kept.set(i, true);
		if (i >= nodeLocationInChain.size()) continue;
		if (nodeLocationInChain[i].first != std::numeric_limits<size_t>::max())
		{
			assert(nodeLocationInChain[i].first < anchorChains.size());
			assert(nodeLocationInChain[i].second < anchorChains[nodeLocationInChain[i].first].nodes.size());
			if ((anchorChains[nodeLocationInChain[i].first].nodes[nodeLocationInChain[i].second] & maskUint64_t) == i)
			{
				if (solvedAnchors.count(nodeLocationInChain[i]) == 0) continue;
			}
			else
			{
				if (solvedBetweenAnchorLocations.count(nodeLocationInChain[i]) == 0) continue;
			}
			// std::cerr << "don't keep " << i << ", inside phased chain" << std::endl;
			kept.set(i, false);
		}
		if (nodeLocationInInterchainTangles[std::make_pair(i, true)] != std::numeric_limits<size_t>::max() && nodeLocationInInterchainTangles[std::make_pair(i, true)] == nodeLocationInInterchainTangles[std::make_pair(i, false)] && solvedInterchainTangles.count(nodeLocationInInterchainTangles[std::make_pair(i, true)]) == 1)
		{
			// std::cerr << "don't keep " << i << ", inside phased tangle " << nodeLocationInInterchainTangles[std::make_pair(i, true)] << std::endl;
			kept.set(i, false);
		}
	}
	return filterUnitigGraph(resultGraph, resultPaths, kept);
}

SparseEdgeContainer getValidChainEdges(const UnitigGraph& unitigGraph, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
	std::vector<double> chainCoverages = getChainCoverages(unitigGraph, anchorChains);
	chainCoverages.resize(anchorChains.size(), 0);
	phmap::flat_hash_map<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>, size_t> edgeCoverage;
	size_t maxChain = 0;
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		if (chainPositionsInReads[i].size() == 0) continue;
		maxChain = std::max(maxChain, chainPositionsInReads[i][0].chain & maskUint64_t);
		for (size_t j = 1; j < chainPositionsInReads[i].size(); j++)
		{
			maxChain = std::max(maxChain, chainPositionsInReads[i][j].chain & maskUint64_t);
			const std::pair<size_t, bool> prevChain = std::make_pair(chainPositionsInReads[i][j-1].chain & maskUint64_t, chainPositionsInReads[i][j-1].chain & firstBitUint64_t);
			const std::pair<size_t, bool> currChain = std::make_pair(chainPositionsInReads[i][j].chain & maskUint64_t, chainPositionsInReads[i][j].chain & firstBitUint64_t);
			auto key = canon(prevChain, currChain);
			edgeCoverage[key] += 1;
		}
	}
	SparseEdgeContainer result;
	result.resize(maxChain+1);
	for (auto pair : edgeCoverage)
	{
		std::cerr << "check valid edge chain " << pair.first.first.first << (pair.first.first.second ? "+" : "-") << " to chain " << pair.first.second.first << (pair.first.second.second ? "+" : "-") << std::endl;
		std::cerr << "coverage " << pair.second << " vs " << (approxOneHapCoverage * 0.5) << " vs " << (chainCoverages[pair.first.first.first] * 0.5) << " vs " << (chainCoverages[pair.first.second.first] * 0.5) << std::endl;
		if (pair.second < approxOneHapCoverage * 0.5 && pair.second < chainCoverages[pair.first.first.first] * 0.5 && pair.second < chainCoverages[pair.first.second.first] * 0.5) continue;
		std::cerr << "include" << std::endl;
		// std::cerr << "chain edge between chains " << (pair.first.first.first) << " and " << (pair.first.second.first) << std::endl;
		result.addEdge(pair.first.first, pair.first.second);
		result.addEdge(reverse(pair.first.second), reverse(pair.first.first));
	}
	return result;
}

void fixFakeSinglePloidyChains(std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, size_t>> chainEdgeCoverage;
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		uint64_t lastChain = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < chainPositionsInReads[i].size(); j++)
		{
			if (anchorChains[chainPositionsInReads[i][j].chain & maskUint64_t].ploidy != 1) continue;
			if (lastChain == std::numeric_limits<size_t>::max())
			{
				lastChain = chainPositionsInReads[i][j].chain;
				continue;
			}
			uint64_t thisChain = chainPositionsInReads[i][j].chain;
			chainEdgeCoverage[lastChain][thisChain] += 1;
			chainEdgeCoverage[thisChain ^ firstBitUint64_t][lastChain ^ firstBitUint64_t] += 1;
			lastChain = thisChain;
		}
	}
	phmap::flat_hash_map<uint64_t, uint64_t> uniqueEdge;
	for (const auto& pair : chainEdgeCoverage)
	{
		uint64_t bestEdge = 0;
		size_t bestEdgeCoverage = 0;
		for (auto pair2 : pair.second)
		{
			if (pair2.second <= bestEdgeCoverage) continue;
			bestEdgeCoverage = pair2.second;
			bestEdge = pair2.first;
		}
		if (bestEdgeCoverage < approxOneHapCoverage * 0.5) continue;
		bool valid = true;
		for (auto pair2 : pair.second)
		{
			if (pair2.second < approxOneHapCoverage * 0.25) continue;
			if (pair2.first == bestEdge) continue;
			valid = false;
			break;
		}
		if (valid) uniqueEdge[pair.first] = bestEdge;
	}
	phmap::flat_hash_set<size_t> actuallyDiploid;
	for (const auto& pair : chainEdgeCoverage)
	{
		if (pair.second.size() < 2) continue;
		if (uniqueEdge.count(pair.first) == 1) continue;
		size_t numSolids = 0;
		bool allSolidsMatched = true;
		for (auto pair2 : pair.second)
		{
			if (pair2.second < approxOneHapCoverage * 0.25) continue;
			if (pair2.second < approxOneHapCoverage * 0.5)
			{
				allSolidsMatched = false;
				continue;
			}
			numSolids += 1;
			if (uniqueEdge.count(pair2.first ^ firstBitUint64_t) == 0)
			{
				allSolidsMatched = false;
				break;
			}
			if (uniqueEdge.at(pair2.first ^ firstBitUint64_t) != (pair.first ^ firstBitUint64_t))
			{
				allSolidsMatched = false;
				break;
			}
		}
		if (numSolids != 2) continue;
		if (!allSolidsMatched) continue;
		actuallyDiploid.insert(pair.first & maskUint64_t);
	}
	for (auto chain : actuallyDiploid)
	{
		std::cerr << "fixed chain " << chain << " from plody " << anchorChains[chain].ploidy << " to " << 2 << std::endl;
		anchorChains[chain].ploidy = 2;
	}
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphPolyploidTransitiveClosure(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage)
{
	std::cerr << "try phase polyploid chains by transitive closure" << std::endl;
	std::cerr << "get anchor chains" << std::endl;
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	fixFakeSinglePloidyChains(anchorChains, chainPositionsInReads, approxOneHapCoverage);
	SparseEdgeContainer validChainEdges = getValidChainEdges(unitigGraph, anchorChains, chainPositionsInReads, approxOneHapCoverage);
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<ReadDiagonalAlleles>> allelesPerReadPerChain;
	std::tie(allelesPerChain, allelesPerReadPerChain) = getRawReadInfoPerChain(anchorChains, readPaths, unitigGraph, validChainEdges, chainPositionsInReads);
	std::vector<PhaseBlock> chainHaplotypes = phasePolyploidsTransitiveClosure(anchorChains, allelesPerReadPerChain, approxOneHapCoverage);
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipPhaseBlocks(unitigGraph, readPaths, anchorChains, allelesPerReadPerChain, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphDiploidMEC(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t numThreads, const size_t graphk, const double approxOneHapCoverage)
{
	std::cerr << "try phase diploid chains with MEC" << std::endl;
	std::cerr << "get anchor chains" << std::endl;
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	fixFakeSinglePloidyChains(anchorChains, chainPositionsInReads, approxOneHapCoverage);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::cerr << "chain " << i << " ploidy " << anchorChains[i].ploidy;
		for (auto node : anchorChains[i].nodes)
		{
			std::cerr << " " << (node & maskUint64_t);
		}
		std::cerr << std::endl;
	}
/*	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		for (auto chain : chainPositionsInReads[i])
		{
			std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
		}
	}*/
	SparseEdgeContainer validChainEdges = getValidChainEdges(unitigGraph, anchorChains, chainPositionsInReads, approxOneHapCoverage);
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<ReadDiagonalAlleles>> allelesPerReadPerChain;
	std::tie(allelesPerChain, allelesPerReadPerChain) = getRawReadInfoPerChain(anchorChains, readPaths, unitigGraph, validChainEdges, chainPositionsInReads);
	std::vector<PhaseBlock> chainHaplotypes = phaseDiploidChains(anchorChains, allelesPerReadPerChain, unitigGraph, numThreads, approxOneHapCoverage);
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipDiploidPhaseBlocks(unitigGraph, readPaths, anchorChains, allelesPerReadPerChain, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}

std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> getHetInfoPerRead(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ReadDiagonalAlleles>>& allelesPerReadPerChain, const size_t readCount, const double approxOneHapCoverage)
{
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> result;
	std::vector<std::vector<bool>> haplotypeInformative;
	haplotypeInformative.resize(anchorChains.size());
	std::vector<std::vector<size_t>> bubbleAlleleCoverage;
	bubbleAlleleCoverage.resize(anchorChains.size());
	std::vector<std::vector<size_t>> bubbleAnchorCoverage;
	bubbleAnchorCoverage.resize(anchorChains.size());
	for (size_t chain = 0; chain < anchorChains.size(); chain++)
	{
		if (anchorChains[chain].ploidy == 1) continue;
		std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> mergedAllelesPerRead;
		std::vector<bool> maybeHaplotypeInformative;
		std::vector<size_t> numAllelesPerSite;
		std::tie(mergedAllelesPerRead, maybeHaplotypeInformative, numAllelesPerSite) = getPhasingTableInformation(chain, allelesPerReadPerChain[chain], anchorChains[chain].nodes, anchorChains[chain].ploidy, anchorChains[chain].nodes.size()+1, approxOneHapCoverage);
		haplotypeInformative[chain] = maybeHaplotypeInformative;
		bubbleAlleleCoverage[chain].resize(maybeHaplotypeInformative.size());
		bubbleAnchorCoverage[chain].resize(maybeHaplotypeInformative.size());
		if (haplotypeInformative[chain].size() != 0)
		{
			assert(haplotypeInformative[chain].size() >= 2);
			haplotypeInformative[chain][0] = false;
			haplotypeInformative[chain].back() = false;
		}
	}
	for (size_t chain = 0; chain < anchorChains.size(); chain++)
	{
		if (anchorChains[chain].ploidy == 1) continue;
		if (haplotypeInformative[chain].size() == 0) continue;
		for (const auto& readInfo : allelesPerReadPerChain[chain])
		{
			for (const auto& t : readInfo.alleles)
			{
				if (!haplotypeInformative[chain][t.site]) continue;
				bubbleAlleleCoverage[chain][t.site] += 1;
			}
		}
	}
	std::vector<std::pair<size_t, size_t>> nodeLocationInChain;
	nodeLocationInChain.resize(unitigGraph.nodeCount(), std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	for (size_t chain = 0; chain < anchorChains.size(); chain++)
	{
		if (haplotypeInformative[chain].size() == 0) continue;
		for (size_t j = 0; j < anchorChains[chain].nodes.size(); j++)
		{
			nodeLocationInChain[anchorChains[chain].nodes[j] & maskUint64_t] = std::make_pair(chain, j);
		}
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::pair<size_t, size_t> lastAnchor { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				std::pair<size_t, size_t> location = nodeLocationInChain[node & maskUint64_t];
				if (location.first == std::numeric_limits<size_t>::max()) continue;
				bool fw = (node & firstBitUint64_t) == (anchorChains[location.first].nodes[location.second] & firstBitUint64_t);
				assert(location.second != std::numeric_limits<size_t>::max());
				if (location.first != lastAnchor.first)
				{
					lastAnchor = location;
					continue;
				}
				if (fw && location.second != lastAnchor.second+1)
				{
					lastAnchor = location;
					continue;
				}
				if (!fw && location.second+1 != lastAnchor.second)
				{
					lastAnchor = location;
					continue;
				}
				if (haplotypeInformative[location.first].size() == 0)
				{
					lastAnchor = location;
					continue;
				}
				if (!haplotypeInformative[location.first][std::min(location.second, lastAnchor.second)+1])
				{
					lastAnchor = location;
					continue;
				}
				bubbleAnchorCoverage[location.first][std::min(location.second, lastAnchor.second)+1] += 1;
				lastAnchor = location;
			}
		}
	}
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (haplotypeInformative[i].size() == 0) continue;
		for (size_t j = 0; j < haplotypeInformative[i].size(); j++)
		{
			if (!haplotypeInformative[i][j]) continue;
			assert(bubbleAnchorCoverage[i][j] >= bubbleAlleleCoverage[i][j]);
			if (bubbleAlleleCoverage[i][j] < bubbleAnchorCoverage[i][j] * 0.95)
			{
				haplotypeInformative[i][j] = false;
			}
		}
	}
	result.resize(readCount);
	for (size_t chain = 0; chain < anchorChains.size(); chain++)
	{
		if (anchorChains[chain].ploidy == 1) continue;
		if (haplotypeInformative[chain].size() == 0) continue;
		for (const auto& readInfo : allelesPerReadPerChain[chain])
		{
			for (const auto& t : readInfo.alleles)
			{
				if (!haplotypeInformative[chain][t.site]) continue;
				// std::cerr << "read " << readInfo.readId << " has het at chain " << chain << (readInfo.fw ? "+" : "-") << " site " << t.site << " (" << (anchorChains[chain].nodes[t.site-1] & maskUint64_t) << " to " << (anchorChains[chain].nodes[t.site] & maskUint64_t) << ") allele " << t.allele << std::endl;
				result[readInfo.readId].emplace_back(t.readStartPos, t.readEndPos, chain, readInfo.fw, t.site, t.allele);
			}
		}
	}
	return result;
}

Context getContext(const int contextWindowStartPos, const size_t nodeStartPos, const size_t nodeEndPos, const std::vector<std::tuple<size_t, size_t, uint64_t>>& processedHetInfos, const std::vector<ChainPosition>& chainPositionsInReads, const size_t contextLength, bool fw)
{
	Context result;
	std::vector<std::pair<int, uint64_t>> preHets;
	std::vector<std::pair<int, uint64_t>> overlapHets;
	std::vector<std::pair<int, uint64_t>> postHets;
	for (auto m : processedHetInfos)
	{
		if ((int)std::get<1>(m) < contextWindowStartPos) continue;
		if ((int)std::get<0>(m) >= contextWindowStartPos+(int)contextLength) continue;
		if ((int)std::get<0>(m) >= (int)nodeEndPos)
		{
			postHets.emplace_back((std::get<0>(m) + std::get<1>(m))*0.5, std::get<2>(m));
		}
		else if ((int)std::get<1>(m) <= (int)nodeStartPos)
		{
			preHets.emplace_back((std::get<0>(m) + std::get<1>(m))*0.5, std::get<2>(m));
		}
		else
		{
			overlapHets.emplace_back((std::get<0>(m) + std::get<1>(m))*0.5, std::get<2>(m));
		}
	}
	std::sort(preHets.begin(), preHets.end());
	std::sort(overlapHets.begin(), overlapHets.end());
	std::sort(postHets.begin(), postHets.end());
	for (auto pair : preHets) result.preHets.emplace_back(pair.second);
	for (auto pair : overlapHets) result.overlapHets.emplace_back(pair.second);
	for (auto pair : postHets) result.postHets.emplace_back(pair.second);
	std::vector<std::pair<int, uint64_t>> preChains;
	std::vector<std::pair<int, uint64_t>> overlapChains;
	std::vector<std::pair<int, uint64_t>> postChains;
	for (auto m : chainPositionsInReads)
	{
		if (m.chainEndPosInRead < contextWindowStartPos) continue;
		if (m.chainStartPosInRead >= contextWindowStartPos+(int)contextLength) continue;
		if (m.chainStartPosInRead >= (int)nodeEndPos)
		{
			postChains.emplace_back((m.chainEndPosInRead + m.chainStartPosInRead)/2, m.chain);
		}
		else if (m.chainEndPosInRead <= (int)nodeStartPos)
		{
			preChains.emplace_back((m.chainEndPosInRead + m.chainStartPosInRead)/2, m.chain);
		}
		else
		{
			overlapChains.emplace_back((m.chainEndPosInRead + m.chainStartPosInRead)/2, m.chain);
		}
	}
	std::sort(preChains.begin(), preChains.end());
	std::sort(overlapChains.begin(), overlapChains.end());
	std::sort(postChains.begin(), postChains.end());
	for (auto pair : preChains) result.preChains.emplace_back(pair.second);
	for (auto pair : overlapChains) result.overlapChains.emplace_back(pair.second);
	for (auto pair : postChains) result.postChains.emplace_back(pair.second);
	if (!fw)
	{
		return result.reverse();
	}
	return result;
}

Context getCheapContext(const int contextWindowStartPos, const std::vector<std::tuple<size_t, size_t, uint64_t>>& processedHetInfos, const std::vector<ChainPosition>& chainPositionsInReads, const size_t contextLength)
{
	Context result;
	std::vector<std::pair<int, uint64_t>> hets;
	for (auto m : processedHetInfos)
	{
		if ((int)std::get<1>(m) < contextWindowStartPos) continue;
		if ((int)std::get<0>(m) >= contextWindowStartPos+(int)contextLength) continue;
		hets.emplace_back((std::get<0>(m) + std::get<1>(m))*0.5, std::get<2>(m));
	}
	std::sort(hets.begin(), hets.end());
	for (auto pair : hets) result.overlapHets.emplace_back(pair.second);
	std::vector<std::pair<int, uint64_t>> chains;
	for (auto m : chainPositionsInReads)
	{
		if (m.chainEndPosInRead < contextWindowStartPos) continue;
		if (m.chainStartPosInRead >= contextWindowStartPos+(int)contextLength) continue;
		chains.emplace_back((m.chainEndPosInRead + m.chainStartPosInRead)/2, m.chain);
	}
	std::sort(chains.begin(), chains.end());
	for (auto pair : chains) result.overlapChains.emplace_back(pair.second);
	return result;
}

template <typename F>
void iterateContextNodes(const std::vector<bool>& checkTheseNodes, const UnitigGraph& unitigGraph, const ReadPathBundle& readPath, const std::vector<size_t>& contextCheckStartPoses, const std::vector<std::tuple<size_t, size_t, uint64_t>>& processedHetInfos, const std::vector<ChainPosition>& chainPositionsInReads, const size_t resolveLength, F callback)
{
	for (size_t j = 0; j < readPath.paths.size(); j++)
	{
		size_t readPos = readPath.paths[j].readStartPos;
		for (size_t k = 0; k < readPath.paths[j].path.size(); k++)
		{
			uint64_t node = readPath.paths[j].path[k];
			readPos += unitigGraph.lengths[node & maskUint64_t];
			if (k == 0) readPos -= readPath.paths[j].pathLeftClipKmers;
			if (!checkTheseNodes[node & maskUint64_t]) continue;
			size_t startPos = readPos;
			if (k == 0) startPos += readPath.paths[j].pathLeftClipKmers;
			size_t endPos = readPos;
			if (k == readPath.paths[j].path.size()-1) endPos -= readPath.paths[j].pathRightClipKmers;
			assert(startPos >= unitigGraph.lengths[node & maskUint64_t]);
			startPos -= unitigGraph.lengths[node & maskUint64_t];
			int minCheck = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t] - (int)resolveLength;
			minCheck = std::max(minCheck, 0);
			int maxCheck = (int)readPos;
			maxCheck = std::min(maxCheck, (int)readPath.readLength - (int)resolveLength);
			phmap::flat_hash_set<Context> contexts;
			phmap::flat_hash_set<int> checkThese;
			checkThese.insert(minCheck);
			checkThese.insert(maxCheck-1);
			for (auto pos : contextCheckStartPoses)
			{
				if ((int)pos < minCheck) continue;
				if ((int)pos >= maxCheck) continue;
				checkThese.insert(pos);
			}
			for (int pos : checkThese)
			{
				Context context = getContext(pos, startPos, endPos, processedHetInfos, chainPositionsInReads, resolveLength, node & firstBitUint64_t);
				contexts.emplace(std::move(context));
			}
			if (contexts.size() == 0)
			{
				contexts.emplace();
			}
			std::vector<Context> contextVec { contexts.begin(), contexts.end() };
			callback(j, k, contextVec);
		}
	}
}

template <typename F>
void iterateContextNodes(const std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t>>>& nodePositionsInReads, const size_t checkThisNode, const std::vector<ReadPathBundle>& readPaths, const UnitigGraph& unitigGraph, const std::vector<std::vector<std::pair<size_t, size_t>>>& cheapContextsFw, const std::vector<std::vector<std::pair<size_t, size_t>>>& cheapContextsBw, const size_t resolveLength, F callback)
{
	for (auto t : nodePositionsInReads[checkThisNode])
	{
		size_t i = std::get<0>(t);
		size_t j = std::get<1>(t);
		size_t k = std::get<2>(t);
		size_t readPos = std::get<3>(t);
		uint64_t node = readPaths[i].paths[j].path[k];
		assert((node & maskUint64_t) == checkThisNode);
		int minCheck = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t] - (int)resolveLength;
		minCheck = std::max(minCheck, 0);
		int maxCheck = (int)readPos;
		maxCheck = std::min(maxCheck, (int)readPaths[i].readLength - (int)resolveLength);
		phmap::flat_hash_set<size_t> contexts;
		if (node & firstBitUint64_t)
		{
			for (auto pair : cheapContextsFw[i])
			{
				if ((int)pair.first < minCheck) continue;
				if ((int)pair.first >= maxCheck) continue;
				contexts.insert(pair.second);
			}
		}
		else
		{
			// for (auto pair : cheapContextsFw[i])
			for (auto pair : cheapContextsBw[i])
			{
				if ((int)pair.first < minCheck) continue;
				if ((int)pair.first >= maxCheck) continue;
				contexts.insert(pair.second);
			}
		}
		assert(contexts.size() >= 0);
		std::vector<size_t> contextVec { contexts.begin(), contexts.end() };
		callback(i, j, k, contextVec);
	}
}

template <typename F>
void iterateContextNodes(const std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t>>>& nodePositionsInReads, const size_t checkThisNode, const std::vector<ReadPathBundle>& readPaths, const UnitigGraph& unitigGraph, const std::vector<std::vector<size_t>>& contextCheckStartPoses, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& processedHetInfos, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const size_t resolveLength, F callback)
{
	for (auto t : nodePositionsInReads[checkThisNode])
	{
		size_t i = std::get<0>(t);
		size_t j = std::get<1>(t);
		size_t k = std::get<2>(t);
		size_t readPos = std::get<3>(t);
		uint64_t node = readPaths[i].paths[j].path[k];
		assert((node & maskUint64_t) == checkThisNode);
		size_t startPos = readPos;
		if (k == 0) startPos += readPaths[i].paths[j].pathLeftClipKmers;
		size_t endPos = readPos;
		if (k == readPaths[i].paths[j].path.size()-1) endPos -= readPaths[i].paths[j].pathRightClipKmers;
		assert(startPos >= unitigGraph.lengths[node & maskUint64_t]);
		startPos -= unitigGraph.lengths[node & maskUint64_t];
		int minCheck = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t] - (int)resolveLength;
		minCheck = std::max(minCheck, 0);
		int maxCheck = (int)readPos;
		maxCheck = std::min(maxCheck, (int)readPaths[i].readLength - (int)resolveLength);
		phmap::flat_hash_set<Context> contexts;
		phmap::flat_hash_set<int> checkThese;
		checkThese.insert(minCheck);
		checkThese.insert(maxCheck-1);
		for (auto pos : contextCheckStartPoses[i])
		{
			if ((int)pos < minCheck) continue;
			if ((int)pos >= maxCheck) continue;
			checkThese.insert(pos);
		}
		for (int pos : checkThese)
		{
			Context context = getContext(pos, startPos, endPos, processedHetInfos[i], chainPositionsInReads[i], resolveLength, node & firstBitUint64_t);
			contexts.emplace(std::move(context));
		}
		if (contexts.size() == 0)
		{
			contexts.emplace();
		}
		std::vector<Context> contextVec { contexts.begin(), contexts.end() };
		callback(i, j, k, contextVec);
	}
}

Context find(phmap::flat_hash_map<Context, Context>& parent, const Context& key)
{
	if (parent.count(key) == 0)
	{
		parent[key] = key;
		return key;
	}
	while (parent.at(parent.at(key)) != parent.at(key))
	{
		parent[key] = parent[parent[key]];
	}
	return parent.at(key);
}

void merge(phmap::flat_hash_map<Context, Context>& parent, const Context& left, const Context& right)
{
	auto leftp = find(parent, left);
	auto rightp = find(parent, right);
	parent[rightp] = leftp;
}

void markNonCheckable(std::vector<bool>& checkTheseNodes, const uint64_t startNode, const uint64_t endNode, const SparseEdgeContainer& edges)
{
	phmap::flat_hash_set<uint64_t> visited;
	std::vector<uint64_t> stack;
	stack.push_back(startNode);
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		visited.insert(top);
		if (top != startNode && top != endNode) stack.push_back(top ^ firstBitUint64_t);
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
	for (auto node : visited)
	{
		checkTheseNodes[node & maskUint64_t] = false;
	}
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipHapmers(const bool cheap, const std::vector<AnchorChain>& anchorChains, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>>& hetInfoPerRead, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const size_t resolveLength, const double approxOneHapCoverage)
{
	const size_t minContextCoverage = 3;
	assert(readPaths.size() == hetInfoPerRead.size());
	assert(readPaths.size() == chainPositionsInReads.size());
	std::vector<bool> checkTheseNodes;
	checkTheseNodes.resize(unitigGraph.nodeCount(), true);
	SparseEdgeContainer graphEdges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 1) continue;
		for (size_t j = 0; j+1 < anchorChains[i].nodes.size(); j++)
		{
			markNonCheckable(checkTheseNodes, anchorChains[i].nodes[j], anchorChains[i].nodes[j+1] ^ firstBitUint64_t, graphEdges);
		}
		for (size_t j = 0; j < anchorChains[i].nodes.size(); j++)
		{
			checkTheseNodes[anchorChains[i].nodes[j] & maskUint64_t] = false;
		}
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.lengths[i] >= resolveLength) checkTheseNodes[i] = false;
	}
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t>>> nodePositionsInReads;
	nodePositionsInReads.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				nodePositionsInReads[node & maskUint64_t].emplace_back(i, j, k, readPos);
			}
		}
	}
	phmap::flat_hash_map<std::tuple<size_t, size_t, size_t>, uint64_t> bubbleAlleleToIndex;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> processedHetInfos;
	processedHetInfos.resize(readPaths.size());
	for (size_t i = 0; i < hetInfoPerRead.size(); i++)
	{
		for (const auto t : hetInfoPerRead[i])
		{
			std::tuple<size_t, size_t, size_t> key { std::get<2>(t), std::get<4>(t), std::get<5>(t) };
			if (bubbleAlleleToIndex.count(key) == 0)
			{
				size_t index = bubbleAlleleToIndex.size();
				bubbleAlleleToIndex[key] = index;
			}
			processedHetInfos[i].emplace_back(std::get<0>(t), std::get<1>(t), bubbleAlleleToIndex.at(key) + (std::get<3>(t) ? firstBitUint64_t : 0));
		}
	}
	std::vector<std::vector<std::pair<size_t, size_t>>> cheapContextsFw;
	std::vector<std::vector<std::pair<size_t, size_t>>> cheapContextsBw;
	std::vector<std::vector<size_t>> contextCheckStartPoses;
	if (cheap)
	{
		cheapContextsFw.resize(readPaths.size());
		cheapContextsBw.resize(readPaths.size());
		for (size_t i = 0; i < readPaths.size(); i++)
		{
			std::vector<size_t> contextCheckStartPosesThisRead;
			for (auto t : hetInfoPerRead[i])
			{
				assert(std::get<0>(t) < readPaths[i].readLength);
				if (std::get<0>(t)+1 >= resolveLength) contextCheckStartPosesThisRead.emplace_back(std::get<0>(t)+1-resolveLength);
				assert(std::get<1>(t) >= 0);
				if ((int)std::get<1>(t)+1 < (int)readPaths[i].readLength - (int)resolveLength) contextCheckStartPosesThisRead.emplace_back(std::get<1>(t)+1);
			}
			for (auto t : chainPositionsInReads[i])
			{
				assert(t.chainStartPosInRead < (int)readPaths[i].readLength);
				if (t.chainStartPosInRead+1 > (int)resolveLength) contextCheckStartPosesThisRead.emplace_back(t.chainStartPosInRead+1-(int)resolveLength);
				assert(t.chainEndPosInRead >= 0);
				if (t.chainEndPosInRead+1 < (int)readPaths[i].readLength - (int)resolveLength) contextCheckStartPosesThisRead.emplace_back(t.chainEndPosInRead+1);
			}
			for (size_t j = 0; j < readPaths[i].paths.size(); j++)
			{
				size_t readPos = readPaths[i].paths[j].readStartPos;
				for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
				{
					uint64_t node = readPaths[i].paths[j].path[k];
					readPos += unitigGraph.lengths[node & maskUint64_t];
					if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
					if ((int)readPos+1 < (int)readPaths[i].readLength - (int)resolveLength) contextCheckStartPosesThisRead.push_back(readPos+1);
					if (readPos > resolveLength + unitigGraph.lengths[node & maskUint64_t]) contextCheckStartPosesThisRead.push_back(readPos - resolveLength - unitigGraph.lengths[node & maskUint64_t]);
				}
			}
			contextCheckStartPosesThisRead.push_back(0);
			contextCheckStartPosesThisRead.push_back(readPaths[i].readLength - resolveLength);
			phmap::flat_hash_set<size_t> uniq { contextCheckStartPosesThisRead.begin(), contextCheckStartPosesThisRead.end() };
			contextCheckStartPosesThisRead.clear();
			contextCheckStartPosesThisRead.insert(contextCheckStartPosesThisRead.end(), uniq.begin(), uniq.end());
			std::sort(contextCheckStartPosesThisRead.begin(), contextCheckStartPosesThisRead.end());
			for (auto pos : contextCheckStartPosesThisRead)
			{
				Context context = getCheapContext(pos, processedHetInfos[i], chainPositionsInReads[i], resolveLength);
				cheapContextsFw[i].emplace_back(pos, context.getHash());
				cheapContextsBw[i].emplace_back(pos, context.reverse().getHash());
				assert(context.getHash() == context.reverse().reverse().getHash());
			}
			std::sort(cheapContextsFw[i].begin(), cheapContextsFw[i].end());
			std::sort(cheapContextsBw[i].begin(), cheapContextsBw[i].end());
		}
	}
	else
	{
		contextCheckStartPoses.resize(readPaths.size());
		for (size_t i = 0; i < readPaths.size(); i++)
		{
			for (auto t : hetInfoPerRead[i])
			{
				assert(std::get<0>(t) < readPaths[i].readLength);
				if (std::get<0>(t)+1 >= resolveLength) contextCheckStartPoses[i].emplace_back(std::get<0>(t)+1-resolveLength);
				assert(std::get<1>(t) >= 0);
				if ((int)std::get<1>(t)+1 < (int)readPaths[i].readLength - (int)resolveLength) contextCheckStartPoses[i].emplace_back(std::get<1>(t)+1);
			}
			for (auto t : chainPositionsInReads[i])
			{
				assert(t.chainStartPosInRead < (int)readPaths[i].readLength);
				if (t.chainStartPosInRead+1 > (int)resolveLength) contextCheckStartPoses[i].emplace_back(t.chainStartPosInRead+1-(int)resolveLength);
				assert(t.chainEndPosInRead >= 0);
				if (t.chainEndPosInRead+1 < (int)readPaths[i].readLength - (int)resolveLength) contextCheckStartPoses[i].emplace_back(t.chainEndPosInRead+1);
			}
			contextCheckStartPoses[i].push_back(0);
			contextCheckStartPoses[i].push_back(readPaths[i].readLength - resolveLength);
			phmap::flat_hash_set<size_t> uniq { contextCheckStartPoses[i].begin(), contextCheckStartPoses[i].end() };
			contextCheckStartPoses[i].clear();
			contextCheckStartPoses[i].insert(contextCheckStartPoses[i].end(), uniq.begin(), uniq.end());
			std::sort(contextCheckStartPoses[i].begin(), contextCheckStartPoses[i].end());
		}
	}
	size_t nextNodeNumber = unitigGraph.nodeCount();
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	for (size_t node = 0; node < unitigGraph.nodeCount(); node++)
	{
		std::cerr << "hap phase check node " << node << "/" << unitigGraph.nodeCount() << " " << (checkTheseNodes[node] ? "yes" : "no") << std::endl;
		if (!checkTheseNodes[node]) continue;
		if (cheap)
		{
			phmap::flat_hash_map<size_t, phmap::flat_hash_set<size_t>> contextsPerRead;
			std::vector<std::vector<size_t>> totalNodeContexts;
			iterateContextNodes(nodePositionsInReads, node, readPaths, unitigGraph, cheapContextsFw, cheapContextsBw, resolveLength, [&contextsPerRead, &totalNodeContexts](const size_t i, const size_t j, const size_t k, const std::vector<size_t>& nodecontexts)
			{
				for (const auto& context : nodecontexts)
				{
					contextsPerRead[i].insert(context);
				}
				totalNodeContexts.emplace_back(nodecontexts);
			});
			phmap::flat_hash_map<size_t, size_t> contextCoverage;
			for (const auto& pair : contextsPerRead)
			{
				for (auto context : pair.second)
				{
					contextCoverage[context] += 1;
				}
			}
			size_t nodeContextsWithOnlyMinor = 0;
			for (const auto& contextGroup : totalNodeContexts)
			{
				bool hasMajor = false;
				for (auto context : contextGroup)
				{
					if (contextCoverage[context] >= minContextCoverage)
					{
						hasMajor = true;
						break;
					}
				}
				if (!hasMajor) nodeContextsWithOnlyMinor += 1;
			}
			std::cerr << "node-contexts with minor only " << nodeContextsWithOnlyMinor << " vs " << totalNodeContexts.size() << std::endl;
			if (nodeContextsWithOnlyMinor > approxOneHapCoverage * 0.25 || nodeContextsWithOnlyMinor > totalNodeContexts.size() * 0.1)
			{
				std::cerr << "too many minor-only node-contexts" << std::endl;
				checkTheseNodes[node] = false;
				continue;
			}
			phmap::flat_hash_map<size_t, size_t> nodeParent;
			iterateContextNodes(nodePositionsInReads, node, readPaths, unitigGraph, cheapContextsFw, cheapContextsBw, resolveLength, [&nodeParent, &contextCoverage, &readPaths, minContextCoverage](const size_t i, const size_t j, const size_t k, const std::vector<size_t>& contexts)
			{
				size_t anyContext;
				for (size_t context : contexts)
				{
					if (contextCoverage[context] < minContextCoverage) continue;
					anyContext = context;
					if (nodeParent.count(context) == 0) nodeParent[context] = context;
				}
				for (size_t context : contexts)
				{
					if (contextCoverage[context] < minContextCoverage) continue;
					merge(nodeParent, context, anyContext);
				}
			});
			if (nodeParent.size() == 1)
			{
				std::cerr << "only one context" << std::endl;
				checkTheseNodes[node] = false;
				continue;
			}
			phmap::flat_hash_map<size_t, size_t> contextToNewNodeNumber;
			phmap::flat_hash_set<size_t> keys;
			for (const auto& pair : nodeParent)
			{
				keys.insert(pair.first);
			}
			phmap::flat_hash_set<size_t> clusters;
			for (const auto& key : keys)
			{
				clusters.insert(find(nodeParent, key));
			}
			size_t splitNum = 0;
			for (const auto& cluster : clusters)
			{
				contextToNewNodeNumber[cluster] = nextNodeNumber;
				resultGraph.lengths.emplace_back(unitigGraph.lengths[node]);
				resultGraph.coverages.emplace_back(0);
				nextNodeNumber += 1;
				splitNum += 1;
			}
			if (splitNum == 1 || splitNum == 0)
			{
				std::cerr << "only " << splitNum << " clusters" << std::endl;
				checkTheseNodes[node] = false;
				continue;
			}
			std::cerr << "split node " << node << " into " << splitNum << " nodes" << std::endl;
			iterateContextNodes(nodePositionsInReads, node, readPaths, unitigGraph, cheapContextsFw, cheapContextsBw, resolveLength, [&nodeParent, &readPaths, &contextCoverage, &contextToNewNodeNumber, &resultPaths, minContextCoverage](const size_t i, const size_t j, const size_t k, const std::vector<size_t>& contexts)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				for (const auto& context : contexts)
				{
					if (contextCoverage[context] < minContextCoverage) continue;
					auto cluster = find(nodeParent, context);
					if (contextToNewNodeNumber.count(cluster) == 1)
					{
						// std::cerr << "replace read " << i << " " << j << " " << k << " (" << (node & maskUint64_t) << ") with " << contextToNewNodeNumber[node & maskUint64_t].at(cluster) << std::endl;
						resultPaths[i].paths[j].path[k] = (node & firstBitUint64_t) + (contextToNewNodeNumber.at(cluster));
						break;
					}
				}
			});
		}
		else
		{
			phmap::flat_hash_map<Context, size_t> contextCoverage;
			std::vector<std::vector<Context>> totalNodeContexts;
			iterateContextNodes(nodePositionsInReads, node, readPaths, unitigGraph, contextCheckStartPoses, processedHetInfos, chainPositionsInReads, resolveLength, [&contextCoverage, &totalNodeContexts](const size_t i, const size_t j, const size_t k, const std::vector<Context>& nodecontexts)
			{
				for (const auto& context : nodecontexts)
				{
					contextCoverage[context] += 1;
				}
				totalNodeContexts.emplace_back(nodecontexts);
			});
			size_t nodeContextsWithOnlyMinor = 0;
			for (const auto& contextGroup : totalNodeContexts)
			{
				bool hasMajor = false;
				for (const auto& context : contextGroup)
				{
					if (contextCoverage[context] >= minContextCoverage)
					{
						hasMajor = true;
						break;
					}
				}
				if (!hasMajor) nodeContextsWithOnlyMinor += 1;
			}
			std::cerr << "node-contexts with minor only " << nodeContextsWithOnlyMinor << " vs " << totalNodeContexts.size() << std::endl;
			if (nodeContextsWithOnlyMinor > approxOneHapCoverage * 0.25 || nodeContextsWithOnlyMinor > totalNodeContexts.size() * 0.1)
			{
				std::cerr << "too many minor-only node-contexts" << std::endl;
				checkTheseNodes[node] = false;
				continue;
			}
			phmap::flat_hash_map<Context, Context> nodeParent;
			iterateContextNodes(nodePositionsInReads, node, readPaths,unitigGraph, contextCheckStartPoses, processedHetInfos, chainPositionsInReads, resolveLength, [&nodeParent, &contextCoverage, &readPaths, minContextCoverage](const size_t i, const size_t j, const size_t k, const std::vector<Context>& contexts)
			{
				Context anyContext;
				for (const auto& context : contexts)
				{
					if (contextCoverage[context] < minContextCoverage) continue;
					anyContext = context;
					find(nodeParent, context);
				}
				for (const auto& context : contexts)
				{
					if (contextCoverage[context] < minContextCoverage) continue;
					merge(nodeParent, context, anyContext);
				}
			});
			phmap::flat_hash_map<Context, size_t> contextToNewNodeNumber;
			phmap::flat_hash_set<Context> keys;
			for (const auto& pair : nodeParent)
			{
				keys.insert(pair.first);
			}
			phmap::flat_hash_set<Context> clusters;
			for (const auto& key : keys)
			{
				clusters.insert(find(nodeParent, key));
			}
			size_t splitNum = 0;
			for (const auto& cluster : clusters)
			{
				contextToNewNodeNumber[cluster] = nextNodeNumber;
				resultGraph.lengths.emplace_back(unitigGraph.lengths[node]);
				resultGraph.coverages.emplace_back(0);
				nextNodeNumber += 1;
				splitNum += 1;
			}
			if (splitNum == 1 || splitNum == 0)
			{
				std::cerr << "only " << splitNum << " clusters" << std::endl;
				checkTheseNodes[node] = false;
				continue;
			}
			std::cerr << "split node " << node << " into " << splitNum << " nodes" << std::endl;
			iterateContextNodes(nodePositionsInReads, node, readPaths, unitigGraph, contextCheckStartPoses, processedHetInfos, chainPositionsInReads, resolveLength, [&nodeParent, &readPaths, &contextCoverage, &contextToNewNodeNumber, &resultPaths, minContextCoverage](const size_t i, const size_t j, const size_t k, const std::vector<Context>& contexts)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				for (const auto& context : contexts)
				{
					if (contextCoverage[context] < minContextCoverage) continue;
					auto cluster = find(nodeParent, context);
					if (contextToNewNodeNumber.count(cluster) == 1)
					{
						// std::cerr << "replace read " << i << " " << j << " " << k << " (" << (node & maskUint64_t) << ") with " << contextToNewNodeNumber[node & maskUint64_t].at(cluster) << std::endl;
						resultPaths[i].paths[j].path[k] = (node & firstBitUint64_t) + (contextToNewNodeNumber.at(cluster));
						break;
					}
				}
			});
		}
	}
	// std::vector<phmap::flat_hash_map<Context, Context>> nodeParent;
	// phmap::flat_hash_map<Context, size_t> contextCoverage;
	// nodeParent.resize(unitigGraph.nodeCount());
	// phmap::flat_hash_map<std::tuple<size_t, size_t, size_t>, uint64_t> bubbleAlleleToIndex;
	// std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> processedHetInfos;
	// processedHetInfos.resize(readPaths.size());
	// for (size_t i = 0; i < hetInfoPerRead.size(); i++)
	// {
	// 	for (const auto t : hetInfoPerRead[i])
	// 	{
	// 		std::tuple<size_t, size_t, size_t> key { std::get<2>(t), std::get<4>(t), std::get<5>(t) };
	// 		if (bubbleAlleleToIndex.count(key) == 0)
	// 		{
	// 			size_t index = bubbleAlleleToIndex.size();
	// 			bubbleAlleleToIndex[key] = index;
	// 		}
	// 		processedHetInfos[i].emplace_back(std::get<0>(t), std::get<1>(t), bubbleAlleleToIndex.at(key) + (std::get<3>(t) ? firstBitUint64_t : 0));
	// 	}
	// }
	// std::vector<std::vector<size_t>> contextCheckStartPoses;
	// contextCheckStartPoses.resize(readPaths.size());
	// for (size_t i = 0; i < readPaths.size(); i++)
	// {
	// 	for (auto t : hetInfoPerRead[i])
	// 	{
	// 		assert(std::get<0>(t) < readPaths[i].readLength);
	// 		if (std::get<0>(t)+1 >= resolveLength) contextCheckStartPoses[i].emplace_back(std::get<0>(t)+1-resolveLength);
	// 		assert(std::get<1>(t) >= 0);
	// 		if ((int)std::get<1>(t)+1 < (int)readPaths[i].readLength - (int)resolveLength) contextCheckStartPoses[i].emplace_back(std::get<1>(t)+1);
	// 	}
	// 	for (auto t : chainPositionsInReads[i])
	// 	{
	// 		assert(t.chainStartPosInRead < (int)readPaths[i].readLength);
	// 		if (t.chainStartPosInRead+1 > (int)resolveLength) contextCheckStartPoses[i].emplace_back(t.chainStartPosInRead+1-(int)resolveLength);
	// 		assert(t.chainEndPosInRead >= 0);
	// 		if (t.chainEndPosInRead+1 < (int)readPaths[i].readLength - (int)resolveLength) contextCheckStartPoses[i].emplace_back(t.chainEndPosInRead+1);
	// 	}
	// 	contextCheckStartPoses[i].push_back(0);
	// 	contextCheckStartPoses[i].push_back(readPaths[i].readLength - resolveLength);
	// 	phmap::flat_hash_set<size_t> uniq { contextCheckStartPoses[i].begin(), contextCheckStartPoses[i].end() };
	// 	contextCheckStartPoses[i].clear();
	// 	contextCheckStartPoses[i].insert(contextCheckStartPoses[i].end(), uniq.begin(), uniq.end());
	// 	std::sort(contextCheckStartPoses[i].begin(), contextCheckStartPoses[i].end());
	// }
	// std::cerr << "counting context coverages" << std::endl;
	// for (size_t i = 0; i < readPaths.size(); i++)
	// {
	// 	phmap::flat_hash_set<Context> contexts;
	// 	iterateContextNodes(checkTheseNodes, unitigGraph, readPaths[i], contextCheckStartPoses[i], processedHetInfos[i], chainPositionsInReads[i], resolveLength, [&contexts](const size_t j, const size_t k, const std::vector<Context>& nodecontexts)
	// 	{
	// 		contexts.insert(nodecontexts.begin(), nodecontexts.end());
	// 	});
	// 	for (const auto& context : contexts)
	// 	{
	// 		contextCoverage[context] += 1;
	// 	}
	// }
	// std::cerr << "merging per-node contexts" << std::endl;
	// for (size_t i = 0; i < readPaths.size(); i++)
	// {
	// 	iterateContextNodes(checkTheseNodes, unitigGraph, readPaths[i], contextCheckStartPoses[i], processedHetInfos[i], chainPositionsInReads[i], resolveLength, [&nodeParent, &contextCoverage, &readPaths, i, minContextCoverage](const size_t j, const size_t k, const std::vector<Context>& contexts)
	// 	{
	// 		uint64_t node = readPaths[i].paths[j].path[k];
	// 		Context anyContext;
	// 		for (const auto& context : contexts)
	// 		{
	// 			if (contextCoverage[context] < minContextCoverage) continue;
	// 			anyContext = context;
	// 			find(nodeParent[node & maskUint64_t], context);
	// 		}
	// 		for (const auto& context : contexts)
	// 		{
	// 			if (contextCoverage[context] < minContextCoverage) continue;
	// 			merge(nodeParent[node & maskUint64_t], context, anyContext);
	// 		}
	// 	});
	// }
	// std::cerr << "getting context clusters" << std::endl;
	// UnitigGraph resultGraph = unitigGraph;
	// std::vector<ReadPathBundle> resultPaths = readPaths;
	// std::vector<phmap::flat_hash_map<Context, size_t>> contextToNewNodeNumber;
	// contextToNewNodeNumber.resize(unitigGraph.nodeCount());
	// size_t nextNodeNumber = unitigGraph.nodeCount();
	// for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	// {
	// 	phmap::flat_hash_set<Context> keys;
	// 	for (const auto& pair : nodeParent[i])
	// 	{
	// 		keys.insert(pair.first);
	// 	}
	// 	phmap::flat_hash_set<Context> clusters;
	// 	for (const auto& key : keys)
	// 	{
	// 		clusters.insert(find(nodeParent[i], key));
	// 	}
	// 	for (const auto& cluster : clusters)
	// 	{
	// 		contextToNewNodeNumber[i][cluster] = nextNodeNumber;
	// 		resultGraph.lengths.emplace_back(unitigGraph.lengths[i]);
	// 		resultGraph.coverages.emplace_back(0);
	// 		nextNodeNumber += 1;
	// 	}
	// }
	// std::cerr << "replacing nodes" << std::endl;
	// for (size_t i = 0; i < readPaths.size(); i++)
	// {
	// 	iterateContextNodes(checkTheseNodes, unitigGraph, readPaths[i], contextCheckStartPoses[i], processedHetInfos[i], chainPositionsInReads[i], resolveLength, [&nodeParent, &readPaths, &contextCoverage, &contextToNewNodeNumber, &resultPaths, i, minContextCoverage](const size_t j, const size_t k, const std::vector<Context>& contexts)
	// 	{
	// 		uint64_t node = readPaths[i].paths[j].path[k];
	// 		for (const auto& context : contexts)
	// 		{
	// 			if (contextCoverage[context] < minContextCoverage) continue;
	// 			auto cluster = find(nodeParent[node & maskUint64_t], context);
	// 			if (contextToNewNodeNumber[node & maskUint64_t].count(cluster) == 1)
	// 			{
	// 				// std::cerr << "replace read " << i << " " << j << " " << k << " (" << (node & maskUint64_t) << ") with " << contextToNewNodeNumber[node & maskUint64_t].at(cluster) << std::endl;
	// 				resultPaths[i].paths[j].path[k] = (node & firstBitUint64_t) + (contextToNewNodeNumber[node & maskUint64_t].at(cluster));
	// 				break;
	// 			}
	// 		}
	// 	});
	// }
	RankBitvector kept;
	kept.resize(resultGraph.nodeCount());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (checkTheseNodes[i])
		{
			kept.set(i, false);
		}
		else
		{
			kept.set(i, true);
		}
	}
	for (size_t i = unitigGraph.nodeCount(); i < resultGraph.nodeCount(); i++)
	{
		kept.set(i, true);
	}
	return filterUnitigGraph(resultGraph, resultPaths, kept);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphHapmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t graphk, const size_t resolveLength)
{
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	// for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	// {
	// 	for (auto chain : chainPositionsInReads[i])
	// 	{
	// 		std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
	// 	}
	// }
	SparseEdgeContainer validChainEdges = getValidChainEdges(unitigGraph, anchorChains, chainPositionsInReads, approxOneHapCoverage);
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<ReadDiagonalAlleles>> allelesPerReadPerChain;
	std::tie(allelesPerChain, allelesPerReadPerChain) = getRawReadInfoPerChain(anchorChains, readPaths, unitigGraph, validChainEdges, chainPositionsInReads);
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> hetInfoPerRead = getHetInfoPerRead(unitigGraph, readPaths, anchorChains, allelesPerReadPerChain, readPaths.size(), approxOneHapCoverage);
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipHapmers(false, anchorChains, unitigGraph, readPaths, hetInfoPerRead, chainPositionsInReads, resolveLength, approxOneHapCoverage);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}

phmap::flat_hash_map<uint64_t, std::vector<std::tuple<uint64_t, size_t, size_t>>> getLocallyUniqueNodeEdges(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<bool>& uniques)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<int>> distances;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<int>> distancesWithoutGaps;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		uint64_t lastLocalUniq = std::numeric_limits<size_t>::max();
		int lastLocalUniqEndPos = 0;
		size_t lastAnchorJ = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (!uniques[node & maskUint64_t]) continue;
				if (lastLocalUniq == std::numeric_limits<size_t>::max())
				{
					lastLocalUniq = node;
					lastLocalUniqEndPos = readPos;
					lastAnchorJ = j;
					continue;
				}
				assert(lastAnchorJ <= j);
				auto key = canon(std::make_pair(lastLocalUniq & maskUint64_t, lastLocalUniq & firstBitUint64_t), std::make_pair(node & maskUint64_t, node & firstBitUint64_t));
				int distance = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t] - lastLocalUniqEndPos;
				if (j != lastAnchorJ && distance < 0)
				{
					lastLocalUniq = node;
					lastLocalUniqEndPos = readPos;
					lastAnchorJ = j;
					continue;
				}
				if (key.first == std::make_pair<size_t, bool>(lastLocalUniq & maskUint64_t, lastLocalUniq & firstBitUint64_t) && key.second == std::make_pair<size_t, bool>(node & maskUint64_t, node & firstBitUint64_t))
				{
					if (j != lastAnchorJ) distances[std::make_pair(lastLocalUniq, node)].emplace_back(distance);
					if (j == lastAnchorJ) distancesWithoutGaps[std::make_pair(lastLocalUniq, node)].emplace_back(distance);
					// distancesWithoutGaps[std::make_pair(lastLocalUniq, node)].emplace_back(distance);
				}
				else
				{
					if (j != lastAnchorJ) distances[std::make_pair(node ^ firstBitUint64_t, lastLocalUniq ^ firstBitUint64_t)].emplace_back(distance);
					if (j == lastAnchorJ) distancesWithoutGaps[std::make_pair(node ^ firstBitUint64_t, lastLocalUniq ^ firstBitUint64_t)].emplace_back(distance);
					// distancesWithoutGaps[std::make_pair(node ^ firstBitUint64_t, lastLocalUniq ^ firstBitUint64_t)].emplace_back(distance);
				}
				lastLocalUniq = node;
				lastLocalUniqEndPos = readPos;
				lastAnchorJ = j;
			}
		}
	}
	phmap::flat_hash_map<uint64_t, std::vector<std::tuple<uint64_t, size_t, size_t>>> result;
	for (auto pair : distancesWithoutGaps)
	{
		std::vector<int> dists = pair.second;
		if (distances.count(pair.first) == 1) dists.insert(dists.end(), distances.at(pair.first).begin(), distances.at(pair.first).end());
		if (dists.size() < 3) continue;
		std::sort(dists.begin(), dists.end());
		result[pair.first.first].emplace_back(pair.first.second, dists[dists.size() / 2], dists.size());
		result[pair.first.second ^ firstBitUint64_t].emplace_back(pair.first.first ^ firstBitUint64_t, dists[dists.size()/2], dists.size());
	}
	for (auto pair : distances)
	{
		if (distancesWithoutGaps.count(pair.first) == 1) continue;
		std::vector<int> dists = pair.second;
		if (result.count(pair.first.first) == 1) continue;
		if (result.count(pair.first.second ^ firstBitUint64_t) == 1) continue;
		if (dists.size() < 3) continue;
		std::sort(dists.begin(), dists.end());
		result[pair.first.first].emplace_back(pair.first.second, dists[dists.size() / 2], dists.size());
		result[pair.first.second ^ firstBitUint64_t].emplace_back(pair.first.first ^ firstBitUint64_t, dists[dists.size()/2], dists.size());
	}
	return result;
}

void makeFakeLocalUniqGraph(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<bool>& uniques, const phmap::flat_hash_map<uint64_t, std::vector<std::tuple<uint64_t, size_t, size_t>>>& edges)
{
	std::ofstream graph { "localuniqgraph.gfa" };
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (!uniques[i]) continue;
		graph << "S\t" << i << "\t*\tLN:i:" << unitigGraph.lengths[i] << "\tll:i:" << unitigGraph.coverages[i] << "\tFC:i:" << (unitigGraph.coverages[i]*unitigGraph.lengths[i]) << std::endl;
	}
	size_t edgenum = 0;
	for (const auto& pair : edges)
	{
		uint64_t from = pair.first;
		for (auto t : pair.second)
		{
			uint64_t target = std::get<0>(t);
			if ((target ^ firstBitUint64_t) < from) continue; // only print canonical
			size_t distance = std::get<1>(t);
			if (distance < 1) distance = 1;
			size_t coverage = std::get<2>(t);
			graph << "S\tedge_" << edgenum << "\t*\tLN:i:" << distance << "\tll:i:" << coverage << "\tFC:i:" << (distance * coverage) << std::endl;
			graph << "L\t" << (from & maskUint64_t) << "\t" << ((from & firstBitUint64_t) ? "+" : "-") << "\tedge_" << edgenum << "\t+\t0M\tec:i:" << coverage << std::endl;
			graph << "L\tedge_" << edgenum << "\t+\t" << (target & maskUint64_t) << "\t" << ((target & firstBitUint64_t) ? "+" : "-") << "\t0M\tec:i:" << coverage << std::endl;
			edgenum += 1;
		}
	}
}

AnchorChain getLocallyUniqueChain(const phmap::flat_hash_map<uint64_t, std::vector<std::tuple<uint64_t, size_t, size_t>>>& edges, const uint64_t startNode, const UnitigGraph& unitigGraph)
{
	uint64_t pos = startNode ^ firstBitUint64_t;
	while (edges.count(pos) == 1)
	{
		if (edges.at(pos).size() != 1) break;
		uint64_t next = std::get<0>(edges.at(pos)[0]);
		if (edges.count(next ^ firstBitUint64_t) == 0) break;
		if (edges.at(next ^ firstBitUint64_t).size() != 1) break;
		assert(std::get<0>(edges.at(next ^ firstBitUint64_t)[0]) == (pos ^ firstBitUint64_t));
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if ((next ^ firstBitUint64_t) == startNode) break;
		pos = next;
	}
	pos ^= firstBitUint64_t;
	const uint64_t chainstart = pos;
	AnchorChain result;
	result.nodes.push_back(pos);
	result.nodeOffsets.push_back(0);
	while (edges.count(pos) == 1)
	{
		if (edges.at(pos).size() != 1) break;
		uint64_t next = std::get<0>(edges.at(pos)[0]);
		if (edges.count(next ^ firstBitUint64_t) == 0) break;
		if (edges.at(next ^ firstBitUint64_t).size() != 1) break;
		assert(std::get<0>(edges.at(next ^ firstBitUint64_t)[0]) == (pos ^ firstBitUint64_t));
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if (next == chainstart) break;
		result.nodes.push_back(next);
		result.nodeOffsets.push_back(result.nodeOffsets.back() + std::get<1>(edges.at(pos)[0]) + unitigGraph.lengths[pos & maskUint64_t]);
		pos = next;
	}
	return result;
}

std::vector<AnchorChain> getUniqueChains(const UnitigGraph& unitigGraph, const phmap::flat_hash_map<uint64_t, std::vector<std::tuple<uint64_t, size_t, size_t>>>& edges, const std::vector<bool>& uniques)
{
	std::vector<AnchorChain> result;
	std::vector<bool> checked;
	checked.resize(uniques.size(), false);
	for (size_t i = 0; i < uniques.size(); i++)
	{
		if (checked[i]) continue;
		if (!uniques[i]) continue;
		AnchorChain chain = getLocallyUniqueChain(edges, i, unitigGraph);
		assert(chain.nodes.size() >= 1);
		result.emplace_back(chain);
		for (uint64_t node : chain.nodes)
		{
			assert(!checked[node & maskUint64_t]);
			checked[node & maskUint64_t] = true;
		}
	}
	return result;
}

void filterOutRepetitiveLocallyUniqueNodes(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, std::vector<bool>& uniques, const size_t minLength)
{
	auto edges = getLocallyUniqueNodeEdges(unitigGraph, readPaths, uniques);
	auto fakeChains = getUniqueChains(unitigGraph, edges, uniques);
	size_t countFiltered = 0;
	for (const auto& chain : fakeChains)
	{
		size_t totalChainLength = 0;
		for (auto node : chain.nodes)
		{
			totalChainLength += unitigGraph.lengths[node & maskUint64_t];
		}
		if (totalChainLength >= minLength) continue;
		if (edges.count(chain.nodes.back()) == 0) continue;
		if (edges.count(chain.nodes[0] ^ firstBitUint64_t) == 0) continue;
		if (edges.at(chain.nodes.back()).size() < 2) continue;
		if (edges.at(chain.nodes[0] ^ firstBitUint64_t).size() < 2) continue;
		for (auto node : chain.nodes)
		{
			assert(uniques[node & maskUint64_t]);
			uniques[node & maskUint64_t] = false;
			countFiltered += 1;
		}
	}
	std::cerr << "filtered out " << countFiltered << " nodes" << std::endl;
	if (countFiltered >= 1) filterOutRepetitiveLocallyUniqueNodes(unitigGraph, readPaths, uniques, minLength);
}

std::pair<std::vector<AnchorChain>, std::vector<bool>> getReadLocalUniqChains(const std::vector<AnchorChain>& anchorChains, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t maxCopyCount, const size_t minLength, const size_t uniqSpanLength)
{
	const size_t notSameNodeDistance = 100;
	const size_t branchMinLength = 500;
	std::vector<bool> uniques;
	uniques.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.coverages[i] < approxOneHapCoverage * 0.25) continue;
		if (unitigGraph.coverages[i] > approxOneHapCoverage * ((double)maxCopyCount + 0.5)) continue;
		uniques[i] = true;
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			phmap::flat_hash_map<size_t, size_t> found;
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				size_t node = readPaths[i].paths[j].path[k] & maskUint64_t;
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (found.count(node) == 1)
				{
					if (found.at(node) + uniqSpanLength >= readPos)
					{
						uniques[node & maskUint64_t] = false;
					}
				}
				found[node] = readPos;
			}
		}
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		phmap::flat_hash_map<uint64_t, std::vector<size_t>> nodePositions;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (!uniques[node & maskUint64_t]) continue;
				nodePositions[node].emplace_back(readPos);
			}
		}
		for (auto& pair : nodePositions)
		{
			if (!uniques[pair.first & maskUint64_t]) continue;
			std::sort(pair.second.begin(), pair.second.end());
			for (size_t j = 1; j < pair.second.size(); j++)
			{
				if (pair.second[j] < pair.second[j-1] + notSameNodeDistance) continue;
				if (pair.second[j] > pair.second[j-1] + uniqSpanLength) continue;
				uniques[pair.first & maskUint64_t] = false;
				break;
			}
		}
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.lengths[i] < uniqSpanLength) continue;
		uniques[i] = true;
	}
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 1) continue;
		size_t chainLength = 0;
		for (auto node : anchorChains[i].nodes)
		{
			chainLength += unitigGraph.lengths[node & maskUint64_t];
		}
		if (chainLength < uniqSpanLength) continue;
		for (auto node : anchorChains[i].nodes)
		{
			uniques[node & maskUint64_t] = true;
		}
	}
	filterOutRepetitiveLocallyUniqueNodes(unitigGraph, readPaths, uniques, branchMinLength);
	auto edges = getLocallyUniqueNodeEdges(unitigGraph, readPaths, uniques);
	makeFakeLocalUniqGraph(unitigGraph, readPaths, uniques, edges);
	auto fakeChainsUnfiltered = getUniqueChains(unitigGraph, edges, uniques);
	std::vector<AnchorChain> fakeChains;
	for (const auto& chain : fakeChainsUnfiltered)
	{
		size_t totalChainLength = 0;
		for (uint64_t node : chain.nodes)
		{
			totalChainLength += unitigGraph.lengths[node & maskUint64_t];
		}
		if (totalChainLength < minLength) continue;
		fakeChains.emplace_back(chain);
		std::cerr << "local-uniq chain length: " << totalChainLength << " kmers nodes:";
		for (uint64_t node : chain.nodes)
		{
			std::cerr << " " << ((node & firstBitUint64_t) ? ">" : "<") << (node & maskUint64_t);
		}
		std::cerr << std::endl;
	}
	std::cerr << fakeChains.size() << " local-uniq chains" << std::endl;
	return std::make_pair(fakeChains, std::move(uniques));
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLocalUniqmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t uniqSpanLength, const size_t resolveLength, const size_t maxCopyCount)
{
	// std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	// std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	// {
	// 	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	// 	fixFakeSinglePloidyChains(anchorChains, chainPositionsInReads, approxOneHapCoverage);
	// }
	// std::vector<bool> localUniqs;
	// std::vector<std::vector<ChainPosition>> localUniqPositions;
	// std::tie(localUniqPositions, localUniqs) = getReadLocalUniqChains(anchorChains, unitigGraph, readPaths, approxOneHapCoverage, maxCopyCount, 100, uniqSpanLength);
	// // for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	// // {
	// // 	for (auto chain : chainPositionsInReads[i])
	// // 	{
	// // 		std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
	// // 	}
	// // }
	// std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> hetInfoPerRead;
	// hetInfoPerRead.resize(readPaths.size());
	// UnitigGraph unzippedGraph;
	// std::vector<ReadPathBundle> unzippedReadPaths;
	// std::tie(unzippedGraph, unzippedReadPaths) = unzipHapmers(false, anchorChains, unitigGraph, readPaths, hetInfoPerRead, localUniqPositions, resolveLength, approxOneHapCoverage);
	// RankBitvector kept;
	// kept.resize(unzippedGraph.nodeCount());
	// for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	// {
	// 	kept.set(i, unzippedGraph.coverages[i] >= 2);
	// }
	// return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}

std::vector<std::vector<std::pair<uint64_t, int>>> getReadExtrapolatedLocalUniqLocations(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& localUniqChains)
{
	const size_t maxClusterDistance = 100;
	phmap::flat_hash_map<uint64_t, std::pair<uint64_t, size_t>> nodeUniqChainLocation;
	for (size_t i = 0; i < localUniqChains.size(); i++)
	{
		for (size_t j = 0; j < localUniqChains[i].nodes.size(); j++)
		{
			uint64_t node = localUniqChains[i].nodes[j];
			nodeUniqChainLocation[node] = std::make_pair(i + firstBitUint64_t, j);
			nodeUniqChainLocation[node ^ firstBitUint64_t] = std::make_pair(i, localUniqChains[i].nodes.size()-1-j);
		}
	}
	std::vector<std::vector<size_t>> localUniqCenterOffsets;
	localUniqCenterOffsets.resize(localUniqChains.size());
	for (size_t i = 0; i < localUniqChains.size(); i++)
	{
		localUniqCenterOffsets[i].push_back(unitigGraph.lengths[localUniqChains[i].nodes[0] & maskUint64_t] / 2);
		for (size_t j = 1; j < localUniqChains[i].nodes.size(); j++)
		{
			size_t offset = localUniqChains[i].nodeOffsets[j] + unitigGraph.lengths[localUniqChains[i].nodes[j] & maskUint64_t] / 2;
			assert(offset > localUniqCenterOffsets[i].back());
			localUniqCenterOffsets[i].push_back(offset);
		}
	}
	std::vector<std::vector<std::pair<uint64_t, int>>> result;
	result.resize(readPaths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		phmap::flat_hash_map<uint64_t, std::vector<int>> nodeApproxPositions;
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (nodeUniqChainLocation.count(node) == 0) continue;
				nodeApproxPositions[node].emplace_back((int)readPos - (int)unitigGraph.lengths[node & maskUint64_t]/2);
			}
		}
		std::vector<std::pair<uint64_t, int>> nodePositions;
		for (auto& pair : nodeApproxPositions)
		{
			std::sort(pair.second.begin(), pair.second.end());
			size_t clusterStart = 0;
			for (size_t j = 1; j < pair.second.size(); j++)
			{
				assert(pair.second[j] >= pair.second[j-1]);
				if (pair.second[j] - pair.second[j-1] < (int)maxClusterDistance) continue;
				nodePositions.emplace_back(pair.first, pair.second[clusterStart+(j-clusterStart)/2]);
				clusterStart = j;
			}
			nodePositions.emplace_back(pair.first, pair.second[clusterStart+(pair.second.size()-clusterStart)/2]);
		}
		std::sort(nodePositions.begin(), nodePositions.end(), [](auto left, auto right) { return left.second < right.second; });
		phmap::flat_hash_map<uint64_t, std::vector<std::pair<size_t, int>>> nodesPerChain;
		for (auto pair : nodePositions)
		{
			uint64_t node = pair.first;
			std::pair<uint64_t, size_t> chainPos = nodeUniqChainLocation.at(node);
			nodesPerChain[chainPos.first].emplace_back(chainPos.second, pair.second);
		}
		std::cerr << "read " << i << " real positions:";
		for (auto pos : nodePositions)
		{
			std::cerr << " pos " << pos.second << " node " << pos.first << " (" << (pos.first & maskUint64_t) << ")";
		}
		std::cerr << std::endl;
		std::vector<std::pair<uint64_t, int>> interpolated;
		for (auto& pair : nodesPerChain)
		{
			std::sort(pair.second.begin(), pair.second.end(), [](auto left, auto right) { return left.second < right.second; });
			if (pair.second[0].first != 0)
			{
				std::cerr << "insert starts chain " << ((pair.first & firstBitUint64_t) ? ">" : "<") << (pair.first & maskUint64_t) << " till " << pair.second[0].first << std::endl;
				for (size_t j = 0; j < pair.second[0].first; j++)
				{
					if (pair.first & firstBitUint64_t)
					{
						uint64_t thisNode = localUniqChains[pair.first & maskUint64_t].nodes[j];
						int distance = (int)localUniqCenterOffsets[pair.first & maskUint64_t][pair.second[0].first] - (int)localUniqCenterOffsets[pair.first & maskUint64_t][j];
						assert(distance > 0);
						interpolated.emplace_back(thisNode, pair.second[0].second - distance);
					}
					else
					{
						uint64_t thisNode = localUniqChains[pair.first & maskUint64_t].nodes[localUniqChains[pair.first & maskUint64_t].nodes.size()-1-j] ^ firstBitUint64_t;
						int distance = (int)localUniqCenterOffsets[pair.first & maskUint64_t][localUniqCenterOffsets[pair.first & maskUint64_t].size() - 1 - j] - (int)localUniqCenterOffsets[pair.first & maskUint64_t][localUniqCenterOffsets[pair.first & maskUint64_t].size() - 1 - pair.second[0].first];
						assert(distance > 0);
						interpolated.emplace_back(thisNode, pair.second[0].second - distance);
					}
				}
			}
			for (size_t j = 1; j < pair.second.size(); j++)
			{
				if (pair.second[j].first == pair.second[j-1].first+1) continue;
				if (pair.second[j].first <= pair.second[j-1].first) continue;
				std::cerr << "insert mids chain " << ((pair.first & firstBitUint64_t) ? ">" : "<") << (pair.first & maskUint64_t) << " " << pair.second[j-1].first+1 << " till " << pair.second[j].first << std::endl;
				if (pair.first & firstBitUint64_t)
				{
					int distanceInChain = localUniqCenterOffsets[pair.first & maskUint64_t][pair.second[j].first] - localUniqCenterOffsets[pair.first & maskUint64_t][pair.second[j-1].first];
					assert(distanceInChain > 0);
					int distanceInRead = pair.second[j].second - pair.second[j-1].second;
					assert(distanceInRead > 0);
					for (size_t k = pair.second[j-1].first+1; k < pair.second[j].first; k++)
					{
						int offsetInChain = localUniqCenterOffsets[pair.first & maskUint64_t][k] - localUniqCenterOffsets[pair.first & maskUint64_t][pair.second[j-1].first];
						assert(offsetInChain > 0);
						interpolated.emplace_back(localUniqChains[pair.first & maskUint64_t].nodes[k], pair.second[j-1].second + (double)distanceInRead/(double)distanceInChain * offsetInChain);
					}
				}
				else
				{
					int distanceInChain = localUniqCenterOffsets[pair.first & maskUint64_t][localUniqCenterOffsets[pair.first & maskUint64_t].size()-1-pair.second[j-1].first] - localUniqCenterOffsets[pair.first & maskUint64_t][localUniqCenterOffsets[pair.first & maskUint64_t].size()-1-pair.second[j].first];
					assert(distanceInChain > 0);
					int distanceInRead = pair.second[j].second - pair.second[j-1].second;
					assert(distanceInRead > 0);
					for (size_t k = pair.second[j-1].first+1; k < pair.second[j].first; k++)
					{
						uint64_t thisNode = localUniqChains[pair.first & maskUint64_t].nodes[localUniqChains[pair.first & maskUint64_t].nodes.size()-1-k];
						int offsetInChain = localUniqCenterOffsets[pair.first & maskUint64_t][localUniqCenterOffsets[pair.first & maskUint64_t].size()-1-k] - localUniqCenterOffsets[pair.first & maskUint64_t][localUniqCenterOffsets[pair.first & maskUint64_t].size()-1-pair.second[j].first];
						assert(offsetInChain > 0);
						interpolated.emplace_back(thisNode, pair.second[j].second - (double)distanceInRead/(double)distanceInChain * offsetInChain);
					}
				}
			}
			if (pair.second.back().first+1 < localUniqChains[pair.first & maskUint64_t].nodes.size()) std::cerr << "insert ends chain " << ((pair.first & firstBitUint64_t) ? ">" : "<") << (pair.first & maskUint64_t) << " starting " << pair.second.back().first+1 << std::endl;
			for (size_t j = pair.second.back().first+1; j < localUniqChains[pair.first & maskUint64_t].nodes.size(); j++)
			{
				if (pair.first & firstBitUint64_t)
				{
					uint64_t thisNode = localUniqChains[pair.first & maskUint64_t].nodes[j];
					int distance = localUniqCenterOffsets[pair.first & maskUint64_t][j] - localUniqCenterOffsets[pair.first & maskUint64_t][pair.second.back().first];
					assert(distance > 0);
					interpolated.emplace_back(thisNode, pair.second.back().second + distance);
				}
				else
				{
					uint64_t thisNode = localUniqChains[pair.first & maskUint64_t].nodes[localUniqChains[pair.first & maskUint64_t].nodes.size()-1-j] ^ firstBitUint64_t;
					int distance = localUniqCenterOffsets[pair.first & maskUint64_t][localUniqCenterOffsets[pair.first & maskUint64_t].size()-1-pair.second.back().first] - localUniqCenterOffsets[pair.first & maskUint64_t][localUniqCenterOffsets[pair.first & maskUint64_t].size()-1-j];
					assert(distance > 0);
					interpolated.emplace_back(thisNode, pair.second.back().second + distance);
				}
			}
		}
		std::cerr << "read " << i << " interpolated positions:";
		for (auto pos : interpolated)
		{
			std::cerr << " pos " << pos.second << " node " << pos.first << " (" << (pos.first & maskUint64_t) << ")";
		}
		std::cerr << std::endl;
		result[i].insert(result[i].end(), nodePositions.begin(), nodePositions.end());
		result[i].insert(result[i].end(), interpolated.begin(), interpolated.end());
		std::sort(result[i].begin(), result[i].end(), [](auto left, auto right) { return left.second < right.second; });
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> splitNodesByChainLocation(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& localUniqChains, const std::vector<bool>& localUniq)
{
	std::cerr << "split nodes by chain location" << std::endl;
	auto localUniqNodeLocations = getReadExtrapolatedLocalUniqLocations(unitigGraph, readPaths, localUniqChains);
	size_t maxClusterDistance = 100;
	assert(localUniqNodeLocations.size() == readPaths.size());
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, int>>> nodeReadPositions;
	nodeReadPositions.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			int readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= (int)readPaths[i].paths[j].pathLeftClipKmers;
				if (localUniq[node & maskUint64_t]) continue;
				nodeReadPositions[node & maskUint64_t].emplace_back(i, j, k, readPos);
			}
		}
	}
	for (size_t iii = 0; iii < nodeReadPositions.size(); iii++)
	{
		if (localUniq[iii]) continue;
		phmap::flat_hash_map<uint64_t, std::vector<size_t>> nodePositionsInChains;
		phmap::flat_hash_map<uint64_t, std::vector<size_t>> fwDistances;
		phmap::flat_hash_map<uint64_t, std::vector<size_t>> bwDistances;
		std::vector<std::tuple<size_t, size_t, size_t, uint64_t, uint64_t, size_t, size_t, size_t>> nodeInfos; // i, j, k, bwnode, fwnode, bwdist, fwdist, readpos
		for (auto t : nodeReadPositions[iii])
		{
			size_t i = std::get<0>(t);
			size_t j = std::get<1>(t);
			size_t k = std::get<2>(t);
			int readPos = std::get<3>(t);
			uint64_t node = readPaths[i].paths[j].path[k];
			readPos -= (int)unitigGraph.lengths[node & maskUint64_t] / 2;
			assert((node & maskUint64_t) == iii);
			uint64_t closestBw = std::numeric_limits<size_t>::max();
			uint64_t closestFw = std::numeric_limits<size_t>::max();
			size_t closestBwDistance = std::numeric_limits<size_t>::max();
			size_t closestFwDistance = std::numeric_limits<size_t>::max();
			for (auto pos : localUniqNodeLocations[i])
			{
				if (pos.second < readPos)
				{
					size_t distance = readPos - pos.second;
					if (distance < closestBwDistance)
					{
						closestBwDistance = distance;
						closestBw = pos.first;
					}
				}
				if (pos.second > readPos)
				{
					size_t distance = pos.second - readPos;
					if (distance < closestFwDistance)
					{
						closestFwDistance = distance;
						closestFw = pos.first;
					}
				}
			}
			if (!(node & firstBitUint64_t))
			{
				std::swap(closestFwDistance, closestBwDistance);
				std::swap(closestBw, closestFw);
				if (closestFw != std::numeric_limits<size_t>::max()) closestFw ^= firstBitUint64_t;
				if (closestBw != std::numeric_limits<size_t>::max()) closestBw ^= firstBitUint64_t;
			}
			nodeInfos.emplace_back(i, j, k, closestBw, closestFw, closestBwDistance, closestFwDistance, readPos);
		}
		for (auto t : nodeInfos)
		{
			if (std::get<5>(t) != std::numeric_limits<size_t>::max())
			{
				bwDistances[std::get<3>(t)].emplace_back(std::get<5>(t));
			}
			if (std::get<6>(t) != std::numeric_limits<size_t>::max())
			{
				fwDistances[std::get<4>(t)].emplace_back(std::get<6>(t));
			}
		}
		phmap::flat_hash_map<uint64_t, std::vector<std::pair<size_t, size_t>>> fwSplitters;
		phmap::flat_hash_map<uint64_t, std::vector<std::pair<size_t, size_t>>> bwSplitters;
		for (auto& pair : fwDistances)
		{
			assert(pair.second.size() >= 1);
			std::sort(pair.second.begin(), pair.second.end());
			assert(pair.second[0] >= 0);
			fwSplitters[pair.first].emplace_back(0, pair.second[0]);
			for (size_t j = 1; j < pair.second.size(); j++)
			{
				assert(pair.second[j] >= pair.second[j-1]);
				if (pair.second[j] - pair.second[j-1] > maxClusterDistance)
				{
					assert((pair.second[j] + pair.second[j-1])/2 > fwSplitters[pair.first].back().first);
					fwSplitters[pair.first].emplace_back((pair.second[j] + pair.second[j-1])/2, pair.second[j]);
				}
			}
		}
		for (auto& pair : bwDistances)
		{
			assert(pair.second.size() >= 1);
			std::sort(pair.second.begin(), pair.second.end());
			assert(pair.second[0] >= 0);
			bwSplitters[pair.first].emplace_back(0, pair.second[0]);
			for (size_t j = 1; j < pair.second.size(); j++)
			{
				assert(pair.second[j] >= pair.second[j-1]);
				if (pair.second[j] - pair.second[j-1] > maxClusterDistance)
				{
					assert((pair.second[j] + pair.second[j-1])/2 > bwSplitters[pair.first].back().first);
					bwSplitters[pair.first].emplace_back((pair.second[j] + pair.second[j-1])/2, pair.second[j]);
				}
			}
		}
		phmap::flat_hash_map<size_t, size_t> parent;
		for (auto t : nodeInfos)
		{
			size_t bwAllele = std::numeric_limits<size_t>::max()-1;
			size_t fwAllele = std::numeric_limits<size_t>::max();
			if (std::get<3>(t) != std::numeric_limits<size_t>::max())
			{
				assert(bwSplitters.count(std::get<3>(t)) == 1);
				for (auto pair : bwSplitters.at(std::get<3>(t)))
				{
					if (pair.first > std::get<5>(t)) break;
					assert((std::get<3>(t) & maskUint64_t) < std::numeric_limits<uint32_t>::max());
					assert(pair.second < std::numeric_limits<uint32_t>::max());
					bwAllele = (std::get<3>(t) & firstBitUint64_t) + ((std::get<3>(t) & maskUint64_t) << 32) + pair.second;
				}
				assert(bwAllele != std::numeric_limits<size_t>::max());
				assert((bwAllele & (firstBitUint64_t >> 1)) == 0);
				bwAllele += firstBitUint64_t >> 1;
			}
			if (std::get<4>(t) != std::numeric_limits<size_t>::max())
			{
				assert(fwSplitters.count(std::get<4>(t)) == 1);
				for (auto pair : fwSplitters.at(std::get<4>(t)))
				{
					if (pair.first > std::get<6>(t)) break;
					assert((std::get<4>(t) & maskUint64_t) < std::numeric_limits<uint32_t>::max());
					assert(pair.second < std::numeric_limits<uint32_t>::max());
					fwAllele = (std::get<4>(t) & firstBitUint64_t) + ((std::get<4>(t) & maskUint64_t) << 32) + pair.second;
				}
				assert(fwAllele != std::numeric_limits<size_t>::max());
			}
			if (parent.count(bwAllele) == 0) parent[bwAllele] = bwAllele;
			if (parent.count(fwAllele) == 0) parent[fwAllele] = fwAllele;
			merge(parent, bwAllele, fwAllele);
		}
		phmap::flat_hash_map<size_t, size_t> clusterToNewNode;
		for (auto pair : parent)
		{
			if (pair.second != pair.first) continue;
			clusterToNewNode[pair.first] = resultGraph.lengths.size();
			resultGraph.lengths.emplace_back(resultGraph.lengths[iii]);
			resultGraph.coverages.emplace_back(0);
		}
		// std::cerr << "node " << iii << " has " << clusterToNewNode.size() << " clusters" << std::endl;
		for (auto t : nodeInfos)
		{
			size_t bwAllele = std::numeric_limits<size_t>::max()-1;
			size_t fwAllele = std::numeric_limits<size_t>::max();
			if (std::get<3>(t) != std::numeric_limits<size_t>::max())
			{
				assert(bwSplitters.count(std::get<3>(t)) == 1);
				for (auto pair : bwSplitters.at(std::get<3>(t)))
				{
					if (pair.first > std::get<5>(t)) break;
					assert((std::get<3>(t) & maskUint64_t) < std::numeric_limits<uint32_t>::max());
					assert(pair.second < std::numeric_limits<uint32_t>::max());
					bwAllele = (std::get<3>(t) & firstBitUint64_t) + ((std::get<3>(t) & maskUint64_t) << 32) + pair.second;
				}
				assert(bwAllele != std::numeric_limits<size_t>::max());
				assert((bwAllele & (firstBitUint64_t >> 1)) == 0);
				bwAllele += firstBitUint64_t >> 1;
			}
			if (std::get<4>(t) != std::numeric_limits<size_t>::max())
			{
				assert(fwSplitters.count(std::get<4>(t)) == 1);
				for (auto pair : fwSplitters.at(std::get<4>(t)))
				{
					if (pair.first > std::get<6>(t)) break;
					assert((std::get<4>(t) & maskUint64_t) < std::numeric_limits<uint32_t>::max());
					assert(pair.second < std::numeric_limits<uint32_t>::max());
					fwAllele = (std::get<4>(t) & firstBitUint64_t) + ((std::get<4>(t) & maskUint64_t) << 32) + pair.second;
				}
				assert(fwAllele != std::numeric_limits<size_t>::max());
			}
			auto cluster = find(parent, bwAllele);
			assert(find(parent, fwAllele) == cluster);
			assert(parent.at(cluster) == cluster);
			assert(clusterToNewNode.count(cluster) == 1);
			size_t replacement = clusterToNewNode.at(cluster);
			// std::cerr << "replace read " << std::get<0>(t) << " node " << iii << " pos " << std::get<7>(t) << " with " << replacement << " (bw was " << ((std::get<3>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<3>(t) & maskUint64_t) << " dist " << std::get<5>(t) << ", fw was " << ((std::get<4>(t) & firstBitUint64_t) ? ">" : "<") << (std::get<4>(t) & maskUint64_t) << " dist " << std::get<6>(t) << ")" << std::endl;
			assert((replacement & firstBitUint64_t) == 0);
			resultPaths[std::get<0>(t)].paths[std::get<1>(t)].path[std::get<2>(t)] = replacement + (resultPaths[std::get<0>(t)].paths[std::get<1>(t)].path[std::get<2>(t)] & firstBitUint64_t);
		}
	}
	RankBitvector kept;
	kept.resize(resultGraph.nodeCount());
	for (size_t i = 0; i < kept.size(); i++)
	{
		kept.set(i, true);
	}
	std::cerr << "unitigify" << std::endl;
	return filterUnitigGraph(resultGraph, resultPaths, kept);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLocalUniqmersLocation(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t uniqSpanLength, const size_t maxCopyCount)
{
	std::cerr << "unzip local uniqmer chain location" << std::endl;
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<bool> localUniqs;
	std::vector<AnchorChain> localUniqChains;
	std::tie(localUniqChains, localUniqs) = getReadLocalUniqChains(anchorChains, unitigGraph, readPaths, approxOneHapCoverage, maxCopyCount, 100, uniqSpanLength);
	// for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	// {
	// 	for (auto chain : chainPositionsInReads[i])
	// 	{
	// 		std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
	// 	}
	// }
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> hetInfoPerRead;
	hetInfoPerRead.resize(readPaths.size());
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = splitNodesByChainLocation(unitigGraph, readPaths, localUniqChains, localUniqs);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphChainmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t resolveLength)
{
	std::cerr << "unzip chainmers" << std::endl;
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	// for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	// {
	// 	for (auto chain : chainPositionsInReads[i])
	// 	{
	// 		std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
	// 	}
	// }
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> hetInfoPerRead;
	hetInfoPerRead.resize(readPaths.size());
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::cerr << "unzip chainmers hapmers" << std::endl;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipHapmers(false, anchorChains, unitigGraph, readPaths, hetInfoPerRead, chainPositionsInReads, resolveLength, approxOneHapCoverage);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}

std::pair<std::vector<bool>, std::vector<bool>> getChainTips(const std::vector<AnchorChain>& anchorChains, const UnitigGraph& graph)
{
	auto edges = getActiveEdges(graph.edgeCoverages, graph.nodeCount());
	std::vector<bool> startsWithTip;
	std::vector<bool> endsInTip;
	startsWithTip.resize(anchorChains.size(), false);
	endsInTip.resize(anchorChains.size(), false);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::pair<size_t, bool> start = std::make_pair(anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] ^ firstBitUint64_t) & firstBitUint64_t);
		std::pair<size_t, bool> end = std::make_pair(anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t);
		if (edges.getEdges(start).size() == 0)
		{
			std::cerr << "chain " << i << " start tip" << std::endl;
			startsWithTip[i] = true;
		}
		if (edges.getEdges(end).size() == 0)
		{
			std::cerr << "chain " << i << " end tip" << std::endl;
			endsInTip[i] = true;
		}
	}
	return std::make_pair(startsWithTip, endsInTip);
}

std::vector<std::tuple<uint64_t, uint64_t, size_t>> getGapFills(const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const std::vector<bool>& startsWithTip, const std::vector<bool>& endsInTip, const size_t minSafeCoverage, const size_t maxSpuriousCoverage)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<size_t>> gapfillLengths;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint64_t>> gapfillEdges;
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		for (size_t j = 1; j < chainPositionsInReads[i].size(); j++)
		{
			uint64_t from = chainPositionsInReads[i][j-1].chain;
			uint64_t to = chainPositionsInReads[i][j].chain;
			if (from & firstBitUint64_t)
			{
				if (!endsInTip[from & maskUint64_t]) continue;
			}
			else
			{
				if (!startsWithTip[from & maskUint64_t]) continue;
			}
			if (to & firstBitUint64_t)
			{
				if (!startsWithTip[to & maskUint64_t]) continue;
			}
			else
			{
				if (!endsInTip[to & maskUint64_t]) continue;
			}
			int gapLength = chainPositionsInReads[i][j].chainStartPosInRead - chainPositionsInReads[i][j-1].chainEndPosInRead;
			std::cerr << "read " << i << " gap support chain " << (from & maskUint64_t) << " to " << (to & maskUint64_t) << " length " << gapLength << std::endl;
			if (gapLength < 50) gapLength = 50;
			auto key = canon(from, to);
			gapfillLengths[key].emplace_back(gapLength);
			gapfillEdges[from].emplace(to);
			gapfillEdges[to ^ firstBitUint64_t].emplace(from ^ firstBitUint64_t);
		}
	}
	phmap::flat_hash_map<uint64_t, uint64_t> hasUniqueConnection;
	for (const auto& pair : gapfillEdges)
	{
		std::cerr << "check unique connection chain " << ((pair.first & firstBitUint64_t) ? ">" : "<") << (pair.first & maskUint64_t) << std::endl;
		uint64_t from = pair.first;
		size_t bestCoverage = 0;
		uint64_t bestEdge = 0;
		for (uint64_t to : pair.second)
		{
			auto key = canon(from, to);
			assert(gapfillLengths.count(key) == 1);
			assert(gapfillLengths.at(key).size() >= 1);
			if (gapfillLengths.at(key).size() > bestCoverage)
			{
				bestEdge = to;
				bestCoverage = gapfillLengths.at(key).size();
			}
		}
		if (bestCoverage < minSafeCoverage)
		{
			std::cerr << "not enough coverage (" << bestCoverage << " vs " << minSafeCoverage << ")" << std::endl;
			continue;
		}
		bool valid = true;
		for (uint64_t to : pair.second)
		{
			if (to == bestEdge) continue;
			auto key = canon(from, to);
			assert(gapfillLengths.count(key) == 1);
			assert(gapfillLengths.at(key).size() >= 1);
			if (gapfillLengths.at(key).size() > maxSpuriousCoverage)
			{
				std::cerr << "has other spurious (" << gapfillLengths.at(key).size() << " vs " << maxSpuriousCoverage << ")" << std::endl;
				valid = false;
			}
		}
		if (valid)
		{
			std::cerr << "unique connection from chain " << ((from & firstBitUint64_t) ? ">" : "<") << (from & maskUint64_t) << " to " << ((bestEdge & firstBitUint64_t) ? ">" : "<") << (bestEdge & maskUint64_t) << " coverage " << bestCoverage << std::endl;
			hasUniqueConnection[from] = bestEdge;
		}
	}
	std::vector<std::tuple<uint64_t, uint64_t, size_t>> result;
	for (const auto& pair : gapfillLengths)
	{
		uint64_t from = pair.first.first;
		uint64_t to = pair.first.second;
		if (hasUniqueConnection.count(from) == 0) continue;
		if (hasUniqueConnection.at(from) != to) continue;
		if (hasUniqueConnection.count(to ^ firstBitUint64_t) == 0) continue;
		if (hasUniqueConnection.at(to ^ firstBitUint64_t) != (from ^ firstBitUint64_t)) continue;
		if (pair.second.size() < minSafeCoverage) continue;
		std::cerr << "fill chain " << ((from & firstBitUint64_t) ? ">" : "<") << (from & maskUint64_t) << " to " << ((to & firstBitUint64_t) ? ">" : "<") << (to & maskUint64_t) << std::endl;
		result.emplace_back(from, to, pair.second[pair.second.size()/2]);
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> applyGapFills(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<std::tuple<uint64_t, uint64_t, size_t>>& fills)
{
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	MostlySparse2DHashmap<uint8_t, size_t> newEdgeCoverages;
	newEdgeCoverages.resize(resultGraph.nodeCount() + fills.size());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edges : unitigGraph.edgeCoverages.getValues(fw))
		{
			newEdgeCoverages.set(fw, edges.first, edges.second);
		}
		for (auto edges : unitigGraph.edgeCoverages.getValues(bw))
		{
			newEdgeCoverages.set(bw, edges.first, edges.second);
		}
	}
	size_t zeroFill = unitigGraph.nodeCount();
	for (size_t i = 0; i < fills.size(); i++)
	{
		uint64_t fromnode;
		uint64_t tonode;
		if (std::get<0>(fills[i]) & firstBitUint64_t)
		{
			fromnode = anchorChains[std::get<0>(fills[i]) & maskUint64_t].nodes.back();
		}
		else
		{
			fromnode = anchorChains[std::get<0>(fills[i]) & maskUint64_t].nodes[0] ^ firstBitUint64_t;
		}
		if (std::get<1>(fills[i]) & firstBitUint64_t)
		{
			tonode = anchorChains[std::get<1>(fills[i]) & maskUint64_t].nodes[0];
		}
		else
		{
			tonode = anchorChains[std::get<1>(fills[i]) & maskUint64_t].nodes.back() ^ firstBitUint64_t;
		}
		std::cerr << "fill from node " << ((fromnode & firstBitUint64_t) ? ">" : "<") << (fromnode & maskUint64_t) << " to node " << ((tonode & firstBitUint64_t) ? ">" : "<") << (tonode & maskUint64_t) << " with length " << std::get<2>(fills[i]) << std::endl;
		assert(resultGraph.lengths.size() == zeroFill + i);
		assert(resultGraph.coverages.size() == zeroFill + i);
		resultGraph.lengths.push_back(std::get<2>(fills[i]));
		resultGraph.coverages.push_back(1);
		newEdgeCoverages.set(std::make_pair(fromnode & maskUint64_t, fromnode & firstBitUint64_t), std::make_pair(zeroFill + i, true), 1);
		newEdgeCoverages.set(std::make_pair(zeroFill + i, true), std::make_pair(tonode & maskUint64_t, tonode & firstBitUint64_t), 1);
		// newEdgeCoverages.set(std::make_pair(tonode & maskUint64_t, (tonode ^ firstBitUint64_t) & firstBitUint64_t), std::make_pair(zeroFill + i, false), 1);
		// newEdgeCoverages.set(std::make_pair(zeroFill + i, false), std::make_pair(fromnode & maskUint64_t, (fromnode ^ firstBitUint64_t) & firstBitUint64_t), 1);
	}
	assert(newEdgeCoverages.size() == resultGraph.lengths.size());
	resultGraph.edgeCoverages = newEdgeCoverages;
	return unitigify(resultGraph, resultPaths);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> connectChainGaps(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t minSafeCoverage, const size_t maxSpuriousCoverage)
{
	std::cerr << "try gap fill" << std::endl;
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::cerr << "chain " << i << " ploidy " << anchorChains[i].ploidy << " size " << anchorChains[i].nodes.size();
		for (auto node : anchorChains[i].nodes)
		{
			std::cerr << " " << (node & maskUint64_t);
		}
		std::cerr << std::endl;
	}
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	std::vector<bool> endsInTip;
	std::vector<bool> startsWithTip;
	std::tie(startsWithTip, endsInTip) = getChainTips(anchorChains, unitigGraph);
	auto fills = getGapFills(anchorChains, chainPositionsInReads, startsWithTip, endsInTip, minSafeCoverage, maxSpuriousCoverage);
	UnitigGraph resultGraph;
	std::vector<ReadPathBundle> resultPaths;
	std::tie(resultGraph, resultPaths) = applyGapFills(unitigGraph, readPaths, anchorChains, fills);
	return std::make_pair(std::move(resultGraph), std::move(resultPaths));
}

void checkUniqueTangle(const SparseEdgeContainer& edges, const std::pair<size_t, bool> start, const std::vector<bool>& anchor, VectorWithDirection<size_t>& result, VectorWithDirection<bool>& checked, const size_t tangleNumber)
{
	assert(result[start] == std::numeric_limits<size_t>::max());
	std::vector<uint64_t> stack;
	stack.push_back(start.first + (start.second ? firstBitUint64_t : 0));
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (checked[std::make_pair(top & maskUint64_t, top & firstBitUint64_t)]) continue;
		std::cerr << "nodetip " << ((top & firstBitUint64_t) ? ">" : "<") << (top & maskUint64_t) << " is in tangle " << tangleNumber << std::endl;
		assert(result[std::make_pair(top & maskUint64_t, top & firstBitUint64_t)] == std::numeric_limits<size_t>::max());
		result[std::make_pair(top & maskUint64_t, top & firstBitUint64_t)] = tangleNumber;
		checked[std::make_pair(top & maskUint64_t, top & firstBitUint64_t)] = true;
		if (!anchor[top & maskUint64_t]) stack.push_back(top ^ firstBitUint64_t);
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			stack.push_back(edge.first + (edge.second ? 0 : firstBitUint64_t));
		}
	}
}

VectorWithDirection<size_t> getNodeTipTangleAssignmentsUniqueChains(const UnitigGraph& unitigGraph, const std::vector<AnchorChain>& anchorChains)
{
	auto edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<bool> anchor;
	anchor.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 1) continue;
		anchor[anchorChains[i].nodes[0] & maskUint64_t] = true;
		anchor[anchorChains[i].nodes.back() & maskUint64_t] = true;
	}
	VectorWithDirection<size_t> result;
	VectorWithDirection<bool> checked;
	checked.resize(unitigGraph.nodeCount(), false);
	result.resize(unitigGraph.nodeCount(), std::numeric_limits<size_t>::max());
	size_t nextTangleNum = 0;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::pair<size_t, bool> bw { anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] & firstBitUint64_t) ^ firstBitUint64_t };
		std::pair<size_t, bool> fw { anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t };
		if (!checked[bw])
		{
			checkUniqueTangle(edges, bw, anchor, result, checked, nextTangleNum);
			nextTangleNum += 1;
		}
		if (!checked[fw])
		{
			checkUniqueTangle(edges, fw, anchor, result, checked, nextTangleNum);
			nextTangleNum += 1;
		}
	}
	return result;
}

phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> getTangleSpanners(const bool removeContained, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, VectorWithDirection<size_t>& nodeTipBelongsToTangle, const size_t minValidCoverage, const size_t maxSpuriousCoverage, const std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> spannerCoverages;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint64_t>> spannerEdges;
	for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	{
		std::vector<ChainPosition> uniqueChains;
		for (size_t j = 0; j < chainPositionsInReads[i].size(); j++)
		{
			if (anchorChains[chainPositionsInReads[i][j].chain & maskUint64_t].ploidy != 1) continue;
			size_t chain = chainPositionsInReads[i][j].chain & maskUint64_t;
			if (removeContained)
			{
				if (nodeTipBelongsToTangle[std::make_pair(anchorChains[chain].nodes.back() & maskUint64_t, anchorChains[chain].nodes.back() & firstBitUint64_t)] == nodeTipBelongsToTangle[std::make_pair(anchorChains[chain].nodes[0] & maskUint64_t, (anchorChains[chain].nodes[0] & firstBitUint64_t) ^ firstBitUint64_t)]) continue;
			}
			uniqueChains.emplace_back(chainPositionsInReads[i][j]);
		}
		for (size_t j = 1; j < uniqueChains.size(); j++)
		{
			uint64_t from = uniqueChains[j-1].chain;
			uint64_t to = uniqueChains[j].chain;
			assert(anchorChains[from & maskUint64_t].ploidy == 1);
			assert(anchorChains[to & maskUint64_t].ploidy == 1);
			int gapLength = uniqueChains[j].chainStartPosInRead - uniqueChains[j-1].chainEndPosInRead;
			if (gapLength < 50) gapLength = 50;
			auto key = canon(from, to);
			spannerCoverages[key] += 1;
			spannerEdges[from].emplace(to);
			spannerEdges[to ^ firstBitUint64_t].emplace(from ^ firstBitUint64_t);
		}
	}
	phmap::flat_hash_map<uint64_t, uint64_t> hasUniqueConnection;
	for (const auto& pair : spannerEdges)
	{
		uint64_t from = pair.first;
		size_t bestCoverage = 0;
		uint64_t bestEdge = 0;
		for (uint64_t to : pair.second)
		{
			auto key = canon(from, to);
			assert(spannerCoverages.count(key) == 1);
			assert(spannerCoverages.at(key) >= 1);
			if (spannerCoverages.at(key) > bestCoverage)
			{
				bestEdge = to;
				bestCoverage = spannerCoverages.at(key);
			}
		}
		if (bestCoverage < minValidCoverage) continue;
		bool valid = true;
		for (uint64_t to : pair.second)
		{
			if (to == bestEdge) continue;
			auto key = canon(from, to);
			assert(spannerCoverages.count(key) == 1);
			assert(spannerCoverages.at(key) >= 1);
			if (spannerCoverages.at(key) > maxSpuriousCoverage) valid = false;
		}
		if (!valid) continue;
		hasUniqueConnection[from] = bestEdge;
		std::cerr << "unique connection from " << ((from & firstBitUint64_t) ? ">" : "<") << (from & maskUint64_t) << " to " << ((bestEdge & firstBitUint64_t) ? ">" : "<") << (bestEdge & maskUint64_t) << std::endl;
	}
	phmap::flat_hash_map<size_t, size_t> tangleMerge;
	for (size_t i = 0; i < nodeTipBelongsToTangle.size(); i++)
	{
		if (nodeTipBelongsToTangle.at(std::make_pair(i, true)) != std::numeric_limits<size_t>::max())
		{
			tangleMerge[nodeTipBelongsToTangle.at(std::make_pair(i, true))] = nodeTipBelongsToTangle.at(std::make_pair(i, true));
		}
		if (nodeTipBelongsToTangle.at(std::make_pair(i, false)) != std::numeric_limits<size_t>::max())
		{
			tangleMerge[nodeTipBelongsToTangle.at(std::make_pair(i, false))] = nodeTipBelongsToTangle.at(std::make_pair(i, false));
		}
	}
	bool tanglesNeedFixing = false;
	for (auto pair : hasUniqueConnection)
	{
		if (hasUniqueConnection.count(pair.second ^ firstBitUint64_t) == 0) continue;
		if (hasUniqueConnection.at(pair.second ^ firstBitUint64_t) != (pair.first ^ firstBitUint64_t)) continue;
		size_t nodeFrom;
		size_t nodeTo;
		if (pair.first & firstBitUint64_t)
		{
			nodeFrom = anchorChains[pair.first & maskUint64_t].nodes.back();
		}
		else
		{
			nodeFrom = anchorChains[pair.first & maskUint64_t].nodes[0] ^ firstBitUint64_t;
		}
		uint64_t otherChain = pair.second ^ firstBitUint64_t;
		if (otherChain & firstBitUint64_t)
		{
			nodeTo = anchorChains[otherChain & maskUint64_t].nodes.back();
		}
		else
		{
			nodeTo = anchorChains[otherChain & maskUint64_t].nodes[0] ^ firstBitUint64_t;
		}
		size_t tangleFrom = nodeTipBelongsToTangle[std::make_pair(nodeFrom & maskUint64_t, nodeFrom & firstBitUint64_t)];
		size_t tangleTo = nodeTipBelongsToTangle[std::make_pair(nodeTo & maskUint64_t, nodeTo & firstBitUint64_t)];
		assert(tangleFrom != std::numeric_limits<size_t>::max());
		assert(tangleTo != std::numeric_limits<size_t>::max());
		if (tangleFrom == tangleTo) continue;
		merge(tangleMerge, tangleFrom, tangleTo);
		tanglesNeedFixing = true;
	}
	if (tanglesNeedFixing)
	{
		for (size_t i = 0; i < nodeTipBelongsToTangle.size(); i++)
		{
			if (nodeTipBelongsToTangle[std::make_pair(i, true)] != std::numeric_limits<size_t>::max()) nodeTipBelongsToTangle[std::make_pair(i, true)] = find(tangleMerge, nodeTipBelongsToTangle[std::make_pair(i, true)]);
			if (nodeTipBelongsToTangle[std::make_pair(i, false)] != std::numeric_limits<size_t>::max()) nodeTipBelongsToTangle[std::make_pair(i, false)] = find(tangleMerge, nodeTipBelongsToTangle[std::make_pair(i, false)]);
		}
	}
	std::vector<std::vector<uint64_t>> chainsPerTangle;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 1) continue;
		std::pair<size_t, bool> bw { anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] & firstBitUint64_t) ^ firstBitUint64_t };
		std::pair<size_t, bool> fw { anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t };
		assert(nodeTipBelongsToTangle[bw] != std::numeric_limits<size_t>::max());
		if (removeContained)
		{
			if (nodeTipBelongsToTangle[bw] == nodeTipBelongsToTangle[fw]) continue;
		}
		while (chainsPerTangle.size() <= nodeTipBelongsToTangle[bw])
		{
			chainsPerTangle.emplace_back();
		}
		chainsPerTangle[nodeTipBelongsToTangle[bw]].emplace_back(i);
		assert(nodeTipBelongsToTangle[fw] != std::numeric_limits<size_t>::max());
		while (chainsPerTangle.size() <= nodeTipBelongsToTangle[fw])
		{
			chainsPerTangle.emplace_back();
		}
		chainsPerTangle[nodeTipBelongsToTangle[fw]].emplace_back(i + firstBitUint64_t);
	}
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> result;
	for (size_t i = 0; i < chainsPerTangle.size(); i++)
	{
		bool allValid = true;
		for (uint64_t tip : chainsPerTangle[i])
		{
			if (hasUniqueConnection.count(tip) == 0)
			{
				allValid = false;
				break;
			}
			uint64_t otherChain = hasUniqueConnection.at(tip) ^ firstBitUint64_t;
			uint64_t otherNode;
			if (otherChain & firstBitUint64_t)
			{
				otherNode = anchorChains[otherChain & maskUint64_t].nodes.back();
			}
			else
			{
				otherNode = anchorChains[otherChain & maskUint64_t].nodes[0] ^ firstBitUint64_t;
			}
			if (hasUniqueConnection.count(hasUniqueConnection.at(tip) ^ firstBitUint64_t) == 0)
			{
				allValid = false;
				break;
			}
			if (hasUniqueConnection.at(hasUniqueConnection.at(tip) ^ firstBitUint64_t) != (tip ^ firstBitUint64_t))
			{
				allValid = false;
				break;
			}
			assert(nodeTipBelongsToTangle[std::make_pair(otherNode & maskUint64_t, otherNode & firstBitUint64_t)] == i);
		}
		if (!allValid) continue;
		for (uint64_t tip : chainsPerTangle[i])
		{
			uint64_t other = hasUniqueConnection.at(tip);
			auto key = canon(tip, other);
			result.emplace(key);
		}
	}
	return result;
}

std::vector<std::tuple<int, int, size_t, std::pair<uint64_t, uint64_t>, bool>> getReplaceableSpans(const bool removeContained, const std::vector<AnchorChain>& anchorChains, const VectorWithDirection<size_t>& nodeTipBelongsToTangle, phmap::flat_hash_set<size_t>& clearedTangles, const ReadPathBundle& readPath, const std::vector<ChainPosition>& chainPositionsInRead, const phmap::flat_hash_set<std::pair<uint64_t, uint64_t>>& validSpans)
{
	std::vector<std::tuple<int, int, size_t, std::pair<uint64_t, uint64_t>, bool>> replaceableSpans;
	std::vector<ChainPosition> uniqueChains;
	for (size_t j = 0; j < chainPositionsInRead.size(); j++)
	{
		if (anchorChains[chainPositionsInRead[j].chain & maskUint64_t].ploidy != 1) continue;
		if (removeContained)
		{
			if (nodeTipBelongsToTangle[std::make_pair(anchorChains[chainPositionsInRead[j].chain].nodes.back() & maskUint64_t, anchorChains[chainPositionsInRead[j].chain].nodes.back() & firstBitUint64_t)] == nodeTipBelongsToTangle[std::make_pair(anchorChains[chainPositionsInRead[j].chain].nodes[0] & maskUint64_t, (anchorChains[chainPositionsInRead[j].chain].nodes[0] & firstBitUint64_t) ^ firstBitUint64_t)]) continue;
		}
		uniqueChains.emplace_back(chainPositionsInRead[j]);
	}
	for (size_t j = 1; j < uniqueChains.size(); j++)
	{
		auto key = canon(uniqueChains[j-1].chain, uniqueChains[j].chain);
		if (validSpans.count(key) == 0) continue;
		bool fw = (key.first == uniqueChains[j-1].chain && key.second == uniqueChains[j].chain);
		uint64_t prevNode;
		uint64_t thisNode;
		if (uniqueChains[j-1].chain & firstBitUint64_t)
		{
			prevNode = anchorChains[uniqueChains[j-1].chain & maskUint64_t].nodes.back();
		}
		else
		{
			prevNode = anchorChains[uniqueChains[j-1].chain & maskUint64_t].nodes[0] ^ firstBitUint64_t;
		}
		if (uniqueChains[j].chain & firstBitUint64_t)
		{
			thisNode = anchorChains[uniqueChains[j].chain & maskUint64_t].nodes[0] ^ firstBitUint64_t;
		}
		else
		{
			thisNode = anchorChains[uniqueChains[j].chain & maskUint64_t].nodes.back();
		}
		assert(nodeTipBelongsToTangle[std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)] == nodeTipBelongsToTangle[std::make_pair(thisNode & maskUint64_t, thisNode & firstBitUint64_t)]);
		size_t tangle = nodeTipBelongsToTangle[std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t)];
		clearedTangles.insert(tangle);
		assert(tangle != std::numeric_limits<size_t>::max());
		replaceableSpans.emplace_back((uniqueChains[j-1].chainEndPosInRead + uniqueChains[j-1].chainStartPosInRead)/2, (uniqueChains[j].chainEndPosInRead + uniqueChains[j].chainStartPosInRead)/2, tangle, key, fw);
	}
	return replaceableSpans;
}

phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<std::pair<double, size_t>>>> getTangleSpannedNodeSeparatorReplacers(UnitigGraph& resultGraph, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, phmap::flat_hash_map<uint64_t, std::vector<double>>>& nodeApproxPositionsWithinTangleSpanners, const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<size_t>>& spannerLengths)
{
	const size_t maxClusterDistance = 100;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<std::pair<double, size_t>>>> result;
	for (auto& pair : nodeApproxPositionsWithinTangleSpanners)
	{
		assert(spannerLengths.count(pair.first) == 1);
		size_t lengthSum = 0;
		for (auto len : spannerLengths.at(pair.first))
		{
			lengthSum += len;
		}
		size_t length = lengthSum / spannerLengths.at(pair.first).size();
		for (auto& pair2 : pair.second)
		{
			assert((pair2.first & firstBitUint64_t) == 0);
			assert(pair2.second.size() >= 1);
			std::sort(pair2.second.begin(), pair2.second.end());
			result[pair2.first][pair.first].emplace_back(pair2.second[0] - 1.0, resultGraph.lengths.size());
			// std::cerr << "breakpoint node " << pair2.first << " spanners " << pair.first.first << " " << pair.first.second << " breakpoint " << pair2.second[0] - 1.0 << std::endl;
			resultGraph.lengths.emplace_back(resultGraph.lengths[pair2.first & maskUint64_t]);
			resultGraph.coverages.emplace_back(0);
			for (size_t i = 1; i < pair2.second.size(); i++)
			{
				if ((pair2.second[i] - pair2.second[i-1]) * length < maxClusterDistance + resultGraph.lengths[pair2.first & maskUint64_t]) continue;
				// std::cerr << "breakpoint node " << pair2.first << " spanners " << pair.first.first << " " << pair.first.second << " breakpoint " << (pair2.second[i] + pair2.second[i-1]) / 2.0 << std::endl;
				result[pair2.first][pair.first].emplace_back((pair2.second[i] + pair2.second[i-1]) / 2.0, resultGraph.lengths.size());
				resultGraph.lengths.emplace_back(resultGraph.lengths[pair2.first & maskUint64_t]);
				resultGraph.coverages.emplace_back(0);
			}
		}
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> applyTangleSpanners(const bool removeContained, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const VectorWithDirection<size_t>& nodeTipBelongsToTangle, const phmap::flat_hash_set<std::pair<uint64_t, uint64_t>>& validSpans, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const std::vector<AnchorChain>& anchorChains)
{
	std::cerr << validSpans.size() << " valid spans" << std::endl;
	for (auto span : validSpans)
	{
		std::cerr << "span from " << ((span.first & firstBitUint64_t) ? ">" : "<") << (span.first & maskUint64_t) << " " << ((span.second & firstBitUint64_t) ? ">" : "<") << (span.second & maskUint64_t) << std::endl;
	}
	assert(readPaths.size() == chainPositionsInReads.size());
	phmap::flat_hash_map<uint64_t, uint64_t> nodeIsAnchorChainTip;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 1) continue;
		assert(nodeIsAnchorChainTip.count(anchorChains[i].nodes[0] ^ firstBitUint64_t) == 0);
		nodeIsAnchorChainTip[anchorChains[i].nodes[0] ^ firstBitUint64_t] = i;
		assert(nodeIsAnchorChainTip.count(anchorChains[i].nodes.back()) == 0);
		nodeIsAnchorChainTip[anchorChains[i].nodes.back()] = i + firstBitUint64_t;
	}
	phmap::flat_hash_set<size_t> mustNotReplace;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 1) continue;
		if (anchorChains[i].nodes.size() != 1) continue;
		if (removeContained)
		{
			if (nodeTipBelongsToTangle[std::make_pair(anchorChains[i].nodes.back() & maskUint64_t, anchorChains[i].nodes.back() & firstBitUint64_t)] == nodeTipBelongsToTangle[std::make_pair(anchorChains[i].nodes[0] & maskUint64_t, (anchorChains[i].nodes[0] & firstBitUint64_t) ^ firstBitUint64_t)]) continue;
		}
		mustNotReplace.insert(anchorChains[i].nodes[0] & maskUint64_t);
	}
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	phmap::flat_hash_set<size_t> clearedTangles;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, phmap::flat_hash_map<uint64_t, std::vector<double>>> nodeApproxPositionsWithinTangleSpanners;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<size_t>> spannerLengths;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<int, int, size_t, std::pair<uint64_t, uint64_t>, bool>> replaceableSpans = getReplaceableSpans(removeContained, anchorChains, nodeTipBelongsToTangle, clearedTangles, readPaths[i], chainPositionsInReads[i], validSpans);
		if (replaceableSpans.size() == 0) continue;
		for (size_t j = 0; j < resultPaths[i].paths.size(); j++)
		{
			size_t readPos = resultPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < resultPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = resultPaths[i].paths[j].path[k];
				readPos += resultGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= resultPaths[i].paths[j].pathLeftClipKmers;
				if (mustNotReplace.count(node & maskUint64_t) == 1) continue;
				size_t tangle = std::numeric_limits<size_t>::max();
				std::pair<uint64_t, uint64_t> replaceKey;
				int spannerStartPos, spannerEndPos;
				bool spannerFw;
				for (auto t : replaceableSpans)
				{
					if (std::get<0>(t) > (int)readPos - (int)resultGraph.lengths[node & maskUint64_t]) continue;
					if (std::get<1>(t) < (int)readPos) continue;
					if (nodeTipBelongsToTangle[std::make_pair(node & maskUint64_t, false)] != std::get<2>(t)) continue;
					if (nodeTipBelongsToTangle[std::make_pair(node & maskUint64_t, true)] != std::get<2>(t)) continue;
					if (tangle == std::numeric_limits<size_t>::max())
					{
						tangle = std::get<2>(t);
						spannerStartPos = std::get<0>(t);
						spannerEndPos = std::get<1>(t);
						spannerFw = std::get<4>(t);
					}
					else
					{
						tangle = std::numeric_limits<size_t>::max()-1;
					}
					replaceKey = std::get<3>(t);
				}
				if (tangle >= std::numeric_limits<size_t>::max()-1) continue;
				assert(spannerEndPos > spannerStartPos);
				int nodeApproxPos = (int)readPos - (int)(resultGraph.lengths[node & maskUint64_t])/2;
				double nodeFractionalPos = (double)(nodeApproxPos - spannerStartPos) / (double)(spannerEndPos - spannerStartPos);
				if (spannerFw)
				{
					// std::cerr << "found read " << i << " node " << (node & maskUint64_t) << " pos " << nodeFractionalPos << " fw spanners " << replaceKey.first << " " << replaceKey.second << std::endl;
					nodeApproxPositionsWithinTangleSpanners[replaceKey][node & maskUint64_t].emplace_back(nodeFractionalPos);
				}
				else
				{
					nodeFractionalPos = 1.0 - nodeFractionalPos;
					// std::cerr << "found read " << i << " node " << (node & maskUint64_t) << " pos " << nodeFractionalPos << " bw spanners " << replaceKey.first << " " << replaceKey.second << std::endl;
					nodeApproxPositionsWithinTangleSpanners[replaceKey][node & maskUint64_t].emplace_back(nodeFractionalPos);
				}
				spannerLengths[replaceKey].emplace_back(spannerEndPos - spannerStartPos);
			}
		}
	}
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<std::pair<double, size_t>>>> replaceNodeWithThisNode = getTangleSpannedNodeSeparatorReplacers(resultGraph, nodeApproxPositionsWithinTangleSpanners, spannerLengths);
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		std::vector<std::tuple<int, int, size_t, std::pair<uint64_t, uint64_t>, bool>> replaceableSpans = getReplaceableSpans(removeContained, anchorChains, nodeTipBelongsToTangle, clearedTangles, readPaths[i], chainPositionsInReads[i], validSpans);
		if (replaceableSpans.size() == 0) continue;
		for (size_t j = 0; j < resultPaths[i].paths.size(); j++)
		{
			size_t readPos = resultPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < resultPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = resultPaths[i].paths[j].path[k];
				readPos += resultGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= resultPaths[i].paths[j].pathLeftClipKmers;
				if (mustNotReplace.count(node & maskUint64_t) == 1) continue;
				size_t tangle = std::numeric_limits<size_t>::max();
				std::pair<uint64_t, uint64_t> replaceKey;
				int spannerStartPos, spannerEndPos;
				bool spannerFw;
				for (auto t : replaceableSpans)
				{
					if (std::get<0>(t) > (int)readPos - (int)resultGraph.lengths[node & maskUint64_t]) continue;
					if (std::get<1>(t) < (int)readPos) continue;
					if (nodeTipBelongsToTangle[std::make_pair(node & maskUint64_t, false)] != std::get<2>(t)) continue;
					if (nodeTipBelongsToTangle[std::make_pair(node & maskUint64_t, true)] != std::get<2>(t)) continue;
					if (tangle == std::numeric_limits<size_t>::max())
					{
						tangle = std::get<2>(t);
						spannerStartPos = std::get<0>(t);
						spannerEndPos = std::get<1>(t);
						spannerFw = std::get<4>(t);
					}
					else
					{
						tangle = std::numeric_limits<size_t>::max()-1;
					}
					replaceKey = std::get<3>(t);
				}
				if (tangle >= std::numeric_limits<size_t>::max()-1) continue;
				assert(spannerEndPos > spannerStartPos);
				int nodeApproxPos = (int)readPos - (int)(resultGraph.lengths[node & maskUint64_t])/2;
				double nodeFractionalPos = (double)(nodeApproxPos - spannerStartPos) / (double)(spannerEndPos - spannerStartPos);
				size_t replacement = std::numeric_limits<size_t>::max();
				if (spannerFw)
				{
					for (auto pair : replaceNodeWithThisNode.at(node & maskUint64_t).at(replaceKey))
					{
						if (pair.first <= nodeFractionalPos) replacement = pair.second;
					}
					assert(replacement != std::numeric_limits<size_t>::max());
				}
				if (!spannerFw)
				{
					nodeFractionalPos = 1.0 - nodeFractionalPos;
					for (auto pair : replaceNodeWithThisNode.at(node & maskUint64_t).at(replaceKey))
					{
						if (pair.first <= nodeFractionalPos)
						{
							replacement = pair.second;
						}
					}
					assert(replacement != std::numeric_limits<size_t>::max());
				}
				assert((replacement & firstBitUint64_t) == 0);
				// std::cerr << "replace read " << i << " node " << (node & maskUint64_t) << " pos " << nodeFractionalPos << " " << (spannerFw ? "fw" : "bw") << " spanners " << replaceKey.first << " " << replaceKey.second << " with " << replacement << std::endl;
				resultPaths[i].paths[j].path[k] = replacement + (node & firstBitUint64_t);
			}
		}
	}
	for (size_t i = 0; i < resultPaths.size(); i++)
	{
		phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> cutPathBeforeHere;
		for (size_t j = 0; j < resultPaths[i].paths.size(); j++)
		{
			size_t readPos = resultPaths[i].paths[j].readStartPos;
			readPos += resultGraph.lengths[resultPaths[i].paths[j].path[0] & maskUint64_t];
			readPos -= resultPaths[i].paths[j].pathLeftClipKmers;
			for (size_t k = 1; k < resultPaths[i].paths[j].path.size(); k++)
			{
				readPos += resultGraph.lengths[resultPaths[i].paths[j].path[k] & maskUint64_t];
				if (nodeIsAnchorChainTip.count(resultPaths[i].paths[j].path[k-1]) == 0) continue;
				if (nodeIsAnchorChainTip.count(resultPaths[i].paths[j].path[k] ^ firstBitUint64_t) == 0) continue;
				if (nodeTipBelongsToTangle[std::make_pair(resultPaths[i].paths[j].path[k-1] & maskUint64_t, resultPaths[i].paths[j].path[k-1] & firstBitUint64_t)] != nodeTipBelongsToTangle[std::make_pair(resultPaths[i].paths[j].path[k] & maskUint64_t, (resultPaths[i].paths[j].path[k] ^ firstBitUint64_t) & firstBitUint64_t)]) continue;
				if (nodeTipBelongsToTangle[std::make_pair(resultPaths[i].paths[j].path[k-1] & maskUint64_t, resultPaths[i].paths[j].path[k-1] & firstBitUint64_t)] == std::numeric_limits<size_t>::max()) continue;
				if (clearedTangles.count(nodeTipBelongsToTangle[std::make_pair(resultPaths[i].paths[j].path[k-1] & maskUint64_t, resultPaths[i].paths[j].path[k-1] & firstBitUint64_t)]) == 0) continue;
				uint64_t chainFrom = nodeIsAnchorChainTip.at(resultPaths[i].paths[j].path[k-1]);
				uint64_t chainTo = nodeIsAnchorChainTip.at(resultPaths[i].paths[j].path[k] ^ firstBitUint64_t) ^ firstBitUint64_t;
				auto key = canon(chainFrom, chainTo);
				if (validSpans.count(key) == 1) continue;
				assert(readPos >= resultGraph.lengths[resultPaths[i].paths[j].path[k] & maskUint64_t]);
				size_t posAtStartOfNode = readPos - resultGraph.lengths[resultPaths[i].paths[j].path[k] & maskUint64_t];
				std::cerr << "should cut read " << i << " pos " << posAtStartOfNode << std::endl;
				cutPathBeforeHere[std::make_pair(j, k)] = posAtStartOfNode;
			}
		}
		if (cutPathBeforeHere.size() == 0) continue;
		bool didcut = false;
		for (size_t j = resultPaths[i].paths.size()-1; j < resultPaths[i].paths.size(); j--)
		{
			for (size_t k = resultPaths[i].paths[j].path.size()-1; k > 0; k--)
			{
				if (cutPathBeforeHere.count(std::make_pair(j, k)) == 0) continue;
				assert(cutPathBeforeHere.at(std::make_pair(j, k)) < readPaths[i].readLength);
				didcut = true;
				std::cerr << "cut read " << i << " pos " << cutPathBeforeHere.at(std::make_pair(j, k)) << std::endl;
				resultPaths[i].paths.emplace_back();
				resultPaths[i].paths.back().readStartPos = cutPathBeforeHere.at(std::make_pair(j, k));
				resultPaths[i].paths.back().pathLeftClipKmers = 0;
				resultPaths[i].paths.back().pathRightClipKmers = resultPaths[i].paths[j].pathRightClipKmers;
				resultPaths[i].paths.back().path.insert(resultPaths[i].paths.back().path.end(), resultPaths[i].paths[j].path.begin()+k, resultPaths[i].paths[j].path.end());
				resultPaths[i].paths[j].path.erase(resultPaths[i].paths[j].path.begin() + k, resultPaths[i].paths[j].path.end());
				resultPaths[i].paths[j].pathRightClipKmers = 0;
			}
		}
		assert(didcut);
	}
	std::cerr << clearedTangles.size() << " cleared tangles" << std::endl;
	RankBitvector kept;
	kept.resize(resultGraph.nodeCount());
	for (size_t i = 0; i < resultGraph.nodeCount(); i++)
	{
		kept.set(i, true);
		if (i >= unitigGraph.nodeCount()) continue;
		if (mustNotReplace.count(i) == 1) continue;
		if (nodeTipBelongsToTangle[std::make_pair(i, false)] != nodeTipBelongsToTangle[std::make_pair(i, true)]) continue;
		if (clearedTangles.count(nodeTipBelongsToTangle[std::make_pair(i, true)]) != 1) continue;
		kept.set(i, false);
	}
	return filterUnitigGraph(resultGraph, resultPaths, kept);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolveSpannedTanglesBipartite(const bool removeContained, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage)
{
	std::cerr << "try bipartite spanning" << std::endl;
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	fixFakeSinglePloidyChains(anchorChains, chainPositionsInReads, approxOneHapCoverage);
	VectorWithDirection<size_t> nodeTipBelongsToTangle = getNodeTipTangleAssignmentsUniqueChains(unitigGraph, anchorChains);
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> validSpans = getTangleSpanners(removeContained, unitigGraph, readPaths, nodeTipBelongsToTangle, 2, 1, anchorChains, chainPositionsInReads);
	UnitigGraph resultGraph;
	std::vector<ReadPathBundle> resultPaths;
	std::tie(resultGraph, resultPaths) = applyTangleSpanners(false, unitigGraph, readPaths, nodeTipBelongsToTangle, validSpans, chainPositionsInReads, anchorChains);
	return std::make_pair(std::move(resultGraph), std::move(resultPaths));
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolveSpannedTangles(const bool removeContained, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage)
{
	std::cerr << "try gap fill" << std::endl;
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	fixFakeSinglePloidyChains(anchorChains, chainPositionsInReads, approxOneHapCoverage);
	VectorWithDirection<size_t> nodeTipBelongsToTangle = getNodeTipTangleAssignmentsUniqueChains(unitigGraph, anchorChains);
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> validSpans = getTangleSpanners(removeContained, unitigGraph, readPaths, nodeTipBelongsToTangle, approxOneHapCoverage * 0.5, 3, anchorChains, chainPositionsInReads);
	UnitigGraph resultGraph;
	std::vector<ReadPathBundle> resultPaths;
	std::tie(resultGraph, resultPaths) = applyTangleSpanners(false, unitigGraph, readPaths, nodeTipBelongsToTangle, validSpans, chainPositionsInReads, anchorChains);
	return std::make_pair(std::move(resultGraph), std::move(resultPaths));
}

bool hasAcyclicConnection(const SparseEdgeContainer& edges, const uint64_t startNode, const uint64_t endNode)
{
	phmap::flat_hash_set<uint64_t> reachableFromFw;
	phmap::flat_hash_set<uint64_t> reachableFromBw;
	std::vector<uint64_t> checkStack;
	checkStack.emplace_back(startNode);
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (reachableFromFw.count(top) == 1) continue;
		reachableFromFw.insert(top);
		if (top == endNode) continue;
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			checkStack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
		}
	}
	checkStack.emplace_back(endNode ^ firstBitUint64_t);
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (reachableFromBw.count(top ^ firstBitUint64_t) == 1) continue;
		reachableFromBw.insert(top ^ firstBitUint64_t);
		if (top == (startNode ^ firstBitUint64_t)) continue;
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			checkStack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
		}
	}
	phmap::flat_hash_set<uint64_t> seenNodes;
	for (auto edge : edges.getEdges(std::make_pair(startNode & maskUint64_t, startNode & firstBitUint64_t)))
	{
		checkStack.emplace_back(edge.first + (edge.second ? firstBitUint64_t : 0));
	}
	seenNodes.emplace(startNode);
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (seenNodes.count(top) == 1) continue;
		bool hasNonvisitedPredecessor = false;
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, (top ^ firstBitUint64_t) & firstBitUint64_t)))
		{
			uint64_t otherNode = edge.first + (edge.second ? 0 : firstBitUint64_t);
			if (reachableFromFw.count(otherNode) == 0) continue;
			if (reachableFromBw.count(otherNode) == 0) continue;
			if (seenNodes.count(otherNode) == 0)
			{
				hasNonvisitedPredecessor = true;
				break;
			}
		}
		if (hasNonvisitedPredecessor) continue;
		if (top == endNode) return true;
		seenNodes.insert(top);
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			checkStack.emplace_back(edge.first + (edge.second ? firstBitUint64_t : 0 ));
		}
	}
	return false;
}

bool splitChainsAtCyclicTangles(const UnitigGraph& unitigGraph, std::vector<AnchorChain>& chains)
{
	bool result = false;
	auto edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	phmap::flat_hash_set<std::pair<size_t, size_t>> splitBeforeTheseNodes;
	for (size_t i = 0; i < chains.size(); i++)
	{
		if (chains[i].ploidy != 1) continue;
		for (size_t j = 1; j < chains[i].nodes.size(); j++)
		{
			if (hasAcyclicConnection(edges, chains[i].nodes[j-1], chains[i].nodes[j])) continue;
			splitBeforeTheseNodes.emplace(i, j);
		}
	}
	for (size_t i = 0; i < chains.size(); i++)
	{
		if (chains[i].ploidy != 1) continue;
		for (size_t j = chains[i].nodes.size()-1; j > 0; j--)
		{
			if (splitBeforeTheseNodes.count(std::make_pair(i, j)) == 0) continue;
			std::cerr << "split chain " << i << " at index " << j << " / " << chains[i].nodes.size() << std::endl;
			result = true;
			chains.emplace_back();
			chains.back().ploidy = 1;
			chains.back().nodes.insert(chains.back().nodes.end(), chains[i].nodes.begin()+j, chains[i].nodes.end());
			chains.back().nodeOffsets.insert(chains.back().nodeOffsets.end(), chains[i].nodeOffsets.begin()+j, chains[i].nodeOffsets.end());
			chains[i].nodes.erase(chains[i].nodes.begin()+j, chains[i].nodes.end());
			chains[i].nodeOffsets.erase(chains[i].nodeOffsets.begin()+j, chains[i].nodeOffsets.end());
		}
	}
	for (size_t i = chains.size()-1; i < chains.size(); i--)
	{
		if (chains[i].nodeOffsets[0] != 0)
		{
			for (size_t j = chains[i].nodeOffsets.size()-1; j < chains[i].nodeOffsets.size(); j--)
			{
				chains[i].nodeOffsets[j] -= chains[i].nodeOffsets[0];
			}
		}
		size_t chainLength = 0;
		for (uint64_t node : chains[i].nodes)
		{
			chainLength += unitigGraph.lengths[node & maskUint64_t];
		}
		if (chainLength < 500)
		{
			result = true;
			std::swap(chains[i], chains.back());
			chains.pop_back();
		}
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolveUniqueChainContainedTangles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage)
{
	std::cerr << "try gap fill" << std::endl;
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	fixFakeSinglePloidyChains(anchorChains, chainPositionsInReads, approxOneHapCoverage);
	bool splitAnything = splitChainsAtCyclicTangles(unitigGraph, anchorChains);
	if (!splitAnything) return std::make_pair(unitigGraph, readPaths);
	chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	VectorWithDirection<size_t> nodeTipBelongsToTangle = getNodeTipTangleAssignmentsUniqueChains(unitigGraph, anchorChains);
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> validSpans = getTangleSpanners(false, unitigGraph, readPaths, nodeTipBelongsToTangle, approxOneHapCoverage * 0.5, 3, anchorChains, chainPositionsInReads);
	UnitigGraph resultGraph;
	std::vector<ReadPathBundle> resultPaths;
	std::tie(resultGraph, resultPaths) = applyTangleSpanners(false, unitigGraph, readPaths, nodeTipBelongsToTangle, validSpans, chainPositionsInReads, anchorChains);
	return std::make_pair(std::move(resultGraph), std::move(resultPaths));
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolveSpannedTangles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage)
{
	auto result = resolveSpannedTangles(false, unitigGraph, readPaths, graphk, approxOneHapCoverage);
	result = resolveSpannedTangles(true, result.first, result.second, graphk, approxOneHapCoverage);
	result = resolveSpannedTanglesBipartite(false, result.first, result.second, graphk, approxOneHapCoverage);
	result = resolveSpannedTanglesBipartite(true, result.first, result.second, graphk, approxOneHapCoverage);
	result = resolveUniqueChainContainedTangles(result.first, result.second, graphk, approxOneHapCoverage);
	return result;
}

void tryPopBubble(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, RankBitvector& kept, const SparseEdgeContainer& edges, uint64_t startNode, uint64_t endNode)
{
	if (!hasAcyclicConnection(edges, startNode, endNode)) return;
	std::vector<uint64_t> checkStack;
	phmap::flat_hash_set<uint64_t> seenFw;
	phmap::flat_hash_set<uint64_t> seenBw;
	checkStack.emplace_back(startNode);
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (seenFw.count(top) == 1) continue;
		seenFw.insert(top);
		if (top == endNode) continue;
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			checkStack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
		}
	}
	checkStack.emplace_back(endNode ^ firstBitUint64_t);
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (seenBw.count(top ^ firstBitUint64_t) == 1) continue;
		seenBw.insert(top ^ firstBitUint64_t);
		if (top == (startNode ^ firstBitUint64_t)) continue;
		for (auto edge : edges.getEdges(std::make_pair(top & maskUint64_t, top & firstBitUint64_t)))
		{
			checkStack.push_back(edge.first + (edge.second ? firstBitUint64_t : 0));
		}
	}
	for (auto node : seenFw)
	{
		kept.set(node & maskUint64_t, false);
	}
	for (auto node : seenBw)
	{
		kept.set(node & maskUint64_t, false);
	}
	kept.set(startNode & maskUint64_t, true);
	kept.set(endNode & maskUint64_t, true);
	uint64_t pos = startNode;
	// todo fix: better algorithm
	while (pos != endNode)
	{
		kept.set(pos & maskUint64_t, true);
		uint64_t maxEdge = pos;
		double maxEdgeCoverage = -1;
		for (auto edge : edges.getEdges(std::make_pair(pos & maskUint64_t, pos & firstBitUint64_t)))
		{
			if (seenFw.count(edge.first + (edge.second ? firstBitUint64_t : 0)) == 0) continue;
			if (seenBw.count(edge.first + (edge.second ? firstBitUint64_t : 0)) == 0) continue;
			double coverage = unitigGraph.coverages[edge.first];
			if (coverage > maxEdgeCoverage)
			{
				maxEdge = edge.first + (edge.second ? firstBitUint64_t : 0);
				maxEdgeCoverage = coverage;
			}
		}
		assert(maxEdge != pos);
		pos = maxEdge;
	}
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> popHaploidBubbles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<AnchorChain>& anchorChains)
{
	RankBitvector keptNodes;
	keptNodes.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < keptNodes.size(); i++)
	{
		keptNodes.set(i, true);
	}
	auto edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		if (anchorChains[i].ploidy != 1) continue;
		if (anchorChains[i].nodes.size() < 2) continue;
		for (size_t j = 1; j < anchorChains[i].nodes.size(); j++)
		{
			tryPopBubble(unitigGraph, readPaths, keptNodes, edges, anchorChains[i].nodes[j-1], anchorChains[i].nodes[j]);
		}
	}
	return filterUnitigGraph(unitigGraph, readPaths, keptNodes);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> popHaploidChainBubbles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage)
{
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, 100+graphk, approxOneHapCoverage);
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	fixFakeSinglePloidyChains(anchorChains, chainPositionsInReads, approxOneHapCoverage);
	UnitigGraph resultGraph;
	std::vector<ReadPathBundle> resultPaths;
	std::tie(resultGraph, resultPaths) = popHaploidBubbles(unitigGraph, readPaths, anchorChains);
	return std::make_pair(std::move(resultGraph), std::move(resultPaths));
}

bool isGoodPhasableFork(const UnitigGraph& unitigGraph, const SparseEdgeContainer& edges, const std::pair<size_t, bool> fork, const double approxOneHapCoverage)
{
	if (edges.getEdges(fork).size() < 2) return false;
	if (unitigGraph.coverages[fork.first] < approxOneHapCoverage * 1.5) return false;
	size_t estimatedMainNodeCoverage = (double)unitigGraph.coverages[fork.first] / approxOneHapCoverage + 0.5;
	size_t totalEstimatedCoverage = 0;
	for (auto edge : edges.getEdges(fork))
	{
		if (edges.getEdges(reverse(edge)).size() != 1) return false;
		if (unitigGraph.coverages[edge.first] < approxOneHapCoverage * 0.5) return false;
		totalEstimatedCoverage += (double)unitigGraph.coverages[edge.first] / approxOneHapCoverage + 0.5;
	}
	if (totalEstimatedCoverage < 2) return false;
	if (totalEstimatedCoverage < estimatedMainNodeCoverage - 1) return false;
	if (totalEstimatedCoverage > estimatedMainNodeCoverage + 1) return false;
	return true;
}

std::tuple<std::vector<bool>, std::vector<bool>, phmap::flat_hash_map<uint64_t, size_t>> getGoodForks(const UnitigGraph& unitigGraph, const double approxOneHapCoverage)
{
	SparseEdgeContainer edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	std::vector<bool> hasGoodForkFw;
	std::vector<bool> hasGoodForkBw;
	phmap::flat_hash_map<uint64_t, size_t> forkAlleleCount;
	hasGoodForkFw.resize(unitigGraph.nodeCount(), false);
	hasGoodForkBw.resize(unitigGraph.nodeCount(), false);
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (isGoodPhasableFork(unitigGraph, edges, std::make_pair(i, true), approxOneHapCoverage))
		{
			hasGoodForkFw[i] = true;
			forkAlleleCount[i + firstBitUint64_t] = edges.getEdges(std::make_pair(i, true)).size();
		}
		if (isGoodPhasableFork(unitigGraph, edges, std::make_pair(i, false), approxOneHapCoverage))
		{
			hasGoodForkBw[i] = true;
			forkAlleleCount[i] = edges.getEdges(std::make_pair(i, false)).size();
		}
	}
	return std::make_tuple(std::move(hasGoodForkBw), std::move(hasGoodForkFw), std::move(forkAlleleCount));
}

std::pair<std::vector<std::vector<HaplotypeInformativeSite>>, std::vector<std::vector<HaplotypeInformativeSite>>> getReadForkEdges(const std::vector<bool>& hasGoodForkBw, const std::vector<bool>& hasGoodForkFw, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths)
{
	std::vector<std::vector<HaplotypeInformativeSite>> readForkEdgesFw;
	std::vector<std::vector<HaplotypeInformativeSite>> readForkEdgesBw;
	readForkEdgesBw.resize(readPaths.size());
	readForkEdgesFw.resize(readPaths.size());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (k == 0) continue;
				uint64_t prevNode = readPaths[i].paths[j].path[k-1];
				if (((prevNode & firstBitUint64_t) && hasGoodForkFw[prevNode & maskUint64_t]) || (((prevNode ^ firstBitUint64_t) & firstBitUint64_t) && hasGoodForkBw[prevNode & maskUint64_t]))
				{
					readForkEdgesFw[i].emplace_back();
					readForkEdgesFw[i].back().readPos = readPos - 1 - unitigGraph.lengths[node & maskUint64_t];
					readForkEdgesFw[i].back().bubble = prevNode;
					readForkEdgesFw[i].back().allele = node;
				}
				if (((node & firstBitUint64_t) && hasGoodForkBw[node & maskUint64_t]) || (((node ^ firstBitUint64_t) & firstBitUint64_t) && hasGoodForkFw[node & maskUint64_t]))
				{
					readForkEdgesBw[i].emplace_back();
					readForkEdgesBw[i].back().readPos = readPos - unitigGraph.lengths[node & maskUint64_t];
					readForkEdgesBw[i].back().bubble = node ^ firstBitUint64_t;
					readForkEdgesBw[i].back().allele = prevNode ^ firstBitUint64_t;
				}
			}
		}
	}
	return std::make_pair(std::move(readForkEdgesBw), std::move(readForkEdgesFw));
}

bool isUsableCluster(const phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t>& allelesInCluster, const double approxOneHapCoverage, const size_t firstMainCoverage, const size_t firstAlleleCount, const size_t secondMainCoverage, const size_t secondAlleleCount, const phmap::flat_hash_map<uint64_t, uint64_t>& forkAlleleCoverage)
{
	assert(firstAlleleCount >= 2);
	assert(firstAlleleCount == secondAlleleCount);
	phmap::flat_hash_set<uint64_t> allelesInFirst;
	phmap::flat_hash_set<uint64_t> allelesInSecond;
	size_t validPairs = 0;
	size_t validPairsCoverage = 0;
	size_t falsePositivePairsCoverage = 0;
	for (auto pair : allelesInCluster)
	{
		if (pair.second >= approxOneHapCoverage * 0.5)
		{
			validPairs += 1;
			validPairsCoverage += pair.second;
			allelesInFirst.emplace(pair.first.first);
			allelesInSecond.emplace(pair.first.second);
		}
		else if (pair.second >= approxOneHapCoverage * 0.25)
		{
			// ambiguous, just don't handle and assume not good
			return false;
		}
		else
		{
			falsePositivePairsCoverage += pair.second;
		}
	}
	if (allelesInFirst.size() < firstAlleleCount) return false;
	if (allelesInSecond.size() < secondAlleleCount) return false;
	if (validPairs != firstAlleleCount) return false;
	if (validPairsCoverage < approxOneHapCoverage * 1.5) return false;
	if (falsePositivePairsCoverage >= approxOneHapCoverage * 0.5) return false;
	if (falsePositivePairsCoverage > validPairsCoverage * 0.1) return false;
	if (validPairsCoverage < std::min(firstMainCoverage, secondMainCoverage)*0.25) return false;
	return true;
}

bool forkIsHaplotypeInformative(const uint64_t left, const uint64_t right, const phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<size_t, std::vector<uint64_t>>>& forkPositionsInReads, const std::vector<std::vector<HaplotypeInformativeSite>>& forkEdgesBw, const std::vector<std::vector<HaplotypeInformativeSite>>& forkEdgesFw, const double approxOneHapCoverage, const phmap::flat_hash_map<uint64_t, size_t>& forkMainCoverage, const phmap::flat_hash_map<uint64_t, size_t>& forkAlleleCoverage, const phmap::flat_hash_map<uint64_t, size_t>& forkAlleleCount)
{
	const size_t clusteringMaxDistance = 100;
	phmap::flat_hash_map<std::tuple<size_t, uint64_t, uint64_t>, size_t> allelePairsHereMap;
	for (const auto& pair : forkPositionsInReads.at(left))
	{
		size_t read = pair.first;
		if (forkPositionsInReads.at(right).count(read) == 0) continue;
		for (const uint64_t secondpos : forkPositionsInReads.at(right).at(read))
		{
			for (const uint64_t firstpos : pair.second)
			{
				if ((firstpos & firstBitUint64_t) == (secondpos & firstBitUint64_t)) continue;
				if (firstpos & firstBitUint64_t)
				{
					assert((secondpos ^ firstBitUint64_t) & firstBitUint64_t);
					size_t i = firstpos & maskUint64_t;
					size_t j = secondpos & maskUint64_t;
					assert(forkEdgesFw[read][i].bubble == left);
					assert(forkEdgesBw[read][j].bubble == right);
					if (forkEdgesFw[read][i].readPos < forkEdgesBw[read][j].readPos) continue;
					allelePairsHereMap[std::make_tuple(forkEdgesFw[read][i].readPos - forkEdgesBw[read][j].readPos, forkEdgesFw[read][i].allele, forkEdgesBw[read][j].allele)] += 1;
				}
				else
				{
					assert(secondpos & firstBitUint64_t);
					size_t i = firstpos & maskUint64_t;
					size_t j = secondpos & maskUint64_t;
					assert(forkEdgesBw[read][i].bubble == left);
					assert(forkEdgesFw[read][j].bubble == right);
					if (forkEdgesFw[read][j].readPos < forkEdgesBw[read][i].readPos) continue;
					allelePairsHereMap[std::make_tuple(forkEdgesFw[read][j].readPos - forkEdgesBw[read][i].readPos, forkEdgesBw[read][i].allele, forkEdgesFw[read][j].allele)] += 1;
				}
			}
		}
	}
	std::vector<std::tuple<size_t, uint64_t, uint64_t, size_t>> allelePairsHere;
	for (auto t : allelePairsHereMap)
	{
		allelePairsHere.emplace_back(std::get<0>(t.first), std::get<1>(t.first), std::get<2>(t.first), t.second);
	}
	std::sort(allelePairsHere.begin(), allelePairsHere.end());
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> allelesInCluster;
	allelesInCluster[std::make_pair(std::get<1>(allelePairsHere[0]), std::get<2>(allelePairsHere[0]))] += std::get<3>(allelePairsHere[0]);
	for (size_t i = 1; i < allelePairsHere.size(); i++)
	{
		if (std::get<0>(allelePairsHere[i]) > std::get<0>(allelePairsHere[i-1])+clusteringMaxDistance)
		{
			if (isUsableCluster(allelesInCluster, approxOneHapCoverage, forkMainCoverage.at(left), forkAlleleCount.at(left), forkMainCoverage.at(right), forkAlleleCount.at(right), forkAlleleCoverage))
			{
				return true;
			}
			allelesInCluster.clear();
		}
		allelesInCluster[std::make_pair(std::get<1>(allelePairsHere[i]), std::get<2>(allelePairsHere[i]))] += std::get<3>(allelePairsHere[i]);
	}
	if (isUsableCluster(allelesInCluster, approxOneHapCoverage, forkMainCoverage.at(left), forkAlleleCount.at(left), forkMainCoverage.at(right), forkAlleleCount.at(right), forkAlleleCoverage))
	{
		return true;
	}
	return false;
}

std::pair<std::vector<bool>, std::vector<bool>> getHaplotypeInformativeForkPairs(const std::vector<std::vector<HaplotypeInformativeSite>>& forkEdgesBw, const std::vector<std::vector<HaplotypeInformativeSite>>& forkEdgesFw, const size_t nodeCount, const phmap::flat_hash_map<uint64_t, size_t>& forkAlleleCount, const double approxOneHapCoverage)
{
	assert(forkEdgesFw.size() == forkEdgesBw.size());
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, phmap::flat_hash_map<std::tuple<size_t, uint64_t, uint64_t>, size_t>> forkCrossers;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<size_t, std::vector<uint64_t>>> forkPositionsInReads;
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> checkablePairsCoverage;
	phmap::flat_hash_map<uint64_t, size_t> forkMainCoverage;
	phmap::flat_hash_map<uint64_t, size_t> forkAlleleCoverage;
	for (size_t i = 0; i < forkEdgesFw.size(); i++)
	{
		for (size_t j = 0; j < forkEdgesFw[i].size(); j++)
		{
			auto site = forkEdgesFw[i][j];
			forkPositionsInReads[site.bubble][i].emplace_back(j + firstBitUint64_t);
			forkMainCoverage[site.bubble] += 1;
			forkAlleleCoverage[site.allele] += 1;
		}
		for (size_t j = 0; j < forkEdgesBw[i].size(); j++)
		{
			auto site = forkEdgesBw[i][j];
			forkPositionsInReads[site.bubble][i].emplace_back(j);
			forkMainCoverage[site.bubble] += 1;
			forkAlleleCoverage[site.allele] += 1;
		}
	}
	std::cerr << "find checkable pairs" << std::endl;
	std::vector<bool> resultFw;
	std::vector<bool> resultBw;
	resultFw.resize(nodeCount, false);
	resultBw.resize(nodeCount, false);
	phmap::flat_hash_set<std::pair<uint64_t, uint64_t>> checkedPairs;
	for (size_t i = 0; i < forkEdgesBw.size(); i++)
	{
		std::cerr << "read " << i << " has " << forkEdgesBw[i].size() << " bw forks, " << forkEdgesFw[i].size() << " fw forks" << std::endl;
		std::vector<std::vector<HaplotypeInformativeSite>> fwSitesPerAlleleCount;
		std::vector<std::vector<HaplotypeInformativeSite>> fwSitesNeedsCheckingPerAlleleCount;
		std::vector<std::vector<HaplotypeInformativeSite>> bwSitesResolvedPerAlleleCount;
		std::vector<std::vector<HaplotypeInformativeSite>> bwSitesNeedsCheckingPerAlleleCount;
		fwSitesPerAlleleCount.resize(5);
		fwSitesNeedsCheckingPerAlleleCount.resize(5);
		bwSitesResolvedPerAlleleCount.resize(5);
		bwSitesNeedsCheckingPerAlleleCount.resize(5);
		for (auto fwsite : forkEdgesFw[i])
		{
			assert(forkAlleleCount.at(fwsite.bubble) >= 2);
			size_t alleleCount = forkAlleleCount.at(fwsite.bubble)-2;
			if (alleleCount >= fwSitesPerAlleleCount.size()) alleleCount = fwSitesPerAlleleCount.size()-1;
			assert(fwSitesPerAlleleCount[alleleCount].size() == 0 || fwsite.readPos > fwSitesPerAlleleCount[alleleCount].back().readPos);
			fwSitesPerAlleleCount[alleleCount].emplace_back(fwsite);
			if ((fwsite.bubble & firstBitUint64_t) && !resultFw[fwsite.bubble & maskUint64_t]) fwSitesNeedsCheckingPerAlleleCount[alleleCount].emplace_back(fwsite);
			if (((fwsite.bubble ^ firstBitUint64_t) & firstBitUint64_t) && !resultBw[fwsite.bubble & maskUint64_t]) fwSitesNeedsCheckingPerAlleleCount[alleleCount].emplace_back(fwsite);
		}
		for (auto bwsite : forkEdgesBw[i])
		{
			assert(forkAlleleCount.at(bwsite.bubble) >= 2);
			size_t alleleCount = forkAlleleCount.at(bwsite.bubble)-2;
			if (alleleCount >= bwSitesResolvedPerAlleleCount.size()) alleleCount = bwSitesResolvedPerAlleleCount.size()-1;
			if (((bwsite.bubble & firstBitUint64_t) && !resultFw[bwsite.bubble & maskUint64_t]) || (((bwsite.bubble ^ firstBitUint64_t) & firstBitUint64_t) && !resultBw[bwsite.bubble & maskUint64_t]))
			{
				bwSitesNeedsCheckingPerAlleleCount[alleleCount].emplace_back(bwsite);
			}
			else
			{
				bwSitesResolvedPerAlleleCount[alleleCount].emplace_back(bwsite);
			}
		}
		for (size_t count = 0; count < fwSitesPerAlleleCount.size(); count++)
		{
			for (size_t bwi = 0; bwi < bwSitesNeedsCheckingPerAlleleCount[count].size(); bwi++)
			{
				auto bwsite = bwSitesNeedsCheckingPerAlleleCount[count][bwi];
				if ((bwsite.bubble & firstBitUint64_t) && resultFw[bwsite.bubble & maskUint64_t]) break;
				if (((bwsite.bubble ^ firstBitUint64_t) & firstBitUint64_t) && resultBw[bwsite.bubble & maskUint64_t]) break;
				size_t fwistart = 0;
				for (size_t fwi = fwistart; fwi < fwSitesPerAlleleCount[count].size(); fwi++)
				{
					auto fwsite = fwSitesPerAlleleCount[count][fwi];
					if (fwsite.readPos < bwsite.readPos)
					{
						fwistart = fwi+1;
						continue;
					}
					if (forkAlleleCount.at(fwsite.bubble) != forkAlleleCount.at(bwsite.bubble)) continue;
					std::pair<uint64_t, uint64_t> key;
					// size_t distance = fwsite.readPos - bwsite.readPos;
					if (fwsite.bubble < bwsite.bubble)
					{
						key = std::make_pair(fwsite.bubble, bwsite.bubble);
						// checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)] += 1;
						// if (checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)] >= approxOneHapCoverage * 1.5 && checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)]-1 < approxOneHapCoverage * 1.5)
						// {
						// 	pairsWithEnoughCoverage += 1;
						// }
						// forkCrossers[std::make_pair(fwsite.bubble, bwsite.bubble)][std::make_tuple(distance, fwsite.allele, bwsite.allele)] += 1;
					}
					else
					{
						key = std::make_pair(bwsite.bubble, fwsite.bubble);
						// checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)] += 1;
						// if (checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)] >= approxOneHapCoverage * 1.5 && checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)]-1 < approxOneHapCoverage * 1.5)
						// {
						// 	pairsWithEnoughCoverage += 1;
						// }
						// forkCrossers[std::make_pair(bwsite.bubble, fwsite.bubble)][std::make_tuple(distance, bwsite.allele, fwsite.allele)] += 1;
					}
					if (checkedPairs.count(key) == 1) continue;
					checkedPairs.insert(key);
					std::cerr << "check " << ((key.first & firstBitUint64_t) ? ">" : "<") << (key.first & maskUint64_t) << " " << ((key.second & firstBitUint64_t) ? ">" : "<") << (key.second & maskUint64_t) << std::endl;
					if (forkIsHaplotypeInformative(key.first, key.second, forkPositionsInReads, forkEdgesBw, forkEdgesFw, approxOneHapCoverage, forkMainCoverage, forkAlleleCoverage, forkAlleleCount))
					{
						std::cerr << "haplotype informative " << ((key.first & firstBitUint64_t) ? ">" : "<") << (key.first & maskUint64_t) << " " << ((key.second & firstBitUint64_t) ? ">" : "<") << (key.second & maskUint64_t) << std::endl;
						if (key.first & firstBitUint64_t)
						{
							resultFw[key.first & maskUint64_t] = true;
						}
						else
						{
							resultBw[key.first & maskUint64_t] = true;
						}
						if (key.second & firstBitUint64_t)
						{
							resultFw[key.second & maskUint64_t] = true;
						}
						else
						{
							resultBw[key.second & maskUint64_t] = true;
						}
					}
				}
			}
			for (size_t bwi = 0; bwi < bwSitesResolvedPerAlleleCount[count].size(); bwi++)
			{
				auto bwsite = bwSitesResolvedPerAlleleCount[count][bwi];
				if ((bwsite.bubble & firstBitUint64_t) && resultFw[bwsite.bubble & maskUint64_t]) break;
				if (((bwsite.bubble ^ firstBitUint64_t) & firstBitUint64_t) && resultBw[bwsite.bubble & maskUint64_t]) break;
				size_t fwistart = 0;
				for (size_t fwi = fwistart; fwi < fwSitesNeedsCheckingPerAlleleCount[count].size(); fwi++)
				{
					auto fwsite = fwSitesNeedsCheckingPerAlleleCount[count][fwi];
					if (fwsite.readPos < bwsite.readPos)
					{
						fwistart = fwi+1;
						continue;
					}
					if (forkAlleleCount.at(fwsite.bubble) != forkAlleleCount.at(bwsite.bubble)) continue;
					std::pair<uint64_t, uint64_t> key;
					// size_t distance = fwsite.readPos - bwsite.readPos;
					if (fwsite.bubble < bwsite.bubble)
					{
						key = std::make_pair(fwsite.bubble, bwsite.bubble);
						// checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)] += 1;
						// if (checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)] >= approxOneHapCoverage * 1.5 && checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)]-1 < approxOneHapCoverage * 1.5)
						// {
						// 	pairsWithEnoughCoverage += 1;
						// }
						// forkCrossers[std::make_pair(fwsite.bubble, bwsite.bubble)][std::make_tuple(distance, fwsite.allele, bwsite.allele)] += 1;
					}
					else
					{
						key = std::make_pair(bwsite.bubble, fwsite.bubble);
						// checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)] += 1;
						// if (checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)] >= approxOneHapCoverage * 1.5 && checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)]-1 < approxOneHapCoverage * 1.5)
						// {
						// 	pairsWithEnoughCoverage += 1;
						// }
						// forkCrossers[std::make_pair(bwsite.bubble, fwsite.bubble)][std::make_tuple(distance, bwsite.allele, fwsite.allele)] += 1;
					}
					if (checkedPairs.count(key) == 1) continue;
					checkedPairs.insert(key);
					std::cerr << "check " << ((key.first & firstBitUint64_t) ? ">" : "<") << (key.first & maskUint64_t) << " " << ((key.second & firstBitUint64_t) ? ">" : "<") << (key.second & maskUint64_t) << std::endl;
					if (forkIsHaplotypeInformative(key.first, key.second, forkPositionsInReads, forkEdgesBw, forkEdgesFw, approxOneHapCoverage, forkMainCoverage, forkAlleleCoverage, forkAlleleCount))
					{
						std::cerr << "haplotype informative " << ((key.first & firstBitUint64_t) ? ">" : "<") << (key.first & maskUint64_t) << " " << ((key.second & firstBitUint64_t) ? ">" : "<") << (key.second & maskUint64_t) << std::endl;
						if (key.first & firstBitUint64_t)
						{
							resultFw[key.first & maskUint64_t] = true;
						}
						else
						{
							resultBw[key.first & maskUint64_t] = true;
						}
						if (key.second & firstBitUint64_t)
						{
							resultFw[key.second & maskUint64_t] = true;
						}
						else
						{
							resultBw[key.second & maskUint64_t] = true;
						}
					}
				}
			}
		}
		// for (auto bwsite : forkEdgesBw[i])
		// {
		// 	for (auto fwsite : forkEdgesFw[i])
		// 	{
		// 		if (fwsite.readPos < bwsite.readPos) continue;
		// 		if (forkAlleleCount.at(fwsite.bubble) != forkAlleleCount.at(bwsite.bubble)) continue;
		// 		// size_t distance = fwsite.readPos - bwsite.readPos;
		// 		if (fwsite.bubble < bwsite.bubble || (fwsite.bubble == bwsite.bubble && fwsite.allele < bwsite.allele))
		// 		{
		// 			checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)] += 1;
		// 			if (checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)] >= approxOneHapCoverage * 1.5 && checkablePairsCoverage[std::make_pair(fwsite.bubble, bwsite.bubble)]-1 < approxOneHapCoverage * 1.5)
		// 			{
		// 				pairsWithEnoughCoverage += 1;
		// 			}
		// 			// forkCrossers[std::make_pair(fwsite.bubble, bwsite.bubble)][std::make_tuple(distance, fwsite.allele, bwsite.allele)] += 1;
		// 		}
		// 		else
		// 		{
		// 			checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)] += 1;
		// 			if (checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)] >= approxOneHapCoverage * 1.5 && checkablePairsCoverage[std::make_pair(bwsite.bubble, fwsite.bubble)]-1 < approxOneHapCoverage * 1.5)
		// 			{
		// 				pairsWithEnoughCoverage += 1;
		// 			}
		// 			// forkCrossers[std::make_pair(bwsite.bubble, fwsite.bubble)][std::make_tuple(distance, bwsite.allele, fwsite.allele)] += 1;
		// 		}
		// 	}
		// }
	}
	// std::cerr << checkablePairsCoverage.size() << " total pairs" << std::endl;
	// std::cerr << pairsWithEnoughCoverage << " pairs with enough coverage" << std::endl;
	// std::cerr << "check if pairs might be haplotype informative" << std::endl;
	// for (auto covpair : checkablePairsCoverage)
	// {
	// 	if (covpair.second < approxOneHapCoverage * 1.5) continue;
	// 	auto pair = covpair.first;
	// 	if (((pair.first & firstBitUint64_t) && resultFw[pair.first & maskUint64_t]) || (((pair.first ^ firstBitUint64_t) & firstBitUint64_t) && resultBw[pair.first & maskUint64_t]))
	// 	{
	// 		if (((pair.second & firstBitUint64_t) && resultFw[pair.second & maskUint64_t]) || (((pair.second ^ firstBitUint64_t) & firstBitUint64_t) && resultBw[pair.second & maskUint64_t]))
	// 		{
	// 			continue;
	// 		}
	// 	}
	// 	phmap::flat_hash_map<size_t, std::vector<uint64_t>> firstPositionsInReads;
	// 	phmap::flat_hash_map<size_t, std::vector<uint64_t>> secondPositionsInReads;
	// 	for (const auto& pair2 : forkPositionsInReads.at(pair.first))
	// 	{
	// 		firstPositionsInReads[pair2.first].emplace_back(pair2.second);
	// 	}
	// 	for (const auto& pair2 : forkPositionsInReads.at(pair.second))
	// 	{
	// 		secondPositionsInReads[pair2.first].emplace_back(pair2.second);
	// 	}
	// 	phmap::flat_hash_map<std::tuple<size_t, uint64_t, uint64_t>, size_t> allelePairsHereMap;
	// 	for (const auto& pair2 : firstPositionsInReads)
	// 	{
	// 		size_t read = pair2.first;
	// 		if (secondPositionsInReads.count(read) == 0) continue;
	// 		for (const uint64_t firstpos : pair2.second)
	// 		{
	// 			for (const uint64_t secondpos : secondPositionsInReads.at(read))
	// 			{
	// 				if ((firstpos & firstBitUint64_t) == (secondpos & firstBitUint64_t)) continue;
	// 				if (firstpos & firstBitUint64_t)
	// 				{
	// 					assert((secondpos ^ firstBitUint64_t) & firstBitUint64_t);
	// 					size_t i = firstpos & maskUint64_t;
	// 					size_t j = secondpos & maskUint64_t;
	// 					if (forkEdgesFw[read][i].readPos < forkEdgesBw[read][j].readPos) continue;
	// 					allelePairsHereMap[std::make_tuple(forkEdgesFw[read][i].readPos - forkEdgesBw[read][j].readPos, forkEdgesFw[read][i].allele, forkEdgesBw[read][j].allele)] += 1;
	// 				}
	// 				else
	// 				{
	// 					assert(secondpos & firstBitUint64_t);
	// 					size_t j = firstpos & maskUint64_t;
	// 					size_t i = secondpos & maskUint64_t;
	// 					if (forkEdgesFw[read][i].readPos < forkEdgesBw[read][j].readPos) continue;
	// 					allelePairsHereMap[std::make_tuple(forkEdgesFw[read][i].readPos - forkEdgesBw[read][j].readPos, forkEdgesFw[read][i].allele, forkEdgesBw[read][j].allele)] += 1;
	// 				}
	// 			}
	// 		}
	// 	}
	// 	std::vector<std::tuple<size_t, uint64_t, uint64_t, size_t>> allelePairsHere;
	// 	for (auto t : allelePairsHereMap)
	// 	{
	// 		allelePairsHere.emplace_back(std::get<0>(t.first), std::get<1>(t.first), std::get<2>(t.first), t.second);
	// 	}
	// 	std::sort(allelePairsHere.begin(), allelePairsHere.end());
	// 	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, size_t> allelesInCluster;
	// 	allelesInCluster[std::make_pair(std::get<1>(allelePairsHere[0]), std::get<2>(allelePairsHere[0]))] += std::get<3>(allelePairsHere[0]);
	// 	for (size_t i = 1; i < allelePairsHere.size(); i++)
	// 	{
	// 		if (std::get<0>(allelePairsHere[i]) > std::get<0>(allelePairsHere[i-1])+clusteringMaxDistance)
	// 		{
	// 			if (isUsableCluster(allelesInCluster, approxOneHapCoverage, forkMainCoverage[pair.first], forkAlleleCount.at(pair.first), forkMainCoverage[pair.second], forkAlleleCount.at(pair.second), forkAlleleCoverage))
	// 			{
	// 				if (pair.first & firstBitUint64_t)
	// 				{
	// 					resultFw[pair.first & maskUint64_t] = true;
	// 				}
	// 				else
	// 				{
	// 					resultBw[pair.first & maskUint64_t] = true;
	// 				}
	// 				if (pair.second & firstBitUint64_t)
	// 				{
	// 					resultFw[pair.second & maskUint64_t] = true;
	// 				}
	// 				else
	// 				{
	// 					resultBw[pair.second & maskUint64_t] = true;
	// 				}
	// 				break;
	// 			}
	// 			allelesInCluster.clear();
	// 		}
	// 		allelesInCluster[std::make_pair(std::get<1>(allelePairsHere[i]), std::get<2>(allelePairsHere[i]))] += std::get<3>(allelePairsHere[i]);
	// 	}
	// 	if (isUsableCluster(allelesInCluster, approxOneHapCoverage, forkMainCoverage[pair.first], forkAlleleCount.at(pair.first), forkMainCoverage[pair.second], forkAlleleCount.at(pair.second), forkAlleleCoverage))
	// 	{
	// 		if (pair.first & firstBitUint64_t)
	// 		{
	// 			resultFw[pair.first & maskUint64_t] = true;
	// 		}
	// 		else
	// 		{
	// 			resultBw[pair.first & maskUint64_t] = true;
	// 		}
	// 		if (pair.second & firstBitUint64_t)
	// 		{
	// 			resultFw[pair.second & maskUint64_t] = true;
	// 		}
	// 		else
	// 		{
	// 			resultBw[pair.second & maskUint64_t] = true;
	// 		}
	// 	}
	// }
	return std::make_pair(std::move(resultBw), std::move(resultFw));
}

std::vector<std::vector<HaplotypeInformativeSite>> getHaplotypeInformativeForks(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	std::vector<bool> hasGoodForkFw;
	std::vector<bool> hasGoodForkBw;
	phmap::flat_hash_map<uint64_t, size_t> forkAlleleCount;
	std::tie(hasGoodForkBw, hasGoodForkFw, forkAlleleCount) = getGoodForks(unitigGraph, approxOneHapCoverage);
	for (size_t i = 0; i < hasGoodForkFw.size(); i++)
	{
		if (hasGoodForkFw[i]) std::cerr << "maybe good fork >" << i << std::endl;
		if (hasGoodForkBw[i]) std::cerr << "maybe good fork <" << i << std::endl;
	}
	std::vector<std::vector<HaplotypeInformativeSite>> readForkEdgesFw;
	std::vector<std::vector<HaplotypeInformativeSite>> readForkEdgesBw;
	std::tie(readForkEdgesBw, readForkEdgesFw) = getReadForkEdges(hasGoodForkBw, hasGoodForkFw, unitigGraph, readPaths);
	std::vector<bool> hasHaplotypeInformativeForkFw;
	std::vector<bool> hasHaplotypeInformativeForkBw;
	std::tie(hasHaplotypeInformativeForkBw, hasHaplotypeInformativeForkFw) = getHaplotypeInformativeForkPairs(readForkEdgesBw, readForkEdgesFw, unitigGraph.nodeCount(), forkAlleleCount, approxOneHapCoverage);
	for (size_t i = 0; i < hasHaplotypeInformativeForkFw.size(); i++)
	{
		if (hasHaplotypeInformativeForkFw[i]) std::cerr << "haplotype informative fork >" << i << std::endl;
		if (hasHaplotypeInformativeForkBw[i]) std::cerr << "haplotype informative fork <" << i << std::endl;
	}
	std::vector<std::vector<HaplotypeInformativeSite>> result;
	result.resize(readPaths.size());
	for (size_t i = 0; i < readForkEdgesFw.size(); i++)
	{
		for (auto site : readForkEdgesFw[i])
		{
			if ((site.bubble & firstBitUint64_t) && !hasHaplotypeInformativeForkFw[site.bubble & maskUint64_t]) continue;
			if (((site.bubble ^ firstBitUint64_t) & firstBitUint64_t) && !hasHaplotypeInformativeForkBw[site.bubble & maskUint64_t]) continue;
			result[i].emplace_back(site);
		}
	}
	for (size_t i = 0; i < readForkEdgesBw.size(); i++)
	{
		for (auto site : readForkEdgesBw[i])
		{
			if ((site.bubble & firstBitUint64_t) && !hasHaplotypeInformativeForkBw[site.bubble & maskUint64_t]) continue;
			if (((site.bubble ^ firstBitUint64_t) & firstBitUint64_t) && !hasHaplotypeInformativeForkFw[site.bubble & maskUint64_t]) continue;
			result[i].emplace_back();
			result[i].back().readPos += 1;
			result[i].back().bubble ^= firstBitUint64_t;
			result[i].back().allele ^= firstBitUint64_t;
		}
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphHaploforks(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t resolveLength, const bool cheapUnzip)
{
	std::vector<bool> hasGoodForkFw;
	std::vector<bool> hasGoodForkBw;
	phmap::flat_hash_map<uint64_t, size_t> forkAlleleCount;
	std::tie(hasGoodForkBw, hasGoodForkFw, forkAlleleCount) = getGoodForks(unitigGraph, approxOneHapCoverage);
	for (size_t i = 0; i < hasGoodForkFw.size(); i++)
	{
		if (hasGoodForkFw[i]) std::cerr << "maybe good fork >" << i << std::endl;
		if (hasGoodForkBw[i]) std::cerr << "maybe good fork <" << i << std::endl;
	}
	std::vector<std::vector<HaplotypeInformativeSite>> readForkEdgesFw;
	std::vector<std::vector<HaplotypeInformativeSite>> readForkEdgesBw;
	std::tie(readForkEdgesBw, readForkEdgesFw) = getReadForkEdges(hasGoodForkBw, hasGoodForkFw, unitigGraph, readPaths);
	std::vector<bool> hasHaplotypeInformativeForkFw;
	std::vector<bool> hasHaplotypeInformativeForkBw;
	std::tie(hasHaplotypeInformativeForkBw, hasHaplotypeInformativeForkFw) = getHaplotypeInformativeForkPairs(readForkEdgesBw, readForkEdgesFw, unitigGraph.nodeCount(), forkAlleleCount, approxOneHapCoverage);
	for (size_t i = 0; i < hasHaplotypeInformativeForkFw.size(); i++)
	{
		if (hasHaplotypeInformativeForkFw[i]) std::cerr << "haplotype informative fork >" << i << std::endl;
		if (hasHaplotypeInformativeForkBw[i]) std::cerr << "haplotype informative fork <" << i << std::endl;
	}
	std::vector<bool> informativeNode;
	informativeNode.resize(unitigGraph.nodeCount(), false);
	auto edges = getActiveEdges(unitigGraph.edgeCoverages, unitigGraph.nodeCount());
	for (size_t i = 0; i < hasHaplotypeInformativeForkFw.size(); i++)
	{
		if (hasHaplotypeInformativeForkFw[i])
		{
			for (auto edge : edges.getEdges(std::make_pair(i, true)))
			{
				informativeNode[edge.first] = true;
			}
		}
		if (hasHaplotypeInformativeForkBw[i])
		{
			for (auto edge : edges.getEdges(std::make_pair(i, false)))
			{
				informativeNode[edge.first] = true;
			}
		}
	}
	std::vector<AnchorChain> fakeChains;
	for (size_t i = 0; i < informativeNode.size(); i++)
	{
		if (!informativeNode[i]) continue;
		std::cerr << "fork informative node " << i << std::endl;
		fakeChains.emplace_back();
		fakeChains.back().ploidy = 0;
		fakeChains.back().nodes.emplace_back(i);
		fakeChains.back().nodeOffsets.emplace_back(0);
	}
	std::vector<std::vector<ChainPosition>> fakePositions = getReadChainPositions(unitigGraph, readPaths, fakeChains);
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> hetInfoPerRead;
	hetInfoPerRead.resize(readPaths.size());
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipHapmers(cheapUnzip, fakeChains, unitigGraph, readPaths, hetInfoPerRead, fakePositions, resolveLength, approxOneHapCoverage);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);

}
