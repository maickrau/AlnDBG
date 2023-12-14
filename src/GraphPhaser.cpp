#include <iostream>
#include <tuple>
#include <vector>
#include <phmap.h>
#include "MBGCommon.h"
#include "GraphPhaser.h"
#include "Common.h"
#include "AnchorFinder.h"

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
			std::vector<uint64_t> currAllele;
			currAllele.emplace_back(chainPositionsInReads[readi][j-1].chain ^ firstBitUint64_t);
			size_t currIndex = getAlleleIndex(allelesPerChain, chainPositionsInReads[readi][j].chain & maskUint64_t, currFw ? 0 : anchorChains[chainPositionsInReads[readi][j].chain & maskUint64_t].nodes.size(), currAllele);
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
	if (intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) >= requiredCoveragePerHap && intersectSize(bestAllelesLeft[1], bestAllelesRight[1]) >= requiredCoveragePerHap && normalCounts > (crossCounts+extraNonMatchers)*9) return true;
	if (intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) >= requiredCoveragePerHap && intersectSize(bestAllelesLeft[1], bestAllelesRight[0]) >= requiredCoveragePerHap && crossCounts > (normalCounts+extraNonMatchers)*9) return true;
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
	if (normalCounts >= requiredCoverage && intersectSize(bestAllelesLeft[0], bestAllelesRight[0]) >= 1 && intersectSize(bestAllelesLeft[1], bestAllelesRight[1]) >= 1 && normalCounts > crossCounts*9) return true;
	if (crossCounts >= requiredCoverage && intersectSize(bestAllelesLeft[0], bestAllelesRight[1]) >= 1 && intersectSize(bestAllelesLeft[1], bestAllelesRight[0]) >= 1 && crossCounts > normalCounts*9) return true;
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
	if (notgood[0])
	{
		size_t minorCount;
		std::tie(bestAlleles[0], minorCount) = getTopNCoveredAlleles(readsPerAllele[0], ploidy);
		if (minorCount == 0) notgood[0] = bestAlleles[0].size() != ploidy;
	}
	if (notgood.back())
	{
		size_t minorCount;
		std::tie(bestAlleles.back(), minorCount) = getTopNCoveredAlleles(readsPerAllele.back(), ploidy);
		if (minorCount == 0) notgood.back() = bestAlleles.back().size() != ploidy;
	}
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
				// std::cerr << "extended haplotype informative to bubble index " << i << "/" << readsPerAllele.size() << std::endl;
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

std::pair<size_t, size_t> find(std::vector<std::vector<std::pair<size_t, size_t>>>& parent, std::pair<size_t, size_t> index)
{
	while (true)
	{
		std::pair<size_t, size_t> oneparent = parent[index.first][index.second];
		std::pair<size_t, size_t> twoparent = parent[oneparent.first][oneparent.second];
		if (oneparent == twoparent) break;
		parent[index.first][index.second] = parent[twoparent.first][twoparent.second];
	}
	return parent[index.first][index.second];
}

void merge(std::vector<std::vector<std::pair<size_t, size_t>>>& parent, std::pair<size_t, size_t> left, std::pair<size_t, size_t> right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent[left.first][left.second] == left);
	assert(parent[right.first][right.second] == right);
	parent[right.first][right.second] = left;
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

std::vector<PhaseBlock> phaseCoreChains(std::vector<AnchorChain>& anchorChains, const std::vector<std::vector<ReadDiagonalAlleles>>& allelesPerReadPerChain, const double approxOneHapCoverage)
{
	assert(anchorChains.size() == allelesPerReadPerChain.size());
	std::vector<PhaseBlock> result;
	for (size_t i = 0; i < anchorChains.size(); i++)
	{
		std::cerr << "check chain " << i << " start " << ((anchorChains[i].nodes[0] & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes[0] & maskUint64_t) << " end " << ((anchorChains[i].nodes.back() & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes.back() & maskUint64_t) << " ploidy " << anchorChains[i].ploidy << std::endl;
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
		if (anchorChains[i].ploidy == 2)
		{
			std::cerr << "try phase diploid chain " << i << " start " << ((anchorChains[i].nodes[0] & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes[0] & maskUint64_t) << " end " << ((anchorChains[i].nodes.back() & firstBitUint64_t) ? ">" : "<") << (anchorChains[i].nodes.back() & maskUint64_t) << std::endl;
			auto phaseBlocks = getChainPhaseBlocksMEC(i, allelesPerReadPerChain[i], anchorChains[i].nodes, anchorChains[i].ploidy, anchorChains[i].nodes.size()+1, approxOneHapCoverage);
			if (phaseBlocks.size() == 0) continue;
			std::cerr << "phased with " << phaseBlocks.size() << " blocks. start: " << (phaseBlocks[0].chainStartPhased ? "yes" : "no") << ", end: " << (phaseBlocks.back().chainEndPhased ? "yes" : "no") << std::endl;
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
						std::cerr << "red flag: allele shared between haps, ploidy " << anchorChains[i].ploidy << ", num distinct alleles: " << foundAlleles.size() << ", skipping this phasing" << std::endl;
						include = false;
					}
				}
			}
			if (include) result.insert(result.end(), phaseBlocks.begin(), phaseBlocks.end());
		}
		else
		{
			std::cerr << "skipped chain with ploidy " << anchorChains[i].ploidy << std::endl;
		}
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

void unzipPhaseBlocks(UnitigGraph& resultGraph, std::vector<ReadPathBundle>& resultPaths, const std::vector<AnchorChain>& anchorChains, const std::vector<PhaseBlock>& chainHaplotypes, const std::vector<std::pair<size_t, size_t>>& nodeLocationInChain, const phmap::flat_hash_set<std::pair<size_t, size_t>>& solvedLocations, const phmap::flat_hash_set<std::pair<size_t, size_t>>& solvedAnchors, const phmap::flat_hash_set<size_t>& solvedInterchainTangles, const VectorWithDirection<size_t>& nodeLocationInInterchainTangles, const size_t tangleCount, const SparseEdgeContainer& validChainEdges, const std::vector<std::vector<std::tuple<size_t, bool, int, size_t>>>& haplotypeDiagonalsPerRead, const phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>>& validHaplotypeConnectionsPerChainEdge)
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
					if (solvedLocations.count(std::make_pair(currChain, std::min(lastOffset, currOffset))) == 1)
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
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> unfilteredResult;
	for (const auto& pair : edgeHapCoverage)
	{
		for (const auto& pair2 : pair.second)
		{
			for (const auto& pair3 : pair2.second)
			{
				std::cerr << "hap check chain " << (pair.first & maskUint64_t) << ((pair.first & firstBitUint64_t) ? "+" : "-") << " to chain " << (pair2.first & maskUint64_t) << ((pair2.first & firstBitUint64_t) ? "+" : "-") << " haps " << pair3.first.first << " " << pair3.first.second << " coverage " << pair3.second << std::endl;
				if (pair3.second < approxOneHapCoverage * 0.5) continue;
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
	haplotypeDiagonalsPerRead.resize(resultPaths.size());
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
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedLocations;
	phmap::flat_hash_set<std::pair<size_t, size_t>> solvedAnchors;
	std::vector<bool> canUnzipStart;
	std::vector<bool> canUnzipEnd;
	phmap::flat_hash_map<uint64_t, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<std::pair<size_t, size_t>>>> validHaplotypeConnectionsPerChainEdge;
	std::tie(validHaplotypeConnectionsPerChainEdge, canUnzipStart, canUnzipEnd) = getValidHaplotypeConnectionsPerChainEdge(unitigGraph, readPaths, anchorChains, haplotypeDiagonalsPerRead, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	// for (size_t i = 0; i < anchorChains.size(); i++)
	// {
	// 	std::cerr << "anchor chain " << i << " (" << (anchorChains[i].nodes[0] & maskUint64_t) << " to " << (anchorChains[i].nodes.back() & maskUint64_t) << ") unzippable at start: " << (canUnzipStart[i] ? "yes" : "no") << " at end: " << (canUnzipEnd[i] ? "yes" : "no") << std::endl;
	// }
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
	for (size_t chain = 0; chain < anchorChains.size(); chain++)
	{
		if (anchorChains[chain].ploidy != 1) continue;
		size_t minSolvedIndex = 0;
		size_t maxSolvedIndex = anchorChains[chain].nodes.size();
		// std::cerr << "haploid chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << " (" << anchorChains[chain].nodes.size() << ")" << std::endl;
		if (!canUnzipStart[chain] && minSolvedIndex == 0) minSolvedIndex = 1;
		if (!canUnzipEnd[chain] && maxSolvedIndex == anchorChains[chain].nodes.size()) maxSolvedIndex = anchorChains[chain].nodes.size();
		// std::cerr << "haploid chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << std::endl;
		if (maxSolvedIndex < minSolvedIndex) continue;
		for (size_t j = minSolvedIndex; j <= maxSolvedIndex; j++)
		{
			assert(solvedLocations.count(std::make_pair(chain, j)) == 0);
			// std::cerr << "solve " << chain << " " << j << " (" << (j > 0 ? (anchorChains[chain].nodes[j-1] & maskUint64_t) : -1) << " to " << (j < anchorChains[chain].nodes.size() ? (anchorChains[chain].nodes[j] & maskUint64_t) : std::numeric_limits<size_t>::max()) << ")" << std::endl;
			solvedLocations.emplace(chain, j);
			if (j != maxSolvedIndex) solvedAnchors.emplace(chain, j);
		}
		if (minSolvedIndex > 0) solvedLocations.emplace(chain, minSolvedIndex-1);
	}
	for (size_t i = 0; i < chainHaplotypes.size(); i++)
	{
		const size_t chain = chainHaplotypes[i].chainNumber;
		size_t minSolvedIndex = chainHaplotypes[i].bubbleIndices[0];
		size_t maxSolvedIndex = chainHaplotypes[i].bubbleIndices.back();
		assert(chainHaplotypes[i].bubbleIndices.size() >= 2);
		assert(anchorChains[chainHaplotypes[i].chainNumber].ploidy >= 2);
		// std::cerr << "block " << i << " chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << " (" << anchorChains[chain].nodes.size() << ")" << std::endl;
		if (!canUnzipStart[chain] && minSolvedIndex == 0) minSolvedIndex = chainHaplotypes[i].bubbleIndices[1];
		if (!canUnzipEnd[chain] && maxSolvedIndex == anchorChains[chain].nodes.size()) maxSolvedIndex = chainHaplotypes[i].bubbleIndices[chainHaplotypes[i].bubbleIndices.size()-2];
		// std::cerr << "block " << i << " chain " << chain << " indices " << minSolvedIndex << " " << maxSolvedIndex << std::endl;
		for (size_t j = minSolvedIndex; j <= maxSolvedIndex; j++)
		{
			assert(solvedLocations.count(std::make_pair(chain, j)) == 0);
			// std::cerr << "solve " << chain << " " << j << " (" << (j > 0 ? (anchorChains[chain].nodes[j-1] & maskUint64_t) : -1) << " to " << (j < anchorChains[chain].nodes.size() ? (anchorChains[chain].nodes[j] & maskUint64_t) : std::numeric_limits<size_t>::max()) << ")" << std::endl;
			solvedLocations.emplace(chain, j);
			if (j != maxSolvedIndex) solvedAnchors.emplace(chain, j);
		}
		if (minSolvedIndex > 0) solvedLocations.emplace(chain, minSolvedIndex-1);
	}
	unzipPhaseBlocks(resultGraph, resultPaths, anchorChains, chainHaplotypes, nodeLocationInChain, solvedLocations, solvedAnchors, solvedInterchainTangles, nodeLocationInInterchainTangles, tangleCount, validChainEdges, haplotypeDiagonalsPerRead, validHaplotypeConnectionsPerChainEdge);
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
				if (solvedLocations.count(nodeLocationInChain[i]) == 0) continue;
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

SparseEdgeContainer getValidChainEdges(const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const double approxOneHapCoverage)
{
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
		if (pair.second < approxOneHapCoverage * 0.5) continue;
		// std::cerr << "chain edge between chains " << (pair.first.first.first) << " and " << (pair.first.second.first) << std::endl;
		result.addEdge(pair.first.first, pair.first.second);
		result.addEdge(reverse(pair.first.second), reverse(pair.first.first));
	}
	return result;
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLinearizable(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage)
{
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	// for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	// {
	// 	for (auto chain : chainPositionsInReads[i])
	// 	{
	// 		std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
	// 	}
	// }
	SparseEdgeContainer validChainEdges = getValidChainEdges(chainPositionsInReads, approxOneHapCoverage);
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<ReadDiagonalAlleles>> allelesPerReadPerChain;
	std::tie(allelesPerChain, allelesPerReadPerChain) = getRawReadInfoPerChain(anchorChains, readPaths, unitigGraph, validChainEdges, chainPositionsInReads);
	std::vector<PhaseBlock> chainHaplotypes = phaseCoreChains(anchorChains, allelesPerReadPerChain, approxOneHapCoverage);
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipPhaseBlocks(unitigGraph, readPaths, anchorChains, allelesPerReadPerChain, chainHaplotypes, validChainEdges, chainPositionsInReads, approxOneHapCoverage);
	return std::make_pair(std::move(unzippedGraph), std::move(unzippedReadPaths));
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

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipHapmers(const std::vector<AnchorChain>& anchorChains, const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>>& hetInfoPerRead, const std::vector<std::vector<ChainPosition>>& chainPositionsInReads, const size_t resolveLength)
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
	std::vector<phmap::flat_hash_map<Context, Context>> nodeParent;
	phmap::flat_hash_map<Context, size_t> contextCoverage;
	nodeParent.resize(unitigGraph.nodeCount());
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
	std::vector<std::vector<size_t>> contextCheckStartPoses;
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
	std::cerr << "counting context coverages" << std::endl;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		phmap::flat_hash_set<Context> contexts;
		iterateContextNodes(checkTheseNodes, unitigGraph, readPaths[i], contextCheckStartPoses[i], processedHetInfos[i], chainPositionsInReads[i], resolveLength, [&contexts](const size_t j, const size_t k, const std::vector<Context>& nodecontexts)
		{
			contexts.insert(nodecontexts.begin(), nodecontexts.end());
		});
		for (const auto& context : contexts)
		{
			contextCoverage[context] += 1;
		}
	}
	std::cerr << "merging per-node contexts" << std::endl;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		iterateContextNodes(checkTheseNodes, unitigGraph, readPaths[i], contextCheckStartPoses[i], processedHetInfos[i], chainPositionsInReads[i], resolveLength, [&nodeParent, &contextCoverage, &readPaths, i, minContextCoverage](const size_t j, const size_t k, const std::vector<Context>& contexts)
		{
			uint64_t node = readPaths[i].paths[j].path[k];
			Context anyContext;
			for (const auto& context : contexts)
			{
				if (contextCoverage[context] < minContextCoverage) continue;
				anyContext = context;
				find(nodeParent[node & maskUint64_t], context);
			}
			for (const auto& context : contexts)
			{
				if (contextCoverage[context] < minContextCoverage) continue;
				merge(nodeParent[node & maskUint64_t], context, anyContext);
			}
		});
	}
	std::cerr << "getting context clusters" << std::endl;
	UnitigGraph resultGraph = unitigGraph;
	std::vector<ReadPathBundle> resultPaths = readPaths;
	std::vector<phmap::flat_hash_map<Context, size_t>> contextToNewNodeNumber;
	contextToNewNodeNumber.resize(unitigGraph.nodeCount());
	size_t nextNodeNumber = unitigGraph.nodeCount();
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		phmap::flat_hash_set<Context> keys;
		for (const auto& pair : nodeParent[i])
		{
			keys.insert(pair.first);
		}
		phmap::flat_hash_set<Context> clusters;
		for (const auto& key : keys)
		{
			clusters.insert(find(nodeParent[i], key));
		}
		for (const auto& cluster : clusters)
		{
			contextToNewNodeNumber[i][cluster] = nextNodeNumber;
			resultGraph.lengths.emplace_back(unitigGraph.lengths[i]);
			resultGraph.coverages.emplace_back(0);
			nextNodeNumber += 1;
		}
	}
	std::cerr << "replacing nodes" << std::endl;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		iterateContextNodes(checkTheseNodes, unitigGraph, readPaths[i], contextCheckStartPoses[i], processedHetInfos[i], chainPositionsInReads[i], resolveLength, [&nodeParent, &readPaths, &contextCoverage, &contextToNewNodeNumber, &resultPaths, i, minContextCoverage](const size_t j, const size_t k, const std::vector<Context>& contexts)
		{
			uint64_t node = readPaths[i].paths[j].path[k];
			for (const auto& context : contexts)
			{
				if (contextCoverage[context] < minContextCoverage) continue;
				auto cluster = find(nodeParent[node & maskUint64_t], context);
				if (contextToNewNodeNumber[node & maskUint64_t].count(cluster) == 1)
				{
					// std::cerr << "replace read " << i << " " << j << " " << k << " (" << (node & maskUint64_t) << ") with " << contextToNewNodeNumber[node & maskUint64_t].at(cluster) << std::endl;
					resultPaths[i].paths[j].path[k] = (node & firstBitUint64_t) + (contextToNewNodeNumber[node & maskUint64_t].at(cluster));
					break;
				}
			}
		});
	}
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

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphHapmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t resolveLength)
{
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	// for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	// {
	// 	for (auto chain : chainPositionsInReads[i])
	// 	{
	// 		std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
	// 	}
	// }
	SparseEdgeContainer validChainEdges = getValidChainEdges(chainPositionsInReads, approxOneHapCoverage);
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<ReadDiagonalAlleles>> allelesPerReadPerChain;
	std::tie(allelesPerChain, allelesPerReadPerChain) = getRawReadInfoPerChain(anchorChains, readPaths, unitigGraph, validChainEdges, chainPositionsInReads);
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> hetInfoPerRead = getHetInfoPerRead(unitigGraph, readPaths, anchorChains, allelesPerReadPerChain, readPaths.size(), approxOneHapCoverage);
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipHapmers(anchorChains, unitigGraph, readPaths, hetInfoPerRead, chainPositionsInReads, resolveLength);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}

phmap::flat_hash_map<uint64_t, std::pair<uint64_t, size_t>> getUniqueLocallyUniqueNodeEdges(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const phmap::flat_hash_set<size_t>& notUniques)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<int>> distances;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			uint64_t lastLocalUniq = std::numeric_limits<size_t>::max();
			int lastLocalUniqEndPos = 0;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (notUniques.count(node & maskUint64_t) == 1) continue;
				if (unitigGraph.coverages[node & maskUint64_t] < approxOneHapCoverage * 0.5) continue;
				if (lastLocalUniq == std::numeric_limits<size_t>::max())
				{
					lastLocalUniq = node;
					lastLocalUniqEndPos = readPos;
					continue;
				}
				auto key = canon(std::make_pair(lastLocalUniq & maskUint64_t, lastLocalUniq & firstBitUint64_t), std::make_pair(node & maskUint64_t, node & firstBitUint64_t));
				int distance = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t] - lastLocalUniqEndPos;
				if (distance < 0)
				{
					lastLocalUniq = node;
					lastLocalUniqEndPos = readPos;
					continue;
				}
				if (key.first == std::make_pair<size_t, bool>(lastLocalUniq & maskUint64_t, lastLocalUniq & firstBitUint64_t) && key.second == std::make_pair<size_t, bool>(node & maskUint64_t, node & firstBitUint64_t))
				{
					distances[std::make_pair(lastLocalUniq, node)].emplace_back(distance);
				}
				else
				{
					distances[std::make_pair(node ^ firstBitUint64_t, lastLocalUniq ^firstBitUint64_t)].emplace_back(distance);
				}
				lastLocalUniq = node;
				lastLocalUniqEndPos = readPos;
			}
		}
	}
	phmap::flat_hash_map<uint64_t, std::pair<uint64_t, size_t>> result;
	phmap::flat_hash_set<uint64_t> hasMultipleOutEdges;
	for (auto pair : distances)
	{
		if (pair.second.size() < 3) continue;
		if (result.count(pair.first.first) == 1)
		{
			hasMultipleOutEdges.insert(pair.first.first);
		}
		else
		{
			result[pair.first.first] = std::make_pair(pair.first.second, pair.second[pair.second.size() / 2]);
		}
		if (result.count(pair.first.second ^ firstBitUint64_t) == 1)
		{
			hasMultipleOutEdges.insert(pair.first.second ^ firstBitUint64_t);
		}
		else
		{
			result[pair.first.second ^ firstBitUint64_t] = std::make_pair(pair.first.first ^ firstBitUint64_t, pair.second[pair.second.size()/2]);
		}
	}
	for (auto key : hasMultipleOutEdges)
	{
		result.erase(key);
	}
	return result;
}

void makeFakeLocalUniqGraph(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const phmap::flat_hash_set<size_t>& notUniques)
{
	phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, std::vector<int>> distances;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			size_t readPos = readPaths[i].paths[j].readStartPos;
			uint64_t lastLocalUniq = std::numeric_limits<size_t>::max();
			int lastLocalUniqEndPos = 0;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
				if (notUniques.count(node & maskUint64_t) == 1) continue;
				if (unitigGraph.coverages[node & maskUint64_t] < approxOneHapCoverage * 0.5) continue;
				if (lastLocalUniq == std::numeric_limits<size_t>::max())
				{
					lastLocalUniq = node;
					lastLocalUniqEndPos = readPos;
					continue;
				}
				auto key = canon(std::make_pair(lastLocalUniq & maskUint64_t, lastLocalUniq & firstBitUint64_t), std::make_pair(node & maskUint64_t, node & firstBitUint64_t));
				int distance = (int)readPos - (int)unitigGraph.lengths[node & maskUint64_t] - lastLocalUniqEndPos;
				if (distance < 0)
				{
					lastLocalUniq = node;
					lastLocalUniqEndPos = readPos;
					continue;
				}
				if (key.first == std::make_pair<size_t, bool>(lastLocalUniq & maskUint64_t, lastLocalUniq & firstBitUint64_t) && key.second == std::make_pair<size_t, bool>(node & maskUint64_t, node & firstBitUint64_t))
				{
					distances[std::make_pair(lastLocalUniq, node)].emplace_back(distance);
				}
				else
				{
					distances[std::make_pair(node ^ firstBitUint64_t, lastLocalUniq ^firstBitUint64_t)].emplace_back(distance);
				}
				lastLocalUniq = node;
				lastLocalUniqEndPos = readPos;
			}
		}
	}
	std::ofstream graph { "localuniqgraph.gfa" };
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (unitigGraph.coverages[i] < approxOneHapCoverage * 0.5) continue;
		if (notUniques.count(i) == 1) continue;
		graph << "S\t" << i << "\t*\tLN:i:" << unitigGraph.lengths[i] << "\tll:i:" << unitigGraph.coverages[i] << "\tFC:i:" << (unitigGraph.coverages[i]*unitigGraph.lengths[i]) << std::endl;
	}
	size_t edgenum = 0;
	for (auto pair : distances)
	{
		if (pair.second.size() < 3) continue;
		int distance = pair.second[pair.second.size()/2];
		if (distance < 1) distance = 1;
		graph << "S\tedge_" << edgenum << "\t*\tLN:i:" << distance << "\tll:i:" << pair.second.size() << "\tFC:i:" << (distance * pair.second.size()) << std::endl;
		graph << "L\t" << (pair.first.first & maskUint64_t) << "\t" << ((pair.first.first & firstBitUint64_t) ? "+" : "-") << "\tedge_" << edgenum << "\t+\t0M\tec:i:" << pair.second.size() << std::endl;
		graph << "L\tedge_" << edgenum << "\t+\t" << (pair.first.second & maskUint64_t) << "\t" << ((pair.first.second & firstBitUint64_t) ? "+" : "-") << "\t0M\tec:i:" << pair.second.size() << std::endl;
		edgenum += 1;
	}
}

AnchorChain getLocallyUniqueChain(const phmap::flat_hash_map<uint64_t, std::pair<uint64_t, size_t>>& edges, const uint64_t startNode, const UnitigGraph& unitigGraph)
{
	uint64_t pos = startNode ^ firstBitUint64_t;
	while (edges.count(pos) == 1)
	{
		uint64_t next = edges.at(pos).first;
		if (edges.count(next ^ firstBitUint64_t) == 0) break;
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
		uint64_t next = edges.at(pos).first;
		if (edges.count(next ^ firstBitUint64_t) == 0) break;
		if ((next & maskUint64_t) == (pos & maskUint64_t)) break;
		if ((next ^ firstBitUint64_t) == chainstart) break;
		result.nodes.push_back(next);
		result.nodeOffsets.push_back(result.nodeOffsets.back() + edges.at(pos).second + unitigGraph.lengths[pos & maskUint64_t]);
		pos = next;
	}
	return result;
}

std::vector<std::vector<ChainPosition>> getReadLocalUniqChains(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t minLength, const size_t resolveLength)
{
	const size_t notSameNodeDistance = 100;
	phmap::flat_hash_set<size_t> notUniques;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			phmap::flat_hash_set<size_t> found;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				size_t node = readPaths[i].paths[j].path[k] & maskUint64_t;
				if (found.count(node) == 1) notUniques.insert(node);
				found.insert(node);
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
				if (notUniques.count(node & maskUint64_t) == 1) continue;
				nodePositions[node].emplace_back(readPos);
			}
		}
		for (auto& pair : nodePositions)
		{
			if (notUniques.count(pair.first & maskUint64_t) == 1) continue;
			std::sort(pair.second.begin(), pair.second.end());
			for (size_t j = 1; j < pair.second.size(); j++)
			{
				if (pair.second[j] < pair.second[j-1] + notSameNodeDistance) continue;
				if (pair.second[j] > pair.second[j-1] + resolveLength) continue;
				notUniques.insert(pair.first);
				break;
			}
		}
	}
	makeFakeLocalUniqGraph(unitigGraph, readPaths, approxOneHapCoverage, notUniques);
	auto edges = getUniqueLocallyUniqueNodeEdges(unitigGraph, readPaths, approxOneHapCoverage, notUniques);
	std::vector<bool> checked;
	checked.resize(unitigGraph.nodeCount(), false);
	std::vector<AnchorChain> fakeChains;
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		if (checked[i]) continue;
		if (notUniques.count(i) == 1) continue;
		if (unitigGraph.coverages[i] < approxOneHapCoverage * 0.5) continue;
		AnchorChain chain = getLocallyUniqueChain(edges, i, unitigGraph);
		assert(chain.nodes.size() >= 1);
		size_t totalChainLength = 0;
		for (uint64_t node : chain.nodes)
		{
			assert(!checked[node & maskUint64_t]);
			checked[node & maskUint64_t] = true;
			totalChainLength += unitigGraph.lengths[node & maskUint64_t];
		}
		if (totalChainLength < minLength) continue;
		fakeChains.emplace_back(chain);
		std::cerr << "local-uniq chain length: " << totalChainLength << " kmers nodes:";
		for (uint64_t node : chain.nodes)
		{
			std::cerr << " " << (node & maskUint64_t);
		}
		std::cerr << std::endl;
	}
	std::cerr << fakeChains.size() << " local-uniq chains" << std::endl;
	return getReadChainPositions(unitigGraph, readPaths, fakeChains);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLocalUniqmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t resolveLength)
{
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> localUniqPositions = getReadLocalUniqChains(unitigGraph, readPaths, approxOneHapCoverage, 100, 2000);
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
	std::tie(unzippedGraph, unzippedReadPaths) = unzipHapmers(anchorChains, unitigGraph, readPaths, hetInfoPerRead, localUniqPositions, resolveLength);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphChainmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t resolveLength)
{
	std::vector<AnchorChain> anchorChains = getAnchorChains(unitigGraph, readPaths, approxOneHapCoverage);
	std::cerr << anchorChains.size() << " anchor chains" << std::endl;
	std::vector<std::vector<ChainPosition>> chainPositionsInReads = getReadChainPositions(unitigGraph, readPaths, anchorChains);
	// for (size_t i = 0; i < chainPositionsInReads.size(); i++)
	// {
	// 	for (auto chain : chainPositionsInReads[i])
	// 	{
	// 		std::cerr << "read " << i << " chain " << (chain.chain & maskUint64_t) << " " << ((chain.chain & firstBitUint64_t) ? "fw" : "bw") << " (" << (anchorChains[chain.chain].nodes[0] & maskUint64_t) << " " << (anchorChains[chain.chain].nodes.back() & maskUint64_t) << ") " << chain.chainStartPosInRead << " " << chain.chainEndPosInRead << std::endl;
	// 	}
	// }
	SparseEdgeContainer validChainEdges = getValidChainEdges(chainPositionsInReads, approxOneHapCoverage);
	std::vector<std::vector<std::vector<std::vector<uint64_t>>>> allelesPerChain;
	std::vector<std::vector<ReadDiagonalAlleles>> allelesPerReadPerChain;
	std::tie(allelesPerChain, allelesPerReadPerChain) = getRawReadInfoPerChain(anchorChains, readPaths, unitigGraph, validChainEdges, chainPositionsInReads);
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, bool, size_t, size_t>>> hetInfoPerRead;
	hetInfoPerRead.resize(readPaths.size());
	UnitigGraph unzippedGraph;
	std::vector<ReadPathBundle> unzippedReadPaths;
	std::tie(unzippedGraph, unzippedReadPaths) = unzipHapmers(anchorChains, unitigGraph, readPaths, hetInfoPerRead, chainPositionsInReads, resolveLength);
	RankBitvector kept;
	kept.resize(unzippedGraph.nodeCount());
	for (size_t i = 0; i < unzippedGraph.nodeCount(); i++)
	{
		kept.set(i, unzippedGraph.coverages[i] >= 2);
	}
	return filterUnitigGraph(unzippedGraph, unzippedReadPaths, kept);
}
