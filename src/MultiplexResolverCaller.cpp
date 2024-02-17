#include "MultiplexResolverCaller.h"
#include "UnitigResolver.h"
#include "UnitigHelper.h"

std::tuple<MBG::HashList, MBG::UnitigGraph, std::vector<MBG::ReadPath>> translateToMBGDataStructures(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk)
{
	MBG::HashList fakeHashList { graphk };
	MBG::UnitigGraph fakeUnitigGraph;
	fakeUnitigGraph.unitigs.resize(unitigGraph.nodeCount());
	fakeUnitigGraph.leftClip.resize(unitigGraph.nodeCount(), 0);
	fakeUnitigGraph.rightClip.resize(unitigGraph.nodeCount(), 0);
	fakeUnitigGraph.unitigCoverage.resize(unitigGraph.nodeCount());
	fakeUnitigGraph.edges.resize(unitigGraph.nodeCount());
	fakeUnitigGraph.edgeCov.resize(unitigGraph.nodeCount());
	fakeUnitigGraph.edgeOvlp.resize(unitigGraph.nodeCount());
	std::vector<RankBitvector> fakeHashPositions;
	fakeHashPositions.resize(unitigGraph.nodeCount());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		fakeHashPositions[i].resize(unitigGraph.lengths[i]);
		fakeHashPositions[i].set(0, true);
		fakeHashPositions[i].set(unitigGraph.lengths[i]-1, true);
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			if (readPaths[i].paths[j].pathLeftClipKmers != 0)
			{
				if (readPaths[i].paths[j].path[0] & firstBitUint64_t)
				{
					fakeHashPositions[readPaths[i].paths[j].path[0] & maskUint64_t].set(readPaths[i].paths[j].pathLeftClipKmers, true);
				}
				else
				{
					fakeHashPositions[readPaths[i].paths[j].path[0] & maskUint64_t].set(fakeHashPositions[readPaths[i].paths[j].path[0] & maskUint64_t].size()-1-readPaths[i].paths[j].pathLeftClipKmers, true);
				}
			}
			if (readPaths[i].paths[j].pathRightClipKmers != 0)
			{
				if (readPaths[i].paths[j].path.back() & firstBitUint64_t)
				{
					fakeHashPositions[readPaths[i].paths[j].path.back() & maskUint64_t].set(fakeHashPositions[readPaths[i].paths[j].path.back() & maskUint64_t].size()-1-readPaths[i].paths[j].pathRightClipKmers, true);
				}
				else
				{
					fakeHashPositions[readPaths[i].paths[j].path.back() & maskUint64_t].set(readPaths[i].paths[j].pathRightClipKmers, true);
				}
			}
		}
	}
	for (size_t i = 0; i < fakeHashPositions.size(); i++)
	{
		size_t lastHashPos = 0;
		for (size_t j = 1; j < fakeHashPositions[i].size(); j++)
		{
			if (fakeHashPositions[i].get(j))
			{
				lastHashPos = j;
				continue;
			}
			if (j == lastHashPos+graphk-1)
			{
				fakeHashPositions[i].set(j, true);
				lastHashPos = j;
				continue;
			}
		}
		fakeHashPositions[i].buildRanks();
	}
	std::vector<size_t> unitigFirstHash;
	std::vector<size_t> unitigLastHash;
	unitigFirstHash.resize(unitigGraph.nodeCount(), std::numeric_limits<size_t>::max());
	unitigLastHash.resize(unitigGraph.nodeCount(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		size_t lastHash = fakeHashList.addFakeNode();
		size_t lastHashPos = 0;
		unitigFirstHash[i] = lastHash;
		fakeUnitigGraph.unitigs[i].emplace_back(lastHash, true);
		fakeUnitigGraph.unitigCoverage[i].emplace_back(0);
		for (size_t j = 1; j < unitigGraph.lengths[i]; j++)
		{
			if (!fakeHashPositions[i].get(j)) continue;
			size_t thisHash = fakeHashList.addFakeNode();
			fakeHashList.addSequenceOverlap(std::make_pair(lastHash, true), std::make_pair(thisHash, true), graphk - (j - lastHashPos));
			fakeHashList.addEdgeCoverage(std::make_pair(lastHash, true), std::make_pair(thisHash, true));
			lastHashPos = j;
			lastHash = thisHash;
			fakeUnitigGraph.unitigs[i].emplace_back(lastHash, true);
			fakeUnitigGraph.unitigCoverage[i].emplace_back(0);
		}
		assert(lastHashPos+1 == unitigGraph.lengths[i]);
		unitigLastHash[i] = lastHash;
	}
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				size_t startIndex = 0;
				size_t endIndex = fakeHashPositions[node & maskUint64_t].getRank(fakeHashPositions[node & maskUint64_t].size()-1)+1;
				if (k == 0)
				{
					if (node & firstBitUint64_t)
					{
						assert(fakeHashPositions[node & maskUint64_t].get(readPaths[i].paths[j].pathLeftClipKmers));
						startIndex = fakeHashPositions[node & maskUint64_t].getRank(readPaths[i].paths[j].pathLeftClipKmers);
					}
					else
					{
						assert(fakeHashPositions[node & maskUint64_t].get(fakeHashPositions[node & maskUint64_t].size()-1-readPaths[i].paths[j].pathLeftClipKmers));
						endIndex = fakeHashPositions[node & maskUint64_t].getRank(fakeHashPositions[node & maskUint64_t].size()-1-readPaths[i].paths[j].pathLeftClipKmers)+1;
					}
				}
				if (k == readPaths[i].paths[j].path.size()-1)
				{
					if (node & firstBitUint64_t)
					{
						assert(fakeHashPositions[node & maskUint64_t].get(fakeHashPositions[node & maskUint64_t].size()-1-readPaths[i].paths[j].pathRightClipKmers));
						endIndex = fakeHashPositions[node & maskUint64_t].getRank(fakeHashPositions[node & maskUint64_t].size()-1-readPaths[i].paths[j].pathRightClipKmers)+1;
					}
					else
					{
						assert(fakeHashPositions[node & maskUint64_t].get(readPaths[i].paths[j].pathRightClipKmers));
						startIndex = fakeHashPositions[node & maskUint64_t].getRank(readPaths[i].paths[j].pathRightClipKmers);
					}
				}
				assert(endIndex > startIndex);
				assert(endIndex <= fakeUnitigGraph.unitigCoverage[node & maskUint64_t].size());
				for (size_t n = startIndex; n < endIndex; n++)
				{
					fakeUnitigGraph.unitigCoverage[node & maskUint64_t][n] += 1;
				}
			}
			for (size_t k = 1; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t thisNode = readPaths[i].paths[j].path[k];
				uint64_t prevNode = readPaths[i].paths[j].path[k-1];
				std::pair<size_t, bool> prevNodePair { prevNode & maskUint64_t, prevNode & firstBitUint64_t };
				std::pair<size_t, bool> thisNodePair { thisNode & maskUint64_t, thisNode & firstBitUint64_t };
				fakeUnitigGraph.setEdgeOverlap(prevNodePair, thisNodePair, graphk-1);
				size_t coverage = 0;
				auto key = MBG::canon(prevNodePair, thisNodePair);
				if (fakeUnitigGraph.edgeCov.hasValue(key.first, key.second))
				{
					coverage = fakeUnitigGraph.edgeCov.get(key.first, key.second);
				}
				fakeUnitigGraph.setEdgeCoverage(prevNodePair, thisNodePair, coverage+1);
			}
		}
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto key : fakeUnitigGraph.edgeCov.getValues(fw))
		{
			fakeUnitigGraph.edges.addEdge(fw, key.first);
			fakeUnitigGraph.edges.addEdge(reverse(key.first), reverse(fw));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto key : fakeUnitigGraph.edgeCov.getValues(bw))
		{
			fakeUnitigGraph.edges.addEdge(bw, key.first);
			fakeUnitigGraph.edges.addEdge(reverse(key.first), reverse(bw));
		}
	}
	for (size_t i = 0; i < unitigGraph.nodeCount(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : fakeUnitigGraph.edges.getEdges(fw))
		{
			std::pair<size_t, bool> from { unitigLastHash[i], true };
			std::pair<size_t, bool> to;
			if (edge.second)
			{
				to = std::make_pair(unitigFirstHash[edge.first], true);
			}
			else
			{
				to = std::make_pair(unitigLastHash[edge.first], false);
			}
			fakeHashList.addSequenceOverlap(from, to, graphk-1);
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : fakeUnitigGraph.edges.getEdges(bw))
		{
			std::pair<size_t, bool> from { unitigFirstHash[i], false };
			std::pair<size_t, bool> to;
			if (edge.second)
			{
				to = std::make_pair(unitigFirstHash[edge.first], true);
			}
			else
			{
				to = std::make_pair(unitigLastHash[edge.first], false);
			}
			fakeHashList.addSequenceOverlap(from, to, graphk-1);
		}
	}
	std::vector<MBG::ReadPath> fakeReadPaths;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		for (size_t j = 0; j < readPaths[i].paths.size(); j++)
		{
			fakeReadPaths.emplace_back();
			fakeReadPaths.back().readName.first = std::to_string(i);
			fakeReadPaths.back().readName.second = j;
			fakeReadPaths.back().leftClip = 0;
			fakeReadPaths.back().rightClip = 0;
			for (size_t n = 0; n < readPaths[i].paths[j].pathLeftClipKmers; n++)
			{
				if (readPaths[i].paths[j].path[0] & firstBitUint64_t)
				{
					if (fakeHashPositions[readPaths[i].paths[j].path[0] & maskUint64_t].get(n)) fakeReadPaths.back().leftClip += 1;
				}
				else
				{
					if (fakeHashPositions[readPaths[i].paths[j].path[0] & maskUint64_t].get(fakeHashPositions[readPaths[i].paths[j].path[0] & maskUint64_t].size()-1-n)) fakeReadPaths.back().leftClip += 1;
				}
			}
			for (size_t n = 0; n < readPaths[i].paths[j].pathRightClipKmers; n++)
			{
				if (readPaths[i].paths[j].path.back() & firstBitUint64_t)
				{
					if (fakeHashPositions[readPaths[i].paths[j].path.back() & maskUint64_t].get(fakeHashPositions[readPaths[i].paths[j].path.back() & maskUint64_t].size()-1-n)) fakeReadPaths.back().rightClip += 1;
				}
				else
				{
					if (fakeHashPositions[readPaths[i].paths[j].path.back() & maskUint64_t].get(n)) fakeReadPaths.back().rightClip += 1;
				}
			}
			fakeReadPaths.back().readLength = readPaths[i].readLength;
			fakeReadPaths.back().readLengthHPC = readPaths[i].readLength;
			size_t readPos = readPaths[i].paths[j].readStartPos;
			for (size_t k = 0; k < readPaths[i].paths[j].path.size(); k++)
			{
				uint64_t node = readPaths[i].paths[j].path[k];
				fakeReadPaths.back().path.emplace_back(node & maskUint64_t, node & firstBitUint64_t);
				size_t nodeStartPos = 0;
				size_t nodeEndPos = unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) nodeStartPos = readPaths[i].paths[j].pathLeftClipKmers;
				if (k+1 == readPaths[i].paths[j].path.size()) nodeEndPos = unitigGraph.lengths[node & maskUint64_t] - readPaths[i].paths[j].pathRightClipKmers;
				assert(nodeEndPos > 0);
				size_t added = 0;
				if (node & firstBitUint64_t)
				{
					assert(fakeHashPositions[node & maskUint64_t].get(nodeStartPos));
					assert(fakeHashPositions[node & maskUint64_t].get(nodeEndPos-1));
					for (size_t n = nodeStartPos; n < nodeEndPos; n++)
					{
						if (!fakeHashPositions[node & maskUint64_t].get(n)) continue;
						fakeReadPaths.back().readPoses.emplace_back(readPos + n - nodeStartPos);
						added += 1;
					}
					assert(added == fakeHashPositions[node & maskUint64_t].getRank(nodeEndPos-1)+1 - fakeHashPositions[node & maskUint64_t].getRank(nodeStartPos));
				}
				else
				{
					assert(fakeHashPositions[node & maskUint64_t].get(fakeHashPositions[node & maskUint64_t].size()-1-nodeStartPos));
					assert(fakeHashPositions[node & maskUint64_t].get(fakeHashPositions[node & maskUint64_t].size()-1-(nodeEndPos-1)));
					for (size_t n = nodeStartPos; n < nodeEndPos; n++)
					{
						if (!fakeHashPositions[node & maskUint64_t].get(fakeHashPositions[node & maskUint64_t].size()-1-n)) continue;
						fakeReadPaths.back().readPoses.emplace_back(readPos + n - nodeStartPos);
						added += 1;
					}
					assert(added == fakeHashPositions[node & maskUint64_t].getRank(fakeHashPositions[node & maskUint64_t].size()-1-nodeStartPos)+1 - fakeHashPositions[node & maskUint64_t].getRank(fakeHashPositions[node & maskUint64_t].size()-nodeEndPos));
				}
				assert(k == 0 || k+1 == readPaths[i].paths[j].path.size() || added == fakeUnitigGraph.unitigs[node & maskUint64_t].size());
				readPos += unitigGraph.lengths[node & maskUint64_t];
				if (k == 0) readPos -= readPaths[i].paths[j].pathLeftClipKmers;
			}
			readPos -= readPaths[i].paths[j].pathRightClipKmers;
			assert(fakeReadPaths.back().readPoses[0] == readPaths[i].paths[j].readStartPos);
			if (!(fakeReadPaths.back().readPoses.back()+1 == readPos))
			{
				std::cerr << fakeReadPaths.back().readPoses.back()+1 << " " << readPos << std::endl;
			}
			assert(fakeReadPaths.back().readPoses.back()+1 == readPos);
		}
	}
	return std::make_tuple(fakeHashList, fakeUnitigGraph, fakeReadPaths);
}

std::tuple<UnitigGraph, std::vector<ReadPathBundle>> translateToAlnDBGDataStructures(const MBG::UnitigGraph& fakeUnitigGraph, const std::vector<MBG::ReadPath>& fakeReadPaths, const MBG::HashList& fakeHashList, const size_t graphk, const std::vector<ReadPathBundle>& oldReads, const UnitigGraph& oldUnitigGraph)
{
	UnitigGraph result;
	std::vector<ReadPathBundle> resultPaths;
	for (size_t i = 0; i < fakeUnitigGraph.unitigs.size(); i++)
	{
		result.lengths.emplace_back(MBG::getUnitigSize(fakeHashList, graphk, fakeUnitigGraph, i) - (graphk-1));
		result.coverages.emplace_back(0);
	}
	result.edgeCoverages.resize(result.lengths.size());
	result.edgeKmerOverlaps.resize(result.lengths.size());
	for (size_t i = 0; i < fakeUnitigGraph.unitigs.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		for (auto edge : fakeUnitigGraph.edges.getEdges(fw))
		{
			size_t overlap = MBG::getUnitigOverlap(fakeHashList, graphk, fakeUnitigGraph, fw, edge);
			assert(overlap >= graphk-1);
			result.edgeKmerOverlaps.set(fw, edge, overlap-(graphk-1));
		}
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : fakeUnitigGraph.edges.getEdges(bw))
		{
			size_t overlap = MBG::getUnitigOverlap(fakeHashList, graphk, fakeUnitigGraph, bw, edge);
			assert(overlap >= graphk-1);
			result.edgeKmerOverlaps.set(bw, edge, overlap-(graphk-1));
		}
	}
	resultPaths.resize(oldReads.size());
	for (size_t i = 0; i < oldReads.size(); i++)
	{
		resultPaths[i].readName = oldReads[i].readName;
		resultPaths[i].readLength = oldReads[i].readLength;
	}
	for (size_t i = 0; i < fakeReadPaths.size(); i++)
	{
		size_t readIndex = std::stoi(fakeReadPaths[i].readName.first);
		size_t oldPathIndex = fakeReadPaths[i].readName.second;
		assert(fakeReadPaths[i].readPoses.size() >= 1);
		assert(fakeReadPaths[i].path.size() >= 1);
		resultPaths[readIndex].paths.emplace_back();
		resultPaths[readIndex].paths.back().readStartPos = fakeReadPaths[i].readPoses[0];
		assert(resultPaths[readIndex].paths.back().readStartPos < resultPaths[readIndex].readLength);
		resultPaths[readIndex].paths.back().pathLeftClipKmers = 0;
		resultPaths[readIndex].paths.back().pathRightClipKmers = 0;
		std::vector<std::pair<size_t, bool>> checkHashSequence;
		for (size_t j = 0; j < fakeReadPaths[i].path.size(); j++)
		{
			std::vector<std::pair<size_t, bool>> add = fakeUnitigGraph.unitigs[fakeReadPaths[i].path[j].id()];
			if (!fakeReadPaths[i].path[j].forward())
			{
				std::reverse(add.begin(), add.end());
				for (auto& node : add) node = MBG::reverse(node);
			}
			size_t overlap = 0;
			if (j > 0)
			{
				overlap = fakeUnitigGraph.edgeOverlap(fakeReadPaths[i].path[j-1], fakeReadPaths[i].path[j]);
				assert(add.size() >= overlap);
				assert(checkHashSequence.size() >= overlap);
				for (size_t k = 0; k < overlap; k++)
				{
					assert(checkHashSequence[checkHashSequence.size()-overlap+k] == add[k]);
				}
			}
			checkHashSequence.insert(checkHashSequence.end(), add.begin()+overlap, add.end());
		}
		assert(checkHashSequence.size() > fakeReadPaths[i].leftClip + fakeReadPaths[i].rightClip);
		if (fakeReadPaths[i].leftClip > 0) checkHashSequence.erase(checkHashSequence.begin(), checkHashSequence.begin()+fakeReadPaths[i].leftClip);
		if (fakeReadPaths[i].rightClip > 0) checkHashSequence.erase(checkHashSequence.end()-fakeReadPaths[i].rightClip, checkHashSequence.end());
		assert(checkHashSequence.size() == fakeReadPaths[i].readPoses.size());
		size_t checkPos = 0;
		for (size_t j = 1; j < checkHashSequence.size(); j++)
		{
			checkPos += graphk;
			checkPos -= fakeHashList.getOverlap(checkHashSequence[j-1], checkHashSequence[j]);
			assert(fakeReadPaths[i].readPoses[j] - fakeReadPaths[i].readPoses[0] == checkPos);
		}
		assert(fakeReadPaths[i].leftClip < fakeUnitigGraph.unitigs[fakeReadPaths[i].path[0].id()].size());
		assert(fakeReadPaths[i].rightClip < fakeUnitigGraph.unitigs[fakeReadPaths[i].path.back().id()].size());
		assert(fakeReadPaths[i].path.size() >= 2 || fakeReadPaths[i].leftClip + fakeReadPaths[i].rightClip < fakeUnitigGraph.unitigs[fakeReadPaths[i].path[0].id()].size());
		for (size_t j = 0; j < fakeReadPaths[i].leftClip; j++)
		{
			resultPaths[readIndex].paths.back().pathLeftClipKmers += graphk;
			if (fakeReadPaths[i].path[0].forward())
			{
				resultPaths[readIndex].paths.back().pathLeftClipKmers -= fakeHashList.getOverlap(fakeUnitigGraph.unitigs[fakeReadPaths[i].path[0].id()][j], fakeUnitigGraph.unitigs[fakeReadPaths[i].path[0].id()][j+1]);
			}
			else
			{
				resultPaths[readIndex].paths.back().pathLeftClipKmers -= fakeHashList.getOverlap(fakeUnitigGraph.unitigs[fakeReadPaths[i].path[0].id()][fakeUnitigGraph.unitigs[fakeReadPaths[i].path[0].id()].size()-2-j], fakeUnitigGraph.unitigs[fakeReadPaths[i].path[0].id()][fakeUnitigGraph.unitigs[fakeReadPaths[i].path[0].id()].size()-1-j]);
			}
		}
		if (fakeReadPaths[i].path[0].forward())
		{
			if (resultPaths[readIndex].paths.back().pathLeftClipKmers >= fakeUnitigGraph.leftClip[fakeReadPaths[i].path[0].id()])
			{
				resultPaths[readIndex].paths.back().pathLeftClipKmers -= fakeUnitigGraph.leftClip[fakeReadPaths[i].path[0].id()];
			}
			else
			{
				resultPaths[readIndex].paths.back().readStartPos += fakeUnitigGraph.leftClip[fakeReadPaths[i].path[0].id()] - resultPaths[readIndex].paths.back().pathLeftClipKmers;
				resultPaths[readIndex].paths.back().pathLeftClipKmers = 0;
			}
		}
		else
		{
			if (resultPaths[readIndex].paths.back().pathLeftClipKmers >= fakeUnitigGraph.rightClip[fakeReadPaths[i].path[0].id()])
			{
				resultPaths[readIndex].paths.back().pathLeftClipKmers -= fakeUnitigGraph.rightClip[fakeReadPaths[i].path[0].id()];
			}
			else
			{
				resultPaths[readIndex].paths.back().readStartPos += fakeUnitigGraph.rightClip[fakeReadPaths[i].path[0].id()] - resultPaths[readIndex].paths.back().pathLeftClipKmers;
				resultPaths[readIndex].paths.back().pathLeftClipKmers = 0;
			}
		}
		for (size_t j = 0; j < fakeReadPaths[i].rightClip; j++)
		{
			resultPaths[readIndex].paths.back().pathRightClipKmers += graphk;
			if (fakeReadPaths[i].path.back().forward())
			{
				resultPaths[readIndex].paths.back().pathRightClipKmers -= fakeHashList.getOverlap(fakeUnitigGraph.unitigs[fakeReadPaths[i].path.back().id()][fakeUnitigGraph.unitigs[fakeReadPaths[i].path.back().id()].size()-2-j], fakeUnitigGraph.unitigs[fakeReadPaths[i].path.back().id()][fakeUnitigGraph.unitigs[fakeReadPaths[i].path.back().id()].size()-1-j]);
			}
			else
			{
				resultPaths[readIndex].paths.back().pathRightClipKmers -= fakeHashList.getOverlap(fakeUnitigGraph.unitigs[fakeReadPaths[i].path.back().id()][j], fakeUnitigGraph.unitigs[fakeReadPaths[i].path.back().id()][j+1]);
			}
		}
		if (fakeReadPaths[i].path.back().forward())
		{
			if (resultPaths[readIndex].paths.back().pathRightClipKmers >= fakeUnitigGraph.rightClip[fakeReadPaths[i].path.back().id()])
			{
				resultPaths[readIndex].paths.back().pathRightClipKmers -= fakeUnitigGraph.rightClip[fakeReadPaths[i].path.back().id()];
			}
			else
			{
				// resultPaths[readIndex].paths.back().readStartPos += fakeUnitigGraph.rightClip[fakeReadPaths[i].path.back().id()] - resultPaths[readIndex].paths.back().pathRightClipKmers;
				resultPaths[readIndex].paths.back().pathRightClipKmers = 0;
			}
		}
		else
		{
			if (resultPaths[readIndex].paths.back().pathRightClipKmers >= fakeUnitigGraph.leftClip[fakeReadPaths[i].path.back().id()])
			{
				resultPaths[readIndex].paths.back().pathRightClipKmers -= fakeUnitigGraph.leftClip[fakeReadPaths[i].path.back().id()];
			}
			else
			{
				// resultPaths[readIndex].paths.back().readStartPos += fakeUnitigGraph.leftClip[fakeReadPaths[i].path.back().id()] - resultPaths[readIndex].paths.back().pathRightClipKmers;
				resultPaths[readIndex].paths.back().pathRightClipKmers = 0;
			}
		}
		for (auto node : fakeReadPaths[i].path)
		{
			resultPaths[readIndex].paths.back().path.emplace_back(node.id() + (node.forward() ? firstBitUint64_t : 0));
		}
		assert(resultPaths[readIndex].paths.back().pathLeftClipKmers < result.lengths[fakeReadPaths[i].path[0].id()]);
		assert(resultPaths[readIndex].paths.back().pathRightClipKmers < result.lengths[fakeReadPaths[i].path.back().id()]);
		assert(resultPaths[readIndex].paths.back().path.size() >= 2 || resultPaths[readIndex].paths.back().pathLeftClipKmers + resultPaths[readIndex].paths.back().pathRightClipKmers < result.lengths[fakeReadPaths[i].path[0].id()]);
		// assert(resultPaths[readIndex].paths.back().pathLeftClipKmers < result.lengths[fakeReadPaths[i].path[0].id()] + fakeUnitigGraph.leftClip[fakeReadPaths[i].path[0].id()] + fakeUnitigGraph.rightClip[fakeReadPaths[i].path[0].id()]);
		// assert(resultPaths[readIndex].paths.back().pathRightClipKmers < result.lengths[fakeReadPaths[i].path.back().id()] + fakeUnitigGraph.leftClip[fakeReadPaths[i].path.back().id()] + fakeUnitigGraph.rightClip[fakeReadPaths[i].path.back().id()]);
		// assert(resultPaths[readIndex].paths.back().path.size() >= 2 || resultPaths[readIndex].paths.back().pathLeftClipKmers + resultPaths[readIndex].paths.back().pathRightClipKmers < result.lengths[fakeReadPaths[i].path[0].id()] + fakeUnitigGraph.leftClip[fakeReadPaths[i].path[0].id()] + fakeUnitigGraph.rightClip[fakeReadPaths[i].path[0].id()]);
		assert(resultPaths[readIndex].paths.back().readStartPos >= oldReads[readIndex].paths[oldPathIndex].readStartPos);
		size_t readEndPos = resultPaths[readIndex].paths.back().readStartPos;
		size_t oldReadEndPos = oldReads[readIndex].paths[oldPathIndex].readStartPos;
		for (size_t k = 0; k < oldReads[readIndex].paths[oldPathIndex].path.size(); k++)
		{
			uint64_t node = oldReads[readIndex].paths[oldPathIndex].path[k];
			oldReadEndPos += oldUnitigGraph.lengths[node & maskUint64_t];
		}
		oldReadEndPos -= oldReads[readIndex].paths[oldPathIndex].pathLeftClipKmers;
		oldReadEndPos -= oldReads[readIndex].paths[oldPathIndex].pathRightClipKmers;
		for (size_t k = 0; k < resultPaths[readIndex].paths.back().path.size(); k++)
		{
			uint64_t node = resultPaths[readIndex].paths.back().path[k];
			readEndPos += result.lengths[node & maskUint64_t];
			if (k > 0)
			{
				uint64_t prevNode = resultPaths[readIndex].paths.back().path[k-1];
				readEndPos -= result.edgeKmerOverlaps.get(std::make_pair(prevNode & maskUint64_t, prevNode & firstBitUint64_t), std::make_pair(node & maskUint64_t, node & firstBitUint64_t));
			}
		}
		readEndPos -= resultPaths[readIndex].paths.back().pathLeftClipKmers;
		readEndPos -= resultPaths[readIndex].paths.back().pathRightClipKmers;
		assert(resultPaths[readIndex].paths.back().readStartPos < readEndPos);
		assert(readEndPos <= oldReadEndPos);
		assert(readEndPos <= resultPaths[readIndex].readLength);
	}
	for (size_t i = 0; i < resultPaths.size(); i++)
	{
		std::sort(resultPaths[i].paths.begin(), resultPaths[i].paths.end(), [](const auto& left, const auto& right) { return left.readStartPos < right.readStartPos; });
	}
	return unitigify(result, resultPaths);
	// return std::make_pair(result, resultPaths);
}

std::pair<UnitigGraph, std::vector<ReadPathBundle>> runMBGMultiplexResolution(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const size_t maxk)
{
	std::cerr << "do multiplexing" << std::endl;
	std::cerr << "switch to MBG data structures" << std::endl;
	auto MBGstructures = translateToMBGDataStructures(unitigGraph, readPaths, graphk);
	MBG::UnitigGraph resolvedGraph;
	std::vector<MBG::ReadPath> resolvedPaths;
	// std::pair<UnitigGraph, std::vector<ReadPath>> resolveUnitigs(const UnitigGraph& initial, const HashList& hashlist, std::vector<ReadPath>& readPaths, const size_t minCoverage, const size_t kmerSize, const size_t maxResolveLength, const size_t maxUnconditionalResolveLength, const bool keepGaps, const bool guesswork, const bool copycountFilterHeuristic, const size_t maxLocalResolve, const bool doCleaning, std::ostream& log);
	std::cerr << "call MBG" << std::endl;
	std::tie(resolvedGraph, resolvedPaths) = MBG::resolveUnitigs(std::get<1>(MBGstructures), std::get<0>(MBGstructures), std::get<2>(MBGstructures), 1, graphk, 100, 0, true, false, false, std::numeric_limits<size_t>::max(), false, std::cerr);
	// std::tie(resolvedGraph, resolvedPaths) = MBG::resolveUnitigs(std::get<1>(MBGstructures), std::get<0>(MBGstructures), std::get<2>(MBGstructures), 1, graphk, 5000, 200, true, false, false, false, false, std::cerr);
	UnitigGraph resultGraph;
	std::vector<ReadPathBundle> resultPaths;
	std::cerr << "switch to AlnDBG data structures" << std::endl;
	std::tie(resultGraph, resultPaths) = translateToAlnDBGDataStructures(resolvedGraph, resolvedPaths, std::get<0>(MBGstructures), graphk, readPaths, unitigGraph);
	std::cerr << "got graph" << std::endl;
	return std::make_pair(resultGraph, resultPaths);
}
