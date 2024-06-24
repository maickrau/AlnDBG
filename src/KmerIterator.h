#ifndef KmerIterator_h
#define KmerIterator_h

#include <vector>
#include <cassert>
#include <tuple>
#include <phmap.h>
#include "TwobitString.h"

template <typename F>
void iterateKmers(const TwobitString& baseSequence, const size_t start, const size_t end, const bool fw, const size_t k, F callback)
{
	const size_t mask = (1ull << (2ull*k))-1;
	size_t kmer = 0;
	if (fw)
	{
		for (size_t m = start; m <= end; m++)
		{
			kmer <<= 2;
			kmer += baseSequence.get(m);
			kmer &= mask;
			if (m < start+k-1) continue;
			callback(kmer, m - (start+k-1));
		}
	}
	else
	{
		for (size_t m = end; m >= start && m < baseSequence.size(); m--)
		{
			kmer <<= 2;
			kmer += 3 - baseSequence.get(m);
			kmer &= mask;
			if (m > end-k+1) continue;
			callback(kmer, end - m + 1 - k);
		}
	}
}

template <typename KmerType, typename F>
phmap::flat_hash_map<size_t, std::vector<std::pair<double, double>>> iterateSolidKmersWithKmerType(const std::vector<TwobitString>& chunkSequences, const size_t kmerSize, const size_t minSolidThreshold, const bool allowTwoAlleleRepeats, const bool needClusters, F callback)
{
	assert(chunkSequences.size() < (size_t)std::numeric_limits<uint32_t>::max());
	assert(chunkSequences.size() >= 1);
	double repetitiveClusteringDistance = 1.0;
	// if one sequence is too small then don't do anything even if theoretically other sequences could have solid kmers
	size_t smallestSequence = chunkSequences[0].size();
	for (size_t i = 0; i < chunkSequences.size(); i++)
	{
		if (chunkSequences[i].size() < kmerSize+2) return phmap::flat_hash_map<size_t, std::vector<std::pair<double, double>>> {};
		smallestSequence = std::min(smallestSequence, chunkSequences[i].size());
	}
	std::vector<size_t> currentPosPerOccurrence;
	// 40 blocks so it advances 2.5% of sequence every block so kmers seen before block cannot cluster with kmers seen after block if no kmers in block
	size_t numBlocks = 40;
	// unless the strings are small, then no blocks
	if (smallestSequence < 2000)
	{
		repetitiveClusteringDistance = 100 * 20.0/(double)smallestSequence;
	}
	if (smallestSequence < 2000 || smallestSequence < 40*kmerSize*2)
	{
		numBlocks = 1;
	}
	for (size_t j = 0; j < chunkSequences.size(); j++)
	{
		assert(chunkSequences[j].size() < (size_t)std::numeric_limits<uint16_t>::max());
		currentPosPerOccurrence.emplace_back(0);
	}
	std::vector<std::vector<std::pair<double, double>>> validClusters;
	std::vector<std::vector<std::tuple<uint16_t, uint16_t, KmerType>>> solidKmers; // pos, cluster, kmer
	solidKmers.resize(chunkSequences.size());
	phmap::flat_hash_map<KmerType, size_t> kmerToIndex;
	std::vector<std::vector<std::tuple<float, uint32_t, uint16_t>>> kmers; // kmer -> (extrapolated-pos, occurrence, occurrence-pos)
	std::vector<size_t> nextCluster;
	std::vector<size_t> kmersDroppingOutOfContext;
	std::vector<KmerType> actualKmer;
	for (size_t offsetBreakpoint = 0; offsetBreakpoint < numBlocks; offsetBreakpoint++)
	{
		phmap::flat_hash_set<size_t> inContext;
		for (size_t j = 0; j < chunkSequences.size(); j++)
		{
			size_t nextOffset = (double)(offsetBreakpoint+1)/(double)numBlocks * (chunkSequences[j].size()-1);
			assert(offsetBreakpoint != (numBlocks-1) || nextOffset == chunkSequences[j].size()-1);
			assert(nextOffset > currentPosPerOccurrence[j] + kmerSize);
			size_t startPos, endPos;
			startPos = currentPosPerOccurrence[j];
			endPos = nextOffset;
			iterateKmers(chunkSequences[j], startPos, endPos, true, kmerSize, [&kmers, &nextCluster, &kmerToIndex, &inContext, &currentPosPerOccurrence, &validClusters, &actualKmer, &chunkSequences, kmerSize, j, offsetBreakpoint, startPos, endPos](const size_t kmer, const size_t pos)
			{
				assert(pos < endPos-startPos);
				assert(pos + currentPosPerOccurrence[j] < chunkSequences[j].size());
				assert((size_t)(pos+currentPosPerOccurrence[j]) < (size_t)std::numeric_limits<uint16_t>::max());
				double extrapolatedPos = 100.0 * (double)(pos+currentPosPerOccurrence[j]+((double)kmerSize/2.0)) / (double)chunkSequences[j].size();
				size_t index = kmerToIndex.size();
				auto found = kmerToIndex.find(kmer);
				if (found != kmerToIndex.end())
				{
					index = found->second;
				}
				else
				{
					kmerToIndex[kmer] = index;
					kmers.emplace_back();
					nextCluster.emplace_back(0);
					validClusters.emplace_back();
					actualKmer.emplace_back(kmer);
				}
				kmers[index].emplace_back(extrapolatedPos, j, (pos+currentPosPerOccurrence[j]));
				inContext.insert(index);
			});
			currentPosPerOccurrence[j] = nextOffset - kmerSize + 2;
		}
		if (offsetBreakpoint == (numBlocks-1))
		{
			kmersDroppingOutOfContext.insert(kmersDroppingOutOfContext.end(), inContext.begin(), inContext.end());
			inContext.clear();
		}
		for (size_t checkindex = kmersDroppingOutOfContext.size()-1; checkindex < kmersDroppingOutOfContext.size(); checkindex--)
		{
			size_t index = kmersDroppingOutOfContext[checkindex];
			assert(index < kmers.size());
			if (inContext.count(index) == 1) continue;
			if (kmers[index].size() >= minSolidThreshold)
			{
				std::sort(kmers[index].begin(), kmers[index].end());
				size_t clusterStart = 0;
				std::vector<bool> occurrencesHere;
				occurrencesHere.resize(chunkSequences.size(), false);
				bool currentlyNonrepetitive = true;
				occurrencesHere[std::get<1>(kmers[index][0])] = true;
				size_t clusterNum = nextCluster[index];
				for (size_t j = 1; j <= kmers[index].size(); j++)
				{
					if (j < kmers[index].size() && std::get<0>(kmers[index][j]) < std::get<0>(kmers[index][j-1])+repetitiveClusteringDistance)
					{
						if (occurrencesHere[std::get<1>(kmers[index][j])])
						{
							currentlyNonrepetitive = false;
						}
						occurrencesHere[std::get<1>(kmers[index][j])] = true;
						continue;
					}
					if (j-clusterStart >= minSolidThreshold && (currentlyNonrepetitive || allowTwoAlleleRepeats))
					{
						double minPos = std::get<0>(kmers[index][clusterStart]);
						double maxPos = std::get<0>(kmers[index][j-1]);
						if (minPos + repetitiveClusteringDistance > maxPos)
						{
							if (!currentlyNonrepetitive && allowTwoAlleleRepeats)
							{
								size_t countOccurrences = 0;
								for (size_t k = 0; k < occurrencesHere.size(); k++)
								{
									countOccurrences += occurrencesHere[k] ? 1 : 0;
								}
								if (countOccurrences >= minSolidThreshold)
								{
									phmap::flat_hash_map<size_t, size_t> kmerCountPerOccurrence;
									for (size_t k = clusterStart; k < j; k++)
									{
										kmerCountPerOccurrence[std::get<1>(kmers[index][k])] += 1;
									}
									phmap::flat_hash_set<size_t> foundCounts;
									for (auto pair2 : kmerCountPerOccurrence)
									{
										foundCounts.insert(pair2.second);
									}
									if (foundCounts.size() == 2)
									{
										size_t firstCount = *foundCounts.begin();
										size_t secondCount = *(++foundCounts.begin());
										assert(firstCount != secondCount);
										for (auto pair2 : kmerCountPerOccurrence)
										{
											if (pair2.second == firstCount)
											{
												// todo pos should be something reasonable not this. but nothing uses poses of repetitive kmers as of comment time
												assert(clusterNum < (size_t)std::numeric_limits<uint16_t>::max());
												solidKmers[pair2.first].emplace_back(std::numeric_limits<uint16_t>::max(), clusterNum, index);
											}
											else
											{
												assert(pair2.second == secondCount);
												// todo pos should be something reasonable not this. but nothing uses poses of repetitive kmers as of comment time
												assert(clusterNum+1 < (size_t)std::numeric_limits<uint16_t>::max());
												solidKmers[pair2.first].emplace_back(std::numeric_limits<uint16_t>::max(), clusterNum+1, index);
											}
										}
										clusterNum += 2;
										validClusters[index].emplace_back(minPos-0.01, maxPos+0.01);
										validClusters[index].emplace_back(minPos-0.01, maxPos+0.01);
									}
								}
							}
							else
							{
								for (size_t k = clusterStart; k < j; k++)
								{
									assert((size_t)std::get<2>(kmers[index][k]) < (size_t)std::numeric_limits<uint16_t>::max());
									assert(std::get<0>(kmers[index][k]) >= minPos);
									assert(std::get<0>(kmers[index][k]) <= maxPos);
									assert(clusterNum < (size_t)std::numeric_limits<uint16_t>::max());
									solidKmers[std::get<1>(kmers[index][k])].emplace_back(std::get<2>(kmers[index][k]), clusterNum, index);
								}
								clusterNum += 1;
								validClusters[index].emplace_back(minPos-0.01, maxPos+0.01);
							}
						}
					}
					currentlyNonrepetitive = true;
					occurrencesHere.assign(occurrencesHere.size(), false);
					if (j < kmers[index].size())
					{
						occurrencesHere[std::get<1>(kmers[index][j])] = true;
					}
					clusterStart = j;
				}
				nextCluster[index] = clusterNum;
			}
			std::swap(kmersDroppingOutOfContext[checkindex], kmersDroppingOutOfContext.back());
			kmersDroppingOutOfContext.pop_back();
			{
				std::remove_reference<decltype(kmers[index])>::type tmp;
				std::swap(tmp, kmers[index]);
			}
		}
		kmersDroppingOutOfContext.insert(kmersDroppingOutOfContext.end(), inContext.begin(), inContext.end());
	}
	for (size_t i = 0; i < kmers.size(); i++)
	{
		assert(kmers[i].size() == 0);
	}
	phmap::flat_hash_map<size_t, std::vector<std::pair<double, double>>> fixedClusters;
	if (needClusters)
	{
		for (size_t i = 0; i < validClusters.size(); i++)
		{
			fixedClusters[actualKmer[i]];
			std::swap(fixedClusters[actualKmer[i]], validClusters[i]);
		}
	}
	for (size_t i = 0; i < solidKmers.size(); i++)
	{
		std::sort(solidKmers[i].begin(), solidKmers[i].end());
		for (auto kmer : solidKmers[i])
		{
			// occurrence ID, chunk startpos in read, chunk endpos in read, chunk, kmer, clusternum, kmerpos in chunk
			callback(i, 0, chunkSequences[i].size(), std::numeric_limits<size_t>::max(), (size_t)actualKmer[std::get<2>(kmer)], (size_t)std::get<1>(kmer), (size_t)std::get<0>(kmer));
		}
	}
	return fixedClusters;
}

template <typename F>
phmap::flat_hash_map<size_t, std::vector<std::pair<double, double>>> iterateSolidKmers(const std::vector<TwobitString>& chunkSequences, const size_t kmerSize, const size_t minSolidThreshold, const bool allowTwoAlleleRepeats, const bool needClusters, F callback)
{
	if (kmerSize < 16)
	{
		return iterateSolidKmersWithKmerType<uint32_t>(chunkSequences, kmerSize, minSolidThreshold, allowTwoAlleleRepeats, needClusters, callback);
	}
	else
	{
		return iterateSolidKmersWithKmerType<size_t>(chunkSequences, kmerSize, minSolidThreshold, allowTwoAlleleRepeats, needClusters, callback);
	}
}

#endif
