#ifndef TransitiveClosure_h
#define TransitiveClosure_h

#include <mutex>
#include <vector>
#include "Common.h"
#include "UnionFind.h"

std::vector<std::vector<int>> getBinCheckOrderNotSymmetric(const int maxDistance);
std::vector<std::vector<int>> getBinCheckOrder(const int maxDistance);

template <typename F, typename F2>
std::vector<size_t> getFastTransitiveClosure(const size_t itemCount, const size_t maxDistance, F distanceFunction, F2 allowedPairwiseDistanceFunction)
{
	std::vector<size_t> parent;
	parent.resize(itemCount);
	for (size_t i = 0; i < itemCount; i++)
	{
		parent[i] = i;
	}
	std::vector<size_t> clusterExample;
	std::vector<std::vector<size_t>> clusterAdditionals;
	std::vector<size_t> clusterMaxDistance;
	for (size_t i = 0; i < itemCount; i++)
	{
		bool found = false;
		for (size_t j = clusterExample.size()-1; j < clusterExample.size(); j--)
		{
			size_t allowedDistance = allowedPairwiseDistanceFunction(i, clusterExample[j]);
			assert(allowedDistance <= maxDistance);
			size_t distance = distanceFunction(i, clusterExample[j], allowedDistance);
			if (distance <= allowedDistance)
			{
				clusterAdditionals[j].emplace_back(i);
				found = true;
				clusterMaxDistance[j] = std::max(distance, clusterMaxDistance[j]);
				merge(parent, i, clusterExample[j]);
				break;
			}
		}
		if (!found)
		{
			clusterExample.emplace_back(i);
			clusterAdditionals.emplace_back();
			clusterMaxDistance.emplace_back(0);
		}
	}
	for (size_t i = 1; i < clusterExample.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (find(parent, clusterExample[i]) == find(parent, clusterExample[j])) continue;
			size_t distance = distanceFunction(clusterExample[i], clusterExample[j], maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j]);
			if (distance <= allowedPairwiseDistanceFunction(clusterExample[i], clusterExample[j]))
			{
				merge(parent, clusterExample[i], clusterExample[j]);
				continue;
			}
			else if (distance <= maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j])
			{
				size_t minPossibleDistance = 0;
				if (distance >= clusterMaxDistance[i] + clusterMaxDistance[j])
				{
					minPossibleDistance = distance - clusterMaxDistance[i] - clusterMaxDistance[j];
				}
				bool found = false;
				for (size_t k : clusterAdditionals[i])
				{
					size_t allowedDistance = allowedPairwiseDistanceFunction(k, clusterExample[j]);
					assert(allowedDistance <= maxDistance);
					if (allowedDistance < minPossibleDistance) continue;
					if (distanceFunction(k, clusterExample[j], allowedDistance) <= allowedDistance)
					{
						found = true;
						merge(parent, k, clusterExample[j]);
						break;
					}
				}
				if (!found)
				{
					for (size_t l : clusterAdditionals[j])
					{
						size_t allowedDistance = allowedPairwiseDistanceFunction(l, clusterExample[i]);
						assert(allowedDistance <= maxDistance);
						if (allowedDistance < minPossibleDistance) continue;
						if (distanceFunction(l, clusterExample[i], allowedDistance) <= allowedDistance)
						{
							found = true;
							merge(parent, l, clusterExample[i]);
							break;
						}
					}
				}
				if (!found)
				{
					for (size_t k : clusterAdditionals[i])
					{
						for (size_t l : clusterAdditionals[j])
						{
							size_t allowedDistance = allowedPairwiseDistanceFunction(k, l);
							assert(allowedDistance <= maxDistance);
							if (allowedDistance < minPossibleDistance) continue;
							if (distanceFunction(k, l, allowedDistance) <= allowedDistance)
							{
								found = true;
								merge(parent, k, l);
								break;
							}
						}
						if (found) break;
					}
				}
			}
		}
	}
	return parent;
}

template <typename F>
std::vector<size_t> getFastTransitiveClosure(const size_t itemCount, const size_t maxDistance, F distanceFunction)
{
	return getFastTransitiveClosure(itemCount, maxDistance, distanceFunction, [maxDistance](const size_t left, const size_t right) { return maxDistance; });
}

template <typename F, typename F2>
std::vector<size_t> getFastTransitiveClosureMultithread(const size_t itemCount, const size_t maxAllowedDistance, const size_t maxDistanceEver, const size_t numThreads, F distanceFunction, F2 allowedPairwiseDistanceFunction)
{
	const size_t maxClusterSize = 100;
	const size_t chunkSize = 100;
//	static const std::vector<std::vector<int>> binCheckOrder = getBinCheckOrder(1);
//	assert(binCheckOrder.size() == 27);
	static const std::vector<std::vector<int>> binCheckOrder = std::vector<std::vector<int>>{{0, 0, 0}};
	assert(binCheckOrder.size() == 1);
//	static const std::vector<std::vector<int>> binCheckOrderDistanceTwoNotsymmetric = getBinCheckOrderNotSymmetric(2);
//	assert(binCheckOrderDistanceTwoNotsymmetric.size() == 1*1*3 + 1*2*5 + 2*5*5);
	static const std::vector<std::vector<int>> binCheckOrderDistanceTwoNotsymmetric = getBinCheckOrderNotSymmetric(1);
	assert(binCheckOrderDistanceTwoNotsymmetric.size() == 1*1*2 + 1*1*3 + 1*3*3);
	std::vector<size_t> parent;
	parent.resize(itemCount);
	for (size_t i = 0; i < itemCount; i++)
	{
		parent[i] = i;
	}
	std::vector<size_t> anchors;
	anchors.emplace_back(0);
	anchors.emplace_back(itemCount/2);
	anchors.emplace_back(itemCount-1);
	std::vector<size_t> clusterExample;
	std::vector<std::vector<size_t>> clusterAdditionals;
	std::vector<size_t> clusterMaxDistance;
	clusterExample.emplace_back(anchors[0]);
	clusterExample.emplace_back(anchors[1]);
	clusterExample.emplace_back(anchors[2]);
	clusterAdditionals.emplace_back();
	clusterAdditionals.emplace_back();
	clusterAdditionals.emplace_back();
	clusterMaxDistance.emplace_back(0);
	clusterMaxDistance.emplace_back(0);
	clusterMaxDistance.emplace_back(0);
	std::vector<std::vector<std::vector<std::vector<size_t>>>> clusterBins;
	size_t clusteringDistance = maxAllowedDistance+1;
	if (clusteringDistance < maxDistanceEver*0.02) clusteringDistance = maxDistanceEver * 0.02;
	clusterBins.resize((maxDistanceEver+clusteringDistance-1)/clusteringDistance);
	for (size_t i = 0; i < clusterBins.size(); i++)
	{
		clusterBins[i].resize(clusterBins.size());
		for (size_t j = 0; j < clusterBins[i].size(); j++)
		{
			clusterBins[i][j].resize(clusterBins.size());
		}
	}
	if (distanceFunction(anchors[0], anchors[1], allowedPairwiseDistanceFunction(anchors[0], anchors[1])) <= allowedPairwiseDistanceFunction(anchors[0], anchors[1]))
	{
		merge(parent, anchors[0], anchors[1]);
	}
	if (distanceFunction(anchors[0], anchors[2], allowedPairwiseDistanceFunction(anchors[0], anchors[2])) <= allowedPairwiseDistanceFunction(anchors[0], anchors[2]))
	{
		merge(parent, anchors[0], anchors[2]);
	}
	if (distanceFunction(anchors[1], anchors[2], allowedPairwiseDistanceFunction(anchors[1], anchors[2])) <= allowedPairwiseDistanceFunction(anchors[1], anchors[2]))
	{
		merge(parent, anchors[1], anchors[2]);
	}
	std::mutex resultMutex;
	iterateMultithreaded(0, (itemCount+chunkSize-1)/chunkSize, numThreads, [&resultMutex, &clusterExample, &clusterBins, &clusterAdditionals, &clusterMaxDistance, &parent, &anchors, distanceFunction, allowedPairwiseDistanceFunction, maxAllowedDistance, clusteringDistance, maxDistanceEver, itemCount](const size_t blockindex)
	{
		for (size_t i = blockindex*chunkSize; i < blockindex*chunkSize+chunkSize && i < itemCount; i++)
		{
			if (i == anchors[0]) continue;
			if (i == anchors[1]) continue;
			if (i == anchors[2]) continue;
			bool found = false;
			std::vector<size_t> binIndices;
			for (size_t j = 0; j < anchors.size(); j++)
			{
				size_t distance = distanceFunction(i, anchors[j], maxDistanceEver);
				if (distance <= allowedPairwiseDistanceFunction(i, anchors[j]))
				{
					std::lock_guard<std::mutex> lock { resultMutex };
					if (clusterAdditionals[j].size() < maxClusterSize)
					{
						found = true;
						if (distance > 0)
						{
							clusterAdditionals[j].emplace_back(i);
							clusterMaxDistance[j] = std::max(distance, clusterMaxDistance[j]);
						}
						merge(parent, i, anchors[j]);
						break;
					}
				}
				binIndices.emplace_back(distance / clusteringDistance);
			}
			if (found) continue;
			for (const std::vector<int>& offset : binCheckOrder)
			{
				int dim1 = (int)binIndices[0] + (int)offset[0];
				int dim2 = (int)binIndices[1] + (int)offset[1];
				int dim3 = (int)binIndices[2] + (int)offset[2];
				if (dim1 < 0 || dim1 >= (int)clusterBins.size()) continue;
				if (dim2 < 0 || dim2 >= (int)clusterBins.size()) continue;
				if (dim3 < 0 || dim3 >= (int)clusterBins.size()) continue;
				assert(dim1 < (int)clusterBins.size());
				assert(dim2 < (int)clusterBins[dim1].size());
				assert(dim3 < (int)clusterBins[dim1][dim2].size());
				size_t clusterIndex = 0;
				while (true)
				{
					size_t cluster = std::numeric_limits<size_t>::max();
					size_t clusterCenterIndex = 0;
					{
						std::lock_guard<std::mutex> lock { resultMutex };
						if (clusterIndex < clusterBins[dim1][dim2][dim3].size())
						{
							cluster = clusterBins[dim1][dim2][dim3][clusterIndex];
							clusterCenterIndex = clusterExample[cluster];
						}
						else
						{
							break;
						}
						clusterIndex += 1;
						if (clusterAdditionals[cluster].size() > maxClusterSize) continue;
					}
					if (cluster == std::numeric_limits<size_t>::max()) break;
					size_t allowedDistance = allowedPairwiseDistanceFunction(i, clusterCenterIndex);
					size_t distance = distanceFunction(i, clusterCenterIndex, allowedDistance);
					if (distance > allowedDistance) continue;
					found = true;
					std::lock_guard<std::mutex> lock { resultMutex };
					if (distance > 0)
					{
						clusterAdditionals[cluster].emplace_back(i);
						clusterMaxDistance[cluster] = std::max(distance, clusterMaxDistance[cluster]);
					}
					merge(parent, i, clusterCenterIndex);
					break;
				}
				if (found) break;
			}
			if (found) continue;
			std::lock_guard<std::mutex> lock { resultMutex };
			size_t cluster = clusterExample.size();
			clusterExample.emplace_back(i);
			clusterAdditionals.emplace_back();
			clusterMaxDistance.emplace_back(0);
			clusterBins[binIndices[0]][binIndices[1]][binIndices[2]].emplace_back(cluster);
		}
	});
	assert(clusterExample.size() >= 1);
	std::vector<size_t> remainingCheckClusterParent;
	for (size_t i = 0; i < clusterExample.size(); i++)
	{
		remainingCheckClusterParent.emplace_back(i);
	}
	iterateMultithreaded(1, clusterExample.size(), numThreads, [&clusterExample, &anchors, &clusterMaxDistance, &resultMutex, &remainingCheckClusterParent, &parent, distanceFunction, allowedPairwiseDistanceFunction, maxAllowedDistance](const size_t i)
	{
		for (size_t j = 0; j < anchors.size() && j < i; j++)
		{
			size_t distance = distanceFunction(clusterExample[j], clusterExample[i], maxAllowedDistance + clusterMaxDistance[j] + clusterMaxDistance[i]);
			if (distance <= allowedPairwiseDistanceFunction(clusterExample[j], clusterExample[i]))
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				merge(parent, clusterExample[j], clusterExample[i]);
			}
			else if (distance <= maxAllowedDistance + clusterMaxDistance[j] + clusterMaxDistance[i])
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				merge(remainingCheckClusterParent, i, j);
			}
		}
	});
	for (const std::vector<int>& offset : binCheckOrderDistanceTwoNotsymmetric)
	{
		for (size_t clusterDimOne = 0; clusterDimOne < clusterBins.size(); clusterDimOne++)
		{
			int otherDimOne = (int)clusterDimOne + offset[0];
			if (otherDimOne < 0) continue;
			if (otherDimOne >= (int)clusterBins.size()) continue;
			for (size_t clusterDimTwo = 0; clusterDimTwo < clusterBins.size(); clusterDimTwo++)
			{
				int otherDimTwo = (int)clusterDimTwo + offset[1];
				if (otherDimTwo >= (int)clusterBins.size()) continue;
				if (otherDimTwo < 0) continue;
				for (size_t clusterDimThree = 0; clusterDimThree < clusterBins.size(); clusterDimThree++)
				{
					int otherDimThree = (int)clusterDimThree + offset[2];
					if (otherDimThree >= (int)clusterBins.size()) continue;
					if (otherDimThree < 0) continue;
					if (clusterBins[clusterDimOne][clusterDimTwo][clusterDimThree].size() == 0) continue;
					if (clusterBins[otherDimOne][otherDimTwo][otherDimThree].size() == 0) continue;
					size_t thisIndex = clusterDimOne * clusterBins.size() * clusterBins.size() + clusterDimTwo * clusterBins.size() + clusterDimThree;
					size_t otherIndex = otherDimOne * clusterBins.size() * clusterBins.size() + otherDimTwo * clusterBins.size() + otherDimThree;
					if (otherIndex < thisIndex) continue;
					iterateMultithreaded(0, clusterBins[clusterDimOne][clusterDimTwo][clusterDimThree].size(), numThreads, [&clusterExample, &parent, &resultMutex, &offset, &clusterBins, &clusterMaxDistance, &remainingCheckClusterParent, maxAllowedDistance, clusterDimOne, clusterDimTwo, clusterDimThree, otherDimOne, otherDimTwo, otherDimThree, allowedPairwiseDistanceFunction, distanceFunction](const size_t i)
					{
						size_t clusteri = clusterBins[clusterDimOne][clusterDimTwo][clusterDimThree][i];
						for (size_t j = 0; j < clusterBins[otherDimOne][otherDimTwo][otherDimThree].size(); j++)
						{
							if (offset[0] == 0 && offset[1] == 0 && offset[2] == 0 && j == i) break;
							size_t clusterj = clusterBins[otherDimOne][otherDimTwo][otherDimThree][j];
							{
								std::lock_guard<std::mutex> lock { resultMutex };
								if (find(parent, clusterExample[clusteri]) == find(parent, clusterExample[clusterj])) continue;
							}
							size_t distance = distanceFunction(clusterExample[clusteri], clusterExample[clusterj], maxAllowedDistance + clusterMaxDistance[clusteri] + clusterMaxDistance[clusterj]);
							size_t allowedDistance = allowedPairwiseDistanceFunction(clusterExample[clusteri], clusterExample[clusterj]);
							if (distance <= allowedDistance)
							{
								std::lock_guard<std::mutex> lock { resultMutex };
								merge(parent, clusterExample[clusteri], clusterExample[clusterj]);
							}
							else if (distance <= maxAllowedDistance + clusterMaxDistance[clusteri] + clusterMaxDistance[clusterj])
							{
								std::lock_guard<std::mutex> lock { resultMutex };
								merge(remainingCheckClusterParent, clusteri, clusterj);
							}
						}
					});
				}
			}
		}
	}
	std::vector<std::vector<size_t>> remainingCheckClusters;
	phmap::flat_hash_map<size_t, size_t> clusterToId;
	for (size_t i = 0; i < remainingCheckClusterParent.size(); i++)
	{
		auto found = find(remainingCheckClusterParent, i);
		if (clusterToId.count(found) == 0)
		{
			clusterToId[found] = remainingCheckClusters.size();
			remainingCheckClusters.emplace_back();
		}
		remainingCheckClusters[clusterToId.at(found)].emplace_back(i);
	}
	for (size_t bigclusteri = 0; bigclusteri < remainingCheckClusters.size(); bigclusteri++)
	{
		assert(remainingCheckClusters[bigclusteri].size() >= 1);
		if (remainingCheckClusters[bigclusteri].size() == 1) continue;
		std::mutex resultMutex;
		size_t nexti = 1;
		size_t nextj = 0;
		iterateMultithreaded(0, remainingCheckClusters[bigclusteri].size(), numThreads, [&resultMutex, &nexti, &nextj, &parent, &remainingCheckClusters, &clusterExample, &clusterAdditionals, &clusterMaxDistance, maxAllowedDistance, distanceFunction, allowedPairwiseDistanceFunction, bigclusteri](const size_t dummyThreadId)
		{
			while (true)
			{
				size_t pickedI = 0;
				size_t pickedJ = 0;
				{
					std::lock_guard<std::mutex> lock { resultMutex };
					if (nexti == remainingCheckClusters[bigclusteri].size()) break;
					pickedI = remainingCheckClusters[bigclusteri][nexti];
					pickedJ = remainingCheckClusters[bigclusteri][nextj];
					assert(nexti < remainingCheckClusters[bigclusteri].size());
					assert(nextj < nexti);
					nextj += 1;
					if (nextj == nexti)
					{
						nexti += 1;
						nextj = 0;
					}
					if (find(parent, clusterExample[pickedI]) == find(parent, clusterExample[pickedJ])) continue;
				}
				size_t distance = distanceFunction(clusterExample[pickedI], clusterExample[pickedJ], maxAllowedDistance + clusterMaxDistance[pickedI] + clusterMaxDistance[pickedJ]);
				if (distance > maxAllowedDistance + clusterMaxDistance[pickedI] + clusterMaxDistance[pickedJ]) continue;
				size_t minPossibleDistance = 0;
				if (distance >= clusterMaxDistance[pickedI] + clusterMaxDistance[pickedJ])
				{
					minPossibleDistance = distance - clusterMaxDistance[pickedI] - clusterMaxDistance[pickedJ];
				}
				bool found = false;
				std::vector<size_t> additionalsInI;
				for (size_t k : clusterAdditionals[pickedI])
				{
					size_t allowedDistance = allowedPairwiseDistanceFunction(k, clusterExample[pickedJ]);
					assert(allowedDistance <= maxAllowedDistance);
					if (allowedDistance < minPossibleDistance) continue;
					size_t maxDistance = clusterMaxDistance[pickedJ] + maxAllowedDistance;
					size_t distance = distanceFunction(k, clusterExample[pickedJ], maxDistance);
					if (distance <= allowedDistance)
					{
						std::lock_guard<std::mutex> lock { resultMutex };
						merge(parent, k, clusterExample[pickedJ]);
						found = true;
						break;
					}
					if (distance <= maxDistance)
					{
						additionalsInI.push_back(k);
					}
				}
				if (!found)
				{
					// check in case another thread transitively merged these
					std::lock_guard<std::mutex> lock { resultMutex };
					if (find(parent, clusterExample[pickedI]) == find(parent, clusterExample[pickedJ])) found = true;
				}
				if (found) continue;
				std::vector<size_t> additionalsInJ;
				for (size_t l : clusterAdditionals[pickedJ])
				{
					size_t allowedDistance = allowedPairwiseDistanceFunction(l, clusterExample[pickedI]);
					assert(allowedDistance <= maxAllowedDistance);
					if (allowedDistance < minPossibleDistance) continue;
					size_t maxDistance = clusterMaxDistance[pickedI] + maxAllowedDistance;
					size_t distance = distanceFunction(l, clusterExample[pickedI], maxDistance);
					if (distance <= allowedDistance)
					{
						std::lock_guard<std::mutex> lock { resultMutex };
						merge(parent, l, clusterExample[pickedI]);
						found = true;
						break;
					}
					if (distance <= maxDistance)
					{
						additionalsInJ.push_back(l);
					}
				}
				if (!found)
				{
					// check in case another thread transitively merged these
					std::lock_guard<std::mutex> lock { resultMutex };
					if (find(parent, clusterExample[pickedI]) == find(parent, clusterExample[pickedJ])) found = true;
				}
				if (found) continue;
				for (size_t k : additionalsInI)
				{
					for (size_t l : additionalsInJ)
					{
						size_t allowedDistance = allowedPairwiseDistanceFunction(k, l);
						assert(allowedDistance <= maxAllowedDistance);
						if (allowedDistance < minPossibleDistance) continue;
						if (distanceFunction(k, l, allowedDistance) > allowedDistance) continue;
						std::lock_guard<std::mutex> lock { resultMutex };
						merge(parent, k, l);
						found = true;
						break;
					}
					if (found) break;
				}
			}
		});
	}
	return parent;
}

template <typename F>
std::vector<size_t> getFastTransitiveClosureMultithread(const size_t itemCount, const size_t maxDistance, const size_t numThreads, F distanceFunction)
{
	return getFastTransitiveClosureMultithread(itemCount, maxDistance, numThreads, distanceFunction, [maxDistance](const size_t left, const size_t right) { return maxDistance; });
}

#endif
