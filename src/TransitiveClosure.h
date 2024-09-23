#ifndef TransitiveClosure_h
#define TransitiveClosure_h

#include <mutex>
#include <vector>
#include "Common.h"
#include "UnionFind.h"

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
std::vector<size_t> getFastTransitiveClosureMultithread(const size_t itemCount, const size_t maxDistance, const size_t numThreads, F distanceFunction, F2 allowedPairwiseDistanceFunction)
{
	const size_t chunkSize = 100;
	std::vector<size_t> parent;
	parent.resize(itemCount);
	for (size_t i = 0; i < itemCount; i++)
	{
		parent[i] = i;
	}
	std::vector<size_t> clusterExample;
	std::vector<std::vector<size_t>> clusterAdditionals;
	std::vector<size_t> clusterMaxDistance;
	std::mutex resultMutex;
	iterateMultithreaded(0, (itemCount+chunkSize-1)/chunkSize, numThreads, [&resultMutex, &clusterExample, &clusterAdditionals, &clusterMaxDistance, &parent, distanceFunction, allowedPairwiseDistanceFunction, maxDistance, itemCount](const size_t blockindex)
	{
		for (size_t i = blockindex*chunkSize; i < blockindex*chunkSize+chunkSize && i < itemCount; i++)
		{
			bool found = false;
			for (size_t j = 0; ; j++)
			{
				size_t clusterCenterIndex = 0;
				{
					std::lock_guard<std::mutex> lock { resultMutex };
					if (j == clusterExample.size()) break;
					clusterCenterIndex = clusterExample[j];
					if (clusterAdditionals[j].size() > 500) continue;
				}
				size_t allowedDistance = allowedPairwiseDistanceFunction(i, clusterCenterIndex);
				assert(allowedDistance <= maxDistance);
				size_t distance = distanceFunction(i, clusterCenterIndex, allowedDistance);
				if (distance <= allowedDistance)
				{
					std::lock_guard<std::mutex> lock { resultMutex };
					clusterAdditionals[j].emplace_back(i);
					found = true;
					clusterMaxDistance[j] = std::max(distance, clusterMaxDistance[j]);
					merge(parent, i, clusterCenterIndex);
					break;
				}
			}
			if (!found)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				clusterExample.emplace_back(i);
				clusterAdditionals.emplace_back();
				clusterMaxDistance.emplace_back(0);
			}
		}
	});
	assert(clusterExample.size() >= 1);
	std::vector<size_t> clusterOrder;
	for (size_t i = 0; i < clusterExample.size(); i++)
	{
		clusterOrder.emplace_back(i);
	}
	std::sort(clusterOrder.begin(), clusterOrder.end(), [&clusterAdditionals](size_t left, size_t right) { return clusterAdditionals[left].size() > clusterAdditionals[right].size(); });
	size_t maxPairs = (clusterExample.size()-1)*clusterExample.size()/2;
	std::vector<std::thread> threads;
	size_t chunki = 0;
	size_t chunkj = 1;
	for (size_t threadi = 0; threadi < numThreads; threadi++)
	{
		threads.emplace_back([&chunki, &chunkj, &resultMutex, &clusterExample, &clusterAdditionals, &clusterMaxDistance, &clusterOrder, &parent, distanceFunction, allowedPairwiseDistanceFunction, maxDistance, maxPairs]()
		{
			while (true)
			{
				size_t i, j;
				{
					std::lock_guard<std::mutex> lock { resultMutex };
					if (chunki+1 == clusterAdditionals.size()) break;
					assert(chunkj < clusterAdditionals.size());
					i = clusterOrder[chunki];
					j = clusterOrder[chunkj];
					chunkj += 1;
					if (chunkj == clusterAdditionals.size())
					{
						chunki += 1;
						if (chunki+1 == clusterAdditionals.size()) break;
						chunkj = chunki+1;
					}
					if (find(parent, clusterExample[i]) == find(parent, clusterExample[j])) continue;
				}
				size_t distance = distanceFunction(clusterExample[i], clusterExample[j], maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j]);
				if (distance > maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j]) continue;
				if (distance <= allowedPairwiseDistanceFunction(clusterExample[i], clusterExample[j]))
				{
					std::lock_guard<std::mutex> lock { resultMutex };
					merge(parent, clusterExample[i], clusterExample[j]);
					continue;
				}
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
					if (distanceFunction(k, clusterExample[j], allowedDistance) > allowedDistance) continue;
					std::lock_guard<std::mutex> lock { resultMutex };
					merge(parent, k, clusterExample[j]);
					found = true;
					break;
				}
				if (found) continue;
				for (size_t l : clusterAdditionals[j])
				{
					size_t allowedDistance = allowedPairwiseDistanceFunction(l, clusterExample[i]);
					assert(allowedDistance <= maxDistance);
					if (allowedDistance < minPossibleDistance) continue;
					if (distanceFunction(l, clusterExample[i], allowedDistance) > allowedDistance) continue;
					std::lock_guard<std::mutex> lock { resultMutex };
					merge(parent, l, clusterExample[i]);
					found = true;
					break;
				}
				if (found) continue;
				for (size_t k : clusterAdditionals[i])
				{
					for (size_t l : clusterAdditionals[j])
					{
						size_t allowedDistance = allowedPairwiseDistanceFunction(k, l);
						assert(allowedDistance <= maxDistance);
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
	for (size_t threadi = 0; threadi < numThreads; threadi++)
	{
		threads[threadi].join();
	}
	return parent;
}

template <typename F>
std::vector<size_t> getFastTransitiveClosureMultithread(const size_t itemCount, const size_t maxDistance, const size_t numThreads, F distanceFunction)
{
	return getFastTransitiveClosureMultithread(itemCount, maxDistance, numThreads, distanceFunction, [maxDistance](const size_t left, const size_t right) { return maxDistance; });
}

#endif
