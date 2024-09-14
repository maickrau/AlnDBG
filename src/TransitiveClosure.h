#ifndef TransitiveClosure_h
#define TransitiveClosure_h

#include <mutex>
#include <vector>
#include "Common.h"
#include "UnionFind.h"

template <typename F>
std::vector<size_t> getFastTransitiveClosure(const size_t itemCount, const size_t maxDistance, F distanceFunction)
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
			size_t distance = distanceFunction(i, clusterExample[j], maxDistance);
			if (distance <= maxDistance)
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
			if (distance <= maxDistance)
			{
				merge(parent, clusterExample[i], clusterExample[j]);
				continue;
			}
			else if (distance <= maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j])
			{
				bool found = false;
				for (size_t k : clusterAdditionals[i])
				{
					if (distanceFunction(k, clusterExample[j], maxDistance) <= maxDistance)
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
						if (distanceFunction(l, clusterExample[i], maxDistance) <= maxDistance)
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
							if (distanceFunction(k, l, maxDistance) <= maxDistance)
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
std::vector<size_t> getFastTransitiveClosureMultithread(const size_t itemCount, const size_t maxDistance, const size_t numThreads, F distanceFunction)
{
	const size_t chunkSize = 500;
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
	iterateMultithreaded(0, (itemCount+chunkSize-1)/chunkSize, numThreads, [&resultMutex, &clusterExample, &clusterAdditionals, &clusterMaxDistance, &parent, distanceFunction, maxDistance, itemCount](const size_t blockindex)
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
				}
				size_t distance = distanceFunction(i, clusterCenterIndex, maxDistance);
				if (distance <= maxDistance)
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
	size_t maxPairs = (clusterExample.size()-1)*clusterExample.size()/2;
	iterateMultithreaded(0, (maxPairs+chunkSize-1)/chunkSize, numThreads, [&resultMutex, &clusterExample, &clusterAdditionals, &clusterMaxDistance, &parent, distanceFunction, maxDistance, maxPairs](const size_t blockindex)
	{
		for (size_t index = blockindex * chunkSize; index < blockindex*chunkSize+chunkSize && index < maxPairs; index++)
		{
			size_t i = ceil(sqrt(1 + 8 * (index+1))/2.0 + 0.5) - 1;
			assert(index >= (i-1)*i/2);
			size_t j = index - (i-1)*i/2;
			assert(j < i);
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				// std::cerr << "check i " << i << " j " << j << " / " << clusterExample.size() << std::endl;
				if (find(parent, clusterExample[i]) == find(parent, clusterExample[j])) continue;
			}
			size_t distance = distanceFunction(clusterExample[i], clusterExample[j], maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j]);
			if (distance > maxDistance + clusterMaxDistance[i] + clusterMaxDistance[j]) continue;
			if (distance <= maxDistance)
			{
				std::lock_guard<std::mutex> lock { resultMutex };
				merge(parent, clusterExample[i], clusterExample[j]);
				continue;
			}
			bool found = false;
			for (size_t k : clusterAdditionals[i])
			{
				if (distanceFunction(k, clusterExample[j], maxDistance) > maxDistance) continue;
				std::lock_guard<std::mutex> lock { resultMutex };
				merge(parent, k, clusterExample[j]);
				found = true;
				break;
			}
			if (found) continue;
			for (size_t l : clusterAdditionals[j])
			{
				if (distanceFunction(l, clusterExample[i], maxDistance) > maxDistance) continue;
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
					if (distanceFunction(k, l, maxDistance) > maxDistance) continue;
					std::lock_guard<std::mutex> lock { resultMutex };
					merge(parent, k, l);
					found = true;
					break;
				}
				if (found) break;
			}
		}
	});
	return parent;
}

#endif
