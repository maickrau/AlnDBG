#ifndef CanonHelper_h
#define CanonHelper_h

#include <vector>
#include <tuple>
#include <algorithm>
#include "Common.h"

std::vector<bool> getCanonicalChunks(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead);
std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> extrapolateCanonInformation(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerReadBeforeModification, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& modifiedChunksPerRead);

template <typename F>
void iterateCanonicalChunksByCoverage(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t numThreads, F callback)
{
	std::vector<std::vector<std::pair<size_t, size_t>>> occurrencesPerChunk;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			auto t = chunksPerRead[i][j];
			if (NonexistantChunk(std::get<2>(t))) continue;
			if ((std::get<2>(t) & maskUint64_t) >= occurrencesPerChunk.size()) occurrencesPerChunk.resize((std::get<2>(t) & maskUint64_t)+1);
			occurrencesPerChunk[(std::get<2>(t) & maskUint64_t)].emplace_back(i, j);
		}
	}
	std::vector<bool> canonical = getCanonicalChunks(chunksPerRead);
	std::vector<size_t> iterationOrder;
	iterationOrder.reserve(occurrencesPerChunk.size());
	for (size_t i = 0; i < occurrencesPerChunk.size(); i++)
	{
		if (!canonical[i]) continue;
		iterationOrder.emplace_back(i);
	}
	std::sort(iterationOrder.begin(), iterationOrder.end(), [&occurrencesPerChunk](size_t left, size_t right) { return occurrencesPerChunk[left].size() > occurrencesPerChunk[right].size(); });
	iterateMultithreaded(0, iterationOrder.size(), numThreads, [&occurrencesPerChunk, &iterationOrder, callback](const size_t iterationIndex)
	{
		callback(iterationOrder[iterationIndex], occurrencesPerChunk);
	});
}

#endif
