#ifndef Common_h
#define Common_h

#include <cstddef>
#include <cstdint>
#include <vector>
#include <thread>
#include "SparseEdgeContainer.h"
#include "MostlySparse2DHashmap.h"

const uint64_t firstBitUint64_t = 1ull << 63ull;
const uint64_t maskUint64_t = firstBitUint64_t-1;

class ReadPathBundle
{
public:
	class ReadPath
	{
	public:
		std::vector<uint64_t> path;
		size_t pathLeftClipKmers;
		size_t pathRightClipKmers;
		size_t readStartPos;
	};
	std::vector<ReadPath> paths;
	size_t readName;
	size_t readLength;
};

SparseEdgeContainer getActiveEdges(const MostlySparse2DHashmap<uint8_t, size_t>& edges, const size_t numNodes);
size_t intersectSize(const std::vector<size_t>& left, const std::vector<size_t>& right);

template <typename F>
void iterateMultithreaded(size_t start, size_t end, size_t numThreads, F callback)
{
	std::vector<std::thread> threads;
	std::atomic<size_t> nextIndex;
	nextIndex = start;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&nextIndex, start, end, callback]()
		{
			while (true)
			{
				size_t gotIndex = nextIndex++;
				if (gotIndex >= end) break;
				callback(gotIndex);
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
}

#endif
