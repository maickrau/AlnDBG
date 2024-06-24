#ifndef Common_h
#define Common_h

#include <cstddef>
#include <cstdint>
#include <vector>
#include <thread>
#include "MostlySparse2DHashmap.h"
#include "SparseEdgeContainer.h"

const uint64_t firstBitUint64_t = 1ull << 63ull;
const uint64_t maskUint64_t = firstBitUint64_t-1;

std::pair<size_t, bool> reverse(std::pair<size_t, bool> pos);
std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> canon(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
std::string revCompRaw(const std::string& raw);
bool NonexistantChunk(const uint64_t chunk);
std::chrono::time_point<std::chrono::steady_clock> getTime();
std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end);
extern std::chrono::time_point<std::chrono::steady_clock> programStartTime;

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
