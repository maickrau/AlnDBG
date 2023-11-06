#ifndef Common_h
#define Common_h

#include <cstddef>
#include <cstdint>
#include <vector>
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
};

SparseEdgeContainer getActiveEdges(const MostlySparse2DHashmap<uint8_t, size_t>& edges, const size_t numNodes);
size_t intersectSize(const std::vector<size_t>& left, const std::vector<size_t>& right);

#endif
