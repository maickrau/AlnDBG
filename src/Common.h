#ifndef Common_h
#define Common_h

#include <cstddef>
#include <cstdint>

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

#endif
