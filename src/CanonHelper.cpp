#include "CanonHelper.h"

std::vector<bool> getCanonicalChunks(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead)
{
	std::vector<bool> canonical;
	assert((chunksPerRead.size() % 2) == 0);
	for (size_t i = 0; i < chunksPerRead.size()/2; i++)
	{
		size_t otheri = i + chunksPerRead.size()/2;
		assert(chunksPerRead[i].size() == chunksPerRead[otheri].size());
		for (size_t j = 0; j < chunksPerRead[i].size(); j++)
		{
			size_t otherj = chunksPerRead[otheri].size()-1-j;
			assert(std::get<1>(chunksPerRead[i][j]) - std::get<0>(chunksPerRead[i][j]) == std::get<1>(chunksPerRead[otheri][otherj]) - std::get<0>(chunksPerRead[otheri][otherj]));
			assert(NonexistantChunk(std::get<2>(chunksPerRead[i][j])) == NonexistantChunk(std::get<2>(chunksPerRead[otheri][otherj])));
			if (NonexistantChunk(std::get<2>(chunksPerRead[i][j]))) continue;
			while (canonical.size() <= (std::get<2>(chunksPerRead[i][j]) & maskUint64_t))
			{
				canonical.emplace_back(true);
			}
			while (canonical.size() <= (std::get<2>(chunksPerRead[otheri][otherj]) & maskUint64_t))
			{
				canonical.emplace_back(true);
			}
			if ((std::get<2>(chunksPerRead[i][j]) & maskUint64_t) > (std::get<2>(chunksPerRead[otheri][otherj]) & maskUint64_t))
			{
				canonical[std::get<2>(chunksPerRead[i][j]) & maskUint64_t] = false;
			}
			if ((std::get<2>(chunksPerRead[i][j]) & maskUint64_t) < (std::get<2>(chunksPerRead[otheri][otherj]) & maskUint64_t))
			{
				canonical[std::get<2>(chunksPerRead[otheri][otherj]) & maskUint64_t] = false;
			}
		}
	}
	return canonical;
}

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> extrapolateCanonInformation(const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerReadBeforeModification, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& modifiedChunksPerRead)
{
	std::vector<bool> canonical = getCanonicalChunks(chunksPerReadBeforeModification);
	phmap::flat_hash_map<size_t, size_t> newChunkToExtrapolatedIndex;
	size_t nextNum = 0;
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> result;
	result = chunksPerReadBeforeModification;
	for (size_t i = 0; i < chunksPerReadBeforeModification.size(); i++)
	{
		assert(modifiedChunksPerRead[i].size() == chunksPerReadBeforeModification[i].size());
		size_t otheri;
		if (i < chunksPerReadBeforeModification.size()/2)
		{
			otheri = i + chunksPerReadBeforeModification.size()/2;
		}
		else
		{
			otheri = i - chunksPerReadBeforeModification.size()/2;
		}
		assert(modifiedChunksPerRead[i].size() == modifiedChunksPerRead[otheri].size());
		assert(modifiedChunksPerRead[otheri].size() == chunksPerReadBeforeModification[otheri].size());
		for (size_t j = 0; j < chunksPerReadBeforeModification[i].size(); j++)
		{
			assert(std::get<0>(chunksPerReadBeforeModification[i][j]) == std::get<0>(modifiedChunksPerRead[i][j]));
			assert(std::get<1>(chunksPerReadBeforeModification[i][j]) == std::get<1>(modifiedChunksPerRead[i][j]));
			size_t otherj = chunksPerReadBeforeModification[otheri].size()-1-j;
			assert(std::get<0>(chunksPerReadBeforeModification[otheri][otherj]) == std::get<0>(modifiedChunksPerRead[otheri][otherj]));
			assert(std::get<1>(chunksPerReadBeforeModification[otheri][otherj]) == std::get<1>(modifiedChunksPerRead[otheri][otherj]));
			if (NonexistantChunk(std::get<2>(chunksPerReadBeforeModification[i][j])))
			{
				assert(NonexistantChunk(std::get<2>(modifiedChunksPerRead[i][j])));
				continue;
			}
			if (NonexistantChunk(std::get<2>(modifiedChunksPerRead[i][j])))
			{
				std::get<2>(result[i][j]) = std::numeric_limits<size_t>::max();
				continue;
			}
			if (!canonical[std::get<2>(chunksPerReadBeforeModification[i][j]) & maskUint64_t])
			{
				assert(canonical[std::get<2>(chunksPerReadBeforeModification[otheri][otherj]) & maskUint64_t]);
				if (newChunkToExtrapolatedIndex.count(std::get<2>(modifiedChunksPerRead[otheri][otherj]) & maskUint64_t) == 1)
				{
					std::get<2>(result[i][j]) = newChunkToExtrapolatedIndex.at(std::get<2>(modifiedChunksPerRead[otheri][otherj]) & maskUint64_t)+1 + firstBitUint64_t;
				}
				else
				{
					size_t index = nextNum;
					newChunkToExtrapolatedIndex[std::get<2>(modifiedChunksPerRead[otheri][otherj]) & maskUint64_t] = index;
					nextNum += 2;
					std::get<2>(result[i][j]) = index+1 + firstBitUint64_t;
				}
				continue;
			}
			if (std::get<2>(chunksPerReadBeforeModification[i][j]) == std::get<2>(chunksPerReadBeforeModification[otheri][otherj]))
			{
				// palindrome
				if (newChunkToExtrapolatedIndex.count(std::get<2>(modifiedChunksPerRead[i][j]) & maskUint64_t) == 1)
				{
					std::get<2>(result[i][j]) = newChunkToExtrapolatedIndex.at(std::get<2>(modifiedChunksPerRead[i][j]) & maskUint64_t) + firstBitUint64_t;
				}
				else
				{
					size_t index = nextNum;
					newChunkToExtrapolatedIndex[std::get<2>(modifiedChunksPerRead[i][j]) & maskUint64_t] = index;
					nextNum += 1;
					std::get<2>(result[i][j]) = index + firstBitUint64_t;
				}
				continue;
			}
			assert(!canonical[std::get<2>(chunksPerReadBeforeModification[otheri][otherj]) & maskUint64_t]);
			if (newChunkToExtrapolatedIndex.count(std::get<2>(modifiedChunksPerRead[i][j]) & maskUint64_t) == 1)
			{
				std::get<2>(result[i][j]) = newChunkToExtrapolatedIndex.at(std::get<2>(modifiedChunksPerRead[i][j]) & maskUint64_t) + firstBitUint64_t;
			}
			else
			{
				size_t index = nextNum;
				newChunkToExtrapolatedIndex[std::get<2>(modifiedChunksPerRead[i][j]) & maskUint64_t] = index;
				nextNum += 2;
				std::get<2>(result[i][j]) = index + firstBitUint64_t;
			}
		}
	}
	return result;
}
