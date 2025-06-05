#include "SequenceHelper.h"
#include "Common.h"

std::string getReadSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t index)
{
	if (index < sequenceIndex.size())
	{
		return sequenceIndex.getSequence(index);
	}
	assert(index < sequenceIndex.size()*2);
	return revCompRaw(sequenceIndex.getSequence(index - sequenceIndex.size()));
}

std::string getChunkSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t read, const size_t index)
{
	if (read < sequenceIndex.size())
	{
		auto t = chunksPerRead[read][index];
		std::string result = sequenceIndex.getSubstring(read, std::get<0>(t), std::get<1>(t)-std::get<0>(t)+1);
		if (std::get<2>(t) & firstBitUint64_t)
		{
		}
		else
		{
			result = revCompRaw(result);
		}
		return result;
	}
	else
	{
		assert(read < sequenceIndex.size()*2);
		auto t = chunksPerRead[read][index];
		std::string result = sequenceIndex.getSubstring(read - sequenceIndex.size(), rawReadLengths[read-sequenceIndex.size()] - 1 - std::get<1>(t), std::get<1>(t) - std::get<0>(t) + 1);
		if (std::get<2>(t) & firstBitUint64_t)
		{
			result = revCompRaw(result);
		}
		return result;
	}
}

std::string getChunkSequenceMaybeMemoized(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t read, const size_t index)
{
	thread_local std::string lastSequence;
	thread_local size_t lastRead = std::numeric_limits<size_t>::max();
	if (read < sequenceIndex.size())
	{
		auto t = chunksPerRead[read][index];
		if (read != lastRead)
		{
			lastSequence = sequenceIndex.getSequence(read);
			lastRead = read;
		}
		std::string result = lastSequence.substr(std::get<0>(t), std::get<1>(t)-std::get<0>(t)+1);
		if (std::get<2>(t) & firstBitUint64_t)
		{
		}
		else
		{
			result = revCompRaw(result);
		}
		return result;
	}
	else
	{
		assert(read < sequenceIndex.size()*2);
		auto t = chunksPerRead[read][index];
		if (read != lastRead)
		{
			lastSequence = sequenceIndex.getSequence(read - sequenceIndex.size());
			lastRead = read;
		}
		std::string result = lastSequence.substr(rawReadLengths[read-sequenceIndex.size()] - 1 - std::get<1>(t), std::get<1>(t) - std::get<0>(t) + 1);
		if (std::get<2>(t) & firstBitUint64_t)
		{
			result = revCompRaw(result);
		}
		return result;
	}
}
