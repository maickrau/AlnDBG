#include <limits>
#include "OverlapMatcher.h"

std::vector<std::pair<size_t, size_t>> getLocallyUniqueKmers(const std::string& readSeq, const size_t kmerSize)
{
	std::vector<std::pair<size_t, size_t>> result;
	iterateLocallyUniqueKmers(readSeq, kmerSize, 100, [&result](const size_t kmer, const size_t pos)
	{
		result.emplace_back(kmer, pos);
	});
	std::sort(result.begin(), result.end());
	return result;
}

ReadOverlapInformation getReadOverlapInformation(const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t readIndex, const size_t kmerSize)
{
	ReadOverlapInformation result;
	std::string readSeq;
	if (readIndex < sequenceIndex.size())
	{
		readSeq = sequenceIndex.getSequence(readIndex);
	}
	else
	{
		readSeq = revCompRaw(sequenceIndex.getSequence(readIndex - sequenceIndex.size()));
	}
	result.readKmers = getLocallyUniqueKmers(readSeq, kmerSize);
	result.readIndex = readIndex;
	result.readLength = readSeq.size();
	return result;
}
