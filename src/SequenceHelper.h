#ifndef SequenceHelper_h
#define SequenceHelper_h

#include <string>
#include <vector>
#include <tuple>
#include "CompressedStringIndex.h"

std::string getReadSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const size_t index);
std::string getChunkSequence(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t read, const size_t index);
std::string getChunkSequenceMaybeMemoized(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t read, const size_t index);

#endif
