#include <thread>
#include <chrono>
#include "ChunkExtractor.h"
#include "Common.h"
#include "fastqloader.h"
#include "KmerIterator.h"

std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> getMinimizerBoundedChunksPerRead(const FastaCompressor::CompressedStringIndex& sequenceIndex, const std::vector<size_t>& rawReadLengths, const size_t numThreads, const size_t k, const size_t windowSize)
{
	std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>> result;
	result.resize(sequenceIndex.size()*2);
	const size_t backwardOffset = sequenceIndex.size();
	iterateMultithreaded(0, sequenceIndex.size(), numThreads, [&result, &sequenceIndex, backwardOffset, k, windowSize](const size_t readIndex)
	{
		std::string readSequence = sequenceIndex.getSequence(readIndex);
		std::vector<size_t> minimizerPositions;
		iterateMinimizers(readSequence, k, windowSize, [&minimizerPositions](const size_t pos)
		{
			assert(minimizerPositions.size() == 0 || pos > minimizerPositions.back());
			minimizerPositions.emplace_back(pos);
		});
		for (size_t i = 1; i < minimizerPositions.size(); i++)
		{
			assert(minimizerPositions[i] > minimizerPositions[i-1]);
			result[readIndex].emplace_back(minimizerPositions[i-1], minimizerPositions[i]+k-1, 0);
		};
		for (size_t i = 0; i < result[readIndex].size(); i++)
		{
			result[readIndex+backwardOffset].emplace_back(result[readIndex][result[readIndex].size()-1-i]);
			std::swap(std::get<0>(result[readIndex+backwardOffset].back()), std::get<1>(result[readIndex+backwardOffset].back()));
			std::get<0>(result[readIndex+backwardOffset].back()) = readSequence.size()-1-std::get<0>(result[readIndex+backwardOffset].back());
			std::get<1>(result[readIndex+backwardOffset].back()) = readSequence.size()-1-std::get<1>(result[readIndex+backwardOffset].back());
		}
	});
	return result;
}

void addMissingPiecesBetweenChunks(std::vector<std::vector<std::tuple<size_t, size_t, uint64_t>>>& chunksPerRead, const size_t kmerSize)
{
	std::vector<std::tuple<size_t, size_t, size_t>> additionalPieces;
	size_t maxChunk = 0;
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		if (chunksPerRead[i].size() == 0) continue;
		maxChunk = std::max(std::get<2>(chunksPerRead[i][0]) & maskUint64_t, maxChunk);
		for (size_t j = 1; j < chunksPerRead[i].size(); j++)
		{
			maxChunk = std::max(std::get<2>(chunksPerRead[i][j]) & maskUint64_t, maxChunk);
			if (std::get<1>(chunksPerRead[i][j-1]) < std::get<0>(chunksPerRead[i][j]))
			{
				size_t start = std::get<1>(chunksPerRead[i][j-1])-kmerSize+1;
				size_t end = std::get<0>(chunksPerRead[i][j])+kmerSize-1;
				assert(start > std::get<0>(chunksPerRead[i][j-1]));
				assert(end > std::get<1>(chunksPerRead[i][j-1]));
				assert(start < std::get<0>(chunksPerRead[i][j]));
				assert(end < std::get<1>(chunksPerRead[i][j]));
				additionalPieces.emplace_back(start, end, i);
			}
		}
	}
	for (auto t : additionalPieces)
	{
		chunksPerRead[std::get<2>(t)].emplace_back(std::get<0>(t), std::get<1>(t), maxChunk+1);
	}
	for (size_t i = 0; i < chunksPerRead.size(); i++)
	{
		std::sort(chunksPerRead[i].begin(), chunksPerRead[i].end());
	}
}

void readFilesAndAddToSequenceIndex(const std::string& filename, FastaCompressor::CompressedStringIndex& sequenceIndex, std::vector<size_t>& readBasepairLengths, const size_t numThreads)
{
	std::vector<std::tuple<size_t, std::string*, std::string*>> sequenceStack;
	std::vector<std::thread> threads;
	std::atomic<bool> allDone;
	allDone = false;
	std::mutex stackMutex;
	size_t nextNum = sequenceIndex.size();
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&allDone, &stackMutex, &sequenceStack, &sequenceIndex]()
		{
			while (true)
			{
				std::tuple<size_t, std::string*, std::string*> item { std::numeric_limits<size_t>::max(), nullptr, nullptr };
				{
					std::lock_guard<std::mutex> lock { stackMutex };
					if (sequenceStack.size() == 0)
					{
						if (allDone) break;
					}
					else
					{
						item = sequenceStack.back();
						sequenceStack.pop_back();
					}
				}
				if (std::get<0>(item) == std::numeric_limits<size_t>::max())
				{
					assert(std::get<1>(item) == nullptr);
					assert(std::get<2>(item) == nullptr);
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
					continue;
				}
				sequenceIndex.addString(std::get<0>(item), *std::get<1>(item), *std::get<2>(item));
				delete std::get<1>(item);
				delete std::get<2>(item);
			}
		});
	}
	FastQ::streamFastqFromFile(filename, false, [&sequenceStack, &stackMutex, &nextNum, &readBasepairLengths](FastQ& read)
	{
		std::tuple<size_t, std::string*, std::string*> item;
		readBasepairLengths.emplace_back(read.sequence.size());
		std::get<0>(item) = nextNum;
		std::get<1>(item) = new std::string;
		std::get<2>(item) = new std::string;
		std::swap(*std::get<1>(item), read.seq_id);
		std::swap(*std::get<2>(item), read.sequence);
		nextNum += 1;
		while (true)
		{
			{
				std::lock_guard<std::mutex> lock { stackMutex };
				if (sequenceStack.size() < 100)
				{
					sequenceStack.emplace_back(item);
					break;
				}
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			continue;
		}
	});
	allDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}
