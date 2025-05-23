#include "KmerIterator.h"

// https://naml.us/post/inverse-of-a-hash-function/
uint64_t hash(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

bool hasHomopolymerWithLengthTwoOrMore(uint64_t kmer, const size_t kmerSize)
{
	kmer ^= kmer >> 2;
	kmer = (~(kmer | (kmer >> 1))) & 0x5555555555555555ull;
	kmer &= (1ull << (2ull * (kmerSize-1)))-1;
	return kmer != 0;
}

bool hasHomopolymerWithLengthThreeOrMore(uint64_t kmer, const size_t kmerSize)
{
	kmer ^= kmer >> 2;
	kmer = (~(kmer | (kmer >> 1))) & 0x5555555555555555ull;
	kmer &= kmer >> 2;
	kmer &= (1ull << (2ull * (kmerSize-2)))-1;
	return kmer != 0;
}

bool hasDinucleotideMotifLengthTwo(uint64_t kmer, const size_t kmerSize)
{
	kmer ^= kmer >> 4;
	kmer = (~(kmer | (kmer >> 1))) & 0x5555555555555555ull;
	kmer &= kmer >> 2;
	kmer &= (1ull << (2ull * (kmerSize-3)))-1;
	return kmer != 0;
}

bool hasTrinucleotideMotifLengthTwo(uint64_t kmer, const size_t kmerSize)
{
	kmer ^= kmer >> 6;
	kmer = (~(kmer | (kmer >> 1))) & 0x5555555555555555ull;
	kmer = kmer & (kmer >> 2) & (kmer >> 4);
	kmer &= (1ull << (2ull * (kmerSize-5)))-1;
	return kmer != 0;
}

bool hasQuadnucleotideMotifLengthTwo(uint64_t kmer, const size_t kmerSize)
{
	kmer ^= kmer >> 8;
	kmer = (~(kmer | (kmer >> 1))) & 0x5555555555555555ull;
	kmer = kmer & (kmer >> 2) & (kmer >> 4) & (kmer >> 6);
	kmer &= (1ull << (2ull * (kmerSize-7)))-1;
	return kmer != 0;
}

bool hasPentanucleotideMotifLengthTwo(uint64_t kmer, const size_t kmerSize)
{
	kmer ^= kmer >> 10;
	kmer = (~(kmer | (kmer >> 1))) & 0x5555555555555555ull;
	kmer = kmer & (kmer >> 2) & (kmer >> 4) & (kmer >> 6) & (kmer >> 8);
	kmer &= (1ull << (2ull * (kmerSize-9)))-1;
	return kmer != 0;
}

uint64_t weightedHash(uint64_t kmer, const size_t kmerSize, const bool goodKmer)
{
	const size_t hashmask = 0x00FFFFFFFFFFFFFFull;
	uint64_t result = hash(kmer) & hashmask;
	if (!goodKmer) result |= 0x8000000000000000ull;
	if (hasHomopolymerWithLengthThreeOrMore(kmer, kmerSize)) result |= 0x4000000000000000ull;
	if (hasDinucleotideMotifLengthTwo(kmer, kmerSize)) result |= 0x2000000000000000ull;
	if (hasTrinucleotideMotifLengthTwo(kmer, kmerSize)) result |= 0x1000000000000000ull;
	if (hasQuadnucleotideMotifLengthTwo(kmer, kmerSize)) result |= 0x0800000000000000ull;
	if (hasPentanucleotideMotifLengthTwo(kmer, kmerSize)) result |= 0x0400000000000000ull;
	if (hasHomopolymerWithLengthTwoOrMore(kmer, kmerSize)) result |= 0x0200000000000000ull;
	return result;
}

uint64_t reverseKmer(uint64_t kmer, const size_t kmerSize)
{
	kmer = ~kmer;
	kmer = ((kmer >> 2) & 0x3333333333333333ull) + ((kmer << 2) & 0xCCCCCCCCCCCCCCCCull);
	kmer = ((kmer >> 4) & 0x0F0F0F0F0F0F0F0Full) + ((kmer << 4) & 0xF0F0F0F0F0F0F0F0ull);
	kmer = ((kmer >> 8) & 0x00FF00FF00FF00FFull) + ((kmer << 8) & 0xFF00FF00FF00FF00ull);
	kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFull) + ((kmer << 16) & 0xFFFF0000FFFF0000ull);
	kmer = ((kmer >> 32) & 0x00000000FFFFFFFFull) + ((kmer << 32) & 0xFFFFFFFF00000000ull);
	kmer >>= (32ull-kmerSize)*2ull;
	return kmer;
}

uint64_t hashFwAndBw(const uint64_t kmer, const size_t kmerSize)
{
	uint64_t fwHash = weightedHash(kmer, kmerSize, true);
	uint64_t bwHash = weightedHash(reverseKmer(kmer, kmerSize), kmerSize, true);
	return std::min(fwHash, bwHash);
}

uint64_t hashFwAndBw(const uint64_t kmer, const size_t kmerSize, const bool goodKmer)
{
	uint64_t fwHash = weightedHash(kmer, kmerSize, goodKmer);
	uint64_t bwHash = weightedHash(reverseKmer(kmer, kmerSize), kmerSize, goodKmer);
	return std::min(fwHash, bwHash);
}

size_t charToInt(const unsigned char c)
{
	switch(c)
	{
	case 'a':
	case 'A':
		return 0;
	case 'C':
	case 'c':
		return 1;
	case 'G':
	case 'g':
		return 2;
	case 'T':
	case 't':
		return 3;
	}
	assert(false);
	return 0;
}
