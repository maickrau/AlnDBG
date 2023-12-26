#include <cassert>
#include <mutex>
#include <thread>
#include <phmap.h>
#include <iostream>
#include "KmerMatcher.h"
#include "RankBitvector.h"

void heapify(std::vector<uint16_t>& vec)
{
	size_t pos = 0;
	while (pos*2+2 < vec.size())
	{
		if (vec[pos] < vec[pos*2+2] && vec[pos] < vec[pos*2+1]) return;
		if (vec[pos*2+2] < vec[pos*2+1])
		{
			std::swap(vec[pos], vec[pos*2+2]);
			pos = pos*2+2;
		}
		else
		{
			assert(vec[pos*2+1] < vec[pos*2+2]);
			std::swap(vec[pos], vec[pos*2+1]);
			pos = pos*2+1;
		}
	}
	if (pos*2+1 < vec.size())
	{
		if (vec[pos*2+1] < vec[pos])
		{
			std::swap(vec[pos], vec[pos*2+1]);
		}
	}
}

class KmerHashtableContainer
{
public:
	void addKmer(uint64_t kmer, size_t pos)
	{
		if (firstKmerPositionInLeft.count(kmer) == 0)
		{
			firstKmerPositionInLeft[kmer] = 0xFFFFFFFFFFFF0000ull + pos;
		}
		else
		{
			auto& val = firstKmerPositionInLeft[kmer];
			if ((val & 0x00000000FFFF0000ull) == 0x00000000FFFF0000ull)
			{
				val &= 0xFFFFFFFF0000FFFFull;
				val += (uint64_t)pos << 16ull;
			}
			else if ((val & 0x0000FFFF00000000ull) == 0x0000FFFF00000000ull)
			{
				val &= 0xFFFF0000FFFFFFFFull;
				val += (uint64_t)pos << 32ull;
			}
			else if ((val & 0xFFFF000000000000ull) == 0xFFFF000000000000ull)
			{
				val &= 0x0000FFFFFFFFFFFFull;
				val += (uint64_t)pos << 48ull;
			}
			else
			{
				extraKmerPositionsInLeft[kmer].push_back(pos);
			}
		}
	}
	void removeKmer(uint64_t kmer, size_t pos)
	{
		auto found = firstKmerPositionInLeft.find(kmer);
		assert(found != firstKmerPositionInLeft.end());
		bool removed = false;
		if ((found->second & 0x000000000000FFFFull) == pos)
		{
			found->second >>= 16;
			found->second |= 0xFFFF000000000000ull;
			if (found->second == 0xFFFFFFFFFFFFFFFFull)
			{
				assert(extraKmerPositionsInLeft.count(kmer) == 0);
				firstKmerPositionInLeft.erase(found);
				return;
			}
			removed = true;
		}
		auto found2 = extraKmerPositionsInLeft.find(kmer);
		if (found2 == extraKmerPositionsInLeft.end()) return;
		assert(found2->second.size() >= 1);
		if (removed) found->second = (found->second & 0x0000FFFFFFFFFFFFull) + ((uint64_t)found2->second[0] << 48ull);
		std::swap(found2->second[0], found2->second.back());
		found2->second.pop_back();
		if (found2->second.size() == 0)
		{
			extraKmerPositionsInLeft.erase(found2);
		}
		else
		{
			heapify(found2->second);
		}
	}
	template <typename F>
	void iterateKmerMatchPositions(const uint64_t kmer, F callback) const
	{
		auto found = firstKmerPositionInLeft.find(kmer);
		if (found == firstKmerPositionInLeft.end()) return;
		callback(found->second & 0x000000000000FFFFull);
		if ((found->second & 0x00000000FFFF0000ull) == 0x00000000FFFF0000ull) return;
		callback((found->second >> 16ull) & 0x000000000000FFFFull);
		if ((found->second & 0x0000FFFF00000000ull) == 0x0000FFFF00000000ull) return;
		callback((found->second >> 32ull) & 0x000000000000FFFFull);
		if ((found->second & 0xFFFF000000000000ull) == 0xFFFF000000000000ull) return;
		callback((found->second >> 48ull) & 0x000000000000FFFFull);
		auto found2 = extraKmerPositionsInLeft.find(kmer);
		if (found2 == extraKmerPositionsInLeft.end()) return;
		for (auto pos : found2->second) callback(pos);
	}
private:
	phmap::flat_hash_map<uint64_t, uint64_t> firstKmerPositionInLeft;
	phmap::flat_hash_map<uint64_t, std::vector<uint16_t>> extraKmerPositionsInLeft;
};

class HashKmerContainer
{
public:
	void initializeKeys(const phmap::flat_hash_set<uint64_t>& hashes)
	{
		for (uint64_t hash : hashes)
		{
			size_t index = hashToIndex.size();
			hashToIndex[hash] = index;
		}
		activeLocationsPerHash.resize(hashToIndex.size());
	}
	void addKmer(const uint64_t kmer, const size_t pos)
	{
		activeLocationsPerHash[hashToIndex.at(kmer)].push_back(pos);
	}
	void removeKmer(const uint64_t kmer, const size_t pos)
	{
		size_t index = hashToIndex.at(kmer);
		assert(activeLocationsPerHash[index].size() >= 1);
		assert(activeLocationsPerHash[index][0] == pos);
		std::swap(activeLocationsPerHash[index][0], activeLocationsPerHash[index].back());
		activeLocationsPerHash[index].pop_back();
		heapify(activeLocationsPerHash[index]);
	}
	void addKmerIndex(const size_t index, const size_t pos)
	{
		activeLocationsPerHash[index].push_back(pos);
	}
	void removeKmerIndex(const size_t index, const size_t pos)
	{
		assert(activeLocationsPerHash[index].size() >= 1);
		assert(activeLocationsPerHash[index][0] == pos);
		std::swap(activeLocationsPerHash[index][0], activeLocationsPerHash[index].back());
		activeLocationsPerHash[index].pop_back();
		heapify(activeLocationsPerHash[index]);
	}
	template <typename F>
	void iterateKmerMatchPositions(const uint64_t kmer, F callback) const
	{
		auto found = hashToIndex.find(kmer);
		if (found == hashToIndex.end()) return;
		for (auto pos : activeLocationsPerHash[found->second])
		{
			callback(pos);
		}
	}
	size_t getIndex(const uint64_t kmer) const
	{
		return hashToIndex.at(kmer);
	}
private:
	std::vector<std::vector<uint16_t>> activeLocationsPerHash;
	phmap::flat_hash_map<uint64_t, size_t> hashToIndex;
};

class NonminimalWeirdoKmerContainer
{
public:
	void initializeKeys(const phmap::flat_hash_set<uint64_t>& hashes)
	{
		std::vector<std::vector<uint64_t>> hashBuckets;
		hashBuckets.resize(numBuckets);
		for (auto hash : hashes)
		{
			hashBuckets[hash % numBuckets].emplace_back(hash);
		}
		size_t totalSize = 0;
		for (size_t i = 0; i < hashBuckets.size(); i++)
		{
			std::sort(hashBuckets[i].begin(), hashBuckets[i].end());
			if (hashBuckets[i].size() == 0)
			{
				countBeforeBucket.emplace_back(totalSize);
				minus.push_back(0);
				continue;
			}
			countBeforeBucket.push_back(totalSize);
			minus.push_back(hashBuckets[i][0]);
			size_t maxDivisor = std::numeric_limits<size_t>::max();
			for (size_t j = 1; j < hashBuckets[i].size(); j++)
			{
				maxDivisor = std::min(maxDivisor, hashBuckets[i][j] - hashBuckets[i][j-1]);
				while ((hashBuckets[i][j] - minus.back()) / maxDivisor == (hashBuckets[i][j-1] - minus.back()) / maxDivisor) maxDivisor -= 1;
			}
			assert(maxDivisor >= numBuckets);
			bucketDivisors.push_back(maxDivisor);
			bucketMax.push_back((hashBuckets[i].back() - minus.back()) / maxDivisor);
			totalSize += bucketMax.back();
			std::cerr << "bucket " << i << " divisor " << bucketDivisors.back() << " size " << bucketMax.back() << std::endl;
			std::cerr << "values:";
			for (auto hash : hashBuckets[i])
			{
				std::cerr << " " << hash;
			}
			std::cerr << std::endl;
		}
		std::cerr << "try build phf with size " << totalSize << std::endl;
		realHashes.resize(totalSize, std::numeric_limits<size_t>::max());
		activeLocationsPerHash.resize(totalSize);
		for (auto hash : hashes)
		{
			size_t bucket = hash % numBuckets;
			assert(hash >= minus[bucket]);
			size_t index = (hash-minus[bucket]) / bucketDivisors[bucket];
			assert(index <= bucketMax[bucket]);
			assert(realHashes[index] == std::numeric_limits<size_t>::max());
			realHashes[index] = hash;
		}
	}
	void addKmer(const uint64_t hash, const size_t pos)
	{
		activeLocationsPerHash[getIndex(hash)].emplace_back(pos);
	}
	void removeKmer(const uint64_t hash, const uint16_t pos)
	{
		size_t index = getIndex(hash);
		assert(activeLocationsPerHash[index].size() >= 1);
		assert(activeLocationsPerHash[index][0] == pos);
		std::swap(activeLocationsPerHash[index][0], activeLocationsPerHash[index].back());
		activeLocationsPerHash[index].pop_back();
		heapify(activeLocationsPerHash[index]);
	}
	template <typename F>
	void iterateKmerMatchPositions(const uint64_t kmer, F callback) const
	{
		size_t index = getIndex(kmer);
		if (index == std::numeric_limits<size_t>::max()) return;
		if (realHashes[index] != kmer) return;
		for (auto pos : activeLocationsPerHash[index])
		{
			callback(pos);
		}
	}
private:
	size_t getIndex(uint64_t hash) const
	{
		size_t bucket = hash % numBuckets;
		if (hash >= minus[bucket]) return std::numeric_limits<size_t>::max();
		size_t offset = (hash - minus[bucket]) / bucketDivisors[bucket];
		if (offset > bucketMax[bucket]) return std::numeric_limits<size_t>::max();
		return countBeforeBucket[bucket] + offset;
	}
	const size_t numBuckets = 1000;
	std::vector<uint64_t> realHashes;
	std::vector<std::vector<uint16_t>> activeLocationsPerHash;
	std::vector<size_t> bucketDivisors;
	std::vector<size_t> countBeforeBucket;
	std::vector<size_t> bucketMax;
	std::vector<size_t> minus;
};

class NonMinimalOneLevelPHFContainer
{
public:
	void initializeKeys(const phmap::flat_hash_set<uint64_t>& hashes)
	{
		valid.resize(hashes.size() * 23);
		std::vector<bool> collision;
		collision.resize(valid.size());
		for (uint64_t hash : hashes)
		{
			size_t pos = hash % valid.size();
			if (collision[pos]) continue;
			if (valid[pos])
			{
				collision[pos] = true;
				valid[pos] = false;
			}
			valid[pos] = true;
		}
		for (uint64_t hash : hashes)
		{
			size_t pos = hash % valid.size();
			if (!valid[pos])
			{
				size_t index = valid.size() + remainers.size();
				remainers[hash] = index;
			}
		}
		realHashes.resize(valid.size() + remainers.size(), std::numeric_limits<size_t>::max());
		activeLocationsPerHash.resize(realHashes.size());
		for (auto hash : hashes)
		{
			realHashes[getIndex(hash)] = hash;
		}
	}
	void addKmer(const uint64_t hash, const size_t pos)
	{
		activeLocationsPerHash[getIndex(hash)].emplace_back(pos);
	}
	void removeKmer(const uint64_t hash, const uint16_t pos)
	{
		size_t index = getIndex(hash);
		assert(activeLocationsPerHash[index].size() >= 1);
		assert(activeLocationsPerHash[index][0] == pos);
		std::swap(activeLocationsPerHash[index][0], activeLocationsPerHash[index].back());
		activeLocationsPerHash[index].pop_back();
		heapify(activeLocationsPerHash[index]);
	}
	void addKmerIndex(const size_t index, const size_t pos)
	{
		activeLocationsPerHash[index].emplace_back(pos);
	}
	void removeKmerIndex(const size_t index, const uint16_t pos)
	{
		assert(activeLocationsPerHash[index].size() >= 1);
		assert(activeLocationsPerHash[index][0] == pos);
		std::swap(activeLocationsPerHash[index][0], activeLocationsPerHash[index].back());
		activeLocationsPerHash[index].pop_back();
		heapify(activeLocationsPerHash[index]);
	}
	template <typename F>
	void iterateKmerMatchPositions(const uint64_t kmer, F callback) const
	{
		size_t index = getIndex(kmer);
		if (index == std::numeric_limits<size_t>::max()) return;
		if (realHashes[index] != kmer) return;
		for (auto pos : activeLocationsPerHash[index])
		{
			callback(pos);
		}
	}
	size_t getIndex(uint64_t kmer) const
	{
		size_t pos = kmer % valid.size();
		if (valid[pos]) return pos;
		auto found = remainers.find(kmer);
		if (found == remainers.end()) return std::numeric_limits<size_t>::max();
		return found->second;
	}
private:
	std::vector<bool> valid;
	std::vector<uint64_t> realHashes;
	std::vector<std::vector<uint16_t>> activeLocationsPerHash;
	phmap::flat_hash_map<uint64_t, size_t> remainers;
};

class MPHFKmerContainer
{
public:
	void initializeKeys(const phmap::flat_hash_set<uint64_t>& hashes)
	{
		activeLocationsPerHash.resize(0);
		realHashes.resize(0);
		buckets.resize(0);
		countBeforeBucket.resize(0);
		buckets.resize(numBuckets);
		countBeforeBucket.resize(numBuckets);
		std::vector<phmap::flat_hash_set<uint64_t>> hashesPerBucket;
		hashesPerBucket.resize(numBuckets);
		for (uint64_t hash : hashes)
		{
			hashesPerBucket[hash % numBuckets].emplace(hash);
		}
		for (size_t i = 0; i < numBuckets; i++)
		{
			if (i > 0)
			{
				countBeforeBucket[i] = countBeforeBucket[i-1] + hashesPerBucket[i-1].size();
			}
			buildBucket(i, hashesPerBucket[i]);
		}
		activeLocationsPerHash.resize(countBeforeBucket.back() + hashesPerBucket.back().size());
		realHashes.resize(activeLocationsPerHash.size());
		for (uint64_t hash : hashes)
		{
			realHashes[getIndex(hash)] = hash;
		}
	}
	void addKmer(const uint64_t hash, const size_t pos)
	{
		activeLocationsPerHash[getIndex(hash)].emplace_back(pos);
	}
	void removeKmer(const uint64_t hash, const uint16_t pos)
	{
		size_t index = getIndex(hash);
		assert(activeLocationsPerHash[index].size() >= 1);
		assert(activeLocationsPerHash[index][0] == pos);
		std::swap(activeLocationsPerHash[index][0], activeLocationsPerHash[index].back());
		activeLocationsPerHash[index].pop_back();
		heapify(activeLocationsPerHash[index]);
	}
	void addKmerIndex(const size_t index, const size_t pos)
	{
		activeLocationsPerHash[index].emplace_back(pos);
	}
	void removeKmerIndex(const size_t index, const uint16_t pos)
	{
		assert(activeLocationsPerHash[index].size() >= 1);
		assert(activeLocationsPerHash[index][0] == pos);
		std::swap(activeLocationsPerHash[index][0], activeLocationsPerHash[index].back());
		activeLocationsPerHash[index].pop_back();
		heapify(activeLocationsPerHash[index]);
	}
	template <typename F>
	void iterateKmerMatchPositions(const uint64_t kmer, F callback) const
	{
		size_t index = getIndex(kmer);
		if (index == std::numeric_limits<size_t>::max()) return;
		if (realHashes[index] != kmer) return;
		for (auto pos : activeLocationsPerHash[index])
		{
			callback(pos);
		}
	}
	size_t getIndex(uint64_t hash) const
	{
		size_t bucket = hash % numBuckets;
		size_t indexInBucket = buckets[bucket].getIndex(hash);
		if (indexInBucket == std::numeric_limits<size_t>::max()) return indexInBucket;
		return countBeforeBucket[bucket] + indexInBucket;
	}
private:
	void buildBucket(const size_t bucket, const phmap::flat_hash_set<uint64_t>& bucketHashes)
	{
		buckets[bucket].build(bucketHashes);
	}
	class Bucket
	{
	public:
		void build(const phmap::flat_hash_set<uint64_t>& keys)
		{
			if (keys.size() == 0)
			{
				return;
			}
			phmap::flat_hash_set<uint64_t> stillNeedsHashing = keys;
			size_t hashedSoFar = 0;
			for (size_t i = 0; i < 1; i++)
			{
//				std::cerr << "trying to build a bucket with " << keys.size() << " hashes, level " << i << ", still needs hashing " << stillNeedsHashing.size() << std::endl;
				valids.emplace_back();
				collisions.emplace_back();
				valids.back().resize((((stillNeedsHashing.size() * 11)+63)/64)*64-1);
				collisions.back().resize(valids.back().size());
				countsBeforeLevel.emplace_back(hashedSoFar);
				phmap::flat_hash_set<uint64_t> collidingHashes;
				for (uint64_t hash : stillNeedsHashing)
				{
					size_t pos = (hash * 17) % valids.back().size();
					if (valids.back().get(pos))
					{
						collisions.back().set(pos, true);
					}
					valids.back().set(pos, true);
				}
				for (uint64_t hash : stillNeedsHashing)
				{
					size_t pos = (hash * 17) % valids.back().size();
					if (collisions.back().get(pos))
					{
						collidingHashes.emplace(hash);
					}
				}
				valids.back().buildRanks();
				collisions.back().buildRanks();
				hashedSoFar += stillNeedsHashing.size() - collidingHashes.size();
				if (collidingHashes.size() == stillNeedsHashing.size())
				{
					valids.pop_back();
					collisions.pop_back();
					countsBeforeLevel.pop_back();
					break;
				}
				stillNeedsHashing = collidingHashes;
				if (stillNeedsHashing.size() == 0) return;
			}
			for (uint64_t hash : stillNeedsHashing)
			{
				lastColliderIndices[hash] = hashedSoFar;
				hashedSoFar += 1;
			}
		}
		size_t getIndex(uint64_t hash) const
		{
			for (size_t i = 0; i < valids.size(); i++)
			{
				size_t pos = (hash * 17) % valids[i].size();
				if (!valids[i].get(pos)) return std::numeric_limits<size_t>::max();
				if (collisions[i].get(pos)) continue;
				return countsBeforeLevel[i] + valids[i].getRank(pos);
			}
			auto found = lastColliderIndices.find(hash);
			if (found == lastColliderIndices.end()) return std::numeric_limits<size_t>::max();
			return found->second;
		}
		std::vector<RankBitvector> valids;
		std::vector<RankBitvector> collisions;
		std::vector<size_t> countsBeforeLevel;
		phmap::flat_hash_map<uint64_t, size_t> lastColliderIndices;
	};
	const size_t numBuckets = 1;
	std::vector<uint64_t> realHashes;
	std::vector<std::vector<uint16_t>> activeLocationsPerHash;
	std::vector<Bucket> buckets;
	std::vector<size_t> countBeforeBucket;
};

template <typename F1, typename F2>
void iterateSyncmers(const std::vector<TwobitString>& readSequences, const size_t k, const size_t w, const size_t maxLen, F1 callback, F2 getchar)
{
	thread_local std::vector<std::tuple<size_t, uint64_t>> smerOrder;
	assert(w < k);
	assert(w >= 3);
	const uint64_t mask = (1ull << (2ull*k)) - 1;
	const size_t s = k-w+1;
	const uint64_t smask = (1ull << (2ull*s)) - 1;
	uint64_t kmer = 0;
	uint64_t smer = 0;
	for (size_t i = 0; i < s; i++)
	{
		uint64_t c = getchar(i);
		smer <<= 2;
		smer += c;
	}
	kmer = smer;
	smerOrder.emplace_back(0, smer);
	for (size_t i = s; i < k; i++)
	{
		uint64_t c = getchar(i);
		kmer <<= 2;
		kmer += c;
		smer <<= 2;
		smer += c;
		smer &= smask;
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > smer) smerOrder.pop_back();
		smerOrder.emplace_back(i-s+1, smer);
	}
	if ((std::get<0>(smerOrder.front()) == 0) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == w-1))
	{
		callback(kmer, 0);
	}
	for (size_t i = k; i < maxLen; i++)
	{
		uint64_t c = getchar(i);
		kmer <<= 2;
		kmer += c;
		kmer &= mask;
		smer <<= 2;
		smer += c;
		smer &= smask;
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > smer) smerOrder.pop_back();
		while (smerOrder.size() > 0 && std::get<0>(smerOrder.front()) <= i-s+1-w) smerOrder.erase(smerOrder.begin());
		smerOrder.emplace_back(i-s+1, smer);
		if ((std::get<0>(smerOrder.front()) == i-s+2-w) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == i-s+1))
		{
			callback(kmer, i-k+1);
		}
	}
	smerOrder.clear();
}

template <typename F>
void iterateSyncmers(const std::vector<TwobitString>& readSequences, const size_t k, const size_t w, const size_t read, const size_t readStart, const size_t readEnd, const bool fw, F callback)
{
	if (fw)
	{
		iterateSyncmers(readSequences, k, w, readEnd-readStart, callback, [&readSequences, read, readStart, readEnd](size_t index){ return readSequences[read].get(readStart+index); });
	}
	else
	{
		iterateSyncmers(readSequences, k, w, readEnd-readStart, callback, [&readSequences, read, readStart, readEnd](size_t index){ return 3-readSequences[read].get(readSequences[read].size() - 1 - (readStart+index)); });
	}
}

void getKmerMatches(const std::vector<TwobitString>& readSequences, MatchGroup& mappingMatch, const size_t k, const size_t d)
{
	thread_local size_t lastLeftRead = std::numeric_limits<size_t>::max();
	thread_local size_t lastLeftStart = std::numeric_limits<size_t>::max();
	thread_local size_t lastLeftEnd = std::numeric_limits<size_t>::max();
	thread_local std::vector<std::pair<size_t, size_t>> leftSyncmers;
	thread_local HashKmerContainer kmerContainer;
	assert(mappingMatch.leftEnd - mappingMatch.leftStart < std::numeric_limits<uint16_t>::max());
	assert(mappingMatch.rightEnd - mappingMatch.rightStart < std::numeric_limits<uint16_t>::max());
	const size_t syncmerw = 13;
	assert(k <= 31);
	assert(k % 2 == 1);
	size_t leftStart = mappingMatch.leftStart;
	size_t leftEnd = mappingMatch.leftEnd;
	size_t rightStart = mappingMatch.rightStart;
	size_t rightEnd = mappingMatch.rightEnd;
	bool rightFw = mappingMatch.rightFw;
	size_t left = mappingMatch.leftRead;
	size_t right = mappingMatch.rightRead;
	if (left != lastLeftRead)
	{
		phmap::flat_hash_set<uint64_t> hashes;
		iterateSyncmers(readSequences, k, syncmerw, left, 0, readSequences[left].size(), true, [&hashes](const size_t kmer, const size_t pos)
		{
			hashes.emplace(kmer);
		});
		kmerContainer.initializeKeys(hashes);
	}
	if (left != lastLeftRead || leftStart != lastLeftStart || leftEnd > lastLeftEnd)
	{
		leftSyncmers.clear();
		iterateSyncmers(readSequences, k, syncmerw, left, leftStart, leftEnd, true, [&leftSyncmers, &kmerContainer](const size_t kmer, const size_t pos)
		{
			leftSyncmers.emplace_back(kmerContainer.getIndex(kmer), pos);
		});
		lastLeftRead = left;
		lastLeftStart = leftStart;
		lastLeftEnd = leftEnd;
	}
	std::vector<std::pair<size_t, size_t>> currentMatchesPerDiagonal;
	size_t diagonalCount = 2*d+1;
	size_t zeroDiagonal = d;
	if (leftEnd-leftStart > rightEnd-rightStart)
	{
		diagonalCount += (leftEnd-leftStart)-(rightEnd-rightStart);
		zeroDiagonal = d + (leftEnd-leftStart)-(rightEnd-rightStart);
	}
	if (rightEnd-rightStart > leftEnd-leftStart) diagonalCount += (rightEnd-rightStart)-(leftEnd-leftStart);
	currentMatchesPerDiagonal.resize(diagonalCount, std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	size_t leftLeadingIndex = 0;
	size_t leftTrailingIndex = 0;
	iterateSyncmers(readSequences, k, syncmerw, right, rightStart, rightEnd, rightFw, [&kmerContainer, &currentMatchesPerDiagonal, diagonalCount, zeroDiagonal, rightFw, left, right, leftStart, rightStart, leftEnd, rightEnd, d, k, &mappingMatch, &leftTrailingIndex, &leftLeadingIndex, &leftSyncmers](const size_t kmer, const size_t rightPos)
	{
		size_t interpolatedLeftPos = (double)(rightPos) / (double)(rightEnd-rightStart) * (double)(leftEnd-leftStart);
		while (leftLeadingIndex < leftSyncmers.size() && leftSyncmers[leftLeadingIndex].second <= interpolatedLeftPos + d && leftSyncmers[leftLeadingIndex].second < leftEnd)
		{
			kmerContainer.addKmerIndex(leftSyncmers[leftLeadingIndex].first, leftSyncmers[leftLeadingIndex].second);
			leftLeadingIndex += 1;
		}
		while (leftTrailingIndex < leftSyncmers.size() && leftSyncmers[leftTrailingIndex].second + d < interpolatedLeftPos && leftSyncmers[leftTrailingIndex].second < leftEnd)
		{
			kmerContainer.removeKmerIndex(leftSyncmers[leftTrailingIndex].first, leftSyncmers[leftTrailingIndex].second);
			leftTrailingIndex += 1;
		}
		assert(rightPos + zeroDiagonal >= interpolatedLeftPos + d);
		assert(rightPos + zeroDiagonal + d >= interpolatedLeftPos);
		size_t minDiagonal = rightPos + zeroDiagonal - interpolatedLeftPos - d;
		size_t maxDiagonal = rightPos + zeroDiagonal - interpolatedLeftPos + d;
		assert(maxDiagonal <= diagonalCount);
		kmerContainer.iterateKmerMatchPositions(kmer, [zeroDiagonal, rightPos, minDiagonal, maxDiagonal, &currentMatchesPerDiagonal, &mappingMatch, rightFw, left, right, leftStart, rightStart, diagonalCount, k](const size_t leftPos)
		{
			if (leftPos > zeroDiagonal + rightPos) return;
			if (zeroDiagonal + rightPos - leftPos >= diagonalCount) return;
			size_t diagonal = zeroDiagonal + rightPos - leftPos;
			if (diagonal < minDiagonal || diagonal > maxDiagonal) return;
			if (currentMatchesPerDiagonal[diagonal].first != std::numeric_limits<size_t>::max() && currentMatchesPerDiagonal[diagonal].second + k > rightPos)
			{
				currentMatchesPerDiagonal[diagonal].second = rightPos+1;
				return;
			}
			if (currentMatchesPerDiagonal[diagonal].first != std::numeric_limits<size_t>::max())
			{
				assert(currentMatchesPerDiagonal[diagonal].second != std::numeric_limits<size_t>::max());
				assert(currentMatchesPerDiagonal[diagonal].second > currentMatchesPerDiagonal[diagonal].first);
				assert(currentMatchesPerDiagonal[diagonal].first + zeroDiagonal >= diagonal);
				size_t leftMatchStart = currentMatchesPerDiagonal[diagonal].first + zeroDiagonal - diagonal;
				size_t rightMatchStart = currentMatchesPerDiagonal[diagonal].first;
				size_t length = currentMatchesPerDiagonal[diagonal].second - currentMatchesPerDiagonal[diagonal].first;
				assert(leftMatchStart < std::numeric_limits<uint16_t>::max());
				assert(rightMatchStart < std::numeric_limits<uint16_t>::max());
				assert(length < std::numeric_limits<uint16_t>::max());
				mappingMatch.matches.emplace_back();
				mappingMatch.matches.back().leftStart = leftMatchStart;
				mappingMatch.matches.back().rightStart = rightMatchStart;
				mappingMatch.matches.back().length = length;
			}
			currentMatchesPerDiagonal[diagonal].first = rightPos;
			currentMatchesPerDiagonal[diagonal].second = rightPos+1;
		});
	});
	while (leftTrailingIndex < leftLeadingIndex)
	{
		kmerContainer.removeKmerIndex(leftSyncmers[leftTrailingIndex].first, leftSyncmers[leftTrailingIndex].second);
		leftTrailingIndex += 1;
	}
	for (size_t diagonal = 0; diagonal < diagonalCount; diagonal++)
	{
		if (currentMatchesPerDiagonal[diagonal].first == std::numeric_limits<size_t>::max()) continue;
		assert(currentMatchesPerDiagonal[diagonal].second != std::numeric_limits<size_t>::max());
		assert(currentMatchesPerDiagonal[diagonal].second > currentMatchesPerDiagonal[diagonal].first);
		assert(currentMatchesPerDiagonal[diagonal].first + zeroDiagonal >= diagonal);
		size_t leftMatchStart = currentMatchesPerDiagonal[diagonal].first + zeroDiagonal - diagonal;
		size_t rightMatchStart = currentMatchesPerDiagonal[diagonal].first;
		size_t length = currentMatchesPerDiagonal[diagonal].second - currentMatchesPerDiagonal[diagonal].first;
		assert(leftMatchStart < std::numeric_limits<uint16_t>::max());
		assert(rightMatchStart < std::numeric_limits<uint16_t>::max());
		assert(length < std::numeric_limits<uint16_t>::max());
		mappingMatch.matches.emplace_back();
		mappingMatch.matches.back().leftStart = leftMatchStart;
		mappingMatch.matches.back().rightStart = rightMatchStart;
		mappingMatch.matches.back().length = length;
	}
}

bool matchContained(MatchGroup::Match smaller, MatchGroup::Match bigger)
{
	if (smaller.leftStart <= bigger.leftStart) return false;
	if (smaller.rightStart <= bigger.rightStart) return false;
	if (smaller.leftStart + smaller.length >= bigger.leftStart + bigger.length) return false;
	if (smaller.rightStart + smaller.length >= bigger.rightStart + bigger.length) return false;
	return true;
}

void removeContainedKmerMatches(MatchGroup& matches)
{
	std::sort(matches.matches.begin(), matches.matches.end(), [](auto left, auto right) { return left.length < right.length; });
	std::vector<bool> remove;
	remove.resize(matches.matches.size(), false);
	for (size_t i = 0; i < matches.matches.size(); i++)
	{
		for (size_t j = matches.matches.size()-1; j > i; j--)
		{
			if (matches.matches[j].length < matches.matches[i].length+2) break;
			if (!matchContained(matches.matches[i], matches.matches[j])) continue;
			remove[i] = true;
			break;
		}
	}
	for (size_t i = matches.matches.size()-1; i < matches.matches.size(); i--)
	{
		if (!remove[i]) continue;
		std::swap(matches.matches[i], matches.matches.back());
		matches.matches.pop_back();
	}
}

void addKmerMatches(const size_t numThreads, const std::vector<TwobitString>& readSequences, std::vector<MatchGroup>& matches, const size_t graphk, const size_t graphd)
{
	std::atomic<size_t> kmerMatchCount;
	kmerMatchCount = 0;
	size_t nextIndex = 0;
	std::mutex indexMutex;
	std::vector<std::thread> threads;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&matches, &readSequences, graphk, graphd, &nextIndex, &indexMutex, &kmerMatchCount](){
			while (true)
			{
				size_t startIndex = 0;
				size_t endIndex = 0;
				{
					std::lock_guard<std::mutex> lock { indexMutex };
					startIndex = nextIndex;
					nextIndex += 1;
					while (nextIndex < matches.size() && matches[nextIndex].leftRead == matches[nextIndex-1].leftRead && matches[nextIndex].leftStart == matches[nextIndex-1].leftStart) nextIndex += 1;
					endIndex = nextIndex;
				}
				if (startIndex >= matches.size()) break;
				for (size_t i = startIndex; i < endIndex; i++)
				{
					getKmerMatches(readSequences, matches[i], graphk, graphd);
					removeContainedKmerMatches(matches[i]);
					std::sort(matches[i].matches.begin(), matches[i].matches.end(), [](auto left, auto right) { return left.leftStart < right.leftStart; });
					kmerMatchCount += matches[i].matches.size();
				}
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	std::cerr << kmerMatchCount << " kmer matches" << std::endl;
}
