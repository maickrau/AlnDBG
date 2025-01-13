#ifndef UnionFind_h
#define UnionFind_h

#include <vector>
#include <phmap.h>

void merge(phmap::flat_hash_map<size_t, size_t>& parent, size_t left, size_t right);
size_t find(phmap::flat_hash_map<size_t, size_t>& parent, size_t val);
void merge(std::vector<std::pair<size_t, bool>>& parent, size_t left, size_t right, bool fw);
void merge(std::vector<std::pair<size_t, bool>>& parent, std::pair<size_t, bool> left, std::pair<size_t, bool> right, bool fw);
void mergeAllowPalindrome(std::vector<std::pair<size_t, bool>>& parent, std::pair<size_t, bool> left, std::pair<size_t, bool> right, bool fw);
std::pair<size_t, bool> find(std::vector<std::pair<size_t, bool>>& parent, size_t val); // bool true is fw, false is reverse orientation
std::pair<size_t, bool> find(std::vector<std::pair<size_t, bool>>& parent, std::pair<size_t, bool> val); // bool true is fw, false is reverse orientation
void merge(std::vector<size_t>& parent, size_t left, size_t right);
size_t find(std::vector<size_t>& parent, size_t val);
void merge(std::vector<std::vector<std::pair<size_t, size_t>>>& parent, std::pair<size_t, size_t> left, std::pair<size_t, size_t> right);
std::pair<size_t, size_t> find(std::vector<std::vector<std::pair<size_t, size_t>>>& parent, std::pair<size_t, size_t> index);

template <typename T>
T find(phmap::flat_hash_map<T, T>& parent, T key)
{
	while (parent.at(key) != parent.at(parent.at(key)))
	{
		parent[key] = parent.at(parent.at(key));
	}
	return parent.at(key);
}

template <typename T>
void merge(phmap::flat_hash_map<T, T>& parent, T left, T right)
{
	auto leftp = find(parent, left);
	auto rightp = find(parent, right);
	assert(parent.at(leftp) == leftp);
	assert(parent.at(rightp) == rightp);
	parent[rightp] = leftp;
}

#endif
