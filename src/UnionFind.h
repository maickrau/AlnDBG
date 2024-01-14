#ifndef UnionFind_h
#define UnionFind_h

#include <vector>
#include <phmap.h>

void merge(phmap::flat_hash_map<size_t, size_t>& parent, size_t left, size_t right);
size_t find(phmap::flat_hash_map<size_t, size_t>& parent, size_t val);
void merge(std::vector<std::pair<size_t, bool>>& parent, size_t left, size_t right, bool fw);
std::pair<size_t, bool> find(std::vector<std::pair<size_t, bool>>& parent, size_t val);
void merge(std::vector<size_t>& parent, size_t left, size_t right);
size_t find(std::vector<size_t>& parent, size_t val);
void merge(std::vector<std::vector<std::pair<size_t, size_t>>>& parent, std::pair<size_t, size_t> left, std::pair<size_t, size_t> right);
std::pair<size_t, size_t> find(std::vector<std::vector<std::pair<size_t, size_t>>>& parent, std::pair<size_t, size_t> index);

#endif
