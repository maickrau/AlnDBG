#include "UnionFind.h"

std::pair<size_t, size_t> find(std::vector<std::vector<std::pair<size_t, size_t>>>& parent, std::pair<size_t, size_t> index)
{
	while (true)
	{
		std::pair<size_t, size_t> oneparent = parent[index.first][index.second];
		std::pair<size_t, size_t> twoparent = parent[oneparent.first][oneparent.second];
		if (oneparent == twoparent) break;
		parent[index.first][index.second] = parent[twoparent.first][twoparent.second];
	}
	return parent[index.first][index.second];
}

void merge(std::vector<std::vector<std::pair<size_t, size_t>>>& parent, std::pair<size_t, size_t> left, std::pair<size_t, size_t> right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent[left.first][left.second] == left);
	assert(parent[right.first][right.second] == right);
	parent[right.first][right.second] = left;
}

size_t find(phmap::flat_hash_map<size_t, size_t>& parent, size_t val)
{
	while (parent.at(val) != parent.at(parent.at(val)))
	{
		parent[val] = parent.at(parent.at(val));
	}
	return parent.at(val);
}

void merge(phmap::flat_hash_map<size_t, size_t>& parent, size_t left, size_t right)
{
	left = find(parent, left);
	right = find(parent, right);
	parent[right] = left;
}

size_t find(std::vector<size_t>& parent, size_t val)
{
	while (parent[val] != parent[parent[val]])
	{
		parent[val] = parent[parent[val]];
	}
	return parent[val];
}

void merge(std::vector<size_t>& parent, size_t left, size_t right)
{
	left = find(parent, left);
	right = find(parent, right);
	parent[right] = left;
}
