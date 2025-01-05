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

std::pair<size_t, bool> find(std::vector<std::pair<size_t, bool>>& parent, size_t val)
{
	assert(val < parent.size());
	std::pair<size_t, bool> result;
	while (parent[val].first != parent[parent[val].first].first)
	{
		if (parent[val].second)
		{
			parent[val] = parent[parent[val].first];
		}
		else
		{
			parent[val] = parent[parent[val].first];
			parent[val].second = !parent[val].second;
		}
	}
	return parent[val];
}

std::pair<size_t, bool> find(std::vector<std::pair<size_t, bool>>& parent, std::pair<size_t, bool> val)
{
	auto result = find(parent, val.first);
	if (!val.second) result.second = !result.second;
	return result;
}

void merge(std::vector<std::pair<size_t, bool>>& parent, size_t left, size_t right, bool fw)
{
	std::pair<size_t, bool> leftp = find(parent, left);
	std::pair<size_t, bool> rightp = find(parent, right);
	if (leftp.first == rightp.first)
	{
		assert(leftp.second ^ rightp.second ^ fw);
		return;
	}
	parent[rightp.first] = std::make_pair(leftp.first, rightp.second ^ leftp.second ^ fw);
}

void merge(std::vector<std::pair<size_t, bool>>& parent, std::pair<size_t, bool> left, std::pair<size_t, bool> right, bool fw)
{
	std::pair<size_t, bool> leftp = find(parent, left);
	std::pair<size_t, bool> rightp = find(parent, right);
	if (leftp.first == rightp.first)
	{
		assert(leftp.second ^ rightp.second ^ fw);
		return;
	}
	parent[rightp.first] = std::make_pair(leftp.first, rightp.second ^ leftp.second ^ fw);
}
