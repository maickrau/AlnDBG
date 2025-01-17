#include <algorithm>
#include "TransitiveClosure.h"

std::vector<std::vector<int>> getBinCheckOrder(const int maxDistance)
{
	std::vector<std::vector<int>> result;
	for (int i = -maxDistance; i <= maxDistance; i++)
	{
		for (int j = -maxDistance; j <= maxDistance; j++)
		{
			for (int k = -maxDistance; k <= maxDistance; k++)
			{
				result.emplace_back();
				result.back().emplace_back(i);
				result.back().emplace_back(j);
				result.back().emplace_back(k);
			}
		}
	}
	assert(result.size() == pow(2*maxDistance+1, 3));
	std::sort(result.begin(), result.end(), [](const std::vector<int>& left, const std::vector<int>& right)
	{
		assert(left.size() == right.size());
		size_t leftDistance = 0;
		size_t rightDistance = 0;
		for (size_t i = 0; i < left.size(); i++)
		{
			leftDistance += left[i] < 0 ? -left[i] : left[i];
			rightDistance += right[i] < 0 ? -right[i] : right[i];
		}
		return leftDistance < rightDistance;
	});
	return result;
}

std::vector<std::vector<int>> getBinCheckOrderNotSymmetric(const int maxDistance)
{
	std::vector<std::vector<int>> result;
	for (int i = 0; i <= maxDistance; i++)
	{
		for (int j = -maxDistance; j <= maxDistance; j++)
		{
			for (int k = -maxDistance; k <= maxDistance; k++)
			{
				result.emplace_back();
				result.back().emplace_back(i);
				result.back().emplace_back(j);
				result.back().emplace_back(k);
			}
		}
	}
	assert(result.size() == pow(2*maxDistance+1, 2)*(maxDistance+1));
	std::sort(result.begin(), result.end(), [](const std::vector<int>& left, const std::vector<int>& right)
	{
		assert(left.size() == right.size());
		size_t leftDistance = 0;
		size_t rightDistance = 0;
		for (size_t i = 0; i < left.size(); i++)
		{
			leftDistance += left[i] < 0 ? -left[i] : left[i];
			rightDistance += right[i] < 0 ? -right[i] : right[i];
		}
		return leftDistance < rightDistance;
	});
	return result;
}