#include <cassert>
#include "TwobitString.h"

TwobitString::TwobitString() :
	realSize(0),
	bits()
{
}

TwobitString::TwobitString(const std::string& str)
{
	resize(str.size());
	for (size_t i = 0; i < str.size(); i++)
	{
		switch(str[i])
		{
			case 'a':
			case 'A':
				set(i, 0);
				break;
			case 'c':
			case 'C':
				set(i, 1);
				break;
			case 'g':
			case 'G':
				set(i, 2);
				break;
			case 't':
			case 'T':
				set(i, 3);
				break;
			default:
				assert(false);
				break;
		}
	}
}

std::string TwobitString::toString() const
{
	std::string result;
	result.resize(size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = "ACGT"[get(i)];
	}
	return result;
}

void TwobitString::resize(size_t size)
{
	realSize = size;
	bits.resize((realSize+3)/4);
}

uint8_t TwobitString::get(size_t i) const
{
	size_t index = i / 4;
	size_t offset = (i % 4) * 2;
	return (bits[index] >> offset) & 3;
}

void TwobitString::set(size_t i, uint8_t v)
{
	assert(v < 4);
	size_t index = i / 4;
	size_t offset = (i % 4) * 2;
	uint8_t removeMask = ~(3 << offset);
	bits[index] &= removeMask;
	uint8_t addMask = (uint8_t)v << offset;
	bits[index] |= addMask;
}

void TwobitString::emplace_back(uint8_t v)
{
	if (realSize == 0) bits.emplace_back(0);
	realSize += 1;
	if (realSize % 4 == 0) bits.emplace_back(0);
	set(realSize-1, v);
}

size_t TwobitString::size() const
{
	return realSize;
}

std::string TwobitString::substr(size_t start) const
{
	std::string result;
	result.reserve(size() - start);
	for (size_t i = start; i < size(); i++)
	{
		result.push_back("ACGT"[get(i)]);
	}
	return result;
}

std::string TwobitString::substr(size_t start, size_t size) const
{
	std::string result;
	result.reserve(size);
	for (size_t i = 0; i < size; i++)
	{
		result.push_back("ACGT"[get(start+i)]);
	}
	return result;
}

bool TwobitString::operator==(const TwobitString& other) const
{
	if (other.size() != size()) return false;
	for (size_t i = 0; i < size(); i++)
	{
		if (get(i) != other.get(i)) return false;
	}
	return true;
}

bool TwobitString::operator<(const TwobitString& other) const
{
	if (size() < other.size()) return true;
	if (size() > other.size()) return false;
	for (size_t i = 0; i < size(); i++)
	{
		if (get(i) < other.get(i)) return true;
		if (get(i) > other.get(i)) return false;
	}
	return false;
}

void TwobitString::pop_back()
{
	assert(realSize >= 1);
	realSize -= 1;
}

TwobitString& TwobitString::operator+=(const std::string& str)
{
	for (size_t i = 0; i < str.size(); i++)
	{
		switch(str[i])
		{
		case 'a':
		case 'A':
			emplace_back(0);
			break;
		case 'c':
		case 'C':
			emplace_back(1);
			break;
		case 'g':
		case 'G':
			emplace_back(2);
			break;
		case 't':
		case 'T':
			emplace_back(3);
			break;
		default:
			assert(false);
		}
	}
	return *this;
}
