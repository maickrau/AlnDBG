#ifndef MatchGroup_h
#define MatchGroup_h

#include <vector>
#include <cstdio>
#include <cstdint>

class MatchGroup
{
public:
	class Match
	{
	public:
		uint16_t leftStart;
		uint16_t rightStart;
		uint16_t length;
	};
	size_t leftRead;
	size_t rightRead;
	bool rightFw;
	size_t leftStart;
	size_t leftEnd;
	size_t rightStart;
	size_t rightEnd;
	std::vector<Match> matches;
};

#endif
