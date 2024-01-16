#ifndef KmerMatcher_h
#define KmerMatcher_h

#include <vector>
#include "MatchGroup.h"
#include "TwobitString.h"

std::vector<MatchGroup> addKmerMatches(const size_t numThreads, const std::vector<TwobitString>& readSequences, const std::vector<MatchGroup>& matches, const size_t graphk, const size_t graphd);

#endif
