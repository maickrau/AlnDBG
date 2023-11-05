#ifndef KmerMatcher_h
#define KmerMatcher_h

#include <vector>
#include "MatchGroup.h"
#include "TwobitString.h"

void addKmerMatches(const size_t numThreads, const std::vector<TwobitString>& readSequences, std::vector<MatchGroup>& matches, const size_t graphk, const size_t graphd);

#endif
