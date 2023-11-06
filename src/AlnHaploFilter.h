#ifndef AlnHaploFilter_h
#define AlnHaploFilter_h

#include <vector>
#include "MatchGroup.h"
#include "TwobitString.h"

std::vector<bool> getValidAlignments(const std::vector<MatchGroup>& matches, const std::vector<size_t>& readLengths, const size_t minHaploCoverage, const size_t maxCrossCoverage);

#endif
