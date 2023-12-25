#ifndef GraphPhaser_h
#define GraphPhaser_h

#include <vector>
#include <tuple>
#include "UnitigGraph.h"
#include "MatchGroup.h"

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLinearizable(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphChainmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t resolveLength);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphHapmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t resolveLength);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLocalUniqmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t uniqSpanLength, const size_t resolveLength);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> connectChainGaps(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage, const size_t minSafeCoverage, const size_t maxSpuriousCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolveSpannedTangles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> popHaploidChainBubbles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const double approxOneHapCoverage);

#endif
