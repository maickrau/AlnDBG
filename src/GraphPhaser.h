#ifndef GraphPhaser_h
#define GraphPhaser_h

#include <vector>
#include <tuple>
#include "UnitigGraph.h"
#include "MatchGroup.h"

class HaplotypeInformativeSite
{
public:
	size_t readPos;
	size_t bubble;
	size_t allele;
};

std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphDiploidMEC(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const size_t numThreads, const double approxOneHapCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphPolyploidTransitiveClosure(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const size_t numThreads, const double approxOneHapCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphChainmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t resolveLength);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphHapmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t resolveLength);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLocalUniqmers(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t uniqSpanLength, const size_t resolveLength, const size_t maxCopyCount);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphLocalUniqmersLocation(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t uniqSpanLength, const size_t maxCopyCount);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> connectChainGaps(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t minSafeCoverage, const size_t maxSpuriousCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> resolveSpannedTangles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> popHaploidChainBubbles(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage);
std::vector<std::vector<HaplotypeInformativeSite>> getHaplotypeInformativeForks(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage);
std::pair<UnitigGraph, std::vector<ReadPathBundle>> unzipGraphHaploforks(const UnitigGraph& unitigGraph, const std::vector<ReadPathBundle>& readPaths, const size_t graphk, const double approxOneHapCoverage, const size_t resolveLength, const bool cheapUnzip);

#endif
