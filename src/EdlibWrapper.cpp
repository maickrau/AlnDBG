#include "edlib.h"
#include "EdlibWrapper.h"

size_t getNumMismatches(const std::string& leftSequence, const std::string& rightSequence, const size_t maxMismatches)
{
	EdlibAlignResult result = edlibAlign(leftSequence.data(), leftSequence.size(), rightSequence.data(), rightSequence.size(), edlibNewAlignConfig(maxMismatches+1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	size_t mismatches = maxMismatches+1;
	if (result.status == EDLIB_STATUS_OK)
	{
		mismatches = result.editDistance;
	}
	edlibFreeAlignResult(result);
	return mismatches;
}
