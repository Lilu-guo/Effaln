#ifndef ALIGNER_REPORT_H_
#define ALIGNER_REPORT_H_
#include "aligner_cache.h"
class Reporter
{
public:
	bool report(const AlignmentCacheIface &cache, const QVal &qv)
	{
		return true;
	}
};
#endif
