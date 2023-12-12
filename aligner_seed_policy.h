#ifndef ALIGNER_SEED_POLICY_H_
#define ALIGNER_SEED_POLICY_H_
#include "scoring.h"
#include "simple_func.h"
#define DEFAULT_SEEDMMS 0
#define DEFAULT_SEEDLEN 22
#define DEFAULT_LOCAL_SEEDLEN 20
#define DEFAULT_IVAL SIMPLE_FUNC_SQRT
#define DEFAULT_IVAL_A 1.15f
#define DEFAULT_IVAL_B 0.0f
#define DEFAULT_UNGAPPED_HITS 6
class SeedAlignmentPolicy
{
public:
	static void parseString(
		const std::string &s,
		bool local,
		bool noisyHpolymer,
		bool ignoreQuals,
		int &bonusMatchType,
		int &bonusMatch,
		int &penMmcType,
		int &penMmcMax,
		int &penMmcMin,
		int &penNType,
		int &penN,
		int &penRdExConst,
		int &penRfExConst,
		int &penRdExLinear,
		int &penRfExLinear,
		SimpleFunc &costMin,
		SimpleFunc &nCeil,
		bool &nCatPair,
		int &multiseedMms,
		int &multiseedLen,
		SimpleFunc &multiseedIval,
		size_t &failStreak,
		size_t &seedRounds);
};
extern int gDefaultSeedLen;
#endif
