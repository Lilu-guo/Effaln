#ifndef UNIQUE_H_
#define UNIQUE_H_
#include <string>
#include "aligner_result.h"
#include "simple_func.h"
#include "util.h"
#include "scoring.h"
typedef int64_t TMapq;
class Uniqueness
{
public:
	static bool bestIsUnique(
		const AlnSetSumm &s,
		const AlnFlags &flags,
		bool mate1,
		size_t rdlen,
		size_t ordlen,
		char *inps)
	{
		assert(!s.empty());
		return !VALID_AL_SCORE(s.bestUnchosenScore(mate1));
	}
};
class Mapq
{
public:
	virtual ~Mapq() {}
	virtual TMapq mapq(
		const AlnSetSumm &s,
		const AlnFlags &flags,
		bool mate1,
		size_t rdlen,
		size_t ordlen,
		char *inps) const = 0;
};
extern const TMapq unp_nosec_perf;
extern const TMapq unp_nosec[11];
extern const TMapq unp_sec_perf[11];
extern const TMapq unp_sec[11][11];
extern const TMapq pair_nosec_perf;
class Mapq3 : public Mapq
{
public:
	Mapq3(
		const SimpleFunc &scoreMin,
		const Scoring &sc) : scoreMin_(scoreMin),
							 sc_(sc)
	{
	}
	virtual ~Mapq3() {}
	virtual TMapq mapq(
		const AlnSetSumm &s,
		const AlnFlags &flags,
		bool mate1,
		size_t rdlen,
		size_t ordlen,
		char *inps)
		const
	{
		if (s.paired())
		{
			return pair_nosec_perf;
		}
		else
		{
			bool hasSecbest = VALID_AL_SCORE(s.bestUnchosenScore(mate1));
			if (!flags.canMax() && !s.exhausted(mate1) && !hasSecbest)
			{
				return 255;
			}
			TAlScore scMax = (TAlScore)sc_.perfectScore(rdlen);
			TAlScore scMin = scoreMin_.f<TAlScore>((float)rdlen);
			assert_geq(scMax, scMin);
			TAlScore best = scMax - s.bestScore(mate1).score();
			size_t best_bin = (size_t)((double)best * (10.0 / (double)(scMax - scMin)) + 0.5);
			assert_geq(best_bin, 0);
			assert_lt(best_bin, 11);
			if (hasSecbest)
			{
				assert_geq(s.bestScore(mate1).score(), s.bestUnchosenScore(mate1).score());
				size_t diff = s.bestScore(mate1).score() - s.bestUnchosenScore(mate1).score();
				size_t diff_bin = (size_t)((double)diff * (10.0 / (double)(scMax - scMin)) + 0.5);
				assert_geq(diff_bin, 0);
				assert_lt(diff_bin, 11);
				if (best == scMax)
				{
					return unp_sec_perf[best_bin];
				}
				else
				{
					return unp_sec[diff_bin][best_bin];
				}
			}
			else
			{
				if (best == scMax)
				{
					return unp_nosec_perf;
				}
				else
				{
					return unp_nosec[best_bin];
				}
			}
		}
	}
protected:
	SimpleFunc scoreMin_;
	const Scoring &sc_;
};
class Mapq2 : public Mapq
{
public:
	Mapq2(
		const SimpleFunc &scoreMin,
		const Scoring &sc) : scoreMin_(scoreMin),
							 sc_(sc)
	{
	}
	virtual ~Mapq2() {}
	virtual TMapq mapq(
		const AlnSetSumm &s,
		const AlnFlags &flags,
		bool mate1,
		size_t rdlen,
		size_t ordlen,
		char *inps)
		const
	{
		bool hasSecbest = s.paired() ? VALID_AL_SCORE(s.bestUnchosenCScore()) : VALID_AL_SCORE(s.bestUnchosenScore(mate1));
		if (!flags.isPrimary() ||
			(!flags.canMax() && !s.exhausted(mate1) && !hasSecbest))
		{
			return 255;
		}
		TAlScore scPer = (TAlScore)sc_.perfectScore(rdlen);
		if (s.paired())
		{
			scPer += (TAlScore)sc_.perfectScore(ordlen);
		}
		TAlScore scMin = scoreMin_.f<TAlScore>((float)rdlen);
		if (s.paired())
		{
			scMin += scoreMin_.f<TAlScore>((float)ordlen);
		}
		TAlScore secbest = scMin - 1;
		TAlScore diff = (scPer - scMin);
		TMapq ret = 0;
		TAlScore best = s.paired() ? s.bestCScore().score() : s.bestScore(mate1).score();
		TAlScore bestOver = best - scMin;
		if (sc_.monotone)
		{
			if (!hasSecbest)
			{
				if (bestOver >= diff * (double)0.8f)
					ret = 42;
				else if (bestOver >= diff * (double)0.7f)
					ret = 40;
				else if (bestOver >= diff * (double)0.6f)
					ret = 24;
				else if (bestOver >= diff * (double)0.5f)
					ret = 23;
				else if (bestOver >= diff * (double)0.4f)
					ret = 8;
				else if (bestOver >= diff * (double)0.3f)
					ret = 3;
				else
					ret = 0;
			}
			else
			{
				secbest = s.paired() ? s.bestUnchosenCScore().score() : s.bestUnchosenScore(mate1).score();
				TAlScore bestdiff = abs(abs(static_cast<long>(best)) - abs(static_cast<long>(secbest)));
				if (bestdiff >= diff * (double)0.9f)
				{
					if (bestOver == diff)
					{
						ret = 39;
					}
					else
					{
						ret = 33;
					}
				}
				else if (bestdiff >= diff * (double)0.8f)
				{
					if (bestOver == diff)
					{
						ret = 38;
					}
					else
					{
						ret = 27;
					}
				}
				else if (bestdiff >= diff * (double)0.7f)
				{
					if (bestOver == diff)
					{
						ret = 37;
					}
					else
					{
						ret = 26;
					}
				}
				else if (bestdiff >= diff * (double)0.6f)
				{
					if (bestOver == diff)
					{
						ret = 36;
					}
					else
					{
						ret = 22;
					}
				}
				else if (bestdiff >= diff * (double)0.5f)
				{
					if (bestOver == diff)
					{
						ret = 35;
					}
					else if (bestOver >= diff * (double)0.84f)
					{
						ret = 25;
					}
					else if (bestOver >= diff * (double)0.68f)
					{
						ret = 16;
					}
					else
					{
						ret = 5;
					}
				}
				else if (bestdiff >= diff * (double)0.4f)
				{
					if (bestOver == diff)
					{
						ret = 34;
					}
					else if (bestOver >= diff * (double)0.84f)
					{
						ret = 21;
					}
					else if (bestOver >= diff * (double)0.68f)
					{
						ret = 14;
					}
					else
					{
						ret = 4;
					}
				}
				else if (bestdiff >= diff * (double)0.3f)
				{
					if (bestOver == diff)
					{
						ret = 32;
					}
					else if (bestOver >= diff * (double)0.88f)
					{
						ret = 18;
					}
					else if (bestOver >= diff * (double)0.67f)
					{
						ret = 15;
					}
					else
					{
						ret = 3;
					}
				}
				else if (bestdiff >= diff * (double)0.2f)
				{
					if (bestOver == diff)
					{
						ret = 31;
					}
					else if (bestOver >= diff * (double)0.88f)
					{
						ret = 17;
					}
					else if (bestOver >= diff * (double)0.67f)
					{
						ret = 11;
					}
					else
					{
						ret = 0;
					}
				}
				else if (bestdiff >= diff * (double)0.1f)
				{
					if (bestOver == diff)
					{
						ret = 30;
					}
					else if (bestOver >= diff * (double)0.88f)
					{
						ret = 12;
					}
					else if (bestOver >= diff * (double)0.67f)
					{
						ret = 7;
					}
					else
					{
						ret = 0;
					}
				}
				else if (bestdiff > 0)
				{
					if (bestOver >= diff * (double)0.67f)
					{
						ret = 6;
					}
					else
					{
						ret = 2;
					}
				}
				else
				{
					assert_eq(bestdiff, 0);
					if (bestOver >= diff * (double)0.67f)
					{
						ret = 1;
					}
					else
					{
						ret = 0;
					}
				}
			}
		}
		else
		{
			if (!hasSecbest)
			{
				if (bestOver >= diff * (double)0.8f)
					ret = 44;
				else if (bestOver >= diff * (double)0.7f)
					ret = 42;
				else if (bestOver >= diff * (double)0.6f)
					ret = 41;
				else if (bestOver >= diff * (double)0.5f)
					ret = 36;
				else if (bestOver >= diff * (double)0.4f)
					ret = 28;
				else if (bestOver >= diff * (double)0.3f)
					ret = 24;
				else
					ret = 22;
			}
			else
			{
				secbest = s.paired() ? s.bestUnchosenCScore().score() : s.bestUnchosenScore(mate1).score();
				TAlScore bestdiff = abs(abs(static_cast<long>(best)) - abs(static_cast<long>(secbest)));
				if (bestdiff >= diff * (double)0.9f)
					ret = 40;
				else if (bestdiff >= diff * (double)0.8f)
					ret = 39;
				else if (bestdiff >= diff * (double)0.7f)
					ret = 38;
				else if (bestdiff >= diff * (double)0.6f)
					ret = 37;
				else if (bestdiff >= diff * (double)0.5f)
				{
					if (bestOver == diff)
						ret = 35;
					else if (bestOver >= diff * (double)0.50f)
						ret = 25;
					else
						ret = 20;
				}
				else if (bestdiff >= diff * (double)0.4f)
				{
					if (bestOver == diff)
						ret = 34;
					else if (bestOver >= diff * (double)0.50f)
						ret = 21;
					else
						ret = 19;
				}
				else if (bestdiff >= diff * (double)0.3f)
				{
					if (bestOver == diff)
						ret = 33;
					else if (bestOver >= diff * (double)0.5f)
						ret = 18;
					else
						ret = 16;
				}
				else if (bestdiff >= diff * (double)0.2f)
				{
					if (bestOver == diff)
						ret = 32;
					else if (bestOver >= diff * (double)0.5f)
						ret = 17;
					else
						ret = 12;
				}
				else if (bestdiff >= diff * (double)0.1f)
				{
					if (bestOver == diff)
						ret = 31;
					else if (bestOver >= diff * (double)0.5f)
						ret = 14;
					else
						ret = 9;
				}
				else if (bestdiff > 0)
				{
					if (bestOver >= diff * (double)0.5f)
						ret = 11;
					else
						ret = 2;
				}
				else
				{
					assert_eq(bestdiff, 0);
					if (bestOver >= diff * (double)0.5f)
						ret = 1;
					else
						ret = 0;
				}
			}
		}
		return ret;
	}
protected:
	SimpleFunc scoreMin_;
	const Scoring &sc_;
};
class eMapq : public Mapq
{
public:
	eMapq(
		const SimpleFunc &scoreMin,
		const Scoring &sc) : scoreMin_(scoreMin),
							 sc_(sc)
	{
	}
	virtual ~eMapq() {}
	virtual TMapq mapq(
		const AlnSetSumm &s,
		const AlnFlags &flags,
		bool mate1,
		size_t rdlen,
		size_t ordlen,
		char *inps)
		const
	{
		bool hasSecbest = VALID_AL_SCORE(s.bestUnchosenScore(mate1));
		if (!flags.canMax() && !s.exhausted(mate1) && !hasSecbest)
		{
			return 255;
		}
		TAlScore scPer = (TAlScore)sc_.perfectScore(rdlen);
		TAlScore scMin = scoreMin_.f<TAlScore>((float)rdlen);
		TAlScore secbest = scMin - 1;
		TAlScore diff = (scPer - scMin);
		float sixth_2 = (float)(scPer - diff * (double)0.1666f * 2);
		float sixth_3 = (float)(scPer - diff * (double)0.1666f * 3);
		TMapq ret = 0;
		TAlScore best = s.bestScore(mate1).score();
		if (!hasSecbest)
		{
			if (best >= sixth_2)
			{
				ret = 37;
			}
			else if (best >= sixth_3)
			{
				ret = 25;
			}
			else
			{
				ret = 10;
			}
		}
		else
		{
			secbest = s.bestUnchosenScore(mate1).score();
			TAlScore bestdiff = abs(abs(static_cast<long>(best)) - abs(static_cast<long>(secbest)));
			if (bestdiff >= diff * 0.1666 * 5)
			{
				ret = 6;
			}
			else if (bestdiff >= diff * 0.1666 * 4)
			{
				ret = 5;
			}
			else if (bestdiff >= diff * 0.1666 * 3)
			{
				ret = 4;
			}
			else if (bestdiff >= diff * 0.1666 * 2)
			{
				ret = 3;
			}
			else if (bestdiff >= diff * 0.1666 * 1)
			{
				ret = 2;
			}
			else
			{
				ret = 1;
			}
		}
		return ret;
	}
protected:
	SimpleFunc scoreMin_;
	const Scoring &sc_;
};
static inline Mapq *new_mapq(
	int version,
	const SimpleFunc &scoreMin,
	const Scoring &sc)
{
	if (version == 3)
	{
		return new Mapq3(scoreMin, sc);
	}
	else if (version == 2)
	{
		return new Mapq2(scoreMin, sc);
	}
	else
	{
		return new eMapq(scoreMin, sc);
	}
}
#endif
