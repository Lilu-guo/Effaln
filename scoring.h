#ifndef SCORING_H_
#define SCORING_H_
#include <limits>
#include "qual.h"
#include "simple_func.h"
#define DEFAULT_MATCH_BONUS_TYPE COST_MODEL_CONSTANT
#define DEFAULT_MATCH_BONUS 0
#define DEFAULT_MATCH_BONUS_TYPE_LOCAL COST_MODEL_CONSTANT
#define DEFAULT_MATCH_BONUS_LOCAL 2
#define DEFAULT_MM_PENALTY_TYPE COST_MODEL_QUAL
#define DEFAULT_MM_PENALTY_TYPE_IGNORE_QUALS COST_MODEL_CONSTANT
#define DEFAULT_MM_PENALTY_MAX 6
#define DEFAULT_MM_PENALTY_MIN 2
#define DEFAULT_N_PENALTY_TYPE COST_MODEL_CONSTANT
#define DEFAULT_N_PENALTY 1
#define DEFAULT_MIN_CONST (-0.4f)
#define DEFAULT_MIN_LINEAR (-0.4f)
#define DEFAULT_MIN_CONST_LOCAL (20.0f)
#define DEFAULT_MIN_LINEAR_LOCAL (8.0f)
#define DEFAULT_N_CEIL_CONST 0.0f
#define DEFAULT_N_CEIL_LINEAR 0.15f
#define DEFAULT_N_CAT_PAIR false
#define DEFAULT_READ_GAP_CONST 5
#define DEFAULT_READ_GAP_LINEAR 3
#define DEFAULT_READ_GAP_CONST_BADHPOLY 3
#define DEFAULT_READ_GAP_LINEAR_BADHPOLY 1
#define DEFAULT_REF_GAP_CONST 5
#define DEFAULT_REF_GAP_LINEAR 3
#define DEFAULT_REF_GAP_CONST_BADHPOLY 3
#define DEFAULT_REF_GAP_LINEAR_BADHPOLY 1
enum
{
	COST_MODEL_ROUNDED_QUAL = 1,
	COST_MODEL_QUAL,
	COST_MODEL_CONSTANT
};
class Scoring
{
	template <typename T>
	void initPens(
		T *pens,
		int type, int consMin, int consMax)
	{
		if (type == COST_MODEL_ROUNDED_QUAL)
		{
			for (int i = 0; i < 256; i++)
			{
				pens[i] = (T)qualRounds[i];
			}
		}
		else if (type == COST_MODEL_QUAL)
		{
			assert_neq(consMin, 0);
			assert_neq(consMax, 0);
			for (int i = 0; i < 256; i++)
			{
				int ii = min(i, 40);
				float frac = (float)ii / 40.0f;
				pens[i] = consMin + (T)(frac * (consMax - consMin));
				assert_gt(pens[i], 0);
			}
		}
		else if (type == COST_MODEL_CONSTANT)
		{
			for (int i = 0; i < 256; i++)
			{
				pens[i] = (T)consMax;
			}
		}
		else
		{
			throw 1;
		}
	}
public:
	Scoring(
		int mat,
		int mmcType, int mmpMax_, int mmpMin_, const SimpleFunc &scoreMin_, const SimpleFunc &nCeil_, int nType, int n, bool ncat, int rdGpConst, int rfGpConst, int rdGpLinear, int rfGpLinear, int gapbar_)
	{
		matchType = COST_MODEL_CONSTANT;
		matchConst = mat;
		mmcostType = mmcType;
		mmpMax = mmpMax_;
		mmpMin = mmpMin_;
		scoreMin = scoreMin_;
		nCeil = nCeil_;
		npenType = nType;
		npen = n;
		ncatpair = ncat;
		rdGapConst = rdGpConst;
		rfGapConst = rfGpConst;
		rdGapLinear = rdGpLinear;
		rfGapLinear = rfGpLinear;
		qualsMatter_ = mmcostType != COST_MODEL_CONSTANT;
		gapbar = gapbar_;
		monotone = matchType == COST_MODEL_CONSTANT && matchConst == 0;
		initPens<int>(mmpens, mmcostType, mmpMin_, mmpMax_);
		initPens<int>(npens, npenType, npen, npen);
		initPens<float>(matchBonuses, matchType, matchConst, matchConst);
		assert(repOk());
	}
	void setMatchBonus(int bonus)
	{
		matchType = COST_MODEL_CONSTANT;
		matchConst = bonus;
		initPens<float>(matchBonuses, matchType, matchConst, matchConst);
		assert(repOk());
	}
	void setMmPen(int mmType_, int mmpMax_, int mmpMin_)
	{
		mmcostType = mmType_;
		mmpMax = mmpMax_;
		mmpMin = mmpMin_;
		initPens<int>(mmpens, mmcostType, mmpMin, mmpMax);
	}
	void setNPen(int nType, int n)
	{
		npenType = nType;
		npen = n;
		initPens<int>(npens, npenType, npen, npen);
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_geq(matchConst, 0);
		assert_gt(rdGapConst, 0);
		assert_gt(rdGapLinear, 0);
		assert_gt(rfGapConst, 0);
		assert_gt(rfGapLinear, 0);
		return true;
	}
#endif
	static float linearFunc(int64_t x, float cnst, float lin)
	{
		return (float)((double)cnst + ((double)lin * x));
	}
	inline int mm(int rdc, int refm, int q) const
	{
		assert_range(0, 255, q);
		return (rdc > 3 || refm > 15) ? npens[q] : mmpens[q];
	}
	inline int score(int rdc, int refm, int q) const
	{
		assert_range(0, 255, q);
		if (rdc > 3 || refm > 15)
		{
			return -npens[q];
		}
		if ((refm & (1 << rdc)) != 0)
		{
			return (int)matchBonuses[q];
		}
		else
		{
			return -mmpens[q];
		}
	}
	inline int score(int rdc, int refm, int q, int &ns) const
	{
		assert_range(0, 255, q);
		if (rdc > 3 || refm > 15)
		{
			ns++;
			return -npens[q];
		}
		if ((refm & (1 << rdc)) != 0)
		{
			return (int)matchBonuses[q];
		}
		else
		{
			return -mmpens[q];
		}
	}
	inline int mm(int rdc, int q) const
	{
		assert_range(0, 255, q);
		return (rdc > 3) ? npens[q] : mmpens[q];
	}
	inline int mm(int q) const
	{
		assert_geq(q, 0);
		return q < 255 ? mmpens[q] : mmpens[255];
	}
	inline int64_t match() const
	{
		return match(30);
	}
	inline int64_t match(int q) const
	{
		assert_geq(q, 0);
		return (int64_t)((q < 255 ? matchBonuses[q] : matchBonuses[255]) + 0.5f);
	}
	inline int64_t perfectScore(size_t rdlen) const
	{
		if (monotone)
		{
			return 0;
		}
		else
		{
			return rdlen * match(30);
		}
	}
	inline bool qualitiesMatter() const { return qualsMatter_; }
	inline int n(int q) const
	{
		assert_geq(q, 0);
		return q < 255 ? npens[q] : npens[255];
	}
	inline int ins(int ext) const
	{
		assert_geq(ext, 0);
		if (ext == 0)
			return readGapOpen();
		return readGapExtend();
	}
	inline int del(int ext) const
	{
		assert_geq(ext, 0);
		if (ext == 0)
			return refGapOpen();
		return refGapExtend();
	}
	bool scoreFilter(
		int64_t minsc,
		size_t rdlen) const;
	int maxReadGaps(
		int64_t minsc,
		size_t rdlen) const;
	int maxRefGaps(
		int64_t minsc,
		size_t rdlen) const;
	bool nFilter(const BTDnaString &rd, size_t &ns) const;
	void nFilterPair(
		const BTDnaString *rd1,
		const BTDnaString *rd2, size_t &ns1, size_t &ns2, bool &filt1, bool &filt2) const;
	inline int readGapOpen() const
	{
		return rdGapConst + rdGapLinear;
	}
	inline int refGapOpen() const
	{
		return rfGapConst + rfGapLinear;
	}
	inline int readGapExtend() const
	{
		return rdGapLinear;
	}
	inline int refGapExtend() const
	{
		return rfGapLinear;
	}
	int matchType;
	int matchConst;
	int mmcostType;
	int mmpMax;
	int mmpMin;
	SimpleFunc scoreMin;
	SimpleFunc nCeil;
	int npenType;
	int npen;
	bool ncatpair;
	int rdGapConst;
	int rfGapConst;
	int rdGapLinear;
	int rfGapLinear;
	int gapbar;
	bool monotone;
	float matchBonuses[256];
	int mmpens[256];
	int npens[256];
	static Scoring base1()
	{
		const double DMAX = std::numeric_limits<double>::max();
		SimpleFunc scoreMin(SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 37.0f, 0.3f);
		SimpleFunc nCeil(SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 2.0f, 0.1f);
		return Scoring(
			1,
			COST_MODEL_CONSTANT, 3, 3, scoreMin, nCeil, COST_MODEL_CONSTANT, 3, false, 11, 11, 4, 4, 5);
	}
protected:
	bool qualsMatter_;
};
#endif
