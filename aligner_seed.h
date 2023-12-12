#ifndef ALIGNER_SEED_H_
#define ALIGNER_SEED_H_
#include <iostream>
#include <utility>
#include <limits>
#include "qual.h"
#include "ds.h"
#include "sstring.h"
#include "alphabet.h"
#include "edit.h"
#include "read.h"
#include "threading.h"
#include "aligner_result.h"
#include "aligner_cache.h"
#include "scoring.h"
#include "mem_ids.h"
#include "simple_func.h"
#include "btypes.h"
struct Constraint
{
	Constraint() { init(); }
	void init()
	{
		edits = mms = ins = dels = penalty = editsCeil = mmsCeil =
			insCeil = delsCeil = penaltyCeil = MAX_I;
		penFunc.reset();
		instantiated = false;
	}
	bool mustMatch()
	{
		assert(instantiated);
		return (mms == 0 && edits == 0) ||
			   penalty == 0 ||
			   (mms == 0 && dels == 0 && ins == 0);
	}
	bool canMismatch(int q, const Scoring &cm)
	{
		assert(instantiated);
		return (mms > 0 || edits > 0) &&
			   penalty >= cm.mm(q);
	}
	bool canN(int q, const Scoring &cm)
	{
		assert(instantiated);
		return (mms > 0 || edits > 0) &&
			   penalty >= cm.n(q);
	}
	bool canMismatch()
	{
		assert(instantiated);
		return (mms > 0 || edits > 0) && penalty > 0;
	}
	bool canN()
	{
		assert(instantiated);
		return (mms > 0 || edits > 0);
	}
	bool canDelete(int ex, const Scoring &cm)
	{
		assert(instantiated);
		return (dels > 0 && edits > 0) &&
			   penalty >= cm.del(ex);
	}
	bool canDelete()
	{
		assert(instantiated);
		return (dels > 0 || edits > 0) &&
			   penalty > 0;
	}
	bool canInsert(int ex, const Scoring &cm)
	{
		assert(instantiated);
		return (ins > 0 || edits > 0) &&
			   penalty >= cm.ins(ex);
	}
	bool canInsert()
	{
		assert(instantiated);
		return (ins > 0 || edits > 0) &&
			   penalty > 0;
	}
	bool canGap()
	{
		assert(instantiated);
		return ((ins > 0 || dels > 0) || edits > 0) && penalty > 0;
	}
	void chargeMismatch(int q, const Scoring &cm)
	{
		assert(instantiated);
		if (mms == 0)
		{
			assert_gt(edits, 0);
			edits--;
		}
		else
			mms--;
		penalty -= cm.mm(q);
		assert_geq(mms, 0);
		assert_geq(edits, 0);
		assert_geq(penalty, 0);
	}
	void chargeN(int q, const Scoring &cm)
	{
		assert(instantiated);
		if (mms == 0)
		{
			assert_gt(edits, 0);
			edits--;
		}
		else
			mms--;
		penalty -= cm.n(q);
		assert_geq(mms, 0);
		assert_geq(edits, 0);
		assert_geq(penalty, 0);
	}
	void chargeDelete(int ex, const Scoring &cm)
	{
		assert(instantiated);
		dels--;
		edits--;
		penalty -= cm.del(ex);
		assert_geq(dels, 0);
		assert_geq(edits, 0);
		assert_geq(penalty, 0);
	}
	void chargeInsert(int ex, const Scoring &cm)
	{
		assert(instantiated);
		ins--;
		edits--;
		penalty -= cm.ins(ex);
		assert_geq(ins, 0);
		assert_geq(edits, 0);
		assert_geq(penalty, 0);
	}
	bool acceptable()
	{
		assert(instantiated);
		return edits <= editsCeil &&
			   mms <= mmsCeil &&
			   ins <= insCeil &&
			   dels <= delsCeil &&
			   penalty <= penaltyCeil;
	}
	static int instantiate(size_t rdlen, const SimpleFunc &func)
	{
		return func.f<int>((double)rdlen);
	}
	void instantiate(size_t rdlen)
	{
		assert(!instantiated);
		if (penFunc.initialized())
		{
			penalty = Constraint::instantiate(rdlen, penFunc);
		}
		instantiated = true;
	}
	int edits;
	int mms;
	int ins;
	int dels;
	int penalty;
	int editsCeil;
	int mmsCeil;
	int insCeil;
	int delsCeil;
	int penaltyCeil;
	SimpleFunc penFunc;
	bool instantiated;
	static Constraint exact();
	static Constraint penaltyBased(int pen);
	static Constraint penaltyFuncBased(const SimpleFunc &func);
	static Constraint mmBased(int mms);
	static Constraint editBased(int edits);
};
enum
{
	SEED_TYPE_EXACT = 1,
	SEED_TYPE_LEFT_TO_RIGHT,
	SEED_TYPE_RIGHT_TO_LEFT,
	SEED_TYPE_INSIDE_OUT
};
struct InstantiatedSeed;
struct Seed
{
	int len;
	int type;
	Constraint *overall;
	Seed() { init(0, 0, NULL); }
	Seed(int ln, int ty, Constraint *oc)
	{
		init(ln, ty, oc);
	}
	void init(int ln, int ty, Constraint *oc)
	{
		len = ln;
		type = ty;
		overall = oc;
	}
	Constraint zones[3];
	bool acceptable()
	{
		assert(overall != NULL);
		return zones[0].acceptable() &&
			   zones[1].acceptable() &&
			   zones[2].acceptable() &&
			   overall->acceptable();
	}
	bool instantiate(
		const Read &read,
		const BTDnaString &seq,
		const BTString &qual,
		const Scoring &pens,
		int depth,
		int seedoffidx,
		int seedtypeidx,
		bool fw,
		InstantiatedSeed &si) const;
	static void mmSeeds(
		int mms,
		int ln,
		EList<Seed> &pols,
		Constraint &oall)
	{
		if (mms == 0)
		{
			zeroMmSeeds(ln, pols, oall);
		}
		else if (mms == 1)
		{
			oneMmSeeds(ln, pols, oall);
		}
		else if (mms == 2)
		{
			twoMmSeeds(ln, pols, oall);
		}
		else
			throw 1;
	}
	static void zeroMmSeeds(int ln, EList<Seed> &, Constraint &);
	static void oneMmSeeds(int ln, EList<Seed> &, Constraint &);
	static void twoMmSeeds(int ln, EList<Seed> &, Constraint &);
};
struct InstantiatedSeed
{
	InstantiatedSeed() : steps(AL_CAT), zones(AL_CAT) {}
	EList<int> steps;
	EList<pair<int, int>> zones;
	BTDnaString *seq;
	BTString *qual;
	Constraint cons[3];
	Constraint overall;
	int maxjump;
	int seedoff;
	int seedoffidx;
	int seedtypeidx;
	bool fw;
	bool nfiltered;
	Seed s;
#ifndef NDEBUG
	bool repOk() const
	{
		return true;
	}
#endif
};
struct EEHit
{
	EEHit() { reset(); }
	void reset()
	{
		top = bot = 0;
		fw = false;
		e1.reset();
		e2.reset();
		score = MIN_I64;
	}
	void init(
		TIndexOffU top_,
		TIndexOffU bot_,
		const Edit *e1_,
		const Edit *e2_,
		bool fw_,
		int64_t score_)
	{
		top = top_;
		bot = bot_;
		if (e1_ != NULL)
		{
			e1 = *e1_;
		}
		else
		{
			e1.reset();
		}
		if (e2_ != NULL)
		{
			e2 = *e2_;
		}
		else
		{
			e2.reset();
		}
		fw = fw_;
		score = score_;
	}
	int mms() const
	{
		if (e2.inited())
			return 2;
		else if (e1.inited())
			return 1;
		else
			return 0;
	}
	int ns() const
	{
		int ns = 0;
		if (e1.inited() && e1.hasN())
		{
			ns++;
			if (e2.inited() && e2.hasN())
			{
				ns++;
			}
		}
		return ns;
	}
	int refns() const
	{
		int ns = 0;
		if (e1.inited() && e1.chr == 'N')
		{
			ns++;
			if (e2.inited() && e2.chr == 'N')
			{
				ns++;
			}
		}
		return ns;
	}
	bool empty() const
	{
		return bot <= top;
	}
	bool operator<(const EEHit &o) const
	{
		return score > o.score;
	}
	TIndexOffU size() const { return bot - top; }
#ifndef NDEBUG
	bool repOk(const Read &rd) const
	{
		assert_gt(bot, top);
		if (e1.inited())
		{
			assert_lt(e1.pos, rd.length());
			if (e2.inited())
			{
				assert_lt(e2.pos, rd.length());
			}
		}
		return true;
	}
#endif
	TIndexOffU top;
	TIndexOffU bot;
	Edit e1;
	Edit e2;
	bool fw;
	int64_t score;
};
class SeedResults
{
public:
	SeedResults() : seqFw_(AL_CAT),
					seqRc_(AL_CAT),
					qualFw_(AL_CAT),
					qualRc_(AL_CAT),
					hitsFw_(AL_CAT),
					hitsRc_(AL_CAT),
					isFw_(AL_CAT),
					isRc_(AL_CAT),
					sortedFw_(AL_CAT),
					sortedRc_(AL_CAT),
					offIdx2off_(AL_CAT),
					rankOffs_(AL_CAT),
					rankFws_(AL_CAT),
					mm1Hit_(AL_CAT)
	{
		clear();
	}
	void nextRead(const Read &read)
	{
		read_ = &read;
	}
	void add(
		const QVal &qv,
		const AlignmentCache &ac,
		uint32_t seedIdx,
		bool seedFw)
	{
		assert(qv.repOk(ac));
		assert(repOk(&ac));
		assert_lt(seedIdx, hitsFw_.size());
		assert_gt(numOffs_, 0);
		if (qv.empty())
			return;
		if (seedFw)
		{
			assert(!hitsFw_[seedIdx].valid());
			hitsFw_[seedIdx] = qv;
			numEltsFw_ += qv.numElts();
			numRangesFw_ += qv.numRanges();
			if (qv.numRanges() > 0)
				nonzFw_++;
		}
		else
		{
			assert(!hitsRc_[seedIdx].valid());
			hitsRc_[seedIdx] = qv;
			numEltsRc_ += qv.numElts();
			numRangesRc_ += qv.numRanges();
			if (qv.numRanges() > 0)
				nonzRc_++;
		}
		numElts_ += qv.numElts();
		numRanges_ += qv.numRanges();
		if (qv.numRanges() > 0)
		{
			nonzTot_++;
			if (qv.numRanges() == 1 && qv.numElts() == 1)
			{
				uniTot_++;
				uniTotS_[seedFw ? 0 : 1]++;
			}
			else
			{
				repTot_++;
				repTotS_[seedFw ? 0 : 1]++;
			}
		}
		assert(repOk(&ac));
	}
	void reset(
		const Read &read,
		const EList<uint32_t> &offIdx2off,
		size_t numOffs)
	{
		assert_gt(numOffs, 0);
		clearSeeds();
		numOffs_ = numOffs;
		seqFw_.resize(numOffs_);
		seqRc_.resize(numOffs_);
		qualFw_.resize(numOffs_);
		qualRc_.resize(numOffs_);
		hitsFw_.resize(numOffs_);
		hitsRc_.resize(numOffs_);
		isFw_.resize(numOffs_);
		isRc_.resize(numOffs_);
		sortedFw_.resize(numOffs_);
		sortedRc_.resize(numOffs_);
		offIdx2off_ = offIdx2off;
		for (size_t i = 0; i < numOffs_; i++)
		{
			sortedFw_[i] = sortedRc_[i] = false;
			hitsFw_[i].reset();
			hitsRc_[i].reset();
			isFw_[i].clear();
			isRc_[i].clear();
		}
		read_ = &read;
		sorted_ = false;
	}
	void clearSeeds()
	{
		sortedFw_.clear();
		sortedRc_.clear();
		rankOffs_.clear();
		rankFws_.clear();
		offIdx2off_.clear();
		hitsFw_.clear();
		hitsRc_.clear();
		isFw_.clear();
		isRc_.clear();
		seqFw_.clear();
		seqRc_.clear();
		nonzTot_ = 0;
		uniTot_ = uniTotS_[0] = uniTotS_[1] = 0;
		repTot_ = repTotS_[0] = repTotS_[1] = 0;
		nonzFw_ = 0;
		nonzRc_ = 0;
		numOffs_ = 0;
		numRanges_ = 0;
		numElts_ = 0;
		numRangesFw_ = 0;
		numEltsFw_ = 0;
		numRangesRc_ = 0;
		numEltsRc_ = 0;
	}
	void clear()
	{
		clearSeeds();
		read_ = NULL;
		exactFwHit_.reset();
		exactRcHit_.reset();
		mm1Hit_.clear();
		mm1Sorted_ = false;
		mm1Elt_ = 0;
		assert(empty());
	}
	void toSeedAlSumm(SeedAlSumm &ssum) const
	{
		ssum.nonzTot = nonzTot_;
		ssum.nonzFw = nonzFw_;
		ssum.nonzRc = nonzRc_;
		ssum.nrangeTot = numRanges_;
		ssum.nrangeFw = numRangesFw_;
		ssum.nrangeRc = numRangesRc_;
		ssum.neltTot = numElts_;
		ssum.neltFw = numEltsFw_;
		ssum.neltRc = numEltsRc_;
		ssum.maxNonzRangeFw = ssum.minNonzRangeFw = 0;
		ssum.maxNonzRangeRc = ssum.minNonzRangeRc = 0;
		ssum.maxNonzEltFw = ssum.minNonzEltFw = 0;
		ssum.maxNonzEltRc = ssum.minNonzEltRc = 0;
		for (size_t i = 0; i < numOffs_; i++)
		{
			if (hitsFw_[i].valid())
			{
				if (ssum.minNonzEltFw == 0 || hitsFw_[i].numElts() < ssum.minNonzEltFw)
				{
					ssum.minNonzEltFw = hitsFw_[i].numElts();
				}
				if (ssum.maxNonzEltFw == 0 || hitsFw_[i].numElts() > ssum.maxNonzEltFw)
				{
					ssum.maxNonzEltFw = hitsFw_[i].numElts();
				}
				if (ssum.minNonzRangeFw == 0 || hitsFw_[i].numRanges() < ssum.minNonzRangeFw)
				{
					ssum.minNonzRangeFw = hitsFw_[i].numRanges();
				}
				if (ssum.maxNonzRangeFw == 0 || hitsFw_[i].numRanges() > ssum.maxNonzRangeFw)
				{
					ssum.maxNonzRangeFw = hitsFw_[i].numRanges();
				}
			}
			if (hitsRc_[i].valid())
			{
				if (ssum.minNonzEltRc == 0 || hitsRc_[i].numElts() < ssum.minNonzEltRc)
				{
					ssum.minNonzEltRc = hitsRc_[i].numElts();
				}
				if (ssum.maxNonzEltRc == 0 || hitsRc_[i].numElts() > ssum.maxNonzEltRc)
				{
					ssum.maxNonzEltRc = hitsRc_[i].numElts();
				}
				if (ssum.minNonzRangeRc == 0 || hitsRc_[i].numRanges() < ssum.minNonzRangeRc)
				{
					ssum.minNonzRangeRc = hitsRc_[i].numRanges();
				}
				if (ssum.maxNonzRangeRc == 0 || hitsRc_[i].numRanges() > ssum.maxNonzRangeRc)
				{
					ssum.maxNonzRangeRc = hitsRc_[i].numRanges();
				}
			}
		}
	}
	float averageHitsPerSeed() const
	{
		return nonzTot_ == 0 ? 0 : (float)numElts_ / (float)nonzTot_;
	}
	size_t numUniqueSeeds() const
	{
		return uniTot_;
	}
	size_t numUniqueSeedsStrand(bool fw) const
	{
		return uniTotS_[fw ? 0 : 1];
	}
	size_t numRepeatSeeds() const
	{
		return repTot_;
	}
	size_t numRepeatSeedsStrand(bool fw) const
	{
		return repTotS_[fw ? 0 : 1];
	}
	float medianHitsPerSeed() const
	{
		EList<size_t> &median = const_cast<EList<size_t> &>(tmpMedian_);
		median.clear();
		for (size_t i = 0; i < numOffs_; i++)
		{
			if (hitsFw_[i].valid() && hitsFw_[i].numElts() > 0)
			{
				median.push_back(hitsFw_[i].numElts());
			}
			if (hitsRc_[i].valid() && hitsRc_[i].numElts() > 0)
			{
				median.push_back(hitsRc_[i].numElts());
			}
		}
		if (tmpMedian_.empty())
		{
			return 0.0f;
		}
		median.sort();
		float med1 = (float)median[tmpMedian_.size() >> 1];
		float med2 = med1;
		if ((median.size() & 1) == 0)
		{
			med2 = (float)median[(tmpMedian_.size() >> 1) - 1];
		}
		return med1 + med2 * 0.5f;
	}
	double uniquenessFactor() const
	{
		double result = 0.0;
		for (size_t i = 0; i < numOffs_; i++)
		{
			if (hitsFw_[i].valid())
			{
				size_t nelt = hitsFw_[i].numElts();
				result += (1.0 / (double)(nelt * nelt));
			}
			if (hitsRc_[i].valid())
			{
				size_t nelt = hitsRc_[i].numElts();
				result += (1.0 / (double)(nelt * nelt));
			}
		}
		return result;
	}
	size_t numRanges() const { return numRanges_; }
	size_t numElts() const { return numElts_; }
	size_t numRangesFw() const { return numRangesFw_; }
	size_t numEltsFw() const { return numEltsFw_; }
	size_t numRangesRc() const { return numRangesRc_; }
	size_t numEltsRc() const { return numEltsRc_; }
	size_t idx2off(size_t off) const
	{
		return offIdx2off_[off];
	}
	bool empty() const { return numRanges() == 0; }
	const QVal &hitsAtOffIdx(bool fw, size_t seedoffidx) const
	{
		assert_lt(seedoffidx, numOffs_);
		assert(repOk(NULL));
		return fw ? hitsFw_[seedoffidx] : hitsRc_[seedoffidx];
	}
	EList<InstantiatedSeed> &instantiatedSeeds(bool fw, size_t seedoffidx)
	{
		assert_lt(seedoffidx, numOffs_);
		assert(repOk(NULL));
		return fw ? isFw_[seedoffidx] : isRc_[seedoffidx];
	}
	size_t numOffs() const { return numOffs_; }
	const Read &read() const { return *read_; }
#ifndef NDEBUG
	bool repOk(
		const AlignmentCache *ac,
		bool requireInited = false) const
	{
		if (requireInited)
		{
			assert(read_ != NULL);
		}
		if (numOffs_ > 0)
		{
			assert_eq(numOffs_, hitsFw_.size());
			assert_eq(numOffs_, hitsRc_.size());
			assert_leq(numRanges_, numElts_);
			assert_leq(nonzTot_, numRanges_);
			size_t nonzs = 0;
			for (int fw = 0; fw <= 1; fw++)
			{
				const EList<QVal> &rrs = (fw ? hitsFw_ : hitsRc_);
				for (size_t i = 0; i < numOffs_; i++)
				{
					if (rrs[i].valid())
					{
						if (rrs[i].numRanges() > 0)
							nonzs++;
						if (ac != NULL)
						{
							assert(rrs[i].repOk(*ac));
						}
					}
				}
			}
			assert_eq(nonzs, nonzTot_);
			assert(!sorted_ || nonzTot_ == rankFws_.size());
			assert(!sorted_ || nonzTot_ == rankOffs_.size());
		}
		return true;
	}
#endif
	void rankSeedHits(RandomSource &rnd, bool all)
	{
		if (all)
		{
			for (uint32_t i = 1; i < numOffs_; i++)
			{
				for (int fwi = 0; fwi <= 1; fwi++)
				{
					bool fw = fwi == 0;
					EList<QVal> &rrs = (fw ? hitsFw_ : hitsRc_);
					if (rrs[i].valid() && rrs[i].numElts() > 0)
					{
						rankOffs_.push_back(i);
						rankFws_.push_back(fw);
					}
				}
			}
			if (hitsFw_[0].valid() && hitsFw_[0].numElts() > 0)
			{
				rankOffs_.push_back(0);
				rankFws_.push_back(true);
			}
			if (hitsRc_[0].valid() && hitsRc_[0].numElts() > 0)
			{
				rankOffs_.push_back(0);
				rankFws_.push_back(false);
			}
		}
		else
		{
			while (rankOffs_.size() < nonzTot_)
			{
				TIndexOffU minsz = MAX_U32;
				uint32_t minidx = 0;
				bool minfw = true;
				bool rb = rnd.nextBool();
				assert(rb == 0 || rb == 1);
				for (int fwi = 0; fwi <= 1; fwi++)
				{
					bool fw = (fwi == (rb ? 1 : 0));
					EList<QVal> &rrs = (fw ? hitsFw_ : hitsRc_);
					EList<bool> &sorted = (fw ? sortedFw_ : sortedRc_);
					uint32_t i = (rnd.nextU32() % (uint32_t)numOffs_);
					for (uint32_t ii = 0; ii < numOffs_; ii++)
					{
						if (rrs[i].valid() &&
							rrs[i].numElts() > 0 &&
							!sorted[i] &&
							rrs[i].numElts() < minsz)
						{
							minsz = rrs[i].numElts();
							minidx = i;
							minfw = (fw == 1);
						}
						if ((++i) == numOffs_)
						{
							i = 0;
						}
					}
				}
				assert_neq(MAX_U32, minsz);
				if (minfw)
				{
					sortedFw_[minidx] = true;
				}
				else
				{
					sortedRc_[minidx] = true;
				}
				rankOffs_.push_back(minidx);
				rankFws_.push_back(minfw);
			}
		}
		assert_eq(rankOffs_.size(), rankFws_.size());
		sorted_ = true;
	}
	size_t nonzeroOffsets() const
	{
		assert(!sorted_ || nonzTot_ == rankFws_.size());
		assert(!sorted_ || nonzTot_ == rankOffs_.size());
		return nonzTot_;
	}
	bool allFwSeedsHit() const
	{
		return nonzFw_ == numOffs();
	}
	bool allRcSeedsHit() const
	{
		return nonzRc_ == numOffs();
	}
	size_t fewestEditsEE(bool fw, int seedlen, int per) const
	{
		assert_gt(seedlen, 0);
		assert_gt(per, 0);
		size_t nonz = fw ? nonzFw_ : nonzRc_;
		if (nonz < numOffs())
		{
			int maxdepth = (seedlen + per - 1) / per;
			int missing = (int)(numOffs() - nonz);
			return (missing + maxdepth - 1) / maxdepth;
		}
		else
		{
			return 0;
		}
	}
	size_t nonzeroOffsetsFw() const
	{
		return nonzFw_;
	}
	size_t nonzeroOffsetsRc() const
	{
		return nonzRc_;
	}
	const QVal &hitsByRank(
		size_t r,
		uint32_t &offidx,
		uint32_t &off,
		bool &fw,
		uint32_t &seedlen)
	{
		assert(sorted_);
		assert_lt(r, nonzTot_);
		if (rankFws_[r])
		{
			fw = true;
			offidx = rankOffs_[r];
			assert_lt(offidx, offIdx2off_.size());
			off = offIdx2off_[offidx];
			seedlen = (uint32_t)seqFw_[rankOffs_[r]].length();
			return hitsFw_[rankOffs_[r]];
		}
		else
		{
			fw = false;
			offidx = rankOffs_[r];
			assert_lt(offidx, offIdx2off_.size());
			off = offIdx2off_[offidx];
			seedlen = (uint32_t)seqRc_[rankOffs_[r]].length();
			return hitsRc_[rankOffs_[r]];
		}
	}
	const BTDnaString &seqByRank(size_t r)
	{
		assert(sorted_);
		assert_lt(r, nonzTot_);
		return rankFws_[r] ? seqFw_[rankOffs_[r]] : seqRc_[rankOffs_[r]];
	}
	const BTString &qualByRank(size_t r)
	{
		assert(sorted_);
		assert_lt(r, nonzTot_);
		return rankFws_[r] ? qualFw_[rankOffs_[r]] : qualRc_[rankOffs_[r]];
	}
	EList<BTDnaString> &seqs(bool fw) { return fw ? seqFw_ : seqRc_; }
	EList<BTString> &quals(bool fw) { return fw ? qualFw_ : qualRc_; }
	EEHit exactFwEEHit() const { return exactFwHit_; }
	EEHit exactRcEEHit() const { return exactRcHit_; }
	const EList<EEHit> &mm1EEHits() const { return mm1Hit_; }
	void sort1mmEe(RandomSource &rnd)
	{
		assert(!mm1Sorted_);
		mm1Hit_.sort();
		size_t streak = 0;
		for (size_t i = 1; i < mm1Hit_.size(); i++)
		{
			if (mm1Hit_[i].score == mm1Hit_[i - 1].score)
			{
				if (streak == 0)
				{
					streak = 1;
				}
				streak++;
			}
			else
			{
				if (streak > 1)
				{
					assert_geq(i, streak);
					mm1Hit_.shufflePortion(i - streak, streak, rnd);
				}
				streak = 0;
			}
		}
		if (streak > 1)
		{
			mm1Hit_.shufflePortion(mm1Hit_.size() - streak, streak, rnd);
		}
		mm1Sorted_ = true;
	}
	void add1mmEe(
		TIndexOffU top,
		TIndexOffU bot,
		const Edit *e1,
		const Edit *e2,
		bool fw,
		int64_t score)
	{
		mm1Hit_.expand();
		mm1Hit_.back().init(top, bot, e1, e2, fw, score);
		mm1Elt_ += (bot - top);
	}
	void addExactEeFw(
		TIndexOffU top,
		TIndexOffU bot,
		const Edit *e1,
		const Edit *e2,
		bool fw,
		int64_t score)
	{
		exactFwHit_.init(top, bot, e1, e2, fw, score);
	}
	void addExactEeRc(
		TIndexOffU top,
		TIndexOffU bot,
		const Edit *e1,
		const Edit *e2,
		bool fw,
		int64_t score)
	{
		exactRcHit_.init(top, bot, e1, e2, fw, score);
	}
	void clearExactE2eHits()
	{
		exactFwHit_.reset();
		exactRcHit_.reset();
	}
	void clear1mmE2eHits()
	{
		mm1Hit_.clear();
		mm1Elt_ = 0;
		mm1Sorted_ = false;
	}
	size_t numE2eHits() const
	{
		return exactFwHit_.size() + exactRcHit_.size() + mm1Elt_;
	}
	size_t numExactE2eHits() const
	{
		return exactFwHit_.size() + exactRcHit_.size();
	}
	size_t num1mmE2eHits() const
	{
		return mm1Elt_;
	}
	size_t readLength() const
	{
		assert(read_ != NULL);
		return read_->length();
	}
protected:
	EList<BTDnaString> seqFw_;
	EList<BTDnaString> seqRc_;
	EList<BTString> qualFw_;
	EList<BTString> qualRc_;
	EList<QVal> hitsFw_;
	EList<QVal> hitsRc_;
	EList<EList<InstantiatedSeed>> isFw_;
	EList<EList<InstantiatedSeed>> isRc_;
	EList<bool> sortedFw_;
	EList<bool> sortedRc_;
	size_t nonzTot_;
	size_t uniTot_;
	size_t uniTotS_[2];
	size_t repTot_;
	size_t repTotS_[2];
	size_t nonzFw_;
	size_t nonzRc_;
	size_t numRanges_;
	size_t numElts_;
	size_t numRangesFw_;
	size_t numEltsFw_;
	size_t numRangesRc_;
	size_t numEltsRc_;
	EList<uint32_t> offIdx2off_;
	EList<uint32_t> rankOffs_;
	EList<bool> rankFws_;
	bool sorted_;
	size_t numOffs_;
	const Read *read_;
	EEHit exactFwHit_;
	EEHit exactRcHit_;
	EList<EEHit> mm1Hit_;
	size_t mm1Elt_;
	bool mm1Sorted_;
	EList<size_t> tmpMedian_;
};
class Ebwt;
struct SideLocus;
struct SeedSearchMetrics
{
	SeedSearchMetrics() : mutex_m()
	{
		reset();
	}
	void merge(const SeedSearchMetrics &m, bool getLock = false)
	{
		seedsearch += m.seedsearch;
		nrange += m.nrange;
		nelt += m.nelt;
		possearch += m.possearch;
		intrahit += m.intrahit;
		interhit += m.interhit;
		filteredseed += m.filteredseed;
		ooms += m.ooms;
		bwops += m.bwops;
		bweds += m.bweds;
		bestmin0 += m.bestmin0;
		bestmin1 += m.bestmin1;
		bestmin2 += m.bestmin2;
	}
	void reset()
	{
		seedsearch =
			nrange =
				nelt =
					possearch =
						intrahit =
							interhit =
								filteredseed =
									ooms =
										bwops =
											bweds =
												bestmin0 =
													bestmin1 =
														bestmin2 = 0;
	}
	uint64_t seedsearch;
	uint64_t nrange;
	uint64_t nelt;
	uint64_t possearch;
	uint64_t intrahit;
	uint64_t interhit;
	uint64_t filteredseed;
	uint64_t ooms;
	uint64_t bwops;
	uint64_t bweds;
	uint64_t bestmin0;
	uint64_t bestmin1;
	uint64_t bestmin2;
	MUTEX_T mutex_m;
};
class SeedAligner
{
public:
	SeedAligner() : edits_(AL_CAT), offIdx2off_(AL_CAT) {}
	void instantiateSeq(
		const Read &read,
		BTDnaString &seq,
		BTString &qual,
		int len,
		int depth,
		bool fw) const;
	std::pair<int, int> instantiateSeeds(
		const EList<Seed> &seeds,
		size_t off,
		int per,
		const Read &read,
		const Scoring &pens,
		bool nofw,
		bool norc,
		AlignmentCacheIface &cache,
		SeedResults &sr,
		SeedSearchMetrics &met,
		std::pair<int, int> &instFw,
		std::pair<int, int> &instRc);
	void searchAllSeeds(
		const EList<Seed> &seeds,
		const Ebwt *ebwtFw,
		const Ebwt *ebwtBw,
		const Read &read,
		const Scoring &pens,
		AlignmentCacheIface &cache,
		SeedResults &hits,
		SeedSearchMetrics &met,
		PerReadMetrics &prm);
	bool sanityPartial(
		const Ebwt *ebwtFw,
		const Ebwt *ebwtBw,
		const BTDnaString &seq,
		size_t dep,
		size_t len,
		bool do1mm,
		TIndexOffU topfw,
		TIndexOffU botfw,
		TIndexOffU topbw,
		TIndexOffU botbw);
	size_t exactSweep(
		const Ebwt &ebwt,
		const Read &read,
		const Scoring &sc,
		bool nofw,
		bool norc,
		size_t mineMax,
		size_t &mineFw,
		size_t &mineRc,
		bool repex,
		SeedResults &hits,
		SeedSearchMetrics &met);
	bool oneMmSearch(
		const Ebwt *ebwtFw,
		const Ebwt *ebwtBw,
		const Read &read,
		const Scoring &sc,
		int64_t minsc,
		bool nofw,
		bool norc,
		bool local,
		bool repex,
		bool rep1mm,
		SeedResults &hits,
		SeedSearchMetrics &met);
protected:
	bool extendAndReportHit(
		TIndexOffU topf,
		TIndexOffU botf,
		TIndexOffU topb,
		TIndexOffU botb,
		uint16_t len,
		DoublyLinkedList<Edit> *prevEdit);
	bool reportHit(
		TIndexOffU topf,
		TIndexOffU botf,
		TIndexOffU topb,
		TIndexOffU botb,
		uint16_t len,
		DoublyLinkedList<Edit> *prevEdit);
	bool searchSeedBi();
	bool searchSeedBi(
		int step,
		int depth,
		TIndexOffU topf,
		TIndexOffU botf,
		TIndexOffU topb,
		TIndexOffU botb,
		SideLocus tloc,
		SideLocus bloc,
		Constraint c0,
		Constraint c1,
		Constraint c2,
		Constraint overall,
		DoublyLinkedList<Edit> *prevEdit);
	inline void nextLocsBi(
		SideLocus &tloc,
		SideLocus &bloc,
		TIndexOffU topf,
		TIndexOffU botf,
		TIndexOffU topb,
		TIndexOffU botb,
		int step);
	const Ebwt *ebwtFw_;
	const Ebwt *ebwtBw_;
	const Scoring *sc_;
	const InstantiatedSeed *s_;
	const Read *read_;
	const BTDnaString *seq_;
	const BTString *qual_;
	size_t off_;
	bool fw_;
	EList<Edit> edits_;
	AlignmentCacheIface *ca_;
	EList<uint32_t> offIdx2off_;
	uint64_t bwops_;
	uint64_t bwedits_;
	BTDnaString tmprfdnastr_;
	ASSERT_ONLY(ESet<BTDnaString> hits_);
	BTDnaString tmpdnastr_;
};
#define INIT_LOCS(top, bot, tloc, bloc, e)                                         \
	{                                                                              \
		if (bot - top == 1)                                                        \
		{                                                                          \
			tloc.initFromRow(top, (e).eh(), (e).ebwt());                           \
			bloc.invalidate();                                                     \
		}                                                                          \
		else                                                                       \
		{                                                                          \
			SideLocus::initFromTopBot(top, bot, (e).eh(), (e).ebwt(), tloc, bloc); \
			assert(bloc.valid());                                                  \
		}                                                                          \
	}
#define SANITY_CHECK_4TUP(t, b, tp, bp)                                                                       \
	{                                                                                                         \
		ASSERT_ONLY(TIndexOffU tot = (b[0] - t[0]) + (b[1] - t[1]) + (b[2] - t[2]) + (b[3] - t[3]));          \
		ASSERT_ONLY(TIndexOffU totp = (bp[0] - tp[0]) + (bp[1] - tp[1]) + (bp[2] - tp[2]) + (bp[3] - tp[3])); \
		assert_eq(tot, totp);                                                                                 \
	}
#endif
