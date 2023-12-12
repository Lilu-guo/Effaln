#ifndef ALIGNER_DRIVER_H_
#define ALIGNER_DRIVER_H_
#include "aligner_seed2.h"
#include "simple_func.h"
#include "aln_sink.h"
class IntervalRootSelector : public DescentRootSelector
{
public:
	IntervalRootSelector(
		double consExp,
		const SimpleFunc &rootIval,
		size_t landing)
	{
		consExp_ = consExp;
		rootIval_ = rootIval;
		landing_ = landing;
	}
	virtual ~IntervalRootSelector() {}
	virtual void select(
		const Read &q,
		const Read *qo,
		bool nofw,
		bool norc,
		EList<DescentConfig> &confs,
		EList<DescentRoot> &roots);
protected:
	double consExp_;
	SimpleFunc rootIval_;
	size_t landing_;
};
class PrioritizedRootSelector : public DescentRootSelector
{
public:
	PrioritizedRootSelector(
		double consExp,
		const SimpleFunc &rootIval,
		size_t landing)
	{
		consExp_ = consExp;
		rootIval_ = rootIval;
		landing_ = landing;
	}
	virtual ~PrioritizedRootSelector() {}
	virtual void select(
		const Read &q,
		const Read *qo,
		bool nofw,
		bool norc,
		EList<DescentConfig> &confs,
		EList<DescentRoot> &roots);
protected:
	double consExp_;
	SimpleFunc rootIval_;
	size_t landing_;
	EHeap<DescentRoot> rootHeap_;
	EList<int> scoresOrig_[2];
	EList<int> scores_[2];
};
enum
{
	ALDRIVER_EXHAUSTED_CANDIDATES = 1,
	ALDRIVER_POLICY_FULFILLED,
	ALDRIVER_EXCEEDED_LIMIT
};
class AlignerDriver
{
public:
	AlignerDriver(
		double consExp,
		bool prioritizeRoots,
		const SimpleFunc &rootIval,
		size_t landing,
		bool veryVerbose,
		const SimpleFunc &totsz,
		const SimpleFunc &totfmops) : alsel_(),
									  dr1_(veryVerbose),
									  dr2_(veryVerbose)
	{
		assert_gt(landing, 0);
		totsz_ = totsz;
		totfmops_ = totfmops;
		if (prioritizeRoots)
		{
			sel_ = new PrioritizedRootSelector(consExp, rootIval, landing);
		}
		else
		{
			sel_ = new IntervalRootSelector(consExp, rootIval, landing);
		}
	}
	virtual ~AlignerDriver()
	{
		delete sel_;
	}
	void initRead(
		const Read &q1,
		bool nofw,
		bool norc,
		TAlScore minsc,
		TAlScore maxpen,
		const Read *q2)
	{
		dr1_.initRead(q1, nofw, norc, minsc, maxpen, q2, sel_);
		red1_.init(q1.length());
		paired_ = false;
		if (q2 != NULL)
		{
			dr2_.initRead(*q2, nofw, norc, minsc, maxpen, &q1, sel_);
			red2_.init(q2->length());
			paired_ = true;
		}
		else
		{
			dr2_.reset();
		}
		size_t totsz = totsz_.f<size_t>(q1.length());
		size_t totfmops = totfmops_.f<size_t>(q1.length());
		stop_.init(
			totsz,
			0,
			true,
			totfmops);
	}
	int go(
		const Scoring &sc,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		const BitPairReference &ref,
		DescentMetrics &met,
		WalkMetrics &wlm,
		PerReadMetrics &prm,
		RandomSource &rnd,
		AlnSinkWrap &sink);
	void reset()
	{
		dr1_.reset();
		dr2_.reset();
		red1_.reset();
		red2_.reset();
	}
	const DescentDriver &dr1() { return dr1_; }
	const DescentDriver &dr2() { return dr2_; }
protected:
	DescentRootSelector *sel_;
	DescentAlignmentSelector alsel_;
	DescentDriver dr1_;
	DescentDriver dr2_;
	DescentStoppingConditions stop_;
	bool paired_;
	SimpleFunc totsz_;
	SimpleFunc totfmops_;
	RedundantAlns red1_;
	RedundantAlns red2_;
	ASSERT_ONLY(SStringExpandable<char> raw_refbuf_);
	ASSERT_ONLY(SStringExpandable<uint32_t> raw_destU32_);
	ASSERT_ONLY(EList<bool> raw_matches_);
	ASSERT_ONLY(BTDnaString tmp_rf_);
	ASSERT_ONLY(BTDnaString tmp_rdseq_);
	ASSERT_ONLY(BTString tmp_qseq_);
};
#endif
