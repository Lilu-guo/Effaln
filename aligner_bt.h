#ifndef ALIGNER_BT_H_
#define ALIGNER_BT_H_
#include <utility>
#include <stdint.h>
#include "aligner_sw_common.h"
#include "aligner_result.h"
#include "scoring.h"
#include "edit.h"
#include "limit.h"
#include "dp_framer.h"
#include "sse_util.h"
enum
{
	BT_NOT_FOUND = 1,
	BT_FOUND,
	BT_REJECTED_N,
	BT_REJECTED_CORE_DIAG
};
class BtBranchProblem
{
public:
	BtBranchProblem() { reset(); }
	void initRef(
		const char *qry,
		const char *qual,
		size_t qrylen,
		const char *ref,
		TRefOff reflen,
		TRefOff treflen,
		TRefId refid,
		TRefOff refoff,
		bool fw,
		const DPRect *rect,
		const Checkpointer *cper,
		const Scoring *sc,
		size_t nceil)
	{
		qry_ = qry;
		qual_ = qual;
		qrylen_ = qrylen;
		ref_ = ref;
		reflen_ = reflen;
		treflen_ = treflen;
		refid_ = refid;
		refoff_ = refoff;
		fw_ = fw;
		rect_ = rect;
		cper_ = cper;
		sc_ = sc;
		nceil_ = nceil;
	}
	void initBt(
		size_t row,
		size_t col,
		bool fill,
		bool usecp,
		TAlScore targ)
	{
		row_ = row;
		col_ = col;
		targ_ = targ;
		fill_ = fill;
		usecp_ = usecp;
		if (fill)
		{
			assert(usecp_);
		}
	}
	void reset()
	{
		qry_ = qual_ = ref_ = NULL;
		cper_ = NULL;
		rect_ = NULL;
		sc_ = NULL;
		qrylen_ = reflen_ = treflen_ = refid_ = refoff_ = row_ = col_ = targ_ = nceil_ = 0;
		fill_ = fw_ = usecp_ = false;
	}
	bool inited() const
	{
		return qry_ != NULL;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_gt(qrylen_, 0);
		assert_gt(reflen_, 0);
		assert_gt(treflen_, 0);
		assert_lt(row_, qrylen_);
		assert_lt((TRefOff)col_, reflen_);
		return true;
	}
#endif
	size_t reflen() const
	{
		return reflen_;
	}
	size_t treflen() const { return treflen_; }
protected:
	const char *qry_;
	const char *qual_;
	size_t qrylen_;
	const char *ref_;
	TRefOff reflen_;
	TRefOff treflen_;
	TRefId refid_;
	TRefOff refoff_;
	bool fw_;
	const DPRect *rect_;
	size_t row_;
	size_t col_;
	TAlScore targ_;
	const Checkpointer *cper_;
	bool fill_;
	bool usecp_;
	const Scoring *sc_;
	size_t nceil_;
	friend class BtBranch;
	friend class BtBranchQ;
	friend class BtBranchTracer;
};
class BtBranch
{
public:
	BtBranch() { reset(); }
	BtBranch(
		const BtBranchProblem &prob,
		size_t parentId,
		TAlScore penalty,
		TAlScore score_en,
		int64_t row,
		int64_t col,
		Edit e,
		int hef,
		bool root,
		bool extend)
	{
		init(prob, parentId, penalty, score_en, row, col, e, hef, root, extend);
	}
	void reset()
	{
		parentId_ = 0;
		score_st_ = score_en_ = len_ = row_ = col_ = 0;
		curtailed_ = false;
		e_.reset();
	}
	void init(
		const BtBranchProblem &prob,
		size_t parentId,
		TAlScore penalty,
		TAlScore score_en,
		int64_t row,
		int64_t col,
		Edit e,
		int hef,
		bool root,
		bool extend);
	bool isSolution(const BtBranchProblem &prob) const
	{
		const bool end2end = prob.sc_->monotone;
		return score_st_ == prob.targ_ && (!end2end || endsInFirstRow());
	}
	bool isValid(const BtBranchProblem &prob) const
	{
		int64_t scoreFloor = prob.sc_->monotone ? MIN_I64 : 0;
		if (score_st_ < scoreFloor)
		{
			return false;
		}
		if (isSolution(prob))
		{
			return true;
		}
		if ((int64_t)len_ > row_)
		{
			return score_st_ == prob.targ_;
		}
		else
		{
			int64_t match = prob.sc_->match();
			int64_t bonusLeft = (row_ + 1 - len_) * match;
			return score_st_ + bonusLeft >= prob.targ_;
		}
	}
	bool overlap(const BtBranchProblem &prob, const BtBranch &bt) const
	{
		assert_lt(row_, (int64_t)prob.qrylen_);
		size_t fromend = prob.qrylen_ - row_ - 1;
		size_t diag = fromend + col_;
		int64_t lo = 0, hi = row_ + 1;
		if (len_ == 0)
		{
			lo = row_;
		}
		else
		{
			lo = row_ - (len_ - 1);
		}
		assert_lt(bt.row_, (int64_t)prob.qrylen_);
		size_t ofromend = prob.qrylen_ - bt.row_ - 1;
		size_t odiag = ofromend + bt.col_;
		if (diag != odiag)
		{
			return false;
		}
		int64_t olo = 0, ohi = bt.row_ + 1;
		if (bt.len_ == 0)
		{
			olo = bt.row_;
		}
		else
		{
			olo = bt.row_ - (bt.len_ - 1);
		}
		int64_t losm = olo, hism = ohi;
		if (hi - lo < ohi - olo)
		{
			swap(lo, losm);
			swap(hi, hism);
		}
		if ((lo <= losm && hi > losm) || (lo < hism && hi >= hism))
		{
			return true;
		}
		return false;
	}
	bool operator<(const BtBranch &o) const
	{
		if (uppermostRow() != o.uppermostRow())
		{
			return uppermostRow() < o.uppermostRow();
		}
		if (score_st_ != o.score_st_)
			return score_st_ > o.score_st_;
		if (row_ != o.row_)
			return row_ < o.row_;
		if (col_ != o.col_)
			return col_ > o.col_;
		if (parentId_ != o.parentId_)
			return parentId_ > o.parentId_;
		assert(false);
		return false;
	}
	bool endsInFirstRow() const
	{
		assert_leq((int64_t)len_, row_ + 1);
		return (int64_t)len_ == row_ + 1;
	}
	size_t uppermostRow() const
	{
		assert_geq(row_ + 1, (int64_t)len_);
		return row_ + 1 - (int64_t)len_;
	}
	size_t leftmostCol() const
	{
		assert_geq(col_ + 1, (int64_t)len_);
		return col_ + 1 - (int64_t)len_;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert(root_ || e_.inited());
		assert_gt(len_, 0);
		assert_geq(col_ + 1, (int64_t)len_);
		assert_geq(row_ + 1, (int64_t)len_);
		return true;
	}
#endif
protected:
	size_t parentId_;
	TAlScore penalty_;
	TAlScore score_st_;
	TAlScore score_en_;
	size_t len_;
	int64_t row_;
	int64_t col_;
	Edit e_;
	bool root_;
	bool curtailed_;
	friend class BtBranchQ;
	friend class BtBranchTracer;
};
class BtBranchTracer
{
public:
	explicit BtBranchTracer() : prob_(), bs_(), seenPaths_(DP_CAT), sawcell_(DP_CAT), doTri_() {}
	void add(size_t id)
	{
		assert(!bs_[id].isSolution(prob_));
		unsorted_.push_back(make_pair(bs_[id].score_st_, id));
	}
	void addSolution(size_t id)
	{
		assert(bs_[id].isSolution(prob_));
		solutions_.push_back(id);
	}
	void examineBranch(
		int64_t row,
		int64_t col,
		const Edit &e,
		TAlScore pen,
		TAlScore sc,
		size_t parentId);
	void addOffshoots(size_t bid);
	size_t best(RandomSource &rnd)
	{
		assert(!empty());
		flushUnsorted();
		assert_gt(sortedSel_ ? sorted1_.size() : sorted2_.size(), cur_);
		size_t id = sortedSel_ ? sorted1_[cur_] : sorted2_[cur_];
		cur_++;
		return id;
	}
	bool empty() const
	{
		return size() == 0;
	}
	size_t size() const
	{
		return unsorted_.size() +
			   (sortedSel_ ? sorted1_.size() : sorted2_.size()) - cur_;
	}
	bool emptySolution() const
	{
		return sizeSolution() == 0;
	}
	size_t sizeSolution() const
	{
		return solutions_.size();
	}
	void flushUnsorted();
#ifndef NDEBUG
	bool repOk() const
	{
		assert_lt(cur_, (sortedSel_ ? sorted1_.size() : sorted2_.size()));
		return true;
	}
#endif
	void initRef(
		const char *rd,
		const char *qu,
		size_t rdlen,
		const char *rf,
		size_t rflen,
		TRefOff trflen,
		TRefId refid,
		TRefOff refoff,
		bool fw,
		const DPRect *rect,
		const Checkpointer *cper,
		const Scoring &sc,
		size_t nceil)
	{
		prob_.initRef(rd, qu, rdlen, rf, rflen, trflen, refid, refoff, fw, rect, cper, &sc, nceil);
		const size_t ndiag = rflen + rdlen - 1;
		seenPaths_.resize(ndiag);
		for (size_t i = 0; i < ndiag; i++)
		{
			seenPaths_[i].clear();
		}
		if (sawcell_.size() < rflen)
		{
			size_t isz = sawcell_.size();
			sawcell_.resize(rflen);
			for (size_t i = isz; i < rflen; i++)
			{
				sawcell_[i].setCat(DP_CAT);
			}
		}
		for (size_t i = 0; i < rflen; i++)
		{
			sawcell_[i].setCat(DP_CAT);
			sawcell_[i].clear();
		}
	}
	void initBt(
		TAlScore escore,
		size_t row,
		size_t col,
		bool fill,
		bool usecp,
		bool doTri,
		RandomSource &rnd)
	{
		prob_.initBt(row, col, fill, usecp, escore);
		Edit e;
		e.reset();
		unsorted_.clear();
		solutions_.clear();
		sorted1_.clear();
		sorted2_.clear();
		cur_ = 0;
		nmm_ = 0;
		nnmm_ = 0;
		nrdop_ = 0;
		nrfop_ = 0;
		nrdex_ = 0;
		nrfex_ = 0;
		nmmPrune_ = 0;
		nnmmPrune_ = 0;
		nrdopPrune_ = 0;
		nrfopPrune_ = 0;
		nrdexPrune_ = 0;
		nrfexPrune_ = 0;
		row_ = row;
		col_ = col;
		doTri_ = doTri;
		bs_.clear();
		if (!prob_.fill_)
		{
			size_t id = bs_.alloc();
			bs_[id].init(
				prob_,
				0,
				0,
				0,
				row,
				col,
				e,
				0,
				true,
				true);
			if (bs_[id].isSolution(prob_))
			{
				addSolution(id);
			}
			else
			{
				add(id);
			}
		}
		else
		{
			int64_t row = row_, col = col_;
			TAlScore targsc = prob_.targ_;
			int hef = 0;
			bool done = false, abort = false;
			size_t depth = 0;
			while (!done && !abort)
			{
				if (doTri_)
				{
					triangleFill(
						row,
						col,
						hef,
						targsc,
						prob_.targ_,
						rnd,
						row,
						col,
						hef,
						targsc,
						done,
						abort);
				}
				else
				{
					squareFill(
						row,
						col,
						hef,
						targsc,
						prob_.targ_,
						rnd,
						row,
						col,
						hef,
						targsc,
						done,
						abort);
				}
				if (depth >= ndep_.size())
				{
					ndep_.resize(depth + 1);
					ndep_[depth] = 1;
				}
				else
				{
					ndep_[depth]++;
				}
				depth++;
				assert((row >= 0 && col >= 0) || done);
			}
		}
		ASSERT_ONLY(seen_.clear());
	}
	bool nextAlignment(
		size_t maxiter,
		SwResult &res,
		size_t &off,
		size_t &nrej,
		size_t &niter,
		RandomSource &rnd);
	bool inited() const
	{
		return prob_.inited();
	}
	bool doTri() const { return doTri_; }
	void triangleFill(
		int64_t rw,
		int64_t cl,
		int hef,
		TAlScore targ,
		TAlScore targ_final,
		RandomSource &rnd,
		int64_t &row_new,
		int64_t &col_new,
		int &hef_new,
		TAlScore &targ_new,
		bool &done,
		bool &abort);
	void squareFill(
		int64_t rw,
		int64_t cl,
		int hef,
		TAlScore targ,
		TAlScore targ_final,
		RandomSource &rnd,
		int64_t &row_new,
		int64_t &col_new,
		int &hef_new,
		TAlScore &targ_new,
		bool &done,
		bool &abort);
protected:
	bool nextAlignmentBacktrace(
		size_t maxiter,
		SwResult &res,
		size_t &off,
		size_t &nrej,
		size_t &niter,
		RandomSource &rnd);
	bool nextAlignmentFill(
		size_t maxiter,
		SwResult &res,
		size_t &off,
		size_t &nrej,
		size_t &niter,
		RandomSource &rnd);
	bool trySolutions(
		bool lookForOlap,
		SwResult &res,
		size_t &off,
		size_t &nrej,
		RandomSource &rnd,
		bool &success);
	int trySolution(
		size_t id,
		bool lookForOlap,
		SwResult &res,
		size_t &off,
		size_t &nrej,
		RandomSource &rnd);
	BtBranchProblem prob_;
	EFactory<BtBranch> bs_;
	ELList<std::pair<size_t, size_t>> seenPaths_;
	ELSet<size_t> sawcell_;
	EList<std::pair<TAlScore, size_t>> unsorted_;
	EList<size_t> sorted1_;
	EList<size_t> sorted2_;
	EList<size_t> solutions_;
	bool sortedSel_;
	size_t cur_;
	size_t nmm_;
	size_t nnmm_;
	size_t nrdop_;
	size_t nrfop_;
	size_t nrdex_;
	size_t nrfex_;
	size_t nmmPrune_;
	size_t nnmmPrune_;
	size_t nrdopPrune_;
	size_t nrfopPrune_;
	size_t nrdexPrune_;
	size_t nrfexPrune_;
	size_t row_;
	size_t col_;
	bool doTri_;
	EList<CpQuad> sq_;
	ELList<CpQuad> tri_;
	EList<size_t> ndep_;
#ifndef NDEBUG
	ESet<size_t> seen_;
#endif
};
#endif
