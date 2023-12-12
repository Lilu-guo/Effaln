#include <limits>
#include <ctype.h>
#include "aligner_seed2.h"
#include "assert_helpers.h"
#include "aln_idx.h"
void DescentDriver::go(
	const Scoring &sc,
	const Ebwt &ebwtFw, const Ebwt &ebwtBw, DescentMetrics &met, PerReadMetrics &prm)
{
	assert(q_.repOk());
	for (size_t i = 0; i < roots_.size(); i++)
	{
		size_t dfsz = df_.size();
		size_t pfsz = pf_.size();
		TDescentId id = df_.alloc();
		Edit e_null;
		assert(!e_null.inited());
		bool succ = df_[id].init(
			q_,
			i, sc, minsc_, maxpen_, id, ebwtFw, ebwtBw, re_, df_, pf_, roots_, confs_, heap_, alsink_, met, prm);
		if (veryVerbose_)
		{
			bool fw = roots_[i].fw;
			tmpedit_.clear();
			df_[id].print(
				&cerr,
				"",
				q_,
				0,
				0,
				fw,
				tmpedit_,
				0,
				tmpedit_.size(),
				tmprfdnastr_);
		}
		if (!succ)
		{
			df_.resize(dfsz);
			pf_.resize(pfsz);
		}
	}
	bool stop = heap_.empty();
	while (!stop)
	{
		TDescentPair p = heap_.pop();
		df_.alloc();
		df_.pop();
		df_[p.second].followBestOutgoing(
			q_,
			ebwtFw, ebwtBw, sc, minsc_, maxpen_, re_, df_, pf_, roots_, confs_, heap_, alsink_, met, prm);
		stop = heap_.empty();
	}
}
int DescentDriver::advance(
	const DescentStoppingConditions &stopc,
	const Scoring &sc, const Ebwt &ebwtFw, const Ebwt &ebwtBw, DescentMetrics &met, PerReadMetrics &prm)
{
	size_t nbwop_i = met.bwops;
	while (rootsInited_ < roots_.size())
	{
		size_t dfsz = df_.size();
		size_t pfsz = pf_.size();
		TDescentId id = df_.alloc();
		Edit e_null;
		assert(!e_null.inited());
		bool succ = df_[id].init(
			q_,
			rootsInited_, sc, minsc_, maxpen_, id, ebwtFw, ebwtBw, re_, df_, pf_, roots_, confs_, heap_, alsink_, met, prm);
		if (!succ)
		{
			df_.resize(dfsz);
			pf_.resize(pfsz);
		}
		rootsInited_++;
		TAlScore best = std::numeric_limits<TAlScore>::max();
		if (!heap_.empty())
		{
			best = heap_.top().first.pen;
		}
		if (stopc.nfound > 0 && alsink_.nelt() > stopc.nfound)
		{
			return DESCENT_DRIVER_ALN;
		}
		if (alsink_.stratumDone(best))
		{
			return DESCENT_DRIVER_STRATA;
		}
		if (stopc.nbwop > 0 && (met.bwops - nbwop_i) > stopc.nbwop)
		{
			return DESCENT_DRIVER_BWOPS;
		}
		if (stopc.totsz > 0 && totalSizeBytes() > stopc.totsz)
		{
			return DESCENT_DRIVER_MEM;
		}
	}
	bool stop = heap_.empty();
	while (!stop)
	{
		TDescentPair p = heap_.pop();
		df_.alloc();
		df_.pop();
		df_[p.second].followBestOutgoing(
			q_,
			ebwtFw, ebwtBw, sc, minsc_, maxpen_, re_, df_, pf_, roots_, confs_, heap_, alsink_, met, prm);
		TAlScore best = std::numeric_limits<TAlScore>::max();
		if (!heap_.empty())
		{
			best = heap_.top().first.pen;
		}
		if (stopc.nfound > 0 && alsink_.nelt() > stopc.nfound)
		{
			return DESCENT_DRIVER_ALN;
		}
		if (alsink_.stratumDone(best))
		{
			return DESCENT_DRIVER_STRATA;
		}
		if (stopc.nbwop > 0 && (met.bwops - nbwop_i) > stopc.nbwop)
		{
			return DESCENT_DRIVER_BWOPS;
		}
		if (stopc.totsz > 0 && totalSizeBytes() > stopc.totsz)
		{
			return DESCENT_DRIVER_MEM;
		}
		stop = heap_.empty();
	}
	return DESCENT_DRIVER_DONE;
}
bool DescentAlignmentSink::reportAlignment(
	const Read &q,
	const Ebwt &ebwtFw, const Ebwt &ebwtBw, TIndexOffU topf, TIndexOffU botf, TIndexOffU topb, TIndexOffU botb, TDescentId id, TRootId rid, const Edit &e, TScore pen, EFactory<Descent> &df, EFactory<DescentPos> &pf, const EList<DescentRoot> &rs, const EList<DescentConfig> &cs)
{
	TDescentId cur = id;
	ASSERT_ONLY(const Descent &desc = df[id]);
	const bool fw = rs[rid].fw;
	ASSERT_ONLY(size_t len = q.length());
	assert(q.repOk());
	assert_lt(desc.al5pf(), len);
	Triple<TIndexOffU, TIndexOffU, size_t> lhs(topf, botf, 0);
	Triple<TIndexOffU, TIndexOffU, size_t> rhs(topb, botb, q.length() - 1);
	if (!lhs_.insert(lhs))
	{
		rhs_.insert(rhs);
		return false;
	}
	if (!rhs_.insert(rhs))
	{
		return false;
	}
	size_t ei = edits_.size();
	df[cur].collectEdits(edits_, &e, df);
	size_t en = edits_.size() - ei;
#ifndef NDEBUG
	{
		for (size_t i = 1; i < en; i++)
		{
			assert_geq(edits_[ei + i].pos, edits_[ei + i - 1].pos);
		}
		size_t trimLf = 0;
		size_t trimRg = 0;
		BTDnaString &rf = tmprfdnastr_;
		rf.clear();
		if (!fw)
		{
			Edit::invertPoss(edits_, len, ei, en, true);
		}
		desc.print(NULL, "", q, trimLf, trimRg, fw, edits_, ei, en, rf);
		if (!fw)
		{
			Edit::invertPoss(edits_, len, ei, en, true);
		}
		ASSERT_ONLY(TIndexOffU toptmp = 0);
		ASSERT_ONLY(TIndexOffU bottmp = 0);
		if (!ebwtFw.contains(rf, &toptmp, &bottmp))
		{
			std::cerr << rf << std::endl;
			assert(false);
		}
	}
#endif
	als_.expand();
	als_.back().init(pen, fw, topf, botf, ei, en);
	nelt_ += (botf - topf);
	if (bestPen_ == std::numeric_limits<TAlScore>::max() || pen < bestPen_)
	{
		bestPen_ = pen;
	}
	if (worstPen_ == std::numeric_limits<TAlScore>::max() || pen > worstPen_)
	{
		worstPen_ = pen;
	}
	return true;
}
bool Descent::init(
	const Read &q,
	TRootId rid, const Scoring &sc, TAlScore minsc, TAlScore maxpen, TReadOff al5pi, TReadOff al5pf, TIndexOffU topf, TIndexOffU botf, TIndexOffU topb, TIndexOffU botb, bool l2r, size_t descid, TDescentId parent, TScore pen, const Edit &e, const Ebwt &ebwtFw, const Ebwt &ebwtBw, DescentRedundancyChecker &re, EFactory<Descent> &df, EFactory<DescentPos> &pf, const EList<DescentRoot> &rs, const EList<DescentConfig> &cs, EHeap<TDescentPair> &heap, DescentAlignmentSink &alsink, DescentMetrics &met, PerReadMetrics &prm)
{
	assert(q.repOk());
	rid_ = rid;
	al5pi_ = al5pi;
	al5pf_ = al5pf;
	l2r_ = l2r;
	topf_ = topf;
	botf_ = botf;
	topb_ = topb;
	botb_ = botb;
	descid_ = descid;
	parent_ = parent;
	pen_ = pen;
	posid_ = std::numeric_limits<size_t>::max();
	len_ = 0;
	out_.clear();
	edit_ = e;
	lastRecalc_ = true;
	gapadd_ = df[parent].gapadd_;
	if (e.inited())
	{
		if (e.isReadGap())
		{
			gapadd_++;
		}
		else if (e.isRefGap())
		{
			gapadd_--;
		}
	}
	bool branches = false, hitEnd = false, done = false;
	TIndexOffU topf_new = 0, botf_new = 0, topb_new = 0, botb_new = 0;
	off5p_i_ = 0;
#ifndef NDEBUG
	size_t depth = al5pf_ - al5pi_ + 1;
	TAlScore maxpend = cs[rid_].cons.get(depth, q.length(), maxpen);
	assert_geq(maxpend, pen_);
#endif
	bool matchSucc = followMatches(
		q,
		sc,
		ebwtFw,
		ebwtBw,
		re,
		df,
		pf,
		rs,
		cs,
		heap,
		alsink,
		met,
		prm,
		branches,
		hitEnd,
		done,
		off5p_i_,
		topf_new,
		botf_new,
		topb_new,
		botb_new);
	bool bounceSucc = false;
	if (matchSucc && hitEnd && !done)
	{
		assert(topf_new > 0 || botf_new > 0);
		bounceSucc = bounce(
			q,
			topf_new,
			botf_new,
			topb_new,
			botb_new,
			ebwtFw,
			ebwtBw,
			sc,
			minsc,
			maxpen, re,
			df,
			pf,
			rs,
			cs,
			heap,
			alsink,
			met,
			prm);
	}
	if (matchSucc)
	{
		recalcOutgoing(q, sc, minsc, maxpen, re, pf, rs, cs, prm);
		if (!empty())
		{
			heap.insert(make_pair(out_.bestPri(), descid));
		}
	}
	return !empty() || bounceSucc;
}
bool Descent::init(
	const Read &q,
	TRootId rid, const Scoring &sc, TAlScore minsc, TAlScore maxpen, size_t descid, const Ebwt &ebwtFw, const Ebwt &ebwtBw, DescentRedundancyChecker &re, EFactory<Descent> &df, EFactory<DescentPos> &pf, const EList<DescentRoot> &rs, const EList<DescentConfig> &cs, EHeap<TDescentPair> &heap, DescentAlignmentSink &alsink, DescentMetrics &met, PerReadMetrics &prm)
{
	rid_ = rid;
	al5pi_ = rs[rid].off5p;
	al5pf_ = rs[rid].off5p;
	assert_lt(al5pi_, q.length());
	assert_lt(al5pf_, q.length());
	l2r_ = rs[rid].l2r;
	topf_ = botf_ = topb_ = botb_ = 0;
	descid_ = descid;
	parent_ = std::numeric_limits<size_t>::max();
	pen_ = 0;
	posid_ = std::numeric_limits<size_t>::max();
	len_ = 0;
	out_.clear();
	edit_.reset();
	lastRecalc_ = true;
	gapadd_ = 0;
	bool branches = false, hitEnd = false, done = false;
	TIndexOffU topf_new = 0, botf_new = 0, topb_new = 0, botb_new = 0;
	off5p_i_ = 0;
	bool matchSucc = followMatches(
		q,
		sc,
		ebwtFw,
		ebwtBw,
		re,
		df,
		pf,
		rs,
		cs,
		heap,
		alsink,
		met,
		prm,
		branches,
		hitEnd,
		done,
		off5p_i_,
		topf_new,
		botf_new,
		topb_new,
		botb_new);
	bool bounceSucc = false;
	if (matchSucc && hitEnd && !done)
	{
		assert(topf_new > 0 || botf_new > 0);
		bounceSucc = bounce(
			q,
			topf_new,
			botf_new,
			topb_new,
			botb_new,
			ebwtFw,
			ebwtBw,
			sc,
			minsc,
			maxpen, re,
			df,
			pf,
			rs,
			cs,
			heap,
			alsink,
			met,
			prm);
	}
	assert(empty());
	if (matchSucc)
	{
		recalcOutgoing(q, sc, minsc, maxpen, re, pf, rs, cs, prm);
		if (!empty())
		{
			heap.insert(make_pair(out_.bestPri(), descid));
		}
	}
	return !empty() || bounceSucc;
}
size_t Descent::recalcOutgoing(
	const Read &q,
	const Scoring &sc, TAlScore minsc, TAlScore maxpen, DescentRedundancyChecker &re, EFactory<DescentPos> &pf, const EList<DescentRoot> &rs, const EList<DescentConfig> &cs, PerReadMetrics &prm)
{
	assert_eq(botf_ - topf_, botb_ - topb_);
	assert(out_.empty());
	assert(repOk(&q));
	bool fw = rs[rid_].fw;
	float rootpri = rs[rid_].pri;
	bool toward3p = (l2r_ == fw);
	size_t off5p = off5p_i_;
	assert_geq(al5pf_, al5pi_);
	size_t off3p = q.length() - off5p - 1;
	size_t depth, extrai = 0, extraf = 0;
	size_t cur5pi = al5pi_, cur5pf = al5pf_;
	if (toward3p)
	{
		cur5pf = off5p;
		depth = off5p - al5pi_;
		if (al5pf_ < q.length() - 1)
		{
			extraf = 1;
		}
	}
	else
	{
		cur5pi = off5p;
		depth = al5pf_ - off5p;
		if (al5pi_ > 0)
		{
			extrai = 1;
		}
	}
	TScore pen_rdg_ex = sc.readGapExtend(), pen_rfg_ex = sc.refGapExtend();
	TScore pen_rdg_op = sc.readGapOpen(), pen_rfg_op = sc.refGapOpen();
	TIndexOffU top = l2r_ ? topb_ : topf_;
	TIndexOffU bot = l2r_ ? botb_ : botf_;
	TIndexOffU topp = l2r_ ? topf_ : topb_;
	TIndexOffU botp = l2r_ ? botf_ : botb_;
	assert_eq(botp - topp, bot - top);
	DescentEdge edge;
	size_t nout = 0;
	size_t d = posid_;
	while (off5p >= al5pi_ - extrai && off5p <= al5pf_ + extraf)
	{
		assert_lt(off5p, q.length());
		assert_lt(off3p, q.length());
		TScore maxpend = cs[rid_].cons.get(depth, q.length(), maxpen);
		assert(depth > 0 || maxpend == 0);
		assert_geq(maxpend, pen_);
		TScore diff = maxpend - pen_;
		const TIndexOffU *t = l2r_ ? pf[d].topb : pf[d].topf;
		const TIndexOffU *b = l2r_ ? pf[d].botb : pf[d].botf;
		const TIndexOffU *tp = l2r_ ? pf[d].topf : pf[d].topb;
		const TIndexOffU *bp = l2r_ ? pf[d].botf : pf[d].botb;
		assert_eq(pf[d].botf - pf[d].topf, pf[d].botb - pf[d].topb);
		std::pair<int, int> p = q.get(off5p, fw);
		int c = p.first;
		assert_range(0, 4, c);
		if (!pf[d].flags.exhausted() && diff > 0)
		{
			int qq = p.second;
			assert_geq(qq, 0);
			TScore pen_mm = sc.mm(c, qq);
			if (pen_mm <= diff)
			{
				for (int j = 0; j < 4; j++)
				{
					if (j == c)
						continue;
					if (b[j] <= t[j])
					{
						continue;
					}
					if (!pf[d].flags.mmExplore(j))
					{
						continue;
					}
					TIndexOffU topf = pf[d].topf[j], botf = pf[d].botf[j];
					ASSERT_ONLY(TIndexOffU topb = pf[d].topb[j], botb = pf[d].botb[j]);
					if (re.contains(fw, l2r_, cur5pi, cur5pf, cur5pf - cur5pi + 1 + gapadd_, topf, botf, pen_ + pen_mm))
					{
						prm.nRedSkip++;
						continue;
					}
					prm.nRedFail++;
					TIndexOffU width = b[j] - t[j];
					Edit edit((uint32_t)off5p, (int)("ACGTN"[j]), (int)("ACGTN"[c]), EDIT_TYPE_MM);
					DescentPriority pri(pen_ + pen_mm, depth, width, rootpri);
					assert(topf != 0 || botf != 0);
					assert(topb != 0 || botb != 0);
					assert_eq(botb - topb, botf - topf);
					edge.init(edit, off5p, pri, d
#ifndef NDEBUG
							  ,
							  d, topf, botf, topb, botb
#endif
					);
					out_.update(edge);
					nout++;
				}
			}
			bool gapsAllowed = (off5p >= (size_t)sc.gapbar && off3p >= (size_t)sc.gapbar);
			if (gapsAllowed)
			{
				assert_gt(depth, 0);
				size_t totwidth = (b[0] - t[0]) +
								  (b[1] - t[1]) +
								  (b[2] - t[2]) +
								  (b[3] - t[3]);
				assert(c > 3 || b[c] - t[c] <= totwidth);
				bool allmatch = c < 4 && (totwidth == (b[c] - t[c]));
				bool rdex = false, rfex = false;
				size_t cur5pi_i = cur5pi, cur5pf_i = cur5pf;
				if (toward3p)
				{
					cur5pf_i--;
				}
				else
				{
					cur5pi_i++;
				}
				if (off5p == off5p_i_ && edit_.inited())
				{
					if (pen_rdg_ex <= diff && edit_.isReadGap())
					{
						rdex = true;
						for (int j = 0; j < 4; j++)
						{
							if (b[j] <= t[j])
							{
								continue;
							}
							if (!pf[d].flags.rdgExplore(j))
							{
								continue;
							}
							TIndexOffU topf = pf[d].topf[j], botf = pf[d].botf[j];
							ASSERT_ONLY(TIndexOffU topb = pf[d].topb[j], botb = pf[d].botb[j]);
							assert(topf != 0 || botf != 0);
							assert(topb != 0 || botb != 0);
							if (re.contains(fw, l2r_, cur5pi_i, cur5pf_i, cur5pf - cur5pi + 1 + gapadd_, topf, botf, pen_ + pen_rdg_ex))
							{
								prm.nRedSkip++;
								continue;
							}
							prm.nRedFail++;
							TIndexOffU width = b[j] - t[j];
							uint32_t off = (uint32_t)off5p + (toward3p ? 0 : 1);
							Edit edit(off, (int)("ACGTN"[j]), '-', EDIT_TYPE_READ_GAP);
							assert(edit.pos2 != std::numeric_limits<uint32_t>::max());
							edit.pos2 = edit_.pos2 + (toward3p ? 1 : -1);
							DescentPriority pri(pen_ + pen_rdg_ex, depth, width, rootpri);
							assert(topf != 0 || botf != 0);
							assert(topb != 0 || botb != 0);
							assert_eq(botb - topb, botf - topf);
							edge.init(edit, off5p, pri, d
#ifndef NDEBUG
									  ,
									  d,
									  topf, botf, topb, botb
#endif
							);
							out_.update(edge);
							nout++;
						}
					}
					if (pen_rfg_ex <= diff && edit_.isRefGap())
					{
						rfex = true;
						if (pf[d].flags.rfgExplore())
						{
							TIndexOffU topf = l2r_ ? topp : top;
							TIndexOffU botf = l2r_ ? botp : bot;
							ASSERT_ONLY(TIndexOffU topb = l2r_ ? top : topp);
							ASSERT_ONLY(TIndexOffU botb = l2r_ ? bot : botp);
							assert(topf != 0 || botf != 0);
							assert(topb != 0 || botb != 0);
							size_t nrefal = cur5pf - cur5pi + gapadd_;
							if (!re.contains(fw, l2r_, cur5pi, cur5pf, nrefal, topf, botf, pen_ + pen_rfg_ex))
							{
								TIndexOffU width = bot - top;
								Edit edit((uint32_t)off5p, '-', (int)("ACGTN"[c]), EDIT_TYPE_REF_GAP);
								DescentPriority pri(pen_ + pen_rfg_ex, depth, width, rootpri);
								assert(topf != 0 || botf != 0);
								assert(topb != 0 || botb != 0);
								edge.init(edit, off5p, pri, d
#ifndef NDEBUG
										  ,
										  (d == posid_) ? std::numeric_limits<size_t>::max() : (d - 1),
										  topf, botf, topb, botb
#endif
								);
								out_.update(edge);
								nout++;
								prm.nRedFail++;
							}
							else
							{
								prm.nRedSkip++;
							}
						}
					}
				}
				if (!allmatch && pen_rdg_op <= diff && !rdex)
				{
					for (int j = 0; j < 4; j++)
					{
						if (b[j] <= t[j])
						{
							continue;
						}
						if (!pf[d].flags.rdgExplore(j))
						{
							continue;
						}
						TIndexOffU topf = pf[d].topf[j], botf = pf[d].botf[j];
						ASSERT_ONLY(TIndexOffU topb = pf[d].topb[j], botb = pf[d].botb[j]);
						assert(topf != 0 || botf != 0);
						assert(topb != 0 || botb != 0);
						if (re.contains(fw, l2r_, cur5pi_i, cur5pf_i, cur5pf - cur5pi + 1 + gapadd_, topf, botf, pen_ + pen_rdg_op))
						{
							prm.nRedSkip++;
							continue;
						}
						prm.nRedFail++;
						TIndexOffU width = b[j] - t[j];
						uint32_t off = (uint32_t)off5p + (toward3p ? 0 : 1);
						Edit edit(off, (int)("ACGTN"[j]), '-', EDIT_TYPE_READ_GAP);
						assert(edit.pos2 != std::numeric_limits<uint32_t>::max());
						DescentPriority pri(pen_ + pen_rdg_op, depth, width, rootpri);
						assert(topf != 0 || botf != 0);
						assert(topb != 0 || botb != 0);
						assert_eq(botb - topb, botf - topf);
						edge.init(edit, off5p, pri, d
#ifndef NDEBUG
								  ,
								  d, topf, botf, topb, botb
#endif
						);
						out_.update(edge);
						nout++;
					}
				}
				if (!allmatch && pen_rfg_op <= diff && !rfex)
				{
					if (pf[d].flags.rfgExplore())
					{
						TIndexOffU topf = l2r_ ? topp : top;
						TIndexOffU botf = l2r_ ? botp : bot;
						ASSERT_ONLY(TIndexOffU topb = l2r_ ? top : topp);
						ASSERT_ONLY(TIndexOffU botb = l2r_ ? bot : botp);
						assert(topf != 0 || botf != 0);
						assert(topb != 0 || botb != 0);
						size_t nrefal = cur5pf - cur5pi + gapadd_;
						if (!re.contains(fw, l2r_, cur5pi, cur5pf, nrefal, topf, botf, pen_ + pen_rfg_op))
						{
							TIndexOffU width = bot - top;
							Edit edit((uint32_t)off5p, '-', (int)("ACGTN"[c]), EDIT_TYPE_REF_GAP);
							DescentPriority pri(pen_ + pen_rfg_op, depth, width, rootpri);
							assert(topf != 0 || botf != 0);
							assert(topb != 0 || botb != 0);
							edge.init(edit, off5p, pri, d
#ifndef NDEBUG
									  ,
									  (d == posid_) ? std::numeric_limits<size_t>::max() : (d - 1),
									  topf, botf, topb, botb
#endif
							);
							out_.update(edge);
							nout++;
							prm.nRedFail++;
						}
						else
						{
							prm.nRedSkip++;
						}
					}
				}
			}
		}
		d++;
		depth++;
		assert_leq(depth, al5pf_ - al5pi_ + 2);
		if (toward3p)
		{
			if (off3p == 0)
			{
				break;
			}
			off5p++;
			off3p--;
			cur5pf++;
		}
		else
		{
			if (off5p == 0)
			{
				break;
			}
			off3p++;
			off5p--;
			cur5pi--;
		}
		if (off5p >= al5pi_ - extrai && off5p <= al5pf_ + extraf)
		{
			assert_range(0, 3, c);
			top = t[c], topp = tp[c];
			bot = b[c], botp = bp[c];
			assert_eq(bot - top, botp - topp);
		}
	}
	lastRecalc_ = (nout <= 5);
	out_.best1.updateFlags(pf);
	out_.best2.updateFlags(pf);
	out_.best3.updateFlags(pf);
	out_.best4.updateFlags(pf);
	out_.best5.updateFlags(pf);
	return nout;
}
void Descent::print(
	std::ostream *os,
	const char *prefix,
	const Read &q,
	size_t trimLf,
	size_t trimRg,
	bool fw,
	const EList<Edit> &edits,
	size_t ei,
	size_t en,
	BTDnaString &rf) const
{
	const BTDnaString &read = fw ? q.patFw : q.patRc;
	size_t eidx = ei;
	if (os != NULL)
	{
		*os << prefix;
	}
	for (size_t i = 0; i < read.length(); i++)
	{
		if (i < trimLf || i >= read.length() - trimRg)
		{
			if (os != NULL)
			{
				*os << (char)tolower(read.toChar(i));
			}
			continue;
		}
		bool del = false, mm = false;
		while (eidx < ei + en && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				if (os != NULL)
				{
					*os << '-';
				}
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				if (os != NULL)
				{
					*os << read.toChar(i);
				}
			}
			else
			{
				mm = true;
				assert(edits[eidx].isMismatch());
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				if (os != NULL)
				{
					*os << (char)edits[eidx].qchr;
				}
			}
			eidx++;
		}
		if (!del && !mm)
		{
			if (os != NULL)
			{
				*os << read.toChar(i);
			}
		}
	}
	if (os != NULL)
	{
		*os << endl;
		*os << prefix;
	}
	eidx = ei;
	for (size_t i = 0; i < read.length(); i++)
	{
		if (i < trimLf || i >= read.length() - trimRg)
		{
			if (os != NULL)
			{
				*os << ' ';
			}
			continue;
		}
		bool del = false, mm = false;
		while (eidx < ei + en && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				if (os != NULL)
				{
					*os << ' ';
				}
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				if (os != NULL)
				{
					*os << ' ';
				}
			}
			else
			{
				mm = true;
				assert(edits[eidx].isMismatch());
				if (os != NULL)
				{
					*os << ' ';
				}
			}
			eidx++;
		}
		if (!del && !mm && os != NULL)
		{
			*os << '|';
		}
	}
	if (os != NULL)
	{
		*os << endl;
		*os << prefix;
	}
	eidx = ei;
	for (size_t i = 0; i < read.length(); i++)
	{
		if (i < trimLf || i >= read.length() - trimRg)
		{
			if (os != NULL)
			{
				*os << ' ';
			}
			continue;
		}
		bool del = false, mm = false;
		while (eidx < ei + en && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				rf.appendChar((char)edits[eidx].chr);
				if (os != NULL)
				{
					*os << (char)edits[eidx].chr;
				}
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				if (os != NULL)
				{
					*os << '-';
				}
			}
			else
			{
				mm = true;
				assert(edits[eidx].isMismatch());
				rf.appendChar((char)edits[eidx].chr);
				if (os != NULL)
				{
					*os << (char)edits[eidx].chr;
				}
			}
			eidx++;
		}
		if (!del && !mm)
		{
			rf.append(read[i]);
			if (os != NULL)
			{
				*os << read.toChar(i);
			}
		}
	}
	if (os != NULL)
	{
		*os << endl;
	}
}
bool Descent::bounce(
	const Read &q,
	TIndexOffU topf, TIndexOffU botf, TIndexOffU topb, TIndexOffU botb, const Ebwt &ebwtFw, const Ebwt &ebwtBw, const Scoring &sc, TAlScore minsc, TAlScore maxpen, DescentRedundancyChecker &re, EFactory<Descent> &df, EFactory<DescentPos> &pf, const EList<DescentRoot> &rs, const EList<DescentConfig> &cs, EHeap<TDescentPair> &heap, DescentAlignmentSink &alsink, DescentMetrics &met, PerReadMetrics &prm)
{
	assert_gt(botf, topf);
	assert(al5pi_ == 0 || al5pf_ == q.length() - 1);
	assert(!(al5pi_ == 0 && al5pf_ == q.length() - 1));
	size_t dfsz = df.size();
	size_t pfsz = pf.size();
	TDescentId id = df.alloc();
	Edit e_null;
	assert(!e_null.inited());
	bool succ = df[id].init(
		q,
		rid_, sc, minsc, maxpen, al5pi_, al5pf_, topf, botf, topb, botb, !l2r_, id, descid_, pen_, e_null, ebwtFw, ebwtBw, re, df, pf, rs, cs, heap, alsink, met, prm);
	if (!succ)
	{
		df.resize(dfsz);
		pf.resize(pfsz);
	}
	return succ;
}
void Descent::followBestOutgoing(
	const Read &q,
	const Ebwt &ebwtFw, const Ebwt &ebwtBw, const Scoring &sc, TAlScore minsc, TAlScore maxpen, DescentRedundancyChecker &re, EFactory<Descent> &df, EFactory<DescentPos> &pf, const EList<DescentRoot> &rs, const EList<DescentConfig> &cs, EHeap<TDescentPair> &heap, DescentAlignmentSink &alsink, DescentMetrics &met, PerReadMetrics &prm)
{
	assert(q.repOk());
	assert(!empty());
	assert(!out_.empty());
	while (!out_.empty())
	{
		DescentPriority best = out_.bestPri();
		DescentEdge e = out_.rotate();
		TReadOff al5pi_new = al5pi_, al5pf_new = al5pf_;
		bool fw = rs[rid_].fw;
		bool toward3p = (l2r_ == fw);
		TReadOff edoff = e.off5p;
		assert_leq(edoff, al5pf_ + 1);
		assert_geq(edoff + 1, al5pi_);
		if (out_.empty())
		{
			if (!lastRecalc_)
			{
				recalcOutgoing(q, sc, minsc, maxpen, re, pf, rs, cs, prm);
				if (empty())
				{
					break;
				}
			}
			else
			{
				assert(empty());
			}
		}
		TReadOff doff;
		int chr = asc2dna[e.e.chr];
		bool hitEnd = false;
		bool done = false;
		if (toward3p)
		{
			al5pf_new = doff = edoff;
			if (e.e.isReadGap())
			{
				assert_gt(al5pf_new, 0);
				al5pf_new--;
			}
			assert_lt(al5pf_new, q.length());
			hitEnd = (al5pf_new == q.length() - 1);
			done = (hitEnd && al5pi_new == 0);
			assert_geq(doff, off5p_i_);
			doff = doff - off5p_i_;
			assert_leq(doff, len_);
		}
		else
		{
			al5pi_new = doff = edoff;
			if (e.e.isReadGap())
			{
				al5pi_new++;
			}
			hitEnd = (al5pi_new == 0);
			done = (hitEnd && al5pf_new == q.length() - 1);
			assert_geq(off5p_i_, doff);
			doff = off5p_i_ - doff;
			assert_leq(doff, len_);
		}
		bool l2r = l2r_;
		if (!done && hitEnd)
		{
			l2r = !l2r;
		}
		size_t dfsz = df.size();
		size_t pfsz = pf.size();
		TIndexOffU topf, botf, topb, botb;
		size_t d = posid_ + doff;
		if (e.e.isRefGap())
		{
			d--;
			if (doff == 0)
			{
				topf = topf_;
				botf = botf_;
				topb = topb_;
				botb = botb_;
				d = std::numeric_limits<size_t>::max();
				assert_eq(botf - topf, botb - topb);
			}
			else
			{
				assert_gt(al5pf_new, 0);
				assert_gt(d, 0);
				chr = pf[d].c;
				assert(pf[d].inited());
				assert_range(0, 3, chr);
				topf = pf[d].topf[chr];
				botf = pf[d].botf[chr];
				topb = pf[d].topb[chr];
				botb = pf[d].botb[chr];
				assert_eq(botf - topf, botb - topb);
			}
		}
		else
		{
			assert(pf[d].inited());
			topf = pf[d].topf[chr];
			botf = pf[d].botf[chr];
			topb = pf[d].topb[chr];
			botb = pf[d].botb[chr];
			assert_eq(botf - topf, botb - topb);
		}
		assert_eq(d, e.d);
		assert_eq(topf, e.topf);
		assert_eq(botf, e.botf);
		assert_eq(topb, e.topb);
		assert_eq(botb, e.botb);
		if (done)
		{
			alsink.reportAlignment(
				q,
				ebwtFw, ebwtBw, topf, botf, topb, botb, descid_, rid_, e.e, best.pen, df, pf, rs, cs);
			assert(alsink.repOk());
			return;
		}
		assert(al5pi_new != 0 || al5pf_new != q.length() - 1);
		TDescentId id = df.alloc();
		bool succ = df[id].init(
			q,
			rid_, sc, minsc, maxpen, al5pi_new, al5pf_new, topf, botf, topb, botb, l2r, id, descid_, best.pen, e.e, ebwtFw, ebwtBw, re, df, pf, rs, cs, heap, alsink, met, prm);
		if (!succ)
		{
			df.resize(dfsz);
			pf.resize(pfsz);
		}
		break;
	}
	if (!empty())
	{
		heap.insert(make_pair(out_.bestPri(), descid_));
	}
}
void Descent::nextLocsBi(
	const Ebwt &ebwtFw,
	const Ebwt &ebwtBw, SideLocus &tloc, SideLocus &bloc, TIndexOffU topf, TIndexOffU botf, TIndexOffU topb, TIndexOffU botb)
{
	assert_gt(botf, 0);
	if (l2r_)
	{
		if (botb - topb == 1)
		{
			tloc.initFromRow(topb, ebwtBw.eh(), ebwtBw.ebwt());
			bloc.invalidate();
		}
		else
		{
			SideLocus::initFromTopBot(
				topb, botb, ebwtBw.eh(), ebwtBw.ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	}
	else
	{
		if (botf - topf == 1)
		{
			tloc.initFromRow(topf, ebwtFw.eh(), ebwtFw.ebwt());
			bloc.invalidate();
		}
		else
		{
			SideLocus::initFromTopBot(
				topf, botf, ebwtFw.eh(), ebwtFw.ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	}
	assert(botf - topf == 1 || bloc.valid());
	assert(botf - topf > 1 || !bloc.valid());
}
bool Descent::followMatches(
	const Read &q,
	const Scoring &sc, const Ebwt &ebwtFw, const Ebwt &ebwtBw, DescentRedundancyChecker &re, EFactory<Descent> &df, EFactory<DescentPos> &pf, const EList<DescentRoot> &rs, const EList<DescentConfig> &cs, EHeap<TDescentPair> &heap, DescentAlignmentSink &alsink, DescentMetrics &met, PerReadMetrics &prm, bool &branches, bool &hitEnd, bool &done, TReadOff &off5p_i, TIndexOffU &topf_bounce, TIndexOffU &botf_bounce, TIndexOffU &topb_bounce, TIndexOffU &botb_bounce)
{
	size_t nobranchDepth = 20;
	bool stopOnN = true;
	assert(q.repOk());
	assert(repOk(&q));
	assert_eq(ebwtFw.eh().ftabChars(), ebwtBw.eh().ftabChars());
#ifndef NDEBUG
	for (int i = 0; i < 4; i++)
	{
		assert_eq(ebwtFw.fchr()[i], ebwtBw.fchr()[i]);
	}
#endif
	SideLocus tloc, bloc;
	TIndexOffU topf = topf_, botf = botf_, topb = topb_, botb = botb_;
	bool fw = rs[rid_].fw;
	bool toward3p;
	size_t off5p;
	assert_lt(al5pi_, q.length());
	assert_lt(al5pf_, q.length());
	while (true)
	{
		toward3p = (l2r_ == fw);
		assert_geq(al5pf_, al5pi_);
		assert(al5pi_ != 0 || al5pf_ != q.length() - 1);
		if (toward3p)
		{
			if (al5pf_ == q.length() - 1)
			{
				l2r_ = !l2r_;
				continue;
			}
			if (al5pi_ == al5pf_ && root())
			{
				off5p = off5p_i = al5pi_;
			}
			else
			{
				off5p = off5p_i = (al5pf_ + 1);
			}
		}
		else
		{
			if (al5pi_ == 0)
			{
				l2r_ = !l2r_;
				continue;
			}
			assert_gt(al5pi_, 0);
			if (al5pi_ == al5pf_ && root())
			{
				off5p = off5p_i = al5pi_;
			}
			else
			{
				off5p = off5p_i = (al5pi_ - 1);
			}
		}
		break;
	}
	size_t off3p = q.length() - off5p - 1;
	assert_lt(off5p, q.length());
	assert_lt(off3p, q.length());
	bool firstPos = true;
	assert_eq(0, len_);
	size_t nalloc = 0;
	branches = false;
	hitEnd = false;
	done = false;
	if (root())
	{
		assert_eq(al5pi_, al5pf_);
		int ftabLen = ebwtFw.eh().ftabChars();
		bool ftabFits = true;
		if (toward3p && ftabLen + off5p > q.length())
		{
			ftabFits = false;
		}
		else if (!toward3p && off5p < (size_t)ftabLen)
		{
			ftabFits = false;
		}
		bool useFtab = ftabLen > 1 && (size_t)ftabLen <= nobranchDepth && ftabFits;
		bool ftabFailed = false;
		if (useFtab)
		{
			prm.nFtabs++;
			size_t off_r2l = fw ? off5p : q.length() - off5p - 1;
			if (l2r_)
			{
			}
			else
			{
				assert_geq((int)off_r2l, ftabLen - 1);
				off_r2l -= (ftabLen - 1);
			}
			bool ret = ebwtFw.ftabLoHi(fw ? q.patFw : q.patRc, off_r2l,
									   false,
									   topf, botf);
			if (!ret)
			{
				ftabFailed = true;
			}
			else
			{
				if (botf - topf == 0)
				{
					return false;
				}
				int c_r2l = fw ? q.patFw[off_r2l] : q.patRc[off_r2l];
				size_t off_l2r = fw ? off5p : q.length() - off5p - 1;
				if (l2r_)
				{
				}
				else
				{
					assert_geq((int)off_l2r, ftabLen - 1);
					off_l2r -= (ftabLen - 1);
				}
				ASSERT_ONLY(bool ret2 =)
				ebwtBw.ftabLoHi(fw ? q.patFw : q.patRc, off_l2r,
								false,
								topb, botb);
				assert(ret == ret2);
				int c_l2r = fw ? q.patFw[off_l2r + ftabLen - 1] : q.patRc[off_l2r + ftabLen - 1];
				assert_eq(botf - topf, botb - topb);
				if (toward3p)
				{
					assert_geq((int)off3p, ftabLen - 1);
					off5p += ftabLen;
					off3p -= ftabLen;
				}
				else
				{
					assert_geq((int)off5p, ftabLen - 1);
					off5p -= ftabLen;
					off3p += ftabLen;
				}
				len_ += ftabLen;
				if (toward3p)
				{
					al5pf_ += (ftabLen - 1);
					if (al5pf_ == q.length() - 1)
					{
						hitEnd = true;
						done = (al5pi_ == 0);
					}
				}
				else
				{
					al5pi_ -= (ftabLen - 1);
					if (al5pi_ == 0)
					{
						hitEnd = true;
						done = (al5pf_ == q.length() - 1);
					}
				}
				size_t id = 0;
				if (firstPos)
				{
					posid_ = pf.alloc();
					pf[posid_].reset();
					firstPos = false;
					for (int i = 1; i < ftabLen; i++)
					{
						id = pf.alloc();
						pf[id].reset();
					}
				}
				else
				{
					for (int i = 0; i < ftabLen; i++)
					{
						id = pf.alloc();
						pf[id].reset();
					}
				}
				assert_eq(botf - topf, botb - topb);
				pf[id].c = l2r_ ? c_l2r : c_r2l;
				pf[id].topf[l2r_ ? c_l2r : c_r2l] = topf;
				pf[id].botf[l2r_ ? c_l2r : c_r2l] = botf;
				pf[id].topb[l2r_ ? c_l2r : c_r2l] = topb;
				pf[id].botb[l2r_ ? c_l2r : c_r2l] = botb;
				assert(pf[id].inited());
				nalloc += ftabLen;
			}
		}
		if (!useFtab || ftabFailed)
		{
			int rdc = q.getc(off5p, fw);
			if (rdc > 3)
			{
				return false;
			}
			assert_range(0, 3, rdc);
			topf = topb = ebwtFw.fchr()[rdc];
			botf = botb = ebwtFw.fchr()[rdc + 1];
			if (botf - topf == 0)
			{
				return false;
			}
			if (toward3p)
			{
				off5p++;
				off3p--;
			}
			else
			{
				off5p--;
				off3p++;
			}
			len_++;
			if (toward3p)
			{
				if (al5pf_ == q.length() - 1)
				{
					hitEnd = true;
					done = (al5pi_ == 0);
				}
			}
			else
			{
				if (al5pi_ == 0)
				{
					hitEnd = true;
					done = (al5pf_ == q.length() - 1);
				}
			}
			size_t id = 0;
			if (firstPos)
			{
				posid_ = id = pf.alloc();
				firstPos = false;
			}
			else
			{
				id = pf.alloc();
			}
			assert_eq(botf - topf, botb - topb);
			pf[id].c = rdc;
			pf[id].topf[rdc] = topf;
			pf[id].botf[rdc] = botf;
			pf[id].topb[rdc] = topb;
			pf[id].botb[rdc] = botb;
			assert(pf[id].inited());
			nalloc++;
		}
		assert_gt(botf, topf);
		assert_eq(botf - topf, botb - topb);
		if (!re.check(fw, l2r_, al5pi_, al5pf_, al5pf_ - al5pi_ + 1 + gapadd_,
					  topf, botf, pen_))
		{
			prm.nRedSkip++;
			return false;
		}
		prm.nRedFail++;
		prm.nRedIns++;
	}
	if (done)
	{
		Edit eempty;
		alsink.reportAlignment(
			q,
			ebwtFw, ebwtBw, topf, botf, topb, botb, descid_, rid_, eempty, pen_, df, pf, rs, cs);
		assert(alsink.repOk());
		return true;
	}
	else if (hitEnd)
	{
		assert(botf > 0 || topf > 0);
		assert_gt(botf, topf);
		topf_bounce = topf;
		botf_bounce = botf;
		topb_bounce = topb;
		botb_bounce = botb;
		return true;
	}
	nextLocsBi(ebwtFw, ebwtBw, tloc, bloc, topf, botf, topb, botb);
	assert(tloc.valid());
	assert(botf - topf == 1 || bloc.valid());
	assert(botf - topf > 1 || !bloc.valid());
	TIndexOffU t[4], b[4];
	TIndexOffU tp[4], bp[4];
	ASSERT_ONLY(TIndexOffU lasttot = botf - topf);
	bool fail = false;
	while (!fail && !hitEnd)
	{
		assert(!done);
		int rdc = q.getc(off5p, fw);
		int rdq = q.getq(off5p);
		assert_range(0, 4, rdc);
		assert_gt(botf, topf);
		assert(botf - topf == 1 || bloc.valid());
		assert(botf - topf > 1 || !bloc.valid());
		assert(tloc.valid());
		TIndexOffU width = botf - topf;
		bool ltr = l2r_;
		const Ebwt &ebwt = ltr ? ebwtBw : ebwtFw;
		t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
		int only = -1;
		size_t nopts = 1;
		if (bloc.valid())
		{
			if (ltr)
			{
				tp[0] = tp[1] = tp[2] = tp[3] = topf;
				bp[0] = bp[1] = bp[2] = bp[3] = botf;
			}
			else
			{
				tp[0] = tp[1] = tp[2] = tp[3] = topb;
				bp[0] = bp[1] = bp[2] = bp[3] = botb;
			}
			met.bwops++;
			met.bwops_bi++;
			prm.nSdFmops++;
			if (prm.doFmString)
			{
				prm.fmString.add(false, pen_, 1);
			}
			ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
			ASSERT_ONLY(TIndexOffU tot = (b[0] - t[0]) + (b[1] - t[1]) + (b[2] - t[2]) + (b[3] - t[3]));
			ASSERT_ONLY(TIndexOffU totp = (bp[0] - tp[0]) + (bp[1] - tp[1]) + (bp[2] - tp[2]) + (bp[3] - tp[3]));
			assert_eq(tot, totp);
			assert_leq(tot, lasttot);
			ASSERT_ONLY(lasttot = tot);
			fail = (rdc > 3 || b[rdc] <= t[rdc]);
			size_t nopts = 0;
			if (b[0] > t[0])
			{
				nopts++;
				only = 0;
			}
			if (b[1] > t[1])
			{
				nopts++;
				only = 1;
			}
			if (b[2] > t[2])
			{
				nopts++;
				only = 2;
			}
			if (b[3] > t[3])
			{
				nopts++;
				only = 3;
			}
			if (!fail && b[rdc] - t[rdc] < width)
			{
				branches = true;
			}
		}
		else
		{
			tp[0] = tp[1] = tp[2] = tp[3] = bp[0] = bp[1] = bp[2] = bp[3] = 0;
			TIndexOffU ntop = ltr ? topb : topf;
			met.bwops++;
			met.bwops_1++;
			prm.nSdFmops++;
			if (prm.doFmString)
			{
				prm.fmString.add(false, pen_, 1);
			}
			int cc = ebwt.mapLF1(ntop, tloc);
			assert_range(-1, 3, cc);
			fail = (cc != rdc);
			if (fail)
			{
				branches = true;
			}
			if (cc >= 0)
			{
				only = cc;
				t[cc] = ntop;
				b[cc] = ntop + 1;
				tp[cc] = ltr ? topf : topb;
				bp[cc] = ltr ? botf : botb;
			}
		}
		int origRdc = rdc;
		if (rdc == 4)
		{
			fail = true;
		}
		else
		{
			topf = ltr ? tp[rdc] : t[rdc];
			botf = ltr ? bp[rdc] : b[rdc];
			topb = ltr ? t[rdc] : tp[rdc];
			botb = ltr ? b[rdc] : bp[rdc];
			assert_eq(botf - topf, botb - topb);
		}
		if (rdc == 4 && !stopOnN && nopts == 1)
		{
			fail = false;
			rdc = only;
			int pen = sc.n(rdq);
			assert_gt(pen, 0);
			pen_ += pen;
		}
		assert_range(0, 4, origRdc);
		assert_range(0, 4, rdc);
		TIndexOffU *tf = ltr ? tp : t;
		TIndexOffU *bf = ltr ? bp : b;
		TIndexOffU *tb = ltr ? t : tp;
		TIndexOffU *bb = ltr ? b : bp;
		if (firstPos)
		{
			posid_ = pf.alloc();
			firstPos = false;
		}
		else
		{
			pf.alloc();
		}
		nalloc++;
		pf[posid_ + len_].reset();
		pf[posid_ + len_].c = origRdc;
		for (size_t i = 0; i < 4; i++)
		{
			pf[posid_ + len_].topf[i] = tf[i];
			pf[posid_ + len_].botf[i] = bf[i];
			pf[posid_ + len_].topb[i] = tb[i];
			pf[posid_ + len_].botb[i] = bb[i];
			assert_eq(pf[posid_ + len_].botf[i] - pf[posid_ + len_].topf[i],
					  pf[posid_ + len_].botb[i] - pf[posid_ + len_].topb[i]);
		}
		if (!fail)
		{
			size_t al5pf = al5pf_, al5pi = al5pi_;
			if (toward3p)
			{
				al5pf++;
			}
			else
			{
				al5pi--;
			}
			fail = !re.check(fw, l2r_, al5pi, al5pf,
							 al5pf - al5pi + 1 + gapadd_, topf, botf, pen_);
			if (fail)
			{
				prm.nRedSkip++;
			}
			else
			{
				prm.nRedFail++;
				prm.nRedIns++;
			}
		}
		if (!fail)
		{
			len_++;
			if (toward3p)
			{
				al5pf_++;
				off5p++;
				off3p--;
				if (al5pf_ == q.length() - 1)
				{
					hitEnd = true;
					done = (al5pi_ == 0);
				}
			}
			else
			{
				assert_gt(al5pi_, 0);
				al5pi_--;
				off5p--;
				off3p++;
				if (al5pi_ == 0)
				{
					hitEnd = true;
					done = (al5pf_ == q.length() - 1);
				}
			}
		}
		if (!fail && !hitEnd)
		{
			nextLocsBi(ebwtFw, ebwtBw, tloc, bloc, tf[rdc], bf[rdc], tb[rdc], bb[rdc]);
		}
	}
	assert_geq(al5pf_, al5pi_);
	assert(!root() || al5pf_ - al5pi_ + 1 == nalloc || al5pf_ - al5pi_ + 2 == nalloc);
	assert_geq(pf.size(), nalloc);
	if (done)
	{
		Edit eempty;
		alsink.reportAlignment(
			q,
			ebwtFw, ebwtBw, topf, botf, topb, botb, descid_, rid_, eempty, pen_, df, pf, rs, cs);
		assert(alsink.repOk());
		return true;
	}
	else if (hitEnd)
	{
		assert(botf > 0 || topf > 0);
		assert_gt(botf, topf);
		topf_bounce = topf;
		botf_bounce = botf;
		topb_bounce = topb;
		botb_bounce = botb;
		return true;
	}
	assert(repOk(&q));
	assert(!hitEnd || topf_bounce > 0 || botf_bounce > 0);
	return true;
}
#ifdef ALIGNER_SEED2_MAIN
#include <string>
#include "sstring.h"
#include "aligner_driver.h"
using namespace std;
bool gReportOverhangs = true;
int main(int argc, char **argv)
{
	EList<string> strs;
	strs.push_back(string("CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA"
						  "NNNNNNNNNN"
						  "CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA"));
	bool packed = false;
	pair<Ebwt *, Ebwt *> ebwts = Ebwt::fromStrings<SString<char>>(
		strs,
		packed,
		0,
		REF_READ_REVERSE,
		Ebwt::default_bigEndian,
		Ebwt::default_lineRate,
		Ebwt::default_offRate,
		Ebwt::default_ftabChars,
		".aligner_seed2.cpp.tmp",
		Ebwt::default_useBlockwise,
		Ebwt::default_bmax,
		Ebwt::default_bmaxMultSqrt,
		Ebwt::default_bmaxDivN,
		Ebwt::default_dcv,
		Ebwt::default_seed,
		false,
		false, false);
	cerr << "main...simple...test to the seed alignment" << endl;
	ebwts.first->loadIntoMemory(0, -1, true, true, true, true, false);
	ebwts.second->loadIntoMemory(0, 1, true, true, true, true, false);
	int testnum = 0;
	for (int rc = 0; rc < 2; rc++)
	{
		for (int i = 0; i < 2; i++)
		{
			cerr << "Test " << (++testnum) << endl;
			cerr << "  Query with length greater than ftab" << endl;
			DescentMetrics mets;
			PerReadMetrics prm;
			DescentDriver dr;
			BTDnaString seq("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			if (rc)
			{
				seq.reverseComp();
				qual.reverse();
			}
			dr.initRead(
				Read("test", seq.toZBuf(), qual.toZBuf()),
				false,
				false, -30, 30);
			DescentConfig conf;
			conf.cons.init(Ebwt::default_ftabChars, 1.0);
			conf.expol = DESC_EX_NONE;
			dr.addRoot(
				conf,
				(i == 0) ? 0 : (seq.length() - 1), (i == 0) ? true : false, rc == 0, 1, 0.0f);
			Scoring sc = Scoring::base1();
			dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
			assert_eq(1, dr.sink().nrange());
			assert_eq(2, dr.sink().nelt());
		}
	}
	for (int i = 0; i < 2; i++)
	{
		cerr << "Test " << (++testnum) << endl;
		cerr << "  Query with length equal to ftab" << endl;
		DescentMetrics mets;
		PerReadMetrics prm;
		DescentDriver dr;
		BTDnaString seq("GCTATATAGC", true);
		BTString qual("ABCDEFGHIa");
		dr.initRead(
			Read("test", seq.toZBuf(), qual.toZBuf()),
			false,
			false, -30, 30);
		DescentConfig conf;
		conf.cons.init(Ebwt::default_ftabChars, 1.0);
		conf.expol = DESC_EX_NONE;
		dr.addRoot(
			conf,
			(i == 0) ? 0 : (seq.length() - 1), (i == 0) ? true : false, true, 1, 0.0f);
		Scoring sc = Scoring::base1();
		dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
		assert_eq(1, dr.sink().nrange());
		assert_eq(2, dr.sink().nelt());
	}
	for (int i = 0; i < 2; i++)
	{
		cerr << "Test " << (++testnum) << endl;
		cerr << "  Query with length less than ftab" << endl;
		DescentMetrics mets;
		PerReadMetrics prm;
		DescentDriver dr;
		BTDnaString seq("GCTATATAG", true);
		BTString qual("ABCDEFGHI");
		dr.initRead(
			Read("test", seq.toZBuf(), qual.toZBuf()),
			false,
			false, -30, 30);
		DescentConfig conf;
		conf.cons.init(Ebwt::default_ftabChars, 1.0);
		conf.expol = DESC_EX_NONE;
		dr.addRoot(
			conf,
			(i == 0) ? 0 : (seq.length() - 1), (i == 0) ? true : false, true, 1, 0.0f);
		Scoring sc = Scoring::base1();
		dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
		assert_eq(1, dr.sink().nrange());
		assert_eq(2, dr.sink().nelt());
	}
	for (int i = 0; i < 2; i++)
	{
		cerr << "Test " << (++testnum) << endl;
		cerr << "  Search root in middle of read" << endl;
		DescentMetrics mets;
		PerReadMetrics prm;
		DescentDriver dr;
		BTDnaString seq("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
		BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
		TIndexOffU top, bot;
		top = bot = 0;
		bool ret = ebwts.first->contains("GCGCTCGCATCATTTTGTGT", &top, &bot);
		cerr << ret << ", " << top << ", " << bot << endl;
		dr.initRead(
			Read("test", seq.toZBuf(), qual.toZBuf()),
			false,
			false, -30, 30);
		DescentConfig conf;
		conf.cons.init(Ebwt::default_ftabChars, 1.0);
		conf.expol = DESC_EX_NONE;
		dr.addRoot(
			conf,
			(i == 0) ? 10 : (seq.length() - 1 - 10), (i == 0) ? true : false, true, 1, 0.0f);
		Scoring sc = Scoring::base1();
		dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
		assert_eq(1, dr.sink().nrange());
		assert_eq(2, dr.sink().nelt());
	}
	delete ebwts.first;
	delete ebwts.second;
	strs.clear();
	strs.push_back(string("CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA"
						  "NNNNNNNNNN"
						  "CATGTCAGCTATATAGCG"));
	ebwts = Ebwt::fromStrings<SString<char>>(
		strs,
		packed,
		0,
		REF_READ_REVERSE,
		Ebwt::default_bigEndian,
		Ebwt::default_lineRate,
		Ebwt::default_offRate,
		Ebwt::default_ftabChars,
		".aligner_seed2.cpp.tmp",
		Ebwt::default_useBlockwise,
		Ebwt::default_bmax,
		Ebwt::default_bmaxMultSqrt,
		Ebwt::default_bmaxDivN,
		Ebwt::default_dcv,
		Ebwt::default_seed,
		false,
		false, false);
	ebwts.first->loadIntoMemory(0, -1, true, true, true, true, false);
	ebwts.second->loadIntoMemory(0, 1, true, true, true, true, false);
	{
		BTDnaString seq("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
		BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
		DescentDriver dr;
		PrioritizedRootSelector sel(
			2.0,
			SimpleFunc(SIMPLE_FUNC_CONST, 10.0, 10.0, 10.0, 10.0),
			10);
		dr.initRead(
			Read("test", seq.toZBuf(), qual.toZBuf()),
			false,
			false, -30, 30, NULL, &sel);
		dr.printRoots(std::cerr);
		assert_eq(12, dr.roots().size());
		assert_eq(652, dr.roots()[0].pri);
		assert_eq(652, dr.roots()[1].pri);
	}
	{
		BTDnaString seq("NCTATATAGCGCGCTCGCATCNTTTTGTGTGCTATATAGCGCGCTCGCATCATTTTGTGTTTAT", true);
		BTString qual("ABCDEFGHIJKLMNOPabcdefghijklmnopABCDEFGHIJKLMNOPabcdefghijklmnop");
		DescentDriver dr;
		PrioritizedRootSelector sel(
			2.0,
			SimpleFunc(SIMPLE_FUNC_CONST, 10.0, 10.0, 10.0, 10.0),
			15);
		dr.initRead(
			Read("test", seq.toZBuf(), qual.toZBuf()),
			false,
			false, -30, 30, NULL, &sel);
		dr.printRoots(std::cerr);
		assert_eq(24, dr.roots().size());
		assert_eq(1230, dr.roots()[0].pri);
		assert_eq(1230, dr.roots()[1].pri);
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString seq("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			for (size_t j = 0; j < seq.length(); j++)
			{
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query with length greater than ftab and matches exactly once" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -30, 30);
				DescentConfig conf;
				conf.cons.init(Ebwt::default_ftabChars, 1.0);
				conf.expol = DESC_EX_NONE;
				dr.addRoot(
					conf,
					j, i == 0, true, 1, 0.0f);
				Scoring sc = Scoring::base1();
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
				assert_eq(1, dr.sink().nrange());
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				assert_eq(1, dr.sink().nelt());
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString seq("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			for (size_t j = 0; j < seq.length(); j++)
			{
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query with length greater than ftab and reverse complement matches exactly once" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -30, 30);
				DescentConfig conf;
				conf.cons.init(Ebwt::default_ftabChars, 1.0);
				conf.expol = DESC_EX_NONE;
				dr.addRoot(
					conf,
					j, i == 0, true, 1, 0.0f);
				dr.addRoot(
					conf,
					j, i == 0, false, 1, 1.0f);
				Scoring sc = Scoring::base1();
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
				assert_eq(1, dr.sink().nrange());
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				assert_eq(1, dr.sink().nelt());
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString orig("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			for (size_t k = 0; k < orig.length(); k++)
			{
				BTDnaString seq = orig;
				seq.set(seq[k] ^ 3, k);
				for (size_t j = 0; j < seq.length(); j++)
				{
					size_t beg = j;
					size_t end = j + Ebwt::default_ftabChars;
					if ((i > 0 && j > 0) || j == seq.length() - 1)
					{
						if (beg < Ebwt::default_ftabChars)
						{
							beg = 0;
						}
						else
						{
							beg -= Ebwt::default_ftabChars;
						}
						end -= Ebwt::default_ftabChars;
					}
					size_t kk = k;
					if (beg <= kk && end > kk)
					{
						continue;
					}
					if ((j > kk) ? (j - kk <= 2) : (kk - j <= 2))
					{
						continue;
					}
					cerr << "Test " << (++testnum) << endl;
					cerr << "  Query with length greater than ftab and matches exactly once with 1mm" << endl;
					DescentMetrics mets;
					PerReadMetrics prm;
					DescentDriver dr;
					dr.initRead(
						Read("test", seq.toZBuf(), qual.toZBuf()),
						false,
						false, -30, 30);
					DescentConfig conf;
					conf.cons.init(0, 1.0);
					conf.expol = DESC_EX_NONE;
					dr.addRoot(
						conf,
						j, i == 0, true, 1, 0.0f);
					Scoring sc = Scoring::base1();
					dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
					assert_eq(1, dr.sink().nrange());
					assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
					assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
					cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
					assert_eq(1, dr.sink().nelt());
					last_topf = dr.sink()[0].topf;
					last_botf = dr.sink()[0].botf;
				}
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString orig("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			for (size_t k = 0; k < orig.length(); k++)
			{
				BTDnaString seq = orig;
				seq.set(4, k);
				for (size_t j = 0; j < seq.length(); j++)
				{
					size_t beg = j;
					size_t end = j + Ebwt::default_ftabChars;
					if ((i > 0 && j > 0) || j == seq.length() - 1)
					{
						if (beg < Ebwt::default_ftabChars)
						{
							beg = 0;
						}
						else
						{
							beg -= Ebwt::default_ftabChars;
						}
						end -= Ebwt::default_ftabChars;
					}
					if (beg <= k && end > k)
					{
						continue;
					}
					if ((j > k) ? (j - k <= 2) : (k - j <= 2))
					{
						continue;
					}
					cerr << "Test " << (++testnum) << endl;
					cerr << "  Query with length greater than ftab and matches exactly once with 1mm" << endl;
					DescentMetrics mets;
					PerReadMetrics prm;
					DescentDriver dr;
					dr.initRead(
						Read("test", seq.toZBuf(), qual.toZBuf()),
						false,
						false, -30, 30);
					DescentConfig conf;
					conf.cons.init(0, 1.0);
					conf.expol = DESC_EX_NONE;
					dr.addRoot(
						conf,
						j, i == 0, true, 1, 0.0f);
					Scoring sc = Scoring::base1();
					dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
					assert_eq(1, dr.sink().nrange());
					assert_eq(sc.n(40), dr.sink()[0].pen);
					assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
					assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
					cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
					assert_eq(1, dr.sink().nelt());
					last_topf = dr.sink()[0].topf;
					last_botf = dr.sink()[0].botf;
				}
			}
		}
	}
	{
		RandomSource rnd(79);
		for (int i = 0; i < 2; i++)
		{
			BTDnaString orig("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			if (i == 1)
			{
				orig.reverseComp();
				qual.reverse();
			}
			for (size_t trials = 0; trials < 100; trials++)
			{
				BTDnaString seq = orig;
				size_t ns = 10;
				for (size_t k = 0; k < ns; k++)
				{
					size_t pos = rnd.nextU32() % seq.length();
					seq.set(4, pos);
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query with a bunch of Ns" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -30, 30);
				DescentConfig conf;
				conf.cons.init(Ebwt::default_ftabChars, 1.0);
				conf.expol = DESC_EX_NONE;
				for (size_t k = 0; k < ns; k++)
				{
					size_t j = rnd.nextU32() % seq.length();
					bool ltr = (rnd.nextU2() == 0) ? true : false;
					bool fw = (rnd.nextU2() == 0) ? true : false;
					dr.addRoot(
						conf,
						j, ltr, fw, 1, 0.0f);
				}
				Scoring sc = Scoring::base1();
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
			}
		}
	}
	{
		RandomSource rnd(77);
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString orig("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			bool fwi = (i == 0);
			if (!fwi)
			{
				orig.reverseComp();
			}
			for (size_t k = 0; k < orig.length(); k++)
			{
				BTDnaString seq = orig;
				seq.set(seq[k] ^ 3, k);
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query with length greater than ftab and matches exactly once with 1mm.  Many search roots." << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -30, 30);
				DescentConfig conf;
				conf.cons.init(0, 1.0);
				conf.expol = DESC_EX_NONE;
				bool onegood = false;
				for (size_t y = 0; y < 10; y++)
				{
					size_t j = rnd.nextU32() % seq.length();
					bool ltr = (rnd.nextU2() == 0) ? true : false;
					bool fw = (rnd.nextU2() == 0) ? true : false;
					dr.addRoot(
						conf,
						(TReadOff)j, ltr, fw, 1, (float)((float)y * 1.0f));
					size_t beg = j;
					size_t end = j + Ebwt::default_ftabChars;
					if (!ltr)
					{
						if (beg < Ebwt::default_ftabChars)
						{
							beg = 0;
						}
						else
						{
							beg -= Ebwt::default_ftabChars;
						}
						end -= Ebwt::default_ftabChars;
					}
					bool good = true;
					if (fw != fwi)
					{
						good = false;
					}
					if (beg <= k && end > k)
					{
						good = false;
					}
					if ((j > k) ? (j - k <= 2) : (k - j <= 2))
					{
						good = false;
					}
					if (good)
					{
						onegood = true;
					}
				}
				if (!onegood)
				{
					continue;
				}
				Scoring sc = Scoring::base1();
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
				assert_eq(1, dr.sink().nrange());
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			for (int k = 0; k < 2; k++)
			{
				BTDnaString seq("GCTATATAGCGCGCTGCATCATTTTGTGT", true);
				BTString qual("ABCDEFGHIabcdefghiABCDEFGHIab");
				if (k == 1)
				{
					seq.reverseComp();
					qual.reverse();
				}
				assert_eq(seq.length(), qual.length());
				for (size_t j = 0; j < seq.length(); j++)
				{
					size_t beg = j;
					if (k == 1)
					{
						beg = seq.length() - beg - 1;
					}
					size_t end = beg + Ebwt::default_ftabChars;
					if ((i > 0 && j > 0) || j == seq.length() - 1)
					{
						if (beg < Ebwt::default_ftabChars)
						{
							beg = 0;
						}
						else
						{
							beg -= Ebwt::default_ftabChars;
						}
						end -= Ebwt::default_ftabChars;
					}
					assert_geq(end, beg);
					if (beg <= 15 && end >= 15)
					{
						continue;
					}
					cerr << "Test " << (++testnum) << endl;
					cerr << "  Query matches once with a read gap of length 1" << endl;
					DescentMetrics mets;
					PerReadMetrics prm;
					DescentDriver dr;
					Read q("test", seq.toZBuf(), qual.toZBuf());
					assert(q.repOk());
					dr.initRead(
						q,
						false, false, -30, 30);
					DescentConfig conf;
					conf.cons.init(0, 0.5);
					conf.expol = DESC_EX_NONE;
					dr.addRoot(
						conf,
						j, i == 0, k == 0, 1, 0.0f);
					Scoring sc = Scoring::base1();
					dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
					assert_eq(1, dr.sink().nrange());
					assert_eq(sc.readGapOpen() + 0 * sc.readGapExtend(), dr.sink()[0].pen);
					assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
					assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
					cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
					assert_eq(1, dr.sink().nelt());
					last_topf = dr.sink()[0].topf;
					last_botf = dr.sink()[0].botf;
				}
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			for (int k = 0; k < 2; k++)
			{
				BTDnaString seq("GCTATATAGCGCGC"
								"CATCATTTTGTGT",
								true);
				BTString qual("ABCDEFGHIabcde"
							  "fghiABCDEFGHI");
				if (k == 1)
				{
					seq.reverseComp();
					qual.reverse();
				}
				for (size_t j = 0; j < seq.length(); j++)
				{
					size_t beg = j;
					if (k == 1)
					{
						beg = seq.length() - beg - 1;
					}
					size_t end = beg + Ebwt::default_ftabChars;
					if ((i > 0 && j > 0) || j == seq.length() - 1)
					{
						if (beg < Ebwt::default_ftabChars)
						{
							beg = 0;
						}
						else
						{
							beg -= Ebwt::default_ftabChars;
						}
						end -= Ebwt::default_ftabChars;
					}
					if (beg <= 14 && end >= 14)
					{
						continue;
					}
					cerr << "Test " << (++testnum) << endl;
					cerr << "  Query matches once with a read gap of length 3" << endl;
					DescentMetrics mets;
					PerReadMetrics prm;
					DescentDriver dr;
					dr.initRead(
						Read("test", seq.toZBuf(), qual.toZBuf()),
						false,
						false, -30, 30);
					DescentConfig conf;
					conf.cons.init(0, 0.2);
					conf.expol = DESC_EX_NONE;
					dr.addRoot(
						conf,
						j, i == 0, k == 0, 1, 0.0f);
					Scoring sc = Scoring::base1();
					sc.setMmPen(COST_MODEL_CONSTANT, 6, 6);
					dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
					assert_eq(1, dr.sink().nrange());
					assert_eq(sc.readGapOpen() + 2 * sc.readGapExtend(), dr.sink()[0].pen);
					assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
					assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
					cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
					assert_eq(1, dr.sink().nelt());
					last_topf = dr.sink()[0].topf;
					last_botf = dr.sink()[0].botf;
				}
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString seq("GCTATATAGCGCGCA"
							"TCGCATCATTTTGTGT",
							true);
			BTString qual("ABCDEFGHIabcdef"
						  "ghiABCDEFGHIabcd");
			for (size_t j = 0; j < seq.length(); j++)
			{
				size_t beg = j;
				size_t end = j + Ebwt::default_ftabChars;
				if ((i > 0 && j > 0) || j == seq.length() - 1)
				{
					if (beg < Ebwt::default_ftabChars)
					{
						beg = 0;
					}
					else
					{
						beg -= Ebwt::default_ftabChars;
					}
					end -= Ebwt::default_ftabChars;
				}
				if (beg <= 14 && end >= 14)
				{
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a reference gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -30, 30);
				DescentConfig conf;
				conf.cons.init(1, 0.5);
				conf.expol = DESC_EX_NONE;
				dr.addRoot(
					conf,
					j, i == 0, true, 1, 0.0f);
				Scoring sc = Scoring::base1();
				sc.setMmPen(COST_MODEL_CONSTANT, 6, 6);
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
				assert_eq(1, dr.sink().nrange());
				assert_eq(sc.refGapOpen() + 0 * sc.refGapExtend(), dr.sink()[0].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString seq("GCTATATAGCGCGCATG"
							"TCGCATCATTTTGTGT",
							true);
			BTString qual("ABCDEFGHIabcdefgh"
						  "iABCDEFGHIabcdef");
			for (size_t j = 0; j < seq.length(); j++)
			{
				size_t beg = j;
				size_t end = j + Ebwt::default_ftabChars;
				if ((i > 0 && j > 0) || j == seq.length() - 1)
				{
					if (beg < Ebwt::default_ftabChars)
					{
						beg = 0;
					}
					else
					{
						beg -= Ebwt::default_ftabChars;
					}
					end -= Ebwt::default_ftabChars;
				}
				if (beg <= 14 && end >= 14)
				{
					continue;
				}
				if (beg <= 15 && end >= 15)
				{
					continue;
				}
				if (beg <= 16 && end >= 16)
				{
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a reference gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -30, 30);
				DescentConfig conf;
				conf.cons.init(1, 0.25);
				conf.expol = DESC_EX_NONE;
				dr.addRoot(
					conf,
					j, i == 0, true, 1, 0.0f);
				Scoring sc = Scoring::base1();
				sc.setMmPen(COST_MODEL_CONSTANT, 6, 6);
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
				assert_eq(1, dr.sink().nrange());
				assert_eq(sc.refGapOpen() + 2 * sc.refGapExtend(), dr.sink()[0].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString seq("CATGTCAGCT"
							"GATATAGCGCGCT"
							"GCATCAATTTGTGTGTAAAC",
							true);
			BTString qual("ABCDEFGHIa"
						  "bcdefghiACDEF"
						  "GHIabcdefghijkABCDEF");
			for (size_t j = 0; j < seq.length(); j++)
			{
				size_t beg = j;
				size_t end = j + Ebwt::default_ftabChars;
				if ((i > 0 && j > 0) || j == seq.length() - 1)
				{
					if (beg < Ebwt::default_ftabChars)
					{
						beg = 0;
					}
					else
					{
						beg -= Ebwt::default_ftabChars;
					}
					end -= Ebwt::default_ftabChars;
				}
				if (beg <= 10 && end >= 10)
				{
					continue;
				}
				if (beg <= 22 && end >= 22)
				{
					continue;
				}
				if (beg <= 30 && end >= 30)
				{
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a read gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -50, 50);
				DescentConfig conf;
				conf.cons.init(1, 0.5);
				conf.expol = DESC_EX_NONE;
				dr.addRoot(
					conf,
					j, i == 0, true, 1, 0.0f);
				Scoring sc = Scoring::base1();
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
				assert_eq(1, dr.sink().nrange());
				assert_eq(sc.readGapOpen() + sc.refGapOpen() + sc.mm((int)'d' - 33), dr.sink()[0].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
	}
	delete ebwts.first;
	delete ebwts.second;
	strs.clear();
	strs.push_back(string("CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAAC"
						  "NNNNNNNNNN"
						  "CATGTCAGCTGATATAGCGCGCTCGCATCATTTTGTGTGTAAAC"
						  "N"
						  "CATGTCAGCTATATAGCGCGCTGCATCATTTTGTGTGTAAAC"
						  "N"
						  "CATGTCAGCTATATAGCGCGCTCGCATCAATTTGTGTGTAAAC"
						  "N"
						  "CATGTCAGCTGATATAGCGCGCTGCATCAATTTGTGTGTAAAC"));
	ebwts = Ebwt::fromStrings<SString<char>>(
		strs,
		packed,
		0,
		REF_READ_REVERSE,
		Ebwt::default_bigEndian,
		Ebwt::default_lineRate,
		Ebwt::default_offRate,
		Ebwt::default_ftabChars,
		".aligner_seed2.cpp.tmp",
		Ebwt::default_useBlockwise,
		Ebwt::default_bmax,
		Ebwt::default_bmaxMultSqrt,
		Ebwt::default_bmaxDivN,
		Ebwt::default_dcv,
		Ebwt::default_seed,
		false,
		false, false);
	ebwts.first->loadIntoMemory(0, -1, true, true, true, true, false);
	ebwts.second->loadIntoMemory(0, 1, true, true, true, true, false);
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString seq("CATGTCAGCT"
							"GATATAGCGCGCT"
							"GCATCAATTTGTGTGTAAAC",
							true);
			BTString qual("ABCDEFGHIa"
						  "bcdefghiACDEF"
						  "GHIabcdefghijkABCDEF");
			for (size_t j = 0; j < seq.length(); j++)
			{
				size_t beg = j;
				size_t end = j + Ebwt::default_ftabChars;
				if ((i > 0 && j > 0) || j == seq.length() - 1)
				{
					if (beg < Ebwt::default_ftabChars)
					{
						beg = 0;
					}
					else
					{
						beg -= Ebwt::default_ftabChars;
					}
					end -= Ebwt::default_ftabChars;
				}
				if (beg <= 10 && end >= 10)
				{
					continue;
				}
				if (beg <= 22 && end >= 22)
				{
					continue;
				}
				if (beg <= 30 && end >= 30)
				{
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a read gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -50, 50);
				DescentConfig conf;
				conf.cons.init(1, 0.5);
				conf.expol = DESC_EX_NONE;
				dr.addRoot(
					conf,
					j, i == 0, true, 1, 0.0f);
				Scoring sc = Scoring::base1();
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
				assert_eq(5, dr.sink().nrange());
				assert_eq(0, dr.sink()[0].pen);
				assert_eq(min(sc.readGapOpen(), sc.refGapOpen()) + sc.mm((int)'d' - 33), dr.sink()[1].pen);
				assert_eq(max(sc.readGapOpen(), sc.refGapOpen()) + sc.mm((int)'d' - 33), dr.sink()[2].pen);
				assert_eq(sc.readGapOpen() + sc.refGapOpen(), dr.sink()[3].pen);
				assert_eq(sc.readGapOpen() + sc.refGapOpen() + sc.mm((int)'d' - 33), dr.sink()[4].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(5, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
	}
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for (int i = 0; i < 2; i++)
		{
			BTDnaString seq("CATGTCAGCT"
							"GATATAGCGCGCT"
							"GCATCAATTTGTGNGTAAAC",
							true);
			BTString qual("ABCDEFGHIa"
						  "bcdefghiACDEF"
						  "GHIabcdefghijkABCDEF");
			for (size_t j = 0; j < seq.length(); j++)
			{
				size_t beg = j;
				size_t end = j + Ebwt::default_ftabChars;
				if ((i > 0 && j > 0) || j == seq.length() - 1)
				{
					if (beg < Ebwt::default_ftabChars)
					{
						beg = 0;
					}
					else
					{
						beg -= Ebwt::default_ftabChars;
					}
					end -= Ebwt::default_ftabChars;
				}
				if (beg <= 10 && end >= 10)
				{
					continue;
				}
				if (beg <= 22 && end >= 22)
				{
					continue;
				}
				if (beg <= 30 && end >= 30)
				{
					continue;
				}
				if (beg <= 36 && end >= 36)
				{
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "Query matches with various patterns of gaps, mismatches and Ns" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				dr.initRead(
					Read("test", seq.toZBuf(), qual.toZBuf()),
					false,
					false, -50, 50);
				DescentConfig conf;
				conf.cons.init(1, 0.5);
				conf.expol = DESC_EX_NONE;
				dr.addRoot(
					conf,
					j, i == 0, true, 1, 0.0f);
				Scoring sc = Scoring::base1();
				sc.setNPen(COST_MODEL_CONSTANT, 1);
				dr.go(sc, *ebwts.first, *ebwts.second, mets, prm);
				assert_eq(5, dr.sink().nrange());
				assert_eq(sc.n(40), dr.sink()[0].pen);
				assert_eq(sc.n(40) + min(sc.readGapOpen(), sc.refGapOpen()) + sc.mm((int)'d' - 33), dr.sink()[1].pen);
				assert_eq(sc.n(40) + max(sc.readGapOpen(), sc.refGapOpen()) + sc.mm((int)'d' - 33), dr.sink()[2].pen);
				assert_eq(sc.n(40) + sc.readGapOpen() + sc.refGapOpen(), dr.sink()[3].pen);
				assert_eq(sc.n(40) + sc.readGapOpen() + sc.refGapOpen() + sc.mm((int)'d' - 33), dr.sink()[4].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(5, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
	}
	delete ebwts.first;
	delete ebwts.second;
	cerr << "DONE" << endl;
}
#endif
