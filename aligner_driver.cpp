#include "aligner_driver.h"
void PrioritizedRootSelector::select(
	const Read &q,
	const Read *qo,
	bool nofw,
	bool norc,
	EList<DescentConfig> &confs,
	EList<DescentRoot> &roots)
{
	assert_gt(landing_, 0);
	const int nPenalty = 150;
	const int endBonus = 150;
	const size_t qlen = q.length();
	int interval = rootIval_.f<int>((double)qlen);
	size_t sizeTarget = qlen - landing_ + 1;
	sizeTarget = (size_t)(ceilf((sizeTarget / (float)interval)));
	sizeTarget *= 4;
	for (int i = 0; i < 2; i++)
	{
		bool fw = (i == 0);
		scoresOrig_[i].resize(qlen);
		scores_[i].resize(qlen);
		for (size_t j = 0; j < qlen; j++)
		{
			size_t off5p = fw ? j : (qlen - j - 1);
			int c = q.getc(off5p, fw);
			int sc = q.getq(off5p) - ((c > 3) ? nPenalty : 0);
			scoresOrig_[i][j] = scores_[i][j] = sc;
		}
	}
	rootHeap_.clear();
	for (int fwi = 0; fwi < 2; fwi++)
	{
		bool fw = (fwi == 0);
		if ((fw && nofw) || (!fw && norc))
		{
			continue;
		}
		int pri = 0;
		size_t revi = qlen;
		for (size_t i = 0; i < qlen; i++)
		{
			revi--;
			pri += scoresOrig_[fwi][i];
			if (i >= landing_)
			{
				pri -= scoresOrig_[fwi][i - landing_];
			}
			if (i >= landing_ - 1 && scoresOrig_[fwi][i] > 0)
			{
				rootHeap_.insert(DescentRoot(
					fw ? i : revi,
					false, fw, landing_, qlen, pri + ((revi == 0) ? endBonus : 0)));
			}
		}
		pri = 0;
		size_t i = qlen - revi;
		for (size_t revi = 0; revi < qlen; revi++)
		{
			i--;
			pri += scoresOrig_[fwi][i];
			if (revi >= landing_)
			{
				pri -= scoresOrig_[fwi][i + landing_];
			}
			if (revi >= landing_ - 1 && scoresOrig_[fwi][i] > 0)
			{
				rootHeap_.insert(DescentRoot(
					fw ? i : revi,
					true, fw, landing_, qlen, pri + ((i == 0) ? endBonus : 0)));
			}
		}
	}
	while (roots.size() < sizeTarget)
	{
		if (rootHeap_.empty())
		{
			break;
		}
		DescentRoot r = rootHeap_.pop();
		const size_t off = r.fw ? r.off5p : (qlen - r.off5p - 1);
		int fwi = r.fw ? 0 : 1;
		int pri = 0;
		if (r.l2r)
		{
			for (size_t i = 0; i < landing_; i++)
			{
				pri += scores_[fwi][off + i];
			}
		}
		else
		{
			for (size_t i = 0; i < landing_; i++)
			{
				pri += scores_[fwi][off - i];
			}
		}
		if ((r.l2r && (off == 0)) || (!r.l2r && (off == qlen - 1)))
		{
			pri += endBonus;
		}
		if (pri == r.pri)
		{
			if (r.l2r)
			{
				for (size_t i = 0; i < landing_; i++)
				{
					float frac = ((float)i / (float)landing_);
					scores_[fwi][off + i] = (int)(scores_[fwi][off + i] * frac);
				}
			}
			else
			{
				for (size_t i = 0; i < landing_; i++)
				{
					float frac = ((float)i / (float)landing_);
					scores_[fwi][off - i] = (int)(scores_[fwi][off - i] * frac);
				}
			}
			confs.expand();
			confs.back().cons.init(landing_, consExp_);
			roots.push_back(r);
		}
		else
		{
			assert_gt(roots.size(), 0);
			r.pri = pri;
			rootHeap_.insert(r);
		}
	}
	assert(!roots.empty());
}
void IntervalRootSelector::select(
	const Read &q,
	const Read *qo,
	bool nofw,
	bool norc,
	EList<DescentConfig> &confs,
	EList<DescentRoot> &roots)
{
	int interval = rootIval_.f<int>((double)q.length());
	if (qo != NULL)
	{
		interval = (int)(interval * 1.2 + 0.5);
	}
	float pri = 0.0f;
	for (int fwi = 0; fwi < 2; fwi++)
	{
		bool fw = (fwi == 0);
		if ((fw && nofw) || (!fw && norc))
		{
			continue;
		}
		{
			bool first = true;
			size_t i = 0;
			while (first || (i + landing_ <= q.length()))
			{
				confs.expand();
				confs.back().cons.init(landing_, consExp_);
				roots.expand();
				roots.back().init(
					i,
					true, fw, 1, q.length(), pri);
				i += interval;
				first = false;
			}
		}
		{
			bool first = true;
			size_t i = 0;
			while (first || (i + landing_ <= q.length()))
			{
				confs.expand();
				confs.back().cons.init(landing_, consExp_);
				roots.expand();
				roots.back().init(
					q.length() - i - 1,
					false, fw, 1, q.length(), pri);
				i += interval;
				first = false;
			}
		}
	}
}
int AlignerDriver::go(
	const Scoring &sc,
	const Ebwt &ebwtFw,
	const Ebwt &ebwtBw,
	const BitPairReference &ref,
	DescentMetrics &met,
	WalkMetrics &wlm,
	PerReadMetrics &prm,
	RandomSource &rnd,
	AlnSinkWrap &sink)
{
	if (paired_)
	{
		bool first1 = rnd.nextBool();
		bool first = true;
		DescentStoppingConditions stopc1 = stop_;
		DescentStoppingConditions stopc2 = stop_;
		size_t totszIncr = (stop_.totsz + 7) / 8;
		stopc1.totsz = totszIncr;
		stopc2.totsz = totszIncr;
		while (stopc1.totsz <= stop_.totsz && stopc2.totsz <= stop_.totsz)
		{
			if (first && first1 && stopc1.totsz <= stop_.totsz)
			{
				dr1_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
				stopc1.totsz += totszIncr;
			}
			if (stopc2.totsz <= stop_.totsz)
			{
				dr2_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
				stopc2.totsz += totszIncr;
			}
			first = false;
		}
	}
	else
	{
		size_t iter = 1;
		while (true)
		{
			int ret = dr1_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
			if (ret == DESCENT_DRIVER_ALN)
			{
				cerr << iter << ". DESCENT_DRIVER_ALN" << endl;
			}
			else if (ret == DESCENT_DRIVER_MEM)
			{
				cerr << iter << ". DESCENT_DRIVER_MEM" << endl;
				break;
			}
			else if (ret == DESCENT_DRIVER_STRATA)
			{
				AlnRes res;
				alsel_.init(
					dr1_.query(),
					dr1_.sink(),
					ebwtFw,
					ref,
					rnd,
					wlm);
				while (!alsel_.done() && !sink.state().doneWithMate(true))
				{
					res.reset();
					bool ret2 = alsel_.next(
						dr1_,
						ebwtFw,
						ref,
						rnd,
						res,
						wlm,
						prm);
					if (ret2)
					{
						assert(res.matchesRef(
							dr1_.query(),
							ref,
							tmp_rf_,
							tmp_rdseq_,
							tmp_qseq_,
							raw_refbuf_,
							raw_destU32_,
							raw_matches_));
						Interval refival(res.refid(), 0, res.fw(), res.reflen());
						assert_gt(res.refExtent(), 0);
						if (gReportOverhangs &&
							!refival.containsIgnoreOrient(res.refival()))
						{
							res.clipOutside(true, 0, res.reflen());
							if (res.refExtent() == 0)
							{
								continue;
							}
						}
						assert(gReportOverhangs ||
							   refival.containsIgnoreOrient(res.refival()));
						if (!refival.overlapsIgnoreOrient(res.refival()))
						{
							continue;
						}
						if (red1_.overlap(res))
						{
							continue;
						}
						red1_.add(res);
						assert(!sink.state().doneWithMate(true));
						assert(!sink.maxed());
						if (sink.report(0, &res, NULL))
						{
							return ALDRIVER_POLICY_FULFILLED;
						}
					}
				}
				dr1_.sink().advanceStratum();
			}
			else if (ret == DESCENT_DRIVER_BWOPS)
			{
				cerr << iter << ". DESCENT_DRIVER_BWOPS" << endl;
				break;
			}
			else if (ret == DESCENT_DRIVER_DONE)
			{
				cerr << iter << ". DESCENT_DRIVER_DONE" << endl;
				break;
			}
			else
			{
				assert(false);
			}
			iter++;
		}
	}
	return ALDRIVER_EXHAUSTED_CANDIDATES;
}
