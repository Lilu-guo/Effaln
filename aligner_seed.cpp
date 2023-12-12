#include "aligner_cache.h"
#include "aligner_seed.h"
#include "search_globals.h"
#include "aln_idx.h"
using namespace std;
Constraint Constraint::exact()
{
	Constraint c;
	c.edits = c.mms = c.ins = c.dels = c.penalty = 0;
	return c;
}
Constraint Constraint::penaltyBased(int pen)
{
	Constraint c;
	c.penalty = pen;
	return c;
}
Constraint Constraint::penaltyFuncBased(const SimpleFunc &f)
{
	Constraint c;
	c.penFunc = f;
	return c;
}
Constraint Constraint::mmBased(int mms)
{
	Constraint c;
	c.mms = mms;
	c.edits = c.dels = c.ins = 0;
	return c;
}
Constraint Constraint::editBased(int edits)
{
	Constraint c;
	c.edits = edits;
	c.dels = c.ins = c.mms = 0;
	return c;
}
bool Seed::instantiate(
	const Read &read,
	const BTDnaString &seq,
	const BTString &qual,
	const Scoring &pens,
	int depth,
	int seedoffidx,
	int seedtypeidx,
	bool fw,
	InstantiatedSeed &is) const
{
	assert(overall != NULL);
	int seedlen = len;
	if ((int)read.length() < seedlen)
	{
		seedlen = (int)read.length();
	}
	assert_gt(seedlen, 0);
	is.steps.resize(seedlen);
	is.zones.resize(seedlen);
	switch (type)
	{
	case SEED_TYPE_EXACT:
	{
		for (int k = 0; k < seedlen; k++)
		{
			is.steps[k] = -(seedlen - k);
			is.zones[k].first = is.zones[k].second = 0;
		}
		break;
	}
	case SEED_TYPE_LEFT_TO_RIGHT:
	{
		for (int k = 0; k < seedlen; k++)
		{
			is.steps[k] = k + 1;
			is.zones[k].first = is.zones[k].second = ((k < (seedlen + 1) / 2) ? 0 : 1);
		}
		is.zones[seedlen - 1].first = is.zones[seedlen - 1].second = -1;
		break;
	}
	case SEED_TYPE_RIGHT_TO_LEFT:
	{
		for (int k = 0; k < seedlen; k++)
		{
			is.steps[k] = -(seedlen - k);
			is.zones[k].first = ((k < seedlen / 2) ? 0 : 1);
			is.zones[k].second = ((k < (seedlen + 1) / 2 + 1) ? 0 : 1);
		}
		is.zones[seedlen - 1].first = is.zones[seedlen - 1].second = -1;
		break;
	}
	case SEED_TYPE_INSIDE_OUT:
	{
		int step = 0;
		for (int k = (seedlen + 3) / 4; k < seedlen - (seedlen / 4); k++)
		{
			is.zones[step].first = is.zones[step].second = 0;
			is.steps[step++] = k + 1;
		}
		for (int k = seedlen - (seedlen / 4); k < seedlen; k++)
		{
			is.zones[step].first = is.zones[step].second = 1;
			is.steps[step++] = k + 1;
		}
		is.zones[step - 1].first = is.zones[step - 1].second = -1;
		for (int k = ((seedlen + 3) / 4) - 1; k >= 0; k--)
		{
			is.zones[step].first = is.zones[step].second = 2;
			is.steps[step++] = -(k + 1);
		}
		assert_eq(2, is.zones[step - 1].first);
		is.zones[step - 1].first = is.zones[step - 1].second = -2;
		assert_eq(seedlen, step);
		break;
	}
	default:
		throw 1;
	}
	for (int i = 0; i < 3; i++)
	{
		is.cons[i] = zones[i];
		is.cons[i].instantiate(read.length());
	}
	is.overall = *overall;
	is.overall.instantiate(read.length());
	bool streak = true;
	is.maxjump = 0;
	bool ret = true;
	bool ltr = (is.steps[0] > 0);
	for (size_t i = 0; i < is.steps.size(); i++)
	{
		assert_neq(0, is.steps[i]);
		int off = is.steps[i];
		off = abs(off) - 1;
		Constraint &cons = is.cons[abs(is.zones[i].first)];
		int c = seq[off];
		assert_range(0, 4, c);
		int q = qual[off];
		if (ltr != (is.steps[i] > 0) ||
			is.zones[i].first != 0 ||
			is.zones[i].second != 0)
		{
			streak = false;
		}
		if (c == 4)
		{
			if (cons.canN(q, pens))
			{
				cons.chargeN(q, pens);
			}
			else
			{
				return false;
			}
		}
		if (streak)
			is.maxjump++;
	}
	is.seedoff = depth;
	is.seedoffidx = seedoffidx;
	is.fw = fw;
	is.s = *this;
	return ret;
}
void Seed::zeroMmSeeds(int ln, EList<Seed> &pols, Constraint &oall)
{
	oall.init();
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_EXACT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::exact();
	pols.back().zones[2] = Constraint::exact();
	pols.back().overall = &oall;
}
void Seed::oneMmSeeds(int ln, EList<Seed> &pols, Constraint &oall)
{
	oall.init();
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[2] = Constraint::exact();
	pols.back().overall = &oall;
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0;
	pols.back().zones[2] = Constraint::exact();
	pols.back().overall = &oall;
}
void Seed::twoMmSeeds(int ln, EList<Seed> &pols, Constraint &oall)
{
	oall.init();
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[2] = Constraint::exact();
	pols.back().overall = &oall;
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[1].mmsCeil = 1;
	pols.back().zones[2] = Constraint::exact();
	pols.back().overall = &oall;
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_INSIDE_OUT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0;
	pols.back().zones[2] = Constraint::mmBased(1);
	pols.back().zones[2].mmsCeil = 0;
	pols.back().overall = &oall;
}
enum
{
	SA_ACTION_TYPE_RESET = 1,
	SA_ACTION_TYPE_SEARCH_SEED,
	SA_ACTION_TYPE_FTAB,
	SA_ACTION_TYPE_FCHR,
	SA_ACTION_TYPE_MATCH,
	SA_ACTION_TYPE_EDIT
};
void SeedAligner::instantiateSeq(
	const Read &read,
	BTDnaString &seq,
	BTString &qual,
	int len,
	int depth,
	bool fw) const
{
	int seedlen = len;
	if ((int)read.length() < seedlen)
		seedlen = (int)read.length();
	seq.resize(len);
	qual.resize(len);
	for (int i = 0; i < len; i++)
	{
		seq.set(read.patFw.windowGetDna(i, fw, depth, len), i);
		qual.set(read.qual.windowGet(i, fw, depth, len), i);
	}
}
pair<int, int> SeedAligner::instantiateSeeds(
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
	pair<int, int> &instFw,
	pair<int, int> &instRc)
{
	assert(!seeds.empty());
	assert_gt(read.length(), 0);
	offIdx2off_.clear();
	int len = seeds[0].len;
#ifndef NDEBUG
	for (size_t i = 1; i < seeds.size(); i++)
	{
		assert_eq(len, seeds[i].len);
	}
#endif
	int nseeds = 1;
	if ((int)read.length() - (int)off > len)
	{
		nseeds += ((int)read.length() - (int)off - len) / per;
	}
	for (int i = 0; i < nseeds; i++)
	{
		offIdx2off_.push_back(per * i + (int)off);
	}
	pair<int, int> ret;
	ret.first = 0;
	ret.second = 0;
	sr.reset(read, offIdx2off_, nseeds);
	assert(sr.repOk(&cache.current(), true));
	for (int fwi = 0; fwi < 2; fwi++)
	{
		bool fw = (fwi == 0);
		if ((fw && nofw) || (!fw && norc))
		{
			continue;
		}
		for (int i = 0; i < nseeds; i++)
		{
			int depth = i * per + (int)off;
			int seedlen = seeds[0].len;
			instantiateSeq(
				read,
				sr.seqs(fw)[i],
				sr.quals(fw)[i],
				std::min<int>((int)seedlen, (int)read.length()),
				depth,
				fw);
			QKey qk(sr.seqs(fw)[i] ASSERT_ONLY(, tmpdnastr_));
			EList<InstantiatedSeed> &iss = sr.instantiatedSeeds(fw, i);
			for (int j = 0; j < (int)seeds.size(); j++)
			{
				iss.expand();
				assert_eq(seedlen, seeds[j].len);
				InstantiatedSeed *is = &iss.back();
				if (seeds[j].instantiate(
						read,
						sr.seqs(fw)[i],
						sr.quals(fw)[i],
						pens,
						depth,
						i,
						j,
						fw,
						*is))
				{
					ret.first++;
					if (fwi == 0)
					{
						instFw.first++;
					}
					else
					{
						instRc.first++;
					}
				}
				else
				{
					met.filteredseed++;
					iss.pop_back();
				}
			}
		}
	}
	return ret;
}
void SeedAligner::searchAllSeeds(
	const EList<Seed> &seeds,
	const Ebwt *ebwtFw,
	const Ebwt *ebwtBw,
	const Read &read,
	const Scoring &pens,
	AlignmentCacheIface &cache,
	SeedResults &sr,
	SeedSearchMetrics &met,
	PerReadMetrics &prm)
{
	assert(!seeds.empty());
	assert(ebwtFw != NULL);
	assert(ebwtFw->isInMemory());
	assert(sr.repOk(&cache.current()));
	ebwtFw_ = ebwtFw;
	ebwtBw_ = ebwtBw;
	sc_ = &pens;
	read_ = &read;
	ca_ = &cache;
	bwops_ = bwedits_ = 0;
	uint64_t possearches = 0, seedsearches = 0, intrahits = 0, interhits = 0, ooms = 0;
	for (int i = 0; i < (int)sr.numOffs(); i++)
	{
		size_t off = sr.idx2off(i);
		for (int fwi = 0; fwi < 2; fwi++)
		{
			bool fw = (fwi == 0);
			assert(sr.repOk(&cache.current()));
			EList<InstantiatedSeed> &iss = sr.instantiatedSeeds(fw, i);
			if (iss.empty())
			{
				continue;
			}
			QVal qv;
			seq_ = &sr.seqs(fw)[i];
			qual_ = &sr.quals(fw)[i];
			off_ = off;
			fw_ = fw;
			int ret = cache.beginAlign(*seq_, *qual_, qv);
			ASSERT_ONLY(hits_.clear());
			if (ret == -1)
			{
				ooms++;
				continue;
			}
			bool abort = false;
			if (ret == 0)
			{
				assert(cache.aligning());
				possearches++;
				for (size_t j = 0; j < iss.size(); j++)
				{
					assert_eq(fw, iss[j].fw);
					assert_eq(i, (int)iss[j].seedoffidx);
					s_ = &iss[j];
					if (!searchSeedBi())
					{
						ooms++;
						abort = true;
						break;
					}
					seedsearches++;
					assert(cache.aligning());
				}
				if (!abort)
				{
					qv = cache.finishAlign();
				}
			}
			else
			{
				assert_eq(1, ret);
				assert(qv.valid());
				intrahits++;
			}
			assert(abort || !cache.aligning());
			if (qv.valid())
			{
				sr.add(
					qv,
					cache.current(),
					i,
					fw);
			}
		}
	}
	prm.nSeedRanges = sr.numRanges();
	prm.nSeedElts = sr.numElts();
	prm.nSeedRangesFw = sr.numRangesFw();
	prm.nSeedRangesRc = sr.numRangesRc();
	prm.nSeedEltsFw = sr.numEltsFw();
	prm.nSeedEltsRc = sr.numEltsRc();
	prm.seedMedian = (uint64_t)(sr.medianHitsPerSeed() + 0.5);
	prm.seedMean = (uint64_t)sr.averageHitsPerSeed();
	prm.nSdFmops += bwops_;
	met.seedsearch += seedsearches;
	met.nrange += sr.numRanges();
	met.nelt += sr.numElts();
	met.possearch += possearches;
	met.intrahit += intrahits;
	met.interhit += interhits;
	met.ooms += ooms;
	met.bwops += bwops_;
	met.bweds += bwedits_;
}
bool SeedAligner::sanityPartial(
	const Ebwt *ebwtFw,
	const Ebwt *ebwtBw,
	const BTDnaString &seq,
	size_t dep,
	size_t len,
	bool do1mm,
	TIndexOffU topfw,
	TIndexOffU botfw,
	TIndexOffU topbw,
	TIndexOffU botbw)
{
	tmpdnastr_.clear();
	for (size_t i = dep; i < len; i++)
	{
		tmpdnastr_.append(seq[i]);
	}
	TIndexOffU top_fw = 0, bot_fw = 0;
	ebwtFw->contains(tmpdnastr_, &top_fw, &bot_fw);
	assert_eq(top_fw, topfw);
	assert_eq(bot_fw, botfw);
	if (do1mm && ebwtBw != NULL)
	{
		tmpdnastr_.reverse();
		TIndexOffU top_bw = 0, bot_bw = 0;
		ebwtBw->contains(tmpdnastr_, &top_bw, &bot_bw);
		assert_eq(top_bw, topbw);
		assert_eq(bot_bw, botbw);
	}
	return true;
}
size_t SeedAligner::exactSweep(
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
	SeedSearchMetrics &met)
{
	assert_gt(mineMax, 0);
	TIndexOffU top = 0, bot = 0;
	SideLocus tloc, bloc;
	const size_t len = read.length();
	size_t nelt = 0;
	for (int fwi = 0; fwi < 2; fwi++)
	{
		bool fw = (fwi == 0);
		if (fw && nofw)
			continue;
		if (!fw && norc)
			continue;
		const BTDnaString &seq = fw ? read.patFw : read.patRc;
		assert(!seq.empty());
		int ftabLen = ebwt.eh().ftabChars();
		size_t dep = 0;
		size_t nedit = 0;
		bool done = false;
		while (dep < len && !done)
		{
			top = bot = 0;
			size_t left = len - dep;
			assert_gt(left, 0);
			bool doFtab = ftabLen > 1 && left >= (size_t)ftabLen;
			if (doFtab)
			{
				for (size_t i = 0; i < (size_t)ftabLen; i++)
				{
					int c = seq[len - dep - 1 - i];
					if (c > 3)
					{
						doFtab = false;
						break;
					}
				}
			}
			if (doFtab)
			{
				ebwt.ftabLoHi(seq, len - dep - ftabLen, false, top, bot);
				dep += (size_t)ftabLen;
			}
			else
			{
				int c = seq[len - dep - 1];
				if (c < 4)
				{
					top = ebwt.fchr()[c];
					bot = ebwt.fchr()[c + 1];
				}
				dep++;
			}
			if (bot <= top)
			{
				nedit++;
				if (nedit >= mineMax)
				{
					if (fw)
					{
						mineFw = nedit;
					}
					else
					{
						mineRc = nedit;
					}
					break;
				}
				continue;
			}
			INIT_LOCS(top, bot, tloc, bloc, ebwt);
			while (dep < len)
			{
				int c = seq[len - dep - 1];
				if (c > 3)
				{
					top = bot = 0;
				}
				else
				{
					if (bloc.valid())
					{
						bwops_ += 2;
						top = ebwt.mapLF(tloc, c);
						bot = ebwt.mapLF(bloc, c);
					}
					else
					{
						bwops_++;
						top = ebwt.mapLF1(top, tloc, c);
						if (top == OFF_MASK)
						{
							top = bot = 0;
						}
						else
						{
							bot = top + 1;
						}
					}
				}
				if (bot <= top)
				{
					nedit++;
					if (nedit >= mineMax)
					{
						if (fw)
						{
							mineFw = nedit;
						}
						else
						{
							mineRc = nedit;
						}
						done = true;
					}
					break;
				}
				INIT_LOCS(top, bot, tloc, bloc, ebwt);
				dep++;
			}
			if (done)
			{
				break;
			}
			if (dep == len)
			{
				if (fw)
				{
					mineFw = nedit;
				}
				else
				{
					mineRc = nedit;
				}
				if (nedit == 0 && bot > top)
				{
					if (repex)
					{
						int64_t score = len * sc.match();
						if (fw)
						{
							hits.addExactEeFw(top, bot, NULL, NULL, fw, score);
							assert(ebwt.contains(seq, NULL, NULL));
						}
						else
						{
							hits.addExactEeRc(top, bot, NULL, NULL, fw, score);
							assert(ebwt.contains(seq, NULL, NULL));
						}
					}
					nelt += (bot - top);
				}
				break;
			}
			dep++;
		}
	}
	return nelt;
}
bool SeedAligner::oneMmSearch(
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
	SeedSearchMetrics &met)
{
	assert(!rep1mm || ebwtBw != NULL);
	const size_t len = read.length();
	int nceil = sc.nCeil.f<int>((double)len);
	size_t ns = read.ns();
	if (ns > 1)
	{
		return false;
	}
	else if (ns == 1 && !rep1mm)
	{
		return false;
	}
	assert_geq(len, 2);
	assert(!rep1mm || ebwtBw->eh().ftabChars() == ebwtFw->eh().ftabChars());
#ifndef NDEBUG
	if (ebwtBw != NULL)
	{
		for (int i = 0; i < 4; i++)
		{
			assert_eq(ebwtBw->fchr()[i], ebwtFw->fchr()[i]);
		}
	}
#endif
	size_t halfFw = len >> 1;
	size_t halfBw = len >> 1;
	if ((len & 1) != 0)
	{
		halfBw++;
	}
	assert_geq(halfFw, 1);
	assert_geq(halfBw, 1);
	SideLocus tloc, bloc;
	TIndexOffU t[4], b[4];
	t[0] = t[1] = t[2] = t[3] = 0;
	b[0] = b[1] = b[2] = b[3] = 0;
	TIndexOffU tp[4], bp[4];
	tp[0] = tp[1] = tp[2] = tp[3] = 0;
	bp[0] = bp[1] = bp[2] = bp[3] = 0;
	TIndexOffU top = 0, bot = 0, topp = 0, botp = 0;
	bool results = false;
	for (int fwi = 0; fwi < 2; fwi++)
	{
		bool fw = (fwi == 0);
		if (fw && nofw)
			continue;
		if (!fw && norc)
			continue;
		int lim = rep1mm ? 2 : 1;
		for (int ebwtfwi = 0; ebwtfwi < lim; ebwtfwi++)
		{
			bool ebwtfw = (ebwtfwi == 0);
			const Ebwt *ebwt = (ebwtfw ? ebwtFw : ebwtBw);
			const Ebwt *ebwtp = (ebwtfw ? ebwtBw : ebwtFw);
			assert(rep1mm || ebwt->fw());
			const BTDnaString &seq =
				(fw ? (ebwtfw ? read.patFw : read.patFwRev) : (ebwtfw ? read.patRc : read.patRcRev));
			assert(!seq.empty());
			const BTString &qual =
				(fw ? (ebwtfw ? read.qual : read.qualRev) : (ebwtfw ? read.qualRev : read.qual));
			int ftabLen = ebwt->eh().ftabChars();
			size_t nea = ebwtfw ? halfFw : halfBw;
			bool skip = false;
			for (size_t dep = 0; dep < nea; dep++)
			{
				if (seq[len - dep - 1] > 3)
				{
					skip = true;
					break;
				}
			}
			if (skip)
			{
				continue;
			}
			size_t dep = 0;
			if (ftabLen > 1 && (size_t)ftabLen <= nea)
			{
				bool rev = !ebwtfw;
				ebwt->ftabLoHi(seq, len - ftabLen, rev, top, bot);
				if (rep1mm)
				{
					ebwtp->ftabLoHi(seq, len - ftabLen, rev, topp, botp);
					assert_eq(bot - top, botp - topp);
				}
				if (bot - top == 0)
				{
					continue;
				}
				int c = seq[len - ftabLen];
				t[c] = top;
				b[c] = bot;
				tp[c] = topp;
				bp[c] = botp;
				dep = ftabLen;
			}
			else
			{
				int c = seq[len - 1];
				assert_range(0, 3, c);
				top = topp = tp[c] = ebwt->fchr()[c];
				bot = botp = bp[c] = ebwt->fchr()[c + 1];
				if (bot - top == 0)
				{
					continue;
				}
				dep = 1;
			}
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
			assert(sanityPartial(ebwt, ebwtp, seq, len - dep, len, rep1mm, top, bot, topp, botp));
			bool do_continue = false;
			for (; dep < nea; dep++)
			{
				assert_lt(dep, len);
				int rdc = seq[len - dep - 1];
				tp[0] = tp[1] = tp[2] = tp[3] = topp;
				bp[0] = bp[1] = bp[2] = bp[3] = botp;
				if (bloc.valid())
				{
					bwops_++;
					t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
					ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
					SANITY_CHECK_4TUP(t, b, tp, bp);
					top = t[rdc];
					bot = b[rdc];
					if (bot <= top)
					{
						do_continue = true;
						break;
					}
					topp = tp[rdc];
					botp = bp[rdc];
					assert(!rep1mm || bot - top == botp - topp);
				}
				else
				{
					assert_eq(bot, top + 1);
					assert(!rep1mm || botp == topp + 1);
					bwops_++;
					top = ebwt->mapLF1(top, tloc, rdc);
					if (top == OFF_MASK)
					{
						do_continue = true;
						break;
					}
					bot = top + 1;
					t[rdc] = top;
					b[rdc] = bot;
					tp[rdc] = topp;
					bp[rdc] = botp;
					assert(!rep1mm || b[rdc] - t[rdc] == bp[rdc] - tp[rdc]);
				}
				INIT_LOCS(top, bot, tloc, bloc, *ebwt);
				assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
			}
			if (do_continue)
			{
				continue;
			}
			for (; dep < len; dep++)
			{
				int rdc = seq[len - dep - 1];
				int quc = qual[len - dep - 1];
				if (rdc > 3 && nceil == 0)
				{
					break;
				}
				tp[0] = tp[1] = tp[2] = tp[3] = topp;
				bp[0] = bp[1] = bp[2] = bp[3] = botp;
				int clo = 0, chi = 3;
				bool match = true;
				if (bloc.valid())
				{
					bwops_++;
					t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
					ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
					SANITY_CHECK_4TUP(t, b, tp, bp);
					match = rdc < 4;
					top = t[rdc];
					bot = b[rdc];
					topp = tp[rdc];
					botp = bp[rdc];
				}
				else
				{
					assert_eq(bot, top + 1);
					assert(!rep1mm || botp == topp + 1);
					bwops_++;
					clo = ebwt->mapLF1(top, tloc);
					match = (clo == rdc);
					assert_range(-1, 3, clo);
					if (clo < 0)
					{
						break;
					}
					else
					{
						t[clo] = top;
						b[clo] = bot = top + 1;
					}
					bp[clo] = botp;
					tp[clo] = topp;
					assert(!rep1mm || bot - top == botp - topp);
					assert(!rep1mm || b[clo] - t[clo] == bp[clo] - tp[clo]);
					chi = clo;
				}
				if (rep1mm && (ns == 0 || rdc > 3))
				{
					for (int j = clo; j <= chi; j++)
					{
						if (j == rdc || b[j] == t[j])
						{
							continue;
						}
						size_t depm = dep + 1;
						TIndexOffU topm = t[j], botm = b[j];
						TIndexOffU topmp = tp[j], botmp = bp[j];
						assert_eq(botm - topm, botmp - topmp);
						TIndexOffU tm[4], bm[4];
						tm[0] = t[0];
						tm[1] = t[1];
						tm[2] = t[2];
						tm[3] = t[3];
						bm[0] = b[0];
						bm[1] = t[1];
						bm[2] = b[2];
						bm[3] = t[3];
						TIndexOffU tmp[4], bmp[4];
						tmp[0] = tp[0];
						tmp[1] = tp[1];
						tmp[2] = tp[2];
						tmp[3] = tp[3];
						bmp[0] = bp[0];
						bmp[1] = tp[1];
						bmp[2] = bp[2];
						bmp[3] = tp[3];
						SideLocus tlocm, blocm;
						INIT_LOCS(topm, botm, tlocm, blocm, *ebwt);
						for (; depm < len; depm++)
						{
							int rdcm = seq[len - depm - 1];
							tmp[0] = tmp[1] = tmp[2] = tmp[3] = topmp;
							bmp[0] = bmp[1] = bmp[2] = bmp[3] = botmp;
							if (blocm.valid())
							{
								bwops_++;
								tm[0] = tm[1] = tm[2] = tm[3] =
									bm[0] = bm[1] = bm[2] = bm[3] = 0;
								ebwt->mapBiLFEx(tlocm, blocm, tm, bm, tmp, bmp);
								SANITY_CHECK_4TUP(tm, bm, tmp, bmp);
								topm = tm[rdcm];
								botm = bm[rdcm];
								topmp = tmp[rdcm];
								botmp = bmp[rdcm];
								if (botm <= topm)
								{
									break;
								}
							}
							else
							{
								assert_eq(botm, topm + 1);
								assert_eq(botmp, topmp + 1);
								bwops_++;
								topm = ebwt->mapLF1(topm, tlocm, rdcm);
								if (topm == OFF_MASK)
								{
									break;
								}
								botm = topm + 1;
							}
							INIT_LOCS(topm, botm, tlocm, blocm, *ebwt);
						}
						if (depm == len)
						{
							size_t off5p = dep;
							size_t offstr = dep;
							if (fw == ebwtfw)
							{
								off5p = len - off5p - 1;
							}
							if (!ebwtfw)
							{
								offstr = len - offstr - 1;
							}
							Edit e((uint32_t)off5p, j, rdc, EDIT_TYPE_MM, false);
							results = true;
							int64_t score = (len - 1) * sc.match();
							int pen = sc.score(rdc, (int)(1 << j), quc - 33);
							score += pen;
							bool valid = true;
							if (local)
							{
								int64_t locscore_fw = 0, locscore_bw = 0;
								for (size_t i = 0; i < len; i++)
								{
									if (i == dep)
									{
										if (locscore_fw + pen <= 0)
										{
											valid = false;
											break;
										}
										locscore_fw += pen;
									}
									else
									{
										locscore_fw += sc.match();
									}
									if (len - i - 1 == dep)
									{
										if (locscore_bw + pen <= 0)
										{
											valid = false;
											break;
										}
										locscore_bw += pen;
									}
									else
									{
										locscore_bw += sc.match();
									}
								}
							}
							if (valid)
							{
								valid = score >= minsc;
							}
							if (valid)
							{
#ifndef NDEBUG
								BTDnaString &rf = tmprfdnastr_;
								rf.clear();
								edits_.clear();
								edits_.push_back(e);
								if (!fw)
									Edit::invertPoss(edits_, len, false);
								Edit::toRef(fw ? read.patFw : read.patRc, edits_, rf);
								if (!fw)
									Edit::invertPoss(edits_, len, false);
								assert_eq(len, rf.length());
								for (size_t i = 0; i < len; i++)
								{
									assert_lt((int)rf[i], 4);
								}
								ASSERT_ONLY(TIndexOffU toptmp = 0);
								ASSERT_ONLY(TIndexOffU bottmp = 0);
								assert(ebwtFw->contains(rf, &toptmp, &bottmp));
#endif
								TIndexOffU toprep = ebwtfw ? topm : topmp;
								TIndexOffU botrep = ebwtfw ? botm : botmp;
								assert_eq(toprep, toptmp);
								assert_eq(botrep, bottmp);
								hits.add1mmEe(toprep, botrep, &e, NULL, fw, score);
							}
						}
					}
				}
				if (bot > top && match)
				{
					assert_lt(rdc, 4);
					if (dep == len - 1)
					{
						if (ebwtfw && repex)
						{
							if (fw)
							{
								results = true;
								int64_t score = len * sc.match();
								hits.addExactEeFw(
									ebwtfw ? top : topp,
									ebwtfw ? bot : botp,
									NULL, NULL, fw, score);
								assert(ebwtFw->contains(seq, NULL, NULL));
							}
							else
							{
								results = true;
								int64_t score = len * sc.match();
								hits.addExactEeRc(
									ebwtfw ? top : topp,
									ebwtfw ? bot : botp,
									NULL, NULL, fw, score);
								assert(ebwtFw->contains(seq, NULL, NULL));
							}
						}
						break;
					}
					else
					{
						INIT_LOCS(top, bot, tloc, bloc, *ebwt);
						assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
					}
				}
				else
				{
					break;
				}
			}
		}
	}
	return results;
}
bool SeedAligner::searchSeedBi()
{
	return searchSeedBi(
		0, 0,
		0, 0, 0, 0,
		SideLocus(), SideLocus(),
		s_->cons[0], s_->cons[1], s_->cons[2], s_->overall,
		NULL);
}
inline void
SeedAligner::nextLocsBi(
	SideLocus &tloc,
	SideLocus &bloc,
	TIndexOffU topf,
	TIndexOffU botf,
	TIndexOffU topb,
	TIndexOffU botb,
	int step
#if 0
	, const SABWOffTrack* prevOt, 
	SABWOffTrack& ot
#endif
)
{
	assert_gt(botf, 0);
	assert(ebwtBw_ == NULL || botb > 0);
	assert_geq(step, 0);
	assert(ebwtBw_ == NULL || botf - topf == botb - topb);
	if (step == (int)s_->steps.size())
		return;
	if (s_->steps[step] > 0)
	{
		if (botb - topb == 1)
		{
			tloc.initFromRow(topb, ebwtBw_->eh(), ebwtBw_->ebwt());
			bloc.invalidate();
		}
		else
		{
			SideLocus::initFromTopBot(
				topb, botb, ebwtBw_->eh(), ebwtBw_->ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	}
	else
	{
		if (botf - topf == 1)
		{
			tloc.initFromRow(topf, ebwtFw_->eh(), ebwtFw_->ebwt());
			bloc.invalidate();
		}
		else
		{
			SideLocus::initFromTopBot(
				topf, botf, ebwtFw_->eh(), ebwtFw_->ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	}
#if 0
	if(botf-topf <= BW_OFF_TRACK_CEIL) {
		if(ot.size() == 0 && prevOt != NULL && prevOt->size() > 0) {
			ot = *prevOt;
		}
		bool ltr = s_->steps[step-1] > 0;
		int adj = abs(s_->steps[step-1])-1;
		const Ebwt* ebwt = ltr ? ebwtBw_ : ebwtFw_;
		ot.update(
			ltr ? topb : topf,    
			ltr ? botb : botf,    
			adj,                  
			ebwt->offs(),         
			ebwt->eh().offRate(), 
			NULL                  
		);
		assert_gt(ot.size(), 0);
	}
#endif
	assert(botf - topf == 1 || bloc.valid());
	assert(botf - topf > 1 || !bloc.valid());
}
bool SeedAligner::extendAndReportHit(
	TIndexOffU topf,
	TIndexOffU botf,
	TIndexOffU topb,
	TIndexOffU botb,
	uint16_t len,
	DoublyLinkedList<Edit> *prevEdit)
{
	size_t nlex = 0, nrex = 0;
	TIndexOffU t[4], b[4];
	TIndexOffU tp[4], bp[4];
	SideLocus tloc, bloc;
	if (off_ > 0)
	{
		const Ebwt *ebwt = ebwtFw_;
		assert(ebwt != NULL);
		const BTDnaString &seq = fw_ ? read_->patFw : read_->patRc;
		TIndexOffU top = topf, bot = botf;
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		tp[0] = tp[1] = tp[2] = tp[3] = topb;
		bp[0] = bp[1] = bp[2] = bp[3] = botb;
		SideLocus tloc, bloc;
		INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		for (size_t ii = off_; ii > 0; ii--)
		{
			size_t i = ii - 1;
			int rdc = seq.get(i);
			if (bloc.valid())
			{
				bwops_++;
				t[0] = t[1] = t[2] = t[3] =
					b[0] = b[1] = b[2] = b[3] = 0;
				ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
				SANITY_CHECK_4TUP(t, b, tp, bp);
				int nonz = -1;
				bool abort = false;
				for (int j = 0; j < 4; j++)
				{
					if (b[i] > t[i])
					{
						if (nonz >= 0)
						{
							abort = true;
							break;
						}
						nonz = j;
						top = t[i];
						bot = b[i];
					}
				}
				if (abort || nonz != rdc)
				{
					break;
				}
			}
			else
			{
				assert_eq(bot, top + 1);
				bwops_++;
				int c = ebwt->mapLF1(top, tloc);
				if (c != rdc)
				{
					break;
				}
				bot = top + 1;
			}
			if (++nlex == 255)
			{
				break;
			}
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		}
	}
	size_t rdlen = read_->length();
	size_t nright = rdlen - off_ - len;
	if (nright > 0 && ebwtBw_ != NULL)
	{
		const Ebwt *ebwt = ebwtBw_;
		assert(ebwt != NULL);
		const BTDnaString &seq = fw_ ? read_->patFw : read_->patRc;
		TIndexOffU top = topb, bot = botb;
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		tp[0] = tp[1] = tp[2] = tp[3] = topb;
		bp[0] = bp[1] = bp[2] = bp[3] = botb;
		INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		for (size_t i = off_ + len; i < rdlen; i++)
		{
			int rdc = seq.get(i);
			if (bloc.valid())
			{
				bwops_++;
				t[0] = t[1] = t[2] = t[3] =
					b[0] = b[1] = b[2] = b[3] = 0;
				ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
				SANITY_CHECK_4TUP(t, b, tp, bp);
				int nonz = -1;
				bool abort = false;
				for (int j = 0; j < 4; j++)
				{
					if (b[i] > t[i])
					{
						if (nonz >= 0)
						{
							abort = true;
							break;
						}
						nonz = j;
						top = t[i];
						bot = b[i];
					}
				}
				if (abort || nonz != rdc)
				{
					break;
				}
			}
			else
			{
				assert_eq(bot, top + 1);
				bwops_++;
				int c = ebwt->mapLF1(top, tloc);
				if (c != rdc)
				{
					break;
				}
				bot = top + 1;
			}
			if (++nrex == 255)
			{
				break;
			}
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		}
	}
	assert_lt(nlex, rdlen);
	assert_leq(nlex, off_);
	assert_lt(nrex, rdlen);
	return reportHit(topf, botf, topb, botb, len, prevEdit);
}
bool SeedAligner::reportHit(
	TIndexOffU topf,
	TIndexOffU botf,
	TIndexOffU topb,
	TIndexOffU botb,
	uint16_t len,
	DoublyLinkedList<Edit> *prevEdit)
{
	BTDnaString &rf = tmprfdnastr_;
	rf.clear();
	edits_.clear();
	if (prevEdit != NULL)
	{
		prevEdit->toList(edits_);
		Edit::sort(edits_);
		assert(Edit::repOk(edits_, *seq_));
		Edit::toRef(*seq_, edits_, rf);
	}
	else
	{
		rf = *seq_;
	}
	assert_eq(hits_.size(), ca_->curNumRanges());
	assert(hits_.insert(rf));
	if (!ca_->addOnTheFly(rf, topf, botf, topb, botb))
	{
		return false;
	}
	assert_eq(hits_.size(), ca_->curNumRanges());
#ifndef NDEBUG
	{
		BTDnaString rfr;
		TIndexOffU tpf, btf, tpb, btb;
		tpf = btf = tpb = btb = 0;
		assert(ebwtFw_->contains(rf, &tpf, &btf));
		if (ebwtBw_ != NULL)
		{
			rfr = rf;
			rfr.reverse();
			assert(ebwtBw_->contains(rfr, &tpb, &btb));
			assert_eq(tpf, topf);
			assert_eq(btf, botf);
			assert_eq(tpb, topb);
			assert_eq(btb, botb);
		}
	}
#endif
	return true;
}
bool SeedAligner::searchSeedBi(
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
	DoublyLinkedList<Edit> *prevEdit
#if 0
	, const SABWOffTrack* prevOt
#endif
)
{
	assert(s_ != NULL);
	const InstantiatedSeed &s = *s_;
	assert_gt(s.steps.size(), 0);
	assert(ebwtBw_ == NULL || ebwtBw_->eh().ftabChars() == ebwtFw_->eh().ftabChars());
#ifndef NDEBUG
	for (int i = 0; i < 4; i++)
	{
		assert(ebwtBw_ == NULL || ebwtBw_->fchr()[i] == ebwtFw_->fchr()[i]);
	}
#endif
	if (step == (int)s.steps.size())
	{
		assert(c0.acceptable());
		assert(c1.acceptable());
		assert(c2.acceptable());
		if (!reportHit(topf, botf, topb, botb, seq_->length(), prevEdit))
		{
			return false;
		}
		return true;
	}
#ifndef NDEBUG
	if (depth > 0)
	{
		assert(botf - topf == 1 || bloc.valid());
		assert(botf - topf > 1 || !bloc.valid());
	}
#endif
	int off;
	TIndexOffU tp[4], bp[4];
	if (step == 0)
	{
		assert(prevEdit == NULL);
		assert(!tloc.valid());
		assert(!bloc.valid());
		off = s.steps[0];
		bool ltr = off > 0;
		off = abs(off) - 1;
		int ftabLen = ebwtFw_->eh().ftabChars();
		if (ftabLen > 1 && ftabLen <= s.maxjump)
		{
			if (!ltr)
			{
				assert_geq(off + 1, ftabLen - 1);
				off = off - ftabLen + 1;
			}
			ebwtFw_->ftabLoHi(*seq_, off, false, topf, botf);
#ifdef NDEBUG
			if (botf - topf == 0)
				return true;
#endif
#ifdef NDEBUG
			if (ebwtBw_ != NULL)
			{
				topb = ebwtBw_->ftabHi(*seq_, off);
				botb = topb + (botf - topf);
			}
#else
			if (ebwtBw_ != NULL)
			{
				ebwtBw_->ftabLoHi(*seq_, off, false, topb, botb);
				assert_eq(botf - topf, botb - topb);
			}
			if (botf - topf == 0)
				return true;
#endif
			step += ftabLen;
		}
		else if (s.maxjump > 0)
		{
			int c = (*seq_)[off];
			assert_range(0, 3, c);
			topf = topb = ebwtFw_->fchr()[c];
			botf = botb = ebwtFw_->fchr()[c + 1];
			if (botf - topf == 0)
				return true;
			step++;
		}
		else
		{
			assert_eq(0, s.maxjump);
			topf = topb = 0;
			botf = botb = ebwtFw_->fchr()[4];
		}
		if (step == (int)s.steps.size())
		{
			assert(c0.acceptable());
			assert(c1.acceptable());
			assert(c2.acceptable());
			if (!reportHit(topf, botf, topb, botb, seq_->length(), prevEdit))
			{
				return false;
			}
			return true;
		}
		nextLocsBi(tloc, bloc, topf, botf, topb, botb, step);
		assert(tloc.valid());
	}
	else
		assert(prevEdit != NULL);
	assert(tloc.valid());
	assert(botf - topf == 1 || bloc.valid());
	assert(botf - topf > 1 || !bloc.valid());
	assert_geq(step, 0);
	TIndexOffU t[4], b[4];
	Constraint *zones[3] = {&c0, &c1, &c2};
	ASSERT_ONLY(TIndexOffU lasttot = botf - topf);
	for (int i = step; i < (int)s.steps.size(); i++)
	{
		assert_gt(botf, topf);
		assert(botf - topf == 1 || bloc.valid());
		assert(botf - topf > 1 || !bloc.valid());
		assert(ebwtBw_ == NULL || botf - topf == botb - topb);
		assert(tloc.valid());
		off = s.steps[i];
		bool ltr = off > 0;
		const Ebwt *ebwt = ltr ? ebwtBw_ : ebwtFw_;
		assert(ebwt != NULL);
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
		t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
		if (bloc.valid())
		{
			bwops_++;
			ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
			ASSERT_ONLY(TIndexOffU tot = (b[0] - t[0]) + (b[1] - t[1]) + (b[2] - t[2]) + (b[3] - t[3]));
			ASSERT_ONLY(TIndexOffU totp = (bp[0] - tp[0]) + (bp[1] - tp[1]) + (bp[2] - tp[2]) + (bp[3] - tp[3]));
			assert_eq(tot, totp);
			assert_leq(tot, lasttot);
			ASSERT_ONLY(lasttot = tot);
		}
		TIndexOffU *tf = ltr ? tp : t, *tb = ltr ? t : tp;
		TIndexOffU *bf = ltr ? bp : b, *bb = ltr ? b : bp;
		off = abs(off) - 1;
		bool leaveZone = s.zones[i].first < 0;
		Constraint &cons = *zones[abs(s.zones[i].first)];
		int c = (*seq_)[off];
		assert_range(0, 4, c);
		int q = (*qual_)[off];
		if (!(cons.mustMatch() && !overall.mustMatch()) || c == 4)
		{
			bool bail = false;
			if (!bloc.valid())
			{
				TIndexOffU ntop = ltr ? topb : topf;
				bwops_++;
				int cc = ebwt->mapLF1(ntop, tloc);
				assert_range(-1, 3, cc);
				if (cc < 0)
					bail = true;
				else
				{
					t[cc] = ntop;
					b[cc] = ntop + 1;
				}
			}
			if (!bail)
			{
				if ((cons.canMismatch(q, *sc_) && overall.canMismatch(q, *sc_)) || c == 4)
				{
					Constraint oldCons = cons, oldOvCons = overall;
					SideLocus oldTloc = tloc, oldBloc = bloc;
					if (c != 4)
					{
						cons.chargeMismatch(q, *sc_);
						overall.chargeMismatch(q, *sc_);
					}
					if (!leaveZone || (cons.acceptable() && overall.acceptable()))
					{
						for (int j = 0; j < 4; j++)
						{
							if (j == c || b[j] == t[j])
								continue;
							nextLocsBi(tloc, bloc, tf[j], bf[j], tb[j], bb[j], i + 1);
							int loff = off;
							if (!ltr)
								loff = (int)(s.steps.size() - loff - 1);
							assert(prevEdit == NULL || prevEdit->next == NULL);
							Edit edit(off, j, c, EDIT_TYPE_MM, false);
							DoublyLinkedList<Edit> editl;
							editl.payload = edit;
							if (prevEdit != NULL)
							{
								prevEdit->next = &editl;
								editl.prev = prevEdit;
							}
							assert(editl.next == NULL);
							bwedits_++;
							if (!searchSeedBi(
									i + 1,
									depth + 1,
									tf[j],
									bf[j],
									tb[j],
									bb[j],
									tloc,
									bloc,
									c0,
									c1,
									c2,
									overall,
									&editl))
							{
								return false;
							}
							if (prevEdit != NULL)
								prevEdit->next = NULL;
						}
					}
					else
					{
					}
					cons = oldCons;
					overall = oldOvCons;
					tloc = oldTloc;
					bloc = oldBloc;
				}
				if (cons.canGap() && overall.canGap())
				{
					throw 1;
				}
			}
		}
		if (c == 4)
		{
			return true;
		}
		if (leaveZone && (!cons.acceptable() || !overall.acceptable()))
		{
			return true;
		}
		if (!bloc.valid())
		{
			assert(ebwtBw_ == NULL || bp[c] == tp[c] + 1);
			TIndexOffU top = ltr ? topb : topf;
			bwops_++;
			t[c] = ebwt->mapLF1(top, tloc, c);
			if (t[c] == OFF_MASK)
			{
				return true;
			}
			assert_geq(t[c], ebwt->fchr()[c]);
			assert_lt(t[c], ebwt->fchr()[c + 1]);
			b[c] = t[c] + 1;
			assert_gt(b[c], 0);
		}
		assert(ebwtBw_ == NULL || bf[c] - tf[c] == bb[c] - tb[c]);
		assert_leq(bf[c] - tf[c], lasttot);
		ASSERT_ONLY(lasttot = bf[c] - tf[c]);
		if (b[c] == t[c])
		{
			return true;
		}
		topf = tf[c];
		botf = bf[c];
		topb = tb[c];
		botb = bb[c];
		if (i + 1 == (int)s.steps.size())
		{
			assert(c0.acceptable());
			assert(c1.acceptable());
			assert(c2.acceptable());
			if (!reportHit(topf, botf, topb, botb, seq_->length(), prevEdit))
			{
				return false;
			}
			return true;
		}
		nextLocsBi(tloc, bloc, tf[c], bf[c], tb[c], bb[c], i + 1);
	}
	return true;
}
#ifdef ALIGNER_SEED_MAIN
#include <getopt.h>
#include <string>
static int parseInt(const char *errmsg, const char *arg)
{
	long l;
	char *endPtr = NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL)
	{
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	throw 1;
	return -1;
}
enum
{
	ARG_NOFW = 256,
	ARG_NORC,
	ARG_MM,
	ARG_SHMEM,
	ARG_TESTS,
	ARG_RANDOM_TESTS,
	ARG_SEED
};
static const char *short_opts = "vCt";
static struct option long_opts[] = {
	{(char *)"verbose", no_argument, 0, 'v'},
	{(char *)"timing", no_argument, 0, 't'},
	{(char *)"nofw", no_argument, 0, ARG_NOFW},
	{(char *)"norc", no_argument, 0, ARG_NORC},
	{(char *)"mm", no_argument, 0, ARG_MM},
	{(char *)"shmem", no_argument, 0, ARG_SHMEM},
	{(char *)"tests", no_argument, 0, ARG_TESTS},
	{(char *)"random", required_argument, 0, ARG_RANDOM_TESTS},
	{(char *)"seed", required_argument, 0, ARG_SEED},
};
bool gNorc = false;
bool gNofw = false;
int gVerbose = 0;
int gGapBarrier = 1;
int gSnpPhred = 30;
bool gReportOverhangs = true;
extern void aligner_seed_tests();
extern void aligner_random_seed_tests(
	int num_tests,
	TIndexOffU qslo,
	TIndexOffU qshi,
	uint32_t seed);
int main(int argc, char **argv)
{
	bool useMm = false;
	bool useShmem = false;
	bool mmSweep = false;
	bool noRefNames = false;
	bool sanity = false;
	bool timing = false;
	int option_index = 0;
	int seed = 777;
	int next_option;
	do
	{
		next_option = getopt_long(
			argc, argv, short_opts, long_opts, &option_index);
		switch (next_option)
		{
		case 'v':
			gVerbose = true;
			break;
		case 't':
			timing = true;
			break;
		case ARG_NOFW:
			gNofw = true;
			break;
		case ARG_NORC:
			gNorc = true;
			break;
		case ARG_MM:
			useMm = true;
			break;
		case ARG_SHMEM:
			useShmem = true;
			break;
		case ARG_SEED:
			seed = parseInt("", optarg);
			break;
		case ARG_TESTS:
		{
			aligner_seed_tests();
			aligner_random_seed_tests(
				100,
				100,
				400,
				18);
			return 0;
		}
		case ARG_RANDOM_TESTS:
		{
			seed = parseInt("", optarg);
			aligner_random_seed_tests(
				100,
				100,
				400,
				seed);
			return 0;
		}
		case -1:
			break;
		default:
		{
			cerr << "Unknown option: " << (char)next_option << endl;
			exit(1);
		}
		}
	} while (next_option != -1);
	char *reffn;
	if (optind >= argc)
	{
		cerr << "No reference; quitting..." << endl;
		return 1;
	}
	reffn = argv[optind++];
	if (optind >= argc)
	{
		cerr << "No reads; quitting..." << endl;
		return 1;
	}
	string ebwtBase(reffn);
	BitPairReference ref(
		ebwtBase,
		false,
		sanity,
		NULL,
		NULL,
		false,
		useMm,
		useShmem,
		mmSweep,
		gVerbose,
		gVerbose);
	Timer *t = new Timer(cerr, "Time loading fw index: ", timing);
	Ebwt ebwtFw(
		ebwtBase,
		false,
		0,
		true,
		-1,
		useMm,
		useShmem,
		mmSweep,
		!noRefNames,
		false,
		true,
		true,
		NULL,
		gVerbose,
		gVerbose,
		false,
		sanity);
	delete t;
	t = new Timer(cerr, "Time loading bw index: ", timing);
	Ebwt ebwtBw(
		ebwtBase + ".rev",
		false,
		1,
		false,
		-1,
		useMm,
		useShmem,
		mmSweep,
		!noRefNames,
		false,
		true,
		false,
		NULL,
		gVerbose,
		gVerbose,
		false,
		sanity);
	delete t;
	for (int i = optind; i < argc; i++)
	{
	}
}
#endif
