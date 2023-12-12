#define TIMER_START()           \
	struct timeval tv_i, tv_f;  \
	struct timezone tz_i, tz_f; \
	size_t total_usecs;         \
	gettimeofday(&tv_i, &tz_i)
#define IF_TIMER_END()                                                         \
	gettimeofday(&tv_f, &tz_f);                                                \
	total_usecs =                                                              \
		(tv_f.tv_sec - tv_i.tv_sec) * 1000000 + (tv_f.tv_usec - tv_i.tv_usec); \
	if (total_usecs > 300000)
extern "C"
{
#include "semiWFA/gap_affine/affine_wavefront_align.h"
	extern mm_allocator_t *mm_allocator_new(const uint64_t segment_size);
	extern affine_wavefronts_t *affine_wavefronts_new_complete(const int pattern_length, const int text_length, affine_penalties_t *const penalties, wavefronts_stats_t *const wavefronts_stats, mm_allocator_t *const mm_allocator);
	extern void affine_wavefronts_align(affine_wavefronts_t *const affine_wavefronts, const char *const pattern, const int pattern_length, const char *const text, const int text_length);
	extern int edit_cigar_score_gap_affine(edit_cigar_t *const edit_cigar, affine_penalties_t *const penalties);
	extern void edit_cigar_print_pretty(FILE *const stream, const char *const pattern, const int pattern_length, const char *const text, const int text_length, edit_cigar_t *const edit_cigar, mm_allocator_t *const mm_allocator);
	extern void affine_wavefronts_delete(affine_wavefronts_t *const affine_wavefronts);
	extern void mm_allocator_delete(mm_allocator_t *const mm_allocator);
}
#include <iostream>
#include "aligner_cache.h"
#include "aligner_sw_driver.h"
#include "pe.h"
#include "dp_framer.h"
#include <stdlib.h>
#include <sys/time.h>
using namespace std;
bool SwDriver::eeSaTups(
	const Read &rd,
	SeedResults &sh,
	const Ebwt &ebwt,
	const BitPairReference &ref,
	RandomSource &rnd,
	WalkMetrics &wlm,
	SwMetrics &swmSeed,
	size_t &nelt_out,
	size_t maxelt,
	bool all)
{
	assert_eq(0, nelt_out);
	gws_.clear();
	rands_.clear();
	satpos_.clear();
	eehits_.clear();
	size_t nobj = 0;
	if (!sh.exactFwEEHit().empty())
		nobj++;
	if (!sh.exactRcEEHit().empty())
		nobj++;
	nobj += sh.mm1EEHits().size();
	nobj = min(nobj, maxelt);
	gws_.ensure(nobj);
	rands_.ensure(nobj);
	satpos_.ensure(nobj);
	eehits_.ensure(nobj);
	size_t tot = sh.exactFwEEHit().size() + sh.exactRcEEHit().size();
	bool succ = false;
	bool firstEe = true;
	bool done = false;
	if (tot > 0)
	{
		bool fwFirst = true;
#ifdef _64BIT_INDEX
		TIndexOffU rn64 = rnd.nextU64();
		TIndexOffU rn = rn64 % (uint64_t)tot;
#else
		TIndexOffU rn32 = rnd.nextU32();
		TIndexOffU rn = rn32 % (uint32_t)tot;
#endif
		if (rn >= sh.exactFwEEHit().size())
		{
			fwFirst = false;
		}
		for (int fwi = 0; fwi < 2 && !done; fwi++)
		{
			bool fw = ((fwi == 0) == fwFirst);
			EEHit hit = fw ? sh.exactFwEEHit() : sh.exactRcEEHit();
			if (hit.empty())
			{
				continue;
			}
			assert(hit.fw == fw);
			if (hit.bot > hit.top)
			{
				TIndexOffU tops[2] = {hit.top, 0};
				TIndexOffU bots[2] = {hit.bot, 0};
				TIndexOffU width = hit.bot - hit.top;
				if (nelt_out + width > maxelt)
				{
					TIndexOffU trim = (TIndexOffU)((nelt_out + width) - maxelt);
#ifdef _64BIT_INDEX
					TIndexOffU rn = rnd.nextU64() % width;
#else
					TIndexOffU rn = rnd.nextU32() % width;
#endif
					TIndexOffU newwidth = width - trim;
					if (hit.top + rn + newwidth > hit.bot)
					{
						tops[0] = hit.top + rn;
						bots[0] = hit.bot;
						tops[1] = hit.top;
						bots[1] = hit.top + newwidth - (bots[0] - tops[0]);
					}
					else
					{
						tops[0] = hit.top + rn;
						bots[0] = tops[0] + newwidth;
					}
					assert_leq(bots[0], hit.bot);
					assert_leq(bots[1], hit.bot);
					assert_geq(bots[0], tops[0]);
					assert_geq(bots[1], tops[1]);
					assert_eq(newwidth, (bots[0] - tops[0]) + (bots[1] - tops[1]));
				}
				for (int i = 0; i < 2 && !done; i++)
				{
					if (bots[i] <= tops[i])
						break;
					TIndexOffU width = bots[i] - tops[i];
					TIndexOffU top = tops[i];
					swmSeed.exranges++;
					swmSeed.exrows += width;
					if (!succ)
					{
						swmSeed.exsucc++;
						succ = true;
					}
					if (firstEe)
					{
						salistEe_.clear();
						pool_.clear();
						firstEe = false;
					}
					TSlice o(salistEe_, (TIndexOffU)salistEe_.size(), width);
					for (TIndexOffU i = 0; i < width; i++)
					{
						if (!salistEe_.add(pool_, OFF_MASK))
						{
							swmSeed.exooms++;
							return false;
						}
					}
					assert(!done);
					eehits_.push_back(hit);
					satpos_.expand();
					satpos_.back().sat.init(SAKey(), top, OFF_MASK, o);
					satpos_.back().sat.key.seq = MAX_U64;
					satpos_.back().sat.key.len = (uint32_t)rd.length();
					satpos_.back().pos.init(fw, 0, 0, (uint32_t)rd.length());
					satpos_.back().origSz = width;
					rands_.expand();
					rands_.back().init(width, all);
					gws_.expand();
					SARangeWithOffs<TSlice> sa;
					sa.topf = satpos_.back().sat.topf;
					sa.len = satpos_.back().sat.key.len;
					sa.offs = satpos_.back().sat.offs;
					gws_.back().init(
						ebwt,
						ref,
						sa,
						rnd,
						wlm);
					assert(gws_.back().repOk(sa));
					nelt_out += width;
					if (nelt_out >= maxelt)
					{
						done = true;
					}
				}
			}
		}
	}
	succ = false;
	if (!done && !sh.mm1EEHits().empty())
	{
		sh.sort1mmEe(rnd);
		size_t sz = sh.mm1EEHits().size();
		for (size_t i = 0; i < sz && !done; i++)
		{
			EEHit hit = sh.mm1EEHits()[i];
			assert(hit.repOk(rd));
			assert(!hit.empty());
			TIndexOffU tops[2] = {hit.top, 0};
			TIndexOffU bots[2] = {hit.bot, 0};
			TIndexOffU width = hit.bot - hit.top;
			if (nelt_out + width > maxelt)
			{
				TIndexOffU trim = (TIndexOffU)((nelt_out + width) - maxelt);
#ifdef _64BIT_INDEX
				TIndexOffU rn = rnd.nextU64() % width;
#else
				TIndexOffU rn = rnd.nextU32() % width;
#endif
				TIndexOffU newwidth = width - trim;
				if (hit.top + rn + newwidth > hit.bot)
				{
					tops[0] = hit.top + rn;
					bots[0] = hit.bot;
					tops[1] = hit.top;
					bots[1] = hit.top + newwidth - (bots[0] - tops[0]);
				}
				else
				{
					tops[0] = hit.top + rn;
					bots[0] = tops[0] + newwidth;
				}
				assert_leq(bots[0], hit.bot);
				assert_leq(bots[1], hit.bot);
				assert_geq(bots[0], tops[0]);
				assert_geq(bots[1], tops[1]);
				assert_eq(newwidth, (bots[0] - tops[0]) + (bots[1] - tops[1]));
			}
			for (int i = 0; i < 2 && !done; i++)
			{
				if (bots[i] <= tops[i])
					break;
				TIndexOffU width = bots[i] - tops[i];
				TIndexOffU top = tops[i];
				swmSeed.mm1ranges++;
				swmSeed.mm1rows += width;
				if (!succ)
				{
					swmSeed.mm1succ++;
					succ = true;
				}
				if (firstEe)
				{
					salistEe_.clear();
					pool_.clear();
					firstEe = false;
				}
				TSlice o(salistEe_, (TIndexOffU)salistEe_.size(), width);
				for (size_t i = 0; i < width; i++)
				{
					if (!salistEe_.add(pool_, OFF_MASK))
					{
						swmSeed.mm1ooms++;
						return false;
					}
				}
				eehits_.push_back(hit);
				satpos_.expand();
				satpos_.back().sat.init(SAKey(), top, OFF_MASK, o);
				satpos_.back().sat.key.seq = MAX_U64;
				satpos_.back().sat.key.len = (uint32_t)rd.length();
				satpos_.back().pos.init(hit.fw, 0, 0, (uint32_t)rd.length());
				satpos_.back().origSz = width;
				rands_.expand();
				rands_.back().init(width, all);
				gws_.expand();
				SARangeWithOffs<TSlice> sa;
				sa.topf = satpos_.back().sat.topf;
				sa.len = satpos_.back().sat.key.len;
				sa.offs = satpos_.back().sat.offs;
				gws_.back().init(
					ebwt,
					ref,
					sa,
					rnd,
					wlm);
				assert(gws_.back().repOk(sa));
				nelt_out += width;
				if (nelt_out >= maxelt)
				{
					done = true;
				}
			}
		}
	}
	return true;
}
void SwDriver::extend(
	const Read &rd,
	const Ebwt &ebwtFw,
	const Ebwt *ebwtBw,
	TIndexOffU topf,
	TIndexOffU botf,
	TIndexOffU topb,
	TIndexOffU botb,
	bool fw,
	size_t off,
	size_t len,
	PerReadMetrics &prm,
	size_t &nlex,
	size_t &nrex)
{
	TIndexOffU t[4], b[4];
	TIndexOffU tp[4], bp[4];
	SideLocus tloc, bloc;
	size_t rdlen = rd.length();
	size_t lim = fw ? off : rdlen - len - off;
#ifndef NDEBUG
	if (false)
	{
	}
#endif
	ASSERT_ONLY(tmp_rdseq_.reverse());
	if (lim > 0)
	{
		const Ebwt *ebwt = &ebwtFw;
		assert(ebwt != NULL);
		const BTDnaString &seq = fw ? rd.patFw : rd.patRc;
		TIndexOffU top = topf, bot = botf;
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		tp[0] = tp[1] = tp[2] = tp[3] = topb;
		bp[0] = bp[1] = bp[2] = bp[3] = botb;
		SideLocus tloc, bloc;
		INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		for (size_t ii = 0; ii < lim; ii++)
		{
			size_t i = 0;
			if (fw)
			{
				i = off - ii - 1;
			}
			else
			{
				i = rdlen - off - len - 1 - ii;
			}
			int rdc = seq.get(i);
			if (bloc.valid())
			{
				prm.nSdFmops++;
				t[0] = t[1] = t[2] = t[3] =
					b[0] = b[1] = b[2] = b[3] = 0;
				ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
				SANITY_CHECK_4TUP(t, b, tp, bp);
				int nonz = -1;
				bool abort = false;
				size_t origSz = bot - top;
				for (int j = 0; j < 4; j++)
				{
					if (b[j] > t[j])
					{
						if (nonz >= 0)
						{
							abort = true;
							break;
						}
						nonz = j;
						top = t[j];
						bot = b[j];
					}
				}
				assert_leq(bot - top, origSz);
				if (abort || (nonz != rdc && rdc <= 3) || bot - top < origSz)
				{
					break;
				}
			}
			else
			{
				assert_eq(bot, top + 1);
				prm.nSdFmops++;
				int c = ebwt->mapLF1(top, tloc);
				if (c != rdc && rdc <= 3)
				{
					break;
				}
				bot = top + 1;
			}
			ASSERT_ONLY(tmp_rdseq_.append(rdc));
			if (++nlex == 255)
			{
				break;
			}
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		}
	}
	ASSERT_ONLY(tmp_rdseq_.reverse());
	lim = fw ? rdlen - len - off : off;
	if (lim > 0 && ebwtBw != NULL)
	{
		const Ebwt *ebwt = ebwtBw;
		assert(ebwt != NULL);
		const BTDnaString &seq = fw ? rd.patFw : rd.patRc;
		TIndexOffU top = topb, bot = botb;
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		tp[0] = tp[1] = tp[2] = tp[3] = topf;
		bp[0] = bp[1] = bp[2] = bp[3] = botf;
		INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		for (size_t ii = 0; ii < lim; ii++)
		{
			size_t i;
			if (fw)
			{
				i = ii + len + off;
			}
			else
			{
				i = rdlen - off + ii;
			}
			int rdc = seq.get(i);
			if (bloc.valid())
			{
				prm.nSdFmops++;
				t[0] = t[1] = t[2] = t[3] =
					b[0] = b[1] = b[2] = b[3] = 0;
				ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
				SANITY_CHECK_4TUP(t, b, tp, bp);
				int nonz = -1;
				bool abort = false;
				size_t origSz = bot - top;
				for (int j = 0; j < 4; j++)
				{
					if (b[j] > t[j])
					{
						if (nonz >= 0)
						{
							abort = true;
							break;
						}
						nonz = j;
						top = t[j];
						bot = b[j];
					}
				}
				assert_leq(bot - top, origSz);
				if (abort || (nonz != rdc && rdc <= 3) || bot - top < origSz)
				{
					break;
				}
			}
			else
			{
				assert_eq(bot, top + 1);
				prm.nSdFmops++;
				int c = ebwt->mapLF1(top, tloc);
				if (c != rdc && rdc <= 3)
				{
					break;
				}
				bot = top + 1;
			}
			ASSERT_ONLY(tmp_rdseq_.append(rdc));
			if (++nrex == 255)
			{
				break;
			}
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		}
	}
#ifndef NDEBUG
	if (false)
	{
	}
#endif
	assert_lt(nlex, rdlen);
	assert_lt(nrex, rdlen);
	return;
}
void SwDriver::prioritizeSATups(
	const Read &read,
	SeedResults &sh,
	const Ebwt &ebwtFw,
	const Ebwt *ebwtBw,
	const BitPairReference &ref,
	int seedmms,
	size_t maxelt,
	bool doExtend,
	bool lensq,
	bool szsq,
	size_t nsm,
	AlignmentCacheIface &ca,
	RandomSource &rnd,
	WalkMetrics &wlm,
	PerReadMetrics &prm,
	size_t &nelt_out,
	bool all)
{
	const size_t nonz = sh.nonzeroOffsets();
	const int matei = (read.mate <= 1 ? 0 : 1);
	satups_.clear();
	gws_.clear();
	rands_.clear();
	rands2_.clear();
	satpos_.clear();
	satpos2_.clear();
	size_t nrange = 0, nelt = 0, nsmall = 0, nsmall_elts = 0;
	bool keepWhole = false;
	EList<SATupleAndPos, 16> &satpos = keepWhole ? satpos_ : satpos2_;
	for (size_t i = 0; i < nonz; i++)
	{
		bool fw = true;
		uint32_t offidx = 0, rdoff = 0, seedlen = 0;
		QVal qv = sh.hitsByRank(i, offidx, rdoff, fw, seedlen);
		assert(qv.valid());
		assert(!qv.empty());
		assert(qv.repOk(ca.current()));
		ca.queryQval(qv, satups_, nrange, nelt);
		for (size_t j = 0; j < satups_.size(); j++)
		{
			const size_t sz = satups_[j].size();
			if (seedmms == 0)
			{
				EList<ExtendRange> &range =
					fw ? seedExRangeFw_[matei] : seedExRangeRc_[matei];
				bool skip = false;
				for (size_t k = 0; k < range.size(); k++)
				{
					size_t p5 = range[k].off;
					size_t len = range[k].len;
					if (p5 <= rdoff && p5 + len >= (rdoff + seedlen))
					{
						if (sz <= range[k].sz)
						{
							skip = true; 
							break;
						}
					}
				}
				if (skip)
				{
					assert_gt(nrange, 0);
					nrange--;
					assert_geq(nelt, sz);
					nelt -= sz;
					continue;
				}
			}
			satpos.expand();
			satpos.back().sat = satups_[j];
			satpos.back().origSz = sz;
			satpos.back().pos.init(fw, offidx, rdoff, seedlen);
			if (sz <= nsm)
			{
				nsmall++;
				nsmall_elts += sz;
			}
			satpos.back().nlex = satpos.back().nrex = 0;
#ifndef NDEBUG
			tmp_rdseq_.clear();
			uint64_t key = satpos.back().sat.key.seq;
			for (size_t k = 0; k < seedlen; k++)
			{
				int c = (int)(key & 3);
				tmp_rdseq_.append(c);
				key >>= 2;
			}
			tmp_rdseq_.reverse();
#endif
			size_t nlex = 0, nrex = 0;
			if (doExtend)
			{
				extend(
					read,
					ebwtFw,
					ebwtBw,
					satpos.back().sat.topf,
					(TIndexOffU)(satpos.back().sat.topf + sz),
					satpos.back().sat.topb,
					(TIndexOffU)(satpos.back().sat.topb + sz),
					fw,
					rdoff,
					seedlen,
					prm,
					nlex,
					nrex);
			}
			satpos.back().nlex = nlex;
			satpos.back().nrex = nrex;
			if (seedmms == 0 && (nlex > 0 || nrex > 0))
			{
				assert_geq(rdoff, (fw ? nlex : nrex));
				size_t p5 = rdoff - (fw ? nlex : nrex);
				EList<ExtendRange> &range =
					fw ? seedExRangeFw_[matei] : seedExRangeRc_[matei];
				range.expand();
				range.back().off = p5;
				range.back().len = seedlen + nlex + nrex;
				range.back().sz = sz;
			}
		}
		satups_.clear();
	}
	assert_leq(nsmall, nrange);
	nelt_out = nelt;
	assert_eq(nrange, satpos.size());
	satpos.sort();
	if (keepWhole)
	{
		gws_.ensure(nrange);
		rands_.ensure(nrange);
		for (size_t i = 0; i < nrange; i++)
		{
			gws_.expand();
			SARangeWithOffs<TSlice> sa;
			sa.topf = satpos_.back().sat.topf;
			sa.len = satpos_.back().sat.key.len;
			sa.offs = satpos_.back().sat.offs;
			gws_.back().init(
				ebwtFw,
				ref,
				sa,
				rnd,
				wlm);
			assert(gws_.back().initialized());
			rands_.expand();
			rands_.back().init(satpos_[i].sat.size(), all);
		}
		return;
	}
	satpos_.ensure(min(maxelt, nelt));
	gws_.ensure(min(maxelt, nelt));
	rands_.ensure(min(maxelt, nelt));
	rands2_.ensure(min(maxelt, nelt));
	size_t nlarge_elts = nelt - nsmall_elts;
	if (maxelt < nelt)
	{
		size_t diff = nelt - maxelt;
		if (diff >= nlarge_elts)
		{
			nlarge_elts = 0;
		}
		else
		{
			nlarge_elts -= diff;
		}
	}
	size_t nelt_added = 0;
	for (size_t j = 0; j < nsmall && nelt_added < maxelt; j++)
	{
		satpos_.expand();
		satpos_.back() = satpos2_[j];
		gws_.expand();
		SARangeWithOffs<TSlice> sa;
		sa.topf = satpos_.back().sat.topf;
		sa.len = satpos_.back().sat.key.len;
		sa.offs = satpos_.back().sat.offs;
		gws_.back().init(
			ebwtFw,
			ref,
			sa,
			rnd,
			wlm);
		assert(gws_.back().initialized());
		rands_.expand();
		rands_.back().init(satpos_.back().sat.size(), all);
		nelt_added += satpos_.back().sat.size();
#ifndef NDEBUG
		for (size_t k = 0; k < satpos_.size() - 1; k++)
		{
			assert(!(satpos_[k] == satpos_.back()));
		}
#endif
	}
	if (nelt_added >= maxelt || nsmall == satpos2_.size())
	{
		nelt_out = nelt_added;
		return;
	}
	rowsamp_.init(satpos2_, nsmall, satpos2_.size(), lensq, szsq);
	rands2_.resize(satpos2_.size());
	for (size_t j = 0; j < satpos2_.size(); j++)
	{
		rands2_[j].reset();
	}
	while (nelt_added < maxelt && nelt_added < nelt)
	{
		size_t ri = rowsamp_.next(rnd) + nsmall;
		assert_geq(ri, nsmall);
		assert_lt(ri, satpos2_.size());
		if (!rands2_[ri].inited())
		{
			rands2_[ri].init(satpos2_[ri].sat.size(), all);
			assert(!rands2_[ri].done());
		}
		assert(!rands2_[ri].done());
		size_t r = rands2_[ri].next(rnd);
		if (rands2_[ri].done())
		{
			rowsamp_.finishedRange(ri - nsmall);
		}
		SATuple sat;
		TSlice o;
		o.init(satpos2_[ri].sat.offs, r, r + 1);
		sat.init(satpos2_[ri].sat.key, (TIndexOffU)(satpos2_[ri].sat.topf + r), OFF_MASK, o);
		satpos_.expand();
		satpos_.back().sat = sat;
		satpos_.back().origSz = satpos2_[ri].origSz;
		satpos_.back().pos = satpos2_[ri].pos;
		gws_.expand();
		SARangeWithOffs<TSlice> sa;
		sa.topf = sat.topf;
		sa.len = sat.key.len;
		sa.offs = sat.offs;
		gws_.back().init(
			ebwtFw,
			ref,
			sa,
			rnd,
			wlm);
		assert(gws_.back().initialized());
		rands_.expand();
		rands_.back().init(1, all);
		nelt_added++;
	}
	nelt_out = nelt_added;
	return;
}
enum
{
	FOUND_NONE = 0,
	FOUND_EE,
	FOUND_UNGAPPED,
};
int SwDriver::extendSeeds(
	Read &rd,
	bool mate1,
	SeedResults &sh,
	const Ebwt &ebwtFw,
	const Ebwt *ebwtBw,
	const BitPairReference &ref,
	SwAligner &swa,
	const Scoring &sc,
	int seedmms,
	int seedlen,
	int seedival,
	TAlScore &minsc,
	int nceil,
	size_t maxhalf,
	bool doUngapped,
	size_t maxIters,
	size_t maxUg,
	size_t maxDp,
	size_t maxUgStreak,
	size_t maxDpStreak,
	bool doExtend,
	bool enable8,
	size_t cminlen,
	size_t cpow2,
	bool doTri,
	int tighten,
	AlignmentCacheIface &ca,
	RandomSource &rnd,
	WalkMetrics &wlm,
	SwMetrics &swmSeed,
	PerReadMetrics &prm,
	AlnSinkWrap *msink,
	bool reportImmediately,
	bool &exhaustive)
{
#ifdef TIME_STATS
	auto start_extendSeeds = std::chrono::system_clock::now();
#endif
	bool all = msink->allHits();
	assert(!reportImmediately || msink != NULL);
	assert(!reportImmediately || !msink->maxed());
	assert_geq(nceil, 0);
	assert_leq((size_t)nceil, rd.length());
	const size_t rdlen = rd.length();
	TAlScore perfectScore = sc.perfectScore(rdlen);
	DynProgFramer dpframe(!gReportOverhangs);
	swa.reset();
	const size_t nsm = 5;
	const size_t nonz = sh.nonzeroOffsets();
	size_t eeHits = sh.numE2eHits();
	bool eeMode = eeHits > 0;
	bool firstEe = true;
	bool firstExtend = true;
	prm.nEeFail = 0;
	prm.nUgFail = 0;
	prm.nDpFail = 0;
	size_t nelt = 0, neltLeft = 0;
	size_t rows = rdlen;
	size_t eltsDone = 0;
	while (true)
	{
		if (eeMode)
		{
			if (firstEe)
			{
				firstEe = false;
#ifdef TIME_STATS
				auto start_extendSeeds_eeSaTups = std::chrono::system_clock::now();
#endif
				eeMode = eeSaTups(
					rd,
					sh,
					ebwtFw,
					ref,
					rnd,
					wlm,
					swmSeed,
					nelt,
					maxIters,
					all);
#ifdef TIME_STATS
				auto end_extendSeeds_eeSaTups = std::chrono::system_clock::now();
				auto elapsed_extendSeeds_eeSaTups = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds_eeSaTups - start_extendSeeds_eeSaTups);
				time_extendSeeds_eeSaTups += elapsed_extendSeeds_eeSaTups.count();
#endif
				assert_eq(gws_.size(), rands_.size());
				assert_eq(gws_.size(), satpos_.size());
			}
			else
			{
				eeMode = false;
			}
		}
		if (!eeMode)
		{
			if (nonz == 0)
			{
				return EXTEND_EXHAUSTED_CANDIDATES;
			}
			if (minsc == perfectScore)
			{
				return EXTEND_PERFECT_SCORE;
			}
			if (firstExtend)
			{
				nelt = 0;
#ifdef TIME_STATS
				auto start_extendSeeds_prioritizeSATups = std::chrono::system_clock::now();
#endif
				prioritizeSATups(
					rd,
					sh,
					ebwtFw,
					ebwtBw,
					ref,
					seedmms,
					maxIters,
					doExtend,
					true,
					true,
					nsm,
					ca,
					rnd,
					wlm,
					prm,
					nelt,
					all);
#ifdef TIME_STATS
				auto end_extendSeeds_prioritizeSATups = std::chrono::system_clock::now();
				auto elapsed_extendSeeds_prioritizeSATups = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds_prioritizeSATups - start_extendSeeds_prioritizeSATups);
				time_extendSeeds_prioritizeSATups += elapsed_extendSeeds_prioritizeSATups.count();
#endif
				assert_eq(gws_.size(), rands_.size());
				assert_eq(gws_.size(), satpos_.size());
				neltLeft = nelt;
				firstExtend = false;
			}
			if (neltLeft == 0)
			{
				break;
			}
		}
		for (size_t i = 0; i < gws_.size(); i++)
		{
			if (eeMode && eehits_[i].score < minsc)
			{
				return EXTEND_PERFECT_SCORE;
			}
			bool is_small = satpos_[i].sat.size() < nsm;
			bool fw = satpos_[i].pos.fw;
			uint32_t rdoff = satpos_[i].pos.rdoff;
			uint32_t seedhitlen = satpos_[i].pos.seedlen;
			if (!fw)
			{
				rdoff = (uint32_t)(rdlen - rdoff - seedhitlen);
			}
			bool first = true;
			size_t riter = 0;
			while (!rands_[i].done() && (first || is_small || eeMode))
			{
				assert(!gws_[i].done());
				riter++;
				if (minsc == perfectScore)
				{
					if (!eeMode || eehits_[i].score < perfectScore)
					{
						return EXTEND_PERFECT_SCORE;
					}
				}
				else if (eeMode && eehits_[i].score < minsc)
				{
					break;
				}
				if (prm.nExDps >= maxDp || prm.nMateDps >= maxDp)
				{
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if (prm.nExUgs >= maxUg || prm.nMateUgs >= maxUg)
				{
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if (prm.nExIters >= maxIters)
				{
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				prm.nExIters++;
				first = false;
				WalkResult wr;
				size_t elt = rands_[i].next(rnd);
				SARangeWithOffs<TSlice> sa;
				sa.topf = satpos_[i].sat.topf;
				sa.len = satpos_[i].sat.key.len;
				sa.offs = satpos_[i].sat.offs;
#ifdef TIME_STATS
				auto start_extendSeeds_advanceElement = std::chrono::system_clock::now();
#endif
				gws_[i].advanceElement((TIndexOffU)elt, ebwtFw, ref, sa, gwstate_, wr, wlm, prm);
#ifdef TIME_STATS
				auto end_extendSeeds_advanceElement = std::chrono::system_clock::now();
				auto elapsed_extendSeeds_advanceElement = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds_advanceElement - start_extendSeeds_advanceElement);
				time_extendSeeds_advanceElement += elapsed_extendSeeds_advanceElement.count();
#endif
				eltsDone++;
				if (!eeMode)
				{
					assert_gt(neltLeft, 0);
					neltLeft--;
				}
				assert_neq(OFF_MASK, wr.toff);
				TIndexOffU tidx = 0, toff = 0, tlen = 0;
				bool straddled = false;
#ifdef TIME_STATS
				auto start_extendSeeds_joinedToTextOff = std::chrono::system_clock::now();
#endif
				ebwtFw.joinedToTextOff(
					wr.elt.len,
					wr.toff,
					tidx,
					toff,
					tlen,
					eeMode,
					straddled);
#ifdef TIME_STATS
				auto end_extendSeeds_joinedToTextOff = std::chrono::system_clock::now();
				auto elapsed_extendSeeds_joinedToTextOff = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds_joinedToTextOff - start_extendSeeds_joinedToTextOff);
				time_extendSeeds_joinedToTextOff += elapsed_extendSeeds_joinedToTextOff.count();
#endif
				if (tidx == OFF_MASK)
				{
					continue;
				}
#ifndef NDEBUG
				if (!eeMode && !straddled)
				{
					uint64_t key = satpos_[i].sat.key.seq;
					for (size_t k = 0; k < wr.elt.len; k++)
					{
						int c = ref.getBase(tidx, toff + wr.elt.len - k - 1);
						assert_leq(c, 3);
						int ck = (int)(key & 3);
						key >>= 2;
						assert_eq(c, ck);
					}
				}
#endif
				int64_t refoff = (int64_t)toff - rdoff;
				Coord refcoord(tidx, refoff, fw);
				if (seenDiags1_.locusPresent(refcoord))
				{
					prm.nRedundants++;
					swmSeed.rshit++;
#ifdef TIME_STATS
					redundants_Nums++;
#endif
					continue;
				}
				int readGaps = 0, refGaps = 0;
				bool ungapped = false;
				if (!eeMode)
				{
					readGaps = sc.maxReadGaps(minsc, rdlen);
					refGaps = sc.maxRefGaps(minsc, rdlen);
					ungapped = (readGaps == 0 && refGaps == 0);
				}
				int state = FOUND_NONE;
				bool found = false;
				if (eeMode)
				{
#ifdef TIME_STATS
					auto start_extendSeeds_resEe = std::chrono::system_clock::now();
#endif
					resEe_.reset();
					resEe_.alres.reset();
					const EEHit &h = eehits_[i];
					assert_leq(h.score, perfectScore);
					resEe_.alres.setScore(AlnScore(h.score,
												   (int)(rdlen - h.mms()),
												   h.mms(), h.ns(), 0));
					resEe_.alres.setShape(
						refcoord.ref(),
						refcoord.off(),
						tlen,
						fw,
						rdlen,
						true,
						0,
						0,
						true,
						0,
						0);
					resEe_.alres.setRefNs(h.refns());
					if (h.mms() > 0)
					{
						assert_eq(1, h.mms());
						assert_lt(h.e1.pos, rd.length());
						resEe_.alres.ned().push_back(h.e1);
					}
					assert(resEe_.repOk(rd));
					state = FOUND_EE;
					found = true;
					Interval refival(refcoord, 1);
					seenDiags1_.add(refival);
#ifdef TIME_STATS
					auto end_extendSeeds_resEe = std::chrono::system_clock::now();
					auto elapsed_extendSeeds_resEe = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds_resEe - start_extendSeeds_resEe);
					time_extendSeeds_resEe += elapsed_extendSeeds_resEe.count();
#endif
				}
				else if (doUngapped && ungapped)
				{
					resUngap_.reset();
#ifdef TIME_STATS
					auto start_extendSeeds_resUngap = std::chrono::system_clock::now();
#endif
					int al = swa.ungappedAlign(
						fw ? rd.patFw : rd.patRc,
						fw ? rd.qual : rd.qualRev,
						refcoord,
						ref,
						tlen,
						sc,
						gReportOverhangs,
						minsc,
						resUngap_);
					Interval refival(refcoord, 1);
					seenDiags1_.add(refival);
					prm.nExUgs++;
#ifdef TIME_STATS
					auto end_extendSeeds_resUngap = std::chrono::system_clock::now();
					auto elapsed_extendSeeds_resUngap = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds_resUngap - start_extendSeeds_resUngap);
					time_extendSeeds_resUngap += elapsed_extendSeeds_resUngap.count();
#endif
					if (al == 0)
					{
						prm.nExUgFails++;
						prm.nUgFail++;
						if (prm.nUgFail >= maxUgStreak)
						{
							return EXTEND_EXCEEDED_SOFT_LIMIT;
						}
						swmSeed.ungapfail++;
						continue;
					}
					else if (al == -1)
					{
						prm.nExUgFails++;
						prm.nUgFail++;
						if (prm.nUgFail >= maxUgStreak)
						{
							return EXTEND_EXCEEDED_SOFT_LIMIT;
						}
						swmSeed.ungapnodec++;
					}
					else
					{
						prm.nExUgSuccs++;
						prm.nUgLastSucc = prm.nExUgs - 1;
						if (prm.nUgFail > prm.nUgFailStreak)
						{
							prm.nUgFailStreak = prm.nUgFail;
						}
						prm.nUgFail = 0;
						found = true;
						state = FOUND_UNGAPPED;
						swmSeed.ungapsucc++;
					}
				}
				int64_t pastedRefoff = (int64_t)wr.toff - rdoff;
				DPRect rect;
				if (state == FOUND_NONE)
				{
					found = dpframe.frameSeedExtensionRect(
						refoff,
						rows,
						tlen,
						readGaps,
						refGaps,
						(size_t)nceil,
						maxhalf,
						rect);
					assert(rect.repOk());
					seenDiags1_.add(Interval(refcoord, 1));
					if (!found)
					{
						continue;
					}
				}
				int64_t leftShift = refoff - rect.refl;
				size_t nwindow = 0;
				if ((int64_t)toff >= rect.refl)
				{
					nwindow = (size_t)(toff - rect.refl);
				}
				pastedRefoff -= leftShift;
				size_t nsInLeftShift = 0;
				if (state == FOUND_NONE)
				{
#ifdef TIME_STATS
					auto start_extendSeeds_init = std::chrono::system_clock::now();
#endif
					if (!swa.initedRead())
					{
						swa.initRead(
							rd.patFw,
							rd.patRc,
							rd.qual,
							rd.qualRev,
							0,
							rdlen,
							sc);
					}
					swa.initRef(
						fw,
						tidx,
						rect,
						ref,
						tlen,
						sc,
						minsc,
						enable8,
						cminlen,
						cpow2,
						doTri,
						true,
						nwindow,
						nsInLeftShift);
#ifdef TIME_STATS
					auto end_extendSeeds_init = std::chrono::system_clock::now();
					auto elapsed_extendSeeds_init = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds_init - start_extendSeeds_init);
					time_extendSeeds_init += elapsed_extendSeeds_init.count();
#endif
					Interval refival(tidx, 0, fw, 0);
					rect.initIval(refival);
					seenDiags1_.add(refival);
					TAlScore bestCell = std::numeric_limits<TAlScore>::min();
#ifdef TIME_STATS
					auto start_extendSeeds_align = std::chrono::system_clock::now();
#endif
					found = swa.align(bestCell, resGap_, refcoord);
#ifdef TIME_STATS
					auto end_extendSeeds_align = std::chrono::system_clock::now();
					auto elapsed_extendSeeds_align = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds_align - start_extendSeeds_align);
					time_extendSeeds_align += elapsed_extendSeeds_align.count();
#endif
					swmSeed.tallyGappedDp(readGaps, refGaps);
					prm.nExDps++;
					if (!found)
					{
						prm.nExDpFails++;
						prm.nDpFail++;
						if (prm.nDpFail >= maxDpStreak)
						{
							return EXTEND_EXCEEDED_SOFT_LIMIT;
						}
						if (bestCell > std::numeric_limits<TAlScore>::min() && bestCell > prm.bestLtMinscMate1)
						{
							prm.bestLtMinscMate1 = bestCell;
						}
						continue;
					}
					else
					{
						prm.nExDpSuccs++;
						prm.nDpLastSucc = prm.nExDps - 1;
						if (prm.nDpFail > prm.nDpFailStreak)
						{
							prm.nDpFailStreak = prm.nDpFail;
						}
						prm.nDpFail = 0;
					}
				}
				bool firstInner = true;
				while (true)
				{
					SwResult *res = NULL;
					if (state == FOUND_EE)
					{
						if (!firstInner)
						{
							break;
						}
						res = &resEe_;
					}
					else if (state == FOUND_UNGAPPED)
					{
						if (!firstInner)
						{
							break;
						}
						res = &resUngap_;
					}
					else
					{
						if (!firstInner)
						{
							break;
						}
						if (swa.done())
						{
							break;
						}
						found = !resGap_.empty();
						if (!found)
						{
							break;
						}
						res = &resGap_;
						int lenSTR=5; 
						string sb=res->alres.mycigar;
						int nI=count(sb.begin(),sb.end(),'I'); 
						int nD=count(sb.begin(),sb.end(),'D');
						int AnI=0, AnD=0; 
						int posD[nD];
						int posD2[nD];
						int posI[nI];
						int posI2[nI];
						int p=0, i=0;
						if(nD>0){
							while((p=sb.find("D",p))!=string::npos){ 
								posD[i]=p;
								posD2[i]=p; 
								p++;
								i++;
							}
							string r="";
							if(res->alres.fw()){
								r=string(rd.patFw.toZBuf());
							}else{
								r=string(rd.patRc.toZBuf());
							}
							for(int i=0; i<nD; i++){
								if(posD[i]>rdoff && posD[i]<rdoff+seedlen){ 
									if(i>0 && posD[i]==posD[i-1]+1){ 
										posD2[i]=posD2[i-1];
										continue;
									}
									int slen=0;
									int lmv=1, plmv=0;
									if(rdoff>0){
										string sl=r.substr(0, posD2[i]); 
										slen=sl.length();
										while(true){
											int k=0;
											for(int i=0; i<lenSTR; i++){
												lmv++; 
												if(slen-lmv==0){ 
													goto end1;
												}
												if(isStr(sl.substr(slen-lmv))){ 
													k++;
													plmv=lmv;
													break; 
												}
											}
											if(k==0){ 
												break;
											}
										}
										end1:;
									}
									int rmv=1, prmv=0;
									if(rdoff<(rdlen-seedlen)){
										string sr=r.substr(posD2[i], rdlen+nD-posD2[i]-1); 
										slen=sr.length();
										while(true){
											int k=0;
											for(int i=0; i<lenSTR; i++){
												rmv++; 
												if(posD2[i]+rmv==rdlen+nD){ 
													goto end2;
												}
												if(isStr(sr.substr(0, rmv))){ 
													k++;
													prmv=rmv;
													break; 
												}
											}
											if(k==0){ 
												break;
											}
										}
										end2:;
									}
									if(plmv>=prmv){ 
										posD2[i]=posD[i]-plmv;
									}else{
										posD2[i]=posD[i]+prmv;
									}
								}                            
							}                                
							for(int i=0; i<nD; i++){ 
								if(posD2[i]<=rdoff){
									AnD++;
								}
							}
						}
						p=0; i=0;
						if(nI>0){
							while((p=sb.find("I",p))!=string::npos){ 
								posI[i]=p;
								posI2[i]=p;
								p++;
								i++;
							}
							string r="";
							if(res->alres.fw()){
								r=string(rd.patFw.toZBuf());
							}else{
								r=string(rd.patRc.toZBuf());
							}
							for(int i=0; i<nI; i++){
								if(posI[i]>rdoff && posI[i]<rdoff+seedlen){ 
									if(i>0 && posI[i]==posI[i-1]+1){ 
										posI2[i]=posI2[i-1];
										continue;
									}
									int slen=0;
									int lmv=1, plmv=0;
									if(rdoff>0){
										string sl=r.substr(0, posI2[i]); 
										slen=sl.length();
										while(true){
											int k=0;
											for(int i=0; i<lenSTR; i++){
												lmv++; 
												if(slen-lmv==0){ 
													goto end3;
												}
												if(isStr(sl.substr(slen-lmv))){ 
													k++;
													plmv=lmv;
													break; 
												}
											}
											if(k==0){ 
												break;
											}
										}
										end3:;
									}
									int rmv=1, prmv=0;
									if(rdoff<(rdlen-seedlen)){
										string sr=r.substr(posI2[i]+i, rdlen-posI2[i]-1); 
										slen=sr.length();
										while(true){
											int k=0;
											for(int i=0; i<lenSTR; i++){
												rmv++; 
												if(posI2[i]+rmv==rdlen){ 
													goto end4;
												}
												if(isStr(sr.substr(0, rmv))){ 
													k++;
													prmv=rmv;
													break; 
												}
											}
											if(k==0){ 
												break;
											}
										}
										end4:;
									}
									if(plmv>=prmv){
										posI2[i]=posI[i]-plmv;
									}else{
										posI2[i]=posI[i]+prmv;
									}
								}                         
							}                                
							for(int i=0; i<nI; i++){ 
								if(posI2[i]<=rdoff){
									AnI++;
								}
							}
						}
						res->alres.refcoord_.off_=res->alres.refcoord_.off_+AnI-AnD;
						int rlen=rdlen;
						string cigar =res->alres.mycigar;
						res->alres.mycigar = swa.cigarformat(cigar, rlen);
					}
					assert(res != NULL);
					firstInner = false;
					assert(res->alres.matchesRef(
						rd,
						ref,
						tmp_rf_,
						tmp_rdseq_,
						tmp_qseq_,
						raw_refbuf_,
						raw_destU32_,
						raw_matches_));
					Interval refival(tidx, 0, fw, tlen);
					assert_gt(res->alres.refExtent(), 0);
					if (gReportOverhangs &&
						!refival.containsIgnoreOrient(res->alres.refival()))
					{
						res->alres.clipOutside(true, 0, tlen);
						if (res->alres.refExtent() == 0)
						{
							continue;
						}
					}
					assert(gReportOverhangs ||
						   refival.containsIgnoreOrient(res->alres.refival()));
					if (!refival.overlapsIgnoreOrient(res->alres.refival()))
					{
						continue;
					}
					if (redAnchor_.overlap(res->alres))
					{
						continue;
					}
					redAnchor_.add(res->alres);
					res->alres.setParams(
						seedmms,
						seedlen,
						seedival,
						minsc);
					if (reportImmediately)
					{
						assert(msink != NULL);
						assert(res->repOk());
						assert(res->alres.matchesRef(
							rd,
							ref,
							tmp_rf_,
							tmp_rdseq_,
							tmp_qseq_,
							raw_refbuf_,
							raw_destU32_,
							raw_matches_));
						assert(!msink->maxed());
						if (msink->report(
								0,
								mate1 ? &res->alres : NULL,
								mate1 ? NULL : &res->alres))
						{
							return EXTEND_POLICY_FULFILLED;
						}
						if (tighten > 0 &&
							msink->Mmode() &&
							msink->hasSecondBestUnp1())
						{
							if (tighten == 1)
							{
								if (msink->bestUnp1() >= minsc)
								{
									minsc = msink->bestUnp1();
									if (minsc < perfectScore &&
										msink->bestUnp1() == msink->secondBestUnp1())
									{
										minsc++;
									}
								}
							}
							else if (tighten == 2)
							{
								if (msink->secondBestUnp1() >= minsc)
								{
									minsc = msink->secondBestUnp1();
									if (minsc < perfectScore)
									{
										minsc++;
									}
								}
							}
							else
							{
								TAlScore diff = msink->bestUnp1() - msink->secondBestUnp1();
								TAlScore bot = msink->secondBestUnp1() + ((diff * 3) / 4);
								if (bot >= minsc)
								{
									minsc = bot;
									if (minsc < perfectScore)
									{
										minsc++;
									}
								}
							}
							assert_leq(minsc, perfectScore);
						}
					}
				}
			}
		}
	}
#ifdef TIME_STATS
	auto end_extendSeeds = std::chrono::system_clock::now();
	auto elapsed_extendSeeds = std::chrono::duration_cast<std::chrono::microseconds>(end_extendSeeds - start_extendSeeds);
	time_extendSeeds += elapsed_extendSeeds.count();
#endif
	return EXTEND_EXHAUSTED_CANDIDATES;
}
int SwDriver::extendSeedsPaired(
	Read &rd,
	Read &ord,
	bool anchor1,
	bool oppFilt,
	SeedResults &sh,
	const Ebwt &ebwtFw,
	const Ebwt *ebwtBw,
	const BitPairReference &ref,
	SwAligner &swa,
	SwAligner &oswa,
	const Scoring &sc,
	const PairedEndPolicy &pepol,
	int seedmms,
	int seedlen,
	int seedival,
	TAlScore &minsc,
	TAlScore &ominsc,
	int nceil,
	int onceil,
	bool nofw,
	bool norc,
	size_t maxhalf,
	bool doUngapped,
	size_t maxIters,
	size_t maxUg,
	size_t maxDp,
	size_t maxEeStreak,
	size_t maxUgStreak,
	size_t maxDpStreak,
	size_t maxMateStreak,
	bool doExtend,
	bool enable8,
	size_t cminlen,
	size_t cpow2,
	bool doTri,
	int tighten,
	AlignmentCacheIface &ca,
	RandomSource &rnd,
	WalkMetrics &wlm,
	SwMetrics &swmSeed,
	SwMetrics &swmMate,
	PerReadMetrics &prm,
	AlnSinkWrap *msink,
	bool swMateImmediately,
	bool reportImmediately,
	bool discord,
	bool mixed,
	bool &exhaustive)
{
	bool all = msink->allHits();
	assert(!reportImmediately || msink != NULL);
	assert(!reportImmediately || !msink->maxed());
	assert(!msink->state().doneWithMate(anchor1));
	assert_geq(nceil, 0);
	assert_geq(onceil, 0);
	assert_leq((size_t)nceil, rd.length());
	assert_leq((size_t)onceil, ord.length());
	const size_t rdlen = rd.length();
	const size_t ordlen = ord.length();
	const TAlScore perfectScore = sc.perfectScore(rdlen);
	const TAlScore operfectScore = sc.perfectScore(ordlen);
	assert_leq(minsc, perfectScore);
	assert(oppFilt || ominsc <= operfectScore);
	TAlScore bestPairScore = perfectScore + operfectScore;
	if (tighten > 0 && msink->Mmode() && msink->hasSecondBestPair())
	{
		TAlScore ps;
		if (tighten == 1)
		{
			ps = msink->bestPair();
		}
		else if (tighten == 2)
		{
			ps = msink->secondBestPair();
		}
		else
		{
			TAlScore diff = msink->bestPair() - msink->secondBestPair();
			ps = msink->secondBestPair() + (diff * 3) / 4;
		}
		if (tighten == 1 && ps < bestPairScore &&
			msink->bestPair() == msink->secondBestPair())
		{
			ps++;
		}
		if (tighten >= 2 && ps < bestPairScore)
		{
			ps++;
		}
		TAlScore nc = ps - operfectScore;
		if (nc > minsc)
		{
			minsc = nc;
		}
		assert_leq(minsc, perfectScore);
	}
	DynProgFramer dpframe(!gReportOverhangs);
	swa.reset();
	oswa.reset();
	const size_t nsm = 5;
	const size_t nonz = sh.nonzeroOffsets();
	size_t eeHits = sh.numE2eHits();
	bool eeMode = eeHits > 0;
	bool firstEe = true;
	bool firstExtend = true;
	prm.nEeFail = 0;
	prm.nUgFail = 0;
	prm.nDpFail = 0;
	size_t nelt = 0, neltLeft = 0;
	const size_t rows = rdlen;
	const size_t orows = ordlen;
	size_t eltsDone = 0;
	while (true)
	{
		if (eeMode)
		{
			if (firstEe)
			{
				firstEe = false;
				eeMode = eeSaTups(
					rd,
					sh,
					ebwtFw,
					ref,
					rnd,
					wlm,
					swmSeed,
					nelt,
					maxIters,
					all);
				assert_eq(gws_.size(), rands_.size());
				assert_eq(gws_.size(), satpos_.size());
				neltLeft = nelt;
				mateStreaks_.resize(gws_.size());
				mateStreaks_.fill(0);
			}
			else
			{
				eeMode = false;
			}
		}
		if (!eeMode)
		{
			if (nonz == 0)
			{
				return EXTEND_EXHAUSTED_CANDIDATES;
			}
			if (msink->Mmode() && minsc == perfectScore)
			{
				return EXTEND_PERFECT_SCORE;
			}
			if (firstExtend)
			{
				nelt = 0;
				prioritizeSATups(
					rd,
					sh,
					ebwtFw,
					ebwtBw,
					ref,
					seedmms,
					maxIters,
					doExtend,
					true,
					true,
					nsm,
					ca,
					rnd,
					wlm,
					prm,
					nelt,
					all);
				assert_eq(gws_.size(), rands_.size());
				assert_eq(gws_.size(), satpos_.size());
				neltLeft = nelt;
				firstExtend = false;
				mateStreaks_.resize(gws_.size());
				mateStreaks_.fill(0);
			}
			if (neltLeft == 0)
			{
				break;
			}
		}
		for (size_t i = 0; i < gws_.size(); i++)
		{
			if (eeMode && eehits_[i].score < minsc)
			{
				return EXTEND_PERFECT_SCORE;
			}
			bool is_small = satpos_[i].sat.size() < nsm;
			bool fw = satpos_[i].pos.fw;
			uint32_t rdoff = satpos_[i].pos.rdoff;
			uint32_t seedhitlen = satpos_[i].pos.seedlen;
			if (!fw)
			{
				rdoff = (uint32_t)(rdlen - rdoff - seedhitlen);
			}
			bool first = true;
			while (!rands_[i].done() && (first || is_small || eeMode))
			{
				if (minsc == perfectScore)
				{
					if (!eeMode || eehits_[i].score < perfectScore)
					{
						return EXTEND_PERFECT_SCORE;
					}
				}
				else if (eeMode && eehits_[i].score < minsc)
				{
					break;
				}
				if (prm.nExDps >= maxDp || prm.nMateDps >= maxDp)
				{
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if (prm.nExUgs >= maxUg || prm.nMateUgs >= maxUg)
				{
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if (prm.nExIters >= maxIters)
				{
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if (eeMode && prm.nEeFail >= maxEeStreak)
				{
					return EXTEND_EXCEEDED_SOFT_LIMIT;
				}
				if (!eeMode && prm.nDpFail >= maxDpStreak)
				{
					return EXTEND_EXCEEDED_SOFT_LIMIT;
				}
				if (!eeMode && prm.nUgFail >= maxUgStreak)
				{
					return EXTEND_EXCEEDED_SOFT_LIMIT;
				}
				if (mateStreaks_[i] >= maxMateStreak)
				{
					rands_[i].setDone();
					assert(rands_[i].done());
					break;
				}
				prm.nExIters++;
				first = false;
				assert(!gws_[i].done());
				WalkResult wr;
				size_t elt = rands_[i].next(rnd);
				SARangeWithOffs<TSlice> sa;
				sa.topf = satpos_[i].sat.topf;
				sa.len = satpos_[i].sat.key.len;
				sa.offs = satpos_[i].sat.offs;
				gws_[i].advanceElement((TIndexOffU)elt, ebwtFw, ref, sa, gwstate_, wr, wlm, prm);
				eltsDone++;
				assert_gt(neltLeft, 0);
				neltLeft--;
				assert_neq(OFF_MASK, wr.toff);
				TIndexOffU tidx = 0, toff = 0, tlen = 0;
				bool straddled = false;
				ebwtFw.joinedToTextOff(
					wr.elt.len,
					wr.toff,
					tidx,
					toff,
					tlen,
					eeMode,
					straddled);
				if (tidx == OFF_MASK)
				{
					continue;
				}
#ifndef NDEBUG
				if (!eeMode && !straddled)
				{
					uint64_t key = satpos_[i].sat.key.seq;
					for (size_t k = 0; k < wr.elt.len; k++)
					{
						int c = ref.getBase(tidx, toff + wr.elt.len - k - 1);
						assert_leq(c, 3);
						int ck = (int)(key & 3);
						key >>= 2;
						assert_eq(c, ck);
					}
				}
#endif
				int64_t refoff = (int64_t)toff - rdoff;
				cerr << "refoff:" << refoff << "  toff:" << toff << " rdoff:" << rdoff << endl;
				EIvalMergeListBinned &seenDiags = anchor1 ? seenDiags1_ : seenDiags2_;
				Coord refcoord(tidx, refoff, fw);
				if (seenDiags.locusPresent(refcoord))
				{
					prm.nRedundants++;
					swmSeed.rshit++;
					continue;
				}
				int readGaps = 0, refGaps = 0;
				bool ungapped = false;
				if (!eeMode)
				{
					readGaps = sc.maxReadGaps(minsc, rdlen);
					refGaps = sc.maxRefGaps(minsc, rdlen);
					ungapped = (readGaps == 0 && refGaps == 0);
				}
				int state = FOUND_NONE;
				bool found = false;
				if (eeMode)
				{
					resEe_.reset();
					resEe_.alres.reset();
					const EEHit &h = eehits_[i];
					assert_leq(h.score, perfectScore);
					resEe_.alres.setScore(AlnScore(
						h.score,
						(int)(rdlen - h.mms()),
						h.mms(), h.ns(), 0));
					resEe_.alres.setShape(
						refcoord.ref(),
						refcoord.off(),
						tlen,
						fw,
						rdlen,
						true,
						0,
						0,
						true,
						0,
						0);
					resEe_.alres.setRefNs(h.refns());
					if (h.mms() > 0)
					{
						assert_eq(1, h.mms());
						assert_lt(h.e1.pos, rd.length());
						resEe_.alres.ned().push_back(h.e1);
					}
					assert(resEe_.repOk(rd));
					state = FOUND_EE;
					found = true;
					Interval refival(refcoord, 1);
					seenDiags.add(refival);
					prm.nExEes++;
					prm.nEeFail++;
					prm.nExEeFails++;
				}
				else if (doUngapped && ungapped)
				{
					resUngap_.reset();
					int al = swa.ungappedAlign(
						fw ? rd.patFw : rd.patRc,
						fw ? rd.qual : rd.qualRev,
						refcoord,
						ref,
						tlen,
						sc,
						gReportOverhangs,
						minsc,
						resUngap_);
					Interval refival(refcoord, 1);
					seenDiags.add(refival);
					prm.nExUgs++;
					prm.nUgFail++;
					prm.nExUgFails++;
					if (al == 0)
					{
						swmSeed.ungapfail++;
						continue;
					}
					else if (al == -1)
					{
						swmSeed.ungapnodec++;
					}
					else
					{
						found = true;
						state = FOUND_UNGAPPED;
						swmSeed.ungapsucc++;
					}
				}
				int64_t pastedRefoff = (int64_t)wr.toff - rdoff;
				DPRect rect;
				if (state == FOUND_NONE)
				{
					found = dpframe.frameSeedExtensionRect(
						refoff,
						rows,
						tlen,
						readGaps,
						refGaps,
						(size_t)nceil,
						maxhalf,
						rect);
					assert(rect.repOk());
					seenDiags.add(Interval(refcoord, 1));
					if (!found)
					{
						continue;
					}
				}
				int64_t leftShift = refoff - rect.refl;
				size_t nwindow = 0;
				if ((int64_t)toff >= rect.refl)
				{
					nwindow = (size_t)(toff - rect.refl);
				}
				pastedRefoff -= leftShift;
				size_t nsInLeftShift = 0;
				if (state == FOUND_NONE)
				{
					if (!swa.initedRead())
					{
						swa.initRead(
							rd.patFw,
							rd.patRc,
							rd.qual,
							rd.qualRev,
							0,
							rdlen,
							sc);
					}
					swa.initRef(
						fw,
						tidx,
						rect,
						ref,
						tlen,
						sc,
						minsc,
						enable8,
						cminlen,
						cpow2,
						doTri,
						true,
						nwindow,
						nsInLeftShift);
					Interval refival(tidx, 0, fw, 0);
					rect.initIval(refival);
					seenDiags.add(refival);
					TAlScore bestCell = std::numeric_limits<TAlScore>::min();
					found = swa.align(bestCell, resGap_, refcoord);
					swmSeed.tallyGappedDp(readGaps, refGaps);
					prm.nExDps++;
					prm.nDpFail++;
					prm.nExDpFails++;
					if (!found)
					{
						TAlScore bestLast = anchor1 ? prm.bestLtMinscMate1 : prm.bestLtMinscMate2;
						if (bestCell > std::numeric_limits<TAlScore>::min() && bestCell > bestLast)
						{
							if (anchor1)
							{
								prm.bestLtMinscMate1 = bestCell;
							}
							else
							{
								prm.bestLtMinscMate2 = bestCell;
							}
						}
						continue;
					}
				}
				bool firstInner = true;
				bool foundConcordant = false;
				while (true)
				{
					assert(found);
					SwResult *res = NULL;
					if (state == FOUND_EE)
					{
						if (!firstInner)
						{
							break;
						}
						res = &resEe_;
						assert(res->repOk(rd));
					}
					else if (state == FOUND_UNGAPPED)
					{
						if (!firstInner)
						{
							break;
						}
						res = &resUngap_;
						assert(res->repOk(rd));
					}
					else
					{
						resGap_.reset();
						assert(resGap_.empty());
						if (swa.done())
						{
							break;
						}
						swa.nextAlignment(resGap_, minsc, rnd);
						found = !resGap_.empty();
						if (!found)
						{
							break;
						}
						res = &resGap_;
						assert(res->repOk(rd));
					}
					assert(res != NULL);
					firstInner = false;
					assert(res->alres.matchesRef(
						rd,
						ref,
						tmp_rf_,
						tmp_rdseq_,
						tmp_qseq_,
						raw_refbuf_,
						raw_destU32_,
						raw_matches_));
					Interval refival(tidx, 0, fw, tlen);
					assert_gt(res->alres.refExtent(), 0);
					if (gReportOverhangs &&
						!refival.containsIgnoreOrient(res->alres.refival()))
					{
						res->alres.clipOutside(true, 0, tlen);
						if (res->alres.refExtent() == 0)
						{
							continue;
						}
					}
					assert(gReportOverhangs ||
						   refival.containsIgnoreOrient(res->alres.refival()));
					if (!refival.overlapsIgnoreOrient(res->alres.refival()))
					{
						continue;
					}
					if (redAnchor_.overlap(res->alres))
					{
						continue;
					}
					redAnchor_.add(res->alres);
					res->alres.setParams(
						seedmms,
						seedlen,
						seedival,
						minsc);
					bool foundMate = false;
					TRefOff off = res->alres.refoff();
					if (msink->state().doneWithMate(!anchor1) &&
						!msink->state().doneWithMate(anchor1))
					{
						swMateImmediately = false;
					}
					if (found && swMateImmediately)
					{
						assert(!msink->state().doneWithMate(!anchor1));
						bool oleft = false, ofw = false;
						int64_t oll = 0, olr = 0, orl = 0, orr = 0;
						assert(!msink->state().done());
						foundMate = !oppFilt;
						TAlScore ominsc_cur = ominsc;
						int oreadGaps = 0, orefGaps = 0;
						if (foundMate)
						{
							ominsc_cur = ominsc;
							if (tighten > 0 && msink->Mmode() && msink->hasSecondBestPair())
							{
								TAlScore ps;
								if (tighten == 1)
								{
									ps = msink->bestPair();
								}
								else if (tighten == 2)
								{
									ps = msink->secondBestPair();
								}
								else
								{
									TAlScore diff = msink->bestPair() - msink->secondBestPair();
									ps = msink->secondBestPair() + (diff * 3) / 4;
								}
								if (tighten == 1 && ps < bestPairScore &&
									msink->bestPair() == msink->secondBestPair())
								{
									ps++;
								}
								if (tighten >= 2 && ps < bestPairScore)
								{
									ps++;
								}
								TAlScore nc = ps - res->alres.score().score();
								if (nc > ominsc_cur)
								{
									ominsc_cur = nc;
									assert_leq(ominsc_cur, operfectScore);
								}
							}
							oreadGaps = sc.maxReadGaps(ominsc_cur, ordlen);
							orefGaps = sc.maxRefGaps(ominsc_cur, ordlen);
							foundMate = pepol.otherMate(
								anchor1,
								fw,
								off,
								orows + oreadGaps,
								tlen,
								anchor1 ? rd.length() : ord.length(),
								anchor1 ? ord.length() : rd.length(),
								oleft,
								oll,
								olr,
								orl,
								orr,
								ofw);
						}
						DPRect orect;
						if (foundMate)
						{
							foundMate = dpframe.frameFindMateRect(
								!oleft,
								oll,
								olr,
								orl,
								orr,
								orows,
								tlen,
								oreadGaps,
								orefGaps,
								(size_t)onceil,
								maxhalf,
								orect);
							assert(!foundMate || orect.refr >= orect.refl);
						}
						if (foundMate)
						{
							oresGap_.reset();
							assert(oresGap_.empty());
							if (!oswa.initedRead())
							{
								oswa.initRead(
									ord.patFw,
									ord.patRc,
									ord.qual,
									ord.qualRev,
									0,
									ordlen,
									sc);
							}
							size_t onsInLeftShift = 0;
							assert_geq(orect.refr, orect.refl);
							oswa.initRef(
								ofw,
								tidx,
								orect,
								ref,
								tlen,
								sc,
								ominsc_cur,
								enable8,
								cminlen,
								cpow2,
								doTri,
								false,
								0,
								onsInLeftShift);
							TAlScore bestCell = std::numeric_limits<TAlScore>::min();
							foundMate = oswa.align(bestCell, oresGap_, refcoord);
							prm.nMateDps++;
							swmMate.tallyGappedDp(oreadGaps, orefGaps);
							if (!foundMate)
							{
								TAlScore bestLast = anchor1 ? prm.bestLtMinscMate2 : prm.bestLtMinscMate1;
								if (bestCell > std::numeric_limits<TAlScore>::min() && bestCell > bestLast)
								{
									if (anchor1)
									{
										prm.bestLtMinscMate2 = bestCell;
									}
									else
									{
										prm.bestLtMinscMate1 = bestCell;
									}
								}
							}
						}
						bool didAnchor = false;
						do
						{
							oresGap_.reset();
							assert(oresGap_.empty());
							if (foundMate && oswa.done())
							{
								foundMate = false;
							}
							else if (foundMate)
							{
								oswa.nextAlignment(oresGap_, ominsc_cur, rnd);
								foundMate = !oresGap_.empty();
								assert(!foundMate || oresGap_.alres.matchesRef(
														 ord,
														 ref,
														 tmp_rf_,
														 tmp_rdseq_,
														 tmp_qseq_,
														 raw_refbuf_,
														 raw_destU32_,
														 raw_matches_));
							}
							if (foundMate)
							{
								if (!redAnchor_.overlap(oresGap_.alres))
								{
									redAnchor_.add(oresGap_.alres);
								}
								assert_eq(ofw, oresGap_.alres.fw());
								oresGap_.alres.setParams(
									seedmms,
									seedlen,
									seedival,
									ominsc);
								assert_gt(oresGap_.alres.refExtent(), 0);
								if (gReportOverhangs &&
									!refival.containsIgnoreOrient(oresGap_.alres.refival()))
								{
									oresGap_.alres.clipOutside(true, 0, tlen);
									foundMate = oresGap_.alres.refExtent() > 0;
								}
								if (foundMate &&
									((!gReportOverhangs &&
									  !refival.containsIgnoreOrient(oresGap_.alres.refival())) ||
									 !refival.overlapsIgnoreOrient(oresGap_.alres.refival())))
								{
									foundMate = false;
								}
							}
							ASSERT_ONLY(TRefId refid);
							TRefOff off1, off2;
							size_t len1, len2;
							bool fw1, fw2;
							int pairCl = PE_ALS_DISCORD;
							if (foundMate)
							{
								ASSERT_ONLY(refid =)
								res->alres.refid();
								assert_eq(refid, oresGap_.alres.refid());
								off1 = anchor1 ? off : oresGap_.alres.refoff();
								off2 = anchor1 ? oresGap_.alres.refoff() : off;
								len1 = anchor1 ? res->alres.refExtent() : oresGap_.alres.refExtent();
								len2 = anchor1 ? oresGap_.alres.refExtent() : res->alres.refExtent();
								fw1 = anchor1 ? res->alres.fw() : oresGap_.alres.fw();
								fw2 = anchor1 ? oresGap_.alres.fw() : res->alres.fw();
								pairCl = pepol.peClassifyPair(
									off1,
									len1,
									fw1,
									off2,
									len2,
									fw2);
							}
							if (msink->state().doneConcordant())
							{
								foundMate = false;
							}
							if (reportImmediately)
							{
								if (foundMate)
								{
									assert(!msink->state().doneConcordant());
									assert(msink != NULL);
									assert(res->repOk());
									assert(oresGap_.repOk());
									assert(!msink->maxed());
									assert(!msink->state().done());
									bool doneUnpaired = false;
									if (!anchor1 || !didAnchor)
									{
										if (anchor1)
										{
											didAnchor = true;
										}
										const AlnRes &r1 = anchor1 ? res->alres : oresGap_.alres;
										if (!redMate1_.overlap(r1))
										{
											redMate1_.add(r1);
											if (msink->report(0, &r1, NULL))
											{
												doneUnpaired = true;
											}
										}
									}
									if (anchor1 || !didAnchor)
									{
										if (!anchor1)
										{
											didAnchor = true;
										}
										const AlnRes &r2 = anchor1 ? oresGap_.alres : res->alres;
										if (!redMate2_.overlap(r2))
										{
											redMate2_.add(r2);
											if (msink->report(0, NULL, &r2))
											{
												doneUnpaired = true;
											}
										}
									}
									bool donePaired = false;
									if (pairCl != PE_ALS_DISCORD)
									{
										foundConcordant = true;
										if (msink->report(
												0,
												anchor1 ? &res->alres : &oresGap_.alres,
												anchor1 ? &oresGap_.alres : &res->alres))
										{
											donePaired = true;
										}
										else
										{
											if (tighten > 0 && msink->Mmode() && msink->hasSecondBestPair())
											{
												TAlScore ps;
												if (tighten == 1)
												{
													ps = msink->bestPair();
												}
												else if (tighten == 2)
												{
													ps = msink->secondBestPair();
												}
												else
												{
													TAlScore diff = msink->bestPair() - msink->secondBestPair();
													ps = msink->secondBestPair() + (diff * 3) / 4;
												}
												if (tighten == 1 && ps < bestPairScore &&
													msink->bestPair() == msink->secondBestPair())
												{
													ps++;
												}
												if (tighten >= 2 && ps < bestPairScore)
												{
													ps++;
												}
												TAlScore nc = ps - operfectScore;
												if (nc > minsc)
												{
													minsc = nc;
													assert_leq(minsc, perfectScore);
													if (minsc > res->alres.score().score())
													{
														break;
													}
												}
												assert_leq(minsc, perfectScore);
											}
										}
									}
									if (donePaired || doneUnpaired)
									{
										return EXTEND_POLICY_FULFILLED;
									}
									if (msink->state().doneWithMate(anchor1))
									{
										return EXTEND_POLICY_FULFILLED;
									}
								}
								else if ((mixed || discord) && !didAnchor)
								{
									didAnchor = true;
									assert(msink != NULL);
									assert(res->repOk());
									assert(res->alres.matchesRef(
										rd,
										ref,
										tmp_rf_,
										tmp_rdseq_,
										tmp_qseq_,
										raw_refbuf_,
										raw_destU32_,
										raw_matches_));
									assert(!msink->maxed());
									assert(!msink->state().done());
									if (!msink->state().doneUnpaired(anchor1))
									{
										const AlnRes &r = res->alres;
										RedundantAlns &red = anchor1 ? redMate1_ : redMate2_;
										const AlnRes *r1 = anchor1 ? &res->alres : NULL;
										const AlnRes *r2 = anchor1 ? NULL : &res->alres;
										if (!red.overlap(r))
										{
											red.add(r);
											if (msink->report(0, r1, r2))
											{
												return EXTEND_POLICY_FULFILLED;
											}
										}
									}
									if (msink->state().doneWithMate(anchor1))
									{
										return EXTEND_POLICY_FULFILLED;
									}
								}
							}
						} while (!oresGap_.empty());
					}
					else if (found)
					{
						assert(!msink->state().doneWithMate(anchor1));
						if (reportImmediately && (mixed || discord))
						{
							assert(msink != NULL);
							assert(res->repOk());
							assert(res->alres.matchesRef(
								rd,
								ref,
								tmp_rf_,
								tmp_rdseq_,
								tmp_qseq_,
								raw_refbuf_,
								raw_destU32_,
								raw_matches_));
							assert(!msink->maxed());
							assert(!msink->state().done());
							if (!msink->state().doneUnpaired(anchor1))
							{
								const AlnRes &r = res->alres;
								RedundantAlns &red = anchor1 ? redMate1_ : redMate2_;
								const AlnRes *r1 = anchor1 ? &res->alres : NULL;
								const AlnRes *r2 = anchor1 ? NULL : &res->alres;
								if (!red.overlap(r))
								{
									red.add(r);
									if (msink->report(0, r1, r2))
									{
										return EXTEND_POLICY_FULFILLED;
									}
								}
							}
							if (msink->state().doneWithMate(anchor1))
							{
								return EXTEND_POLICY_FULFILLED;
							}
						}
					}
				}
				if (foundConcordant)
				{
					prm.nMateDpSuccs++;
					mateStreaks_[i] = 0;
					if (state == FOUND_UNGAPPED)
					{
						assert_gt(prm.nUgFail, 0);
						assert_gt(prm.nExUgFails, 0);
						prm.nExUgFails--;
						prm.nExUgSuccs++;
						prm.nUgLastSucc = prm.nExUgs - 1;
						if (prm.nUgFail > prm.nUgFailStreak)
						{
							prm.nUgFailStreak = prm.nUgFail;
						}
						prm.nUgFail = 0;
					}
					else if (state == FOUND_EE)
					{
						assert_gt(prm.nEeFail, 0);
						assert_gt(prm.nExEeFails, 0);
						prm.nExEeFails--;
						prm.nExEeSuccs++;
						prm.nEeLastSucc = prm.nExEes - 1;
						if (prm.nEeFail > prm.nEeFailStreak)
						{
							prm.nEeFailStreak = prm.nEeFail;
						}
						prm.nEeFail = 0;
					}
					else
					{
						assert_gt(prm.nDpFail, 0);
						assert_gt(prm.nExDpFails, 0);
						prm.nExDpFails--;
						prm.nExDpSuccs++;
						prm.nDpLastSucc = prm.nExDps - 1;
						if (prm.nDpFail > prm.nDpFailStreak)
						{
							prm.nDpFailStreak = prm.nDpFail;
						}
						prm.nDpFail = 0;
					}
				}
				else
				{
					prm.nMateDpFails++;
					mateStreaks_[i]++;
				}
			}
		}
	}
	return EXTEND_EXHAUSTED_CANDIDATES;
}
