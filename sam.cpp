#include <string>
#include <sys/time.h>
#include "sam.h"
#include "filebuf.h"
using namespace std;
void SamConfig::printRefName(
	BTString &o,
	const std::string &name) const
{
	size_t namelen = name.length();
	for (size_t i = 0; i < namelen; i++)
	{
		if (isspace(name[i]))
		{
			return;
		}
		o.append(name[i]);
	}
}
void SamConfig::printRefNameFromIndex(BTString &o, size_t i) const
{
	printRefName(o, refnames_[i]);
}
void SamConfig::printHeader(
	BTString &o,
	const string &rgid,
	const string &rgs,
	bool printHd,
	bool printSq,
	bool printPg) const
{
	if (printHd)
		printHdLine(o, "1.0");
	if (printSq)
		printSqLines(o);
	if (!rgid.empty())
	{
		o.append("@RG");
		o.append(rgid.c_str());
		o.append(rgs.c_str());
		o.append('\n');
	}
	if (printPg)
		printPgLine(o);
}
void SamConfig::printHdLine(BTString &o, const char *samver) const
{
	o.append("@HD\tVN:");
	o.append(samver);
	o.append("\tSO:unsorted\n");
}
void SamConfig::printSqLines(BTString &o) const
{
	char buf[1024];
	for (size_t i = 0; i < refnames_.size(); i++)
	{
		o.append("@SQ\tSN:");
		printRefName(o, refnames_[i]);
		o.append("\tLN:");
		itoa10<size_t>(reflens_[i], buf);
		o.append(buf);
		o.append('\n');
	}
}
void SamConfig::printPgLine(BTString &o) const
{
	o.append("@PG\tID:");
	o.append(pg_pn_.c_str());
	o.append("\tPN:");
	o.append(pg_pn_.c_str());
	o.append("\tVN:");
	o.append(pg_vn_.c_str());
	o.append("\tCL:\"");
	o.append(pg_cl_.c_str());
	o.append('"');
	o.append('\n');
}
#define WRITE_SEP()         \
	{                       \
		if (!first)         \
			o.append('\t'); \
		first = false;      \
	}
void SamConfig::printAlignedOptFlags(
	BTString &o,
	bool first,
	const Read &rd,
	const Read *rdo,
	AlnRes &res,
	StackedAln &staln,
	const AlnFlags &flags,
	const AlnSetSumm &summ,
	const SeedAlSumm &ssm,
	const PerReadMetrics &prm,
	const Scoring &sc,
	const char *mapqInp)
	const
{
	char buf[1024];
	assert(summ.bestScore(rd.mate < 2).valid());
	if (print_as_)
	{
		itoa10<TAlScore>(res.score().score(), buf);
		WRITE_SEP();
		o.append("AS:i:");
		o.append(buf);
	}
	if (print_xs_)
	{
		AlnScore sco;
		if (flags.partOfPair())
		{
			sco = summ.bestUnchosenPScore(rd.mate < 2);
		}
		else
		{
			sco = summ.bestUnchosenUScore();
		}
		if (sco.valid())
		{
			itoa10<TAlScore>(sco.score(), buf);
			WRITE_SEP();
			o.append("XS:i:");
			o.append(buf);
		}
	}
	if (print_xn_)
	{
		itoa10<size_t>(res.refNs(), buf);
		WRITE_SEP();
		o.append("XN:i:");
		o.append(buf);
	}
	if (print_x0_)
	{
	}
	if (print_x1_)
	{
	}
	size_t num_mm = 0;
	size_t num_go = 0;
	size_t num_gx = 0;
	for (size_t i = 0; i < res.ned().size(); i++)
	{
		if (res.ned()[i].isMismatch())
		{
			num_mm++;
		}
		else if (res.ned()[i].isReadGap())
		{
			num_go++;
			num_gx++;
			while (i < res.ned().size() - 1 &&
				   res.ned()[i + 1].pos == res.ned()[i].pos &&
				   res.ned()[i + 1].isReadGap())
			{
				i++;
				num_gx++;
			}
		}
		else if (res.ned()[i].isRefGap())
		{
			num_go++;
			num_gx++;
			while (i < res.ned().size() - 1 &&
				   res.ned()[i + 1].pos == res.ned()[i].pos + 1 &&
				   res.ned()[i + 1].isRefGap())
			{
				i++;
				num_gx++;
			}
		}
	}
	if (print_xm_)
	{
		itoa10<size_t>(num_mm, buf);
		WRITE_SEP();
		o.append("XM:i:");
		o.append(buf);
	}
	if (print_xo_)
	{
		itoa10<size_t>(num_go, buf);
		WRITE_SEP();
		o.append("XO:i:");
		o.append(buf);
	}
	if (print_xg_)
	{
		itoa10<size_t>(num_gx, buf);
		WRITE_SEP();
		o.append("XG:i:");
		o.append(buf);
	}
	if (print_nm_)
	{
		itoa10<size_t>(res.ned().size(), buf);
		WRITE_SEP();
		o.append("NM:i:");
		o.append(buf);
	}
	if (print_md_)
	{
		WRITE_SEP();
		o.append("MD:Z:");
		staln.buildMdz();
		staln.writeMdz(
			&o,
			NULL);
	}
	if (print_ys_ && summ.paired())
	{
		assert(res.oscore().valid());
		itoa10<TAlScore>(res.oscore().score(), buf);
		WRITE_SEP();
		o.append("YS:i:");
		o.append(buf);
	}
	if (print_yn_)
	{
		TAlScore mn = sc.scoreMin.f<TAlScore>(rd.length());
		itoa10<TAlScore>(mn, buf);
		WRITE_SEP();
		o.append("YN:i:");
		o.append(buf);
		TAlScore pe = sc.perfectScore(rd.length());
		itoa10<TAlScore>(pe, buf);
		WRITE_SEP();
		o.append("Yn:i:");
		o.append(buf);
		if (summ.paired())
		{
			assert(rdo != NULL);
			TAlScore mn = sc.scoreMin.f<TAlScore>(rdo->length());
			itoa10<TAlScore>(mn, buf);
			WRITE_SEP();
			o.append("ZN:i:");
			o.append(buf);
			TAlScore pe = sc.perfectScore(rdo->length());
			itoa10<TAlScore>(pe, buf);
			WRITE_SEP();
			o.append("Zn:i:");
			o.append(buf);
		}
	}
	if (print_xss_)
	{
		bool one = true;
		if (flags.partOfPair() && !flags.readMate1())
		{
			one = false;
		}
		TAlScore bst = one ? prm.bestLtMinscMate1 : prm.bestLtMinscMate2;
		if (bst > std::numeric_limits<TAlScore>::min())
		{
			itoa10<TAlScore>(bst, buf);
			WRITE_SEP();
			o.append("Xs:i:");
			o.append(buf);
		}
		if (flags.partOfPair())
		{
			bst = one ? prm.bestLtMinscMate2 : prm.bestLtMinscMate1;
			if (bst > std::numeric_limits<TAlScore>::min())
			{
				itoa10<TAlScore>(bst, buf);
				WRITE_SEP();
				o.append("Ys:i:");
				o.append(buf);
			}
		}
	}
	if (print_zs_)
	{
		itoa10<uint32_t>(rd.seed, buf);
		WRITE_SEP();
		o.append("ZS:i:");
		o.append(buf);
	}
	if (print_yt_)
	{
		WRITE_SEP();
		flags.printYT(o);
	}
	if (print_yp_ && flags.partOfPair() && flags.canMax())
	{
		WRITE_SEP();
		flags.printYP(o);
	}
	if (print_ym_ && flags.canMax() && (flags.isMixedMode() || !flags.partOfPair()))
	{
		WRITE_SEP();
		flags.printYM(o);
	}
	if (print_yf_ && flags.filtered())
	{
		first = flags.printYF(o, first) && first;
	}
	if (print_yi_)
	{
		if (mapqInp[0] != '\0')
		{
			WRITE_SEP();
			o.append("YI:Z:");
			o.append(mapqInp);
		}
	}
	if (flags.partOfPair() && print_zp_)
	{
		if (summ.bestCScore().valid())
		{
			WRITE_SEP();
			o.append("ZP:i:");
			itoa10<TAlScore>(summ.bestCScore().score(), buf);
			o.append(buf);
		}
		if (summ.bestUnchosenCScore().valid())
		{
			WRITE_SEP();
			o.append("Zp:i:");
			itoa10<TAlScore>(summ.bestUnchosenCScore().score(), buf);
			o.append(buf);
		}
	}
	if (print_zu_)
	{
		AlnScore best = summ.bestScore(rd.mate <= 1);
		AlnScore secbest = summ.bestUnchosenPScore(rd.mate <= 1);
		WRITE_SEP();
		o.append("ZU:i:");
		if (best.valid())
		{
			itoa10<TAlScore>(best.score(), buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
		WRITE_SEP();
		o.append("Zu:i:");
		if (secbest.valid())
		{
			itoa10<TAlScore>(secbest.score(), buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
	}
	if (!rgs_.empty())
	{
		WRITE_SEP();
		o.append(rgs_.c_str());
	}
	if (print_xt_)
	{
		WRITE_SEP();
		struct timeval tv_end;
		struct timezone tz_end;
		gettimeofday(&tv_end, &tz_end);
		size_t total_usecs =
			(tv_end.tv_sec - prm.tv_beg.tv_sec) * 1000000 +
			(tv_end.tv_usec - prm.tv_beg.tv_usec);
		itoa10<size_t>(total_usecs, buf);
		o.append("XT:i:");
		o.append(buf);
	}
	if (print_xd_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nExDps, buf);
		o.append("XD:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nMateDps, buf);
		o.append("Xd:i:");
		o.append(buf);
	}
	if (print_xu_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nExUgs, buf);
		o.append("XU:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nMateUgs, buf);
		o.append("Xu:i:");
		o.append(buf);
	}
	if (print_ye_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nDpFail, buf);
		o.append("YE:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nUgFail, buf);
		o.append("Ye:i:");
		o.append(buf);
	}
	if (print_yl_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nDpFailStreak, buf);
		o.append("YL:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nUgFailStreak, buf);
		o.append("Yl:i:");
		o.append(buf);
	}
	if (print_yu_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nDpLastSucc, buf);
		o.append("YU:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nUgLastSucc, buf);
		o.append("Yu:i:");
		o.append(buf);
	}
	if (print_xp_)
	{
		WRITE_SEP();
		o.append("XP:B:I,");
		itoa10<uint64_t>(prm.nSeedElts, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.nSeedEltsFw, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.nSeedEltsRc, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.seedMean, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.seedMedian, buf);
		o.append(buf);
	}
	if (print_yr_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nRedundants, buf);
		o.append("YR:i:");
		o.append(buf);
	}
	if (print_zb_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nFtabs, buf);
		o.append("ZB:i:");
		o.append(buf);
	}
	if (print_zr_)
	{
		WRITE_SEP();
		o.append("ZR:Z:");
		itoa10<uint64_t>(prm.nRedSkip, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.nRedFail, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.nRedIns, buf);
		o.append(buf);
	}
	if (print_zf_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nSdFmops, buf);
		o.append("ZF:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nExFmops, buf);
		o.append("Zf:i:");
		o.append(buf);
	}
	if (print_zm_)
	{
		WRITE_SEP();
		o.append("ZM:Z:");
		prm.fmString.print(o, buf);
	}
	if (print_zi_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nExIters, buf);
		o.append("ZI:i:");
		o.append(buf);
	}
	if (print_xr_)
	{
		o.append("\n");
		printOptFieldNewlineEscapedZ(o, rd.readOrigBuf);
	}
	if (print_zt_)
	{
		WRITE_SEP();
		const bool paired = flags.partOfPair();
		const TAlScore MN = std::numeric_limits<TAlScore>::min();
		TAlScore secondBest[2] = {MN, MN};
		TAlScore thirdBest[2] = {MN, MN};
		const int ED_MAX = std::numeric_limits<int>::max();
		AlnScore best[2] = {res.score(), res.oscore()};
		TAlScore diffEd[2] = {ED_MAX, ED_MAX};
		for (int self = 0; self < (paired ? 2 : 1); self++)
		{
			AlnScore sco;
			bool mate1 = rd.mate < 2;
			if (self > 0)
				mate1 = !mate1;
			if (flags.partOfPair())
			{
				sco = summ.bestUnchosenPScore(mate1);
			}
			else
			{
				sco = summ.bestUnchosenUScore();
			}
			if (sco.valid())
			{
				secondBest[self] = sco.score();
			}
			thirdBest[self] = mate1 ? prm.bestLtMinscMate1 : prm.bestLtMinscMate2;
			if (flags.partOfPair())
			{
				if (summ.bestUnchosenPDist(mate1).valid())
				{
					diffEd[self] = best[self].basesAligned() - summ.bestUnchosenPDist(mate1).basesAligned();
				}
			}
			else
			{
				if (summ.bestUnchosenUDist().valid())
				{
					diffEd[self] = best[self].basesAligned() - summ.bestUnchosenUDist().basesAligned();
				}
			}
		}
		TAlScore diff[2] = {MN, MN};
		for (int self = 0; self < 2; self++)
		{
			const TAlScore mx = max(secondBest[self], thirdBest[self]);
			if (best[self].score() > MN && mx > MN)
			{
				diff[self] = best[self].score() - mx;
			}
		}
		TAlScore best_conc = MN, diff_conc = MN;
		int diffEd_conc = ED_MAX;
		if (paired && summ.bestCScore().valid())
		{
			best_conc = summ.bestCScore().score();
			if (summ.bestUnchosenCScore().valid())
			{
				diff_conc = best_conc - summ.bestUnchosenCScore().score();
			}
			if (summ.bestUnchosenCDist().valid())
			{
				diffEd_conc = summ.bestCDist().basesAligned() - summ.bestUnchosenCDist().basesAligned();
			}
		}
		o.append("ZT:Z:");
		itoa10<TAlScore>((int)best[0].score(), buf);
		o.append(buf);
		o.append(",");
		if (diff[0] > MN)
		{
			itoa10<TAlScore>((int)diff[0], buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
		o.append(",");
		if (diffEd[0] != ED_MAX)
		{
			itoa10<TAlScore>((int)diffEd[0], buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
		o.append(",");
		if (best[1].score() > MN)
		{
			itoa10<TAlScore>((int)best[1].score(), buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
		o.append(",");
		if (diff[1] > MN)
		{
			itoa10<TAlScore>((int)diff[1], buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
		o.append(",");
		if (best_conc > MN)
		{
			itoa10<TAlScore>((int)best_conc, buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
		o.append(",");
		if (diff_conc > MN)
		{
			itoa10<TAlScore>((int)diff_conc, buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
		o.append(",");
		if (diffEd_conc != ED_MAX)
		{
			itoa10<TAlScore>((int)diffEd_conc, buf);
			o.append(buf);
		}
		else
		{
			o.append("NA");
		}
		int mate = (rd.mate < 2 ? 0 : 1);
		o.append(",");
		itoa10<TAlScore>((int)((prm.seedsPerNucMS[2 * mate] + prm.seedsPerNucMS[2 * mate + 1]) * 1000), buf);
		o.append(buf);
		o.append(",");
		itoa10<TAlScore>((int)((prm.seedPctUniqueMS[2 * mate] + prm.seedPctUniqueMS[2 * mate + 1]) * 1000), buf);
		o.append(buf);
		o.append(",");
		itoa10<TAlScore>((int)((prm.seedPctRepMS[2 * mate] + prm.seedPctRepMS[2 * mate + 1]) * 1000), buf);
		o.append(buf);
		o.append(",");
		itoa10<TAlScore>((int)((prm.seedHitAvgMS[2 * mate] + prm.seedHitAvgMS[2 * mate + 1]) + 0.5f), buf);
		o.append(buf);
		int fw = res.fw() ? 0 : 1;
		o.append(",");
		itoa10<TAlScore>((int)(prm.seedsPerNucMS[2 * mate + fw] * 1000), buf);
		o.append(buf);
		o.append(",");
		itoa10<TAlScore>((int)(prm.seedPctUniqueMS[2 * mate + fw] * 1000), buf);
		o.append(buf);
		o.append(",");
		itoa10<TAlScore>((int)(prm.seedPctRepMS[2 * mate + fw] * 1000), buf);
		o.append(buf);
		o.append(",");
		itoa10<TAlScore>((int)(prm.seedHitAvgMS[2 * mate + fw] + 0.5f), buf);
		o.append(buf);
	}
}
void SamConfig::printEmptyOptFlags(
	BTString &o,
	bool first,
	const Read &rd,
	const AlnFlags &flags,
	const AlnSetSumm &summ,
	const SeedAlSumm &ssm,
	const PerReadMetrics &prm,
	const Scoring &sc)
	const
{
	char buf[1024];
	if (print_yn_)
	{
		TAlScore mn = sc.scoreMin.f<TAlScore>(rd.length());
		itoa10<TAlScore>(mn, buf);
		WRITE_SEP();
		o.append("YN:i:");
		o.append(buf);
		TAlScore pe = sc.perfectScore(rd.length());
		itoa10<TAlScore>(pe, buf);
		WRITE_SEP();
		o.append("Yn:i:");
		o.append(buf);
	}
	if (print_zs_)
	{
		itoa10<uint32_t>(rd.seed, buf);
		WRITE_SEP();
		o.append("ZS:i:");
		o.append(buf);
	}
	if (print_yt_)
	{
		WRITE_SEP();
		flags.printYT(o);
	}
	if (print_yp_ && flags.partOfPair() && flags.canMax())
	{
		WRITE_SEP();
		flags.printYP(o);
	}
	if (print_ym_ && flags.canMax() && (flags.isMixedMode() || !flags.partOfPair()))
	{
		WRITE_SEP();
		flags.printYM(o);
	}
	if (print_yf_ && flags.filtered())
	{
		first = flags.printYF(o, first) && first;
	}
	if (!rgs_.empty())
	{
		WRITE_SEP();
		o.append(rgs_.c_str());
	}
	if (print_xt_)
	{
		WRITE_SEP();
		struct timeval tv_end;
		struct timezone tz_end;
		gettimeofday(&tv_end, &tz_end);
		size_t total_usecs =
			(tv_end.tv_sec - prm.tv_beg.tv_sec) * 1000000 +
			(tv_end.tv_usec - prm.tv_beg.tv_usec);
		itoa10<size_t>(total_usecs, buf);
		o.append("XT:i:");
		o.append(buf);
	}
	if (print_xd_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nExDps, buf);
		o.append("XD:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nMateDps, buf);
		o.append("Xd:i:");
		o.append(buf);
	}
	if (print_xu_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nExUgs, buf);
		o.append("XU:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nMateUgs, buf);
		o.append("Xu:i:");
		o.append(buf);
	}
	if (print_ye_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nDpFail, buf);
		o.append("YE:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nUgFail, buf);
		o.append("Ye:i:");
		o.append(buf);
	}
	if (print_yl_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nDpFailStreak, buf);
		o.append("YL:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nUgFailStreak, buf);
		o.append("Yl:i:");
		o.append(buf);
	}
	if (print_yu_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nDpLastSucc, buf);
		o.append("YU:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nUgLastSucc, buf);
		o.append("Yu:i:");
		o.append(buf);
	}
	if (print_xp_)
	{
		WRITE_SEP();
		o.append("XP:B:I,");
		itoa10<uint64_t>(prm.nSeedElts, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.nSeedEltsFw, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.nSeedEltsRc, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.seedMean, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.seedMedian, buf);
		o.append(buf);
	}
	if (print_yr_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nRedundants, buf);
		o.append("YR:i:");
		o.append(buf);
	}
	if (print_zb_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nFtabs, buf);
		o.append("ZB:i:");
		o.append(buf);
	}
	if (print_zr_)
	{
		WRITE_SEP();
		o.append("ZR:Z:");
		itoa10<uint64_t>(prm.nRedSkip, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.nRedFail, buf);
		o.append(buf);
		o.append(',');
		itoa10<uint64_t>(prm.nRedIns, buf);
		o.append(buf);
	}
	if (print_zf_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nSdFmops, buf);
		o.append("ZF:i:");
		o.append(buf);
		WRITE_SEP();
		itoa10<uint64_t>(prm.nExFmops, buf);
		o.append("Zf:i:");
		o.append(buf);
	}
	if (print_zm_)
	{
		WRITE_SEP();
		o.append("ZM:Z:");
		prm.fmString.print(o, buf);
	}
	if (print_zi_)
	{
		WRITE_SEP();
		itoa10<uint64_t>(prm.nExIters, buf);
		o.append("ZI:i:");
		o.append(buf);
	}
	if (print_xr_)
	{
		o.append("\n");
		printOptFieldNewlineEscapedZ(o, rd.readOrigBuf);
	}
}
void SamConfig::printPreservedOptFlags(BTString &o, const Read &rd) const
{
	if (rd.preservedOptFlags.length() != 0)
	{
		char buf[1024];
		const char *b = rd.preservedOptFlags.buf();
		int i = 0, len = rd.preservedOptFlags.length();
		while (i < len)
		{
			o.append('\t');
			char tag[2], val_type;
			memcpy(tag, b + i, 2 * sizeof(char));
			o.append(tag, 2);
			i += 2 * sizeof(char);
			memcpy(&val_type, b + i, 1);
			o.append(':');
			if (val_type == 'c' || val_type == 'C' || val_type == 'i' || val_type == 'I' || val_type == 's' || val_type == 'S')
			{
				o.append('i');
			}
			else
			{
				o.append(val_type);
			}
			o.append(':');
			i += sizeof(char);
			switch (val_type)
			{
			case 'A':
				char A_val;
				memcpy(&A_val, b + i, sizeof(A_val));
				i += sizeof(A_val);
				itoa10<char>(A_val, buf);
				o.append(buf);
				break;
			case 'c':
				int8_t c_val;
				memcpy(&c_val, b + i, sizeof(c_val));
				i += sizeof(c_val);
				itoa10<int8_t>(c_val, buf);
				o.append(buf);
				break;
			case 'C':
				uint8_t C_val;
				memcpy(&C_val, b + i, sizeof(C_val));
				i += sizeof(C_val);
				itoa10<uint8_t>(C_val, buf);
				o.append(buf);
				break;
			case 's':
				int16_t s_val;
				memcpy(&s_val, b + i, sizeof(s_val));
				i += sizeof(s_val);
				itoa10<int16_t>(s_val, buf);
				o.append(buf);
				break;
			case 'S':
				uint16_t S_val;
				memcpy(&S_val, b + i, sizeof(S_val));
				i += sizeof(S_val);
				itoa10<uint16_t>(S_val, buf);
				o.append(buf);
				break;
			case 'i':
				int32_t i_val;
				memcpy(&i_val, b + i, sizeof(i_val));
				i += sizeof(i_val);
				itoa10<int32_t>(i_val, buf);
				o.append(buf);
				break;
			case 'I':
				uint32_t I_val;
				memcpy(&I_val, b + i, sizeof(I_val));
				i += sizeof(I_val);
				itoa10<uint32_t>(I_val, buf);
				o.append(buf);
				break;
			case 'Z':
				char c;
				memcpy(&c, b + i, sizeof(char));
				while (c != '\0')
				{
					o.append(c);
					i++;
					memcpy(&c, b + i, sizeof(char));
				}
				i++;
				break;
			default:
				break;
			}
		}
	}
}
