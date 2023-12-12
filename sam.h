#ifndef SAM_H_
#define SAM_H_
#include <string>
#include "ds.h"
#include "read.h"
#include "util.h"
#include "aligner_result.h"
#include "scoring.h"
enum
{
	SAM_FLAG_PAIRED = 1,
	SAM_FLAG_MAPPED_PAIRED = 2,
	SAM_FLAG_UNMAPPED = 4,
	SAM_FLAG_MATE_UNMAPPED = 8,
	SAM_FLAG_QUERY_STRAND = 16,
	SAM_FLAG_MATE_STRAND = 32,
	SAM_FLAG_FIRST_IN_PAIR = 64,
	SAM_FLAG_SECOND_IN_PAIR = 128,
	SAM_FLAG_NOT_PRIMARY = 256,
	SAM_FLAG_FAILS_CHECKS = 512,
	SAM_FLAG_DUPLICATE = 1024
};
class AlnRes;
class AlnFlags;
class AlnSetSumm;
class SamConfig
{
	typedef EList<std::string> StrList;
	typedef EList<size_t> LenList;
public:
	SamConfig(
		const StrList &refnames,
		const LenList &reflens,
		bool truncQname,
		bool appendComment,
		bool omitsec,
		bool noUnal,
		const std::string &pg_pn,
		const std::string &pg_vn,
		const std::string &pg_cl,
		const std::string &rgs,
		bool print_as,
		bool print_xs,
		bool print_xss,
		bool print_yn,
		bool print_xn,
		bool print_x0,
		bool print_x1,
		bool print_xm,
		bool print_xo,
		bool print_xg,
		bool print_nm,
		bool print_md,
		bool print_yf,
		bool print_yi,
		bool print_ym,
		bool print_yp,
		bool print_yt,
		bool print_ys,
		bool print_zs,
		bool print_xr,
		bool print_xt,
		bool print_xd,
		bool print_xu,
		bool print_ye,
		bool print_yl,
		bool print_yu,
		bool print_xp,
		bool print_yr,
		bool print_zb,
		bool print_zr,
		bool print_zf,
		bool print_zm,
		bool print_zi,
		bool print_zp,
		bool print_zu,
		bool print_zt) : truncQname_(truncQname),
						 appendComment_(appendComment),
						 omitsec_(omitsec),
						 noUnal_(noUnal),
						 pg_pn_(pg_pn),
						 pg_vn_(pg_vn),
						 pg_cl_(pg_cl),
						 rgs_(rgs),
						 refnames_(refnames),
						 reflens_(reflens),
						 print_as_(print_as),
						 print_xs_(print_xs),
						 print_xss_(print_xss),
						 print_yn_(print_yn),
						 print_xn_(print_xn),
						 print_x0_(print_x0),
						 print_x1_(print_x1),
						 print_xm_(print_xm),
						 print_xo_(print_xo),
						 print_xg_(print_xg),
						 print_nm_(print_nm),
						 print_md_(print_md),
						 print_yf_(print_yf),
						 print_yi_(print_yi),
						 print_ym_(print_ym),
						 print_yp_(print_yp),
						 print_yt_(print_yt),
						 print_ys_(print_ys),
						 print_zs_(print_zs),
						 print_xr_(print_xr),
						 print_xt_(print_xt),
						 print_xd_(print_xd),
						 print_xu_(print_xu),
						 print_ye_(print_ye),
						 print_yl_(print_yl),
						 print_yu_(print_yu),
						 print_xp_(print_xp),
						 print_yr_(print_yr),
						 print_zb_(print_zb),
						 print_zr_(print_zr),
						 print_zf_(print_zf),
						 print_zm_(print_zm),
						 print_zi_(print_zi),
						 print_zp_(print_zp),
						 print_zu_(print_zu),
						 print_zt_(print_zt)
	{
		assert_eq(refnames_.size(), reflens_.size());
	}
	void printRefName(
		BTString &o,
		const std::string &name)
		const;
	template <typename T>
	void printOptFieldEscapedZ(BTString &o, const T &s) const
	{
		size_t len = s.length();
		for (size_t i = 0; i < len; i++)
		{
			if (s[i] < 33 || s[i] > 126 || s[i] == ':' || s[i] == '%')
			{
				o.append('%');
				int ms = s[i] >> 4;
				int ls = s[i] & 15;
				assert_range(0, 15, ms);
				assert_range(0, 15, ls);
				o.append("0123456789ABCDEF"[ms]);
				o.append("0123456789ABCDEF"[ls]);
			}
			else
			{
				o.append(s[i]);
			}
		}
	}
	template <typename T>
	void printOptFieldNewlineEscapedZ(BTString &o, const T &s) const
	{
		size_t len = s.length();
		for (size_t i = 0; i < len; i++)
		{
			if (s[i] == 10 || s[i] == 13 || s[i] == '%')
			{
				o.append('%');
				int ms = s[i] >> 4;
				int ls = s[i] & 15;
				assert_range(0, 15, ms);
				assert_range(0, 15, ls);
				o.append("0123456789ABCDEF"[ms]);
				o.append("0123456789ABCDEF"[ls]);
			}
			else
			{
				o.append(s[i]);
			}
		}
	}
	template <typename TStr>
	void printReadName(
		BTString &o,
		const TStr &name,
		bool omitSlashMate)
		const
	{
		size_t namelen = name.length();
		if (truncQname_ && namelen > 255)
		{
			namelen = 255;
		}
		for (size_t i = 0; i < namelen; i++)
		{
			if (truncQname_ && isspace(name[i]))
			{
				break;
			}
			o.append(name[i]);
		}
		size_t olen = o.length();
		if (omitSlashMate &&
			olen >= 2 &&
			o[olen - 2] == '/' &&
			(o[olen - 1] == '1' || o[olen - 1] == '2' || o[olen - 1] == '3'))
		{
			o.resize(olen - 2);
		}
	}
	void printRefNameFromIndex(
		BTString &o,
		size_t i)
		const;
	void printHeader(
		BTString &o,
		const std::string &rgid,
		const std::string &rgs,
		bool printHd,
		bool printSq,
		bool printPg)
		const;
	void printHdLine(BTString &o, const char *samver) const;
	void printSqLines(BTString &o) const;
	void printPgLine(BTString &o) const;
	void printAlignedOptFlags(
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
		const;
	void printEmptyOptFlags(
		BTString &o,
		bool first,
		const Read &rd,
		const AlnFlags &flags,
		const AlnSetSumm &summ,
		const SeedAlSumm &ssm,
		const PerReadMetrics &prm,
		const Scoring &sc)
		const;
	void printPreservedOptFlags(BTString &o, const Read &rd) const;
	template <typename TStr>
	void printComment(BTString &o, TStr &name) const
	{
		if (appendComment_)
		{
			size_t i;
			for (i = 0; i < name.length() && !isspace(name[i]); i++)
				;
			o.append('\t');
			o.append(name.toZBuf() + i + 1);
		}
	}
	bool omitSecondarySeqQual() const
	{
		return omitsec_;
	}
	bool omitUnalignedReads() const
	{
		return noUnal_;
	}
protected:
	bool truncQname_;
	bool appendComment_;
	bool omitsec_;
	bool noUnal_;
	std::string pg_id_;
	std::string pg_pn_;
	std::string pg_vn_;
	std::string pg_cl_;
	std::string rgs_;
	const StrList &refnames_;
	const LenList &reflens_;
	bool print_as_;
	bool print_xs_;
	bool print_xss_;
	bool print_yn_;
	bool print_xn_;
	bool print_x0_;
	bool print_x1_;
	bool print_xm_;
	bool print_xo_;
	bool print_xg_;
	bool print_nm_;
	bool print_md_;
	bool print_yf_;
	bool print_yi_;
	bool print_ym_;
	bool print_yp_;
	bool print_yt_;
	bool print_ys_;
	bool print_zs_;
	bool print_xr_;
	bool print_xt_;
	bool print_xd_;
	bool print_xu_;
	bool print_ye_;
	bool print_yl_;
	bool print_yu_;
	bool print_xp_;
	bool print_yr_;
	bool print_zb_;
	bool print_zr_;
	bool print_zf_;
	bool print_zm_;
	bool print_zi_;
	bool print_zp_;
	bool print_zu_;
	bool print_zt_;
};
#endif
