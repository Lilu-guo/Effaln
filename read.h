#ifndef READ_H_
#define READ_H_
#include <stdint.h>
#include <sys/time.h>
#include "ds.h"
#include "sstring.h"
#include "filebuf.h"
#include "util.h"
typedef uint64_t TReadId;
typedef size_t TReadOff;
typedef int64_t TAlScore;
class HitSet;
struct Read
{
	typedef SStringExpandable<char, 1024, 2, 1024> TBuf;
	Read() { reset(); }
	Read(const char *nm, const char *seq, const char *ql) { init(nm, seq, ql); }
	void reset()
	{
		rdid = 0;
		trimmed5 = trimmed3 = 0;
		readOrigBuf.clear();
		patFw.clear();
		patRc.clear();
		qual.clear();
		patFwRev.clear();
		patRcRev.clear();
		qualRev.clear();
		name.clear();
		preservedOptFlags.clear();
		filter = '?';
		seed = 0;
		parsed = false;
		ns_ = 0;
	}
	void finalize()
	{
		for (size_t i = 0; i < patFw.length(); i++)
		{
			if ((int)patFw[i] > 3)
			{
				ns_++;
			}
		}
		constructRevComps();
		constructReverses();
	}
	void init(
		const char *nm,
		const char *seq,
		const char *ql)
	{
		reset();
		patFw.installChars(seq);
		qual.install(ql);
		for (size_t i = 0; i < patFw.length(); i++)
		{
			if ((int)patFw[i] > 3)
			{
				ns_++;
			}
		}
		constructRevComps();
		constructReverses();
		if (nm != NULL)
			name.install(nm);
	}
	bool empty() const
	{
		return patFw.empty();
	}
	size_t length() const
	{
		return patFw.length();
	}
	size_t ns() const
	{
		return ns_;
	}
	void constructRevComps()
	{
		patRc.installReverseComp(patFw);
	}
	void constructReverses()
	{
		patFwRev.installReverse(patFw);
		patRcRev.installReverse(patRc);
		qualRev.installReverse(qual);
	}
	void fixMateName(int i)
	{
		assert(i == 1 || i == 2);
		size_t namelen = name.length();
		bool append = false;
		if (namelen < 2)
		{
			append = true;
		}
		else
		{
			if (i == 1)
			{
				append =
					name[namelen - 2] != '/' ||
					name[namelen - 1] != '1';
			}
			else
			{
				append =
					name[namelen - 2] != '/' ||
					name[namelen - 1] != '2';
			}
		}
		if (append)
		{
			name.append('/');
			name.append("012"[i]);
		}
	}
	void dump(std::ostream &os) const
	{
		using namespace std;
		os << name << ' ';
		os << patFw;
		os << ' ';
		os << qual.toZBuf() << " ";
	}
	static bool same(
		const BTDnaString &seq1,
		const BTString &qual1,
		const BTDnaString &seq2,
		const BTString &qual2,
		bool qualitiesMatter)
	{
		if (seq1.length() != seq2.length())
		{
			return false;
		}
		for (size_t i = 0; i < seq1.length(); i++)
		{
			if (seq1[i] != seq2[i])
				return false;
		}
		if (qualitiesMatter)
		{
			if (qual1.length() != qual2.length())
			{
				return false;
			}
			for (size_t i = 0; i < qual1.length(); i++)
			{
				if (qual1[i] != qual2[i])
					return false;
			}
		}
		return true;
	}
	std::pair<int, int> get(TReadOff off5p, bool fw) const
	{
		assert_lt(off5p, length());
		int c = (int)patFw[off5p];
		int q = qual[off5p];
		assert_geq(q, 33);
		return make_pair((!fw && c < 4) ? (c ^ 3) : c, q - 33);
	}
	int getc(TReadOff off5p, bool fw) const
	{
		assert_lt(off5p, length());
		int c = (int)patFw[off5p];
		return (!fw && c < 4) ? (c ^ 3) : c;
	}
	int getq(TReadOff off5p) const
	{
		assert_lt(off5p, length());
		int q = qual[off5p];
		assert_geq(q, 33);
		return q - 33;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		if (patFw.empty())
			return true;
		assert_eq(qual.length(), patFw.length());
		return true;
	}
#endif
	BTDnaString patFw;
	BTDnaString patRc;
	BTString qual;
	BTDnaString patFwRev;
	BTDnaString patRcRev;
	BTString qualRev;
	TBuf readOrigBuf;
	BTString name;
	BTString preservedOptFlags;
	TReadId rdid;
	int mate;
	uint32_t seed;
	bool parsed;
	size_t ns_;
	char filter;
	int trimmed5;
	int trimmed3;
	HitSet *hitset;
};
struct FmStringOp
{
	bool alignment;
	TAlScore pen;
	size_t n;
};
struct FmString
{
	void add(bool alignment, TAlScore pen, size_t nops)
	{
		if (ops.empty() || ops.back().pen != pen)
		{
			ops.expand();
			ops.back().alignment = alignment;
			ops.back().pen = pen;
			ops.back().n = 0;
		}
		ops.back().n++;
	}
	void reset()
	{
		pen = std::numeric_limits<TAlScore>::max();
		ops.clear();
	}
	void print(BTString &o, char *buf) const
	{
		for (size_t i = 0; i < ops.size(); i++)
		{
			if (i > 0)
			{
				o.append(';');
			}
			if (ops[i].alignment)
			{
				o.append("A,");
				itoa10(ops[i].pen, buf);
				o.append(buf);
			}
			else
			{
				o.append("F,");
				itoa10(ops[i].pen, buf);
				o.append(buf);
				o.append(',');
				itoa10(ops[i].n, buf);
				o.append(buf);
			}
		}
	}
	TAlScore pen;
	EList<FmStringOp> ops;
};
struct PerReadMetrics
{
	PerReadMetrics() { reset(); }
	void reset()
	{
		nExIters =
			nExDps = nExDpSuccs = nExDpFails =
				nMateDps = nMateDpSuccs = nMateDpFails =
					nExUgs = nExUgSuccs = nExUgFails =
						nMateUgs = nMateUgSuccs = nMateUgFails =
							nExEes = nExEeSuccs = nExEeFails =
								nRedundants =
									nEeFmops = nSdFmops = nExFmops =
										nDpFail = nDpFailStreak = nDpLastSucc =
											nUgFail = nUgFailStreak = nUgLastSucc =
												nEeFail = nEeFailStreak = nEeLastSucc =
													nFilt = 0;
		nFtabs = 0;
		nRedSkip = 0;
		nRedFail = 0;
		nRedIns = 0;
		doFmString = false;
		nSeedRanges = nSeedElts = 0;
		nSeedRangesFw = nSeedEltsFw = 0;
		nSeedRangesRc = nSeedEltsRc = 0;
		seedMedian = seedMean = 0;
		bestLtMinscMate1 =
			bestLtMinscMate2 = std::numeric_limits<TAlScore>::min();
		seedPctUnique = seedPctRep = seedsPerNuc = seedHitAvg = 0.0f;
		fmString.reset();
	}
	struct timeval tv_beg;
	struct timezone tz_beg;
	uint64_t nExIters;
	uint64_t nExDps;
	uint64_t nExDpSuccs;
	uint64_t nExDpFails;
	uint64_t nExUgs;
	uint64_t nExUgSuccs;
	uint64_t nExUgFails;
	uint64_t nExEes;
	uint64_t nExEeSuccs;
	uint64_t nExEeFails;
	uint64_t nMateDps;
	uint64_t nMateDpSuccs;
	uint64_t nMateDpFails;
	uint64_t nMateUgs;
	uint64_t nMateUgSuccs;
	uint64_t nMateUgFails;
	uint64_t nRedundants;
	uint64_t nSeedRanges;
	uint64_t nSeedElts;
	uint64_t nSeedRangesFw;
	uint64_t nSeedEltsFw;
	uint64_t nSeedRangesRc;
	uint64_t nSeedEltsRc;
	uint64_t seedMedian;
	uint64_t seedMean;
	uint64_t nEeFmops;
	uint64_t nSdFmops;
	uint64_t nExFmops;
	uint64_t nFtabs;
	uint64_t nRedSkip;
	uint64_t nRedFail;
	uint64_t nRedIns;
	uint64_t nDpFail;
	uint64_t nDpFailStreak;
	uint64_t nDpLastSucc;
	uint64_t nUgFail;
	uint64_t nUgFailStreak;
	uint64_t nUgLastSucc;
	uint64_t nEeFail;
	uint64_t nEeFailStreak;
	uint64_t nEeLastSucc;
	uint64_t nFilt;
	TAlScore bestLtMinscMate1;
	TAlScore bestLtMinscMate2;
	float seedPctUnique;
	float seedPctUniqueMS[4];
	float seedPctRep;
	float seedPctRepMS[4];
	float seedHitAvg;
	float seedHitAvgMS[4];
	float seedsPerNuc;
	float seedsPerNucMS[4];
	bool doFmString;
	FmString fmString;
};
#endif
