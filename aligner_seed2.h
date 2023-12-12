#ifndef ALIGNER_SEED2_H_
#define ALIGNER_SEED2_H_
#include <stdint.h>
#include <math.h>
#include <utility>
#include <limits>
#include "assert_helpers.h"
#include "random_util.h"
#include "aligner_result.h"
#include "aln_idx.h"
#include "simple_func.h"
#include "scoring.h"
#include "edit.h"
#include "read.h"
#include "ds.h"
#include "group_walk.h"
#include "btypes.h"
typedef size_t TReadOff;
typedef int64_t TScore;
typedef float TRootPri;
typedef size_t TDescentId;
typedef size_t TRootId;
enum
{
	DESC_EX_NONE = 1,
	DESC_EX_FROM_1ST_BRANCH = 2,
	DESC_EX_EACH_EDGE = 3
};
struct DescentMetrics
{
	DescentMetrics() { reset(); }
	void reset()
	{
		bwops = bwops_1 = bwops_bi = recalc = branch = branch_mm =
			branch_del = branch_ins = heap_max = descent_max = descentpos_max =
				nex = 0;
	}
	uint64_t bwops;
	uint64_t bwops_1;
	uint64_t bwops_bi;
	uint64_t recalc;
	uint64_t branch;
	uint64_t branch_mm;
	uint64_t branch_del;
	uint64_t branch_ins;
	uint64_t heap_max;
	uint64_t descent_max;
	uint64_t descentpos_max;
	uint64_t nex;
};
struct DescentPriority
{
	DescentPriority() { reset(); }
	DescentPriority(
		TScore pen_,
		size_t depth_,
		TIndexOffU width_,
		float rootpri_)
	{
		pen = pen_;
		depth = depth_;
		width = width_;
		rootpri = rootpri_;
	}
	void init(TScore pen_, size_t depth_, TIndexOffU width_, float rootpri_)
	{
		pen = pen_;
		depth = depth_;
		width = width_;
		rootpri = rootpri_;
	}
	void reset()
	{
		width = 0;
	}
	bool inited() const
	{
		return width > 0;
	}
	bool operator<(const DescentPriority &o) const
	{
		assert(inited());
		assert(o.inited());
		if (pen < o.pen)
			return true;
		if (pen > o.pen)
			return false;
		if (depth > o.depth)
			return true;
		if (depth < o.depth)
			return false;
		if (width < o.width)
			return true;
		if (width > o.width)
			return false;
		if (rootpri > o.rootpri)
			return true;
		return false;
	}
	bool operator<=(const DescentPriority &o) const
	{
		assert(inited());
		assert(o.inited());
		if (pen < o.pen)
			return true;
		if (pen > o.pen)
			return false;
		if (depth > o.depth)
			return true;
		if (depth < o.depth)
			return false;
		if (width < o.depth)
			return true;
		if (width > o.width)
			return false;
		if (rootpri > o.rootpri)
			return true;
		return true;
	}
	bool operator==(const DescentPriority &o) const
	{
		assert(inited());
		assert(o.inited());
		return pen == o.pen && depth == o.depth && width == o.width && rootpri == o.rootpri;
	}
	TScore pen;
	size_t depth;
	TIndexOffU width;
	float rootpri;
};
static inline std::ostream &operator<<(
	std::ostream &os,
	const DescentPriority &o)
{
	os << "[" << o.pen << ", " << o.depth << ", " << o.width << ", " << o.rootpri << "]";
	return os;
}
static inline std::ostream &operator<<(
	std::ostream &os,
	const std::pair<DescentPriority, TDescentId> &o)
{
	os << "{[" << o.first.pen << ", " << o.first.depth << ", "
	   << o.first.width << ", " << o.first.rootpri << "], " << o.second << "}";
	return os;
}
typedef std::pair<DescentPriority, TDescentId> TDescentPair;
struct DescentConstraints
{
	DescentConstraints() { reset(); }
	DescentConstraints(size_t nzero, double exp)
	{
		init(nzero, exp);
	}
	void init(size_t nzero_, double exp_)
	{
		nzero = nzero_ > 0 ? nzero_ : 1;
		exp = exp_;
#ifndef NDEBUG
		for (size_t i = 1; i < nzero_ + 5; i++)
		{
			assert_geq(get(i, nzero_ + 10, 100), get(i - 1, nzero_ + 10, 100));
		}
#endif
	}
	void reset()
	{
		nzero = 0;
		exp = -1.0f;
	}
	bool inited() const
	{
		return exp >= 0.0f;
	}
	inline TScore get(TReadOff off, TReadOff rdlen, TAlScore maxpen) const
	{
		if (off < nzero || nzero >= rdlen)
		{
			return 0;
		}
		double frac = (double)(off - nzero) / (rdlen - nzero);
		if (fabs(exp - 1.0f) > 0.00001)
		{
			if (fabs(exp - 2.0f) < 0.00001)
			{
				frac *= frac;
			}
			else
			{
				frac = pow(frac, exp);
			}
		}
		return (TAlScore)(frac * maxpen + 0.5f);
	}
	size_t nzero;
	double exp;
};
struct DescentConfig
{
	DescentConfig() { reset(); }
	void reset() { expol = 0; }
	bool inited() const { return expol != 0; }
	DescentConstraints cons;
	int expol;
};
struct DescentRedundancyKey
{
	DescentRedundancyKey() { reset(); }
	DescentRedundancyKey(
		TReadOff al5pf_,
		size_t rflen_,
		TIndexOffU topf_,
		TIndexOffU botf_)
	{
		init(al5pf_, rflen_, topf_, botf_);
	}
	void reset()
	{
		al5pf = 0;
		rflen = 0;
		topf = botf = 0;
	}
	bool inited() const { return rflen > 0; }
	void init(
		TReadOff al5pf_,
		size_t rflen_,
		TIndexOffU topf_,
		TIndexOffU botf_)
	{
		al5pf = al5pf_;
		rflen = rflen_;
		topf = topf_;
		botf = botf_;
	}
	bool operator==(const DescentRedundancyKey &o) const
	{
		return al5pf == o.al5pf && rflen == o.rflen && topf == o.topf && botf == o.botf;
	}
	bool operator<(const DescentRedundancyKey &o) const
	{
		if (al5pf < o.al5pf)
			return true;
		if (al5pf > o.al5pf)
			return false;
		if (rflen < o.rflen)
			return true;
		if (rflen > o.rflen)
			return false;
		if (topf < o.topf)
			return true;
		if (topf > o.topf)
			return false;
		return botf < o.botf;
	}
	TReadOff al5pf;
	size_t rflen;
	TIndexOffU topf;
	TIndexOffU botf;
};
class DescentRedundancyChecker
{
public:
	DescentRedundancyChecker() { reset(); }
	void clear() { reset(); }
	void reset()
	{
		bits_.reset();
		inited_ = false;
		totsz_ = 0;
		totcap_ = 0;
	}
	const static int NPARTS = 8;
	const static int PART_MASK = 7;
	const static int NBITS = (1 << 16);
	void init(TReadOff rdlen)
	{
		reset();
		bits_.resize(NBITS);
		maplist_fl_.resize(NPARTS);
		maplist_fr_.resize(NPARTS);
		maplist_rl_.resize(NPARTS);
		maplist_rr_.resize(NPARTS);
		for (int i = 0; i < NPARTS; i++)
		{
			maplist_fl_[i].resize(rdlen);
			maplist_fr_[i].resize(rdlen);
			maplist_rl_[i].resize(rdlen);
			maplist_rr_[i].resize(rdlen);
			totcap_ += maplist_fl_[i].totalCapacityBytes();
			totcap_ += maplist_fr_[i].totalCapacityBytes();
			totcap_ += maplist_rl_[i].totalCapacityBytes();
			totcap_ += maplist_rr_[i].totalCapacityBytes();
			for (size_t j = 0; j < rdlen; j++)
			{
				maplist_fl_[i][j].clear();
				maplist_fr_[i][j].clear();
				maplist_rl_[i][j].clear();
				maplist_rr_[i][j].clear();
				totcap_ += maplist_fl_[i][j].totalCapacityBytes();
				totcap_ += maplist_fr_[i][j].totalCapacityBytes();
				totcap_ += maplist_rl_[i][j].totalCapacityBytes();
				totcap_ += maplist_rr_[i][j].totalCapacityBytes();
			}
		}
		inited_ = true;
	}
	bool inited() const
	{
		return inited_;
	}
	bool check(
		bool fw,
		bool l2r,
		TReadOff al5pi,
		TReadOff al5pf,
		size_t rflen,
		TIndexOffU topf,
		TIndexOffU botf,
		TScore pen)
	{
		assert(inited_);
		assert(topf > 0 || botf > 0);
		DescentRedundancyKey k(al5pf, rflen, topf, botf);
		size_t i = std::numeric_limits<size_t>::max();
		size_t mask = topf & PART_MASK;
		EMap<DescentRedundancyKey, TScore> &map =
			(fw ? (l2r ? maplist_fl_[mask][al5pi] : maplist_fr_[mask][al5pi]) : (l2r ? maplist_rl_[mask][al5pi] : maplist_rr_[mask][al5pi]));
		size_t key = (topf & 255) | ((botf & 255) << 8);
		if (bits_.test(key) && map.containsEx(k, i))
		{
			assert_lt(i, map.size());
			assert_geq(pen, map[i].second);
			return false;
		}
		assert(!map.containsEx(k, i));
		size_t oldsz = map.totalSizeBytes();
		size_t oldcap = map.totalCapacityBytes();
		map.insert(make_pair(k, pen));
		bits_.set(key);
		totsz_ += (map.totalSizeBytes() - oldsz);
		totcap_ += (map.totalCapacityBytes() - oldcap);
		return true;
	}
	bool contains(
		bool fw,
		bool l2r,
		TReadOff al5pi,
		TReadOff al5pf,
		size_t rflen,
		TIndexOffU topf,
		TIndexOffU botf,
		TScore pen)
	{
		assert(inited_);
		size_t key = (topf & 255) | ((botf & 255) << 8);
		if (!bits_.test(key))
		{
			return false;
		}
		DescentRedundancyKey k(al5pf, rflen, topf, botf);
		size_t mask = topf & PART_MASK;
		EMap<DescentRedundancyKey, TScore> &map =
			(fw ? (l2r ? maplist_fl_[mask][al5pi] : maplist_fr_[mask][al5pi]) : (l2r ? maplist_rl_[mask][al5pi] : maplist_rr_[mask][al5pi]));
		return map.contains(k);
	}
	size_t totalSizeBytes() const
	{
		return totsz_;
	}
	size_t totalCapacityBytes() const
	{
		return totcap_;
	}
protected:
	bool inited_;
	size_t totsz_;
	size_t totcap_;
	ELList<EMap<DescentRedundancyKey, TScore>, NPARTS, 100> maplist_fl_;
	ELList<EMap<DescentRedundancyKey, TScore>, NPARTS, 100> maplist_rl_;
	ELList<EMap<DescentRedundancyKey, TScore>, NPARTS, 100> maplist_fr_;
	ELList<EMap<DescentRedundancyKey, TScore>, NPARTS, 100> maplist_rr_;
	EBitList<128> bits_;
};
struct DescentRoot
{
	DescentRoot() { reset(); }
	DescentRoot(
		size_t off5p_,
		bool l2r_,
		bool fw_,
		size_t landing_,
		size_t len,
		float pri_)
	{
		init(off5p_, l2r_, fw_, landing_, len, pri_);
	}
	void init(
		size_t off5p_,
		bool l2r_,
		bool fw_,
		size_t landing_,
		size_t len,
		float pri_)
	{
		off5p = off5p_;
		l2r = l2r_;
		fw = fw_;
		landing = landing_;
		pri = pri_;
		assert_lt(off5p, len);
	}
	void reset()
	{
		off5p = std::numeric_limits<size_t>::max();
	}
	bool inited() const
	{
		return off5p == std::numeric_limits<size_t>::max();
	}
	bool operator==(const DescentRoot &o) const
	{
		return pri == o.pri && off5p == o.off5p && l2r == o.l2r &&
			   fw == o.fw && landing == o.landing;
	}
	bool operator<(const DescentRoot &o) const
	{
		if (pri > o.pri)
			return true;
		if (pri < o.pri)
			return false;
		if (off5p < o.off5p)
			return true;
		if (off5p > o.off5p)
			return false;
		if (fw != o.fw)
			return fw;
		if (l2r != o.l2r)
			return l2r;
		if (landing < o.landing)
			return true;
		if (landing > o.landing)
			return false;
		return false;
	}
	bool operator<=(const DescentRoot &o) const
	{
		return (*this) < o || (*this) == o;
	}
	TReadOff off5p;
	bool l2r;
	bool fw;
	size_t landing;
	float pri;
};
struct DescentPosFlags
{
	DescentPosFlags() { reset(); }
	void reset()
	{
		mm_a = mm_c = mm_g = mm_t = rdg_a = rdg_c = rdg_g = rdg_t = rfg = 1;
		reserved = 0;
	}
	bool exhausted() const
	{
		return ((uint16_t *)this)[0] == 0;
	}
	bool mmExplore(int c)
	{
		assert_range(0, 3, c);
		if (c == 0)
		{
			return mm_a;
		}
		else if (c == 1)
		{
			return mm_c;
		}
		else if (c == 2)
		{
			return mm_g;
		}
		else
		{
			return mm_t;
		}
	}
	bool mmSet(int c)
	{
		assert_range(0, 3, c);
		if (c == 0)
		{
			bool ret = mm_a;
			mm_a = 0;
			return ret;
		}
		else if (c == 1)
		{
			bool ret = mm_c;
			mm_c = 0;
			return ret;
		}
		else if (c == 2)
		{
			bool ret = mm_g;
			mm_g = 0;
			return ret;
		}
		else
		{
			bool ret = mm_t;
			mm_t = 0;
			return ret;
		}
	}
	bool rdgExplore(int c)
	{
		assert_range(0, 3, c);
		if (c == 0)
		{
			return rdg_a;
		}
		else if (c == 1)
		{
			return rdg_c;
		}
		else if (c == 2)
		{
			return rdg_g;
		}
		else
		{
			return rdg_t;
		}
	}
	bool rdgSet(int c)
	{
		assert_range(0, 3, c);
		if (c == 0)
		{
			bool ret = rdg_a;
			rdg_a = 0;
			return ret;
		}
		else if (c == 1)
		{
			bool ret = rdg_c;
			rdg_c = 0;
			return ret;
		}
		else if (c == 2)
		{
			bool ret = rdg_g;
			rdg_g = 0;
			return ret;
		}
		else
		{
			bool ret = rdg_t;
			rdg_t = 0;
			return ret;
		}
	}
	bool rfgExplore()
	{
		return rfg;
	}
	bool rfgSet()
	{
		bool ret = rfg;
		rfg = 0;
		return ret;
	}
	uint16_t mm_a : 1;
	uint16_t mm_c : 1;
	uint16_t mm_g : 1;
	uint16_t mm_t : 1;
	uint16_t rdg_a : 1;
	uint16_t rdg_c : 1;
	uint16_t rdg_g : 1;
	uint16_t rdg_t : 1;
	uint16_t rfg : 1;
	uint16_t reserved : 7;
};
struct DescentPos
{
	void reset()
	{
		topf[0] = topf[1] = topf[2] = topf[3] = 0;
		botf[0] = botf[1] = botf[2] = botf[3] = 0;
		topb[0] = topb[1] = topb[2] = topb[3] = 0;
		botb[0] = botb[1] = botb[2] = botb[3] = 0;
		c = -1;
		flags.reset();
	}
	bool inited() const
	{
		return c >= 0;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_range(0, 3, (int)c);
		return true;
	}
#endif
	TIndexOffU topf[4];
	TIndexOffU botf[4];
	TIndexOffU topb[4];
	TIndexOffU botb[4];
	char c;
	DescentPosFlags flags;
};
struct DescentEdge
{
	DescentEdge() { reset(); }
	DescentEdge(
		Edit e_,
		TReadOff off5p_,
		DescentPriority pri_,
		size_t posFlag_,
		TReadOff nex_
#ifndef NDEBUG
		,
		size_t d_,
		TIndexOffU topf_,
		TIndexOffU botf_,
		TIndexOffU topb_,
		TIndexOffU botb_
#endif
	)
	{
		init(e_, off5p_, pri_, posFlag_
#ifndef NDEBUG
			 ,
			 d_, topf_, botf_, topb_, botb_
#endif
		);
	}
	bool inited() const { return e.inited(); }
	void reset() { e.reset(); }
	void init(
		Edit e_,
		TReadOff off5p_,
		DescentPriority pri_,
		size_t posFlag_
#ifndef NDEBUG
		,
		size_t d_,
		TIndexOffU topf_,
		TIndexOffU botf_,
		TIndexOffU topb_,
		TIndexOffU botb_
#endif
	)
	{
		e = e_;
		off5p = off5p_;
		pri = pri_;
		posFlag = posFlag_;
#ifndef NDEBUG
		d = d_;
		topf = topf_;
		botf = botf_;
		topb = topb_;
		botb = botb_;
#endif
	}
	void updateFlags(EFactory<DescentPos> &pf)
	{
		if (inited())
		{
			if (e.isReadGap())
			{
				assert_neq('-', e.chr);
				pf[posFlag].flags.rdgSet(asc2dna[e.chr]);
			}
			else if (e.isRefGap())
			{
				pf[posFlag].flags.rfgSet();
			}
			else
			{
				assert_neq('-', e.chr);
				pf[posFlag].flags.mmSet(asc2dna[e.chr]);
			}
		}
	}
	bool operator<(const DescentEdge &o) const
	{
		if (inited() && !o.inited())
		{
			return true;
		}
		else if (!inited())
		{
			return false;
		}
		return pri < o.pri;
	}
	DescentPriority pri;
	TReadOff nex;
	size_t posFlag;
#ifndef NDEBUG
	size_t d;
	TIndexOffU topf, botf, topb, botb;
#endif
	Edit e;
	TReadOff off5p;
};
class DescentOutgoing
{
public:
	DescentEdge rotate()
	{
		DescentEdge tmp = best1;
		assert(!(best2 < tmp));
		best1 = best2;
		assert(!(best3 < best2));
		best2 = best3;
		assert(!(best4 < best3));
		best3 = best4;
		assert(!(best5 < best4));
		best4 = best5;
		best5.reset();
		return tmp;
	}
	void update(DescentEdge e)
	{
		if (!best1.inited())
		{
			best1 = e;
		}
		else if (e < best1)
		{
			best5 = best4;
			best4 = best3;
			best3 = best2;
			best2 = best1;
			best1 = e;
		}
		else if (!best2.inited())
		{
			best2 = e;
		}
		else if (e < best2)
		{
			best5 = best4;
			best4 = best3;
			best3 = best2;
			best2 = e;
		}
		else if (!best3.inited())
		{
			best3 = e;
		}
		else if (e < best3)
		{
			best5 = best4;
			best4 = best3;
			best3 = e;
		}
		else if (!best4.inited())
		{
			best4 = e;
		}
		else if (e < best4)
		{
			best5 = best4;
			best4 = e;
		}
		else if (!best5.inited() || e < best5)
		{
			best5 = e;
		}
	}
	void clear()
	{
		best1.reset();
		best2.reset();
		best3.reset();
		best4.reset();
		best5.reset();
	}
	bool empty() const
	{
		return !best1.inited();
	}
	DescentPriority bestPri() const
	{
		assert(!empty());
		return best1.pri;
	}
	DescentEdge best1;
	DescentEdge best2;
	DescentEdge best3;
	DescentEdge best4;
	DescentEdge best5;
};
class DescentAlignmentSink;
class Descent
{
public:
	Descent() { reset(); }
	bool init(
		const Read &q,
		TRootId rid,
		const Scoring &sc,
		TAlScore minsc,
		TAlScore maxpen,
		TReadOff al5pi,
		TReadOff al5pf,
		TIndexOffU topf,
		TIndexOffU botf,
		TIndexOffU topb,
		TIndexOffU botb,
		bool l2r,
		size_t descid,
		TDescentId parent,
		TScore pen,
		const Edit &e,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		DescentRedundancyChecker &re,
		EFactory<Descent> &df,
		EFactory<DescentPos> &pf,
		const EList<DescentRoot> &rs,
		const EList<DescentConfig> &cs,
		EHeap<TDescentPair> &heap,
		DescentAlignmentSink &alsink,
		DescentMetrics &met,
		PerReadMetrics &prm);
	bool init(
		const Read &q,
		TRootId rid,
		const Scoring &sc,
		TAlScore minsc,
		TAlScore maxpen,
		size_t descid,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		DescentRedundancyChecker &re,
		EFactory<Descent> &df,
		EFactory<DescentPos> &pf,
		const EList<DescentRoot> &rs,
		const EList<DescentConfig> &cs,
		EHeap<TDescentPair> &heap,
		DescentAlignmentSink &alsink,
		DescentMetrics &met,
		PerReadMetrics &prm);
	bool inited() const
	{
		return descid_ != std::numeric_limits<size_t>::max();
	}
	void reset()
	{
		lastRecalc_ = true;
		descid_ = std::numeric_limits<size_t>::max();
	}
	bool root() const
	{
		return parent_ == std::numeric_limits<TDescentId>::max();
	}
	const Edit &edit() const
	{
		return edit_;
	}
	TDescentId parent() const
	{
		return parent_;
	}
	void followBestOutgoing(
		const Read &q,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		const Scoring &sc,
		TAlScore minsc,
		TAlScore maxpen,
		DescentRedundancyChecker &re,
		EFactory<Descent> &df,
		EFactory<DescentPos> &pf,
		const EList<DescentRoot> &rs,
		const EList<DescentConfig> &cs,
		EHeap<TDescentPair> &heap,
		DescentAlignmentSink &alsink,
		DescentMetrics &met,
		PerReadMetrics &prm);
	bool empty() const { return lastRecalc_ && out_.empty(); }
#ifndef NDEBUG
	bool repOk(const Read *q) const
	{
		assert(!root() || !edit_.inited());
		assert_eq(botf_ - topf_, botb_ - topb_);
		if (q != NULL)
		{
			assert_leq(len_, q->length());
		}
		return true;
	}
#endif
	size_t al5pi() const
	{
		return al5pi_;
	}
	size_t al5pf() const { return al5pf_; }
	bool l2r() const { return l2r_; }
	void print(
		std::ostream *os,
		const char *prefix,
		const Read &q,
		size_t trimLf,
		size_t trimRg,
		bool fw,
		const EList<Edit> &edits,
		size_t ei,
		size_t en,
		BTDnaString &rf) const;
	void collectEdits(
		EList<Edit> &edits,
		const Edit *e,
		EFactory<Descent> &df)
	{
		size_t nuninited = 0;
		size_t ei = edits.size();
		size_t en = 0;
		if (e != NULL && e->inited())
		{
			edits.push_back(*e);
			en++;
		}
		size_t cur = descid_;
		while (cur != std::numeric_limits<TDescentId>::max())
		{
			if (!df[cur].edit().inited())
			{
				nuninited++;
				assert_leq(nuninited, 2);
			}
			else
			{
				edits.push_back(df[cur].edit());
				en++;
			}
			cur = df[cur].parent();
		}
		edits.sortPortion(ei, en);
	}
	TIndexOffU topf() const { return topf_; }
	TIndexOffU botf() const { return botf_; }
protected:
	bool bounce(
		const Read &q,
		TIndexOffU topf,
		TIndexOffU botf,
		TIndexOffU topb,
		TIndexOffU botb,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		const Scoring &sc,
		TAlScore minsc,
		TAlScore maxpen,
		DescentRedundancyChecker &re,
		EFactory<Descent> &df,
		EFactory<DescentPos> &pf,
		const EList<DescentRoot> &rs,
		const EList<DescentConfig> &cs,
		EHeap<TDescentPair> &heap,
		DescentAlignmentSink &alsink,
		DescentMetrics &met,
		PerReadMetrics &prm);
	void nextLocsBi(
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		SideLocus &tloc,
		SideLocus &bloc,
		TIndexOffU topf,
		TIndexOffU botf,
		TIndexOffU topb,
		TIndexOffU botb);
	bool followMatches(
		const Read &q,
		const Scoring &sc,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		DescentRedundancyChecker &re,
		EFactory<Descent> &df,
		EFactory<DescentPos> &pf,
		const EList<DescentRoot> &rs,
		const EList<DescentConfig> &cs,
		EHeap<TDescentPair> &heap,
		DescentAlignmentSink &alsink,
		DescentMetrics &met,
		PerReadMetrics &prm,
		bool &branches,
		bool &hitEnd,
		bool &done,
		TReadOff &off5p_i,
		TIndexOffU &topf_bounce,
		TIndexOffU &botf_bounce,
		TIndexOffU &topb_bounce,
		TIndexOffU &botb_bounce);
	size_t recalcOutgoing(
		const Read &q,
		const Scoring &sc,
		TAlScore minsc,
		TAlScore maxpen,
		DescentRedundancyChecker &re,
		EFactory<DescentPos> &pf,
		const EList<DescentRoot> &rs,
		const EList<DescentConfig> &cs,
		PerReadMetrics &prm);
	TRootId rid_;
	TReadOff al5pi_;
	TReadOff al5pf_;
	bool l2r_;
	int gapadd_;
	TReadOff off5p_i_;
	TIndexOffU topf_, botf_;
	TIndexOffU topb_, botb_;
	size_t descid_;
	TDescentId parent_;
	TScore pen_;
	size_t posid_;
	size_t len_;
	DescentOutgoing out_;
	Edit edit_;
	bool lastRecalc_;
};
struct DescentAlignment
{
	DescentAlignment() { reset(); }
	void reset()
	{
		topf = botf = 0;
		pen = 0;
		fw = false;
		ei = en = 0;
	}
	void init(
		TScore pen_,
		bool fw_,
		TIndexOffU topf_,
		TIndexOffU botf_,
		size_t ei_,
		size_t en_)
	{
		assert_gt(botf_, topf_);
		pen = pen_;
		fw = fw_;
		topf = topf_;
		botf = botf_;
		ei = ei_;
		en = en_;
	}
	bool inited() const
	{
		return botf > topf;
	}
	bool perfect() const
	{
		return pen == 0;
	}
	size_t size() const
	{
		return botf - topf;
	}
	TScore pen;
	bool fw;
	TIndexOffU topf;
	TIndexOffU botf;
	size_t ei;
	size_t en;
};
struct DescentPartialResolvedAlignment
{
	DescentPartialResolvedAlignment() { reset(); }
	void reset()
	{
		topf = botf = 0;
		pen = 0;
		fw = false;
		ei = en = 0;
		refcoord.reset();
	}
	void init(
		TScore pen_,
		bool fw_,
		TIndexOffU topf_,
		TIndexOffU botf_,
		size_t ei_,
		size_t en_,
		const Coord &refcoord_)
	{
		assert_gt(botf_, topf_);
		pen = pen_;
		fw = fw_;
		topf = topf_;
		botf = botf_;
		ei = ei_;
		en = en_;
		refcoord = refcoord_;
	}
	bool inited() const
	{
		return botf > topf;
	}
	size_t size() const
	{
		return botf - topf;
	}
	TScore pen;
	bool fw;
	TIndexOffU topf;
	TIndexOffU botf;
	size_t ei;
	size_t en;
	Coord refcoord;
};
class DescentAlignmentSink
{
public:
	bool reportAlignment(
		const Read &q,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		TIndexOffU topf,
		TIndexOffU botf,
		TIndexOffU topb,
		TIndexOffU botb,
		TDescentId id,
		TRootId rid,
		const Edit &e,
		TScore pen,
		EFactory<Descent> &df,
		EFactory<DescentPos> &pf,
		const EList<DescentRoot> &rs,
		const EList<DescentConfig> &cs);
	void reset()
	{
		edits_.clear();
		als_.clear();
		lhs_.clear();
		rhs_.clear();
		nelt_ = 0;
		bestPen_ = worstPen_ = std::numeric_limits<TAlScore>::max();
	}
	size_t totalSizeBytes() const
	{
		return edits_.totalSizeBytes() +
			   als_.totalSizeBytes() +
			   lhs_.totalSizeBytes() +
			   rhs_.totalSizeBytes() +
			   sizeof(size_t);
	}
	size_t totalCapacityBytes() const
	{
		return edits_.totalCapacityBytes() +
			   als_.totalCapacityBytes() +
			   lhs_.totalCapacityBytes() +
			   rhs_.totalCapacityBytes() +
			   sizeof(size_t);
	}
	size_t nrange() const
	{
		return als_.size();
	}
	size_t nelt() const
	{
		return nelt_;
	}
	void elt(size_t i, DescentAlignment &al, size_t &ri, size_t &off) const
	{
		assert_lt(i, nelt());
		for (size_t j = 0; j < als_.size(); j++)
		{
			if (i < als_[j].size())
			{
				al = als_[j];
				ri = j;
				off = i;
				return;
			}
			i -= als_[j].size();
		}
		assert(false);
	}
	const DescentAlignment &operator[](size_t i) const
	{
		return als_[i];
	}
	bool stratumDone(TAlScore bestPen) const
	{
		if (nelt_ > 0 && bestPen > worstPen_)
		{
			return true;
		}
		return false;
	}
	void advanceStratum()
	{
		assert_gt(nelt_, 0);
		edits_.clear();
		als_.clear();
		nelt_ = 0;
		bestPen_ = worstPen_ = std::numeric_limits<TAlScore>::max();
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_geq(nelt_, als_.size());
		for (size_t i = 1; i < als_.size(); i++)
		{
			assert_geq(als_[i].pen, als_[i - 1].pen);
		}
		assert(bestPen_ == std::numeric_limits<TAlScore>::max() || worstPen_ >= bestPen_);
		return true;
	}
#endif
	TAlScore bestPenalty() const
	{
		return bestPen_;
	}
	TAlScore worstPenalty() const { return worstPen_; }
	size_t editsSize() const { return edits_.size(); }
	size_t alsSize() const { return als_.size(); }
	size_t lhsSize() const { return lhs_.size(); }
	size_t rhsSize() const { return rhs_.size(); }
	const EList<Edit> &edits() const { return edits_; }
protected:
	EList<Edit> edits_;
	EList<DescentAlignment> als_;
	ESet<Triple<TIndexOffU, TIndexOffU, size_t>> lhs_;
	ESet<Triple<TIndexOffU, TIndexOffU, size_t>> rhs_;
	size_t nelt_;
	TAlScore bestPen_;
	TAlScore worstPen_;
#ifndef NDEBUG
	BTDnaString tmprfdnastr_;
#endif
};
class DescentPartialResolvedAlignmentSink
{
public:
	void reset()
	{
		edits_.clear();
		als_.clear();
	}
	size_t totalSizeBytes() const
	{
		return edits_.totalSizeBytes() +
			   als_.totalSizeBytes() +
			   sizeof(size_t);
	}
	size_t totalCapacityBytes() const
	{
		return edits_.totalCapacityBytes() +
			   als_.totalCapacityBytes() +
			   sizeof(size_t);
	}
	const DescentPartialResolvedAlignment &operator[](size_t i) const
	{
		return als_[i];
	}
	size_t editsSize() const { return edits_.size(); }
	size_t alsSize() const { return als_.size(); }
	const EList<Edit> &edits() const { return edits_; }
protected:
	EList<Edit> edits_;
	EList<DescentPartialResolvedAlignment> als_;
};
class DescentRootSelector
{
public:
	virtual ~DescentRootSelector() {}
	virtual void select(
		const Read &q,
		const Read *qo,
		bool nofw,
		bool norc,
		EList<DescentConfig> &confs,
		EList<DescentRoot> &roots) = 0;
};
struct DescentStoppingConditions
{
	DescentStoppingConditions() { reset(); }
	DescentStoppingConditions(
		size_t totsz_,
		size_t nfound_,
		bool stra_,
		size_t nbwop_)
	{
		init(totsz_, nfound_, stra_, nbwop_);
	}
	void reset()
	{
		totsz = nfound = nbwop = std::numeric_limits<size_t>::max();
		stra = false;
		assert(!inited());
	}
	void init(
		size_t totsz_,
		size_t nfound_,
		bool stra_,
		size_t nbwop_)
	{
		totsz = totsz_;
		nfound = nfound_;
		stra = stra_;
		nbwop = nbwop_;
		assert(inited());
	}
	bool inited() const
	{
		return totsz != std::numeric_limits<size_t>::max();
	}
	size_t totsz;
	size_t nfound;
	bool stra;
	size_t nbwop;
};
enum
{
	DESCENT_DRIVER_ALN = 1,
	DESCENT_DRIVER_STRATA = 2,
	DESCENT_DRIVER_MEM = 4,
	DESCENT_DRIVER_BWOPS = 8,
	DESCENT_DRIVER_DONE = 16
};
class DescentDriver
{
public:
	DescentDriver(bool veryVerbose = false) : veryVerbose_(veryVerbose)
	{
		reset();
	}
	void initRead(
		const Read &q,
		bool nofw,
		bool norc,
		TAlScore minsc,
		TAlScore maxpen,
		const Read *qmate = NULL,
		DescentRootSelector *sel = NULL)
	{
		reset();
		q_ = q;
		minsc_ = minsc;
		maxpen_ = maxpen;
		if (sel != NULL)
		{
			sel->select(
				q_,
				qmate,
				nofw,
				norc,
				confs_,
				roots_);
		}
		re_.init(q.length());
	}
	void addRoot(
		const DescentConfig &conf,
		TReadOff off,
		bool l2r,
		bool fw,
		size_t landing,
		float pri)
	{
		confs_.push_back(conf);
		assert_lt(off, q_.length());
		if (l2r && off == q_.length() - 1)
		{
			l2r = !l2r;
		}
		else if (!l2r && off == 0)
		{
			l2r = !l2r;
		}
		roots_.push_back(DescentRoot(off, l2r, fw, landing, q_.length(), pri));
	}
	void clearRoots()
	{
		confs_.clear();
		roots_.clear();
	}
	void printRoots(std::ostream &os)
	{
		std::ostringstream fwstr, rcstr;
		fwstr << q_.patFw << std::endl
			  << q_.qual << std::endl;
		rcstr << q_.patRc << std::endl
			  << q_.qualRev << std::endl;
		for (size_t i = 0; i < roots_.size(); i++)
		{
			if (roots_[i].fw)
			{
				for (size_t j = 0; j < roots_[i].off5p; j++)
				{
					fwstr << " ";
				}
				fwstr << (roots_[i].l2r ? ">" : "<");
				fwstr << " " << i << ":";
				fwstr << roots_[i].pri;
				fwstr << "\n";
			}
			else
			{
				size_t off = q_.length() - roots_[i].off5p - 1;
				for (size_t j = 0; j < off; j++)
				{
					rcstr << " ";
				}
				rcstr << (roots_[i].l2r ? ">" : "<");
				rcstr << " " << i << ":";
				rcstr << roots_[i].pri;
				rcstr << "\n";
			}
		}
		os << fwstr.str() << rcstr.str();
	}
	void resetRead()
	{
		df_.clear();
		assert_leq(df_.totalSizeBytes(), 100);
		pf_.clear();
		assert_leq(pf_.totalSizeBytes(), 100);
		heap_.clear();
		assert_leq(heap_.totalSizeBytes(), 100);
		roots_.clear();
		assert_leq(roots_.totalSizeBytes(), 100);
		confs_.clear();
		assert_leq(confs_.totalSizeBytes(), 100);
		alsink_.reset();
		assert_leq(alsink_.totalSizeBytes(), 100);
		re_.reset();
		assert_leq(re_.totalSizeBytes(), 100);
		rootsInited_ = 0;
		curPen_ = 0;
	}
	void reset()
	{
		resetRead();
	}
	void go(
		const Scoring &sc,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		DescentMetrics &met,
		PerReadMetrics &prm);
	int advance(
		const DescentStoppingConditions &stopc,
		const Scoring &sc,
		const Ebwt &ebwtFw,
		const Ebwt &ebwtBw,
		DescentMetrics &met,
		PerReadMetrics &prm);
#ifndef NDEBUG
	bool repOk() const
	{
		return true;
	}
#endif
	size_t numAlignments() const
	{
		return alsink_.nelt();
	}
	const DescentAlignmentSink &sink() const
	{
		return alsink_;
	}
	DescentAlignmentSink &sink()
	{
		return alsink_;
	}
	size_t totalSizeBytes() const
	{
		return df_.totalSizeBytes() +
			   pf_.totalSizeBytes() +
			   heap_.totalSizeBytes() +
			   roots_.totalSizeBytes() +
			   confs_.totalSizeBytes() +
			   alsink_.totalSizeBytes() +
			   re_.totalSizeBytes();
	}
	size_t totalCapacityBytes() const
	{
		return df_.totalCapacityBytes() +
			   pf_.totalCapacityBytes() +
			   heap_.totalCapacityBytes() +
			   roots_.totalCapacityBytes() +
			   confs_.totalCapacityBytes() +
			   alsink_.totalCapacityBytes() +
			   re_.totalCapacityBytes();
	}
	const Read &query() const
	{
		return q_;
	}
	TAlScore minScore() const
	{
		return minsc_;
	}
	const EList<DescentRoot> &roots() { return roots_; }
	void nextPartial()
	{
	}
protected:
	Read q_;
	TAlScore minsc_;
	TAlScore maxpen_;
	EFactory<Descent> df_;
	EFactory<DescentPos> pf_;
	EList<DescentRoot> roots_;
	EList<DescentConfig> confs_;
	size_t rootsInited_;
	EHeap<TDescentPair> heap_;
	DescentAlignmentSink alsink_;
	DescentRedundancyChecker re_;
	TAlScore curPen_;
	bool veryVerbose_;
	EList<Edit> tmpedit_;
	BTDnaString tmprfdnastr_;
};
class DescentAlignmentSelector
{
public:
	DescentAlignmentSelector() : gwstate_(GW_CAT) { reset(); }
	void init(
		const Read &q,
		const DescentAlignmentSink &sink,
		const Ebwt &ebwtFw,
		const BitPairReference &ref,
		RandomSource &rnd,
		WalkMetrics &met)
	{
		rnd_.init(
			sink.nelt(),
			true);
		offs_.resize(sink.nelt());
		offs_.fill(std::numeric_limits<TIndexOffU>::max());
		sas_.resize(sink.nrange());
		gws_.resize(sink.nrange());
		size_t ei = 0;
		for (size_t i = 0; i < sas_.size(); i++)
		{
			size_t en = sink[i].botf - sink[i].topf;
			sas_[i].init(sink[i].topf, EListSlice<TIndexOffU, 16>(offs_, ei, en));
			gws_[i].init(ebwtFw, ref, sas_[i], rnd, met);
			ei += en;
		}
	}
	void reset()
	{
		rnd_.reset();
	}
	bool inited() const
	{
		return rnd_.size() > 0;
	}
	bool next(
		const DescentDriver &dr,
		const Ebwt &ebwtFw,
		const BitPairReference &ref,
		RandomSource &rnd,
		AlnRes &rs,
		WalkMetrics &met,
		PerReadMetrics &prm)
	{
		size_t ri = (size_t)rnd_.next(rnd);
		size_t off = 0;
		DescentAlignment al;
		size_t rangei = 0;
		dr.sink().elt(ri, al, rangei, off);
		assert_lt(off, al.size());
		Coord refcoord;
		WalkResult wr;
		TIndexOffU tidx = 0, toff = 0, tlen = 0;
		gws_[rangei].advanceElement(
			(TIndexOffU)off,
			ebwtFw,
			ref,
			sas_[rangei],
			gwstate_,
			wr,
			met,
			prm);
		assert_neq(OFF_MASK, wr.toff);
		bool straddled = false;
		ebwtFw.joinedToTextOff(
			wr.elt.len,
			wr.toff,
			tidx,
			toff,
			tlen,
			true,
			straddled);
		if (tidx == OFF_MASK)
		{
			return false;
		}
		refcoord.init(tidx, (int64_t)toff, dr.sink()[rangei].fw);
		const EList<Edit> &edits = dr.sink().edits();
		size_t ns = 0, ngap = 0, nrefn = 0;
		for (size_t i = al.ei; i < al.ei + al.en; i++)
		{
			if (edits[i].qchr == 'N' || edits[i].chr == 'N')
				ns++;
			if (edits[i].chr == 'N')
				nrefn++;
			if (edits[i].isGap())
				ngap++;
		}
		AlnScore asc(
			-dr.sink().bestPenalty(),
			dr.query().length() - edits.size(),
			(int)edits.size(),
			ns,
			ngap);
		rs.init(
			dr.query().length(),
			asc,
			&dr.sink().edits(),
			al.ei,
			al.en,
			NULL,
			0,
			0,
			refcoord,
			tlen,
			-1,
			-1,
			-1,
			dr.minScore(),
			false,
			0,
			0,
			false,
			0,
			0);
		rs.setRefNs(nrefn);
		return true;
	}
	bool done() const
	{
		return rnd_.done();
	}
	size_t totalSizeBytes() const
	{
		return rnd_.totalSizeBytes() +
			   offs_.totalSizeBytes() +
			   sas_.totalSizeBytes() +
			   gws_.totalSizeBytes();
	}
	size_t totalCapacityBytes() const
	{
		return rnd_.totalCapacityBytes() +
			   offs_.totalCapacityBytes() +
			   sas_.totalCapacityBytes() +
			   gws_.totalCapacityBytes();
	}
protected:
	Random1toN rnd_;
	EList<TIndexOffU, 16> offs_;
	EList<SARangeWithOffs<EListSlice<TIndexOffU, 16>>> sas_;
	EList<GroupWalk2S<EListSlice<TIndexOffU, 16>, 16>> gws_;
	GroupWalkState gwstate_;
};
class DescentPartialAlignmentSelector
{
public:
	static const size_t BIG_RANGE = 5;
	static const size_t NRANGE_AT_A_TIME = 3;
	DescentPartialAlignmentSelector() : gwstate_(GW_CAT) { reset(); }
	void init(
		const Read &q,
		const EHeap<TDescentPair> &heap,
		EFactory<Descent> &df,
		EFactory<DescentPos> &pf,
		TAlScore depthBonus,
		size_t nbatch,
		const Ebwt &ebwtFw,
		const BitPairReference &ref,
		RandomSource &rnd,
		WalkMetrics &met)
	{
		if (depthBonus > 0)
		{
			heap_.clear();
			for (size_t i = 0; i < heap.size(); i++)
			{
				TDescentPair p = heap[i];
				p.first.pen += depthBonus * p.first.depth;
				heap_.insert(p);
			}
		}
		else
			heap_ = heap;
		assert(!heap_.empty());
		nextRanges(df, pf, ebwtFw, ref, rnd, met);
		assert(!rangeExhausted());
	}
	bool empty() const
	{
		return heap_.empty();
	}
	void reset()
	{
		heap_.clear();
		offs_.clear();
		sas_.clear();
		gws_.clear();
	}
	bool inited() const
	{
		return !heap_.empty();
	}
	bool nextPartial(
		const DescentDriver &dr,
		const Ebwt &ebwtFw,
		const BitPairReference &ref,
		RandomSource &rnd,
		AlnRes &rs,
		WalkMetrics &met,
		PerReadMetrics &prm)
	{
		size_t ri = (size_t)rnd_.next(rnd);
		size_t off = 0;
		DescentAlignment al;
		size_t rangei = 0;
		dr.sink().elt(ri, al, rangei, off);
		assert_lt(off, al.size());
		Coord refcoord;
		WalkResult wr;
		TIndexOffU tidx = 0, toff = 0, tlen = 0;
		gws_[rangei].advanceElement(
			(TIndexOffU)off,
			ebwtFw,
			ref,
			sas_[rangei],
			gwstate_,
			wr,
			met,
			prm);
		assert_neq(OFF_MASK, wr.toff);
		bool straddled = false;
		ebwtFw.joinedToTextOff(
			wr.elt.len,
			wr.toff,
			tidx,
			toff,
			tlen,
			true,
			straddled);
		if (tidx == OFF_MASK)
		{
			return false;
		}
		refcoord.init(tidx, (int64_t)toff, dr.sink()[rangei].fw);
		const EList<Edit> &edits = dr.sink().edits();
		size_t ns = 0, ngap = 0, nrefn = 0;
		for (size_t i = al.ei; i < al.ei + al.en; i++)
		{
			if (edits[i].qchr == 'N' || edits[i].chr == 'N')
				ns++;
			if (edits[i].chr == 'N')
				nrefn++;
			if (edits[i].isGap())
				ngap++;
		}
		return true;
	}
	size_t totalSizeBytes() const
	{
		return heap_.totalSizeBytes() +
			   rnd_.totalSizeBytes() +
			   offs_.totalSizeBytes() +
			   sas_.totalSizeBytes() +
			   gws_.totalSizeBytes();
	}
	size_t totalCapacityBytes() const
	{
		return heap_.totalCapacityBytes() +
			   rnd_.totalCapacityBytes() +
			   offs_.totalCapacityBytes() +
			   sas_.totalCapacityBytes() +
			   gws_.totalCapacityBytes();
	}
protected:
	bool rangeExhausted() const
	{
		return rnd_.left() == 0;
	}
	void nextRanges(
		EFactory<Descent> &df,
		EFactory<DescentPos> &pf,
		const Ebwt &ebwtFw,
		const BitPairReference &ref,
		RandomSource &rnd,
		WalkMetrics &met)
	{
		assert(!heap_.empty());
		TDescentPair p = heap_.pop();
		TIndexOffU topf = df[p.second].topf(), botf = df[p.second].botf();
		assert_gt(botf, topf);
		offs_.resize(botf - topf);
		offs_.fill(std::numeric_limits<TIndexOffU>::max());
		rnd_.init(botf - topf, true);
		sas_.resize(1);
		gws_.resize(1);
		sas_[0].init(topf, EListSlice<TIndexOffU, 16>(offs_, 0, botf - topf));
		gws_[0].init(ebwtFw, ref, sas_[0], rnd, met);
	}
	DescentPartialResolvedAlignmentSink palsink_;
	EHeap<TDescentPair> heap_;
	Random1toN rnd_;
	EList<TIndexOffU, 16> offs_;
	EList<SARangeWithOffs<EListSlice<TIndexOffU, 16>>> sas_;
	EList<GroupWalk2S<EListSlice<TIndexOffU, 16>, 16>> gws_;
	GroupWalkState gwstate_;
};
#endif
