#ifndef ALIGNER_RESULT_H_
#define ALIGNER_RESULT_H_
#include <utility>
#include <limits>
#include "mem_ids.h"
#include "ref_coord.h"
#include "read.h"
#include "filebuf.h"
#include "ds.h"
#include "edit.h"
#include "limit.h"
typedef int64_t TAlScore;
#define VALID_AL_SCORE(x) ((x).score_ > MIN_I64)
#define VALID_SCORE(x) ((x) > MIN_I64)
#define INVALIDATE_SCORE(x) ((x) = MIN_I64)
class AlnScore
{
public:
	inline AlnScore()
	{
		reset();
		invalidate();
		assert(!valid());
	}
	inline AlnScore(TAlScore score, int basesAligned, int edits, TAlScore ns, TAlScore gaps)
	{
		score_ = score;
		basesAligned_ = basesAligned;
		edits_ = edits;
		ns_ = ns;
		gaps_ = gaps;
		assert(valid());
	}
	void reset()
	{
		score_ = basesAligned_ = edits_ = ns_ = gaps_ = 0;
	}
	inline static AlnScore INVALID()
	{
		AlnScore s;
		s.invalidate();
		assert(!s.valid());
		return s;
	}
	inline bool valid() const
	{
		return score_ != MIN_I64;
	}
	inline void invalidate()
	{
		score_ = MIN_I64;
		edits_ = basesAligned_ = std::numeric_limits<int>::min();
		ns_ = gaps_ = 0;
		assert(!valid());
	}
	inline void incNs(int nceil)
	{
		if (++ns_ > nceil)
		{
			invalidate();
		}
		assert_lt(ns_, 0x7fffffff);
	}
	inline bool operator>(const AlnScore &o) const
	{
		if (!VALID_AL_SCORE(o))
		{
			if (!VALID_AL_SCORE(*this))
			{
				return false;
			}
			else
			{
				return true;
			}
		}
		else if (!VALID_AL_SCORE(*this))
		{
			return false;
		}
		return score_ > o.score_;
	}
	inline AlnScore &operator=(const AlnScore &o)
	{
		gaps_ = o.gaps_;
		basesAligned_ = o.basesAligned_;
		ns_ = o.ns_;
		edits_ = o.edits_;
		score_ = o.score_;
		assert_lt(ns_, 0x7fffffff);
		return *this;
	}
	inline bool operator==(const AlnScore &o) const
	{
		return VALID_AL_SCORE(*this) && VALID_AL_SCORE(o) && score_ == o.score_;
	}
	inline bool operator!=(const AlnScore &o) const
	{
		return !(*this == o);
	}
	inline bool operator>=(const AlnScore &o) const
	{
		if (!VALID_AL_SCORE(o))
		{
			if (!VALID_AL_SCORE(*this))
			{
				return false;
			}
			else
			{
				return true;
			}
		}
		else if (!VALID_AL_SCORE(*this))
		{
			return false;
		}
		return score_ >= o.score_;
	}
	inline bool operator<(const AlnScore &o) const
	{
		return !operator>=(o);
	}
	inline AlnScore operator+(const AlnScore &o) const
	{
		if (!VALID_AL_SCORE(*this))
			return *this;
		AlnScore s;
		s.gaps_ = gaps_ + o.gaps_;
		s.basesAligned_ = basesAligned_ + o.basesAligned_;
		s.ns_ = ns_ + o.ns_;
		s.edits_ = edits_ + o.edits_;
		s.score_ = score_ + o.score_;
		assert_lt(s.ns_, 0x7fffffff);
		return s;
	}
	inline AlnScore operator+=(const AlnScore &o)
	{
		if (VALID_AL_SCORE(*this))
		{
			gaps_ += o.gaps_;
			basesAligned_ += o.basesAligned_;
			score_ += o.score_;
			edits_ += o.edits_;
			ns_ += o.ns_;
		}
		return (*this);
	}
	TAlScore score() const { return score_; }
	TAlScore penalty() const { return -score_; }
	TAlScore gaps() const { return gaps_; }
	TAlScore ns() const { return ns_; }
	int basesAligned() const { return basesAligned_; }
	int nedit() const { return edits_; }
	TAlScore score_;
	int basesAligned_;
	int edits_;
	TAlScore ns_;
	TAlScore gaps_;
};
enum
{
	ALN_FLAG_PAIR_CONCORD_MATE1 = 1,
	ALN_FLAG_PAIR_CONCORD_MATE2,
	ALN_FLAG_PAIR_DISCORD_MATE1,
	ALN_FLAG_PAIR_DISCORD_MATE2,
	ALN_FLAG_PAIR_UNPAIRED_MATE1,
	ALN_FLAG_PAIR_UNPAIRED_MATE2,
	ALN_FLAG_PAIR_UNPAIRED
};
class AlnFlags
{
public:
	AlnFlags()
	{
		init(
			ALN_FLAG_PAIR_UNPAIRED,
			false,
			false,
			false, false, false, false, false, false, false, false, false, false, false);
	}
	AlnFlags(
		int pairing,
		bool canMax,
		bool maxed,
		bool maxedPair,
		bool nfilt,
		bool scfilt,
		bool lenfilt,
		bool qcfilt,
		bool mixedMode,
		bool primary,
		bool oppAligned,
		bool oppFw,
		bool scUnMapped,
		bool xeq)
	{
		init(pairing, canMax, maxed, maxedPair, nfilt, scfilt,
			 lenfilt, qcfilt, mixedMode, primary, oppAligned,
			 oppFw, scUnMapped, xeq);
	}
	void init(
		int pairing,
		bool canMax,
		bool maxed,
		bool maxedPair,
		bool nfilt,
		bool scfilt,
		bool lenfilt,
		bool qcfilt,
		bool mixedMode,
		bool primary,
		bool oppAligned,
		bool oppFw,
		bool scUnMapped,
		bool xeq)
	{
		assert_gt(pairing, 0);
		assert_leq(pairing, ALN_FLAG_PAIR_UNPAIRED);
		pairing_ = pairing;
		canMax_ = canMax;
		maxed_ = maxed;
		maxedPair_ = maxedPair;
		nfilt_ = nfilt;
		scfilt_ = scfilt;
		lenfilt_ = lenfilt;
		qcfilt_ = qcfilt;
		mixedMode_ = mixedMode;
		primary_ = primary;
		oppAligned_ = oppAligned;
		oppFw_ = oppFw;
		scUnMapped_ = scUnMapped;
		xeq_ = xeq;
	}
	bool partOfPair() const
	{
		assert_gt(pairing_, 0);
		return pairing_ < ALN_FLAG_PAIR_UNPAIRED;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert(partOfPair() || !maxedPair_);
		return true;
	}
#endif
	bool printYF(BTString &o, bool first) const;
	void printYM(BTString &o) const;
	void printYP(BTString &o) const;
	void printYT(BTString &o) const;
	inline int pairing() const { return pairing_; }
	inline bool maxed() const { return maxed_; }
	inline bool maxedPair() const { return maxedPair_; }
	inline bool isPrimary() const
	{
		return primary_;
	}
	void setPrimary(bool primary)
	{
		primary_ = primary;
	}
	inline bool isMixedMode() const
	{
		return mixedMode_;
	}
	inline bool canMax() const
	{
		return canMax_;
	}
	bool filtered() const
	{
		return !nfilt_ || !scfilt_ || !lenfilt_ || !qcfilt_;
	}
	bool readMate1() const
	{
		return pairing_ == ALN_FLAG_PAIR_CONCORD_MATE1 ||
			   pairing_ == ALN_FLAG_PAIR_DISCORD_MATE1 ||
			   pairing_ == ALN_FLAG_PAIR_UNPAIRED_MATE1;
	}
	bool readMate2() const
	{
		return pairing_ == ALN_FLAG_PAIR_CONCORD_MATE2 ||
			   pairing_ == ALN_FLAG_PAIR_DISCORD_MATE2 ||
			   pairing_ == ALN_FLAG_PAIR_UNPAIRED_MATE2;
	}
	bool alignedConcordant() const
	{
		return pairing_ == ALN_FLAG_PAIR_CONCORD_MATE1 ||
			   pairing_ == ALN_FLAG_PAIR_CONCORD_MATE2;
	}
	bool alignedDiscordant() const
	{
		return pairing_ == ALN_FLAG_PAIR_DISCORD_MATE1 ||
			   pairing_ == ALN_FLAG_PAIR_DISCORD_MATE2;
	}
	bool alignedPaired() const
	{
		return alignedConcordant() && alignedDiscordant();
	}
	bool alignedUnpaired() const
	{
		return pairing_ == ALN_FLAG_PAIR_UNPAIRED;
	}
	bool alignedUnpairedMate() const
	{
		return pairing_ == ALN_FLAG_PAIR_UNPAIRED_MATE1 ||
			   pairing_ == ALN_FLAG_PAIR_UNPAIRED_MATE2;
	}
	bool mateAligned() const
	{
		return oppAligned_;
	}
	bool isOppFw() const
	{
		return oppFw_;
	}
	bool scUnMapped() const
	{
		return scUnMapped_;
	}
	bool xeq() const
	{
		return xeq_;
	}
protected:
	int pairing_;
	bool canMax_;
	bool maxed_;
	bool maxedPair_;
	bool nfilt_;
	bool scfilt_;
	bool lenfilt_;
	bool qcfilt_;
	bool mixedMode_;
	bool primary_;
	bool oppAligned_;
	bool oppFw_;
	bool scUnMapped_;
	bool xeq_;
};
static inline ostream &operator<<(ostream &os, const AlnScore &o)
{
	os << o.score();
	return os;
}
class BitPairReference;
enum
{
	ALN_RES_TYPE_UNPAIRED = 1,
	ALN_RES_TYPE_UNPAIRED_MATE1,
	ALN_RES_TYPE_UNPAIRED_MATE2,
	ALN_RES_TYPE_MATE1,
	ALN_RES_TYPE_MATE2
};
struct SeedAlSumm
{
	SeedAlSumm() { reset(); }
	void reset()
	{
		nonzTot = nonzFw = nonzRc = 0;
		nrangeTot = nrangeFw = nrangeRc = 0;
		neltTot = neltFw = neltRc = 0;
		minNonzRangeFw = minNonzRangeRc = 0;
		maxNonzRangeFw = maxNonzRangeRc = 0;
		minNonzEltFw = minNonzEltRc = 0;
		maxNonzEltFw = maxNonzEltRc = 0;
	}
	size_t nonzTot;
	size_t nonzFw;
	size_t nonzRc;
	size_t nrangeTot;
	size_t nrangeFw;
	size_t nrangeRc;
	size_t neltTot;
	size_t neltFw;
	size_t neltRc;
	size_t minNonzRangeFw;
	size_t minNonzRangeRc;
	size_t maxNonzRangeFw;
	size_t maxNonzRangeRc;
	size_t minNonzEltFw;
	size_t minNonzEltRc;
	size_t maxNonzEltFw;
	size_t maxNonzEltRc;
};
class StackedAln
{
public:
	StackedAln() : stackRef_(RES_CAT),
				   stackRel_(RES_CAT),
				   stackRead_(RES_CAT),
				   cigOp_(RES_CAT),
				   cigRun_(RES_CAT),
				   mdzOp_(RES_CAT),
				   mdzChr_(RES_CAT),
				   mdzRun_(RES_CAT)
	{
		reset();
	}
	void reset()
	{
		inited_ = false;
		trimLS_ = trimLH_ = trimRS_ = trimRH_ = 0;
		stackRef_.clear();
		stackRel_.clear();
		stackRead_.clear();
		cigDistMm_ = cigCalc_ = false;
		cigOp_.clear();
		cigRun_.clear();
		mdzCalc_ = false;
		mdzOp_.clear();
		mdzChr_.clear();
		mdzRun_.clear();
	}
	bool inited() const { return inited_; }
	void init(
		const BTDnaString &s,
		const EList<Edit> &ed,
		size_t trimLS,
		size_t trimLH,
		size_t trimRS,
		size_t trimRH);
	void leftAlign(bool pastMms);
	bool buildCigar(bool xeq);
	bool buildMdz();
	void writeCigar(BTString *o, char *oc) const;
	void writeMdz(BTString *o, char *oc) const;
#ifndef NDEBUG
	bool repOk() const
	{
		if (inited_)
		{
			assert_eq(stackRef_.size(), stackRead_.size());
			assert_eq(stackRef_.size(), stackRel_.size());
		}
		return true;
	}
#endif
protected:
	bool inited_;
	size_t trimLS_;
	size_t trimLH_;
	size_t trimRS_;
	size_t trimRH_;
	EList<char> stackRef_;
	EList<char> stackRel_;
	EList<char> stackRead_;
	bool cigDistMm_;
	bool cigCalc_;
	EList<char> cigOp_;
	EList<size_t> cigRun_;
	bool mdzCalc_;
	EList<char> mdzOp_;
	EList<char> mdzChr_;
	EList<size_t> mdzRun_;
};
class AlnRes
{
public:
	AlnRes() : ned_(RES_CAT),
			   aed_(RES_CAT)
	{
		reset();
	}
	void reset();
	void reverseEdits()
	{
		ned_.reverse();
		aed_.reverse();
	}
	void invertEdits()
	{
		assert(shapeSet_);
		assert_gt(rdlen_, 0);
		assert_gt(rdrows_, 0);
		Edit::invertPoss(ned_, rdexrows_, false);
		Edit::invertPoss(aed_, rdexrows_, false);
	}
	bool empty() const
	{
		if (!VALID_AL_SCORE(score_))
		{
			assert(ned_.empty());
			assert(aed_.empty());
			assert(!refcoord_.inited());
			assert(!refival_.inited());
			return true;
		}
		else
		{
			return false;
		}
	}
	inline TRefId refid() const
	{
		assert(shapeSet_);
		return refcoord_.ref();
	}
	inline int orient() const
	{
		assert(shapeSet_);
		return refcoord_.orient();
	}
	inline TRefOff refoff() const
	{
		assert(shapeSet_);
		return refcoord_.off();
	}
	inline void getCoords(
		Coord &st,
		Coord &en)
		const
	{
		assert(shapeSet_);
		st.init(refcoord_);
		en.init(refcoord_);
		en.adjustOff(refExtent() - 1);
	}
	inline void getExtendedCoords(
		Coord &st,
		Coord &en,
		const AlnFlags &flags)
		const
	{
		getCoords(st, en);
		if (!flags.scUnMapped())
		{
			int64_t trim_st = (fw() ? trim5p_ : trim3p_);
			int64_t trim_en = (fw() ? trim3p_ : trim5p_);
			trim_st += (fw() ? pretrim5p_ : pretrim3p_);
			trim_en += (fw() ? pretrim3p_ : pretrim5p_);
			st.adjustOff(-trim_st);
			en.adjustOff(trim_en);
		}
	}
	void setShape(
		TRefId id,
		TRefOff off,
		TRefOff reflen, bool fw, size_t rdlen, bool pretrimSoft, size_t pretrim5p, size_t pretrim3p, bool trimSoft, size_t trim5p, size_t trim3p);
	bool within(
		TRefId id,
		TRefOff off,
		bool fw,
		size_t extent) const
	{
		if (refcoord_.ref() == id &&
			refcoord_.off() >= off &&
			refcoord_.off() + refExtent() <= off + extent &&
			refcoord_.fw() == fw)
		{
			return true;
		}
		return false;
	}
	void setScore(AlnScore score)
	{
		score_ = score;
	}
	void setNucs(bool fw, int nup, int ndn)
	{
		nuc5p_ = fw ? nup : ndn;
		nuc3p_ = fw ? ndn : nup;
	}
	const Coord &refcoord() const
	{
		return refcoord_;
	}
	const Interval &refival() const
	{
		return refival_;
	}
	Coord &refcoord()
	{
		return refcoord_;
	}
	inline bool fw() const
	{
		return refcoord_.fw();
	}
	AlnScore score() const { return score_; }
	AlnScore oscore() const { return oscore_; }
	EList<Edit> &ned() { return ned_; }
	EList<Edit> &aed() { return aed_; }
	const EList<Edit> &ned() const { return ned_; }
	const EList<Edit> &aed() const { return aed_; }
	size_t readExtent() const { return rdextent_; }
	size_t readExtentRows() const { return rdexrows_; }
	size_t readLength() const { return rdlen_; }
	size_t refExtent() const
	{
		return rfextent_;
	}
	TRefOff reflen() const
	{
		return reflen_;
	}
	size_t refNucExtent() const
	{
		return rfextent_;
	}
	void printSeq(
		const Read &rd,
		const BTDnaString *dns,
		BTString &o) const;
	void printQuals(
		const Read &rd,
		const BTString *dqs,
		BTString &o) const;
	void printStacked(
		const Read &rd,
		std::ostream &o) const
	{
		printStacked(refcoord_.fw() ? rd.patFw : rd.patRc, o);
	}
	void printStacked(
		const BTDnaString &seq,
		std::ostream &o) const
	{
		Edit::printQAlign(o, seq, ned_);
		o << "^" << std::endl;
		o << "(" << refcoord_.ref() << "," << refcoord_.off() << ")" << std::endl;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert(refcoord_.repOk());
		if (shapeSet_)
		{
			assert_lt(refoff(), reflen_);
		}
		assert(refival_.repOk());
		assert(VALID_AL_SCORE(score_) || ned_.empty());
		assert(VALID_AL_SCORE(score_) || aed_.empty());
		assert(empty() || refcoord_.inited());
		assert(empty() || refival_.inited());
		assert_geq(rdexrows_, rdextent_);
		assert(empty() || rdextent_ > 0);
		assert(empty() || rfextent_ > 0);
		return true;
	}
	bool repOk(const Read &rd) const
	{
		assert(Edit::repOk(ned_, refcoord_.fw() ? rd.patFw : rd.patRc,
						   refcoord_.fw(), trimmed5p(true), trimmed3p(true)));
		return repOk();
	}
#endif
#ifndef NDEBUG
	bool matchesRef(
		const Read &rd,
		const BitPairReference &ref,
		BTDnaString &rf,
		BTDnaString &rdseq,
		BTString &qseq,
		SStringExpandable<char> &raw_refbuf,
		SStringExpandable<uint32_t> &destU32,
		EList<bool> &matches);
#endif
	void setParams(
		int seedmms,
		int seedlen,
		int seedival,
		int64_t minsc)
	{
		seedmms_ = seedmms;
		seedlen_ = seedlen;
		seedival_ = seedival;
		minsc_ = minsc;
	}
	int seedmms() const { return seedmms_; }
	int seedlen() const { return seedlen_; }
	int seedival() const { return seedival_; }
	int64_t minScore() const { return minsc_; }
	inline bool trimmedRow5p(size_t i) const
	{
		return i < trim5p_ || rdrows_ - i - 1 < trim3p_;
	}
	inline bool trimmedPos5p(size_t i) const
	{
		return i < trim5p_ || rdlen_ - i - 1 < trim3p_;
	}
	inline bool alignedRow5p(size_t i) const
	{
		return !trimmedRow5p(i);
	}
	inline bool alignedPos5p(size_t i) const
	{
		return !trimmedPos5p(i);
	}
	bool overlap(AlnRes &res);
	inline bool readUnpaired() const
	{
		assert_gt(type_, 0);
		return type_ == ALN_RES_TYPE_UNPAIRED;
	}
	inline bool alignedUnpaired() const
	{
		assert_gt(type_, 0);
		return type_ == ALN_RES_TYPE_UNPAIRED ||
			   type_ == ALN_RES_TYPE_UNPAIRED_MATE1 ||
			   type_ == ALN_RES_TYPE_UNPAIRED_MATE2;
	}
	inline bool alignedPaired() const
	{
		assert_gt(type_, 0);
		return type_ == ALN_RES_TYPE_MATE1 ||
			   type_ == ALN_RES_TYPE_MATE2;
	}
	inline bool readMate1() const
	{
		assert_gt(type_, 0);
		return type_ == ALN_RES_TYPE_MATE1 ||
			   type_ == ALN_RES_TYPE_UNPAIRED_MATE1;
	}
	inline bool alignedMate1() const
	{
		assert_gt(type_, 0);
		return type_ == ALN_RES_TYPE_MATE1;
	}
	inline bool readMate2() const
	{
		assert_gt(type_, 0);
		return type_ == ALN_RES_TYPE_MATE2 ||
			   type_ == ALN_RES_TYPE_UNPAIRED_MATE2;
	}
	inline bool alignedMate2() const
	{
		assert_gt(type_, 0);
		return type_ == ALN_RES_TYPE_MATE2;
	}
	bool isFraglenSet() const
	{
		return fraglenSet_;
	}
	void setMateParams(
		int type,
		const AlnRes *omate,
		const AlnFlags &flags)
	{
		assert_gt(type, 0);
		type_ = type;
		fraglen_ = 0;
		if (omate != NULL)
		{
			oscore_ = omate->score_;
			bool sameChr = true;
			if ((sameChr && refcoord_.ref() == omate->refcoord_.ref()) ||
				flags.alignedConcordant())
			{
				setFragmentLength(*omate, flags);
			}
			else
			{
				assert(!isFraglenSet());
			}
		}
	}
	int64_t setFragmentLength(const AlnRes &omate, const AlnFlags &flags)
	{
		Coord st, en;
		Coord ost, oen;
		assert_eq(refid(), omate.refid());
		getExtendedCoords(st, en, flags);
		omate.getExtendedCoords(ost, oen, flags);
		bool imUpstream;
		if (st.off() == ost.off())
		{
			if (st.fw() && ost.fw() && readMate1())
			{
				imUpstream = true;
			}
			else if (st.fw() && !ost.fw())
			{
				imUpstream = true;
			}
			else
			{
				imUpstream = false;
			}
		}
		else if (st.off() < ost.off())
		{
			imUpstream = true;
		}
		else
		{
			imUpstream = false;
		}
		TRefOff up = std::min(st.off(), ost.off());
		TRefOff dn = std::max(en.off(), oen.off());
		assert_geq(dn, up);
		fraglen_ = 1 + dn - up;
		if (!imUpstream)
		{
			fraglen_ = -fraglen_;
		}
		fraglenSet_ = true;
		return fraglen_;
	}
	int64_t fragmentLength() const
	{
		assert_gt(type_, 0);
		assert(fraglenSet_);
		return fraglen_;
	}
	void init(
		size_t rdlen,
		AlnScore score,
		const EList<Edit> *ned, size_t ned_i, size_t ned_n, const EList<Edit> *aed, size_t aed_i, size_t aed_n, Coord refcoord, TRefOff reflen, int seedmms = -1, int seedlen = -1, int seedival = -1, int64_t minsc = -1, int nuc5p = -1, int nuc3p = -1,
		bool pretrimSoft = false,
		size_t pretrim5p = 0,
		size_t pretrim3p = 0,
		bool trimSoft = true,
		size_t trim5p = 0,
		size_t trim3p = 0);
	size_t trimmed5p(bool soft) const
	{
		size_t trim = 0;
		if (pretrimSoft_ == soft)
			trim += pretrim5p_;
		if (trimSoft_ == soft)
			trim += trim5p_;
		return trim;
	}
	size_t trimmed3p(bool soft) const
	{
		size_t trim = 0;
		if (pretrimSoft_ == soft)
			trim += pretrim3p_;
		if (trimSoft_ == soft)
			trim += trim3p_;
		return trim;
	}
	size_t trimmedLeft(bool soft) const
	{
		return fw() ? trimmed5p(soft) : trimmed3p(soft);
	}
	size_t trimmedRight(bool soft) const
	{
		return fw() ? trimmed3p(soft) : trimmed5p(soft);
	}
	void setRefNs(size_t refns)
	{
		refns_ = refns;
	}
	size_t refNs() const { return refns_; }
	void clipOutside(bool soft, TRefOff refi, TRefOff reff);
	void clipLeft(size_t rd_amt, size_t rf_amt);
	void clipRight(size_t rd_amt, size_t rf_amt);
	ASSERT_ONLY(BTDnaString drd);
	bool operator<(const AlnRes &o) const
	{
		return score_ > o.score_;
	}
	bool operator==(const AlnRes &o) const
	{
		return shapeSet_ == o.shapeSet_ &&
			   rdlen_ == o.rdlen_ &&
			   rdrows_ == o.rdrows_ &&
			   score_ == o.score_ &&
			   ned_ == o.ned_ &&
			   aed_ == o.aed_ &&
			   refcoord_ == o.refcoord_ &&
			   reflen_ == o.reflen_ &&
			   refival_ == o.refival_ &&
			   rdextent_ == o.rdextent_ &&
			   rdexrows_ == o.rdexrows_ &&
			   rfextent_ == o.rfextent_ &&
			   seedmms_ == o.seedmms_ &&
			   seedlen_ == o.seedlen_ &&
			   seedival_ == o.seedival_ &&
			   minsc_ == o.minsc_ &&
			   nuc5p_ == o.nuc5p_ &&
			   nuc3p_ == o.nuc3p_ &&
			   refns_ == o.refns_ &&
			   type_ == o.type_ &&
			   fraglen_ == o.fraglen_ &&
			   pretrimSoft_ == o.pretrimSoft_ &&
			   pretrim5p_ == o.pretrim5p_ &&
			   pretrim3p_ == o.pretrim3p_ &&
			   trimSoft_ == o.trimSoft_ &&
			   trim5p_ == o.trim5p_ &&
			   trim3p_ == o.trim3p_;
	}
	void initStacked(const Read &rd, StackedAln &st) const
	{
		size_t trimLS = trimmed5p(true);
		size_t trimLH = trimmed5p(false);
		size_t trimRS = trimmed3p(true);
		size_t trimRH = trimmed3p(false);
		size_t len_trimmed = rd.length() - trimLS - trimRS;
		if (!fw())
		{
			Edit::invertPoss(const_cast<EList<Edit> &>(ned_), len_trimmed, false);
			swap(trimLS, trimRS);
			swap(trimLH, trimRH);
		}
		st.init(
			fw() ? rd.patFw : rd.patRc,
			ned_, trimLS, trimLH, trimRS, trimRH);
		if (!fw())
		{
			Edit::invertPoss(const_cast<EList<Edit> &>(ned_), len_trimmed, false);
		}
	}
	void calcRefExtent()
	{
		assert_gt(rdextent_, 0);
		rfextent_ = rdextent_;
		for (size_t i = 0; i < ned_.size(); i++)
		{
			if (ned_[i].isRefGap())
				rfextent_--;
			if (ned_[i].isReadGap())
				rfextent_++;
		}
	}
	bool shapeSet_;
	size_t rdlen_;
	size_t rdrows_;
	AlnScore score_;
	AlnScore oscore_;
	EList<Edit> ned_;
	EList<Edit> aed_;
	Coord refcoord_;
	TRefOff reflen_;
	Interval refival_;
	size_t rdextent_;
	size_t rdexrows_;
	size_t rfextent_;
	int seedmms_;
	int seedlen_;
	int seedival_;
	int64_t minsc_;
	int nuc5p_;
	int nuc3p_;
	size_t refns_;
	int type_;
	bool fraglenSet_;
	int64_t fraglen_;
	string mycigar;
	bool pretrimSoft_;
	size_t pretrim5p_;
	size_t pretrim3p_;
	bool trimSoft_;
	size_t trim5p_;
	size_t trim3p_;
};
struct RedundantCell
{
	RedundantCell()
	{
		rfid = 0;
		fw = true;
		rfoff = 0;
		rdoff = 0;
	}
	RedundantCell(
		TRefId rfid_,
		bool fw_,
		TRefOff rfoff_,
		size_t rdoff_)
	{
		init(rfid_, fw_, rfoff_, rdoff_);
	}
	void init(
		TRefId rfid_,
		bool fw_,
		TRefOff rfoff_,
		size_t rdoff_)
	{
		rfid = rfid_;
		fw = fw_;
		rfoff = rfoff_;
		rdoff = rdoff_;
	}
	inline bool operator<(const RedundantCell &c) const
	{
		if (rfid < c.rfid)
			return true;
		if (rfid > c.rfid)
			return false;
		if (!fw && c.fw)
			return true;
		if (fw && !c.fw)
			return false;
		if (rfoff < c.rfoff)
			return true;
		if (rfoff > c.rfoff)
			return false;
		return rdoff < c.rdoff;
	}
	inline bool operator>(const RedundantCell &c) const
	{
		if (rfid > c.rfid)
			return true;
		if (rfid < c.rfid)
			return false;
		if (fw && !c.fw)
			return true;
		if (!fw && c.fw)
			return false;
		if (rfoff > c.rfoff)
			return true;
		if (rfoff < c.rfoff)
			return false;
		return rdoff > c.rdoff;
	}
	inline bool operator==(const RedundantCell &c) const
	{
		return rfid == c.rfid &&
			   fw == c.fw &&
			   rfoff == c.rfoff &&
			   rdoff == c.rdoff;
	}
	TRefId rfid;
	bool fw;
	TRefOff rfoff;
	size_t rdoff;
};
class RedundantAlns
{
public:
	RedundantAlns(int cat = DP_CAT) : cells_(cat) {}
	void reset() { cells_.clear(); }
	void init(size_t npos)
	{
		cells_.resize(npos);
		for (size_t i = 0; i < npos; i++)
		{
			cells_[i].clear();
		}
	}
	void add(const AlnRes &res);
	bool overlap(const AlnRes &res);
protected:
	EList<ESet<RedundantCell>> cells_;
};
typedef uint64_t TNumAlns;
class AlnSetSumm
{
public:
	AlnSetSumm() { reset(); }
	explicit AlnSetSumm(
		const Read *rd1,
		const Read *rd2,
		const EList<AlnRes> *rs1,
		const EList<AlnRes> *rs2,
		const EList<AlnRes> *rs1u,
		const EList<AlnRes> *rs2u,
		bool exhausted1,
		bool exhausted2,
		TRefId orefid,
		TRefOff orefoff)
	{
		init(rd1, rd2, rs1, rs2, rs1u, rs2u, exhausted1, exhausted2,
			 orefid, orefoff);
	}
	explicit AlnSetSumm(
		TNumAlns other1,
		TNumAlns other2,
		bool paired,
		bool exhausted1,
		bool exhausted2,
		TRefId orefid,
		TRefOff orefoff)
	{
		init(
			other1,
			other2,
			paired,
			exhausted1,
			exhausted2,
			orefid,
			orefoff);
	}
	void reset()
	{
		bestUScore_.invalidate();
		bestP1Score_.invalidate();
		bestP2Score_.invalidate();
		bestCScore_.invalidate();
		bestUDist_.invalidate();
		bestP1Dist_.invalidate();
		bestP2Dist_.invalidate();
		bestCDist_.invalidate();
		bestUnchosenUScore_.invalidate();
		bestUnchosenP1Score_.invalidate();
		bestUnchosenP2Score_.invalidate();
		bestUnchosenCScore_.invalidate();
		bestUnchosenUDist_.invalidate();
		bestUnchosenP1Dist_.invalidate();
		bestUnchosenP2Dist_.invalidate();
		bestUnchosenCDist_.invalidate();
		other1_ = other2_ = 0;
		paired_ = false;
		exhausted1_ = exhausted2_ = false;
		orefid_ = -1;
		orefoff_ = -1;
	}
	void init(
		const Read *rd1,
		const Read *rd2,
		const EList<AlnRes> *rs1,
		const EList<AlnRes> *rs2,
		const EList<AlnRes> *rs1u,
		const EList<AlnRes> *rs2u,
		bool exhausted1,
		bool exhausted2,
		TRefId orefid,
		TRefOff orefoff);
	void init(
		TNumAlns other1,
		TNumAlns other2,
		bool paired,
		bool exhausted1,
		bool exhausted2,
		TRefId orefid,
		TRefOff orefoff)
	{
		other1_ = other1;
		other2_ = other2;
		paired_ = paired;
		exhausted1_ = exhausted1;
		exhausted2_ = exhausted2;
		orefid_ = orefid;
		orefoff_ = orefoff;
		assert(repOk());
	}
	bool empty() const
	{
		assert(repOk());
		return !VALID_AL_SCORE(bestScore(true));
	}
#ifndef NDEBUG
	bool repOk() const
	{
		return true;
	}
#endif
	TNumAlns other1() const
	{
		return other1_;
	}
	TNumAlns other2() const { return other2_; }
	bool paired() const { return paired_; }
	bool exhausted1() const { return exhausted1_; }
	bool exhausted2() const { return exhausted2_; }
	TRefId orefid() const { return orefid_; }
	TRefOff orefoff() const { return orefoff_; }
	AlnScore bestUScore() const { return bestUScore_; }
	AlnScore bestP1Score() const { return bestP1Score_; }
	AlnScore bestP2Score() const { return bestP2Score_; }
	AlnScore bestCScore() const { return bestCScore_; }
	AlnScore bestUDist() const { return bestUDist_; }
	AlnScore bestP1Dist() const { return bestP1Dist_; }
	AlnScore bestP2Dist() const { return bestP2Dist_; }
	AlnScore bestCDist() const { return bestCDist_; }
	AlnScore bestUnchosenUScore() const { return bestUnchosenUScore_; }
	AlnScore bestUnchosenP1Score() const { return bestUnchosenP1Score_; }
	AlnScore bestUnchosenP2Score() const { return bestUnchosenP2Score_; }
	AlnScore bestUnchosenCScore() const { return bestUnchosenCScore_; }
	AlnScore bestUnchosenUDist() const { return bestUnchosenUDist_; }
	AlnScore bestUnchosenP1Dist() const { return bestUnchosenP1Dist_; }
	AlnScore bestUnchosenP2Dist() const { return bestUnchosenP2Dist_; }
	AlnScore bestUnchosenCDist() const { return bestUnchosenCDist_; }
	AlnScore bestUnchosenPScore(bool mate1) const
	{
		return mate1 ? bestUnchosenP1Score_ : bestUnchosenP2Score_;
	}
	AlnScore bestUnchosenPDist(bool mate1) const
	{
		return mate1 ? bestUnchosenP1Dist_ : bestUnchosenP2Dist_;
	}
	AlnScore bestUnchosenScore(bool mate1) const
	{
		return paired_ ? (mate1 ? bestUnchosenP1Score_ : bestUnchosenP2Score_) : bestUnchosenUScore();
	}
	AlnScore bestUnchosenDist(bool mate1) const
	{
		return paired_ ? (mate1 ? bestUnchosenP1Dist_ : bestUnchosenP2Dist_) : bestUnchosenUDist();
	}
	bool exhausted(bool mate1) const
	{
		return mate1 ? exhausted1_ : exhausted2_;
	}
	AlnScore bestScore(bool mate1) const
	{
		return paired_ ? (mate1 ? bestP1Score_ : bestP2Score_) : bestUScore_;
	}
	AlnScore bestDist(bool mate1) const
	{
		return paired_ ? (mate1 ? bestP1Dist_ : bestP2Dist_) : bestUDist_;
	}
	void setBest(
		AlnScore bestUScore,
		AlnScore bestUDist,
		AlnScore bestP1Score,
		AlnScore bestP1Dist,
		AlnScore bestP2Score,
		AlnScore bestP2Dist,
		AlnScore bestCScore,
		AlnScore bestCDist,
		AlnScore bestUnchosenUScore,
		AlnScore bestUnchosenUDist,
		AlnScore bestUnchosenP1Score,
		AlnScore bestUnchosenP1Dist,
		AlnScore bestUnchosenP2Score,
		AlnScore bestUnchosenP2Dist,
		AlnScore bestUnchosenCScore,
		AlnScore bestUnchosenCDist)
	{
		assert(bestUScore.valid() == bestUDist.valid());
		assert(bestP1Score.valid() == bestP1Dist.valid());
		assert(bestP2Score.valid() == bestP2Dist.valid());
		assert(bestCScore.valid() == bestCDist.valid());
		assert(bestUnchosenUScore.valid() == bestUnchosenUDist.valid());
		assert(bestUnchosenP1Score.valid() == bestUnchosenP1Dist.valid());
		assert(bestUnchosenP2Score.valid() == bestUnchosenP2Dist.valid());
		assert(bestUnchosenCScore.valid() == bestUnchosenCDist.valid());
		bestUScore_ = bestUScore;
		bestUDist_ = bestUDist;
		bestP1Score_ = bestP1Score;
		bestP1Dist_ = bestP1Dist;
		bestP2Score_ = bestP2Score;
		bestP2Dist_ = bestP2Dist;
		bestCScore_ = bestCScore;
		bestCDist_ = bestCDist;
		bestUnchosenUScore_ = bestUnchosenUScore;
		bestUnchosenUDist_ = bestUnchosenUDist;
		bestUnchosenP1Score_ = bestUnchosenP1Score;
		bestUnchosenP1Dist_ = bestUnchosenP1Dist;
		bestUnchosenP2Score_ = bestUnchosenP2Score;
		bestUnchosenP2Dist_ = bestUnchosenP2Dist;
		bestUnchosenCScore_ = bestUnchosenCScore;
		bestUnchosenCDist_ = bestUnchosenCDist;
	}
protected:
	TNumAlns other1_;
	TNumAlns other2_;
	bool paired_;
	bool exhausted1_;
	bool exhausted2_;
	TRefId orefid_;
	TRefOff orefoff_;
	AlnScore bestUScore_;
	AlnScore bestUDist_;
	AlnScore bestP1Score_;
	AlnScore bestP1Dist_;
	AlnScore bestP2Score_;
	AlnScore bestP2Dist_;
	AlnScore bestCScore_;
	AlnScore bestCDist_;
	AlnScore bestUnchosenUScore_;
	AlnScore bestUnchosenUDist_;
	AlnScore bestUnchosenP1Score_;
	AlnScore bestUnchosenP1Dist_;
	AlnScore bestUnchosenP2Score_;
	AlnScore bestUnchosenP2Dist_;
	AlnScore bestUnchosenCScore_;
	AlnScore bestUnchosenCDist_;
};
#endif
