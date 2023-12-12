#ifndef ALN_SINK_H_
#define ALN_SINK_H_
#include <limits>
#include "read.h"
#include "unique.h"
#include "sam.h"
#include "ds.h"
#include "simple_func.h"
#include "outq.h"
#include <utility>
class SeedResults;
enum
{
	OUTPUT_SAM = 1
};
struct ReportingMetrics
{
	ReportingMetrics() : mutex_m()
	{
		reset();
	}
	void reset()
	{
		init(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	}
	void init(
		uint64_t nread_,
		uint64_t npaired_,
		uint64_t nunpaired_,
		uint64_t nconcord_uni_,
		uint64_t nconcord_uni1_,
		uint64_t nconcord_uni2_,
		uint64_t nconcord_rep_,
		uint64_t nconcord_0_,
		uint64_t ndiscord_,
		uint64_t nunp_0_uni_,
		uint64_t nunp_0_uni1_,
		uint64_t nunp_0_uni2_,
		uint64_t nunp_0_rep_,
		uint64_t nunp_0_0_,
		uint64_t nunp_rep_uni_,
		uint64_t nunp_rep_uni1_,
		uint64_t nunp_rep_uni2_,
		uint64_t nunp_rep_rep_,
		uint64_t nunp_rep_0_,
		uint64_t nunp_uni_,
		uint64_t nunp_uni1_,
		uint64_t nunp_uni2_,
		uint64_t nunp_rep_,
		uint64_t nunp_0_,
		uint64_t sum_best1_,
		uint64_t sum_best2_,
		uint64_t sum_best_)
	{
		nread = nread_;
		npaired = npaired_;
		nunpaired = nunpaired_;
		nconcord_uni = nconcord_uni_;
		nconcord_uni1 = nconcord_uni1_;
		nconcord_uni2 = nconcord_uni2_;
		nconcord_rep = nconcord_rep_;
		nconcord_0 = nconcord_0_;
		ndiscord = ndiscord_;
		nunp_0_uni = nunp_0_uni_;
		nunp_0_uni1 = nunp_0_uni1_;
		nunp_0_uni2 = nunp_0_uni2_;
		nunp_0_rep = nunp_0_rep_;
		nunp_0_0 = nunp_0_0_;
		nunp_rep_uni = nunp_rep_uni_;
		nunp_rep_uni1 = nunp_rep_uni1_;
		nunp_rep_uni2 = nunp_rep_uni2_;
		nunp_rep_rep = nunp_rep_rep_;
		nunp_rep_0 = nunp_rep_0_;
		nunp_uni = nunp_uni_;
		nunp_uni1 = nunp_uni1_;
		nunp_uni2 = nunp_uni2_;
		nunp_rep = nunp_rep_;
		nunp_0 = nunp_0_;
		sum_best1 = sum_best1_;
		sum_best2 = sum_best2_;
		sum_best = sum_best_;
	}
	void merge(const ReportingMetrics &met)
	{
		ThreadSafe ts(mutex_m);
		nread += met.nread;
		npaired += met.npaired;
		nunpaired += met.nunpaired;
		nconcord_uni += met.nconcord_uni;
		nconcord_uni1 += met.nconcord_uni1;
		nconcord_uni2 += met.nconcord_uni2;
		nconcord_rep += met.nconcord_rep;
		nconcord_0 += met.nconcord_0;
		ndiscord += met.ndiscord;
		nunp_0_uni += met.nunp_0_uni;
		nunp_0_uni1 += met.nunp_0_uni1;
		nunp_0_uni2 += met.nunp_0_uni2;
		nunp_0_rep += met.nunp_0_rep;
		nunp_0_0 += met.nunp_0_0;
		nunp_rep_uni += met.nunp_rep_uni;
		nunp_rep_uni1 += met.nunp_rep_uni1;
		nunp_rep_uni2 += met.nunp_rep_uni2;
		nunp_rep_rep += met.nunp_rep_rep;
		nunp_rep_0 += met.nunp_rep_0;
		nunp_uni += met.nunp_uni;
		nunp_uni1 += met.nunp_uni1;
		nunp_uni2 += met.nunp_uni2;
		nunp_rep += met.nunp_rep;
		nunp_0 += met.nunp_0;
		sum_best1 += met.sum_best1;
		sum_best2 += met.sum_best2;
		sum_best += met.sum_best;
	}
	uint64_t nread;
	uint64_t npaired;
	uint64_t nunpaired;
	uint64_t nconcord_uni;
	uint64_t nconcord_uni1;
	uint64_t nconcord_uni2;
	uint64_t nconcord_rep;
	uint64_t nconcord_0;
	uint64_t ndiscord;
	uint64_t nunp_0_uni;
	uint64_t nunp_0_uni1;
	uint64_t nunp_0_uni2;
	uint64_t nunp_0_rep;
	uint64_t nunp_0_0;
	uint64_t nunp_rep_uni;
	uint64_t nunp_rep_uni1;
	uint64_t nunp_rep_uni2;
	uint64_t nunp_rep_rep;
	uint64_t nunp_rep_0;
	uint64_t nunp_uni;
	uint64_t nunp_uni1;
	uint64_t nunp_uni2;
	uint64_t nunp_rep;
	uint64_t nunp_0;
	uint64_t sum_best1;
	uint64_t sum_best2;
	uint64_t sum_best;
	MUTEX_T mutex_m;
};
typedef int64_t THitInt;
struct ReportingParams
{
	explicit ReportingParams(
		THitInt khits_,
		THitInt mhits_,
		THitInt pengap_,
		bool msample_,
		bool discord_,
		bool mixed_)
	{
		init(khits_, mhits_, pengap_, msample_, discord_, mixed_);
	}
	void init(
		THitInt khits_,
		THitInt mhits_,
		THitInt pengap_,
		bool msample_,
		bool discord_,
		bool mixed_)
	{
		khits = khits_;
		mhits = ((mhits_ == 0) ? std::numeric_limits<THitInt>::max() : mhits_);
		pengap = pengap_;
		msample = msample_;
		discord = discord_;
		mixed = mixed_;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_geq(khits, 1);
		assert_geq(mhits, 1);
		return true;
	}
#endif
	inline bool mhitsSet() const
	{
		return mhits < std::numeric_limits<THitInt>::max();
	}
	inline THitInt mult() const
	{
		if (mhitsSet())
		{
			return mhits + 1;
		}
		return khits;
	}
	void boostThreshold(SimpleFunc &func)
	{
		THitInt mul = mult();
		assert_gt(mul, 0);
		if (mul == std::numeric_limits<THitInt>::max())
		{
			func.setMin(std::numeric_limits<double>::max());
		}
		else if (mul > 1)
		{
			func.mult(mul);
		}
	}
	bool allHits() const
	{
		return khits == std::numeric_limits<THitInt>::max();
	}
	THitInt khits;
	THitInt mhits, pengap;
	bool msample;
	bool discord;
	bool mixed;
};
class ReportingState
{
public:
	enum
	{
		NO_READ = 1,
		CONCORDANT_PAIRS,
		DISCORDANT_PAIRS,
		UNPAIRED,
		DONE
	};
	enum
	{
		EXIT_DID_NOT_EXIT = 1,
		EXIT_DID_NOT_ENTER,
		EXIT_SHORT_CIRCUIT_k,
		EXIT_SHORT_CIRCUIT_M,
		EXIT_SHORT_CIRCUIT_TRUMPED,
		EXIT_CONVERTED_TO_DISCORDANT,
		EXIT_NO_ALIGNMENTS,
		EXIT_WITH_ALIGNMENTS
	};
	ReportingState(const ReportingParams &p) : p_(p) { reset(); }
	void reset()
	{
		state_ = ReportingState::NO_READ;
		paired_ = false;
		nconcord_ = 0;
		ndiscord_ = 0;
		nunpair1_ = 0;
		nunpair2_ = 0;
		doneConcord_ = false;
		doneDiscord_ = false;
		doneUnpair_ = false;
		doneUnpair1_ = false;
		doneUnpair2_ = false;
		exitConcord_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitDiscord_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair1_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair2_ = ReportingState::EXIT_DID_NOT_ENTER;
		done_ = false;
	}
	bool inited() const { return state_ != ReportingState::NO_READ; }
	void nextRead(bool paired);
	bool foundConcordant();
	bool foundUnpaired(bool mate1);
	void finish();
	void getReport(
		uint64_t &nconcordAln,
		uint64_t &ndiscordAln,
		uint64_t &nunpair1Aln,
		uint64_t &nunpair2Aln,
		bool &pairMax,
		bool &unpair1Max,
		bool &unpair2Max)
		const;
	inline int state() const { return state_; }
	inline bool doneConcordant() const { return doneConcord_; }
	inline bool doneDiscordant() const { return doneDiscord_; }
	inline bool doneUnpaired(bool mate1) const
	{
		return mate1 ? doneUnpair1_ : doneUnpair2_;
	}
	inline bool doneWithMate(bool mate1) const
	{
		bool doneUnpair = mate1 ? doneUnpair1_ : doneUnpair2_;
		uint64_t nun = mate1 ? nunpair1_ : nunpair2_;
		if (!doneUnpair || !doneConcord_)
		{
			return false;
		}
		if (!doneDiscord_ && nun == 0)
		{
			return false;
		}
		return true;
	}
	inline bool doneUnpaired() const { return doneUnpair_; }
	inline bool done() const { return done_; }
	inline uint64_t numConcordant() const { return nconcord_; }
	inline uint64_t numDiscordant() const { return ndiscord_; }
	inline uint64_t numUnpaired1() const { return nunpair1_; }
	inline uint64_t numUnpaired2() const { return nunpair2_; }
	inline int exitConcordant() const { return exitConcord_; }
	inline int exitDiscordant() const { return exitDiscord_; }
	inline int exitUnpaired1() const { return exitUnpair1_; }
	inline int exitUnpaired2() const { return exitUnpair2_; }
#ifndef NDEBUG
	bool repOk() const
	{
		assert(p_.discord || doneDiscord_);
		assert(p_.mixed || !paired_ || doneUnpair_);
		assert(doneUnpair_ || !doneUnpair1_ || !doneUnpair2_);
		if (p_.mhitsSet())
		{
			assert_leq(numConcordant(), (uint64_t)p_.mhits + 1);
			assert_leq(numDiscordant(), (uint64_t)p_.mhits + 1);
			assert(paired_ || numUnpaired1() <= (uint64_t)p_.mhits + 1);
			assert(paired_ || numUnpaired2() <= (uint64_t)p_.mhits + 1);
		}
		assert(done() || !doneWithMate(true) || !doneWithMate(false));
		return true;
	}
#endif
	const ReportingParams &params() const
	{
		return p_;
	}
	void convertUnpairedToDiscordant()
	{
		assert_eq(1, numUnpaired1());
		assert_eq(1, numUnpaired2());
		assert_eq(0, numDiscordant());
		exitUnpair1_ = exitUnpair2_ = ReportingState::EXIT_CONVERTED_TO_DISCORDANT;
		nunpair1_ = nunpair2_ = 0;
		ndiscord_ = 1;
		assert_eq(1, numDiscordant());
	}
	inline void areDone(
		uint64_t cnt,
		bool &done,
		int &exit) const;
	inline void updateDone()
	{
		doneUnpair_ = doneUnpair1_ && doneUnpair2_;
		done_ = doneUnpair_ && doneDiscord_ && doneConcord_;
	}
	const ReportingParams &p_;
	int state_;
	bool paired_;
	uint64_t nconcord_;
	uint64_t ndiscord_;
	uint64_t nunpair1_;
	uint64_t nunpair2_;
	bool doneConcord_;
	bool doneDiscord_;
	bool doneUnpair_;
	bool doneUnpair1_;
	bool doneUnpair2_;
	int exitConcord_;
	int exitDiscord_;
	int exitUnpair1_;
	int exitUnpair2_;
	bool done_;
};
class AlnSink
{
	typedef EList<std::string> StrList;
public:
	explicit AlnSink(
		OutputQueue &oq,
		const StrList &refnames,
		bool quiet) : oq_(oq),
					  refnames_(refnames),
					  quiet_(quiet)
	{
	}
	virtual ~AlnSink() {}
	virtual void append(
		BTString &o,
		StackedAln &staln,
		size_t threadId,
		const Read *rd1,
		const Read *rd2,
		const TReadId rdid,
		AlnRes *rs1,
		AlnRes *rs2,
		const AlnSetSumm &summ,
		const SeedAlSumm &ssm1,
		const SeedAlSumm &ssm2,
		const AlnFlags *flags1,
		const AlnFlags *flags2,
		const PerReadMetrics &prm,
		const Mapq &mapq,
		const Scoring &sc,
		bool report2) = 0;
	virtual void reportHits(
		BTString &o,
		StackedAln &staln,
		size_t threadId,
		const Read *rd1,
		const Read *rd2,
		const TReadId rdid,
		const EList<size_t> &select1,
		const EList<size_t> *select2,
		EList<AlnRes> *rs1,
		EList<AlnRes> *rs2,
		bool maxed,
		const AlnSetSumm &summ,
		const SeedAlSumm &ssm1,
		const SeedAlSumm &ssm2,
		const AlnFlags *flags1,
		const AlnFlags *flags2,
		const PerReadMetrics &prm,
		const Mapq &mapq,
		const Scoring &sc,
		bool reportBoth,
		bool getLock = true)
	{
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);
		AlnFlags flagscp1, flagscp2;
		if (flags1 != NULL)
		{
			flagscp1 = *flags1;
			flags1 = &flagscp1;
			flagscp1.setPrimary(true);
		}
		if (flags2 != NULL)
		{
			flagscp2 = *flags2;
			flags2 = &flagscp2;
			flagscp2.setPrimary(true);
		}
		if (select2 != NULL)
		{
			assert(rd1 != NULL);
			assert(flags1 != NULL);
			assert(rd2 != NULL);
			assert(flags2 != NULL);
			assert_gt(select1.size(), 0);
			assert_gt(select2->size(), 0);
			AlnRes *r1pri = ((rs1 != NULL) ? &rs1->get(select1[0]) : NULL);
			AlnRes *r2pri = ((rs2 != NULL) ? &rs2->get((*select2)[0]) : NULL);
			append(o, staln, threadId, rd1, rd2, rdid, r1pri, r2pri, summ,
				   ssm1, ssm2, flags1, flags2, prm, mapq, sc, false);
			flagscp1.setPrimary(false);
			flagscp2.setPrimary(false);
			for (size_t i = 1; i < select1.size(); i++)
			{
				AlnRes *r1 = ((rs1 != NULL) ? &rs1->get(select1[i]) : NULL);
				append(o, staln, threadId, rd1, rd2, rdid, r1, r2pri, summ,
					   ssm1, ssm2, flags1, flags2, prm, mapq, sc, false);
			}
			if (reportBoth)
			{
				for (size_t i = 1; i < select2->size(); i++)
				{
					AlnRes *r2 = ((rs2 != NULL) ? &rs2->get((*select2)[i]) : NULL);
					append(o, staln, threadId, rd2, rd1, rdid, r2, r1pri, summ,
						   ssm2, ssm1, flags2, flags1, prm, mapq, sc, false);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < select1.size(); i++)
			{
				AlnRes *r1 = ((rs1 != NULL) ? &rs1->get(select1[i]) : NULL);
				AlnRes *r2 = ((rs2 != NULL) ? &rs2->get(select1[i]) : NULL);
				append(o, staln, threadId, rd1, rd2, rdid, r1, r2, summ,
					   ssm1, ssm2, flags1, flags2, prm, mapq, sc, true);
				if (flags1 != NULL)
				{
					flagscp1.setPrimary(false);
				}
				if (flags2 != NULL)
				{
					flagscp2.setPrimary(false);
				}
			}
		}
	}
	virtual void reportUnaligned(
		BTString &o,
		StackedAln &staln,
		size_t threadId,
		const Read *rd1,
		const Read *rd2,
		const TReadId rdid,
		const AlnSetSumm &summ,
		const SeedAlSumm &ssm1,
		const SeedAlSumm &ssm2,
		const AlnFlags *flags1,
		const AlnFlags *flags2,
		const PerReadMetrics &prm,
		const Mapq &mapq,
		const Scoring &sc,
		bool report2,
		bool getLock = true)
	{
		append(o, staln, threadId, rd1, rd2, rdid, NULL, NULL, summ,
			   ssm1, ssm2, flags1, flags2, prm, mapq, sc, report2);
	}
	void printAlSumm(
		const ReportingMetrics &met,
		size_t repThresh,
		bool discord,
		bool mixed,
		bool hadoopOut);
	void finish(
		size_t repThresh,
		bool discord,
		bool mixed,
		bool hadoopOut)
	{
		if (!quiet_)
		{
			printAlSumm(
				met_,
				repThresh,
				discord,
				mixed,
				hadoopOut);
		}
	}
#ifndef NDEBUG
	bool repOk() const
	{
		return true;
	}
#endif
	void reportSeedSummary(
		BTString &o,
		const Read &rd,
		TReadId rdid,
		size_t threadId,
		const SeedResults &rs,
		bool getLock = true);
	void reportEmptySeedSummary(
		BTString &o,
		const Read &rd,
		TReadId rdid,
		size_t threadId,
		bool getLock = true);
	virtual void appendSeedSummary(
		BTString &o,
		const Read &rd,
		const TReadId rdid,
		size_t seedsTried,
		size_t nonzero,
		size_t ranges,
		size_t elts,
		size_t seedsTriedFw,
		size_t nonzeroFw,
		size_t rangesFw,
		size_t eltsFw,
		size_t seedsTriedRc,
		size_t nonzeroRc,
		size_t rangesRc,
		size_t eltsRc);
	void mergeMetrics(const ReportingMetrics &met)
	{
		met_.merge(met);
	}
	OutputQueue &outq()
	{
		return oq_;
	}
	OutputQueue &oq_;
	const StrList &refnames_;
	bool quiet_;
	ReportingMetrics met_;
};
class AlnSinkWrap
{
public:
	AlnSinkWrap(
		AlnSink &g,
		const ReportingParams &rp,
		Mapq &mapq,
		size_t threadId) : g_(g),
						   rp_(rp),
						   threadid_(threadId),
						   mapq_(mapq),
						   init_(false),
						   maxed1_(false),
						   maxed2_(false),
						   maxedOverall_(false),
						   bestPair_(std::numeric_limits<TAlScore>::min()),
						   best2Pair_(std::numeric_limits<TAlScore>::min()),
						   bestUnp1_(std::numeric_limits<TAlScore>::min()),
						   best2Unp1_(std::numeric_limits<TAlScore>::min()),
						   bestUnp2_(std::numeric_limits<TAlScore>::min()),
						   best2Unp2_(std::numeric_limits<TAlScore>::min()),
						   rd1_(NULL),
						   rd2_(NULL),
						   rdid_(std::numeric_limits<TReadId>::max()),
						   rs1_(),
						   rs2_(),
						   rs1u_(),
						   rs2u_(),
						   select1_(),
						   select2_(),
						   st_(rp)
	{
		assert(rp_.repOk());
	}
	int nextRead(
		const Read *rd1,
		const Read *rd2,
		TReadId rdid,
		bool qualitiesMatter);
	void finishRead(
		const SeedResults *sr1,
		const SeedResults *sr2,
		bool exhaust1,
		bool exhaust2,
		bool nfilt1,
		bool nfilt2,
		bool scfilt1,
		bool scfilt2,
		bool lenfilt1,
		bool lenfilt2,
		bool qcfilt1,
		bool qcfilt2,
		RandomSource &rnd,
		ReportingMetrics &met,
		const PerReadMetrics &prm,
		const Scoring &sc,
		bool suppressSeedSummary = true,
		bool suppressAlignments = false,
		bool scUnMapped = false,
		bool xeq = false);
	bool report(
		int stage,
		const AlnRes *rs1,
		const AlnRes *rs2);
#ifndef NDEBUG
	bool repOk() const
	{
		assert_eq(rs2_.size(), rs1_.size());
		if (rp_.mhitsSet())
		{
			assert_gt(rp_.mhits, 0);
			assert_leq((int)rs1_.size(), rp_.mhits + 1);
			assert_leq((int)rs2_.size(), rp_.mhits + 1);
			assert(readIsPair() || (int)rs1u_.size() <= rp_.mhits + 1);
			assert(readIsPair() || (int)rs2u_.size() <= rp_.mhits + 1);
		}
		if (init_)
		{
			assert(rd1_ != NULL);
			assert_neq(std::numeric_limits<TReadId>::max(), rdid_);
		}
		assert_eq(st_.numConcordant() + st_.numDiscordant(), rs1_.size());
		assert(st_.repOk());
		return true;
	}
#endif
	bool empty() const
	{
		return rs1_.empty() && rs1u_.empty() && rs2u_.empty();
	}
	bool maxed() const
	{
		return maxedOverall_;
	}
	bool readIsPair() const
	{
		return rd1_ != NULL && rd2_ != NULL;
	}
	bool inited() const { return init_; }
	const ReportingState &state() const { return st_; }
	bool Mmode() const
	{
		return rp_.mhitsSet();
	}
	bool allHits() const
	{
		return rp_.allHits();
	}
	bool hasSecondBestUnp1() const
	{
		return best2Unp1_ != std::numeric_limits<TAlScore>::min();
	}
	bool hasSecondBestUnp2() const
	{
		return best2Unp2_ != std::numeric_limits<TAlScore>::min();
	}
	bool hasSecondBestPair() const
	{
		return best2Pair_ != std::numeric_limits<TAlScore>::min();
	}
	TAlScore bestUnp1() const
	{
		return bestUnp1_;
	}
	TAlScore secondBestUnp1() const
	{
		return best2Unp1_;
	}
	TAlScore bestUnp2() const
	{
		return bestUnp2_;
	}
	TAlScore secondBestUnp2() const
	{
		return best2Unp2_;
	}
	TAlScore bestPair() const
	{
		return bestPair_;
	}
	TAlScore secondBestPair() const
	{
		return best2Pair_;
	}
	bool sameRead(
		const Read *rd1,
		const Read *rd2,
		bool qualitiesMatter);
	bool prepareDiscordants();
	size_t selectAlnsToReport(
		const EList<AlnRes> &rs,
		uint64_t num,
		EList<size_t> &select,
		RandomSource &rnd)
		const;
	size_t selectByScore(
		const EList<AlnRes> *rs1,
		const EList<AlnRes> *rs2,
		uint64_t num,
		EList<size_t> &select,
		const EList<AlnRes> *rs1u,
		const EList<AlnRes> *rs2u,
		AlnScore &bestUScore,
		AlnScore &bestUDist,
		AlnScore &bestP1Score,
		AlnScore &bestP1Dist,
		AlnScore &bestP2Score,
		AlnScore &bestP2Dist,
		AlnScore &bestCScore,
		AlnScore &bestCDist,
		AlnScore &bestUnchosenUScore,
		AlnScore &bestUnchosenUDist,
		AlnScore &bestUnchosenP1Score,
		AlnScore &bestUnchosenP1Dist,
		AlnScore &bestUnchosenP2Score,
		AlnScore &bestUnchosenP2Dist,
		AlnScore &bestUnchosenCScore,
		AlnScore &bestUnchosenCDist,
		RandomSource &rnd)
		const;
	AlnSink &g_;
	ReportingParams rp_;
	size_t threadid_;
	Mapq &mapq_;
	bool init_;
	bool maxed1_;
	bool maxed2_;
	bool maxedOverall_;
	TAlScore bestPair_;
	TAlScore best2Pair_;
	TAlScore bestUnp1_;
	TAlScore best2Unp1_;
	TAlScore bestUnp2_;
	TAlScore best2Unp2_;
	const Read *rd1_;
	const Read *rd2_;
	TReadId rdid_;
	EList<AlnRes> rs1_;
	EList<AlnRes> rs2_;
	EList<AlnRes> rs1u_;
	EList<AlnRes> rs2u_;
	EList<size_t> select1_;
	EList<size_t> select2_;
	ReportingState st_;
	EList<std::pair<AlnScore, size_t>> selectBuf_;
	BTString obuf_;
	StackedAln staln_;
};
class AlnSinkSam : public AlnSink
{
	typedef EList<std::string> StrList;
public:
	AlnSinkSam(
		OutputQueue &oq,
		const SamConfig &samc,
		const StrList &refnames,
		bool quiet) : AlnSink(oq,
							  refnames,
							  quiet),
					  samc_(samc)
	{
	}
	virtual ~AlnSinkSam() {}
	virtual void append(
		BTString &o,
		StackedAln &staln,
		size_t threadId,
		const Read *rd1,
		const Read *rd2,
		const TReadId rdid,
		AlnRes *rs1,
		AlnRes *rs2,
		const AlnSetSumm &summ,
		const SeedAlSumm &ssm1,
		const SeedAlSumm &ssm2,
		const AlnFlags *flags1,
		const AlnFlags *flags2,
		const PerReadMetrics &prm,
		const Mapq &mapq,
		const Scoring &sc,
		bool report2)
	{
		assert(rd1 != NULL || rd2 != NULL);
		if (rd1 != NULL)
		{
			assert(flags1 != NULL);
			appendMate(o, staln, *rd1, rd2, rdid, rs1, rs2, summ, ssm1, ssm2, *flags1, prm, mapq, sc);
		}
		if (rd2 != NULL && report2)
		{
			assert(flags2 != NULL);
			appendMate(o, staln, *rd2, rd1, rdid, rs2, rs1, summ, ssm2, ssm1, *flags2, prm, mapq, sc);
		}
	}
protected:
	void appendMate(
		BTString &o,
		StackedAln &staln,
		const Read &rd,
		const Read *rdo,
		const TReadId rdid,
		AlnRes *rs,
		AlnRes *rso,
		const AlnSetSumm &summ,
		const SeedAlSumm &ssm,
		const SeedAlSumm &ssmo,
		const AlnFlags &flags,
		const PerReadMetrics &prm,
		const Mapq &mapq,
		const Scoring &sc);
	const SamConfig &samc_;
	BTDnaString dseq_;
	BTString dqual_;
};
#endif
