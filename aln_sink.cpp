#include <iomanip>
#include <limits>
#include "aln_sink.h"
#include "aligner_seed.h"
#include "util.h"
using namespace std;
void ReportingState::nextRead(bool paired)
{
	paired_ = paired;
	if (paired)
	{
		state_ = CONCORDANT_PAIRS;
		doneConcord_ = false;
		doneDiscord_ = p_.discord ? false : true;
		doneUnpair1_ = p_.mixed ? false : true;
		doneUnpair2_ = p_.mixed ? false : true;
		exitConcord_ = ReportingState::EXIT_DID_NOT_EXIT;
		exitDiscord_ = p_.discord ? ReportingState::EXIT_DID_NOT_EXIT : ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair1_ = p_.mixed ? ReportingState::EXIT_DID_NOT_EXIT : ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair2_ = p_.mixed ? ReportingState::EXIT_DID_NOT_EXIT : ReportingState::EXIT_DID_NOT_ENTER;
	}
	else
	{
		state_ = UNPAIRED;
		doneConcord_ = true;
		doneDiscord_ = true;
		doneUnpair1_ = false;
		doneUnpair2_ = true;
		exitConcord_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitDiscord_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair1_ = ReportingState::EXIT_DID_NOT_EXIT;
		exitUnpair2_ = ReportingState::EXIT_DID_NOT_ENTER;
	}
	doneUnpair_ = doneUnpair1_ && doneUnpair2_;
	done_ = false;
	nconcord_ = ndiscord_ = nunpair1_ = nunpair2_ = 0;
}
bool ReportingState::foundConcordant()
{
	assert(paired_);
	assert_geq(state_, ReportingState::CONCORDANT_PAIRS);
	assert(!doneConcord_);
	nconcord_++;
	areDone(nconcord_, doneConcord_, exitConcord_);
	doneDiscord_ = true;
	exitDiscord_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
	if (doneConcord_)
	{
		assert_neq(ReportingState::EXIT_NO_ALIGNMENTS, exitConcord_);
		if (exitConcord_ != ReportingState::EXIT_SHORT_CIRCUIT_M)
		{
			if (!doneUnpair1_)
			{
				doneUnpair1_ = true;
				exitUnpair1_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
			}
			if (!doneUnpair2_)
			{
				doneUnpair2_ = true;
				exitUnpair2_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
			}
		}
	}
	updateDone();
	return done();
}
bool ReportingState::foundUnpaired(bool mate1)
{
	assert_gt(state_, ReportingState::NO_READ);
	if (mate1)
	{
		nunpair1_++;
		if (!doneUnpair1_)
		{
			areDone(nunpair1_, doneUnpair1_, exitUnpair1_);
			if (doneUnpair1_)
			{
				doneUnpair_ = doneUnpair1_ && doneUnpair2_;
				updateDone();
			}
		}
		if (nunpair1_ > 1)
		{
			doneDiscord_ = true;
			exitDiscord_ = ReportingState::EXIT_NO_ALIGNMENTS;
		}
	}
	else
	{
		nunpair2_++;
		if (!doneUnpair2_)
		{
			areDone(nunpair2_, doneUnpair2_, exitUnpair2_);
			if (doneUnpair2_)
			{
				doneUnpair_ = doneUnpair1_ && doneUnpair2_;
				updateDone();
			}
		}
		if (nunpair2_ > 1)
		{
			doneDiscord_ = true;
			exitDiscord_ = ReportingState::EXIT_NO_ALIGNMENTS;
		}
	}
	return done();
}
void ReportingState::finish()
{
	if (!doneConcord_)
	{
		doneConcord_ = true;
		exitConcord_ =
			((nconcord_ > 0) ? ReportingState::EXIT_WITH_ALIGNMENTS : ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitConcord_, EXIT_DID_NOT_EXIT);
	if (!doneUnpair1_)
	{
		doneUnpair1_ = true;
		exitUnpair1_ =
			((nunpair1_ > 0) ? ReportingState::EXIT_WITH_ALIGNMENTS : ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitUnpair1_, EXIT_DID_NOT_EXIT);
	if (!doneUnpair2_)
	{
		doneUnpair2_ = true;
		exitUnpair2_ =
			((nunpair2_ > 0) ? ReportingState::EXIT_WITH_ALIGNMENTS : ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitUnpair2_, EXIT_DID_NOT_EXIT);
	if (!doneDiscord_)
	{
		assert_eq(0, ndiscord_);
		if (nconcord_ == 0 && nunpair1_ == 1 && nunpair2_ == 1)
		{
			convertUnpairedToDiscordant();
		}
		doneDiscord_ = true;
		exitDiscord_ =
			((ndiscord_ > 0) ? ReportingState::EXIT_WITH_ALIGNMENTS : ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert(!paired_ || exitDiscord_ > ReportingState::EXIT_DID_NOT_EXIT);
	doneUnpair_ = done_ = true;
	assert(done());
}
void ReportingState::getReport(
	uint64_t &nconcordAln,
	uint64_t &ndiscordAln,
	uint64_t &nunpair1Aln,
	uint64_t &nunpair2Aln,
	bool &pairMax,
	bool &unpair1Max,
	bool &unpair2Max)
	const
{
	nconcordAln = ndiscordAln = nunpair1Aln = nunpair2Aln = 0;
	pairMax = unpair1Max = unpair2Max = false;
	assert_gt(p_.khits, 0);
	assert_gt(p_.mhits, 0);
	if (paired_)
	{
		if (exitConcord_ == ReportingState::EXIT_SHORT_CIRCUIT_k)
		{
			assert_geq(nconcord_, (uint64_t)p_.khits);
			nconcordAln = p_.khits;
			return;
		}
		else if (exitConcord_ == ReportingState::EXIT_SHORT_CIRCUIT_M)
		{
			assert(p_.msample);
			assert_gt(nconcord_, 0);
			pairMax = true;
			if (p_.mixed)
			{
				unpair1Max = nunpair1_ > (uint64_t)p_.mhits;
				unpair2Max = nunpair2_ > (uint64_t)p_.mhits;
			}
			nconcordAln = 1;
			return;
		}
		else if (exitConcord_ == ReportingState::EXIT_WITH_ALIGNMENTS)
		{
			assert_gt(nconcord_, 0);
			nconcordAln = min<uint64_t>(nconcord_, p_.khits);
			return;
		}
		assert(!p_.mhitsSet() || nconcord_ <= (uint64_t)p_.mhits + 1);
		if (exitDiscord_ == ReportingState::EXIT_WITH_ALIGNMENTS)
		{
			assert(p_.discord);
			ndiscordAln = 1;
			return;
		}
	}
	assert_neq(ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED, exitUnpair1_);
	assert_neq(ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED, exitUnpair2_);
	if ((paired_ && !p_.mixed) || nunpair1_ + nunpair2_ == 0)
	{
		return;
	}
	if (exitUnpair1_ == ReportingState::EXIT_SHORT_CIRCUIT_k)
	{
		assert_geq(nunpair1_, (uint64_t)p_.khits);
		nunpair1Aln = p_.khits;
	}
	else if (exitUnpair1_ == ReportingState::EXIT_SHORT_CIRCUIT_M)
	{
		assert(p_.msample);
		assert_gt(nunpair1_, 0);
		unpair1Max = true;
		nunpair1Aln = 1;
	}
	else if (exitUnpair1_ == ReportingState::EXIT_WITH_ALIGNMENTS)
	{
		assert_gt(nunpair1_, 0);
		nunpair1Aln = min<uint64_t>(nunpair1_, (uint64_t)p_.khits);
	}
	assert(!p_.mhitsSet() || paired_ || nunpair1_ <= (uint64_t)p_.mhits + 1);
	if (exitUnpair2_ == ReportingState::EXIT_SHORT_CIRCUIT_k)
	{
		nunpair2Aln = p_.khits;
	}
	else if (exitUnpair2_ == ReportingState::EXIT_SHORT_CIRCUIT_M)
	{
		assert(p_.msample);
		assert_gt(nunpair2_, 0);
		unpair2Max = true;
		nunpair2Aln = 1;
	}
	else if (exitUnpair2_ == ReportingState::EXIT_WITH_ALIGNMENTS)
	{
		assert_gt(nunpair2_, 0);
		nunpair2Aln = min<uint64_t>(nunpair2_, (uint64_t)p_.khits);
	}
	assert(!p_.mhitsSet() || paired_ || nunpair2_ <= (uint64_t)p_.mhits + 1);
}
inline void ReportingState::areDone(
	uint64_t cnt,
	bool &done,
	int &exit) const
{
	assert(!done);
	assert_gt(p_.khits, 0);
	assert_gt(p_.mhits, 0);
	if (cnt >= (uint64_t)p_.khits && !p_.mhitsSet())
	{
		done = true;
		exit = ReportingState::EXIT_SHORT_CIRCUIT_k;
	}
	else if (p_.mhitsSet() && cnt > (uint64_t)p_.mhits)
	{
		done = true;
		assert(p_.msample);
		exit = ReportingState::EXIT_SHORT_CIRCUIT_M;
	}
}
static std::ostream &printPct(
	std::ostream &os,
	uint64_t num,
	uint64_t denom)
{
	double pct = 0.0f;
	if (denom != 0)
	{
		pct = 100.0 * (double)num / (double)denom;
	}
	os << fixed << setprecision(2) << pct << '%';
	return os;
}
void AlnSink::printAlSumm(
	const ReportingMetrics &met,
	size_t repThresh,
	bool discord,
	bool mixed,
	bool hadoopOut)
{
	bool canRep = repThresh != MAX_SIZE_T;
	if (hadoopOut)
	{
		cerr << "reporter:counter,Reads processed," << met.nread << endl;
	}
	uint64_t totread = met.nread;
	if (totread > 0)
	{
		cerr << "" << met.nread << " reads; of these:" << endl;
	}
	else
	{
		assert_eq(0, met.npaired);
		assert_eq(0, met.nunpaired);
		cerr << "" << totread << " reads" << endl;
	}
	uint64_t totpair = met.npaired;
	if (totpair > 0)
	{
		cerr << "  " << totpair << " (";
		printPct(cerr, totpair, totread);
		cerr << ") were paired; of these:" << endl;
		cerr << "    " << met.nconcord_0 << " (";
		printPct(cerr, met.nconcord_0, met.npaired);
		cerr << ") aligned concordantly 0 times" << endl;
		if (canRep)
		{
			assert_eq(met.nconcord_uni, met.nconcord_uni1 + met.nconcord_uni2);
			cerr << "    " << met.nconcord_uni1 << " (";
			printPct(cerr, met.nconcord_uni1, met.npaired);
			cerr << ") aligned concordantly exactly 1 time" << endl;
			cerr << "    " << met.nconcord_uni2 + met.nconcord_rep << " (";
			printPct(cerr, met.nconcord_uni2 + met.nconcord_rep, met.npaired);
			cerr << ") aligned concordantly >1 times" << endl;
		}
		else
		{
			assert_eq(met.nconcord_uni, met.nconcord_uni1 + met.nconcord_uni2);
			cerr << "    " << met.nconcord_uni1 << " (";
			printPct(cerr, met.nconcord_uni1, met.npaired);
			cerr << ") aligned concordantly exactly 1 time" << endl;
			cerr << "    " << met.nconcord_uni2 << " (";
			printPct(cerr, met.nconcord_uni2, met.npaired);
			cerr << ") aligned concordantly >1 times" << endl;
		}
		if (discord)
		{
			cerr << "    ----" << endl;
			cerr << "    " << met.nconcord_0
				 << " pairs aligned concordantly 0 times; of these:" << endl;
			cerr << "      " << met.ndiscord << " (";
			printPct(cerr, met.ndiscord, met.nconcord_0);
			cerr << ") aligned discordantly 1 time" << endl;
		}
		uint64_t ncondiscord_0 = met.nconcord_0 - met.ndiscord;
		if (mixed)
		{
			cerr << "    ----" << endl;
			cerr << "    " << ncondiscord_0
				 << " pairs aligned 0 times concordantly or discordantly; of these:" << endl;
			cerr << "      " << (ncondiscord_0 * 2) << " mates make up the pairs; of these:" << endl;
			cerr << "        " << met.nunp_0_0 << " "
				 << "(";
			printPct(cerr, met.nunp_0_0, ncondiscord_0 * 2);
			cerr << ") aligned 0 times" << endl;
			if (canRep)
			{
				assert_eq(met.nunp_0_uni, met.nunp_0_uni1 + met.nunp_0_uni2);
				cerr << "        " << met.nunp_0_uni1 << " (";
				printPct(cerr, met.nunp_0_uni1, ncondiscord_0 * 2);
				cerr << ") aligned exactly 1 time" << endl;
				cerr << "        " << met.nunp_0_uni2 + met.nunp_0_rep << " (";
				printPct(cerr, met.nunp_0_uni2 + met.nunp_0_rep, ncondiscord_0 * 2);
				cerr << ") aligned >1 times" << endl;
			}
			else
			{
				assert_eq(met.nunp_0_uni, met.nunp_0_uni1 + met.nunp_0_uni2);
				cerr << "        " << met.nunp_0_uni1 << " (";
				printPct(cerr, met.nunp_0_uni1, ncondiscord_0 * 2);
				cerr << ") aligned exactly 1 time" << endl;
				cerr << "        " << met.nunp_0_uni2 << " (";
				printPct(cerr, met.nunp_0_uni2, ncondiscord_0 * 2);
				cerr << ") aligned >1 times" << endl;
			}
		}
	}
	uint64_t totunpair = met.nunpaired;
	if (totunpair > 0)
	{
		cerr << "  " << totunpair << " (";
		printPct(cerr, totunpair, totread);
		cerr << ") were unpaired; of these:" << endl;
		cerr << "    " << met.nunp_0 << " (";
		printPct(cerr, met.nunp_0, met.nunpaired);
		cerr << ") aligned 0 times" << endl;
		if (hadoopOut)
		{
			cerr << "reporter:counter,Unpaired reads with 0 alignments,"
				 << met.nunpaired << endl;
		}
		if (canRep)
		{
			assert_eq(met.nunp_uni, met.nunp_uni1 + met.nunp_uni2);
			cerr << "    " << met.nunp_uni1 << " (";
			printPct(cerr, met.nunp_uni1, met.nunpaired);
			cerr << ") aligned exactly 1 time" << endl;
			cerr << "    " << met.nunp_uni2 + met.nunp_rep << " (";
			printPct(cerr, met.nunp_uni2 + met.nunp_rep, met.nunpaired);
			cerr << ") aligned >1 times" << endl;
		}
		else
		{
			assert_eq(met.nunp_uni, met.nunp_uni1 + met.nunp_uni2);
			cerr << "    " << met.nunp_uni1 << " (";
			printPct(cerr, met.nunp_uni1, met.nunpaired);
			cerr << ") aligned exactly 1 time" << endl;
			cerr << "    " << met.nunp_uni2 << " (";
			printPct(cerr, met.nunp_uni2, met.nunpaired);
			cerr << ") aligned >1 times" << endl;
		}
	}
	uint64_t tot_al_cand = totunpair + totpair * 2;
	uint64_t tot_al =
		(met.nconcord_uni + met.nconcord_rep) * 2 +
		(met.ndiscord) * 2 +
		met.nunp_0_uni +
		met.nunp_0_rep +
		met.nunp_uni +
		met.nunp_rep;
	assert_leq(tot_al, tot_al_cand);
	printPct(cerr, tot_al, tot_al_cand);
	cerr << " overall alignment rate" << endl;
}
bool AlnSinkWrap::sameRead(
	const Read *rd1,
	const Read *rd2,
	bool qualitiesMatter)
{
	bool same = false;
	if (rd1_ != NULL || rd2_ != NULL)
	{
		if ((rd1_ == NULL) == (rd1 == NULL) &&
			(rd2_ == NULL) == (rd2 == NULL))
		{
			bool m1same = (rd1 == NULL && rd1_ == NULL);
			if (!m1same)
			{
				assert(rd1 != NULL);
				assert(rd1_ != NULL);
				m1same = Read::same(
					rd1->patFw,
					rd1->qual,
					rd1_->patFw,
					rd1_->qual,
					qualitiesMatter);
			}
			if (m1same)
			{
				bool m2same = (rd2 == NULL && rd2_ == NULL);
				if (!m2same)
				{
					m2same = Read::same(
						rd2->patFw,
						rd2->qual,
						rd2_->patFw,
						rd2_->qual,
						qualitiesMatter);
				}
				same = m2same;
			}
		}
	}
	return same;
}
int AlnSinkWrap::nextRead(
	const Read *rd1,
	const Read *rd2,
	TReadId rdid,
	bool qualitiesMatter)
{
	assert(!init_);
	assert(rd1 != NULL || rd2 != NULL);
	init_ = true;
	if (rd1 != NULL)
	{
		rd1_ = rd1;
	}
	else
		rd1_ = NULL;
	if (rd2 != NULL)
	{
		rd2_ = rd2;
	}
	else
		rd2_ = NULL;
	rdid_ = rdid;
	maxed1_ = false;
	maxed2_ = false;
	maxedOverall_ = false;
	bestPair_ = best2Pair_ =
		bestUnp1_ = best2Unp1_ =
			bestUnp2_ = best2Unp2_ = std::numeric_limits<THitInt>::min();
	rs1_.clear();
	rs2_.clear();
	rs1u_.clear();
	rs2u_.clear();
	st_.nextRead(readIsPair());
	assert(empty());
	assert(!maxed());
	return 0;
}
void AlnSinkWrap::finishRead(
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
	bool suppressSeedSummary,
	bool suppressAlignments,
	bool scUnMapped,
	bool xeq)
{
	obuf_.clear();
	OutputQueueMark qqm(g_.outq(), obuf_, rdid_, threadid_);
	assert(init_);
	if (!suppressSeedSummary)
	{
		if (sr1 != NULL)
		{
			assert(rd1_ != NULL);
			g_.reportSeedSummary(obuf_, *rd1_, rdid_, threadid_, *sr1, true);
		}
		else if (rd1_ != NULL)
		{
			g_.reportEmptySeedSummary(obuf_, *rd1_, rdid_, true);
		}
		if (sr2 != NULL)
		{
			assert(rd2_ != NULL);
			g_.reportSeedSummary(obuf_, *rd2_, rdid_, threadid_, *sr2, true);
		}
		else if (rd2_ != NULL)
		{
			g_.reportEmptySeedSummary(obuf_, *rd2_, rdid_, true);
		}
	}
	if (!suppressAlignments)
	{
		st_.finish();
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		st_.getReport(
			nconcord,
			ndiscord,
			nunpair1,
			nunpair2,
			pairMax,
			unpair1Max,
			unpair2Max);
		assert_leq(nconcord, rs1_.size());
		assert_leq(nunpair1, rs1u_.size());
		assert_leq(nunpair2, rs2u_.size());
		assert_leq(ndiscord, 1);
		assert_gt(rp_.khits, 0);
		assert_gt(rp_.mhits, 0);
		assert(!pairMax || rs1_.size() >= (uint64_t)rp_.mhits);
		assert(!unpair1Max || rs1u_.size() >= (uint64_t)rp_.mhits);
		assert(!unpair2Max || rs2u_.size() >= (uint64_t)rp_.mhits);
		met.nread++;
		if (readIsPair())
		{
			met.npaired++;
		}
		else
		{
			met.nunpaired++;
		}
		if (nconcord > 0)
		{
			AlnSetSumm concordSumm(
				rd1_, rd2_, &rs1_, &rs2_, &rs1u_, &rs2u_,
				exhaust1, exhaust2, -1, -1);
			AlnScore bestUScore, bestP1Score, bestP2Score, bestCScore;
			AlnScore bestUDist, bestP1Dist, bestP2Dist, bestCDist;
			AlnScore bestUnchosenUScore, bestUnchosenP1Score, bestUnchosenP2Score, bestUnchosenCScore;
			AlnScore bestUnchosenUDist, bestUnchosenP1Dist, bestUnchosenP2Dist, bestUnchosenCDist;
			size_t off = selectByScore(
				&rs1_, &rs2_,
				nconcord, select1_,
				&rs1u_, &rs2u_,
				bestUScore,
				bestUDist,
				bestP1Score,
				bestP1Dist,
				bestP2Score,
				bestP2Dist,
				bestCScore,
				bestCDist,
				bestUnchosenUScore,
				bestUnchosenUDist,
				bestUnchosenP1Score,
				bestUnchosenP1Dist,
				bestUnchosenP2Score,
				bestUnchosenP2Dist,
				bestUnchosenCScore,
				bestUnchosenCDist,
				rnd);
			concordSumm.setBest(
				bestUScore,
				bestUDist,
				bestP1Score,
				bestP1Dist,
				bestP2Score,
				bestP2Dist,
				bestCScore,
				bestCDist,
				bestUnchosenUScore,
				bestUnchosenUDist,
				bestUnchosenP1Score,
				bestUnchosenP1Dist,
				bestUnchosenP2Score,
				bestUnchosenP2Dist,
				bestUnchosenCScore,
				bestUnchosenCDist);
			assert(concordSumm.bestScore(true).valid());
			assert(concordSumm.bestScore(false).valid());
			assert_lt(off, rs1_.size());
			const AlnRes *rs1 = &rs1_[off];
			const AlnRes *rs2 = &rs2_[off];
			AlnFlags flags1(
				ALN_FLAG_PAIR_CONCORD_MATE1,
				st_.params().mhitsSet(),
				unpair1Max,
				pairMax,
				nfilt1,
				scfilt1,
				lenfilt1,
				qcfilt1,
				st_.params().mixed,
				true,
				true,
				rs2->fw(),
				scUnMapped,
				xeq);
			AlnFlags flags2(
				ALN_FLAG_PAIR_CONCORD_MATE2,
				st_.params().mhitsSet(),
				unpair2Max,
				pairMax,
				nfilt2,
				scfilt2,
				lenfilt2,
				qcfilt2,
				st_.params().mixed,
				false,
				true,
				rs1->fw(),
				scUnMapped,
				xeq);
			SeedAlSumm ssm1, ssm2;
			sr1->toSeedAlSumm(ssm1);
			sr2->toSeedAlSumm(ssm2);
			for (size_t i = 0; i < rs1_.size(); i++)
			{
				rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1);
				rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2);
				assert_eq(abs(rs1_[i].fragmentLength()), abs(rs2_[i].fragmentLength()));
			}
			assert(!select1_.empty());
			g_.reportHits(
				obuf_,
				staln_,
				threadid_,
				rd1_,
				rd2_,
				rdid_,
				select1_,
				NULL,
				&rs1_,
				&rs2_,
				pairMax,
				concordSumm,
				ssm1,
				ssm2,
				&flags1,
				&flags2,
				prm,
				mapq_,
				sc,
				false);
			if (pairMax)
			{
				met.nconcord_rep++;
			}
			else
			{
				met.nconcord_uni++;
				assert(!rs1_.empty());
				if (rs1_.size() == 1)
				{
					met.nconcord_uni1++;
				}
				else
				{
					met.nconcord_uni2++;
				}
			}
			init_ = false;
			return;
		}
		else if (ndiscord > 0)
		{
			ASSERT_ONLY(bool ret =)
			prepareDiscordants();
			assert(ret);
			assert_eq(1, rs1_.size());
			assert_eq(1, rs2_.size());
			AlnSetSumm discordSumm(
				rd1_, rd2_, &rs1_, &rs2_, &rs1u_, &rs2u_,
				exhaust1, exhaust2, -1, -1);
			const AlnRes *rs1 = &rs1_[0];
			const AlnRes *rs2 = &rs2_[0];
			AlnFlags flags1(
				ALN_FLAG_PAIR_DISCORD_MATE1,
				st_.params().mhitsSet(),
				false,
				pairMax,
				nfilt1,
				scfilt1,
				lenfilt1,
				qcfilt1,
				st_.params().mixed,
				true,
				true,
				rs2->fw(),
				scUnMapped,
				xeq);
			AlnFlags flags2(
				ALN_FLAG_PAIR_DISCORD_MATE2,
				st_.params().mhitsSet(),
				false,
				pairMax,
				nfilt2,
				scfilt2,
				lenfilt2,
				qcfilt2,
				st_.params().mixed,
				false,
				true,
				rs1->fw(),
				scUnMapped,
				xeq);
			SeedAlSumm ssm1, ssm2;
			sr1->toSeedAlSumm(ssm1);
			sr2->toSeedAlSumm(ssm2);
			for (size_t i = 0; i < rs1_.size(); i++)
			{
				rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1);
				rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2);
				assert(rs1_[i].isFraglenSet() == rs2_[i].isFraglenSet());
				assert(!rs1_[i].isFraglenSet() || abs(rs1_[i].fragmentLength()) == abs(rs2_[i].fragmentLength()));
			}
			AlnScore bestUScore, bestP1Score, bestP2Score, bestCScore;
			AlnScore bestUDist, bestP1Dist, bestP2Dist, bestCDist;
			AlnScore bestUnchosenUScore, bestUnchosenP1Score, bestUnchosenP2Score, bestUnchosenCScore;
			AlnScore bestUnchosenUDist, bestUnchosenP1Dist, bestUnchosenP2Dist, bestUnchosenCDist;
			ASSERT_ONLY(size_t off =)
			selectByScore(
				&rs1_, &rs2_,
				ndiscord, select1_,
				&rs1u_, &rs2u_,
				bestUScore,
				bestUDist,
				bestP1Score,
				bestP1Dist,
				bestP2Score,
				bestP2Dist,
				bestCScore,
				bestCDist,
				bestUnchosenUScore,
				bestUnchosenUDist,
				bestUnchosenP1Score,
				bestUnchosenP1Dist,
				bestUnchosenP2Score,
				bestUnchosenP2Dist,
				bestUnchosenCScore,
				bestUnchosenCDist,
				rnd);
			discordSumm.setBest(
				bestUScore,
				bestUDist,
				bestP1Score,
				bestP1Dist,
				bestP2Score,
				bestP2Dist,
				bestCScore,
				bestCDist,
				bestUnchosenUScore,
				bestUnchosenUDist,
				bestUnchosenP1Score,
				bestUnchosenP1Dist,
				bestUnchosenP2Score,
				bestUnchosenP2Dist,
				bestUnchosenCScore,
				bestUnchosenCDist);
			assert_eq(0, off);
			assert(!select1_.empty());
			g_.reportHits(
				obuf_,
				staln_,
				threadid_,
				rd1_,
				rd2_,
				rdid_,
				select1_,
				NULL,
				&rs1_,
				&rs2_,
				pairMax,
				discordSumm,
				ssm1,
				ssm2,
				&flags1,
				&flags2,
				prm,
				mapq_,
				sc,
				false);
			met.nconcord_0++;
			met.ndiscord++;
			init_ = false;
			return;
		}
		assert(!pairMax);
		if (readIsPair())
		{
			met.nconcord_0++;
		}
		if (rd1_ != NULL)
		{
			if (nunpair1 > 0)
			{
				if (readIsPair())
				{
					if (unpair1Max)
						met.nunp_0_rep++;
					else
					{
						met.nunp_0_uni++;
						assert(!rs1u_.empty());
						if (rs1u_.size() == 1)
						{
							met.nunp_0_uni1++;
						}
						else
						{
							met.nunp_0_uni2++;
						}
					}
				}
				else
				{
					if (unpair1Max)
						met.nunp_rep++;
					else
					{
						met.nunp_uni++;
						assert(!rs1u_.empty());
						if (rs1u_.size() == 1)
						{
							met.nunp_uni1++;
						}
						else
						{
							met.nunp_uni2++;
						}
					}
				}
			}
			else if (unpair1Max)
			{
				if (readIsPair())
					met.nunp_0_rep++;
				else
					met.nunp_rep++;
			}
			else
			{
				if (readIsPair())
					met.nunp_0_0++;
				else
					met.nunp_0++;
			}
		}
		if (rd2_ != NULL)
		{
			if (nunpair2 > 0)
			{
				if (readIsPair())
				{
					if (unpair2Max)
						met.nunp_0_rep++;
					else
					{
						assert(!rs2u_.empty());
						met.nunp_0_uni++;
						if (rs2u_.size() == 1)
						{
							met.nunp_0_uni1++;
						}
						else
						{
							met.nunp_0_uni2++;
						}
					}
				}
				else
				{
					if (unpair2Max)
						met.nunp_rep++;
					else
					{
						assert(!rs2u_.empty());
						met.nunp_uni++;
						if (rs2u_.size() == 1)
						{
							met.nunp_uni1++;
						}
						else
						{
							met.nunp_uni2++;
						}
					}
				}
			}
			else if (unpair2Max)
			{
				if (readIsPair())
					met.nunp_0_rep++;
				else
					met.nunp_rep++;
			}
			else
			{
				if (readIsPair())
					met.nunp_0_0++;
				else
					met.nunp_0++;
			}
		}
		const AlnRes *repRs1 = NULL, *repRs2 = NULL;
		AlnSetSumm summ1, summ2;
		AlnFlags flags1, flags2;
		TRefId refid = -1;
		TRefOff refoff = -1;
		bool rep1 = rd1_ != NULL && nunpair1 > 0;
		bool rep2 = rd2_ != NULL && nunpair2 > 0;
		if (rep1)
		{
			summ1.init(
				rd1_, NULL, NULL, NULL, &rs1u_, NULL,
				exhaust1, exhaust2, -1, -1);
			AlnScore bestUScore, bestP1Score, bestP2Score, bestCScore;
			AlnScore bestUDist, bestP1Dist, bestP2Dist, bestCDist;
			AlnScore bestUnchosenUScore, bestUnchosenP1Score, bestUnchosenP2Score, bestUnchosenCScore;
			AlnScore bestUnchosenUDist, bestUnchosenP1Dist, bestUnchosenP2Dist, bestUnchosenCDist;
			size_t off = selectByScore(
				&rs1u_, NULL, nunpair1, select1_, NULL, NULL,
				bestUScore,
				bestUDist,
				bestP1Score,
				bestP1Dist,
				bestP2Score,
				bestP2Dist,
				bestCScore,
				bestCDist,
				bestUnchosenUScore,
				bestUnchosenUDist,
				bestUnchosenP1Score,
				bestUnchosenP1Dist,
				bestUnchosenP2Score,
				bestUnchosenP2Dist,
				bestUnchosenCScore,
				bestUnchosenCDist,
				rnd);
			summ1.setBest(
				bestUScore,
				bestUDist,
				bestP1Score,
				bestP1Dist,
				bestP2Score,
				bestP2Dist,
				bestCScore,
				bestCDist,
				bestUnchosenUScore,
				bestUnchosenUDist,
				bestUnchosenP1Score,
				bestUnchosenP1Dist,
				bestUnchosenP2Score,
				bestUnchosenP2Dist,
				bestUnchosenCScore,
				bestUnchosenCDist);
			repRs1 = &rs1u_[off];
		}
		else if (rd1_ != NULL)
		{
			assert(!unpair1Max);
		}
		if (rep2)
		{
			summ2.init(
				NULL, rd2_, NULL, NULL, NULL, &rs2u_,
				exhaust1, exhaust2, -1, -1);
			AlnScore bestUScore, bestP1Score, bestP2Score, bestCScore;
			AlnScore bestUDist, bestP1Dist, bestP2Dist, bestCDist;
			AlnScore bestUnchosenUScore, bestUnchosenP1Score, bestUnchosenP2Score, bestUnchosenCScore;
			AlnScore bestUnchosenUDist, bestUnchosenP1Dist, bestUnchosenP2Dist, bestUnchosenCDist;
			size_t off = selectByScore(
				&rs2u_, NULL, nunpair2, select2_, NULL, NULL,
				bestUScore,
				bestUDist,
				bestP1Score,
				bestP1Dist,
				bestP2Score,
				bestP2Dist,
				bestCScore,
				bestCDist,
				bestUnchosenUScore,
				bestUnchosenUDist,
				bestUnchosenP1Score,
				bestUnchosenP1Dist,
				bestUnchosenP2Score,
				bestUnchosenP2Dist,
				bestUnchosenCScore,
				bestUnchosenCDist,
				rnd);
			summ2.setBest(
				bestUScore,
				bestUDist,
				bestP1Score,
				bestP1Dist,
				bestP2Score,
				bestP2Dist,
				bestCScore,
				bestCDist,
				bestUnchosenUScore,
				bestUnchosenUDist,
				bestUnchosenP1Score,
				bestUnchosenP1Dist,
				bestUnchosenP2Score,
				bestUnchosenP2Dist,
				bestUnchosenCScore,
				bestUnchosenCDist);
			repRs2 = &rs2u_[off];
		}
		else if (rd2_ != NULL)
		{
			assert(!unpair2Max);
		}
		if (rep1)
		{
			flags1.init(
				readIsPair() ? ALN_FLAG_PAIR_UNPAIRED_MATE1 : ALN_FLAG_PAIR_UNPAIRED,
				st_.params().mhitsSet(),
				unpair1Max,
				pairMax,
				nfilt1,
				scfilt1,
				lenfilt1,
				qcfilt1,
				st_.params().mixed,
				true,
				repRs2 != NULL,
				repRs2 == NULL || repRs2->fw(),
				scUnMapped,
				xeq);
			for (size_t i = 0; i < rs1u_.size(); i++)
			{
				rs1u_[i].setMateParams(ALN_RES_TYPE_UNPAIRED_MATE1, NULL, flags1);
			}
		}
		if (rep2)
		{
			flags2.init(
				readIsPair() ? ALN_FLAG_PAIR_UNPAIRED_MATE2 : ALN_FLAG_PAIR_UNPAIRED,
				st_.params().mhitsSet(),
				unpair2Max,
				pairMax,
				nfilt2,
				scfilt2,
				lenfilt2,
				qcfilt2,
				st_.params().mixed,
				true,
				repRs1 != NULL,
				repRs1 == NULL || repRs1->fw(),
				scUnMapped,
				xeq);
			for (size_t i = 0; i < rs2u_.size(); i++)
			{
				rs2u_[i].setMateParams(ALN_RES_TYPE_UNPAIRED_MATE2, NULL, flags2);
			}
		}
		if (rep1)
		{
			SeedAlSumm ssm1, ssm2;
			if (sr1 != NULL)
				sr1->toSeedAlSumm(ssm1);
			if (sr2 != NULL)
				sr2->toSeedAlSumm(ssm2);
			assert(!select1_.empty());
			g_.reportHits(
				obuf_,
				staln_,
				threadid_,
				rd1_,
				repRs2 != NULL ? rd2_ : NULL,
				rdid_,
				select1_,
				repRs2 != NULL ? &select2_ : NULL,
				&rs1u_,
				repRs2 != NULL ? &rs2u_ : NULL,
				unpair1Max,
				summ1,
				ssm1,
				ssm2,
				&flags1,
				repRs2 != NULL ? &flags2 : NULL,
				prm,
				mapq_,
				sc,
				false);
			assert_lt(select1_[0], rs1u_.size());
			refid = rs1u_[select1_[0]].refid();
			refoff = rs1u_[select1_[0]].refoff();
		}
		if (rep2)
		{
			SeedAlSumm ssm1, ssm2;
			if (sr1 != NULL)
				sr1->toSeedAlSumm(ssm1);
			if (sr2 != NULL)
				sr2->toSeedAlSumm(ssm2);
			assert(!select2_.empty());
			g_.reportHits(
				obuf_,
				staln_,
				threadid_,
				rd2_,
				repRs1 != NULL ? rd1_ : NULL,
				rdid_,
				select2_,
				repRs1 != NULL ? &select1_ : NULL,
				&rs2u_,
				repRs1 != NULL ? &rs1u_ : NULL,
				unpair2Max,
				summ2,
				ssm1,
				ssm2,
				&flags2,
				repRs1 != NULL ? &flags1 : NULL,
				prm,
				mapq_,
				sc,
				false);
			assert_lt(select2_[0], rs2u_.size());
			refid = rs2u_[select2_[0]].refid();
			refoff = rs2u_[select2_[0]].refoff();
		}
		if (rd1_ != NULL && nunpair1 == 0)
		{
			if (nunpair2 > 0)
			{
				assert_neq(-1, refid);
				summ1.init(
					rd1_, NULL, NULL, NULL, NULL, NULL,
					exhaust1, exhaust2, refid, refoff);
			}
			else
			{
				summ1.init(
					rd1_, NULL, NULL, NULL, NULL, NULL,
					exhaust1, exhaust2, -1, -1);
			}
			SeedAlSumm ssm1, ssm2;
			if (sr1 != NULL)
				sr1->toSeedAlSumm(ssm1);
			if (sr2 != NULL)
				sr2->toSeedAlSumm(ssm2);
			flags1.init(
				readIsPair() ? ALN_FLAG_PAIR_UNPAIRED_MATE1 : ALN_FLAG_PAIR_UNPAIRED,
				st_.params().mhitsSet(),
				false,
				false,
				nfilt1,
				scfilt1,
				lenfilt1,
				qcfilt1,
				st_.params().mixed,
				true,
				repRs2 != NULL,
				(repRs2 != NULL) ? repRs2->fw() : false,
				scUnMapped,
				xeq);
			g_.reportUnaligned(
				obuf_,
				staln_,
				threadid_,
				rd1_,
				NULL,
				rdid_,
				summ1,
				ssm1,
				ssm2,
				&flags1,
				NULL,
				prm,
				mapq_,
				sc,
				true);
		}
		if (rd2_ != NULL && nunpair2 == 0)
		{
			if (nunpair1 > 0)
			{
				assert_neq(-1, refid);
				summ2.init(
					NULL, rd2_, NULL, NULL, NULL, NULL,
					exhaust1, exhaust2, refid, refoff);
			}
			else
			{
				summ2.init(
					NULL, rd2_, NULL, NULL, NULL, NULL,
					exhaust1, exhaust2, -1, -1);
			}
			SeedAlSumm ssm1, ssm2;
			if (sr1 != NULL)
				sr1->toSeedAlSumm(ssm1);
			if (sr2 != NULL)
				sr2->toSeedAlSumm(ssm2);
			flags2.init(
				readIsPair() ? ALN_FLAG_PAIR_UNPAIRED_MATE2 : ALN_FLAG_PAIR_UNPAIRED,
				st_.params().mhitsSet(),
				false,
				false,
				nfilt2,
				scfilt2,
				lenfilt2,
				qcfilt2,
				st_.params().mixed,
				true,
				repRs1 != NULL,
				(repRs1 != NULL) ? repRs1->fw() : false,
				scUnMapped,
				xeq);
			g_.reportUnaligned(
				obuf_,
				staln_,
				threadid_,
				rd2_,
				NULL,
				rdid_,
				summ2,
				ssm1,
				ssm2,
				&flags2,
				NULL,
				prm,
				mapq_,
				sc,
				true);
		}
	}
	init_ = false;
	return;
}
bool AlnSinkWrap::report(
	int stage,
	const AlnRes *rs1,
	const AlnRes *rs2)
{
	assert(init_);
	assert(rs1 != NULL || rs2 != NULL);
	assert(rs1 == NULL || !rs1->empty());
	assert(rs2 == NULL || !rs2->empty());
	assert(rs1 == NULL || rs1->repOk());
	assert(rs2 == NULL || rs2->repOk());
	bool paired = (rs1 != NULL && rs2 != NULL);
	bool one = (rs1 != NULL);
	const AlnRes *rsa = one ? rs1 : rs2;
	const AlnRes *rsb = one ? rs2 : rs1;
	if (paired)
	{
		assert(readIsPair());
		st_.foundConcordant();
		rs1_.push_back(*rs1);
		rs2_.push_back(*rs2);
	}
	else
	{
		st_.foundUnpaired(one);
		if (one)
		{
			rs1u_.push_back(*rs1);
		}
		else
		{
			rs2u_.push_back(*rs2);
		}
	}
	TAlScore score = rsa->score().score();
	if (rsb != NULL)
		score += rsb->score().score();
	if (paired)
	{
		if (score > bestPair_)
		{
			best2Pair_ = bestPair_;
			bestPair_ = score;
		}
		else if (score > best2Pair_)
		{
			best2Pair_ = score;
		}
	}
	else
	{
		if (one)
		{
			if (score > bestUnp1_)
			{
				best2Unp1_ = bestUnp1_;
				bestUnp1_ = score;
			}
			else if (score > best2Unp1_)
			{
				best2Unp1_ = score;
			}
		}
		else
		{
			if (score > bestUnp2_)
			{
				best2Unp2_ = bestUnp2_;
				bestUnp2_ = score;
			}
			else if (score > best2Unp2_)
			{
				best2Unp2_ = score;
			}
		}
	}
	return st_.done();
}
bool AlnSinkWrap::prepareDiscordants()
{
	if (rs1u_.size() == 1 && rs2u_.size() == 1)
	{
		assert(rs1_.empty());
		assert(rs2_.empty());
		rs1_.push_back(rs1u_[0]);
		rs2_.push_back(rs2u_[0]);
		return true;
	}
	return false;
}
size_t AlnSinkWrap::selectByScore(
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
	const
{
	assert(init_);
	assert(repOk());
	assert_gt(num, 0);
	assert(rs1 != NULL);
	assert(rs2 == NULL || rs1u != NULL);
	assert(rs2 == NULL || rs2u != NULL);
	bestUScore.invalidate();
	bestUDist.invalidate();
	bestUnchosenUScore.invalidate();
	bestUnchosenUDist.invalidate();
	bestCScore.invalidate();
	bestP1Score.invalidate();
	bestP2Score.invalidate();
	bestCDist.invalidate();
	bestP1Dist.invalidate();
	bestP2Dist.invalidate();
	bestUnchosenCScore.invalidate();
	bestUnchosenP1Score.invalidate();
	bestUnchosenP2Score.invalidate();
	bestUnchosenCDist.invalidate();
	bestUnchosenP1Dist.invalidate();
	bestUnchosenP2Dist.invalidate();
	size_t sz = rs1->size();
	assert_leq(num, sz);
	if (sz < num)
	{
		num = sz;
	}
	if (sz == 0)
	{
		return 0;
	}
	select.resize((size_t)num);
	EList<std::pair<AlnScore, size_t>> &buf =
		const_cast<EList<std::pair<AlnScore, size_t>> &>(selectBuf_);
	buf.resize(sz);
	for (size_t i = 0; i < sz; i++)
	{
		buf[i].first = (*rs1)[i].score();
		if (rs2 != NULL)
		{
			buf[i].first += (*rs2)[i].score();
		}
		buf[i].second = i;
	}
	buf.sort();
	buf.reverse();
	size_t streak = 0;
	for (size_t i = 1; i < buf.size(); i++)
	{
		if (buf[i].first == buf[i - 1].first)
		{
			if (streak == 0)
			{
				streak = 1;
			}
			streak++;
		}
		else
		{
			if (streak > 1)
			{
				assert_geq(i, streak);
				buf.shufflePortion(i - streak, streak, rnd);
			}
			streak = 0;
		}
	}
	if (streak > 1)
	{
		buf.shufflePortion(buf.size() - streak, streak, rnd);
	}
	for (size_t i = 0; i < num; i++)
	{
		select[i] = buf[i].second;
	}
	if (rs2 == NULL)
	{
		bestUScore = bestUDist = (*rs1)[select[0]].score();
	}
	if (rs2 != NULL)
	{
		bestCScore = bestCDist = (*rs1)[select[0]].score() + (*rs2)[select[0]].score();
		bestP1Score = bestP1Dist = (*rs1)[select[0]].score();
		bestP2Score = bestP2Dist = (*rs2)[select[0]].score();
		for (size_t i = 0; i < rs1u->size(); i++)
		{
			if ((*rs1u)[i].refcoord() == (*rs1)[select[0]].refcoord())
			{
				continue;
			}
			if ((*rs1u)[i].score() > bestUnchosenP1Score)
			{
				bestUnchosenP1Score = (*rs1u)[i].score();
			}
			if ((*rs1u)[i].score().basesAligned() > bestUnchosenP1Dist.basesAligned())
			{
				bestUnchosenP1Dist = (*rs1u)[i].score();
			}
		}
		for (size_t i = 0; i < rs2u->size(); i++)
		{
			if ((*rs2u)[i].refcoord() == (*rs2)[select[0]].refcoord())
			{
				continue;
			}
			if ((*rs2u)[i].score() > bestUnchosenP2Score)
			{
				bestUnchosenP2Score = (*rs2u)[i].score();
			}
			if ((*rs2u)[i].score().basesAligned() > bestUnchosenP2Dist.basesAligned())
			{
				bestUnchosenP2Dist = (*rs2u)[i].score();
			}
		}
		if (buf.size() > 1)
		{
			bestUnchosenCScore = buf[1].first;
			for (size_t i = 1; i < buf.size(); i++)
			{
				AlnScore dist = (*rs1)[buf[i].second].score() +
								(*rs2)[buf[i].second].score();
				if (dist.basesAligned() > bestUnchosenCDist.basesAligned())
				{
					bestUnchosenCDist = dist;
				}
			}
		}
	}
	else if (buf.size() > 1)
	{
		bestUnchosenUScore = (*rs1)[buf[1].second].score();
		for (size_t i = 1; i < buf.size(); i++)
		{
			if ((*rs1)[buf[1].second].score().basesAligned() > bestUnchosenUDist.basesAligned())
			{
				bestUnchosenUDist = (*rs1)[buf[1].second].score();
			}
		}
	}
	return selectBuf_[0].second;
}
size_t AlnSinkWrap::selectAlnsToReport(
	const EList<AlnRes> &rs,
	uint64_t num,
	EList<size_t> &select,
	RandomSource &rnd)
	const
{
	assert(init_);
	assert(repOk());
	assert_gt(num, 0);
	size_t sz = rs.size();
	if (sz < num)
	{
		num = sz;
	}
	if (sz < 1)
	{
		return 0;
	}
	select.resize((size_t)num);
	if (sz == 1)
	{
		assert_eq(1, num);
		select[0] = 0;
		return 0;
	}
	uint32_t off = rnd.nextU32() % (uint32_t)sz;
	uint32_t offOrig = off;
	for (size_t i = 0; i < num; i++)
	{
		select[i] = off;
		off++;
		if (off == sz)
		{
			off = 0;
		}
	}
	return offOrig;
}
#define NOT_SUPPRESSED !suppress_[field++]
#define BEGIN_FIELD             \
	{                           \
		if (firstfield)         \
			firstfield = false; \
		else                    \
			o.append('\t');     \
	}
#define WRITE_TAB               \
	{                           \
		if (firstfield)         \
			firstfield = false; \
		else                    \
			o.append('\t');     \
	}
#define WRITE_NUM(o, x) \
	{                   \
		itoa10(x, buf); \
		o.append(buf);  \
	}
void AlnSink::reportSeedSummary(
	BTString &o,
	const Read &rd,
	TReadId rdid,
	size_t threadId,
	const SeedResults &rs,
	bool getLock)
{
	appendSeedSummary(
		o,
		rd,
		rdid,
		rs.numOffs() * 2,
		rs.nonzeroOffsets(),
		rs.numRanges(),
		rs.numElts(),
		rs.numOffs(),
		rs.nonzeroOffsetsFw(),
		rs.numRangesFw(),
		rs.numEltsFw(),
		rs.numOffs(),
		rs.nonzeroOffsetsRc(),
		rs.numRangesRc(),
		rs.numEltsRc());
}
void AlnSink::reportEmptySeedSummary(
	BTString &o,
	const Read &rd,
	TReadId rdid,
	size_t threadId,
	bool getLock)
{
	appendSeedSummary(
		o,
		rd,
		rdid,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0);
}
template <typename T>
static inline void printUptoWs(
	BTString &s,
	const T &str,
	bool chopws)
{
	size_t len = str.length();
	for (size_t i = 0; i < len; i++)
	{
		if (!chopws || (str[i] != ' ' && str[i] != '\t'))
		{
			s.append(str[i]);
		}
		else
		{
			break;
		}
	}
}
void AlnSink::appendSeedSummary(
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
	size_t eltsRc)
{
	char buf[1024];
	bool firstfield = true;
	BEGIN_FIELD;
	printUptoWs(o, rd.name, true);
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTried);
	BEGIN_FIELD;
	WRITE_NUM(o, nonzero);
	BEGIN_FIELD;
	WRITE_NUM(o, ranges);
	BEGIN_FIELD;
	WRITE_NUM(o, elts);
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTriedFw);
	BEGIN_FIELD;
	WRITE_NUM(o, nonzeroFw);
	BEGIN_FIELD;
	WRITE_NUM(o, rangesFw);
	BEGIN_FIELD;
	WRITE_NUM(o, eltsFw);
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTriedRc);
	BEGIN_FIELD;
	WRITE_NUM(o, nonzeroRc);
	BEGIN_FIELD;
	WRITE_NUM(o, rangesRc);
	BEGIN_FIELD;
	WRITE_NUM(o, eltsRc);
	o.append('\n');
}
void AlnSinkSam::appendMate(
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
	const Mapq &mapqCalc,
	const Scoring &sc)
{
	if (rs == NULL && samc_.omitUnalignedReads())
	{
		return;
	}
	char buf[1024];
	char mapqInps[1024];
	if (rs != NULL)
	{
		staln.reset();
		rs->initStacked(rd, staln);
		staln.leftAlign(false);
	}
	int offAdj = 0;
	samc_.printReadName(o, rd.name, flags.partOfPair());
	o.append('\t');
	int fl = 0;
	if (flags.partOfPair())
	{
		fl |= SAM_FLAG_PAIRED;
		if (flags.alignedConcordant())
		{
			fl |= SAM_FLAG_MAPPED_PAIRED;
		}
		if (!flags.mateAligned())
		{
			fl |= SAM_FLAG_MATE_UNMAPPED;
		}
		fl |= (flags.readMate1() ? SAM_FLAG_FIRST_IN_PAIR : SAM_FLAG_SECOND_IN_PAIR);
		if (flags.mateAligned())
		{
			bool oppFw = (rso != NULL) ? rso->fw() : flags.isOppFw();
			if (!oppFw)
			{
				fl |= SAM_FLAG_MATE_STRAND;
			}
		}
	}
	if (!flags.isPrimary())
	{
		fl |= SAM_FLAG_NOT_PRIMARY;
	}
	if (rs != NULL && !rs->fw())
	{
		fl |= SAM_FLAG_QUERY_STRAND;
	}
	if (rs == NULL)
	{
		fl |= SAM_FLAG_UNMAPPED;
	}
	itoa10<int>(fl, buf);
	o.append(buf);
	o.append('\t');
	if (rs != NULL)
	{
		samc_.printRefNameFromIndex(o, (size_t)rs->refid());
		o.append('\t');
	}
	else
	{
		if (summ.orefid() != -1)
		{
			assert(flags.partOfPair());
			samc_.printRefNameFromIndex(o, (size_t)summ.orefid());
		}
		else
		{
			o.append('*');
		}
		o.append('\t');
	}
	if (rs != NULL)
	{	
		if(rs->refoff() + 1 + offAdj<0){
			itoa10<int64_t>(0, buf);
		}else{
			itoa10<int64_t>(rs->refoff() + 1 + offAdj, buf);
		}
		o.append(buf);
		o.append('\t');
	}
	else
	{
		if (summ.orefid() != -1)
		{
			assert(flags.partOfPair());
			itoa10<int64_t>(summ.orefoff() + 1 + offAdj, buf);
			o.append(buf);
		}
		else
		{
			o.append('0');
		}
		o.append('\t');
	}
	mapqInps[0] = '\0';
	if (rs != NULL)
	{
		itoa10<TMapq>(mapqCalc.mapq(
						  summ, flags, rd.mate < 2, rd.length(),
						  rdo == NULL ? 0 : rdo->length(), mapqInps),
					  buf);
		o.append(buf);
		o.append('\t');
	}
	else
	{
		o.append("0\t");
	}
	if (rs != NULL)
	{	
		if(rs->mycigar.empty()){
			int len=rd.length();
			o.append((to_string(len)+"M").c_str());
		}else{
			o.append((rs->mycigar).c_str());
		}
		o.append('\t');
	}
	else
	{
		o.append("*\t");
	}
	if (rs != NULL && flags.partOfPair())
	{
		if (rso != NULL && rs->refid() != rso->refid())
		{
			samc_.printRefNameFromIndex(o, (size_t)rso->refid());
			o.append('\t');
		}
		else
		{
			o.append("=\t");
		}
	}
	else if (summ.orefid() != -1)
	{
		o.append("=\t");
	}
	else
	{
		o.append("*\t");
	}
	if (rs != NULL && flags.partOfPair())
	{
		if (rso != NULL)
		{
			itoa10<int64_t>(rso->refoff() + 1, buf);
			o.append(buf);
			o.append('\t');
		}
		else
		{
			itoa10<int64_t>(rs->refoff() + 1, buf);
			o.append(buf);
			o.append('\t');
		}
	}
	else if (summ.orefid() != -1)
	{
		itoa10<int64_t>(summ.orefoff() + 1, buf);
		o.append(buf);
		o.append('\t');
	}
	else
	{
		o.append("0\t");
	}
	if (rs != NULL && rs->isFraglenSet())
	{
		itoa10<int64_t>(rs->fragmentLength(), buf);
		o.append(buf);
		o.append('\t');
	}
	else
	{
		o.append("0\t");
	}
	if (!flags.isPrimary() && samc_.omitSecondarySeqQual())
	{
		o.append('*');
	}
	else
	{
		if (rd.patFw.length() == 0)
		{
			o.append('*');
		}
		else
		{
			if (rs == NULL || rs->fw())
			{
				o.append(rd.patFw.toZBuf());
			}
			else
			{
				o.append(rd.patRc.toZBuf());
			}
		}
	}
	o.append('\t');
	if (!flags.isPrimary() && samc_.omitSecondarySeqQual())
	{
		o.append('*');
	}
	else
	{
		if (rd.qual.length() == 0)
		{
			o.append('*');
		}
		else
		{
			if (rs == NULL || rs->fw())
			{
				o.append(rd.qual.toZBuf());
			}
			else
			{
				o.append(rd.qualRev.toZBuf());
			}
		}
	}
	o.append('\t');
	o.append('\n');
}
#ifdef ALN_SINK_MAIN
#include <iostream>
bool testDones(
	const ReportingState &st,
	bool done1,
	bool done2,
	bool done3,
	bool done4,
	bool done5,
	bool done6)
{
	assert(st.doneConcordant() == done1);
	assert(st.doneDiscordant() == done2);
	assert(st.doneUnpaired(true) == done3);
	assert(st.doneUnpaired(false) == done4);
	assert(st.doneUnpaired() == done5);
	assert(st.done() == done6);
	assert(st.repOk());
	return true;
}
int main(void)
{
	cerr << "Case 1 (simple unpaired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,
			0,
			0,
			false,
			false,
			false);
		ReportingState st(rp);
		st.nextRead(false);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, true, true, true, true));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(2, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;
	cerr << "Case 2 (simple unpaired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,
			3,
			0,
			false,
			false,
			false);
		ReportingState st(rp);
		st.nextRead(false);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;
	cerr << "Case 3 (simple paired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,
			3,
			0,
			false,
			false,
			false);
		ReportingState st(rp);
		st.nextRead(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;
	cerr << "Case 4 (simple paired 2) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,
			3,
			0,
			false,
			true,
			true);
		ReportingState st(rp);
		st.nextRead(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(unpair1Max);
		assert(unpair2Max);
	}
	cerr << "PASSED" << endl;
	cerr << "Case 5 (potential discordant after concordant) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,
			3,
			0,
			false,
			true,
			true);
		ReportingState st(rp);
		st.nextRead(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		st.foundUnpaired(false);
		st.foundConcordant();
		assert(testDones(st, false, true, false, false, false, false));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(1, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(1, st.numUnpaired1());
		assert_eq(1, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(1, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;
	cerr << "Case 6 (true discordant) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,
			3,
			0,
			false,
			true,
			true);
		ReportingState st(rp);
		st.nextRead(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		st.foundUnpaired(false);
		assert(testDones(st, false, false, false, false, false, false));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(1, st.numDiscordant());
		assert_eq(0, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(1, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;
	cerr << "Case 7 (unaligned pair & uniquely aligned mate, mixed-mode) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,
			1,
			0,
			false,
			true,
			true);
		ReportingState st(rp);
		st.nextRead(true);
		st.foundUnpaired(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, false, false, false));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;
	cerr << "Case 8 (unaligned pair & uniquely aligned mate, NOT mixed-mode) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,
			1,
			0,
			false,
			true,
			false);
		ReportingState st(rp);
		st.nextRead(true);
		st.foundUnpaired(true);
		assert(testDones(st, false, false, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;
	cerr << "Case 9 (repetitive pair, only one mate repetitive) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,
			1,
			0,
			true,
			true,
			true);
		ReportingState st(rp);
		st.nextRead(true);
		st.foundConcordant();
		assert(st.repOk());
		st.foundUnpaired(true);
		assert(st.repOk());
		st.foundUnpaired(false);
		assert(st.repOk());
		assert(testDones(st, false, true, false, false, false, false));
		assert(st.repOk());
		st.foundConcordant();
		assert(st.repOk());
		st.foundUnpaired(true);
		assert(st.repOk());
		assert(testDones(st, true, true, true, false, false, false));
		assert_eq(2, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(1, st.numUnpaired2());
		st.foundUnpaired(false);
		assert(st.repOk());
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(2, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(2, st.numUnpaired2());
		st.finish();
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
					 pairMax, unpair1Max, unpair2Max);
		assert_eq(1, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(unpair1Max);
		assert(unpair2Max);
	}
	cerr << "PASSED" << endl;
}
#endif
