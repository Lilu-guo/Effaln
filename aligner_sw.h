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
extern unsigned long long int countAlign_Static;
extern unsigned long long int countGather_Static;
extern unsigned long long int align_time;
extern unsigned long long int gather_time;
extern unsigned long long int time_multiseedSearchWorker;
extern unsigned long long int time_effaln;
extern unsigned long long int time_exact;
extern unsigned long long int time_1mm;
extern unsigned long long int time_seed;
extern unsigned long long int time_extendSeeds;
extern unsigned long long int time_seed_extendSeeds;
extern unsigned long long int time_seed_searchAllSeeds;
extern unsigned long long int time_seed_mmSeeds;
extern unsigned long long int time_seed_instantiateSeeds;
extern unsigned long long int time_extendSeeds_eeSaTups;
extern unsigned long long int time_extendSeeds_prioritizeSATups;
extern unsigned long long int time_extendSeeds_advanceElement;
extern unsigned long long int time_extendSeeds_joinedToTextOff;
extern unsigned long long int redundants_Nums;
extern unsigned long long int time_extendSeeds_resEe;
extern unsigned long long int time_extendSeeds_resUngap;
extern unsigned long long int time_extendSeeds_init;
extern unsigned long long int time_extendSeeds_align;
extern unsigned long long int bestFilter_Nums;
extern unsigned long long int time_wfa0;
extern unsigned long long int time_wfa1;
extern unsigned long long int seed_prune;
extern uint8_t code[256];
extern char rcsymbol[6];
#ifndef ALIGNER_SW_H_
#define ALIGNER_SW_H_
#define INLINE_CUPS
#include <stdint.h>
#include <iostream>
#include <limits>
#include "threading.h"
#include "sse_wrap.h"
#include "aligner_sw_common.h"
#include "aligner_sw_nuc.h"
#include "ds.h"
#include "aligner_seed.h"
#include "reference.h"
#include "random_source.h"
#include "mem_ids.h"
#include "aligner_result.h"
#include "mask.h"
#include "dp_framer.h"
#include "aligner_swsse.h"
#include "aligner_bt.h"
#include "embedding.h"
#define QUAL2(d, f) sc_->mm((int)(*rd_)[rdi_ + d], \
							(int)rf_[rfi_ + f],    \
							(int)(*qu_)[rdi_ + d] - 33)
#define QUAL(d) sc_->mm((int)(*rd_)[rdi_ + d], \
						(int)(*qu_)[rdi_ + d] - 33)
#define N_SNP_PEN(c) (((int)rf_[rfi_ + c] > 15) ? sc_->n(30) : sc_->penSnp)
class SwAligner
{
	typedef std::pair<size_t, size_t> SizeTPair;
	enum
	{
		STATE_UNINIT,
		STATE_INITED,
		STATE_ALIGNED,
	};
	const static size_t ALPHA_SIZE = 5;
public:
	explicit SwAligner(std::ostream *dpLog, bool firstRead = true) : sseU8fw_(DP_CAT),
																	 sseU8rc_(DP_CAT),
																	 sseI16fw_(DP_CAT),
																	 sseI16rc_(DP_CAT),
																	 state_(STATE_UNINIT),
																	 initedRead_(false),
																	 readSse16_(false),
																	 initedRef_(false),
																	 rfwbuf_(DP_CAT),
																	 btnstack_(DP_CAT),
																	 btcells_(DP_CAT),
																	 btdiag_(),
																	 btncand_(DP_CAT),
																	 btncanddone_(DP_CAT),
																	 btncanddoneSucc_(0),
																	 btncanddoneFail_(0),
																	 cper_(),
																	 cperMinlen_(),
																	 cperPerPow2_(),
																	 cperEf_(),
																	 cperTri_(),
																	 colstop_(0),
																	 lastsolcol_(0),
																	 cural_(0),
																	 dpLog_(dpLog),
																	 firstRead_(firstRead)
																		 ASSERT_ONLY(, cand_tmp_(DP_CAT))
	{
		int loop = CGK_LOOP;
		mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		affine_penalties = {
			.match = 0,
			.mismatch = 3,
			.gap_opening = 5,
			.gap_extension = 3,
		};
		embedding = new Embedding();
		bianma = new char *[3];
		for (int i = 0; i < 3; i++)
		{
			bianma[i] = new char[MAX_ELEN + 1];
		}
		embed_rd = false;
		for (int i = 0; i < loop; i++)
		{
			embeddedQ[i] = new char *[2];
			for (int j = 0; j < 2; j++)
			{
				embeddedQ[i][j] = new char[MAX_ELEN + 1];
			}
		}
	}
	void initRead(
		const BTDnaString &rdfw,
		const BTDnaString &rdrc,
		const BTString &qufw,
		const BTString &qurc,
		size_t rdi,
		size_t rdf,
		const Scoring &sc);
	void initRef(
		bool fw,
		TRefId refidx,
		const DPRect &rect,
		char *rf,
		size_t rfi,
		size_t rff,
		TRefOff reflen,
		const Scoring &sc,
		TAlScore minsc,
		bool enable8,
		size_t cminlen,
		size_t cpow2,
		bool doTri,
		bool extend);
	void initRef(
		bool fw,
		TRefId refidx,
		const DPRect &rect,
		const BitPairReference &refs,
		TRefOff reflen,
		const Scoring &sc,
		TAlScore minsc,
		bool enable8,
		size_t cminlen,
		size_t cpow2,
		bool doTri,
		bool extend,
		size_t upto,
		size_t &nsUpto);
	int ungappedAlign(
		const BTDnaString &rd,
		const BTString &qu,
		const Coord &coord,
		const BitPairReference &refs,
		size_t reflen,
		const Scoring &sc,
		bool ohang,
		TAlScore minsc,
		SwResult &res);
	bool align(TAlScore &best, SwResult &res, const Coord &coord);
	bool nextAlignment(
		SwResult &res,
		TAlScore minsc,
		RandomSource &rnd);
	void printResultStacked(
		const SwResult &res,
		std::ostream &os)
	{
		res.alres.printStacked(*rd_, os);
	}
	bool done() const
	{
		assert(initedRead() && initedRef());
		return cural_ == btncand_.size();
	}
	inline bool initedRef() const { return initedRef_; }
	inline bool initedRead() const { return initedRead_; }
	inline void reset()
	{
		initedRef_ = initedRead_ = false;
		embed_rd = false;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_gt(dpRows(), 0);
		for (size_t i = 0; i < btncand_.size(); i++)
		{
			assert(btncand_[i].repOk());
			assert_geq(btncand_[i].score, minsc_);
		}
		return true;
	}
#endif
	size_t numAlignmentsReported() const
	{
		return cural_;
	}
	void merge(
		SSEMetrics &sseU8ExtendMet,
		SSEMetrics &sseU8MateMet,
		SSEMetrics &sseI16ExtendMet,
		SSEMetrics &sseI16MateMet,
		uint64_t &nbtfiltst,
		uint64_t &nbtfiltsc,
		uint64_t &nbtfiltdo)
	{
		sseU8ExtendMet.merge(sseU8ExtendMet_);
		sseU8MateMet.merge(sseU8MateMet_);
		sseI16ExtendMet.merge(sseI16ExtendMet_);
		sseI16MateMet.merge(sseI16MateMet_);
		nbtfiltst += nbtfiltst_;
		nbtfiltsc += nbtfiltsc_;
		nbtfiltdo += nbtfiltdo_;
	}
	void resetCounters()
	{
		sseU8ExtendMet_.reset();
		sseU8MateMet_.reset();
		sseI16ExtendMet_.reset();
		sseI16MateMet_.reset();
		nbtfiltst_ = nbtfiltsc_ = nbtfiltdo_ = 0;
	}
	size_t size() const
	{
		return dpRows() * (rff_ - rfi_);
	}
	inline size_t dpRows() const
	{
		assert(initedRead_);
		return rdf_ - rdi_;
	}
	TAlScore alignNucleotidesEnd2EndSseU8(
		int &flag, bool debug);
	TAlScore alignNucleotidesLocalSseU8(
		int &flag, bool debug);
	TAlScore alignNucleotidesEnd2EndSseI16(
		int &flag, bool debug);
	TAlScore alignNucleotidesLocalSseI16(
		int &flag, bool debug);
	TAlScore alignGatherEE8(
		int &flag, bool debug);
	TAlScore alignGatherLoc8(
		int &flag, bool debug);
	TAlScore alignGatherEE16(
		int &flag, bool debug);
	TAlScore alignGatherLoc16(
		int &flag, bool debug);
	void buildQueryProfileEnd2EndSseU8(bool fw);
	void buildQueryProfileLocalSseU8(bool fw);
	void buildQueryProfileEnd2EndSseI16(bool fw);
	void buildQueryProfileLocalSseI16(bool fw);
	bool gatherCellsNucleotidesLocalSseU8(TAlScore best);
	bool gatherCellsNucleotidesEnd2EndSseU8(TAlScore best);
	bool gatherCellsNucleotidesLocalSseI16(TAlScore best);
	bool gatherCellsNucleotidesEnd2EndSseI16(TAlScore best);
	bool backtraceNucleotidesLocalSseU8(
		TAlScore escore,
		SwResult &res,
		size_t &off,
		size_t &nbts,
		size_t row,
		size_t col,
		RandomSource &rand);
	bool backtraceNucleotidesLocalSseI16(
		TAlScore escore,
		SwResult &res,
		size_t &off,
		size_t &nbts,
		size_t row,
		size_t col,
		RandomSource &rand);
	bool backtraceNucleotidesEnd2EndSseU8(
		TAlScore escore,
		SwResult &res,
		size_t &off,
		size_t &nbts,
		size_t row,
		size_t col,
		RandomSource &rand);
	bool backtraceNucleotidesEnd2EndSseI16(
		TAlScore escore,
		SwResult &res,
		size_t &off,
		size_t &nbts,
		size_t row,
		size_t col,
		RandomSource &rand);
	bool backtrace(
		TAlScore escore,
		bool fill,
		bool usecp,
		SwResult &res,
		size_t &off,
		size_t row,
		size_t col,
		size_t maxiter,
		size_t &niter,
		RandomSource &rnd)
	{
		bter_.initBt(
			escore,
			row,
			col,
			fill,
			usecp,
			cperTri_,
			rnd);
		assert(bter_.inited());
		size_t nrej = 0;
		if (bter_.emptySolution())
		{
			return false;
		}
		else
		{
			return bter_.nextAlignment(maxiter, res, off, nrej, niter, rnd);
		}
	}
	string cigarformat(string rawCigar, int& rdlen){
		int len=rawCigar.length();
		int nD=count(rawCigar.begin(),rawCigar.end(),'D'); 
		rdlen=rdlen+nD; 
		if(len<rdlen){ 
			for(int i=0; i<(rdlen-len); i++){
				rawCigar.push_back('M');
			}
		}
		string res="";
		char last=0;
		int n=0;
		for(int i=0; i<rdlen; i++){
			char c=rawCigar[i];
			if(c=='X') c='M';
			if(c!=last && i!=0){
				res.append(to_string(n));
				res.push_back(last);
				n=0; 
			}
			n++;
			if(i==rdlen-1 && n!=0){
				res.append(to_string(n));
				res.push_back(c);
			}
			last=c;
		}
		return res;
	}
	Embedding *embedding;
	char **embeddedQ[CGK_LOOP];
	char **bianma;
	char *rev;
	char *rev_str;
	bool embed_rd;
	mm_allocator_t *mm_allocator;
	affine_penalties_t affine_penalties;
	const BTDnaString *rd_;
	const BTString *qu_;
	const BTDnaString *rdfw_;
	const BTDnaString *rdrc_;
	const BTString *qufw_;
	const BTString *qurc_;
	TReadOff rdi_;
	TReadOff rdf_;
	bool fw_;
	TRefId refidx_;
	TRefOff reflen_;
	const DPRect *rect_;
	char *rf_;
	TRefOff rfi_;
	TRefOff rff_;
	size_t rdgap_;
	size_t rfgap_;
	bool enable8_;
	bool extend_;
	const Scoring *sc_;
	TAlScore minsc_;
	int nceil_;
	bool sse8succ_;
	bool sse16succ_;
	SSEData sseU8fw_;
	SSEData sseU8rc_;
	SSEData sseI16fw_;
	SSEData sseI16rc_;
	bool sseU8fwBuilt_;
	bool sseU8rcBuilt_;
	bool sseI16fwBuilt_;
	bool sseI16rcBuilt_;
	SSEMetrics sseU8ExtendMet_;
	SSEMetrics sseU8MateMet_;
	SSEMetrics sseI16ExtendMet_;
	SSEMetrics sseI16MateMet_;
	int state_;
	bool initedRead_;
	bool readSse16_;
	bool initedRef_;
	EList<uint32_t> rfwbuf_;
	EList<DpNucFrame> btnstack_;
	EList<SizeTPair> btcells_;
	NBest<DpBtCandidate> btdiag_;
	EList<DpBtCandidate> btncand_;
	EList<DpBtCandidate> btncanddone_;
	size_t btncanddoneSucc_;
	size_t btncanddoneFail_;
	BtBranchTracer bter_;
	Checkpointer cper_;
	size_t cperMinlen_;
	size_t cperPerPow2_;
	bool cperEf_;
	bool cperTri_;
	size_t colstop_;
	size_t lastsolcol_;
	size_t cural_;
	uint64_t nbtfiltst_;
	uint64_t nbtfiltsc_;
	uint64_t nbtfiltdo_;
	std::ostream *dpLog_;
	bool firstRead_;
	ASSERT_ONLY(SStringExpandable<uint32_t> tmp_destU32_);
	ASSERT_ONLY(BTDnaString tmp_editstr_, tmp_refstr_);
	ASSERT_ONLY(EList<DpBtCandidate> cand_tmp_);
};
#endif
