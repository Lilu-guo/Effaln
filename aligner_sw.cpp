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
#ifndef myNDEBUG
#define jj printf
#else
#define jj
#endif
#include <limits>
#include "aligner_sw.h"
#include "aligner_result.h"
#include "search_globals.h"
#include "scoring.h"
#include "mask.h"
#include "embedding.h"
#include "aln_search.h"
static void parse(char seq[], char fwd[])
{
	unsigned len = strlen(seq);
	for (size_t i = 0; i < len; i++)
	{
		uint8_t c = *(code + seq[i]);
		fwd[i] = c;
	}
	*(fwd + len) = '\0';
}
void SwAligner::initRead(
	const BTDnaString &rdfw,
	const BTDnaString &rdrc,
	const BTString &qufw,
	const BTString &qurc,
	size_t rdi,
	size_t rdf,
	const Scoring &sc)
{
	assert_gt(rdf, rdi);
	int nceil = sc.nCeil.f<int>((double)rdfw.length());
	rdfw_ = &rdfw;
	rdrc_ = &rdrc;
	qufw_ = &qufw;
	qurc_ = &qurc;
	rdi_ = rdi;
	rdf_ = rdf;
	sc_ = &sc;
	nceil_ = nceil;
	readSse16_ = false;
	initedRead_ = true;
#ifndef NO_SSE
	sseU8fwBuilt_ = false;
	sseU8rcBuilt_ = false;
	sseI16fwBuilt_ = false;
	sseI16rcBuilt_ = false;
#endif
	if (dpLog_ != NULL)
	{
		if (!firstRead_)
		{
			(*dpLog_) << '\n';
		}
		(*dpLog_) << rdfw.toZBuf() << '\t' << qufw.toZBuf();
	}
	firstRead_ = false;
}
void SwAligner::initRef(
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
	bool extend)
{
	size_t readGaps = sc.maxReadGaps(minsc, rdfw_->length());
	size_t refGaps = sc.maxRefGaps(minsc, rdfw_->length());
	assert_geq(readGaps, 0);
	assert_geq(refGaps, 0);
	assert_gt(rff, rfi);
	rdgap_ = readGaps;
	rfgap_ = refGaps;
	state_ = STATE_INITED;
	fw_ = fw;
	rd_ = fw ? rdfw_ : rdrc_;
	qu_ = fw ? qufw_ : qurc_;
	refidx_ = refidx;
	rf_ = rf;
	rfi_ = rfi;
	rff_ = rff;
	reflen_ = reflen;
	rect_ = &rect;
	minsc_ = minsc;
	cural_ = 0;
	initedRef_ = true;
	enable8_ = enable8;
	extend_ = extend;
	cperMinlen_ = cminlen;
	cperPerPow2_ = cpow2;
	cperEf_ = true;
	cperTri_ = doTri;
	bter_.initRef(
		fw_ ? rdfw_->buf() : rdrc_->buf(),
		fw_ ? qufw_->buf() : qurc_->buf(),
		rd_->length(),
		rf_ + rfi_,
		rff_ - rfi_,
		reflen,
		refidx_,
		rfi_,
		fw_,
		rect_,
		&cper_,
		*sc_,
		nceil_);
	if (dpLog_ != NULL)
	{
		(*dpLog_) << '\t';
		(*dpLog_) << refidx_ << ',';
		(*dpLog_) << reflen_ << ',';
		(*dpLog_) << minsc_ << ',';
		(*dpLog_) << (fw ? '+' : '-') << ',';
		rect_->write(*dpLog_);
		(*dpLog_) << ',';
		for (TRefOff i = rfi_; i < rff_; i++)
		{
			(*dpLog_) << mask2dna[(int)rf[i]];
		}
	}
}
void SwAligner::initRef(
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
	size_t &nsUpto)
{
	TRefOff rfi = rect.refl;
	TRefOff rff = rect.refr + 1;
	assert_gt(rff, rfi);
	rff++;
	const size_t rflen = (size_t)(rff - rfi);
	size_t leftNs =
		(rfi >= 0 ? 0 : (size_t)std::abs(static_cast<long>(rfi)));
	leftNs = min(leftNs, rflen);
	size_t rightNs =
		(rff <= (TRefOff)reflen ? 0 : (size_t)std::abs(static_cast<long>(rff - reflen)));
	rightNs = min(rightNs, rflen);
	assert_geq(rflen, leftNs + rightNs);
	const size_t rflenInner = rflen - (leftNs + rightNs);
#ifndef NDEBUG
	bool haveRfbuf2 = false;
	EList<char> rfbuf2(rflen);
	if ((rand() % 10) == 0)
	{
		TRefOff rfii = rfi;
		for (size_t i = 0; i < rflen; i++)
		{
			if (rfii < 0 || (TRefOff)rfii >= reflen)
			{
				rfbuf2.push_back(4);
			}
			else
			{
				rfbuf2.push_back(refs.getBase(refidx, (size_t)rfii));
			}
			rfii++;
		}
		haveRfbuf2 = true;
	}
#endif
	rfwbuf_.resize((rflen + 16) / 4);
	int offset = refs.getStretch(
		rfwbuf_.ptr(),
		refidx,
		(rfi < 0) ? 0 : (size_t)rfi,
		rflenInner
			ASSERT_ONLY(, tmp_destU32_));
	assert_leq(offset, 16);
	rf_ = (char *)rfwbuf_.ptr() + offset;
	if (leftNs > 0)
	{
		for (size_t i = rflenInner; i > 0; i--)
		{
			rf_[i + leftNs - 1] = rf_[i - 1];
		}
		for (size_t i = 0; i < leftNs; i++)
		{
			rf_[i] = 4;
		}
	}
	if (rightNs > 0)
	{
		for (size_t i = 0; i < rightNs; i++)
		{
			rf_[i + leftNs + rflenInner] = 4;
		}
	}
#ifndef NDEBUG
	for (size_t i = 0; i < rflen; i++)
	{
		assert(!haveRfbuf2 || rf_[i] == rfbuf2[i]);
		assert_range(0, 4, (int)rf_[i]);
	}
#endif
	nsUpto = 0;
	for (size_t i = 0; i < rflen; i++)
	{
		if (i < upto && rf_[i] > 3)
		{
			nsUpto++;
		}
		rf_[i] = (1 << rf_[i]);
	}
	rff--;
	initRef(
		fw,
		refidx,
		rect,
		rf_,
		0,
		(size_t)(rff - rfi),
		reflen,
		sc,
		minsc,
		enable8,
		cminlen,
		cpow2,
		doTri,
		extend);
}
int SwAligner::ungappedAlign(
	const BTDnaString &rd,
	const BTString &qu,
	const Coord &coord,
	const BitPairReference &refs,
	size_t reflen,
	const Scoring &sc,
	bool ohang,
	TAlScore minsc,
	SwResult &res)
{
	const size_t len = rd.length();
	int nceil = sc.nCeil.f<int>((double)len);
	int ns = 0;
	TRefOff rfi = coord.off();
	TRefOff rff = rfi + (TRefOff)len;
	TRefId refidx = coord.ref();
	assert_gt(rff, rfi);
	size_t leftNs = 0;
	if (rfi < 0)
	{
		if (ohang)
		{
			leftNs = (size_t)(-rfi);
		}
		else
		{
			return 0;
		}
	}
	size_t rightNs = 0;
	if (rff > (TRefOff)reflen)
	{
		if (ohang)
		{
			rightNs = (size_t)(rff - (TRefOff)reflen);
		}
		else
		{
			return 0;
		}
	}
	if ((leftNs + rightNs) > (size_t)nceil)
	{
		return 0;
	}
	assert_geq(len, leftNs + rightNs);
	const size_t rflenInner = len - (leftNs + rightNs);
#ifndef NDEBUG
	bool haveRfbuf2 = false;
	EList<char> rfbuf2(len);
	if ((rand() % 10) == 0)
	{
		TRefOff rfii = rfi;
		for (size_t i = 0; i < len; i++)
		{
			if (rfii < 0 || (size_t)rfii >= reflen)
			{
				rfbuf2.push_back(4);
			}
			else
			{
				rfbuf2.push_back(refs.getBase(refidx, (size_t)rfii));
			}
			rfii++;
		}
		haveRfbuf2 = true;
	}
#endif
	rfwbuf_.resize((len + 16) / 4);
	int offset = refs.getStretch(
		rfwbuf_.ptr(),
		refidx,
		(rfi < 0) ? 0 : (size_t)rfi,
		rflenInner
			ASSERT_ONLY(, tmp_destU32_));
	assert_leq(offset, 16);
	rf_ = (char *)rfwbuf_.ptr() + offset;
	if (leftNs > 0)
	{
		for (size_t i = rflenInner; i > 0; i--)
		{
			rf_[i + leftNs - 1] = rf_[i - 1];
		}
		for (size_t i = 0; i < leftNs; i++)
		{
			rf_[i] = 4;
		}
	}
	if (rightNs > 0)
	{
		for (size_t i = 0; i < rightNs; i++)
		{
			rf_[i + leftNs + rflenInner] = 4;
		}
	}
#ifndef NDEBUG
	for (size_t i = 0; i < len; i++)
	{
		assert(!haveRfbuf2 || rf_[i] == rfbuf2[i]);
		assert_range(0, 4, (int)rf_[i]);
	}
#endif
	TAlScore score = 0;
	res.alres.reset();
	size_t rowi = 0;
	size_t rowf = len - 1;
	if (sc.monotone)
	{
		for (size_t i = 0; i < len; i++)
		{
			assert_geq(qu[i], 33);
			score += sc.score(rd[i], (int)(1 << rf_[i]), qu[i] - 33, ns);
			assert_leq(score, 0);
			if (score < minsc || ns > nceil)
			{
				return 0;
			}
		}
	}
	else
	{
		TAlScore floorsc = 0;
		TAlScore scoreMax = floorsc;
		size_t lastfloor = 0;
		rowi = MAX_SIZE_T;
		size_t sols = 0;
		for (size_t i = 0; i < len; i++)
		{
			score += sc.score(rd[i], (int)(1 << rf_[i]), qu[i] - 33, ns);
			if (score >= minsc && score >= scoreMax)
			{
				scoreMax = score;
				rowf = i;
				if (rowi != lastfloor)
				{
					rowi = lastfloor;
					sols++;
				}
			}
			if (score <= floorsc)
			{
				score = floorsc;
				lastfloor = i + 1;
			}
		}
		if (ns > nceil || scoreMax < minsc)
		{
			return 0;
		}
		if (sols > 1)
		{
			return -1;
		}
		score = scoreMax;
	}
	assert_geq(rowf, rowi);
	EList<Edit> &ned = res.alres.ned();
	size_t refns = 0;
	ASSERT_ONLY(BTDnaString refstr);
	for (size_t i = rowi; i <= rowf; i++)
	{
		ASSERT_ONLY(refstr.append((int)rf_[i]));
		if (rf_[i] > 3 || rd[i] != rf_[i])
		{
			Edit e((int)i,
				   mask2dna[1 << (int)rf_[i]],
				   "ACGTN"[(int)rd[i]],
				   EDIT_TYPE_MM);
			ned.push_back(e);
			if (rf_[i] > 3)
			{
				refns++;
			}
		}
	}
	res.alres.setScore(AlnScore(score,
								(int)(rd.length() - ned.size()),
								(int)ned.size(), ns, 0));
	assert(Edit::repOk(ned, rd));
	bool fw = coord.fw();
	assert_leq(rowf, len - 1);
	size_t trimEnd = (len - 1) - rowf;
	res.alres.setShape(
		coord.ref(),
		coord.off() + rowi,
		reflen,
		fw,
		len,
		true,
		0,
		0,
		true,
		fw ? rowi : trimEnd,
		fw ? trimEnd : rowi);
	res.alres.setRefNs(refns);
	assert(res.repOk());
#ifndef NDEBUG
	BTDnaString editstr;
	Edit::toRef(rd, ned, editstr, true, rowi, trimEnd);
	if (refstr != editstr)
	{
		cerr << "Decoded nucleotides and edits don't match reference:" << endl;
		cerr << "           score: " << res.alres.score().score() << endl;
		cerr << "           edits: ";
		Edit::print(cerr, ned);
		cerr << endl;
		cerr << "    decoded nucs: " << rd << endl;
		cerr << "     edited nucs: " << editstr << endl;
		cerr << "  reference nucs: " << refstr << endl;
		assert(0);
	}
#endif
	if (!fw)
	{
		res.alres.invertEdits();
	}
	return 1;
}
unsigned long long int countAlign_Static = 0;
unsigned long long int countGather_Static = 0;
unsigned long long int align_time = 0;
unsigned long long int gather_time = 0;
unsigned long long int time_multiseedSearchWorker = 0;
unsigned long long int time_effaln = 0;
unsigned long long int time_exact = 0;
unsigned long long int time_1mm = 0;
unsigned long long int time_seed = 0;
unsigned long long int time_extendSeeds = 0;
unsigned long long int time_seed_extendSeeds = 0;
unsigned long long int time_seed_searchAllSeeds = 0;
unsigned long long int time_seed_mmSeeds = 0;
unsigned long long int time_seed_instantiateSeeds = 0;
unsigned long long int time_extendSeeds_eeSaTups = 0;
unsigned long long int time_extendSeeds_prioritizeSATups = 0;
unsigned long long int time_extendSeeds_advanceElement = 0;
unsigned long long int time_extendSeeds_joinedToTextOff = 0;
unsigned long long int redundants_Nums = 0;
unsigned long long int time_extendSeeds_resEe = 0;
unsigned long long int time_extendSeeds_resUngap = 0;
unsigned long long int time_extendSeeds_init = 0;
unsigned long long int time_extendSeeds_align = 0;
unsigned long long int bestFilter_Nums = 0;
unsigned long long int time_wfa0 = 0;
unsigned long long int time_wfa1 = 0;
unsigned long long int seed_prune = 0;
bool SwAligner::align(
	TAlScore &best,
	SwResult &res,
	const Coord &coord)
{
#ifdef TIME_STATS
	countAlign_Static++;
#endif
#ifdef TIME_STATS
	auto start = std::chrono::system_clock::now();
#endif
	assert(initedRef() && initedRead());
	assert_eq(STATE_INITED, state_);
	state_ = STATE_ALIGNED;
	btncand_.clear();
	btncanddone_.clear();
	btncanddoneSucc_ = btncanddoneFail_ = 0;
	best = std::numeric_limits<TAlScore>::min();
	sse8succ_ = sse16succ_ = false;
	int flag = 0;
	size_t rdlen = rdf_ - rdi_;
	bool checkpointed = rdlen >= cperMinlen_;
	bool gathered = false;
	if (sc_->monotone)
	{
		if (enable8_ && !readSse16_) 
		{
			if (checkpointed)
			{
			}
			else
			{
#ifdef TIME_STATS
				auto start_wfa0 = std::chrono::system_clock::now();
#endif
				char pattern[MAX_ELEN], pat[MAX_ELEN];
				for (int i = rfi_; i < rff_; i++)
				{
					pattern[i] = mask2iupac[(int)rf_[i]];
				}
				int a = 0;
				if ((rff_ - rfi_) > rdlen)
				{
					a = ((rff_ - rfi_) - rdlen) * 0.5;
				}
				strncpy(pat, pattern + a, rdlen);
				parse(pat, bianma[2]);
				int thred = rdlen * CGK_THRED;
				int loop = CGK_LOOP;
				if (!embed_rd)
				{
					char *text = const_cast<char *>(rdfw_->toZBuf());
					parse(text, bianma[0]);
					text = const_cast<char *>(rdrc_->toZBuf());
					parse(text, bianma[1]);
					for (int i = 0; i < loop; i++)
					{
						embedding->cgk_embed((const char **)bianma, rdlen, thred, 0, 0, embeddedQ[i][0], i);
						embedding->cgk_embed((const char **)bianma, rdlen, thred, 1, 0, embeddedQ[i][1], i);
					}
					embed_rd = true;
				}
				int last_hm = 99999, tmp = 99999;
				for (int i = 0; i < loop; i++)
				{
					tmp = embedding->cgk_embed((const char **)bianma, rdlen, thred, 2, 0, fw_ ? embeddedQ[i][0] : embeddedQ[i][1], i);
					if (last_hm > tmp)
						last_hm = tmp;
				}
#ifdef TIME_STATS
				auto end_wfa0 = std::chrono::system_clock::now();
				auto elapsed_wfa0 = std::chrono::duration_cast<std::chrono::microseconds>(end_wfa0 - start_wfa0);
				time_wfa0 += elapsed_wfa0.count();
#endif
				if (last_hm > (rdlen * CGK_THRED))
				{
					return false;
				}
#ifdef TIME_STATS
				auto start_wfa1 = std::chrono::system_clock::now();
#endif
				char *text = const_cast<char *>(rd_->toZBuf());
				affine_wavefronts_t *affine_wavefronts = affine_wavefronts_new_complete(strlen(pattern), strlen(text), &affine_penalties, NULL, mm_allocator);
				affine_wavefronts_align(affine_wavefronts, pattern, strlen(pattern), text, strlen(text));
				best = edit_cigar_score_gap_affine(&affine_wavefronts->edit_cigar, &affine_penalties);
				edit_cigar_t *const edit_cigar = &affine_wavefronts->edit_cigar;
				char *const operations = edit_cigar->operations;
				char cur;
				int num = 0;
				string sb;
				for (int i = edit_cigar->begin_offset; i < edit_cigar->end_offset; ++i)
				{
					cur = operations[i];
					sb.push_back(cur);
				}
				res.alres.mycigar = sb;
				res.alres.setScore(AlnScore(best-1, rdlen, 0, 0, 0));
				res.alres.setShape(refidx_, coord.off(), reflen_, fw_, rdf_ - rdi_, true, 0, 0, true, 0, 0);
				res.alres.setRefNs(0);
				affine_wavefronts_delete(affine_wavefronts);
#ifdef TIME_STATS
				auto end_wfa1 = std::chrono::system_clock::now();
				auto elapsed_wfa1 = std::chrono::duration_cast<std::chrono::microseconds>(end_wfa1 - start_wfa1);
				time_wfa1 += elapsed_wfa1.count();
#endif
#ifndef NDEBUG
#endif
			}
			sse8succ_ = (flag == 0);
#ifndef NDEBUG
			{
			}
#endif
		}
		else
		{
		}
	}
	else
	{
	}
#ifndef NDEBUG
	if (!checkpointed && (rand() & 15) == 0 && sse8succ_ && sse16succ_)
	{
	}
#endif
	assert(repOk());
	cural_ = 0;
	if (best == MIN_I64 || best < minsc_)
	{
		if (dpLog_ != NULL)
		{
			(*dpLog_) << ",0,0";
		}
#ifdef TIME_STATS
		bestFilter_Nums++;
#endif
		return false;
	}
	if (!gathered)
	{
		assert(sse8succ_ || sse16succ_);
		if (sc_->monotone)
		{
			if (sse8succ_)
			{
				btncand_.clear();
				btncand_.expand();
				btncand_.back().init(33, 33, best);
			}
		}
		else
		{
		}
	}
	if (dpLog_ != NULL)
	{
		(*dpLog_) << ",1," << best;
	}
#ifdef TIME_STATS
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	align_time += elapsed.count();
#endif
	return !btncand_.empty();
}
bool SwAligner::nextAlignment(
	SwResult &res,
	TAlScore minsc,
	RandomSource &rnd)
{
	return true;
}
#ifdef MAIN_ALIGNER_SW
#include <sstream>
#include <utility>
#include <getopt.h>
#include "scoring.h"
#include "aligner_seed_policy.h"
int gGapBarrier;
int gSnpPhred;
static int bonusMatchType;
static int bonusMatch;
static int penMmcType;
static int penMmc;
static int penNType;
static int penN;
static bool nPairCat;
static int penRdExConst;
static int penRfExConst;
static int penRdExLinear;
static int penRfExLinear;
static float costMinConst;
static float costMinLinear;
static float costFloorConst;
static float costFloorLinear;
static float nCeilConst;
static float nCeilLinear;
static bool nCatPair;
static int multiseedMms;
static int multiseedLen;
static int multiseedIvalType;
static float multiseedIvalA;
static float multiseedIvalB;
static float posmin;
static float posfrac;
static float rowmult;
enum
{
	ARG_TESTS = 256
};
static const char *short_opts = "s:m:r:d:i:";
static struct option long_opts[] = {
	{(char *)"snppen", required_argument, 0, 's'},
	{(char *)"misspen", required_argument, 0, 'm'},
	{(char *)"seed", required_argument, 0, 'r'},
	{(char *)"align-policy", no_argument, 0, 'A'},
	{(char *)"test", no_argument, 0, ARG_TESTS},
};
template <typename T>
T parse(const char *s)
{
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}
static EList<bool> stbuf, enbuf;
static BTDnaString btread;
static BTString btqual;
static BTString btref;
static BTString btref2;
static BTDnaString readrc;
static BTString qualrc;
static void doTestCase(
	SwAligner &al,
	const BTDnaString &read,
	const BTString &qual,
	const BTString &refin,
	TRefOff off,
	EList<bool> *en,
	const Scoring &sc,
	TAlScore minsc,
	SwResult &res,
	bool nsInclusive,
	bool filterns,
	uint32_t seed)
{
	RandomSource rnd(seed);
	btref2 = refin;
	assert_eq(read.length(), qual.length());
	size_t nrow = read.length();
	TRefOff rfi, rff;
	size_t maxgaps;
	size_t padi, padf;
	{
		int readGaps = sc.maxReadGaps(minsc, read.length());
		int refGaps = sc.maxRefGaps(minsc, read.length());
		assert_geq(readGaps, 0);
		assert_geq(refGaps, 0);
		int maxGaps = max(readGaps, refGaps);
		padi = 2 * maxGaps;
		padf = maxGaps;
		maxgaps = (size_t)maxGaps;
	}
	size_t nceil = (size_t)sc.nCeil.f((double)read.length());
	size_t width = 1 + padi + padf;
	rfi = off;
	off = 0;
	if (rfi < padi)
	{
		size_t beginpad = (size_t)(padi - rfi);
		for (size_t i = 0; i < beginpad; i++)
		{
			btref2.insert('N', 0);
			off--;
		}
		rfi = 0;
	}
	else
	{
		rfi -= padi;
	}
	assert_geq(rfi, 0);
	while (rfi + nrow + padi + padf > btref2.length())
	{
		btref2.append('N');
	}
	rff = rfi + nrow + padi + padf;
	for (size_t i = 0; i < btref2.length(); i++)
	{
		if (toupper(btref2[i]) == 'N' && !nsInclusive)
		{
			btref2.set(16, i);
		}
		else
		{
			int num = 0;
			int alts[] = {4, 4, 4, 4};
			decodeNuc(toupper(btref2[i]), num, alts);
			assert_leq(num, 4);
			assert_gt(num, 0);
			btref2.set(0, i);
			for (int j = 0; j < num; j++)
			{
				btref2.set(btref2[i] | (1 << alts[j]), i);
			}
		}
	}
	bool fw = true;
	uint32_t refidx = 0;
	size_t solwidth = width;
	if (maxgaps >= solwidth)
	{
		solwidth = 0;
	}
	else
	{
		solwidth -= maxgaps;
	}
	if (en == NULL)
	{
		enbuf.resize(solwidth);
		enbuf.fill(true);
		en = &enbuf;
	}
	assert_geq(rfi, 0);
	assert_gt(rff, rfi);
	readrc = read;
	qualrc = qual;
	al.initRead(
		read,
		readrc,
		qual,
		qualrc,
		0,
		read.length(),
		floorsc);
	al.initRef(
		fw,
		refidx,
		off,
		btref2.wbuf(),
		rfi,
		rff,
		width,
		solwidth,
		sc,
		minsc,
		maxgaps,
		0,
		en);
	if (filterns)
	{
		al.filter((int)nceil);
	}
	al.align(rnd);
}
static void doTestCase2(
	SwAligner &al,
	const char *read,
	const char *qual,
	const char *refin,
	TRefOff off,
	const Scoring &sc,
	float costMinConst,
	float costMinLinear,
	SwResult &res,
	bool nsInclusive = false,
	bool filterns = false,
	uint32_t seed = 0)
{
	btread.install(read, true);
	TAlScore minsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costMinConst,
		costMinLinear));
	TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costFloorConst,
		costFloorLinear));
	btqual.install(qual);
	btref.install(refin);
	doTestCase(
		al,
		btread,
		btqual,
		btref,
		off,
		NULL,
		sc,
		minsc,
		floorsc,
		res,
		nsInclusive,
		filterns,
		seed);
}
static void doTestCase3(
	SwAligner &al,
	const char *read,
	const char *qual,
	const char *refin,
	TRefOff off,
	Scoring &sc,
	float costMinConst,
	float costMinLinear,
	float nCeilConst,
	float nCeilLinear,
	SwResult &res,
	bool nsInclusive = false,
	bool filterns = false,
	uint32_t seed = 0)
{
	btread.install(read, true);
	TAlScore minsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costMinConst,
		costMinLinear));
	TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costFloorConst,
		costFloorLinear));
	btqual.install(qual);
	btref.install(refin);
	sc.nCeil.setType(SIMPLE_FUNC_LINEAR);
	sc.nCeil.setConst(costMinConst);
	sc.nCeil.setCoeff(costMinLinear);
	doTestCase(
		al,
		btread,
		btqual,
		btref,
		off,
		NULL,
		sc,
		minsc,
		floorsc,
		res,
		nsInclusive,
		filterns,
		seed);
}
static void doTestCase4(
	SwAligner &al,
	const char *read,
	const char *qual,
	const char *refin,
	TRefOff off,
	EList<bool> &en,
	Scoring &sc,
	float costMinConst,
	float costMinLinear,
	float nCeilConst,
	float nCeilLinear,
	SwResult &res,
	bool nsInclusive = false,
	bool filterns = false,
	uint32_t seed = 0)
{
	btread.install(read, true);
	TAlScore minsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costMinConst,
		costMinLinear));
	TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costFloorConst,
		costFloorLinear));
	btqual.install(qual);
	btref.install(refin);
	sc.nCeil.setType(SIMPLE_FUNC_LINEAR);
	sc.nCeil.setConst(costMinConst);
	sc.nCeil.setCoeff(costMinLinear);
	doTestCase(
		al,
		btread,
		btqual,
		btref,
		off,
		&en,
		sc,
		minsc,
		floorsc,
		res,
		nsInclusive,
		filterns,
		seed);
}
static void doTests()
{
	bonusMatchType = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch = DEFAULT_MATCH_BONUS;
	penMmcType = DEFAULT_MM_PENALTY_TYPE;
	penMmc = DEFAULT_MM_PENALTY;
	penSnp = DEFAULT_SNP_PENALTY;
	penNType = DEFAULT_N_PENALTY_TYPE;
	penN = DEFAULT_N_PENALTY;
	nPairCat = DEFAULT_N_CAT_PAIR;
	penRdExConst = DEFAULT_READ_GAP_CONST;
	penRfExConst = DEFAULT_REF_GAP_CONST;
	penRdExLinear = DEFAULT_READ_GAP_LINEAR;
	penRfExLinear = DEFAULT_REF_GAP_LINEAR;
	costMinConst = DEFAULT_MIN_CONST;
	costMinLinear = DEFAULT_MIN_LINEAR;
	costFloorConst = DEFAULT_FLOOR_CONST;
	costFloorLinear = DEFAULT_FLOOR_LINEAR;
	nCeilConst = 1.0f;
	nCeilLinear = 0.1f;
	multiseedMms = DEFAULT_SEEDMMS;
	multiseedLen = DEFAULT_SEEDLEN;
	Scoring sc(
		bonusMatch,
		penMmcType,
		30,
		30,
		costMinConst,
		costMinLinear,
		costFloorConst,
		costFloorLinear,
		nCeilConst,
		nCeilLinear,
		penNType,
		penN,
		nPairCat,
		25,
		25,
		15,
		15,
		1,
		-1,
		false);
	Scoring sc2(
		bonusMatch,
		COST_MODEL_QUAL,
		30,
		30,
		costMinConst,
		costMinLinear,
		costFloorConst,
		costFloorLinear,
		1.0f,
		1.0f,
		penNType,
		penN,
		nPairCat,
		25,
		25,
		15,
		15,
		1,
		-1,
		false);
	SwResult res;
	cerr << "Running tests..." << endl;
	int tests = 1;
	bool nIncl = false;
	bool nfilter = false;
	SwAligner al;
	RandomSource rnd(73);
	for (int i = 0; i < 3; i++)
	{
		cerr << "  Test " << tests++ << " (nuc space, offset "
			 << (i * 4) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		sc.rdGapLinear = 15;
		sc.rfGapLinear = 15;
		doTestCase2(
			al,
			"ACGTACGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1mm allowed by minsc)...";
		sc.setMmPen(COST_MODEL_CONSTANT, 30);
		doTestCase2(
			al,
			"ACGTTCGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1mm allowed by minsc, check qual 1)...";
		doTestCase2(
			al,
			"ACGTTCGT",
			"ABCDEFGH",
			"ACGTACGTACGTACGT",
			i * 4,
			sc2,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		size_t lo, hi;
		if (i == 0)
		{
			lo = 0;
			hi = 1;
		}
		else if (i == 1)
		{
			lo = 1;
			hi = 2;
		}
		else
		{
			lo = 2;
			hi = 3;
		}
		for (size_t j = lo; j < hi; j++)
		{
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(j * 4, res.alres.refoff());
			assert_eq(8, res.alres.refExtent());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), -36);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
		}
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1mm allowed by minsc, check qual 2)...";
		doTestCase2(
			al,
			"ACGAACGT",
			"ABCDEFGH",
			"ACGTACGTACGTACGT",
			i * 4,
			sc2,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -35);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1mm allowed by minsc, check qual )...";
		assert(res.empty());
		doTestCase2(
			al,
			"TCGTACGT",
			"ABCDEFGH",
			"ACGTACGTACGTACGT",
			i * 4,
			sc2,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -32);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1mm at the beginning, allowed by minsc)...";
		doTestCase2(
			al,
			"CCGTACGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(res.alres.score().ns(), 0);
		assert_eq(1, res.alres.ned().size());
		assert_eq(0, res.alres.aed().size());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1 n in read, allowed)...";
		doTestCase3(
			al,
			"ACGTNCGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			1.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(res.alres.score().ns(), 1);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 2 n in read, allowed)...";
		doTestCase3(
			al,
			"ACGNNCGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			2.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 2 n in read, 1 at beginning, allowed)...";
		doTestCase2(
			al,
			"NCGTNCGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1 n in ref, allowed)...";
		doTestCase2(
			al,
			"ACGTACGT",
			"IIIIIIII",
			"ACGTNCGTACGTANGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(res.alres.score().ns(), 1);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1mm disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTTCGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-10.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", read gap allowed by minsc)...";
		assert(res.empty());
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", read gap disallowed by minsc)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", ref gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", read gap disallowed by gap barrier)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		sc.gapbar = 4;
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		sc.gapbar = 1;
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", ref gap allowed by minsc, gapbar=3)...";
		sc.gapbar = 3;
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		sc.gapbar = 1;
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", ref gap allowed by minsc, gapbar=4)...";
		sc.gapbar = 4;
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		sc.gapbar = 1;
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", ref gap disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(al.done());
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", ref gap disallowed by gap barrier)...";
		sc.gapbar = 5;
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		sc.gapbar = 1;
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(al.done());
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1 read gap, ref gaps disallowed by minsc)...";
		sc.rfGapConst = 35;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 20;
		sc.rdGapLinear = 10;
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -35);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", gaps disallowed by minsc)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 10;
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		sc.rfGapConst = 25;
		sc.rdGapConst = 35;
		sc.rfGapLinear = 12;
		sc.rdGapLinear = 22;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1 ref gap, read gaps disallowed by minsc)...";
		assert(res.empty());
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -37);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;
		sc.rfGapConst = 20;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 15;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1 read gap, 2 ref gaps allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;
		sc.rfGapConst = 25;
		sc.rdGapConst = 11;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 10;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1 ref gap, 2 read gaps allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4) << ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(al.done());
		cerr << "PASSED" << endl;
		sc.rfGapConst = 15;
		sc.rdGapConst = 15;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 10;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 2 ref gaps, 2 read gaps allowed by minsc)...";
		doTestCase3(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-40.0f,
			0.0f,
			1.0,
			0.0,
			res,
			nIncl,
			true);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		if (!res.empty())
		{
		}
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -25);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(i == 2 || res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		sc.rfGapConst = 10;
		sc.rdGapConst = 10;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 10;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1 ref gap, 1 read gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		sc.rfGapConst = 15;
		sc.rdGapConst = 15;
		sc.rfGapLinear = 5;
		sc.rdGapLinear = 5;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 2 ref gaps, 2 read gaps allowed by minsc)...";
		doTestCase3(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-35.0f,
			0.0f,
			1.0f,
			0.0f,
			res,
			nIncl,
			true);
		if (i == 0)
		{
			lo = 0;
			hi = 1;
		}
		else if (i == 1)
		{
			lo = 1;
			hi = 2;
		}
		else
		{
			lo = 2;
			hi = 3;
		}
		for (size_t j = lo; j < hi; j++)
		{
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert(res.alres.refoff() == 0 ||
				   res.alres.refoff() == 4 ||
				   res.alres.refoff() == 8);
			assert_eq(8, res.alres.refExtent());
			assert_eq(res.alres.score().gaps(), 1);
			assert_eq(res.alres.score().score(), -20);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
		}
		al.nextAlignment(res, rnd);
		cerr << "PASSED" << endl;
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 4;
		sc.rdGapLinear = 4;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1 ref gap, 1 read gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -29);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", short read)...";
		doTestCase2(
			al,
			"A",
			"I",
			"AAAAAAAAAAAA",
			i * 4,
			sc,
			-30.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		if (i == 0)
		{
			cerr << "  Test " << tests++
				 << " (nuc space, offset 0, short read & ref)...";
			doTestCase2(
				al,
				"A",
				"I",
				"A",
				0,
				sc,
				-30.0f,
				0.0f,
				res,
				nIncl,
				nfilter);
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
			cerr << "PASSED" << endl;
		}
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", short read, many allowed gaps)...";
		doTestCase2(
			al,
			"A",
			"I",
			"AAAAAAAAAAAA",
			i * 4,
			sc,
			-150.0f,
			0.0f,
			res,
			nIncl,
			nfilter);
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		if (i == 0)
		{
			cerr << "  Test " << tests++
				 << " (nuc space, offset 0, short read & ref, "
				 << "many allowed gaps)...";
			doTestCase2(
				al,
				"A",
				"I",
				"A",
				0,
				sc,
				-150.0f,
				0.0f,
				res,
				nIncl,
				nfilter);
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
			cerr << "PASSED" << endl;
		}
	}
	cerr << "  Test " << tests++ << " (N ceiling 1)...";
	sc.mmcostType = COST_MODEL_CONSTANT;
	sc.mmcost = 10;
	sc.snp = 30;
	sc.nCeilConst = 0.0f;
	sc.nCeilLinear = 0.0f;
	sc.rfGapConst = 10;
	sc.rdGapLinear = 10;
	sc.rfGapConst = 10;
	sc.rfGapLinear = 10;
	sc.setNPen(COST_MODEL_CONSTANT, 2);
	sc.gapbar = 1;
	doTestCase2(
		al,
		"ACGTACGT",
		"IIIIIIII",
		"NCGTACGT",
		0,
		sc,
		-25.0f,
		0.0f,
		res,
		false,
		true,
		0);
	al.nextAlignment(res, rnd);
	assert(res.empty());
	cerr << "PASSED" << endl;
	res.reset();
	cerr << "  Test " << tests++ << " (N ceiling 2)...";
	doTestCase3(
		al,
		"ACGTACGT",
		"IIIIIIII",
		"NCGTACGT",
		0,
		sc,
		-25.0f,
		0.0f,
		1.0f,
		0.0f,
		res,
		false,
		false,
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(0, res.alres.score().gaps());
	assert_eq(-2, res.alres.score().score());
	assert_eq(1, res.alres.score().ns());
	cerr << "PASSED" << endl;
	res.reset();
	for (size_t i = 0; i < 2; i++)
	{
		cerr << "  Test " << tests++ << " (N ceiling 2 with st_ override)...";
		EList<bool> en;
		en.resize(3);
		en.fill(true);
		if (i == 1)
		{
			en[1] = false;
		}
		sc.rfGapConst = 10;
		sc.rdGapLinear = 10;
		sc.rfGapConst = 10;
		sc.rfGapLinear = 10;
		doTestCase4(
			al,
			"ACGTACGT",
			"IIIIIIII",
			"NCGTACGT",
			0,
			en,
			sc,
			-25.0f,
			0.0f,
			1.0f,
			0.0f,
			res,
			false,
			false,
			0);
		al.nextAlignment(res, rnd);
		if (i > 0)
		{
			assert(res.empty());
		}
		else
		{
			assert(!res.empty());
		}
		cerr << "PASSED" << endl;
		res.reset();
	}
	cerr << "  Test " << tests++ << " (N ceiling 3)...";
	sc.nCeilConst = 1.0f;
	sc.nCeilLinear = 0.0f;
	doTestCase2(
		al,
		"ACGTACGT",
		"IIIIIIII",
		"NCGTACGT",
		0,
		sc,
		-25.0f,
		0.0f,
		res,
		false,
		true,
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(0, res.alres.score().gaps());
	assert_eq(-2, res.alres.score().score());
	assert_eq(1, res.alres.score().ns());
	cerr << "PASSED" << endl;
	res.reset();
	cerr << "  Test " << tests++ << " (redundant alignment elimination 1)...";
	sc.nCeilConst = 1.0f;
	sc.nCeilLinear = 0.0f;
	sc.rfGapConst = 25;
	sc.rdGapLinear = 15;
	sc.rfGapConst = 25;
	sc.rfGapLinear = 15;
	doTestCase2(
		al,
		"AGGCTATGCCTCTGACGCGATATCGGCGCCCACTTCAGAGCTAACCG",
		"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
		"TTTTTTTTAGGCTATGCCTCTGACGCGATATCGGCGCCCACTTCAGAGCTAACCGTTTTTTT",
		8,
		sc,
		-25.0f,
		-5.0f,
		res,
		false,
		true,
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(8, res.alres.refoff());
	assert_eq(47, res.alres.refExtent());
	assert_eq(0, res.alres.score().gaps());
	assert_eq(0, res.alres.score().score());
	assert_eq(0, res.alres.score().ns());
	res.reset();
	al.nextAlignment(res, rnd);
	assert(res.empty());
	assert(al.done());
	cerr << "PASSED" << endl;
	res.reset();
}
static void doLocalTests()
{
	bonusMatchType = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch = DEFAULT_MATCH_BONUS_LOCAL;
	penMmcType = DEFAULT_MM_PENALTY_TYPE;
	penMmc = DEFAULT_MM_PENALTY;
	penSnp = DEFAULT_SNP_PENALTY;
	penNType = DEFAULT_N_PENALTY_TYPE;
	penN = DEFAULT_N_PENALTY;
	nPairCat = DEFAULT_N_CAT_PAIR;
	penRdExConst = DEFAULT_READ_GAP_CONST;
	penRfExConst = DEFAULT_REF_GAP_CONST;
	penRdExLinear = DEFAULT_READ_GAP_LINEAR;
	penRfExLinear = DEFAULT_REF_GAP_LINEAR;
	costMinConst = DEFAULT_MIN_CONST_LOCAL;
	costMinLinear = DEFAULT_MIN_LINEAR_LOCAL;
	costFloorConst = DEFAULT_FLOOR_CONST_LOCAL;
	costFloorLinear = DEFAULT_FLOOR_LINEAR_LOCAL;
	nCeilConst = 1.0f;
	nCeilLinear = 0.1f;
	multiseedMms = DEFAULT_SEEDMMS;
	multiseedLen = DEFAULT_SEEDLEN;
	Scoring sc(
		10,
		penMmcType,
		30,
		penSnp,
		costMinConst,
		costMinLinear,
		costFloorConst,
		costFloorLinear,
		nCeilConst,
		nCeilLinear,
		penNType,
		penN,
		nPairCat,
		25,
		25,
		15,
		15,
		1,
		-1,
		false);
	SwResult res;
	cerr << "Running local tests..." << endl;
	int tests = 1;
	bool nIncl = false;
	bool nfilter = false;
	SwAligner al;
	RandomSource rnd(73);
	for (int i = 0; i < 3; i++)
	{
		cerr << "  Test " << tests++ << " (short nuc space, offset "
			 << (i * 4) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGT",
			"IIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			0.0f,
			8.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(4, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 40);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (short nuc space, offset "
			 << (i * 4) << ", 1mm)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"CCGT",
			"IIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			0.0f,
			7.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4 + 1, res.alres.refoff());
		assert_eq(3, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 30);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (short nuc space, offset "
			 << (i * 4) << ", 1mm)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGA",
			"IIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			0.0f,
			7.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(3, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 30);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		if (i == 0)
		{
			cerr << "  Test " << tests++ << " (short nuc space, offset "
				 << (i * 4) << ", 1mm, match bonus=20)...";
			sc.rdGapConst = 40;
			sc.rfGapConst = 40;
			sc.setMatchBonus(20);
			doTestCase2(
				al,
				"TTGT",
				"IIII",
				"TTGA",
				i * 4,
				sc,
				25.0f,
				0.0f,
				res,
				nIncl,
				nfilter);
			assert(!al.done());
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(i * 4, res.alres.refoff());
			assert_eq(3, res.alres.refExtent());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 60);
			assert_eq(res.alres.score().ns(), 0);
			assert(res.alres.ned().empty());
			assert(res.alres.aed().empty());
			res.reset();
			al.nextAlignment(res, rnd);
			assert(res.empty());
			assert(al.done());
			res.reset();
			sc.setMatchBonus(10);
			cerr << "PASSED" << endl;
		}
		cerr << "  Test " << tests++ << " (nuc space, offset "
			 << (i * 4) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGTACGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			0.0f,
			8.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 80);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (long nuc space, offset "
			 << (i * 8) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGTACGTACGTACGTACGTA",
			"IIIIIIIIIIIIIIIIIIIII",
			"ACGTACGTACGTACGTACGTACGTACGTACGTACGTA",
			i * 8,
			sc,
			0.0f,
			8.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 8, res.alres.refoff());
		assert_eq(21, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 210);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		res.reset();
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i * 4)
			 << ", 1mm allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTTCGT",
			"IIIIIIII",
			"ACGTACGTACGTACGT",
			i * 4,
			sc,
			0.0f,
			5.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		cerr << "PASSED" << endl;
		cerr << "  Test " << tests++ << " (long nuc space, offset "
			 << (i * 8) << ", 6mm allowed by minsc)...";
		sc.rdGapConst = 50;
		sc.rfGapConst = 50;
		sc.rdGapLinear = 45;
		sc.rfGapLinear = 45;
		doTestCase2(
			al,
			"ACGTACGATGCATCGTACGTA",
			"IIIIIIIIIIIIIIIIIIIII",
			"ACGTACGTACGTACGTACGTACGTACGTACGTACGTA",
			i * 8,
			sc,
			0.0f,
			1.0f,
			res,
			nIncl,
			nfilter);
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i * 8 + 13, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 80);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		res.reset();
		cerr << "PASSED" << endl;
	}
}
int main(int argc, char **argv)
{
	int option_index = 0;
	int next_option;
	unsigned seed = 0;
	gGapBarrier = 1;
	gSnpPhred = 30;
	bool nsInclusive = false;
	bool nfilter = false;
	bonusMatchType = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch = DEFAULT_MATCH_BONUS;
	penMmcType = DEFAULT_MM_PENALTY_TYPE;
	penMmc = DEFAULT_MM_PENALTY;
	penSnp = DEFAULT_SNP_PENALTY;
	penNType = DEFAULT_N_PENALTY_TYPE;
	penN = DEFAULT_N_PENALTY;
	penRdExConst = DEFAULT_READ_GAP_CONST;
	penRfExConst = DEFAULT_REF_GAP_CONST;
	penRdExLinear = DEFAULT_READ_GAP_LINEAR;
	penRfExLinear = DEFAULT_REF_GAP_LINEAR;
	costMinConst = DEFAULT_MIN_CONST;
	costMinLinear = DEFAULT_MIN_LINEAR;
	costFloorConst = DEFAULT_FLOOR_CONST;
	costFloorLinear = DEFAULT_FLOOR_LINEAR;
	nCeilConst = 1.0f;
	nCeilLinear = 1.0f;
	nCatPair = false;
	multiseedMms = DEFAULT_SEEDMMS;
	multiseedLen = DEFAULT_SEEDLEN;
	multiseedIvalType = DEFAULT_IVAL;
	multiseedIvalA = DEFAULT_IVAL_A;
	multiseedIvalB = DEFAULT_IVAL_B;
	mhits = 1;
	do
	{
		next_option = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		switch (next_option)
		{
		case 's':
			gSnpPhred = parse<int>(optarg);
			break;
		case 'r':
			seed = parse<unsigned>(optarg);
			break;
		case ARG_TESTS:
		{
			doTests();
			cout << "PASSED end-to-ends" << endl;
			doLocalTests();
			cout << "PASSED locals" << endl;
			return 0;
		}
		case 'A':
		{
			bool localAlign = false;
			bool noisyHpolymer = false;
			bool ignoreQuals = false;
			SeedAlignmentPolicy::parseString(
				optarg,
				localAlign,
				noisyHpolymer,
				ignoreQuals,
				bonusMatchType,
				bonusMatch,
				penMmcType,
				penMmc,
				penNType,
				penN,
				penRdExConst,
				penRfExConst,
				penRdExLinear,
				penRfExLinear,
				costMinConst,
				costMinLinear,
				costFloorConst,
				costFloorLinear,
				nCeilConst,
				nCeilLinear,
				nCatPair,
				multiseedMms,
				multiseedLen,
				multiseedIvalType,
				multiseedIvalA,
				multiseedIvalB,
				posmin);
			break;
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
	srand(seed);
	if (argc - optind < 4)
	{
		cerr << "Need at least 4 arguments" << endl;
		exit(1);
	}
	BTDnaString read;
	BTString ref, qual;
	read.installChars(argv[optind]);
	qual.install(argv[optind + 1]);
	assert_eq(read.length(), qual.length());
	ref.install(argv[optind + 2]);
	size_t off = parse<size_t>(argv[optind + 3]);
	Scoring sc(
		false,
		false,
		bonusMatchType,
		bonusMatch,
		penMmcType,
		penMmc,
		costMinConst,
		costMinLinear,
		costFloorConst,
		costFloorLinear,
		nCeilConst,
		nCeilLinear,
		penNType,
		penN,
		nCatPair,
		penRdExConst,
		penRfExConst,
		penRdExLinear,
		penRfExLinear);
	TAlScore minsc = Scoring::linearFunc(
		read.length(),
		costMinConst,
		costMinLinear);
	TAlScore floorsc = Scoring::linearFunc(
		read.length(),
		costFloorConst,
		costFloorLinear);
	SwResult res;
	SwAligner al;
	doTestCase(
		al,
		read,
		qual,
		ref,
		off,
		NULL,
		sc,
		minsc,
		res,
		nsInclusive,
		nfilter,
		seed);
}
#endif
