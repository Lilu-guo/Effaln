#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <stdexcept>
#include <getopt.h>
#include <math.h>
#include <utility>
#include <limits>
#include <time.h>
#include <dirent.h>
#ifndef _WIN32
#include <signal.h>
#endif
#include "alphabet.h"
#include "assert_helpers.h"
#include "endian_swap.h"
#include "aln_idx.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "aln_sink.h"
#include "pat.h"
#include "threading.h"
#include "ds.h"
#include "aligner_metrics.h"
#include "sam.h"
#include "aligner_seed.h"
#include "aligner_seed_policy.h"
#include "aligner_driver.h"
#include "aligner_sw.h"
#include "aligner_sw_driver.h"
#include "aligner_cache.h"
#include "util.h"
#include "pe.h"
#include "simple_func.h"
#include "presets.h"
#include "opts.h"
#include "outq.h"
#include "aligner_seed2.h"
#include "aln_search.h"
#ifdef WITH_TBB
#include <thread>
#endif
#if __cplusplus <= 199711L
#define unique_ptr auto_ptr
#endif
using namespace std;
static int FNAME_SIZE;
#ifdef WITH_TBB
static std::atomic<int> thread_counter;
#else
static int thread_counter;
static MUTEX_T thread_counter_mutex;
#endif
static EList<string> mates1;
static EList<string> mates2;
static EList<string> mates12;
static string adjIdxBase;
int gVerbose;
static bool startVerbose;
int gQuiet;
static int sanityCheck;
static int format;
static bool interleaved;
static string origString;
static int seed;
static int timing;
static int metricsIval;
static string metricsFile;
static bool metricsStderr;
static bool metricsPerRead;
static bool allHits;
static bool showVersion;
static int ipause;
static uint64_t qUpto;
static int gTrim5;
static int gTrim3;
static pair<short, size_t> trimTo;
static int offRate;
static bool solexaQuals;
static bool phred64Quals;
static bool integerQuals;
static int nthreads;
static int thread_ceiling;
static string thread_stealing_dir;
static bool thread_stealing;
static int outType;
static bool noRefNames;
static uint32_t khits;
static uint32_t mhits;
static int partitionSz;
static int readsPerBatch;
static bool fileParallel;
static bool useShmem;
static bool useMm;
static bool mmSweep;
int gMinInsert;
int gMaxInsert;
bool gMate1fw;
bool gMate2fw;
bool gFlippedMatesOK;
bool gDovetailMatesOK;
bool gContainMatesOK;
bool gOlapMatesOK;
bool gExpandToFrag;
bool gReportDiscordant;
bool gReportMixed;
static uint32_t cacheLimit;
static uint32_t cacheSize;
static uint32_t skipReads;
bool gNofw;
bool gNorc;
static uint32_t fastaContLen;
static uint32_t fastaContFreq;
static bool hadoopOut;
static bool fullRef;
static bool samTruncQname;
static bool samAppendComment;
static bool samOmitSecSeqQual;
static bool samNoUnal;
static bool samNoHead;
static bool samNoSQ;
static bool sam_print_as;
static bool sam_print_xs;
static bool sam_print_xss;
static bool sam_print_yn;
static bool sam_print_xn;
static bool sam_print_x0;
static bool sam_print_x1;
static bool sam_print_xm;
static bool sam_print_xo;
static bool sam_print_xg;
static bool sam_print_nm;
static bool sam_print_md;
static bool sam_print_yf;
static bool sam_print_yi;
static bool sam_print_ym;
static bool sam_print_yp;
static bool sam_print_yt;
static bool sam_print_ys;
static bool sam_print_zs;
static bool sam_print_xr;
static bool sam_print_xt;
static bool sam_print_xd;
static bool sam_print_xu;
static bool sam_print_yl;
static bool sam_print_ye;
static bool sam_print_yu;
static bool sam_print_xp;
static bool sam_print_yr;
static bool sam_print_zb;
static bool sam_print_zr;
static bool sam_print_zf;
static bool sam_print_zm;
static bool sam_print_zi;
static bool sam_print_zp;
static bool sam_print_zu;
static bool sam_print_zt;
static bool preserve_tags;
static bool align_paired_reads;
static bool bwaSwLike;
static bool gSeedLenIsSet;
static float bwaSwLikeC;
static float bwaSwLikeT;
static bool qcFilter;
bool gReportOverhangs;
static string rgid;
static string rgs;
static string rgs_optflag;
static bool msample;
int gGapBarrier;
int gDefaultSeedLen;
static EList<string> qualities;
static EList<string> qualities1;
static EList<string> qualities2;
static string polstr;
static bool msNoCache;
static int bonusMatchType;
static int bonusMatch;
static int penMmcType;
static int penMmcMax;
static int penMmcMin;
static int penNType;
static int penN;
static bool penNCatPair;
static bool localAlign;
static bool noisyHpolymer;
static int penRdGapConst;
static int penRfGapConst;
static int penRdGapLinear;
static int penRfGapLinear;
static SimpleFunc scoreMin;
static SimpleFunc nCeil;
static SimpleFunc msIval;
static double descConsExp;
static bool descPrioritizeRoots;
static size_t descLanding;
static SimpleFunc descentTotSz;
static SimpleFunc descentTotFmops;
static int multiseedMms;
static int multiseedLen;
static size_t multiseedOff;
static uint32_t seedCacheLocalMB;
static uint32_t seedCacheCurrentMB;
static uint32_t exactCacheCurrentMB;
static size_t maxhalf;
static bool seedSumm;
static bool scUnMapped;
static bool doUngapped;
static bool xeq;
static size_t maxIters;
static size_t maxUg;
static size_t maxDp;
static size_t maxItersIncr;
static size_t maxEeStreak;
static size_t maxUgStreak;
static size_t maxDpStreak;
static size_t maxStreakIncr;
static size_t maxMateStreak;
static bool doExtend;
static bool enable8;
static size_t cminlen;
static size_t cpow2;
static bool doTri;
static string defaultPreset;
static bool ignoreQuals;
static string wrapper;
static EList<string> queries;
static string outfile;
static int mapqv;
static int tighten;
static bool doExactUpFront;
static bool do1mmUpFront;
static size_t do1mmMinLen;
static int seedBoostThresh;
static size_t nSeedRounds;
static bool reorder;
static float sampleFrac;
static bool arbitraryRandom;
static bool nextVer;
static string logDps;
static string logDpsOpp;
static string bt2index;
static EList<pair<int, string>> extra_opts;
static size_t extra_opts_cur;
uint8_t code[256];
char rcsymbol[6] = "TGCAN";
#ifdef USE_SRA
static EList<string> sra_accs;
#endif
#define DMAX std::numeric_limits<double>::max()
static void resetOptions()
{
	mates1.clear();
	mates2.clear();
	mates12.clear();
	adjIdxBase = "";
	gVerbose = 0;
	startVerbose = 0;
	gQuiet = false;
	sanityCheck = 0;
	format = FASTQ;
	interleaved = false;
	origString = "";
	seed = 0;
	timing = 0;
	metricsIval = 1;
	metricsFile = "";
	metricsStderr = false;
	metricsPerRead = false;
	allHits = false;
	showVersion = false;
	ipause = 0;
	qUpto = 0xffffffffffffffff;
	gTrim5 = 0;
	gTrim3 = 0;
	trimTo = pair<short, size_t>(5, 0);
	offRate = -1;
	solexaQuals = false;
	phred64Quals = false;
	integerQuals = false;
	nthreads = 1;
	thread_ceiling = 0;
	thread_stealing_dir = "";
	thread_stealing = false;
	FNAME_SIZE = 4096;
	outType = OUTPUT_SAM;
	noRefNames = false;
	khits = 1;
	mhits = 50;
	partitionSz = 0;
	readsPerBatch = 16;
	fileParallel = false;
	useShmem = false;
	useMm = false;
	mmSweep = false;
	gMinInsert = 0;
	gMaxInsert = 500;
	gMate1fw = true;
	gMate2fw = false;
	gFlippedMatesOK = false;
	gDovetailMatesOK = false;
	gContainMatesOK = true;
	gOlapMatesOK = true;
	gExpandToFrag = true;
	gReportDiscordant = true;
	gReportMixed = true;
	cacheLimit = 5;
	cacheSize = 0;
	skipReads = 0;
	gNofw = false;
	gNorc = false;
	fastaContLen = 0;
	fastaContFreq = 0;
	hadoopOut = false;
	fullRef = false;
	samTruncQname = true;
	samAppendComment = false;
	samOmitSecSeqQual = false;
	samNoUnal = false;
	samNoHead = false;
	samNoSQ = false;
	sam_print_as = true;
	sam_print_xs = true;
	sam_print_xss = false;
	sam_print_yn = false;
	sam_print_xn = true;
	sam_print_x0 = true;
	sam_print_x1 = true;
	sam_print_xm = true;
	sam_print_xo = true;
	sam_print_xg = true;
	sam_print_nm = true;
	sam_print_md = true;
	sam_print_yf = true;
	sam_print_yi = false;
	sam_print_ym = false;
	sam_print_yp = false;
	sam_print_yt = true;
	sam_print_ys = true;
	sam_print_zs = false;
	sam_print_xr = false;
	sam_print_xt = false;
	sam_print_xd = false;
	sam_print_xu = false;
	sam_print_yl = false;
	sam_print_ye = false;
	sam_print_yu = false;
	sam_print_xp = false;
	sam_print_yr = false;
	sam_print_zb = false;
	sam_print_zr = false;
	sam_print_zf = false;
	sam_print_zm = false;
	sam_print_zi = false;
	sam_print_zp = false;
	sam_print_zu = false;
	sam_print_zt = false;
	preserve_tags = false;
	align_paired_reads = false;
	bwaSwLike = false;
	gSeedLenIsSet = false;
	bwaSwLikeC = 5.5f;
	bwaSwLikeT = 20.0f;
	gDefaultSeedLen = DEFAULT_SEEDLEN;
	qcFilter = false;
	rgid = "";
	rgs = "";
	rgs_optflag = "";
	msample = true;
	gGapBarrier = 4;
	qualities.clear();
	qualities1.clear();
	qualities2.clear();
	polstr.clear();
	msNoCache = true;
	bonusMatchType = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch = DEFAULT_MATCH_BONUS;
	penMmcType = DEFAULT_MM_PENALTY_TYPE;
	penMmcMax = DEFAULT_MM_PENALTY_MAX;
	penMmcMin = DEFAULT_MM_PENALTY_MIN;
	penNType = DEFAULT_N_PENALTY_TYPE;
	penN = DEFAULT_N_PENALTY;
	penNCatPair = DEFAULT_N_CAT_PAIR;
	localAlign = false;
	noisyHpolymer = false;
	penRdGapConst = DEFAULT_READ_GAP_CONST;
	penRfGapConst = DEFAULT_REF_GAP_CONST;
	penRdGapLinear = DEFAULT_READ_GAP_LINEAR;
	penRfGapLinear = DEFAULT_REF_GAP_LINEAR;
	scoreMin.init(SIMPLE_FUNC_LINEAR, DEFAULT_MIN_CONST, DEFAULT_MIN_LINEAR);
	nCeil.init(SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 2.0f, 0.1f);
	msIval.init(SIMPLE_FUNC_LINEAR, 1.0f, DMAX, DEFAULT_IVAL_B, DEFAULT_IVAL_A);
	descConsExp = 2.0;
	descPrioritizeRoots = false;
	descLanding = 20;
	descentTotSz.init(SIMPLE_FUNC_LINEAR, 1024.0, DMAX, 0.0, 1024.0);
	descentTotFmops.init(SIMPLE_FUNC_LINEAR, 100.0, DMAX, 0.0, 10.0);
	multiseedMms = DEFAULT_SEEDMMS;
	multiseedLen = gDefaultSeedLen;
	multiseedOff = 0;
	seedCacheLocalMB = 32;
	seedCacheCurrentMB = 20;
	exactCacheCurrentMB = 20;
	maxhalf = 15;
	seedSumm = false;
	scUnMapped = false;
	xeq = false;
	doUngapped = true;
	maxIters = 400;
	maxUg = 300;
	maxDp = 300;
	maxItersIncr = 20;
	maxEeStreak = 15;
	maxUgStreak = 15;
	maxDpStreak = 15;
	maxStreakIncr = 10;
	maxMateStreak = 10;
	doExtend = true;
	enable8 = true;
	cminlen = 2000;
	cpow2 = 4;
	doTri = false;
	defaultPreset = "sensitive%LOCAL%";
	extra_opts.clear();
	extra_opts_cur = 0;
	bt2index.clear();
	ignoreQuals = false;
	wrapper.clear();
	queries.clear();
	outfile.clear();
	mapqv = 2;
	tighten = 3;
	doExactUpFront = true;
	do1mmUpFront = true;
	seedBoostThresh = 300;
	nSeedRounds = 2;
	do1mmMinLen = 60;
	reorder = false;
	sampleFrac = 1.1f;
	arbitraryRandom = false;
	nextVer = false;
	logDps.clear();
	logDpsOpp.clear();
#ifdef USE_SRA
	sra_accs.clear();
#endif
}
static const char *short_options = "bfF:qbzhcu:rv:s:aP:t3:5:w:p:k:M:1:2:I:X:CQ:N:i:L:U:x:S:g:O:D:R:";
static struct option long_options[] = {
	{(char *)"verbose", no_argument, 0, ARG_VERBOSE},
	{(char *)"startverbose", no_argument, 0, ARG_STARTVERBOSE},
	{(char *)"quiet", no_argument, 0, ARG_QUIET},
	{(char *)"sanity", no_argument, 0, ARG_SANITY},
	{(char *)"pause", no_argument, &ipause, 1},
	{(char *)"orig", required_argument, 0, ARG_ORIG},
	{(char *)"all", no_argument, 0, 'a'},
	{(char *)"solexa-quals", no_argument, 0, ARG_SOLEXA_QUALS},
	{(char *)"integer-quals", no_argument, 0, ARG_INTEGER_QUALS},
	{(char *)"int-quals", no_argument, 0, ARG_INTEGER_QUALS},
	{(char *)"metrics", required_argument, 0, ARG_METRIC_IVAL},
	{(char *)"metrics-file", required_argument, 0, ARG_METRIC_FILE},
	{(char *)"metrics-stderr", no_argument, 0, ARG_METRIC_STDERR},
	{(char *)"metrics-per-read", no_argument, 0, ARG_METRIC_PER_READ},
	{(char *)"met-read", no_argument, 0, ARG_METRIC_PER_READ},
	{(char *)"met", required_argument, 0, ARG_METRIC_IVAL},
	{(char *)"met-file", required_argument, 0, ARG_METRIC_FILE},
	{(char *)"met-stderr", no_argument, 0, ARG_METRIC_STDERR},
	{(char *)"time", no_argument, 0, 't'},
	{(char *)"trim3", required_argument, 0, '3'},
	{(char *)"trim5", required_argument, 0, '5'},
	{(char *)"seed", required_argument, 0, ARG_SEED},
	{(char *)"qupto", required_argument, 0, 'u'},
	{(char *)"upto", required_argument, 0, 'u'},
	{(char *)"version", no_argument, 0, ARG_VERSION},
	{(char *)"reads-per-batch", required_argument, 0, ARG_READS_PER_BATCH},
	{(char *)"filepar", no_argument, 0, ARG_FILEPAR},
	{(char *)"help", no_argument, 0, 'h'},
	{(char *)"threads", required_argument, 0, 'p'},
	{(char *)"khits", required_argument, 0, 'k'},
	{(char *)"minins", required_argument, 0, 'I'},
	{(char *)"maxins", required_argument, 0, 'X'},
	{(char *)"quals", required_argument, 0, 'Q'},
	{(char *)"Q1", required_argument, 0, ARG_QUALS1},
	{(char *)"Q2", required_argument, 0, ARG_QUALS2},
	{(char *)"refidx", no_argument, 0, ARG_REFIDX},
	{(char *)"partition", required_argument, 0, ARG_PARTITION},
	{(char *)"ff", no_argument, 0, ARG_FF},
	{(char *)"fr", no_argument, 0, ARG_FR},
	{(char *)"rf", no_argument, 0, ARG_RF},
	{(char *)"cachelim", required_argument, 0, ARG_CACHE_LIM},
	{(char *)"cachesz", required_argument, 0, ARG_CACHE_SZ},
	{(char *)"nofw", no_argument, 0, ARG_NO_FW},
	{(char *)"norc", no_argument, 0, ARG_NO_RC},
	{(char *)"skip", required_argument, 0, 's'},
	{(char *)"12", required_argument, 0, ARG_ONETWO},
	{(char *)"tab5", required_argument, 0, ARG_TAB5},
	{(char *)"tab6", required_argument, 0, ARG_TAB6},
	{(char *)"interleaved", required_argument, 0, ARG_INTERLEAVED},
	{(char *)"phred33-quals", no_argument, 0, ARG_PHRED33},
	{(char *)"phred64-quals", no_argument, 0, ARG_PHRED64},
	{(char *)"phred33", no_argument, 0, ARG_PHRED33},
	{(char *)"phred64", no_argument, 0, ARG_PHRED64},
	{(char *)"solexa1.3-quals", no_argument, 0, ARG_PHRED64},
	{(char *)"mm", no_argument, 0, ARG_MM},
	{(char *)"shmem", no_argument, 0, ARG_SHMEM},
	{(char *)"mmsweep", no_argument, 0, ARG_MMSWEEP},
	{(char *)"hadoopout", no_argument, 0, ARG_HADOOPOUT},
	{(char *)"fullref", no_argument, 0, ARG_FULLREF},
	{(char *)"usage", no_argument, 0, ARG_USAGE},
	{(char *)"sam-no-qname-trunc", no_argument, 0, ARG_SAM_NO_QNAME_TRUNC},
	{(char *)"sam-omit-sec-seq", no_argument, 0, ARG_SAM_OMIT_SEC_SEQ},
	{(char *)"omit-sec-seq", no_argument, 0, ARG_SAM_OMIT_SEC_SEQ},
	{(char *)"sam-no-head", no_argument, 0, ARG_SAM_NOHEAD},
	{(char *)"sam-nohead", no_argument, 0, ARG_SAM_NOHEAD},
	{(char *)"sam-noHD", no_argument, 0, ARG_SAM_NOHEAD},
	{(char *)"sam-no-hd", no_argument, 0, ARG_SAM_NOHEAD},
	{(char *)"sam-nosq", no_argument, 0, ARG_SAM_NOSQ},
	{(char *)"sam-no-sq", no_argument, 0, ARG_SAM_NOSQ},
	{(char *)"sam-noSQ", no_argument, 0, ARG_SAM_NOSQ},
	{(char *)"no-head", no_argument, 0, ARG_SAM_NOHEAD},
	{(char *)"no-hd", no_argument, 0, ARG_SAM_NOHEAD},
	{(char *)"no-sq", no_argument, 0, ARG_SAM_NOSQ},
	{(char *)"no-HD", no_argument, 0, ARG_SAM_NOHEAD},
	{(char *)"no-SQ", no_argument, 0, ARG_SAM_NOSQ},
	{(char *)"no-unal", no_argument, 0, ARG_SAM_NO_UNAL},
	{(char *)"sam-RG", required_argument, 0, ARG_SAM_RG},
	{(char *)"sam-rg", required_argument, 0, ARG_SAM_RG},
	{(char *)"sam-rg-id", required_argument, 0, ARG_SAM_RGID},
	{(char *)"RG", required_argument, 0, ARG_SAM_RG},
	{(char *)"rg", required_argument, 0, ARG_SAM_RG},
	{(char *)"rg-id", required_argument, 0, ARG_SAM_RGID},
	{(char *)"snpphred", required_argument, 0, ARG_SNPPHRED},
	{(char *)"snpfrac", required_argument, 0, ARG_SNPFRAC},
	{(char *)"gbar", required_argument, 0, ARG_GAP_BAR},
	{(char *)"qseq", no_argument, 0, ARG_QSEQ},
	{(char *)"policy", required_argument, 0, ARG_ALIGN_POLICY},
	{(char *)"preset", required_argument, 0, 'P'},
	{(char *)"seed-summ", no_argument, 0, ARG_SEED_SUMM},
	{(char *)"seed-summary", no_argument, 0, ARG_SEED_SUMM},
	{(char *)"overhang", no_argument, 0, ARG_OVERHANG},
	{(char *)"no-cache", no_argument, 0, ARG_NO_CACHE},
	{(char *)"cache", no_argument, 0, ARG_USE_CACHE},
	{(char *)"454", no_argument, 0, ARG_NOISY_HPOLY},
	{(char *)"ion-torrent", no_argument, 0, ARG_NOISY_HPOLY},
	{(char *)"no-mixed", no_argument, 0, ARG_NO_MIXED},
	{(char *)"no-discordant", no_argument, 0, ARG_NO_DISCORDANT},
	{(char *)"local", no_argument, 0, ARG_LOCAL},
	{(char *)"end-to-end", no_argument, 0, ARG_END_TO_END},
	{(char *)"ungapped", no_argument, 0, ARG_UNGAPPED},
	{(char *)"no-ungapped", no_argument, 0, ARG_UNGAPPED_NO},
	{(char *)"sse8", no_argument, 0, ARG_SSE8},
	{(char *)"no-sse8", no_argument, 0, ARG_SSE8_NO},
	{(char *)"scan-narrowed", no_argument, 0, ARG_SCAN_NARROWED},
	{(char *)"qc-filter", no_argument, 0, ARG_QC_FILTER},
	{(char *)"bwa-sw-like", no_argument, 0, ARG_BWA_SW_LIKE},
	{(char *)"multiseed", required_argument, 0, ARG_MULTISEED_IVAL},
	{(char *)"ma", required_argument, 0, ARG_SCORE_MA},
	{(char *)"mp", required_argument, 0, ARG_SCORE_MMP},
	{(char *)"np", required_argument, 0, ARG_SCORE_NP},
	{(char *)"rdg", required_argument, 0, ARG_SCORE_RDG},
	{(char *)"rfg", required_argument, 0, ARG_SCORE_RFG},
	{(char *)"score-min", required_argument, 0, ARG_SCORE_MIN},
	{(char *)"min-score", required_argument, 0, ARG_SCORE_MIN},
	{(char *)"n-ceil", required_argument, 0, ARG_N_CEIL},
	{(char *)"dpad", required_argument, 0, ARG_DPAD},
	{(char *)"mapq-print-inputs", no_argument, 0, ARG_SAM_PRINT_YI},
	{(char *)"very-fast", no_argument, 0, ARG_PRESET_VERY_FAST},
	{(char *)"fast", no_argument, 0, ARG_PRESET_FAST},
	{(char *)"sensitive", no_argument, 0, ARG_PRESET_SENSITIVE},
	{(char *)"very-sensitive", no_argument, 0, ARG_PRESET_VERY_SENSITIVE},
	{(char *)"very-fast-local", no_argument, 0, ARG_PRESET_VERY_FAST_LOCAL},
	{(char *)"fast-local", no_argument, 0, ARG_PRESET_FAST_LOCAL},
	{(char *)"sensitive-local", no_argument, 0, ARG_PRESET_SENSITIVE_LOCAL},
	{(char *)"very-sensitive-local", no_argument, 0, ARG_PRESET_VERY_SENSITIVE_LOCAL},
	{(char *)"seedlen", required_argument, 0, 'L'},
	{(char *)"seedmms", required_argument, 0, 'N'},
	{(char *)"seedival", required_argument, 0, 'i'},
	{(char *)"ignore-quals", no_argument, 0, ARG_IGNORE_QUALS},
	{(char *)"index", required_argument, 0, 'x'},
	{(char *)"arg-desc", no_argument, 0, ARG_DESC},
	{(char *)"wrapper", required_argument, 0, ARG_WRAPPER},
	{(char *)"unpaired", required_argument, 0, 'U'},
	{(char *)"output", required_argument, 0, 'S'},
	{(char *)"mapq-v", required_argument, 0, ARG_MAPQ_V},
	{(char *)"dovetail", no_argument, 0, ARG_DOVETAIL},
	{(char *)"no-dovetail", no_argument, 0, ARG_NO_DOVETAIL},
	{(char *)"contain", no_argument, 0, ARG_CONTAIN},
	{(char *)"no-contain", no_argument, 0, ARG_NO_CONTAIN},
	{(char *)"overlap", no_argument, 0, ARG_OVERLAP},
	{(char *)"no-overlap", no_argument, 0, ARG_NO_OVERLAP},
	{(char *)"tighten", required_argument, 0, ARG_TIGHTEN},
	{(char *)"exact-upfront", no_argument, 0, ARG_EXACT_UPFRONT},
	{(char *)"1mm-upfront", no_argument, 0, ARG_1MM_UPFRONT},
	{(char *)"no-exact-upfront", no_argument, 0, ARG_EXACT_UPFRONT_NO},
	{(char *)"no-1mm-upfront", no_argument, 0, ARG_1MM_UPFRONT_NO},
	{(char *)"1mm-minlen", required_argument, 0, ARG_1MM_MINLEN},
	{(char *)"seed-off", required_argument, 0, 'O'},
	{(char *)"seed-boost", required_argument, 0, ARG_SEED_BOOST_THRESH},
	{(char *)"read-times", no_argument, 0, ARG_READ_TIMES},
	{(char *)"show-rand-seed", no_argument, 0, ARG_SHOW_RAND_SEED},
	{(char *)"dp-fail-streak", required_argument, 0, ARG_DP_FAIL_STREAK_THRESH},
	{(char *)"ee-fail-streak", required_argument, 0, ARG_EE_FAIL_STREAK_THRESH},
	{(char *)"ug-fail-streak", required_argument, 0, ARG_UG_FAIL_STREAK_THRESH},
	{(char *)"fail-streak", required_argument, 0, 'D'},
	{(char *)"dp-fails", required_argument, 0, ARG_DP_FAIL_THRESH},
	{(char *)"ug-fails", required_argument, 0, ARG_UG_FAIL_THRESH},
	{(char *)"extends", required_argument, 0, ARG_EXTEND_ITERS},
	{(char *)"no-extend", no_argument, 0, ARG_NO_EXTEND},
	{(char *)"mapq-extra", no_argument, 0, ARG_MAPQ_EX},
	{(char *)"seed-rounds", required_argument, 0, 'R'},
	{(char *)"reorder", no_argument, 0, ARG_REORDER},
	{(char *)"passthrough", no_argument, 0, ARG_READ_PASSTHRU},
	{(char *)"sample", required_argument, 0, ARG_SAMPLE},
	{(char *)"cp-min", required_argument, 0, ARG_CP_MIN},
	{(char *)"cp-ival", required_argument, 0, ARG_CP_IVAL},
	{(char *)"tri", no_argument, 0, ARG_TRI},
	{(char *)"nondeterministic", no_argument, 0, ARG_NON_DETERMINISTIC},
	{(char *)"non-deterministic", no_argument, 0, ARG_NON_DETERMINISTIC},
	{(char *)"local-seed-cache-sz", required_argument, 0, ARG_LOCAL_SEED_CACHE_SZ},
	{(char *)"seed-cache-sz", required_argument, 0, ARG_CURRENT_SEED_CACHE_SZ},
	{(char *)"no-unal", no_argument, 0, ARG_SAM_NO_UNAL},
	{(char *)"test-25", no_argument, 0, ARG_TEST_25},
	{(char *)"desc-kb", required_argument, 0, ARG_DESC_KB},
	{(char *)"desc-landing", required_argument, 0, ARG_DESC_LANDING},
	{(char *)"desc-exp", required_argument, 0, ARG_DESC_EXP},
	{(char *)"desc-prioritize", no_argument, 0, ARG_DESC_PRIORITIZE},
	{(char *)"desc-fmops", required_argument, 0, ARG_DESC_FMOPS},
	{(char *)"log-dp", required_argument, 0, ARG_LOG_DP},
	{(char *)"log-dp-opp", required_argument, 0, ARG_LOG_DP_OPP},
	{(char *)"soft-clipped-unmapped-tlen", no_argument, 0, ARG_SC_UNMAPPED},
	{(char *)"xeq", no_argument, 0, ARG_XEQ},
	{(char *)"thread-ceiling", required_argument, 0, ARG_THREAD_CEILING},
	{(char *)"thread-piddir", required_argument, 0, ARG_THREAD_PIDDIR},
	{(char *)"trim-to", required_argument, 0, ARG_TRIM_TO},
	{(char *)"preserve-tags", no_argument, 0, ARG_PRESERVE_TAGS},
	{(char *)"align-paired-reads", no_argument, 0, ARG_ALIGN_PAIRED_READS},
#ifdef USE_SRA
	{(char *)"sra-acc", required_argument, 0, ARG_SRA_ACC},
#endif
	{(char *)"sam-append-comment", no_argument, 0, ARG_SAM_APPEND_COMMENT},
	{(char *)0, 0, 0, 0}};
static void printArgDesc(ostream &out)
{
	size_t i = 0;
	while (long_options[i].name != 0)
	{
		out << long_options[i].name << "\t"
			<< (long_options[i].has_arg == no_argument ? 0 : 1)
			<< endl;
		i++;
	}
	size_t solen = strlen(short_options);
	for (i = 0; i < solen; i++)
	{
		if (i == solen - 1)
		{
			assert_neq(':', short_options[i]);
			cout << (char)short_options[i] << "\t" << 0 << endl;
		}
		else
		{
			if (short_options[i + 1] == ':')
			{
				cout << (char)short_options[i] << "\t" << 1 << endl;
				i++;
			}
			else
			{
				cout << (char)short_options[i] << "\t" << 0 << endl;
			}
		}
	}
}
static int parseInt(int lower, int upper, const char *errmsg, const char *arg)
{
	long l;
	char *endPtr = NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL)
	{
		if (l < lower || l > upper)
		{
			cerr << errmsg << endl;
			throw 1;
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	throw 1;
	return -1;
}
static int parseInt(int lower, const char *errmsg, const char *arg)
{
	return parseInt(lower, std::numeric_limits<int>::max(), errmsg, arg);
}
template <typename T>
T parse(const char *s)
{
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}
template <typename T>
pair<T, T> parsePair(const char *str, char delim)
{
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	pair<T, T> ret;
	ret.first = parse<T>(ss[0].c_str());
	ret.second = parse<T>(ss[1].c_str());
	return ret;
}
template <typename T>
void parseTuple(const char *str, char delim, EList<T> &ret)
{
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	for (size_t i = 0; i < ss.size(); i++)
	{
		ret.push_back(parse<T>(ss[i].c_str()));
	}
}
static string applyPreset(const string &sorig, Presets &presets)
{
	string s = sorig;
	size_t found = s.find("%LOCAL%");
	if (found != string::npos)
	{
		s.replace(found, strlen("%LOCAL%"), localAlign ? "-local" : "");
	}
	if (gVerbose)
	{
		cerr << "Applying preset: '" << s.c_str() << "' using preset menu '"
			 << presets.name() << "'" << endl;
	}
	string pol;
	presets.apply(s, pol, extra_opts);
	return pol;
}
static bool saw_M;
static bool saw_a;
static bool saw_k;
static bool saw_trim3;
static bool saw_trim5;
static bool saw_trim_to;
static bool saw_bam;
static bool saw_preserve_tags;
static bool saw_align_paired_reads;
static EList<string> presetList;
static void parseOption(int next_option, const char *arg)
{
	switch (next_option)
	{
	case ARG_TEST_25:
		nextVer = true;
		break;
	case ARG_DESC_KB:
		descentTotSz = SimpleFunc::parse(arg, 0.0, 1024.0, 1024.0, DMAX);
		break;
	case ARG_DESC_FMOPS:
		descentTotFmops = SimpleFunc::parse(arg, 0.0, 10.0, 100.0, DMAX);
		break;
	case ARG_LOG_DP:
		logDps = arg;
		break;
	case ARG_LOG_DP_OPP:
		logDpsOpp = arg;
		break;
	case ARG_DESC_LANDING:
	{
		descLanding = parse<int>(arg);
		if (descLanding < 1)
		{
			cerr << "Error: --desc-landing must be greater than or equal to 1" << endl;
			throw 1;
		}
		break;
	}
	case ARG_DESC_EXP:
	{
		descConsExp = parse<double>(arg);
		if (descConsExp < 0.0)
		{
			cerr << "Error: --desc-exp must be greater than or equal to 0" << endl;
			throw 1;
		}
		break;
	}
	case ARG_DESC_PRIORITIZE:
		descPrioritizeRoots = true;
		break;
	case '1':
		tokenize(arg, ",", mates1);
		break;
	case '2':
		tokenize(arg, ",", mates2);
		break;
	case ARG_ONETWO:
		tokenize(arg, ",", mates12);
		format = TAB_MATE5;
		break;
	case ARG_TAB5:
		tokenize(arg, ",", mates12);
		format = TAB_MATE5;
		break;
	case ARG_TAB6:
		tokenize(arg, ",", mates12);
		format = TAB_MATE6;
		break;
	case ARG_INTERLEAVED:
	{
		tokenize(arg, ",", mates12);
		interleaved = true;
		break;
	}
	case 'b':
	{
		format = BAM;
		saw_bam = true;
		break;
	}
	case 'f':
		format = FASTA;
		break;
	case 'F':
	{
		format = FASTA_CONT;
		pair<uint32_t, uint32_t> p = parsePair<uint32_t>(arg, ',');
		fastaContLen = p.first;
		fastaContFreq = p.second;
		break;
	}
	case ARG_BWA_SW_LIKE:
	{
		bwaSwLikeC = 5.5f;
		bwaSwLikeT = 30;
		bwaSwLike = true;
		localAlign = true;
		polstr += ";MA=1;MMP=C3;RDG=5,2;RFG=5,2";
		break;
	}
	case 'q':
		format = FASTQ;
		break;
	case 'r':
		format = RAW;
		break;
	case 'c':
		format = CMDLINE;
		break;
	case ARG_QSEQ:
		format = QSEQ;
		break;
	case 'I':
		gMinInsert = parseInt(0, "-I arg must be positive", arg);
		break;
	case 'X':
		gMaxInsert = parseInt(1, "-X arg must be at least 1", arg);
		break;
	case ARG_NO_DISCORDANT:
		gReportDiscordant = false;
		break;
	case ARG_NO_MIXED:
		gReportMixed = false;
		break;
	case 's':
		skipReads = (uint32_t)parseInt(0, "-s arg must be positive", arg);
		break;
	case ARG_FF:
		gMate1fw = true;
		gMate2fw = true;
		break;
	case ARG_RF:
		gMate1fw = false;
		gMate2fw = true;
		break;
	case ARG_FR:
		gMate1fw = true;
		gMate2fw = false;
		break;
	case ARG_SHMEM:
		useShmem = true;
		break;
	case ARG_SEED_SUMM:
		seedSumm = true;
		break;
	case ARG_SC_UNMAPPED:
		scUnMapped = true;
		break;
	case ARG_XEQ:
		xeq = true;
		break;
	case ARG_PRESERVE_TAGS:
	{
		preserve_tags = true;
		saw_preserve_tags = true;
		break;
	}
	case ARG_ALIGN_PAIRED_READS:
	{
		align_paired_reads = true;
		saw_align_paired_reads = true;
		break;
	}
	case ARG_MM:
	{
#ifdef _MM
		useMm = true;
		break;
#else
		cerr << "Memory-mapped I/O mode is disabled." << endl;
		throw 1;
#endif
	}
	case ARG_MMSWEEP:
		mmSweep = true;
		break;
	case ARG_HADOOPOUT:
		hadoopOut = true;
		break;
	case ARG_SOLEXA_QUALS:
		solexaQuals = true;
		break;
	case ARG_INTEGER_QUALS:
		integerQuals = true;
		break;
	case ARG_PHRED64:
		phred64Quals = true;
		break;
	case ARG_PHRED33:
		solexaQuals = false;
		phred64Quals = false;
		break;
	case ARG_OVERHANG:
		gReportOverhangs = true;
		break;
	case ARG_NO_CACHE:
		msNoCache = true;
		break;
	case ARG_USE_CACHE:
		msNoCache = false;
		break;
	case ARG_LOCAL_SEED_CACHE_SZ:
		seedCacheLocalMB = (uint32_t)parseInt(1, "--local-seed-cache-sz arg must be at least 1", arg);
		break;
	case ARG_CURRENT_SEED_CACHE_SZ:
		seedCacheCurrentMB = (uint32_t)parseInt(1, "--seed-cache-sz arg must be at least 1", arg);
		break;
	case ARG_REFIDX:
		noRefNames = true;
		break;
	case ARG_FULLREF:
		fullRef = true;
		break;
	case ARG_GAP_BAR:
		gGapBarrier = parseInt(1, "--gbar must be no less than 1", arg);
		break;
	case ARG_SEED:
		seed = parseInt(0, "--seed arg must be at least 0", arg);
		break;
	case ARG_NON_DETERMINISTIC:
		arbitraryRandom = true;
		break;
	case 'u':
		qUpto = (uint32_t)parseInt(1, "-u/--qupto arg must be at least 1", arg);
		break;
	case 'Q':
		tokenize(arg, ",", qualities);
		integerQuals = true;
		break;
	case ARG_QUALS1:
		tokenize(arg, ",", qualities1);
		integerQuals = true;
		break;
	case ARG_QUALS2:
		tokenize(arg, ",", qualities2);
		integerQuals = true;
		break;
	case ARG_CACHE_LIM:
		cacheLimit = (uint32_t)parseInt(1, "--cachelim arg must be at least 1", arg);
		break;
	case ARG_CACHE_SZ:
		cacheSize = (uint32_t)parseInt(1, "--cachesz arg must be at least 1", arg);
		cacheSize *= (1024 * 1024);
		break;
	case ARG_WRAPPER:
		wrapper = arg;
		break;
	case 'p':
		nthreads = parseInt(1, "-p/--threads arg must be at least 1", arg);
		break;
	case ARG_THREAD_CEILING:
		thread_ceiling = parseInt(0, "--thread-ceiling must be at least 0", arg);
		break;
	case ARG_THREAD_PIDDIR:
		thread_stealing_dir = arg;
		break;
	case ARG_FILEPAR:
		fileParallel = true;
		break;
	case '3':
		gTrim3 = parseInt(0, "-3/--trim3 arg must be at least 0", arg);
		break;
	case '5':
		gTrim5 = parseInt(0, "-5/--trim5 arg must be at least 0", arg);
		break;
	case ARG_TRIM_TO:
	{
		if (strlen(arg) > 1 && arg[1] != ':')
		{
			trimTo.first = 3;
			trimTo.second = parseInt(0, "--trim-to: the number of bases to trim must be at least 0", arg);
			break;
		}
		pair<int, int> res = parsePair<int>(arg, ':');
		if (res.first != 3 && res.first != 5)
		{
			cerr << "--trim-to: trim position must be either 3 or 5" << endl;
			throw 1;
		}
		if (res.second < 0)
		{
			cerr << "--trim-to: the number bases to trim must be at least 0" << endl;
			throw 1;
		}
		trimTo = static_cast<pair<short, size_t>>(res);
		break;
	}
	case 'M':
	{
		msample = true;
		mhits = parse<uint32_t>(arg);
		if (saw_a || saw_k)
		{
			cerr << "Warning: -M, -k and -a are mutually exclusive. "
				 << "-M will override" << endl;
			khits = 1;
		}
		assert_eq(1, khits);
		saw_M = true;
		cerr << "Warning: -M is deprecated.  Use -D and -R to adjust "
			 << "effort instead." << endl;
		break;
	}
	case ARG_EXTEND_ITERS:
	{
		maxIters = parse<size_t>(arg);
		break;
	}
	case ARG_NO_EXTEND:
	{
		doExtend = false;
		break;
	}
	case 'R':
	{
		polstr += ";ROUNDS=";
		polstr += arg;
		break;
	}
	case 'D':
	{
		polstr += ";DPS=";
		polstr += arg;
		break;
	}
	case ARG_DP_MATE_STREAK_THRESH:
	{
		maxMateStreak = parse<size_t>(arg);
		break;
	}
	case ARG_DP_FAIL_STREAK_THRESH:
	{
		maxDpStreak = parse<size_t>(arg);
		break;
	}
	case ARG_EE_FAIL_STREAK_THRESH:
	{
		maxEeStreak = parse<size_t>(arg);
		break;
	}
	case ARG_UG_FAIL_STREAK_THRESH:
	{
		maxUgStreak = parse<size_t>(arg);
		break;
	}
	case ARG_DP_FAIL_THRESH:
	{
		maxDp = parse<size_t>(arg);
		break;
	}
	case ARG_UG_FAIL_THRESH:
	{
		maxUg = parse<size_t>(arg);
		break;
	}
	case ARG_SEED_BOOST_THRESH:
	{
		seedBoostThresh = parse<int>(arg);
		break;
	}
	case 'a':
	{
		msample = false;
		allHits = true;
		mhits = 0;
		if (saw_M || saw_k)
		{
			cerr << "Warning: -M, -k and -a are mutually exclusive. "
				 << "-a will override" << endl;
		}
		saw_a = true;
		break;
	}
	case 'k':
	{
		msample = false;
		khits = (uint32_t)parseInt(1, "-k arg must be at least 1", arg);
		mhits = 0;
		if (saw_M || saw_a)
		{
			cerr << "Warning: -M, -k and -a are mutually exclusive. "
				 << "-k will override" << endl;
		}
		saw_k = true;
		break;
	}
	case ARG_VERBOSE:
		gVerbose = 1;
		break;
	case ARG_STARTVERBOSE:
		startVerbose = true;
		break;
	case ARG_QUIET:
		gQuiet = true;
		break;
	case ARG_SANITY:
		sanityCheck = true;
		break;
	case 't':
		timing = true;
		break;
	case ARG_METRIC_IVAL:
	{
		metricsIval = parseInt(1, "--metrics arg must be at least 1", arg);
		break;
	}
	case ARG_METRIC_FILE:
		metricsFile = arg;
		break;
	case ARG_METRIC_STDERR:
		metricsStderr = true;
		break;
	case ARG_METRIC_PER_READ:
		metricsPerRead = true;
		break;
	case ARG_NO_FW:
		gNofw = true;
		break;
	case ARG_NO_RC:
		gNorc = true;
		break;
	case ARG_SAM_NO_QNAME_TRUNC:
		samTruncQname = false;
		break;
	case ARG_SAM_APPEND_COMMENT:
		samAppendComment = true;
		break;
	case ARG_SAM_OMIT_SEC_SEQ:
		samOmitSecSeqQual = true;
		break;
	case ARG_SAM_NO_UNAL:
		samNoUnal = true;
		break;
	case ARG_SAM_NOHEAD:
		samNoHead = true;
		break;
	case ARG_SAM_NOSQ:
		samNoSQ = true;
		break;
	case ARG_SAM_PRINT_YI:
		sam_print_yi = true;
		break;
	case ARG_REORDER:
		reorder = true;
		break;
	case ARG_MAPQ_EX:
	{
		sam_print_zt = true;
		break;
	}
	case ARG_SHOW_RAND_SEED:
	{
		sam_print_zs = true;
		break;
	}
	case ARG_SAMPLE:
		sampleFrac = parse<float>(arg);
		break;
	case ARG_CP_MIN:
		cminlen = parse<size_t>(arg);
		break;
	case ARG_CP_IVAL:
		cpow2 = parse<size_t>(arg);
		break;
	case ARG_TRI:
		doTri = true;
		break;
	case ARG_READ_PASSTHRU:
	{
		sam_print_xr = true;
		break;
	}
	case ARG_READ_TIMES:
	{
		sam_print_xt = true;
		sam_print_xd = true;
		sam_print_xu = true;
		sam_print_yl = true;
		sam_print_ye = true;
		sam_print_yu = true;
		sam_print_yr = true;
		sam_print_zb = true;
		sam_print_zr = true;
		sam_print_zf = true;
		sam_print_zm = true;
		sam_print_zi = true;
		break;
	}
	case ARG_SAM_RG:
	{
		string argstr = arg;
		if (argstr.substr(0, 3) == "ID:")
		{
			rgid = "\t";
			rgid += argstr;
			rgs_optflag = "RG:Z:" + argstr.substr(3);
		}
		else
		{
			rgs += '\t';
			rgs += argstr;
		}
		break;
	}
	case ARG_SAM_RGID:
	{
		string argstr = arg;
		rgid = "\t";
		rgid = "\tID:" + argstr;
		rgs_optflag = "RG:Z:" + argstr;
		break;
	}
	case ARG_PARTITION:
		partitionSz = parse<int>(arg);
		break;
	case ARG_READS_PER_BATCH:
		readsPerBatch = parseInt(1, "--reads-per-batch arg must be at least 1", arg);
		break;
	case ARG_DPAD:
		maxhalf = parseInt(0, "--dpad must be no less than 0", arg);
		break;
	case ARG_ORIG:
		if (arg == NULL || strlen(arg) == 0)
		{
			cerr << "--orig arg must be followed by a string" << endl;
			throw 1;
		}
		origString = arg;
		break;
	case ARG_LOCAL:
	{
		localAlign = true;
		gDefaultSeedLen = DEFAULT_LOCAL_SEEDLEN;
		break;
	}
	case ARG_END_TO_END:
		localAlign = false;
		break;
	case ARG_SSE8:
		enable8 = true;
		break;
	case ARG_SSE8_NO:
		enable8 = false;
		break;
	case ARG_UNGAPPED:
		doUngapped = true;
		break;
	case ARG_UNGAPPED_NO:
		doUngapped = false;
		break;
	case ARG_NO_DOVETAIL:
		gDovetailMatesOK = false;
		break;
	case ARG_NO_CONTAIN:
		gContainMatesOK = false;
		break;
	case ARG_NO_OVERLAP:
		gOlapMatesOK = false;
		break;
	case ARG_DOVETAIL:
		gDovetailMatesOK = true;
		break;
	case ARG_CONTAIN:
		gContainMatesOK = true;
		break;
	case ARG_OVERLAP:
		gOlapMatesOK = true;
		break;
	case ARG_QC_FILTER:
		qcFilter = true;
		break;
	case ARG_IGNORE_QUALS:
		ignoreQuals = true;
		break;
	case ARG_MAPQ_V:
		mapqv = parse<int>(arg);
		break;
	case ARG_TIGHTEN:
		tighten = parse<int>(arg);
		break;
	case ARG_EXACT_UPFRONT:
		doExactUpFront = true;
		break;
	case ARG_1MM_UPFRONT:
		do1mmUpFront = true;
		break;
	case ARG_EXACT_UPFRONT_NO:
		doExactUpFront = false;
		break;
	case ARG_1MM_UPFRONT_NO:
		do1mmUpFront = false;
		break;
	case ARG_1MM_MINLEN:
		do1mmMinLen = parse<size_t>(arg);
		break;
	case ARG_NOISY_HPOLY:
		noisyHpolymer = true;
		break;
	case 'x':
		bt2index = arg;
		break;
	case ARG_PRESET_VERY_FAST_LOCAL:
		localAlign = true;
	case ARG_PRESET_VERY_FAST:
	{
		presetList.push_back("very-fast%LOCAL%");
		break;
	}
	case ARG_PRESET_FAST_LOCAL:
		localAlign = true;
	case ARG_PRESET_FAST:
	{
		presetList.push_back("fast%LOCAL%");
		break;
	}
	case ARG_PRESET_SENSITIVE_LOCAL:
		localAlign = true;
	case ARG_PRESET_SENSITIVE:
	{
		presetList.push_back("sensitive%LOCAL%");
		break;
	}
	case ARG_PRESET_VERY_SENSITIVE_LOCAL:
		localAlign = true;
	case ARG_PRESET_VERY_SENSITIVE:
	{
		presetList.push_back("very-sensitive%LOCAL%");
		break;
	}
	case 'P':
	{
		presetList.push_back(arg);
		break;
	}
	case ARG_ALIGN_POLICY:
	{
		if (strlen(arg) > 0)
		{
			polstr += ";";
			polstr += arg;
		}
		break;
	}
	case 'N':
	{
		int64_t len = parse<size_t>(arg);
		if (len < 0 || len > 1)
		{
			cerr << "Error: -N argument must be within the interval [0,1]; was " << arg << endl;
			throw 1;
		}
		polstr += ";SEED=";
		polstr += arg;
		break;
	}
	case 'L':
	{
		int64_t len = parse<size_t>(arg);
		if (len < 1 || len > 32)
		{
			cerr << "Error: -L argument must be within the interval [1,32]; was " << arg << endl;
			throw 1;
		}
		polstr += ";SEEDLEN=";
		polstr += arg;
		break;
	}
	case 'O':
		multiseedOff = parse<size_t>(arg);
		break;
	case 'i':
	{
		EList<string> args;
		tokenize(arg, ",", args);
		if (args.size() > 3 || args.size() == 0)
		{
			cerr << "Error: expected 3 or fewer comma-separated "
				 << "arguments to -i option, got "
				 << args.size() << endl;
			throw 1;
		}
		polstr += (";IVAL=" + args[0]);
		if (args.size() > 1)
		{
			polstr += ("," + args[1]);
		}
		if (args.size() > 2)
		{
			polstr += ("," + args[2]);
		}
		break;
	}
	case ARG_MULTISEED_IVAL:
	{
		polstr += ";";
		EList<string> args;
		tokenize(arg, ",", args);
		if (args.size() > 5 || args.size() == 0)
		{
			cerr << "Error: expected 5 or fewer comma-separated "
				 << "arguments to --multiseed option, got "
				 << args.size() << endl;
			throw 1;
		}
		polstr += "SEED=";
		polstr += (args[0]);
		if (args.size() > 1)
			polstr += (";SEEDLEN=" + args[1]);
		if (args.size() > 2)
			polstr += (";IVAL=" + args[2]);
		if (args.size() > 3)
			polstr += ("," + args[3]);
		if (args.size() > 4)
			polstr += ("," + args[4]);
		break;
	}
	case ARG_N_CEIL:
	{
		EList<string> args;
		tokenize(arg, ",", args);
		if (args.size() > 3)
		{
			cerr << "Error: expected 3 or fewer comma-separated "
				 << "arguments to --n-ceil option, got "
				 << args.size() << endl;
			throw 1;
		}
		if (args.size() == 0)
		{
			cerr << "Error: expected at least one argument to --n-ceil option" << endl;
			throw 1;
		}
		polstr += ";NCEIL=";
		if (args.size() == 3)
		{
			polstr += (args[0] + "," + args[1] + "," + args[2]);
		}
		else
		{
			polstr += ("L," + args[0]);
			if (args.size() > 1)
			{
				polstr += ("," + (args[1]));
			}
		}
		break;
	}
	case ARG_SCORE_MA:
		polstr += ";MA=";
		polstr += arg;
		break;
	case ARG_SCORE_MMP:
	{
		EList<string> args;
		tokenize(arg, ",", args);
		if (args.size() > 2 || args.size() == 0)
		{
			cerr << "Error: expected 1 or 2 comma-separated "
				 << "arguments to --mmp option, got " << args.size() << endl;
			throw 1;
		}
		if (args.size() >= 1)
		{
			polstr += ";MMP=Q,";
			polstr += args[0];
			if (args.size() >= 2)
			{
				polstr += ",";
				polstr += args[1];
			}
		}
		break;
	}
	case ARG_SCORE_NP:
		polstr += ";NP=C";
		polstr += arg;
		break;
	case ARG_SCORE_RDG:
		polstr += ";RDG=";
		polstr += arg;
		break;
	case ARG_SCORE_RFG:
		polstr += ";RFG=";
		polstr += arg;
		break;
	case ARG_SCORE_MIN:
	{
		polstr += ";";
		EList<string> args;
		tokenize(arg, ",", args);
		if (args.size() > 3 || args.size() == 0)
		{
			cerr << "Error: expected 3 or fewer comma-separated "
				 << "arguments to --n-ceil option, got "
				 << args.size() << endl;
			throw 1;
		}
		polstr += ("MIN=" + args[0]);
		if (args.size() > 1)
		{
			polstr += ("," + args[1]);
		}
		if (args.size() > 2)
		{
			polstr += ("," + args[2]);
		}
		break;
	}
	case ARG_DESC:
		printArgDesc(cout);
		throw 0;
	case 'S':
		outfile = arg;
		break;
	case 'U':
	{
		EList<string> args;
		tokenize(arg, ",", args);
		for (size_t i = 0; i < args.size(); i++)
		{
			queries.push_back(args[i]);
		}
		break;
	}
#ifdef USE_SRA
	case ARG_SRA_ACC:
	{
		tokenize(arg, ",", sra_accs);
		format = SRA_FASTA;
		break;
	}
#endif
	case ARG_VERSION:
		showVersion = 1;
		break;
	default:
		throw 1;
	}
}
static void parseOptions(int argc, const char **argv)
{
	int option_index = 0;
	int next_option;
	saw_M = false;
	saw_a = false;
	saw_k = false;
	saw_trim3 = false;
	saw_trim5 = false;
	saw_trim_to = false;
	saw_bam = false;
	saw_preserve_tags = false;
	saw_align_paired_reads = false;
	presetList.clear();
	if (startVerbose)
	{
		cerr << "Parsing options: ";
		logTime(cerr, true);
	}
	while (true)
	{
		next_option = getopt_long(
			argc, const_cast<char **>(argv),
			short_options, long_options, &option_index);
		const char *arg = optarg;
		if (next_option == EOF)
		{
			if (extra_opts_cur < extra_opts.size())
			{
				next_option = extra_opts[extra_opts_cur].first;
				arg = extra_opts[extra_opts_cur].second.c_str();
				extra_opts_cur++;
			}
			else
			{
				break;
			}
		}
		parseOption(next_option, arg);
	}
	if (!localAlign && scUnMapped)
	{
		cerr << "ERROR: --soft-clipped-unmapped-tlen can only be set for local alignments." << endl;
		exit(1);
	}
	if ((saw_trim3 || saw_trim5) && saw_trim_to)
	{
		cerr << "ERROR: --trim5/--trim3 and --trim-to are mutually exclusive "
			 << "options." << endl;
		exit(1);
	}
	if (!saw_bam && saw_preserve_tags)
	{
		cerr << "--preserve_tags can only be used when aligning BAM reads." << endl;
		exit(1);
	}
	if (!saw_bam && saw_align_paired_reads)
	{
		cerr << "--align-paired-reads can only be used when aligning BAM reads." << endl;
		exit(1);
	}
	unique_ptr<Presets> presets(new PresetsV0());
	if (presetList.empty())
		polstr = applyPreset(defaultPreset, *presets.get()) + polstr;
	else
	{
		for (size_t i = presetList.size(); i != 0; i--)
			polstr = applyPreset(presetList[i - 1], *presets.get()) + polstr;
	}
	for (size_t i = 0; i < presetList.size(); i++)
	{
		polstr += applyPreset(presetList[i], *presets.get());
	}
	for (size_t i = 0; i < extra_opts.size(); i++)
	{
		next_option = extra_opts[extra_opts_cur].first;
		const char *arg = extra_opts[extra_opts_cur].second.c_str();
		parseOption(next_option, arg);
	}
	while (!polstr.empty() && polstr[0] == ';')
	{
		polstr = polstr.substr(1);
	}
	if (gVerbose)
	{
		cerr << "Final policy string: '" << polstr.c_str() << "'" << endl;
	}
	size_t failStreakTmp = 0;
	SeedAlignmentPolicy::parseString(
		polstr,
		localAlign,
		noisyHpolymer,
		ignoreQuals,
		bonusMatchType,
		bonusMatch,
		penMmcType,
		penMmcMax,
		penMmcMin,
		penNType,
		penN,
		penRdGapConst,
		penRfGapConst,
		penRdGapLinear,
		penRfGapLinear,
		scoreMin,
		nCeil,
		penNCatPair,
		multiseedMms,
		multiseedLen,
		msIval,
		failStreakTmp,
		nSeedRounds);
	if (failStreakTmp > 0)
	{
		maxEeStreak = failStreakTmp;
		maxUgStreak = failStreakTmp;
		maxDpStreak = failStreakTmp;
	}
	if (saw_a || saw_k)
	{
		msample = false;
		mhits = 0;
	}
	else
	{
		assert_gt(mhits, 0);
		msample = true;
	}
	if (mates1.size() != mates2.size())
	{
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << mates2.size() << endl
			 << "mate files/sequences were specified with -2.  The same number of mate files/" << endl
			 << "sequences must be specified with -1 and -2." << endl;
		throw 1;
	}
	if (interleaved && (format != FASTA && format != FASTQ))
	{
		cerr << "Error: --interleaved only works in combination with FASTA (-f) and FASTQ (-q) formats." << endl;
		throw 1;
	}
	if (samAppendComment && (format != FASTA && format != FASTQ))
	{
		cerr << "Error --sam-append-comment only works with FASTA (-f) and FASTQ (-q) formats. " << endl;
		throw 1;
	}
	if (qualities.size() && format != FASTA)
	{
		cerr << "Error: one or more quality files were specified with -Q but -f was not" << endl
			 << "enabled.  -Q works only in combination with -f and -C." << endl;
		throw 1;
	}
	if (qualities1.size() && format != FASTA)
	{
		cerr << "Error: one or more quality files were specified with --Q1 but -f was not" << endl
			 << "enabled.  --Q1 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if (qualities2.size() && format != FASTA)
	{
		cerr << "Error: one or more quality files were specified with --Q2 but -f was not" << endl
			 << "enabled.  --Q2 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if (qualities1.size() > 0 && mates1.size() != qualities1.size())
	{
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << qualities1.size() << endl
			 << "quality files were specified with --Q1.  The same number of mate and quality" << endl
			 << "files must sequences must be specified with -1 and --Q1." << endl;
		throw 1;
	}
	if (qualities2.size() > 0 && mates2.size() != qualities2.size())
	{
		cerr << "Error: " << mates2.size() << " mate files/sequences were specified with -2, but " << qualities2.size() << endl
			 << "quality files were specified with --Q2.  The same number of mate and quality" << endl
			 << "files must sequences must be specified with -2 and --Q2." << endl;
		throw 1;
	}
	if (!rgs.empty() && rgid.empty())
	{
		cerr << "Warning: --rg was specified without --rg-id also "
			 << "being specified.  @RG line is not printed unless --rg-id "
			 << "is specified." << endl;
	}
	if (format != CMDLINE)
	{
		for (size_t i = 0; i < mates1.size(); i++)
		{
			for (size_t j = 0; j < mates2.size(); j++)
			{
				if (mates1[i] == mates2[j] && !gQuiet)
				{
					cerr << "Warning: Same mate file \"" << mates1[i].c_str() << "\" appears as argument to both -1 and -2" << endl;
				}
			}
		}
	}
	if (qUpto + skipReads > qUpto)
	{
		qUpto += skipReads;
	}
	if (useShmem && useMm && !gQuiet)
	{
		cerr << "Warning: --shmem overrides --mm..." << endl;
		useMm = false;
	}
	if (gGapBarrier < 1)
	{
		cerr << "Warning: --gbar was set less than 1 (=" << gGapBarrier
			 << "); setting to 1 instead" << endl;
		gGapBarrier = 1;
	}
	if (bonusMatch > 0 && !scoreMin.alwaysPositive())
	{
		cerr << "Error: the match penalty is greater than 0 (" << bonusMatch
			 << ") but the --score-min function can be less than or equal to "
			 << "zero.  Either let the match penalty be 0 or make --score-min "
			 << "always positive." << endl;
		throw 1;
	}
	if (multiseedMms >= multiseedLen)
	{
		assert_gt(multiseedLen, 0);
		cerr << "Warning: seed mismatches (" << multiseedMms
			 << ") is less than seed length (" << multiseedLen
			 << "); setting mismatches to " << (multiseedMms - 1)
			 << " instead" << endl;
		multiseedMms = multiseedLen - 1;
	}
	sam_print_zm = sam_print_zm && nextVer;
#ifndef NDEBUG
	if (!gQuiet)
	{
		cerr << "Warning: Running in debug mode." << endl;
	}
#endif
}
static const char *argv0 = NULL;
static PatternSourcePerThreadFactory *
createPatsrcFactory(
	PatternComposer &patcomp,
	const PatternParams &pp,
	int tid)
{
	PatternSourcePerThreadFactory *patsrcFact;
	patsrcFact = new PatternSourcePerThreadFactory(patcomp, pp, tid);
	assert(patsrcFact != NULL);
	return patsrcFact;
}
#define PTHREAD_ATTRS (PTHREAD_CREATE_JOINABLE | PTHREAD_CREATE_DETACHED)
static PatternComposer *multiseed_patsrc;
static PatternParams multiseed_pp;
static Ebwt *multiseed_ebwtFw;
static Ebwt *multiseed_ebwtBw;
static Scoring *multiseed_sc;
static BitPairReference *multiseed_refs;
static AlignmentCache *multiseed_ca;
static AlnSink *multiseed_msink;
static OutFileBuf *multiseed_metricsOfb;
struct OuterLoopMetrics
{
	OuterLoopMetrics()
	{
		reset();
	}
	void reset()
	{
		reads = bases = srreads = srbases =
			freads = fbases = ureads = ubases = 0;
	}
	void merge(const OuterLoopMetrics &m)
	{
		reads += m.reads;
		bases += m.bases;
		srreads += m.srreads;
		srbases += m.srbases;
		freads += m.freads;
		fbases += m.fbases;
		ureads += m.ureads;
		ubases += m.ubases;
	}
	uint64_t reads;
	uint64_t bases;
	uint64_t srreads;
	uint64_t srbases;
	uint64_t freads;
	uint64_t fbases;
	uint64_t ureads;
	uint64_t ubases;
	MUTEX_T mutex_m;
};
struct PerfMetrics
{
	PerfMetrics() : first(true) { reset(); }
	void reset()
	{
		olm.reset();
		sdm.reset();
		wlm.reset();
		swmSeed.reset();
		swmMate.reset();
		rpm.reset();
		dpSse8Seed.reset();
		dpSse8Mate.reset();
		dpSse16Seed.reset();
		dpSse16Mate.reset();
		nbtfiltst = 0;
		nbtfiltsc = 0;
		nbtfiltdo = 0;
		olmu.reset();
		sdmu.reset();
		wlmu.reset();
		swmuSeed.reset();
		swmuMate.reset();
		rpmu.reset();
		dpSse8uSeed.reset();
		dpSse8uMate.reset();
		dpSse16uSeed.reset();
		dpSse16uMate.reset();
		nbtfiltst_u = 0;
		nbtfiltsc_u = 0;
		nbtfiltdo_u = 0;
	}
	void merge(
		const OuterLoopMetrics *ol,
		const SeedSearchMetrics *sd,
		const WalkMetrics *wl,
		const SwMetrics *swSeed,
		const SwMetrics *swMate,
		const ReportingMetrics *rm,
		const SSEMetrics *dpSse8Ex,
		const SSEMetrics *dpSse8Ma,
		const SSEMetrics *dpSse16Ex,
		const SSEMetrics *dpSse16Ma,
		uint64_t nbtfiltst_,
		uint64_t nbtfiltsc_,
		uint64_t nbtfiltdo_)
	{
		ThreadSafe ts(mutex_m);
		if (ol != NULL)
		{
			olmu.merge(*ol);
		}
		if (sd != NULL)
		{
			sdmu.merge(*sd);
		}
		if (wl != NULL)
		{
			wlmu.merge(*wl);
		}
		if (swSeed != NULL)
		{
			swmuSeed.merge(*swSeed);
		}
		if (swMate != NULL)
		{
			swmuMate.merge(*swMate);
		}
		if (rm != NULL)
		{
			rpmu.merge(*rm);
		}
		if (dpSse8Ex != NULL)
		{
			dpSse8uSeed.merge(*dpSse8Ex);
		}
		if (dpSse8Ma != NULL)
		{
			dpSse8uMate.merge(*dpSse8Ma);
		}
		if (dpSse16Ex != NULL)
		{
			dpSse16uSeed.merge(*dpSse16Ex);
		}
		if (dpSse16Ma != NULL)
		{
			dpSse16uMate.merge(*dpSse16Ma);
		}
		nbtfiltst_u += nbtfiltst_;
		nbtfiltsc_u += nbtfiltsc_;
		nbtfiltdo_u += nbtfiltdo_;
	}
	void reportInterval(
		OutFileBuf *o,
		bool metricsStderr,
		bool total,
		const BTString *name)
	{
		ThreadSafe ts(mutex_m);
		ostringstream stderrSs;
		time_t curtime = time(0);
		char buf[1024];
		if (first)
		{
			const char *str =
				"Time"
				"\t"
				"Read"
				"\t"
				"Base"
				"\t"
				"SameRead"
				"\t"
				"SameReadBase"
				"\t"
				"UnfilteredRead"
				"\t"
				"UnfilteredBase"
				"\t"
				"Paired"
				"\t"
				"Unpaired"
				"\t"
				"AlConUni"
				"\t"
				"AlConRep"
				"\t"
				"AlConFail"
				"\t"
				"AlDis"
				"\t"
				"AlConFailUni"
				"\t"
				"AlConFailRep"
				"\t"
				"AlConFailFail"
				"\t"
				"AlConRepUni"
				"\t"
				"AlConRepRep"
				"\t"
				"AlConRepFail"
				"\t"
				"AlUnpUni"
				"\t"
				"AlUnpRep"
				"\t"
				"AlUnpFail"
				"\t"
				"SeedSearch"
				"\t"
				"NRange"
				"\t"
				"NElt"
				"\t"
				"IntraSCacheHit"
				"\t"
				"InterSCacheHit"
				"\t"
				"OutOfMemory"
				"\t"
				"AlBWOp"
				"\t"
				"AlBWBranch"
				"\t"
				"ResBWOp"
				"\t"
				"ResBWBranch"
				"\t"
				"ResResolve"
				"\t"
				"ResReport"
				"\t"
				"RedundantSHit"
				"\t"
				"BestMinEdit0"
				"\t"
				"BestMinEdit1"
				"\t"
				"BestMinEdit2"
				"\t"
				"ExactAttempts"
				"\t"
				"ExactSucc"
				"\t"
				"ExactRanges"
				"\t"
				"ExactRows"
				"\t"
				"ExactOOMs"
				"\t"
				"1mmAttempts"
				"\t"
				"1mmSucc"
				"\t"
				"1mmRanges"
				"\t"
				"1mmRows"
				"\t"
				"1mmOOMs"
				"\t"
				"UngappedSucc"
				"\t"
				"UngappedFail"
				"\t"
				"UngappedNoDec"
				"\t"
				"DPExLt10Gaps"
				"\t"
				"DPExLt5Gaps"
				"\t"
				"DPExLt3Gaps"
				"\t"
				"DPMateLt10Gaps"
				"\t"
				"DPMateLt5Gaps"
				"\t"
				"DPMateLt3Gaps"
				"\t"
				"DP16ExDps"
				"\t"
				"DP16ExDpSat"
				"\t"
				"DP16ExDpFail"
				"\t"
				"DP16ExDpSucc"
				"\t"
				"DP16ExCol"
				"\t"
				"DP16ExCell"
				"\t"
				"DP16ExInner"
				"\t"
				"DP16ExFixup"
				"\t"
				"DP16ExGathSol"
				"\t"
				"DP16ExBt"
				"\t"
				"DP16ExBtFail"
				"\t"
				"DP16ExBtSucc"
				"\t"
				"DP16ExBtCell"
				"\t"
				"DP16ExCoreRej"
				"\t"
				"DP16ExNRej"
				"\t"
				"DP8ExDps"
				"\t"
				"DP8ExDpSat"
				"\t"
				"DP8ExDpFail"
				"\t"
				"DP8ExDpSucc"
				"\t"
				"DP8ExCol"
				"\t"
				"DP8ExCell"
				"\t"
				"DP8ExInner"
				"\t"
				"DP8ExFixup"
				"\t"
				"DP8ExGathSol"
				"\t"
				"DP8ExBt"
				"\t"
				"DP8ExBtFail"
				"\t"
				"DP8ExBtSucc"
				"\t"
				"DP8ExBtCell"
				"\t"
				"DP8ExCoreRej"
				"\t"
				"DP8ExNRej"
				"\t"
				"DP16MateDps"
				"\t"
				"DP16MateDpSat"
				"\t"
				"DP16MateDpFail"
				"\t"
				"DP16MateDpSucc"
				"\t"
				"DP16MateCol"
				"\t"
				"DP16MateCell"
				"\t"
				"DP16MateInner"
				"\t"
				"DP16MateFixup"
				"\t"
				"DP16MateGathSol"
				"\t"
				"DP16MateBt"
				"\t"
				"DP16MateBtFail"
				"\t"
				"DP16MateBtSucc"
				"\t"
				"DP16MateBtCell"
				"\t"
				"DP16MateCoreRej"
				"\t"
				"DP16MateNRej"
				"\t"
				"DP8MateDps"
				"\t"
				"DP8MateDpSat"
				"\t"
				"DP8MateDpFail"
				"\t"
				"DP8MateDpSucc"
				"\t"
				"DP8MateCol"
				"\t"
				"DP8MateCell"
				"\t"
				"DP8MateInner"
				"\t"
				"DP8MateFixup"
				"\t"
				"DP8MateGathSol"
				"\t"
				"DP8MateBt"
				"\t"
				"DP8MateBtFail"
				"\t"
				"DP8MateBtSucc"
				"\t"
				"DP8MateBtCell"
				"\t"
				"DP8MateCoreRej"
				"\t"
				"DP8MateNRej"
				"\t"
				"DPBtFiltStart"
				"\t"
				"DPBtFiltScore"
				"\t"
				"DpBtFiltDom"
				"\t"
#ifdef USE_MEM_TALLY
				"MemPeak"
				"\t"
				"UncatMemPeak"
				"\t"
				"EbwtMemPeak"
				"\t"
				"CacheMemPeak"
				"\t"
				"ResolveMemPeak"
				"\t"
				"AlignMemPeak"
				"\t"
				"DPMemPeak"
				"\t"
				"MiscMemPeak"
				"\t"
				"DebugMemPeak"
				"\t"
#endif
				"\n";
			if (name != NULL)
			{
				if (o != NULL)
					o->writeChars("Name\t");
				if (metricsStderr)
					stderrSs << "Name\t";
			}
			if (o != NULL)
				o->writeChars(str);
			if (metricsStderr)
				stderrSs << str;
			first = false;
		}
		if (total)
			mergeIncrementals();
		if (name != NULL)
		{
			if (o != NULL)
			{
				o->writeChars(name->toZBuf());
				o->write('\t');
			}
			if (metricsStderr)
			{
				stderrSs << (*name) << '\t';
			}
		}
		itoa10<time_t>(curtime, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		const OuterLoopMetrics &ol = total ? olm : olmu;
		itoa10<uint64_t>(ol.reads, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(ol.bases, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(ol.srreads, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(ol.srbases, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(ol.ureads, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(ol.ubases, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		const ReportingMetrics &rp = total ? rpm : rpmu;
		itoa10<uint64_t>(rp.npaired, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunpaired, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nconcord_uni, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nconcord_rep, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nconcord_0, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.ndiscord, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_0_uni, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_0_rep, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_0_0, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_rep_uni, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_rep_rep, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_rep_0, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_uni, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_rep, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(rp.nunp_0, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		const SeedSearchMetrics &sd = total ? sdm : sdmu;
		itoa10<uint64_t>(sd.seedsearch, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(sd.nrange, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(sd.nelt, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(sd.intrahit, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(sd.interhit, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(sd.ooms, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(sd.bwops, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(sd.bweds, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		const WalkMetrics &wl = total ? wlm : wlmu;
		itoa10<uint64_t>(wl.bwops, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(wl.branches, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(wl.resolves, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(wl.reports, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.rshit : swmuSeed.rshit, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? sdm.bestmin0 : sdmu.bestmin0, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? sdm.bestmin1 : sdmu.bestmin1, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? sdm.bestmin2 : sdmu.bestmin2, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.exatts : swmuSeed.exatts, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.exsucc : swmuSeed.exsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.exranges : swmuSeed.exranges, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.exrows : swmuSeed.exrows, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.exooms : swmuSeed.exooms, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.mm1atts : swmuSeed.mm1atts, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.mm1succ : swmuSeed.mm1succ, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.mm1ranges : swmuSeed.mm1ranges, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.mm1rows : swmuSeed.mm1rows, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.mm1ooms : swmuSeed.mm1ooms, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.ungapsucc : swmuSeed.ungapsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.ungapfail : swmuSeed.ungapfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.ungapnodec : swmuSeed.ungapnodec, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.sws10 : swmuSeed.sws10, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.sws5 : swmuSeed.sws5, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmSeed.sws3 : swmuSeed.sws3, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmMate.sws10 : swmuMate.sws10, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmMate.sws5 : swmuMate.sws5, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? swmMate.sws3 : swmuMate.sws3, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		const SSEMetrics &dpSse16s = total ? dpSse16Seed : dpSse16uSeed;
		itoa10<uint64_t>(dpSse16s.dp, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.dpsat, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.dpfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.dpsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.col, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.cell, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.inner, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.fixup, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.gathsol, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.bt, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.btfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.btsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.btcell, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.corerej, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16s.nrej, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		const SSEMetrics &dpSse8s = total ? dpSse8Seed : dpSse8uSeed;
		itoa10<uint64_t>(dpSse8s.dp, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.dpsat, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.dpfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.dpsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.col, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.cell, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.inner, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.fixup, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.gathsol, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.bt, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.btfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.btsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.btcell, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.corerej, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8s.nrej, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		const SSEMetrics &dpSse16m = total ? dpSse16Mate : dpSse16uMate;
		itoa10<uint64_t>(dpSse16m.dp, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.dpsat, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.dpfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.dpsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.col, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.cell, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.inner, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.fixup, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.gathsol, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.bt, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.btfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.btsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.btcell, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.corerej, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse16m.nrej, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		const SSEMetrics &dpSse8m = total ? dpSse8Mate : dpSse8uMate;
		itoa10<uint64_t>(dpSse8m.dp, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.dpsat, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.dpfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.dpsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.col, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.cell, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.inner, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.fixup, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.gathsol, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.bt, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.btfail, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.btsucc, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.btcell, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.corerej, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(dpSse8m.nrej, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? nbtfiltst : nbtfiltst_u, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? nbtfiltsc : nbtfiltsc_u, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<uint64_t>(total ? nbtfiltdo : nbtfiltdo_u, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
#ifdef USE_MEM_TALLY
		itoa10<size_t>(gMemTally.peak() >> 20, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<size_t>(gMemTally.peak(0) >> 20, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<size_t>(gMemTally.peak(EBWT_CAT) >> 20, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<size_t>(gMemTally.peak(CA_CAT) >> 20, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<size_t>(gMemTally.peak(GW_CAT) >> 20, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<size_t>(gMemTally.peak(AL_CAT) >> 20, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<size_t>(gMemTally.peak(DP_CAT) >> 20, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<size_t>(gMemTally.peak(MISC_CAT) >> 20, buf);
		if (metricsStderr)
			stderrSs << buf << '\t';
		if (o != NULL)
		{
			o->writeChars(buf);
			o->write('\t');
		}
		itoa10<size_t>(gMemTally.peak(DEBUG_CAT) >> 20, buf);
		if (metricsStderr)
			stderrSs << buf;
		if (o != NULL)
		{
			o->writeChars(buf);
		}
#endif
		if (o != NULL)
		{
			o->write('\n');
		}
		if (metricsStderr)
			cerr << stderrSs.str().c_str() << endl;
		if (!total)
			mergeIncrementals();
	}
	void mergeIncrementals()
	{
		olm.merge(olmu);
		sdm.merge(sdmu);
		wlm.merge(wlmu);
		swmSeed.merge(swmuSeed);
		swmMate.merge(swmuMate);
		dpSse8Seed.merge(dpSse8uSeed);
		dpSse8Mate.merge(dpSse8uMate);
		dpSse16Seed.merge(dpSse16uSeed);
		dpSse16Mate.merge(dpSse16uMate);
		nbtfiltst_u += nbtfiltst;
		nbtfiltsc_u += nbtfiltsc;
		nbtfiltdo_u += nbtfiltdo;
		olmu.reset();
		sdmu.reset();
		wlmu.reset();
		swmuSeed.reset();
		swmuMate.reset();
		rpmu.reset();
		dpSse8uSeed.reset();
		dpSse8uMate.reset();
		dpSse16uSeed.reset();
		dpSse16uMate.reset();
		nbtfiltst_u = 0;
		nbtfiltsc_u = 0;
		nbtfiltdo_u = 0;
	}
	OuterLoopMetrics olm;
	SeedSearchMetrics sdm;
	WalkMetrics wlm;
	SwMetrics swmSeed;
	SwMetrics swmMate;
	ReportingMetrics rpm;
	SSEMetrics dpSse8Seed;
	SSEMetrics dpSse8Mate;
	SSEMetrics dpSse16Seed;
	SSEMetrics dpSse16Mate;
	uint64_t nbtfiltst;
	uint64_t nbtfiltsc;
	uint64_t nbtfiltdo;
	OuterLoopMetrics olmu;
	SeedSearchMetrics sdmu;
	WalkMetrics wlmu;
	SwMetrics swmuSeed;
	SwMetrics swmuMate;
	ReportingMetrics rpmu;
	SSEMetrics dpSse8uSeed;
	SSEMetrics dpSse8uMate;
	SSEMetrics dpSse16uSeed;
	SSEMetrics dpSse16uMate;
	uint64_t nbtfiltst_u;
	uint64_t nbtfiltsc_u;
	uint64_t nbtfiltdo_u;
	MUTEX_T mutex_m;
	bool first;
	time_t lastElapsed;
};
static PerfMetrics metrics;
#define ROTL(n, x) (((x) << (n)) | ((x) >> (32 - n)))
#define ROTR(n, x) (((x) >> (n)) | ((x) << (32 - n)))
static inline void printMmsSkipMsg(
	const PatternSourcePerThread &ps,
	bool paired,
	bool mate1,
	int seedmms)
{
	ostringstream os;
	if (paired)
	{
		os << "Warning: skipping mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "' because length (" << (mate1 ? ps.read_a().patFw.length() : ps.read_b().patFw.length())
		   << ") <= # seed mismatches (" << seedmms << ")" << endl;
	}
	else
	{
		os << "Warning: skipping read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "' because length (" << (mate1 ? ps.read_a().patFw.length() : ps.read_b().patFw.length())
		   << ") <= # seed mismatches (" << seedmms << ")" << endl;
	}
	cerr << os.str().c_str();
}
static inline void printLenSkipMsg(
	const PatternSourcePerThread &ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if (paired)
	{
		os << "Warning: skipping mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "' because it was < 2 characters long" << endl;
	}
	else
	{
		os << "Warning: skipping read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "' because it was < 2 characters long" << endl;
	}
	cerr << os.str().c_str();
}
static inline void printLocalScoreMsg(
	const PatternSourcePerThread &ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if (paired)
	{
		os << "Warning: minimum score function gave negative number in "
		   << "--local mode for mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "; setting to 0 instead" << endl;
	}
	else
	{
		os << "Warning: minimum score function gave negative number in "
		   << "--local mode for read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "; setting to 0 instead" << endl;
	}
	cerr << os.str().c_str();
}
static inline void printEEScoreMsg(
	const PatternSourcePerThread &ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if (paired)
	{
		os << "Warning: minimum score function gave positive number in "
		   << "--end-to-end mode for mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "; setting to 0 instead" << endl;
	}
	else
	{
		os << "Warning: minimum score function gave positive number in "
		   << "--end-to-end mode for read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "; setting to 0 instead" << endl;
	}
	cerr << os.str().c_str();
}
static void setupMinScores(
	const PatternSourcePerThread &ps,
	bool paired,
	bool localAlign,
	const Scoring &sc,
	const size_t *rdlens,
	TAlScore *minsc,
	TAlScore *maxpen)
{
	if (bwaSwLike)
	{
		float a = (float)sc.match(30);
		float T = bwaSwLikeT, c = bwaSwLikeC;
		minsc[0] = (TAlScore)max<float>(a * T, a * c * log(rdlens[0]));
		if (paired)
		{
			minsc[1] = (TAlScore)max<float>(a * T, a * c * log(rdlens[1]));
		}
	}
	else
	{
		minsc[0] = scoreMin.f<TAlScore>(rdlens[0]);
		if (paired)
			minsc[1] = scoreMin.f<TAlScore>(rdlens[1]);
		if (localAlign)
		{
			if (minsc[0] < 0)
			{
				if (!gQuiet)
					printLocalScoreMsg(ps, paired, true);
				minsc[0] = 0;
			}
			if (paired && minsc[1] < 0)
			{
				if (!gQuiet)
					printLocalScoreMsg(ps, paired, false);
				minsc[1] = 0;
			}
		}
		else
		{
			if (minsc[0] > 0)
			{
				if (!gQuiet)
					printEEScoreMsg(ps, paired, true);
				minsc[0] = 0;
			}
			if (paired && minsc[1] > 0)
			{
				if (!gQuiet)
					printEEScoreMsg(ps, paired, false);
				minsc[1] = 0;
			}
		}
	}
	if (localAlign)
	{
		TAlScore perfect0 = sc.perfectScore(rdlens[0]);
		assert_geq(perfect0, minsc[0]);
		maxpen[0] = perfect0 - minsc[0];
		if (paired)
		{
			TAlScore perfect1 = sc.perfectScore(rdlens[1]);
			assert_geq(perfect1, minsc[1]);
			maxpen[1] = perfect1 - minsc[1];
		}
		else
		{
			maxpen[1] = std::numeric_limits<TAlScore>::min();
		}
	}
	else
	{
		assert_leq(minsc[0], 0);
		maxpen[0] = -minsc[0];
		if (paired)
		{
			assert_leq(minsc[1], 0);
			maxpen[1] = -minsc[1];
		}
		else
		{
			maxpen[1] = std::numeric_limits<TAlScore>::min();
		}
	}
}
#define MERGE_METRICS(met)       \
	{                            \
		msink.mergeMetrics(rpm); \
		met.merge(               \
			&olm,                \
			&sdm,                \
			&wlm,                \
			&swmSeed,            \
			&swmMate,            \
			&rpm,                \
			&sseU8ExtendMet,     \
			&sseU8MateMet,       \
			&sseI16ExtendMet,    \
			&sseI16MateMet,      \
			nbtfiltst,           \
			nbtfiltsc,           \
			nbtfiltdo);          \
		olm.reset();             \
		sdm.reset();             \
		wlm.reset();             \
		swmSeed.reset();         \
		swmMate.reset();         \
		rpm.reset();             \
		sseU8ExtendMet.reset();  \
		sseU8MateMet.reset();    \
		sseI16ExtendMet.reset(); \
		sseI16MateMet.reset();   \
	}
#define MERGE_SW(x)          \
	{                        \
		x.merge(             \
			sseU8ExtendMet,  \
			sseU8MateMet,    \
			sseI16ExtendMet, \
			sseI16MateMet,   \
			nbtfiltst,       \
			nbtfiltsc,       \
			nbtfiltdo);      \
		x.resetCounters();   \
	}
#ifdef PER_THREAD_TIMING
void get_cpu_and_node(int &cpu, int &node)
{
	unsigned long a, d, c;
	__asm__ volatile("rdtscp"
					 : "=a"(a), "=d"(d), "=c"(c));
	node = (c & 0xFFF000) >> 12;
	cpu = c & 0xFFF;
}
#endif
class ThreadCounter
{
public:
	ThreadCounter()
	{
#ifdef WITH_TBB
		thread_counter.fetch_add(1);
#else
		ThreadSafe ts(thread_counter_mutex);
		thread_counter++;
#endif
	}
	~ThreadCounter()
	{
#ifdef WITH_TBB
		thread_counter.fetch_sub(1);
#else
		ThreadSafe ts(thread_counter_mutex);
		thread_counter--;
#endif
	}
};
void make_code(void)
{
	for (size_t i = 0; i < 256; i++)
		code[i] = 4;
	code['A'] = code['a'] = 0;
	code['C'] = code['c'] = 1;
	code['G'] = code['g'] = 2;
	code['T'] = code['t'] = 3;
	code['N'] = code['n'] = 0;
}
#ifdef WITH_TBB
static void multiseedSearchWorker(void *vp)
{
	thread_tracking_pair *p = (thread_tracking_pair *)vp;
	int tid = p->tid;
#else
static void multiseedSearchWorker(void *vp)
{
	int tid = *((int *)vp);
#endif
	assert(multiseed_ebwtFw != NULL);
	assert(multiseedMms == 0 || multiseed_ebwtBw != NULL);
	PatternComposer &patsrc = *multiseed_patsrc;
	PatternParams pp = multiseed_pp;
	const Ebwt &ebwtFw = *multiseed_ebwtFw;
	const Ebwt *ebwtBw = multiseed_ebwtBw;
	const Scoring &sc = *multiseed_sc;
	const BitPairReference &ref = *multiseed_refs;
	AlnSink &msink = *multiseed_msink;
	OutFileBuf *metricsOfb = multiseed_metricsOfb;
#ifdef TIME_STATS
	auto start_multiseedSearchWorker = std::chrono::system_clock::now();
#endif
	{
		make_code();
#ifdef PER_THREAD_TIMING
		uint64_t ncpu_changeovers = 0;
		uint64_t nnuma_changeovers = 0;
		int current_cpu = 0, current_node = 0;
		get_cpu_and_node(current_cpu, current_node);
		std::stringstream ss;
		std::string msg;
		ss << "thread: " << tid << " time: ";
		msg = ss.str();
		Timer timer(std::cout, msg.c_str());
#endif
		unique_ptr<PatternSourcePerThreadFactory> patsrcFact(createPatsrcFactory(patsrc, pp, tid));
		unique_ptr<PatternSourcePerThread> ps(patsrcFact->create());
		PtrWrap<AlignmentCache> scLocal;
		if (!msNoCache)
		{
			scLocal.init(new AlignmentCache(seedCacheLocalMB * 1024 * 1024, false));
		}
		AlignmentCache scCurrent(seedCacheCurrentMB * 1024 * 1024, false);
		AlignmentCacheIface ca(
			&scCurrent,
			scLocal.get(),
			msNoCache ? NULL : multiseed_ca);
		ReportingParams rp(
			(allHits ? std::numeric_limits<THitInt>::max() : khits),
			mhits,
			0,
			msample,
			gReportDiscordant,
			gReportMixed);
		unique_ptr<Mapq> bmapq(new_mapq(mapqv, scoreMin, sc));
		AlnSinkWrap msinkwrap(
			msink,
			rp,
			*bmapq,
			(size_t)tid);
		ofstream *dpLog = NULL, *dpLogOpp = NULL;
		if (!logDps.empty())
		{
			dpLog = new ofstream(logDps.c_str(), ofstream::out);
			dpLog->sync_with_stdio(false);
		}
		if (!logDpsOpp.empty())
		{
			dpLogOpp = new ofstream(logDpsOpp.c_str(), ofstream::out);
			dpLogOpp->sync_with_stdio(false);
		}
		SeedAligner al;
		SwDriver sd(exactCacheCurrentMB * 1024 * 1024);
		SwAligner sw(dpLog), osw(dpLogOpp);
		SeedResults shs[2];
		OuterLoopMetrics olm;
		SeedSearchMetrics sdm;
		WalkMetrics wlm;
		SwMetrics swmSeed, swmMate;
		ReportingMetrics rpm;
		RandomSource rnd, rndArb;
		SSEMetrics sseU8ExtendMet;
		SSEMetrics sseU8MateMet;
		SSEMetrics sseI16ExtendMet;
		SSEMetrics sseI16MateMet;
		uint64_t nbtfiltst = 0;
		uint64_t nbtfiltsc = 0;
		uint64_t nbtfiltdo = 0;
		ASSERT_ONLY(BTDnaString tmp);
		int pepolFlag;
		if (gMate1fw && gMate2fw)
		{
			pepolFlag = PE_POLICY_FF;
		}
		else if (gMate1fw && !gMate2fw)
		{
			pepolFlag = PE_POLICY_FR;
		}
		else if (!gMate1fw && gMate2fw)
		{
			pepolFlag = PE_POLICY_RF;
		}
		else
		{
			pepolFlag = PE_POLICY_RR;
		}
		assert_geq(gMaxInsert, gMinInsert);
		assert_geq(gMinInsert, 0);
		PairedEndPolicy pepol(
			pepolFlag,
			gMaxInsert,
			gMinInsert,
			localAlign,
			gFlippedMatesOK,
			gDovetailMatesOK,
			gContainMatesOK,
			gOlapMatesOK,
			gExpandToFrag);
		PerfMetrics metricsPt;
		BTString nametmp;
		EList<Seed> seeds1, seeds2;
		EList<Seed> *seeds[2] = {&seeds1, &seeds2};
		PerReadMetrics prm;
		time_t iTime = time(0);
		bool exhaustive[2] = {false, false};
		bool filt[2] = {true, true};
		bool nfilt[2] = {true, true};
		bool scfilt[2] = {true, true};
		bool lenfilt[2] = {true, true};
		bool qcfilt[2] = {true, true};
		rndArb.init((uint32_t)time(0));
		int mergei = 0;
		int mergeival = 16;
		bool done = false;
		while (!done)
		{
			pair<bool, bool> ret = ps->nextReadPair();
			bool success = ret.first;
			done = ret.second;
			if (!success && done)
			{
				break;
			}
			else if (!success)
			{
				continue;
			}
			TReadId rdid = ps->read_a().rdid;
			bool sample = true;
			if (arbitraryRandom)
			{
				ps->read_a().seed = rndArb.nextU32();
				ps->read_b().seed = rndArb.nextU32();
			}
			if (sampleFrac < 1.0f)
			{
				rnd.init(ROTL(ps->read_a().seed, 2));
				sample = rnd.nextFloat() < sampleFrac;
			}
			if (rdid >= skipReads && rdid < qUpto && sample)
			{
				bool retry = true;
				if (metricsIval > 0 &&
					(metricsOfb != NULL || metricsStderr) &&
					!metricsPerRead &&
					++mergei == mergeival)
				{
					MERGE_METRICS(metrics);
					mergei = 0;
					if (tid == 0)
					{
						time_t curTime = time(0);
						if (curTime - iTime >= metricsIval)
						{
							metrics.reportInterval(metricsOfb, metricsStderr, false, NULL);
							iTime = curTime;
						}
					}
				}
				prm.reset();
				prm.doFmString = false;
				if (sam_print_xt)
				{
					gettimeofday(&prm.tv_beg, &prm.tz_beg);
				}
#ifdef PER_THREAD_TIMING
				int cpu = 0, node = 0;
				get_cpu_and_node(cpu, node);
				if (cpu != current_cpu)
				{
					ncpu_changeovers++;
					current_cpu = cpu;
				}
				if (node != current_node)
				{
					nnuma_changeovers++;
					current_node = node;
				}
#endif
				while (retry)
				{
					retry = false;
					ca.nextRead();
					olm.reads++;
					assert(!ca.aligning());
					bool paired = !ps->read_b().empty();
					const size_t rdlen1 = ps->read_a().length();
					const size_t rdlen2 = paired ? ps->read_b().length() : 0;
					olm.bases += (rdlen1 + rdlen2);
					msinkwrap.nextRead(
						&ps->read_a(),
						paired ? &ps->read_b() : NULL,
						rdid,
						sc.qualitiesMatter());
					assert(msinkwrap.inited());
					size_t rdlens[2] = {rdlen1, rdlen2};
					size_t rdrows[2] = {rdlen1, rdlen2};
					TAlScore minsc[2];
					minsc[0] = minsc[1] = std::numeric_limits<TAlScore>::max();
					if (bwaSwLike)
					{
						float a = (float)sc.match(30);
						float T = bwaSwLikeT, c = bwaSwLikeC;
						minsc[0] = (TAlScore)max<float>(a * T, a * c * log(rdlens[0]));
						if (paired)
						{
							minsc[1] = (TAlScore)max<float>(a * T, a * c * log(rdlens[1]));
						}
					}
					else
					{
						minsc[0] = scoreMin.f<TAlScore>(rdlens[0]);
						if (paired)
							minsc[1] = scoreMin.f<TAlScore>(rdlens[1]);
						if (localAlign)
						{
							if (minsc[0] < 0)
							{
								if (!gQuiet)
									printLocalScoreMsg(*ps, paired, true);
								minsc[0] = 0;
							}
							if (paired && minsc[1] < 0)
							{
								if (!gQuiet)
									printLocalScoreMsg(*ps, paired, false);
								minsc[1] = 0;
							}
						}
						else
						{
							if (minsc[0] > 0)
							{
								if (!gQuiet)
									printEEScoreMsg(*ps, paired, true);
								minsc[0] = 0;
							}
							if (paired && minsc[1] > 0)
							{
								if (!gQuiet)
									printEEScoreMsg(*ps, paired, false);
								minsc[1] = 0;
							}
						}
					}
					size_t readns[2] = {0, 0};
					sc.nFilterPair(
						&ps->read_a().patFw,
						paired ? &ps->read_b().patFw : NULL,
						readns[0],
						readns[1],
						nfilt[0],
						nfilt[1]);
					scfilt[0] = sc.scoreFilter(minsc[0], rdlens[0]);
					scfilt[1] = sc.scoreFilter(minsc[1], rdlens[1]);
					lenfilt[0] = lenfilt[1] = true;
					if (rdlens[0] <= (size_t)multiseedMms || rdlens[0] < 2)
					{
						if (!gQuiet)
							printMmsSkipMsg(*ps, paired, true, multiseedMms);
						lenfilt[0] = false;
					}
					if ((rdlens[1] <= (size_t)multiseedMms || rdlens[1] < 2) && paired)
					{
						if (!gQuiet)
							printMmsSkipMsg(*ps, paired, false, multiseedMms);
						lenfilt[1] = false;
					}
					if (rdlens[0] < 2)
					{
						if (!gQuiet)
							printLenSkipMsg(*ps, paired, true);
						lenfilt[0] = false;
					}
					if (rdlens[1] < 2 && paired)
					{
						if (!gQuiet)
							printLenSkipMsg(*ps, paired, false);
						lenfilt[1] = false;
					}
					qcfilt[0] = qcfilt[1] = true;
					if (qcFilter)
					{
						qcfilt[0] = (ps->read_a().filter != '0');
						qcfilt[1] = (ps->read_b().filter != '0');
					}
					filt[0] = (nfilt[0] && scfilt[0] && lenfilt[0] && qcfilt[0]);
					filt[1] = (nfilt[1] && scfilt[1] && lenfilt[1] && qcfilt[1]);
					prm.nFilt += (filt[0] ? 0 : 1) + (filt[1] ? 0 : 1);
					Read *rds[2] = {&ps->read_a(), &ps->read_b()};
					assert(msinkwrap.empty());
					sd.nextRead(paired, rdrows[0], rdrows[1]);
					size_t minedfw[2] = {0, 0};
					size_t minedrc[2] = {0, 0};
					bool nofw[2] = {false, false};
					bool norc[2] = {false, false};
					nofw[0] = paired ? (gMate1fw ? gNofw : gNorc) : gNofw;
					norc[0] = paired ? (gMate1fw ? gNorc : gNofw) : gNorc;
					nofw[1] = paired ? (gMate2fw ? gNofw : gNorc) : gNofw;
					norc[1] = paired ? (gMate2fw ? gNorc : gNofw) : gNorc;
					int nceil[2] = {0, 0};
					nceil[0] = nCeil.f<int>((double)rdlens[0]);
					nceil[0] = min(nceil[0], (int)rdlens[0]);
					if (paired)
					{
						nceil[1] = nCeil.f<int>((double)rdlens[1]);
						nceil[1] = min(nceil[1], (int)rdlens[1]);
					}
					exhaustive[0] = exhaustive[1] = false;
					size_t matemap[2] = {0, 1};
					bool pairPostFilt = filt[0] && filt[1];
					if (pairPostFilt)
					{
						rnd.init(ps->read_a().seed ^ ps->read_b().seed);
					}
					else
					{
						rnd.init(ps->read_a().seed);
					}
					int interval[2] = {0, 0};
					for (size_t mate = 0; mate < (paired ? 2 : 1); mate++)
					{
						interval[mate] = msIval.f<int>((double)rdlens[mate]);
						if (filt[0] && filt[1])
						{
							interval[mate] = (int)(interval[mate] * 1.2 + 0.5);
						}
						interval[mate] = max(interval[mate], 1);
					}
					size_t streak[2] = {maxDpStreak, maxDpStreak};
					size_t mtStreak[2] = {maxMateStreak, maxMateStreak};
					size_t mxDp[2] = {maxDp, maxDp};
					size_t mxUg[2] = {maxUg, maxUg};
					size_t mxIter[2] = {maxIters, maxIters};
					if (allHits)
					{
						streak[0] = streak[1] = std::numeric_limits<size_t>::max();
						mtStreak[0] = mtStreak[1] = std::numeric_limits<size_t>::max();
						mxDp[0] = mxDp[1] = std::numeric_limits<size_t>::max();
						mxUg[0] = mxUg[1] = std::numeric_limits<size_t>::max();
						mxIter[0] = mxIter[1] = std::numeric_limits<size_t>::max();
					}
					else if (khits > 1)
					{
						for (size_t mate = 0; mate < 2; mate++)
						{
							streak[mate] += (khits - 1) * maxStreakIncr;
							mtStreak[mate] += (khits - 1) * maxStreakIncr;
							mxDp[mate] += (khits - 1) * maxItersIncr;
							mxUg[mate] += (khits - 1) * maxItersIncr;
							mxIter[mate] += (khits - 1) * maxItersIncr;
						}
					}
					if (filt[0] && filt[1])
					{
						streak[0] = (size_t)ceil((double)streak[0] / 2.0);
						streak[1] = (size_t)ceil((double)streak[1] / 2.0);
						assert_gt(streak[1], 0);
					}
					assert_gt(streak[0], 0);
					size_t nrounds[2] = {nSeedRounds, nSeedRounds};
					if (filt[0] && filt[1])
					{
						nrounds[0] = (size_t)ceil((double)nrounds[0] / 2.0);
						nrounds[1] = (size_t)ceil((double)nrounds[1] / 2.0);
						assert_gt(nrounds[1], 0);
					}
					assert_gt(nrounds[0], 0);
					for (size_t mate = 0; mate < (paired ? 2 : 1); mate++)
					{
						if (!filt[mate])
						{
							olm.freads++;
							olm.fbases += rdlens[mate];
						}
						else
						{
							shs[mate].clear();
							shs[mate].nextRead(mate == 0 ? ps->read_a() : ps->read_b());
							assert(shs[mate].empty());
							olm.ureads++;
							olm.ubases += rdlens[mate];
						}
					}
					size_t eePeEeltLimit = std::numeric_limits<size_t>::max();
					bool done[2] = {!filt[0], !filt[1]};
					size_t nelt[2] = {0, 0};
#ifdef TIME_STATS
					auto start_exact = std::chrono::system_clock::now();
#endif
					if (doExactUpFront)
					{
						for (size_t matei = 0; matei < (paired ? 2 : 1); matei++)
						{
							size_t mate = matemap[matei];
							if (!filt[mate] || done[mate] || msinkwrap.state().doneWithMate(mate == 0))
							{
								continue;
							}
							swmSeed.exatts++;
#ifdef TIME_STATS
							exact_Nums++;
#endif
							nelt[mate] = al.exactSweep(
								ebwtFw,
								*rds[mate],
								sc,
								nofw[mate],
								norc[mate],
								2,
								minedfw[mate],
								minedrc[mate],
								true,
								shs[mate],
								sdm);
							size_t bestmin = min(minedfw[mate], minedrc[mate]);
							if (bestmin == 0)
							{
								sdm.bestmin0++;
							}
							else if (bestmin == 1)
							{
								sdm.bestmin1++;
							}
							else
							{
								assert_eq(2, bestmin);
								sdm.bestmin2++;
							}
						}
						matemap[0] = 0;
						matemap[1] = 1;
						if (nelt[0] > 0 && nelt[1] > 0 && nelt[0] > nelt[1])
						{
							matemap[0] = 1;
							matemap[1] = 0;
						}
						for (size_t matei = 0; matei < (seedSumm ? 0 : 2); matei++)
						{
							size_t mate = matemap[matei];
							if (nelt[mate] == 0 || nelt[mate] > eePeEeltLimit)
							{
								shs[mate].clearExactE2eHits();
								continue;
							}
							if (msinkwrap.state().doneWithMate(mate == 0))
							{
								shs[mate].clearExactE2eHits();
								done[mate] = true;
								continue;
							}
							assert(filt[mate]);
							assert(matei == 0 || paired);
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							int ret = 0;
							if (paired)
							{
								ret = sd.extendSeedsPaired(
									*rds[mate],
									*rds[mate ^ 1],
									mate == 0,
									!filt[mate ^ 1],
									shs[mate],
									ebwtFw,
									ebwtBw,
									ref,
									sw,
									osw,
									sc,
									pepol,
									-1,
									0,
									0,
									minsc[mate],
									minsc[mate ^ 1],
									nceil[mate],
									nceil[mate ^ 1],
									nofw[mate],
									norc[mate],
									maxhalf,
									doUngapped,
									mxIter[mate],
									mxUg[mate],
									mxDp[mate],
									streak[mate],
									streak[mate],
									streak[mate],
									mtStreak[mate],
									doExtend,
									enable8,
									cminlen,
									cpow2,
									doTri,
									tighten,
									ca,
									rnd,
									wlm,
									swmSeed,
									swmMate,
									prm,
									&msinkwrap,
									true,
									true,
									gReportDiscordant,
									gReportMixed,
									exhaustive[mate]);
							}
							else
							{
								ret = sd.extendSeeds(
									*rds[mate],
									mate == 0,
									shs[mate],
									ebwtFw,
									ebwtBw,
									ref,
									sw,
									sc,
									-1,
									0,
									0,
									minsc[mate],
									nceil[mate],
									maxhalf,
									doUngapped,
									mxIter[mate],
									mxUg[mate],
									mxDp[mate],
									streak[mate],
									streak[mate],
									doExtend,
									enable8,
									cminlen,
									cpow2,
									doTri,
									tighten,
									ca,
									rnd,
									wlm,
									swmSeed,
									prm,
									&msinkwrap,
									true,
									exhaustive[mate]);
							}
							assert_gt(ret, 0);
							MERGE_SW(sw);
							MERGE_SW(osw);
							shs[mate].clearExactE2eHits();
							if (ret == EXTEND_EXHAUSTED_CANDIDATES)
							{
							}
							else if (ret == EXTEND_POLICY_FULFILLED)
							{
								if (msinkwrap.state().doneWithMate(mate == 0))
								{
									done[mate] = true;
								}
								if (msinkwrap.state().doneWithMate(mate == 1))
								{
									done[mate ^ 1] = true;
								}
							}
							else if (ret == EXTEND_PERFECT_SCORE)
							{
								done[mate] = true;
							}
							else if (ret == EXTEND_EXCEEDED_HARD_LIMIT)
							{
								done[mate] = true;
							}
							else if (ret == EXTEND_EXCEEDED_SOFT_LIMIT)
							{
							}
							else
							{
								cerr << "Bad return value: " << ret << endl;
								throw 1;
							}
							if (!done[mate])
							{
								TAlScore perfectScore = sc.perfectScore(rdlens[mate]);
								if (!done[mate] && minsc[mate] == perfectScore)
								{
									done[mate] = true;
								}
							}
						}
					}
#ifdef TIME_STATS
					auto end_exact = std::chrono::system_clock::now();
					auto elapsed_exact = std::chrono::duration_cast<std::chrono::microseconds>(end_exact - start_exact);
					time_exact += elapsed_exact.count();
#endif
#ifdef TIME_STATS
					auto start_1mm = std::chrono::system_clock::now();
#endif
					if (do1mmUpFront && !seedSumm)
					{
						for (size_t matei = 0; matei < (paired ? 2 : 1); matei++)
						{
							size_t mate = matemap[matei];
							if (!filt[mate] || done[mate] || nelt[mate] > eePeEeltLimit)
							{
								shs[mate].clear1mmE2eHits();
								nelt[mate] = 0;
								continue;
							}
							nelt[mate] = 0;
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							assert(shs[mate].empty());
							assert(shs[mate].repOk(&ca.current()));
							bool yfw = minedfw[mate] <= 1 && !nofw[mate];
							bool yrc = minedrc[mate] <= 1 && !norc[mate];
							if (yfw || yrc)
							{
								swmSeed.mm1atts++;
#ifdef TIME_STATS
								onemm_Nums++;
#endif
								al.oneMmSearch(
									&ebwtFw,
									ebwtBw,
									*rds[mate],
									sc,
									minsc[mate],
									!yfw,
									!yrc,
									localAlign,
									false,
									true,
									shs[mate],
									sdm);
								nelt[mate] = shs[mate].num1mmE2eHits();
							}
						}
						matemap[0] = 0;
						matemap[1] = 1;
						if (nelt[0] > 0 && nelt[1] > 0 && nelt[0] > nelt[1])
						{
							matemap[0] = 1;
							matemap[1] = 0;
						}
						for (size_t matei = 0; matei < (seedSumm ? 0 : 2); matei++)
						{
							size_t mate = matemap[matei];
							if (nelt[mate] == 0 || nelt[mate] > eePeEeltLimit)
							{
								continue;
							}
							if (msinkwrap.state().doneWithMate(mate == 0))
							{
								done[mate] = true;
								continue;
							}
							int ret = 0;
							if (paired)
							{
								ret = sd.extendSeedsPaired(
									*rds[mate],
									*rds[mate ^ 1],
									mate == 0,
									!filt[mate ^ 1],
									shs[mate],
									ebwtFw,
									ebwtBw,
									ref,
									sw,
									osw,
									sc,
									pepol,
									-1,
									0,
									0,
									minsc[mate],
									minsc[mate ^ 1],
									nceil[mate],
									nceil[mate ^ 1],
									nofw[mate],
									norc[mate],
									maxhalf,
									doUngapped,
									mxIter[mate],
									mxUg[mate],
									mxDp[mate],
									streak[mate],
									streak[mate],
									streak[mate],
									mtStreak[mate],
									doExtend,
									enable8,
									cminlen,
									cpow2,
									doTri,
									tighten,
									ca,
									rnd,
									wlm,
									swmSeed,
									swmMate,
									prm,
									&msinkwrap,
									true,
									true,
									gReportDiscordant,
									gReportMixed,
									exhaustive[mate]);
							}
							else
							{
								ret = sd.extendSeeds(
									*rds[mate],
									mate == 0,
									shs[mate],
									ebwtFw,
									ebwtBw,
									ref,
									sw,
									sc,
									-1,
									0,
									0,
									minsc[mate],
									nceil[mate],
									maxhalf,
									doUngapped,
									mxIter[mate],
									mxUg[mate],
									mxDp[mate],
									streak[mate],
									streak[mate],
									doExtend,
									enable8,
									cminlen,
									cpow2,
									doTri,
									tighten,
									ca,
									rnd,
									wlm,
									swmSeed,
									prm,
									&msinkwrap,
									true,
									exhaustive[mate]);
							}
							assert_gt(ret, 0);
							MERGE_SW(sw);
							MERGE_SW(osw);
							shs[mate].clear1mmE2eHits();
							if (ret == EXTEND_EXHAUSTED_CANDIDATES)
							{
							}
							else if (ret == EXTEND_POLICY_FULFILLED)
							{
								if (msinkwrap.state().doneWithMate(mate == 0))
								{
									done[mate] = true;
								}
								if (msinkwrap.state().doneWithMate(mate == 1))
								{
									done[mate ^ 1] = true;
								}
							}
							else if (ret == EXTEND_PERFECT_SCORE)
							{
								done[mate] = true;
							}
							else if (ret == EXTEND_EXCEEDED_HARD_LIMIT)
							{
								done[mate] = true;
							}
							else if (ret == EXTEND_EXCEEDED_SOFT_LIMIT)
							{
							}
							else
							{
								cerr << "Bad return value: " << ret << endl;
								throw 1;
							}
							if (!done[mate])
							{
								TAlScore perfectScore = sc.perfectScore(rdlens[mate]);
								if (!done[mate] && minsc[mate] == perfectScore)
								{
									done[mate] = true;
								}
							}
						}
					}
#ifdef TIME_STATS
					auto end_1mm = std::chrono::system_clock::now();
					auto elapsed_1mm = std::chrono::duration_cast<std::chrono::microseconds>(end_1mm - start_1mm);
					time_1mm += elapsed_1mm.count();
#endif
					int seedlens[2] = {multiseedLen, multiseedLen};
					nrounds[0] = min<size_t>(nrounds[0], interval[0]);
					nrounds[1] = min<size_t>(nrounds[1], interval[1]);
					Constraint gc = Constraint::penaltyFuncBased(scoreMin);
					size_t seedsTried = 0;
					size_t seedsTriedMS[] = {0, 0, 0, 0};
					size_t nUniqueSeeds = 0, nRepeatSeeds = 0, seedHitTot = 0;
					size_t nUniqueSeedsMS[] = {0, 0, 0, 0};
					size_t nRepeatSeedsMS[] = {0, 0, 0, 0};
					size_t seedHitTotMS[] = {0, 0, 0, 0};
#ifdef TIME_STATS
					auto start_seed = std::chrono::system_clock::now();
#endif
					for (size_t roundi = 0; roundi < nSeedRounds; roundi++)
					{
						ca.nextRead();
						shs[0].clearSeeds();
						shs[1].clearSeeds();
						assert(shs[0].empty());
						assert(shs[1].empty());
						assert(shs[0].repOk(&ca.current()));
						assert(shs[1].repOk(&ca.current()));
						for (size_t matei = 0; matei < (paired ? 2 : 1); matei++)
						{
							size_t mate = matemap[matei];
							if (done[mate] || msinkwrap.state().doneWithMate(mate == 0))
							{
								done[mate] = true;
								continue;
							}
							if (roundi >= nrounds[mate])
							{
								continue;
							}
							if (interval[mate] <= (int)roundi)
							{
								continue;
							}
							size_t offset = (interval[mate] * roundi) / nrounds[mate];
							assert(roundi == 0 || offset > 0);
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							assert(shs[mate].repOk(&ca.current()));
							swmSeed.sdatts++;
							seeds[mate]->clear();
#ifdef TIME_STATS
							auto start_seed_mmSeeds = std::chrono::system_clock::now();
#endif
							Seed::mmSeeds(
								multiseedMms,
								seedlens[mate],
								*seeds[mate],
								gc);
#ifdef TIME_STATS
							auto end_seed_mmSeeds = std::chrono::system_clock::now();
							auto elapsed_seed_mmSeeds = std::chrono::duration_cast<std::chrono::microseconds>(end_seed_mmSeeds - start_seed_mmSeeds);
							time_seed_mmSeeds += elapsed_seed_mmSeeds.count();
#endif
							if (offset > 0 && (*seeds[mate])[0].len + offset > rds[mate]->length())
							{
								continue;
							}
							std::pair<int, int> instFw, instRc;
#ifdef TIME_STATS
							auto start_seed_instantiateSeeds = std::chrono::system_clock::now();
#endif
							std::pair<int, int> inst = al.instantiateSeeds(
								*seeds[mate],
								offset,
								interval[mate],
								*rds[mate],
								sc,
								nofw[mate],
								norc[mate],
								ca,
								shs[mate],
								sdm,
								instFw,
								instRc);
#ifdef TIME_STATS
							auto end_seed_instantiateSeeds = std::chrono::system_clock::now();
							auto elapsed_seed_instantiateSeeds = std::chrono::duration_cast<std::chrono::microseconds>(end_seed_instantiateSeeds - start_seed_instantiateSeeds);
							time_seed_instantiateSeeds += elapsed_seed_instantiateSeeds.count();
#endif
							assert(shs[mate].repOk(&ca.current()));
							if (inst.first + inst.second == 0)
							{
								assert(shs[mate].empty());
								done[mate] = true;
								break;
							}
							seedsTried += (inst.first + inst.second);
							seedsTriedMS[mate * 2 + 0] = instFw.first + instFw.second;
							seedsTriedMS[mate * 2 + 1] = instRc.first + instRc.second;
#ifdef TIME_STATS
							auto start_seed_searchAllSeeds = std::chrono::system_clock::now();
#endif
							al.searchAllSeeds(
								*seeds[mate],
								&ebwtFw,
								ebwtBw,
								*rds[mate],
								sc,
								ca,
								shs[mate],
								sdm,
								prm);
#ifdef TIME_STATS
							auto end_seed_searchAllSeeds = std::chrono::system_clock::now();
							auto elapsed_seed_searchAllSeeds = std::chrono::duration_cast<std::chrono::microseconds>(end_seed_searchAllSeeds - start_seed_searchAllSeeds);
							time_seed_searchAllSeeds += elapsed_seed_searchAllSeeds.count();
#endif
							assert(shs[mate].repOk(&ca.current()));
							if (shs[mate].empty())
							{
								done[mate] = true;
								break;
							}
						}
						for (size_t mate = 0; mate < 2; mate++)
						{
							if (!shs[mate].empty())
							{
								nUniqueSeeds += shs[mate].numUniqueSeeds();
								nUniqueSeedsMS[mate * 2 + 0] += shs[mate].numUniqueSeedsStrand(true);
								nUniqueSeedsMS[mate * 2 + 1] += shs[mate].numUniqueSeedsStrand(false);
								nRepeatSeeds += shs[mate].numRepeatSeeds();
								nRepeatSeedsMS[mate * 2 + 0] += shs[mate].numRepeatSeedsStrand(true);
								nRepeatSeedsMS[mate * 2 + 1] += shs[mate].numRepeatSeedsStrand(false);
								seedHitTot += shs[mate].numElts();
								seedHitTotMS[mate * 2 + 0] += shs[mate].numEltsFw();
								seedHitTotMS[mate * 2 + 1] += shs[mate].numEltsRc();
							}
						}
						double uniqFactor[2] = {0.0f, 0.0f};
						for (size_t i = 0; i < 2; i++)
						{
							if (!shs[i].empty())
							{
								swmSeed.sdsucc++;
								uniqFactor[i] = shs[i].uniquenessFactor();
							}
						}
						matemap[0] = 0;
						matemap[1] = 1;
						if (!shs[0].empty() && !shs[1].empty() && uniqFactor[1] > uniqFactor[0])
						{
							matemap[0] = 1;
							matemap[1] = 0;
						}
						for (size_t matei = 0; matei < (paired ? 2 : 1); matei++)
						{
							size_t mate = matemap[matei];
							if (done[mate] || msinkwrap.state().doneWithMate(mate == 0))
							{
								done[mate] = true;
								continue;
							}
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							assert(shs[mate].repOk(&ca.current()));
							if (!seedSumm)
							{
								if (shs[mate].empty())
								{
									continue;
								}
								shs[mate].rankSeedHits(rnd, msinkwrap.allHits());
#ifdef TIME_STATS
								allseed_Nums++;
#endif
								int ret = 0;
								if (paired)
								{
									ret = sd.extendSeedsPaired(
										*rds[mate],
										*rds[mate ^ 1],
										mate == 0,
										!filt[mate ^ 1],
										shs[mate],
										ebwtFw,
										ebwtBw,
										ref,
										sw,
										osw,
										sc,
										pepol,
										multiseedMms,
										seedlens[mate],
										interval[mate],
										minsc[mate],
										minsc[mate ^ 1],
										nceil[mate],
										nceil[mate ^ 1],
										nofw[mate],
										norc[mate],
										maxhalf,
										doUngapped,
										mxIter[mate],
										mxUg[mate],
										mxDp[mate],
										streak[mate],
										streak[mate],
										streak[mate],
										mtStreak[mate],
										doExtend,
										enable8,
										cminlen,
										cpow2,
										doTri,
										tighten,
										ca,
										rnd,
										wlm,
										swmSeed,
										swmMate,
										prm,
										&msinkwrap,
										true,
										true,
										gReportDiscordant,
										gReportMixed,
										exhaustive[mate]);
								}
								else
								{
#ifdef TIME_STATS
									auto start_seed_extendSeeds = std::chrono::system_clock::now();
#endif
									ret = sd.extendSeeds(
										*rds[mate],
										mate == 0,
										shs[mate],
										ebwtFw,
										ebwtBw,
										ref,
										sw,
										sc,
										multiseedMms,
										seedlens[mate],
										interval[mate],
										minsc[mate],
										nceil[mate],
										maxhalf,
										doUngapped,
										mxIter[mate],
										mxUg[mate],
										mxDp[mate],
										streak[mate],
										streak[mate],
										doExtend,
										enable8,
										cminlen,
										cpow2,
										doTri,
										tighten,
										ca,
										rnd,
										wlm,
										swmSeed,
										prm,
										&msinkwrap,
										true,
										exhaustive[mate]);
#ifdef TIME_STATS
									auto end_seed_extendSeeds = std::chrono::system_clock::now();
									auto elapsed_seed_extendSeeds = std::chrono::duration_cast<std::chrono::microseconds>(end_seed_extendSeeds - start_seed_extendSeeds);
									time_seed_extendSeeds += elapsed_seed_extendSeeds.count();
#endif
								}
								assert_gt(ret, 0);
								MERGE_SW(sw);
								MERGE_SW(osw);
								if (ret == EXTEND_EXHAUSTED_CANDIDATES)
								{
								}
								else if (ret == EXTEND_POLICY_FULFILLED)
								{
									if (msinkwrap.state().doneWithMate(mate == 0))
									{
										done[mate] = true;
									}
									if (msinkwrap.state().doneWithMate(mate == 1))
									{
										done[mate ^ 1] = true;
									}
								}
								else if (ret == EXTEND_PERFECT_SCORE)
								{
									done[mate] = true;
								}
								else if (ret == EXTEND_EXCEEDED_HARD_LIMIT)
								{
									done[mate] = true;
								}
								else if (ret == EXTEND_EXCEEDED_SOFT_LIMIT)
								{
								}
								else
								{
									cerr << "Bad return value: " << ret << endl;
									throw 1;
								}
							}
						}
						for (size_t mate = 0; mate < 2; mate++)
						{
							if (!done[mate] && shs[mate].averageHitsPerSeed() < seedBoostThresh)
							{
								done[mate] = true;
							}
						}
					}
#ifdef TIME_STATS
					auto end_seed = std::chrono::system_clock::now();
					auto elapsed_seed = std::chrono::duration_cast<std::chrono::microseconds>(end_seed - start_seed);
					time_seed += elapsed_seed.count();
#endif
					if (seedsTried > 0)
					{
						prm.seedPctUnique = (float)nUniqueSeeds / seedsTried;
						prm.seedPctRep = (float)nRepeatSeeds / seedsTried;
						prm.seedHitAvg = (float)seedHitTot / seedsTried;
					}
					else
					{
						prm.seedPctUnique = -1.0f;
						prm.seedPctRep = -1.0f;
						prm.seedHitAvg = -1.0f;
					}
					for (int i = 0; i < 4; i++)
					{
						if (seedsTriedMS[i] > 0)
						{
							prm.seedPctUniqueMS[i] = (float)nUniqueSeedsMS[i] / seedsTriedMS[i];
							prm.seedPctRepMS[i] = (float)nRepeatSeedsMS[i] / seedsTriedMS[i];
							prm.seedHitAvgMS[i] = (float)seedHitTotMS[i] / seedsTriedMS[i];
						}
						else
						{
							prm.seedPctUniqueMS[i] = -1.0f;
							prm.seedPctRepMS[i] = -1.0f;
							prm.seedHitAvgMS[i] = -1.0f;
						}
					}
					size_t totnucs = 0;
					for (size_t mate = 0; mate < (paired ? 2 : 1); mate++)
					{
						if (filt[mate])
						{
							size_t len = rdlens[mate];
							if (!nofw[mate] && !norc[mate])
							{
								len *= 2;
							}
							totnucs += len;
						}
					}
					prm.seedsPerNuc = totnucs > 0 ? ((float)seedsTried / totnucs) : -1;
					for (int i = 0; i < 4; i++)
					{
						prm.seedsPerNucMS[i] = totnucs > 0 ? ((float)seedsTriedMS[i] / totnucs) : -1;
					}
					for (size_t i = 0; i < 2; i++)
					{
						assert_leq(prm.nExIters, mxIter[i]);
						assert_leq(prm.nExDps, mxDp[i]);
						assert_leq(prm.nMateDps, mxDp[i]);
						assert_leq(prm.nExUgs, mxUg[i]);
						assert_leq(prm.nMateUgs, mxUg[i]);
						assert_leq(prm.nDpFail, streak[i]);
						assert_leq(prm.nUgFail, streak[i]);
						assert_leq(prm.nEeFail, streak[i]);
					}
					msinkwrap.finishRead(
						&shs[0],
						&shs[1],
						exhaustive[0],
						exhaustive[1],
						nfilt[0],
						nfilt[1],
						scfilt[0],
						scfilt[1],
						lenfilt[0],
						lenfilt[1],
						qcfilt[0],
						qcfilt[1],
						rnd,
						rpm,
						prm,
						sc,
						!seedSumm,
						seedSumm,
						scUnMapped,
						xeq);
					assert(!retry || msinkwrap.empty());
				}
			}
			else if (rdid >= qUpto)
			{
				break;
			}
			if (metricsPerRead)
			{
				MERGE_METRICS(metricsPt);
				nametmp = ps->read_a().name;
				metricsPt.reportInterval(
					metricsOfb, metricsStderr, true, &nametmp);
				metricsPt.reset();
			}
		}
		MERGE_METRICS(metrics);
		if (dpLog != NULL)
			dpLog->close();
		if (dpLogOpp != NULL)
			dpLogOpp->close();
#ifdef PER_THREAD_TIMING
		ss.str("");
		ss.clear();
		ss << "thread: " << tid << " cpu_changeovers: " << ncpu_changeovers << std::endl
		   << "thread: " << tid << " node_changeovers: " << nnuma_changeovers << std::endl;
		std::cout << ss.str();
#endif
	}
#ifdef WITH_TBB
	p->done->fetch_add(1);
#endif
#ifdef TIME_STATS
	auto end_multiseedSearchWorker = std::chrono::system_clock::now();
	auto elapsed_multiseedSearchWorker = std::chrono::duration_cast<std::chrono::microseconds>(end_multiseedSearchWorker - start_multiseedSearchWorker);
	time_multiseedSearchWorker += elapsed_multiseedSearchWorker.count();
#endif
	return;
}
#ifdef WITH_TBB
static void multiseedSearchWorker_2p5(void *vp)
{
	thread_tracking_pair *p = (thread_tracking_pair *)vp;
	int tid = p->tid;
#else
static void multiseedSearchWorker_2p5(void *vp)
{
	int tid = *((int *)vp);
#endif
	assert(multiseed_ebwtFw != NULL);
	assert(multiseedMms == 0 || multiseed_ebwtBw != NULL);
	PatternComposer &patsrc = *multiseed_patsrc;
	PatternParams pp = multiseed_pp;
	const Ebwt &ebwtFw = *multiseed_ebwtFw;
	const Ebwt &ebwtBw = *multiseed_ebwtBw;
	const Scoring &sc = *multiseed_sc;
	const BitPairReference &ref = *multiseed_refs;
	AlnSink &msink = *multiseed_msink;
	OutFileBuf *metricsOfb = multiseed_metricsOfb;
	ThreadCounter tc;
	unique_ptr<PatternSourcePerThreadFactory> patsrcFact(createPatsrcFactory(patsrc, pp, tid));
	unique_ptr<PatternSourcePerThread> ps(patsrcFact->create());
	ReportingParams rp(
		(allHits ? std::numeric_limits<THitInt>::max() : khits),
		mhits,
		0,
		msample,
		gReportDiscordant,
		gReportMixed);
	unique_ptr<Mapq> bmapq(new_mapq(mapqv, scoreMin, sc));
	AlnSinkWrap msinkwrap(
		msink,
		rp,
		*bmapq,
		(size_t)tid);
	OuterLoopMetrics olm;
	SeedSearchMetrics sdm;
	WalkMetrics wlm;
	SwMetrics swmSeed, swmMate;
	DescentMetrics descm;
	ReportingMetrics rpm;
	RandomSource rnd, rndArb;
	SSEMetrics sseU8ExtendMet;
	SSEMetrics sseU8MateMet;
	SSEMetrics sseI16ExtendMet;
	SSEMetrics sseI16MateMet;
	uint64_t nbtfiltst = 0;
	uint64_t nbtfiltsc = 0;
	uint64_t nbtfiltdo = 0;
	ASSERT_ONLY(BTDnaString tmp);
	int pepolFlag;
	if (gMate1fw && gMate2fw)
	{
		pepolFlag = PE_POLICY_FF;
	}
	else if (gMate1fw && !gMate2fw)
	{
		pepolFlag = PE_POLICY_FR;
	}
	else if (!gMate1fw && gMate2fw)
	{
		pepolFlag = PE_POLICY_RF;
	}
	else
	{
		pepolFlag = PE_POLICY_RR;
	}
	assert_geq(gMaxInsert, gMinInsert);
	assert_geq(gMinInsert, 0);
	PairedEndPolicy pepol(
		pepolFlag,
		gMaxInsert,
		gMinInsert,
		localAlign,
		gFlippedMatesOK,
		gDovetailMatesOK,
		gContainMatesOK,
		gOlapMatesOK,
		gExpandToFrag);
	AlignerDriver ald(
		descConsExp,
		descPrioritizeRoots,
		msIval,
		descLanding,
		gVerbose,
		descentTotSz,
		descentTotFmops);
	PerfMetrics metricsPt;
	BTString nametmp;
	PerReadMetrics prm;
	time_t iTime = time(0);
	bool exhaustive[2] = {false, false};
	bool filt[2] = {true, true};
	bool nfilt[2] = {true, true};
	bool scfilt[2] = {true, true};
	bool lenfilt[2] = {true, true};
	bool qcfilt[2] = {true, true};
	rndArb.init((uint32_t)time(0));
	int mergei = 0;
	int mergeival = 16;
	while (true)
	{
		pair<bool, bool> ret = ps->nextReadPair();
		bool success = ret.first;
		bool done = ret.second;
		if (!success && done)
		{
			break;
		}
		else if (!success)
		{
			continue;
		}
		TReadId rdid = ps->read_a().rdid;
		bool sample = true;
		if (arbitraryRandom)
		{
			ps->read_a().seed = rndArb.nextU32();
			ps->read_b().seed = rndArb.nextU32();
		}
		if (sampleFrac < 1.0f)
		{
			rnd.init(ROTL(ps->read_a().seed, 2));
			sample = rnd.nextFloat() < sampleFrac;
		}
		if (rdid >= skipReads && rdid < qUpto && sample)
		{
			if (metricsIval > 0 &&
				(metricsOfb != NULL || metricsStderr) &&
				!metricsPerRead &&
				++mergei == mergeival)
			{
				MERGE_METRICS(metrics);
				mergei = 0;
				if (tid == 0)
				{
					time_t curTime = time(0);
					if (curTime - iTime >= metricsIval)
					{
						metrics.reportInterval(metricsOfb, metricsStderr, false, NULL);
						iTime = curTime;
					}
				}
			}
			prm.reset();
			prm.doFmString = sam_print_zm;
			if (sam_print_xt)
			{
				gettimeofday(&prm.tv_beg, &prm.tz_beg);
			}
			olm.reads++;
			bool paired = !ps->read_b().empty();
			const size_t rdlen1 = ps->read_a().length();
			const size_t rdlen2 = paired ? ps->read_b().length() : 0;
			olm.bases += (rdlen1 + rdlen2);
			rnd.init(ROTL(ps->read_a().seed, 5));
			msinkwrap.nextRead(
				&ps->read_a(),
				paired ? &ps->read_b() : NULL,
				rdid,
				sc.qualitiesMatter());
			assert(msinkwrap.inited());
			size_t rdlens[2] = {rdlen1, rdlen2};
			TAlScore minsc[2], maxpen[2];
			minsc[0] = minsc[1] = std::numeric_limits<TAlScore>::max();
			setupMinScores(*ps, paired, localAlign, sc, rdlens, minsc, maxpen);
			size_t readns[2] = {0, 0};
			sc.nFilterPair(
				&ps->read_a().patFw,
				paired ? &ps->read_b().patFw : NULL,
				readns[0],
				readns[1],
				nfilt[0],
				nfilt[1]);
			scfilt[0] = sc.scoreFilter(minsc[0], rdlens[0]);
			scfilt[1] = sc.scoreFilter(minsc[1], rdlens[1]);
			lenfilt[0] = lenfilt[1] = true;
			if (rdlens[0] <= (size_t)multiseedMms || rdlens[0] < 2)
			{
				if (!gQuiet)
					printMmsSkipMsg(*ps, paired, true, multiseedMms);
				lenfilt[0] = false;
			}
			if ((rdlens[1] <= (size_t)multiseedMms || rdlens[1] < 2) && paired)
			{
				if (!gQuiet)
					printMmsSkipMsg(*ps, paired, false, multiseedMms);
				lenfilt[1] = false;
			}
			if (rdlens[0] < 2)
			{
				if (!gQuiet)
					printLenSkipMsg(*ps, paired, true);
				lenfilt[0] = false;
			}
			if (rdlens[1] < 2 && paired)
			{
				if (!gQuiet)
					printLenSkipMsg(*ps, paired, false);
				lenfilt[1] = false;
			}
			qcfilt[0] = qcfilt[1] = true;
			if (qcFilter)
			{
				qcfilt[0] = (ps->read_a().filter != '0');
				qcfilt[1] = (ps->read_b().filter != '0');
			}
			filt[0] = (nfilt[0] && scfilt[0] && lenfilt[0] && qcfilt[0]);
			filt[1] = (nfilt[1] && scfilt[1] && lenfilt[1] && qcfilt[1]);
			prm.nFilt += (filt[0] ? 0 : 1) + (filt[1] ? 0 : 1);
			Read *rds[2] = {&ps->read_a(), &ps->read_b()};
			assert(msinkwrap.empty());
			bool nofw[2] = {false, false};
			bool norc[2] = {false, false};
			nofw[0] = paired ? (gMate1fw ? gNofw : gNorc) : gNofw;
			norc[0] = paired ? (gMate1fw ? gNorc : gNofw) : gNorc;
			nofw[1] = paired ? (gMate2fw ? gNofw : gNorc) : gNofw;
			norc[1] = paired ? (gMate2fw ? gNorc : gNofw) : gNorc;
			int nceil[2] = {0, 0};
			nceil[0] = nCeil.f<int>((double)rdlens[0]);
			nceil[0] = min(nceil[0], (int)rdlens[0]);
			if (paired)
			{
				nceil[1] = nCeil.f<int>((double)rdlens[1]);
				nceil[1] = min(nceil[1], (int)rdlens[1]);
			}
			exhaustive[0] = exhaustive[1] = false;
			bool pairPostFilt = filt[0] && filt[1];
			if (pairPostFilt)
			{
				rnd.init(ROTL((rds[0]->seed ^ rds[1]->seed), 10));
			}
			size_t streak[2] = {maxDpStreak, maxDpStreak};
			size_t mtStreak[2] = {maxMateStreak, maxMateStreak};
			size_t mxDp[2] = {maxDp, maxDp};
			size_t mxUg[2] = {maxUg, maxUg};
			size_t mxIter[2] = {maxIters, maxIters};
			if (allHits)
			{
				streak[0] = streak[1] = std::numeric_limits<size_t>::max();
				mtStreak[0] = mtStreak[1] = std::numeric_limits<size_t>::max();
				mxDp[0] = mxDp[1] = std::numeric_limits<size_t>::max();
				mxUg[0] = mxUg[1] = std::numeric_limits<size_t>::max();
				mxIter[0] = mxIter[1] = std::numeric_limits<size_t>::max();
			}
			else if (khits > 1)
			{
				for (size_t mate = 0; mate < 2; mate++)
				{
					streak[mate] += (khits - 1) * maxStreakIncr;
					mtStreak[mate] += (khits - 1) * maxStreakIncr;
					mxDp[mate] += (khits - 1) * maxItersIncr;
					mxUg[mate] += (khits - 1) * maxItersIncr;
					mxIter[mate] += (khits - 1) * maxItersIncr;
				}
			}
			if (filt[0] && filt[1])
			{
				streak[0] = (size_t)ceil((double)streak[0] / 2.0);
				streak[1] = (size_t)ceil((double)streak[1] / 2.0);
				assert_gt(streak[1], 0);
			}
			assert_gt(streak[0], 0);
			for (size_t mate = 0; mate < (paired ? 2 : 1); mate++)
			{
				if (!filt[mate])
				{
					olm.freads++;
					olm.fbases += rdlens[mate];
				}
				else
				{
					olm.ureads++;
					olm.ubases += rdlens[mate];
				}
			}
			if (filt[0])
			{
				ald.initRead(ps->read_a(), nofw[0], norc[0], minsc[0], maxpen[0], filt[1] ? &ps->read_b() : NULL);
			}
			else if (filt[1])
			{
				ald.initRead(ps->read_b(), nofw[1], norc[1], minsc[1], maxpen[1], NULL);
			}
			if (filt[0] || filt[1])
			{
				ald.go(sc, ebwtFw, ebwtBw, ref, descm, wlm, prm, rnd, msinkwrap);
			}
			uint32_t sd = rds[0]->seed ^ rds[1]->seed;
			rnd.init(ROTL(sd, 20));
			msinkwrap.finishRead(
				NULL,
				NULL,
				exhaustive[0],
				exhaustive[1],
				nfilt[0],
				nfilt[1],
				scfilt[0],
				scfilt[1],
				lenfilt[0],
				lenfilt[1],
				qcfilt[0],
				qcfilt[1],
				rnd,
				rpm,
				prm,
				sc,
				!seedSumm,
				seedSumm,
				scUnMapped,
				xeq);
		}
		else if (rdid >= qUpto)
		{
			break;
		}
		if (metricsPerRead)
		{
			MERGE_METRICS(metricsPt);
			nametmp = ps->read_a().name;
			metricsPt.reportInterval(
				metricsOfb, metricsStderr, true, &nametmp);
			metricsPt.reset();
		}
	}
	MERGE_METRICS(metrics);
#ifdef WITH_TBB
	p->done->fetch_add(1);
#endif
	return;
}
#ifndef _WIN32
static void errno_message()
{
	int errnum = errno;
	cerr << "errno is " << errnum << endl;
	perror("perror error: ");
}
void del_pid(const char *dirname, int pid)
{
	char *fname = (char *)calloc(FNAME_SIZE, sizeof(char));
	if (fname == NULL)
	{
		errno_message();
		cerr << "del_pid: could not allocate buffer" << endl;
		throw 1;
	}
	snprintf(fname, FNAME_SIZE, "%s/%d", dirname, pid);
	if (unlink(fname) != 0)
	{
		if (errno != ENOENT)
		{
			errno_message();
			cerr << "del_pid: could not delete PID file " << fname << endl;
			free(fname);
			throw 1;
		}
		else
		{
		}
	}
	free(fname);
}
static void write_pid(const char *dirname, int pid)
{
	struct stat dinfo;
	if (stat(dirname, &dinfo) != 0)
	{
		if (mkdir(dirname, 0755) != 0)
		{
			if (errno != EEXIST)
			{
				errno_message();
				cerr << "write_pid: could not create PID directory " << dirname << endl;
				throw 1;
			}
		}
	}
	char *fname = (char *)calloc(FNAME_SIZE, sizeof(char));
	if (fname == NULL)
	{
		errno_message();
		cerr << "write_pid: could not allocate buffer" << endl;
		throw 1;
	}
	snprintf(fname, FNAME_SIZE, "%s/%d", dirname, pid);
	FILE *f = fopen(fname, "w");
	if (f == NULL)
	{
		errno_message();
		cerr << "write_pid: could not open PID file " << fname << endl;
		throw 1;
	}
	if (fclose(f) != 0)
	{
		errno_message();
		cerr << "write_pid: could not close PID file " << fname << endl;
		throw 1;
	}
	free(fname);
}
static int read_dir(const char *dirname, int *num_pids)
{
	DIR *dir;
	struct dirent *ent;
	char *fname = (char *)calloc(FNAME_SIZE, sizeof(char));
	if (fname == NULL)
	{
		errno_message();
		cerr << "read_dir: could not allocate buffer" << endl;
		throw 1;
	}
	dir = opendir(dirname);
	if (dir == NULL)
	{
		errno_message();
		cerr << "read_dir: could not open directory " << dirname << endl;
		free(fname);
		throw 1;
	}
	int lowest_pid = -1;
	while (true)
	{
		errno = 0;
		ent = readdir(dir);
		if (ent == NULL)
		{
			if (errno != 0)
			{
				errno_message();
				cerr << "read_dir: could not read directory " << dirname << endl;
				free(fname);
				throw 1;
			}
			break;
		}
		if (ent->d_name[0] == '.')
		{
			continue;
		}
		int pid = atoi(ent->d_name);
		if (kill(pid, 0) != 0)
		{
			if (errno == ESRCH)
			{
				del_pid(dirname, pid);
				continue;
			}
			else
			{
				errno_message();
				cerr << "read_dir: could not interrogate pid " << pid << endl;
				free(fname);
				throw 1;
			}
		}
		(*num_pids)++;
		if (pid < lowest_pid || lowest_pid == -1)
		{
			lowest_pid = pid;
		}
	}
	if (closedir(dir) != 0)
	{
		errno_message();
		cerr << "read_dir: could not close directory " << dir << endl;
		free(fname);
		throw 1;
	}
	free(fname);
	return lowest_pid;
}
template <typename T>
static void steal_threads(int pid, int orig_nthreads, EList<int> &tids, EList<T *> &threads)
{
	int ncpu = thread_ceiling;
	if (thread_ceiling <= nthreads)
	{
		return;
	}
	int num_pids = 0;
	int lowest_pid = read_dir(thread_stealing_dir.c_str(), &num_pids);
	if (lowest_pid != pid)
	{
		return;
	}
	int in_use = ((num_pids - 1) * orig_nthreads) + nthreads;
	if (in_use < ncpu)
	{
		nthreads++;
		tids.push_back(nthreads);
		threads.push_back(new T(multiseedSearchWorker, (void *)&tids.back()));
		cerr << "pid " << pid << " started new worker # " << nthreads << endl;
	}
}
template <typename T>
static void thread_monitor(int pid, int orig_threads, EList<int> &tids, EList<T *> &threads)
{
	for (int j = 0; j < 10; j++)
	{
		sleep(1);
	}
	int steal_ctr = 1;
	while (thread_counter > 0)
	{
		steal_threads(pid, orig_threads, tids, threads);
		steal_ctr++;
		for (int j = 0; j < 10; j++)
		{
			sleep(1);
		}
	}
}
#endif
static void multiseedSearch(
	Scoring &sc,
	const PatternParams &pp,
	PatternComposer &patsrc,
	AlnSink &msink,
	Ebwt &ebwtFw,
	Ebwt *ebwtBw,
	OutFileBuf *metricsOfb)
{
	multiseed_patsrc = &patsrc;
	multiseed_pp = pp;
	multiseed_msink = &msink;
	multiseed_ebwtFw = &ebwtFw;
	multiseed_ebwtBw = ebwtBw;
	multiseed_sc = &sc;
	multiseed_metricsOfb = metricsOfb;
	Timer *_t = new Timer(cerr, "Time loading reference: ", timing);
	unique_ptr<BitPairReference> refs(
		new BitPairReference(
			adjIdxBase,
			false,
			sanityCheck,
			NULL,
			NULL,
			false,
			useMm,
			useShmem,
			mmSweep,
			gVerbose,
			startVerbose));
	delete _t;
	if (!refs->loaded())
		throw 1;
	multiseed_refs = refs.get();
#ifndef _WIN32
	sigset_t set;
	sigemptyset(&set);
	sigaddset(&set, SIGPIPE);
	pthread_sigmask(SIG_BLOCK, &set, NULL);
#endif
	EList<int> tids;
#ifdef WITH_TBB
	EList<std::thread *> threads(nthreads);
	EList<thread_tracking_pair> tps;
	tps.resize(std::max(nthreads, thread_ceiling));
#else
	EList<tthread::thread *> threads;
#endif
	threads.reserveExact(std::max(nthreads, thread_ceiling));
	tids.reserveExact(std::max(nthreads, thread_ceiling));
	{
		assert(!ebwtFw.isInMemory());
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(
			0,
			-1,
			true,
			true,
			true,
			!noRefNames,
			startVerbose);
	}
	if (multiseedMms > 0 || do1mmUpFront)
	{
		assert(!ebwtBw->isInMemory());
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw->loadIntoMemory(
			0,
			1,
			false,
			true,
			false,
			!noRefNames,
			startVerbose);
	}
#ifdef WITH_TBB
	std::atomic<int> all_threads_done;
	all_threads_done = 0;
#endif
	{
		Timer _t(cerr, "Multiseed full-index search: ", timing);
#ifndef _WIN32
		int pid = 0;
		if (thread_stealing)
		{
			pid = getpid();
			write_pid(thread_stealing_dir.c_str(), pid);
			thread_counter = 0;
		}
#endif
		for (int i = 0; i < nthreads; i++)
		{
			tids.push_back(i);
#ifdef WITH_TBB
			tps[i].tid = i;
			tps[i].done = &all_threads_done;
			if (nextVer)
			{
				threads.push_back(new std::thread(multiseedSearchWorker_2p5, (void *)&tps[i]));
			}
			else
			{
				threads.push_back(new std::thread(multiseedSearchWorker, (void *)&tps[i]));
			}
			threads[i]->detach();
			SLEEP(10);
#else
			if (nextVer)
			{
				threads.push_back(new tthread::thread(multiseedSearchWorker_2p5, (void *)&tids.back()));
			}
			else
			{
				threads.push_back(new tthread::thread(multiseedSearchWorker, (void *)&tids.back()));
			}
#endif
		}
#ifndef _WIN32
		if (thread_stealing)
		{
			int orig_threads = nthreads;
			thread_monitor(pid, orig_threads, tids, threads);
		}
#endif
#ifdef WITH_TBB
		while (all_threads_done < nthreads)
		{
			SLEEP(10);
		}
#else
		for (int i = 0; i < nthreads; i++)
		{
			threads[i]->join();
		}
#endif
		for (int i = 0; i < nthreads; ++i)
		{
			delete threads[i];
		}
#ifndef _WIN32
		if (thread_stealing)
		{
			del_pid(thread_stealing_dir.c_str(), pid);
		}
#endif
	}
	if (!metricsPerRead && (metricsOfb != NULL || metricsStderr))
	{
		metrics.reportInterval(metricsOfb, metricsStderr, true, NULL);
	}
}
static string argstr;
template <typename TStr>
static void driver(
	const char *type,
	const string &bt2indexBase,
	const string &outfile)
{
	if (gVerbose || startVerbose)
	{
		cerr << "Entered driver(): ";
		logTime(cerr, true);
	}
	EList<SString<char>> names, os;
	EList<size_t> nameLens, seqLens;
	if (!origString.empty())
	{
		EList<string> origFiles;
		tokenize(origString, ",", origFiles);
		parseFastas(origFiles, names, nameLens, os, seqLens);
	}
	PatternParams pp(
		format,
		interleaved,
		fileParallel,
		seed,
		readsPerBatch,
		solexaQuals,
		phred64Quals,
		integerQuals,
		gTrim5,
		gTrim3,
		trimTo,
		fastaContLen,
		fastaContFreq,
		skipReads,
		qUpto,
		nthreads,
		outType != OUTPUT_SAM,
		preserve_tags,
		align_paired_reads);
	if (gVerbose || startVerbose)
	{
		cerr << "Creating PatternSource: ";
		logTime(cerr, true);
	}
	PatternComposer *patsrc = PatternComposer::setupPatternComposer(
		queries,
		mates1,
		mates2,
		mates12,
		qualities,
		qualities1,
		qualities2,
#ifdef USE_SRA
		sra_accs,
#endif
		pp,
		gVerbose || startVerbose);
	if (gVerbose || startVerbose)
	{
		cerr << "Opening hit output file: ";
		logTime(cerr, true);
	}
	OutFileBuf *fout;
	if (!outfile.empty())
	{
		fout = new OutFileBuf(outfile.c_str(), false);
	}
	else
	{
		fout = new OutFileBuf();
	}
	if (gVerbose || startVerbose)
	{
		cerr << "About to initialize fw Ebwt: ";
		logTime(cerr, true);
	}
	adjIdxBase = adjustEbwtBase(argv0, bt2indexBase, gVerbose);
	Ebwt ebwt(
		adjIdxBase,
		0,
		-1,
		true,
		offRate,
		0,
		useMm,
		useShmem,
		mmSweep,
		!noRefNames,
		true,
		true,
		true,
		gVerbose,
		startVerbose,
		false,
		sanityCheck);
	if (sanityCheck && !os.empty())
	{
		assert_eq(os.size(), ebwt.nPat());
		for (size_t i = 0; i < os.size(); i++)
		{
			assert_eq(os[i].length(), ebwt.plen()[i]);
		}
	}
	if (sanityCheck && !os.empty())
	{
		ebwt.loadIntoMemory(
			0,
			-1,
			true,
			true,
			true,
			!noRefNames,
			startVerbose);
		ebwt.checkOrigs(os, false, false);
		ebwt.evictFromMemory();
	}
	OutputQueue oq(
		*fout,
		reorder && (nthreads > 1 || thread_stealing),
		nthreads,
		nthreads > 1 || thread_stealing,
		readsPerBatch,
		skipReads);
	{
		Timer _t(cerr, "Time searching: ", timing);
		if (bonusMatch > 0 && !localAlign)
		{
			cerr << "Warning: Match bonus always = 0 in --end-to-end mode; ignoring user setting" << endl;
			bonusMatch = 0;
		}
		Scoring sc(
			bonusMatch,
			penMmcType,
			penMmcMax,
			penMmcMin,
			scoreMin,
			nCeil,
			penNType,
			penN,
			penNCatPair,
			penRdGapConst,
			penRfGapConst,
			penRdGapLinear,
			penRfGapLinear,
			gGapBarrier);
		EList<size_t> reflens;
		for (size_t i = 0; i < ebwt.nPat(); i++)
		{
			reflens.push_back(ebwt.plen()[i]);
		}
		EList<string> refnames;
		readEbwtRefnames(adjIdxBase, refnames);
		SamConfig samc(
			refnames,
			reflens,
			samTruncQname,
			samAppendComment,
			samOmitSecSeqQual,
			samNoUnal,
			string("effaln"),
			string("1.0"),
			argstr,
			rgs_optflag,
			sam_print_as,
			sam_print_xs,
			sam_print_xss,
			sam_print_yn,
			sam_print_xn,
			sam_print_x0,
			sam_print_x1,
			sam_print_xm,
			sam_print_xo,
			sam_print_xg,
			sam_print_nm,
			sam_print_md,
			sam_print_yf,
			sam_print_yi,
			sam_print_ym,
			sam_print_yp,
			sam_print_yt,
			sam_print_ys,
			sam_print_zs,
			sam_print_xr,
			sam_print_xt,
			sam_print_xd,
			sam_print_xu,
			sam_print_yl,
			sam_print_ye,
			sam_print_yu,
			sam_print_xp,
			sam_print_yr,
			sam_print_zb,
			sam_print_zr,
			sam_print_zf,
			sam_print_zm,
			sam_print_zi,
			sam_print_zp,
			sam_print_zu,
			sam_print_zt);
		AlnSink *mssink = NULL;
		switch (outType)
		{
		case OUTPUT_SAM:
		{
			mssink = new AlnSinkSam(
				oq,
				samc,
				refnames,
				gQuiet);
			if (!samNoHead)
			{
				bool printHd = true, printSq = true;
				BTString buf;
				samc.printHeader(buf, rgid, rgs, printHd, !samNoSQ, printSq);
				fout->writeString(buf);
			}
			break;
		}
		default:
			cerr << "Invalid output type: " << outType << endl;
			throw 1;
		}
		if (gVerbose || startVerbose)
		{
			cerr << "Dispatching to search driver: ";
			logTime(cerr, true);
		}
		OutFileBuf *metricsOfb = NULL;
		if (!metricsFile.empty() && metricsIval > 0)
		{
			metricsOfb = new OutFileBuf(metricsFile);
		}
		assert(patsrc != NULL);
		assert(mssink != NULL);
		if (multiseedMms > 0 || do1mmUpFront)
		{
			if (gVerbose || startVerbose)
			{
				cerr << "About to initialize rev Ebwt: ";
				logTime(cerr, true);
			}
			Ebwt ebwtBw = Ebwt(
				adjIdxBase + ".rev",
				0,
				1,
				false,
				offRate,
				0,
				useMm,
				useShmem,
				mmSweep,
				!noRefNames,
				true,
				true,
				true,
				gVerbose,
				startVerbose,
				false,
				sanityCheck);
			multiseedSearch(
				sc,
				pp,
				*patsrc,
				*mssink,
				ebwt,
				&ebwtBw,
				metricsOfb);
		}
		else
		{
			multiseedSearch(
				sc,
				pp,
				*patsrc,
				*mssink,
				ebwt,
				NULL,
				metricsOfb);
		}
		if (ebwt.isInMemory())
		{
			ebwt.evictFromMemory();
		}
		if (!gQuiet && !seedSumm)
		{
			size_t repThresh = mhits;
			if (repThresh == 0)
			{
				repThresh = std::numeric_limits<size_t>::max();
			}
			mssink->finish(
				repThresh,
				gReportDiscordant,
				gReportMixed,
				hadoopOut);
		}
		oq.flush(true);
		assert_eq(oq.numStarted(), oq.numFinished());
		assert_eq(oq.numStarted(), oq.numFlushed());
		delete patsrc;
		delete mssink;
		delete metricsOfb;
		if (fout != NULL)
		{
			delete fout;
		}
	}
}
extern "C"
{
	int effaln(int argc, const char **argv)
	{
#ifdef TIME_STATS
		auto start_effaln = std::chrono::system_clock::now();
#endif
		try
		{
#ifdef WITH_TBB
#ifdef WITH_AFFINITY
			pinning_observer pinner(2);
			pinner.observe(true);
#endif
#endif
			opterr = optind = 1;
			resetOptions();
			for (int i = 0; i < argc; i++)
			{
				argstr += argv[i];
				if (i < argc - 1)
					argstr += " ";
			}
			if (startVerbose)
			{
				cerr << "Entered main(): ";
				logTime(cerr, true);
			}
			parseOptions(argc, argv);
			argv0 = argv[0];
			if (showVersion)
			{
				cout << argv0 << " version 1.0" << endl;
				if (sizeof(void *) == 4)
				{
					cout << "32-bit" << endl;
				}
				else if (sizeof(void *) == 8)
				{
					cout << "64-bit" << endl;
				}
				else
				{
					cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void *) << endl;
				}
				cout << "Built on " << BUILD_HOST << endl;
				cout << BUILD_TIME << endl;
				cout << "Compiler: " << COMPILER_VERSION << endl;
				cout << "Options: " << COMPILER_OPTIONS << endl;
				cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
					 << sizeof(int)
					 << ", " << sizeof(long) << ", " << sizeof(long long)
					 << ", " << sizeof(void *) << ", " << sizeof(size_t)
					 << ", " << sizeof(off_t) << "}" << endl;
				return 0;
			}
			{
				Timer _t(cerr, "Overall time: ", timing);
				if (startVerbose)
				{
					cerr << "Parsing index and read arguments: ";
					logTime(cerr, true);
				}
				if (bt2index.empty())
				{
					cerr << "No index, query, or output file specified!" << endl;
					return 1;
				}
#ifndef _WIN32
				thread_stealing = thread_ceiling > nthreads;
#endif
				if (thread_stealing && thread_stealing_dir.empty())
				{
					cerr << "When --thread-ceiling is specified, must also specify --thread-piddir" << endl;
					return 1;
				}
				bool got_reads = !queries.empty() || !mates1.empty() || !mates12.empty();
#ifdef USE_SRA
				got_reads = got_reads || !sra_accs.empty();
#endif
				if (optind >= argc)
				{
					if (!got_reads)
					{
						cerr << "***" << endl
#ifdef USE_SRA
							 << "Error: Must specify at least one read input with -U/-1/-2/--sra-acc" << endl;
#else
							 << "Error: Must specify at least one read input with -U/-1/-2" << endl;
#endif
						return 1;
					}
				}
				else if (!got_reads)
				{
					tokenize(argv[optind++], ",", queries);
					if (queries.empty())
					{
						cerr << "Tokenized query file list was empty!" << endl;
						return 1;
					}
				}
				if (optind < argc && outfile.empty())
				{
					outfile = argv[optind++];
					cerr << "Warning: Output file '" << outfile.c_str()
						 << "' was specified without -S." << endl;
				}
				if (optind < argc)
				{
					cerr << "Extra parameter(s) specified: ";
					for (int i = optind; i < argc; i++)
					{
						cerr << "\"" << argv[i] << "\"";
						if (i < argc - 1)
							cerr << ", ";
					}
					cerr << endl;
					if (mates1.size() > 0)
					{
						cerr << "Note that if <mates> files are specified using -1/-2, a <singles> file cannot also be specified." << endl;
					}
					throw 1;
				}
				if (gVerbose)
				{
					cout << "Input " + gEbwt_ext + " file: \"" << bt2index.c_str() << "\"" << endl;
					cout << "Query inputs (DNA, " << file_format_names[format].c_str() << "):" << endl;
					for (size_t i = 0; i < queries.size(); i++)
					{
						cout << "  " << queries[i].c_str() << endl;
					}
					cout << "Quality inputs:" << endl;
					for (size_t i = 0; i < qualities.size(); i++)
					{
						cout << "  " << qualities[i].c_str() << endl;
					}
					cout << "Output file: \"" << outfile.c_str() << "\"" << endl;
					cout << "Local endianness: " << (currentlyBigEndian() ? "big" : "little") << endl;
					cout << "Sanity checking: " << (sanityCheck ? "enabled" : "disabled") << endl;
#ifdef NDEBUG
					cout << "Assertions: disabled" << endl;
#else
					cout << "Assertions: enabled" << endl;
#endif
				}
				if (ipause)
				{
					cout << "Press key to continue..." << endl;
					getchar();
				}
				driver<SString<char>>("DNA", bt2index, outfile);
#ifdef TIME_STATS
				auto end_effaln = std::chrono::system_clock::now();
				auto elapsed_effaln = std::chrono::duration_cast<std::chrono::microseconds>(end_effaln - start_effaln);
				time_effaln += elapsed_effaln.count();
#endif
#ifdef TIME_STATS
				cerr << endl
					 << "----------------------------------------" << endl;
				cerr << "LF-call-Nums:               " << setw(12) << setiosflags(ios::right) << countLF_Static << endl;
				cerr << "Align-call-Nums:            " << setw(12) << setiosflags(ios::right) << countAlign_Static << endl;
				cerr << "Gather-call-Nums:           " << setw(12) << setiosflags(ios::right) << countGather_Static << endl;
				cerr << "align_time:                 " << setw(12) << setiosflags(ios::right) << align_time << endl;
				cerr << "gather_time:                " << setw(12) << setiosflags(ios::right) << gather_time << endl;
				cerr << "0_time_effaln:              " << setw(12) << setiosflags(ios::right) << time_effaln << endl;
				cerr << "0_time_multiseedSearchWorke:" << setw(12) << setiosflags(ios::right) << time_multiseedSearchWorker << endl;
				cerr << "1_time_exact:               " << setw(12) << setiosflags(ios::right) << time_exact << endl;
				cerr << "1_time_1mm:                 " << setw(12) << setiosflags(ios::right) << time_1mm << endl;
				cerr << "1_time_seed:                " << setw(12) << setiosflags(ios::right) << time_seed << endl;
				cerr << "2_time_seed_mmSeeds:        " << setw(12) << setiosflags(ios::right) << time_seed_mmSeeds << endl;
				cerr << "2_time_seed_instantiateSeed:" << setw(12) << setiosflags(ios::right) << time_seed_instantiateSeeds << endl;
				cerr << "2_time_seed_searchAllSeeds: " << setw(12) << setiosflags(ios::right) << time_seed_searchAllSeeds << endl;
				cerr << "2_time_seed_extendSeeds:    " << setw(12) << setiosflags(ios::right) << time_seed_extendSeeds << endl;
				cerr << "3_time_extendSeeds_eeSaTups:" << setw(12) << setiosflags(ios::right) << time_extendSeeds_eeSaTups << endl;
				cerr << "3_exendSeed_prioritizeSATup:" << setw(12) << setiosflags(ios::right) << time_extendSeeds_prioritizeSATups << endl;
				cerr << "3_extendSeed_advanceElement:" << setw(12) << setiosflags(ios::right) << time_extendSeeds_advanceElement << endl;
				cerr << "3_extendSeed_joinedToTxtOff:" << setw(12) << setiosflags(ios::right) << time_extendSeeds_joinedToTextOff << endl;
				cerr << "3_redundants_Nums:          " << setw(12) << setiosflags(ios::right) << redundants_Nums << endl;
				cerr << "3_time_extendSeeds_resEe:   " << setw(12) << setiosflags(ios::right) << time_extendSeeds_resEe << endl;
				cerr << "3_time_extendSeeds_resUngap:" << setw(12) << setiosflags(ios::right) << time_extendSeeds_resUngap << endl;
				cerr << "3_time_extendSeeds_init:    " << setw(12) << setiosflags(ios::right) << time_extendSeeds_init << endl;
				cerr << "3_time_extendSeeds_align:   " << setw(12) << setiosflags(ios::right) << time_extendSeeds_align << endl;
				cerr << "4_time_cgk:                 " << setw(12) << setiosflags(ios::right) << time_wfa0 << endl;
				cerr << "4_time_wfa:                 " << setw(12) << setiosflags(ios::right) << time_wfa1 << endl;
				cerr << "3_time_extendSeeds_nxtAlign:" << setw(12) << setiosflags(ios::right) << 0 << endl;
				cerr << "exact_Nums:                 " << setw(12) << setiosflags(ios::right) << exact_Nums << endl;
				cerr << "onemm_Nums:                 " << setw(12) << setiosflags(ios::right) << onemm_Nums << endl;
				cerr << "allseed_Nums:               " << setw(12) << setiosflags(ios::right) << allseed_Nums << endl;
				cerr << "bestFilter_Nums:            " << setw(12) << setiosflags(ios::right) << bestFilter_Nums << endl;
				cerr << "----------------------------------------" << endl;
#endif
			}
#ifdef WITH_TBB
#ifdef WITH_AFFINITY
			pinner.observe(false);
#endif
#endif
			return 0;
		}
		catch (std::exception &e)
		{
			cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
			cerr << "Command: ";
			for (int i = 0; i < argc; i++)
				cerr << argv[i] << " ";
			cerr << endl;
			return 1;
		}
		catch (int e)
		{
			if (e != 0)
			{
				cerr << "Error: Encountered internal exception (#" << e << ")" << endl;
				cerr << "Command: ";
				for (int i = 0; i < argc; i++)
					cerr << argv[i] << " ";
				cerr << endl;
			}
			return e;
		}
	}
}
