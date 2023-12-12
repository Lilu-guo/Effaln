#ifndef OPTS_H_
#define OPTS_H_
enum
{
	ARG_ORIG = 256,
	ARG_SEED,
	ARG_SOLEXA_QUALS,
	ARG_VERBOSE,
	ARG_STARTVERBOSE,
	ARG_QUIET,
	ARG_METRIC_IVAL,
	ARG_METRIC_FILE,
	ARG_METRIC_STDERR,
	ARG_METRIC_PER_READ,
	ARG_REFIDX,
	ARG_SANITY,
	ARG_PARTITION,
	ARG_INTEGER_QUALS,
	ARG_FILEPAR,
	ARG_READS_PER_BATCH,
	ARG_SHMEM,
	ARG_MM,
	ARG_MMSWEEP,
	ARG_FF,
	ARG_FR,
	ARG_RF,
	ARG_NO_MIXED,
	ARG_NO_DISCORDANT,
	ARG_CACHE_LIM,
	ARG_CACHE_SZ,
	ARG_NO_FW,
	ARG_NO_RC,
	ARG_SKIP,
	ARG_ONETWO,
	ARG_PHRED64,
	ARG_PHRED33,
	ARG_HADOOPOUT,
	ARG_FULLREF,
	ARG_USAGE,
	ARG_SNPPHRED,
	ARG_SNPFRAC,
	ARG_SAM_NO_QNAME_TRUNC,
	ARG_SAM_OMIT_SEC_SEQ,
	ARG_SAM_NOHEAD,
	ARG_SAM_NOSQ,
	ARG_SAM_RG,
	ARG_SAM_RGID,
	ARG_GAP_BAR,
	ARG_QUALS1,
	ARG_QUALS2,
	ARG_QSEQ,
	ARG_SEED_SUMM,
	ARG_SC_UNMAPPED,
	ARG_OVERHANG,
	ARG_NO_CACHE,
	ARG_USE_CACHE,
	ARG_NOISY_HPOLY,
	ARG_LOCAL,
	ARG_END_TO_END,
	ARG_SCAN_NARROWED,
	ARG_QC_FILTER,
	ARG_BWA_SW_LIKE,
	ARG_MULTISEED_IVAL,
	ARG_SCORE_MIN,
	ARG_SCORE_MA,
	ARG_SCORE_MMP,
	ARG_SCORE_NP,
	ARG_SCORE_RDG,
	ARG_SCORE_RFG,
	ARG_N_CEIL,
	ARG_DPAD,
	ARG_SAM_PRINT_YI,
	ARG_ALIGN_POLICY,
	ARG_PRESET_VERY_FAST,
	ARG_PRESET_FAST,
	ARG_PRESET_SENSITIVE,
	ARG_PRESET_VERY_SENSITIVE,
	ARG_PRESET_VERY_FAST_LOCAL,
	ARG_PRESET_FAST_LOCAL,
	ARG_PRESET_SENSITIVE_LOCAL,
	ARG_PRESET_VERY_SENSITIVE_LOCAL,
	ARG_IGNORE_QUALS,
	ARG_DESC,
	ARG_TAB5,
	ARG_TAB6,
	ARG_WRAPPER,
	ARG_DOVETAIL,
	ARG_NO_DOVETAIL,
	ARG_CONTAIN,
	ARG_NO_CONTAIN,
	ARG_OVERLAP,
	ARG_NO_OVERLAP,
	ARG_MAPQ_V,
	ARG_SSE8,
	ARG_SSE8_NO,
	ARG_UNGAPPED,
	ARG_UNGAPPED_NO,
	ARG_TIGHTEN,
	ARG_UNGAP_THRESH,
	ARG_EXACT_UPFRONT,
	ARG_1MM_UPFRONT,
	ARG_EXACT_UPFRONT_NO,
	ARG_1MM_UPFRONT_NO,
	ARG_1MM_MINLEN,
	ARG_VERSION,
	ARG_SEED_OFF,
	ARG_SEED_BOOST_THRESH,
	ARG_READ_TIMES,
	ARG_EXTEND_ITERS,
	ARG_DP_MATE_STREAK_THRESH,
	ARG_DP_FAIL_STREAK_THRESH,
	ARG_UG_FAIL_STREAK_THRESH,
	ARG_EE_FAIL_STREAK_THRESH,
	ARG_DP_FAIL_THRESH,
	ARG_UG_FAIL_THRESH,
	ARG_MAPQ_EX,
	ARG_NO_EXTEND,
	ARG_REORDER,
	ARG_SHOW_RAND_SEED,
	ARG_READ_PASSTHRU,
	ARG_SAMPLE,
	ARG_CP_MIN,
	ARG_CP_IVAL,
	ARG_TRI,
	ARG_LOCAL_SEED_CACHE_SZ,
	ARG_CURRENT_SEED_CACHE_SZ,
	ARG_SAM_NO_UNAL,
	ARG_NON_DETERMINISTIC,
	ARG_TEST_25,
	ARG_DESC_KB,
	ARG_DESC_LANDING,
	ARG_DESC_EXP,
	ARG_DESC_PRIORITIZE,
	ARG_DESC_FMOPS,
	ARG_LOG_DP,
	ARG_LOG_DP_OPP,
	ARG_XEQ,
	ARG_THREAD_CEILING,
	ARG_THREAD_PIDDIR,
	ARG_INTERLEAVED,
	ARG_TRIM_TO,
	ARG_PRESERVE_TAGS,
	ARG_ALIGN_PAIRED_READS,
	ARG_SRA_ACC,
	ARG_SAM_APPEND_COMMENT,
};
#endif