#pragma once
#define MOD ((1UL << 29) - 1)
#define MOD ((1UL << 29) - 1)
#define BATCH_SIZE 15000
#define MAX_INDEL 15
#define SC_MCH 2
#define SC_MIS 8
#define GAPO 12
#define GAPE 2
#define SC_AMBI 1
#define END_BONUS 10
#define Z_DROP 100
#define BANDWIDTH 151
#define SHORT_READ_MODE 0
#define EMBED_PAD 4
#define NUM_STR 1
#define NUM_CHAR 5
#define MAX_ELEN 2000 
#define RBITS_PER_STRING (MAX_ELEN * NUM_CHAR)
#define TOTAL_RBITS (RBITS_PER_STRING * NUM_STR)
#define BITPOS(STR_ID, OFFSET, CHAR_ID) (STR_ID * RBITS_PER_STRING + OFFSET + CHAR_ID) 
#define MIS_PENALTY -1
#define MAX_LEN 1000 
#define KSW_EZ_SCORE_ONLY 0x01
#define KSW_EZ_RIGHT 0x02
#define KSW_EZ_GENERIC_SC 0x04
#define KSW_EZ_APPROX_MAX 0x08
#define KSW_EZ_APPROX_DROP 0x10
#define KSW_EZ_EXTZ_ONLY 0x40
#define KSW_EZ_REV_CIGAR 0x80
#ifndef kroundup32
#define kroundup32(x) (--(x), (x) |= (x) >> 1, (x) |= (x) >> 2, (x) |= (x) >> 4, (x) |= (x) >> 8, (x) |= (x) >> 16, ++(x))
#endif