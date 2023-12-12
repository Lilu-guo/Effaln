#ifndef SSE_WRAP_H_
#define SSE_WRAP_H_
#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
#include "simde/x86/sse2.h"
#else
#include <emmintrin.h>
#endif
#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
typedef simde__m128i __m128i;
#define _mm_adds_epi16(x, y) simde_mm_adds_epi16(x, y)
#define _mm_adds_epu8(x, y) simde_mm_adds_epu8(x, y)
#define _mm_cmpeq_epi16(x, y) simde_mm_cmpeq_epi16(x, y)
#define _mm_cmpeq_epi8(x, y) simde_mm_cmpeq_epi8(x, y)
#define _mm_cmpgt_epi16(x, y) simde_mm_cmpgt_epi16(x, y)
#define _mm_cmpgt_epi8(x, y) simde_mm_cmpgt_epi8(x, y)
#define _mm_cmplt_epi16(x, y) simde_mm_cmplt_epi16(x, y)
#define _mm_cmplt_epu8(x, y) simde_mm_cmplt_epu8(x, y)
#define _mm_extract_epi16(x, y) simde_mm_extract_epi16(x, y)
#define _mm_insert_epi16(x, y, z) simde_mm_insert_epi16(x, y, z)
#define _mm_load_si128(x) simde_mm_load_si128(x)
#define _mm_max_epi16(x, y) simde_mm_max_epi16(x, y)
#define _mm_max_epu8(x, y) simde_mm_max_epu8(x, y)
#define _mm_movemask_epi8(x) simde_mm_movemask_epi8(x)
#define _mm_or_si128(x, y) simde_mm_or_si128(x, y)
#define _mm_setzero_si128() simde_mm_setzero_si128()
#define _mm_shuffle_epi32(x, y) simde_mm_shuffle_epi32(x, y)
#define _mm_shufflelo_epi16(x, y) simde_mm_shufflelo_epi16(x, y)
#define _mm_slli_epi16(x, y) simde_mm_slli_epi16(x, y)
#define _mm_slli_si128(x, y) simde_mm_slli_si128(x, y)
#define _mm_srli_epi16(x, y) simde_mm_srli_epi16(x, y)
#define _mm_srli_epu8(x, y) simde_mm_srli_epu8(x, y)
#define _mm_srli_si128(x, y) simde_mm_srli_si128(x, y)
#define _mm_store_si128(x, y) simde_mm_store_si128(x, y)
#define _mm_subs_epi16(x, y) simde_mm_subs_epi16(x, y)
#define _mm_subs_epu8(x, y) simde_mm_subs_epu8(x, y)
#define _mm_xor_si128(x, y) simde_mm_xor_si128(x, y)
#endif
#endif
