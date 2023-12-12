#ifndef AFFINE_WAVEFRONT_H_
#define AFFINE_WAVEFRONT_H_

#include "../utils/commons.h"
#include "../system/profiler_counter.h"
#include "../system/profiler_timer.h"

#include "affine_table.h"
#include "affine_wavefront_penalties.h"
#include "affine_wavefront_reduction.h"
#include "wavefront_stats.h"

#define AFFINE_WAVEFRONT_OFFSET_NULL (-10)
#define AFFINE_WAVEFRONT_K_NULL (INT_MAX / 2)

#define AFFINE_WAVEFRONT_V(k, offset) ((offset) - (k))
#define AFFINE_WAVEFRONT_H(k, offset) (offset)

#define AFFINE_WAVEFRONT_DIAGONAL(h, v) ((v) - (h))
#define AFFINE_WAVEFRONT_OFFSET(h, v) (h)

#define AFFINE_WAVEFRONT_W32

#ifdef AFFINE_WAVEFRONT_W8
typedef int8_t awf_offset_t;
#else
#ifdef AFFINE_WAVEFRONT_W16
    typedef int16_t awf_offset_t;
#else 
    typedef int32_t awf_offset_t;
#endif
#endif

typedef struct
{
  bool null;
  int lo;
  int hi;
  int lo_base;
  int hi_base;
  awf_offset_t *offsets; 
#ifdef AFFINE_WAVEFRONT_DEBUG
  awf_offset_t *offsets_base;
#endif
} affine_wavefront_t;

typedef struct
{
  int pattern_length;
  int text_length;
  int num_wavefronts;
  int max_penalty;
  int max_k;
  int min_k;
  affine_wavefront_t **mwavefronts; 
  affine_wavefront_t **iwavefronts;
  affine_wavefront_t **dwavefronts;
  affine_wavefront_t wavefront_null;
  affine_wavefronts_reduction_t reduction;
  affine_wavefronts_penalties_t penalties;
  edit_cigar_t edit_cigar;
  mm_allocator_t *mm_allocator;
  affine_wavefront_t *wavefronts_mem;
  affine_wavefront_t *wavefronts_current;
  wavefronts_stats_t *wavefronts_stats;
#ifdef AFFINE_WAVEFRONT_DEBUG
      affine_table_t gap_affine_table;
#endif
} affine_wavefronts_t;

typedef struct
{

  affine_wavefront_t *in_mwavefront_sub;
  affine_wavefront_t *in_mwavefront_gap;
  affine_wavefront_t *in_iwavefront_ext;
  affine_wavefront_t *in_dwavefront_ext;

  affine_wavefront_t *out_mwavefront;
  affine_wavefront_t *out_iwavefront;
  affine_wavefront_t *out_dwavefront;
} affine_wavefront_set;

void affine_wavefronts_clear(
    affine_wavefronts_t *const affine_wavefronts);
void affine_wavefronts_delete(
    affine_wavefronts_t *const affine_wavefronts);

affine_wavefronts_t *affine_wavefronts_new_complete(
    const int pattern_length,
    const int text_length,
    affine_penalties_t *const penalties,
    wavefronts_stats_t *const wavefronts_stats,
    mm_allocator_t *const mm_allocator);
affine_wavefronts_t *affine_wavefronts_new_reduced(
    const int pattern_length,
    const int text_length,
    affine_penalties_t *const penalties,
    const int min_wavefront_length,
    const int max_distance_threshold,
    wavefronts_stats_t *const wavefronts_stats,
    mm_allocator_t *const mm_allocator);

affine_wavefront_t *affine_wavefronts_allocate_wavefront(
    affine_wavefronts_t *const affine_wavefronts,
    const int lo_base,
    const int hi_base);

#endif
