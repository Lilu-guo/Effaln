#include "gap_affine/affine_wavefront_display.h"
#include "gap_affine/affine_wavefront_extend.h"
#include "gap_affine/affine_wavefront_reduction.h"
#include "gap_affine/affine_wavefront_utils.h"
#include "utils/string_padded.h"

void affine_wavefronts_reduce_wavefront_offsets(
    affine_wavefronts_t *const affine_wavefronts,
    affine_wavefront_t *const wavefront,
    const int pattern_length,
    const int text_length,
    const int min_distance,
    const int max_distance_threshold,
    const int alignment_k)
{
  const awf_offset_t *const offsets = wavefront->offsets;
  int k;
  const int top_limit = MIN(alignment_k - 1, wavefront->hi);
  for (k = wavefront->lo; k < top_limit; ++k)
  {
    const int distance = affine_wavefronts_compute_distance(pattern_length, text_length, offsets[k], k);
    if (distance - min_distance <= max_distance_threshold)
      break;
    ++(wavefront->lo);
  }
  const int botton_limit = MAX(alignment_k + 1, wavefront->lo);
  for (k = wavefront->hi; k > botton_limit; --k)
  {
    const int distance = affine_wavefronts_compute_distance(pattern_length, text_length, offsets[k], k);
    if (distance - min_distance <= max_distance_threshold)
      break;
    --(wavefront->hi);
  }
  if (wavefront->lo > wavefront->hi)
  {
    wavefront->null = true;
  }
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_reduced_cells,
                              (wavefront->hi_base - wavefront->hi) + (wavefront->lo - wavefront->lo_base));
}

void affine_wavefronts_reduce_wavefronts(
    affine_wavefronts_t *const affine_wavefronts,
    const int pattern_length,
    const int text_length,
    const int score)
{
  const int min_wavefront_length = affine_wavefronts->reduction.min_wavefront_length;
  const int max_distance_threshold = affine_wavefronts->reduction.max_distance_threshold;
  const int alignment_k = abs(AFFINE_WAVEFRONT_DIAGONAL(text_length, pattern_length));

  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront == NULL)
    return;
  if ((mwavefront->hi - mwavefront->lo + 1) < min_wavefront_length)
    return;
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_reduction, 1);
  const awf_offset_t *const offsets = mwavefront->offsets;
  int min_distance = MAX(pattern_length, text_length);
  int k;
  for (k = mwavefront->lo; k <= mwavefront->hi; ++k)
  {
    const int distance = affine_wavefronts_compute_distance(pattern_length, text_length, offsets[k], k);
    min_distance = MIN(min_distance, distance);
  }
  affine_wavefronts_reduce_wavefront_offsets(
      affine_wavefronts, mwavefront, pattern_length, text_length,
      min_distance, max_distance_threshold, alignment_k);
  affine_wavefront_t *const iwavefront = affine_wavefronts->iwavefronts[score];
  if (iwavefront != NULL)
  {
    if (mwavefront->lo > iwavefront->lo)
      iwavefront->lo = mwavefront->lo;
    if (mwavefront->hi < iwavefront->hi)
      iwavefront->hi = mwavefront->hi;
    if (iwavefront->lo > iwavefront->hi)
      iwavefront->null = true;
  }
  affine_wavefront_t *const dwavefront = affine_wavefronts->dwavefronts[score];
  if (dwavefront != NULL)
  {
    if (mwavefront->lo > dwavefront->lo)
      dwavefront->lo = mwavefront->lo;
    if (mwavefront->hi < dwavefront->hi)
      dwavefront->hi = mwavefront->hi;
    if (dwavefront->lo > dwavefront->hi)
      dwavefront->null = true;
  }
}

void affine_wavefronts_extend_mwavefront_epiloge(affine_wavefronts_t *const affine_wavefronts, const int score, const int pattern_length,
                                                 const int text_length)
{
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_extensions,
                              affine_wavefronts->mwavefronts[score]->hi - affine_wavefronts->mwavefronts[score]->lo + 1);
#ifdef AFFINE_WAVEFRONT_DEBUG
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  int k;
  for (k = mwavefront->lo; k <= mwavefront->hi; ++k)
  {
    if (mwavefront->offsets[k] >= 0)
    {
      awf_offset_t offset;
      for (offset = mwavefront->offsets_base[k]; offset <= mwavefront->offsets[k]; ++offset)
      {
        affine_wavefronts_set_edit_table(affine_wavefronts, pattern_length, text_length, k, offset, score);
      }
    }
  }
#endif
}

bool affine_wavefronts_extend_mwavefront_compute_packed(affine_wavefronts_t *const affine_wavefronts,
                                                        const char *const pattern, const int pattern_length,
                                                        const char *const text, const int text_length,
                                                        const int score,
                                                        int *alignment_k)
{
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront == NULL)
    return false;
  awf_offset_t *const offsets = mwavefront->offsets;
  int k;
  for (k = mwavefront->lo; k <= mwavefront->hi; ++k)
  {
    const awf_offset_t offset = offsets[k];
    const int v = AFFINE_WAVEFRONT_V(k, offset);
    const int h = AFFINE_WAVEFRONT_H(k, offset);
    uint64_t *pattern_blocks = (uint64_t *)(pattern + v);
    uint64_t *text_blocks = (uint64_t *)(text + h);
    uint64_t pattern_block = *pattern_blocks;
    uint64_t text_block = *text_blocks;
    uint64_t cmp = pattern_block ^ text_block;
    while (__builtin_expect(!cmp, 0))
    {
      offsets[k] += 8;
      ++pattern_blocks;
      ++text_blocks;
      pattern_block = *pattern_blocks;
      text_block = *text_blocks;
      cmp = pattern_block ^ text_block;
      WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_extend_inner_loop, 1);
    }
    const int equal_right_bits = __builtin_ctzl(cmp);
    const int equal_chars = DIV_FLOOR(equal_right_bits, 8);
    offsets[k] += equal_chars;
    if (offsets[k] >= text_length)
    {
      *alignment_k = k;
      return true;
    }
  }
  affine_wavefronts_extend_mwavefront_epiloge(affine_wavefronts, score, pattern_length, text_length);
  return false;
}

void affine_wavefronts_extend_mwavefront_compute(
    affine_wavefronts_t *const affine_wavefronts,
    const char *const pattern,
    const int pattern_length,
    const char *const text,
    const int text_length,
    const int score)
{
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront == NULL)
    return;
  awf_offset_t *const offsets = mwavefront->offsets;
  int k;
  for (k = mwavefront->lo; k <= mwavefront->hi; ++k)
  {
    const awf_offset_t offset = offsets[k];
    int v = AFFINE_WAVEFRONT_V(k, offset);
    int h = AFFINE_WAVEFRONT_H(k, offset);
    while (pattern[v++] == text[h++])
    {
      ++(offsets[k]);
    }
  }
  affine_wavefronts_extend_mwavefront_epiloge(affine_wavefronts, score, pattern_length, text_length);
}

bool affine_wavefronts_extend_wavefront_packed(
    affine_wavefronts_t *const affine_wavefronts,
    const char *const pattern,
    const int pattern_length,
    const char *const text,
    const int text_length,
    const int score,
    int *alignment_k)
{
  bool done = false;
  done = affine_wavefronts_extend_mwavefront_compute_packed(
      affine_wavefronts, pattern, pattern_length,
      text, text_length, score, alignment_k);

  if (affine_wavefronts->reduction.reduction_strategy == wavefronts_reduction_dynamic)
  {
    affine_wavefronts_reduce_wavefronts(
        affine_wavefronts, pattern_length,
        text_length, score);
  }
  return done;
}
