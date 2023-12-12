#include "affine_wavefront_align.h"
#include "gap_affine/affine_wavefront_backtrace.h"
#include "gap_affine/affine_wavefront_display.h"
#include "gap_affine/affine_wavefront_extend.h"
#include "gap_affine/affine_wavefront_utils.h"
#include "utils/string_padded.h"
#include "stdio.h"

void affine_wavefronts_fetch_wavefronts(
    affine_wavefronts_t *const affine_wavefronts,
    affine_wavefront_set *const wavefront_set, const int score)
{
  const affine_penalties_t *const wavefront_penalties = &(affine_wavefronts->penalties.wavefront_penalties);
  const int mismatch_score = score - wavefront_penalties->mismatch;
  const int gap_open_score = score - wavefront_penalties->gap_opening - wavefront_penalties->gap_extension;
  const int gap_extend_score = score - wavefront_penalties->gap_extension;
  wavefront_set->in_mwavefront_sub = affine_wavefronts_get_source_mwavefront(affine_wavefronts, mismatch_score);
  wavefront_set->in_mwavefront_gap = affine_wavefronts_get_source_mwavefront(affine_wavefronts, gap_open_score);
  wavefront_set->in_iwavefront_ext = affine_wavefronts_get_source_iwavefront(affine_wavefronts, gap_extend_score);
  wavefront_set->in_dwavefront_ext = affine_wavefronts_get_source_dwavefront(affine_wavefronts, gap_extend_score);
}

void affine_wavefronts_allocate_wavefronts(
    affine_wavefronts_t *const affine_wavefronts,
    affine_wavefront_set *const wavefront_set,
    const int score,
    const int lo_effective, const int hi_effective)
{

  wavefront_set->out_mwavefront =
      affine_wavefronts_allocate_wavefront(affine_wavefronts, lo_effective, hi_effective);
  affine_wavefronts->mwavefronts[score] = wavefront_set->out_mwavefront;
  if (!wavefront_set->in_mwavefront_gap->null || !wavefront_set->in_iwavefront_ext->null)
  {
    wavefront_set->out_iwavefront =
        affine_wavefronts_allocate_wavefront(affine_wavefronts, lo_effective, hi_effective);
    affine_wavefronts->iwavefronts[score] = wavefront_set->out_iwavefront;
  }
  else
  {
    wavefront_set->out_iwavefront = NULL;
  }

  if (!wavefront_set->in_mwavefront_gap->null || !wavefront_set->in_dwavefront_ext->null)
  {
    wavefront_set->out_dwavefront =
        affine_wavefronts_allocate_wavefront(affine_wavefronts, lo_effective, hi_effective);
    affine_wavefronts->dwavefronts[score] = wavefront_set->out_dwavefront;
  }
  else
  {
    wavefront_set->out_dwavefront = NULL;
  }
}

void affine_wavefronts_compute_limits(
    affine_wavefronts_t *const affine_wavefronts,
    const affine_wavefront_set *const wavefront_set,
    const int score, int *const lo_effective,
    int *const hi_effective)
{

  int lo = wavefront_set->in_mwavefront_sub->lo;
  if (lo > wavefront_set->in_mwavefront_gap->lo)
    lo = wavefront_set->in_mwavefront_gap->lo;
  if (lo > wavefront_set->in_iwavefront_ext->lo)
    lo = wavefront_set->in_iwavefront_ext->lo;
  if (lo > wavefront_set->in_dwavefront_ext->lo)
    lo = wavefront_set->in_dwavefront_ext->lo;
  --lo;

  int hi = wavefront_set->in_mwavefront_sub->hi;
  if (hi < wavefront_set->in_mwavefront_gap->hi)
    hi = wavefront_set->in_mwavefront_gap->hi;
  if (hi < wavefront_set->in_iwavefront_ext->hi)
    hi = wavefront_set->in_iwavefront_ext->hi;
  if (hi < wavefront_set->in_dwavefront_ext->hi)
    hi = wavefront_set->in_dwavefront_ext->hi;
  ++hi;

  *hi_effective = hi;
  *lo_effective = lo;
}

#define AFFINE_WAVEFRONT_DECLARE(wavefront, prefix)                \
  const awf_offset_t *const prefix##_offsets = wavefront->offsets; \
  const int prefix##_hi = wavefront->hi;                           \
  const int prefix##_lo = wavefront->lo
#define AFFINE_WAVEFRONT_COND_FETCH(prefix, index, value) \
  (prefix##_lo <= (index) && (index) <= prefix##_hi) ? (value) : AFFINE_WAVEFRONT_OFFSET_NULL

void affine_wavefronts_compute_offsets_idm(
    affine_wavefronts_t *const affine_wavefronts,
    const affine_wavefront_set *const wavefront_set, const int lo,
    const int hi)
{
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub, m_sub);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap, m_gap);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_iwavefront_ext, i_ext);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_dwavefront_ext, d_ext);
  awf_offset_t *const out_ioffsets = wavefront_set->out_iwavefront->offsets;
  awf_offset_t *const out_doffsets = wavefront_set->out_dwavefront->offsets;
  awf_offset_t *const out_moffsets = wavefront_set->out_mwavefront->offsets;
  int min_hi = wavefront_set->in_mwavefront_sub->hi;
  if (!wavefront_set->in_mwavefront_gap->null && min_hi > wavefront_set->in_mwavefront_gap->hi - 1)
    min_hi = wavefront_set->in_mwavefront_gap->hi - 1;
  if (!wavefront_set->in_iwavefront_ext->null && min_hi > wavefront_set->in_iwavefront_ext->hi + 1)
    min_hi = wavefront_set->in_iwavefront_ext->hi + 1;
  if (!wavefront_set->in_dwavefront_ext->null && min_hi > wavefront_set->in_dwavefront_ext->hi - 1)
    min_hi = wavefront_set->in_dwavefront_ext->hi - 1;
  int max_lo = wavefront_set->in_mwavefront_sub->lo;
  if (!wavefront_set->in_mwavefront_gap->null && max_lo < wavefront_set->in_mwavefront_gap->lo + 1)
    max_lo = wavefront_set->in_mwavefront_gap->lo + 1;
  if (!wavefront_set->in_iwavefront_ext->null && max_lo < wavefront_set->in_iwavefront_ext->lo + 1)
    max_lo = wavefront_set->in_iwavefront_ext->lo + 1;
  if (!wavefront_set->in_dwavefront_ext->null && max_lo < wavefront_set->in_dwavefront_ext->lo - 1)
    max_lo = wavefront_set->in_dwavefront_ext->lo - 1;

  int k;
  for (k = lo; k < max_lo; ++k)
  {
    const awf_offset_t ins_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap, k - 1, m_gap_offsets[k - 1]);
    const awf_offset_t ins_i = AFFINE_WAVEFRONT_COND_FETCH(i_ext, k - 1, i_ext_offsets[k - 1]);
    const awf_offset_t ins = MAX(ins_g, ins_i) + 1;
    out_ioffsets[k] = ins;
    const awf_offset_t del_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap, k + 1, m_gap_offsets[k + 1]);
    const awf_offset_t del_d = AFFINE_WAVEFRONT_COND_FETCH(d_ext, k + 1, d_ext_offsets[k + 1]);
    const awf_offset_t del = MAX(del_g, del_d);
    out_doffsets[k] = del;
    const awf_offset_t sub = AFFINE_WAVEFRONT_COND_FETCH(m_sub, k, m_sub_offsets[k] + 1);
    out_moffsets[k] = MAX(del, MAX(sub, ins));
  }

#if defined(__clang__) 
#pragma clang loop vectorize(enable)
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC ivdep 
#else
#pragma ivdep
#endif

  for (k = max_lo; k <= min_hi; ++k)
  {
    const awf_offset_t m_gapi_value = m_gap_offsets[k - 1];
    const awf_offset_t i_ext_value = i_ext_offsets[k - 1];
    const awf_offset_t ins = MAX(m_gapi_value, i_ext_value) + 1;
    out_ioffsets[k] = ins;
    const awf_offset_t m_gapd_value = m_gap_offsets[k + 1];
    const awf_offset_t d_ext_value = d_ext_offsets[k + 1];
    const awf_offset_t del = MAX(m_gapd_value, d_ext_value);
    out_doffsets[k] = del;
    const awf_offset_t sub = m_sub_offsets[k] + 1;
    out_moffsets[k] = MAX(del, MAX(sub, ins));
  }

  for (k = min_hi + 1; k <= hi; ++k)
  {
    const awf_offset_t ins_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap, k - 1, m_gap_offsets[k - 1]);
    const awf_offset_t ins_i = AFFINE_WAVEFRONT_COND_FETCH(i_ext, k - 1, i_ext_offsets[k - 1]);
    const awf_offset_t ins = MAX(ins_g, ins_i) + 1;
    out_ioffsets[k] = ins;
    const awf_offset_t del_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap, k + 1, m_gap_offsets[k + 1]);
    const awf_offset_t del_d = AFFINE_WAVEFRONT_COND_FETCH(d_ext, k + 1, d_ext_offsets[k + 1]);
    const awf_offset_t del = MAX(del_g, del_d);
    out_doffsets[k] = del;
    const awf_offset_t sub = AFFINE_WAVEFRONT_COND_FETCH(m_sub, k, m_sub_offsets[k] + 1);
    out_moffsets[k] = MAX(del, MAX(sub, ins));
  }
}

void affine_wavefronts_compute_offsets_im(
    affine_wavefronts_t *const affine_wavefronts,
    const affine_wavefront_set *const wavefront_set,
    const int lo,
    const int hi)
{
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub, m_sub);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap, m_gap);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_iwavefront_ext, i_ext);
  awf_offset_t *const out_ioffsets = wavefront_set->out_iwavefront->offsets;
  awf_offset_t *const out_moffsets = wavefront_set->out_mwavefront->offsets;
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
#pragma GCC ivdep
#else
#pragma ivdep
#endif
  for (k = lo; k <= hi; ++k)
  {
    const awf_offset_t ins_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap, k - 1, m_gap_offsets[k - 1]);
    const awf_offset_t ins_i = AFFINE_WAVEFRONT_COND_FETCH(i_ext, k - 1, i_ext_offsets[k - 1]);
    const awf_offset_t ins = MAX(ins_g, ins_i) + 1;
    out_ioffsets[k] = ins;
    const awf_offset_t sub = AFFINE_WAVEFRONT_COND_FETCH(m_sub, k, m_sub_offsets[k] + 1);
    out_moffsets[k] = MAX(ins, sub);
  }
}
void affine_wavefronts_compute_offsets_dm(
    affine_wavefronts_t *const affine_wavefronts,
    const affine_wavefront_set *const wavefront_set,
    const int lo,
    const int hi)
{
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub, m_sub);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap, m_gap);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_dwavefront_ext, d_ext);
  awf_offset_t *const out_doffsets = wavefront_set->out_dwavefront->offsets;
  awf_offset_t *const out_moffsets = wavefront_set->out_mwavefront->offsets;
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
#pragma GCC ivdep
#else
#pragma ivdep
#endif
  for (k = lo; k <= hi; ++k)
  {
    const awf_offset_t del_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap, k + 1, m_gap_offsets[k + 1]);
    const awf_offset_t del_d = AFFINE_WAVEFRONT_COND_FETCH(d_ext, k + 1, d_ext_offsets[k + 1]);
    const awf_offset_t del = MAX(del_g, del_d);
    out_doffsets[k] = del;
    const awf_offset_t sub = AFFINE_WAVEFRONT_COND_FETCH(m_sub, k, m_sub_offsets[k] + 1);
    out_moffsets[k] = MAX(del, sub);
  }
}

void affine_wavefronts_compute_offsets_m(
    affine_wavefronts_t *const affine_wavefronts,
    const affine_wavefront_set *const wavefront_set,
    const int lo,
    const int hi)
{
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub, m_sub);
  awf_offset_t *const out_moffsets = wavefront_set->out_mwavefront->offsets;
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
#pragma GCC ivdep
#else
#pragma ivdep
#endif
  for (k = lo; k <= hi; ++k)
  {
    out_moffsets[k] = AFFINE_WAVEFRONT_COND_FETCH(m_sub, k, m_sub_offsets[k] + 1);
  }
}

void affine_wavefronts_compute_wavefront(
    affine_wavefronts_t *const affine_wavefronts,
    const char *const pattern,
    const int pattern_length,
    const char *const text,
    const int text_length,
    const int score)
{
  affine_wavefront_set wavefront_set;
  affine_wavefronts_fetch_wavefronts(affine_wavefronts, &wavefront_set, score);
  if (wavefront_set.in_mwavefront_sub->null &&
      wavefront_set.in_mwavefront_gap->null &&
      wavefront_set.in_iwavefront_ext->null &&
      wavefront_set.in_dwavefront_ext->null)
  {
    WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_steps_null, 1);
    return;
  }
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_null_used, (wavefront_set.in_mwavefront_sub->null ? 1 : 0));
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_null_used, (wavefront_set.in_mwavefront_gap->null ? 1 : 0));
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_null_used, (wavefront_set.in_iwavefront_ext->null ? 1 : 0));
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_null_used, (wavefront_set.in_dwavefront_ext->null ? 1 : 0));

  int hi, lo;
  affine_wavefronts_compute_limits(affine_wavefronts, &wavefront_set, score, &lo, &hi);
  affine_wavefronts_allocate_wavefronts(affine_wavefronts, &wavefront_set, score, lo, hi);

  const int kernel = ((wavefront_set.out_iwavefront != NULL) << 1) | (wavefront_set.out_dwavefront != NULL); 
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_compute_kernel[kernel], 1);

  switch (kernel)
  {
  case 3:
    affine_wavefronts_compute_offsets_idm(affine_wavefronts, &wavefront_set, lo, hi);
    break;
  case 2:
    affine_wavefronts_compute_offsets_im(affine_wavefronts, &wavefront_set, lo, hi);
    break;
  case 1:
    affine_wavefronts_compute_offsets_dm(affine_wavefronts, &wavefront_set, lo, hi);
    break;
  case 0:
    affine_wavefronts_compute_offsets_m(affine_wavefronts, &wavefront_set, lo, hi);
    break;
  }

  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts, wf_operations, hi - lo + 1);
#ifdef AFFINE_WAVEFRONT_DEBUG
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront != NULL)
  {
    int k;
    for (k = mwavefront->lo; k <= mwavefront->hi; ++k)
    {
      mwavefront->offsets_base[k] = mwavefront->offsets[k];
    }
  }
#endif
}

void affine_wavefronts_align(
    affine_wavefronts_t *const affine_wavefronts, const char *const pattern,
    const int pattern_length,
    const char *const text,
    const int text_length)
{
  strings_padded_t *const strings_padded =
      strings_padded_new_rhomb(pattern, pattern_length, text, text_length,
                               AFFINE_WAVEFRONT_PADDING, affine_wavefronts->mm_allocator);

  affine_wavefront_initialize(affine_wavefronts);

  int score = 0;
  int alignment_k = 0;
  int done = false;
  while (true)
  {
    done = affine_wavefronts_extend_wavefront_packed(
        affine_wavefronts, strings_padded->pattern_padded, pattern_length,
        strings_padded->text_padded, text_length, score, &alignment_k);
    if (done)
    {
      affine_wavefronts_backtrace(
          affine_wavefronts, strings_padded->pattern_padded, pattern_length,
          strings_padded->text_padded, text_length, score, alignment_k);
      break;
    }
    ++score;
    affine_wavefronts_compute_wavefront(affine_wavefronts, strings_padded->pattern_padded, pattern_length,
                                        strings_padded->text_padded, text_length, score);
  }

  strings_padded_delete(strings_padded);
}
