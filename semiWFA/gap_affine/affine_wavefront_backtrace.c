#include "gap_affine/affine_wavefront_backtrace.h"

bool affine_wavefronts_valid_location(
    const int k,
    const awf_offset_t offset,
    const int pattern_length,
    const int text_length)
{
  const int v = AFFINE_WAVEFRONT_V(k, offset);
  const int h = AFFINE_WAVEFRONT_H(k, offset);
  return (v > 0 && v <= pattern_length &&
          h > 0 && h <= text_length);
}

void affine_wavefronts_offset_add_trailing_gap(
    edit_cigar_t *const edit_cigar,
    const int k, const int alignment_k)
{
  char *const operations = edit_cigar->operations;
  int op_sentinel = edit_cigar->begin_offset;
  int i;
  if (k < alignment_k)
  {
    for (i = k; i < alignment_k; ++i)
      operations[op_sentinel--] = 'I';
  }
  else if (k > alignment_k)
  {
    for (i = alignment_k; i < k; ++i)
      operations[op_sentinel--] = 'D';
  }
  edit_cigar->begin_offset = op_sentinel;
}

awf_offset_t backtrace_wavefront_trace_deletion_open_offset(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k + 1 &&
      k + 1 <= mwavefront->hi_base)
  {
    return mwavefront->offsets[k + 1];
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}

awf_offset_t backtrace_wavefront_trace_deletion_extend_offset(
    affine_wavefronts_t *const affine_wavefronts,
    const int score, const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t *const dwavefront = affine_wavefronts->dwavefronts[score];
  if (dwavefront != NULL &&
      dwavefront->lo_base <= k + 1 &&
      k + 1 <= dwavefront->hi_base)
  {
    return dwavefront->offsets[k + 1];
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}

awf_offset_t backtrace_wavefront_trace_insertion_open_offset(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k - 1 &&
      k - 1 <= mwavefront->hi_base)
  {
    return mwavefront->offsets[k - 1] + 1;
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}

awf_offset_t backtrace_wavefront_trace_insertion_extend_offset(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t *const iwavefront = affine_wavefronts->iwavefronts[score];
  if (iwavefront != NULL &&
      iwavefront->lo_base <= k - 1 &&
      k - 1 <= iwavefront->hi_base)
  {
    return iwavefront->offsets[k - 1] + 1;
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}

awf_offset_t backtrace_wavefront_trace_mismatch_offset(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k &&
      k <= mwavefront->hi_base)
  {
    return mwavefront->offsets[k] + 1;
  }
  else
  {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}

bool backtrace_wavefront_trace_deletion(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return false;
  affine_wavefront_t *const dwavefront = affine_wavefronts->dwavefronts[score];
  return (dwavefront != NULL &&
          dwavefront->lo_base <= k &&
          k <= dwavefront->hi_base &&
          offset == dwavefront->offsets[k]);
}
bool backtrace_wavefront_trace_deletion_open(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return false;
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  return (mwavefront != NULL &&
          mwavefront->lo_base <= k + 1 &&
          k + 1 <= mwavefront->hi_base &&
          offset == mwavefront->offsets[k + 1]);
}
bool backtrace_wavefront_trace_deletion_extend(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return false;
  affine_wavefront_t *const dwavefront = affine_wavefronts->dwavefronts[score];
  return (dwavefront != NULL &&
          dwavefront->lo_base <= k + 1 &&
          k + 1 <= dwavefront->hi_base &&
          offset == dwavefront->offsets[k + 1]);
}
bool backtrace_wavefront_trace_insertion(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return false;
  affine_wavefront_t *const iwavefront = affine_wavefronts->iwavefronts[score];
  return (iwavefront != NULL &&
          iwavefront->lo_base <= k &&
          k <= iwavefront->hi_base &&
          offset == iwavefront->offsets[k]);
}
bool backtrace_wavefront_trace_insertion_open(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return false;
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  return (mwavefront != NULL &&
          mwavefront->lo_base <= k - 1 &&
          k - 1 <= mwavefront->hi_base &&
          offset == mwavefront->offsets[k - 1] + 1);
}
bool backtrace_wavefront_trace_insertion_extend(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return false;
  affine_wavefront_t *const iwavefront = affine_wavefronts->iwavefronts[score];
  return (iwavefront != NULL &&
          iwavefront->lo_base <= k - 1 &&
          k - 1 <= iwavefront->hi_base &&
          offset == iwavefront->offsets[k - 1] + 1);
}
bool backtrace_wavefront_trace_mismatch(
    affine_wavefronts_t *const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset)
{
  if (score < 0)
    return false;
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  return (mwavefront != NULL &&
          mwavefront->lo_base <= k &&
          k <= mwavefront->hi_base &&
          offset == mwavefront->offsets[k] + 1);
}

int affine_wavefronts_backtrace_compute_max_matches(
    affine_wavefronts_t *const affine_wavefronts,
    const char *const pattern,
    const char *const text,
    const int k,
    awf_offset_t offset)
{
  int v = AFFINE_WAVEFRONT_V(k, offset);
  int h = AFFINE_WAVEFRONT_H(k, offset);
  int num_matches = 0;
  while (v > 0 && h > 0 && pattern[--v] == text[--h])
    ++num_matches;
  return num_matches;
}

void affine_wavefronts_backtrace_matches__check(
    affine_wavefronts_t *const affine_wavefronts,
    const char *const pattern,
    const char *const text,
    const int k,
    awf_offset_t offset,
    const bool valid_location, const int num_matches, edit_cigar_t *const edit_cigar)
{
  int i;
  for (i = 0; i < num_matches; ++i)
  {
#ifdef AFFINE_WAVEFRONT_DEBUG

    const int v = AFFINE_WAVEFRONT_V(k, offset);
    const int h = AFFINE_WAVEFRONT_H(k, offset);
    if (!valid_location)
    {
      fprintf(stderr, "Backtrace error: Match outside DP-Table\n");
      exit(1);
    }
    else if (pattern[v - 1] != text[h - 1])
    {
      fprintf(stderr, "Backtrace error: Not a match traceback\n");
      exit(1);
    }
#endif

    edit_cigar->operations[(edit_cigar->begin_offset)--] = 'M';
    --offset;
  }
}

void affine_wavefronts_backtrace_matches(
    edit_cigar_t *const edit_cigar,
    const int num_matches)
{
  int i;
  for (i = 0; i < num_matches; ++i)
  {
    edit_cigar->operations[(edit_cigar->begin_offset)--] = 'M';
  }
}

void affine_wavefronts_backtrace(
    affine_wavefronts_t *const affine_wavefronts,
    char *const pattern,
    const int pattern_length,
    char *const text,
    const int text_length,
    const int alignment_score,
    int alignment_k)
{
  const affine_penalties_t *const wavefront_penalties =
      &(affine_wavefronts->penalties.wavefront_penalties);
  edit_cigar_t *const cigar = &affine_wavefronts->edit_cigar;
  int score = alignment_score;
  int k = alignment_k;
  awf_offset_t offset = affine_wavefronts->mwavefronts[alignment_score]->offsets[k];

  bool valid_location = affine_wavefronts_valid_location(k, offset, pattern_length, text_length);

  backtrace_wavefront_type backtrace_type = backtrace_wavefront_M;
  int v = AFFINE_WAVEFRONT_V(k, offset);
  int h = AFFINE_WAVEFRONT_H(k, offset);
  while (v > 0 && h > 0 && score > 0)
  {

    if (!valid_location)
    {
      valid_location = affine_wavefronts_valid_location(k, offset, pattern_length, text_length);
      if (valid_location)
      {
        affine_wavefronts_offset_add_trailing_gap(cigar, k, alignment_k);
      }
    }

    const int gap_open_score = score - wavefront_penalties->gap_opening - wavefront_penalties->gap_extension;
    const int gap_extend_score = score - wavefront_penalties->gap_extension;
    const int mismatch_score = score - wavefront_penalties->mismatch;
    const awf_offset_t del_ext = (backtrace_type == backtrace_wavefront_I) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_deletion_extend_offset(affine_wavefronts, gap_extend_score, k, offset);
    const awf_offset_t del_open = (backtrace_type == backtrace_wavefront_I) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_deletion_open_offset(affine_wavefronts, gap_open_score, k, offset);
    const awf_offset_t ins_ext = (backtrace_type == backtrace_wavefront_D) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_insertion_extend_offset(affine_wavefronts, gap_extend_score, k, offset);
    const awf_offset_t ins_open = (backtrace_type == backtrace_wavefront_D) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_insertion_open_offset(affine_wavefronts, gap_open_score, k, offset);
    const awf_offset_t misms = (backtrace_type != backtrace_wavefront_M) ? AFFINE_WAVEFRONT_OFFSET_NULL : backtrace_wavefront_trace_mismatch_offset(affine_wavefronts, mismatch_score, k, offset);
    const awf_offset_t max_del = MAX(del_ext, del_open);
    const awf_offset_t max_ins = MAX(ins_ext, ins_open);
    const awf_offset_t max_all = MAX(misms, MAX(max_ins, max_del));
    if (backtrace_type == backtrace_wavefront_M)
    {
      const int num_matches = offset - max_all;
      affine_wavefronts_backtrace_matches__check(affine_wavefronts,
                                                 pattern, text, k, offset, valid_location, num_matches, cigar);
      offset = max_all;
    }

    if (max_all == del_ext)
    {
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'D';
      score = gap_extend_score;
      ++k;
      backtrace_type = backtrace_wavefront_D;
    }
    else if (max_all == del_open)
    {
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'D';
      score = gap_open_score;
      ++k;
      backtrace_type = backtrace_wavefront_M;
    }
    else if (max_all == ins_ext)
    {
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'I';
      score = gap_extend_score;
      --k;
      --offset;
      backtrace_type = backtrace_wavefront_I;
    }
    else if (max_all == ins_open)
    {
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'I';
      score = gap_open_score;
      --k;
      --offset;
      backtrace_type = backtrace_wavefront_M;
    }
    else if (max_all == misms)
    {
      if (valid_location)
        cigar->operations[(cigar->begin_offset)--] = 'X';
      score = mismatch_score;
      --offset;
    }
    else
    {
      fprintf(stderr, "Backtrace error: No link found during backtrace\n");
      exit(1);
    }

    v = AFFINE_WAVEFRONT_V(k, offset);
    h = AFFINE_WAVEFRONT_H(k, offset);
  }
  if (score == 0)
  {
    affine_wavefronts_backtrace_matches__check(affine_wavefronts,
                                               pattern, text, k, offset, valid_location, offset, cigar);
  }
  else
  {
    while (v > 0)
    {
      cigar->operations[(cigar->begin_offset)--] = 'D';
      --v;
    };
    while (h > 0)
    {
      cigar->operations[(cigar->begin_offset)--] = 'I';
      --h;
    };
  }

  ++(cigar->begin_offset);
}
