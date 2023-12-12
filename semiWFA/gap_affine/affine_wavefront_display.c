#include "gap_affine/affine_wavefront_display.h"
#include "gap_affine/affine_wavefront_utils.h"
void affine_wavefronts_set_edit_table(
    affine_wavefronts_t *const affine_wavefronts,
    const int pattern_length,
    const int text_length,
    const int k,
    const awf_offset_t offset, const int score)
{
#ifdef AFFINE_WAVEFRONT_DEBUG
  if (offset < 0)
    return;
  const int h = AFFINE_WAVEFRONT_H(k, offset);
  const int v = AFFINE_WAVEFRONT_V(k, offset);
  if (0 <= v && v <= pattern_length && 0 <= h && h <= text_length)
  {
    if (affine_wavefronts->gap_affine_table.columns[h][v].M == -1)
    {
      affine_wavefronts->gap_affine_table.columns[h][v].M = score;
    }
  }
#endif
}
#define AFFINE_WAVEFRONTS_PRINT_ELEMENT(wavefront, k)                \
                                                                     \
  if (wavefront != NULL && wavefront->lo <= k && k <= wavefront->hi) \
  {                                                                  \
    if (wavefront->offsets[k] >= 0)                                  \
    {                                                                \
      fprintf(stream, "[%2d]", (int)wavefront->offsets[k]);          \
    }                                                                \
    else                                                             \
    {                                                                \
      fprintf(stream, "[  ]");                                       \
    }                                                                \
  }                                                                  \
  else                                                               \
  {                                                                  \
    fprintf(stream, "    ");                                         \
  }
void affine_wavefronts_print_wavefront(
    FILE *const stream,
    affine_wavefronts_t *const affine_wavefronts,
    const int current_score)
{
  fprintf(stream, ">[SCORE=%3d]\n", current_score);
  fprintf(stream, "        ");
  if (affine_wavefronts->iwavefronts != NULL)
  {
    fprintf(stream, " [M] [I] [D] ");
  }
  else
  {
    fprintf(stream, " [M] ");
  }
  fprintf(stream, "\n");
  int k, max_k = 0, min_k = 0;
  for (k = affine_wavefronts->max_k; k >= affine_wavefronts->min_k; k--)
  {
    affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[current_score];
    if (mwavefront != NULL)
    {
      max_k = MAX(max_k, mwavefront->hi);
      min_k = MIN(min_k, mwavefront->lo);
    }
    if (affine_wavefronts->iwavefronts != NULL)
    {
      affine_wavefront_t *const iwavefront = affine_wavefronts->iwavefronts[current_score];
      affine_wavefront_t *const dwavefront = affine_wavefronts->dwavefronts[current_score];
      if (iwavefront != NULL)
      {
        max_k = MAX(max_k, iwavefront->hi);
        min_k = MIN(min_k, iwavefront->lo);
      }
      if (dwavefront != NULL)
      {
        max_k = MAX(max_k, dwavefront->hi);
        min_k = MIN(min_k, dwavefront->lo);
      }
    }
  }
  for (k = max_k; k >= min_k; k--)
  {
    fprintf(stream, "[k=%3d] ", k);
    affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[current_score];
    if (affine_wavefronts->iwavefronts != NULL)
    {
      affine_wavefront_t *const iwavefront = affine_wavefronts->iwavefronts[current_score];
      affine_wavefront_t *const dwavefront = affine_wavefronts->dwavefronts[current_score];
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(mwavefront, k);
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(iwavefront, k);
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(dwavefront, k);
    }
    else
    {
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(mwavefront, k);
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}
void affine_wavefronts_print_wavefronts(
    FILE *const stream,
    affine_wavefronts_t *const affine_wavefronts,
    const int current_score)
{
  int k, s;
  fprintf(stream, ">[SCORE=%3d]\n", current_score);
  fprintf(stream, "        ");
  for (s = 0; s <= current_score; ++s)
  {
    if (affine_wavefronts->iwavefronts != NULL)
    {
      fprintf(stream, " [M] [I] [D] ");
    }
    else
    {
      fprintf(stream, " [M] ");
    }
  }
  fprintf(stream, "\n");
  int max_k = 0, min_k = 0;
  for (k = affine_wavefronts->max_k; k >= affine_wavefronts->min_k; k--)
  {
    for (s = 0; s <= current_score; ++s)
    {
      affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[s];
      if (mwavefront != NULL)
      {
        max_k = MAX(max_k, mwavefront->hi);
        min_k = MIN(min_k, mwavefront->lo);
      }
      affine_wavefront_t *const iwavefront = affine_wavefronts->iwavefronts[s];
      if (iwavefront != NULL)
      {
        max_k = MAX(max_k, iwavefront->hi);
        min_k = MIN(min_k, iwavefront->lo);
      }
      affine_wavefront_t *const dwavefront = affine_wavefronts->dwavefronts[s];
      if (dwavefront != NULL)
      {
        max_k = MAX(max_k, dwavefront->hi);
        min_k = MIN(min_k, dwavefront->lo);
      }
    }
  }
  for (k = max_k; k >= min_k; k--)
  {
    fprintf(stream, "[k=%3d] ", k);
    for (s = 0; s <= current_score; ++s)
    {
      affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[s];
      affine_wavefront_t *const iwavefront = affine_wavefronts->iwavefronts[s];
      affine_wavefront_t *const dwavefront = affine_wavefronts->dwavefronts[s];
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(mwavefront, k);
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(iwavefront, k);
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(dwavefront, k);
      fprintf(stream, " ");
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "SCORE   ");
  for (s = 0; s <= current_score; ++s)
  {
    fprintf(stream, "     %2d      ", s);
  }
  fprintf(stream, "\n");
  fprintf(stream, "\n");
}
void affine_wavefronts_print_wavefronts_pretty(
    FILE *const stream,
    affine_wavefronts_t *const affine_wavefronts,
    const int current_score)
{
  fprintf(stream, ">[SCORE=%3d]\n", current_score);
  int k, s;
  awf_offset_t offset;
  for (k = affine_wavefronts->max_k; k >= affine_wavefronts->min_k; k--)
  {
    fprintf(stream, "[k=%3d] ", k);
    awf_offset_t max_offset = affine_wavefronts_diagonal_length(affine_wavefronts, k);
    if (max_offset <= 0)
    {
      fprintf(stream, "\n");
      continue;
    }
    offset = 0;
    if (k > 0)
    {
      max_offset += k;
      for (; offset < k; ++offset)
        fprintf(stream, "   ");
    }
    int last_offset = -1;
    for (s = 0; s <= current_score; ++s)
    {
      affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[s];
      if (mwavefront == NULL)
        continue;
      if (k < mwavefront->lo || mwavefront->hi < k)
        continue;
      awf_offset_t offsets_sub = mwavefront->offsets[k];
      awf_offset_t offsets_sub_base = offsets_sub;
#ifdef AFFINE_WAVEFRONT_DEBUG
      offsets_sub_base = mwavefront->offsets_base[k];
#endif
      if (offsets_sub < 0)
        continue;
      if (last_offset == -1)
      {
        while (offset < offsets_sub_base - 1)
        {
          if (offset == max_offset)
            fprintf(stream, " | ");
          fprintf(stream, " * ");
          ++offset;
        }
      }
      while (offset < offsets_sub)
      {
        if (offset == max_offset)
          fprintf(stream, " | ");
        fprintf(stream, "%2d ", s);
        ++offset;
      }
      last_offset = offsets_sub;
    }
    for (; offset <= max_offset; ++offset)
    {
      if (offset == max_offset)
        fprintf(stream, " | ");
      fprintf(stream, " * ");
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}
void affine_wavefronts_debug_step(affine_wavefronts_t *const affine_wavefronts,
                                  const char *const pattern,
                                  const char *const text,
                                  const int score)
{
#ifdef AFFINE_WAVEFRONT_DEBUG
  affine_wavefronts_print_wavefronts_pretty(stderr, affine_wavefronts, score);
  affine_wavefronts_print_wavefronts(stderr, affine_wavefronts, score);
  affine_wavefronts_print_wavefront(stderr, affine_wavefronts, score);
#endif
}
