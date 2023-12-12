#include "gap_affine/affine_wavefront_utils.h"
#include <stdio.h>
affine_wavefront_t *affine_wavefronts_get_source_mwavefront(
    affine_wavefronts_t *const affine_wavefronts,
    const int score)
{
  return (score < 0 || affine_wavefronts->mwavefronts[score] == NULL) ? &affine_wavefronts->wavefront_null : affine_wavefronts->mwavefronts[score];
}
affine_wavefront_t *affine_wavefronts_get_source_iwavefront(
    affine_wavefronts_t *const affine_wavefronts,
    const int score)
{
  return (score < 0 || affine_wavefronts->iwavefronts[score] == NULL) ? &affine_wavefronts->wavefront_null : affine_wavefronts->iwavefronts[score];
}
affine_wavefront_t *affine_wavefronts_get_source_dwavefront(
    affine_wavefronts_t *const affine_wavefronts,
    const int score)
{
  return (score < 0 || affine_wavefronts->dwavefronts[score] == NULL) ? &affine_wavefronts->wavefront_null : affine_wavefronts->dwavefronts[score];
}
int affine_wavefronts_diagonal_length(
    affine_wavefronts_t *const affine_wavefronts,
    const int k)
{
  if (k >= 0)
  {
    return MIN(affine_wavefronts->text_length - k, affine_wavefronts->pattern_length);
  }
  else
  {
    return MIN(affine_wavefronts->pattern_length + k, affine_wavefronts->text_length);
  }
}
int affine_wavefronts_compute_distance(
    const int pattern_length,
    const int text_length,
    const awf_offset_t offset,
    const int k)
{
  const int v = AFFINE_WAVEFRONT_V(k, offset);
  const int h = AFFINE_WAVEFRONT_H(k, offset);
  const int left_v = pattern_length - v;
  const int left_h = text_length - h;
  return MAX(left_v, left_h);
}
void affine_wavefront_initialize(
    affine_wavefronts_t *const affine_wavefronts)
{
  int pattern_length = affine_wavefronts->pattern_length;
  affine_wavefronts->mwavefronts[0] = affine_wavefronts_allocate_wavefront(affine_wavefronts, -pattern_length, 0);
  affine_wavefronts->mwavefronts[0]->offsets[0] = 0;
  for (int k = 1; k <= pattern_length; k++)
  {
    affine_wavefronts->mwavefronts[0]->offsets[-k] = 0;
  }
}
bool affine_wavefront_end_reached(
    affine_wavefronts_t *const affine_wavefronts,
    const int pattern_length,
    const int text_length,
    const int score,
    int *alignment_k)
{
  int k = 0;
  const int alignment_offset = AFFINE_WAVEFRONT_OFFSET(text_length, pattern_length);
  affine_wavefront_t *const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront != NULL)
  {
    awf_offset_t *const offsets = mwavefront->offsets;
    for (k = mwavefront->lo; k <= mwavefront->hi; k++)
    {
      if (offsets[k] >= alignment_offset)
      {
        *alignment_k = k;
        return true;
      }
    }
  }
  return false;
}
