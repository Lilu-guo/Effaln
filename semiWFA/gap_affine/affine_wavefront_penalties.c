#include "gap_affine/affine_wavefront_penalties.h"
void affine_wavefronts_penalties_init(
    affine_wavefronts_penalties_t* const wavefronts_penalties,
    affine_penalties_t* const penalties,
    const wavefronts_penalties_strategy penalties_strategy) {
  wavefronts_penalties->base_penalties = *penalties;
  wavefronts_penalties->penalties_strategy =
      (penalties->match==0) ? wavefronts_penalties_match_zero : penalties_strategy;
  switch (wavefronts_penalties->penalties_strategy) {
    case wavefronts_penalties_match_zero:
    case wavefronts_penalties_force_zero_match:
      affine_penalties_mzero(penalties,&(wavefronts_penalties->wavefront_penalties));
      break;
    case wavefronts_penalties_shifted_penalties:
      affine_penalties_shift(penalties,&(wavefronts_penalties->wavefront_penalties),false);
      break;
    case wavefronts_penalties_odd_pair_penalties:
      affine_penalties_shift(penalties,&(wavefronts_penalties->wavefront_penalties),true);
      break;
    default:
      break;
  }
}
void affine_penalties_mzero(
    affine_penalties_t* const base_penalties,
    affine_penalties_t* const shifted_penalties) {
  if (base_penalties->match > 0) {
    fprintf(stderr,"Match score must be negative or zero (M=%d)\n",base_penalties->match);
    exit(1);
  }
  if (base_penalties->mismatch <= 0 ||
      base_penalties->gap_opening <= 0 ||
      base_penalties->gap_extension <= 0) {
    fprintf(stderr,"Mismatch/Gap scores must be strictly positive (X=%d,O=%d,E=%d)\n",
        base_penalties->mismatch,base_penalties->gap_opening,base_penalties->gap_extension);
    exit(1);
  }
  *shifted_penalties = *base_penalties;
  shifted_penalties->match = 0;
}
void affine_penalties_shift(
    affine_penalties_t* const base_penalties,
    affine_penalties_t* const shifted_penalties,
    const bool pair_odd_heuristic) {
  if (base_penalties->match > 0) {
    fprintf(stderr,"Match score must be negative (M=%d)\n",base_penalties->match);
    exit(1);
  }
  if (base_penalties->mismatch <= 0 ||
      base_penalties->gap_opening <= 0 ||
      base_penalties->gap_extension <= 0) {
    fprintf(stderr,"Mismatch/Gap scores must be strictly positive (X=%d,O=%d,E=%d)\n",
        base_penalties->mismatch,base_penalties->gap_opening,base_penalties->gap_extension);
    exit(1);
  }
  *shifted_penalties = *base_penalties;
  shifted_penalties->match = 0;
  shifted_penalties->mismatch -= base_penalties->match;
  shifted_penalties->gap_opening -= base_penalties->match;
  shifted_penalties->gap_extension -= base_penalties->match;
  if (pair_odd_heuristic) {
    const bool is_mismatch_pair = ((shifted_penalties->mismatch%2)==0);
    const bool is_gap_opening_pair = ((shifted_penalties->gap_opening%2)==0);
    const bool is_gap_extension_pair = ((shifted_penalties->gap_extension%2)==0);
    const int total_odd = !is_mismatch_pair + !is_gap_opening_pair + !is_gap_extension_pair;
    const int total_pair = is_mismatch_pair + is_gap_opening_pair + is_gap_extension_pair;
    if (total_odd > total_pair) {
      if (is_mismatch_pair) ++(shifted_penalties->mismatch);
      if (is_gap_opening_pair) ++(shifted_penalties->gap_opening);
      if (is_gap_extension_pair) ++(shifted_penalties->gap_extension);
    } else {
      if (!is_mismatch_pair) ++(shifted_penalties->mismatch);
      if (!is_gap_opening_pair) ++(shifted_penalties->gap_opening);
      if (!is_gap_extension_pair) ++(shifted_penalties->gap_extension);
    }
  }
}
