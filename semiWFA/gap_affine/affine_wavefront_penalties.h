#ifndef AFFINE_WAVEFRONT_PENALTIES_H_
#define AFFINE_WAVEFRONT_PENALTIES_H_
#include "../utils/commons.h"
#include "affine_penalties.h"
typedef enum {
  wavefronts_penalties_match_zero,
  wavefronts_penalties_force_zero_match,
  wavefronts_penalties_shifted_penalties,
  wavefronts_penalties_odd_pair_penalties
} wavefronts_penalties_strategy;
typedef struct {
  affine_penalties_t base_penalties;                
  affine_penalties_t wavefront_penalties;             wavefronts_penalties_strategy penalties_strategy; } affine_wavefronts_penalties_t;
void affine_wavefronts_penalties_init(
    affine_wavefronts_penalties_t* const wavefronts_penalties,
    affine_penalties_t* const penalties,
    const wavefronts_penalties_strategy penalties_strategy);
void affine_penalties_mzero(
    affine_penalties_t* const base_penalties,
    affine_penalties_t* const shifted_penalties);
void affine_penalties_shift(
    affine_penalties_t* const base_penalties,
    affine_penalties_t* const shifted_penalties,
    const bool pair_odd_heuristic);
#endif 
