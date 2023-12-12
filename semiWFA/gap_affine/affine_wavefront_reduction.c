#include "gap_affine/affine_wavefront_reduction.h"
void affine_wavefronts_reduction_set_none(
    affine_wavefronts_reduction_t* const wavefronts_reduction) {
  wavefronts_reduction->reduction_strategy = wavefronts_reduction_none;
}
void affine_wavefronts_reduction_set_dynamic(
    affine_wavefronts_reduction_t* const wavefronts_reduction,
    const int min_wavefront_length,
    const int max_distance_threshold) {
  wavefronts_reduction->reduction_strategy = wavefronts_reduction_dynamic;
  wavefronts_reduction->min_wavefront_length = min_wavefront_length;
  wavefronts_reduction->max_distance_threshold = max_distance_threshold;
}
