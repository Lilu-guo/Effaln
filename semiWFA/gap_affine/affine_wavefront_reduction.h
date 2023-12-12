#ifndef AFFINE_WAVEFRONT_REDUCTION_H_
#define AFFINE_WAVEFRONT_REDUCTION_H_
#include "../utils/commons.h"
typedef enum {
  wavefronts_reduction_none,
  wavefronts_reduction_dynamic,
} wavefront_reduction_type;
typedef struct {
  wavefront_reduction_type reduction_strategy;     
  int min_wavefront_length;                          int max_distance_threshold;                      } affine_wavefronts_reduction_t;
void affine_wavefronts_reduction_set_none(
    affine_wavefronts_reduction_t* const wavefronts_reduction);
void affine_wavefronts_reduction_set_dynamic(
    affine_wavefronts_reduction_t* const wavefronts_reduction,
    const int min_wavefront_length,
    const int max_distance_threshold);
#endif 
