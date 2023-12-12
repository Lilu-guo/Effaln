#ifndef AFFINE_WAVEFRONT_UTILS_H_
#define AFFINE_WAVEFRONT_UTILS_H_
#include "gap_affine/affine_wavefront.h"
affine_wavefront_t* affine_wavefronts_get_source_mwavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const int score);
affine_wavefront_t* affine_wavefronts_get_source_iwavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const int score);
affine_wavefront_t* affine_wavefronts_get_source_dwavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const int score);
int affine_wavefronts_diagonal_length(
    affine_wavefronts_t* const affine_wavefronts,
    const int k);
int affine_wavefronts_compute_distance(
    const int pattern_length,
    const int text_length,
    const awf_offset_t offset,
    const int k);
void affine_wavefront_initialize(
    affine_wavefronts_t* const affine_wavefronts);
bool affine_wavefront_end_reached(
    affine_wavefronts_t* const affine_wavefronts,
    const int pattern_length,
    const int text_length,
    const int score,
    int* alignment_k); 
#endif 
