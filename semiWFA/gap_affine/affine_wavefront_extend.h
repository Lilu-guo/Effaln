#ifndef AFFINE_WAVEFRONT_EXTEND_H_
#define AFFINE_WAVEFRONT_EXTEND_H_
#include "gap_affine/affine_wavefront.h"
#define AFFINE_WAVEFRONT_PADDING  10 
bool affine_wavefronts_extend_wavefront_packed(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score,
    int* alignment_k);
#endif 
