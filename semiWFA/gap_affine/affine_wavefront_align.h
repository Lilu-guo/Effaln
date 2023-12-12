#ifndef AFFINE_WAVEFRONT_ALIGN_H_
#define AFFINE_WAVEFRONT_ALIGN_H_
#include "affine_wavefront.h"
#include "../utils/commons.h"
void affine_wavefronts_align(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);
#endif 
