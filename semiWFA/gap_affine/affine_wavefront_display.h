#ifndef AFFINE_WAVEFRONT_DISPLAY_H_
#define AFFINE_WAVEFRONT_DISPLAY_H_
#include "gap_affine/affine_table.h"
#include "gap_affine/affine_wavefront.h"
#include "utils/commons.h"
void affine_wavefronts_set_edit_table(
    affine_wavefronts_t* const affine_wavefronts,
    const int pattern_length,
    const int text_length,
    const int k,
    const awf_offset_t offset,
    const int score);
void affine_wavefronts_print_wavefront(
    FILE* const stream,
    affine_wavefronts_t* const affine_wavefronts,
    const int current_score);
void affine_wavefronts_print_wavefronts(
    FILE* const stream,
    affine_wavefronts_t* const affine_wavefronts,
    const int current_score);
void affine_wavefronts_print_wavefronts_pretty(
    FILE* const stream,
    affine_wavefronts_t* const affine_wavefronts,
    const int current_score);
void affine_wavefronts_debug_step(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const char* const text,
    const int score);
#endif 
