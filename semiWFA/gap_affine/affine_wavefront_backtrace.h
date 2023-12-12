#ifndef AFFINE_WAVEFRONT_BACKTRACE_H_
#define AFFINE_WAVEFRONT_BACKTRACE_H_
#include "gap_affine/affine_wavefront.h"
typedef struct {
  char* pattern;
  int pattern_length;
  char* text;
  int text_length;
} alignment_sequences_t;
typedef enum {
  backtrace_wavefront_M = 0,
  backtrace_wavefront_I = 1,
  backtrace_wavefront_D = 2
} backtrace_wavefront_type;
void affine_wavefronts_backtrace(
    affine_wavefronts_t* const affine_wavefronts,
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int alignment_score,
    int alignment_k);
#endif 
