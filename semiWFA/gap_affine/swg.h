#ifndef SWG_H_
#define SWG_H_
#include "utils/commons.h"
#include "edit/edit_table.h"
#include "gap_affine/affine_table.h"
void swg_compute(
    affine_table_t* const affine_table,
    affine_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);
void swg_compute_banded(
    affine_table_t* const affine_table,
    affine_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int bandwidth);
#endif 
