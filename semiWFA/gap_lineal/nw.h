#ifndef NW_H_
#define NW_H_
#include "utils/commons.h"
#include "edit/edit_table.h"
void nw_compute(
    edit_table_t* const edit_table,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    lineal_penalties_t* const penalties);
#endif 
