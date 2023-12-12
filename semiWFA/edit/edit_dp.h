#ifndef EDIT_DP_H_
#define EDIT_DP_H_
#include "utils/commons.h"
#include "edit/edit_table.h"
void edit_dp_compute(
    edit_table_t *const edit_table,
    const char *const pattern,
    const int pattern_length,
    const char *const text,
    const int text_length);
void edit_dp_compute_banded(
    edit_table_t *const edit_table,
    const char *const pattern,
    const int pattern_length,
    const char *const text,
    const int text_length,
    const int bandwidth);
#endif
