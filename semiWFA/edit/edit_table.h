#ifndef EDIT_TABLE_H_
#define EDIT_TABLE_H_
#include "utils/commons.h"
#include "system/mm_allocator.h"
#include "edit/edit_cigar.h"
#define SCORE_MAX (10000000)
typedef struct
{
    int **columns;
    int num_rows;
    int num_columns;
    edit_cigar_t edit_cigar;
} edit_table_t;
void edit_table_allocate(
    edit_table_t *const edit_table,
    const int pattern_length,
    const int text_length,
    mm_allocator_t *const mm_allocator);
void edit_table_free(
    edit_table_t *const edit_table,
    mm_allocator_t *const mm_allocator);
void edit_table_print(
    FILE *const stream,
    const edit_table_t *const edit_table,
    const char *const pattern,
    const char *const text);
#endif
