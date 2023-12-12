#ifndef AFFINE_TABLE_H_
#define AFFINE_TABLE_H_
#include "../utils/commons.h"
#include "../edit/edit_cigar.h"
typedef struct {
  int M;   int I;   int D; } affine_cell_t;
typedef struct {
    affine_cell_t** columns;   int num_rows;   int num_columns;     edit_cigar_t edit_cigar; } affine_table_t;
void affine_table_allocate(
    affine_table_t* const table,
    const int pattern_length,
    const int text_length,
    mm_allocator_t* const mm_allocator);
void affine_table_free(
    affine_table_t* const table,
    mm_allocator_t* const mm_allocator);
void affine_table_print(
    FILE* const stream,
    const affine_table_t* const table,
    const char* const pattern,
    const char* const text);
void affine_table_print_extended(
    FILE* const stream,
    const affine_table_t* const table,
    const char* const pattern,
    const char* const text);
#endif 
