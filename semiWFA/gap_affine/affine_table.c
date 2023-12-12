#include "edit/edit_table.h"
#include "gap_affine/affine_table.h"
void affine_table_allocate(
    affine_table_t *const table,
    const int pattern_length,
    const int text_length,
    mm_allocator_t *const mm_allocator)
{
  int h;
  table->num_rows = pattern_length + 1;
  const int num_columns = text_length + 1;
  table->num_columns = num_columns;
  table->columns = mm_allocator_malloc(mm_allocator, (text_length + 1) * sizeof(affine_cell_t *));
  for (h = 0; h < num_columns; ++h)
  {
    table->columns[h] = mm_allocator_calloc(mm_allocator, pattern_length + 1, affine_cell_t, false);
  }
  edit_cigar_allocate(&table->edit_cigar, pattern_length, text_length, mm_allocator);
}
void affine_table_free(
    affine_table_t *const table,
    mm_allocator_t *const mm_allocator)
{
  const int num_columns = table->num_columns;
  int h;
  for (h = 0; h < num_columns; ++h)
  {
    mm_allocator_free(mm_allocator, table->columns[h]);
  }
  mm_allocator_free(mm_allocator, table->columns);
  edit_cigar_free(&table->edit_cigar, mm_allocator);
}
#define AFFINE_TABLE_PRINT_VALUE(value) \
  if (value >= 0 && value < SCORE_MAX)  \
  {                                     \
    fprintf(stream, "%2d", value);      \
  }                                     \
  else                                  \
  {                                     \
    fprintf(stream, " *");              \
  }
void affine_table_print(
    FILE *const stream,
    const affine_table_t *const table, const char *const pattern,
    const char *const text)
{
  affine_cell_t **const dp = table->columns;
  int i, j;
  fprintf(stream, "     ");
  for (i = 0; i < table->num_columns - 1; ++i)
  {
    fprintf(stream, " %c ", text[i]);
  }
  fprintf(stream, "\n ");
  for (i = 0; i < table->num_columns; ++i)
  {
    fprintf(stream, " ");
    AFFINE_TABLE_PRINT_VALUE(dp[i][0].M);
  }
  fprintf(stream, "\n");
  for (i = 1; i < table->num_rows; ++i)
  {
    fprintf(stream, "%c", pattern[i - 1]);
    for (j = 0; j < table->num_columns; ++j)
    {
      fprintf(stream, " ");
      AFFINE_TABLE_PRINT_VALUE(dp[j][i].M);
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}
void affine_table_print_extended(
    FILE *const stream,
    const affine_table_t *const table,
    const char *const pattern,
    const char *const text)
{
  affine_cell_t **const dp = table->columns;
  int i, j;
  fprintf(stream, "         ");
  for (i = 0; i < table->num_columns - 1; ++i)
  {
    fprintf(stream, "     %c     ", text[i]);
  }
  fprintf(stream, "\n ");
  for (i = 0; i < table->num_columns; ++i)
  {
    fprintf(stream, " ");
    AFFINE_TABLE_PRINT_VALUE(dp[i][0].M);
    fprintf(stream, "{");
    AFFINE_TABLE_PRINT_VALUE(dp[i][0].I);
    fprintf(stream, ",");
    AFFINE_TABLE_PRINT_VALUE(dp[i][0].D);
    fprintf(stream, "} ");
  }
  fprintf(stream, "\n");
  for (i = 1; i < table->num_rows; ++i)
  {
    fprintf(stream, "%c", pattern[i - 1]);
    for (j = 0; j < table->num_columns; ++j)
    {
      fprintf(stream, " ");
      AFFINE_TABLE_PRINT_VALUE(dp[j][i].M);
      fprintf(stream, "{");
      AFFINE_TABLE_PRINT_VALUE(dp[j][i].I);
      fprintf(stream, ",");
      AFFINE_TABLE_PRINT_VALUE(dp[j][i].D);
      fprintf(stream, "} ");
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}
