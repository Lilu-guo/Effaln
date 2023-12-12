#include "edit/edit_table.h"
#include "edit/edit_cigar.h"
void edit_table_allocate(
    edit_table_t *const edit_table,
    const int pattern_length,
    const int text_length,
    mm_allocator_t *const mm_allocator)
{
  int h;
  edit_table->num_rows = pattern_length + 1;
  const int num_columns = text_length + 1;
  edit_table->num_columns = num_columns;
  edit_table->columns = mm_allocator_malloc(mm_allocator, (text_length + 1) * sizeof(int *));
  for (h = 0; h < num_columns; ++h)
  {
    edit_table->columns[h] = mm_allocator_calloc(mm_allocator, pattern_length + 1, int, false);
  }
  edit_cigar_allocate(&edit_table->edit_cigar, pattern_length, text_length, mm_allocator);
}
void edit_table_free(
    edit_table_t *const edit_table,
    mm_allocator_t *const mm_allocator)
{
  const int num_columns = edit_table->num_columns;
  int h;
  for (h = 0; h < num_columns; ++h)
  {
    mm_allocator_free(mm_allocator, edit_table->columns[h]);
  }
  mm_allocator_free(mm_allocator, edit_table->columns);
  edit_cigar_free(&edit_table->edit_cigar, mm_allocator);
}
void edit_table_print(
    FILE *const stream,
    const edit_table_t *const edit_table,
    const char *const pattern,
    const char *const text)
{
  printf(">>>>>>");
  int **const dp = edit_table->columns;
  int i, j;
  fprintf(stream, "       ");
  for (i = 0; i < edit_table->num_columns - 1; ++i)
  {
    fprintf(stream, "  %c  ", text[i]);
  }
  fprintf(stream, "\n ");
  for (i = 0; i < edit_table->num_columns; ++i)
  {
    fprintf(stream, " %3d ", i);
  }
  fprintf(stream, "\n ");
  for (i = 0; i < edit_table->num_columns; ++i)
  {
    if (-1 < dp[i][0] && dp[i][0] < 10000)
    {
      fprintf(stream, " %3d ", dp[i][0]);
    }
    else
    {
      fprintf(stream, "   * ");
    }
  }
  fprintf(stream, "\n");
  for (i = 1; i < edit_table->num_rows; ++i)
  {
    fprintf(stream, "%c", pattern[i - 1]);
    for (j = 0; j < edit_table->num_columns; ++j)
    {
      if (-1 < dp[j][i] && dp[j][i] < 10000)
      {
        fprintf(stream, " %3d ", dp[j][i]);
      }
      else
      {
        fprintf(stream, "   * ");
      }
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}
