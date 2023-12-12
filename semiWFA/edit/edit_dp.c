#include "edit_dp.h"
void edit_dp_traceback(
    edit_table_t *const edit_table)
{
  int **const dp = edit_table->columns;
  char *const operations = edit_table->edit_cigar.operations;
  int op_sentinel = edit_table->edit_cigar.end_offset - 1;
  int h, v;
  h = edit_table->num_columns - 1;
  v = edit_table->num_rows - 1;
  while (h > 0 && v > 0)
  {
    if (dp[h][v] == dp[h][v - 1] + 1)
    {
      operations[op_sentinel--] = 'D';
      --v;
    }
    else if (dp[h][v] == dp[h - 1][v] + 1)
    {
      operations[op_sentinel--] = 'I';
      --h;
    }
    else if (dp[h][v] == dp[h - 1][v - 1])
    {
      operations[op_sentinel--] = 'M';
      --h;
      --v;
    }
    else if (dp[h][v] == dp[h - 1][v - 1] + 1)
    {
      operations[op_sentinel--] = 'X';
      --h;
      --v;
    }
    else
    {
      fprintf(stderr, "Edit backtrace. No backtrace operation found");
      exit(1);
    }
  }
  while (h > 0)
  {
    operations[op_sentinel--] = 'I';
    --h;
  }
  while (v > 0)
  {
    operations[op_sentinel--] = 'D';
    --v;
  }
  edit_table->edit_cigar.begin_offset = op_sentinel + 1;
}
void edit_dp_compute(
    edit_table_t *const edit_table,
    const char *const pattern,
    const int pattern_length,
    const char *const text,
    const int text_length)
{
  int **dp = edit_table->columns;
  int h, v;
  for (v = 0; v <= pattern_length; ++v)
    dp[0][v] = v;
  for (h = 0; h <= text_length; ++h)
    dp[h][0] = h;
  for (h = 1; h <= text_length; ++h)
  {
    for (v = 1; v <= pattern_length; ++v)
    {
      int min = dp[h - 1][v - 1] + (text[h - 1] != pattern[v - 1]);
      min = MIN(min, dp[h - 1][v] + 1);
      min = MIN(min, dp[h][v - 1] + 1);
      dp[h][v] = min;
    }
  }
  edit_dp_traceback(edit_table);
  edit_table_print(stderr, edit_table, pattern, text);
}
void edit_dp_compute_banded(
    edit_table_t *const edit_table,
    const char *const pattern,
    const int pattern_length,
    const char *const text,
    const int text_length,
    const int bandwidth)
{
  const int k_end = ABS(text_length - pattern_length) + 1;
  const int effective_bandwidth = MAX(k_end, bandwidth);
  int **dp = edit_table->columns;
  int h, v;
  dp[0][0] = 0;
  for (v = 1; v <= effective_bandwidth; ++v)
    dp[0][v] = v;
  for (h = 1; h <= text_length; ++h)
  {
    const bool lo_band = (h <= effective_bandwidth);
    const int lo = (lo_band) ? 1 : h - effective_bandwidth;
    dp[h][lo - 1] = (lo_band) ? h : INT16_MAX;
    const int hi = MIN(pattern_length, effective_bandwidth + h - 1);
    if (h > 1)
      dp[h - 1][hi] = INT16_MAX;
    for (v = lo; v <= hi; ++v)
    {
      const int sub = dp[h - 1][v - 1] + (text[h - 1] != pattern[v - 1]);
      const int ins = dp[h - 1][v];
      const int del = dp[h][v - 1];
      dp[h][v] = MIN(MIN(ins, del) + 1, sub);
    }
  }
  edit_dp_traceback(edit_table);
  edit_table_print(stderr, edit_table, pattern, text);
}
