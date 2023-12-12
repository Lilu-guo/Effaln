#include "gap_lineal/nw.h"
void nw_traceback(
    edit_table_t* const edit_table,
    lineal_penalties_t* const penalties) {
  int** const dp = edit_table->columns;
  char* const operations = edit_table->edit_cigar.operations;
  int op_sentinel = edit_table->edit_cigar.end_offset-1;
  int h, v;
  h = edit_table->num_columns-1;
  v = edit_table->num_rows-1;
  while (h>0 && v>0) {
    if (dp[h][v] == dp[h][v-1]+penalties->deletion) {
      operations[op_sentinel--] = 'D';
      --v;
    } else if (dp[h][v] == dp[h-1][v]+penalties->insertion) {
      operations[op_sentinel--] = 'I';
      --h;
    } else {
      operations[op_sentinel--] =
          (dp[h][v] == dp[h-1][v-1]+penalties->mismatch) ? 'X' : 'M';
      --h;
      --v;
    }
  }
  while (h>0) {operations[op_sentinel--] = 'I'; --h;}
  while (v>0) {operations[op_sentinel--] = 'D'; --v;}
  edit_table->edit_cigar.begin_offset = op_sentinel+1;
}
void nw_compute(
    edit_table_t* const edit_table,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    lineal_penalties_t* const penalties) {
  int** dp = edit_table->columns;
  int h, v;
  dp[0][0] = 0;
  for (v=1;v<=pattern_length;++v) {
    dp[0][v] = dp[0][v-1] + penalties->deletion;
  }
  for (h=1;h<=text_length;++h) {
    dp[h][0] = dp[h-1][0] + penalties->insertion;
  }
  for (h=1;h<=text_length;++h) {
    for (v=1;v<=pattern_length;++v) {
      int min = dp[h-1][v-1] + ((pattern[v-1]==text[h-1]) ? 0 : penalties->mismatch); 
      min = MIN(min,dp[h-1][v]+penalties->insertion);       min = MIN(min,dp[h][v-1]+penalties->deletion);       dp[h][v] = min;
    }
  }
  nw_traceback(edit_table,penalties);
}
