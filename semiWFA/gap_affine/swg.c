#include "gap_affine/swg.h"
typedef enum {
  swg_M_layer,
  swg_I_layer,
  swg_D_layer
} swg_layer_type;
void swg_traceback(
    affine_table_t* const affine_table,
    affine_penalties_t* const penalties) {
  affine_cell_t** const dp = affine_table->columns;
  char* const operations = affine_table->edit_cigar.operations;
  int op_sentinel = affine_table->edit_cigar.end_offset-1;
  int h, v;
  h = affine_table->num_columns-1;
  v = affine_table->num_rows-1;
  swg_layer_type swg_layer = swg_M_layer;
  while (h>0 && v>0) {
    switch (swg_layer) {
      case swg_D_layer:
        operations[op_sentinel--] = 'D';
        if (dp[h][v].D == dp[h][v-1].M + penalties->gap_opening + penalties->gap_extension) {
          swg_layer = swg_M_layer;
        }
        --v;
        break;
      case swg_I_layer:
        operations[op_sentinel--] = 'I';
        if (dp[h][v].I == dp[h-1][v].M + penalties->gap_opening + penalties->gap_extension) {
          swg_layer = swg_M_layer;
        }
        --h;
        break;
      case swg_M_layer:
        if (dp[h][v].M == dp[h][v].D) {
          swg_layer = swg_D_layer;
        } else if (dp[h][v].M == dp[h][v].I) {
          swg_layer = swg_I_layer;
        } else if (dp[h][v].M == dp[h-1][v-1].M + penalties->match) {
          operations[op_sentinel--] = 'M';
          --h; --v;
        } else if (dp[h][v].M == dp[h-1][v-1].M + penalties->mismatch) {
          operations[op_sentinel--] = 'X';
          --h; --v;
        } else {
          fprintf(stderr,"SWG backtrace. No backtrace operation found");
          exit(1);
        }
        break;
    }
  }
  while (h>0) {operations[op_sentinel--] = 'I'; --h;}
  while (v>0) {operations[op_sentinel--] = 'D'; --v;}
  affine_table->edit_cigar.begin_offset = op_sentinel+1;
}
void swg_compute(
    affine_table_t* const affine_table,
    affine_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  affine_cell_t** const dp = affine_table->columns;
  int h, v;
  dp[0][0].D = SCORE_MAX;
  dp[0][0].I = SCORE_MAX;
  dp[0][0].M = 0;
  for (v=1;v<=pattern_length;++v) { 
    dp[0][v].D = penalties->gap_opening + v*penalties->gap_extension;
    dp[0][v].I = SCORE_MAX;
    dp[0][v].M = dp[0][v].D;
  }
  for (h=1;h<=text_length;++h) { 
    dp[h][0].D = SCORE_MAX;
    dp[h][0].I = penalties->gap_opening + h*penalties->gap_extension;
    dp[h][0].M = dp[h][0].I;
  }
  for (h=1;h<=text_length;++h) {
    for (v=1;v<=pattern_length;++v) {
      const int del_new = dp[h][v-1].M + penalties->gap_opening + penalties->gap_extension;
      const int del_ext = dp[h][v-1].D + penalties->gap_extension;
      const int del = MIN(del_new,del_ext);
      dp[h][v].D = del;
      const int ins_new = dp[h-1][v].M + penalties->gap_opening + penalties->gap_extension;
      const int ins_ext = dp[h-1][v].I + penalties->gap_extension;
      const int ins = MIN(ins_new,ins_ext);
      dp[h][v].I = ins;
      const int m_match = dp[h-1][v-1].M + ((pattern[v-1]==text[h-1]) ? penalties->match : penalties->mismatch);
      dp[h][v].M = MIN(m_match,MIN(ins,del));
    }
  }
  swg_traceback(affine_table,penalties);
  }
void swg_compute_banded(
    affine_table_t* const affine_table,
    affine_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int bandwidth) {
  const int k_end = ABS(text_length-pattern_length)+1;
  int effective_bandwidth = bandwidth;
  if (effective_bandwidth < k_end) effective_bandwidth = k_end;
  if (effective_bandwidth > pattern_length) effective_bandwidth = pattern_length;
  if (effective_bandwidth > text_length) effective_bandwidth = text_length;
  affine_cell_t** const dp = affine_table->columns;
  int h, v;
  dp[0][0].D = SCORE_MAX;
  dp[0][0].I = SCORE_MAX;
  dp[0][0].M = 0;
  for (v=1;v<=effective_bandwidth;++v) {
    dp[0][v].D = penalties->gap_opening + v*penalties->gap_extension;
    dp[0][v].I = SCORE_MAX;
    dp[0][v].M = dp[0][v].D;
  }
  for (h=1;h<=text_length;++h) {
    int lo;
    if (h <= effective_bandwidth) {
      lo = 1;
      dp[h][lo-1].D = SCORE_MAX;
      dp[h][lo-1].I = penalties->gap_opening + h*penalties->gap_extension;
      dp[h][lo-1].M = dp[h][lo-1].I;
    } else {
      lo = h - effective_bandwidth;
      dp[h][lo-1].D = SCORE_MAX;
      dp[h][lo-1].I = SCORE_MAX;
      dp[h][lo-1].M = SCORE_MAX;
    }
    int hi = h + effective_bandwidth - 1;
    if (hi > pattern_length) {
      hi = pattern_length;
    } else if (h > 1) {
      dp[h-1][hi].D = SCORE_MAX;
      dp[h-1][hi].I = SCORE_MAX;
      dp[h-1][hi].M = SCORE_MAX;
    }
    for (v=lo;v<=hi;++v) {
      const int del_new = dp[h][v-1].M + penalties->gap_opening + penalties->gap_extension;
      const int del_ext = dp[h][v-1].D + penalties->gap_extension;
      const int del = MIN(del_new,del_ext);
      dp[h][v].D = del;
      const int ins_new = dp[h-1][v].M + penalties->gap_opening + penalties->gap_extension;
      const int ins_ext = dp[h-1][v].I + penalties->gap_extension;
      const int ins = MIN(ins_new,ins_ext);
      dp[h][v].I = ins;
      const int m_match = dp[h-1][v-1].M + ((pattern[v-1]==text[h-1]) ? penalties->match : penalties->mismatch);
      dp[h][v].M = MIN(m_match,MIN(ins,del));
    }
  }
  swg_traceback(affine_table,penalties);
}
