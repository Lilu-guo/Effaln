#ifndef WAVEFRONT_STATS_H_
#define WAVEFRONT_STATS_H_
#include "../utils/commons.h"
#include "../system/profiler_counter.h"
#include "../system/profiler_timer.h"
#ifdef AFFINE_WAVEFRONT_STATS
  #define WAVEFRONT_STATS_COUNTER_ADD(wavefronts,counter,amount) \
    if (wavefronts->wf_stats!=NULL) counter_add(&(wavefronts->wavefronts_stats->counter),(amount))
  #define WAVEFRONT_STATS_TIMER_START(wavefronts,timer) \
    if (wavefronts->wf_stats!=NULL) timer_start(&(wavefronts->wavefronts_stats->timer))
  #define WAVEFRONT_STATS_TIMER_STOP(wavefronts,timer) \
    if (wavefronts->wf_stats!=NULL) timer_stop(&(wavefronts->wavefronts_stats->timer))
#else
  #define WAVEFRONT_STATS_COUNTER_ADD(wf,counter,amount)
  #define WAVEFRONT_STATS_TIMER_START(wf,timer)
  #define WAVEFRONT_STATS_TIMER_STOP(wf,timer)
#endif
typedef struct {
  profiler_counter_t wf_score;              
  profiler_counter_t wf_steps;                profiler_counter_t wf_steps_null;           profiler_counter_t wf_steps_extra;          profiler_counter_t wf_operations;           profiler_counter_t wf_extensions;           profiler_counter_t wf_reduction;            profiler_counter_t wf_reduced_cells;        profiler_counter_t wf_null_used;            profiler_counter_t wf_extend_inner_loop;    profiler_counter_t wf_compute_kernel[4];    profiler_timer_t wf_time_backtrace;         profiler_counter_t wf_backtrace_paths;      profiler_counter_t wf_backtrace_alg;      } wavefronts_stats_t;
void wavefronts_stats_clear(wavefronts_stats_t* const wavefronts_stats);
void wavefronts_stats_print(
    FILE* const stream,
    wavefronts_stats_t* const wavefronts_stats);
#endif 
