#include "gap_affine/wavefront_stats.h"
void wavefronts_stats_clear(wavefronts_stats_t* const wavefronts_stats) {
  counter_reset(&(wavefronts_stats->wf_score));
  counter_reset(&(wavefronts_stats->wf_steps));
  counter_reset(&(wavefronts_stats->wf_steps_null));
  counter_reset(&(wavefronts_stats->wf_steps_extra));
  counter_reset(&(wavefronts_stats->wf_operations));
  counter_reset(&(wavefronts_stats->wf_extensions));
  counter_reset(&(wavefronts_stats->wf_reduction));
  counter_reset(&(wavefronts_stats->wf_reduced_cells));
  counter_reset(&(wavefronts_stats->wf_null_used));
  counter_reset(&(wavefronts_stats->wf_extend_inner_loop));
  int i;
  for (i=0;i<4;++i) {
    counter_reset(&(wavefronts_stats->wf_compute_kernel[i]));
  }
  timer_reset(&(wavefronts_stats->wf_time_backtrace));
  counter_reset(&(wavefronts_stats->wf_backtrace_paths));
  counter_reset(&(wavefronts_stats->wf_backtrace_alg));
}
void wavefronts_stats_print(
    FILE* const stream,
    wavefronts_stats_t* const wavefronts_stats) {
  fprintf(stream,"[Wavefront.Computations]\n");
  fprintf(stream," => Score                  ");
  counter_print(stream,&wavefronts_stats->wf_score,NULL,"score     ",true);
  fprintf(stream," => Steps                  ");
  counter_print(stream,&wavefronts_stats->wf_steps,NULL,"steps     ",true);
  fprintf(stream,"   => Steps.Null           ");
  counter_print(stream,&wavefronts_stats->wf_steps_null,
                       &wavefronts_stats->wf_steps,"steps     ",true);
  fprintf(stream,"   => WF.Null.Used         ");
  counter_print(stream,&wavefronts_stats->wf_null_used,NULL,"wavefront ",true);
  fprintf(stream,"   => Steps.Extra          ");
  counter_print(stream,&wavefronts_stats->wf_steps_extra,
                       &wavefronts_stats->wf_steps,"steps     ",true);
  fprintf(stream," => Wavefront.Operations   ");
  counter_print(stream,&wavefronts_stats->wf_operations,NULL,"operations",true);
  int i;
  for (i=0;i<4;++i) {
    fprintf(stream,"   => WF.Compute.Kernel.%02d  ",i);
    counter_print(stream,&(wavefronts_stats->wf_compute_kernel[i]),NULL,"operations",true);
  }
  fprintf(stream," => Extend.Operations      ");
  counter_print(stream,&wavefronts_stats->wf_extensions,NULL,"extensions",true);
  fprintf(stream,"   => Inner.Loop.Iters     ");
  counter_print(stream,&wavefronts_stats->wf_extend_inner_loop,
                       &wavefronts_stats->wf_extensions,"iterations",true);
  fprintf(stream," => Reduction.Calls        ");
  counter_print(stream,&wavefronts_stats->wf_reduction,NULL,"calls     ",true);
  fprintf(stream,"   => Reduction.Cells      ");
  counter_print(stream,&wavefronts_stats->wf_reduced_cells,NULL,"cells     ",true);
  fprintf(stream," => Time.Backtrace         ");
  timer_print(stream,&wavefronts_stats->wf_time_backtrace,NULL);
  fprintf(stream,"   => Backtrace.Paths      ");
  counter_print(stream,&wavefronts_stats->wf_backtrace_paths,NULL,"paths     ",true);
  fprintf(stream,"   => Backtrace.Alignments ");
  counter_print(stream,&wavefronts_stats->wf_backtrace_alg,NULL,"alignments",true);
}
