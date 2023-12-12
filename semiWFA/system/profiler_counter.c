#include "profiler_counter.h"
void counter_reset(
    profiler_counter_t* const counter) {
  memset(counter,0,sizeof(profiler_counter_t));
}
void counter_add(
    profiler_counter_t* const counter,
    const uint64_t amount) {
  counter->total += amount;
  ++(counter->samples);
  if (counter->samples == 1) {
    counter->min = amount;
    counter->max = amount;
    counter->m_oldM = amount;
    counter->m_newM = amount;
    counter->m_oldS = 0.0;
  } else {
    counter->min = MIN(counter->min,amount);
    counter->max = MAX(counter->max,amount);
    counter->m_newM = counter->m_oldM + ((double)amount-counter->m_oldM)/(double)counter->samples;
    counter->m_newS = counter->m_oldS + ((double)amount-counter->m_oldM)*((double)amount-counter->m_newM);
    counter->m_oldM = counter->m_newM;
    counter->m_oldS = counter->m_newS;
  }
}
uint64_t counter_get_total(const profiler_counter_t* const counter) {
  return counter->total;
}
uint64_t counter_get_num_samples(const profiler_counter_t* const counter) {
  return counter->samples;
}
uint64_t counter_get_min(const profiler_counter_t* const counter) {
  return counter->min;
}
uint64_t counter_get_max(const profiler_counter_t* const counter) {
  return counter->max;
}
double counter_get_mean(const profiler_counter_t* const counter) {
  return (double)counter->total/(double)counter->samples;
}
double counter_get_variance(const profiler_counter_t* const counter) {
  return ((counter->samples > 1) ? counter->m_newS/(double)(counter->samples - 1) : 0.0);
}
double counter_get_stddev(const profiler_counter_t* const counter) {
  return sqrt(counter_get_variance(counter));
}
void counter_combine_sum(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src) {
  counter_dst->total += counter_src->total;
  counter_dst->samples += counter_src->samples;
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
}
void counter_combine_max(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src) {
  if (counter_dst->total < counter_src->total) {
    counter_dst->total = counter_src->total;
    counter_dst->samples = counter_src->samples;
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
  }
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
}
void counter_combine_min(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src) {
  if (counter_dst->samples==0 || counter_dst->total==0 || counter_dst->total > counter_src->total) {
    counter_dst->total = counter_src->total;
    counter_dst->samples = counter_src->samples;
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
  }
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
}
void counter_combine_mean(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src) {
  counter_dst->total = (counter_dst->total+counter_src->total)/2;
  counter_dst->samples = (counter_dst->samples+counter_src->samples)/2;
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
}
void counter_print_stats(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const profiler_counter_t* const ref_counter,
    const char* const units) {
  const uint64_t num_samples = counter_get_num_samples(counter);
  if (num_samples >= METRIC_FACTOR_1G) {
    fprintf(stream," (samples=%"PRIu64"G",num_samples/METRIC_FACTOR_1G);
  } else if (num_samples >= METRIC_FACTOR_1M) {
    fprintf(stream," (samples=%"PRIu64"M",num_samples/METRIC_FACTOR_1M);
  } else if (num_samples >= METRIC_FACTOR_1K) {
    fprintf(stream," (samples=%"PRIu64"K",num_samples/METRIC_FACTOR_1K);
  } else {
    fprintf(stream," (samples=%"PRIu64"",num_samples);
    if (num_samples==0) {
      fprintf(stream,",--n/a--)}\n");
      return;
    }
  }
  const double mean = counter_get_mean(counter);
  if (mean >= METRIC_FACTOR_1G) {
    fprintf(stream,"{mean%.2fG",mean/METRIC_FACTOR_1G);
  } else if (mean >= METRIC_FACTOR_1M) {
    fprintf(stream,"{mean%.2fM",mean/METRIC_FACTOR_1M);
  } else if (mean >= METRIC_FACTOR_1K) {
    fprintf(stream,"{mean%.2fK",mean/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,"{mean%.2f",mean);
  }
  const uint64_t min = counter_get_min(counter);
  if (min >= METRIC_FACTOR_1G) {
    fprintf(stream,",min%.2fG",(double)min/METRIC_FACTOR_1G);
  } else if (min >= METRIC_FACTOR_1M) {
    fprintf(stream,",min%.2fM",(double)min/METRIC_FACTOR_1M);
  } else if (min >= METRIC_FACTOR_1K) {
    fprintf(stream,",min%.2fK",(double)min/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,",min%.2f",(double)min);
  }
  const uint64_t max = counter_get_max(counter);
  if (max >= METRIC_FACTOR_1G) {
    fprintf(stream,",Max%.2fG",(double)max/METRIC_FACTOR_1G);
  } else if (max >= METRIC_FACTOR_1M) {
    fprintf(stream,",Max%.2fM",(double)max/METRIC_FACTOR_1M);
  } else if (max >= METRIC_FACTOR_1K) {
    fprintf(stream,",Max%.2fK",(double)max/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,",Max%.2f",(double)max);
  }
  const uint64_t var = counter_get_variance(counter);
  if (var >= METRIC_FACTOR_1G) {
    fprintf(stream,",Var%.2fG",(double)var/METRIC_FACTOR_1G);
  } else if (var >= METRIC_FACTOR_1M) {
    fprintf(stream,",Var%.2fM",(double)var/METRIC_FACTOR_1M);
  } else if (var >= METRIC_FACTOR_1K) {
    fprintf(stream,",Var%.2fK",(double)var/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,",Var%.2f",(double)var);
  }
  const uint64_t stdDev = counter_get_stddev(counter);
  if (stdDev >= METRIC_FACTOR_1G) {
    fprintf(stream,",StdDev%.2fG)}\n",(double)stdDev/METRIC_FACTOR_1G);
  } else if (stdDev >= METRIC_FACTOR_1M) {
    fprintf(stream,",StdDev%.2fM)}\n",(double)stdDev/METRIC_FACTOR_1M);
  } else if (stdDev >= METRIC_FACTOR_1K) {
    fprintf(stream,",StdDev%.2fK)}\n",(double)stdDev/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,",StdDev%.2f)}\n",(double)stdDev);
  }
}
void counter_print(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const profiler_counter_t* const ref_counter,
    const char* const units,
    const bool full_report) {
  const uint64_t total = counter_get_total(counter);
  if (total >= METRIC_FACTOR_1G) {
    fprintf(stream,"%7.2f G%s",(double)total/METRIC_FACTOR_1G,units);
  } else if (total >= METRIC_FACTOR_1M) {
    fprintf(stream,"%7.2f M%s",(double)total/METRIC_FACTOR_1M,units);
  } else if (total >= METRIC_FACTOR_1K) {
    fprintf(stream,"%7.2f K%s",(double)total/METRIC_FACTOR_1K,units);
  } else {
    fprintf(stream,"%7.2f %s ",(double)total,units);
  }
  if (ref_counter!=NULL) {
    if (total==0) {
        fprintf(stream," (  0.00 %%)");
    } else {
      const uint64_t total_ref = counter_get_total(ref_counter);
      if (total_ref==0) {
        fprintf(stream," (  n/a  %%)");
      } else {
        const double percentage = (double)(total*100)/(double)total_ref;
        fprintf(stream," (%6.02f %%)",percentage);
      }
    }
  } else {
    fprintf(stream,"           ");
  }
  if (!full_report) {
    fprintf(stream,"\n");
    return;
  } else {
    counter_print_stats(stream,counter,ref_counter,units);
  }
}
void sampler_print(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const profiler_counter_t* const ref_counter,
    const char* const units) {
  fprintf(stream,"\t\t\t\t");
  counter_print_stats(stream,counter,ref_counter,units);
}
void percentage_print(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const char* const units) {
  const double mean = counter_get_mean(counter);
  fprintf(stream,"%7.2f %%%s\t\t",mean,units);
  const uint64_t num_samples = counter_get_num_samples(counter);
  if (num_samples >= METRIC_FACTOR_1G) {
    fprintf(stream," (samples=%"PRIu64"G",num_samples/METRIC_FACTOR_1G);
  } else if (num_samples >= METRIC_FACTOR_1M) {
    fprintf(stream," (samples=%"PRIu64"M",num_samples/METRIC_FACTOR_1M);
  } else if (num_samples >= METRIC_FACTOR_1K) {
    fprintf(stream," (samples=%"PRIu64"K",num_samples/METRIC_FACTOR_1K);
  } else {
    fprintf(stream," (samples=%"PRIu64"",num_samples);
  }
  if (num_samples == 0) {
    fprintf(stream,")\n");
    return;
  }
  fprintf(stream,",min%.2f%%,Max%.2f%%",
      (double)counter_get_min(counter),(double)counter_get_max(counter));
  fprintf(stream,",Var%.2f,StdDev%.2f)\n",
      counter_get_variance(counter),counter_get_stddev(counter));
}
void rcounter_start(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference) {
  rcounter->accumulated = 0;
  rcounter_continue(rcounter,reference);
}
void rcounter_stop(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference) {
  rcounter_pause(rcounter,reference);
  counter_add(&rcounter->counter,rcounter->accumulated);
}
void rcounter_pause(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference) {
  rcounter->accumulated += reference-rcounter->begin_count;
}
void rcounter_continue(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference) {
  rcounter->begin_count = reference;
}
void rcounter_reset(
    profiler_rcounter_t* const rcounter) {
  counter_reset(&rcounter->counter);
}
uint64_t rcounter_get_total(profiler_rcounter_t* const rcounter) {
  return counter_get_total(&rcounter->counter);
}
uint64_t rcounter_get_num_samples(profiler_rcounter_t* const rcounter) {
  return counter_get_num_samples(&rcounter->counter);
}
uint64_t rcounter_get_min(profiler_rcounter_t* const rcounter) {
  return counter_get_min(&rcounter->counter);
}
uint64_t rcounter_get_max(profiler_rcounter_t* const rcounter) {
  return counter_get_max(&rcounter->counter);
}
uint64_t rcounter_get_mean(profiler_rcounter_t* const rcounter) {
  return counter_get_mean(&rcounter->counter);
}
uint64_t rcounter_get_variance(profiler_rcounter_t* const rcounter) {
  return counter_get_variance(&rcounter->counter);
}
uint64_t rcounter_get_stddev(profiler_rcounter_t* const rcounter) {
  return counter_get_stddev(&rcounter->counter);
}
