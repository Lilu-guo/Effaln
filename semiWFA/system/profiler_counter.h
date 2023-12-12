#ifndef PROFILER_COUNTER_H_
#define PROFILER_COUNTER_H_
#include "../utils/commons.h"
typedef struct {
  uint64_t total;
  uint64_t samples;
  uint64_t min;
  uint64_t max;
  double m_oldM;
  double m_newM;
  double m_oldS;
  double m_newS;
} profiler_counter_t;
void counter_reset(
    profiler_counter_t* const counter);
void counter_add(
    profiler_counter_t* const counter,
    const uint64_t amount);
uint64_t counter_get_total(const profiler_counter_t* const counter);
uint64_t counter_get_num_samples(const profiler_counter_t* const counter);
uint64_t counter_get_min(const profiler_counter_t* const counter);
uint64_t counter_get_max(const profiler_counter_t* const counter);
double counter_get_mean(const profiler_counter_t* const counter);
double counter_get_variance(const profiler_counter_t* const counter);
double counter_get_stddev(const profiler_counter_t* const counter);
void counter_combine_sum(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src);
void counter_combine_max(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src);
void counter_combine_min(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src);
void counter_combine_mean(
    profiler_counter_t* const counter_dst,
    profiler_counter_t* const counter_src);
void counter_print(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const profiler_counter_t* const ref_counter,
    const char* const units,
    const bool full_report);
void sampler_print(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const profiler_counter_t* const ref_counter,
    const char* const units);
void percentage_print(
    FILE* const stream,
    const profiler_counter_t* const counter,
    const char* const units);
typedef struct {
  uint64_t begin_count;       
  profiler_counter_t counter;   uint64_t accumulated;       } profiler_rcounter_t;
void rcounter_start(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference);
void rcounter_stop(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference);
void rcounter_pause(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference);
void rcounter_continue(
    profiler_rcounter_t* const rcounter,
    const uint64_t reference);
void rcounter_reset(
    profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_total(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_num_samples(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_min(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_max(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_mean(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_variance(profiler_rcounter_t* const rcounter);
uint64_t rcounter_get_stddev(profiler_rcounter_t* const rcounter);
#define PRIcounter "lu(#%"PRIu64",m%"PRIu64",M%"PRIu64",{%.2f})"
#define PRIcounterVal(counter) \
  counter_get_total(counter), \
  counter_get_num_samples(counter), \
  counter_get_min(counter), \
  counter_get_max(counter), \
  counter_get_mean(counter)
#define PRIcounterX "lu(#%"PRIu64",m%"PRIu64",M%"PRIu64",{%.2f,%.2f,%.2f})"
#define PRIcounterXVal(counter) \
  counter_get_total(counter), \
  counter_get_num_samples(counter), \
  counter_get_min(counter), \
  counter_get_max(counter), \
  counter_get_mean(counter), \
  counter_get_variance(counter), \
  counter_get_stddev(counter)
#endif 
