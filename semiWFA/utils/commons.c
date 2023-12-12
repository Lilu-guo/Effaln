#include "commons.h"
uint64_t rand_iid(const uint64_t min,const uint64_t max) {
  int n_rand = rand(); 
  const uint64_t range = max - min;
  const uint64_t rem = RAND_MAX % range;
  const uint64_t sample = RAND_MAX / range;
  if (n_rand < RAND_MAX - rem) {
    return min + n_rand/sample;
  } else {
    return rand_iid(min,max);
  }
}
