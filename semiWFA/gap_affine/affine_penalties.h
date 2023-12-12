#ifndef AFFINE_PENALTIES_H_
#define AFFINE_PENALTIES_H_
#include "../utils/commons.h"
typedef struct {
  int match;              
  int mismatch;             int gap_opening;          int gap_extension;      } affine_penalties_t;
#endif 
