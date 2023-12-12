#ifndef STRING_PADDED_H
#define STRING_PADDED_H
#include "utils/commons.h"
#include "system/mm_allocator.h"
typedef struct {
  char* pattern_padded_buffer;
  char* pattern_padded;
  char* text_padded_buffer;
  char* text_padded;
  mm_allocator_t* mm_allocator;
} strings_padded_t;
strings_padded_t* strings_padded_new(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int padding_length,
    mm_allocator_t* const mm_allocator);
strings_padded_t* strings_padded_new_rhomb(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int padding_length,
    mm_allocator_t* const mm_allocator);
void strings_padded_delete(
    strings_padded_t* const strings_padded);
#endif 
