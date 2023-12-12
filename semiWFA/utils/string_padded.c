#include "utils/string_padded.h"
#include "system/mm_allocator.h"
void strings_padded_add_padding(
    const char* const buffer, 
    const int buffer_length, 
    const int begin_padding_length, 
    const int end_padding_length, 
    const char padding_value, 
    char** const buffer_padded, 
    char** const buffer_padded_begin, 
    mm_allocator_t* const mm_allocator) {
  const int buffer_padded_length = begin_padding_length + buffer_length + end_padding_length; 
  *buffer_padded = mm_allocator_malloc(mm_allocator,buffer_padded_length);
  memset(*buffer_padded,padding_value,begin_padding_length); 
  *buffer_padded_begin = *buffer_padded + begin_padding_length; 
  memcpy(*buffer_padded_begin,buffer,buffer_length); 
  memset(*buffer_padded_begin+buffer_length,padding_value,end_padding_length); 
}
strings_padded_t* strings_padded_new(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int padding_length,
    mm_allocator_t* const mm_allocator) {
  strings_padded_t* const strings_padded =
      mm_allocator_alloc(mm_allocator,strings_padded_t);
  strings_padded->mm_allocator = mm_allocator;
  const int pattern_begin_padding_length = 0;
  const int pattern_end_padding_length = padding_length;
  const int text_begin_padding_length = 0;
  const int text_end_padding_length = padding_length;
  strings_padded_add_padding(
      pattern,pattern_length,
      pattern_begin_padding_length,pattern_end_padding_length,'X',
      &(strings_padded->pattern_padded_buffer),
      &(strings_padded->pattern_padded),mm_allocator);
  strings_padded_add_padding(
      text,text_length,
      text_begin_padding_length,text_end_padding_length,'Y',
      &(strings_padded->text_padded_buffer),
      &(strings_padded->text_padded),mm_allocator);
  return strings_padded;
}
strings_padded_t* strings_padded_new_rhomb(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int padding_length, 
    mm_allocator_t* const mm_allocator) {
  strings_padded_t* const strings_padded =
      mm_allocator_alloc(mm_allocator,strings_padded_t); 
  strings_padded->mm_allocator = mm_allocator;
  const int pattern_begin_padding_length = text_length + padding_length; 
  const int pattern_end_padding_length = pattern_length + text_length + padding_length; 
  const int text_begin_padding_length = padding_length; 
  const int text_end_padding_length = text_length + padding_length; 
  strings_padded_add_padding( 
      pattern,pattern_length, 
      pattern_begin_padding_length,pattern_end_padding_length,'X', 
      &(strings_padded->pattern_padded_buffer),
      &(strings_padded->pattern_padded),mm_allocator);
  strings_padded_add_padding( 
      text,text_length,
      text_begin_padding_length,text_end_padding_length,'Y',
      &(strings_padded->text_padded_buffer),
      &(strings_padded->text_padded),mm_allocator);
  return strings_padded; 
}
void strings_padded_delete(strings_padded_t* const strings_padded) {
  mm_allocator_free(strings_padded->mm_allocator,strings_padded->pattern_padded_buffer);
  mm_allocator_free(strings_padded->mm_allocator,strings_padded->text_padded_buffer);
  mm_allocator_free(strings_padded->mm_allocator,strings_padded);
}
