#include "vector.h"
#define VECTOR_EXPAND_FACTOR (3.0/2.0)
vector_t* vector_new_(const uint64_t num_initial_elements,const uint64_t element_size) {
  vector_t* const vector_buffer = malloc(sizeof(vector_t));
  vector_buffer->element_size = element_size;
  vector_buffer->elements_allocated = num_initial_elements;
  vector_buffer->memory = malloc(num_initial_elements*element_size);
  if (!vector_buffer->memory) {
    fprintf(stderr,"Could not create new vector (%"PRIu64" bytes requested)",
        num_initial_elements*element_size);
    exit(1);
  }
  vector_buffer->used = 0;
  return vector_buffer;
}
void vector_reserve(vector_t* const vector,const uint64_t num_elements,const bool zero_mem) {
  if (vector->elements_allocated < num_elements) {
    const uint64_t proposed=(float)vector->elements_allocated*VECTOR_EXPAND_FACTOR;
    vector->elements_allocated = num_elements>proposed?num_elements:proposed;
    vector->memory = realloc(vector->memory,vector->elements_allocated*vector->element_size);
    if (!vector->memory) {
      fprintf(stderr,"Could not reserve vector (%"PRIu64" bytes requested)",
          vector->elements_allocated*vector->element_size);
      exit(1);
    }
  }
  if (zero_mem) {
    memset(vector->memory+vector->used*vector->element_size,0,
        (vector->elements_allocated-vector->used)*vector->element_size);
  }
}
void vector_resize__clear(vector_t* const vector,const uint64_t num_elements) {
  if (vector->elements_allocated < num_elements) {
    const uint64_t proposed = (float)vector->elements_allocated*VECTOR_EXPAND_FACTOR;
    vector->elements_allocated = (num_elements>proposed)?num_elements:proposed;
    free(vector->memory);
    vector->memory = malloc(vector->elements_allocated*vector->element_size);
    if (!vector->memory) {
      fprintf(stderr,"Could not reserve vector (%"PRIu64" bytes requested)",
          vector->elements_allocated*vector->element_size);
      exit(1);
    }
  }
  vector->used=0;
}
void vector_cast__clear_(vector_t* const vector,const uint64_t element_size) {
  vector->elements_allocated = (vector->elements_allocated*vector->element_size)/element_size;
  vector->element_size = element_size;
  vector->used = 0;
}
void vector_delete(vector_t* const vector) {
  free(vector->memory);
  free(vector);
}
#ifdef VECTOR_DEBUG
void* vector_get_mem_element(vector_t* const vector,const uint64_t position,const uint64_t element_size) {
  if (position >= (vector)->used) {
    fprintf(stderr,"Vector position out-of-range [0,%"PRIu64")",(vector)->used);
    exit(1);
  }
  return vector->memory + (position*element_size);
}
#endif
void vector_copy(vector_t* const vector_to,vector_t* const vector_from) {
  vector_cast__clear_(vector_to,vector_from->element_size);
  vector_reserve(vector_to,vector_from->used,false);
  vector_set_used(vector_to,vector_from->used);
  memcpy(vector_to->memory,vector_from->memory,vector_from->used*vector_from->element_size);
}
vector_t* vector_dup(vector_t* const vector_src) {
  vector_t* const vector_cpy = vector_new_(vector_src->used,vector_src->element_size);
  vector_set_used(vector_cpy,vector_src->used);
  memcpy(vector_cpy->memory,vector_src->memory,vector_src->used*vector_src->element_size);
  return vector_cpy;
}
