#ifndef MM_ALLOCATOR_H_
#define MM_ALLOCATOR_H_
#include "../utils/vector.h"
typedef struct {
  uint64_t request_ticker;          
  uint64_t segment_size;          
  vector_t* segments;             
  vector_t* segments_free;        
  uint64_t current_segment_idx;     
  vector_t* malloc_requests;    
  } mm_allocator_t;
mm_allocator_t* mm_allocator_new(
    const uint64_t segment_size);
void mm_allocator_clear(
    mm_allocator_t* const mm_allocator);
void mm_allocator_delete(
    mm_allocator_t* const mm_allocator);
void* mm_allocator_allocate(
    mm_allocator_t* const mm_allocator,
    uint64_t num_bytes,
    const bool zero_mem
#ifdef MM_ALLOCATOR_LOG
    ,const char* func_name,
    uint64_t line_no
#endif
    );
#ifdef MM_ALLOCATOR_LOG
#define mm_allocator_alloc(mm_allocator,type) \
  ((type*)mm_allocator_allocate(mm_allocator,sizeof(type),false,__func__,(uint64_t)__LINE__))
#define mm_allocator_malloc(mm_allocator,num_bytes) \
  (mm_allocator_allocate(mm_allocator,(num_bytes),false,__func__,(uint64_t)__LINE__))
#define mm_allocator_calloc(mm_allocator,num_elements,type,clear_mem) \
  ((type*)mm_allocator_allocate(mm_allocator,(num_elements)*sizeof(type),clear_mem,__func__,(uint64_t)__LINE__))
#else
#define mm_allocator_alloc(mm_allocator,type) \
  ((type*)mm_allocator_allocate(mm_allocator,sizeof(type),false))
#define mm_allocator_malloc(mm_allocator,num_bytes) \
  (mm_allocator_allocate(mm_allocator,(num_bytes),false))
#define mm_allocator_calloc(mm_allocator,num_elements,type,clear_mem) \
  ((type*)mm_allocator_allocate(mm_allocator,(num_elements)*sizeof(type),clear_mem))
#endif
#define mm_allocator_uint64(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint64_t))
#define mm_allocator_uint32(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint32_t))
#define mm_allocator_uint16(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint16_t))
#define mm_allocator_uint8(mm_allocator)  mm_allocator_malloc(mm_allocator,sizeof(uint8_t))
void mm_allocator_free(
    mm_allocator_t* const mm_allocator,
    void* const memory);
void mm_allocator_get_occupation(
    mm_allocator_t* const mm_allocator,
    uint64_t* const bytes_used,
    uint64_t* const bytes_free_available,
    uint64_t* const bytes_free_fragmented);
void mm_allocator_print(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool display_requests);
#endif 
