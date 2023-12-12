#include "mm_allocator.h"
#define MM_ALLOCATOR_SEGMENT_INITIAL_REQUESTS   10000
#define MM_ALLOCATOR_INITIAL_SEGMENTS              10
#define MM_ALLOCATOR_INITIAL_MALLOC_REQUESTS       10
#define MM_ALLOCATOR_INITIAL_STATES                10
#define MM_ALLOCATOR_FREED_FLAG                 0x80000000ul
#define MM_ALLOCATOR_REQUEST_IS_FREE(request)  ((request)->size & MM_ALLOCATOR_FREED_FLAG)
#define MM_ALLOCATOR_REQUEST_SET_FREE(request) ((request)->size |= MM_ALLOCATOR_FREED_FLAG)
#define MM_ALLOCATOR_REQUEST_SIZE(request)     ((request)->size & ~(MM_ALLOCATOR_FREED_FLAG))
typedef struct {
  uint32_t offset;
  uint32_t size;
#ifdef MM_ALLOCATOR_LOG
  uint64_t timestamp;
  char* func_name;
  uint64_t line_no;
#endif
} mm_allocator_request_t;
typedef struct {
  uint64_t segment_idx;             uint64_t segment_size;          void* memory;                   uint64_t used;                    vector_t* requests;           } mm_allocator_segment_t;
typedef struct {
  uint32_t segment_idx;
  uint32_t request_idx;
} mm_allocator_reference_t;
mm_allocator_segment_t* mm_allocator_segment_new(
    mm_allocator_t* const mm_allocator) {
  mm_allocator_segment_t* const segment = malloc(sizeof(mm_allocator_segment_t));
  const uint64_t segment_idx = vector_get_used(mm_allocator->segments);
  segment->segment_idx = segment_idx;
  segment->segment_size = mm_allocator->segment_size;
  segment->memory = malloc(mm_allocator->segment_size);
  segment->used = 0;
  segment->requests = vector_new(MM_ALLOCATOR_SEGMENT_INITIAL_REQUESTS,mm_allocator_request_t);
  vector_insert(mm_allocator->segments,segment,mm_allocator_segment_t*);
  return segment;
}
void mm_allocator_segment_clear(
    mm_allocator_segment_t* const segment) {
  segment->used = 0;
  vector_clear(segment->requests);
}
void mm_allocator_segment_delete(
    mm_allocator_segment_t* const segment) {
  vector_delete(segment->requests);
  free(segment->memory);
  free(segment);
}
mm_allocator_request_t* mm_allocator_segment_get_request(
    mm_allocator_segment_t* const segment,
    const uint64_t request_idx) {
  return vector_get_elm(segment->requests,request_idx,mm_allocator_request_t);
}
uint64_t mm_allocator_segment_get_num_requests(
    mm_allocator_segment_t* const segment) {
  return vector_get_used(segment->requests);
}
mm_allocator_t* mm_allocator_new(
    const uint64_t segment_size) {
  mm_allocator_t* const mm_allocator = malloc(sizeof(mm_allocator_t));   mm_allocator->request_ticker = 0;
  mm_allocator->segment_size = segment_size;   mm_allocator->current_segment_idx = 0;
  mm_allocator->segments = vector_new(MM_ALLOCATOR_INITIAL_SEGMENTS,mm_allocator_segment_t*);
  mm_allocator->segments_free = vector_new(MM_ALLOCATOR_INITIAL_SEGMENTS,mm_allocator_segment_t*);
  mm_allocator_segment_new(mm_allocator);
  mm_allocator->malloc_requests = vector_new(MM_ALLOCATOR_INITIAL_MALLOC_REQUESTS,void*);
  return mm_allocator;
}
void mm_allocator_clear(
    mm_allocator_t* const mm_allocator) {
  vector_clear(mm_allocator->segments_free);
  const uint64_t num_segments = vector_get_used(mm_allocator->segments);
  mm_allocator_segment_t** const segments = 
      vector_get_mem(mm_allocator->segments,mm_allocator_segment_t*);
  uint64_t i;
  for (i=1;i<num_segments;++i) {
    mm_allocator_segment_clear(segments[i]); 
    vector_insert(mm_allocator->segments_free,segments[i],mm_allocator_segment_t*);   }
  mm_allocator->current_segment_idx = 0;
  VECTOR_ITERATE(mm_allocator->malloc_requests,malloc_request,m,void*) {
    free(*malloc_request); 
  }
  vector_clear(mm_allocator->malloc_requests);
}
void mm_allocator_delete(
    mm_allocator_t* const mm_allocator) {
  VECTOR_ITERATE(mm_allocator->segments,segment_ptr,p,mm_allocator_segment_t*) {
    mm_allocator_segment_delete(*segment_ptr);
  }
  vector_delete(mm_allocator->segments);
  vector_delete(mm_allocator->segments_free);
  VECTOR_ITERATE(mm_allocator->malloc_requests,malloc_request,m,void*) {
    free(*malloc_request); 
  }
  vector_delete(mm_allocator->malloc_requests);
  free(mm_allocator);
}
mm_allocator_segment_t* mm_allocator_get_segment(
    mm_allocator_t* const mm_allocator,
    const uint64_t segment_idx) {
  return *(vector_get_elm(mm_allocator->segments,segment_idx,mm_allocator_segment_t*));
}
uint64_t mm_allocator_get_num_segments(
    mm_allocator_t* const mm_allocator) {
  return vector_get_used(mm_allocator->segments);
}
mm_allocator_segment_t* mm_allocator_fetch_segment(
    mm_allocator_t* const mm_allocator,
    const uint64_t num_bytes) {
  mm_allocator_segment_t* const curr_segment =
      mm_allocator_get_segment(mm_allocator,mm_allocator->current_segment_idx);
  if (curr_segment->used + num_bytes <= curr_segment->segment_size) {
    return curr_segment;
  }
  if (num_bytes > curr_segment->segment_size) {
    return NULL; 
  }
  const uint64_t free_segments = vector_get_used(mm_allocator->segments_free);
  if (free_segments > 0) {
    mm_allocator_segment_t* const segment =
        *vector_get_elm(mm_allocator->segments_free,free_segments-1,mm_allocator_segment_t*);
    vector_dec_used(mm_allocator->segments_free);
    mm_allocator->current_segment_idx = segment->segment_idx;
    return segment;
  }
  mm_allocator_segment_t* const segment = mm_allocator_segment_new(mm_allocator);
  mm_allocator->current_segment_idx = segment->segment_idx;
  return segment;
}
void* mm_allocator_allocate(
    mm_allocator_t* const mm_allocator,
    uint64_t num_bytes,
    const bool zero_mem
#ifdef MM_ALLOCATOR_LOG
    ,const char* func_name,
    uint64_t line_no
#endif
    ) {
  if (num_bytes == 0) {
    fprintf(stderr,"MM-Allocator error. Zero bytes request\n");
    exit(1);
  }
  num_bytes += sizeof(mm_allocator_reference_t);
  if (num_bytes%16 != 0) { 
    num_bytes += 16 - (num_bytes%16);   }
#ifdef MM_ALLOCATOR_MALLOC
  mm_allocator_segment_t* const segment = NULL; 
#else
  mm_allocator_segment_t* const segment = mm_allocator_fetch_segment(mm_allocator,num_bytes);
#endif
  if (segment != NULL) {
    void* const memory = segment->memory + segment->used;
    if (zero_mem) memset(memory,0,num_bytes); 
        mm_allocator_reference_t* const mm_reference = memory;
    mm_reference->segment_idx = segment->segment_idx;
    mm_reference->request_idx = mm_allocator_segment_get_num_requests(segment);
    mm_allocator_request_t* request;
    vector_alloc_new(segment->requests,mm_allocator_request_t,request);
    request->offset = segment->used;
    request->size = num_bytes;
#ifdef MM_ALLOCATOR_LOG
    request->timestamp = (mm_allocator->request_ticker)++;
    request->func_name = (char*)func_name;
    request->line_no = line_no;
#endif
    segment->used += num_bytes;
    return memory + sizeof(mm_allocator_reference_t);
  } else {
    void* const memory = malloc(num_bytes);
    if (zero_mem) memset(memory,0,num_bytes); 
    vector_insert(mm_allocator->malloc_requests,memory,void*);
    mm_allocator_reference_t* const mm_reference = memory;
    mm_reference->segment_idx = UINT32_MAX;
    return memory + sizeof(mm_allocator_reference_t);
  }
}
void mm_allocator_free_malloc_request(
    mm_allocator_t* const mm_allocator,
    void* memory) {
  const uint64_t num_malloc_requests = vector_get_used(mm_allocator->malloc_requests);
  void** const malloc_requests = vector_get_mem(mm_allocator->malloc_requests,void*);
  uint64_t i;
  for (i=0;i<num_malloc_requests;++i) { 
    if (malloc_requests[i] == memory) {
      free(memory);
      for (;i<num_malloc_requests-1;++i) {
        malloc_requests[i] = malloc_requests[i+1];
      }
      vector_dec_used(mm_allocator->malloc_requests);
      return;
    }
  }
  fprintf(stderr,"MM-Allocator error. Invalid address freed (request not found)\n");
  exit(1);
}
void mm_allocator_free_allocator_request(
    mm_allocator_t* const mm_allocator,
    mm_allocator_segment_t* const segment,
    const uint32_t request_idx,
    mm_allocator_request_t* const request) {
  if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
    fprintf(stderr,"MM-Allocator error: double free\n");
    exit(1);
  }
  MM_ALLOCATOR_REQUEST_SET_FREE(request);
  uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
  if (request_idx == num_requests-1) { 
    --num_requests;
    mm_allocator_request_t* request =
        vector_get_mem(segment->requests,mm_allocator_request_t) + (num_requests-1);
    while (num_requests>0 && MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
      --num_requests; 
      --request;
    }
    if (num_requests > 0) {
      segment->used = request->offset + request->size;
      vector_set_used(segment->requests,num_requests);
    } else {
      mm_allocator_segment_clear(segment);             if (segment->segment_idx != mm_allocator->current_segment_idx) {
        vector_insert(mm_allocator->segments_free,segment,mm_allocator_segment_t*);
      }
    }
  }
}
void mm_allocator_free(
    mm_allocator_t* const mm_allocator,
    void* const memory) {
  void* const effective_memory = memory - sizeof(mm_allocator_reference_t);
  mm_allocator_reference_t* const mm_reference = effective_memory;
  if (mm_reference->segment_idx == UINT32_MAX) {
    mm_allocator_free_malloc_request(mm_allocator,effective_memory);
  } else {
    mm_allocator_segment_t* const segment =
        mm_allocator_get_segment(mm_allocator,mm_reference->segment_idx);
    mm_allocator_request_t* const request =
        mm_allocator_segment_get_request(segment,mm_reference->request_idx);
    mm_allocator_free_allocator_request(
        mm_allocator,segment,mm_reference->request_idx,request);
  }
}
void mm_allocator_get_occupation(
    mm_allocator_t* const mm_allocator,
    uint64_t* const bytes_used,
    uint64_t* const bytes_free_available,
    uint64_t* const bytes_free_fragmented) {
  *bytes_used = 0;
  *bytes_free_available = 0;
  *bytes_free_fragmented = 0;
  const uint64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  int64_t segment_idx, request_idx;
  for (segment_idx=0;segment_idx<num_segments;++segment_idx) {
    mm_allocator_segment_t* const segment = mm_allocator_get_segment(mm_allocator,segment_idx);
    const uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
    bool free_memory = true;
    for (request_idx=num_requests-1;request_idx>=0;--request_idx) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,request_idx);
      const uint64_t size = MM_ALLOCATOR_REQUEST_SIZE(request);
      if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
        if (free_memory) {
          *bytes_free_available += size;
        } else {
          *bytes_free_fragmented += size;
        }
      } else {
        free_memory = false;
        *bytes_used += size;
      }
    }
    if (num_requests > 0) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,num_requests-1);
      *bytes_free_available += segment->used - (request->offset+request->size);
    }
  }
}
void mm_allocator_print_request(
    FILE* const stream,
    mm_allocator_request_t* const request,
    const uint64_t segment_idx,
    const uint64_t request_idx) {
      fprintf(stream,"    [#%03" PRIu64 "/%05" PRIu64 "\t%s\t@%08u\t(%" PRIu64 " Bytes)"
#ifdef MM_ALLOCATOR_LOG
          "\t%s:%" PRIu64 "\t{ts=%" PRIu64 "}"
#endif
          "\n",
          segment_idx,
          request_idx,
          MM_ALLOCATOR_REQUEST_IS_FREE(request) ? "Free]     " : "Allocated]",
          request->offset,
          (uint64_t)MM_ALLOCATOR_REQUEST_SIZE(request)
#ifdef MM_ALLOCATOR_LOG
          ,request->func_name,
          request->line_no,
          request->timestamp
#endif
      );
}
void mm_allocator_print_requests(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool compact_free) {
  uint64_t segment_idx, request_idx;
  uint64_t free_block = 0;
  fprintf(stream,"  => Memory.requests\n");
  const uint64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  for (segment_idx=0;segment_idx<num_segments;++segment_idx) {
    mm_allocator_segment_t* const segment = mm_allocator_get_segment(mm_allocator,segment_idx);
    const uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
    for (request_idx=0;request_idx<num_requests;++request_idx) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,request_idx);
      if (compact_free) {
        if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
          free_block += MM_ALLOCATOR_REQUEST_SIZE(request);
        } else {
          if (free_block > 0) {
            fprintf(stream,"    [n/a\tFree]      \t(%" PRIu64 " Bytes)\n",free_block);
            free_block = 0;
          }
          mm_allocator_print_request(stream,request,segment_idx,request_idx);
        }
      } else {
        mm_allocator_print_request(stream,request,segment_idx,request_idx);
      }
    }
    if (request_idx==0) {
      fprintf(stream,"    -- No requests --\n");
    }
  }
}
void mm_allocator_print(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool display_requests) {
  fprintf(stream,"MM-Allocator.report\n");
  const uint64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  const uint64_t segment_size = mm_allocator_get_segment(mm_allocator,0)->segment_size;
  fprintf(stream,"  => Segments.allocated %" PRIu64 "\n",num_segments);
  fprintf(stream,"    => Segments.size %" PRIu64 " MB\n",segment_size/(1024*1024));
  fprintf(stream,"    => Memory.available %" PRIu64 " MB\n",num_segments*(segment_size/(1024*1024)));
  uint64_t bytes_used, bytes_free_available, bytes_free_fragmented;
  mm_allocator_get_occupation(mm_allocator,&bytes_used,&bytes_free_available,&bytes_free_fragmented);
  fprintf(stream,"  => Memory.used %" PRIu64 "\n",bytes_used);
  fprintf(stream,"  => Memory.free %" PRIu64 "\n",bytes_free_available+bytes_free_fragmented);
  fprintf(stream,"    => Memory.free.available  %" PRIu64 "\n",bytes_free_available);
  fprintf(stream,"    => Memory.free.fragmented %" PRIu64 "\n",bytes_free_fragmented);
  if (display_requests) {
    mm_allocator_print_requests(stream,mm_allocator,false);
  }
}
