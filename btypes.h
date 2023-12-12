#ifndef _INDEX_TYPES_H
#define _INDEX_TYPES_H
#ifdef _64BIT_INDEX
#define OFF_MASK 0xffffffffffffffff
#define OFF_LEN_MASK 0xc000000000000000
#define LS_SIZE 0x100000000000000
#define OFF_SIZE 8
typedef uint64_t TIndexOffU;
typedef int64_t TIndexOff;
#else
#define OFF_MASK 0xffffffff
#define OFF_LEN_MASK 0xc0000000
#define LS_SIZE 0x10000000
#define OFF_SIZE 4
typedef uint32_t TIndexOffU;
typedef int TIndexOff;
#endif
extern const std::string gEbwt_ext;
#endif
