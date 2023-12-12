#ifndef BITPACK_H_
#define BITPACK_H_
#include <stdint.h>
#include "assert_helpers.h"
static inline void pack_2b_in_8b(const int two, uint8_t &eight, const int off)
{
	assert_lt(two, 4);
	assert_lt(off, 4);
	eight |= (two << (off * 2));
}
static inline int unpack_2b_from_8b(const uint8_t eight, const int off)
{
	assert_lt(off, 4);
	return ((eight >> (off * 2)) & 0x3);
}
static inline void pack_1b_in_8b(const int one, uint8_t &eight, const int off)
{
	assert_lt(one, 2);
	assert_lt(off, 8);
	eight |= (one << off);
}
static inline int unpack_1b_from_8b(const uint8_t eight, const int off)
{
	assert_lt(off, 8);
	return ((eight >> off) & 0x1);
}
static inline void pack_2b_in_32b(const int two, uint32_t &thirty2, const int off)
{
	assert_lt(two, 4);
	assert_lt(off, 16);
	thirty2 |= (two << (off * 2));
}
static inline int unpack_2b_from_32b(const uint32_t thirty2, const int off)
{
	assert_lt(off, 16);
	return ((thirty2 >> (off * 2)) & 0x3);
}
#endif
