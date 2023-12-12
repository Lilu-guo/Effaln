#ifndef ENDIAN_SWAP_H
#define ENDIAN_SWAP_H
#include <stdint.h>
#include <inttypes.h>
#include "assert_helpers.h"
#ifdef _64BIT_INDEX
#define endianSwapU(x) endianSwapU64(x)
#define endianSwapI(x) endianSwapI64(x)
#else
#define endianSwapU(x) endianSwapU32(x)
#define endianSwapI(x) endianSwapI32(x)
#endif
static inline bool currentlyBigEndian()
{
	static uint8_t endianCheck[] = {1, 0, 0, 0};
	return *((uint32_t *)endianCheck) != 1;
}
static inline uint32_t endianSwapU32(uint32_t u)
{
	uint32_t tmp = 0;
	tmp |= ((u >> 24) & (0xff << 0));
	tmp |= ((u >> 8) & (0xff << 8));
	tmp |= ((u << 8) & (0xff << 16));
	tmp |= ((u << 24) & (0xff << 24));
	return tmp;
}
static inline uint64_t endianSwapU64(uint64_t u)
{
	uint64_t tmp = 0;
	tmp |= ((u >> 56) & (0xffull << 0));
	tmp |= ((u >> 40) & (0xffull << 8));
	tmp |= ((u >> 24) & (0xffull << 16));
	tmp |= ((u >> 8) & (0xffull << 24));
	tmp |= ((u << 8) & (0xffull << 32));
	tmp |= ((u << 24) & (0xffull << 40));
	tmp |= ((u << 40) & (0xffull << 48));
	tmp |= ((u << 56) & (0xffull << 56));
	return tmp;
}
static inline int32_t endianSwapI32(int32_t i)
{
	int32_t tmp = 0;
	tmp |= ((i >> 24) & (0xff << 0));
	tmp |= ((i >> 8) & (0xff << 8));
	tmp |= ((i << 8) & (0xff << 16));
	tmp |= ((i << 24) & (0xff << 24));
	return tmp;
}
static inline int64_t endianSwapI64(int64_t u)
{
	int64_t tmp = 0;
	tmp |= ((u >> 56) & (0xffull << 0));
	tmp |= ((u >> 40) & (0xffull << 8));
	tmp |= ((u >> 24) & (0xffull << 16));
	tmp |= ((u >> 8) & (0xffull << 24));
	tmp |= ((u << 8) & (0xffull << 32));
	tmp |= ((u << 24) & (0xffull << 40));
	tmp |= ((u << 40) & (0xffull << 48));
	tmp |= ((u << 56) & (0xffull << 56));
	return tmp;
}
template <typename T>
static inline T endianizeU(T u, bool toBig)
{
	if (toBig == currentlyBigEndian())
	{
		return u;
	}
	if (sizeof(T) == 4)
	{
		return (T)endianSwapU32((uint32_t)u);
	}
	else if (sizeof(T) == 8)
	{
		return (T)endianSwapU64((uint64_t)u);
	}
	else
	{
		assert(false);
	}
}
template <typename T>
static inline T endianizeI(T i, bool toBig)
{
	if (toBig == currentlyBigEndian())
	{
		return i;
	}
	if (sizeof(T) == 4)
	{
		return endianSwapI32((int32_t)i);
	}
	else if (sizeof(T) == 8)
	{
		return endianSwapI64((int64_t)i);
	}
	else
	{
		assert(false);
	}
}
#endif
