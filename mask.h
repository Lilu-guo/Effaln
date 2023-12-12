#ifndef MASK_H_
#define MASK_H_
#include <iostream>
#include "random_source.h"
extern int alts5[32];
extern int firsts5[32];
static inline int matchesEx(int i, int j)
{
	if (j >= 16 || i > 3)
	{
		return -1;
	}
	return (((1 << i) & j) != 0) ? 1 : 0;
}
static inline bool matches(int i, int j)
{
	return ((1 << i) & j) != 0;
}
static inline int randFromMask(RandomSource &rnd, int mask)
{
	assert_gt(mask, 0);
	if (alts5[mask] == 1)
	{
		return firsts5[mask];
	}
	assert_gt(mask, 0);
	assert_lt(mask, 32);
	int r = rnd.nextU32() % alts5[mask];
	assert_geq(r, 0);
	assert_lt(r, alts5[mask]);
	for (int i = 0; i < 5; i++)
	{
		if ((mask & (1 << i)) != 0)
		{
			if (r == 0)
				return i;
			r--;
		}
	}
	std::cerr << "Shouldn't get here" << std::endl;
	throw 1;
	return -1;
}
#endif
