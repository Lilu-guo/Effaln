#ifndef RANDOM_GEN_H_
#define RANDOM_GEN_H_
#include <stdint.h>
#include "assert_helpers.h"
#ifndef MERSENNE_TWISTER
class RandomSource
{
public:
	static const uint32_t DEFUALT_A = 1664525;
	static const uint32_t DEFUALT_C = 1013904223;
	RandomSource() : a(DEFUALT_A), c(DEFUALT_C), inited_(false) {}
	RandomSource(uint32_t _last) : a(DEFUALT_A), c(DEFUALT_C), last(_last), inited_(true) {}
	RandomSource(uint32_t _a, uint32_t _c) : a(_a), c(_c), inited_(false) {}
	void init(uint32_t seed = 0)
	{
		last = seed;
		inited_ = true;
		lastOff = 30;
	}
	uint32_t nextU32()
	{
		assert(inited_);
		uint32_t ret;
		last = a * last + c;
		ret = last >> 16;
		last = a * last + c;
		ret ^= last;
		lastOff = 0;
		return ret;
	}
	uint64_t nextU64()
	{
		assert(inited_);
		uint64_t first = nextU32();
		first = first << 32;
		uint64_t second = nextU32();
		return first | second;
	}
	size_t nextSizeT()
	{
		if (sizeof(size_t) == 4)
		{
			return nextU32();
		}
		else
		{
			return nextU64();
		}
	}
	uint32_t nextU32Range(uint32_t lo, uint32_t hi)
	{
		uint32_t ret = lo;
		if (hi > lo)
		{
			ret += (nextU32() % (hi - lo + 1));
		}
		return ret;
	}
	uint32_t nextU2()
	{
		assert(inited_);
		if (lastOff > 30)
		{
			nextU32();
		}
		uint32_t ret = (last >> lastOff) & 3;
		lastOff += 2;
		return ret;
	}
	bool nextBool()
	{
		assert(inited_);
		if (lastOff > 31)
		{
			nextU32();
		}
		uint32_t ret = (last >> lastOff) & 1;
		lastOff++;
		return ret;
	}
	uint32_t nextFromProbs(
		const float *weights,
		size_t numWeights)
	{
		float f = nextFloat();
		float tot = 0.0f;
		for (uint32_t i = 0; i < numWeights; i++)
		{
			tot += weights[i];
			if (f < tot)
				return i;
		}
		return (uint32_t)(numWeights - 1);
	}
	float nextFloat()
	{
		assert(inited_);
		return (float)nextU32() / (float)0xffffffff;
	}
	static uint32_t nextU32(uint32_t last,
							uint32_t a = DEFUALT_A,
							uint32_t c = DEFUALT_C)
	{
		return (a * last) + c;
	}
	uint32_t currentA() const { return a; }
	uint32_t currentC() const { return c; }
	uint32_t currentLast() const { return last; }
private:
	uint32_t a;
	uint32_t c;
	uint32_t last;
	uint32_t lastOff;
	bool inited_;
};
#else
class RandomSource
{
public:
	RandomSource()
	{
		reset();
	}
	RandomSource(uint32_t s)
	{
		init(s);
	}
	RandomSource(const uint32_t *array, int size)
	{
		init(array, size);
	}
	void reset()
	{
		state_[0] = 0;
		p_ = 0;
		inited_ = false;
	}
	virtual ~RandomSource() {}
	d init(const uint32_t *, int size);
	bool nextBool()
	{
		return (nextU32() & 1) == 0;
	}
	inline uint32_t nextU32()
	{
		assert(inited_);
		if (p_ == n)
		{
			gen_state();
		}
		uint32_t x = state_[p_++];
		x ^= (x >> 11);
		x ^= (x << 7) & 0x9D2C5680UL;
		x ^= (x << 15) & 0xEFC60000UL;
		x ^= (x >> 18);
		return x;
	}
	float nextFloat()
	{
		assert(inited_);
		return (float)nextU32() / (float)0xffffffff;
	}
protected:
	static const int n = 624, m = 397;
	uint32_t state_[n];
	int p_;
	bool inited_;
	uint32_t twiddle(uint32_t u, uint32_t v)
	{
		return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1) ^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
	}
	void gen_state();
};
#endif
#endif
