#include "ds.h"
MemoryTally gMemTally;
void MemoryTally::add(int cat, uint64_t amt)
{
	ThreadSafe ts(mutex_m);
	tots_[cat] += amt;
	tot_ += amt;
	if (tots_[cat] > peaks_[cat])
	{
		peaks_[cat] = tots_[cat];
	}
	if (tot_ > peak_)
	{
		peak_ = tot_;
	}
}
void MemoryTally::del(int cat, uint64_t amt)
{
	ThreadSafe ts(mutex_m);
	assert_geq(tots_[cat], amt);
	assert_geq(tot_, amt);
	tots_[cat] -= amt;
	tot_ -= amt;
}
#ifdef MAIN_DS
#include <limits>
#include "random_source.h"
using namespace std;
int main(void)
{
	cerr << "Test EHeap 1...";
	{
		EHeap<float> h;
		h.insert(0.5f);
		h.insert(0.6f);
		h.insert(0.25f);
		h.insert(0.75f);
		h.insert(0.1f);
		h.insert(0.9f);
		h.insert(0.4f);
		assert_eq(7, h.size());
		if (h.pop() != 0.1f)
		{
			throw 1;
		}
		assert_eq(6, h.size());
		if (h.pop() != 0.25f)
		{
			throw 1;
		}
		assert_eq(5, h.size());
		if (h.pop() != 0.4f)
		{
			throw 1;
		}
		assert_eq(4, h.size());
		if (h.pop() != 0.5f)
		{
			throw 1;
		}
		assert_eq(3, h.size());
		if (h.pop() != 0.6f)
		{
			throw 1;
		}
		assert_eq(2, h.size());
		if (h.pop() != 0.75f)
		{
			throw 1;
		}
		assert_eq(1, h.size());
		if (h.pop() != 0.9f)
		{
			throw 1;
		}
		assert_eq(0, h.size());
		assert(h.empty());
	}
	cerr << "PASSED" << endl;
	cerr << "Test EHeap 2...";
	{
		EHeap<size_t> h;
		RandomSource rnd(12);
		size_t lim = 2000;
		while (h.size() < lim)
		{
			h.insert(rnd.nextU32());
		}
		size_t last = std::numeric_limits<size_t>::max();
		bool first = true;
		while (!h.empty())
		{
			size_t p = h.pop();
			assert(first || p >= last);
			last = p;
			first = false;
		}
	}
	cerr << "PASSED" << endl;
	cerr << "Test EBitList 1...";
	{
		EBitList<128> l;
		assert_eq(0, l.size());
		assert_eq(std::numeric_limits<size_t>::max(), l.max());
		assert(!l.test(0));
		assert(!l.test(1));
		assert(!l.test(10));
		for (int i = 0; i < 3; i++)
		{
			l.set(10);
			assert(!l.test(0));
			assert(!l.test(1));
			assert(!l.test(9));
			assert(l.test(10));
			assert(!l.test(11));
		}
		assert_eq(10, l.max());
		l.clear();
		assert(!l.test(10));
		assert_eq(std::numeric_limits<size_t>::max(), l.max());
		RandomSource rnd(12);
		size_t lim = 2000;
		for (size_t i = 0; i < lim; i++)
		{
			uint32_t ri = rnd.nextU32() % 10000;
			l.set(ri);
			assert(l.test(ri));
		}
	}
	cerr << "PASSED" << endl;
}
#endif
