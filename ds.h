#ifndef DS_H_
#define DS_H_
#include <chrono>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <stdint.h>
#include <string.h>
#include <limits>
#include "assert_helpers.h"
#include "threading.h"
#include "random_source.h"
#include "btypes.h"
class MemoryTally
{
public:
	MemoryTally() : tot_(0), peak_(0)
	{
		memset(tots_, 0, 256 * sizeof(uint64_t));
		memset(peaks_, 0, 256 * sizeof(uint64_t));
	}
	void add(int cat, uint64_t amt);
	void del(int cat, uint64_t amt);
	uint64_t total() { return tot_; }
	uint64_t total(int cat) { return tots_[cat]; }
	uint64_t peak() { return peak_; }
	uint64_t peak(int cat) { return peaks_[cat]; }
#ifndef NDEBUG
	bool repOk() const
	{
		uint64_t tot = 0;
		for (int i = 0; i < 256; i++)
		{
			assert_leq(tots_[i], peaks_[i]);
			tot += tots_[i];
		}
		assert_eq(tot, tot_);
		return true;
	}
#endif
protected:
	MUTEX_T mutex_m;
	uint64_t tots_[256];
	uint64_t tot_;
	uint64_t peaks_[256];
	uint64_t peak_;
};
#ifdef USE_MEM_TALLY
extern MemoryTally gMemTally;
#endif
template <typename T>
class AutoArray
{
public:
	AutoArray(size_t sz, int cat = 0) : cat_(cat)
	{
		t_ = NULL;
		t_ = new T[sz];
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#else
		(void)cat_;
#endif
		memset(t_, 0, sz * sizeof(T));
		sz_ = sz;
	}
	~AutoArray()
	{
		if (t_ != NULL)
		{
			delete[] t_;
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, sz_);
#endif
		}
	}
	T &operator[](size_t sz)
	{
		return t_[sz];
	}
	const T &operator[](size_t sz) const
	{
		return t_[sz];
	}
	size_t size() const { return sz_; }
private:
	int cat_;
	T *t_;
	size_t sz_;
};
template <typename T>
class PtrWrap
{
public:
	explicit PtrWrap(
		T *p,
		bool freeable = true,
		int cat = 0) : cat_(cat),
					   p_(NULL)
	{
		init(p, freeable);
	}
	explicit PtrWrap(int cat = 0) : cat_(cat),
									p_(NULL)
	{
		reset();
	}
	void reset()
	{
		free();
		init(NULL);
	}
	~PtrWrap() { free(); }
	void init(T *p, bool freeable = true)
	{
		assert(p_ == NULL);
		p_ = p;
		freeable_ = freeable;
#ifdef USE_MEM_TALLY
		if (p != NULL && freeable_)
		{
			gMemTally.add(cat_, sizeof(T));
		}
#else
		(void)cat_;
#endif
	}
	void free()
	{
		if (p_ != NULL)
		{
			if (freeable_)
			{
				delete p_;
#ifdef USE_MEM_TALLY
				gMemTally.del(cat_, sizeof(T));
#endif
			}
			p_ = NULL;
		}
	}
	inline T *get() { return p_; }
	inline const T *get() const { return p_; }
private:
	int cat_;
	T *p_;
	bool freeable_;
};
template <typename T>
class APtrWrap
{
public:
	explicit APtrWrap(
		T *p,
		size_t sz,
		bool freeable = true,
		int cat = 0) : cat_(cat),
					   p_(NULL)
	{
		init(p, sz, freeable);
	}
	explicit APtrWrap(int cat = 0) : cat_(cat),
									 p_(NULL)
	{
		reset();
	}
	void reset()
	{
		free();
		init(NULL, 0);
	}
	~APtrWrap() { free(); }
	void init(T *p, size_t sz, bool freeable = true)
	{
		assert(p_ == NULL);
		p_ = p;
		sz_ = sz;
		freeable_ = freeable;
#ifdef USE_MEM_TALLY
		if (p != NULL && freeable_)
		{
			gMemTally.add(cat_, sizeof(T) * sz_);
		}
#else
		(void)cat_;
#endif
	}
	void free()
	{
		if (p_ != NULL)
		{
			if (freeable_)
			{
				delete[] p_;
#ifdef USE_MEM_TALLY
				gMemTally.del(cat_, sizeof(T) * sz_);
#endif
			}
			p_ = NULL;
		}
	}
	inline T *get() { return p_; }
	inline const T *get() const { return p_; }
private:
	int cat_;
	T *p_;
	bool freeable_;
	size_t sz_;
};
template <typename T, int S = 128>
class EList
{
public:
	explicit EList() : cat_(0), allocCat_(-1), list_(NULL), sz_(S), cur_(0)
	{
#ifndef USE_MEM_TALLY
		(void)cat_;
#endif
	}
	explicit EList(int cat) : cat_(cat), allocCat_(-1), list_(NULL), sz_(S), cur_(0)
	{
		assert_geq(cat, 0);
	}
	explicit EList(size_t isz, int cat = 0) : cat_(cat), allocCat_(-1), list_(NULL), sz_(isz), cur_(0)
	{
		assert_geq(cat, 0);
	}
	EList(const EList<T, S> &o) : cat_(0), allocCat_(-1), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
	}
	explicit EList(const EList<T, S> &o, int cat) : cat_(cat), allocCat_(-1), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
		assert_geq(cat, 0);
	}
	~EList() { free(); }
	EList<T, S> &operator=(const EList<T, S> &o)
	{
		assert_eq(cat_, o.cat());
		if (o.cur_ == 0)
		{
			cur_ = 0;
			return *this;
		}
		if (list_ == NULL)
		{
			lazyInit();
		}
		if (sz_ < o.cur_)
			expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for (size_t i = 0; i < cur_; i++)
		{
			list_[i] = o.list_[i];
		}
		return *this;
	}
	void xfer(EList<T, S> &o)
	{
		assert_eq(cat_, o.cat());
		free();
		allocCat_ = cat_;
		list_ = o.list_;
		sz_ = o.sz_;
		cur_ = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
		o.allocCat_ = -1;
	}
	inline size_t size() const { return cur_; }
	inline size_t capacity() const { return sz_; }
	size_t totalSizeBytes() const
	{
		return 2 * sizeof(int) +
			   2 * sizeof(size_t) +
			   cur_ * sizeof(T);
	}
	size_t totalCapacityBytes() const
	{
		return 2 * sizeof(int) +
			   2 * sizeof(size_t) +
			   sz_ * sizeof(T);
	}
	inline void ensure(size_t thresh)
	{
		if (list_ == NULL)
			lazyInit();
		expandCopy(cur_ + thresh);
	}
	inline void reserveExact(size_t newsz)
	{
		if (list_ == NULL)
			lazyInitExact(newsz);
		expandCopyExact(newsz);
	}
	inline bool empty() const { return cur_ == 0; }
	inline bool null() const { return list_ == NULL; }
	void push_back(const T &el)
	{
		if (list_ == NULL)
			lazyInit();
		if (cur_ == sz_)
			expandCopy(sz_ + 1);
		list_[cur_++] = el;
	}
	void expand()
	{
		if (list_ == NULL)
			lazyInit();
		if (cur_ == sz_)
			expandCopy(sz_ + 1);
		cur_++;
	}
	void fill(size_t begin, size_t end, const T &v)
	{
		assert_leq(begin, end);
		assert_leq(end, cur_);
		for (size_t i = begin; i < end; i++)
		{
			list_[i] = v;
		}
	}
	void fill(const T &v)
	{
		for (size_t i = 0; i < cur_; i++)
		{
			list_[i] = v;
		}
	}
	void fillZero(size_t begin, size_t end)
	{
		assert_leq(begin, end);
		memset(&list_[begin], 0, sizeof(T) * (end - begin));
	}
	void fillZero()
	{
		memset(list_, 0, sizeof(T) * cur_);
	}
	void resizeNoCopy(size_t sz)
	{
		if (sz > 0 && list_ == NULL)
			lazyInit();
		if (sz <= cur_)
		{
			cur_ = sz;
			return;
		}
		if (sz_ < sz)
			expandNoCopy(sz);
		cur_ = sz;
	}
	void resize(size_t sz)
	{
		if (sz > 0 && list_ == NULL)
			lazyInit();
		if (sz <= cur_)
		{
			cur_ = sz;
			return;
		}
		if (sz_ < sz)
		{
			expandCopy(sz);
		}
		cur_ = sz;
	}
	void resizeExact(size_t sz)
	{
		if (sz > 0 && list_ == NULL)
			lazyInitExact(sz);
		if (sz <= cur_)
		{
			cur_ = sz;
			return;
		}
		if (sz_ < sz)
			expandCopyExact(sz);
		cur_ = sz;
	}
	void erase(size_t idx)
	{
		assert_lt(idx, cur_);
		for (size_t i = idx; i < cur_ - 1; i++)
		{
			list_[i] = list_[i + 1];
		}
		cur_--;
	}
	void erase(size_t idx, size_t len)
	{
		assert_geq(len, 0);
		if (len == 0)
		{
			return;
		}
		assert_lt(idx, cur_);
		for (size_t i = idx; i < cur_ - len; i++)
		{
			list_[i] = list_[i + len];
		}
		cur_ -= len;
	}
	void insert(const T &el, size_t idx)
	{
		if (list_ == NULL)
			lazyInit();
		assert_leq(idx, cur_);
		if (cur_ == sz_)
			expandCopy(sz_ + 1);
		for (size_t i = cur_; i > idx; i--)
		{
			list_[i] = list_[i - 1];
		}
		list_[idx] = el;
		cur_++;
	}
	void insert(const EList<T> &l, size_t idx)
	{
		if (list_ == NULL)
			lazyInit();
		assert_lt(idx, cur_);
		if (l.cur_ == 0)
			return;
		if (cur_ + l.cur_ > sz_)
			expandCopy(cur_ + l.cur_);
		for (size_t i = cur_ + l.cur_ - 1; i > idx + (l.cur_ - 1); i--)
		{
			list_[i] = list_[i - l.cur_];
		}
		for (size_t i = 0; i < l.cur_; i++)
		{
			list_[i + idx] = l.list_[i];
		}
		cur_ += l.cur_;
	}
	void pop_back()
	{
		assert_gt(cur_, 0);
		cur_--;
	}
	void clear()
	{
		cur_ = 0;
	}
	inline T &back()
	{
		assert_gt(cur_, 0);
		return list_[cur_ - 1];
	}
	void reverse()
	{
		if (cur_ > 1)
		{
			size_t n = cur_ >> 1;
			for (size_t i = 0; i < n; i++)
			{
				T tmp = list_[i];
				list_[i] = list_[cur_ - i - 1];
				list_[cur_ - i - 1] = tmp;
			}
		}
	}
	inline const T &back() const
	{
		assert_gt(cur_, 0);
		return list_[cur_ - 1];
	}
	inline T &front()
	{
		assert_gt(cur_, 0);
		return list_[0];
	}
	inline const T &front() const { return front(); }
	bool operator==(const EList<T, S> &o) const
	{
		if (size() != o.size())
		{
			return false;
		}
		for (size_t i = 0; i < size(); i++)
		{
			if (!(get(i) == o.get(i)))
			{
				return false;
			}
		}
		return true;
	}
	bool isSuperset(const EList<T, S> &o) const
	{
		if (o.size() > size())
		{
			return false;
		}
		for (size_t i = 0; i < o.size(); i++)
		{
			bool inthis = false;
			for (size_t j = 0; j < size(); j++)
			{
				if (o[i] == (*this)[j])
				{
					inthis = true;
					break;
				}
			}
			if (!inthis)
			{
				return false;
			}
		}
		return true;
	}
	inline T &operator[](size_t i)
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline const T &operator[](size_t i) const
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline T &get(size_t i)
	{
		return operator[](i);
	}
	inline const T &get(size_t i) const
	{
		return operator[](i);
	}
	T &getSlow(size_t i)
	{
		return operator[](i);
	}
	const T &getSlow(size_t i) const
	{
		return operator[](i);
	}
	void sortPortion(size_t begin, size_t num)
	{
		assert_leq(begin + num, cur_);
		if (num < 2)
			return;
		std::stable_sort(list_ + begin, list_ + begin + num);
	}
	void shufflePortion(size_t begin, size_t num, RandomSource &rnd)
	{
		assert_leq(begin + num, cur_);
		if (num < 2)
			return;
		size_t left = num;
		for (size_t i = begin; i < begin + num - 1; i++)
		{
			size_t rndi = rnd.nextSizeT() % left;
			if (rndi > 0)
			{
				std::swap(list_[i], list_[i + rndi]);
			}
			left--;
		}
	}
	void sort()
	{
		sortPortion(0, cur_);
	}
	bool sorted() const
	{
		for (size_t i = 1; i < cur_; i++)
		{
			if (!(list_[i - 1] < list_[i]))
			{
				return false;
			}
		}
		return true;
	}
	void remove(size_t idx)
	{
		assert_lt(idx, cur_);
		assert_gt(cur_, 0);
		for (size_t i = idx; i < cur_ - 1; i++)
		{
			list_[i] = list_[i + 1];
		}
		cur_--;
	}
	T *ptr() { return list_; }
	const T *ptr() const { return list_; }
	void setCat(int cat)
	{
		assert(null());
		assert_gt(cat, 0);
		cat_ = cat;
	}
	int cat() const { return cat_; }
	size_t bsearchLoBound(const T &el) const
	{
		size_t hi = cur_;
		size_t lo = 0;
		while (true)
		{
			if (lo == hi)
			{
				return lo;
			}
			size_t mid = lo + ((hi - lo) >> 1);
			assert_neq(mid, hi);
			if (list_[mid] < el)
			{
				if (lo == mid)
				{
					return hi;
				}
				lo = mid;
			}
			else
			{
				hi = mid;
			}
		}
	}
private:
	void lazyInit()
	{
		assert(list_ == NULL);
		list_ = alloc(sz_);
	}
	void lazyInitExact(size_t sz)
	{
		assert_gt(sz, 0);
		assert(list_ == NULL);
		sz_ = sz;
		list_ = alloc(sz);
	}
	T *alloc(size_t sz)
	{
		T *tmp = new T[sz];
		assert(tmp != NULL);
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#endif
		allocCat_ = cat_;
		return tmp;
	}
	void free()
	{
		if (list_ != NULL)
		{
			assert_neq(-1, allocCat_);
			assert_eq(allocCat_, cat_);
			delete[] list_;
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, sz_);
#endif
			list_ = NULL;
			sz_ = cur_ = 0;
		}
	}
	void expandCopy(size_t thresh)
	{
		if (thresh <= sz_)
			return;
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
			newsz *= 2;
		expandCopyExact(newsz);
	}
	void expandCopyExact(size_t newsz)
	{
		if (newsz <= sz_)
			return;
		T *tmp = alloc(newsz);
		assert(tmp != NULL);
		size_t cur = cur_;
		if (list_ != NULL)
		{
			for (size_t i = 0; i < cur_; i++)
			{
				tmp[i] = list_[i];
			}
			free();
		}
		list_ = tmp;
		sz_ = newsz;
		cur_ = cur;
	}
	void expandNoCopy(size_t thresh)
	{
		assert(list_ != NULL);
		if (thresh <= sz_)
			return;
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
			newsz *= 2;
		expandNoCopyExact(newsz);
	}
	void expandNoCopyExact(size_t newsz)
	{
		assert(list_ != NULL);
		assert_gt(newsz, 0);
		free();
		T *tmp = alloc(newsz);
		assert(tmp != NULL);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}
	int cat_;
	int allocCat_;
	T *list_;
	size_t sz_;
	size_t cur_;
};
template <typename T, int S1 = 128, int S2 = 128>
class ELList
{
public:
	explicit ELList(int cat = 0) : cat_(cat), list_(NULL), sz_(S2), cur_(0)
	{
#ifndef USE_MEM_TALLY
		(void)cat_;
#endif
		assert_geq(cat, 0);
	}
	explicit ELList(size_t isz, int cat = 0) : cat_(cat), list_(NULL), sz_(isz), cur_(0)
	{
		assert_gt(isz, 0);
		assert_geq(cat, 0);
	}
	ELList(const ELList<T, S1, S2> &o) : cat_(0), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
	}
	explicit ELList(const ELList<T, S1, S2> &o, int cat) : cat_(cat), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
		assert_geq(cat, 0);
	}
	~ELList() { free(); }
	ELList<T, S1, S2> &operator=(const ELList<T, S1, S2> &o)
	{
		assert_eq(cat_, o.cat());
		if (list_ == NULL)
		{
			lazyInit();
		}
		if (o.cur_ == 0)
		{
			cur_ = 0;
			return *this;
		}
		if (sz_ < o.cur_)
			expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for (size_t i = 0; i < cur_; i++)
		{
			assert_eq(list_[i].cat(), o.list_[i].cat());
			list_[i] = o.list_[i];
		}
		return *this;
	}
	void xfer(ELList<T, S1, S2> &o)
	{
		assert_eq(cat_, o.cat());
		list_ = o.list_;
		sz_ = o.sz_;
		cur_ = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
	}
	inline size_t size() const { return cur_; }
	inline bool empty() const { return cur_ == 0; }
	inline bool null() const { return list_ == NULL; }
	void expand()
	{
		if (list_ == NULL)
			lazyInit();
		if (cur_ == sz_)
			expandCopy(sz_ + 1);
		cur_++;
	}
	void resize(size_t sz)
	{
		if (sz > 0 && list_ == NULL)
			lazyInit();
		if (sz <= cur_)
		{
			cur_ = sz;
			return;
		}
		if (sz_ < sz)
		{
			expandCopy(sz);
		}
		cur_ = sz;
	}
	void clear()
	{
		cur_ = 0;
	}
	inline EList<T, S1> &back()
	{
		assert_gt(cur_, 0);
		return list_[cur_ - 1];
	}
	inline const EList<T, S1> &back() const
	{
		assert_gt(cur_, 0);
		return list_[cur_ - 1];
	}
	inline EList<T, S1> &front()
	{
		assert_gt(cur_, 0);
		return list_[0];
	}
	inline const EList<T, S1> &front() const { return front(); }
	inline EList<T, S1> &operator[](size_t i)
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline const EList<T, S1> &operator[](size_t i) const
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline EList<T, S1> &get(size_t i)
	{
		return operator[](i);
	}
	inline const EList<T, S1> &get(size_t i) const
	{
		return operator[](i);
	}
	EList<T, S1> &getSlow(size_t i)
	{
		return operator[](i);
	}
	const EList<T, S1> &getSlow(size_t i) const
	{
		return operator[](i);
	}
	EList<T, S1> *ptr() { return list_; }
	void setCat(int cat)
	{
		assert_gt(cat, 0);
		cat_ = cat;
		if (cat_ != 0)
		{
			for (size_t i = 0; i < sz_; i++)
			{
				assert(list_[i].null());
				list_[i].setCat(cat_);
			}
		}
	}
	int cat() const { return cat_; }
protected:
	void lazyInit()
	{
		assert(list_ == NULL);
		list_ = alloc(sz_);
	}
	EList<T, S1> *alloc(size_t sz)
	{
		assert_gt(sz, 0);
		EList<T, S1> *tmp = new EList<T, S1>[sz];
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#endif
		if (cat_ != 0)
		{
			for (size_t i = 0; i < sz; i++)
			{
				assert(tmp[i].ptr() == NULL);
				tmp[i].setCat(cat_);
			}
		}
		return tmp;
	}
	void free()
	{
		if (list_ != NULL)
		{
			delete[] list_;
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, sz_);
#endif
			list_ = NULL;
		}
	}
	void expandCopy(size_t thresh)
	{
		assert(list_ != NULL);
		if (thresh <= sz_)
			return;
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
			newsz *= 2;
		EList<T, S1> *tmp = alloc(newsz);
		if (list_ != NULL)
		{
			for (size_t i = 0; i < cur_; i++)
			{
				assert_eq(cat_, tmp[i].cat());
				tmp[i].xfer(list_[i]);
				assert_eq(cat_, tmp[i].cat());
			}
			free();
		}
		list_ = tmp;
		sz_ = newsz;
	}
	void expandNoCopy(size_t thresh)
	{
		assert(list_ != NULL);
		if (thresh <= sz_)
			return;
		free();
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
			newsz *= 2;
		EList<T, S1> *tmp = alloc(newsz);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}
	int cat_;
	EList<T, S1> *list_;
	size_t sz_;
	size_t cur_;
};
template <typename T, int S1 = 128, int S2 = 128, int S3 = 128>
class ELLList
{
public:
	explicit ELLList(int cat = 0) : cat_(cat), list_(NULL), sz_(S3), cur_(0)
	{
		assert_geq(cat, 0);
	}
	explicit ELLList(size_t isz, int cat = 0) : cat_(cat), list_(NULL), sz_(isz), cur_(0)
	{
		assert_geq(cat, 0);
		assert_gt(isz, 0);
	}
	ELLList(const ELLList<T, S1, S2, S3> &o) : cat_(0), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
	}
	explicit ELLList(const ELLList<T, S1, S2, S3> &o, int cat) : cat_(cat), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
		assert_geq(cat, 0);
	}
	~ELLList() { free(); }
	ELLList<T, S1, S2, S3> &operator=(const ELLList<T, S1, S2, S3> &o)
	{
		assert_eq(cat_, o.cat());
		if (list_ == NULL)
			lazyInit();
		if (o.cur_ == 0)
		{
			cur_ = 0;
			return *this;
		}
		if (sz_ < o.cur_)
			expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for (size_t i = 0; i < cur_; i++)
		{
			assert_eq(list_[i].cat(), o.list_[i].cat());
			list_[i] = o.list_[i];
		}
		return *this;
	}
	void xfer(ELLList<T, S1, S2, S3> &o)
	{
		assert_eq(cat_, o.cat());
		list_ = o.list_;
		sz_ = o.sz_;
		cur_ = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
	}
	inline size_t size() const { return cur_; }
	inline bool empty() const { return cur_ == 0; }
	inline bool null() const { return list_ == NULL; }
	void expand()
	{
		if (list_ == NULL)
			lazyInit();
		if (cur_ == sz_)
			expandCopy(sz_ + 1);
		cur_++;
	}
	void resize(size_t sz)
	{
		if (sz > 0 && list_ == NULL)
			lazyInit();
		if (sz <= cur_)
		{
			cur_ = sz;
			return;
		}
		if (sz_ < sz)
			expandCopy(sz);
		cur_ = sz;
	}
	void clear()
	{
		cur_ = 0;
	}
	inline ELList<T, S1, S2> &back()
	{
		assert_gt(cur_, 0);
		return list_[cur_ - 1];
	}
	inline const ELList<T, S1, S2> &back() const
	{
		assert_gt(cur_, 0);
		return list_[cur_ - 1];
	}
	inline ELList<T, S1, S2> &front()
	{
		assert_gt(cur_, 0);
		return list_[0];
	}
	inline const ELList<T, S1, S2> &front() const { return front(); }
	inline ELList<T, S1, S2> &operator[](size_t i)
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline const ELList<T, S1, S2> &operator[](size_t i) const
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline ELList<T, S1, S2> &get(size_t i)
	{
		return operator[](i);
	}
	inline const ELList<T, S1, S2> &get(size_t i) const
	{
		return operator[](i);
	}
	ELList<T, S1, S2> &getSlow(size_t i)
	{
		return operator[](i);
	}
	const ELList<T, S1, S2> &getSlow(size_t i) const
	{
		return operator[](i);
	}
	ELList<T, S1, S2> *ptr() { return list_; }
	void setCat(int cat)
	{
		assert_gt(cat, 0);
		cat_ = cat;
		if (cat_ != 0)
		{
			for (size_t i = 0; i < sz_; i++)
			{
				assert(list_[i].null());
				list_[i].setCat(cat_);
			}
		}
	}
	int cat() const { return cat_; }
protected:
	void lazyInit()
	{
		assert(null());
		list_ = alloc(sz_);
	}
	ELList<T, S1, S2> *alloc(size_t sz)
	{
		assert_gt(sz, 0);
		ELList<T, S1, S2> *tmp = new ELList<T, S1, S2>[sz];
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#endif
		if (cat_ != 0)
		{
			for (size_t i = 0; i < sz; i++)
			{
				assert(tmp[i].ptr() == NULL);
				tmp[i].setCat(cat_);
			}
		}
		return tmp;
	}
	void free()
	{
		if (list_ != NULL)
		{
			delete[] list_;
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, sz_);
#endif
			list_ = NULL;
		}
	}
	void expandCopy(size_t thresh)
	{
		assert(list_ != NULL);
		if (thresh <= sz_)
			return;
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
			newsz *= 2;
		ELList<T, S1, S2> *tmp = alloc(newsz);
		if (list_ != NULL)
		{
			for (size_t i = 0; i < cur_; i++)
			{
				assert_eq(cat_, tmp[i].cat());
				tmp[i].xfer(list_[i]);
				assert_eq(cat_, tmp[i].cat());
			}
			free();
		}
		list_ = tmp;
		sz_ = newsz;
	}
	void expandNoCopy(size_t thresh)
	{
		assert(list_ != NULL);
		if (thresh <= sz_)
			return;
		free();
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
			newsz *= 2;
		ELList<T, S1, S2> *tmp = alloc(newsz);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}
	int cat_;
	ELList<T, S1, S2> *list_;
	size_t sz_;
	size_t cur_;
};
template <typename T>
class ESet
{
public:
	ESet(int cat = 0) : cat_(cat),
						list_(NULL),
						sz_(0),
						cur_(0)
	{
#ifndef USE_MEM_TALLY
		(void)cat_;
#endif
		if (sz_ > 0)
		{
			list_ = alloc(sz_);
		}
	}
	ESet(size_t isz, int cat = 0) : cat_(cat),
									list_(NULL),
									sz_(isz),
									cur_(0)
	{
		assert_gt(isz, 0);
		if (sz_ > 0)
		{
			list_ = alloc(sz_);
		}
	}
	ESet(const ESet<T> &o, int cat = 0) : cat_(cat), list_(NULL)
	{
		assert_eq(cat_, o.cat());
		*this = o;
	}
	~ESet() { free(); }
	ESet &operator=(const ESet<T> &o)
	{
		assert_eq(cat_, o.cat());
		sz_ = o.sz_;
		cur_ = o.cur_;
		free();
		if (sz_ > 0)
		{
			list_ = alloc(sz_);
			memcpy(list_, o.list_, cur_ * sizeof(T));
		}
		else
		{
			list_ = NULL;
		}
		return *this;
	}
	size_t size() const { return cur_; }
	size_t totalSizeBytes() const
	{
		return sizeof(int) + cur_ * sizeof(T) + 2 * sizeof(size_t);
	}
	size_t totalCapacityBytes() const
	{
		return sizeof(int) + sz_ * sizeof(T) + 2 * sizeof(size_t);
	}
	bool empty() const { return cur_ == 0; }
	bool null() const { return list_ == NULL; }
	bool insert(const T &el)
	{
		size_t i = 0;
		if (cur_ == 0)
		{
			insert(el, 0);
			return true;
		}
		if (cur_ < 16)
		{
			i = scanLoBound(el);
		}
		else
		{
			i = bsearchLoBound(el);
		}
		if (i < cur_ && list_[i] == el)
			return false;
		insert(el, i);
		return true;
	}
	bool contains(const T &el) const
	{
		if (cur_ == 0)
		{
			return false;
		}
		else if (cur_ == 1)
		{
			return el == list_[0];
		}
		size_t i;
		if (cur_ < 16)
		{
			i = scanLoBound(el);
		}
		else
		{
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i] == el;
	}
	void remove(const T &el)
	{
		size_t i;
		if (cur_ < 16)
		{
			i = scanLoBound(el);
		}
		else
		{
			i = bsearchLoBound(el);
		}
		assert(i != cur_ && list_[i] == el);
		erase(i);
	}
	void resize(size_t sz)
	{
		if (sz <= cur_)
			return;
		if (sz_ < sz)
			expandCopy(sz);
	}
	void clear() { cur_ = 0; }
	int cat() const { return cat_; }
	void setCat(int cat)
	{
		cat_ = cat;
	}
	void xfer(ESet<T> &o)
	{
		assert_eq(cat_, o.cat());
		free();
		list_ = o.list_;
		sz_ = o.sz_;
		cur_ = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
	}
	T *ptr() { return list_; }
	const T *ptr() const { return list_; }
private:
	T *alloc(size_t sz)
	{
		assert_gt(sz, 0);
		T *tmp = new T[sz];
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#endif
		return tmp;
	}
	void free()
	{
		if (list_ != NULL)
		{
			delete[] list_;
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, sz_);
#endif
			list_ = NULL;
		}
	}
	size_t scanLoBound(const T &el) const
	{
		for (size_t i = 0; i < cur_; i++)
		{
			if (!(list_[i] < el))
			{
				return i;
			}
		}
		return cur_;
	}
	size_t bsearchLoBound(const T &el) const
	{
		size_t hi = cur_;
		size_t lo = 0;
		while (true)
		{
			if (lo == hi)
			{
#ifndef NDEBUG
				if ((rand() % 10) == 0)
				{
					assert_eq(lo, scanLoBound(el));
				}
#endif
				return lo;
			}
			size_t mid = lo + ((hi - lo) >> 1);
			assert_neq(mid, hi);
			if (list_[mid] < el)
			{
				if (lo == mid)
				{
#ifndef NDEBUG
					if ((rand() % 10) == 0)
					{
						assert_eq(hi, scanLoBound(el));
					}
#endif
					return hi;
				}
				lo = mid;
			}
			else
			{
				hi = mid;
			}
		}
	}
	bool sorted() const
	{
		if (cur_ <= 1)
			return true;
#ifndef NDEBUG
		if ((rand() % 20) == 0)
		{
			for (size_t i = 0; i < cur_ - 1; i++)
			{
				assert(list_[i] < list_[i + 1]);
			}
		}
#endif
		return true;
	}
	void insert(const T &el, size_t idx)
	{
		assert_leq(idx, cur_);
		if (cur_ == sz_)
		{
			expandCopy(sz_ + 1);
			assert(sorted());
		}
		for (size_t i = cur_; i > idx; i--)
		{
			list_[i] = list_[i - 1];
		}
		list_[idx] = el;
		cur_++;
		assert(sorted());
	}
	void erase(size_t idx)
	{
		assert_lt(idx, cur_);
		for (size_t i = idx; i < cur_ - 1; i++)
		{
			list_[i] = list_[i + 1];
		}
		cur_--;
		assert(sorted());
	}
	void expandCopy(size_t thresh)
	{
		if (thresh <= sz_)
			return;
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
		{
			newsz *= 2;
		}
		T *tmp = alloc(newsz);
		for (size_t i = 0; i < cur_; i++)
		{
			tmp[i] = list_[i];
		}
		free();
		list_ = tmp;
		sz_ = newsz;
	}
	int cat_;
	T *list_;
	size_t sz_;
	size_t cur_;
};
template <typename T, int S = 128>
class ELSet
{
public:
	explicit ELSet(int cat = 0) : cat_(cat), list_(NULL), sz_(S), cur_(0)
	{
#ifndef USE_MEM_TALLY
		(void)cat_;
#endif
		assert_geq(cat, 0);
	}
	explicit ELSet(size_t isz, int cat = 0) : cat_(cat), list_(NULL), sz_(isz), cur_(0)
	{
		assert_gt(isz, 0);
		assert_geq(cat, 0);
	}
	ELSet(const ELSet<T, S> &o) : cat_(0), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
	}
	explicit ELSet(const ELSet<T, S> &o, int cat) : cat_(cat), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
		assert_geq(cat, 0);
	}
	~ELSet() { free(); }
	ELSet<T, S> &operator=(const ELSet<T, S> &o)
	{
		assert_eq(cat_, o.cat());
		if (list_ == NULL)
		{
			lazyInit();
		}
		if (o.cur_ == 0)
		{
			cur_ = 0;
			return *this;
		}
		if (sz_ < o.cur_)
			expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for (size_t i = 0; i < cur_; i++)
		{
			assert_eq(list_[i].cat(), o.list_[i].cat());
			list_[i] = o.list_[i];
		}
		return *this;
	}
	void xfer(ELSet<T, S> &o)
	{
		assert_eq(cat_, o.cat());
		list_ = o.list_;
		sz_ = o.sz_;
		cur_ = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
	}
	inline size_t size() const { return cur_; }
	inline bool empty() const { return cur_ == 0; }
	inline bool null() const { return list_ == NULL; }
	void expand()
	{
		if (list_ == NULL)
			lazyInit();
		if (cur_ == sz_)
			expandCopy(sz_ + 1);
		cur_++;
	}
	void resize(size_t sz)
	{
		if (sz > 0 && list_ == NULL)
			lazyInit();
		if (sz <= cur_)
		{
			cur_ = sz;
			return;
		}
		if (sz_ < sz)
		{
			expandCopy(sz);
		}
		cur_ = sz;
	}
	void clear()
	{
		cur_ = 0;
	}
	inline ESet<T> &back()
	{
		assert_gt(cur_, 0);
		return list_[cur_ - 1];
	}
	inline const ESet<T> &back() const
	{
		assert_gt(cur_, 0);
		return list_[cur_ - 1];
	}
	inline ESet<T> &front()
	{
		assert_gt(cur_, 0);
		return list_[0];
	}
	inline const ESet<T> &front() const { return front(); }
	inline ESet<T> &operator[](size_t i)
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline const ESet<T> &operator[](size_t i) const
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline ESet<T> &get(size_t i)
	{
		return operator[](i);
	}
	inline const ESet<T> &get(size_t i) const
	{
		return operator[](i);
	}
	ESet<T> &getSlow(size_t i)
	{
		return operator[](i);
	}
	const ESet<T> &getSlow(size_t i) const
	{
		return operator[](i);
	}
	ESet<T> *ptr() { return list_; }
	const ESet<T> *ptr() const { return list_; }
	void setCat(int cat)
	{
		assert_gt(cat, 0);
		cat_ = cat;
		if (cat_ != 0)
		{
			for (size_t i = 0; i < sz_; i++)
			{
				assert(list_[i].null());
				list_[i].setCat(cat_);
			}
		}
	}
	int cat() const { return cat_; }
protected:
	void lazyInit()
	{
		assert(list_ == NULL);
		list_ = alloc(sz_);
	}
	ESet<T> *alloc(size_t sz)
	{
		assert_gt(sz, 0);
		ESet<T> *tmp = new ESet<T>[sz];
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#endif
		if (cat_ != 0)
		{
			for (size_t i = 0; i < sz; i++)
			{
				assert(tmp[i].ptr() == NULL);
				tmp[i].setCat(cat_);
			}
		}
		return tmp;
	}
	void free()
	{
		if (list_ != NULL)
		{
			delete[] list_;
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, sz_);
#endif
			list_ = NULL;
		}
	}
	void expandCopy(size_t thresh)
	{
		assert(list_ != NULL);
		if (thresh <= sz_)
			return;
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
			newsz *= 2;
		ESet<T> *tmp = alloc(newsz);
		if (list_ != NULL)
		{
			for (size_t i = 0; i < cur_; i++)
			{
				assert_eq(cat_, tmp[i].cat());
				tmp[i].xfer(list_[i]);
				assert_eq(cat_, tmp[i].cat());
			}
			free();
		}
		list_ = tmp;
		sz_ = newsz;
	}
	void expandNoCopy(size_t thresh)
	{
		assert(list_ != NULL);
		if (thresh <= sz_)
			return;
		free();
		size_t newsz = (sz_ * 2) + 1;
		while (newsz < thresh)
			newsz *= 2;
		ESet<T> *tmp = alloc(newsz);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}
	int cat_;
	ESet<T> *list_;
	size_t sz_;
	size_t cur_;
};
template <typename K, typename V>
class EMap
{
public:
	EMap(int cat = 0) : cat_(cat),
						list_(NULL),
						sz_(128),
						cur_(0)
	{
#ifndef USE_MEM_TALLY
		(void)cat_;
#endif
		list_ = alloc(sz_);
	}
	EMap(size_t isz, int cat = 0) : cat_(cat),
									list_(NULL),
									sz_(isz),
									cur_(0)
	{
		assert_gt(isz, 0);
		list_ = alloc(sz_);
	}
	EMap(const EMap<K, V> &o) : list_(NULL)
	{
		*this = o;
	}
	~EMap() { free(); }
	EMap &operator=(const EMap<K, V> &o)
	{
		sz_ = o.sz_;
		cur_ = o.cur_;
		free();
		list_ = alloc(sz_);
		for (size_t i = 0; i < cur_; i++)
			list_[i] = o.list_[i];
		return *this;
	}
	size_t size() const { return cur_; }
	size_t totalSizeBytes() const
	{
		return sizeof(int) +
			   2 * sizeof(size_t) +
			   cur_ * sizeof(std::pair<K, V>);
	}
	size_t totalCapacityBytes() const
	{
		return sizeof(int) +
			   2 * sizeof(size_t) +
			   sz_ * sizeof(std::pair<K, V>);
	}
	bool empty() const { return cur_ == 0; }
	bool insert(const std::pair<K, V> &el)
	{
		size_t i = 0;
		if (cur_ == 0)
		{
			insert(el, 0);
			return true;
		}
		if (cur_ < 16)
		{
			i = scanLoBound(el.first);
		}
		else
		{
			i = bsearchLoBound(el.first);
		}
		if (list_[i] == el)
			return false;
		insert(el, i);
		return true;
	}
	bool contains(const K &el) const
	{
		if (cur_ == 0)
			return false;
		else if (cur_ == 1)
			return el == list_[0].first;
		size_t i;
		if (cur_ < 16)
		{
			i = scanLoBound(el);
		}
		else
		{
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i].first == el;
	}
	bool containsEx(const K &el, size_t &i) const
	{
		if (cur_ == 0)
			return false;
		else if (cur_ == 1)
		{
			i = 0;
			return el == list_[0].first;
		}
		if (cur_ < 16)
		{
			i = scanLoBound(el);
		}
		else
		{
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i].first == el;
	}
	void remove(const K &el)
	{
		size_t i;
		if (cur_ < 16)
		{
			i = scanLoBound(el);
		}
		else
		{
			i = bsearchLoBound(el);
		}
		assert(i != cur_ && list_[i].first == el);
		erase(i);
	}
	void resize(size_t sz)
	{
		if (sz <= cur_)
			return;
		if (sz_ < sz)
			expandCopy(sz);
	}
	const std::pair<K, V> &get(size_t i) const
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	const std::pair<K, V> &operator[](size_t i) const
	{
		return get(i);
	}
	void clear() { cur_ = 0; }
private:
	std::pair<K, V> *alloc(size_t sz)
	{
		assert_gt(sz, 0);
		std::pair<K, V> *tmp = new std::pair<K, V>[sz];
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#endif
		return tmp;
	}
	void free()
	{
		if (list_ != NULL)
		{
			delete[] list_;
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, sz_);
#endif
			list_ = NULL;
		}
	}
	size_t scanLoBound(const K &el) const
	{
		for (size_t i = 0; i < cur_; i++)
		{
			if (!(list_[i].first < el))
			{
				return i;
			}
		}
		return cur_;
	}
	size_t bsearchLoBound(const K &el) const
	{
		size_t hi = cur_;
		size_t lo = 0;
		while (true)
		{
			if (lo == hi)
			{
#ifndef NDEBUG
				if ((rand() % 10) == 0)
				{
					assert_eq(lo, scanLoBound(el));
				}
#endif
				return lo;
			}
			size_t mid = lo + ((hi - lo) >> 1);
			assert_neq(mid, hi);
			if (list_[mid].first < el)
			{
				if (lo == mid)
				{
#ifndef NDEBUG
					if ((rand() % 10) == 0)
					{
						assert_eq(hi, scanLoBound(el));
					}
#endif
					return hi;
				}
				lo = mid;
			}
			else
			{
				hi = mid;
			}
		}
	}
	bool sorted() const
	{
		if (cur_ <= 1)
			return true;
#ifndef NDEBUG
		for (size_t i = 0; i < cur_ - 1; i++)
		{
			assert(!(list_[i] == list_[i + 1]));
			assert(list_[i] < list_[i + 1]);
		}
#endif
		return true;
	}
	void insert(const std::pair<K, V> &el, size_t idx)
	{
		assert_leq(idx, cur_);
		if (cur_ == sz_)
		{
			expandCopy(sz_ + 1);
		}
		for (size_t i = cur_; i > idx; i--)
		{
			list_[i] = list_[i - 1];
		}
		list_[idx] = el;
		assert(idx == cur_ || list_[idx] < list_[idx + 1]);
		cur_++;
		assert(sorted());
	}
	void erase(size_t idx)
	{
		assert_lt(idx, cur_);
		for (size_t i = idx; i < cur_ - 1; i++)
		{
			list_[i] = list_[i + 1];
		}
		cur_--;
		assert(sorted());
	}
	void expandCopy(size_t thresh)
	{
		if (thresh <= sz_)
			return;
		size_t newsz = sz_ * 2;
		while (newsz < thresh)
			newsz *= 2;
		std::pair<K, V> *tmp = alloc(newsz);
		for (size_t i = 0; i < cur_; i++)
		{
			tmp[i] = list_[i];
		}
		free();
		list_ = tmp;
		sz_ = newsz;
	}
	int cat_;
	std::pair<K, V> *list_;
	size_t sz_;
	size_t cur_;
};
template <typename T, int S = 128>
class EFactory
{
public:
	explicit EFactory(size_t isz, int cat = 0) : l_(isz, cat) {}
	explicit EFactory(int cat = 0) : l_(cat) {}
	void clear()
	{
		l_.clear();
	}
	size_t alloc()
	{
		l_.expand();
		return l_.size() - 1;
	}
	size_t size() const
	{
		return l_.size();
	}
	size_t totalSizeBytes() const
	{
		return l_.totalSizeBytes();
	}
	size_t totalCapacityBytes() const
	{
		return l_.totalCapacityBytes();
	}
	void resize(size_t sz)
	{
		l_.resize(sz);
	}
	bool empty() const
	{
		return size() == 0;
	}
	void pop()
	{
		l_.resize(l_.size() - 1);
	}
	T &operator[](size_t off)
	{
		return l_[off];
	}
	const T &operator[](size_t off) const
	{
		return l_[off];
	}
protected:
	EList<T, S> l_;
};
template <int S = 128>
class EBitList
{
public:
	explicit EBitList(size_t isz, int cat = 0) : l_(isz, cat) { reset(); }
	explicit EBitList(int cat = 0) : l_(cat) { reset(); }
	void clear()
	{
		reset();
	}
	void reset()
	{
		l_.clear();
		max_ = std::numeric_limits<size_t>::max();
	}
	void set(size_t off)
	{
		resize(off);
		l_[off >> 3] |= (1 << (off & 7));
		if (off > max_ || max_ == std::numeric_limits<size_t>::max())
		{
			max_ = off;
		}
	}
	bool test(size_t off) const
	{
		if ((size_t)(off >> 3) >= l_.size())
		{
			return false;
		}
		return (l_[off >> 3] & (1 << (off & 7))) != 0;
	}
	size_t size() const
	{
		return l_.size();
	}
	void resize(size_t off)
	{
		if ((size_t)(off >> 3) >= l_.size())
		{
			size_t oldsz = l_.size();
			l_.resize((off >> 3) + 1);
			for (size_t i = oldsz; i < l_.size(); i++)
			{
				l_[i] = 0;
			}
		}
	}
	size_t max() const
	{
		return max_;
	}
protected:
	EList<uint8_t, S> l_;
	size_t max_;
};
template <typename T, int S = 128>
class EHeap
{
public:
	void insert(T o)
	{
		size_t pos = l_.size();
		l_.push_back(o);
		while (pos > 0)
		{
			size_t parent = (pos - 1) >> 1;
			if (l_[pos] < l_[parent])
			{
				T tmp(l_[pos]);
				l_[pos] = l_[parent];
				l_[parent] = tmp;
				pos = parent;
			}
			else
				break;
		}
		assert(repOk());
	}
	T top()
	{
		assert_gt(l_.size(), 0);
		return l_[0];
	}
	T pop()
	{
		assert_gt(l_.size(), 0);
		T ret = l_[0];
		l_[0] = l_[l_.size() - 1];
		l_.resize(l_.size() - 1);
		size_t cur = 0;
		while (true)
		{
			size_t c1 = ((cur + 1) << 1) - 1;
			size_t c2 = c1 + 1;
			if (c2 < l_.size())
			{
				if (l_[c1] < l_[cur] && l_[c1] <= l_[c2])
				{
					T tmp(l_[c1]);
					l_[c1] = l_[cur];
					l_[cur] = tmp;
					cur = c1;
				}
				else if (l_[c2] < l_[cur])
				{
					T tmp(l_[c2]);
					l_[c2] = l_[cur];
					l_[cur] = tmp;
					cur = c2;
				}
				else
				{
					break;
				}
			}
			else if (c1 < l_.size())
			{
				if (l_[c1] < l_[cur])
				{
					T tmp(l_[c1]);
					l_[c1] = l_[cur];
					l_[cur] = tmp;
					cur = c1;
				}
				else
				{
					break;
				}
			}
			else
			{
				break;
			}
		}
		assert(repOk());
		return ret;
	}
	size_t size() const
	{
		return l_.size();
	}
	size_t totalSizeBytes() const
	{
		return l_.totalSizeBytes();
	}
	size_t totalCapacityBytes() const
	{
		return l_.totalCapacityBytes();
	}
	bool empty() const
	{
		return l_.empty();
	}
	const T &operator[](size_t i) const
	{
		return l_[i];
	}
#ifndef NDEBUG
	bool repOk() const
	{
		if (empty())
			return true;
		return repOkNode(0);
	}
	bool repOkNode(size_t cur) const
	{
		size_t c1 = ((cur + 1) << 1) - 1;
		size_t c2 = c1 + 1;
		if (c1 < l_.size())
		{
			assert(l_[cur] <= l_[c1]);
		}
		if (c2 < l_.size())
		{
			assert(l_[cur] <= l_[c2]);
		}
		if (c2 < l_.size())
		{
			return repOkNode(c1) && repOkNode(c2);
		}
		else if (c1 < l_.size())
		{
			return repOkNode(c1);
		}
		return true;
	}
#endif
	void clear()
	{
		l_.clear();
	}
protected:
	EList<T, S> l_;
};
class Pool
{
public:
	Pool(
		uint64_t bytes,
		uint32_t pagesz,
		int cat = 0) : cat_(cat),
					   cur_(0),
					   pagesz_(pagesz),
					   pages_(cat)
	{
		size_t super_page_num = ((bytes + pagesz - 1) / pagesz + 1);
		super_pages = new uint8_t[pagesz * super_page_num];
		for (size_t i = 0; i < super_page_num; i++)
		{
			pages_.push_back(&super_pages[i * pagesz]);
#ifdef USE_MEM_TALLY
			gMemTally.add(cat, pagesz);
#else
			(void)cat_;
			(void)pagesz_;
#endif
			assert(pages_.back() != NULL);
		}
		assert(repOk());
	}
	~Pool()
	{
		delete[] super_pages;
		for (size_t i = 0; i < pages_.size(); i++)
		{
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, pagesz_);
#else
			(void)cat_;
#endif
		}
	}
	uint8_t *alloc()
	{
		assert(repOk());
		if (cur_ == pages_.size())
			return NULL;
		return pages_[cur_++];
	}
	void clear()
	{
		cur_ = 0;
		assert(repOk());
	}
	void free()
	{
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_leq(cur_, pages_.size());
		assert(!pages_.empty());
		assert_gt(pagesz_, 0);
		return true;
	}
#endif
private:
	uint8_t *super_pages;
	int cat_;
	size_t cur_;
	const size_t pagesz_;
	EList<uint8_t *> pages_;
};
template <typename T, int S>
class PList
{
#define PLIST_PER_PAGE (S / sizeof(T))
public:
	PList(int cat = 0) : cur_(0),
						 curPage_(0),
						 pages_(cat) {}
	bool add(Pool &p, const T &o)
	{
		assert(repOk());
		if (!ensure(p, 1))
			return false;
		if (cur_ == PLIST_PER_PAGE)
		{
			cur_ = 0;
			curPage_++;
		}
		assert_lt(curPage_, pages_.size());
		assert(repOk());
		assert_lt(cur_, PLIST_PER_PAGE);
		pages_[curPage_][cur_++] = o;
		return true;
	}
	bool add(Pool &p, const EList<T> &os)
	{
		if (!ensure(p, os.size()))
			return false;
		for (size_t i = 0; i < os.size(); i++)
		{
			if (cur_ == PLIST_PER_PAGE)
			{
				cur_ = 0;
				curPage_++;
			}
			assert_lt(curPage_, pages_.size());
			assert(repOk());
			assert_lt(cur_, PLIST_PER_PAGE);
			pages_[curPage_][cur_++] = os[i];
		}
		return true;
	}
	bool copy(
		Pool &p,
		const PList<T, S> &src,
		size_t i,
		size_t len)
	{
		if (!ensure(p, src.size()))
			return false;
		for (size_t i = 0; i < src.size(); i++)
		{
			if (cur_ == PLIST_PER_PAGE)
			{
				cur_ = 0;
				curPage_++;
			}
			assert_lt(curPage_, pages_.size());
			assert(repOk());
			assert_lt(cur_, PLIST_PER_PAGE);
			pages_[curPage_][cur_++] = src[i];
		}
		return true;
	}
	bool addFill(Pool &p, size_t num, const T &o)
	{
		if (!ensure(p, num))
			return false;
		for (size_t i = 0; i < num; i++)
		{
			if (cur_ == PLIST_PER_PAGE)
			{
				cur_ = 0;
				curPage_++;
			}
			assert_lt(curPage_, pages_.size());
			assert(repOk());
			assert_lt(cur_, PLIST_PER_PAGE);
			pages_[curPage_][cur_++] = o;
		}
		return true;
	}
	void clear()
	{
		pages_.clear();
		cur_ = curPage_ = 0;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert(pages_.size() == 0 || curPage_ < pages_.size());
		assert_leq(cur_, PLIST_PER_PAGE);
		return true;
	}
#endif
	size_t size() const
	{
		return curPage_ * PLIST_PER_PAGE + cur_;
	}
	bool empty() const
	{
		return size() == 0;
	}
	inline const T &getConst(size_t i) const
	{
		assert_lt(i, size());
		size_t page = i / PLIST_PER_PAGE;
		size_t elt = i % PLIST_PER_PAGE;
		return pages_[page][elt];
	}
	inline T &get(size_t i)
	{
		assert_lt(i, size());
		size_t page = i / PLIST_PER_PAGE;
		size_t elt = i % PLIST_PER_PAGE;
		assert_lt(page, pages_.size());
		assert(page < pages_.size() - 1 || elt < cur_);
		return pages_[page][elt];
	}
	inline T &back()
	{
		size_t page = (size() - 1) / PLIST_PER_PAGE;
		size_t elt = (size() - 1) % PLIST_PER_PAGE;
		assert_lt(page, pages_.size());
		assert(page < pages_.size() - 1 || elt < cur_);
		return pages_[page][elt];
	}
	inline const T &back() const
	{
		size_t page = (size() - 1) / PLIST_PER_PAGE;
		size_t elt = (size() - 1) % PLIST_PER_PAGE;
		assert_lt(page, pages_.size());
		assert(page < pages_.size() - 1 || elt < cur_);
		return pages_[page][elt];
	}
	T &last()
	{
		assert(!pages_.empty());
		assert_gt(PLIST_PER_PAGE, 0);
		if (cur_ == 0)
		{
			assert_gt(pages_.size(), 1);
			return pages_[pages_.size() - 2][PLIST_PER_PAGE - 1];
		}
		else
		{
			return pages_.back()[cur_ - 1];
		}
	}
	bool ensure(Pool &p, size_t num)
	{
		assert(repOk());
		if (num == 0)
			return true;
		if (pages_.size() == 0)
		{
			if (expand(p) == NULL)
			{
				return false;
			}
			assert_eq(1, pages_.size());
		}
		size_t cur = cur_;
		size_t curPage = curPage_;
		while (cur + num > PLIST_PER_PAGE)
		{
			assert_lt(curPage, pages_.size());
			if (curPage == pages_.size() - 1 && expand(p) == NULL)
			{
				return false;
			}
			num -= (PLIST_PER_PAGE - cur);
			cur = 0;
			curPage++;
		}
		return true;
	}
protected:
	T *expand(Pool &p)
	{
		T *newpage = (T *)p.alloc();
		if (newpage == NULL)
		{
			return NULL;
		}
		pages_.push_back(newpage);
		return pages_.back();
	}
	size_t cur_;
	size_t curPage_;
	EList<T *> pages_;
};
template <typename T, int S>
class EListSlice
{
public:
	EListSlice() : i_(0),
				   len_(0),
				   list_()
	{
	}
	EListSlice(
		EList<T, S> &list,
		size_t i,
		size_t len) : i_(i),
					  len_(len),
					  list_(&list)
	{
	}
	void init(const EListSlice<T, S> &sl, size_t first, size_t last)
	{
		assert_gt(last, first);
		assert_leq(last - first, sl.len_);
		i_ = sl.i_ + first;
		len_ = last - first;
		list_ = sl.list_;
	}
	void reset()
	{
		i_ = len_ = 0;
		list_ = NULL;
	}
	inline const T &get(size_t i) const
	{
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i + i_);
	}
	inline T &get(size_t i)
	{
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i + i_);
	}
	inline T &operator[](size_t i)
	{
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i + i_);
	}
	inline const T &operator[](size_t i) const
	{
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i + i_);
	}
	bool valid() const
	{
		return len_ != 0;
	}
	size_t size() const
	{
		return len_;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_leq(i_ + len_, list_->size());
		return true;
	}
#endif
	bool operator==(const EListSlice &sl) const
	{
		return i_ == sl.i_ && len_ == sl.len_ && list_ == sl.list_;
	}
	bool operator!=(const EListSlice &sl) const
	{
		return !(*this == sl);
	}
	void setLength(size_t nlen)
	{
		len_ = nlen;
	}
protected:
	size_t i_;
	size_t len_;
	EList<T, S> *list_;
};
template <typename T, int S>
class PListSlice
{
public:
	PListSlice() : i_(0),
				   len_(0),
				   list_()
	{
	}
	PListSlice(
		PList<T, S> &list,
		TIndexOffU i,
		TIndexOffU len) : i_(i),
						  len_(len),
						  list_(&list)
	{
	}
	void init(const PListSlice<T, S> &sl, size_t first, size_t last)
	{
		assert_gt(last, first);
		assert_leq(last - first, sl.len_);
		i_ = (TIndexOffU)(sl.i_ + first);
		len_ = (TIndexOffU)(last - first);
		list_ = sl.list_;
	}
	void reset()
	{
		i_ = len_ = 0;
		list_ = NULL;
	}
	inline const T &get(size_t i) const
	{
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i + i_);
	}
	inline T &get(size_t i)
	{
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i + i_);
	}
	inline T &operator[](size_t i)
	{
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i + i_);
	}
	inline const T &operator[](size_t i) const
	{
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i + i_);
	}
	bool valid() const
	{
		return len_ != 0;
	}
	size_t size() const
	{
		return len_;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_leq(i_ + len_, list_->size());
		return true;
	}
#endif
	bool operator==(const PListSlice &sl) const
	{
		return i_ == sl.i_ && len_ == sl.len_ && list_ == sl.list_;
	}
	bool operator!=(const PListSlice &sl) const
	{
		return !(*this == sl);
	}
	void setLength(size_t nlen)
	{
		len_ = (TIndexOffU)nlen;
	}
protected:
	TIndexOffU i_;
	TIndexOffU len_;
	PList<T, S> *list_;
};
template <typename K, typename P>
class RedBlackNode
{
	typedef RedBlackNode<K, P> TNode;
public:
	TNode *parent;
	TNode *left;
	TNode *right;
	bool red;
	K key;
	P payload; 
	RedBlackNode *grandparent()
	{
		return parent != NULL ? parent->parent : NULL;
	}
	RedBlackNode *uncle()
	{
		if (parent == NULL)
			return NULL;
		if (parent->parent == NULL)
			return NULL;
		return (parent->parent->left == parent) ? parent->parent->right : parent->parent->left;
	}
	bool isLeftChild() const
	{
		assert(parent != NULL);
		return parent->left == this;
	}
	bool isRightChild() const
	{
		assert(parent != NULL);
		return parent->right == this;
	}
	void replaceChild(RedBlackNode *ol, RedBlackNode *nw)
	{
		if (left == ol)
		{
			left = nw;
		}
		else
		{
			assert(right == ol);
			right = nw;
		}
	}
	int numChildren() const
	{
		return ((left != NULL) ? 1 : 0) + ((right != NULL) ? 1 : 0);
	}
#ifndef NDEBUG
	bool repOk() const
	{
		if (parent != NULL)
		{
			assert(parent->left == this || parent->right == this);
		}
		return true;
	}
#endif
	bool operator<(const TNode &o) const
	{
		return key < o.key;
	}
	bool operator>(const TNode &o) const { return key > o.key; }
	bool operator==(const TNode &o) const { return key == o.key; }
	bool operator<(const K &okey) const { return key < okey; }
	bool operator>(const K &okey) const { return key > okey; }
	bool operator==(const K &okey) const { return key == okey; }
};
template <typename K, typename P>
class RedBlack
{
	typedef RedBlackNode<K, P> TNode;
public:
	RedBlack(size_t pageSz, int cat = 0) : perPage_(pageSz / sizeof(TNode)), pages_(cat) { clear(); }
	inline TNode *lookup(const K &key) const
	{
		TNode *cur = root_;
		while (cur != NULL)
		{
			if ((*cur) == key)
				return cur;
			if ((*cur) < key)
			{
				cur = cur->right;
			}
			else
			{
				cur = cur->left;
			}
		}
		return NULL;
	}
	TNode *add(
		Pool &p,
		const K &key, bool *added)
	{
		TNode *cur = root_;
		assert(root_ == NULL || !root_->red);
		TNode *parent = NULL;
		bool leftChild = true;
		while (cur != NULL)
		{
			if ((*cur) == key)
			{
				break;
			}
			parent = cur;
			if ((*cur) < key)
			{
				if ((cur = cur->right) == NULL)
				{
					leftChild = false;
				}
			}
			else
			{
				if ((cur = cur->left) == NULL)
				{
					leftChild = true;
				}
			}
		}
		if (cur != NULL)
		{
			if (added != NULL)
				*added = false;
		}
		else
		{
			assert(root_ == NULL || !root_->red);
			if (!addNode(p, cur))
			{
				return NULL;
			}
			assert(cur != NULL);
			assert(cur != root_);
			assert(cur != parent);
			cur->key = key;
			cur->left = cur->right = NULL;
			cur->red = true;
			keys_++;
			if (added != NULL)
				*added = true;
			addNode(cur, parent, leftChild);
		}
		return cur;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert(curPage_ == 0 || curPage_ < pages_.size());
		assert_leq(cur_, perPage_);
		assert(root_ == NULL || !root_->red);
		return true;
	}
#endif
	void clear()
	{
		cur_ = curPage_ = 0;
		root_ = NULL;
		keys_ = 0;
		intenseRepOkCnt_ = 0;
		pages_.clear();
	}
	size_t size() const
	{
		return keys_;
	}
	bool empty() const
	{
		return keys_ == 0;
	}
	bool addNode(Pool &p, TNode *&node)
	{
		assert_leq(cur_, perPage_);
		assert(repOk());
		assert(this != NULL);
		if (pages_.size() == 0)
		{
			if (addPage(p) == NULL)
			{
				node = NULL;
				return false;
			}
			assert_eq(1, pages_.size());
		}
		if (cur_ == perPage_)
		{
			assert_lt(curPage_, pages_.size());
			if (curPage_ == pages_.size() - 1 && addPage(p) == NULL)
			{
				return false;
			}
			cur_ = 0;
			curPage_++;
		}
		assert_lt(cur_, perPage_);
		assert_lt(curPage_, pages_.size());
		node = &pages_[curPage_][cur_];
		assert(node != NULL);
		cur_++;
		return true;
	}
protected:
#ifndef NDEBUG
	bool redBlackRepOk(TNode *n)
	{
		if (n == NULL)
			return true;
		if (++intenseRepOkCnt_ < 500)
			return true;
		intenseRepOkCnt_ = 0;
		int minNodes = -1;
		int maxNodes = -1;
		int blackConst = -1;
		size_t nodesTot = 0;
		redBlackRepOk(
			n,
			1,
			n->red ? 0 : 1, 
			blackConst,
			minNodes,
			maxNodes,
			nodesTot);
		if (n == root_)
		{
			assert_eq(nodesTot, keys_);
		}
		assert_gt(minNodes, 0);
		assert_gt(maxNodes, 0);
		assert_leq(maxNodes, 2 * minNodes);
		return true;
	}
	bool redBlackRepOk(
		TNode *n,
		int nodes,
		int black,
		int &blackConst,
		int &minNodes,
		int &maxNodes,
		size_t &nodesTot) const
	{
		assert_gt(black, 0);
		nodesTot++;
		if (n->left == NULL)
		{
			if (blackConst == -1)
				blackConst = black;
			assert_eq(black, blackConst);
			if (nodes + 1 > maxNodes)
				maxNodes = nodes + 1;
			if (nodes + 1 < minNodes || minNodes == -1)
				minNodes = nodes + 1;
		}
		else
		{
			if (n->red)
				assert(!n->left->red);
			redBlackRepOk(
				n->left,
				nodes + 1, black + (n->left->red ? 0 : 1), blackConst, minNodes, maxNodes, nodesTot);
		}
		if (n->right == NULL)
		{
			if (blackConst == -1)
				blackConst = black;
			assert_eq(black, blackConst);
			if (nodes + 1 > maxNodes)
				maxNodes = nodes + 1;
			if (nodes + 1 < minNodes || minNodes == -1)
				minNodes = nodes + 1;
		}
		else
		{
			if (n->red)
				assert(!n->right->red);
			redBlackRepOk(
				n->right,
				nodes + 1, black + (n->right->red ? 0 : 1), blackConst, minNodes, maxNodes, nodesTot);
		}
		return true;
	}
#endif
	void leftRotate(TNode *n)
	{
		TNode *r = n->right;
		assert(n->repOk());
		assert(r->repOk());
		n->right = r->left;
		if (n->right != NULL)
		{
			n->right->parent = n;
			assert(n->right->repOk());
		}
		r->parent = n->parent;
		n->parent = r;
		r->left = n;
		if (r->parent != NULL)
		{
			r->parent->replaceChild(n, r);
		}
		if (root_ == n)
			root_ = r;
		assert(!root_->red);
		assert(n->repOk());
		assert(r->repOk());
	}
	void rightRotate(TNode *n)
	{
		TNode *r = n->left;
		assert(n->repOk());
		assert(r->repOk());
		n->left = r->right;
		if (n->left != NULL)
		{
			n->left->parent = n;
			assert(n->left->repOk());
		}
		r->parent = n->parent;
		n->parent = r;
		r->right = n;
		if (r->parent != NULL)
		{
			r->parent->replaceChild(n, r);
		}
		if (root_ == n)
			root_ = r;
		assert(!root_->red);
		assert(n->repOk());
		assert(r->repOk());
	}
	void addNode(TNode *n, TNode *parent, bool leftChild)
	{
		assert(n != NULL);
		if (parent == NULL)
		{
			root_ = n;
			root_->red = false;
			n->parent = NULL;
			assert(redBlackRepOk(root_));
			assert(n->repOk());
		}
		else
		{
			assert(!root_->red);
			if (leftChild)
			{
				assert(parent->left == NULL);
				parent->left = n;
			}
			else
			{
				assert(parent->right == NULL);
				parent->right = n;
			}
			n->parent = parent;
			int thru = 0;
			while (true)
			{
				thru++;
				parent = n->parent;
				if (parent != NULL)
					assert(parent->repOk());
				if (parent == NULL && n->red)
				{
					n->red = false;
				}
				if (parent == NULL || !parent->red)
				{
					assert(redBlackRepOk(root_));
					break;
				}
				TNode *uncle = n->uncle();
				TNode *gparent = n->grandparent();
				assert(gparent != NULL);
				bool uncleRed = (uncle != NULL ? uncle->red : false);
				if (uncleRed)
				{
					assert(uncle != NULL);
					parent->red = uncle->red = false;
					gparent->red = true;
					n = gparent;
					continue;
				}
				else
				{
					if (parent->isLeftChild())
					{
						if (!n->isLeftChild())
						{
							n = parent;
							leftRotate(n);
						}
						n = n->parent;
						n->red = false;
						n->parent->red = true;
						rightRotate(n->parent);
						assert(redBlackRepOk(n));
						assert(redBlackRepOk(root_));
					}
					else
					{
						if (!n->isRightChild())
						{
							n = parent;
							rightRotate(n);
						}
						n = n->parent;
						n->red = false;
						n->parent->red = true;
						leftRotate(n->parent);
						assert(redBlackRepOk(n));
						assert(redBlackRepOk(root_));
					}
				}
				break;
			}
		}
		assert(redBlackRepOk(root_));
	}
	TNode *addPage(Pool &p)
	{
		TNode *n = (TNode *)p.alloc();
		if (n != NULL)
		{
			pages_.push_back(n);
		}
		return n;
	}
	size_t keys_;
	size_t cur_;
	size_t curPage_;
	const size_t perPage_;
	TNode *root_;
	EList<TNode *> pages_;
	int intenseRepOkCnt_;
};
template <typename T>
struct DoublyLinkedList
{
	DoublyLinkedList() : payload(), prev(NULL), next(NULL) {}
	void toList(EList<T> &l)
	{
		DoublyLinkedList<T> *cur = this;
		while (cur != NULL)
		{
			l.push_back(cur->payload);
			cur = cur->next;
		}
		cur = prev;
		while (cur != NULL)
		{
			l.push_back(cur->payload);
			cur = cur->prev;
		}
	}
	T payload;
	DoublyLinkedList<T> *prev;
	DoublyLinkedList<T> *next;
};
template <typename T1, typename T2>
struct Pair
{
	T1 a;
	T2 b;
	Pair() : a(), b() {}
	Pair(
		const T1 &a_,
		const T2 &b_)
	{
		a = a_;
		b = b_;
	}
	bool operator==(const Pair &o) const
	{
		return a == o.a && b == o.b;
	}
	bool operator<(const Pair &o) const
	{
		if (a < o.a)
			return true;
		if (a > o.a)
			return false;
		if (b < o.b)
			return true;
		return false;
	}
};
template <typename T1, typename T2, typename T3>
struct Triple
{
	T1 a;
	T2 b;
	T3 c;
	Triple() : a(), b(), c() {}
	Triple(
		const T1 &a_,
		const T2 &b_,
		const T3 &c_)
	{
		a = a_;
		b = b_;
		c = c_;
	}
	bool operator==(const Triple &o) const
	{
		return a == o.a && b == o.b && c == o.c;
	}
	bool operator<(const Triple &o) const
	{
		if (a < o.a)
			return true;
		if (a > o.a)
			return false;
		if (b < o.b)
			return true;
		if (b > o.b)
			return false;
		if (c < o.c)
			return true;
		return false;
	}
};
template <typename T1, typename T2, typename T3, typename T4>
struct Quad
{
	Quad() : a(), b(), c(), d() {}
	Quad(
		const T1 &a_,
		const T2 &b_,
		const T3 &c_,
		const T4 &d_)
	{
		a = a_;
		b = b_;
		c = c_;
		d = d_;
	}
	Quad(
		const T1 &a_,
		const T1 &b_,
		const T1 &c_,
		const T1 &d_)
	{
		init(a_, b_, c_, d_);
	}
	void init(
		const T1 &a_,
		const T1 &b_,
		const T1 &c_,
		const T1 &d_)
	{
		a = a_;
		b = b_;
		c = c_;
		d = d_;
	}
	bool operator==(const Quad &o) const
	{
		return a == o.a && b == o.b && c == o.c && d == o.d;
	}
	bool operator<(const Quad &o) const
	{
		if (a < o.a)
			return true;
		if (a > o.a)
			return false;
		if (b < o.b)
			return true;
		if (b > o.b)
			return false;
		if (c < o.c)
			return true;
		if (c > o.c)
			return false;
		if (d < o.d)
			return true;
		return false;
	}
	T1 a;
	T2 b;
	T3 c;
	T4 d;
};
#endif
