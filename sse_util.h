#ifndef SSE_UTIL_H_
#define SSE_UTIL_H_
#include "assert_helpers.h"
#include "ds.h"
#include "limit.h"
#include <iostream>
#include "sse_wrap.h"
class EList_m128i
{
public:
	explicit EList_m128i(int cat = 0) : cat_(cat), last_alloc_(NULL), list_(NULL), sz_(0), cur_(0)
	{
		assert_geq(cat, 0);
	}
	~EList_m128i() { free(); }
	inline size_t size() const { return cur_; }
	inline size_t capacity() const { return sz_; }
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
	void zero()
	{
		if (cur_ > 0)
		{
			memset(list_, 0, cur_ * sizeof(__m128i));
		}
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
		{
			expandNoCopy(sz);
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
	void clear()
	{
		cur_ = 0;
	}
	inline __m128i &operator[](size_t i)
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline __m128i operator[](size_t i) const
	{
		assert_lt(i, cur_);
		return list_[i];
	}
	inline __m128i &get(size_t i)
	{
		return operator[](i);
	}
	inline __m128i get(size_t i) const
	{
		return operator[](i);
	}
	__m128i *ptr() { return list_; }
	const __m128i *ptr() const { return list_; }
	int cat() const { return cat_; }
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
	__m128i *alloc(size_t sz)
	{
		__m128i *last_alloc_;
		try
		{
			last_alloc_ = new __m128i[sz + 2];
		}
		catch (std::bad_alloc &e)
		{
			std::cerr << "Error: Out of memory allocating " << sz << " __m128i's for DP matrix: '" << e.what() << "'" << std::endl;
			throw e;
		}
		this->last_alloc_ = last_alloc_;
		__m128i *tmp = last_alloc_;
		size_t tmpint = (size_t)tmp;
		if ((tmpint & 0xf) != 0)
		{
			tmpint += 15;
			tmpint &= (~0xf);
			tmp = reinterpret_cast<__m128i *>(tmpint);
		}
		assert_eq(0, (tmpint & 0xf));
		assert(tmp != NULL);
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#endif
		return tmp;
	}
	void free()
	{
		if (list_ != NULL)
		{
			delete[] last_alloc_;
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
		__m128i *prev_last_alloc = last_alloc_;
		__m128i *tmp = alloc(newsz);
		assert(tmp != NULL);
		size_t cur = cur_;
		if (list_ != NULL)
		{
			for (size_t i = 0; i < cur_; i++)
			{
				tmp[i] = list_[i];
			}
			__m128i *current_last_alloc = last_alloc_;
			last_alloc_ = prev_last_alloc;
			free();
			last_alloc_ = current_last_alloc;
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
		__m128i *tmp = alloc(newsz);
		assert(tmp != NULL);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}
	int cat_;
	__m128i *last_alloc_;
	__m128i *list_;
	size_t sz_;
	size_t cur_;
};
struct CpQuad
{
	CpQuad() { reset(); }
	void reset() { sc[0] = sc[1] = sc[2] = sc[3] = 0; }
	bool operator==(const CpQuad &o) const
	{
		return sc[0] == o.sc[0] &&
			   sc[1] == o.sc[1] &&
			   sc[2] == o.sc[2] &&
			   sc[3] == o.sc[3];
	}
	int16_t sc[4];
};
class Checkpointer
{
public:
	Checkpointer() { reset(); }
	void init(
		size_t nrow,
		size_t ncol, size_t perpow2, int64_t perfectScore, bool is8, bool doTri, bool local, bool debug)
	{
		assert_gt(perpow2, 0);
		nrow_ = nrow;
		ncol_ = ncol;
		perpow2_ = perpow2;
		per_ = 1 << perpow2;
		lomask_ = ~(0xffffffff << perpow2);
		perf_ = perfectScore;
		local_ = local;
		ndiag_ = (ncol + nrow - 1 + 1) / per_;
		locol_ = MAX_SIZE_T;
		hicol_ = MIN_SIZE_T;
		debug_ = true;
		commitMap_.clear();
		firstCommit_ = true;
		size_t perword = (is8 ? 16 : 8);
		is8_ = is8;
		niter_ = ((nrow_ + perword - 1) / perword);
		if (doTri)
		{
			qdiag1s_.resize(ndiag_ * nrow_);
			qdiag2s_.resize(ndiag_ * nrow_);
		}
		else
		{
			qrows_.resize((nrow_ / per_) * ncol_);
			qcols_.resize((ncol_ / per_) * (niter_ << 2));
		}
		if (debug_)
		{
			qcolsD_.resize(ncol_ * (niter_ << 2));
		}
	}
	bool debug() const { return debug_; }
	int64_t debugCell(size_t row, size_t col, int hef) const
	{
		assert(debug_);
		const __m128i *ptr = qcolsD_.ptr() + hef;
		ptr += ((col * niter_) << 2);
		size_t mod = row % niter_;
		size_t div = row / niter_;
		ptr += (mod << 2);
		int16_t sc = (is8_ ? ((uint8_t *)ptr)[div] : ((int16_t *)ptr)[div]);
		int64_t asc = MIN_I64;
		if (is8_)
		{
			if (local_)
			{
				asc = sc;
			}
			else
			{
				if (sc == 0)
					asc = MIN_I64;
				else
					asc = sc - 0xff;
			}
		}
		else
		{
			if (local_)
			{
				asc = sc + 0x8000;
			}
			else
			{
				if (sc != MIN_I16)
					asc = sc - 0x7fff;
			}
		}
		return asc;
	}
	bool isCheckpointed(size_t row, size_t col) const
	{
		assert_leq(col, hicol_);
		assert_geq(col, locol_);
		size_t mod = (row + col) & lomask_;
		assert_lt(mod, per_);
		return mod >= per_ - 2;
	}
	inline int64_t scoreTriangle(size_t row, size_t col, int hef) const
	{
		assert(isCheckpointed(row, col));
		bool diag1 = ((row + col) & lomask_) == per_ - 2;
		size_t off = (row + col) >> perpow2_;
		if (diag1)
		{
			if (qdiag1s_[off * nrow_ + row].sc[hef] == MIN_I16)
			{
				return MIN_I64;
			}
			else
			{
				return qdiag1s_[off * nrow_ + row].sc[hef];
			}
		}
		else
		{
			if (qdiag2s_[off * nrow_ + row].sc[hef] == MIN_I16)
			{
				return MIN_I64;
			}
			else
			{
				return qdiag2s_[off * nrow_ + row].sc[hef];
			}
		}
	}
	inline int64_t scoreSquare(size_t row, size_t col, int hef) const
	{
		if ((row & lomask_) == lomask_ && hef != 1)
		{
			int64_t sc = qrows_[(row >> perpow2_) * ncol_ + col].sc[hef];
			if (sc == MIN_I16)
				return MIN_I64;
			return sc;
		}
		hef--;
		if (hef == -1)
			hef = 2;
		assert_eq(lomask_, (col & lomask_));
		const __m128i *ptr = qcols_.ptr() + hef;
		ptr += (((col >> perpow2_) * niter_) << 2);
		size_t mod = row % niter_;
		size_t div = row / niter_;
		ptr += (mod << 2);
		int16_t sc = (is8_ ? ((uint8_t *)ptr)[div] : ((int16_t *)ptr)[div]);
		int64_t asc = MIN_I64;
		if (is8_)
		{
			if (local_)
			{
				asc = sc;
			}
			else
			{
				if (sc == 0)
					asc = MIN_I64;
				else
					asc = sc - 0xff;
			}
		}
		else
		{
			if (local_)
			{
				asc = sc + 0x8000;
			}
			else
			{
				if (sc != MIN_I16)
					asc = sc - 0x7fff;
			}
		}
		return asc;
	}
	void commitCol(__m128i *pvH, __m128i *pvE, __m128i *pvF, size_t coli);
	void reset()
	{
		perpow2_ = per_ = lomask_ = nrow_ = ncol_ = 0;
		local_ = false;
		niter_ = ndiag_ = locol_ = hicol_ = 0;
		perf_ = 0;
		firstCommit_ = true;
		is8_ = debug_ = false;
	}
	bool inited() const
	{
		return nrow_ > 0;
	}
	size_t per() const { return per_; }
	size_t perpow2() const { return perpow2_; }
	size_t lomask() const { return lomask_; }
	size_t locol() const { return locol_; }
	size_t hicol() const { return hicol_; }
	size_t nrow() const { return nrow_; }
	size_t ncol() const { return ncol_; }
	const CpQuad *qdiag1sPtr() const { return qdiag1s_.ptr(); }
	const CpQuad *qdiag2sPtr() const { return qdiag2s_.ptr(); }
	size_t perpow2_;
	size_t per_;
	size_t lomask_;
	size_t nrow_;
	size_t ncol_;
	int64_t perf_;
	bool local_;
	size_t ndiag_;
	size_t locol_;
	size_t hicol_;
	EList<size_t> commitMap_;
	bool firstCommit_;
	EList<CpQuad> qdiag1s_;
	EList<CpQuad> qdiag2s_;
	EList<CpQuad> qrows_;
	bool is8_;
	size_t niter_;
	EList_m128i qcols_;
	bool debug_;
	EList_m128i qcolsD_;
};
#endif
