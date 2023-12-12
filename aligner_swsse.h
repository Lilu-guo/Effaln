#ifndef ALIGNER_SWSSE_H_
#define ALIGNER_SWSSE_H_
#include "ds.h"
#include "mem_ids.h"
#include "random_source.h"
#include "scoring.h"
#include "mask.h"
#include "sse_util.h"
#include <strings.h>
struct SSEMetrics
{
	SSEMetrics() : mutex_m() { reset(); }
	void clear() { reset(); }
	void reset()
	{
		dp = dpsat = dpfail = dpsucc =
			col = cell = inner = fixup =
				gathsol = bt = btfail = btsucc = btcell =
					corerej = nrej = 0;
	}
	void merge(const SSEMetrics &o)
	{
		dp += o.dp;
		dpsat += o.dpsat;
		dpfail += o.dpfail;
		dpsucc += o.dpsucc;
		col += o.col;
		cell += o.cell;
		inner += o.inner;
		fixup += o.fixup;
		gathsol += o.gathsol;
		bt += o.bt;
		btfail += o.btfail;
		btsucc += o.btsucc;
		btcell += o.btcell;
		corerej += o.corerej;
		nrej += o.nrej;
	}
	uint64_t dp;
	uint64_t dpsat;
	uint64_t dpfail;
	uint64_t dpsucc;
	uint64_t col;
	uint64_t cell;
	uint64_t inner;
	uint64_t fixup;
	uint64_t gathsol;
	uint64_t bt;
	uint64_t btfail;
	uint64_t btsucc;
	uint64_t btcell;
	uint64_t corerej;
	uint64_t nrej;
	MUTEX_T mutex_m;
};
struct SSEMatrix
{
	const static size_t E = 0;
	const static size_t F = 1;
	const static size_t H = 2;
	const static size_t TMP = 3;
	SSEMatrix(int cat = 0) : nvecPerCell_(4), matbuf_(cat) {}
	inline __m128i *ptr()
	{
		assert(inited_);
		return matbuf_.ptr();
	}
	inline __m128i *evec(size_t row, size_t col)
	{
		assert_lt(row, nvecrow_);
		assert_lt(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + E;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}
	inline __m128i *evecUnsafe(size_t row, size_t col)
	{
		assert_lt(row, nvecrow_);
		assert_leq(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + E;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}
	inline __m128i *fvec(size_t row, size_t col)
	{
		assert_lt(row, nvecrow_);
		assert_lt(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + F;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}
	inline __m128i *hvec(size_t row, size_t col)
	{
		assert_lt(row, nvecrow_);
		assert_lt(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + H;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}
	inline __m128i *tmpvec(size_t row, size_t col)
	{
		assert_lt(row, nvecrow_);
		assert_lt(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + TMP;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}
	inline __m128i *tmpvecUnsafe(size_t row, size_t col)
	{
		assert_lt(row, nvecrow_);
		assert_leq(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + TMP;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}
	void init(
		size_t nrow,
		size_t ncol,
		size_t wperv);
	inline size_t colstride() const { return colstride_; }
	inline size_t rowstride() const { return rowstride_; }
	int eltSlow(size_t row, size_t col, size_t mat) const;
	inline int elt(size_t row, size_t col, size_t mat) const
	{
		assert(inited_);
		assert_lt(row, nrow_);
		assert_lt(col, ncol_);
		assert_lt(mat, 3);
		size_t rowelt = row / nvecrow_;
		size_t rowvec = row % nvecrow_;
		size_t eltvec = (col * colstride_) + (rowvec * rowstride_) + mat;
		assert_lt(eltvec, matbuf_.size());
		if (wperv_ == 16)
		{
			return (int)((uint8_t *)(matbuf_.ptr() + eltvec))[rowelt];
		}
		else
		{
			assert_eq(8, wperv_);
			return (int)((int16_t *)(matbuf_.ptr() + eltvec))[rowelt];
		}
	}
	inline int eelt(size_t row, size_t col) const
	{
		return elt(row, col, E);
	}
	inline int felt(size_t row, size_t col) const
	{
		return elt(row, col, F);
	}
	inline int helt(size_t row, size_t col) const
	{
		return elt(row, col, H);
	}
	inline bool reportedThrough(
		size_t row,
		size_t col) const
	{
		return (masks_[row][col] & (1 << 0)) != 0;
	}
	inline void setReportedThrough(
		size_t row,
		size_t col)
	{
		masks_[row][col] |= (1 << 0);
	}
	bool isHMaskSet(
		size_t row,
		size_t col) const;
	void hMaskSet(
		size_t row,
		size_t col, int mask);
	bool isEMaskSet(
		size_t row,
		size_t col) const;
	void eMaskSet(
		size_t row,
		size_t col, int mask);
	bool isFMaskSet(
		size_t row,
		size_t col) const;
	void fMaskSet(
		size_t row,
		size_t col, int mask);
	void analyzeCell(
		size_t row,
		size_t col,
		size_t ct,
		int refc,
		int readc,
		int readq,
		const Scoring &sc,
		RandomSource &rand,
		bool &empty,
		int &cur,
		bool &branch,
		bool &canMoveThru,
		bool &reportedThru);
	void initMasks();
	size_t nrow() const
	{
		return nrow_;
	}
	size_t ncol() const
	{
		return ncol_;
	}
	void resetRow(size_t i)
	{
		assert(!reset_[i]);
		masks_[i].resizeNoCopy(ncol_);
		masks_[i].fillZero();
		reset_[i] = true;
	}
	bool inited_;
	size_t nrow_;
	size_t ncol_;
	size_t nvecrow_;
	size_t nveccol_;
	size_t wperv_;
	size_t vecshift_;
	size_t nvecPerCol_;
	size_t nvecPerCell_;
	size_t colstride_;
	size_t rowstride_;
	EList_m128i matbuf_;
	ELList<uint16_t> masks_;
	EList<bool> reset_;
};
struct SSEData
{
	SSEData(int cat = 0) : profbuf_(cat), mat_(cat) {}
	EList_m128i profbuf_;
	EList_m128i vecbuf_;
	size_t qprofStride_;
	size_t gbarStride_;
	SSEMatrix mat_;
	size_t maxPen_;
	size_t maxBonus_;
	size_t lastIter_;
	size_t lastWord_;
	int bias_;
};
inline bool SSEMatrix::isHMaskSet(
	size_t row,
	size_t col) const
{
	return (masks_[row][col] & (1 << 1)) != 0;
}
inline void SSEMatrix::hMaskSet(
	size_t row,
	size_t col, int mask)
{
	assert_lt(mask, 32);
	masks_[row][col] &= ~(31 << 1);
	masks_[row][col] |= (1 << 1 | mask << 2);
}
inline bool SSEMatrix::isEMaskSet(
	size_t row,
	size_t col) const
{
	return (masks_[row][col] & (1 << 7)) != 0;
}
inline void SSEMatrix::eMaskSet(
	size_t row,
	size_t col, int mask)
{
	assert_lt(mask, 4);
	masks_[row][col] &= ~(7 << 7);
	masks_[row][col] |= (1 << 7 | mask << 8);
}
inline bool SSEMatrix::isFMaskSet(
	size_t row,
	size_t col) const
{
	return (masks_[row][col] & (1 << 10)) != 0;
}
inline void SSEMatrix::fMaskSet(
	size_t row,
	size_t col, int mask)
{
	assert_lt(mask, 4);
	masks_[row][col] &= ~(7 << 10);
	masks_[row][col] |= (1 << 10 | mask << 11);
}
#define ROWSTRIDE_2COL 4
#define ROWSTRIDE 4
#endif
