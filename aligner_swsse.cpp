#include <string.h>
#include "aligner_sw_common.h"
#include "aligner_swsse.h"
void SSEMatrix::init(
	size_t nrow,
	size_t ncol,
	size_t wperv)
{
	nrow_ = nrow;
	ncol_ = ncol;
	wperv_ = wperv;
	nvecPerCol_ = (nrow + (wperv - 1)) / wperv;
	try
	{
		matbuf_.resizeNoCopy((ncol + 1) * nvecPerCell_ * nvecPerCol_);
	}
	catch (exception &e)
	{
		cerr << "Tried to allocate DP matrix with " << (ncol + 1)
			 << " columns, " << nvecPerCol_
			 << " vectors per column, and and " << nvecPerCell_
			 << " vectors per cell" << endl;
		throw e;
	}
	assert(wperv_ == 8 || wperv_ == 16);
	vecshift_ = (wperv_ == 8) ? 3 : 4;
	nvecrow_ = (nrow + (wperv_ - 1)) >> vecshift_;
	nveccol_ = ncol;
	colstride_ = nvecPerCol_ * nvecPerCell_;
	rowstride_ = nvecPerCell_;
	inited_ = true;
}
void SSEMatrix::initMasks()
{
	assert_gt(nrow_, 0);
	assert_gt(ncol_, 0);
	masks_.resize(nrow_);
	reset_.resizeNoCopy(nrow_);
	reset_.fill(false);
}
int SSEMatrix::eltSlow(size_t row, size_t col, size_t mat) const
{
	assert_lt(row, nrow_);
	assert_lt(col, ncol_);
	assert_leq(mat, 3);
	size_t rowelt = row / nvecrow_;
	size_t rowvec = row % nvecrow_;
	size_t eltvec = (col * colstride_) + (rowvec * rowstride_) + mat;
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
