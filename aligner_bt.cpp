#include "aligner_bt.h"
#include "mask.h"
using namespace std;
#define CHECK_ROW_COL(rowc, colc)                                  \
	if (rowc >= 0 && colc >= 0)                                    \
	{                                                              \
		if (!sawcell_[colc].insert(rowc))                          \
		{                                                          \
			abort = true;                                          \
			return;                                                \
		}                                                          \
		assert(local || prob_.cper_->debugCell(rowc, colc, hefc)); \
	}
void BtBranchTracer::triangleFill(
	int64_t rw,
	int64_t cl,
	int hef,
	TAlScore targ,
	TAlScore targ_final,
	RandomSource &rnd,
	int64_t &row_new,
	int64_t &col_new,
	int &hef_new,
	TAlScore &targ_new,
	bool &done,
	bool &abort)
{
	assert_geq(rw, 0);
	assert_geq(cl, 0);
	assert_range(0, 2, hef);
	assert_lt(rw, (int64_t)prob_.qrylen_);
	assert_lt(cl, (int64_t)prob_.reflen_);
	assert(prob_.usecp_ && prob_.fill_);
	int64_t row = rw, col = cl;
	const int64_t colmin = 0;
	const int64_t rowmin = 0;
	const int64_t colmax = prob_.reflen_ - 1;
	const int64_t rowmax = prob_.qrylen_ - 1;
	assert_leq(prob_.reflen_, (TRefOff)sawcell_.size());
	assert_leq(col, (int64_t)prob_.cper_->hicol());
	assert_geq(col, (int64_t)prob_.cper_->locol());
	assert_geq(prob_.cper_->per(), 2);
	size_t mod = (row + col) & prob_.cper_->lomask();
	assert_lt(mod, prob_.cper_->per());
	size_t depth = mod + 1;
	assert_leq(depth, prob_.cper_->per());
	size_t breadth = depth;
	tri_.resize(depth);
	for (size_t i = 0; i < depth; i++)
	{
		tri_[i].resize(breadth - i);
	}
	bool upperleft = false;
	size_t off = (row + col) >> prob_.cper_->perpow2();
	if (off == 0)
	{
		upperleft = true;
	}
	else
	{
		off--;
	}
	const TAlScore sc_rdo = prob_.sc_->readGapOpen();
	const TAlScore sc_rde = prob_.sc_->readGapExtend();
	const TAlScore sc_rfo = prob_.sc_->refGapOpen();
	const TAlScore sc_rfe = prob_.sc_->refGapExtend();
	const bool local = !prob_.sc_->monotone;
	int64_t row_lo = row - (int64_t)mod;
	const CpQuad *prev2 = NULL, *prev1 = NULL;
	if (!upperleft)
	{
		prev2 = prob_.cper_->qdiag1sPtr() + (off * prob_.cper_->nrow() + row_lo - 1);
		prev1 = prob_.cper_->qdiag2sPtr() + (off * prob_.cper_->nrow() + row_lo - 1);
#ifndef NDEBUG
		if (row >= (int64_t)mod)
		{
			size_t rowc = row - mod, colc = col;
			if (rowc > 0 && prob_.cper_->isCheckpointed(rowc - 1, colc))
			{
				TAlScore al = prev1[0].sc[0];
				if (al == MIN_I16)
					al = MIN_I64;
				assert_eq(prob_.cper_->scoreTriangle(rowc - 1, colc, 0), al);
			}
			if (rowc > 0 && colc > 0 && prob_.cper_->isCheckpointed(rowc - 1, colc - 1))
			{
				TAlScore al = prev2[0].sc[0];
				if (al == MIN_I16)
					al = MIN_I64;
				assert_eq(prob_.cper_->scoreTriangle(rowc - 1, colc - 1, 0), al);
			}
		}
#endif
	}
	for (size_t i = 0; i < depth; i++)
	{
		CpQuad *cur = tri_[i].ptr();
		CpQuad *curc = cur;
		size_t doff = mod - i;
		int64_t rowc = row - doff;
		int64_t colc = col;
		size_t neval = 0;
		ASSERT_ONLY(const CpQuad *last = NULL);
		for (size_t j = 0; j < breadth; j++)
		{
			if (rowc >= rowmin && rowc <= rowmax &&
				colc >= colmin && colc <= colmax)
			{
				neval++;
				int64_t fromend = prob_.qrylen_ - rowc - 1;
				bool allowGaps = fromend >= prob_.sc_->gapbar && rowc >= prob_.sc_->gapbar;
				int qc = prob_.qry_[rowc];
				int qq = prob_.qual_[rowc];
				assert_geq(qq, 33);
				int rc = prob_.ref_[colc];
				assert_range(0, 16, rc);
				int16_t sc_diag = prob_.sc_->score(qc, rc, qq - 33);
				int16_t sc_h_up = MIN_I16;
				int16_t sc_f_up = MIN_I16;
				int16_t sc_h_lf = MIN_I16;
				int16_t sc_e_lf = MIN_I16;
				if (allowGaps)
				{
					if (rowc > 0)
					{
						assert(local || prev1[j + 0].sc[2] < 0);
						if (prev1[j + 0].sc[0] > MIN_I16)
						{
							sc_h_up = prev1[j + 0].sc[0] - sc_rfo;
							if (local)
								sc_h_up = max<int16_t>(sc_h_up, 0);
						}
						if (prev1[j + 0].sc[2] > MIN_I16)
						{
							sc_f_up = prev1[j + 0].sc[2] - sc_rfe;
							if (local)
								sc_f_up = max<int16_t>(sc_f_up, 0);
						}
#ifndef NDEBUG
						TAlScore hup = prev1[j + 0].sc[0];
						TAlScore fup = prev1[j + 0].sc[2];
						if (hup == MIN_I16)
							hup = MIN_I64;
						if (fup == MIN_I16)
							fup = MIN_I64;
						if (local)
						{
							hup = max<int16_t>(hup, 0);
							fup = max<int16_t>(fup, 0);
						}
						if (prob_.cper_->isCheckpointed(rowc - 1, colc))
						{
							assert_eq(hup, prob_.cper_->scoreTriangle(rowc - 1, colc, 0));
							assert_eq(fup, prob_.cper_->scoreTriangle(rowc - 1, colc, 2));
						}
#endif
					}
					if (colc > 0)
					{
						assert(local || prev1[j + 1].sc[1] < 0);
						if (prev1[j + 1].sc[0] > MIN_I16)
						{
							sc_h_lf = prev1[j + 1].sc[0] - sc_rdo;
							if (local)
								sc_h_lf = max<int16_t>(sc_h_lf, 0);
						}
						if (prev1[j + 1].sc[1] > MIN_I16)
						{
							sc_e_lf = prev1[j + 1].sc[1] - sc_rde;
							if (local)
								sc_e_lf = max<int16_t>(sc_e_lf, 0);
						}
#ifndef NDEBUG
						TAlScore hlf = prev1[j + 1].sc[0];
						TAlScore elf = prev1[j + 1].sc[1];
						if (hlf == MIN_I16)
							hlf = MIN_I64;
						if (elf == MIN_I16)
							elf = MIN_I64;
						if (local)
						{
							hlf = max<int16_t>(hlf, 0);
							elf = max<int16_t>(elf, 0);
						}
						if (prob_.cper_->isCheckpointed(rowc, colc - 1))
						{
							assert_eq(hlf, prob_.cper_->scoreTriangle(rowc, colc - 1, 0));
							assert_eq(elf, prob_.cper_->scoreTriangle(rowc, colc - 1, 1));
						}
#endif
					}
				}
				assert(rowc <= 1 || colc <= 0 || prev2 != NULL);
				int16_t sc_h_dg = ((rowc > 0 && colc > 0) ? prev2[j + 0].sc[0] : 0);
				if (colc == 0 && rowc > 0 && !local)
				{
					sc_h_dg = MIN_I16;
				}
				if (sc_h_dg > MIN_I16)
				{
					sc_h_dg += sc_diag;
				}
				if (local)
					sc_h_dg = max<int16_t>(sc_h_dg, 0);
				int mask = 0;
				int16_t sc_best = sc_h_dg;
				if (sc_h_dg > MIN_I64)
				{
					mask = 1;
				}
				if (colc > 0 && sc_h_lf >= sc_best && sc_h_lf > MIN_I64)
				{
					if (sc_h_lf > sc_best)
						mask = 0;
					mask |= 2;
					sc_best = sc_h_lf;
				}
				if (colc > 0 && sc_e_lf >= sc_best && sc_e_lf > MIN_I64)
				{
					if (sc_e_lf > sc_best)
						mask = 0;
					mask |= 4;
					sc_best = sc_e_lf;
				}
				if (rowc > 0 && sc_h_up >= sc_best && sc_h_up > MIN_I64)
				{
					if (sc_h_up > sc_best)
						mask = 0;
					mask |= 8;
					sc_best = sc_h_up;
				}
				if (rowc > 0 && sc_f_up >= sc_best && sc_f_up > MIN_I64)
				{
					if (sc_f_up > sc_best)
						mask = 0;
					mask |= 16;
					sc_best = sc_f_up;
				}
				int16_t sc_e_best = sc_h_lf;
				if (colc > 0)
				{
					if (sc_h_lf >= sc_e_lf && sc_h_lf > MIN_I64)
					{
						if (sc_h_lf == sc_e_lf)
						{
							mask |= 64;
						}
						mask |= 32;
					}
					else if (sc_e_lf > MIN_I64)
					{
						sc_e_best = sc_e_lf;
						mask |= 64;
					}
				}
				if (sc_e_best > sc_best)
				{
					sc_best = sc_e_best;
					mask &= ~31;
				}
				int16_t sc_f_best = sc_h_up;
				if (rowc > 0)
				{
					if (sc_h_up >= sc_f_up && sc_h_up > MIN_I64)
					{
						if (sc_h_up == sc_f_up)
						{
							mask |= 256;
						}
						mask |= 128;
					}
					else if (sc_f_up > MIN_I64)
					{
						sc_f_best = sc_f_up;
						mask |= 256;
					}
				}
				if (sc_f_best > sc_best)
				{
					sc_best = sc_f_best;
					mask &= ~127;
				}
				assert(!prob_.sc_->monotone || sc_best <= 0);
				assert(!prob_.sc_->monotone || sc_e_best <= 0);
				assert(!prob_.sc_->monotone || sc_f_best <= 0);
				curc->sc[0] = sc_best;
				assert(local || sc_e_best < 0);
				assert(local || sc_f_best < 0);
				assert(!local || sc_e_best >= 0 || sc_e_best == MIN_I16);
				assert(!local || sc_f_best >= 0 || sc_f_best == MIN_I16);
				curc->sc[1] = sc_e_best;
				curc->sc[2] = sc_f_best;
				curc->sc[3] = mask;
				ASSERT_ONLY(last = curc);
#ifndef NDEBUG
				if (prob_.cper_->isCheckpointed(rowc, colc))
				{
					if (local)
					{
						sc_e_best = max<int16_t>(sc_e_best, 0);
						sc_f_best = max<int16_t>(sc_f_best, 0);
					}
					TAlScore sc_best64 = sc_best;
					if (sc_best == MIN_I16)
						sc_best64 = MIN_I64;
					TAlScore sc_e_best64 = sc_e_best;
					if (sc_e_best == MIN_I16)
						sc_e_best64 = MIN_I64;
					TAlScore sc_f_best64 = sc_f_best;
					if (sc_f_best == MIN_I16)
						sc_f_best64 = MIN_I64;
					assert_eq(prob_.cper_->scoreTriangle(rowc, colc, 0), sc_best64);
					assert_eq(prob_.cper_->scoreTriangle(rowc, colc, 1), sc_e_best64);
					assert_eq(prob_.cper_->scoreTriangle(rowc, colc, 2), sc_f_best64);
				}
#endif
			}
			assert_lt(rowc, (int64_t)prob_.qrylen_);
			rowc++;
			colc--;
			curc++;
		}
		if (i == depth - 1)
		{
			assert(last != NULL);
			assert_eq(1, neval);
			assert_neq(0, last->sc[3]);
			assert_eq(targ, last->sc[hef]);
		}
		else
		{
			breadth--;
			prev2 = prev1 + 1;
			prev1 = cur;
		}
	}
	int64_t rowc = row, colc = col;
	size_t curid;
	int hefc = hef;
	if (bs_.empty())
	{
		CHECK_ROW_COL(rowc, colc);
		curid = bs_.alloc();
		assert_eq(0, curid);
		Edit e;
		bs_[curid].init(
			prob_,
			0,
			0,
			0,
			rowc,
			colc,
			e,
			0,
			true,
			false);
		bs_[curid].len_ = 0;
	}
	else
	{
		curid = bs_.size() - 1;
	}
	size_t idx_orig = (row + col) >> prob_.cper_->perpow2();
	while (true)
	{
		size_t mod = (rowc + colc) & prob_.cper_->lomask();
		assert_lt(mod, prob_.cper_->per());
		CpQuad *cur = tri_[mod].ptr();
		int64_t row_off = rowc - row_lo - mod;
		assert(!local || cur[row_off].sc[0] > 0);
		assert_geq(row_off, 0);
		int mask = cur[row_off].sc[3];
		assert_gt(mask, 0);
		int sel = -1;
		if (hefc == 0)
		{
			if ((mask & 1) != 0)
			{
				sel = 0;
			}
			else if ((mask & 8) != 0)
			{
				sel = 3;
			}
			else if ((mask & 16) != 0)
			{
				sel = 4;
			}
			else if ((mask & 2) != 0)
			{
				sel = 1;
			}
			else if ((mask & 4) != 0)
			{
				sel = 2;
			}
		}
		else if (hefc == 1)
		{
			if ((mask & 32) != 0)
			{
				sel = 5;
			}
			else if ((mask & 64) != 0)
			{
				sel = 6;
			}
		}
		else
		{
			assert_eq(2, hefc);
			if ((mask & 128) != 0)
			{
				sel = 7;
			}
			else if ((mask & 256) != 0)
			{
				sel = 8;
			}
		}
		assert_geq(sel, 0);
		int qc = prob_.qry_[rowc], qq = prob_.qual_[rowc];
		int rc = prob_.ref_[colc];
		assert_range(0, 16, rc);
		if (sel == 0)
		{
			assert_geq(rowc, 0);
			assert_geq(colc, 0);
			TAlScore scd = prob_.sc_->score(qc, rc, qq - 33);
			if ((rc & (1 << qc)) == 0)
			{
				size_t id = curid;
				bool empty = (bs_[curid].len_ == 0 && curid == 0);
				if (!empty)
				{
					id = bs_.alloc();
				}
				Edit e((int)rowc, mask2dna[rc], "ACGTN"[qc], EDIT_TYPE_MM);
				assert_lt(scd, 0);
				TAlScore score_en = bs_[curid].score_st_ + scd;
				bs_[id].init(
					prob_,
					curid,
					-scd,
					score_en,
					rowc,
					colc,
					e,
					hefc,
					empty,
					false);
				curid = id;
			}
			else
			{
				bs_[curid].score_st_ += prob_.sc_->match();
				bs_[curid].len_++;
				assert_leq((int64_t)bs_[curid].len_, bs_[curid].row_ + 1);
			}
			rowc--;
			colc--;
			assert(local || bs_[curid].score_st_ >= targ_final);
			hefc = 0;
		}
		else if ((sel >= 1 && sel <= 2) || (sel >= 5 && sel <= 6))
		{
			assert_gt(colc, 0);
			size_t id = bs_.alloc();
			Edit e((int)rowc + 1, mask2dna[rc], '-', EDIT_TYPE_READ_GAP);
			TAlScore gapp = prob_.sc_->readGapOpen();
			if (bs_[curid].len_ == 0 && bs_[curid].e_.inited() && bs_[curid].e_.isReadGap())
			{
				gapp = prob_.sc_->readGapExtend();
			}
			TAlScore score_en = bs_[curid].score_st_ - gapp;
			bs_[id].init(
				prob_,
				curid,
				gapp,
				score_en,
				rowc,
				colc - 1,
				e,
				hefc,
				false,
				false);
			colc--;
			curid = id;
			assert(local || bs_[curid].score_st_ >= targ_final);
			if (sel == 1 || sel == 5)
			{
				hefc = 0;
			}
			else
			{
				hefc = 1;
			}
		}
		else
		{
			assert_gt(rowc, 0);
			size_t id = bs_.alloc();
			Edit e((int)rowc, '-', "ACGTN"[qc], EDIT_TYPE_REF_GAP);
			TAlScore gapp = prob_.sc_->refGapOpen();
			if (bs_[curid].len_ == 0 && bs_[curid].e_.inited() && bs_[curid].e_.isRefGap())
			{
				gapp = prob_.sc_->refGapExtend();
			}
			TAlScore score_en = bs_[curid].score_st_ - gapp;
			bs_[id].init(
				prob_,
				curid,
				gapp,
				score_en,
				rowc - 1,
				colc,
				e,
				hefc,
				false,
				false);
			rowc--;
			curid = id;
			if (sel == 3 || sel == 7)
			{
				hefc = 0;
			}
			else
			{
				hefc = 2;
			}
		}
		CHECK_ROW_COL(rowc, colc);
		size_t mod_new = (rowc + colc) & prob_.cper_->lomask();
		size_t idx = (rowc + colc) >> prob_.cper_->perpow2();
		assert_lt(mod_new, prob_.cper_->per());
		int64_t row_off_new = rowc - row_lo - mod_new;
		CpQuad *cur_new = NULL;
		if (colc >= 0 && rowc >= 0 && idx == idx_orig)
		{
			cur_new = tri_[mod_new].ptr();
		}
		bool hit_new_tri = (idx < idx_orig && colc >= 0 && rowc >= 0);
		if (colc < 0 || rowc < 0 ||
			(cur_new != NULL && (local && cur_new[row_off_new].sc[0] == 0)))
		{
			done = true;
			assert(bs_[curid].isSolution(prob_));
			addSolution(curid);
#ifndef NDEBUG
			size_t cur = curid;
			if (!bs_[cur].root_)
			{
				size_t next = bs_[cur].parentId_;
				while (!bs_[next].root_)
				{
					assert_neq(cur, next);
					if (bs_[next].len_ != 0 || bs_[cur].len_ == 0)
					{
						assert(!bs_[cur].overlap(prob_, bs_[next]));
					}
					cur = next;
					next = bs_[cur].parentId_;
				}
			}
#endif
			return;
		}
		if (hit_new_tri)
		{
			assert(rowc < 0 || colc < 0 || prob_.cper_->isCheckpointed(rowc, colc));
			row_new = rowc;
			col_new = colc;
			hef_new = hefc;
			done = false;
			if (rowc < 0 || colc < 0)
			{
				assert(local);
				targ_new = 0;
			}
			else
			{
				targ_new = prob_.cper_->scoreTriangle(rowc, colc, hefc);
			}
			if (local && targ_new == 0)
			{
				done = true;
				assert(bs_[curid].isSolution(prob_));
				addSolution(curid);
			}
			assert((row_new >= 0 && col_new >= 0) || done);
			return;
		}
	}
	assert(false);
}
#ifndef NDEBUG
#define DEBUG_CHECK(ss, row, col, hef)                            \
	{                                                             \
		if (prob_.cper_->debug() && row >= 0 && col >= 0)         \
		{                                                         \
			TAlScore s = ss;                                      \
			if (s == MIN_I16)                                     \
				s = MIN_I64;                                      \
			if (local && s < 0)                                   \
				s = 0;                                            \
			TAlScore deb = prob_.cper_->debugCell(row, col, hef); \
			if (local && deb < 0)                                 \
				deb = 0;                                          \
			assert_eq(s, deb);                                    \
		}                                                         \
	}
#else
#define DEBUG_CHECK(ss, row, col, hef)
#endif
void BtBranchTracer::squareFill(
	int64_t rw,
	int64_t cl,
	int hef,
	TAlScore targ,
	TAlScore targ_final,
	RandomSource &rnd,
	int64_t &row_new,
	int64_t &col_new,
	int &hef_new,
	TAlScore &targ_new,
	bool &done,
	bool &abort)
{
	assert_geq(rw, 0);
	assert_geq(cl, 0);
	assert_range(0, 2, hef);
	assert_lt(rw, (int64_t)prob_.qrylen_);
	assert_lt(cl, (int64_t)prob_.reflen_);
	assert(prob_.usecp_ && prob_.fill_);
	const bool is8_ = prob_.cper_->is8_;
	int64_t row = rw, col = cl;
	assert_leq(prob_.reflen_, (TRefOff)sawcell_.size());
	assert_leq(col, (int64_t)prob_.cper_->hicol());
	assert_geq(col, (int64_t)prob_.cper_->locol());
	assert_geq(prob_.cper_->per(), 2);
	size_t xmod = col & prob_.cper_->lomask();
	size_t ymod = row & prob_.cper_->lomask();
	size_t xdiv = col >> prob_.cper_->perpow2();
	size_t ydiv = row >> prob_.cper_->perpow2();
	size_t sq_ncol = xmod + 1, sq_nrow = ymod + 1;
	sq_.resize(sq_ncol * sq_nrow);
	bool upper = ydiv == 0;
	bool left = xdiv == 0;
	const TAlScore sc_rdo = prob_.sc_->readGapOpen();
	const TAlScore sc_rde = prob_.sc_->readGapExtend();
	const TAlScore sc_rfo = prob_.sc_->refGapOpen();
	const TAlScore sc_rfe = prob_.sc_->refGapExtend();
	const bool local = !prob_.sc_->monotone;
	const CpQuad *qup = NULL;
	const __m128i *qlf = NULL;
	size_t per = prob_.cper_->per_;
	ASSERT_ONLY(size_t nrow = prob_.cper_->nrow());
	size_t ncol = prob_.cper_->ncol();
	assert_eq(prob_.qrylen_, nrow);
	assert_eq(prob_.reflen_, (TRefOff)ncol);
	size_t niter = prob_.cper_->niter_;
	if (!upper)
	{
		qup = prob_.cper_->qrows_.ptr() + (ncol * (ydiv - 1)) + xdiv * per;
	}
	if (!left)
	{
		size_t off = (niter << 2) * (xdiv - 1);
		qlf = prob_.cper_->qcols_.ptr() + off;
	}
	size_t xedge = xdiv * per;
	size_t yedge = ydiv * per;
	size_t xi = xedge, yi = yedge;
	size_t ii = 0;
	size_t m128mod = yi % prob_.cper_->niter_;
	size_t m128div = yi / prob_.cper_->niter_;
	int16_t sc_h_dg_lastrow = MIN_I16;
	for (size_t i = 0; i <= ymod; i++, yi++)
	{
		assert_lt(yi, nrow);
		xi = xedge;
		size_t fromend = prob_.qrylen_ - yi - 1;
		bool allowGaps = fromend >= (size_t)prob_.sc_->gapbar && yi >= (size_t)prob_.sc_->gapbar;
		int qc = prob_.qry_[yi], qq = prob_.qual_[yi];
		assert_geq(qq, 33);
		int16_t sc_h_lf_last = MIN_I16;
		int16_t sc_e_lf_last = MIN_I16;
		for (size_t j = 0; j <= xmod; j++, xi++)
		{
			assert_lt(xi, ncol);
			int rc = prob_.ref_[xi];
			assert_range(0, 16, rc);
			int16_t sc_diag = prob_.sc_->score(qc, rc, qq - 33);
			int16_t sc_h_up = MIN_I16, sc_f_up = MIN_I16,
					sc_h_lf = MIN_I16, sc_e_lf = MIN_I16,
					sc_h_dg = MIN_I16;
			int16_t sc_h_up_c = MIN_I16, sc_f_up_c = MIN_I16,
					sc_h_lf_c = MIN_I16, sc_e_lf_c = MIN_I16,
					sc_h_dg_c = MIN_I16;
			if (yi == 0)
			{
				sc_h_dg = 0;
			}
			else if (xi == 0)
			{
				if (local)
				{
					sc_h_dg = 0;
				}
			}
			else if (i == 0 && j == 0)
			{
				sc_h_dg = qup[-1].sc[0];
			}
			else if (j == 0)
			{
				sc_h_dg = sc_h_dg_lastrow;
			}
			else
			{
				sc_h_dg = qup[j - 1].sc[0];
			}
			if (yi > 0 && xi > 0)
				DEBUG_CHECK(sc_h_dg, yi - 1, xi - 1, 2);
			if (j == 0 && xi > 0)
			{
				if (is8_)
				{
					size_t vecoff = (m128mod << 6) + m128div;
					sc_e_lf = ((uint8_t *)(qlf + 0))[vecoff];
					sc_h_lf = ((uint8_t *)(qlf + 2))[vecoff];
					if (local)
					{
					}
					else
					{
						if (sc_h_lf == 0)
							sc_h_lf = MIN_I16;
						else
							sc_h_lf -= 0xff;
						if (sc_e_lf == 0)
							sc_e_lf = MIN_I16;
						else
							sc_e_lf -= 0xff;
					}
				}
				else
				{
					size_t vecoff = (m128mod << 5) + m128div;
					sc_e_lf = ((int16_t *)(qlf + 0))[vecoff];
					sc_h_lf = ((int16_t *)(qlf + 2))[vecoff];
					if (local)
					{
						sc_h_lf += 0x8000;
						assert_geq(sc_h_lf, 0);
						sc_e_lf += 0x8000;
						assert_geq(sc_e_lf, 0);
					}
					else
					{
						if (sc_h_lf != MIN_I16)
							sc_h_lf -= 0x7fff;
						if (sc_e_lf != MIN_I16)
							sc_e_lf -= 0x7fff;
					}
				}
				DEBUG_CHECK(sc_e_lf, yi, xi - 1, 0);
				DEBUG_CHECK(sc_h_lf, yi, xi - 1, 2);
				sc_h_dg_lastrow = sc_h_lf;
			}
			if (allowGaps)
			{
				if (j == 0 && xi > 0)
				{
					sc_h_lf_c = sc_h_lf;
					sc_e_lf_c = sc_e_lf;
					if (sc_h_lf_c != MIN_I16)
						sc_h_lf_c -= sc_rdo;
					if (sc_e_lf_c != MIN_I16)
						sc_e_lf_c -= sc_rde;
					assert_leq(sc_h_lf_c, prob_.cper_->perf_);
					assert_leq(sc_e_lf_c, prob_.cper_->perf_);
				}
				else if (xi > 0)
				{
					if (sc_h_lf_last != MIN_I16)
					{
						sc_h_lf = sc_h_lf_last;
						sc_h_lf_c = sc_h_lf - sc_rdo;
					}
					if (sc_e_lf_last != MIN_I16)
					{
						sc_e_lf = sc_e_lf_last;
						sc_e_lf_c = sc_e_lf - sc_rde;
					}
				}
				if (yi > 0)
				{
					assert(qup != NULL);
					assert(local || qup[j].sc[2] < 0);
					if (qup[j].sc[0] > MIN_I16)
					{
						DEBUG_CHECK(qup[j].sc[0], yi - 1, xi, 2);
						sc_h_up = qup[j].sc[0];
						sc_h_up_c = sc_h_up - sc_rfo;
					}
					if (qup[j].sc[2] > MIN_I16)
					{
						DEBUG_CHECK(qup[j].sc[2], yi - 1, xi, 1);
						sc_f_up = qup[j].sc[2];
						sc_f_up_c = sc_f_up - sc_rfe;
					}
				}
				if (local)
				{
					sc_h_up_c = max<int16_t>(sc_h_up_c, 0);
					sc_f_up_c = max<int16_t>(sc_f_up_c, 0);
					sc_h_lf_c = max<int16_t>(sc_h_lf_c, 0);
					sc_e_lf_c = max<int16_t>(sc_e_lf_c, 0);
				}
			}
			if (sc_h_dg > MIN_I16)
			{
				sc_h_dg_c = sc_h_dg + sc_diag;
			}
			if (local)
				sc_h_dg_c = max<int16_t>(sc_h_dg_c, 0);
			int mask = 0;
			int16_t sc_best = sc_h_dg_c;
			if (sc_h_dg_c > MIN_I64)
			{
				mask = 1;
			}
			if (xi > 0 && sc_h_lf_c >= sc_best && sc_h_lf_c > MIN_I64)
			{
				if (sc_h_lf_c > sc_best)
					mask = 0;
				mask |= 2;
				sc_best = sc_h_lf_c;
			}
			if (xi > 0 && sc_e_lf_c >= sc_best && sc_e_lf_c > MIN_I64)
			{
				if (sc_e_lf_c > sc_best)
					mask = 0;
				mask |= 4;
				sc_best = sc_e_lf_c;
			}
			if (yi > 0 && sc_h_up_c >= sc_best && sc_h_up_c > MIN_I64)
			{
				if (sc_h_up_c > sc_best)
					mask = 0;
				mask |= 8;
				sc_best = sc_h_up_c;
			}
			if (yi > 0 && sc_f_up_c >= sc_best && sc_f_up_c > MIN_I64)
			{
				if (sc_f_up_c > sc_best)
					mask = 0;
				mask |= 16;
				sc_best = sc_f_up_c;
			}
			int16_t sc_e_best = sc_h_lf_c;
			if (xi > 0)
			{
				if (sc_h_lf_c >= sc_e_lf_c && sc_h_lf_c > MIN_I64)
				{
					if (sc_h_lf_c == sc_e_lf_c)
					{
						mask |= 64;
					}
					mask |= 32;
				}
				else if (sc_e_lf_c > MIN_I64)
				{
					sc_e_best = sc_e_lf_c;
					mask |= 64;
				}
			}
			if (sc_e_best > sc_best)
			{
				sc_best = sc_e_best;
				mask &= ~31;
			}
			int16_t sc_f_best = sc_h_up_c;
			if (yi > 0)
			{
				if (sc_h_up_c >= sc_f_up_c && sc_h_up_c > MIN_I64)
				{
					if (sc_h_up_c == sc_f_up_c)
					{
						mask |= 256;
					}
					mask |= 128;
				}
				else if (sc_f_up_c > MIN_I64)
				{
					sc_f_best = sc_f_up_c;
					mask |= 256;
				}
			}
			if (sc_f_best > sc_best)
			{
				sc_best = sc_f_best;
				mask &= ~127;
			}
			assert(local || sc_best <= 0);
			sq_[ii + j].sc[0] = sc_best;
			assert(local || sc_e_best < 0);
			assert(local || sc_f_best < 0);
			assert(!local || sc_e_best >= 0 || sc_e_best == MIN_I16);
			assert(!local || sc_f_best >= 0 || sc_f_best == MIN_I16);
			sq_[ii + j].sc[1] = sc_e_best;
			sq_[ii + j].sc[2] = sc_f_best;
			sq_[ii + j].sc[3] = mask;
			DEBUG_CHECK(sq_[ii + j].sc[0], yi, xi, 2);
			DEBUG_CHECK(sq_[ii + j].sc[1], yi, xi, 0);
			DEBUG_CHECK(sq_[ii + j].sc[2], yi, xi, 1);
			sc_h_lf_last = sc_best;
			sc_e_lf_last = sc_e_best;
		}
		m128mod++;
		if (m128mod == prob_.cper_->niter_)
		{
			m128mod = 0;
			m128div++;
		}
		ii += sq_ncol;
		qup = sq_.ptr() + sq_ncol * i;
	}
	assert_eq(targ, sq_[ymod * sq_ncol + xmod].sc[hef]);
	int64_t rowc = row, colc = col;
	size_t curid;
	int hefc = hef;
	if (bs_.empty())
	{
		CHECK_ROW_COL(rowc, colc);
		curid = bs_.alloc();
		assert_eq(0, curid);
		Edit e;
		bs_[curid].init(
			prob_,
			0,
			0,
			0,
			rowc,
			colc,
			e,
			0,
			true,
			false);
		bs_[curid].len_ = 0;
	}
	else
	{
		curid = bs_.size() - 1;
	}
	size_t ymodTimesNcol = ymod * sq_ncol;
	while (true)
	{
		assert_eq(ymodTimesNcol, ymod * sq_ncol);
		CpQuad *cur = sq_.ptr() + ymodTimesNcol + xmod;
		int mask = cur->sc[3];
		assert_gt(mask, 0);
		int sel = -1;
		if (hefc == 0)
		{
			if ((mask & 1) != 0)
			{
				sel = 0;
			}
			else if ((mask & 8) != 0)
			{
				sel = 3;
			}
			else if ((mask & 16) != 0)
			{
				sel = 4;
			}
			else if ((mask & 2) != 0)
			{
				sel = 1;
			}
			else if ((mask & 4) != 0)
			{
				sel = 2;
			}
		}
		else if (hefc == 1)
		{
			if ((mask & 32) != 0)
			{
				sel = 5;
			}
			else if ((mask & 64) != 0)
			{
				sel = 6;
			}
		}
		else
		{
			assert_eq(2, hefc);
			if ((mask & 128) != 0)
			{
				sel = 7;
			}
			else if ((mask & 256) != 0)
			{
				sel = 8;
			}
		}
		assert_geq(sel, 0);
		int qc = prob_.qry_[rowc], qq = prob_.qual_[rowc];
		int rc = prob_.ref_[colc];
		assert_range(0, 16, rc);
		bool xexit = false, yexit = false;
		if (sel == 0)
		{
			assert_geq(rowc, 0);
			assert_geq(colc, 0);
			TAlScore scd = prob_.sc_->score(qc, rc, qq - 33);
			if ((rc & (1 << qc)) == 0)
			{
				size_t id = curid;
				bool empty = (bs_[curid].len_ == 0 && curid == 0);
				if (!empty)
				{
					id = bs_.alloc();
				}
				Edit e((int)rowc, mask2dna[rc], "ACGTN"[qc], EDIT_TYPE_MM);
				assert_lt(scd, 0);
				TAlScore score_en = bs_[curid].score_st_ + scd;
				bs_[id].init(
					prob_,
					curid,
					-scd, score_en, rowc, colc, e, hefc, empty, false);
				curid = id;
			}
			else
			{
				bs_[curid].score_st_ += prob_.sc_->match();
				bs_[curid].len_++;
				assert_leq((int64_t)bs_[curid].len_, bs_[curid].row_ + 1);
			}
			if (xmod == 0)
				xexit = true;
			if (ymod == 0)
				yexit = true;
			rowc--;
			ymod--;
			ymodTimesNcol -= sq_ncol;
			colc--;
			xmod--;
			assert(local || bs_[curid].score_st_ >= targ_final);
			hefc = 0;
		}
		else if ((sel >= 1 && sel <= 2) || (sel >= 5 && sel <= 6))
		{
			assert_gt(colc, 0);
			size_t id = bs_.alloc();
			Edit e((int)rowc + 1, mask2dna[rc], '-', EDIT_TYPE_READ_GAP);
			TAlScore gapp = prob_.sc_->readGapOpen();
			if (bs_[curid].len_ == 0 && bs_[curid].e_.inited() && bs_[curid].e_.isReadGap())
			{
				gapp = prob_.sc_->readGapExtend();
			}
			TAlScore score_en = bs_[curid].score_st_ - gapp;
			bs_[id].init(
				prob_,
				curid,
				gapp, score_en, rowc, colc - 1, e, hefc, false, false);
			if (xmod == 0)
				xexit = true;
			colc--;
			xmod--;
			curid = id;
			assert(local || bs_[curid].score_st_ >= targ_final);
			if (sel == 1 || sel == 5)
			{
				hefc = 0;
			}
			else
			{
				hefc = 1;
			}
		}
		else
		{
			assert_gt(rowc, 0);
			size_t id = bs_.alloc();
			Edit e((int)rowc, '-', "ACGTN"[qc], EDIT_TYPE_REF_GAP);
			TAlScore gapp = prob_.sc_->refGapOpen();
			if (bs_[curid].len_ == 0 && bs_[curid].e_.inited() && bs_[curid].e_.isRefGap())
			{
				gapp = prob_.sc_->refGapExtend();
			}
			TAlScore score_en = bs_[curid].score_st_ - gapp;
			bs_[id].init(
				prob_,
				curid,
				gapp, score_en, rowc - 1, colc, e, hefc, false, false);
			if (ymod == 0)
				yexit = true;
			rowc--;
			ymod--;
			ymodTimesNcol -= sq_ncol;
			curid = id;
			assert(local || bs_[curid].score_st_ >= targ_final);
			if (sel == 3 || sel == 7)
			{
				hefc = 0;
			}
			else
			{
				hefc = 2;
			}
		}
		CHECK_ROW_COL(rowc, colc);
		CpQuad *cur_new = NULL;
		if (!xexit && !yexit)
		{
			cur_new = sq_.ptr() + ymodTimesNcol + xmod;
		}
		if (colc < 0 || rowc < 0 ||
			(cur_new != NULL && local && cur_new->sc[0] == 0))
		{
			done = true;
			assert(bs_[curid].isSolution(prob_));
			addSolution(curid);
#ifndef NDEBUG
			size_t cur = curid;
			if (!bs_[cur].root_)
			{
				size_t next = bs_[cur].parentId_;
				while (!bs_[next].root_)
				{
					assert_neq(cur, next);
					if (bs_[next].len_ != 0 || bs_[cur].len_ == 0)
					{
						assert(!bs_[cur].overlap(prob_, bs_[next]));
					}
					cur = next;
					next = bs_[cur].parentId_;
				}
			}
#endif
			return;
		}
		assert(!xexit || hefc == 0 || hefc == 1);
		assert(!yexit || hefc == 0 || hefc == 2);
		if (xexit || yexit)
		{
			row_new = rowc;
			col_new = colc;
			hef_new = hefc;
			done = false;
			if (rowc < 0 || colc < 0)
			{
				assert(local);
				targ_new = 0;
			}
			else
			{
				targ_new = prob_.cper_->scoreSquare(rowc, colc, hefc);
				assert(local || targ_new >= targ);
				assert(local || targ_new >= targ_final);
			}
			if (local && targ_new == 0)
			{
				assert_eq(0, hefc);
				done = true;
				assert(bs_[curid].isSolution(prob_));
				addSolution(curid);
			}
			assert((row_new >= 0 && col_new >= 0) || done);
			return;
		}
	}
	assert(false);
}
void BtBranch::init(
	const BtBranchProblem &prob,
	size_t parentId,
	TAlScore penalty,
	TAlScore score_en,
	int64_t row,
	int64_t col,
	Edit e,
	int hef,
	bool root,
	bool extend)
{
	score_en_ = score_en;
	penalty_ = penalty;
	score_st_ = score_en_;
	row_ = row;
	col_ = col;
	parentId_ = parentId;
	e_ = e;
	root_ = root;
	assert(!root_ || parentId == 0);
	assert_lt(row, (int64_t)prob.qrylen_);
	assert_lt(col, (int64_t)prob.reflen_);
	int64_t rowc = row;
	int64_t colc = col;
	len_ = 0;
	if (e.inited() && e.isMismatch())
	{
		rowc--;
		colc--;
		len_ = 1;
	}
	int64_t match = prob.sc_->match();
	bool cp = prob.usecp_;
	size_t iters = 0;
	curtailed_ = false;
	if (extend)
	{
		while (rowc >= 0 && colc >= 0)
		{
			int rfm = prob.ref_[colc];
			assert_range(0, 16, rfm);
			int rdc = prob.qry_[rowc];
			bool matches = (rfm & (1 << rdc)) != 0;
			if (!matches)
			{
				break;
			}
			score_st_ += match;
			if (cp && rowc - 1 >= 0 && colc - 1 >= 0 &&
				prob.cper_->isCheckpointed(rowc - 1, colc - 1))
			{
				int16_t cpsc;
				cpsc = prob.cper_->scoreTriangle(rowc - 1, colc - 1, hef);
				if (cpsc + score_st_ < prob.targ_)
				{
					curtailed_ = true;
					break;
				}
			}
			iters++;
			rowc--;
			colc--;
		}
	}
	assert_geq(rowc, -1);
	assert_geq(colc, -1);
	len_ = (int64_t)row - rowc;
	assert_leq((int64_t)len_, row_ + 1);
	assert_leq((int64_t)len_, col_ + 1);
	assert_leq((int64_t)score_st_, (int64_t)prob.qrylen_ * match);
}
void BtBranchTracer::examineBranch(
	int64_t row,
	int64_t col,
	const Edit &e,
	TAlScore pen,
	TAlScore sc,
	size_t parentId)
{
	size_t id = bs_.alloc();
	bs_[id].init(prob_, parentId, pen, sc, row, col, e, 0, false, true);
	if (bs_[id].isSolution(prob_))
	{
		assert(bs_[id].isValid(prob_));
		addSolution(id);
	}
	else
	{
		if (bs_[id].isValid(prob_))
		{
			add(id);
		}
		else
		{
			bs_.pop();
		}
	}
}
void BtBranchTracer::addOffshoots(size_t bid)
{
	BtBranch &b = bs_[bid];
	TAlScore sc = b.score_en_;
	int64_t match = prob_.sc_->match();
	int64_t scoreFloor = prob_.sc_->monotone ? MIN_I64 : 0;
	bool cp = prob_.usecp_;
	ASSERT_ONLY(TAlScore perfectScore = prob_.sc_->perfectScore(prob_.qrylen_));
	assert_leq(prob_.targ_, perfectScore);
	for (size_t i = 0; i < b.len_; i++)
	{
		assert_leq((int64_t)i, b.row_ + 1);
		assert_leq((int64_t)i, b.col_ + 1);
		int64_t row = b.row_ - i, col = b.col_ - i;
		int64_t bonusLeft = (row + 1) * match;
		int64_t fromend = prob_.qrylen_ - row - 1;
		bool allowGaps = fromend >= prob_.sc_->gapbar && row >= prob_.sc_->gapbar;
		if (allowGaps && row >= 0 && col >= 0)
		{
			if (col > 0)
			{
				bool extend = b.e_.inited() && b.e_.isReadGap() && i == 0;
				TAlScore rdgapPen = extend ? prob_.sc_->readGapExtend() : prob_.sc_->readGapOpen();
				bool prune = false;
				assert_gt(rdgapPen, 0);
				if (cp && prob_.cper_->isCheckpointed(row, col - 1))
				{
					int16_t cpsc = (int16_t)prob_.cper_->scoreTriangle(row, col - 1, 0);
					assert_leq(cpsc, perfectScore);
					assert_geq(prob_.sc_->readGapOpen(), prob_.sc_->readGapExtend());
					TAlScore bonus = prob_.sc_->readGapOpen() - prob_.sc_->readGapExtend();
					assert_geq(bonus, 0);
					if (cpsc + bonus + sc - rdgapPen < prob_.targ_)
					{
						prune = true;
					}
				}
				if (prune)
				{
					if (extend)
					{
						nrdexPrune_++;
					}
					else
					{
						nrdopPrune_++;
					}
				}
				else if (sc - rdgapPen >= scoreFloor && sc - rdgapPen + bonusLeft >= prob_.targ_)
				{
					Edit e((int)row + 1, mask2dna[(int)prob_.ref_[col]], '-', EDIT_TYPE_READ_GAP);
					assert(e.isReadGap());
					examineBranch(row, col - 1, e, rdgapPen, sc - rdgapPen, bid);
					if (extend)
					{
						nrdex_++;
					}
					else
					{
						nrdop_++;
					}
				}
			}
			if (row > 0)
			{
				bool extend = b.e_.inited() && b.e_.isRefGap() && i == 0;
				TAlScore rfgapPen = (b.e_.inited() && b.e_.isRefGap()) ? prob_.sc_->refGapExtend() : prob_.sc_->refGapOpen();
				bool prune = false;
				assert_gt(rfgapPen, 0);
				if (cp && prob_.cper_->isCheckpointed(row - 1, col))
				{
					int16_t cpsc = (int16_t)prob_.cper_->scoreTriangle(row - 1, col, 0);
					assert_leq(cpsc, perfectScore);
					assert_geq(prob_.sc_->refGapOpen(), prob_.sc_->refGapExtend());
					TAlScore bonus = prob_.sc_->refGapOpen() - prob_.sc_->refGapExtend();
					assert_geq(bonus, 0);
					if (cpsc + bonus + sc - rfgapPen < prob_.targ_)
					{
						prune = true;
					}
				}
				if (prune)
				{
					if (extend)
					{
						nrfexPrune_++;
					}
					else
					{
						nrfopPrune_++;
					}
				}
				else if (sc - rfgapPen >= scoreFloor && sc - rfgapPen + bonusLeft >= prob_.targ_)
				{
					Edit e((int)row, '-', "ACGTN"[(int)prob_.qry_[row]], EDIT_TYPE_REF_GAP);
					assert(e.isRefGap());
					examineBranch(row - 1, col, e, rfgapPen, sc - rfgapPen, bid);
					if (extend)
					{
						nrfex_++;
					}
					else
					{
						nrfop_++;
					}
				}
			}
		}
		if (i == b.len_ && !b.curtailed_ && row >= 0 && col >= 0)
		{
			int rfm = prob_.ref_[col];
			assert_lt(row, (int64_t)prob_.qrylen_);
			int rdc = prob_.qry_[row];
			int rdq = prob_.qual_[row];
			int scdiff = prob_.sc_->score(rdc, rfm, rdq - 33);
			assert_lt(scdiff, 0);
			bool prune = false;
			if (cp && row > 0 && col > 0 && prob_.cper_->isCheckpointed(row - 1, col - 1))
			{
				int16_t cpsc = prob_.cper_->scoreTriangle(row - 1, col - 1, 0);
				assert_leq(cpsc, perfectScore);
				assert_leq(cpsc + scdiff + sc, perfectScore);
				if (cpsc + scdiff + sc < prob_.targ_)
				{
					prune = true;
				}
			}
			if (prune)
			{
				nmm_++;
			}
			else
			{
				if (sc + scdiff >= scoreFloor && sc + scdiff + bonusLeft >= prob_.targ_)
				{
					Edit e((int)row, mask2dna[rfm], "ACGTN"[rdc], EDIT_TYPE_MM);
					bool nmm = (mask2dna[rfm] == 'N' || rdc > 4);
					assert_neq(e.chr, e.qchr);
					assert_lt(scdiff, 0);
					examineBranch(row - 1, col - 1, e, -scdiff, sc + scdiff, bid);
					if (nmm)
					{
						nnmm_++;
					}
					else
					{
						nmm_++;
					}
				}
			}
		}
		sc += match;
	}
}
void BtBranchTracer::flushUnsorted()
{
	if (unsorted_.empty())
	{
		return;
	}
	unsorted_.sort();
	unsorted_.reverse();
#ifndef NDEBUG
	for (size_t i = 1; i < unsorted_.size(); i++)
	{
		assert_leq(bs_[unsorted_[i].second].score_st_, bs_[unsorted_[i - 1].second].score_st_);
	}
#endif
	EList<size_t> *src2 = sortedSel_ ? &sorted1_ : &sorted2_;
	EList<size_t> *dest = sortedSel_ ? &sorted2_ : &sorted1_;
	dest->clear();
	size_t cur1 = 0, cur2 = cur_;
	while (cur1 < unsorted_.size() || cur2 < src2->size())
	{
		bool take1 = true;
		if (cur1 == unsorted_.size())
		{
			take1 = false;
		}
		else if (cur2 == src2->size())
		{
			take1 = true;
		}
		else
		{
			assert_neq(unsorted_[cur1].second, (*src2)[cur2]);
			take1 = bs_[unsorted_[cur1].second] < bs_[(*src2)[cur2]];
		}
		if (take1)
		{
			dest->push_back(unsorted_[cur1++].second);
		}
		else
		{
			dest->push_back((*src2)[cur2++]);
		}
	}
	assert_eq(cur1, unsorted_.size());
	assert_eq(cur2, src2->size());
	sortedSel_ = !sortedSel_;
	cur_ = 0;
	unsorted_.clear();
}
bool BtBranchTracer::trySolutions(
	bool lookForOlap,
	SwResult &res,
	size_t &off,
	size_t &nrej,
	RandomSource &rnd,
	bool &success)
{
	if (solutions_.size() > 0)
	{
		for (size_t i = 0; i < solutions_.size(); i++)
		{
			int ret = trySolution(solutions_[i], lookForOlap, res, off, nrej, rnd);
			if (ret == BT_FOUND)
			{
				success = true;
				return true;
			}
		}
		solutions_.clear();
		success = false;
		return true;
	}
	return false;
}
int BtBranchTracer::trySolution(
	size_t id,
	bool lookForOlap,
	SwResult &res,
	size_t &off,
	size_t &nrej,
	RandomSource &rnd)
{
	AlnScore score;
	BtBranch *br = &bs_[id];
	EList<Edit> &ned = res.alres.ned();
	const BtBranch *cur = br, *prev = NULL;
	size_t ns = 0, nrefns = 0;
	size_t ngap = 0;
	while (true)
	{
		if (cur->e_.inited())
		{
			if (cur->e_.isMismatch())
			{
				if (cur->e_.qchr == 'N' || cur->e_.chr == 'N')
				{
					ns++;
				}
			}
			else if (cur->e_.isGap())
			{
				ngap++;
			}
			if (cur->e_.chr == 'N')
			{
				nrefns++;
			}
			ned.push_back(cur->e_);
		}
		if (cur->root_)
		{
			break;
		}
		cur = &bs_[cur->parentId_];
	}
	if (ns > prob_.nceil_)
	{
		res.reset();
		assert(res.alres.ned().empty());
		nrej++;
		return BT_REJECTED_N;
	}
	cur = br;
	bool rejSeen = false;
	bool rejCore = true;
	while (true)
	{
		int64_t row = cur->row_, col = cur->col_;
		assert_lt(row, (int64_t)prob_.qrylen_);
		size_t fromend = prob_.qrylen_ - row - 1;
		size_t diag = fromend + col;
		int64_t diagi = col - row;
		diagi += prob_.rect_->triml;
		assert_lt(diag, seenPaths_.size());
		if (diagi >= 0)
		{
			size_t diag = (size_t)diagi;
			if (diag >= prob_.rect_->corel &&
				diag <= prob_.rect_->corer)
			{
				rejCore = false;
			}
		}
		if (lookForOlap)
		{
			int64_t newlo, newhi;
			if (cur->len_ == 0)
			{
				if (prev != NULL && prev->len_ > 0)
				{
					newhi = newlo = 0;
				}
				else
				{
					newlo = row;
					newhi = row + 1;
				}
			}
			else
			{
				newlo = row - (cur->len_ - 1);
				newhi = row + 1;
			}
			assert_geq(newlo, 0);
			assert_geq(newhi, 0);
			if (newhi > newlo)
			{
				bool added = false;
				const size_t sz = seenPaths_[diag].size();
				for (size_t i = 0; i < sz; i++)
				{
					size_t lo = seenPaths_[diag][i].first;
					size_t hi = seenPaths_[diag][i].second;
					assert_lt(lo, hi);
					size_t lo_sm = newlo, hi_sm = newhi;
					if (hi - lo < hi_sm - lo_sm)
					{
						swap(lo, lo_sm);
						swap(hi, hi_sm);
					}
					if ((lo <= lo_sm && hi > lo_sm) ||
						(lo < hi_sm && hi >= hi_sm))
					{
						rejSeen = true;
						seenPaths_[diag][i].first = min(lo, lo_sm);
						seenPaths_[diag][i].second = max(hi, hi_sm);
#ifndef NDEBUG
						for (int64_t ii = seenPaths_[diag][i].first;
							 ii < (int64_t)seenPaths_[diag][i].second;
							 ii++)
						{
						}
#endif
						added = true;
						break;
					}
					else if (hi == lo_sm || lo == hi_sm)
					{
						seenPaths_[diag][i].first = min(lo, lo_sm);
						seenPaths_[diag][i].second = max(hi, hi_sm);
#ifndef NDEBUG
						for (int64_t ii = seenPaths_[diag][i].first;
							 ii < (int64_t)seenPaths_[diag][i].second;
							 ii++)
						{
						}
#endif
						added = true;
					}
				}
				if (!added)
				{
					seenPaths_[diag].push_back(make_pair(newlo, newhi));
				}
			}
		}
		if (cur->root_)
		{
			assert_eq(0, cur->parentId_);
			break;
		}
		prev = cur;
		cur = &bs_[cur->parentId_];
	}
	if (rejSeen)
	{
		res.reset();
		assert(res.alres.ned().empty());
		nrej++;
		return BT_NOT_FOUND;
	}
	if (rejCore)
	{
		res.reset();
		assert(res.alres.ned().empty());
		nrej++;
		return BT_REJECTED_CORE_DIAG;
	}
	off = br->leftmostCol();
	size_t trimBeg = br->uppermostRow();
	size_t trimEnd = prob_.qrylen_ - prob_.row_ - 1;
	score.score_ = prob_.targ_;
	score.basesAligned_ = (int)(prob_.qrylen_ - trimBeg - trimEnd - ned.size());
	score.edits_ = (int)ned.size();
	score.ns_ = ns;
	score.gaps_ = ngap;
	res.alres.setScore(score);
	res.alres.setRefNs(nrefns);
	assert_leq(trimBeg, prob_.qrylen_);
	assert_leq(trimEnd, prob_.qrylen_);
	TRefOff refoff = off + prob_.refoff_ + prob_.rect_->refl;
	res.alres.setShape(
		prob_.refid_,
		refoff, prob_.treflen(), prob_.fw_, prob_.qrylen_, true, 0, 0, true, prob_.fw_ ? trimBeg : trimEnd, prob_.fw_ ? trimEnd : trimBeg);
	return BT_FOUND;
}
bool BtBranchTracer::nextAlignmentBacktrace(
	size_t maxiter,
	SwResult &res,
	size_t &off,
	size_t &nrej,
	size_t &niter,
	RandomSource &rnd)
{
	assert(!empty() || !emptySolution());
	assert(prob_.inited());
	bool result = false;
	niter = 0;
	while (!empty())
	{
		if (trySolutions(true, res, off, nrej, rnd, result))
		{
			return result;
		}
		if (niter++ >= maxiter)
		{
			break;
		}
		size_t brid = best(rnd);
		assert(!seen_.contains(brid));
		ASSERT_ONLY(seen_.insert(brid));
#if 0
		BtBranch *br = &bs_[brid];
		cerr << brid
		     << ": targ:" << prob_.targ_
		     << ", sc:" << br->score_st_
		     << ", row:" << br->uppermostRow()
			 << ", nmm:" << nmm_
			 << ", nnmm:" << nnmm_
			 << ", nrdop:" << nrdop_
			 << ", nrfop:" << nrfop_
			 << ", nrdex:" << nrdex_
			 << ", nrfex:" << nrfex_
			 << ", nrdop_pr: " << nrdopPrune_
			 << ", nrfop_pr: " << nrfopPrune_
			 << ", nrdex_pr: " << nrdexPrune_
			 << ", nrfex_pr: " << nrfexPrune_
			 << endl;
#endif
		addOffshoots(brid);
	}
	if (trySolutions(true, res, off, nrej, rnd, result))
	{
		return result;
	}
	return false;
}
bool BtBranchTracer::nextAlignmentFill(
	size_t maxiter,
	SwResult &res,
	size_t &off,
	size_t &nrej,
	size_t &niter,
	RandomSource &rnd)
{
	assert(prob_.inited());
	assert(!emptySolution());
	bool result = false;
	if (trySolutions(false, res, off, nrej, rnd, result))
	{
		return result;
	}
	return false;
}
bool BtBranchTracer::nextAlignment(
	size_t maxiter,
	SwResult &res,
	size_t &off,
	size_t &nrej,
	size_t &niter,
	RandomSource &rnd)
{
	if (prob_.fill_)
	{
		return nextAlignmentFill(
			maxiter,
			res,
			off,
			nrej,
			niter,
			rnd);
	}
	else
	{
		return nextAlignmentBacktrace(
			maxiter,
			res,
			off,
			nrej,
			niter,
			rnd);
	}
}
#ifdef MAIN_ALIGNER_BT
#include <iostream>
int main(int argc, char **argv)
{
	size_t off = 0;
	RandomSource rnd(77);
	BtBranchTracer tr;
	Scoring sc = Scoring::base1();
	SwResult res;
	tr.init(
		"ACGTACGT",
		"IIIIIIII", 8, "ACGTACGT", 8, 0, 0, true, sc, 0, 8, 7, 7, rnd);
	size_t nrej = 0;
	tr.nextAlignment(
		res,
		off,
		nrej,
		rnd);
}
#endif
