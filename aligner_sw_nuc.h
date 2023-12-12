#ifndef ALIGNER_SW_NUC_H_
#define ALIGNER_SW_NUC_H_
#include <stdint.h>
#include "aligner_sw_common.h"
#include "aligner_result.h"
struct DpNucFrame
{
	void init(
		size_t nedsz_,
		size_t aedsz_,
		size_t celsz_,
		size_t row_,
		size_t col_,
		size_t gaps_,
		size_t readGaps_,
		size_t refGaps_,
		AlnScore score_,
		int ct_)
	{
		nedsz = nedsz_;
		aedsz = aedsz_;
		celsz = celsz_;
		row = row_;
		col = col_;
		gaps = gaps_;
		readGaps = readGaps_;
		refGaps = refGaps_;
		score = score_;
		ct = ct_;
	}
	size_t nedsz;
	size_t aedsz;
	size_t celsz;
	size_t row;
	size_t col;
	size_t gaps;
	size_t readGaps;
	size_t refGaps;
	AlnScore score;
	int ct;
};
enum
{
	BT_CAND_FATE_SUCCEEDED = 1,
	BT_CAND_FATE_FAILED,
	BT_CAND_FATE_FILT_START,
	BT_CAND_FATE_FILT_DOMINATED,
	BT_CAND_FATE_FILT_SCORE
};
struct DpBtCandidate
{
	DpBtCandidate() { reset(); }
	DpBtCandidate(size_t row_, size_t col_, TAlScore score_)
	{
		init(row_, col_, score_);
	}
	void reset() { init(0, 0, 0); }
	void init(size_t row_, size_t col_, TAlScore score_)
	{
		row = row_;
		col = col_;
		score = score_;
		fate = 0;
	}
	inline bool dominatedBy(const DpBtCandidate &o)
	{
		const size_t SQ = 40;
		size_t rowhi = row;
		size_t rowlo = o.row;
		if (rowhi < rowlo)
			swap(rowhi, rowlo);
		size_t colhi = col;
		size_t collo = o.col;
		if (colhi < collo)
			swap(colhi, collo);
		return (colhi - collo) <= SQ &&
			   (rowhi - rowlo) <= SQ;
	}
	bool operator>(const DpBtCandidate &o) const
	{
		if (score < o.score)
			return true;
		if (score > o.score)
			return false;
		if (row < o.row)
			return true;
		if (row > o.row)
			return false;
		if (col < o.col)
			return true;
		if (col > o.col)
			return false;
		return false;
	}
	bool operator<(const DpBtCandidate &o) const
	{
		if (score > o.score)
			return true;
		if (score < o.score)
			return false;
		if (row > o.row)
			return true;
		if (row < o.row)
			return false;
		if (col > o.col)
			return true;
		if (col < o.col)
			return false;
		return false;
	}
	bool operator==(const DpBtCandidate &o) const
	{
		return row == o.row &&
			   col == o.col &&
			   score == o.score;
	}
	bool operator>=(const DpBtCandidate &o) const { return !((*this) < o); }
	bool operator<=(const DpBtCandidate &o) const { return !((*this) > o); }
#ifndef NDEBUG
	bool repOk() const
	{
		assert(VALID_SCORE(score));
		return true;
	}
#endif
	size_t row;
	size_t col;
	TAlScore score;
	int fate;
};
template <typename T>
class NBest
{
public:
	NBest<T>() { nelt_ = nbest_ = n_ = 0; }
	bool inited() const { return nelt_ > 0; }
	void init(size_t nelt, size_t nbest)
	{
		nelt_ = nelt;
		nbest_ = nbest;
		elts_.resize(nelt * nbest);
		ncur_.resize(nelt);
		ncur_.fill(0);
		n_ = 0;
	}
	bool add(size_t elt, const T &o)
	{
		assert_lt(elt, nelt_);
		const size_t ncur = ncur_[elt];
		assert_leq(ncur, nbest_);
		n_++;
		for (size_t i = 0; i < nbest_ && i <= ncur; i++)
		{
			if (o > elts_[nbest_ * elt + i] || i >= ncur)
			{
				for (int j = (int)ncur; j > (int)i; j--)
				{
					if (j < (int)nbest_)
					{
						elts_[nbest_ * elt + j] = elts_[nbest_ * elt + j - 1];
					}
				}
				elts_[nbest_ * elt + i] = o;
				if (ncur < nbest_)
				{
					ncur_[elt]++;
				}
				return true;
			}
		}
		return false;
	}
	bool empty() const
	{
		return n_ == 0;
	}
	template <typename TList>
	void dump(TList &l) const
	{
		if (empty())
			return;
		for (size_t i = 0; i < nelt_; i++)
		{
			assert_leq(ncur_[i], nbest_);
			for (size_t j = 0; j < ncur_[i]; j++)
			{
				l.push_back(elts_[i * nbest_ + j]);
			}
		}
	}
protected:
	size_t nelt_;
	size_t nbest_;
	EList<T> elts_;
	EList<size_t> ncur_;
	size_t n_;
};
#endif
