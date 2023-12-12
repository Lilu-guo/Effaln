#ifndef IVAL_LIST_H_
#define IVAL_LIST_H_
#include "ds.h"
#include "ref_coord.h"
#include <utility>
class EIvalMergeList
{
public:
	static const size_t DEFAULT_UNSORTED_SZ = 16;
	explicit EIvalMergeList(int cat = 0) : sorted_(cat),
										   sortedLhs_(cat),
										   unsorted_(cat),
										   unsortedSz_(DEFAULT_UNSORTED_SZ)
	{
	}
	explicit EIvalMergeList(size_t unsortedSz, int cat = 0) : sorted_(cat),
															  sortedLhs_(cat),
															  unsorted_(cat),
															  unsortedSz_(unsortedSz)
	{
	}
	void setUnsortedSize(size_t usz)
	{
		unsortedSz_ = usz;
	}
	void add(const Interval &i)
	{
		assert_leq(unsorted_.size(), unsortedSz_);
		if (unsorted_.size() < unsortedSz_)
		{
			unsorted_.push_back(i);
		}
		if (unsorted_.size() == unsortedSz_)
		{
			flush();
		}
	}
	void flush()
	{
		for (size_t i = 0; i < unsorted_.size(); i++)
		{
			sorted_.push_back(unsorted_[i]);
		}
		sorted_.sort();
		merge();
		sortedLhs_.clear();
		for (size_t i = 0; i < sorted_.size(); i++)
		{
			sortedLhs_.push_back(sorted_[i].upstream());
		}
		assert(sortedLhs_.sorted());
		unsorted_.clear();
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_eq(sorted_.size(), sortedLhs_.size());
		return true;
	}
#endif
	void reset()
	{
		clear();
	}
	void clear()
	{
		sorted_.clear();
		sortedLhs_.clear();
		unsorted_.clear();
	}
	bool locusPresent(const Coord &loc) const
	{
		return locusPresentUnsorted(loc) || locusPresentSorted(loc);
	}
	size_t size() const
	{
		return sorted_.size() + unsorted_.size();
	}
	bool empty() const
	{
		return sorted_.empty() && unsorted_.empty();
	}
protected:
	void merge()
	{
		size_t nmerged = 0;
		for (size_t i = 1; i < sorted_.size(); i++)
		{
			if (sorted_[i - 1].downstream() >= sorted_[i].upstream())
			{
				nmerged++;
				assert_leq(sorted_[i - 1].upstream(), sorted_[i].upstream());
				Coord up = std::min(sorted_[i - 1].upstream(), sorted_[i].upstream());
				Coord dn = std::max(sorted_[i - 1].downstream(), sorted_[i].downstream());
				sorted_[i].setUpstream(up);
				sorted_[i].setLength(dn.off() - up.off());
				sorted_[i - 1].reset();
			}
		}
		sorted_.sort();
		assert_lt(nmerged, sorted_.size());
		sorted_.resize(sorted_.size() - nmerged);
#ifndef NDEBUG
		for (size_t i = 0; i < sorted_.size(); i++)
		{
			assert(sorted_[i].inited());
		}
#endif
	}
	bool locusPresentSorted(const Coord &loc) const
	{
		assert(repOk());
		if (sorted_.empty())
		{
			return false;
		}
		size_t beg = sortedLhs_.bsearchLoBound(loc);
		if (beg == sortedLhs_.size() || sortedLhs_[beg] > loc)
		{
			if (beg == 0)
			{
				return false;
			}
			return sorted_[beg - 1].contains(loc);
		}
		else
		{
			assert_eq(loc, sortedLhs_[beg]);
			return true;
		}
	}
	bool locusPresentUnsorted(const Coord &loc) const
	{
		for (size_t i = 0; i < unsorted_.size(); i++)
		{
			if (unsorted_[i].contains(loc))
			{
				return true;
			}
		}
		return false;
	}
	EList<Interval> sorted_;
	EList<Coord> sortedLhs_;
	EList<Interval> unsorted_;
	size_t unsortedSz_;
};
class EIvalMergeListBinned
{
public:
	static const size_t NBIN = 7;
	explicit EIvalMergeListBinned(int cat = 0) : bins_(1 << NBIN, cat)
	{
		bins_.resize(1 << NBIN);
	}
	explicit EIvalMergeListBinned(
		size_t unsortedSz,
		int cat = 0) : bins_(1 << NBIN, cat)
	{
		bins_.resize(1 << NBIN);
		for (size_t i = 0; i < (1 << NBIN); i++)
		{
			bins_[i].setUnsortedSize(unsortedSz);
		}
	}
	void add(const Interval &i)
	{
		size_t bin = i.ref() & ~(0xffffffff << NBIN);
		assert_lt(bin, bins_.size());
		bins_[bin].add(i);
	}
#ifndef NDEBUG
	bool repOk() const
	{
		for (size_t i = 0; i < bins_.size(); i++)
		{
			assert(bins_[i].repOk());
		}
		return true;
	}
#endif
	void reset()
	{
		clear();
	}
	void clear()
	{
		for (size_t i = 0; i < bins_.size(); i++)
		{
			bins_[i].clear();
		}
	}
	bool locusPresent(const Coord &loc) const
	{
		size_t bin = loc.ref() & ~(0xffffffff << NBIN);
		assert_lt(bin, bins_.size());
		return bins_[bin].locusPresent(loc);
	}
	size_t size() const
	{
		size_t sz = 0;
		for (size_t i = 0; i < bins_.size(); i++)
		{
			sz += bins_[i].size();
		}
		return sz;
	}
	bool empty() const
	{
		return size() == 0;
	}
protected:
	EList<EIvalMergeList> bins_;
};
#endif
