#ifndef RANDOM_UTIL_H_
#define RANDOM_UTIL_H_
#include <algorithm>
#include "random_source.h"
#include "ds.h"
class Random1toN
{
	typedef size_t T;
public:
	static const size_t SWAPLIST_THRESH;
	static const size_t CONVERSION_THRESH;
	static const float CONVERSION_FRAC;
	Random1toN(int cat = 0) : sz_(0), n_(0), cur_(0),
							  list_(SWAPLIST_THRESH, cat), seen_(CONVERSION_THRESH, cat),
							  thresh_(0) {}
	Random1toN(size_t n, int cat = 0) : sz_(0), n_(n), cur_(0),
										list_(SWAPLIST_THRESH, cat), seen_(CONVERSION_THRESH, cat),
										thresh_(0) {}
	void init(size_t n, bool withoutReplacement)
	{
		sz_ = n_ = n;
		converted_ = false;
		swaplist_ = n < SWAPLIST_THRESH || withoutReplacement;
		cur_ = 0;
		list_.clear();
		seen_.clear();
		thresh_ = std::max(CONVERSION_THRESH, (size_t)(CONVERSION_FRAC * n));
	}
	void reset()
	{
		sz_ = n_ = cur_ = 0;
		swaplist_ = converted_ = false;
		list_.clear();
		seen_.clear();
		thresh_ = 0;
	}
	T next(RandomSource &rnd)
	{
		assert(!done());
		if (cur_ == 0 && !converted_)
		{
			if (n_ == 1)
			{
				cur_ = 1;
				return 0;
			}
			if (swaplist_)
			{
				list_.resize(n_);
				for (size_t i = 0; i < n_; i++)
				{
					list_[i] = (T)i;
				}
			}
		}
		if (swaplist_)
		{
			size_t r = cur_ + (rnd.nextU32() % (n_ - cur_));
			if (r != cur_)
			{
				std::swap(list_[cur_], list_[r]);
			}
			return list_[cur_++];
		}
		else
		{
			assert(!converted_);
			bool again = true;
			T rn = 0;
			size_t seenSz = seen_.size();
			while (again)
			{
				rn = rnd.nextU32() % (T)n_;
				again = false;
				for (size_t i = 0; i < seenSz; i++)
				{
					if (seen_[i] == rn)
					{
						again = true;
						break;
					}
				}
			}
			seen_.push_back(rn);
			cur_++;
			assert_leq(cur_, n_);
			assert_gt(thresh_, 0);
			if (seen_.size() >= thresh_ && cur_ < n_)
			{
				assert(!seen_.empty());
				seen_.sort();
				list_.resize(n_ - cur_);
				size_t prev = 0;
				size_t cur = 0;
				for (size_t i = 0; i <= seenSz; i++)
				{
					for (size_t j = prev; j < seen_[i]; j++)
					{
						list_[cur++] = (T)j;
					}
					prev = seen_[i] + 1;
				}
				for (size_t j = prev; j < n_; j++)
				{
					list_[cur++] = (T)j;
				}
				assert_eq(cur, n_ - cur_);
				seen_.clear();
				cur_ = 0;
				n_ = list_.size();
				converted_ = true;
				swaplist_ = true;
			}
			return rn;
		}
	}
	bool inited() const { return n_ > 0; }
	void setDone()
	{
		assert(inited());
		cur_ = n_;
	}
	bool done() const { return inited() && cur_ >= n_; }
	size_t size() const { return n_; }
	size_t left() const { return n_ - cur_; }
	size_t totalSizeBytes() const
	{
		return list_.totalSizeBytes() +
			   seen_.totalSizeBytes();
	}
	size_t totalCapacityBytes() const
	{
		return list_.totalCapacityBytes() +
			   seen_.totalCapacityBytes();
	}
protected:
	size_t sz_;
	size_t n_;
	bool swaplist_;
	bool converted_;
	size_t cur_;
	EList<T> list_;
	EList<T> seen_;
	size_t thresh_;
};
#endif
