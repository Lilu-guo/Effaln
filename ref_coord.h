#ifndef REF_COORD_H_
#define REF_COORD_H_
#include <stdint.h>
#include <iostream>
#include <limits>
#include "assert_helpers.h"
typedef int64_t TRefId;
typedef int64_t TRefOff;
class Coord
{
public:
	Coord() { reset(); }
	Coord(const Coord &c) { init(c); }
	Coord(TRefId rf, TRefOff of, bool fw) { init(rf, of, fw); }
	void init(TRefId rf, TRefOff of, bool fw)
	{
		ref_ = rf;
		off_ = of;
		orient_ = (fw ? 1 : 0);
	}
	void init(const Coord &c)
	{
		ref_ = c.ref_;
		off_ = c.off_;
		orient_ = c.orient_;
	}
	bool operator==(const Coord &o) const
	{
		assert(inited());
		assert(o.inited());
		return ref_ == o.ref_ && off_ == o.off_ && fw() == o.fw();
	}
	bool operator<(const Coord &o) const
	{
		if (ref_ < o.ref_)
			return true;
		if (ref_ > o.ref_)
			return false;
		if (orient_ < o.orient_)
			return true;
		if (orient_ > o.orient_)
			return false;
		if (off_ < o.off_)
			return true;
		if (off_ > o.off_)
			return false;
		return false;
	}
	bool operator>=(const Coord &o) const
	{
		return !((*this) < o);
	}
	bool operator>(const Coord &o) const
	{
		if (ref_ > o.ref_)
			return true;
		if (ref_ < o.ref_)
			return false;
		if (orient_ > o.orient_)
			return true;
		if (orient_ < o.orient_)
			return false;
		if (off_ > o.off_)
			return true;
		if (off_ < o.off_)
			return false;
		return false;
	}
	bool operator<=(const Coord &o) const
	{
		return !((*this) > o);
	}
	void reset()
	{
		ref_ = std::numeric_limits<TRefId>::max();
		off_ = std::numeric_limits<TRefOff>::max();
		orient_ = -1;
	}
	bool inited() const
	{
		if (ref_ != std::numeric_limits<TRefId>::max() &&
			off_ != std::numeric_limits<TRefOff>::max())
		{
			assert(orient_ == 0 || orient_ == 1);
			return true;
		}
		return false;
	}
	bool fw() const
	{
		assert(inited());
		assert(orient_ == 0 || orient_ == 1);
		return orient_ == 1;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		if (ref_ != std::numeric_limits<TRefId>::max() &&
			off_ != std::numeric_limits<TRefOff>::max())
		{
			assert(orient_ == 0 || orient_ == 1);
		}
		return true;
	}
#endif
	bool within(int64_t len, int64_t inbegin, int64_t inend) const
	{
		return off_ >= inbegin && off_ + len <= inend;
	}
	inline TRefId ref() const { return ref_; }
	inline TRefOff off() const { return off_; }
	inline int orient() const { return orient_; }
	inline void setRef(TRefId id) { ref_ = id; }
	inline void setOff(TRefOff off) { off_ = off; }
	inline void adjustOff(TRefOff off) { off_ += off; }
	TRefId ref_;
	TRefOff off_;
	int orient_;
};
std::ostream &operator<<(std::ostream &out, const Coord &c);
class Interval
{
public:
	Interval() { reset(); }
	explicit Interval(const Coord &upstream, TRefOff len)
	{
		init(upstream, len);
	}
	explicit Interval(TRefId rf, TRefOff of, bool fw, TRefOff len)
	{
		init(rf, of, fw, len);
	}
	void init(const Coord &upstream, TRefOff len)
	{
		upstream_ = upstream;
		len_ = len;
	}
	void init(TRefId rf, TRefOff of, bool fw, TRefOff len)
	{
		upstream_.init(rf, of, fw);
		len_ = len;
	}
	void setOff(TRefOff of)
	{
		upstream_.setOff(of);
	}
	void setLen(TRefOff len)
	{
		len_ = len;
	}
	void reset()
	{
		upstream_.reset();
		len_ = 0;
	}
	bool inited() const
	{
		if (upstream_.inited())
		{
			assert_gt(len_, 0);
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator==(const Interval &o) const
	{
		return upstream_ == o.upstream_ &&
			   len_ == o.len_;
	}
	bool operator<(const Interval &o) const
	{
		if (upstream_ < o.upstream_)
			return true;
		if (upstream_ > o.upstream_)
			return false;
		if (len_ < o.len_)
			return true;
		return false;
	}
	bool operator>=(const Interval &o) const
	{
		return !((*this) < o);
	}
	bool operator>(const Interval &o) const
	{
		if (upstream_ > o.upstream_)
			return true;
		if (upstream_ < o.upstream_)
			return false;
		if (len_ > o.len_)
			return true;
		return false;
	}
	bool operator<=(const Interval &o) const
	{
		return !((*this) > o);
	}
	void setUpstream(const Coord &c)
	{
		upstream_ = c;
	}
	void setLength(TRefOff l)
	{
		len_ = l;
	}
	inline TRefId ref() const { return upstream_.ref(); }
	inline TRefOff off() const { return upstream_.off(); }
	inline TRefOff dnoff() const { return upstream_.off() + len_; }
	inline int orient() const { return upstream_.orient(); }
	inline Coord downstream() const
	{
		return Coord(
			upstream_.ref(),
			upstream_.off() + len_,
			upstream_.orient());
	}
	inline bool contains(const Coord &c) const
	{
		return c.ref() == ref() &&
			   c.orient() == orient() &&
			   c.off() >= off() &&
			   c.off() < dnoff();
	}
	inline bool containsIgnoreOrient(const Coord &c) const
	{
		return c.ref() == ref() &&
			   c.off() >= off() &&
			   c.off() < dnoff();
	}
	inline bool contains(const Interval &c) const
	{
		return c.ref() == ref() &&
			   c.orient() == orient() &&
			   c.off() >= off() &&
			   c.dnoff() <= dnoff();
	}
	inline bool containsIgnoreOrient(const Interval &c) const
	{
		return c.ref() == ref() &&
			   c.off() >= off() &&
			   c.dnoff() <= dnoff();
	}
	inline bool overlaps(const Interval &c) const
	{
		return c.ref() == upstream_.ref() &&
			   c.orient() == upstream_.orient() &&
			   ((off() <= c.off() && dnoff() > c.off()) ||
				(off() <= c.dnoff() && dnoff() > c.dnoff()) ||
				(c.off() <= off() && c.dnoff() > off()) ||
				(c.off() <= dnoff() && c.dnoff() > dnoff()));
	}
	inline bool overlapsIgnoreOrient(const Interval &c) const
	{
		return c.ref() == upstream_.ref() &&
			   ((off() <= c.off() && dnoff() > c.off()) ||
				(off() <= c.dnoff() && dnoff() > c.dnoff()) ||
				(c.off() <= off() && c.dnoff() > off()) ||
				(c.off() <= dnoff() && c.dnoff() > dnoff()));
	}
	inline const Coord &upstream() const { return upstream_; }
	inline TRefOff len() const { return len_; }
#ifndef NDEBUG
	bool repOk() const
	{
		assert(upstream_.repOk());
		assert_geq(len_, 0);
		return true;
	}
#endif
	inline void adjustOff(TRefOff off)
	{
		upstream_.adjustOff(off);
	}
protected:
	Coord upstream_;
	TRefOff len_;
};
std::ostream &operator<<(std::ostream &out, const Interval &c);
#endif
