#ifndef DP_FRAMER_H_
#define DP_FRAMER_H_
#include <stdint.h>
#include "ds.h"
#include "ref_coord.h"
struct DPRect
{
	DPRect(int cat = 0)
	{
		refl = refr = triml = trimr = corel = corer = 0;
	}
	int64_t refl;
	int64_t refr;
	int64_t refl_pretrim;
	int64_t refr_pretrim;
	size_t triml;
	size_t trimr;
	size_t corel;
	size_t corer;
	size_t maxgap;
	void write(std::ostream &os) const
	{
		os << refl << ',' << refr << ',' << refl_pretrim << ','
		   << refr_pretrim << ',' << triml << ',' << trimr << ','
		   << corel << ',' << corer << ',' << maxgap;
	}
	bool entirelyTrimmed() const
	{
		bool tr = refr < refl;
		ASSERT_ONLY(size_t width = (size_t)(refr_pretrim - refl_pretrim + 1));
		assert(tr == (width <= triml + trimr));
		return tr;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_geq(corer, corel);
		return true;
	}
#endif
	void initIval(Interval &iv)
	{
		iv.setOff(refl_pretrim + (int64_t)corel);
		iv.setLen(corer - corel + 1);
	}
};
class DynProgFramer
{
public:
	DynProgFramer(bool trimToRef) : trimToRef_(trimToRef) {}
	bool frameSeedExtensionRect(
		int64_t off,
		size_t rdlen, int64_t reflen, size_t maxrdgap, size_t maxrfgap, int64_t maxns, size_t maxhalf, DPRect &rect);
	bool frameFindMateRect(
		bool anchorLeft,
		int64_t ll, int64_t lr, int64_t rl, int64_t rr, size_t rdlen, int64_t reflen, size_t maxrdgap, size_t maxrfgap, int64_t maxns, size_t maxhalf, DPRect &rect) const
	{
		if (anchorLeft)
		{
			return frameFindMateAnchorLeftRect(
				ll,
				lr,
				rl,
				rr,
				rdlen,
				reflen,
				maxrdgap,
				maxrfgap,
				maxns,
				maxhalf,
				rect);
		}
		else
		{
			return frameFindMateAnchorRightRect(
				ll,
				lr,
				rl,
				rr,
				rdlen,
				reflen,
				maxrdgap,
				maxrfgap,
				maxns,
				maxhalf,
				rect);
		}
	}
	bool frameFindMateAnchorLeftRect(
		int64_t ll,
		int64_t lr, int64_t rl, int64_t rr, size_t rdlen, int64_t reflen, size_t maxrdgap, size_t maxrfgap, int64_t maxns, size_t maxhalf, DPRect &rect) const;
	bool frameFindMateAnchorRightRect(
		int64_t ll,
		int64_t lr, int64_t rl, int64_t rr, size_t rdlen, int64_t reflen, size_t maxrdgap, size_t maxrfgap, int64_t maxns, size_t maxhalf, DPRect &rect) const;
protected:
	void trimToRef(
		size_t reflen,
		int64_t &refl, int64_t &refr, size_t &trimup, size_t &trimdn)
	{
		if (refl < 0)
		{
			trimup = (size_t)(-refl);
		}
		if (refr >= (int64_t)reflen)
		{
			trimdn = (size_t)(refr - reflen + 1);
		}
	}
	bool trimToRef_;
};
#endif
