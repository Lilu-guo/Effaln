#ifndef GROUP_WALK_H_
#define GROUP_WALK_H_
#include <stdint.h>
#include <limits>
#include "ds.h"
#include "aln_idx.h"
#include "read.h"
#include "reference.h"
#include "mem_ids.h"
template <typename T>
class SARangeWithOffs
{
public:
	SARangeWithOffs() { reset(); };
	SARangeWithOffs(TIndexOffU tf, size_t len, const T &o)
	{
		init(tf, len, o);
	}
	void init(TIndexOffU tf, const T &o)
	{
		topf = tf, offs = o;
	}
	void reset() { topf = std::numeric_limits<TIndexOffU>::max(); }
	bool inited() const
	{
		return topf != std::numeric_limits<TIndexOffU>::max();
	}
	size_t size() const { return offs.size(); }
	TIndexOffU topf;
	size_t len;
	T offs;
};
struct GroupWalkState
{
	GroupWalkState(int cat) : map(cat)
	{
		masks[0].setCat(cat);
		masks[1].setCat(cat);
		masks[2].setCat(cat);
		masks[3].setCat(cat);
	}
	EList<bool> masks[4];
	EList<TIndexOffU, 16> map;
};
struct WalkMetrics
{
	WalkMetrics()
	{
		reset();
	}
	void merge(const WalkMetrics &m, bool getLock = false)
	{
		if (getLock)
		{
			ThreadSafe ts(mutex_m);
			mergeImpl(m);
		}
		else
		{
			mergeImpl(m);
		}
	}
	void reset()
	{
		bwops = branches = resolves = refresolves = reports = 0;
	}
	uint64_t bwops;
	uint64_t branches;
	uint64_t resolves;
	uint64_t refresolves;
	uint64_t reports;
	MUTEX_T mutex_m;
private:
	void mergeImpl(const WalkMetrics &m)
	{
		bwops += m.bwops;
		branches += m.branches;
		resolves += m.resolves;
		refresolves += m.refresolves;
		reports += m.reports;
	}
};
struct GWElt
{
	GWElt() { reset(); }
	void reset()
	{
		offidx = range = elt = len = OFF_MASK;
		fw = false;
	}
	void init(
		TIndexOffU oi,
		bool f,
		TIndexOffU r,
		TIndexOffU e,
		TIndexOffU l)
	{
		offidx = oi;
		fw = f;
		range = r;
		elt = e;
		len = l;
	}
	bool operator==(const GWElt &o) const
	{
		return offidx == o.offidx &&
			   fw == o.fw &&
			   range == o.range &&
			   elt == o.elt &&
			   len == o.len;
	}
	bool operator!=(const GWElt &o) const
	{
		return !(*this == o);
	}
	TIndexOffU offidx;
	bool fw;
	TIndexOffU range;
	TIndexOffU elt;
	TIndexOffU len;
};
struct WalkResult
{
	WalkResult() { reset(); }
	void reset()
	{
		elt.reset();
		bwrow = toff = OFF_MASK;
	}
	void init(
		TIndexOffU oi,
		bool f,
		TIndexOffU r,
		TIndexOffU e,
		TIndexOffU bwr,
		TIndexOffU len,
		TIndexOffU to)
	{
		elt.init(oi, f, r, e, len);
		bwrow = bwr;
		toff = to;
	}
	GWElt elt;
	TIndexOffU bwrow;
	TIndexOffU toff;
};
template <typename T>
class GWHit
{
public:
	GWHit() : fmap(0, GW_CAT),
			  offidx(OFF_MASK),
			  fw(false),
			  range(OFF_MASK),
			  len(OFF_MASK),
			  reported_(0, GW_CAT),
			  nrep_(0)
	{
		assert(repOkBasic());
	}
	void init(
		SARangeWithOffs<T> &sa,
		TIndexOffU oi,
		bool f,
		TIndexOffU r)
	{
		nrep_ = 0;
		offidx = oi;
		fw = f;
		range = r;
		len = (TIndexOffU)sa.len;
		reported_.resize(sa.offs.size());
		reported_.fill(false);
		fmap.resize(sa.offs.size());
		fmap.fill(make_pair(OFF_MASK, OFF_MASK));
	}
	void reset()
	{
		reported_.clear();
		fmap.clear();
		nrep_ = 0;
		offidx = OFF_MASK;
		fw = false;
		range = OFF_MASK;
		len = OFF_MASK;
	}
#ifndef NDEBUG
	bool repOk(const SARangeWithOffs<T> &sa) const
	{
		assert_eq(reported_.size(), sa.offs.size());
		assert_eq(fmap.size(), sa.offs.size());
		size_t nrep = 0;
		for (size_t i = 0; i < fmap.size(); i++)
		{
			if (reported_[i])
				nrep++;
			if (sa.offs[i] != OFF_MASK)
			{
				continue;
			}
			for (size_t j = i + 1; j < fmap.size(); j++)
			{
				if (sa.offs[j] != OFF_MASK)
				{
					continue;
				}
				assert(fmap[i] != fmap[j]);
			}
		}
		assert_eq(nrep_, nrep);
		return true;
	}
	bool repOkBasic()
	{
		return true;
	}
#endif
	void setReported(size_t i)
	{
		assert(!reported_[i]);
		assert_lt(i, reported_.size());
		reported_[i] = true;
		nrep_++;
	}
	bool reported(size_t i) const
	{
		assert_lt(i, reported_.size());
		return reported_[i];
	}
	bool done() const
	{
		assert_leq(nrep_, reported_.size());
		return nrep_ == reported_.size();
	}
	EList<std::pair<TIndexOffU, TIndexOffU>, 16> fmap;
	TIndexOffU offidx;
	bool fw;
	TIndexOffU range;
	TIndexOffU len;
protected:
	EList<bool, 16> reported_;
	size_t nrep_;
};
template <typename T>
class GWState
{
public:
	GWState() : map_(0, GW_CAT)
	{
		reset();
		assert(repOkBasic());
	}
	template <int S>
	pair<TIndexOff, TIndexOff> init(
		const Ebwt &ebwt,
		const BitPairReference &ref,
		SARangeWithOffs<T> &sa,
		EList<GWState, S> &sts,
		GWHit<T> &hit,
		TIndexOffU range,
		bool reportList,
		EList<WalkResult, 16> *res,
		TIndexOffU tp,
		TIndexOffU bt,
		TIndexOffU st,
		WalkMetrics &met)
	{
		assert_gt(bt, tp);
		assert_lt(range, sts.size());
		top = tp;
		bot = bt;
		step = st;
		assert(!inited_);
		ASSERT_ONLY(inited_ = true);
		ASSERT_ONLY(lastStep_ = step - 1);
		return init(ebwt, ref, sa, sts, hit, range, reportList, res, met);
	}
	template <int S>
	pair<TIndexOff, TIndexOff> init(
		const Ebwt &ebwt,
		const BitPairReference &ref,
		SARangeWithOffs<T> &sa,
		EList<GWState, S> &st,
		GWHit<T> &hit,
		TIndexOffU range,
		bool reportList,
		EList<WalkResult, 16> *res,
		WalkMetrics &met)
	{
		assert(inited_);
		assert_eq(step, lastStep_ + 1);
		ASSERT_ONLY(lastStep_++);
		assert_leq((TIndexOffU)step, ebwt.eh().len());
		assert_lt(range, st.size());
		pair<TIndexOff, TIndexOff> ret = make_pair(0, 0);
		TIndexOffU trimBegin = 0, trimEnd = 0;
		bool empty = true;
		for (size_t i = mapi_; i < map_.size(); i++)
		{
			bool resolved = (off(i, sa) != OFF_MASK);
			if (!resolved)
			{
				TIndexOffU bwrow = (TIndexOff)(top - mapi_ + i);
				TIndexOffU toff = ebwt.tryOffset(bwrow, i);
				ASSERT_ONLY(TIndexOffU origBwRow = sa.topf + map(i));
				assert_eq(bwrow, ebwt.walkLeft(origBwRow, step));
				if (toff != OFF_MASK)
				{
					assert_eq(toff, ebwt.getOffset(bwrow));
					met.resolves++;
					toff += step;
					assert_eq(toff, ebwt.getOffset(origBwRow));
					setOff(i, toff, sa, met);
					if (!reportList)
						ret.first++;
#if 0
					TIndexOffU tidx = OFF_MASK, tof, tlen;
					bool straddled = false;
					ebwt.joinedToTextOff(
						hit.len, 
						toff,    
						tidx,    
						tof,     
						tlen,    
						true,    
						straddled);
					if(tidx != OFF_MASK &&
					   hit.satup->key.seq != std::numeric_limits<uint64_t>::max())
					{
						uint64_t key = hit.satup->key.seq;
						for(int64_t j = tof + hit.len-1; j >= tof; j--) {
							int c = ref.getBase(tidx, j);
							assert_range(0, 3, c);
							if(c != (int)(key & 3)) {
								SString<char> jref;
								ebwt.restore(jref);
								uint64_t key2 = hit.satup->key.seq;
								for(int64_t k = toff + hit.len-1; k >= toff; k--) {
									int c = jref[k];
									assert_range(0, 3, c);
									assert_eq(c, (int)(key2 & 3));
									key2 >>= 2;
								}
								assert(false);
							}
							key >>= 2;
						}
					}
#endif
				}
			}
			if (off(i, sa) != OFF_MASK)
			{
				if (reportList && !hit.reported(map(i)))
				{
					TIndexOffU toff = off(i, sa);
					assert(res != NULL);
					res->expand();
					TIndexOffU origBwRow = sa.topf + map(i);
					res->back().init(
						hit.offidx,
						hit.fw,
						hit.range,
						map(i),
						origBwRow,
						hit.len,
						toff);
					hit.setReported(map(i));
					met.reports++;
				}
				if (empty)
				{
					trimBegin++;
				}
				else
				{
					trimEnd++;
				}
			}
			else
			{
				ret.second++;
				trimEnd = 0;
				empty = false;
				assert_geq(i, mapi_);
				TIndexOffU bmap = map(i);
				hit.fmap[bmap].first = range;
				hit.fmap[bmap].second = (TIndexOffU)i;
#ifndef NDEBUG
				for (size_t j = 0; j < bmap; j++)
				{
					if (sa.offs[j] == OFF_MASK &&
						hit.fmap[j].first == range)
					{
						assert_neq(i, hit.fmap[j].second);
					}
				}
#endif
			}
		}
		assert_geq(trimBegin, 0);
		mapi_ += trimBegin;
		top += trimBegin;
		if (trimEnd > 0)
		{
			map_.resize(map_.size() - trimEnd);
			bot -= trimEnd;
		}
		if (empty)
		{
			assert(done());
#ifndef NDEBUG
			for (size_t i = mapi_; i < map_.size(); i++)
			{
				assert_neq(OFF_MASK, off(i, sa));
			}
			for (size_t i = 0; i < hit.fmap.size(); i++)
			{
				if (sa.offs[i] == OFF_MASK)
				{
					assert_neq(range, hit.fmap[i].first);
				}
			}
#endif
			return ret;
		}
		else
		{
			assert(!done());
		}
		assert_neq(top, ebwt._zOff);
		assert_neq(bot - 1, ebwt._zOff);
		if (ebwt._zOff > top && ebwt._zOff < bot - 1)
		{
			TIndexOffU oldbot = bot;
			bot = ebwt._zOff;
			st.expand();
			st.back().reset();
			TIndexOffU ztop = ebwt._zOff + 1;
			st.back().initMap(oldbot - ztop);
			assert_eq(map_.size(), oldbot - top + mapi_);
			for (size_t i = ztop; i < oldbot; i++)
			{
				st.back().map_[i - ztop] = map(i - top + mapi_);
			}
			map_.resize(bot - top + mapi_);
			st.back().init(
				ebwt,
				ref,
				sa,
				st,
				hit,
				(TIndexOffU)st.size() - 1,
				reportList,
				res,
				ztop,
				oldbot,
				step,
				met);
		}
		assert_gt(bot, top);
		if (bot - top > 1)
		{
			SideLocus::initFromTopBot(top, bot, ebwt.eh(), ebwt.ebwt(), tloc, bloc);
			assert(tloc.valid());
			assert(tloc.repOk(ebwt.eh()));
			assert(bloc.valid());
			assert(bloc.repOk(ebwt.eh()));
		}
		else
		{
			tloc.initFromRow(top, ebwt.eh(), ebwt.ebwt());
			assert(tloc.valid());
			assert(tloc.repOk(ebwt.eh()));
			bloc.invalidate();
		}
		return ret;
	}
#ifndef NDEBUG
	bool repOk(
		const Ebwt &ebwt,
		GWHit<T> &hit,
		TIndexOffU range) const
	{
		assert(done() || bot > top);
		assert(doneResolving(hit) || (tloc.valid() && tloc.repOk(ebwt.eh())));
		assert(doneResolving(hit) || bot == top + 1 || (bloc.valid() && bloc.repOk(ebwt.eh())));
		assert_eq(map_.size() - mapi_, bot - top);
		TIndexOff left = 0;
		for (size_t i = mapi_; i < map_.size(); i++)
		{
			ASSERT_ONLY(TIndexOffU row = (TIndexOffU)(top + i - mapi_));
			ASSERT_ONLY(TIndexOffU origRow = hit.satup->topf + map(i));
			assert(step == 0 || row != origRow);
			assert_eq(row, ebwt.walkLeft(origRow, step));
			assert_lt(map_[i], hit.satup->offs.size());
			if (off(i, hit) == OFF_MASK)
				left++;
		}
		assert(repOkMapRepeats());
		assert(repOkMapInclusive(hit, range));
		return true;
	}
	bool repOkBasic()
	{
		assert_geq(bot, top);
		return true;
	}
	bool repOkMapInclusive(GWHit<T> &hit, TIndexOffU range) const
	{
		for (size_t i = 0; i < hit.fmap.size(); i++)
		{
			if (hit.satup->offs[i] == OFF_MASK)
			{
				if (range == hit.fmap[i].first)
				{
					ASSERT_ONLY(bool found = false);
					for (size_t j = mapi_; j < map_.size(); j++)
					{
						if (map(j) == i)
						{
							ASSERT_ONLY(found = true);
							break;
						}
					}
					assert(found);
				}
			}
		}
		return true;
	}
	bool repOkMapRepeats() const
	{
		for (size_t i = mapi_; i < map_.size(); i++)
		{
			for (size_t j = i + 1; j < map_.size(); j++)
			{
				assert_neq(map_[i], map_[j]);
			}
		}
		return true;
	}
#endif
	TIndexOffU off(
		size_t i,
		const SARangeWithOffs<T> &sa)
	{
		assert_geq(i, mapi_);
		assert_lt(i, map_.size());
		assert_lt(map_[i], sa.offs.size());
		return sa.offs.get(map_[i]);
	}
	TIndexOffU map(size_t i) const
	{
		assert_geq(i, mapi_);
		assert_lt(i, map_.size());
		return map_[i];
	}
	TIndexOffU mapi() const
	{
		return mapi_;
	}
	size_t size() const
	{
		return map_.size() - mapi_;
	}
	bool done() const
	{
		return size() == 0;
	}
	void setOff(
		size_t i,
		TIndexOffU off,
		SARangeWithOffs<T> &sa,
		WalkMetrics &met)
	{
		assert_lt(i + mapi_, map_.size());
		assert_lt(map_[i + mapi_], sa.offs.size());
		size_t saoff = map_[i + mapi_];
		sa.offs[saoff] = off;
		assert_eq(off, sa.offs[saoff]);
	}
	template <int S>
	pair<TIndexOff, TIndexOff> advance(
		const Ebwt &ebwt,
		const BitPairReference &ref,
		SARangeWithOffs<T> &sa,
		GWHit<T> &hit,
		TIndexOffU range,
		bool reportList,
		EList<WalkResult, 16> *res,
		EList<GWState, S> &st,
		GroupWalkState &gws,
		WalkMetrics &met,
		PerReadMetrics &prm)
	{
		ASSERT_ONLY(TIndexOffU origTop = top);
		ASSERT_ONLY(TIndexOffU origBot = bot);
		assert_geq(step, 0);
		assert_eq(step, lastStep_);
		assert_geq(st.capacity(), st.size() + 4);
		assert(tloc.valid());
		assert(tloc.repOk(ebwt.eh()));
		assert_eq(bot - top, map_.size() - mapi_);
		pair<TIndexOff, TIndexOff> ret = make_pair(0, 0);
		assert_eq(top, tloc.toBWRow());
		if (bloc.valid())
		{
			assert_lt(top + 1, bot);
			TIndexOffU upto[4], in[4];
			upto[0] = in[0] = upto[1] = in[1] =
				upto[2] = in[2] = upto[3] = in[3] = 0;
			assert_eq(bot, bloc.toBWRow());
			met.bwops++;
			prm.nExFmops++;
			assert(bot <= ebwt._zOff || top > ebwt._zOff);
			ebwt.mapLFRange(tloc, bloc, bot - top, upto, in, gws.masks);
#ifndef NDEBUG
			for (int i = 0; i < 4; i++)
			{
				assert_eq(bot - top, gws.masks[i].size());
			}
#endif
			bool first = true;
			ASSERT_ONLY(TIndexOffU sum = 0);
			TIndexOffU newtop = 0, newbot = 0;
			gws.map.clear();
			for (int i = 0; i < 4; i++)
			{
				if (in[i] > 0)
				{
					if (first)
					{
						first = false;
						newtop = upto[i];
						newbot = newtop + in[i];
						assert_leq(newbot - newtop, bot - top);
						for (size_t j = 0; j < gws.masks[i].size(); j++)
						{
							assert_lt(j + mapi_, map_.size());
							if (gws.masks[i][j])
							{
								gws.map.push_back(map_[j + mapi_]);
								assert(gws.map.size() <= 1 || gws.map.back() != gws.map[gws.map.size() - 2]);
#ifndef NDEBUG
								assert_lt(gws.map.back(), sa.size());
								if (sa.offs[gws.map.back()] == OFF_MASK)
								{
									assert_eq(newtop + gws.map.size() - 1,
											  ebwt.walkLeft(sa.topf + gws.map.back(), step + 1));
								}
#endif
							}
						}
						assert_eq(newbot - newtop, gws.map.size());
					}
					else
					{
						st.expand();
						st.back().reset();
						TIndexOffU ntop = upto[i];
						TIndexOffU nbot = ntop + in[i];
						assert_lt(nbot - ntop, bot - top);
						st.back().mapi_ = 0;
						st.back().map_.clear();
						met.branches++;
						for (size_t j = 0; j < gws.masks[i].size(); j++)
						{
							if (gws.masks[i][j])
								st.back().map_.push_back(map_[j + mapi_]);
						}
						pair<TIndexOff, TIndexOff> rret =
							st.back().init(
								ebwt,
								ref,
								sa,
								st,
								hit,
								(TIndexOffU)st.size() - 1,
								reportList,
								res,
								ntop,
								nbot,
								step + 1,
								met);
						ret.first += rret.first;
						ret.second += rret.second;
					}
					ASSERT_ONLY(sum += in[i]);
				}
			}
			mapi_ = 0;
			assert_eq(bot - top, sum);
			assert_gt(newbot, newtop);
			assert_leq(newbot - newtop, bot - top);
			assert(top != newtop || bot != newbot);
			top = newtop;
			bot = newbot;
			if (!gws.map.empty())
			{
				map_ = gws.map;
			}
			assert_eq(bot - top, map_.size());
		}
		else
		{
			assert_eq(bot, top + 1);
			assert_eq(1, map_.size() - mapi_);
			ASSERT_ONLY(TIndexOffU oldtop = top);
			met.bwops++;
			prm.nExFmops++;
			ebwt.mapLF1(top, tloc);
			assert_neq(top, oldtop);
			bot = top + 1;
			if (mapi_ > 0)
			{
				map_[0] = map_[mapi_];
				mapi_ = 0;
			}
			map_.resize(1);
		}
		assert(top != origTop || bot != origBot);
		step++;
		assert_gt(step, 0);
		assert_leq((TIndexOffU)step, ebwt.eh().len());
		pair<TIndexOff, TIndexOff> rret =
			init<S>(
				ebwt,
				ref,
				sa,
				st,
				hit,
				range,
				reportList,
				res,
				met);
		ret.first += rret.first;
		ret.second += rret.second;
		return ret;
	}
	void reset()
	{
		top = bot = step = mapi_ = 0;
		ASSERT_ONLY(lastStep_ = -1);
		ASSERT_ONLY(inited_ = false);
		tloc.invalidate();
		bloc.invalidate();
		map_.clear();
	}
	void initMap(size_t newsz)
	{
		mapi_ = 0;
		map_.resize(newsz);
		for (size_t i = 0; i < newsz; i++)
		{
			map_[i] = (TIndexOffU)i;
		}
	}
	bool doneReporting(const GWHit<T> &hit) const
	{
		for (size_t i = mapi_; i < map_.size(); i++)
		{
			if (!hit.reported(map(i)))
				return false;
		}
		return true;
	}
	bool doneResolving(const SARangeWithOffs<T> &sa) const
	{
		for (size_t i = mapi_; i < map_.size(); i++)
		{
			if (sa.offs[map(i)] == OFF_MASK)
				return false;
		}
		return true;
	}
	SideLocus tloc;
	SideLocus bloc;
	TIndexOffU top;
	TIndexOffU bot;
	TIndexOff step;
protected:
	ASSERT_ONLY(bool inited_);
	ASSERT_ONLY(TIndexOff lastStep_);
	EList<TIndexOffU, 16> map_;
	TIndexOffU mapi_;
};
template <typename T, int S>
class GroupWalk2S
{
public:
	typedef EList<GWState<T>, S> TStateV;
	GroupWalk2S() : st_(8, GW_CAT)
	{
		reset();
	}
	void reset()
	{
		elt_ = rep_ = 0;
		ASSERT_ONLY(inited_ = false);
	}
	void init(
		const Ebwt &ebwtFw,
		const BitPairReference &ref,
		SARangeWithOffs<T> &sa,
		RandomSource &rnd,
		WalkMetrics &met)
	{
		reset();
#ifndef NDEBUG
		inited_ = true;
#endif
		hit_.init(sa, 0, false, 0);
		st_.resize(1);
		st_.back().reset();
		assert(st_.back().repOkBasic());
		TIndexOffU top = sa.topf;
		TIndexOffU bot = (TIndexOffU)(top + sa.size());
		st_.back().initMap(bot - top);
		st_.ensure(4);
		st_.back().init(
			ebwtFw,
			ref,
			sa,
			st_,
			hit_,
			0,
			false,
			NULL,
			top,
			bot,
			0,
			met);
		elt_ += sa.size();
		assert(hit_.repOk(sa));
	}
	void resolveAll(WalkMetrics &met, PerReadMetrics &prm)
	{
		WalkResult res;
		for (size_t i = 0; i < elt_; i++)
		{
			advanceElement((TIndexOffU)i, res, met, prm);
		}
	}
	bool advanceElement(
		TIndexOffU elt,
		const Ebwt &ebwtFw,
		const BitPairReference &ref,
		SARangeWithOffs<T> &sa,
		GroupWalkState &gws,
		WalkResult &res,
		WalkMetrics &met,
		PerReadMetrics &prm)
	{
		assert(inited_);
		assert(!done());
		assert(hit_.repOk(sa));
		assert_lt(elt, sa.size());
		while (sa.offs[elt] == OFF_MASK)
		{
			size_t range = hit_.fmap[elt].first;
			st_.ensure(4);
			GWState<T> &st = st_[range];
			assert(!st.doneResolving(sa));
			st.advance(
				ebwtFw,
				ref,
				sa,
				hit_,
				(TIndexOffU)range,
				false,
				NULL,
				st_,
				gws,
				met,
				prm);
			assert(sa.offs[elt] != OFF_MASK ||
				   !st_[hit_.fmap[elt].first].doneResolving(sa));
		}
		assert_neq(OFF_MASK, sa.offs[elt]);
		if (!hit_.reported(elt))
		{
			hit_.setReported(elt);
		}
		met.reports++;
		res.init(
			0,
			false,
			0,
			elt,
			sa.topf + elt,
			(TIndexOffU)sa.len,
			sa.offs[elt]);
		rep_++;
		return true;
	}
	bool done() const { return rep_ == elt_; }
#ifndef NDEBUG
	bool repOk(const SARangeWithOffs<T> &sa) const
	{
		assert(hit_.repOk(sa));
		assert_leq(rep_, elt_);
		size_t resolved = 0, reported = 0;
		const size_t sz = sa.size();
		for (size_t m = 0; m < sz; m++)
		{
			if (sa.offs[m] != OFF_MASK)
			{
				resolved++;
			}
			else
			{
				assert(!hit_.reported(m));
			}
			if (hit_.reported(m))
			{
				reported++;
			}
			assert_geq(resolved, reported);
		}
		assert_geq(resolved, reported);
		assert_eq(rep_, reported);
		assert_eq(elt_, sz);
		return true;
	}
#endif
	size_t numElts() const
	{
		return elt_;
	}
	size_t totalSizeBytes() const
	{
		return 2 * sizeof(size_t) + st_.totalSizeBytes() + sizeof(GWHit<T>);
	}
	size_t totalCapacityBytes() const
	{
		return 2 * sizeof(size_t) + st_.totalCapacityBytes() + sizeof(GWHit<T>);
	}
#ifndef NDEBUG
	bool initialized() const
	{
		return inited_;
	}
#endif
protected:
	ASSERT_ONLY(bool inited_);
	size_t elt_;
	size_t rep_;
	TStateV st_;
	GWHit<T> hit_;
};
#endif