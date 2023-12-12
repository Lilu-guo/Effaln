#ifndef ALIGNER_CACHE_H_
#define ALIGNER_CACHE_H_
#include <iostream>
#include "ds.h"
#include "read.h"
#include "threading.h"
#include "mem_ids.h"
#include "simple_func.h"
#include "btypes.h"
#define CACHE_PAGE_SZ (16 * 1024)
typedef PListSlice<TIndexOffU, CACHE_PAGE_SZ> TSlice;
struct QKey
{
	QKey() { reset(); }
	QKey(const BTDnaString &s ASSERT_ONLY(, BTDnaString &tmp))
	{
		init(s ASSERT_ONLY(, tmp));
	}
	bool init(
		const BTDnaString &s
			ASSERT_ONLY(, BTDnaString &tmp))
	{
		seq = 0;
		len = (uint32_t)s.length();
		ASSERT_ONLY(tmp.clear());
		if (len > 32)
		{
			len = 0xffffffff;
			return false;
		}
		else
		{
			for (size_t i = 0; i < 32 && i < s.length(); i++)
			{
				int c = (int)s.get(i);
				assert_range(0, 4, c);
				if (c == 4)
				{
					len = 0xffffffff;
					return false;
				}
				seq = (seq << 2) | s.get(i);
			}
			ASSERT_ONLY(toString(tmp));
			assert(sstr_eq(tmp, s));
			assert_leq(len, 32);
			return true;
		}
	}
	void toString(BTDnaString &s)
	{
		s.resize(len);
		uint64_t sq = seq;
		for (int i = (len)-1; i >= 0; i--)
		{
			s.set((uint32_t)(sq & 3), i);
			sq >>= 2;
		}
	}
	bool cacheable() const { return len != 0xffffffff; }
	void reset()
	{
		seq = 0;
		len = 0xffffffff;
	}
	bool operator<(const QKey &o) const
	{
		return seq < o.seq || (seq == o.seq && len < o.len);
	}
	bool operator>(const QKey &o) const
	{
		return !(*this < o || *this == o);
	}
	bool operator==(const QKey &o) const
	{
		return seq == o.seq && len == o.len;
	}
	bool operator!=(const QKey &o) const
	{
		return !(*this == o);
	}
#ifndef NDEBUG
	bool repOk() const
	{
		return len != 0xffffffff;
	}
#endif
	uint64_t seq;
	uint32_t len;
};
class AlignmentCache;
class QVal
{
public:
	QVal() { reset(); }
	TIndexOffU offset() const { return i_; }
	TIndexOffU numRanges() const
	{
		assert(valid());
		return rangen_;
	}
	TIndexOffU numElts() const
	{
		assert(valid());
		return eltn_;
	}
	bool empty() const
	{
		assert(valid());
		return numRanges() == 0;
	}
	bool valid() const { return rangen_ != OFF_MASK; }
	void reset()
	{
		i_ = 0;
		rangen_ = eltn_ = OFF_MASK;
	}
	void init(TIndexOffU i, TIndexOffU ranges, TIndexOffU elts)
	{
		i_ = i;
		rangen_ = ranges;
		eltn_ = elts;
	}
	void addRange(TIndexOffU numElts)
	{
		rangen_++;
		eltn_ += numElts;
	}
#ifndef NDEBUG
	bool repOk(const AlignmentCache &ac) const;
#endif
protected:
	TIndexOffU i_;
	TIndexOffU rangen_;
	TIndexOffU eltn_;
};
typedef QKey SAKey;
struct SAVal
{
	SAVal() : topf(), topb(), i(), len(OFF_MASK) {}
	bool valid() { return len != OFF_MASK; }
#ifndef NDEBUG
	bool repOk(const AlignmentCache &ac) const;
#endif
	void init(
		TIndexOffU tf,
		TIndexOffU tb,
		TIndexOffU ii,
		TIndexOffU ln)
	{
		topf = tf;
		topb = tb;
		i = ii;
		len = ln;
	}
	TIndexOffU topf;
	TIndexOffU topb;
	TIndexOffU i;
	TIndexOffU len;
};
class SATuple
{
public:
	SATuple() { reset(); };
	SATuple(SAKey k, TIndexOffU tf, TIndexOffU tb, TSlice o)
	{
		init(k, tf, tb, o);
	}
	void init(SAKey k, TIndexOffU tf, TIndexOffU tb, TSlice o)
	{
		key = k;
		topf = tf;
		topb = tb;
		offs = o;
	}
	void init(const SATuple &src, size_t first, size_t last)
	{
		assert_neq(OFF_MASK, src.topb);
		key = src.key;
		topf = (TIndexOffU)(src.topf + first);
		topb = OFF_MASK;
		offs.init(src.offs, first, last);
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert(offs.repOk());
		return true;
	}
#endif
	bool operator<(const SATuple &o) const
	{
		if (offs.size() < o.offs.size())
		{
			return true;
		}
		if (offs.size() > o.offs.size())
		{
			return false;
		}
		return topf < o.topf;
	}
	bool operator>(const SATuple &o) const
	{
		if (offs.size() < o.offs.size())
		{
			return false;
		}
		if (offs.size() > o.offs.size())
		{
			return true;
		}
		return topf > o.topf;
	}
	bool operator==(const SATuple &o) const
	{
		return key == o.key && topf == o.topf && topb == o.topb && offs == o.offs;
	}
	void reset()
	{
		topf = topb = OFF_MASK;
		offs.reset();
	}
	void setLength(size_t nlen)
	{
		assert_leq(nlen, offs.size());
		offs.setLength(nlen);
	}
	size_t size() const { return offs.size(); }
	SAKey key;
	TIndexOffU topf;
	TIndexOffU topb;
	TSlice offs;
};
class AlignmentCache
{
	typedef RedBlackNode<QKey, QVal> QNode;
	typedef RedBlackNode<SAKey, SAVal> SANode;
	typedef PList<SAKey, CACHE_PAGE_SZ> TQList;
	typedef PList<TIndexOffU, CACHE_PAGE_SZ> TSAList;
public:
	AlignmentCache(
		uint64_t bytes,
		bool shared) : pool_(bytes, CACHE_PAGE_SZ, CA_CAT),
					   qmap_(CACHE_PAGE_SZ, CA_CAT),
					   qlist_(CA_CAT),
					   samap_(CACHE_PAGE_SZ, CA_CAT),
					   salist_(CA_CAT),
					   shared_(shared),
					   mutex_m(),
					   version_(0) {}
	template <int S>
	void queryQval(
		const QVal &qv,
		EList<SATuple, S> &satups,
		size_t &nrange,
		size_t &nelt,
		bool getLock = true)
	{
		if (shared_ && getLock)
		{
			ThreadSafe ts(mutex_m);
			queryQvalImpl(qv, satups, nrange, nelt);
		}
		else
		{
			queryQvalImpl(qv, satups, nrange, nelt);
		}
	}
	bool empty() const
	{
		bool ret = qmap_.empty();
		assert(!ret || qlist_.empty());
		assert(!ret || samap_.empty());
		assert(!ret || salist_.empty());
		return ret;
	}
	QVal *add(
		const QKey &qk,
		bool *added,
		bool getLock = true)
	{
		if (shared_ && getLock)
		{
			ThreadSafe ts(mutex_m);
			return addImpl(qk, added);
		}
		else
		{
			return addImpl(qk, added);
		}
	}
	bool addOnTheFly(
		QVal &qv, const SAKey &sak, TIndexOffU topf, TIndexOffU botf, TIndexOffU topb, TIndexOffU botb, bool getLock = true);
	void clear()
	{
		ThreadSafe ts(mutex_m);
		pool_.clear();
		qmap_.clear();
		qlist_.clear();
		samap_.clear();
		salist_.clear();
		version_++;
	}
	size_t qNumKeys() const { return qmap_.size(); }
	size_t saNumKeys() const { return samap_.size(); }
	size_t qSize() const { return qlist_.size(); }
	size_t saSize() const { return salist_.size(); }
	Pool &pool() { return pool_; }
	MUTEX_T &lock()
	{
		return mutex_m;
	}
	bool shared() const { return shared_; }
	uint32_t version() const { return version_; }
protected:
	Pool pool_;
	RedBlack<QKey, QVal> qmap_;
	TQList qlist_;
	RedBlack<SAKey, SAVal> samap_;
	TSAList salist_;
	bool shared_;
	MUTEX_T mutex_m;
	uint32_t version_;
private:
	template <int S>
	void queryQvalImpl(
		const QVal &qv,
		EList<SATuple, S> &satups,
		size_t &nrange,
		size_t &nelt)
	{
		assert(qv.repOk(*this));
		const size_t refi = qv.offset();
		const size_t reff = refi + qv.numRanges();
		for (size_t i = refi; i < reff; i++)
		{
			SAKey sak = qlist_.get(i);
			assert(i == refi || qlist_.get(i) != qlist_.get(i - 1));
			SANode *n = samap_.lookup(sak);
			assert(n != NULL);
			const SAVal &sav = n->payload;
			assert(sav.repOk(*this));
			if (sav.len > 0)
			{
				nrange++;
				satups.expand();
				satups.back().init(sak, sav.topf, sav.topb, TSlice(salist_, sav.i, sav.len));
				nelt += sav.len;
#ifndef NDEBUG
				if (i > refi)
				{
					const SATuple b1 = satups.back();
					const SATuple b2 = satups[satups.size() - 2];
					assert(b1.key != b2.key || b1.topf != b2.topf || b1.offs != b2.offs);
				}
#endif
			}
		}
	}
	bool addOnTheFlyImpl(
		QVal &qv, const SAKey &sak, TIndexOffU topf, TIndexOffU botf, TIndexOffU topb, TIndexOffU botb);
	QVal *addImpl(
		const QKey &qk,
		bool *added)
	{
		assert(qk.cacheable());
		QNode *n = qmap_.add(pool(), qk, added);
		return (n != NULL ? &n->payload : NULL);
	}
};
class AlignmentCacheIface
{
public:
	AlignmentCacheIface(
		AlignmentCache *current,
		AlignmentCache *local,
		AlignmentCache *shared) : qk_(),
								  qv_(NULL),
								  cacheable_(false),
								  rangen_(0),
								  eltsn_(0),
								  current_(current),
								  local_(local),
								  shared_(shared)
	{
		assert(current_ != NULL);
	}
#if 0
	QVal* queryCopy(const QKey& qk, bool getLock = true) {
		assert(qk.cacheable());
		AlignmentCache* caches[3] = { current_, local_, shared_ };
		for(int i = 0; i < 3; i++) {
			if(caches[i] == NULL) continue;
			QVal* qv = caches[i]->query(qk, getLock);
			if(qv != NULL) {
				if(i == 0) return qv;
				if(!current_->copy(qk, *qv, *caches[i], getLock)) {
															return NULL;
				}
				QVal* curqv = current_->query(qk, getLock);
				assert(curqv != NULL);
				return curqv;
			}
		}
		return NULL;
	}
	inline QVal* query(
		const QKey& qk,
		AlignmentCache** which,
		bool getLock = true)
	{
		assert(qk.cacheable());
		AlignmentCache* caches[3] = { current_, local_, shared_ };
		for(int i = 0; i < 3; i++) {
			if(caches[i] == NULL) continue;
			QVal* qv = caches[i]->query(qk, getLock);
			if(qv != NULL) {
				if(which != NULL) *which = caches[i];
				return qv;
			}
		}
		return NULL;
	}
#endif
	int beginAlign(
		const BTDnaString &seq,
		const BTString &qual,
		QVal &qv, bool getLock = true)
	{
		assert(repOk());
		qk_.init(seq ASSERT_ONLY(, tmpdnastr_));
		if (qk_.cacheable())
		{
			qv_ = current_->add(qk_, &cacheable_, getLock);
		}
		else
		{
			qv_ = &qvbuf_;
		}
		if (qv_ == NULL)
		{
			resetRead();
			return -1;
		}
		qv_->reset();
		return 0;
	}
	ASSERT_ONLY(BTDnaString tmpdnastr_);
	QVal finishAlign(bool getLock = true)
	{
		if (!qv_->valid())
		{
			qv_->init(0, 0, 0);
		}
		QVal *qv = qv_;
#if 0
		if(qk_.cacheable()) {
			AlignmentCache* caches[3] = { current_, local_, shared_ };
			ASSERT_ONLY(AlignmentCache* which);
			ASSERT_ONLY(QVal* qv2 = query(qk_, &which, true));
			assert(qv2 == qv);
			assert(which == current_);
			for(int i = 1; i < 3; i++) {
				if(caches[i] != NULL) {
																				caches[i]->clearCopy(qk_, *qv_, *current_, getLock);
					break;
				}
			}
		}
#endif
		resetRead();
		assert(repOk());
		return *qv;
	}
	void nextRead()
	{
		current_->clear();
		resetRead();
		assert(!aligning());
	}
	bool aligning() const
	{
		return qv_ != NULL;
	}
	void clear()
	{
		if (current_ != NULL)
			current_->clear();
		if (local_ != NULL)
			local_->clear();
		if (shared_ != NULL)
			shared_->clear();
	}
	bool addOnTheFly(
		const BTDnaString &rfseq, TIndexOffU topf, TIndexOffU botf, TIndexOffU topb, TIndexOffU botb, bool getLock = true)
	{
		assert(aligning());
		assert(repOk());
		ASSERT_ONLY(BTDnaString tmp);
		SAKey sak(rfseq ASSERT_ONLY(, tmp));
		if (current_->addOnTheFly((*qv_), sak, topf, botf, topb, botb, getLock))
		{
			rangen_++;
			eltsn_ += (botf - topf);
			return true;
		}
		return false;
	}
	template <int S>
	void queryQval(
		const QVal &qv,
		EList<SATuple, S> &satups,
		size_t &nrange,
		size_t &nelt,
		bool getLock = true)
	{
		current_->queryQval(qv, satups, nrange, nelt, getLock);
	}
	const AlignmentCache *currentCache() const { return current_; }
	size_t curNumRanges() const { return rangen_; }
	size_t curNumElts() const { return eltsn_; }
#ifndef NDEBUG
	bool repOk() const
	{
		assert(current_ != NULL);
		assert_geq(eltsn_, rangen_);
		if (qv_ == NULL)
		{
			assert_eq(0, rangen_);
			assert_eq(0, eltsn_);
		}
		return true;
	}
#endif
	const AlignmentCache &current()
	{
		return *current_;
	}
protected:
	void resetRead()
	{
		cacheable_ = false;
		rangen_ = eltsn_ = 0;
		qv_ = NULL;
	}
	QKey qk_;
	QVal *qv_;
	QVal qvbuf_;
	bool cacheable_;
	size_t rangen_;
	size_t eltsn_;
	AlignmentCache *current_;
	AlignmentCache *local_;
	AlignmentCache *shared_;
};
#endif
