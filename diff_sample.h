#ifndef DIFF_SAMPLE_H_
#define DIFF_SAMPLE_H_
#ifdef WITH_TBB
#include <tbb/tbb.h>
#include <tbb/task_group.h>
#endif
#include <stdint.h>
#include <string.h>
#include "assert_helpers.h"
#include "multikey_qsort.h"
#include "timer.h"
#include "ds.h"
#include "mem_ids.h"
#include "ls.h"
#include "btypes.h"
using namespace std;
#ifndef VMSG_NL
#define VMSG_NL(...)                \
	if (this->verbose())            \
	{                               \
		stringstream tmp;           \
		tmp << __VA_ARGS__ << endl; \
		this->verbose(tmp.str());   \
	}
#endif
#ifndef VMSG
#define VMSG(...)                 \
	if (this->verbose())          \
	{                             \
		stringstream tmp;         \
		tmp << __VA_ARGS__;       \
		this->verbose(tmp.str()); \
	}
#endif
struct sampleEntry
{
	uint32_t maxV;
	uint32_t numSamples;
	uint32_t samples[128];
};
extern struct sampleEntry clDCs[16];
extern bool clDCs_calced;
template <typename T>
static bool dcRepOk(T v, EList<T> &ds)
{
	AutoArray<bool> covered(v, EBWT_CAT);
	for (T i = 1; i < v; i++)
	{
		covered[i] = false;
	}
	for (T di = T(); di < ds.size(); di++)
	{
		for (T dj = di + 1; dj < ds.size(); dj++)
		{
			assert_lt(ds[di], ds[dj]);
			T d1 = (ds[dj] - ds[di]);
			T d2 = (ds[di] + v - ds[dj]);
			assert_lt(d1, v);
			assert_lt(d2, v);
			covered[d1] = true;
			covered[d2] = true;
		}
	}
	bool ok = true;
	for (T i = 1; i < v; i++)
	{
		if (covered[i] == false)
		{
			ok = false;
			break;
		}
	}
	return ok;
}
template <typename T>
static bool increasing(T *ts, size_t limit)
{
	for (size_t i = 0; i < limit - 1; i++)
	{
		if (ts[i + 1] <= ts[i])
			return false;
	}
	return true;
}
template <typename T>
static inline bool hasDifference(T *ds, T d, T v, T diff)
{
	for (T di = T(); di < d; di++)
	{
		for (T dj = di + 1; dj < d; dj++)
		{
			assert_lt(ds[di], ds[dj]);
			T d1 = (ds[dj] - ds[di]);
			T d2 = (ds[di] + v - ds[dj]);
			assert_lt(d1, v);
			assert_lt(d2, v);
			if (d1 == diff || d2 == diff)
				return true;
		}
	}
	return false;
}
template <typename T>
void calcExhaustiveDC(T i, bool verbose = false, bool sanityCheck = false)
{
	T v = i;
	AutoArray<bool> diffs(v, EBWT_CAT);
	T ld = (T)ceil(sqrt(v));
	T ud = v / 2;
	bool ok = true;
	T *ds = NULL;
	T d;
	for (d = ld; d <= ud + 1; d++)
	{
		AutoArray<T> ds(d, EBWT_CAT);
		for (T j = 0; j < d; j++)
		{
			ds[j] = j;
		}
		assert(increasing(ds, d));
		while (true)
		{
			for (T t = 1; t < v; t++)
			{
				diffs[t] = false;
			}
			T diffCnt = 0;
			for (T di = 0; di < d; di++)
			{
				for (T dj = di + 1; dj < d; dj++)
				{
					assert_lt(ds[di], ds[dj]);
					T d1 = (ds[dj] - ds[di]);
					T d2 = (ds[di] + v - ds[dj]);
					assert_lt(d1, v);
					assert_lt(d2, v);
					assert_gt(d1, 0);
					assert_gt(d2, 0);
					if (!diffs[d1])
					{
						diffCnt++;
						diffs[d1] = true;
					}
					if (!diffs[d2])
					{
						diffCnt++;
						diffs[d2] = true;
					}
				}
			}
			ok = diffCnt == v - 1;
			if (ok)
			{
				break;
			}
			else
			{
				assert(increasing(ds, d));
				bool advanced = false;
				bool keepGoing = false;
				do
				{
					keepGoing = false;
					for (T bd = d - 1; bd > 1; bd--)
					{
						T dif = (d - 1) - bd;
						if (ds[bd] < v - 1 - dif)
						{
							ds[bd]++;
							assert_neq(0, ds[bd]);
							for (T bdi = bd + 1; bdi < d; bdi++)
							{
								assert_eq(0, ds[bdi]);
								ds[bdi] = ds[bdi - 1] + 1;
								assert_gt(ds[bdi], ds[bdi - 1]);
							}
							assert(increasing(ds, d));
							advanced = true;
							break;
						}
						else
						{
							ds[bd] = 0;
						}
					}
				} while (keepGoing);
				if (!advanced)
					break;
				assert(increasing(ds, d));
			}
		}
		if (ok)
		{
			break;
		}
	}
	assert(ok);
	cout << "Did exhaustive v=" << v << " |D|=" << d << endl;
	cout << "  ";
	for (T i = 0; i < d; i++)
	{
		cout << ds[i];
		if (i < d - 1)
			cout << ",";
	}
	cout << endl;
}
template <typename T>
void calcColbournAndLingDCs(bool verbose = false, bool sanityCheck = false)
{
	for (T r = 0; r < 16; r++)
	{
		T maxv = 24 * r * r + 36 * r + 13;
		T numsamp = 6 * r + 4;
		clDCs[r].maxV = maxv;
		clDCs[r].numSamples = numsamp;
		memset(clDCs[r].samples, 0, 4 * 128);
		T i;
		for (i = 1; i < r + 1; i++)
		{
			clDCs[r].samples[i] = clDCs[r].samples[i - 1] + 1;
		}
		clDCs[r].samples[r + 1] = clDCs[r].samples[r] + r + 1;
		for (i = r + 2; i < r + 2 + r; i++)
		{
			clDCs[r].samples[i] = clDCs[r].samples[i - 1] + 2 * r + 1;
		}
		for (i = r + 2 + r; i < r + 2 + r + 2 * r + 1; i++)
		{
			clDCs[r].samples[i] = clDCs[r].samples[i - 1] + 4 * r + 3;
		}
		for (i = r + 2 + r + 2 * r + 1; i < r + 2 + r + 2 * r + 1 + r + 1; i++)
		{
			clDCs[r].samples[i] = clDCs[r].samples[i - 1] + 2 * r + 2;
		}
		for (i = r + 2 + r + 2 * r + 1 + r + 1; i < r + 2 + r + 2 * r + 1 + r + 1 + r; i++)
		{
			clDCs[r].samples[i] = clDCs[r].samples[i - 1] + 1;
		}
		assert_eq(i, numsamp);
		assert_lt(i, 128);
		if (sanityCheck)
		{
			AutoArray<bool> diffs(maxv, EBWT_CAT);
			for (T i = 0; i < numsamp; i++)
			{
				for (T j = i + 1; j < numsamp; j++)
				{
					T d1 = (clDCs[r].samples[j] - clDCs[r].samples[i]);
					T d2 = (clDCs[r].samples[i] + maxv - clDCs[r].samples[j]);
					assert_lt(d1, maxv);
					assert_lt(d2, maxv);
					diffs[d1] = true;
					diffs[d2] = true;
				}
			}
			for (T i = 1; i < maxv; i++)
			{
				if (diffs[i] == false)
					cout << r << ", " << i << endl;
				assert(diffs[i] == true);
			}
		}
	}
	clDCs_calced = true;
}
extern uint32_t dc0to64[65][10];
template <typename T>
static EList<T> getDiffCover(
	T v,
	bool verbose = false,
	bool sanityCheck = false)
{
	assert_gt(v, 2);
	EList<T> ret;
	ret.clear();
	if (v <= 64 && dc0to64[v][0] == 0xffffffff)
	{
		if (verbose)
			cout << "v in hardcoded area, but hardcoded entry was all-fs" << endl;
		return ret;
	}
	else if (v <= 64)
	{
		ret.push_back(0);
		for (size_t i = 0; i < 10; i++)
		{
			if (dc0to64[v][i] == 0)
				break;
			ret.push_back(dc0to64[v][i]);
		}
		if (sanityCheck)
			assert(dcRepOk(v, ret));
		return ret;
	}
	if (!clDCs_calced)
	{
		calcColbournAndLingDCs<uint32_t>(verbose, sanityCheck);
		assert(clDCs_calced);
	}
	for (size_t i = 0; i < 16; i++)
	{
		if (v <= clDCs[i].maxV)
		{
			for (size_t j = 0; j < clDCs[i].numSamples; j++)
			{
				T s = clDCs[i].samples[j];
				if (s >= v)
				{
					s %= v;
					for (size_t k = 0; k < ret.size(); k++)
					{
						if (s == ret[k])
							break;
						if (s < ret[k])
						{
							ret.insert(s, k);
							break;
						}
					}
				}
				else
				{
					ret.push_back(s % v);
				}
			}
			if (sanityCheck)
				assert(dcRepOk(v, ret));
			return ret;
		}
	}
	cerr << "Error: Could not find a difference cover sample for v=" << v << endl;
	throw 1;
}
template <typename T>
static EList<T> getDeltaMap(T v, const EList<T> &dc)
{
	EList<T> amap;
	size_t amapEnts = 1;
	amap.resizeExact((size_t)v);
	amap.fill(0xffffffff);
	amap[0] = 0;
	for (size_t i = 0; i < dc.size(); i++)
	{
		for (size_t j = i + 1; j < dc.size(); j++)
		{
			assert_gt(dc[j], dc[i]);
			T diffLeft = dc[j] - dc[i];
			T diffRight = dc[i] + v - dc[j];
			assert_lt(diffLeft, v);
			assert_lt(diffRight, v);
			if (amap[diffLeft] == 0xffffffff)
			{
				amap[diffLeft] = dc[i];
				amapEnts++;
			}
			if (amap[diffRight] == 0xffffffff)
			{
				amap[diffRight] = dc[j];
				amapEnts++;
			}
		}
	}
	return amap;
}
template <typename T>
static unsigned int popCount(T i)
{
	unsigned int cnt = 0;
	for (size_t j = 0; j < sizeof(T) * 8; j++)
	{
		if (i & 1)
			cnt++;
		i >>= 1;
	}
	return cnt;
}
template <typename T>
static unsigned int myLog2(T i)
{
	assert_eq(1, popCount(i));
	for (size_t j = 0; j < sizeof(T) * 8; j++)
	{
		if (i & 1)
			return (int)j;
		i >>= 1;
	}
	assert(false);
	return 0xffffffff;
}
template <typename TStr>
class DifferenceCoverSample
{
public:
	DifferenceCoverSample(const TStr &__text,
						  uint32_t __v,
						  bool __verbose = false,
						  bool __sanity = false,
						  ostream &__logger = cout) : _text(__text),
													  _v(__v),
													  _verbose(__verbose),
													  _sanity(__sanity),
													  _ds(getDiffCover(_v, _verbose, _sanity)),
													  _dmap(getDeltaMap(_v, _ds)),
													  _d((uint32_t)_ds.size()),
													  _doffs(),
													  _isaPrime(),
													  _dInv(),
													  _log2v(myLog2(_v)),
													  _vmask(OFF_MASK << _log2v),
													  _logger(__logger)
	{
		assert_gt(_d, 0);
		assert_eq(1, popCount(_v));
		_dInv.resizeExact((size_t)v());
		_dInv.fill(0xffffffff);
		uint32_t lim = (uint32_t)_ds.size();
		for (uint32_t i = 0; i < lim; i++)
		{
			_dInv[_ds[i]] = i;
		}
	}
	static size_t simulateAllocs(const TStr &text, uint32_t v)
	{
		EList<uint32_t> ds(getDiffCover(v, false, false));
		size_t len = text.length();
		size_t sPrimeSz = (len / v) * ds.size();
		AutoArray<TIndexOffU> aa(sPrimeSz * 3 + (1024 * 1024), EBWT_CAT);
		return sPrimeSz * 4;
	}
	uint32_t v() const { return _v; }
	uint32_t log2v() const { return _log2v; }
	uint32_t vmask() const { return _vmask; }
	uint32_t modv(TIndexOffU i) const { return (uint32_t)(i & ~_vmask); }
	TIndexOffU divv(TIndexOffU i) const { return i >> _log2v; }
	uint32_t d() const { return _d; }
	bool verbose() const { return _verbose; }
	bool sanityCheck() const { return _sanity; }
	const TStr &text() const { return _text; }
	const EList<uint32_t> &ds() const { return _ds; }
	const EList<uint32_t> &dmap() const { return _dmap; }
	ostream &log() const { return _logger; }
	void build(int nthreads);
	uint32_t tieBreakOff(TIndexOffU i, TIndexOffU j) const;
	int64_t breakTie(TIndexOffU i, TIndexOffU j) const;
	bool isCovered(TIndexOffU i) const;
	TIndexOffU rank(TIndexOffU i) const;
	void print(ostream &out)
	{
		for (size_t i = 0; i < _text.length(); i++)
		{
			if (isCovered(i))
			{
				out << rank(i);
			}
			else
			{
				out << "-";
			}
			if (i < _text.length() - 1)
			{
				out << ",";
			}
		}
		out << endl;
	}
private:
	void doBuiltSanityCheck() const;
	void buildSPrime(EList<TIndexOffU> &sPrime, size_t padding);
	bool built() const
	{
		return _isaPrime.size() > 0;
	}
	void verbose(const string &s) const
	{
		if (this->verbose())
		{
			this->log() << s.c_str();
			this->log().flush();
		}
	}
	const TStr &_text;
	uint32_t _v;
	bool _verbose;
	bool _sanity;
	EList<uint32_t> _ds;
	EList<uint32_t> _dmap;
	uint32_t _d;
	EList<TIndexOffU> _doffs;
	EList<TIndexOffU> _isaPrime;
	EList<uint32_t> _dInv;
	uint32_t _log2v;
	TIndexOffU _vmask;
	ostream &_logger;
};
template <typename TStr>
void DifferenceCoverSample<TStr>::doBuiltSanityCheck() const
{
	uint32_t v = this->v();
	assert(built());
	VMSG_NL("  Doing sanity check");
	TIndexOffU added = 0;
	EList<TIndexOffU> sorted;
	sorted.resizeExact(_isaPrime.size());
	sorted.fill(OFF_MASK);
	for (size_t di = 0; di < this->d(); di++)
	{
		uint32_t d = _ds[di];
		size_t i = 0;
		for (size_t doi = _doffs[di]; doi < _doffs[di + 1]; doi++, i++)
		{
			assert_eq(OFF_MASK, sorted[_isaPrime[doi]]);
			sorted[_isaPrime[doi]] = (TIndexOffU)(v * i + d);
			added++;
		}
	}
	assert_eq(added, _isaPrime.size());
#ifndef NDEBUG
	for (size_t i = 0; i < sorted.size() - 1; i++)
	{
		assert(sstr_suf_lt(this->text(), sorted[i], this->text(), sorted[i + 1], false));
	}
#endif
}
template <typename TStr>
void DifferenceCoverSample<TStr>::buildSPrime(
	EList<TIndexOffU> &sPrime,
	size_t padding)
{
	const TStr &t = this->text();
	const EList<uint32_t> &ds = this->ds();
	TIndexOffU tlen = (TIndexOffU)t.length();
	uint32_t v = this->v();
	uint32_t d = this->d();
	assert_gt(v, 2);
	assert_lt(d, v);
	TIndexOffU tlenDivV = this->divv(tlen);
	uint32_t tlenModV = this->modv(tlen);
	TIndexOffU sPrimeSz = 0;
	assert(_doffs.empty());
	_doffs.resizeExact((size_t)d + 1);
	for (uint32_t di = 0; di < d; di++)
	{
		TIndexOffU sz = tlenDivV + ((ds[di] <= tlenModV) ? 1 : 0);
		assert_geq(sz, 0);
		_doffs[di] = sPrimeSz;
		sPrimeSz += sz;
	}
	_doffs[d] = sPrimeSz;
#ifndef NDEBUG
	if (tlenDivV > 0)
	{
		for (size_t i = 0; i < d; i++)
		{
			assert_gt(_doffs[i + 1], _doffs[i]);
			TIndexOffU diff = _doffs[i + 1] - _doffs[i];
			assert(diff == tlenDivV || diff == tlenDivV + 1);
		}
	}
#endif
	assert_eq(_doffs.size(), d + 1);
	sPrime.resizeExact((size_t)sPrimeSz + padding);
	sPrime.fill(OFF_MASK);
	TIndexOffU added = 0;
	TIndexOffU i = 0;
	for (TIndexOffU ti = 0; ti <= tlen; ti += v)
	{
		for (uint32_t di = 0; di < d; di++)
		{
			TIndexOffU tti = ti + ds[di];
			if (tti > tlen)
				break;
			TIndexOffU spi = _doffs[di] + i;
			assert_lt(spi, _doffs[di + 1]);
			assert_leq(tti, tlen);
			assert_lt(spi, sPrimeSz);
			assert_eq(OFF_MASK, sPrime[spi]);
			sPrime[spi] = tti;
			added++;
		}
		i++;
	}
	assert_eq(added, sPrimeSz);
}
template <typename TStr>
static inline bool suffixSameUpTo(
	const TStr &host,
	TIndexOffU suf1,
	TIndexOffU suf2,
	TIndexOffU v)
{
	for (TIndexOffU i = 0; i < v; i++)
	{
		bool endSuf1 = suf1 + i >= host.length();
		bool endSuf2 = suf2 + i >= host.length();
		if ((endSuf1 && !endSuf2) || (!endSuf1 && endSuf2))
			return false;
		if (endSuf1 && endSuf2)
			return true;
		if (host[suf1 + i] != host[suf2 + i])
			return false;
	}
	return true;
}
template <typename TStr>
struct VSortingParam
{
	DifferenceCoverSample<TStr> *dcs;
	TIndexOffU *sPrimeArr;
	size_t sPrimeSz;
	TIndexOffU *sPrimeOrderArr;
	size_t depth;
	const EList<size_t> *boundaries;
	size_t *cur;
	MUTEX_T *mutex;
};
template <typename TStr>
#ifdef WITH_TBB
class VSorting_worker
{
	void *vp;
public:
	VSorting_worker(const VSorting_worker &W) : vp(W.vp){};
	VSorting_worker(void *vp_) : vp(vp_){};
	void operator()() const {
#else
static void VSorting_worker(void *vp)
{
#endif
		VSortingParam<TStr> *param = (VSortingParam<TStr> *)vp;
	DifferenceCoverSample<TStr> *dcs = param->dcs;
	const TStr &host = dcs->text();
	const size_t hlen = host.length();
	uint32_t v = dcs->v();
	while (true)
	{
		size_t cur = 0;
		{
			ThreadSafe ts(*param->mutex);
			cur = *(param->cur);
			(*param->cur)++;
		}
		if (cur >= param->boundaries->size())
			return;
		size_t begin = (cur == 0 ? 0 : (*param->boundaries)[cur - 1]);
		size_t end = (*param->boundaries)[cur];
		assert_leq(begin, end);
		if (end - begin <= 1)
			continue;
		mkeyQSortSuf2(
			host,
			hlen,
			param->sPrimeArr,
			param->sPrimeSz,
			param->sPrimeOrderArr,
			4,
			begin,
			end,
			param->depth,
			v);
	}
}
#ifdef WITH_TBB
}
;
#endif
template <typename TStr>
void DifferenceCoverSample<TStr>::build(int nthreads)
{
	VMSG_NL("Building DifferenceCoverSample");
	const TStr &t = this->text();
	uint32_t v = this->v();
	assert_gt(v, 2);
	EList<TIndexOffU> sPrime;
	size_t padding = 1;
	VMSG_NL("  Building sPrime");
	buildSPrime(sPrime, padding);
	size_t sPrimeSz = sPrime.size() - padding;
	assert_gt(sPrime.size(), padding);
	assert_leq(sPrime.size(), t.length() + padding + 1);
	TIndexOffU nextRank = 0;
	{
		VMSG_NL("  Building sPrimeOrder");
		EList<TIndexOffU> sPrimeOrder;
		sPrimeOrder.resizeExact(sPrimeSz);
		for (TIndexOffU i = 0; i < sPrimeSz; i++)
		{
			sPrimeOrder[i] = i;
		}
		{
			Timer timer(cout, "  V-Sorting samples time: ", this->verbose());
			VMSG_NL("  V-Sorting samples");
			TIndexOffU *sPrimeArr = (TIndexOffU *)sPrime.ptr();
			assert_eq(sPrimeArr[0], sPrime[0]);
			assert_eq(sPrimeArr[sPrimeSz - 1], sPrime[sPrimeSz - 1]);
			TIndexOffU *sPrimeOrderArr = (TIndexOffU *)sPrimeOrder.ptr();
			assert_eq(sPrimeOrderArr[0], sPrimeOrder[0]);
			assert_eq(sPrimeOrderArr[sPrimeSz - 1], sPrimeOrder[sPrimeSz - 1]);
			if (nthreads == 1)
			{
				mkeyQSortSuf2(t, sPrimeArr, sPrimeSz, sPrimeOrderArr, 4,
							  this->verbose(), this->sanityCheck(), v);
			}
			else
			{
				int query_depth = 0;
				int tmp_nthreads = nthreads;
				while (tmp_nthreads > 0)
				{
					query_depth++;
					tmp_nthreads >>= 1;
				}
				EList<size_t> boundaries;
				TIndexOffU *sOrig = NULL;
				if (this->sanityCheck())
				{
					sOrig = new TIndexOffU[sPrimeSz];
					memcpy(sOrig, sPrimeArr, OFF_SIZE * sPrimeSz);
				}
				mkeyQSortSuf2(t, sPrimeArr, sPrimeSz, sPrimeOrderArr, 4,
							  this->verbose(), false, query_depth, &boundaries);
				if (boundaries.size() > 0)
				{
#ifdef WITH_TBB
					tbb::task_group tbb_grp;
#else
					AutoArray<tthread::thread *> threads(nthreads);
#endif
					EList<VSortingParam<TStr>> tparams;
					size_t cur = 0;
					MUTEX_T mutex;
					tparams.resize(nthreads);
					for (int tid = 0; tid < nthreads; tid++)
					{
						tparams[tid].dcs = this;
						tparams[tid].sPrimeArr = sPrimeArr;
						tparams[tid].sPrimeSz = sPrimeSz;
						tparams[tid].sPrimeOrderArr = sPrimeOrderArr;
						tparams[tid].depth = query_depth;
						tparams[tid].boundaries = &boundaries;
						tparams[tid].cur = &cur;
						tparams[tid].mutex = &mutex;
#ifdef WITH_TBB
						tbb_grp.run(VSorting_worker<TStr>(((void *)&tparams[tid])));
					}
					tbb_grp.wait();
#else
						threads[tid] = new tthread::thread(VSorting_worker<TStr>, (void *)&tparams[tid]);
					}
					for (int tid = 0; tid < nthreads; tid++)
					{
						threads[tid]->join();
					}
#endif
				}
				if (this->sanityCheck())
				{
					sanityCheckOrderedSufs(t, t.length(), sPrimeArr, sPrimeSz, v);
					for (size_t i = 0; i < sPrimeSz; i++)
					{
						assert_eq(sPrimeArr[i], sOrig[sPrimeOrderArr[i]]);
					}
					delete[] sOrig;
				}
			}
			assert_eq(sPrimeArr[0], sPrime[0]);
			assert_eq(sPrimeArr[sPrimeSz - 1], sPrime[sPrimeSz - 1]);
			assert_eq(sPrimeOrderArr[0], sPrimeOrder[0]);
			assert_eq(sPrimeOrderArr[sPrimeSz - 1], sPrimeOrder[sPrimeSz - 1]);
		}
		VMSG_NL("  Allocating rank array");
		_isaPrime.resizeExact(sPrime.size());
		ASSERT_ONLY(_isaPrime.fill(OFF_MASK));
		assert_gt(_isaPrime.size(), 0);
		{
			Timer timer(cout, "  Ranking v-sort output time: ", this->verbose());
			VMSG_NL("  Ranking v-sort output");
			for (size_t i = 0; i < sPrimeSz - 1; i++)
			{
				_isaPrime[sPrimeOrder[i]] = nextRank;
				if (!suffixSameUpTo(t, sPrime[i], sPrime[i + 1], v))
					nextRank++;
			}
			_isaPrime[sPrimeOrder[sPrimeSz - 1]] = nextRank;
#ifndef NDEBUG
			for (size_t i = 0; i < sPrimeSz; i++)
			{
				assert_neq(OFF_MASK, _isaPrime[i]);
				assert_lt(_isaPrime[i], sPrimeSz);
			}
#endif
		}
	}
	_isaPrime[_isaPrime.size() - 1] = (TIndexOffU)sPrimeSz;
	sPrime[sPrime.size() - 1] = (TIndexOffU)sPrimeSz;
	{
		Timer timer(cout, "  Invoking Larsson-Sadakane on ranks time: ", this->verbose());
		VMSG_NL("  Invoking Larsson-Sadakane on ranks");
		if (sPrime.size() >= LS_SIZE)
		{
			cerr << "Error; sPrime array has so many elements that it can't be converted to a signed array without overflow." << endl;
			throw 1;
		}
		LarssonSadakane<TIndexOff> ls;
		ls.suffixsort(
			(TIndexOff *)_isaPrime.ptr(),
			(TIndexOff *)sPrime.ptr(),
			(TIndexOff)sPrimeSz,
			(TIndexOff)sPrime.size(),
			0);
	}
	_isaPrime.resizeExact(sPrimeSz);
	for (size_t i = 0; i < _isaPrime.size(); i++)
	{
		_isaPrime[i]--;
	}
#ifndef NDEBUG
	for (size_t i = 0; i < sPrimeSz - 1; i++)
	{
		assert_lt(_isaPrime[i], sPrimeSz);
		assert(i == 0 || _isaPrime[i] != _isaPrime[i - 1]);
	}
#endif
	VMSG_NL("  Sanity-checking and returning");
	if (this->sanityCheck())
		doBuiltSanityCheck();
}
template <typename TStr>
bool DifferenceCoverSample<TStr>::isCovered(TIndexOffU i) const
{
	assert(built());
	uint32_t modi = this->modv(i);
	assert_lt(modi, _dInv.size());
	return _dInv[modi] != 0xffffffff;
}
template <typename TStr>
TIndexOffU DifferenceCoverSample<TStr>::rank(TIndexOffU i) const
{
	assert(built());
	assert_lt(i, this->text().length());
	uint32_t imodv = this->modv(i);
	assert_neq(0xffffffff, _dInv[imodv]);
	TIndexOffU ioff = this->divv(i);
	assert_lt(ioff, _doffs[_dInv[imodv] + 1] - _doffs[_dInv[imodv]]);
	TIndexOffU isaIIdx = _doffs[_dInv[imodv]] + ioff;
	assert_lt(isaIIdx, _isaPrime.size());
	TIndexOffU isaPrimeI = _isaPrime[isaIIdx];
	assert_leq(isaPrimeI, _isaPrime.size());
	return isaPrimeI;
}
template <typename TStr>
int64_t DifferenceCoverSample<TStr>::breakTie(TIndexOffU i, TIndexOffU j) const
{
	assert(built());
	assert_neq(i, j);
	assert_lt(i, this->text().length());
	assert_lt(j, this->text().length());
	uint32_t imodv = this->modv(i);
	uint32_t jmodv = this->modv(j);
	assert_neq(0xffffffff, _dInv[imodv]);
	assert_neq(0xffffffff, _dInv[jmodv]);
	uint32_t dimodv = _dInv[imodv];
	uint32_t djmodv = _dInv[jmodv];
	TIndexOffU ioff = this->divv(i);
	TIndexOffU joff = this->divv(j);
	assert_lt(dimodv + 1, _doffs.size());
	assert_lt(djmodv + 1, _doffs.size());
	assert_lt(ioff, _doffs[dimodv + 1] - _doffs[dimodv]);
	assert_lt(joff, _doffs[djmodv + 1] - _doffs[djmodv]);
	TIndexOffU isaIIdx = _doffs[dimodv] + ioff;
	TIndexOffU isaJIdx = _doffs[djmodv] + joff;
	assert_lt(isaIIdx, _isaPrime.size());
	assert_lt(isaJIdx, _isaPrime.size());
	assert_neq(isaIIdx, isaJIdx);
	TIndexOffU isaPrimeI = _isaPrime[isaIIdx];
	TIndexOffU isaPrimeJ = _isaPrime[isaJIdx];
	assert_neq(isaPrimeI, isaPrimeJ);
	assert_leq(isaPrimeI, _isaPrime.size());
	assert_leq(isaPrimeJ, _isaPrime.size());
	return (int64_t)isaPrimeI - (int64_t)isaPrimeJ;
}
template <typename TStr>
uint32_t DifferenceCoverSample<TStr>::tieBreakOff(TIndexOffU i, TIndexOffU j) const
{
	const TStr &t = this->text();
	const EList<uint32_t> &dmap = this->dmap();
	assert(built());
	if (t[i] != t[j])
		return 0xffffffff;
	uint32_t v = this->v();
	assert_neq(i, j);
	assert_lt(i, t.length());
	assert_lt(j, t.length());
	uint32_t imod = this->modv(i);
	uint32_t jmod = this->modv(j);
	uint32_t diffLeft = (jmod >= imod) ? (jmod - imod) : (jmod + v - imod);
	uint32_t diffRight = (imod >= jmod) ? (imod - jmod) : (imod + v - jmod);
	assert_lt(diffLeft, dmap.size());
	assert_lt(diffRight, dmap.size());
	uint32_t destLeft = dmap[diffLeft];
	uint32_t destRight = dmap[diffRight];
	assert(isCovered(destLeft));
	assert(isCovered(destLeft + diffLeft));
	assert(isCovered(destRight));
	assert(isCovered(destRight + diffRight));
	assert_lt(destLeft, v);
	assert_lt(destRight, v);
	uint32_t deltaLeft = (destLeft >= imod) ? (destLeft - imod) : (destLeft + v - imod);
	if (deltaLeft == v)
		deltaLeft = 0;
	uint32_t deltaRight = (destRight >= jmod) ? (destRight - jmod) : (destRight + v - jmod);
	if (deltaRight == v)
		deltaRight = 0;
	assert_lt(deltaLeft, v);
	assert_lt(deltaRight, v);
	assert(isCovered(i + deltaLeft));
	assert(isCovered(j + deltaLeft));
	assert(isCovered(i + deltaRight));
	assert(isCovered(j + deltaRight));
	return min(deltaLeft, deltaRight);
}
#endif
