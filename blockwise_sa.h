#ifndef BLOCKWISE_SA_H_
#define BLOCKWISE_SA_H_
#ifdef WITH_TBB
#include <tbb/tbb.h>
#include <tbb/task_group.h>
#endif
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <memory>
#include <stdexcept>
#include "assert_helpers.h"
#include "diff_sample.h"
#include "multikey_qsort.h"
#include "random_source.h"
#include "binary_sa_search.h"
#include "zbox.h"
#include "alphabet.h"
#include "timer.h"
#include "ds.h"
#include "mem_ids.h"
#include "word_io.h"
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
template <typename TStr>
class BlockwiseSA
{
public:
    BlockwiseSA(const TStr &__text,
                TIndexOffU __bucketSz,
                bool __sanityCheck = false,
                bool __passMemExc = false,
                bool __verbose = false,
                ostream &__logger = cout) : _text(__text),
                                            _bucketSz(max<TIndexOffU>(__bucketSz, 2u)),
                                            _sanityCheck(__sanityCheck),
                                            _passMemExc(__passMemExc),
                                            _verbose(__verbose),
                                            _itrBucket(EBWTB_CAT),
                                            _itrBucketPos(OFF_MASK),
                                            _itrPushedBackSuffix(OFF_MASK),
                                            _logger(__logger)
    {
    }
    virtual ~BlockwiseSA() {}
    virtual TIndexOffU nextSuffix() = 0;
    bool hasMoreSuffixes()
    {
        if (_itrPushedBackSuffix != OFF_MASK)
            return true;
        try
        {
            _itrPushedBackSuffix = nextSuffix();
        }
        catch (out_of_range &e)
        {
            assert_eq(OFF_MASK, _itrPushedBackSuffix);
            return false;
        }
        return true;
    }
    void resetSuffixItr()
    {
        _itrBucket.clear();
        _itrBucketPos = OFF_MASK;
        _itrPushedBackSuffix = OFF_MASK;
        reset();
        assert(suffixItrIsReset());
    }
    bool suffixItrIsReset()
    {
        return _itrBucket.size() == 0 &&
               _itrBucketPos == OFF_MASK &&
               _itrPushedBackSuffix == OFF_MASK &&
               isReset();
    }
    const TStr &text() const { return _text; }
    TIndexOffU bucketSz() const { return _bucketSz; }
    bool sanityCheck() const { return _sanityCheck; }
    bool verbose() const { return _verbose; }
    ostream &log() const { return _logger; }
    size_t size() const { return _text.length() + 1; }
protected:
    virtual void reset() = 0;
    virtual bool isReset() = 0;
    virtual void nextBlock(int cur_block, int tid = 0) = 0;
    virtual bool hasMoreBlocks() const = 0;
    void verbose(const string &s) const
    {
        if (this->verbose())
        {
            this->log() << s.c_str();
            this->log().flush();
        }
    }
    const TStr &_text;
    const TIndexOffU _bucketSz;
    const bool _sanityCheck;
    const bool _passMemExc;
    const bool _verbose;
    EList<TIndexOffU> _itrBucket;
    TIndexOffU _itrBucketPos;
    TIndexOffU _itrPushedBackSuffix;
    ostream &_logger;
};
template <typename TStr>
class InorderBlockwiseSA : public BlockwiseSA<TStr>
{
public:
    InorderBlockwiseSA(const TStr &__text,
                       TIndexOffU __bucketSz,
                       bool __sanityCheck = false,
                       bool __passMemExc = false,
                       bool __verbose = false,
                       ostream &__logger = cout) : BlockwiseSA<TStr>(__text, __bucketSz, __sanityCheck, __passMemExc, __verbose, __logger)
    {
    }
};
template <typename TStr>
class KarkkainenBlockwiseSA : public InorderBlockwiseSA<TStr>
{
public:
    typedef DifferenceCoverSample<TStr> TDC;
    KarkkainenBlockwiseSA(const TStr &__text,
                          TIndexOffU __bucketSz,
                          int __nthreads,
                          uint32_t __dcV,
                          uint32_t __seed = 0,
                          bool __sanityCheck = false,
                          bool __passMemExc = false,
                          bool __verbose = false,
                          string base_fname = "",
                          ostream &__logger = cout) : InorderBlockwiseSA<TStr>(__text, __bucketSz, __sanityCheck, __passMemExc, __verbose, __logger),
                                                      _sampleSuffs(EBWTB_CAT),
                                                      _nthreads(__nthreads),
                                                      _itrBucketIdx(0),
                                                      _cur(0),
                                                      _dcV(__dcV),
                                                      _dc(EBWTB_CAT),
                                                      _built(false),
                                                      _base_fname(base_fname),
                                                      _bigEndian(currentlyBigEndian()),
#ifdef WITH_TBB
                                                      thread_group_started(false),
#endif
                                                      _done(NULL)
    {
        _randomSrc.init(__seed);
        reset();
    }
    ~KarkkainenBlockwiseSA() throw()
    {
#ifdef WITH_TBB
        tbb_grp.wait();
#else
        if (_threads.size() > 0)
        {
            for (size_t tid = 0; tid < _threads.size(); tid++)
            {
                _threads[tid]->join();
                delete _threads[tid];
            }
        }
#endif
        if (_done != NULL)
            delete[] _done;
    }
    static size_t simulateAllocs(const TStr &text, TIndexOffU bucketSz)
    {
        size_t len = text.length();
        size_t bsz = bucketSz;
        size_t sssz = len / max<TIndexOffU>(bucketSz - 1, 1);
        AutoArray<TIndexOffU> tmp(bsz + sssz + (1024 * 1024), EBWT_CAT);
        return bsz;
    }
    virtual TIndexOffU nextSuffix()
    {
        if (this->_nthreads > 1)
        {
#ifdef WITH_TBB
            if (!thread_group_started)
            {
#else
            if (_threads.size() == 0)
            {
#endif
                _done = new volatile bool[_sampleSuffs.size() + 1];
                for (size_t i = 0; i < _sampleSuffs.size() + 1; i++)
                {
                    _done[i] = false;
                }
                _itrBuckets.resize(this->_nthreads);
                _tparams.resize(this->_nthreads);
                for (int tid = 0; tid < this->_nthreads; tid++)
                {
                    _tparams[tid].first = this;
                    _tparams[tid].second = tid;
#ifdef WITH_TBB
                    tbb_grp.run(nextBlock_Worker((void *)&_tparams[tid]));
                }
                thread_group_started = true;
#else
                    _threads.push_back(new tthread::thread(nextBlock_Worker, (void *)&_tparams[tid]));
                }
                assert_eq(_threads.size(), (size_t)this->_nthreads);
#endif
            }
        }
        if (this->_itrPushedBackSuffix != OFF_MASK)
        {
            TIndexOffU tmp = this->_itrPushedBackSuffix;
            this->_itrPushedBackSuffix = OFF_MASK;
            return tmp;
        }
        while (this->_itrBucketPos >= this->_itrBucket.size() ||
               this->_itrBucket.size() == 0)
        {
            if (!hasMoreBlocks())
            {
                throw out_of_range("No more suffixes");
            }
            if (this->_nthreads == 1)
            {
                nextBlock((int)_cur);
                _cur++;
            }
            else
            {
                while (!_done[this->_itrBucketIdx])
                {
                    SLEEP(1);
                }
                std::ostringstream number;
                number << this->_itrBucketIdx;
                const string fname = _base_fname + "." + number.str() + ".sa";
                ifstream sa_file(fname.c_str(), ios::binary);
                if (!sa_file.good())
                {
                    cerr << "Could not open file for reading a suffix array: \"" << fname << "\"" << endl;
                    throw 1;
                }
                size_t numSAs = readU<TIndexOffU>(sa_file, false);
                this->_itrBucket.resizeExact(numSAs);
                for (size_t i = 0; i < numSAs; i++)
                {
                    this->_itrBucket[i] = readU<TIndexOffU>(sa_file, false);
                }
                sa_file.close();
                std::remove(fname.c_str());
            }
            this->_itrBucketIdx++;
            this->_itrBucketPos = 0;
        }
        return this->_itrBucket[this->_itrBucketPos++];
    }
    virtual void nextBlock(int cur_block, int tid = 0);
    virtual void qsort(EList<TIndexOffU> &bucket);
    virtual bool hasMoreBlocks() const
    {
        return this->_itrBucketIdx <= _sampleSuffs.size();
    }
    uint32_t dcV() const { return _dcV; }
#ifdef WITH_TBB
    class nextBlock_Worker
    {
        void *vp;
    public:
        nextBlock_Worker(const nextBlock_Worker &W) : vp(W.vp){};
        nextBlock_Worker(void *vp_) : vp(vp_){};
        void operator()() const {
#else
    static void nextBlock_Worker(void *vp)
    {
#endif
            pair<KarkkainenBlockwiseSA *, int> param = *(pair<KarkkainenBlockwiseSA *, int> *)vp;
        KarkkainenBlockwiseSA *sa = param.first;
        int tid = param.second;
        while (true)
        {
            size_t cur = 0;
            {
                ThreadSafe ts(sa->_mutex);
                cur = sa->_cur;
                if (cur > sa->_sampleSuffs.size())
                    break;
                sa->_cur++;
            }
            sa->nextBlock((int)cur, tid);
            std::ostringstream number;
            number << cur;
            const string fname = sa->_base_fname + "." + number.str() + ".sa";
            ofstream sa_file(fname.c_str(), ios::binary);
            if (!sa_file.good())
            {
                cerr << "Could not open file for writing a reference graph: \"" << fname << "\"" << endl;
                throw 1;
            }
            const EList<TIndexOffU> &bucket = sa->_itrBuckets[tid];
            writeU<TIndexOffU>(sa_file, (TIndexOffU)bucket.size(), sa->_bigEndian);
            for (size_t i = 0; i < bucket.size(); i++)
            {
                writeU<TIndexOffU>(sa_file, bucket[i], sa->_bigEndian);
            }
            sa_file.close();
            sa->_itrBuckets[tid].clear();
            sa->_done[cur] = true;
        }
    }
#ifdef WITH_TBB
};
#endif
protected:
virtual void reset()
{
    if (!_built)
    {
        build();
    }
    assert(_built);
    _cur = 0;
}
virtual bool isReset()
{
    return _cur == 0;
}
private:
void build()
{
    assert(_dc.get() == NULL);
    if (_dcV != 0)
    {
        _dc.init(new TDC(this->text(), _dcV, this->verbose(), this->sanityCheck()));
        _dc.get()->build(this->_nthreads);
    }
    if (this->bucketSz() <= this->text().length())
    {
        VMSG_NL("Building samples");
        buildSamples();
    }
    else
    {
        VMSG_NL("Skipping building samples since text length " << this->text().length() << " is less than bucket size: " << this->bucketSz());
    }
    _built = true;
}
inline bool tieBreakingLcp(TIndexOffU aOff,
                           TIndexOffU bOff,
                           TIndexOffU &lcp,
                           bool &lcpIsSoft);
inline bool suffixCmp(TIndexOffU cmp,
                      TIndexOffU i,
                      int64_t &j,
                      int64_t &k,
                      bool &kSoft,
                      const EList<TIndexOffU> &z);
void buildSamples();
EList<TIndexOffU> _sampleSuffs;
int _nthreads;
TIndexOffU _itrBucketIdx;
TIndexOffU _cur;
const uint32_t _dcV;
PtrWrap<TDC> _dc;
bool _built;
RandomSource _randomSrc;
MUTEX_T _mutex;
string _base_fname;
bool _bigEndian;
#ifdef WITH_TBB
tbb::task_group tbb_grp;
bool thread_group_started;
#else
    EList<tthread::thread *> _threads;
#endif
EList<pair<KarkkainenBlockwiseSA *, int>> _tparams;
ELList<TIndexOffU> _itrBuckets;
volatile bool *_done;
}
;
template <typename TStr>
inline void KarkkainenBlockwiseSA<TStr>::qsort(EList<TIndexOffU> &bucket)
{
    const TStr &t = this->text();
    TIndexOffU *s = bucket.ptr();
    size_t slen = bucket.size();
    TIndexOffU len = (TIndexOffU)t.length();
    if (_dc.get() != NULL)
    {
        VMSG_NL("  (Using difference cover)");
        const uint8_t *host = (const uint8_t *)t.buf();
        assert(_dc.get() != NULL);
        mkeyQSortSufDcU8(t, host, len, s, slen, *_dc.get(), 4,
                         this->verbose(), this->sanityCheck());
    }
    else
    {
        VMSG_NL("  (Not using difference cover)");
        mkeyQSortSuf(t, s, slen, 4,
                     this->verbose(), this->sanityCheck());
    }
}
template <>
inline void KarkkainenBlockwiseSA<S2bDnaString>::qsort(
    EList<TIndexOffU> &bucket)
{
    const S2bDnaString &t = this->text();
    TIndexOffU *s = bucket.ptr();
    size_t slen = bucket.size();
    size_t len = t.length();
    if (_dc.get() != NULL)
    {
        VMSG_NL("  (Using difference cover)");
        mkeyQSortSufDcU8(t, t, len, s, slen, *_dc.get(), 4,
                         this->verbose(), this->sanityCheck());
    }
    else
    {
        VMSG_NL("  (Not using difference cover)");
        mkeyQSortSuf(t, s, slen, 4,
                     this->verbose(), this->sanityCheck());
    }
}
template <typename TStr>
struct BinarySortingParam
{
    const TStr *t;
    const EList<TIndexOffU> *sampleSuffs;
    EList<TIndexOffU> bucketSzs;
    EList<TIndexOffU> bucketReps;
    size_t begin;
    size_t end;
};
template <typename TStr>
#ifdef WITH_TBB
class BinarySorting_worker
{
    void *vp;
public:
    BinarySorting_worker(const BinarySorting_worker &W) : vp(W.vp){};
    BinarySorting_worker(void *vp_) : vp(vp_){};
    void operator()() const {
#else
static void BinarySorting_worker(void *vp)
{
#endif
        BinarySortingParam<TStr> *param = (BinarySortingParam<TStr> *)vp;
    const TStr &t = *(param->t);
    size_t len = t.length();
    const EList<TIndexOffU> &sampleSuffs = *(param->sampleSuffs);
    EList<TIndexOffU> &bucketSzs = param->bucketSzs;
    EList<TIndexOffU> &bucketReps = param->bucketReps;
    ASSERT_ONLY(size_t numBuckets = bucketSzs.size());
    size_t begin = param->begin;
    size_t end = param->end;
    for (TIndexOffU i = (TIndexOffU)begin; i < end && i < len; i++)
    {
        TIndexOffU r = binarySASearch(t, i, sampleSuffs);
        if (r == std::numeric_limits<TIndexOffU>::max())
            continue;
        assert_lt(r, numBuckets);
        bucketSzs[r]++;
        assert_lt(bucketSzs[r], len);
        if (bucketReps[r] == OFF_MASK || (i & 100) == 0)
        {
            bucketReps[r] = i;
        }
    }
}
#ifdef WITH_TBB
}
;
#endif
template <typename TStr>
void KarkkainenBlockwiseSA<TStr>::buildSamples()
{
    const TStr &t = this->text();
    TIndexOffU bsz = this->bucketSz() - 1;
    size_t len = this->text().length();
    _sampleSuffs.clear();
    TIndexOffU numSamples = (TIndexOffU)((len / bsz) + 1) << 1;
    assert_gt(numSamples, 0);
    VMSG_NL("Reserving space for " << numSamples << " sample suffixes");
    if (this->_passMemExc)
    {
        _sampleSuffs.resizeExact(numSamples);
        VMSG_NL("Generating random suffixes");
        for (size_t i = 0; i < numSamples; i++)
        {
#ifdef _64BIT_INDEX
            _sampleSuffs[i] = (TIndexOffU)(_randomSrc.nextU64() % len);
#else
            _sampleSuffs[i] = (TIndexOffU)(_randomSrc.nextU32() % len);
#endif
        }
    }
    else
    {
        try
        {
            _sampleSuffs.resizeExact(numSamples);
            VMSG_NL("Generating random suffixes");
            for (size_t i = 0; i < numSamples; i++)
            {
#ifdef _64BIT_INDEX
                _sampleSuffs[i] = (TIndexOffU)(_randomSrc.nextU64() % len);
#else
                _sampleSuffs[i] = (TIndexOffU)(_randomSrc.nextU32() % len);
#endif
            }
        }
        catch (bad_alloc &e)
        {
            if (this->_passMemExc)
            {
                throw e;
            }
            else
            {
                cerr << "Could not allocate sample suffix container of " << (numSamples * OFF_SIZE) << " bytes." << endl
                     << "Please try using a smaller number of blocks by specifying a larger --bmax or" << endl
                     << "a smaller --bmaxdivn" << endl;
                throw 1;
            }
        }
    }
    {
        Timer timer(cout, "QSorting sample offsets, eliminating duplicates time: ", this->verbose());
        VMSG_NL("QSorting " << _sampleSuffs.size() << " sample offsets, eliminating duplicates");
        _sampleSuffs.sort();
        size_t sslen = _sampleSuffs.size();
        for (size_t i = 0; i < sslen - 1; i++)
        {
            if (_sampleSuffs[i] == _sampleSuffs[i + 1])
            {
                _sampleSuffs.erase(i--);
                sslen--;
            }
        }
    }
    {
        Timer timer(cout, "  Multikey QSorting samples time: ", this->verbose());
        VMSG_NL("Multikey QSorting " << _sampleSuffs.size() << " samples");
        this->qsort(_sampleSuffs);
    }
    VMSG_NL("Calculating bucket sizes");
    int limit = 5;
    while (--limit >= 0)
    {
        TIndexOffU numBuckets = (TIndexOffU)_sampleSuffs.size() + 1;
#ifdef WITH_TBB
        tbb::task_group tbb_grp;
#else
        AutoArray<tthread::thread *> threads(this->_nthreads);
#endif
        EList<BinarySortingParam<TStr>> tparams;
        tparams.resize(this->_nthreads);
        for (int tid = 0; tid < this->_nthreads; tid++)
        {
            try
            {
                tparams[tid].bucketSzs.resizeExact(numBuckets);
                tparams[tid].bucketReps.resizeExact(numBuckets);
                tparams[tid].bucketSzs.fillZero();
                tparams[tid].bucketReps.fill(OFF_MASK);
            }
            catch (bad_alloc &e)
            {
                if (this->_passMemExc)
                {
                    throw e;
                }
                else
                {
                    cerr << "Could not allocate sizes, representatives (" << ((numBuckets * 8) >> 10) << " KB) for blocks." << endl
                         << "Please try using a smaller number of blocks by specifying a larger --bmax or a" << endl
                         << "smaller --bmaxdivn." << endl;
                    throw 1;
                }
            }
            tparams[tid].t = &t;
            tparams[tid].sampleSuffs = &_sampleSuffs;
            tparams[tid].begin = (tid == 0 ? 0 : len / this->_nthreads * tid);
            tparams[tid].end = (tid + 1 == this->_nthreads ? len : len / this->_nthreads * (tid + 1));
            if (this->_nthreads == 1)
            {
                BinarySorting_worker<TStr>((void *)&tparams[tid]);
            }
            else
            {
#ifdef WITH_TBB
                tbb_grp.run(BinarySorting_worker<TStr>(((void *)&tparams[tid])));
            }
        }
        tbb_grp.wait();
#else
                threads[tid] = new tthread::thread(BinarySorting_worker<TStr>, (void *)&tparams[tid]);
            }
        }
        if (this->_nthreads > 1)
        {
            for (int tid = 0; tid < this->_nthreads; tid++)
            {
                threads[tid]->join();
            }
        }
#endif
        EList<TIndexOffU> &bucketSzs = tparams[0].bucketSzs;
        EList<TIndexOffU> &bucketReps = tparams[0].bucketReps;
        for (int tid = 1; tid < this->_nthreads; tid++)
        {
            for (size_t j = 0; j < numBuckets; j++)
            {
                bucketSzs[j] += tparams[tid].bucketSzs[j];
                if (bucketReps[j] == OFF_MASK)
                {
                    bucketReps[j] = tparams[tid].bucketReps[j];
                }
            }
        }
        TIndexOff added = 0;
        TIndexOff merged = 0;
        assert_eq(bucketSzs.size(), numBuckets);
        assert_eq(bucketReps.size(), numBuckets);
        {
            Timer timer(cout, "  Splitting and merging time: ", this->verbose());
            VMSG_NL("Splitting and merging");
            for (TIndexOffU i = 0; i < numBuckets; i++)
            {
                TIndexOffU mergedSz = bsz + 1;
                assert(bucketSzs[(size_t)i] == 0 || bucketReps[(size_t)i] != OFF_MASK);
                if (i < numBuckets - 1)
                {
                    mergedSz = bucketSzs[(size_t)i] + bucketSzs[(size_t)i + 1] + 1;
                }
                if (mergedSz <= bsz)
                {
                    bucketSzs[(size_t)i + 1] += (bucketSzs[(size_t)i] + 1);
                    bucketReps[(size_t)i + 1] = _sampleSuffs[(size_t)i + added];
                    _sampleSuffs.erase((size_t)i + added);
                    bucketSzs.erase((size_t)i);
                    bucketReps.erase((size_t)i);
                    i--;
                    numBuckets--;
                    merged++;
                    assert_eq(numBuckets, _sampleSuffs.size() + 1 - added);
                    assert_eq(numBuckets, bucketSzs.size());
                }
                else if (bucketSzs[(size_t)i] > bsz)
                {
                    _sampleSuffs.insert(bucketReps[(size_t)i], (TIndexOffU)(i + (added++)));
                }
            }
        }
        if (added == 0)
        {
            break;
        }
        VMSG_NL("Split " << added << ", merged " << merged << "; iterating...");
    }
    VMSG_NL("Avg bucket size: " << ((double)(len - _sampleSuffs.size()) / (_sampleSuffs.size() + 1)) << " (target: " << bsz << ")");
}
template <typename T>
inline static TIndexOffU suffixLcp(const T &t, TIndexOffU aOff, TIndexOffU bOff)
{
    TIndexOffU c = 0;
    size_t len = t.length();
    assert_leq(aOff, len);
    assert_leq(bOff, len);
    while (aOff + c < len && bOff + c < len && t[aOff + c] == t[bOff + c])
        c++;
    return c;
}
template <typename TStr>
inline bool KarkkainenBlockwiseSA<TStr>::tieBreakingLcp(TIndexOffU aOff,
                                                        TIndexOffU bOff,
                                                        TIndexOffU &lcp,
                                                        bool &lcpIsSoft)
{
    const TStr &t = this->text();
    TIndexOffU c = 0;
    TIndexOffU tlen = (TIndexOffU)t.length();
    assert_leq(aOff, tlen);
    assert_leq(bOff, tlen);
    assert(_dc.get() != NULL);
    uint32_t dcDist = _dc.get()->tieBreakOff(aOff, bOff);
    lcpIsSoft = false;
    while (c < dcDist && c < tlen - aOff && c < tlen - bOff && t[aOff + c] == t[bOff + c])
        c++;
    lcp = c;
    if (c == tlen - aOff)
    {
        return false;
    }
    else if (c == tlen - bOff)
    {
        return true;
    }
    else if (c == dcDist)
    {
        lcpIsSoft = true;
        assert_neq(dcDist, 0xffffffff);
        return _dc.get()->breakTie(aOff + c, bOff + c) < 0;
    }
    else
    {
        assert_neq(t[aOff + c], t[bOff + c]);
        return t[aOff + c] < t[bOff + c];
    }
}
template <typename T>
static TIndexOffU lookupSuffixZ(
    const T &t,
    TIndexOffU zOff,
    TIndexOffU off,
    const EList<TIndexOffU> &z)
{
    if (zOff < z.size())
    {
        TIndexOffU ret = z[zOff];
        assert_eq(ret, suffixLcp(t, off + zOff, off));
        return ret;
    }
    assert_leq(off + zOff, t.length());
    return suffixLcp(t, off + zOff, off);
}
template <typename TStr>
inline bool KarkkainenBlockwiseSA<TStr>::suffixCmp(
    TIndexOffU cmp,
    TIndexOffU i,
    int64_t &j,
    int64_t &k,
    bool &kSoft,
    const EList<TIndexOffU> &z)
{
    const TStr &t = this->text();
    TIndexOffU len = (TIndexOffU)t.length();
    TIndexOffU l;
    if ((int64_t)i > k)
    {
        k = i;
        l = 0;
        kSoft = false;
    }
    else
    {
        assert_gt((int64_t)i, j);
        TIndexOffU zIdx = (TIndexOffU)(i - j);
        assert_leq(zIdx, len - cmp);
        if (zIdx < _dcV || _dc.get() == NULL)
        {
            l = lookupSuffixZ(t, zIdx, cmp, z);
            if (i + l > len)
            {
                l = len - i;
            }
            assert_leq(i + l, len);
        }
        else
        {
            bool ret = tieBreakingLcp(i, cmp, l, kSoft);
            if (this->sanityCheck())
            {
                if (ret)
                    assert(sstr_suf_lt(t, i, t, cmp, false));
                else
                    assert(sstr_suf_gt(t, i, t, cmp, false));
            }
            j = i;
            k = i + l;
            if (this->sanityCheck())
            {
                if (kSoft)
                {
                    assert_leq(l, suffixLcp(t, i, cmp));
                }
                else
                {
                    assert_eq(l, suffixLcp(t, i, cmp));
                }
            }
            return ret;
        }
    }
    if ((int64_t)(i + l) == k)
    {
        while (l < len - cmp && k < (int64_t)len && t[(size_t)(cmp + l)] == t[(size_t)k])
        {
            k++;
            l++;
        }
        j = i;
        kSoft = false;
        assert_eq(l, suffixLcp(t, i, cmp));
    }
    else if ((int64_t)(i + l) > k)
    {
        l = (TIndexOffU)(k - i);
        j = i;
        if (kSoft)
        {
            while (l < len - cmp && k < (int64_t)len && t[(size_t)(cmp + l)] == t[(size_t)k])
            {
                k++;
                l++;
            }
            kSoft = false;
            assert_eq(l, suffixLcp(t, i, cmp));
        }
        else
            assert_eq(l, suffixLcp(t, i, cmp));
    }
    if (this->sanityCheck())
    {
        if (!kSoft)
        {
            assert_eq(l, suffixLcp(t, i, cmp));
        }
        else
        {
            assert_leq(l, suffixLcp(t, i, cmp));
        }
    }
    assert_leq(l + i, len);
    assert_leq(l, len - cmp);
    assert(l != len - cmp || i + l != len);
    if (l + i != len && (l == len - cmp ||
                         t[i + l] < t[cmp + l]))
    {
#ifndef NDEBUG
        if (this->sanityCheck())
        {
            assert(sstr_suf_lt(t, i, t, cmp, false));
        }
#endif
        return true;
    }
    else
    {
#ifndef NDEBUG
        if (this->sanityCheck())
        {
            assert(sstr_suf_gt(t, i, t, cmp, false));
        }
#endif
        return false;
    }
}
template <typename TStr>
void KarkkainenBlockwiseSA<TStr>::nextBlock(int cur_block, int tid)
{
#ifndef NDEBUG
    if (this->_nthreads > 1)
    {
        assert_lt(tid, (int)this->_itrBuckets.size());
    }
#endif
    EList<TIndexOffU> &bucket = (this->_nthreads > 1 ? this->_itrBuckets[tid] : this->_itrBucket);
    {
        ThreadSafe ts(_mutex);
        VMSG_NL("Getting block " << (cur_block + 1) << " of " << _sampleSuffs.size() + 1);
    }
    assert(_built);
    assert_gt(_dcV, 3);
    assert_leq(cur_block, (int)_sampleSuffs.size());
    const TStr &t = this->text();
    TIndexOffU len = (TIndexOffU)t.length();
    bucket.clear();
    TIndexOffU lo = OFF_MASK, hi = OFF_MASK;
    if (_sampleSuffs.size() == 0)
    {
        {
            ThreadSafe ts(_mutex);
            VMSG_NL("  No samples; assembling all-inclusive block");
        }
        assert_eq(0, cur_block);
        try
        {
            if (bucket.capacity() < this->bucketSz())
            {
                bucket.reserveExact(len + 1);
            }
            bucket.resize(len);
            for (TIndexOffU i = 0; i < len; i++)
            {
                bucket[i] = i;
            }
        }
        catch (bad_alloc &e)
        {
            if (this->_passMemExc)
            {
                throw e;
            }
            else
            {
                cerr << "Could not allocate a master suffix-array block of " << ((len + 1) * 4) << " bytes" << endl
                     << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
                     << "a larger --bmaxdivn" << endl;
                throw 1;
            }
        }
    }
    else
    {
        try
        {
            {
                ThreadSafe ts(_mutex);
                VMSG_NL("  Reserving size (" << this->bucketSz() << ") for bucket " << (cur_block + 1));
            }
            if (bucket.size() < this->bucketSz() + 100)
            {
                bucket.reserveExact(this->bucketSz() + 100);
            }
        }
        catch (bad_alloc &e)
        {
            if (this->_passMemExc)
            {
                throw e;
            }
            else
            {
                cerr << "Could not allocate a suffix-array block of " << ((this->bucketSz() + 1) * 4) << " bytes" << endl;
                cerr << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
                     << "a larger --bmaxdivn" << endl;
                throw 1;
            }
        }
        EList<TIndexOffU> zLo(EBWTB_CAT), zHi(EBWTB_CAT);
        assert_geq(cur_block, 0);
        assert_leq((size_t)cur_block, _sampleSuffs.size());
        bool first = (cur_block == 0);
        bool last = ((size_t)cur_block == _sampleSuffs.size());
        try
        {
            {
                ThreadSafe ts(_mutex);
                VMSG_NL("  Calculating Z arrays for bucket " << (cur_block + 1));
            }
            if (!last)
            {
                assert_lt(cur_block, (int)_sampleSuffs.size());
                hi = _sampleSuffs[cur_block];
                zHi.resizeExact(_dcV);
                zHi.fillZero();
                assert_eq(zHi[0], 0);
                calcZ(t, hi, zHi, this->verbose(), this->sanityCheck());
            }
            if (!first)
            {
                assert_gt(cur_block, 0);
                assert_leq(cur_block, (int)_sampleSuffs.size());
                lo = _sampleSuffs[cur_block - 1];
                zLo.resizeExact(_dcV);
                zLo.fillZero();
                assert_gt(_dcV, 3);
                assert_eq(zLo[0], 0);
                calcZ(t, lo, zLo, this->verbose(), this->sanityCheck());
            }
        }
        catch (bad_alloc &e)
        {
            if (this->_passMemExc)
            {
                throw e;
            }
            else
            {
                cerr << "Could not allocate a z-array of " << (_dcV * 4) << " bytes" << endl;
                cerr << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
                     << "a larger --bmaxdivn" << endl;
                throw 1;
            }
        }
        int64_t kHi = -1, kLo = -1;
        int64_t jHi = -1, jLo = -1;
        bool kHiSoft = false, kLoSoft = false;
        assert_eq(0, bucket.size());
        {
            {
                ThreadSafe ts(_mutex);
                VMSG_NL("  Entering block accumulator loop for bucket " << (cur_block + 1) << ":");
            }
            TIndexOffU lenDiv10 = (len + 9) / 10;
            for (TIndexOffU iten = 0, ten = 0; iten < len; iten += lenDiv10, ten++)
            {
                TIndexOffU itenNext = iten + lenDiv10;
                {
                    ThreadSafe ts(_mutex);
                    if (ten > 0)
                        VMSG_NL("  bucket " << (cur_block + 1) << ": " << (ten * 10) << "%");
                }
                for (TIndexOffU i = iten; i < itenNext && i < len; i++)
                {
                    assert_lt(jLo, (int64_t)i);
                    assert_lt(jHi, (int64_t)i);
                    if (i == hi || i == lo)
                        continue;
                    if (hi != OFF_MASK && !suffixCmp(hi, i, jHi, kHi, kHiSoft, zHi))
                    {
                        continue;
                    }
                    if (lo != OFF_MASK && suffixCmp(lo, i, jLo, kLo, kLoSoft, zLo))
                    {
                        continue;
                    }
                    assert_lt(i, len);
                    try
                    {
                        bucket.push_back(i);
                    }
                    catch (bad_alloc &e)
                    {
                        cerr << "Could not append element to block of " << ((bucket.size()) * OFF_SIZE) << " bytes" << endl;
                        if (this->_passMemExc)
                        {
                            throw e;
                        }
                        else
                        {
                            cerr << "Please try using a larger number of blocks by specifying a smaller --bmax or" << endl
                                 << "a larger --bmaxdivn" << endl;
                            throw 1;
                        }
                    }
                }
            }
            {
                ThreadSafe ts(_mutex);
                VMSG_NL("  bucket " << (cur_block + 1) << ": 100%");
            }
        }
    }
    if (bucket.size() > 0)
    {
        Timer timer(cout, "  Sorting block time: ", this->verbose());
        {
            ThreadSafe ts(_mutex);
            VMSG_NL("  Sorting block of length " << bucket.size() << " for bucket " << (cur_block + 1));
        }
        this->qsort(bucket);
    }
    if (hi != OFF_MASK)
    {
        bucket.push_back(hi);
    }
    else
    {
        bucket.push_back(len);
    }
    {
        ThreadSafe ts(_mutex);
        VMSG_NL("Returning block of " << bucket.size() << " for bucket " << (cur_block + 1));
    }
}
#endif
