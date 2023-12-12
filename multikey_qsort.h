#ifndef MULTIKEY_QSORT_H_
#define MULTIKEY_QSORT_H_
#include <iostream>
#include "sequence_io.h"
#include "alphabet.h"
#include "assert_helpers.h"
#include "diff_sample.h"
#include "sstring.h"
#include "btypes.h"
using namespace std;
template <typename TStr, typename TPos>
static inline void swap(TStr &s, size_t slen, TPos a, TPos b)
{
	assert_lt(a, slen);
	assert_lt(b, slen);
	swap(s[a], s[b]);
}
template <typename TVal, typename TPos>
static inline void swap(TVal *s, size_t slen, TPos a, TPos b)
{
	assert_lt(a, slen);
	assert_lt(b, slen);
	swap(s[a], s[b]);
}
#define SWAP(s, a, b)         \
	{                         \
		assert_geq(a, begin); \
		assert_geq(b, begin); \
		assert_lt(a, end);    \
		assert_lt(b, end);    \
		swap(s, slen, a, b);  \
	}
#define SWAP2(s, s2, a, b)    \
	{                         \
		SWAP(s, a, b);        \
		swap(s2, slen, a, b); \
	}
#define SWAP1(s, s2, a, b) \
	{                      \
		SWAP(s, a, b);     \
	}
#define VECSWAP(s, i, j, n)                        \
	{                                              \
		if (n > 0)                                 \
		{                                          \
			vecswap(s, slen, i, j, n, begin, end); \
		}                                          \
	}
#define VECSWAP2(s, s2, i, j, n)                        \
	{                                                   \
		if (n > 0)                                      \
		{                                               \
			vecswap2(s, slen, s2, i, j, n, begin, end); \
		}                                               \
	}
template <typename TStr, typename TPos>
static inline void vecswap(TStr &s, size_t slen, TPos i, TPos j, TPos n, TPos begin, TPos end)
{
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while (n-- > 0)
	{
		assert_geq(n, 0);
		TPos a = i + n;
		TPos b = j + n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
	}
}
template <typename TVal, typename TPos>
static inline void vecswap(TVal *s, size_t slen, TPos i, TPos j, TPos n, TPos begin, TPos end)
{
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while (n-- > 0)
	{
		assert_geq(n, 0);
		TPos a = i + n;
		TPos b = j + n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
	}
}
template <typename TStr, typename TPos>
static inline void vecswap2(
	TStr &s,
	size_t slen,
	TStr &s2,
	TPos i,
	TPos j,
	TPos n,
	TPos begin,
	TPos end)
{
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while (n-- > 0)
	{
		assert_geq(n, 0);
		TPos a = i + n;
		TPos b = j + n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
		swap(s2, slen, a, b);
	}
}
template <typename TVal, typename TPos>
static inline void vecswap2(TVal *s, size_t slen, TVal *s2, TPos i, TPos j, TPos n, TPos begin, TPos end)
{
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while (n-- > 0)
	{
		assert_geq(n, 0);
		TPos a = i + n;
		TPos b = j + n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
		swap(s2, slen, a, b);
	}
}
#define CHAR_AT(ss, aa) ((length(s[ss]) > aa) ? (int)(s[ss][aa]) : hi)
#define CHAR_AT_SUF(si, off) \
	(((off + s[si]) < hlen) ? ((int)(host[off + s[si]])) : (hi))
#define CHAR_AT_SUF_U8(si, off) char_at_suf_u8(host, hlen, s, si, off, hi)
#define CHOOSE_AND_SWAP_RANDOM_PIVOT(sw, ch) \
	{                                        \
                                             \
		a = (rand() % n) + begin;            \
		assert_lt(a, end);                   \
		assert_geq(a, begin);                \
		sw(s, s2, begin, a);                 \
	}
#define CHOOSE_AND_SWAP_SMART_PIVOT(sw, ch)                                  \
	{                                                                        \
		a = begin;                                                           \
                                                                             \
		if (n >= 5)                                                          \
		{                                                                    \
			if (ch(begin + 1, depth) == 1 || ch(begin + 1, depth) == 2)      \
				a = begin + 1;                                               \
			else if (ch(begin + 2, depth) == 1 || ch(begin + 2, depth) == 2) \
				a = begin + 2;                                               \
			else if (ch(begin + 3, depth) == 1 || ch(begin + 3, depth) == 2) \
				a = begin + 3;                                               \
			else if (ch(begin + 4, depth) == 1 || ch(begin + 4, depth) == 2) \
				a = begin + 4;                                               \
			if (a != begin)                                                  \
				sw(s, s2, begin, a);                                         \
		}                                                                    \
	}
#define CHOOSE_AND_SWAP_PIVOT CHOOSE_AND_SWAP_SMART_PIVOT
#ifndef NDEBUG
template <typename THost>
bool assertPartitionedSuf(
	const THost &host,
	TIndexOffU *s,
	size_t slen,
	int hi,
	int pivot,
	size_t begin,
	size_t end,
	size_t depth)
{
	size_t hlen = host.length();
	int state = 0;
	for (size_t i = begin; i < end; i++)
	{
		switch (state)
		{
		case 0:
			if (CHAR_AT_SUF(i, depth) < pivot)
			{
				state = 1;
				break;
			}
			else if (CHAR_AT_SUF(i, depth) > pivot)
			{
				state = 2;
				break;
			}
			assert_eq(CHAR_AT_SUF(i, depth), pivot);
			break;
		case 1:
			if (CHAR_AT_SUF(i, depth) > pivot)
			{
				state = 2;
				break;
			}
			else if (CHAR_AT_SUF(i, depth) == pivot)
			{
				state = 3;
				break;
			}
			assert_lt(CHAR_AT_SUF(i, depth), pivot);
			break;
		case 2:
			if (CHAR_AT_SUF(i, depth) == pivot)
			{
				state = 3;
				break;
			}
			assert_gt(CHAR_AT_SUF(i, depth), pivot);
			break;
		case 3:
			assert_eq(CHAR_AT_SUF(i, depth), pivot);
			break;
		}
	}
	return true;
}
template <typename THost>
bool assertPartitionedSuf2(
	const THost &host,
	TIndexOffU *s,
	size_t slen,
	int hi,
	int pivot,
	size_t begin,
	size_t end,
	size_t depth)
{
	size_t hlen = host.length();
	int state = 0;
	for (size_t i = begin; i < end; i++)
	{
		switch (state)
		{
		case 0:
			if (CHAR_AT_SUF(i, depth) == pivot)
			{
				state = 1;
				break;
			}
			else if (CHAR_AT_SUF(i, depth) > pivot)
			{
				state = 2;
				break;
			}
			assert_lt(CHAR_AT_SUF(i, depth), pivot);
			break;
		case 1:
			if (CHAR_AT_SUF(i, depth) > pivot)
			{
				state = 2;
				break;
			}
			assert_eq(CHAR_AT_SUF(i, depth), pivot);
			break;
		case 2:
			assert_gt(CHAR_AT_SUF(i, depth), pivot);
			break;
		}
	}
	return true;
}
#endif
static inline void sanityCheckInputSufs(TIndexOffU *s, size_t slen)
{
	assert_gt(slen, 0);
	for (size_t i = 0; i < slen; i++)
	{
		for (size_t j = i + 1; j < slen; j++)
		{
			assert_neq(s[i], s[j]);
		}
	}
}
template <typename T>
void sanityCheckOrderedSufs(
	const T &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	size_t upto,
	size_t lower = 0,
	size_t upper = OFF_MASK)
{
	assert_lt(s[0], hlen);
	upper = min<size_t>(upper, slen - 1);
	for (size_t i = lower; i < upper; i++)
	{
		if (s[i + 1] >= hlen)
			continue;
#ifndef NDEBUG
		if (upto == OFF_MASK)
		{
			assert(sstr_suf_lt(host, s[i], hlen, host, s[i + 1], hlen, false));
		}
		else
		{
			if (sstr_suf_upto_lt(host, s[i], host, s[i + 1], upto, false))
			{
			}
		}
#endif
	}
}
template <typename T>
void mkeyQSortSuf(
	const T &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	int hi,
	size_t begin,
	size_t end,
	size_t depth,
	size_t upto = OFF_MASK)
{
#define MQS_RECURSE_SUF(nbegin, nend, ndepth)                                  \
	{                                                                          \
		assert(nbegin > begin || nend < end || ndepth > depth);                \
		if (ndepth < upto)                                                     \
		{                                                                      \
			mkeyQSortSuf(host, hlen, s, slen, hi, nbegin, nend, ndepth, upto); \
		}                                                                      \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t a, b, c, d, r;
	size_t n = end - begin;
	if (n <= 1)
		return;
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF);
	int v = CHAR_AT_SUF(begin, depth);
#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for (size_t i = begin; i < end; i++)
		{
			if (depth < (hlen - s[i]))
			{
				stillInBounds = true;
				break;
			}
			else
			{
			}
		}
		assert(stillInBounds);
	}
#endif
	a = b = begin;
	c = d = end - 1;
	while (true)
	{
		int bc = 0;
		while (b <= c && v >= (bc = CHAR_AT_SUF(b, depth)))
		{
			if (v == bc)
			{
				SWAP(s, a, b);
				a++;
			}
			b++;
		}
		int cc = 0;
		while (b <= c && v <= (cc = CHAR_AT_SUF(c, depth)))
		{
			if (v == cc)
			{
				SWAP(s, c, d);
				d--;
			}
			c--;
		}
		if (b > c)
			break;
		SWAP(s, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end - 1);
	assert_lt(d - c, n);
	assert_lt(b - a, n);
	assert(assertPartitionedSuf(host, s, slen, hi, v, begin, end, depth));
	r = min(a - begin, b - a);
	VECSWAP(s, begin, b - r, r);
	r = min(d - c, end - d - 1);
	VECSWAP(s, b, end - r, r);
	assert(assertPartitionedSuf2(host, s, slen, hi, v, begin, end, depth));
	r = b - a;
	if (r > 0)
	{
		MQS_RECURSE_SUF(begin, begin + r, depth);
	}
	if (v != hi)
	{
		MQS_RECURSE_SUF(begin + r, begin + r + (a - begin) + (end - d - 1), depth + 1);
	}
	r = d - c;
	if (r > 0 && v < hi - 1)
	{
		MQS_RECURSE_SUF(end - r, end, depth);
	}
}
template <typename T>
void mkeyQSortSuf(
	const T &host,
	TIndexOffU *s,
	size_t slen,
	int hi,
	bool verbose = false,
	bool sanityCheck = false,
	size_t upto = OFF_MASK)
{
	size_t hlen = host.length();
	assert_gt(slen, 0);
	if (sanityCheck)
		sanityCheckInputSufs(s, slen);
	mkeyQSortSuf(host, hlen, s, slen, hi, (size_t)0, slen, (size_t)0, upto);
	if (sanityCheck)
		sanityCheckOrderedSufs(host, hlen, s, slen, upto);
}
struct QSortRange
{
	size_t begin;
	size_t end;
	size_t depth;
};
template <typename T>
void mkeyQSortSuf2(
	const T &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	TIndexOffU *s2,
	int hi,
	size_t _begin,
	size_t _end,
	size_t _depth,
	size_t upto = OFF_MASK,
	EList<size_t> *boundaries = NULL)
{
	ELList<QSortRange, 3, 1024> block_list;
	while (true)
	{
		size_t begin = 0, end = 0, depth = 0;
		if (block_list.size() == 0)
		{
			begin = _begin;
			end = _end;
			depth = _depth;
		}
		else
		{
			if (block_list.back().size() > 0)
			{
				begin = block_list.back()[0].begin;
				end = block_list.back()[0].end;
				depth = block_list.back()[0].depth;
				block_list.back().erase(0);
			}
			else
			{
				block_list.resize(block_list.size() - 1);
				if (block_list.size() == 0)
				{
					break;
				}
			}
		}
		if (depth == upto)
		{
			if (boundaries != NULL)
			{
				(*boundaries).push_back(end);
			}
			continue;
		}
		assert_leq(begin, slen);
		assert_leq(end, slen);
		size_t a, b, c, d, r;
		size_t n = end - begin;
		if (n <= 1)
		{
			if (n == 1 && boundaries != NULL)
			{
				boundaries->push_back(end);
			}
			continue;
		}
		CHOOSE_AND_SWAP_PIVOT(SWAP2, CHAR_AT_SUF);
		int v = CHAR_AT_SUF(begin, depth);
#ifndef NDEBUG
		{
			bool stillInBounds = false;
			for (size_t i = begin; i < end; i++)
			{
				if (depth < (hlen - s[i]))
				{
					stillInBounds = true;
					break;
				}
				else
				{
				}
			}
			assert(stillInBounds);
		}
#endif
		a = b = begin;
		c = d = end - 1;
		while (true)
		{
			int bc = 0;
			while (b <= c && v >= (bc = CHAR_AT_SUF(b, depth)))
			{
				if (v == bc)
				{
					SWAP2(s, s2, a, b);
					a++;
				}
				b++;
			}
			int cc = 0;
			while (b <= c && v <= (cc = CHAR_AT_SUF(c, depth)))
			{
				if (v == cc)
				{
					SWAP2(s, s2, c, d);
					d--;
				}
				c--;
			}
			if (b > c)
				break;
			SWAP2(s, s2, b, c);
			b++;
			c--;
		}
		assert(a > begin || c < end - 1);
		assert_lt(d - c, n);
		assert_lt(b - a, n);
		assert(assertPartitionedSuf(host, s, slen, hi, v, begin, end, depth));
		r = min(a - begin, b - a);
		VECSWAP2(s, s2, begin, b - r, r);
		r = min(d - c, end - d - 1);
		VECSWAP2(s, s2, b, end - r, r);
		assert(assertPartitionedSuf2(host, s, slen, hi, v, begin, end, depth));
		r = b - a;
		block_list.expand();
		block_list.back().clear();
		if (r > 0)
		{
			block_list.back().expand();
			block_list.back().back().begin = begin;
			block_list.back().back().end = begin + r;
			block_list.back().back().depth = depth;
		}
		block_list.back().expand();
		block_list.back().back().begin = begin + r;
		block_list.back().back().end = begin + r + (a - begin) + (end - d - 1);
		block_list.back().back().depth = depth + 1;
		r = d - c;
		if (r > 0)
		{
			block_list.back().expand();
			block_list.back().back().begin = end - r;
			block_list.back().back().end = end;
			block_list.back().back().depth = depth;
		}
	}
}
template <typename T>
void mkeyQSortSuf2(
	const T &host,
	TIndexOffU *s,
	size_t slen,
	TIndexOffU *s2,
	int hi,
	bool verbose = false,
	bool sanityCheck = false,
	size_t upto = OFF_MASK,
	EList<size_t> *boundaries = NULL)
{
	size_t hlen = host.length();
	if (sanityCheck)
		sanityCheckInputSufs(s, slen);
	TIndexOffU *sOrig = NULL;
	if (sanityCheck)
	{
		sOrig = new TIndexOffU[slen];
		memcpy(sOrig, s, OFF_SIZE * slen);
	}
	mkeyQSortSuf2(host, hlen, s, slen, s2, hi, (size_t)0, slen, (size_t)0, upto, boundaries);
	if (sanityCheck)
	{
		sanityCheckOrderedSufs(host, hlen, s, slen, upto);
		for (size_t i = 0; i < slen; i++)
		{
			assert_eq(s[i], sOrig[s2[i]]);
		}
		delete[] sOrig;
	}
}
template <typename T>
class DifferenceCoverSample;
template <typename T1, typename T2>
inline bool sufDcLt(
	const T1 &host,
	const T2 &s1,
	const T2 &s2,
	const DifferenceCoverSample<T1> &dc,
	bool sanityCheck = false)
{
	size_t diff = dc.tieBreakOff(s1, s2);
	ASSERT_ONLY(size_t hlen = host.length());
	assert_lt(diff, dc.v());
	assert_lt(diff, hlen - s1);
	assert_lt(diff, hlen - s2);
	if (sanityCheck)
	{
		for (size_t i = 0; i < diff; i++)
		{
			assert_eq(host[s1 + i], host[s2 + i]);
		}
	}
	bool ret = dc.breakTie(s1 + diff, s2 + diff) < 0;
#ifndef NDEBUG
	if (sanityCheck && ret != sstr_suf_lt(host, s1, hlen, host, s2, hlen, false))
	{
		assert(false);
	}
#endif
	return ret;
}
template <typename T>
inline void qsortSufDc(
	const T &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	const DifferenceCoverSample<T> &dc,
	size_t begin,
	size_t end,
	bool sanityCheck = false)
{
	assert_leq(end, slen);
	assert_lt(begin, slen);
	assert_gt(end, begin);
	size_t n = end - begin;
	if (n <= 1)
		return;
	size_t a = (rand() % n) + begin;
	assert_lt(a, end);
	assert_geq(a, begin);
	SWAP(s, end - 1, a);
	size_t cur = 0;
	for (size_t i = begin; i < end - 1; i++)
	{
		if (sufDcLt(host, s[i], s[end - 1], dc, sanityCheck))
		{
			if (sanityCheck)
				assert(dollarLt(suffix(host, s[i]), suffix(host, s[end - 1])));
			assert_lt(begin + cur, end - 1);
			SWAP(s, i, begin + cur);
			cur++;
		}
	}
	assert_lt(cur, end - begin);
	SWAP(s, end - 1, begin + cur);
	if (begin + cur > begin)
		qsortSufDc(host, hlen, s, slen, dc, begin, begin + cur);
	if (end > begin + cur + 1)
		qsortSufDc(host, hlen, s, slen, dc, begin + cur + 1, end);
}
template <typename T1, typename T2>
void mkeyQSortSufDcU8(
	const T1 &host1,
	const T2 &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	const DifferenceCoverSample<T1> &dc,
	int hi,
	bool verbose = false,
	bool sanityCheck = false)
{
	if (sanityCheck)
		sanityCheckInputSufs(s, slen);
	mkeyQSortSufDcU8(host1, host, hlen, s, slen, dc, hi, 0, slen, 0, sanityCheck);
	if (sanityCheck)
		sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK);
}
template <typename T1, typename T2>
inline bool sufDcLtU8(
	const T1 &host1,
	const T2 &host,
	size_t hlen,
	size_t s1,
	size_t s2,
	const DifferenceCoverSample<T1> &dc,
	bool sanityCheck = false)
{
	hlen += 0;
	size_t diff = dc.tieBreakOff((TIndexOffU)s1, (TIndexOffU)s2);
	assert_lt(diff, dc.v());
	assert_lt(diff, hlen - s1);
	assert_lt(diff, hlen - s2);
	if (sanityCheck)
	{
		for (size_t i = 0; i < diff; i++)
		{
			assert_eq(host[s1 + i], host1[s2 + i]);
		}
	}
	bool ret = dc.breakTie((TIndexOffU)(s1 + diff), (TIndexOffU)(s2 + diff)) < 0;
#ifndef NDEBUG
	bool ret2 = sstr_suf_lt(host1, s1, hlen, host, s2, hlen, false);
	assert(!sanityCheck || ret == ret2);
#endif
	return ret;
}
template <typename T1, typename T2>
inline void qsortSufDcU8(
	const T1 &host1,
	const T2 &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	const DifferenceCoverSample<T1> &dc,
	size_t begin,
	size_t end,
	bool sanityCheck = false)
{
	assert_leq(end, slen);
	assert_lt(begin, slen);
	assert_gt(end, begin);
	size_t n = end - begin;
	if (n <= 1)
		return;
	size_t a = (rand() % n) + begin;
	assert_lt(a, end);
	assert_geq(a, begin);
	SWAP(s, end - 1, a);
	size_t cur = 0;
	for (size_t i = begin; i < end - 1; i++)
	{
		if (sufDcLtU8(host1, host, hlen, s[i], s[end - 1], dc, sanityCheck))
		{
#ifndef NDEBUG
			if (sanityCheck)
			{
				assert(sstr_suf_lt(host1, s[i], hlen, host1, s[end - 1], hlen, false));
			}
			assert_lt(begin + cur, end - 1);
#endif
			SWAP(s, i, begin + cur);
			cur++;
		}
	}
	assert_lt(cur, end - begin);
	SWAP(s, end - 1, begin + cur);
	if (begin + cur > begin)
		qsortSufDcU8(host1, host, hlen, s, slen, dc, begin, begin + cur);
	if (end > begin + cur + 1)
		qsortSufDcU8(host1, host, hlen, s, slen, dc, begin + cur + 1, end);
}
#define BUCKET_SORT_CUTOFF (4 * 1024 * 1024)
#define SELECTION_SORT_CUTOFF 6
extern TIndexOffU bkts[4][4 * 1024 * 1024];
template <typename TStr>
inline uint8_t get_uint8(const TStr &t, size_t off)
{
	return t[off];
}
template <>
inline uint8_t get_uint8<S2bDnaString>(const S2bDnaString &t, size_t off)
{
	return (uint8_t)t[off];
}
template <typename TStr>
static inline int char_at_suf_u8(
	const TStr &host,
	size_t hlen,
	TIndexOffU *s,
	size_t si,
	size_t off,
	uint8_t hi)
{
	return ((off + s[si]) < hlen) ? get_uint8(host, off + s[si]) : (hi);
}
template <typename T1, typename T2>
static void selectionSortSufDcU8(
	const T1 &host1,
	const T2 &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	const DifferenceCoverSample<T1> &dc,
	uint8_t hi,
	size_t begin,
	size_t end,
	size_t depth,
	bool sanityCheck = false)
{
#define ASSERT_SUF_LT(l, r)                                        \
	if (sanityCheck &&                                             \
		!sstr_suf_lt(host1, s[l], hlen, host1, s[r], hlen, false)) \
	{                                                              \
		assert(false);                                             \
	}
	assert_gt(end, begin + 1);
	assert_leq(end - begin, SELECTION_SORT_CUTOFF);
	assert_eq(hi, 4);
	size_t v = dc.v();
	if (end == begin + 2)
	{
		size_t off = dc.tieBreakOff(s[begin], s[begin + 1]);
		if (off + s[begin] >= hlen ||
			off + s[begin + 1] >= hlen)
		{
			off = OFF_MASK;
		}
		if (off != OFF_MASK)
		{
			if (off < depth)
			{
				qsortSufDcU8<T1, T2>(host1, host, hlen, s, slen, dc,
									 begin, end, sanityCheck);
				if (sanityCheck)
				{
					sanityCheckOrderedSufs(host1, hlen, s, slen,
										   OFF_MASK, begin, end);
				}
				return;
			}
			v = off - depth + 1;
		}
	}
	assert_leq(v, dc.v());
	size_t lim = v;
	assert_geq(lim, 0);
	for (size_t i = begin; i < end - 1; i++)
	{
		size_t targ = i;
		size_t targoff = depth + s[i];
		for (size_t j = i + 1; j < end; j++)
		{
			assert_neq(j, targ);
			size_t joff = depth + s[j];
			size_t k;
			for (k = 0; k <= lim; k++)
			{
				assert_neq(j, targ);
				uint8_t jc = (k + joff < hlen) ? get_uint8(host, k + joff) : hi;
				uint8_t tc = (k + targoff < hlen) ? get_uint8(host, k + targoff) : hi;
				assert(jc != hi || tc != hi);
				if (jc > tc)
				{
					ASSERT_SUF_LT(targ, j);
					break;
				}
				else if (jc < tc)
				{
					ASSERT_SUF_LT(j, targ);
					targ = j;
					targoff = joff;
					break;
				}
				else if (k == lim)
				{
					assert_leq(k + joff + 1, hlen);
					assert_leq(k + targoff + 1, hlen);
					if (k + joff + 1 == hlen)
					{
						assert_neq(k + targoff + 1, hlen);
						ASSERT_SUF_LT(targ, j);
						break;
					}
					else if (k + targoff + 1 == hlen)
					{
						ASSERT_SUF_LT(j, targ);
						targ = j;
						targoff = joff;
						break;
					}
				}
				else
				{
				}
			}
			if (k == lim + 1)
			{
				assert_neq(j, targ);
				if (sufDcLtU8(host1, host, hlen, s[j], s[targ], dc, sanityCheck))
				{
					assert(!sufDcLtU8(host1, host, hlen, s[targ], s[j], dc, sanityCheck));
					ASSERT_SUF_LT(j, targ);
					targ = j;
					targoff = joff;
				}
				else
				{
					assert(sufDcLtU8(host1, host, hlen, s[targ], s[j], dc, sanityCheck));
					ASSERT_SUF_LT(targ, j);
				}
			}
		}
		if (i != targ)
		{
			ASSERT_SUF_LT(targ, i);
			TIndexOffU tmp = s[i];
			s[i] = s[targ];
			s[targ] = tmp;
		}
		for (size_t j = i + 1; j < end; j++)
		{
			ASSERT_SUF_LT(i, j);
		}
	}
	if (sanityCheck)
	{
		sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
	}
}
template <typename T1, typename T2>
static void bucketSortSufDcU8(
	const T1 &host1,
	const T2 &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	const DifferenceCoverSample<T1> &dc,
	uint8_t hi,
	size_t _begin,
	size_t _end,
	size_t _depth,
	bool sanityCheck = false)
{
	TIndexOffU *bkts[4];
	for (size_t i = 0; i < 4; i++)
	{
		bkts[i] = new TIndexOffU[4 * 1024 * 1024];
	}
	ELList<size_t, 5, 1024> block_list;
	bool first = true;
	while (true)
	{
		size_t begin = 0, end = 0;
		if (first)
		{
			begin = _begin;
			end = _end;
			first = false;
		}
		else
		{
			if (block_list.size() == 0)
			{
				break;
			}
			if (block_list.back().size() > 1)
			{
				end = block_list.back().back();
				block_list.back().pop_back();
				begin = block_list.back().back();
			}
			else
			{
				block_list.resize(block_list.size() - 1);
				if (block_list.size() == 0)
				{
					break;
				}
			}
		}
		size_t depth = block_list.size() + _depth;
		assert_leq(end - begin, BUCKET_SORT_CUTOFF);
		assert_eq(hi, 4);
		if (end <= begin + 1)
		{
			continue;
		}
		if (depth > dc.v())
		{
			qsortSufDcU8<T1, T2>(host1, host, hlen, s, slen, dc, begin, end, sanityCheck);
			continue;
		}
		if (end - begin <= SELECTION_SORT_CUTOFF)
		{
			selectionSortSufDcU8(host1, host, hlen, s, slen, dc, hi,
								 begin, end, depth, sanityCheck);
			if (sanityCheck)
			{
				sanityCheckOrderedSufs(host1, hlen, s, slen,
									   OFF_MASK, begin, end);
			}
			continue;
		}
		size_t cnts[] = {0, 0, 0, 0, 0};
		for (size_t i = begin; i < end; i++)
		{
			size_t off = depth + s[i];
			uint8_t c = (off < hlen) ? get_uint8(host, off) : hi;
			assert_leq(c, 4);
			if (c == 0)
			{
				s[begin + cnts[0]++] = s[i];
			}
			else
			{
				bkts[c - 1][cnts[c]++] = s[i];
			}
		}
		assert_eq(cnts[0] + cnts[1] + cnts[2] + cnts[3] + cnts[4], end - begin);
		size_t cur = begin + cnts[0];
		if (cnts[1] > 0)
		{
			memcpy(&s[cur], bkts[0], cnts[1] << (OFF_SIZE / 4 + 1));
			cur += cnts[1];
		}
		if (cnts[2] > 0)
		{
			memcpy(&s[cur], bkts[1], cnts[2] << (OFF_SIZE / 4 + 1));
			cur += cnts[2];
		}
		if (cnts[3] > 0)
		{
			memcpy(&s[cur], bkts[2], cnts[3] << (OFF_SIZE / 4 + 1));
			cur += cnts[3];
		}
		if (cnts[4] > 0)
		{
			memcpy(&s[cur], bkts[3], cnts[4] << (OFF_SIZE / 4 + 1));
		}
		block_list.expand();
		block_list.back().clear();
		block_list.back().push_back(begin);
		for (size_t i = 0; i < 4; i++)
		{
			if (cnts[i] > 0)
			{
				block_list.back().push_back(block_list.back().back() + cnts[i]);
			}
		}
	}
	for (size_t i = 0; i < 4; i++)
	{
		delete[] bkts[i];
	}
}
template <typename T1, typename T2>
void mkeyQSortSufDcU8(
	const T1 &host1,
	const T2 &host,
	size_t hlen,
	TIndexOffU *s,
	size_t slen,
	const DifferenceCoverSample<T1> &dc,
	int hi,
	size_t begin,
	size_t end,
	size_t depth,
	bool sanityCheck = false)
{
#define MQS_RECURSE_SUF_DC_U8(nbegin, nend, ndepth)                                              \
	{                                                                                            \
		assert(nbegin > begin || nend < end || ndepth > depth);                                  \
		mkeyQSortSufDcU8(host1, host, hlen, s, slen, dc, hi, nbegin, nend, ndepth, sanityCheck); \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t n = end - begin;
	if (n <= 1)
		return;
	if (depth > dc.v())
	{
		qsortSufDcU8<T1, T2>(host1, host, hlen, s, slen, dc, begin, end, sanityCheck);
		if (sanityCheck)
		{
			sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
		}
		return;
	}
	if (n <= BUCKET_SORT_CUTOFF)
	{
		bucketSortSufDcU8(host1, host, hlen, s, slen, dc,
						  (uint8_t)hi, begin, end, depth, sanityCheck);
		if (sanityCheck)
		{
			sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
		}
		return;
	}
	size_t a, b, c, d, r;
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF_U8);
	int v = CHAR_AT_SUF_U8(begin, depth);
#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for (size_t i = begin; i < end; i++)
		{
			if (depth < (hlen - s[i]))
			{
				stillInBounds = true;
				break;
			}
			else
			{
			}
		}
		assert(stillInBounds);
	}
#endif
	a = b = begin;
	c = d = end - 1;
	while (true)
	{
		int bc = 0;
		while (b <= c && v >= (bc = CHAR_AT_SUF_U8(b, depth)))
		{
			if (v == bc)
			{
				SWAP(s, a, b);
				a++;
			}
			b++;
		}
		int cc = 0;
		while (b <= c && v <= (cc = CHAR_AT_SUF_U8(c, depth)))
		{
			if (v == cc)
			{
				SWAP(s, c, d);
				d--;
			}
			c--;
		}
		if (b > c)
			break;
		SWAP(s, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end - 1);
	assert_lt(d - c, n);
	assert_lt(b - a, n);
	r = min(a - begin, b - a);
	VECSWAP(s, begin, b - r, r);
	r = min(d - c, end - d - 1);
	VECSWAP(s, b, end - r, r);
	r = b - a;
	if (r > 0)
	{
		MQS_RECURSE_SUF_DC_U8(begin, begin + r, depth);
	}
	if (v != hi)
	{
		MQS_RECURSE_SUF_DC_U8(begin + r, begin + r + (a - begin) + (end - d - 1), depth + 1);
	}
	r = d - c;
	if (r > 0 && v < hi - 1)
	{
		MQS_RECURSE_SUF_DC_U8(end - r, end, depth);
	}
}
#endif
