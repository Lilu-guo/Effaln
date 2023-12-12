#ifndef BINARY_SA_SEARCH_H_
#define BINARY_SA_SEARCH_H_
#include <stdint.h>
#include <iostream>
#include <limits>
#include "alphabet.h"
#include "assert_helpers.h"
#include "ds.h"
#include "btypes.h"
template <typename TStr, typename TSufElt>
inline TIndexOffU binarySASearch(
	const TStr &host,
	TIndexOffU qry,
	const EList<TSufElt> &sa)
{
	TIndexOffU lLcp = 0, rLcp = 0;
	TIndexOffU l = 0, r = (TIndexOffU)sa.size() + 1;
	TIndexOffU hostLen = (TIndexOffU)host.length();
	while (true)
	{
		assert_gt(r, l);
		TIndexOffU m = (l + r) >> 1;
		if (m == l)
		{
			if (m > 0 && sa[m - 1] == qry)
			{
				return std::numeric_limits<TIndexOffU>::max();
			}
			assert_leq(m, sa.size());
			return m;
		}
		assert_gt(m, 0);
		TIndexOffU suf = sa[m - 1];
		if (suf == qry)
		{
			return std::numeric_limits<TIndexOffU>::max();
		}
		TIndexOffU lcp = min(lLcp, rLcp);
#ifndef NDEBUG
		if (sstr_suf_upto_neq(host, qry, host, suf, lcp))
		{
			assert(0);
		}
#endif
		while (suf + lcp < hostLen && qry + lcp < hostLen && host[suf + lcp] == host[qry + lcp])
		{
			lcp++;
		}
		bool fell = (suf + lcp == hostLen || qry + lcp == hostLen);
		if ((fell && qry + lcp == hostLen) || (!fell && host[suf + lcp] < host[qry + lcp]))
		{
			l = m;
			lLcp = max(lLcp, lcp);
		}
		else if ((fell && suf + lcp == hostLen) || (!fell && host[suf + lcp] > host[qry + lcp]))
		{
			r = m;
			rLcp = max(rLcp, lcp);
		}
		else
		{
			assert(false);
		}
	}
	assert(false);
	return std::numeric_limits<TIndexOffU>::max();
}
#endif
