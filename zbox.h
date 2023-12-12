#ifndef ZBOX_H_
#define ZBOX_H_
#include "btypes.h"
template <typename T>
void calcZ(const T &s,
		   TIndexOffU off,
		   EList<TIndexOffU> &z,
		   bool verbose = false,
		   bool sanityCheck = false)
{
	size_t lCur = 0, rCur = 0;
	size_t zlen = z.size();
	size_t slen = s.length();
	assert_gt(zlen, 0);
	assert_eq(z[0], 0);
	for (size_t k = 1; k < zlen && k + off < slen; k++)
	{
		assert_lt(lCur, k);
		assert(z[lCur] == 0 || z[lCur] == rCur - lCur + 1);
		if (k > rCur)
		{
			size_t ki = k;
			while (off + ki < s.length() && s[off + ki] == s[off + ki - k])
				ki++;
			z[k] = (TIndexOffU)(ki - k);
			assert_lt(off + z[k], slen);
			if (z[k] > 0)
			{
				lCur = k;
				rCur = k + z[k] - 1;
			}
		}
		else
		{
			size_t betaLen = rCur - k + 1;
			size_t kPrime = k - lCur;
			assert_eq(s[off + k], s[off + kPrime]);
			if (z[kPrime] < betaLen)
			{
				z[k] = z[kPrime];
				assert_lt(off + z[k], slen);
			}
			else if (z[kPrime] > 0)
			{
				int q = 0;
				while (off + q + rCur + 1 < s.length() && s[off + q + rCur + 1] == s[off + betaLen + q])
					q++;
				z[k] = (TIndexOffU)(betaLen + q);
				assert_lt(off + z[k], slen);
				rCur = rCur + q;
				assert_geq(k, lCur);
				lCur = k;
			}
			else
			{
				z[k] = 0;
				assert_lt(off + z[k], slen);
			}
		}
	}
#ifndef NDEBUG
	if (sanityCheck)
	{
		assert_eq(0, z[0]);
		for (size_t i = 1; i < z.size(); i++)
		{
			size_t j;
			for (j = i; off + j < s.length(); j++)
			{
				if (s[off + j] != s[off + j - i])
					break;
			}
			assert_eq(j - i, z[i]);
		}
	}
#endif
}
#endif
