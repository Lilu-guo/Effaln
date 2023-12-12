#include <iostream>
#include "scoring.h"
using namespace std;
bool Scoring::scoreFilter(
	int64_t minsc,
	size_t rdlen) const
{
	int64_t sc = (int64_t)(rdlen * match(30));
	return sc >= minsc;
}
int Scoring::maxReadGaps(
	int64_t minsc,
	size_t rdlen) const
{
	int64_t sc = (int64_t)(rdlen * match(30));
	assert_geq(sc, minsc);
	bool first = true;
	int num = 0;
	while (sc >= minsc)
	{
		if (first)
		{
			first = false;
			sc -= readGapOpen();
		}
		else
		{
			sc -= readGapExtend();
		}
		num++;
	}
	assert_gt(num, 0);
	return num - 1;
}
int Scoring::maxRefGaps(
	int64_t minsc,
	size_t rdlen) const
{
	int64_t sc = (int64_t)(rdlen * match(30));
	assert_geq(sc, minsc);
	bool first = true;
	int num = 0;
	while (sc >= minsc)
	{
		sc -= match(30);
		if (first)
		{
			first = false;
			sc -= refGapOpen();
		}
		else
		{
			sc -= refGapExtend();
		}
		num++;
	}
	assert_gt(num, 0);
	return num - 1;
}
bool Scoring::nFilter(const BTDnaString &rd, size_t &ns) const
{
	size_t rdlen = rd.length();
	size_t maxns = nCeil.f<size_t>((double)rdlen);
	assert_geq(rd.length(), 0);
	for (size_t i = 0; i < rdlen; i++)
	{
		if (rd[i] == 4)
		{
			ns++;
			if (ns > maxns)
			{
				return false;
			}
		}
	}
	return true;
}
void Scoring::nFilterPair(
	const BTDnaString *rd1,
	const BTDnaString *rd2, size_t &ns1, size_t &ns2, bool &filt1, bool &filt2) const
{
	filt1 = filt2 = false;
	if (rd1 != NULL && rd2 != NULL && ncatpair)
	{
		size_t rdlen1 = rd1->length();
		size_t rdlen2 = rd2->length();
		size_t maxns = nCeil.f<size_t>((double)(rdlen1 + rdlen2));
		for (size_t i = 0; i < rdlen1; i++)
		{
			if ((*rd1)[i] == 4)
				ns1++;
			if (ns1 > maxns)
			{
				return;
			}
		}
		for (size_t i = 0; i < rdlen2; i++)
		{
			if ((*rd2)[i] == 4)
				ns2++;
			if (ns2 > maxns)
			{
				return;
			}
		}
		filt1 = filt2 = true;
	}
	else
	{
		if (rd1 != NULL)
			filt1 = nFilter(*rd1, ns1);
		if (rd2 != NULL)
			filt2 = nFilter(*rd2, ns2);
	}
}
#ifdef SCORING_MAIN
int main()
{
	{
		cout << "Case 1: Simple 1 ... ";
		Scoring sc = Scoring::base1();
		assert_eq(COST_MODEL_CONSTANT, sc.matchType);
		assert_eq(0, sc.maxRefGaps(0, 10));
		assert_eq(0, sc.maxRefGaps(0, 11));
		assert_eq(0, sc.maxRefGaps(0, 12));
		assert_eq(0, sc.maxRefGaps(0, 13));
		assert_eq(0, sc.maxRefGaps(0, 14));
		assert_eq(0, sc.maxRefGaps(0, 15));
		assert_eq(1, sc.maxRefGaps(0, 16));
		assert_eq(1, sc.maxRefGaps(0, 17));
		assert_eq(1, sc.maxRefGaps(0, 18));
		assert_eq(1, sc.maxRefGaps(0, 19));
		assert_eq(1, sc.maxRefGaps(0, 20));
		assert_eq(2, sc.maxRefGaps(0, 21));
		assert_eq(0, sc.maxReadGaps(0, 10));
		assert_eq(0, sc.maxReadGaps(0, 11));
		assert_eq(0, sc.maxReadGaps(0, 12));
		assert_eq(0, sc.maxReadGaps(0, 13));
		assert_eq(0, sc.maxReadGaps(0, 14));
		assert_eq(1, sc.maxReadGaps(0, 15));
		assert_eq(1, sc.maxReadGaps(0, 16));
		assert_eq(1, sc.maxReadGaps(0, 17));
		assert_eq(1, sc.maxReadGaps(0, 18));
		assert_eq(2, sc.maxReadGaps(0, 19));
		assert_eq(2, sc.maxReadGaps(0, 20));
		assert_eq(2, sc.maxReadGaps(0, 21));
		assert_eq(1, sc.nCeil(1));
		assert_eq(2, sc.nCeil(3));
		assert_eq(2, sc.nCeil(5));
		assert_eq(2, sc.nCeil(7));
		assert_eq(2, sc.nCeil(9));
		assert_eq(3, sc.nCeil(10));
		for (int i = 0; i < 30; i++)
		{
			assert_eq(3, sc.n(i));
			assert_eq(3, sc.mm(i));
		}
		assert_eq(5, sc.gapbar);
		cout << "PASSED" << endl;
	}
	{
		cout << "Case 2: Simple 2 ... ";
		Scoring sc(
			4,
			COST_MODEL_QUAL, 0, -3.0f, -3.0f, DEFAULT_FLOOR_CONST, DEFAULT_FLOOR_LINEAR, 3.0f, 0.4f, COST_MODEL_QUAL, 0, true, 25, 25, 10, 10, 5, -1, false);
		assert_eq(COST_MODEL_CONSTANT, sc.matchType);
		assert_eq(4, sc.matchConst);
		assert_eq(COST_MODEL_QUAL, sc.mmcostType);
		assert_eq(COST_MODEL_QUAL, sc.npenType);
		assert_eq(0, sc.maxRefGaps(0, 8));
		assert_eq(0, sc.maxRefGaps(0, 9));
		assert_eq(1, sc.maxRefGaps(0, 10));
		assert_eq(1, sc.maxRefGaps(0, 11));
		assert_eq(1, sc.maxRefGaps(0, 12));
		assert_eq(1, sc.maxRefGaps(0, 13));
		assert_eq(2, sc.maxRefGaps(0, 14));
		assert_eq(0, sc.maxReadGaps(0, 8));
		assert_eq(1, sc.maxReadGaps(0, 9));
		assert_eq(1, sc.maxReadGaps(0, 10));
		assert_eq(1, sc.maxReadGaps(0, 11));
		assert_eq(2, sc.maxReadGaps(0, 12));
		assert_eq(2, sc.maxReadGaps(0, 13));
		assert_eq(3, sc.maxReadGaps(0, 14));
		assert_eq(1, sc.nCeil(1));
		assert_eq(2, sc.nCeil(2));
		assert_eq(3, sc.nCeil(3));
		assert_eq(4, sc.nCeil(4));
		assert_eq(5, sc.nCeil(5));
		assert_eq(5, sc.nCeil(6));
		assert_eq(5, sc.nCeil(7));
		assert_eq(6, sc.nCeil(8));
		assert_eq(6, sc.nCeil(9));
		for (int i = 0; i < 256; i++)
		{
			assert_eq(i, sc.n(i));
			assert_eq(i, sc.mm(i));
		}
		assert_eq(5, sc.gapbar);
		cout << "PASSED" << endl;
	}
}
#endif
