#include "dp_framer.h"
using namespace std;
bool DynProgFramer::frameSeedExtensionRect(
	int64_t off,
	size_t rdlen, int64_t reflen, size_t maxrdgap, size_t maxrfgap, int64_t maxns, size_t maxhalf, DPRect &rect)
{
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	size_t maxgap = max(maxrdgap, maxrfgap);
	maxgap = min(maxgap, maxhalf);
	int64_t refl = off - 2 * maxgap;
	int64_t refr = off + (rdlen - 1) + 2 * maxgap;
	size_t triml = 0, trimr = 0;
	if (trimToRef_)
	{
		maxns = 0;
	}
	else if (maxns == (int64_t)rdlen)
	{
		maxns--;
	}
	if (refr >= reflen + maxns)
	{
		trimr = (size_t)(refr - (reflen + maxns - 1));
	}
	if (refl < -maxns)
	{
		triml = (size_t)(-refl) - (size_t)maxns;
	}
	rect.refl_pretrim = refl;
	rect.refr_pretrim = refr;
	rect.refl = refl + triml;
	rect.refr = refr - trimr;
	rect.triml = triml;
	rect.trimr = trimr;
	rect.maxgap = maxgap;
	rect.corel = maxgap;
	rect.corer = rect.corel + 2 * maxgap;
	assert(rect.repOk());
	return !rect.entirelyTrimmed();
}
bool DynProgFramer::frameFindMateAnchorLeftRect(
	int64_t ll,
	int64_t lr, int64_t rl, int64_t rr, size_t rdlen, int64_t reflen, size_t maxrdgap, size_t maxrfgap, int64_t maxns, size_t maxhalf, DPRect &rect) const
{
	assert_geq(lr, ll);
	assert_geq(rr, rl);
	assert_geq(rr, lr);
	assert_geq(rl, ll);
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	size_t triml = 0, trimr = 0;
	size_t maxgap = max(maxrdgap, maxrfgap);
	maxgap = max(maxgap, maxhalf);
	int64_t pad_left = maxgap;
	int64_t pad_right = maxgap;
	int64_t en_left = rl;
	int64_t en_right = rr;
	int64_t st_left = en_left - (rdlen - 1);
	ASSERT_ONLY(int64_t st_right = en_right - (rdlen - 1));
	int64_t en_right_pad = en_right + pad_right;
	ASSERT_ONLY(int64_t en_left_pad = en_left - pad_left);
	ASSERT_ONLY(int64_t st_right_pad = st_right + pad_right);
	int64_t st_left_pad = st_left - pad_left;
	assert_leq(st_left, en_left);
	assert_geq(en_right, st_right);
	assert_leq(st_left_pad, en_left_pad);
	assert_geq(en_right_pad, st_right_pad);
	int64_t refl = st_left_pad;
	int64_t refr = en_right_pad;
	if (trimToRef_)
	{
		maxns = 0;
	}
	else if (maxns == (int64_t)rdlen)
	{
		maxns--;
	}
	if (refr >= reflen + maxns)
	{
		trimr = (size_t)(refr - (reflen + maxns - 1));
	}
	if (refl < -maxns)
	{
		triml = (size_t)(-refl) - (size_t)maxns;
	}
	size_t width = (size_t)(refr - refl + 1);
	rect.refl_pretrim = refl;
	rect.refr_pretrim = refr;
	rect.refl = refl + triml;
	rect.refr = refr - trimr;
	rect.triml = triml;
	rect.trimr = trimr;
	rect.maxgap = maxgap;
	rect.corel = maxgap;
	rect.corer = width - maxgap - 1;
	assert(rect.repOk());
	return !rect.entirelyTrimmed();
}
bool DynProgFramer::frameFindMateAnchorRightRect(
	int64_t ll,
	int64_t lr, int64_t rl, int64_t rr, size_t rdlen, int64_t reflen, size_t maxrdgap, size_t maxrfgap, int64_t maxns, size_t maxhalf, DPRect &rect) const
{
	assert_geq(lr, ll);
	assert_geq(rr, rl);
	assert_geq(rr, lr);
	assert_geq(rl, ll);
	assert_gt(rdlen, 0);
	assert_gt(reflen, 0);
	size_t triml = 0, trimr = 0;
	size_t maxgap = max(maxrdgap, maxrfgap);
	maxgap = max(maxgap, maxhalf);
	int64_t pad_left = maxgap;
	int64_t pad_right = maxgap;
	int64_t st_left = ll;
	int64_t st_right = lr;
	ASSERT_ONLY(int64_t en_left = st_left + (rdlen - 1));
	int64_t en_right = st_right + (rdlen - 1);
	int64_t en_right_pad = en_right + pad_right;
	ASSERT_ONLY(int64_t en_left_pad = en_left - pad_left);
	ASSERT_ONLY(int64_t st_right_pad = st_right + pad_right);
	int64_t st_left_pad = st_left - pad_left;
	assert_leq(st_left, en_left);
	assert_geq(en_right, st_right);
	assert_leq(st_left_pad, en_left_pad);
	assert_geq(en_right_pad, st_right_pad);
	int64_t refl = st_left_pad;
	int64_t refr = en_right_pad;
	if (trimToRef_)
	{
		maxns = 0;
	}
	else if (maxns == (int64_t)rdlen)
	{
		maxns--;
	}
	if (refr >= reflen + maxns)
	{
		trimr = (size_t)(refr - (reflen + maxns - 1));
	}
	if (refl < -maxns)
	{
		triml = (size_t)(-refl) - (size_t)maxns;
	}
	size_t width = (size_t)(refr - refl + 1);
	rect.refl_pretrim = refl;
	rect.refr_pretrim = refr;
	rect.refl = refl + triml;
	rect.refr = refr - trimr;
	rect.triml = triml;
	rect.trimr = trimr;
	rect.maxgap = maxgap;
	rect.corel = maxgap;
	rect.corer = width - maxgap - 1;
	assert(rect.repOk());
	return !rect.entirelyTrimmed();
}
#ifdef MAIN_DP_FRAMER
#include <iostream>
static void testCaseFindMateAnchorLeft(
	const char *testName,
	bool trimToRef,
	int64_t ll,
	int64_t lr,
	int64_t rl,
	int64_t rr,
	size_t rdlen,
	size_t reflen,
	size_t maxrdgap,
	size_t maxrfgap,
	size_t ex_width,
	size_t ex_solwidth,
	size_t ex_trimup,
	size_t ex_trimdn,
	int64_t ex_refl,
	int64_t ex_refr,
	const char *ex_st,
	const char *ex_en)
{
	cerr << testName << "...";
	DynProgFramer fr(trimToRef);
	size_t width, solwidth;
	int64_t refl, refr;
	EList<bool> st, en;
	size_t trimup, trimdn;
	size_t maxhalf = 1000; 
	size_t maxgaps = 0;
	fr.frameFindMateAnchorLeft(
		ll,
		lr, rl, rr, rdlen, reflen, maxrdgap, maxrfgap, maxns, maxhalf, width, maxgaps, trimup, trimdn, refl, refr, st, en);
	assert_eq(ex_width, width);
	assert_eq(ex_solwidth, solwidth);
	assert_eq(ex_trimup, trimup);
	assert_eq(ex_trimdn, trimdn);
	assert_eq(ex_refl, refl);
	assert_eq(ex_refr, refr);
	for (size_t i = 0; i < width; i++)
	{
		assert_eq((ex_st[i] == '1'), st[i]);
		assert_eq((ex_en[i] == '1'), en[i]);
	}
	cerr << "PASSED" << endl;
}
static void testCaseFindMateAnchorRight(
	const char *testName,
	bool trimToRef,
	int64_t ll,
	int64_t lr,
	int64_t rl,
	int64_t rr,
	size_t rdlen,
	size_t reflen,
	size_t maxrdgap,
	size_t maxrfgap,
	size_t ex_width,
	size_t ex_solwidth,
	size_t ex_trimup,
	size_t ex_trimdn,
	int64_t ex_refl,
	int64_t ex_refr,
	const char *ex_st,
	const char *ex_en)
{
	cerr << testName << "...";
	DynProgFramer fr(trimToRef);
	size_t width, solwidth;
	size_t maxgaps;
	int64_t refl, refr;
	EList<bool> st, en;
	size_t trimup, trimdn;
	size_t maxhalf = 1000; 
	fr.frameFindMateAnchorRight(
		ll,
		lr, rl, rr, rdlen, reflen, maxrdgap, maxrfgap, maxns, maxhalf, width, maxgaps, trimup, trimdn, refl, refr, st, en);
	assert_eq(ex_width, width);
	assert_eq(ex_trimup, trimup);
	assert_eq(ex_trimdn, trimdn);
	assert_eq(ex_refl, refl);
	assert_eq(ex_refr, refr);
	for (size_t i = 0; i < width; i++)
	{
		assert_eq((ex_st[i] == '1'), st[i]);
		assert_eq((ex_en[i] == '1'), en[i]);
	}
	cerr << "PASSED" << endl;
}
int main(void)
{
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft1",
		false,
		3, 15, 10, 16, 5, 30, 3, 3, 13, 0, 0, 3, 19, "1111111111111", "0001111111000");
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft2",
		false,
		9, 14, 10, 15, 5, 30, 2, 2, 7, 3, 0, 7, 17, "0011111", "1111100");
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft3",
		true,
		9, 14, 10, 15, 5, 17, 2, 2, 7, 3, 0, 7, 17, "0011111", "1111100");
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft4",
		true,
		9, 14, 10, 15, 5, 15, 2, 2, 6, 3, 1, 7, 16, "001111", "111100");
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft5",
		true,
		1, 7, 2, 7, 5, 9, 2, 2, 7, 3, 0, -1, 9, "0011111", "1111100");
	testCaseFindMateAnchorLeft(
		"FindMateAnchorLeft6",
		false,
		8, 8, 10, 15, 5, 30, 4, 2, 6, 4, 2, 6, 15, "001000", "111111");
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight1",
		false,
		10, 16, 11, 23, 5, 30, 3, 3, 13, 0, 0, 7, 23, "0001111111000", "1111111111111");
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight2",
		false,
		6, 11, 13, 18, 5, 30, 2, 2, 7, 3, 0, 7, 17, "1111100", "0011111");
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight3",
		true,
		0, 5, 7, 11, 5, 30, 2, 2, 7, 3, 0, 1, 11, "1111100", "0011111");
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight4",
		true,
		-3, 2, 4, 10, 5, 30, 2, 2, 5, 5, 0, 0, 8, "11100", "11111");
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight5",
		true,
		-3, 2, 4, 10, 5, 7, 2, 2, 3, 5, 2, 0, 6, "111", "111");
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight6",
		false,
		6, 11, 14, 14, 5, 30, 4, 2, 6, 2, 4, 6, 15, "111111", "000010");
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight7",
		false,
		6, 11, 14, 14, 5, 30, 2, 4, 4, 6, 2, 8, 15, "1111", "0010");
	testCaseFindMateAnchorRight(
		"FindMateAnchorRight8",
		true,
		-37, 13, -37, 52, 10, 53, 0, 0, 14, 37, 0, 0, 22, "11111111111111", "11111111111111");
}
#endif
