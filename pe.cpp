#include "assert_helpers.h"
#include "pe.h"
using namespace std;
int PairedEndPolicy::peClassifyPair(
	int64_t off1,
	size_t len1, bool fw1, int64_t off2, size_t len2, bool fw2) const
{
	assert_gt(len1, 0);
	assert_gt(len2, 0);
	size_t maxfrag = maxfrag_;
	if (len1 > maxfrag && expandToFit_)
		maxfrag = len1;
	if (len2 > maxfrag && expandToFit_)
		maxfrag = len2;
	size_t minfrag = minfrag_;
	if (minfrag < 1)
	{
		minfrag = 1;
	}
	bool oneLeft = false;
	if (pol_ == PE_POLICY_FF)
	{
		if (fw1 != fw2)
		{
			return PE_ALS_DISCORD;
		}
		oneLeft = fw1;
	}
	else if (pol_ == PE_POLICY_RR)
	{
		if (fw1 != fw2)
		{
			return PE_ALS_DISCORD;
		}
		oneLeft = !fw1;
	}
	else if (pol_ == PE_POLICY_FR)
	{
		if (fw1 == fw2)
		{
			return PE_ALS_DISCORD;
		}
		oneLeft = fw1;
	}
	else if (pol_ == PE_POLICY_RF)
	{
		if (fw1 == fw2)
		{
			return PE_ALS_DISCORD;
		}
		oneLeft = !fw1;
	}
	int64_t fraglo = min<int64_t>(off1, off2);
	int64_t fraghi = max<int64_t>(off1 + len1, off2 + len2);
	assert_gt(fraghi, fraglo);
	size_t frag = (size_t)(fraghi - fraglo);
	if (frag > maxfrag || frag < minfrag)
	{
		return PE_ALS_DISCORD;
	}
	int64_t lo1 = off1;
	int64_t hi1 = off1 + len1 - 1;
	int64_t lo2 = off2;
	int64_t hi2 = off2 + len2 - 1;
	bool containment = false;
	if ((lo1 >= lo2 && hi1 <= hi2) ||
		(lo2 >= lo1 && hi2 <= hi1))
	{
		containment = true;
	}
	int type = PE_ALS_NORMAL;
	bool olap = false;
	if ((lo1 <= lo2 && hi1 >= lo2) ||
		(lo1 <= hi2 && hi1 >= hi2) ||
		containment)
	{
		olap = true;
		if (!olapOk_)
			return PE_ALS_DISCORD;
		type = PE_ALS_OVERLAP;
	}
	if (!olap)
	{
		if ((oneLeft && lo2 < lo1) || (!oneLeft && lo1 < lo2))
		{
			return PE_ALS_DISCORD;
		}
	}
	if (containment)
	{
		if (!containOk_)
			return PE_ALS_DISCORD;
		type = PE_ALS_CONTAIN;
	}
	if ((oneLeft && (hi1 > hi2 || lo2 < lo1)) ||
		(!oneLeft && (hi2 > hi1 || lo1 < lo2)))
	{
		if (!dovetailOk_)
			return PE_ALS_DISCORD;
		type = PE_ALS_DOVETAIL;
	}
	return type;
}
bool PairedEndPolicy::otherMate(
	bool is1,
	bool fw, int64_t off, int64_t maxalcols, size_t reflen, size_t len1, size_t len2, bool &oleft, int64_t &oll, int64_t &olr, int64_t &orl, int64_t &orr, bool &ofw) const
{
	assert_gt(len1, 0);
	assert_gt(len2, 0);
	assert_gt(maxfrag_, 0);
	assert_geq(minfrag_, 0);
	assert_geq(maxfrag_, minfrag_);
	assert(maxalcols == -1 || maxalcols > 0);
	pePolicyMateDir(pol_, is1, fw, oleft, ofw);
	size_t alen = is1 ? len1 : len2;
	size_t maxfrag = maxfrag_;
	size_t minfrag = minfrag_;
	if (minfrag < 1)
	{
		minfrag = 1;
	}
	if (len1 > maxfrag && expandToFit_)
		maxfrag = len1;
	if (len2 > maxfrag && expandToFit_)
		maxfrag = len2;
	if (!expandToFit_ && (len1 > maxfrag || len2 > maxfrag))
	{
		return false;
	}
	if (oleft)
	{
		oll = off + alen - maxfrag;
		olr = off + alen - minfrag;
		assert_geq(olr, oll);
		orl = oll;
		orr = off + maxfrag - 1;
		assert_geq(olr, oll);
		if (!olapOk_)
		{
			orr = min<int64_t>(orr, off - 1);
			if (orr < olr)
				olr = orr;
			assert_leq(oll, olr);
			assert_leq(orl, orr);
			assert_geq(orr, olr);
		}
		else if (!dovetailOk_)
		{
			orr = min<int64_t>(orr, off + alen - 1);
			assert_leq(oll, olr);
			assert_leq(orl, orr);
		}
		else if (!flippingOk_ && maxalcols != -1)
		{
			orr = min<int64_t>(orr, off + alen - 1 + (maxalcols - 1));
			assert_leq(oll, olr);
			assert_leq(orl, orr);
		}
		assert_geq(olr, oll);
		assert_geq(orr, orl);
		assert_geq(orr, olr);
		assert_geq(orl, oll);
	}
	else
	{
		orr = off + (maxfrag - 1);
		orl = off + (minfrag - 1);
		assert_geq(orr, orl);
		oll = off + alen - maxfrag;
		olr = orr;
		assert_geq(olr, oll);
		if (!olapOk_)
		{
			oll = max<int64_t>(oll, off + alen);
			if (oll > orl)
				orl = oll;
			assert_leq(oll, olr);
			assert_leq(orl, orr);
			assert_geq(orl, oll);
		}
		else if (!dovetailOk_)
		{
			oll = max<int64_t>(oll, off);
			assert_leq(oll, olr);
			assert_leq(orl, orr);
		}
		else if (!flippingOk_ && maxalcols != -1)
		{
			oll = max<int64_t>(oll, off - maxalcols + 1);
			assert_leq(oll, olr);
			assert_leq(orl, orr);
		}
		assert_geq(olr, oll);
		assert_geq(orr, orl);
		assert_geq(orr, olr);
		assert_geq(orl, oll);
	}
	return true;
}
#ifdef MAIN_PE
#include <string>
#include <sstream>
void testCaseClassify(
	const string &name,
	int pol,
	size_t maxfrag,
	size_t minfrag,
	bool local,
	bool flip,
	bool dove,
	bool cont,
	bool olap,
	bool expand,
	int64_t off1,
	size_t len1,
	bool fw1,
	int64_t off2,
	size_t len2,
	bool fw2,
	int expect_class)
{
	PairedEndPolicy pepol(
		pol,
		maxfrag,
		minfrag,
		local,
		flip,
		dove,
		cont,
		olap,
		expand);
	int ret = pepol.peClassifyPair(
		off1,
		len1, fw1, off2, len2, fw2);
	assert_eq(expect_class, ret);
	cout << "peClassifyPair: " << name << "...PASSED" << endl;
}
void testCaseOtherMate(
	const string &name,
	int pol,
	size_t maxfrag,
	size_t minfrag,
	bool local,
	bool flip,
	bool dove,
	bool cont,
	bool olap,
	bool expand,
	bool is1,
	bool fw,
	int64_t off,
	int64_t maxalcols,
	size_t reflen,
	size_t len1,
	size_t len2,
	bool expect_ret,
	bool expect_oleft,
	int64_t expect_oll,
	int64_t expect_olr,
	int64_t expect_orl,
	int64_t expect_orr,
	bool expect_ofw)
{
	PairedEndPolicy pepol(
		pol,
		maxfrag,
		minfrag,
		local,
		flip,
		dove,
		cont,
		olap,
		expand);
	int64_t oll = 0, olr = 0;
	int64_t orl = 0, orr = 0;
	bool oleft = false, ofw = false;
	bool ret = pepol.otherMate(
		is1,
		fw,
		off,
		maxalcols,
		reflen,
		len1,
		len2,
		oleft,
		oll,
		olr,
		orl,
		orr,
		ofw);
	assert(ret == expect_ret);
	if (ret)
	{
		assert_eq(expect_oleft, oleft);
		assert_eq(expect_oll, oll);
		assert_eq(expect_olr, olr);
		assert_eq(expect_orl, orl);
		assert_eq(expect_orr, orr);
		assert_eq(expect_ofw, ofw);
	}
	cout << "otherMate: " << name << "...PASSED" << endl;
}
int main(int argc, char **argv)
{
	{
		int policies[] = {PE_POLICY_FF, PE_POLICY_RR, PE_POLICY_FR, PE_POLICY_RF, PE_POLICY_FF, PE_POLICY_RR, PE_POLICY_FR, PE_POLICY_RF};
		bool is1[] = {true, true, true, true, false, false, false, false};
		bool fw[] = {true, false, true, false, false, true, true, false};
		bool oleft[] = {false, false, false, false, false, false, false, false};
		bool ofw[] = {true, false, false, true, false, true, false, true};
		for (int i = 0; i < 8; i++)
		{
			ostringstream oss;
			oss << "Simple";
			oss << i;
			testCaseOtherMate(
				oss.str(),
				policies[i],
				30, 20, false, true, true, true, true, true, is1[i], fw[i], 100, -1, 200, 10, 10, true, oleft[i], 80, 129, 119, 129, ofw[i]);
		}
	}
	{
		int policies[] = {PE_POLICY_FF, PE_POLICY_RR, PE_POLICY_FR, PE_POLICY_RF, PE_POLICY_FF, PE_POLICY_RR, PE_POLICY_FR, PE_POLICY_RF};
		bool is1[] = {false, false, false, false, true, true, true, true};
		bool fw[] = {true, false, false, true, false, true, false, true};
		bool oleft[] = {true, true, true, true, true, true, true, true};
		bool ofw[] = {true, false, true, false, false, true, true, false};
		for (int i = 0; i < 8; i++)
		{
			ostringstream oss;
			oss << "Simple";
			oss << (i + 8);
			testCaseOtherMate(
				oss.str(),
				policies[i],
				30, 20, false, true, true, true, true, true, is1[i], fw[i], 120, -1, 200, 10, 10, true, oleft[i], 100, 110, 100, 149, ofw[i]);
		}
	}
	testCaseOtherMate(
		"MinFragEqMax1",
		PE_POLICY_FR,
		30, 30, false, true, true, true, true, true, false, false, 120, -1, 200, 10, 10, true, true, 100, 100, 100, 149, true);
	testCaseOtherMate(
		"MinFragEqMax2",
		PE_POLICY_FR,
		30, 30, false, true, true, true, true, true, true, true, 100, -1, 200, 10, 10, true, false, 80, 129, 129, 129, false);
	testCaseOtherMate(
		"MinFragEqMax4NoDove1",
		PE_POLICY_FR,
		30, 25, false, true, false, true, true, true, true, true, 100, -1, 200, 10, 10, true, false, 100, 129, 124, 129, false);
	testCaseOtherMate(
		"MinFragEqMax4NoCont1",
		PE_POLICY_FR,
		30, 25, false, true, false, false, true, true, true, true, 100, -1, 200, 10, 10, true, false, 100, 129, 124, 129, false);
	testCaseOtherMate(
		"MinFragEqMax4NoOlap1",
		PE_POLICY_FR,
		30, 25, false, true, false, false, false, true, true, true, 100, -1, 200, 10, 10, true, false, 110, 129, 124, 129, false);
	testCaseOtherMate(
		"MinFragEqMax4NoDove2",
		PE_POLICY_FR,
		30, 25, false, true, false, true, true, true, false, false, 120, -1, 200, 10, 10, true, true, 100, 105, 100, 129, true);
	testCaseOtherMate(
		"MinFragEqMax4NoOlap2",
		PE_POLICY_FR,
		30, 25, false, true, false, false, false, true, false, false, 120, -1, 200, 10, 10, true, true, 100, 105, 100, 119, true);
	{
		int olls[] = {110};
		int olrs[] = {299};
		int orls[] = {149};
		int orrs[] = {299};
		for (int i = 0; i < 1; i++)
		{
			ostringstream oss;
			oss << "Overhang1_";
			oss << (i + 1);
			testCaseOtherMate(
				oss.str(),
				PE_POLICY_FR,
				200, 50, false, true, true, true, false, true, true, true, 100, -1, 200, 10, 10, true, false, olls[i], olrs[i], orls[i], orrs[i], false);
		}
	}
	{
		int olls[] = {-100};
		int olrs[] = {50};
		int orls[] = {-100};
		int orrs[] = {89};
		for (int i = 0; i < 1; i++)
		{
			ostringstream oss;
			oss << "Overhang2_";
			oss << (i + 1);
			testCaseOtherMate(
				oss.str(),
				PE_POLICY_FR,
				200, 50, false, true, true, true, false, true, true, false, 90, -1, 200, 10, 10, true, true, olls[i], olrs[i], orls[i], orrs[i], true);
		}
	}
	{
		int mate2offs[] = {150, 149, 149, 100, 99, 299, 1, 250, 250};
		int mate2lens[] = {50, 50, 51, 100, 101, 1, 50, 50, 51};
		int peExpects[] = {PE_ALS_NORMAL, PE_ALS_DISCORD, PE_ALS_OVERLAP, PE_ALS_CONTAIN, PE_ALS_DOVETAIL, PE_ALS_NORMAL, PE_ALS_DISCORD, PE_ALS_NORMAL, PE_ALS_DISCORD};
		for (int i = 0; i < 9; i++)
		{
			ostringstream oss;
			oss << "Simple1_";
			oss << (i);
			testCaseClassify(
				oss.str(),
				PE_POLICY_FR,
				200, 100, false, true, true, true, true, true, 100, 50, true, mate2offs[i], mate2lens[i], false, peExpects[i]);
		}
	}
	{
		int mate1offs[] = {200, 201, 200, 200, 200, 100, 400, 100, 99};
		int mate1lens[] = {50, 49, 51, 100, 101, 1, 50, 50, 51};
		int peExpects[] = {PE_ALS_NORMAL, PE_ALS_DISCORD, PE_ALS_OVERLAP, PE_ALS_CONTAIN, PE_ALS_DOVETAIL, PE_ALS_NORMAL, PE_ALS_DISCORD, PE_ALS_NORMAL, PE_ALS_DISCORD};
		for (int i = 0; i < 9; i++)
		{
			ostringstream oss;
			oss << "Simple2_";
			oss << (i);
			testCaseClassify(
				oss.str(),
				PE_POLICY_FR,
				200, 100, false, true, true, true, true, true, mate1offs[i], mate1lens[i], true, 250, 50, false, peExpects[i]);
		}
	}
	testCaseOtherMate(
		"Regression1",
		PE_POLICY_FF,
		50, 0, false, true, true, true, true, true, true, false, 3, -1, 53, 10, 10, true, true, -37, 13, -37, 52, false);
}
#endif
