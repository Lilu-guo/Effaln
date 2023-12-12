#ifndef PE_H_
#define PE_H_
#include <iostream>
#include <stdint.h>
enum
{
	PE_POLICY_FF = 1,
	PE_POLICY_RR,
	PE_POLICY_FR,
	PE_POLICY_RF
};
enum
{
	PE_ALS_NORMAL = 1,
	PE_ALS_OVERLAP,
	PE_ALS_CONTAIN,
	PE_ALS_DOVETAIL,
	PE_ALS_DISCORD
};
static inline bool pePolicyCompat(
	int policy,
	bool oneLeft, bool oneWat, bool twoWat)
{
	switch (policy)
	{
	case PE_POLICY_FF:
		return oneWat == twoWat && oneWat == oneLeft;
	case PE_POLICY_RR:
		return oneWat == twoWat && oneWat != oneLeft;
	case PE_POLICY_FR:
		return oneWat != twoWat && oneWat == oneLeft;
	case PE_POLICY_RF:
		return oneWat != twoWat && oneWat != oneLeft;
	default:
	{
		std::cerr << "Bad PE_POLICY: " << policy << std::endl;
		throw 1;
	}
	}
	throw 1;
}
static inline void pePolicyMateDir(
	int policy,
	bool is1, bool fw, bool &left, bool &mfw)
{
	switch (policy)
	{
	case PE_POLICY_FF:
	{
		left = (is1 != fw);
		mfw = fw;
		break;
	}
	case PE_POLICY_RR:
	{
		left = (is1 == fw);
		mfw = fw;
		break;
	}
	case PE_POLICY_FR:
	{
		left = !fw;
		mfw = !fw;
		break;
	}
	case PE_POLICY_RF:
	{
		left = fw;
		mfw = !fw;
		break;
	}
	default:
	{
		std::cerr << "Error: No such PE_POLICY: " << policy << std::endl;
		throw 1;
	}
	}
	return;
}
class PairedEndPolicy
{
public:
	PairedEndPolicy() { reset(); }
	PairedEndPolicy(
		int pol,
		size_t maxfrag,
		size_t minfrag,
		bool local,
		bool flippingOk,
		bool dovetailOk,
		bool containOk,
		bool olapOk,
		bool expandToFit)
	{
		init(
			pol,
			maxfrag,
			minfrag,
			local,
			flippingOk,
			dovetailOk,
			containOk,
			olapOk,
			expandToFit);
	}
	void reset()
	{
		init(-1, 0xffffffff, 0xffffffff, false, false, false, false, false, false);
	}
	void init(
		int pol,
		size_t maxfrag,
		size_t minfrag,
		bool local,
		bool flippingOk,
		bool dovetailOk,
		bool containOk,
		bool olapOk,
		bool expandToFit)
	{
		pol_ = pol;
		maxfrag_ = maxfrag;
		minfrag_ = minfrag;
		local_ = local;
		flippingOk_ = flippingOk;
		dovetailOk_ = dovetailOk;
		containOk_ = containOk;
		olapOk_ = olapOk;
		expandToFit_ = expandToFit;
	}
	bool otherMate(
		bool is1,
		bool fw, int64_t off, int64_t maxalcols, size_t reflen, size_t len1, size_t len2, bool &oleft, int64_t &oll, int64_t &olr, int64_t &orl, int64_t &orr, bool &ofw) const;
	int peClassifyPair(
		int64_t off1,
		size_t len1, bool fw1, int64_t off2, size_t len2, bool fw2) const;
	int policy() const { return pol_; }
	size_t maxFragLen() const { return maxfrag_; }
	size_t minFragLen() const { return minfrag_; }
protected:
	bool local_;
	int pol_;
	bool flippingOk_;
	bool dovetailOk_;
	bool containOk_;
	bool olapOk_;
	bool expandToFit_;
	size_t maxfrag_;
	size_t minfrag_;
};
#endif
