#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include "ds.h"
#include "aligner_seed_policy.h"
#include "mem_ids.h"
using namespace std;
static int parseFuncType(const std::string &otype)
{
	string type = otype;
	if (type == "C" || type == "Constant")
	{
		return SIMPLE_FUNC_CONST;
	}
	else if (type == "L" || type == "Linear")
	{
		return SIMPLE_FUNC_LINEAR;
	}
	else if (type == "S" || type == "Sqrt")
	{
		return SIMPLE_FUNC_SQRT;
	}
	else if (type == "G" || type == "Log")
	{
		return SIMPLE_FUNC_LOG;
	}
	std::cerr << "Error: Bad function type '" << otype.c_str()
			  << "'.  Should be C (constant), L (linear), "
			  << "S (square root) or G (natural log)." << std::endl;
	throw 1;
}
#define PARSE_FUNC(fv)                           \
	{                                            \
		if (ctoks.size() >= 1)                   \
		{                                        \
			fv.setType(parseFuncType(ctoks[0])); \
		}                                        \
		if (ctoks.size() >= 2)                   \
		{                                        \
			double co;                           \
			istringstream tmpss(ctoks[1]);       \
			tmpss >> co;                         \
			fv.setConst(co);                     \
		}                                        \
		if (ctoks.size() >= 3)                   \
		{                                        \
			double ce;                           \
			istringstream tmpss(ctoks[2]);       \
			tmpss >> ce;                         \
			fv.setCoeff(ce);                     \
		}                                        \
		if (ctoks.size() >= 4)                   \
		{                                        \
			double mn;                           \
			istringstream tmpss(ctoks[3]);       \
			tmpss >> mn;                         \
			fv.setMin(mn);                       \
		}                                        \
		if (ctoks.size() >= 5)                   \
		{                                        \
			double mx;                           \
			istringstream tmpss(ctoks[4]);       \
			tmpss >> mx;                         \
			fv.setMin(mx);                       \
		}                                        \
	}
void SeedAlignmentPolicy::parseString(
	const std::string &s,
	bool local,
	bool noisyHpolymer,
	bool ignoreQuals,
	int &bonusMatchType,
	int &bonusMatch,
	int &penMmcType,
	int &penMmcMax,
	int &penMmcMin,
	int &penNType,
	int &penN,
	int &penRdExConst,
	int &penRfExConst,
	int &penRdExLinear,
	int &penRfExLinear,
	SimpleFunc &costMin,
	SimpleFunc &nCeil,
	bool &nCatPair,
	int &multiseedMms,
	int &multiseedLen,
	SimpleFunc &multiseedIval,
	size_t &failStreak,
	size_t &seedRounds)
{
	bonusMatchType = local ? DEFAULT_MATCH_BONUS_TYPE_LOCAL : DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch = local ? DEFAULT_MATCH_BONUS_LOCAL : DEFAULT_MATCH_BONUS;
	penMmcType = ignoreQuals ? DEFAULT_MM_PENALTY_TYPE_IGNORE_QUALS : DEFAULT_MM_PENALTY_TYPE;
	penMmcMax = DEFAULT_MM_PENALTY_MAX;
	penMmcMin = DEFAULT_MM_PENALTY_MIN;
	penNType = DEFAULT_N_PENALTY_TYPE;
	penN = DEFAULT_N_PENALTY;
	const double DMAX = std::numeric_limits<double>::max();
	costMin.init(
		local ? SIMPLE_FUNC_LOG : SIMPLE_FUNC_LINEAR,
		local ? DEFAULT_MIN_CONST_LOCAL : DEFAULT_MIN_CONST,
		local ? DEFAULT_MIN_LINEAR_LOCAL : DEFAULT_MIN_LINEAR);
	nCeil.init(
		SIMPLE_FUNC_LINEAR, 0.0f, DMAX,
		DEFAULT_N_CEIL_CONST, DEFAULT_N_CEIL_LINEAR);
	multiseedIval.init(
		DEFAULT_IVAL, 1.0f, DMAX,
		DEFAULT_IVAL_B, DEFAULT_IVAL_A);
	nCatPair = DEFAULT_N_CAT_PAIR;
	if (!noisyHpolymer)
	{
		penRdExConst = DEFAULT_READ_GAP_CONST;
		penRdExLinear = DEFAULT_READ_GAP_LINEAR;
		penRfExConst = DEFAULT_REF_GAP_CONST;
		penRfExLinear = DEFAULT_REF_GAP_LINEAR;
	}
	else
	{
		penRdExConst = DEFAULT_READ_GAP_CONST_BADHPOLY;
		penRdExLinear = DEFAULT_READ_GAP_LINEAR_BADHPOLY;
		penRfExConst = DEFAULT_REF_GAP_CONST_BADHPOLY;
		penRfExLinear = DEFAULT_REF_GAP_LINEAR_BADHPOLY;
	}
	multiseedMms = DEFAULT_SEEDMMS;
	multiseedLen = gDefaultSeedLen;
	EList<string> toks(MISC_CAT);
	string tok;
	istringstream ss(s);
	int setting = 0;
	while (getline(ss, tok, ';'))
	{
		setting++;
		EList<string> etoks(MISC_CAT);
		string etok;
		istringstream ess(tok);
		while (getline(ess, etok, '='))
		{
			etoks.push_back(etok);
		}
		if (etoks.size() != 2)
		{
			cerr << "Error parsing alignment policy setting " << setting
				 << "; must be bisected by = sign" << endl
				 << "Policy: " << s.c_str() << endl;
			assert(false);
			throw 1;
		}
		string tag = etoks[0], val = etoks[1];
		EList<string> ctoks(MISC_CAT);
		string ctok;
		istringstream css(val);
		while (getline(css, ctok, ','))
		{
			ctoks.push_back(ctok);
		}
		if (ctoks.size() == 0)
		{
			cerr << "Error parsing alignment policy setting " << setting
				 << "; RHS must have at least 1 token" << endl
				 << "Policy: " << s.c_str() << endl;
			assert(false);
			throw 1;
		}
		for (size_t i = 0; i < ctoks.size(); i++)
		{
			if (ctoks[i].length() == 0)
			{
				cerr << "Error parsing alignment policy setting " << setting
					 << "; token " << i + 1 << " on RHS had length=0" << endl
					 << "Policy: " << s.c_str() << endl;
				assert(false);
				throw 1;
			}
		}
		if (tag == "MA")
		{
			if (ctoks.size() != 1)
			{
				cerr << "Error parsing alignment policy setting " << setting
					 << "; RHS must have 1 token" << endl
					 << "Policy: " << s.c_str() << endl;
				assert(false);
				throw 1;
			}
			string tmp = ctoks[0];
			istringstream tmpss(tmp);
			tmpss >> bonusMatch;
		}
		else if (tag == "MMP")
		{
			if (ctoks.size() > 3)
			{
				cerr << "Error parsing alignment policy setting "
					 << "'" << tag.c_str() << "'"
					 << "; RHS must have at most 3 tokens" << endl
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false);
				throw 1;
			}
			if (ctoks[0][0] == 'C')
			{
				string tmp = ctoks[0].substr(1);
				istringstream tmpss(tmp);
				tmpss >> penMmcMax;
				penMmcMin = penMmcMax;
				penMmcType = COST_MODEL_CONSTANT;
			}
			else if (ctoks[0][0] == 'Q')
			{
				if (ctoks.size() >= 2)
				{
					string tmp = ctoks[1];
					istringstream tmpss(tmp);
					tmpss >> penMmcMax;
				}
				else
				{
					penMmcMax = DEFAULT_MM_PENALTY_MAX;
				}
				if (ctoks.size() >= 3)
				{
					string tmp = ctoks[2];
					istringstream tmpss(tmp);
					tmpss >> penMmcMin;
				}
				else
				{
					penMmcMin = DEFAULT_MM_PENALTY_MIN;
				}
				if (penMmcMin > penMmcMax)
				{
					cerr << "Error: Maximum mismatch penalty (" << penMmcMax
						 << ") is less than minimum penalty (" << penMmcMin
						 << endl;
					throw 1;
				}
				penMmcType = COST_MODEL_QUAL;
				if (ignoreQuals)
				{
					if (gVerbose)
						cerr << "Changing MMP=Q," << penMmcMax << " to ";
					penMmcMin = penMmcMax;
					penMmcType = COST_MODEL_CONSTANT;
					if (gVerbose)
					{
						cerr << "MMP=C," << penMmcMax
							 << " because of --ignore-quals"
							 << endl;
					}
				}
			}
			else if (ctoks[0][0] == 'R')
			{
				penMmcType = COST_MODEL_ROUNDED_QUAL;
			}
			else
			{
				cerr << "Error parsing alignment policy setting "
					 << "'" << tag.c_str() << "'"
					 << "; RHS must start with C, Q or R" << endl
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false);
				throw 1;
			}
		}
		else if (tag == "NP")
		{
			if (ctoks.size() != 1)
			{
				cerr << "Error parsing alignment policy setting "
					 << "'" << tag.c_str() << "'"
					 << "; RHS must have 1 token" << endl
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false);
				throw 1;
			}
			if (ctoks[0][0] == 'C')
			{
				string tmp = ctoks[0].substr(1);
				istringstream tmpss(tmp);
				tmpss >> penN;
				penNType = COST_MODEL_CONSTANT;
			}
			else if (ctoks[0][0] == 'Q')
			{
				penNType = COST_MODEL_QUAL;
			}
			else if (ctoks[0][0] == 'R')
			{
				penNType = COST_MODEL_ROUNDED_QUAL;
			}
			else
			{
				cerr << "Error parsing alignment policy setting "
					 << "'" << tag.c_str() << "'"
					 << "; RHS must start with C, Q or R" << endl
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false);
				throw 1;
			}
		}
		else if (tag == "RDG")
		{
			if (ctoks.size() >= 1)
			{
				istringstream tmpss(ctoks[0]);
				tmpss >> penRdExConst;
			}
			else
			{
				penRdExConst = noisyHpolymer ? DEFAULT_READ_GAP_CONST_BADHPOLY : DEFAULT_READ_GAP_CONST;
			}
			if (ctoks.size() >= 2)
			{
				istringstream tmpss(ctoks[1]);
				tmpss >> penRdExLinear;
			}
			else
			{
				penRdExLinear = noisyHpolymer ? DEFAULT_READ_GAP_LINEAR_BADHPOLY : DEFAULT_READ_GAP_LINEAR;
			}
		}
		else if (tag == "RFG")
		{
			if (ctoks.size() >= 1)
			{
				istringstream tmpss(ctoks[0]);
				tmpss >> penRfExConst;
			}
			else
			{
				penRfExConst = noisyHpolymer ? DEFAULT_REF_GAP_CONST_BADHPOLY : DEFAULT_REF_GAP_CONST;
			}
			if (ctoks.size() >= 2)
			{
				istringstream tmpss(ctoks[1]);
				tmpss >> penRfExLinear;
			}
			else
			{
				penRfExLinear = noisyHpolymer ? DEFAULT_REF_GAP_LINEAR_BADHPOLY : DEFAULT_REF_GAP_LINEAR;
			}
		}
		else if (tag == "MIN")
		{
			PARSE_FUNC(costMin);
		}
		else if (tag == "NCEIL")
		{
			PARSE_FUNC(nCeil);
		}
		else if (tag == "SEED")
		{
			if (ctoks.size() > 1)
			{
				cerr << "Error parsing alignment policy setting "
					 << "'" << tag.c_str() << "'; RHS must have 1 token, "
					 << "had " << ctoks.size() << ".  "
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false);
				throw 1;
			}
			if (ctoks.size() >= 1)
			{
				istringstream tmpss(ctoks[0]);
				tmpss >> multiseedMms;
				if (multiseedMms > 1)
				{
					cerr << "Error: -N was set to " << multiseedMms << ", but cannot be set greater than 1" << endl;
					throw 1;
				}
				if (multiseedMms < 0)
				{
					cerr << "Error: -N was set to a number less than 0 (" << multiseedMms << ")" << endl;
					throw 1;
				}
			}
		}
		else if (tag == "SEEDLEN")
		{
			if (ctoks.size() > 1)
			{
				cerr << "Error parsing alignment policy setting "
					 << "'" << tag.c_str() << "'; RHS must have 1 token, "
					 << "had " << ctoks.size() << ".  "
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false);
				throw 1;
			}
			if (ctoks.size() >= 1)
			{
				istringstream tmpss(ctoks[0]);
				tmpss >> multiseedLen;
			}
		}
		else if (tag == "DPS")
		{
			if (ctoks.size() > 1)
			{
				cerr << "Error parsing alignment policy setting "
					 << "'" << tag.c_str() << "'; RHS must have 1 token, "
					 << "had " << ctoks.size() << ".  "
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false);
				throw 1;
			}
			if (ctoks.size() >= 1)
			{
				istringstream tmpss(ctoks[0]);
				tmpss >> failStreak;
			}
		}
		else if (tag == "ROUNDS")
		{
			if (ctoks.size() > 1)
			{
				cerr << "Error parsing alignment policy setting "
					 << "'" << tag.c_str() << "'; RHS must have 1 token, "
					 << "had " << ctoks.size() << ".  "
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false);
				throw 1;
			}
			if (ctoks.size() >= 1)
			{
				istringstream tmpss(ctoks[0]);
				tmpss >> seedRounds;
			}
		}
		else if (tag == "IVAL")
		{
			PARSE_FUNC(multiseedIval);
		}
		else
		{
			cerr << "Unexpected alignment policy setting "
				 << "'" << tag.c_str() << "'" << endl
				 << "Policy: '" << s.c_str() << "'" << endl;
			assert(false);
			throw 1;
		}
	}
}
#ifdef ALIGNER_SEED_POLICY_MAIN
int main()
{
	int bonusMatchType;
	int bonusMatch;
	int penMmcType;
	int penMmc;
	int penNType;
	int penN;
	int penRdExConst;
	int penRfExConst;
	int penRdExLinear;
	int penRfExLinear;
	SimpleFunc costMin;
	SimpleFunc costFloor;
	SimpleFunc nCeil;
	bool nCatPair;
	int multiseedMms;
	int multiseedLen;
	SimpleFunc msIval;
	SimpleFunc posfrac;
	SimpleFunc rowmult;
	uint32_t mhits;
	{
		cout << "Case 1: Defaults 1 ... ";
		const char *pol = "";
		SeedAlignmentPolicy::parseString(
			string(pol),
			false,
			false,
			false,
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmc,
			penNType,
			penN,
			penRdExConst,
			penRfExConst,
			penRdExLinear,
			penRfExLinear,
			costMin,
			costFloor,
			nCeil,
			nCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			mhits);
		assert_eq(DEFAULT_MATCH_BONUS_TYPE, bonusMatchType);
		assert_eq(DEFAULT_MATCH_BONUS, bonusMatch);
		assert_eq(DEFAULT_MM_PENALTY_TYPE, penMmcType);
		assert_eq(DEFAULT_MM_PENALTY_MAX, penMmcMax);
		assert_eq(DEFAULT_MM_PENALTY_MIN, penMmcMin);
		assert_eq(DEFAULT_N_PENALTY_TYPE, penNType);
		assert_eq(DEFAULT_N_PENALTY, penN);
		assert_eq(DEFAULT_MIN_CONST, costMin.getConst());
		assert_eq(DEFAULT_MIN_LINEAR, costMin.getCoeff());
		assert_eq(DEFAULT_FLOOR_CONST, costFloor.getConst());
		assert_eq(DEFAULT_FLOOR_LINEAR, costFloor.getCoeff());
		assert_eq(DEFAULT_N_CEIL_CONST, nCeil.getConst());
		assert_eq(DEFAULT_N_CAT_PAIR, nCatPair);
		assert_eq(DEFAULT_READ_GAP_CONST, penRdExConst);
		assert_eq(DEFAULT_READ_GAP_LINEAR, penRdExLinear);
		assert_eq(DEFAULT_REF_GAP_CONST, penRfExConst);
		assert_eq(DEFAULT_REF_GAP_LINEAR, penRfExLinear);
		assert_eq(DEFAULT_SEEDMMS, multiseedMms);
		assert_eq(DEFAULT_IVAL, msIval.getType());
		assert_eq(DEFAULT_IVAL_A, msIval.getCoeff());
		assert_eq(DEFAULT_IVAL_B, msIval.getConst());
		cout << "PASSED" << endl;
	}
	{
		cout << "Case 2: Defaults 2 ... ";
		const char *pol = "";
		SeedAlignmentPolicy::parseString(
			string(pol),
			false,
			true,
			false,
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmc,
			penNType,
			penN,
			penRdExConst,
			penRfExConst,
			penRdExLinear,
			penRfExLinear,
			costMin,
			costFloor,
			nCeil,
			nCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			mhits);
		assert_eq(DEFAULT_MATCH_BONUS_TYPE, bonusMatchType);
		assert_eq(DEFAULT_MATCH_BONUS, bonusMatch);
		assert_eq(DEFAULT_MM_PENALTY_TYPE, penMmcType);
		assert_eq(DEFAULT_MM_PENALTY_MAX, penMmc);
		assert_eq(DEFAULT_MM_PENALTY_MIN, penMmc);
		assert_eq(DEFAULT_N_PENALTY_TYPE, penNType);
		assert_eq(DEFAULT_N_PENALTY, penN);
		assert_eq(DEFAULT_MIN_CONST, costMin.getConst());
		assert_eq(DEFAULT_MIN_LINEAR, costMin.getCoeff());
		assert_eq(DEFAULT_FLOOR_CONST, costFloor.getConst());
		assert_eq(DEFAULT_FLOOR_LINEAR, costFloor.getCoeff());
		assert_eq(DEFAULT_N_CEIL_CONST, nCeil.getConst());
		assert_eq(DEFAULT_N_CAT_PAIR, nCatPair);
		assert_eq(DEFAULT_READ_GAP_CONST_BADHPOLY, penRdExConst);
		assert_eq(DEFAULT_READ_GAP_LINEAR_BADHPOLY, penRdExLinear);
		assert_eq(DEFAULT_REF_GAP_CONST_BADHPOLY, penRfExConst);
		assert_eq(DEFAULT_REF_GAP_LINEAR_BADHPOLY, penRfExLinear);
		assert_eq(DEFAULT_SEEDMMS, multiseedMms);
		assert_eq(DEFAULT_IVAL, msIval.getType());
		assert_eq(DEFAULT_IVAL_A, msIval.getCoeff());
		assert_eq(DEFAULT_IVAL_B, msIval.getConst());
		cout << "PASSED" << endl;
	}
	{
		cout << "Case 3: Defaults 3 ... ";
		const char *pol = "";
		SeedAlignmentPolicy::parseString(
			string(pol),
			true,
			false,
			false,
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmc,
			penNType,
			penN,
			penRdExConst,
			penRfExConst,
			penRdExLinear,
			penRfExLinear,
			costMin,
			costFloor,
			nCeil,
			nCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			mhits);
		assert_eq(DEFAULT_MATCH_BONUS_TYPE_LOCAL, bonusMatchType);
		assert_eq(DEFAULT_MATCH_BONUS_LOCAL, bonusMatch);
		assert_eq(DEFAULT_MM_PENALTY_TYPE, penMmcType);
		assert_eq(DEFAULT_MM_PENALTY_MAX, penMmcMax);
		assert_eq(DEFAULT_MM_PENALTY_MIN, penMmcMin);
		assert_eq(DEFAULT_N_PENALTY_TYPE, penNType);
		assert_eq(DEFAULT_N_PENALTY, penN);
		assert_eq(DEFAULT_MIN_CONST_LOCAL, costMin.getConst());
		assert_eq(DEFAULT_MIN_LINEAR_LOCAL, costMin.getCoeff());
		assert_eq(DEFAULT_FLOOR_CONST_LOCAL, costFloor.getConst());
		assert_eq(DEFAULT_FLOOR_LINEAR_LOCAL, costFloor.getCoeff());
		assert_eq(DEFAULT_N_CEIL_CONST, nCeil.getConst());
		assert_eq(DEFAULT_N_CEIL_LINEAR, nCeil.getCoeff());
		assert_eq(DEFAULT_N_CAT_PAIR, nCatPair);
		assert_eq(DEFAULT_READ_GAP_CONST, penRdExConst);
		assert_eq(DEFAULT_READ_GAP_LINEAR, penRdExLinear);
		assert_eq(DEFAULT_REF_GAP_CONST, penRfExConst);
		assert_eq(DEFAULT_REF_GAP_LINEAR, penRfExLinear);
		assert_eq(DEFAULT_SEEDMMS, multiseedMms);
		assert_eq(DEFAULT_IVAL, msIval.getType());
		assert_eq(DEFAULT_IVAL_A, msIval.getCoeff());
		assert_eq(DEFAULT_IVAL_B, msIval.getConst());
		cout << "PASSED" << endl;
	}
	{
		cout << "Case 4: Simple string 1 ... ";
		const char *pol = "MMP=C44;MA=4;RFG=24,12;FL=C,8;RDG=2;NP=C4;MIN=C,7";
		SeedAlignmentPolicy::parseString(
			string(pol),
			true,
			false,
			false,
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmc,
			penNType,
			penN,
			penRdExConst,
			penRfExConst,
			penRdExLinear,
			penRfExLinear,
			costMin,
			costFloor,
			nCeil,
			nCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			mhits);
		assert_eq(COST_MODEL_CONSTANT, bonusMatchType);
		assert_eq(4, bonusMatch);
		assert_eq(COST_MODEL_CONSTANT, penMmcType);
		assert_eq(44, penMmc);
		assert_eq(COST_MODEL_CONSTANT, penNType);
		assert_eq(4.0f, penN);
		assert_eq(7, costMin.getConst());
		assert_eq(DEFAULT_MIN_LINEAR_LOCAL, costMin.getCoeff());
		assert_eq(8, costFloor.getConst());
		assert_eq(DEFAULT_FLOOR_LINEAR_LOCAL, costFloor.getCoeff());
		assert_eq(DEFAULT_N_CEIL_CONST, nCeil.getConst());
		assert_eq(DEFAULT_N_CEIL_LINEAR, nCeil.getCoeff());
		assert_eq(DEFAULT_N_CAT_PAIR, nCatPair);
		assert_eq(2.0f, penRdExConst);
		assert_eq(DEFAULT_READ_GAP_LINEAR, penRdExLinear);
		assert_eq(24.0f, penRfExConst);
		assert_eq(12.0f, penRfExLinear);
		assert_eq(DEFAULT_SEEDMMS, multiseedMms);
		assert_eq(DEFAULT_IVAL, msIval.getType());
		assert_eq(DEFAULT_IVAL_A, msIval.getCoeff());
		assert_eq(DEFAULT_IVAL_B, msIval.getConst());
		cout << "PASSED" << endl;
	}
}
#endif
