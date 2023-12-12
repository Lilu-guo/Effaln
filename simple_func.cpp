#include <iostream>
#include "simple_func.h"
#include "ds.h"
#include "mem_ids.h"
int SimpleFunc::parseType(const std::string &otype)
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
SimpleFunc SimpleFunc::parse(
	const std::string &s,
	double defaultConst,
	double defaultLinear,
	double defaultMin,
	double defaultMax)
{
	EList<string> ctoks(MISC_CAT);
	string ctok;
	istringstream css(s);
	SimpleFunc fv;
	while (getline(css, ctok, ','))
	{
		ctoks.push_back(ctok);
	}
	if (ctoks.size() >= 1)
	{
		fv.setType(parseType(ctoks[0]));
	}
	if (ctoks.size() >= 2)
	{
		double co;
		istringstream tmpss(ctoks[1]);
		tmpss >> co;
		fv.setConst(co);
	}
	else
	{
		fv.setConst(defaultConst);
	}
	if (ctoks.size() >= 3)
	{
		double ce;
		istringstream tmpss(ctoks[2]);
		tmpss >> ce;
		fv.setCoeff(ce);
	}
	else
	{
		fv.setCoeff(defaultLinear);
	}
	if (ctoks.size() >= 4)
	{
		double mn;
		istringstream tmpss(ctoks[3]);
		tmpss >> mn;
		fv.setMin(mn);
	}
	else
	{
		fv.setMin(defaultMin);
	}
	if (ctoks.size() >= 5)
	{
		double mx;
		istringstream tmpss(ctoks[4]);
		tmpss >> mx;
		fv.setMax(mx);
	}
	else
	{
		fv.setMax(defaultMax);
	}
	return fv;
}
