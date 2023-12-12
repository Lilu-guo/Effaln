#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "aln_idx.h"
using namespace std;
#ifdef _64BIT_INDEX
const std::string gEbwt_ext("indexL");
#else
const std::string gEbwt_ext("index");
#endif
string gLastIOErrMsg;
unsigned long long int countLF_Static = 0;
unsigned long long int exact_Nums = 0;
unsigned long long int onemm_Nums = 0;
unsigned long long int allseed_Nums = 0;
void Ebwt::joinedToTextOff(
	TIndexOffU qlen,
	TIndexOffU off,
	TIndexOffU &tidx,
	TIndexOffU &textoff,
	TIndexOffU &tlen,
	bool rejectStraddle,
	bool &straddled) const
{
	assert(rstarts() != NULL);
	TIndexOffU top = 0;
	TIndexOffU bot = _nFrag;
	TIndexOffU elt = OFF_MASK;
	while (true)
	{
		ASSERT_ONLY(TIndexOffU oldelt = elt);
		elt = top + ((bot - top) >> 1);
		assert_neq(oldelt, elt);
		TIndexOffU lower = rstarts()[elt * 3];
		TIndexOffU upper;
		if (elt == _nFrag - 1)
		{
			upper = _eh._len;
		}
		else
		{
			upper = rstarts()[((elt + 1) * 3)];
		}
		assert_gt(upper, lower);
		TIndexOffU fraglen = upper - lower;
		if (lower <= off)
		{
			if (upper > off)
			{
				if (off + qlen > upper)
				{
					straddled = true;
					if (rejectStraddle)
					{
						tidx = OFF_MASK;
						assert_lt(elt, _nFrag - 1);
						return;
					}
				}
				tidx = rstarts()[(elt * 3) + 1];
				assert_lt(tidx, this->_nPat);
				assert_leq(fraglen, this->plen()[tidx]);
				TIndexOffU fragoff = off - rstarts()[(elt * 3)];
				if (!this->fw_)
				{
					fragoff = fraglen - fragoff - 1;
					fragoff -= (qlen - 1);
				}
				textoff = fragoff + rstarts()[(elt * 3) + 2];
				assert_lt(textoff, this->plen()[tidx]);
				break;
			}
			else
			{
				top = elt;
			}
		}
		else
		{
			bot = elt;
		}
	}
	tlen = this->plen()[tidx];
}
TIndexOffU Ebwt::walkLeft(TIndexOffU row, TIndexOffU steps) const
{
	assert(offs() != NULL);
	assert_neq(OFF_MASK, row);
	SideLocus l;
	if (steps > 0)
		l.initFromRow(row, _eh, ebwt());
	while (steps > 0)
	{
		if (row == _zOff)
			return OFF_MASK;
		TIndexOffU newrow = this->mapLF(l ASSERT_ONLY(, false));
		assert_neq(OFF_MASK, newrow);
		assert_neq(newrow, row);
		row = newrow;
		steps--;
		if (steps > 0)
			l.initFromRow(row, _eh, ebwt());
	}
	return row;
}
TIndexOffU Ebwt::getOffset(TIndexOffU row) const
{
	assert(offs() != NULL);
	assert_neq(OFF_MASK, row);
	if (row == _zOff)
		return 0;
	if ((row & _eh._offMask) == row)
		return this->offs()[row >> _eh._offRate];
	TIndexOffU jumps = 0;
	SideLocus l;
	l.initFromRow(row, _eh, ebwt());
	while (true)
	{
		TIndexOffU newrow = this->mapLF(l ASSERT_ONLY(, false));
		jumps++;
		assert_neq(OFF_MASK, newrow);
		assert_neq(newrow, row);
		row = newrow;
		if (row == _zOff)
		{
			cout << "碰到$符的位置..." << jumps << endl;
			return jumps;
		}
		else if ((row & _eh._offMask) == row)
		{
			return jumps + this->offs()[row >> _eh._offRate];
		}
		l.initFromRow(row, _eh, ebwt());
	}
}
TIndexOffU Ebwt::getOffset(
	TIndexOffU elt,
	bool fw,
	TIndexOffU hitlen) const
{
	TIndexOffU off = getOffset(elt);
	assert_neq(OFF_MASK, off);
	if (!fw)
	{
		assert_lt(off, _eh._len);
		off = _eh._len - off - 1;
		assert_geq(off, hitlen - 1);
		off -= (hitlen - 1);
		assert_lt(off, _eh._len);
	}
	return off;
}
bool Ebwt::contains(
	const BTDnaString &str,
	TIndexOffU *otop,
	TIndexOffU *obot) const
{
	assert(isInMemory());
	SideLocus tloc, bloc;
	if (str.empty())
	{
		if (otop != NULL && obot != NULL)
			*otop = *obot = 0;
		return true;
	}
	int c = str[str.length() - 1];
	assert_range(0, 4, c);
	TIndexOffU top = 0, bot = 0;
	if (c < 4)
	{
		top = fchr()[c];
		bot = fchr()[c + 1];
	}
	else
	{
		bool set = false;
		for (int i = 0; i < 4; i++)
		{
			if (fchr()[c] < fchr()[c + 1])
			{
				if (set)
				{
					return false;
				}
				else
				{
					set = true;
					top = fchr()[c];
					bot = fchr()[c + 1];
				}
			}
		}
	}
	assert_geq(bot, top);
	tloc.initFromRow(top, eh(), ebwt());
	bloc.initFromRow(bot, eh(), ebwt());
	ASSERT_ONLY(TIndexOffU lastDiff = bot - top);
	for (TIndexOff i = (TIndexOff)str.length() - 2; i >= 0; i--)
	{
		c = str[i];
		assert_range(0, 4, c);
		if (c <= 3)
		{
			top = mapLF(tloc, c);
			bot = mapLF(bloc, c);
		}
		else
		{
			TIndexOffU sz = bot - top;
			int c1 = mapLF1(top, tloc ASSERT_ONLY(, false));
			bot = mapLF(bloc, c1);
			assert_leq(bot - top, sz);
			if (bot - top < sz)
			{
				return false;
			}
		}
		assert_geq(bot, top);
		assert_leq(bot - top, lastDiff);
		ASSERT_ONLY(lastDiff = bot - top);
		if (i > 0)
		{
			tloc.initFromRow(top, eh(), ebwt());
			bloc.initFromRow(bot, eh(), ebwt());
		}
	}
	if (otop != NULL && obot != NULL)
	{
		*otop = top;
		*obot = bot;
	}
	return bot > top;
}
string adjustEbwtBase(const string &cmdline,
					  const string &ebwtFileBase,
					  bool verbose = false)
{
	string str = ebwtFileBase;
	ifstream in;
	if (verbose)
		cout << "Trying " << str.c_str() << endl;
	in.open((str + ".1." + gEbwt_ext).c_str(), ios_base::in | ios::binary);
	if (!in.is_open())
	{
		if (verbose)
			cout << "  didn't work" << endl;
		in.close();
		if (getenv("_INDEXES") != NULL)
		{
			str = string(getenv("_INDEXES")) + "/" + ebwtFileBase;
			if (verbose)
				cout << "Trying " << str.c_str() << endl;
			in.open((str + ".1." + gEbwt_ext).c_str(), ios_base::in | ios::binary);
			if (!in.is_open())
			{
				if (verbose)
					cout << "  didn't work" << endl;
				in.close();
			}
			else
			{
				if (verbose)
					cout << "  worked" << endl;
			}
		}
	}
	if (!in.is_open())
	{
		cerr << "Could not locate a Effaln index corresponding to basename \"" << ebwtFileBase.c_str() << "\"" << endl;
		throw 1;
	}
	return str;
}
