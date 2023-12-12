#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "aln_idx.h"
void Ebwt::sanityCheckUpToSide(TIndexOff upToSide) const
{
	assert(isInMemory());
	TIndexOffU occ[] = {0, 0, 0, 0};
	ASSERT_ONLY(TIndexOffU occ_save[] = {0, 0, 0, 0});
	TIndexOffU cur = 0;
	const EbwtParams &eh = this->_eh;
	bool fw = false;
	while (cur < (TIndexOffU)(upToSide * eh._sideSz))
	{
		assert_leq(cur + eh._sideSz, eh._ebwtTotLen);
		for (uint32_t i = 0; i < eh._sideBwtSz; i++)
		{
			uint8_t by = this->ebwt()[cur + (fw ? i : eh._sideBwtSz - i - 1)];
			for (int j = 0; j < 4; j++)
			{
				int twoBit = unpack_2b_from_8b(by, fw ? j : 3 - j);
				occ[twoBit]++;
			}
			assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3]) % 4);
		}
		assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3]) % eh._sideBwtLen);
		ASSERT_ONLY(const TIndexOffU *u_ebwt = reinterpret_cast<const TIndexOffU *>(&ebwt()[cur + eh._sideBwtSz]));
		ASSERT_ONLY(TIndexOffU as = u_ebwt[0]);
		ASSERT_ONLY(TIndexOffU cs = u_ebwt[1]);
		ASSERT_ONLY(TIndexOffU gs = u_ebwt[2]);
		ASSERT_ONLY(TIndexOffU ts = u_ebwt[3]);
		assert(as == occ_save[0] || as == occ_save[0] - 1);
		assert_eq(cs, occ_save[1]);
		assert_eq(gs, occ_save[2]);
		assert_eq(ts, occ_save[3]);
#ifndef NDEBUG
		occ_save[0] = occ[0];
		occ_save[1] = occ[1];
		occ_save[2] = occ[2];
		occ_save[3] = occ[3];
#endif
		cur += eh._sideSz;
	}
}
void Ebwt::sanityCheckAll(int reverse) const
{
	const EbwtParams &eh = this->_eh;
	assert(isInMemory());
	for (TIndexOffU i = 1; i < eh._ftabLen; i++)
	{
		assert_geq(this->ftabHi(i), this->ftabLo(i - 1));
		assert_geq(this->ftabLo(i), this->ftabHi(i - 1));
		assert_leq(this->ftabHi(i), eh._bwtLen + 1);
	}
	assert_eq(this->ftabHi(eh._ftabLen - 1), eh._bwtLen);
	TIndexOff seenLen = (eh._bwtLen + 31) >> ((TIndexOffU)5);
	TIndexOff *seen;
	try
	{
		seen = new TIndexOff[seenLen];
	}
	catch (bad_alloc &e)
	{
		cerr << "Out of memory allocating seen[] at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}
	memset(seen, 0, OFF_SIZE * seenLen);
	TIndexOffU offsLen = eh._offsLen;
	for (TIndexOffU i = 0; i < offsLen; i++)
	{
		assert_lt(this->offs()[i], eh._bwtLen);
		TIndexOff w = this->offs()[i] >> 5;
		TIndexOff r = this->offs()[i] & 31;
		assert_eq(0, (seen[w] >> r) & 1);
		seen[w] |= (1 << r);
	}
	delete[] seen;
	assert_gt(this->_nPat, 0);
	for (TIndexOffU i = 0; i < this->_nPat; i++)
	{
		assert_geq(this->plen()[i], 0);
	}
	if (this->rstarts() != NULL)
	{
		for (TIndexOffU i = 0; i < this->_nFrag - 1; i++)
		{
			assert_gt(this->rstarts()[(i + 1) * 3], this->rstarts()[i * 3]);
			if (reverse == REF_READ_REVERSE)
			{
				assert(this->rstarts()[(i * 3) + 1] >= this->rstarts()[((i + 1) * 3) + 1]);
			}
			else
			{
				assert(this->rstarts()[(i * 3) + 1] <= this->rstarts()[((i + 1) * 3) + 1]);
			}
		}
	}
	sanityCheckUpToSide(eh._numSides);
	VMSG_NL("Ebwt::sanityCheck passed");
}
void Ebwt::restore(SString<char> &s) const
{
	assert(isInMemory());
	s.resize(this->_eh._len);
	TIndexOffU jumps = 0;
	TIndexOffU i = this->_eh._len;
	SideLocus l(i, this->_eh, this->ebwt());
	while (i != _zOff)
	{
		assert_lt(jumps, this->_eh._len);
		TIndexOffU newi = mapLF(l ASSERT_ONLY(, false));
		assert_neq(newi, i);
		s[this->_eh._len - jumps - 1] = rowL(l);
		i = newi;
		l.initFromRow(i, this->_eh, this->ebwt());
		jumps++;
	}
	assert_eq(jumps, this->_eh._len);
}
void Ebwt::checkOrigs(
	const EList<SString<char>> &os,
	bool color,
	bool mirror) const
{
	SString<char> rest;
	restore(rest);
	TIndexOffU restOff = 0;
	size_t i = 0, j = 0;
	if (mirror)
	{
		return;
	}
	while (i < os.size())
	{
		size_t olen = os[i].length();
		int lastorig = -1;
		for (; j < olen; j++)
		{
			size_t joff = j;
			if (mirror)
				joff = olen - j - 1;
			if ((int)os[i][joff] == 4)
			{
				lastorig = -1;
				if (!mirror)
				{
					while (j < olen && (int)os[i][j] == 4)
						j++;
				}
				else
				{
					while (j < olen && (int)os[i][olen - j - 1] == 4)
						j++;
				}
				j--;
				continue;
			}
			if (lastorig == -1 && color)
			{
				lastorig = os[i][joff];
				continue;
			}
			if (color)
			{
				assert_neq(-1, lastorig);
				assert_eq(dinuc2color[(int)os[i][joff]][lastorig], rest[restOff]);
			}
			else
			{
				assert_eq(os[i][joff], rest[restOff]);
			}
			lastorig = (int)os[i][joff];
			restOff++;
		}
		if (j == os[i].length())
		{
			i++;
			j = 0;
		}
		else
		{
		}
	}
}
