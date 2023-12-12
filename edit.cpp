#include <iostream>
#include "edit.h"
using namespace std;
ostream &operator<<(ostream &os, const Edit &e)
{
	os << e.pos << ":" << (char)e.chr << ">" << (char)e.qchr;
	return os;
}
void Edit::print(ostream &os, const EList<Edit> &edits, char delim)
{
	for (size_t i = 0; i < edits.size(); i++)
	{
		os << edits[i];
		if (i < edits.size() - 1)
			os << delim;
	}
}
void Edit::invertPoss(
	EList<Edit> &edits,
	size_t sz,
	size_t ei,
	size_t en,
	bool sort)
{
	size_t ii = 0;
	for (size_t i = ei; i < ei + en / 2; i++)
	{
		Edit tmp = edits[i];
		edits[i] = edits[ei + en - ii - 1];
		edits[ei + en - ii - 1] = tmp;
		ii++;
	}
	for (size_t i = ei; i < ei + en; i++)
	{
		assert(edits[i].pos < sz ||
			   (edits[i].isReadGap() && edits[i].pos == sz));
		edits[i].pos =
			(uint32_t)(sz - edits[i].pos - (edits[i].isReadGap() ? 0 : 1));
		if (edits[i].isReadGap())
		{
			int64_t pos2diff = (int64_t)(uint64_t)edits[i].pos2 - (int64_t)((uint64_t)std::numeric_limits<uint32_t>::max() >> 1);
			int64_t pos2new = (int64_t)(uint64_t)edits[i].pos2 - 2 * pos2diff;
			assert(pos2diff == 0 || (uint32_t)pos2new != (std::numeric_limits<uint32_t>::max() >> 1));
			edits[i].pos2 = (uint32_t)pos2new;
		}
	}
	if (sort)
	{
		edits.sortPortion(ei, en);
#ifndef NDEBUG
		for (size_t i = ei + 1; i < ei + en; i++)
		{
			assert_geq(edits[i].pos, edits[i - 1].pos);
		}
#endif
	}
}
void Edit::printQAlign(
	std::ostream &os,
	const BTDnaString &read,
	const EList<Edit> &edits)
{
	printQAlign(os, "", read, edits);
}
void Edit::printQAlignNoCheck(
	std::ostream &os,
	const BTDnaString &read,
	const EList<Edit> &edits)
{
	printQAlignNoCheck(os, "", read, edits);
}
void Edit::printQAlign(
	std::ostream &os,
	const char *prefix,
	const BTDnaString &read,
	const EList<Edit> &edits)
{
	size_t eidx = 0;
	os << prefix;
	for (size_t i = 0; i < read.length(); i++)
	{
		bool del = false, mm = false;
		while (eidx < edits.size() && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				os << '-';
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				os << read.toChar(i);
			}
			else
			{
				mm = true;
				assert(edits[eidx].isMismatch());
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				os << (char)edits[eidx].qchr;
			}
			eidx++;
		}
		if (!del && !mm)
			os << read.toChar(i);
	}
	os << endl;
	os << prefix;
	eidx = 0;
	for (size_t i = 0; i < read.length(); i++)
	{
		bool del = false, mm = false;
		while (eidx < edits.size() && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				os << ' ';
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				os << ' ';
			}
			else
			{
				mm = true;
				assert(edits[eidx].isMismatch());
				os << ' ';
			}
			eidx++;
		}
		if (!del && !mm)
			os << '|';
	}
	os << endl;
	os << prefix;
	eidx = 0;
	for (size_t i = 0; i < read.length(); i++)
	{
		bool del = false, mm = false;
		while (eidx < edits.size() && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				os << (char)edits[eidx].chr;
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				os << '-';
			}
			else
			{
				mm = true;
				assert(edits[eidx].isMismatch());
				os << (char)edits[eidx].chr;
			}
			eidx++;
		}
		if (!del && !mm)
			os << read.toChar(i);
	}
	os << endl;
}
void Edit::printQAlignNoCheck(
	std::ostream &os,
	const char *prefix,
	const BTDnaString &read,
	const EList<Edit> &edits)
{
	size_t eidx = 0;
	os << prefix;
	for (size_t i = 0; i < read.length(); i++)
	{
		bool del = false, mm = false;
		while (eidx < edits.size() && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				os << '-';
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				os << read.toChar(i);
			}
			else
			{
				mm = true;
				os << (char)edits[eidx].qchr;
			}
			eidx++;
		}
		if (!del && !mm)
			os << read.toChar(i);
	}
	os << endl;
	os << prefix;
	eidx = 0;
	for (size_t i = 0; i < read.length(); i++)
	{
		bool del = false, mm = false;
		while (eidx < edits.size() && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				os << ' ';
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				os << ' ';
			}
			else
			{
				mm = true;
				os << ' ';
			}
			eidx++;
		}
		if (!del && !mm)
			os << '|';
	}
	os << endl;
	os << prefix;
	eidx = 0;
	for (size_t i = 0; i < read.length(); i++)
	{
		bool del = false, mm = false;
		while (eidx < edits.size() && edits[eidx].pos == i)
		{
			if (edits[eidx].isReadGap())
			{
				os << (char)edits[eidx].chr;
			}
			else if (edits[eidx].isRefGap())
			{
				del = true;
				os << '-';
			}
			else
			{
				mm = true;
				os << (char)edits[eidx].chr;
			}
			eidx++;
		}
		if (!del && !mm)
			os << read.toChar(i);
	}
	os << endl;
}
void Edit::sort(EList<Edit> &edits)
{
	edits.sort();
}
void Edit::toRef(
	const BTDnaString &read,
	const EList<Edit> &edits,
	BTDnaString &ref,
	bool fw,
	size_t trim5,
	size_t trim3)
{
	size_t eidx = 0;
	const size_t rdlen = read.length();
	size_t trimBeg = fw ? trim5 : trim3;
	size_t trimEnd = fw ? trim3 : trim5;
	assert(Edit::repOk(edits, read, fw, trim5, trim3));
	if (!fw)
	{
		invertPoss(const_cast<EList<Edit> &>(edits), read.length() - trimBeg - trimEnd, false);
	}
	for (size_t i = 0; i < rdlen; i++)
	{
		ASSERT_ONLY(int c = read[i]);
		assert_range(0, 4, c);
		bool del = false, mm = false;
		bool append = i >= trimBeg && rdlen - i - 1 >= trimEnd;
		bool appendIns = i >= trimBeg && rdlen - i >= trimEnd;
		while (eidx < edits.size() && edits[eidx].pos + trimBeg == i)
		{
			if (edits[eidx].isReadGap())
			{
				if (appendIns)
				{
					ref.appendChar((char)edits[eidx].chr);
				}
			}
			else if (edits[eidx].isRefGap())
			{
				assert_eq("ACGTN"[c], edits[eidx].qchr);
				del = true;
			}
			else
			{
				mm = true;
				assert(edits[eidx].isMismatch());
				assert(edits[eidx].qchr != edits[eidx].chr || edits[eidx].qchr == 'N');
				assert_eq("ACGTN"[c], edits[eidx].qchr);
				if (append)
				{
					ref.appendChar((char)edits[eidx].chr);
				}
			}
			eidx++;
		}
		if (!del && !mm)
		{
			if (append)
			{
				ref.append(read[i]);
			}
		}
	}
	if (trimEnd == 0)
	{
		while (eidx < edits.size())
		{
			assert_gt(rdlen, edits[eidx].pos);
			if (edits[eidx].isReadGap())
			{
				ref.appendChar((char)edits[eidx].chr);
			}
			eidx++;
		}
	}
	if (!fw)
	{
		invertPoss(const_cast<EList<Edit> &>(edits), read.length() - trimBeg - trimEnd, false);
	}
}
#ifndef NDEBUG
bool Edit::repOk() const
{
	assert(inited());
	assert(qchr != chr || qchr == 'N');
	assert(isRefGap() || chr != '-');
	assert(isReadGap() || qchr != '-');
	assert(!isMismatch() || (qchr != '-' && chr != '-'));
	return true;
}
bool Edit::repOk(
	const EList<Edit> &edits,
	const BTDnaString &s,
	bool fw,
	size_t trimBeg,
	size_t trimEnd)
{
	if (!fw)
	{
		invertPoss(const_cast<EList<Edit> &>(edits), s.length() - trimBeg - trimEnd, false);
		swap(trimBeg, trimEnd);
	}
	for (size_t i = 0; i < edits.size(); i++)
	{
		const Edit &e = edits[i];
		size_t pos = e.pos;
		if (i > 0)
		{
			assert_geq(pos, edits[i - 1].pos);
		}
		bool del = false, mm = false;
		while (i < edits.size() && edits[i].pos == pos)
		{
			const Edit &ee = edits[i];
			assert_lt(ee.pos, s.length());
			if (ee.qchr != '-')
			{
				assert(ee.isRefGap() || ee.isMismatch());
				assert_eq((int)ee.qchr, s.toChar(ee.pos + trimBeg));
			}
			if (ee.isMismatch())
			{
				assert(!mm);
				mm = true;
				assert(!del);
			}
			else if (ee.isReadGap())
			{
				assert(!mm);
			}
			else if (ee.isRefGap())
			{
				assert(!mm);
				assert(!del);
				del = true;
			}
			i++;
		}
	}
	if (!fw)
	{
		invertPoss(const_cast<EList<Edit> &>(edits), s.length() - trimBeg - trimEnd, false);
	}
	return true;
}
#endif
void Edit::merge(EList<Edit> &dst, const EList<Edit> &src)
{
	size_t di = 0, si = 0;
	while (di < dst.size())
	{
		if (src[si].pos < dst[di].pos)
		{
			dst.insert(src[si], di);
			si++;
			di++;
		}
		else if (src[si].pos == dst[di].pos)
		{
			assert(src[si].isReadGap() != dst[di].isReadGap());
			if (src[si].isReadGap())
			{
				dst.insert(src[si], di);
				si++;
				di++;
			}
			else if (dst[di].isReadGap())
			{
				di++;
			}
		}
	}
	while (si < src.size())
		dst.push_back(src[si++]);
}
void Edit::clipLo(EList<Edit> &ed, size_t len, size_t amt)
{
	size_t nrm = 0;
	for (size_t i = 0; i < ed.size(); i++)
	{
		assert_lt(ed[i].pos, len);
		if (ed[i].pos < amt)
		{
			nrm++;
		}
		else
		{
			ed[i].pos -= (uint32_t)amt;
		}
	}
	ed.erase(0, nrm);
}
void Edit::clipHi(EList<Edit> &ed, size_t len, size_t amt)
{
	assert_leq(amt, len);
	size_t max = len - amt;
	size_t nrm = 0;
	for (size_t i = 0; i < ed.size(); i++)
	{
		size_t ii = ed.size() - i - 1;
		assert_lt(ed[ii].pos, len);
		if (ed[ii].pos > max)
		{
			nrm++;
		}
		else if (ed[ii].pos == max && !ed[ii].isReadGap())
		{
			nrm++;
		}
		else
		{
			break;
		}
	}
	ed.resize(ed.size() - nrm);
}
