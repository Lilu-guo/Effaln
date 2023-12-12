#ifndef REF_READ_H_
#define REF_READ_H_
#include <iostream>
#include <cassert>
#include <string>
#include <ctype.h>
#include <fstream>
#include <stdexcept>
#include "alphabet.h"
#include "assert_helpers.h"
#include "filebuf.h"
#include "word_io.h"
#include "ds.h"
#include "endian_swap.h"
using namespace std;
class RefTooLongException : public exception
{
public:
	RefTooLongException()
	{
#ifdef _64BIT_INDEX
		msg = "Error: Reference sequence has more than 2^64-1 characters!  ";
#else
		msg = "Error: Reference sequence has more than 2^32-1 characters!  ";
#endif
	}
	~RefTooLongException() throw() {}
	const char *what() const throw()
	{
		return msg.c_str();
	}
protected:
	string msg;
};
struct RefRecord
{
	RefRecord() : off(), len(), first() {}
	RefRecord(TIndexOffU _off, TIndexOffU _len, bool _first) : off(_off), len(_len), first(_first)
	{
	}
	RefRecord(FILE *in, bool swap)
	{
		assert(in != NULL);
		if (!fread(&off, OFF_SIZE, 1, in))
		{
			cerr << "Error reading RefRecord offset from FILE" << endl;
			throw 1;
		}
		if (swap)
			off = endianSwapU(off);
		if (!fread(&len, OFF_SIZE, 1, in))
		{
			cerr << "Error reading RefRecord offset from FILE" << endl;
			throw 1;
		}
		if (swap)
			len = endianSwapU(len);
		first = fgetc(in) ? true : false;
	}
	void write(std::ostream &out, bool be)
	{
		writeU<TIndexOffU>(out, off, be);
		writeU<TIndexOffU>(out, len, be);
		out.put(first ? 1 : 0);
	}
	TIndexOffU off;
	TIndexOffU len;
	bool first;
};
enum
{
	REF_READ_FORWARD = 0,
	REF_READ_REVERSE,
	REF_READ_REVERSE_EACH
};
struct RefReadInParams
{
	RefReadInParams(bool col, int r, bool nsToA, bool bisulf) : color(col), reverse(r), nsToAs(nsToA), bisulfite(bisulf) {}
	bool color;
	int reverse;
	bool nsToAs;
	bool bisulfite;
};
extern RefRecord
fastaRefReadSize(
	FileBuf &in,
	const RefReadInParams &rparms,
	bool first,
	BitpairOutFileBuf *bpout = NULL);
extern std::pair<size_t, size_t>
fastaRefReadSizes(
	EList<FileBuf *> &in,
	EList<RefRecord> &recs,
	const RefReadInParams &rparms,
	BitpairOutFileBuf *bpout,
	TIndexOff &numSeqs);
extern void
reverseRefRecords(
	const EList<RefRecord> &src,
	EList<RefRecord> &dst,
	bool recursive = false,
	bool verbose = false);
template <typename TStr>
static RefRecord fastaRefReadAppend(
	FileBuf &in,
	bool first,
	TStr &dst, TIndexOffU &dstoff, RefReadInParams &rparms, string *name = NULL)
{
	int c;
	static int lastc = '>';
	if (first)
	{
		c = in.getPastWhitespace();
		if (c != '>')
		{
			cerr << "Reference file does not seem to be a FASTA file" << endl;
			throw 1;
		}
		lastc = c;
	}
	assert_neq(-1, lastc);
	size_t len = 0;
	size_t off = 0;
	first = true;
	size_t ilen = dstoff;
	int lc = -1;
	c = lastc;
	if (c == '>' || c == '#')
	{
		do
		{
			while (c == '#')
			{
				if ((c = in.getPastNewline()) == -1)
				{
					lastc = -1;
					goto bail;
				}
			}
			assert_eq('>', c);
			while (true)
			{
				c = in.get();
				if (c == -1)
				{
					lastc = -1;
					goto bail;
				}
				if (c == '\n' || c == '\r')
				{
					while (c == '\r' || c == '\n')
						c = in.get();
					if (c == -1)
					{
						lastc = -1;
						goto bail;
					}
					break;
				}
				if (name)
					name->push_back(c);
			}
			if (c == '>')
			{
				if (name)
					name->clear();
			}
		} while (c == '>' || c == '#');
	}
	else
	{
		ASSERT_ONLY(int cc = toupper(c));
		assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
		first = false;
	}
	while (true)
	{
		int cat = asc2dnacat[c];
		if (rparms.nsToAs && cat >= 2)
		{
			c = 'A';
		}
		int cc = toupper(c);
		if (rparms.bisulfite && cc == 'C')
			c = cc = 'T';
		if (cat == 1)
		{
			if (rparms.color)
			{
				if (lc != -1)
				{
					break;
				}
				lc = asc2dna[(int)c];
				if (off > 0)
					off++;
			}
			else
			{
				break;
			}
		}
		else if (cat >= 2)
		{
			if (lc != -1 && off == 0)
			{
				off++;
			}
			lc = -1;
			off++;
		}
		else if (c == '>')
		{
			lastc = '>';
			goto bail;
		}
		c = in.get();
		if (c == -1)
		{
			lastc = -1;
			goto bail;
		}
	}
	if (first && rparms.color && off > 0)
	{
		off--;
	}
	assert(!rparms.color || lc != -1);
	assert_eq(1, asc2dnacat[c]);
	while (true)
	{
		int cat = asc2dnacat[c];
		assert_neq(2, cat);
		if (cat == 1)
		{
			if (!rparms.color || lc != -1)
				len++;
			if (rparms.color)
			{
				dst.set((char)dinuc2color[asc2dna[(int)c]][lc], dstoff++);
			}
			else if (!rparms.color)
			{
				dst.set(asc2dna[c], dstoff++);
			}
			assert_lt((int)dst[dstoff - 1], 4);
			lc = asc2dna[(int)c];
		}
		c = in.get();
		if (rparms.nsToAs && asc2dnacat[c] >= 2)
			c = 'A';
		if (c == -1 || c == '>' || c == '#' || asc2dnacat[c] >= 2)
		{
			lastc = c;
			break;
		}
		if (rparms.bisulfite && toupper(c) == 'C')
			c = 'T';
	}
bail:
	if (rparms.reverse == REF_READ_REVERSE_EACH)
	{
		size_t nlen = dstoff;
		dst.reverseWindow(ilen, nlen);
	}
	return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
}
#endif
