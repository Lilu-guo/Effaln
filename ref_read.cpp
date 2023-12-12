#include "ref_read.h"
RefRecord fastaRefReadSize(
	FileBuf &in,
	const RefReadInParams &rparms,
	bool first,
	BitpairOutFileBuf *bpout)
{
	int c;
	static int lastc = '>';
	TIndexOffU len = 0;
	size_t off = 0;
	if (first)
	{
		assert(!in.eof());
		lastc = '>';
		c = in.getPastWhitespace();
		if (in.eof())
		{
			cerr << "Warning: Empty input file" << endl;
			lastc = -1;
			return RefRecord(0, 0, true);
		}
		assert(c == '>');
	}
	first = true;
	if (lastc == '>')
	{
		do
		{
			if ((c = in.getPastNewline()) == -1)
			{
				cerr << "Warning: Encountered empty reference sequence" << endl;
				lastc = -1;
				return RefRecord(0, 0, true);
			}
			if (c == '>')
			{
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
		} while (c == '>');
	}
	else
	{
		first = false;
		off = 1;
		if ((c = in.get()) == -1)
		{
			lastc = -1;
			return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
		}
	}
	int lc = -1;
	while (true)
	{
		int cat = asc2dnacat[c];
		if (rparms.nsToAs && cat >= 2)
			c = 'A';
		if (cat == 1)
		{
			break;
		}
		else if (cat >= 2)
		{
			if (lc != -1 && off == 0)
				off++;
			lc = -1;
			off++;
		}
		else if (c == '>')
		{
			if (off > 0 && lastc == '>')
			{
				cerr << "Warning: Encountered reference sequence with only gaps" << endl;
			}
			else if (lastc == '>')
			{
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
			lastc = '>';
			return RefRecord((TIndexOffU)off, 0, first);
		}
		c = in.get();
		if (c == -1)
		{
			if (off > 0 && lastc == '>')
			{
				cerr << "Warning: Encountered reference sequence with only gaps" << endl;
			}
			else if (lastc == '>')
			{
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
			lastc = -1;
			return RefRecord((TIndexOffU)off, 0, first);
		}
	}
	assert_eq(1, asc2dnacat[c]);
	while (c != -1 && c != '>')
	{
		if (rparms.nsToAs && asc2dnacat[c] >= 2)
			c = 'A';
		uint8_t cat = asc2dnacat[c];
		int cc = toupper(c);
		if (rparms.bisulfite && cc == 'C')
			c = cc = 'T';
		if (cat == 1)
		{
			assert(cc == 'A' || cc == 'C' || cc == 'G' || cc == 'T');
			if ((TIndexOffU)(len + 1) < len)
			{
				throw RefTooLongException();
			}
			len++;
			if (bpout != NULL)
			{
				bpout->write(asc2dna[c]);
			}
			lc = asc2dna[(int)c];
		}
		else if (cat >= 2)
		{
			lastc = c;
			assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
			return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
		}
		else
		{
#ifndef NDEBUG
			if (!isspace(c))
			{
				cerr << "Unexpected character in sequence: ";
				if (isprint(c))
				{
					cerr << ((char)c) << endl;
				}
				else
				{
					cerr << "(" << c << ")" << endl;
				}
			}
#endif
		}
		c = in.get();
	}
	lastc = c;
	return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
}
#if 0
static void
printRecords(ostream& os, const EList<RefRecord>& l) {
	for(size_t i = 0; i < l.size(); i++) {
		os << l[i].first << ", " << l[i].off << ", " << l[i].len << endl;
	}
}
#endif
void reverseRefRecords(
	const EList<RefRecord> &src,
	EList<RefRecord> &dst,
	bool recursive,
	bool verbose)
{
	dst.clear();
	{
		EList<RefRecord> cur;
		for (int i = (int)src.size() - 1; i >= 0; i--)
		{
			bool first = (i == (int)src.size() - 1 || src[i + 1].first);
			if (src[i].len || (first && src[i].off == 0))
			{
				cur.push_back(RefRecord(0, src[i].len, first));
				first = false;
			}
			if (src[i].off)
				cur.push_back(RefRecord(src[i].off, 0, first));
		}
		for (int i = 0; i < (int)cur.size(); i++)
		{
			assert(cur[i].off == 0 || cur[i].len == 0);
			if (i < (int)cur.size() - 1 && cur[i].off != 0 && !cur[i + 1].first)
			{
				dst.push_back(RefRecord(cur[i].off, cur[i + 1].len, cur[i].first));
				i++;
			}
			else
			{
				dst.push_back(cur[i]);
			}
		}
	}
#ifndef NDEBUG
	size_t srcnfirst = 0, dstnfirst = 0;
	for (size_t i = 0; i < src.size(); i++)
	{
		if (src[i].first)
		{
			srcnfirst++;
		}
	}
	for (size_t i = 0; i < dst.size(); i++)
	{
		if (dst[i].first)
		{
			dstnfirst++;
		}
	}
	assert_eq(srcnfirst, dstnfirst);
	if (!recursive)
	{
		EList<RefRecord> tmp;
		reverseRefRecords(dst, tmp, true);
		assert_eq(tmp.size(), src.size());
		for (size_t i = 0; i < src.size(); i++)
		{
			assert_eq(src[i].len, tmp[i].len);
			assert_eq(src[i].off, tmp[i].off);
			assert_eq(src[i].first, tmp[i].first);
		}
	}
#endif
}
std::pair<size_t, size_t>
fastaRefReadSizes(
	EList<FileBuf *> &in,
	EList<RefRecord> &recs,
	const RefReadInParams &rparms,
	BitpairOutFileBuf *bpout,
	TIndexOff &numSeqs)
{
	TIndexOffU unambigTot = 0;
	size_t bothTot = 0;
	assert_gt(in.size(), 0);
	for (size_t i = 0; i < in.size(); i++)
	{
		bool first = true;
		assert(!in[i]->eof());
		while (!in[i]->eof())
		{
			RefRecord rec;
			try
			{
				rec = fastaRefReadSize(*in[i], rparms, first, bpout);
				if ((unambigTot + rec.len) < unambigTot)
				{
					throw RefTooLongException();
				}
			}
			catch (RefTooLongException &e)
			{
				cerr << e.what() << endl;
				throw 1;
			}
			first = false;
			if (rec.len == 0 && rec.first)
			{
				continue;
			}
			else if (rec.first)
			{
				numSeqs++;
			}
			unambigTot += rec.len;
			bothTot += rec.len;
			bothTot += rec.off;
			if (rec.len == 0 && rec.off == 0 && !rec.first)
				continue;
			recs.push_back(rec);
		}
		in[i]->reset();
		assert(!in[i]->eof());
#ifndef NDEBUG
		int c = in[i]->get();
		assert_eq('>', c);
		in[i]->reset();
		assert(!in[i]->eof());
#endif
	}
	assert_geq(bothTot, 0);
	assert_geq(unambigTot, 0);
	return make_pair(
		unambigTot,
		bothTot);
}
