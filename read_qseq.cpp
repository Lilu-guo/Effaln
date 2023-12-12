#include "pat.h"
static int parseName(
	Read::TBuf &buf,
	size_t &cur, Read &r, int upto)
{
	const size_t buflen = buf.length();
	int c;
	while (cur < buflen)
	{
		c = buf[cur++];
		assert(c != '\r' && c != '\n');
		if (c == upto)
		{
			break;
		}
		r.name.append(c);
	}
	if (cur >= buflen)
	{
		return -1;
	}
	return (int)r.name.length();
}
pair<bool, int> QseqPatternSource::nextBatchFromFile(
	PerThreadReadBuf &pt,
	bool batch_a, unsigned readi)
{
	int c = getc_wrapper();
	while (c >= 0 && (c == '\n' || c == '\r'))
	{
		c = getc_wrapper();
	}
	EList<Read> &readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	for (; readi < pt.max_buf_ && c >= 0; readi++)
	{
		readbuf[readi].readOrigBuf.clear();
		while (c >= 0 && c != '\n' && c != '\r')
		{
			readbuf[readi].readOrigBuf.append(c);
			c = getc_wrapper();
		}
		while (c >= 0 && (c == '\n' || c == '\r'))
		{
			c = getc_wrapper();
		}
	}
	if (c != EOF)
	{
		ungetc_wrapper(c);
	}
	return make_pair(c < 0, readi);
}
bool QseqPatternSource::parse(Read &r, Read &rb, TReadId rdid) const
{
	assert(r.empty());
	assert(!r.readOrigBuf.empty());
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = r.readOrigBuf.length();
	assert(r.name.empty());
	if (parseName(r.readOrigBuf, cur, r, '\t') == -1)
	{
		return false;
	}
	r.name.append('_');
	if (parseName(r.readOrigBuf, cur, r, '\t') == -1)
	{
		return false;
	}
	r.name.append('_');
	if (parseName(r.readOrigBuf, cur, r, '\t') == -1)
	{
		return false;
	}
	r.name.append('_');
	if (parseName(r.readOrigBuf, cur, r, '\t') == -1)
	{
		return false;
	}
	r.name.append('_');
	if (parseName(r.readOrigBuf, cur, r, '\t') == -1)
	{
		return false;
	}
	r.name.append('_');
	if (parseName(r.readOrigBuf, cur, r, '\t') == -1)
	{
		return false;
	}
	r.name.append('_');
	if (parseName(r.readOrigBuf, cur, r, '\t') == -1)
	{
		return false;
	}
	r.name.append('/');
	if (parseName(r.readOrigBuf, cur, r, '\t') == -1)
	{
		return false;
	}
	if (cur >= buflen)
	{
		return false;
	}
	c = r.readOrigBuf[cur++];
	assert(c != '\r' && c != '\n');
	if (c == '\t')
	{
		c = r.readOrigBuf[cur++];
		assert(c != '\r' && c != '\n');
		assert_eq('\t', c);
		cerr << "Warning: skipping empty QSEQ read with name '" << r.name << "'" << endl;
	}
	else
	{
		int nchar = 0;
		while (c != '\t')
		{
			if (c == '.')
			{
				c = 'N';
			}
			if (isalpha(c))
			{
				assert_in(toupper(c), "ACGTN");
				if (++nchar > pp_.trim5)
				{
					assert_neq(0, asc2dnacat[c]);
					r.patFw.append(asc2dna[c]);
				}
			}
			if (cur >= buflen)
			{
				break;
			}
			c = r.readOrigBuf[cur++];
		}
		if (cur >= buflen)
		{
			return false;
		}
		r.trimmed5 = (int)(nchar - r.patFw.length());
		r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
		assert(r.qual.empty());
		int nqual = 0;
		if (pp_.intQuals)
		{
			int cur_int = 0;
			while (c != '\t')
			{
				cur_int *= 10;
				cur_int += (int)(c - '0');
				c = r.readOrigBuf[cur++];
				assert(c != '\r' && c != '\n');
				if (c == ' ' || c == '\t')
				{
					char cadd = intToPhred33(cur_int, pp_.solexa64);
					cur_int = 0;
					assert_geq(cadd, 33);
					if (++nqual > pp_.trim5)
					{
						r.qual.append(cadd);
					}
				}
			}
		}
		else
		{
			while (cur < buflen)
			{
				c = r.readOrigBuf[cur++];
				assert(c != '\r' && c != '\n');
				if (c == ' ')
				{
					wrongQualityFormat(r.name);
					return false;
				}
				else if (c == '\t')
				{
					break;
				}
				c = charToPhred33(c, pp_.solexa64, pp_.phred64);
				if (++nqual > r.trimmed5)
				{
					r.qual.append(c);
				}
			}
			r.qual.trimEnd(r.trimmed3);
			if (r.qual.length() < r.patFw.length())
			{
				tooFewQualities(r.name);
				return false;
			}
			else if (r.qual.length() > r.patFw.length())
			{
				tooManyQualities(r.name);
				return false;
			}
		}
	}
	assert_eq('\t', c);
	if (cur >= buflen)
	{
		return false;
	}
	int filt = r.readOrigBuf[cur++];
	r.filter = filt;
	if (filt != '0' && filt != '1')
	{
		cerr << "Error: Bad value '" << filt
			 << "' for qseq filter flag" << endl;
		throw 1;
	}
	assert_eq(cur, buflen);
	r.parsed = true;
	if (!rb.parsed && !rb.readOrigBuf.empty())
	{
		return parse(rb, r, rdid);
	}
	return true;
}
