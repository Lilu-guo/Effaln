#ifndef FILEBUF_H_
#define FILEBUF_H_
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdexcept>
#include "assert_helpers.h"
#include <errno.h>
#include <stdlib.h>
#include <zlib.h>
static inline bool isnewline(int c)
{
	return c == '\r' || c == '\n';
}
static inline bool isspace_notnl(int c)
{
	return isspace(c) && !isnewline(c);
}
class FileBuf
{
public:
	FileBuf()
	{
		init();
	}
	FileBuf(FILE *in)
	{
		init();
		_in = in;
		assert(_in != NULL);
	}
	FileBuf(gzFile in)
	{
		init();
		_zIn = in;
		assert(_zIn != NULL);
	}
	FileBuf(std::ifstream *inf)
	{
		init();
		_inf = inf;
		assert(_inf != NULL);
	}
	FileBuf(std::istream *ins)
	{
		init();
		_ins = ins;
		assert(_ins != NULL);
	}
	~FileBuf()
	{
		close();
	}
	bool isOpen()
	{
		return _in != NULL || _inf != NULL || _ins != NULL;
	}
	void close()
	{
		if (_in != NULL && _in != stdin)
		{
			fclose(_in);
		}
		else if (_inf != NULL)
		{
			_inf->close();
		}
		else if (_zIn != NULL)
		{
			gzclose(_zIn);
		}
		else
		{
		}
	}
	int get()
	{
		assert(_in != NULL || _zIn != NULL || _inf != NULL || _ins != NULL);
		int c = peek();
		if (c != -1)
		{
			_cur++;
			if (_lastn_cur < LASTN_BUF_SZ)
				_lastn_buf[_lastn_cur++] = c;
		}
		return c;
	}
	bool eof()
	{
		return (_cur == _buf_sz) && _done;
	}
	void newFile(FILE *in)
	{
		_in = in;
		_zIn = NULL;
		_inf = NULL;
		_ins = NULL;
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}
	void newFile(gzFile in)
	{
		_in = NULL;
		_zIn = in;
		_inf = NULL;
		_ins = NULL;
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}
	void newFile(std::ifstream *__inf)
	{
		_in = NULL;
		_zIn = NULL;
		_inf = __inf;
		_ins = NULL;
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}
	void newFile(std::istream *__ins)
	{
		_in = NULL;
		_zIn = NULL;
		_inf = NULL;
		_ins = __ins;
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}
	void reset()
	{
		if (_inf != NULL)
		{
			_inf->clear();
			_inf->seekg(0, std::ios::beg);
		}
		else if (_ins != NULL)
		{
			_ins->clear();
			_ins->seekg(0, std::ios::beg);
		}
		else if (_zIn != NULL)
		{
			gzrewind(_zIn);
		}
		else
		{
			rewind(_in);
		}
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}
	int peek()
	{
		assert(_in != NULL || _zIn != NULL || _inf != NULL || _ins != NULL);
		assert_leq(_cur, _buf_sz);
		if (_cur == _buf_sz)
		{
			if (_done)
			{
				return -1;
			}
			else
			{
				if (_inf != NULL)
				{
					_inf->read((char *)_buf, BUF_SZ);
					_buf_sz = _inf->gcount();
				}
				else if (_zIn != NULL)
				{
					_buf_sz = gzread(_zIn, (void *)_buf, BUF_SZ);
				}
				else if (_ins != NULL)
				{
					_ins->read((char *)_buf, BUF_SZ);
					_buf_sz = _ins->gcount();
				}
				else
				{
					assert(_in != NULL);
					_buf_sz = fread(_buf, 1, BUF_SZ, _in);
				}
				_cur = 0;
				if (_buf_sz == 0)
				{
					_done = true;
					return -1;
				}
				else if (_buf_sz < BUF_SZ)
				{
					_done = true;
				}
			}
		}
		return (int)_buf[_cur];
	}
	size_t gets(char *buf, size_t len)
	{
		size_t stored = 0;
		while (true)
		{
			int c = get();
			if (c == -1)
			{
				buf[stored] = '\0';
				return stored;
			}
			if (stored == len - 1 || isnewline(c))
			{
				buf[stored] = '\0';
				int pc = peek();
				while (isnewline(pc))
				{
					get();
					pc = peek();
				}
				return stored;
			}
			buf[stored++] = (char)c;
		}
	}
	size_t get(char *buf, size_t len)
	{
		size_t stored = 0;
		for (size_t i = 0; i < len; i++)
		{
			int c = get();
			if (c == -1)
				return i;
			buf[stored++] = (char)c;
		}
		return len;
	}
	static const size_t LASTN_BUF_SZ = 8 * 1024;
	int getPastWhitespace()
	{
		int c;
		while (isspace(c = get()) && c != -1)
			;
		return c;
	}
	int getPastNewline()
	{
		int c = get();
		while (!isnewline(c) && c != -1)
			c = get();
		while (isnewline(c))
			c = get();
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}
	int peekPastNewline()
	{
		int c = peek();
		while (!isnewline(c) && c != -1)
			c = get();
		while (isnewline(c))
			c = get();
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}
	int peekUptoNewline()
	{
		int c = peek();
		while (!isnewline(c) && c != -1)
		{
			get();
			c = peek();
		}
		while (isnewline(c))
		{
			get();
			c = peek();
		}
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}
	template <typename TNameStr, typename TSeqStr>
	void parseFastaRecord(
		TNameStr &name,
		TSeqStr &seq,
		bool gotCaret = false)
	{
		int c;
		if (!gotCaret)
		{
			c = peek();
			while (isspace_notnl(c) || c == '>')
			{
				get();
				c = peek();
			}
		}
		else
		{
			c = peek();
			while (isspace_notnl(c))
			{
				get();
				c = peek();
			}
		}
		size_t namecur = 0, seqcur = 0;
		while (!isnewline(c) && c != -1)
		{
			name[namecur++] = c;
			get();
			c = peek();
		}
		while (true)
		{
			while (isspace(c))
			{
				get();
				c = peek();
			}
			if (c == '>' || c == -1)
				break;
			seq[seqcur++] = c;
			get();
			c = peek();
		}
	}
	void parseFastaRecordLength(
		size_t &nameLen,
		size_t &seqLen,
		bool gotCaret = false)
	{
		int c;
		nameLen = seqLen = 0;
		if (!gotCaret)
		{
			c = peek();
			while (isspace_notnl(c) || c == '>')
			{
				get();
				c = peek();
			}
		}
		else
		{
			c = peek();
			while (isspace_notnl(c))
			{
				get();
				c = peek();
			}
		}
		while (!isnewline(c) && c != -1)
		{
			nameLen++;
			get();
			c = peek();
		}
		while (true)
		{
			while (isspace(c))
			{
				get();
				c = peek();
			}
			if (c == '>' || c == -1)
				break;
			seqLen++;
			get();
			c = peek();
		}
	}
	void resetLastN()
	{
		_lastn_cur = 0;
	}
	size_t copyLastN(char *buf)
	{
		memcpy(buf, _lastn_buf, _lastn_cur);
		return _lastn_cur;
	}
	const char *lastN() const
	{
		return _lastn_buf;
	}
	size_t lastNLen() const
	{
		return _lastn_cur;
	}
private:
	void init()
	{
		_in = NULL;
		_zIn = NULL;
		_inf = NULL;
		_ins = NULL;
		_cur = _buf_sz = BUF_SZ;
		_done = false;
		_lastn_cur = 0;
	}
	static const size_t BUF_SZ = 256 * 1024;
	FILE *_in;
	gzFile _zIn;
	std::ifstream *_inf;
	std::istream *_ins;
	size_t _cur;
	size_t _buf_sz;
	bool _done;
	uint8_t _buf[BUF_SZ];
	size_t _lastn_cur;
	char _lastn_buf[LASTN_BUF_SZ];
};
class BitpairOutFileBuf
{
public:
	BitpairOutFileBuf(const char *in) : bpPtr_(0), cur_(0)
	{
		assert(in != NULL);
		out_ = fopen(in, "wb");
		if (out_ == NULL)
		{
			std::cerr << "Error: Could not open bitpair-output file " << in << std::endl;
			throw 1;
		}
		memset(buf_, 0, BUF_SZ);
	}
	void write(int bp)
	{
		assert_lt(bp, 4);
		assert_geq(bp, 0);
		buf_[cur_] |= (bp << bpPtr_);
		if (bpPtr_ == 6)
		{
			bpPtr_ = 0;
			cur_++;
			if (cur_ == BUF_SZ)
			{
				if (!fwrite((const void *)buf_, BUF_SZ, 1, out_))
				{
					std::cerr << "Error writing to the reference index file (.4.ebwt)" << std::endl;
					throw 1;
				}
				cur_ = 0;
			}
			buf_[cur_] = 0;
		}
		else
		{
			bpPtr_ += 2;
		}
	}
	void close()
	{
		if (cur_ > 0 || bpPtr_ > 0)
		{
			if (bpPtr_ == 0)
				cur_--;
			if (!fwrite((const void *)buf_, cur_ + 1, 1, out_))
			{
				std::cerr << "Error writing to the reference index file (.4.ebwt)" << std::endl;
				throw 1;
			}
		}
		fclose(out_);
	}
private:
	static const size_t BUF_SZ = 128 * 1024;
	FILE *out_;
	int bpPtr_;
	size_t cur_;
	char buf_[BUF_SZ];
};
class OutFileBuf
{
public:
	OutFileBuf(const std::string &out, bool binary = false) : name_(out.c_str()), cur_(0), closed_(false)
	{
		out_ = fopen(out.c_str(), binary ? "wb" : "w");
		if (out_ == NULL)
		{
			std::cerr << "Error: Could not open alignment output file " << out.c_str() << std::endl;
			throw 1;
		}
		if (setvbuf(out_, NULL, _IOFBF, 10 * 1024 * 1024))
			std::cerr << "Warning: Could not allocate the proper buffer size for output file stream. " << std::endl;
	}
	OutFileBuf(const char *out, bool binary = false) : name_(out), cur_(0), closed_(false)
	{
		assert(out != NULL);
		out_ = fopen(out, binary ? "wb" : "w");
		if (out_ == NULL)
		{
			std::cerr << "Error: Could not open alignment output file " << out << std::endl;
			throw 1;
		}
	}
	OutFileBuf() : name_("cout"), cur_(0), closed_(false)
	{
		out_ = stdout;
	}
	~OutFileBuf() { close(); }
	void setFile(const char *out, bool binary = false)
	{
		assert(out != NULL);
		out_ = fopen(out, binary ? "wb" : "w");
		if (out_ == NULL)
		{
			std::cerr << "Error: Could not open alignment output file " << out << std::endl;
			throw 1;
		}
		reset();
	}
	void write(char c)
	{
		assert(!closed_);
		if (cur_ == BUF_SZ)
			flush();
		buf_[cur_++] = c;
	}
	void writeString(const std::string &s)
	{
		assert(!closed_);
		size_t slen = s.length();
		if (cur_ + slen > BUF_SZ)
		{
			if (cur_ > 0)
				flush();
			if (slen >= BUF_SZ)
			{
				if (slen != fwrite(s.c_str(), 1, slen, out_))
				{
					std::cerr << "Error: outputting data" << std::endl;
					throw 1;
				}
			}
			else
			{
				memcpy(&buf_[cur_], s.data(), slen);
				assert_eq(0, cur_);
				cur_ = slen;
			}
		}
		else
		{
			memcpy(&buf_[cur_], s.data(), slen);
			cur_ += slen;
		}
		assert_leq(cur_, BUF_SZ);
	}
	template <typename T>
	void writeString(const T &s)
	{
		assert(!closed_);
		size_t slen = s.length();
		if (cur_ + slen > BUF_SZ)
		{
			if (cur_ > 0)
				flush();
			if (slen >= BUF_SZ)
			{
				if (slen != fwrite(s.toZBuf(), 1, slen, out_))
				{
					std::cerr << "Error outputting data" << std::endl;
					throw 1;
				}
			}
			else
			{
				memcpy(&buf_[cur_], s.toZBuf(), slen);
				assert_eq(0, cur_);
				cur_ = slen;
			}
		}
		else
		{
			memcpy(&buf_[cur_], s.toZBuf(), slen);
			cur_ += slen;
		}
		assert_leq(cur_, BUF_SZ);
	}
	void writeChars(const char *s, size_t len)
	{
		assert(!closed_);
		if (cur_ + len > BUF_SZ)
		{
			if (cur_ > 0)
				flush();
			if (len >= BUF_SZ)
			{
				if (fwrite(s, len, 1, out_) != 1)
				{
					std::cerr << "Error outputting data" << std::endl;
					throw 1;
				}
			}
			else
			{
				memcpy(&buf_[cur_], s, len);
				assert_eq(0, cur_);
				cur_ = len;
			}
		}
		else
		{
			memcpy(&buf_[cur_], s, len);
			cur_ += len;
		}
		assert_leq(cur_, BUF_SZ);
	}
	void writeChars(const char *s)
	{
		writeChars(s, strlen(s));
	}
	void close()
	{
		if (closed_)
			return;
		if (cur_ > 0)
			flush();
		closed_ = true;
		if (out_ != stdout)
		{
			fclose(out_);
		}
	}
	void reset()
	{
		cur_ = 0;
		closed_ = false;
	}
	void flush()
	{
		if (cur_ != fwrite((const void *)buf_, 1, cur_, out_))
		{
			if (errno == EPIPE)
			{
				exit(EXIT_SUCCESS);
			}
			std::cerr << "Error while flushing and closing output" << std::endl;
			throw 1;
		}
		cur_ = 0;
	}
	bool closed() const
	{
		return closed_;
	}
	const char *name()
	{
		return name_;
	}
private:
	static const size_t BUF_SZ = 16 * 1024;
	const char *name_;
	FILE *out_;
	size_t cur_;
	char buf_[BUF_SZ];
	bool closed_;
};
#endif
