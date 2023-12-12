#ifndef SSTRING_H_
#define SSTRING_H_
#include <string.h>
#include <iostream>
#include "assert_helpers.h"
#include "alphabet.h"
#include "random_source.h"
template <typename T>
class Class_sstr_len
{
public:
	static inline size_t sstr_len(const T &s)
	{
		return s.length();
	}
};
template <unsigned N>
class Class_sstr_len<const char[N]>
{
public:
	static inline size_t sstr_len(const char s[N])
	{
		return strlen(s);
	}
};
template <>
class Class_sstr_len<const char *>
{
public:
	static inline size_t sstr_len(const char *s)
	{
		return strlen(s);
	}
};
template <>
class Class_sstr_len<const unsigned char *>
{
public:
	static inline size_t sstr_len(const unsigned char *s)
	{
		return strlen((const char *)s);
	}
};
template <typename T1, typename T2>
static inline bool sstr_eq(const T1 &s1, const T2 &s2)
{
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	if (len1 != len2)
		return false;
	for (size_t i = 0; i < len1; i++)
	{
		if (s1[i] != s2[i])
			return false;
	}
	return true;
}
template <typename T1, typename T2>
static inline bool sstr_neq(const T1 &s1, const T2 &s2)
{
	return !sstr_eq(s1, s2);
}
template <typename T1, typename T2>
static inline bool sstr_suf_upto_eq(
	const T1 &s1, size_t suf1,
	const T2 &s2, size_t suf2,
	size_t upto,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	if (len1 > upto)
		len1 = upto;
	if (len2 > upto)
		len2 = upto;
	if (len1 != len2)
		return false;
	for (size_t i = 0; i < len1; i++)
	{
		if (s1[suf1 + i] != s2[suf2 + i])
		{
			return false;
		}
	}
	return true;
}
template <typename T1, typename T2>
static inline bool sstr_suf_upto_neq(
	const T1 &s1, size_t suf1,
	const T2 &s2, size_t suf2,
	size_t upto,
	bool endlt = true)
{
	return !sstr_suf_upto_eq(s1, suf1, s2, suf2, upto, endlt);
}
template <typename T1, typename T2>
static inline bool sstr_lt(const T1 &s1, const T2 &s2, bool endlt = true)
{
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[i] < s2[i])
		{
			return true;
		}
		else if (s1[i] > s2[i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return false;
	return (len1 < len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_suf_lt(
	const T1 &s1, size_t suf1,
	const T2 &s2, size_t suf2,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[suf1 + i] < s2[suf2 + i])
		{
			return true;
		}
		else if (s1[suf1 + i] > s2[suf2 + i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return false;
	return (len1 < len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_suf_lt(
	const T1 &s1, size_t suf1, size_t len1,
	const T2 &s2, size_t suf2, size_t len2,
	bool endlt = true)
{
	assert_leq(suf1, len1);
	assert_leq(suf2, len2);
	size_t left1 = len1 - suf1;
	size_t left2 = len2 - suf2;
	size_t minleft = (left1 < left2 ? left1 : left2);
	for (size_t i = 0; i < minleft; i++)
	{
		if (s1[suf1 + i] < s2[suf2 + i])
		{
			return true;
		}
		else if (s1[suf1 + i] > s2[suf2 + i])
		{
			return false;
		}
	}
	if (left1 == left2)
		return false;
	return (left1 < left2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_suf_upto_lt(
	const T1 &s1, size_t suf1,
	const T2 &s2, size_t suf2,
	size_t upto,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	if (len1 > upto)
		len1 = upto;
	if (len2 > upto)
		len2 = upto;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[suf1 + i] < s2[suf2 + i])
		{
			return true;
		}
		else if (s1[suf1 + i] > s2[suf2 + i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return false;
	return (len1 < len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_pre_lt(
	const T1 &s1, size_t pre1,
	const T2 &s2, size_t pre2,
	bool endlt = true)
{
	assert_leq(pre1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(pre2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = pre1;
	size_t len2 = pre2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[i] < s2[i])
		{
			return true;
		}
		else if (s1[i] > s2[i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return false;
	return (len1 < len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_leq(const T1 &s1, const T2 &s2, bool endlt = true)
{
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[i] < s2[i])
		{
			return true;
		}
		else if (s1[i] > s2[i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return true;
	return (len1 < len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_suf_leq(
	const T1 &s1, size_t suf1,
	const T2 &s2, size_t suf2,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[suf1 + i] < s2[suf2 + i])
		{
			return true;
		}
		else if (s1[suf1 + i] > s2[suf2 + i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return true;
	return (len1 < len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_pre_leq(
	const T1 &s1, size_t pre1,
	const T2 &s2, size_t pre2,
	bool endlt = true)
{
	assert_leq(pre1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(pre2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = pre1;
	size_t len2 = pre2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[i] < s2[i])
		{
			return true;
		}
		else if (s1[i] > s2[i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return true;
	return (len1 < len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_gt(const T1 &s1, const T2 &s2, bool endlt = true)
{
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[i] > s2[i])
		{
			return true;
		}
		else if (s1[i] < s2[i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return false;
	return (len1 > len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_suf_gt(
	const T1 &s1, size_t suf1,
	const T2 &s2, size_t suf2,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[suf1 + i] > s2[suf2 + i])
		{
			return true;
		}
		else if (s1[suf1 + i] < s2[suf2 + i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return false;
	return (len1 > len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_pre_gt(
	const T1 &s1, size_t pre1,
	const T2 &s2, size_t pre2,
	bool endlt = true)
{
	assert_leq(pre1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(pre2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = pre1;
	size_t len2 = pre2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[i] > s2[i])
		{
			return true;
		}
		else if (s1[i] < s2[i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return false;
	return (len1 > len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_geq(const T1 &s1, const T2 &s2, bool endlt = true)
{
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1);
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2);
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[i] > s2[i])
		{
			return true;
		}
		else if (s1[i] < s2[i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return true;
	return (len1 > len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_suf_geq(
	const T1 &s1, size_t suf1,
	const T2 &s2, size_t suf2,
	bool endlt = true)
{
	assert_leq(suf1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(suf2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = Class_sstr_len<T1>::sstr_len(s1) - suf1;
	size_t len2 = Class_sstr_len<T2>::sstr_len(s2) - suf2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[suf1 + i] > s2[suf2 + i])
		{
			return true;
		}
		else if (s1[suf1 + i] < s2[suf2 + i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return true;
	return (len1 > len2) == endlt;
}
template <typename T1, typename T2>
static inline bool sstr_pre_geq(
	const T1 &s1, size_t pre1,
	const T2 &s2, size_t pre2,
	bool endlt = true)
{
	assert_leq(pre1, Class_sstr_len<T1>::sstr_len(s1));
	assert_leq(pre2, Class_sstr_len<T2>::sstr_len(s2));
	size_t len1 = pre1;
	size_t len2 = pre2;
	size_t minlen = (len1 < len2 ? len1 : len2);
	for (size_t i = 0; i < minlen; i++)
	{
		if (s1[i] > s2[i])
		{
			return true;
		}
		else if (s1[i] < s2[i])
		{
			return false;
		}
	}
	if (len1 == len2)
		return true;
	return (len1 > len2) == endlt;
}
template <typename T>
static inline const char *sstr_to_cstr(const T &s)
{
	return s.toZBuf();
}
template <>
inline const char *sstr_to_cstr<std::basic_string<char>>(
	const std::basic_string<char> &s)
{
	return s.c_str();
}
template <typename T>
class SString
{
public:
	explicit SString() : cs_(NULL),
						 printcs_(NULL),
						 len_(0)
	{
	}
	explicit SString(size_t sz) : cs_(NULL),
								  printcs_(NULL),
								  len_(0)
	{
		resize(sz);
	}
	SString(const SString<T> &o) : cs_(NULL),
								   printcs_(NULL),
								   len_(0)
	{
		*this = o;
	}
	explicit SString(const std::basic_string<T> &str) : cs_(NULL),
														printcs_(NULL),
														len_(0)
	{
		install(str.c_str(), str.length());
	}
	explicit SString(const T *b, size_t sz) : cs_(NULL),
											  printcs_(NULL),
											  len_(0)
	{
		install(b, sz);
	}
	explicit SString(const T *b) : cs_(NULL),
								   printcs_(NULL),
								   len_(0)
	{
		install(b, strlen(b));
	}
	virtual ~SString()
	{
		if (cs_ != NULL)
		{
			delete[] cs_;
			cs_ = NULL;
		}
		if (printcs_ != NULL)
		{
			delete[] printcs_;
			printcs_ = NULL;
		}
		len_ = 0;
	}
	SString<T> &operator=(const SString<T> &o)
	{
		install(o.cs_, o.len_);
		return *this;
	}
	SString<T> &operator=(const std::basic_string<T> &o)
	{
		install(o);
		return *this;
	}
	void resize(size_t sz)
	{
		if (cs_ != NULL)
		{
			delete cs_;
			cs_ = NULL;
		}
		if (printcs_ != NULL)
		{
			delete printcs_;
			printcs_ = NULL;
		}
		if (sz != 0)
		{
			cs_ = new T[sz + 1];
		}
		len_ = sz;
	}
	T windowGet(
		size_t i,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_lt(i, len);
		assert_leq(len, len_ - depth);
		return fw ? cs_[depth + i] : cs_[depth + len - i - 1];
	}
	void windowGet(
		T &ret,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_leq(len, len_ - depth);
		ret.resize(len);
		for (size_t i = 0; i < len; i++)
		{
			ret.set(fw ? cs_[depth + i] : cs_[depth + len - i - 1], i);
		}
	}
	inline void set(int c, size_t idx)
	{
		assert_lt(idx, len_);
		cs_[idx] = c;
	}
	inline const T &operator[](size_t i) const
	{
		assert_lt(i, len_);
		return cs_[i];
	}
	inline T &operator[](size_t i)
	{
		assert_lt(i, len_);
		return cs_[i];
	}
	inline const T &get(size_t i) const
	{
		assert_lt(i, len_);
		return cs_[i];
	}
	virtual void install(const T *b, size_t sz)
	{
		if (sz == 0)
			return;
		resize(sz);
		memcpy(cs_, b, sz * sizeof(T));
	}
	virtual void install(const std::basic_string<T> &b)
	{
		size_t sz = b.length();
		if (sz == 0)
			return;
		resize(sz);
		memcpy(cs_, b.c_str(), sz * sizeof(T));
	}
	void install(const T *b)
	{
		install(b, strlen(b));
	}
	void installReverse(const char *b, size_t sz)
	{
		if (sz == 0)
			return;
		resize(sz);
		for (size_t i = 0; i < sz; i++)
		{
			cs_[i] = b[sz - i - 1];
		}
		len_ = sz;
	}
	void installReverse(const SString<T> &b)
	{
		installReverse(b.cs_, b.len_);
	}
	bool operator==(const SString<T> &o)
	{
		return sstr_eq(*this, o);
	}
	bool operator!=(const SString<T> &o)
	{
		return sstr_neq(*this, o);
	}
	bool operator<(const SString<T> &o)
	{
		return sstr_lt(*this, o);
	}
	bool operator>(const SString<T> &o)
	{
		return sstr_gt(*this, o);
	}
	bool operator<=(const SString<T> &o)
	{
		return sstr_leq(*this, o);
	}
	bool operator>=(const SString<T> &o)
	{
		return sstr_geq(*this, o);
	}
	void reverse()
	{
		for (size_t i = 0; i < (len_ >> 1); i++)
		{
			T tmp = get(i);
			set(get(len_ - i - 1), i);
			set(tmp, len_ - i - 1);
		}
	}
	void reverseWindow(size_t off, size_t len)
	{
		assert_leq(off, len_);
		assert_leq(off + len, len_);
		size_t mid = len >> 1;
		for (size_t i = 0; i < mid; i++)
		{
			T tmp = get(off + i);
			set(get(off + len - i - 1), off + i);
			set(tmp, off + len - i - 1);
		}
	}
	void fill(size_t len, const T &el)
	{
		assert_leq(len, len_);
		for (size_t i = 0; i < len; i++)
		{
			set(el, i);
		}
	}
	void fill(const T &el)
	{
		fill(len_, el);
	}
	inline size_t length() const { return len_; }
	void clear() { len_ = 0; }
	inline bool empty() const { return len_ == 0; }
	const char *toZBufXForm(const char *xform) const
	{
		ASSERT_ONLY(size_t xformElts = strlen(xform));
		if (printcs_ == NULL)
		{
			const_cast<char *&>(printcs_) = new char[len_ + 1];
		}
		char *printcs = const_cast<char *>(printcs_);
		assert(printcs != NULL);
		for (size_t i = 0; i < len_; i++)
		{
			assert_lt(cs_[i], (int)xformElts);
			printcs[i] = xform[cs_[i]];
		}
		printcs[len_] = 0;
		return printcs_;
	}
	virtual const T *toZBuf() const
	{
		const_cast<T *>(cs_)[len_] = 0;
		return cs_;
	}
	const T *buf() const { return cs_; }
	T *wbuf() { return cs_; }
protected:
	T *cs_;
	char *printcs_;
	size_t len_;
};
class S2bDnaString
{
public:
	explicit S2bDnaString() : cs_(NULL),
							  printcs_(NULL),
							  len_(0)
	{
	}
	explicit S2bDnaString(size_t sz) : cs_(NULL),
									   printcs_(NULL),
									   len_(0)
	{
		resize(sz);
	}
	S2bDnaString(const S2bDnaString &o) : cs_(NULL),
										  printcs_(NULL),
										  len_(0)
	{
		*this = o;
	}
	explicit S2bDnaString(
		const std::basic_string<char> &str,
		bool chars = false,
		bool colors = false) : cs_(NULL),
							   printcs_(NULL),
							   len_(0)
	{
		if (chars)
		{
			if (colors)
			{
				installColors(str.c_str(), str.length());
			}
			else
			{
				installChars(str.c_str(), str.length());
			}
		}
		else
		{
			install(str.c_str(), str.length());
		}
	}
	explicit S2bDnaString(
		const char *b,
		size_t sz,
		bool chars = false,
		bool colors = false) : cs_(NULL),
							   printcs_(NULL),
							   len_(0)
	{
		if (chars)
		{
			if (colors)
			{
				installColors(b, sz);
			}
			else
			{
				installChars(b, sz);
			}
		}
		else
		{
			install(b, sz);
		}
	}
	explicit S2bDnaString(
		const char *b,
		bool chars = false,
		bool colors = false) : cs_(NULL),
							   printcs_(NULL),
							   len_(0)
	{
		if (chars)
		{
			if (colors)
			{
				installColors(b, strlen(b));
			}
			else
			{
				installChars(b, strlen(b));
			}
		}
		else
		{
			install(b, strlen(b));
		}
	}
	virtual ~S2bDnaString()
	{
		if (cs_ != NULL)
		{
			delete[] cs_;
			cs_ = NULL;
		}
		if (printcs_ != NULL)
		{
			delete[] printcs_;
			printcs_ = NULL;
		}
		len_ = 0;
	}
	template <typename T>
	S2bDnaString &operator=(const T &o)
	{
		install(o.c_str(), o.length());
		return *this;
	}
	template <typename T>
	S2bDnaString &operator=(const std::basic_string<char> &o)
	{
		install(o);
		return *this;
	}
	void resize(size_t sz)
	{
		if (cs_ != NULL)
		{
			delete cs_;
			cs_ = NULL;
		}
		if (printcs_ != NULL)
		{
			delete printcs_;
			printcs_ = NULL;
		}
		len_ = sz;
		if (sz != 0)
		{
			cs_ = new uint32_t[nwords()];
		}
	}
	char toChar(size_t idx) const
	{
		int c = (int)get(idx);
		assert_range(0, 3, c);
		return "ACGT"[c];
	}
	char toColor(size_t idx) const
	{
		int c = (int)get(idx);
		assert_range(0, 3, c);
		return "0123"[c];
	}
	char windowGet(
		size_t i,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_lt(i, len);
		assert_leq(len, len_ - depth);
		return fw ? get(depth + i) : get(depth + len - i - 1);
	}
	template <typename T>
	void windowGet(
		T &ret,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_leq(len, len_ - depth);
		ret.resize(len);
		for (size_t i = 0; i < len; i++)
		{
			ret.set((fw ? get(depth + i) : get(depth + len - i - 1)), i);
		}
	}
	size_t nwords() const
	{
		return (len_ + 15) >> 4;
	}
	void set(int c, size_t idx)
	{
		assert_lt(idx, len_);
		assert_range(0, 3, c);
		size_t word = idx >> 4;
		size_t bpoff = (idx & 15) << 1;
		cs_[word] = cs_[word] & ~(uint32_t)(3 << bpoff);
		cs_[word] = cs_[word] | (uint32_t)(c << bpoff);
	}
	void setChar(int c, size_t idx)
	{
		assert_in(toupper(c), "ACGT");
		int bp = asc2dna[c];
		set(bp, idx);
	}
	void setColor(int c, size_t idx)
	{
		assert_in(toupper(c), "0123");
		int co = asc2col[c];
		set(co, idx);
	}
	void setWord(uint32_t w, size_t i)
	{
		assert_lt(i, nwords());
		cs_[i] = w;
	}
	char operator[](size_t i) const
	{
		assert_lt(i, len_);
		return get(i);
	}
	char get(size_t i) const
	{
		assert_lt(i, len_);
		size_t word = i >> 4;
		size_t bpoff = (i & 15) << 1;
		return (char)((cs_[word] >> bpoff) & 3);
	}
	void install(const uint32_t *b, size_t sz)
	{
		if (sz == 0)
			return;
		resize(sz);
		memcpy(cs_, b, sizeof(uint32_t) * nwords());
	}
	void install(const char *b, size_t sz)
	{
		if (sz == 0)
			return;
		resize(sz);
		size_t wordi = 0;
		for (size_t i = 0; i < sz; i += 16)
		{
			uint32_t word = 0;
			for (int j = 0; j < 16 && (size_t)(i + j) < sz; j++)
			{
				uint32_t bp = (int)b[i + j];
				uint32_t shift = (uint32_t)j << 1;
				assert_range(0, 3, (int)bp);
				word |= (bp << shift);
			}
			cs_[wordi++] = word;
		}
	}
	void installChars(const char *b, size_t sz)
	{
		if (sz == 0)
			return;
		resize(sz);
		size_t wordi = 0;
		for (size_t i = 0; i < sz; i += 16)
		{
			uint32_t word = 0;
			for (int j = 0; j < 16 && (size_t)(i + j) < sz; j++)
			{
				char c = b[i + j];
				assert_in(toupper(c), "ACGT");
				int bp = asc2dna[(int)c];
				assert_range(0, 3, (int)bp);
				uint32_t shift = (uint32_t)j << 1;
				word |= (bp << shift);
			}
			cs_[wordi++] = word;
		}
	}
	void installColors(const char *b, size_t sz)
	{
		if (sz == 0)
			return;
		resize(sz);
		size_t wordi = 0;
		for (size_t i = 0; i < sz; i += 16)
		{
			uint32_t word = 0;
			for (int j = 0; j < 16 && (size_t)(i + j) < sz; j++)
			{
				char c = b[i + j];
				assert_in(c, "0123");
				int bp = asc2col[(int)c];
				assert_range(0, 3, (int)bp);
				uint32_t shift = (uint32_t)j << 1;
				word |= (bp << shift);
			}
			cs_[wordi++] = word;
		}
	}
	void install(const char *b)
	{
		install(b, strlen(b));
	}
	void installChars(const char *b)
	{
		installChars(b, strlen(b));
	}
	void installColors(const char *b)
	{
		installColors(b, strlen(b));
	}
	void install(const std::basic_string<char> &b)
	{
		install(b.c_str(), b.length());
	}
	void installChars(const std::basic_string<char> &b)
	{
		installChars(b.c_str(), b.length());
	}
	void installColors(const std::basic_string<char> &b)
	{
		installColors(b.c_str(), b.length());
	}
	void installReverse(const char *b, size_t sz)
	{
		resize(sz);
		if (sz == 0)
			return;
		size_t wordi = 0;
		size_t bpi = 0;
		cs_[0] = 0;
		for (size_t i = sz; i > 0; i--)
		{
			assert_range(0, 3, (int)b[i - 1]);
			cs_[wordi] |= ((int)b[i - 1] << (bpi << 1));
			if (bpi == 15)
			{
				wordi++;
				cs_[wordi] = 0;
				bpi = 0;
			}
			else
				bpi++;
		}
	}
	void installReverse(const char *b)
	{
		installReverse(b, strlen(b));
	}
	void installReverseChars(const char *b, size_t sz)
	{
		resize(sz);
		if (sz == 0)
			return;
		size_t wordi = 0;
		size_t bpi = 0;
		cs_[0] = 0;
		for (size_t i = sz; i > 0; i--)
		{
			char c = b[i - 1];
			assert_in(toupper(c), "ACGT");
			int bp = asc2dna[(int)c];
			assert_range(0, 3, bp);
			cs_[wordi] |= (bp << (bpi << 1));
			if (bpi == 15)
			{
				wordi++;
				cs_[wordi] = 0;
				bpi = 0;
			}
			else
				bpi++;
		}
	}
	void installReverseChars(const char *b)
	{
		installReverseChars(b, strlen(b));
	}
	void installReverseColors(const char *b, size_t sz)
	{
		resize(sz);
		if (sz == 0)
			return;
		size_t wordi = 0;
		size_t bpi = 0;
		cs_[0] = 0;
		for (size_t i = sz; i > 0; i--)
		{
			char c = b[i - 1];
			assert_in(c, "0123");
			int bp = asc2col[(int)c];
			assert_range(0, 3, bp);
			cs_[wordi] |= (bp << (bpi << 1));
			if (bpi == 15)
			{
				wordi++;
				cs_[wordi] = 0;
				bpi = 0;
			}
			else
				bpi++;
		}
	}
	void installReverseColors(const char *b)
	{
		installReverseColors(b, strlen(b));
	}
	void installReverse(const S2bDnaString &b)
	{
		resize(b.len_);
		if (b.len_ == 0)
			return;
		size_t wordi = 0;
		size_t bpi = 0;
		size_t wordb = b.nwords() - 1;
		size_t bpb = (b.len_ - 1) & 15;
		cs_[0] = 0;
		for (size_t i = b.len_; i > 0; i--)
		{
			int bbp = (int)((b[wordb] >> (bpb << 1)) & 3);
			assert_range(0, 3, bbp);
			cs_[wordi] |= (bbp << (bpi << 1));
			if (bpi == 15)
			{
				wordi++;
				cs_[wordi] = 0;
				bpi = 0;
			}
			else
				bpi++;
			if (bpb == 0)
			{
				wordb--;
				bpi = 15;
			}
			else
				bpi--;
		}
	}
	bool operator==(const S2bDnaString &o)
	{
		return sstr_eq(*this, o);
	}
	bool operator!=(const S2bDnaString &o)
	{
		return sstr_neq(*this, o);
	}
	bool operator<(const S2bDnaString &o)
	{
		return sstr_lt(*this, o);
	}
	bool operator>(const S2bDnaString &o)
	{
		return sstr_gt(*this, o);
	}
	bool operator<=(const S2bDnaString &o)
	{
		return sstr_leq(*this, o);
	}
	bool operator>=(const S2bDnaString &o)
	{
		return sstr_geq(*this, o);
	}
	void reverse()
	{
		if (len_ <= 1)
			return;
		size_t wordf = nwords() - 1;
		size_t bpf = (len_ - 1) & 15;
		size_t wordi = 0;
		size_t bpi = 0;
		while (wordf > wordi || (wordf == wordi && bpf > bpi))
		{
			int f = (cs_[wordf] >> (bpf << 1)) & 3;
			int i = (cs_[wordi] >> (bpi << 1)) & 3;
			cs_[wordf] &= ~(uint32_t)(3 << (bpf << 1));
			cs_[wordi] &= ~(uint32_t)(3 << (bpi << 1));
			cs_[wordf] |= (uint32_t)(i << (bpf << 1));
			cs_[wordi] |= (uint32_t)(f << (bpi << 1));
			if (bpf == 0)
			{
				bpf = 15;
				wordf--;
			}
			else
				bpf--;
			if (bpi == 15)
			{
				bpi = 0;
				wordi++;
			}
			else
				bpi++;
		}
	}
	void reverseWindow(size_t off, size_t len)
	{
		assert_leq(off, len_);
		assert_leq(off + len, len_);
		if (len <= 1)
			return;
		size_t wordf = (off + len - 1) >> 4;
		size_t bpf = (off + len - 1) & 15;
		size_t wordi = (off) >> 4;
		size_t bpi = (off)&15;
		while (wordf > wordi || (wordf == wordi && bpf > bpi))
		{
			int f = (cs_[wordf] >> (bpf << 1)) & 3;
			int i = (cs_[wordi] >> (bpi << 1)) & 3;
			cs_[wordf] &= ~(uint32_t)(3 << (bpf << 1));
			cs_[wordi] &= ~(uint32_t)(3 << (bpi << 1));
			cs_[wordf] |= (uint32_t)(i << (bpf << 1));
			cs_[wordi] |= (uint32_t)(f << (bpi << 1));
			if (bpf == 0)
			{
				bpf = 15;
				wordf--;
			}
			else
				bpf--;
			if (bpi == 15)
			{
				bpi = 0;
				wordi++;
			}
			else
				bpi++;
		}
	}
	void fill(size_t len, char el)
	{
		assert_leq(len, len_);
		assert_range(0, 3, (int)el);
		size_t word = 0;
		if (len > 32)
		{
			uint32_t bl = (uint32_t)el;
			bl |= (bl << 2);
			bl |= (bl << 4);
			bl |= (bl << 8);
			bl |= (bl << 16);
			size_t blen = len >> 4;
			for (; word < blen; word++)
			{
				cs_[word] = bl;
			}
			len = len & 15;
		}
		size_t bp = 0;
		for (size_t i = 0; i < len; i++)
		{
			cs_[word] &= ~(uint32_t)(3 << (bp << 1));
			cs_[word] |= (uint32_t)(el << (bp << 1));
			if (bp == 15)
			{
				word++;
				bp = 0;
			}
			else
				bp++;
		}
	}
	void fill(char el)
	{
		fill(len_, el);
	}
	char windowGetDna(
		size_t i,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_lt(i, len);
		assert_leq(len, len_ - depth);
		if (fw)
		{
			return get(depth + i);
		}
		else
		{
			return compDna(get(depth + len - i - 1));
		}
	}
	template <typename T>
	void windowGetDna(
		T &buf,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_leq(len, len_ - depth);
		buf.resize(len);
		for (size_t i = 0; i < len; i++)
		{
			buf.set(
				(fw ? get(depth + i) : compDna(get(depth + len - i - 1))), i);
		}
	}
	inline size_t length() const { return len_; }
	void clear() { len_ = 0; }
	inline bool empty() const { return len_ == 0; }
	const uint32_t *buf() const { return cs_; }
	uint32_t *wbuf() { return cs_; }
	const char *toZBuf() const
	{
		if (printcs_ == NULL)
		{
			const_cast<char *&>(printcs_) = new char[len_ + 1];
		}
		char *printcs = const_cast<char *>(printcs_);
		size_t word = 0, bp = 0;
		for (size_t i = 0; i < len_; i++)
		{
			int c = (cs_[word] >> (bp << 1)) & 3;
			printcs[i] = "ACGT"[c];
			if (bp == 15)
			{
				word++;
				bp = 0;
			}
			else
				bp++;
		}
		printcs[len_] = '\0';
		return printcs_;
	}
protected:
	uint32_t *cs_;
	char *printcs_;
	size_t len_;
};
template <typename T, int S = 1024, int M = 2, int I = 0>
class SStringExpandable
{
public:
	explicit SStringExpandable() : cs_(NULL),
								   printcs_(NULL),
								   len_(0),
								   sz_(0)
	{
		if (I > 0)
		{
			expandNoCopy(I);
		}
	}
	explicit SStringExpandable(size_t sz) : cs_(NULL),
											printcs_(NULL),
											len_(0),
											sz_(0)
	{
		expandNoCopy(sz);
	}
	SStringExpandable(const SStringExpandable<T, S> &o) : cs_(NULL),
														  printcs_(NULL),
														  len_(0),
														  sz_(0)
	{
		*this = o;
	}
	explicit SStringExpandable(const std::basic_string<T> &str) : cs_(NULL),
																  printcs_(NULL),
																  len_(0),
																  sz_(0)
	{
		install(str.c_str(), str.length());
	}
	explicit SStringExpandable(const T *b, size_t sz) : cs_(NULL),
														printcs_(NULL),
														len_(0),
														sz_(0)
	{
		install(b, sz);
	}
	explicit SStringExpandable(const T *b) : cs_(NULL),
											 printcs_(NULL),
											 len_(0),
											 sz_(0)
	{
		install(b, strlen(b));
	}
	virtual ~SStringExpandable()
	{
		if (cs_ != NULL)
		{
			delete[] cs_;
			cs_ = NULL;
		}
		if (printcs_ != NULL)
		{
			delete[] printcs_;
			printcs_ = NULL;
		}
		sz_ = len_ = 0;
	}
	T windowGet(
		size_t i,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_lt(i, len);
		assert_leq(len, len_ - depth);
		return fw ? cs_[depth + i] : cs_[depth + len - i - 1];
	}
	void windowGet(
		T &ret,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_leq(len, len_ - depth);
		for (size_t i = 0; i < len; i++)
		{
			ret.append(fw ? cs_[depth + i] : cs_[depth + len - i - 1]);
		}
	}
	SStringExpandable<T, S, M, I> &operator=(const SStringExpandable<T, S, M, I> &o)
	{
		install(o.cs_, o.len_);
		return *this;
	}
	SStringExpandable<T, S> &operator=(const std::basic_string<T> &o)
	{
		install(o.c_str(), o.length());
		return *this;
	}
	void insert(const T &c, size_t idx)
	{
		assert_lt(idx, len_);
		if (sz_ < len_ + 1)
			expandCopy((len_ + 1 + S) * M);
		len_++;
		for (size_t i = len_; i > idx + 1; i--)
		{
			cs_[i - 1] = cs_[i - 2];
		}
		cs_[idx] = c;
	}
	void set(int c, size_t idx)
	{
		assert_lt(idx, len_);
		cs_[idx] = c;
	}
	void append(const T &c)
	{
		if (sz_ < len_ + 1)
			expandCopy((len_ + 1 + S) * M);
		cs_[len_++] = c;
	}
	void remove(size_t idx)
	{
		assert_lt(idx, len_);
		assert_gt(len_, 0);
		for (size_t i = idx; i < len_ - 1; i++)
		{
			cs_[i] = cs_[i + 1];
		}
		len_--;
	}
	const T &operator[](size_t i) const
	{
		assert_lt(i, len_);
		return cs_[i];
	}
	T &operator[](size_t i)
	{
		assert_lt(i, len_);
		return cs_[i];
	}
	const T &get(size_t i) const
	{
		assert_lt(i, len_);
		return cs_[i];
	}
	virtual void install(const T *b, size_t sz)
	{
		if (sz_ < sz)
			expandNoCopy((sz + S) * M);
		memcpy(cs_, b, sz * sizeof(T));
		len_ = sz;
	}
	void install(const T *b) { install(b, strlen(b)); }
	void installReverse(const char *b, size_t sz)
	{
		if (sz_ < sz)
			expandNoCopy((sz + S) * M);
		for (size_t i = 0; i < sz; i++)
		{
			cs_[i] = b[sz - i - 1];
		}
		len_ = sz;
	}
	void installReverse(const SStringExpandable<T, S> &b)
	{
		if (sz_ < b.len_)
			expandNoCopy((b.len_ + S) * M);
		for (size_t i = 0; i < b.len_; i++)
		{
			cs_[i] = b.cs_[b.len_ - i - 1];
		}
		len_ = b.len_;
	}
	bool operator==(const SStringExpandable<T, S> &o)
	{
		return sstr_eq(*this, o);
	}
	bool operator!=(const SStringExpandable<T, S> &o)
	{
		return sstr_neq(*this, o);
	}
	bool operator<(const SStringExpandable<T, S> &o)
	{
		return sstr_lt(*this, o);
	}
	bool operator>(const SStringExpandable<T, S> &o)
	{
		return sstr_gt(*this, o);
	}
	bool operator<=(const SStringExpandable<T, S> &o)
	{
		return sstr_leq(*this, o);
	}
	bool operator>=(const SStringExpandable<T, S> &o)
	{
		return sstr_geq(*this, o);
	}
	void reverse()
	{
		for (size_t i = 0; i < (len_ >> 1); i++)
		{
			T tmp = get(i);
			set(get(len_ - i - 1), i);
			set(tmp, len_ - i - 1);
		}
	}
	void reverseWindow(size_t off, size_t len)
	{
		assert_leq(off, len_);
		assert_leq(off + len, len_);
		size_t mid = len >> 1;
		for (size_t i = 0; i < mid; i++)
		{
			T tmp = get(off + i);
			set(get(off + len - i - 1), off + i);
			set(tmp, off + len - i - 1);
		}
	}
	void resize(size_t len)
	{
		if (sz_ < len)
			expandCopy((len + S) * M);
		len_ = len;
	}
	void resize(size_t len, const T &el)
	{
		if (sz_ < len)
			expandCopy((len + S) * M);
		if (len > len_)
		{
			for (size_t i = len_; i < len; i++)
			{
				cs_[i] = el;
			}
		}
		len_ = len;
	}
	void fill(size_t len, const T &el)
	{
		assert_leq(len, len_);
		for (size_t i = 0; i < len; i++)
		{
			cs_[i] = el;
		}
	}
	void fill(const T &el)
	{
		fill(len_, el);
	}
	void trimBegin(size_t len)
	{
		assert_leq(len, len_);
		if (len == len_)
		{
			len_ = 0;
			return;
		}
		for (size_t i = 0; i < len_ - len; i++)
		{
			cs_[i] = cs_[i + len];
		}
		len_ -= len;
	}
	size_t trimEnd(size_t len)
	{
		if (len >= len_)
		{
			size_t ret = len_;
			len_ = 0;
			return ret;
		}
		len_ -= len;
		return len;
	}
	void append(const T *b, size_t sz)
	{
		if (sz_ < len_ + sz)
			expandCopy((len_ + sz + S) * M);
		memcpy(cs_ + len_, b, sz * sizeof(T));
		len_ += sz;
	}
	void append(const T *b)
	{
		append(b, strlen(b));
	}
	size_t length() const { return len_; }
	void clear() { len_ = 0; }
	bool empty() const { return len_ == 0; }
	const char *toZBufXForm(const char *xform) const
	{
		ASSERT_ONLY(size_t xformElts = strlen(xform));
		if (empty())
		{
			const_cast<char &>(zero_) = 0;
			return &zero_;
		}
		char *printcs = const_cast<char *>(printcs_);
		for (size_t i = 0; i < len_; i++)
		{
			assert_lt(cs_[i], (int)xformElts);
			printcs[i] = xform[(int)cs_[i]];
		}
		printcs[len_] = 0;
		return printcs_;
	}
	virtual const T *toZBuf() const
	{
		if (empty())
		{
			const_cast<T &>(zeroT_) = 0;
			return &zeroT_;
		}
		assert_leq(len_, sz_);
		const_cast<T *>(cs_)[len_] = 0;
		return cs_;
	}
	bool eq(const char *str) const
	{
		const char *self = toZBuf();
		return strcmp(str, self) == 0;
	}
	const T *buf() const { return cs_; }
	T *wbuf() { return cs_; }
protected:
	void expandCopy(size_t sz)
	{
		if (sz_ >= sz)
			return;
		T *tmp = new T[sz + 1];
		char *ptmp = new char[sz + 1];
		if (cs_ != NULL)
		{
			memcpy(tmp, cs_, sizeof(T) * len_);
			delete[] cs_;
		}
		if (printcs_ != NULL)
		{
			memcpy(ptmp, printcs_, sizeof(char) * len_);
			delete[] printcs_;
		}
		cs_ = tmp;
		printcs_ = ptmp;
		sz_ = sz;
	}
	void expandNoCopy(size_t sz)
	{
		if (sz_ >= sz)
			return;
		if (cs_ != NULL)
			delete[] cs_;
		if (printcs_ != NULL)
			delete[] printcs_;
		cs_ = new T[sz + 1];
		printcs_ = new char[sz + 1];
		sz_ = sz;
	}
	T *cs_;
	char *printcs_;
	char zero_;
	T zeroT_;
	size_t len_;
	size_t sz_;
};
template <typename T, int S>
class SStringFixed
{
public:
	explicit SStringFixed() : len_(0) {}
	SStringFixed(const SStringFixed<T, S> &o)
	{
		*this = o;
	}
	explicit SStringFixed(const std::basic_string<T> &str)
	{
		install(str.c_str(), str.length());
	}
	explicit SStringFixed(const T *b, size_t sz)
	{
		install(b, sz);
	}
	explicit SStringFixed(const T *b)
	{
		install(b, strlen(b));
	}
	virtual ~SStringFixed() {}
	inline const T &operator[](size_t i) const
	{
		return get(i);
	}
	inline T &operator[](size_t i)
	{
		return get(i);
	}
	inline const T &get(size_t i) const
	{
		assert_lt(i, len_);
		return cs_[i];
	}
	inline T &get(size_t i)
	{
		assert_lt(i, len_);
		return cs_[i];
	}
	T windowGet(
		size_t i,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_lt(i, len);
		assert_leq(len, len_ - depth);
		return fw ? cs_[depth + i] : cs_[depth + len - i - 1];
	}
	void windowGet(
		T &ret,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = len_;
		assert_leq(len, len_ - depth);
		for (size_t i = 0; i < len; i++)
		{
			ret.append(fw ? cs_[depth + i] : cs_[depth + len - i - 1]);
		}
	}
	SStringFixed<T, S> &operator=(const SStringFixed<T, S> &o)
	{
		install(o.cs_, o.len_);
		return *this;
	}
	SStringFixed<T, S> &operator=(const std::basic_string<T> &o)
	{
		install(o);
		return *this;
	}
	void insert(const T &c, size_t idx)
	{
		assert_lt(len_, S);
		assert_lt(idx, len_);
		for (int i = len_; i > idx; i--)
		{
			cs_[i] = cs_[i - 1];
		}
		cs_[idx] = c;
		len_++;
	}
	void set(int c, size_t idx)
	{
		assert_lt(idx, len_);
		cs_[idx] = c;
	}
	void append(const T &c)
	{
		assert_lt(len_, S);
		cs_[len_++] = c;
	}
	void remove(size_t idx)
	{
		assert_lt(idx, len_);
		assert_gt(len_, 0);
		for (size_t i = idx; i < len_ - 1; i++)
		{
			cs_[i] = cs_[i + 1];
		}
		len_--;
	}
	virtual void install(const T *b, size_t sz)
	{
		assert_leq(sz, S);
		memcpy(cs_, b, sz * sizeof(T));
		len_ = sz;
	}
	void install(const T *b) { install(b, strlen(b)); }
	void installReverse(const char *b, size_t sz)
	{
		assert_leq(sz, S);
		for (size_t i = 0; i < sz; i++)
		{
			cs_[i] = b[sz - i - 1];
		}
		len_ = sz;
	}
	void installReverse(const SStringFixed<T, S> &b)
	{
		assert_leq(b.len_, S);
		for (size_t i = 0; i < b.len_; i++)
		{
			cs_[i] = b.cs_[b.len_ - i - 1];
		}
		len_ = b.len_;
	}
	bool operator==(const SStringFixed<T, S> &o)
	{
		return sstr_eq(*this, o);
	}
	bool operator!=(const SStringFixed<T, S> &o)
	{
		return sstr_neq(*this, o);
	}
	bool operator<(const SStringFixed<T, S> &o)
	{
		return sstr_lt(*this, o);
	}
	bool operator>(const SStringFixed<T, S> &o)
	{
		return sstr_gt(*this, o);
	}
	bool operator<=(const SStringFixed<T, S> &o)
	{
		return sstr_leq(*this, o);
	}
	bool operator>=(const SStringFixed<T, S> &o)
	{
		return sstr_geq(*this, o);
	}
	void reverse()
	{
		for (size_t i = 0; i < (len_ >> 1); i++)
		{
			T tmp = get(i);
			set(get(len_ - i - 1), i);
			set(tmp, len_ - i - 1);
		}
	}
	void reverseWindow(size_t off, size_t len)
	{
		assert_leq(off, len_);
		assert_leq(off + len, len_);
		size_t mid = len >> 1;
		for (size_t i = 0; i < mid; i++)
		{
			T tmp = get(off + i);
			set(get(off + len - i - 1), off + i);
			set(tmp, off + len - i - 1);
		}
	}
	void resize(size_t len)
	{
		assert_lt(len, S);
		len_ = len;
	}
	void resize(size_t len, const T &el)
	{
		assert_lt(len, S);
		if (len > len_)
		{
			for (size_t i = len_; i < len; i++)
			{
				cs_[i] = el;
			}
		}
		len_ = len;
	}
	void fill(size_t len, const T &el)
	{
		assert_leq(len, len_);
		for (size_t i = 0; i < len; i++)
		{
			cs_[i] = el;
		}
	}
	void fill(const T &el)
	{
		fill(len_, el);
	}
	void trimBegin(size_t len)
	{
		assert_leq(len, len_);
		if (len == len_)
		{
			len_ = 0;
			return;
		}
		for (size_t i = 0; i < len_ - len; i++)
		{
			cs_[i] = cs_[i + len];
		}
		len_ -= len;
	}
	void trimEnd(size_t len)
	{
		if (len >= len_)
			len_ = 0;
		else
			len_ -= len;
	}
	void append(const T *b, size_t sz)
	{
		assert_leq(sz + len_, S);
		memcpy(cs_ + len_, b, sz * sizeof(T));
		len_ += sz;
	}
	void append(const T *b)
	{
		append(b, strlen(b));
	}
	size_t length() const { return len_; }
	void clear() { len_ = 0; }
	bool empty() const { return len_ == 0; }
	virtual const T *toZBuf() const
	{
		const_cast<T *>(cs_)[len_] = 0;
		return cs_;
	}
	bool eq(const char *str) const
	{
		const char *self = toZBuf();
		return strcmp(str, self) == 0;
	}
	const char *toZBufXForm(const char *xform) const
	{
		ASSERT_ONLY(size_t xformElts = strlen(xform));
		char *printcs = const_cast<char *>(printcs_);
		for (size_t i = 0; i < len_; i++)
		{
			assert_lt(cs_[i], (int)xformElts);
			printcs[i] = xform[cs_[i]];
		}
		printcs[len_] = 0;
		return printcs_;
	}
	const T *buf() const { return cs_; }
	T *wbuf() { return cs_; }
protected:
	T cs_[S + 1];
	char printcs_[S + 1];
	size_t len_;
};
template <typename T, int S, int M>
std::ostream &operator<<(std::ostream &os, const SStringExpandable<T, S, M> &str)
{
	os << str.toZBuf();
	return os;
}
template <typename T, int S>
std::ostream &operator<<(std::ostream &os, const SStringFixed<T, S> &str)
{
	os << str.toZBuf();
	return os;
}
extern uint8_t asc2dna[];
extern uint8_t asc2col[];
template <int S>
class SDnaStringFixed : public SStringFixed<char, S>
{
public:
	explicit SDnaStringFixed() : SStringFixed<char, S>() {}
	SDnaStringFixed(const SDnaStringFixed<S> &o) : SStringFixed<char, S>(o) {}
	explicit SDnaStringFixed(const std::basic_string<char> &str) : SStringFixed<char, S>(str) {}
	explicit SDnaStringFixed(const char *b, size_t sz) : SStringFixed<char, S>(b, sz) {}
	explicit SDnaStringFixed(
		const char *b,
		bool chars = false,
		bool colors = false) : SStringFixed<char, S>()
	{
		if (chars)
		{
			if (colors)
			{
				installColors(b, strlen(b));
			}
			else
			{
				installChars(b, strlen(b));
			}
		}
		else
		{
			install(b, strlen(b));
		}
	}
	virtual ~SDnaStringFixed() {}
	void installReverseComp(const char *b, size_t sz)
	{
		assert_leq(sz, S);
		for (size_t i = 0; i < sz; i++)
		{
			this->cs_[i] = (b[sz - i - 1] == 4 ? 4 : b[sz - i - 1] ^ 3);
		}
		this->len_ = sz;
	}
	void installReverseComp(const SDnaStringFixed<S> &b)
	{
		assert_leq(b.len_, S);
		for (size_t i = 0; i < b.len_; i++)
		{
			this->cs_[i] = (b.cs_[b.len_ - i - 1] == 4 ? 4 : b.cs_[b.len_ - i - 1] ^ 3);
		}
		this->len_ = b.len_;
	}
	void reverseComp()
	{
		for (size_t i = 0; i < (this->len_ >> 1); i++)
		{
			char tmp1 = (this->cs_[i] == 4 ? 4 : this->cs_[i] ^ 3);
			char tmp2 = (this->cs_[this->len_ - i - 1] == 4 ? 4 : this->cs_[this->len_ - i - 1] ^ 3);
			this->cs_[i] = tmp2;
			this->cs_[this->len_ - i - 1] = tmp1;
		}
		if ((this->len_ & 1) != 0)
		{
			char tmp = this->cs_[this->len_ >> 1];
			tmp = (tmp == 4 ? 4 : tmp ^ 3);
			this->cs_[this->len_ >> 1] = tmp;
		}
	}
	virtual void install(const char *b, size_t sz)
	{
		assert_leq(sz, S);
		memcpy(this->cs_, b, sz);
#ifndef NDEBUG
		for (size_t i = 0; i < sz; i++)
		{
			assert_leq(this->cs_[i], 4);
			assert_geq(this->cs_[i], 0);
		}
#endif
		this->len_ = sz;
	}
	virtual void installChars(const char *b, size_t sz)
	{
		assert_leq(sz, S);
		for (size_t i = 0; i < sz; i++)
		{
			assert_in(toupper(b[i]), "ACGTN-");
			this->cs_[i] = asc2dna[(int)b[i]];
			assert_geq(this->cs_[i], 0);
			assert_leq(this->cs_[i], 4);
		}
		this->len_ = sz;
	}
	virtual void installColors(const char *b, size_t sz)
	{
		assert_leq(sz, S);
		for (size_t i = 0; i < sz; i++)
		{
			assert_in(b[i], "0123.");
			this->cs_[i] = asc2col[(int)b[i]];
			assert_geq(this->cs_[i], 0);
			assert_leq(this->cs_[i], 4);
		}
		this->len_ = sz;
	}
	virtual void installChars(const std::basic_string<char> &str)
	{
		installChars(str.c_str(), str.length());
	}
	virtual void installColors(const std::basic_string<char> &str)
	{
		installColors(str.c_str(), str.length());
	}
	void set(int c, size_t idx)
	{
		assert_lt(idx, this->len_);
		assert_leq(c, 4);
		assert_geq(c, 0);
		this->cs_[idx] = c;
	}
	void append(const char &c)
	{
		assert_lt(this->len_, S);
		assert_leq(c, 4);
		assert_geq(c, 0);
		this->cs_[this->len_++] = c;
	}
	void setChar(char c, size_t idx)
	{
		assert_lt(idx, this->len_);
		assert_in(toupper(c), "ACGTN");
		this->cs_[idx] = asc2dna[(int)c];
	}
	void appendChar(char c)
	{
		assert_lt(this->len_, S);
		assert_in(toupper(c), "ACGTN");
		this->cs_[this->len_++] = asc2dna[(int)c];
	}
	char toChar(size_t idx) const
	{
		assert_geq((int)this->cs_[idx], 0);
		assert_leq((int)this->cs_[idx], 4);
		return "ACGTN"[(int)this->cs_[idx]];
	}
	const char &operator[](size_t i) const
	{
		return this->get(i);
	}
	const char &get(size_t i) const
	{
		assert_lt(i, this->len_);
		assert_leq(this->cs_[i], 4);
		assert_geq(this->cs_[i], 0);
		return this->cs_[i];
	}
	char windowGetDna(
		size_t i,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = this->len_;
		assert_lt(i, len);
		assert_leq(len, this->len_ - depth);
		if (fw)
			return this->cs_[depth + i];
		else
			return compDna(this->cs_[depth + len - i - 1]);
	}
	void windowGetDna(
		SDnaStringFixed<S> &buf,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = this->len_;
		assert_leq(len, this->len_ - depth);
		for (size_t i = 0; i < len; i++)
		{
			buf.append(fw ? this->cs_[depth + i] : compDna(this->cs_[depth + len - i - 1]));
		}
	}
	virtual const char *toZBuf() const { return this->toZBufXForm("ACGTN"); }
};
template <int S = 1024, int M = 2>
class SDnaStringExpandable : public SStringExpandable<char, S, M>
{
public:
	explicit SDnaStringExpandable() : SStringExpandable<char, S, M>() {}
	SDnaStringExpandable(const SDnaStringExpandable<S, M> &o) : SStringExpandable<char, S, M>(o) {}
	explicit SDnaStringExpandable(
		const std::basic_string<char> &str,
		bool chars = false,
		bool colors = false) : SStringExpandable<char, S, M>()
	{
		if (chars)
		{
			if (colors)
			{
				installColors(str);
			}
			else
			{
				installChars(str);
			}
		}
		else
		{
			install(str);
		}
	}
	explicit SDnaStringExpandable(
		const char *b,
		size_t sz,
		bool chars = false,
		bool colors = false) : SStringExpandable<char, S, M>()
	{
		if (chars)
		{
			if (colors)
			{
				installColors(b, sz);
			}
			else
			{
				installChars(b, sz);
			}
		}
		else
		{
			install(b, sz);
		}
	}
	explicit SDnaStringExpandable(
		const char *b,
		bool chars = false,
		bool colors = false) : SStringExpandable<char, S, M>()
	{
		install(b, chars, colors);
	}
	virtual ~SDnaStringExpandable() {}
	void installReverseComp(const char *b, size_t sz)
	{
		if (this->sz_ < sz)
			this->expandCopy((sz + S) * M);
		for (size_t i = 0; i < sz; i++)
		{
			this->cs_[i] = (b[sz - i - 1] == 4 ? 4 : b[sz - i - 1] ^ 3);
		}
		this->len_ = sz;
	}
	void installReverseComp(const SDnaStringExpandable<S, M> &b)
	{
		if (this->sz_ < b.len_)
			this->expandCopy((b.len_ + S) * M);
		for (size_t i = 0; i < b.len_; i++)
		{
			this->cs_[i] = (b.cs_[b.len_ - i - 1] == 4 ? 4 : b.cs_[b.len_ - i - 1] ^ 3);
		}
		this->len_ = b.len_;
	}
	void reverseComp()
	{
		for (size_t i = 0; i < (this->len_ >> 1); i++)
		{
			char tmp1 = (this->cs_[i] == 4 ? 4 : this->cs_[i] ^ 3);
			char tmp2 = (this->cs_[this->len_ - i - 1] == 4 ? 4 : this->cs_[this->len_ - i - 1] ^ 3);
			this->cs_[i] = tmp2;
			this->cs_[this->len_ - i - 1] = tmp1;
		}
		if ((this->len_ & 1) != 0)
		{
			char tmp = this->cs_[this->len_ >> 1];
			tmp = (tmp == 4 ? 4 : tmp ^ 3);
			this->cs_[this->len_ >> 1] = tmp;
		}
	}
	virtual void install(
		const char *b,
		bool chars = false,
		bool colors = false)
	{
		if (chars)
		{
			if (colors)
			{
				installColors(b, strlen(b));
			}
			else
			{
				installChars(b, strlen(b));
			}
		}
		else
		{
			install(b, strlen(b));
		}
	}
	virtual void install(const char *b, size_t sz)
	{
		if (this->sz_ < sz)
			this->expandCopy((sz + S) * M);
		memcpy(this->cs_, b, sz);
#ifndef NDEBUG
		for (size_t i = 0; i < sz; i++)
		{
			assert_range(0, 4, (int)this->cs_[i]);
		}
#endif
		this->len_ = sz;
	}
	virtual void installChars(const char *b, size_t sz)
	{
		if (this->sz_ < sz)
			this->expandCopy((sz + S) * M);
		for (size_t i = 0; i < sz; i++)
		{
			assert_in(toupper(b[i]), "ACGTN-");
			this->cs_[i] = asc2dna[(int)b[i]];
			assert_range(0, 4, (int)this->cs_[i]);
		}
		this->len_ = sz;
	}
	virtual void installColors(const char *b, size_t sz)
	{
		if (this->sz_ < sz)
			this->expandCopy((sz + S) * M);
		for (size_t i = 0; i < sz; i++)
		{
			assert_in(b[i], "0123.");
			this->cs_[i] = asc2col[(int)b[i]];
			assert_range(0, 4, (int)this->cs_[i]);
		}
		this->len_ = sz;
	}
	virtual void installChars(const std::basic_string<char> &str)
	{
		installChars(str.c_str(), str.length());
	}
	virtual void installColors(const std::basic_string<char> &str)
	{
		installColors(str.c_str(), str.length());
	}
	void set(int c, size_t idx)
	{
		assert_lt(idx, this->len_);
		assert_range(0, 4, c);
		this->cs_[idx] = c;
	}
	void append(const char &c)
	{
		if (this->sz_ < this->len_ + 1)
		{
			this->expandCopy((this->len_ + 1 + S) * M);
		}
		assert_range(0, 4, (int)c);
		this->cs_[this->len_++] = c;
	}
	void setChar(char c, size_t idx)
	{
		assert_lt(idx, this->len_);
		assert_in(toupper(c), "ACGTN");
		this->cs_[idx] = asc2dna[(int)c];
	}
	void appendChar(char c)
	{
		if (this->sz_ < this->len_ + 1)
		{
			this->expandCopy((this->len_ + 1 + S) * M);
		}
		assert_in(toupper(c), "ACGTN");
		this->cs_[this->len_++] = asc2dna[(int)c];
	}
	char toChar(size_t idx) const
	{
		assert_range(0, 4, (int)this->cs_[idx]);
		return "ACGTN"[(int)this->cs_[idx]];
	}
	inline const char &operator[](size_t i) const
	{
		return this->get(i);
	}
	inline const char &get(size_t i) const
	{
		assert_lt(i, this->len_);
		assert_range(0, 4, (int)this->cs_[i]);
		return this->cs_[i];
	}
	char windowGetDna(
		size_t i,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = this->len_;
		assert_lt(i, len);
		assert_leq(len, this->len_ - depth);
		if (fw)
			return this->cs_[depth + i];
		else
			return compDna(this->cs_[depth + len - i - 1]);
	}
	void windowGetDna(
		SDnaStringExpandable<S, M> &buf,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = this->len_;
		assert_leq(len, this->len_ - depth);
		for (size_t i = 0; i < len; i++)
		{
			buf.append(fw ? this->cs_[depth + i] : compDna(this->cs_[depth + len - i - 1]));
		}
	}
	virtual const char *toZBuf() const { return this->toZBufXForm("ACGTN"); }
};
template <int S = 16, int M = 2>
class SDnaMaskString : public SStringExpandable<char, S, M>
{
public:
	explicit SDnaMaskString() : SStringExpandable<char, S, M>() {}
	SDnaMaskString(const SDnaMaskString<S, M> &o) : SStringExpandable<char, S, M>(o) {}
	explicit SDnaMaskString(const std::basic_string<char> &str) : SStringExpandable<char, S, M>(str) {}
	explicit SDnaMaskString(const char *b, size_t sz) : SStringExpandable<char, S, M>(b, sz) {}
	explicit SDnaMaskString(const char *b, bool chars = false) : SStringExpandable<char, S, M>()
	{
		if (chars)
		{
			installChars(b, strlen(b));
		}
		else
		{
			install(b, strlen(b));
		}
	}
	virtual ~SDnaMaskString() {}
	void installReverseComp(const char *b, size_t sz)
	{
		while (this->sz_ < sz)
		{
			this->expandNoCopy((sz + S) * M);
		}
		for (size_t i = 0; i < sz; i++)
		{
			this->cs_[i] = maskcomp[(int)b[sz - i - 1]];
		}
		this->len_ = sz;
	}
	void installReverseComp(const SDnaMaskString<S, M> &b)
	{
		while (this->sz_ < b.len_)
		{
			this->expandNoCopy((b.len_ + S) * M);
		}
		for (size_t i = 0; i < b.len_; i++)
		{
			this->cs_[i] = maskcomp[(int)b.cs_[b.len_ - i - 1]];
		}
		this->len_ = b.len_;
	}
	void reverseComp()
	{
		for (size_t i = 0; i < (this->len_ >> 1); i++)
		{
			char tmp1 = maskcomp[(int)this->cs_[i]];
			char tmp2 = maskcomp[(int)this->cs_[this->len_ - i - 1]];
			this->cs_[i] = tmp2;
			this->cs_[this->len_ - i - 1] = tmp1;
		}
		if ((this->len_ & 1) != 0)
		{
			char tmp = this->cs_[this->len_ >> 1];
			tmp = maskcomp[(int)tmp];
			this->cs_[this->len_ >> 1] = tmp;
		}
	}
	virtual void install(const char *b, size_t sz)
	{
		while (this->sz_ < sz)
		{
			this->expandNoCopy((sz + S) * M);
		}
		memcpy(this->cs_, b, sz);
#ifndef NDEBUG
		for (size_t i = 0; i < sz; i++)
		{
			assert_range((int)this->cs_[i], 0, 15);
		}
#endif
		this->len_ = sz;
	}
	virtual void installChars(const char *b, size_t sz)
	{
		while (this->sz_ < sz)
		{
			this->expandNoCopy((sz + S) * M);
		}
		for (size_t i = 0; i < sz; i++)
		{
			assert_in(b[i], iupacs);
			this->cs_[i] = asc2dnamask[(int)b[i]];
			assert_range((int)this->cs_[i], 0, 15);
		}
		this->len_ = sz;
	}
	virtual void installChars(const std::basic_string<char> &str)
	{
		installChars(str.c_str(), str.length());
	}
	void set(int c, size_t idx)
	{
		assert_lt(idx, this->len_);
		assert_range(c, 0, 15);
		this->cs_[idx] = c;
	}
	void append(const char &c)
	{
		while (this->sz_ < this->len_ + 1)
		{
			this->expandNoCopy((this->len_ + 1 + S) * M);
		}
		assert_range((int)c, 0, 15);
		this->cs_[this->len_++] = c;
	}
	void setChar(char c, size_t idx)
	{
		assert_lt(idx, this->len_);
		assert_in(toupper(c), iupacs);
		this->cs_[idx] = asc2dnamask[(int)c];
	}
	void appendChar(char c)
	{
		while (this->sz_ < this->len_ + 1)
		{
			expandNoCopy((this->len_ + 1 + S) * M);
		}
		assert_in(toupper(c), iupacs);
		this->cs_[this->len_++] = asc2dnamask[(int)c];
	}
	char toChar(size_t idx) const
	{
		assert_range((int)this->cs_[idx], 0, 15);
		return mask2iupac[(int)this->cs_[idx]];
	}
	const char &operator[](size_t i) const
	{
		return this->get(i);
	}
	char &operator[](size_t i)
	{
		return this->get(i);
	}
	const char &get(size_t i) const
	{
		assert_lt(i, this->len_);
		assert_range((int)this->cs_[i], 0, 15);
		return this->cs_[i];
	}
	char &get(size_t i)
	{
		assert_lt(i, this->len_);
		assert_range((int)this->cs_[i], 0, 15);
		return this->cs_[i];
	}
	char windowGetDna(
		size_t i,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = this->len_;
		assert_lt(i, len);
		assert_leq(len, this->len_ - depth);
		if (fw)
			return this->cs_[depth + i];
		else
			return maskcomp[this->cs_[depth + len - i - 1]];
	}
	void windowGetDna(
		SDnaStringFixed<S> &buf,
		bool fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if (len == 0)
			len = this->len_;
		assert_leq(len, this->len_ - depth);
		for (size_t i = 0; i < len; i++)
		{
			buf.append(fw ? this->cs_[depth + i] : maskcomp[this->cs_[depth + len - i - 1]]);
		}
	}
	template <typename T>
	void randSubstr(
		RandomSource &rnd,
		T &dst, size_t len, bool watson = true, bool crick = true)
	{
		assert(watson || crick);
		assert_geq(this->len_, len);
		size_t poss = this->len_ - len + 1;
		assert_gt(poss, 0);
		uint32_t rndoff = (uint32_t)(rnd.nextU32() % poss);
		bool fw;
		if (watson && !crick)
			fw = true;
		else if (!watson && crick)
			fw = false;
		else
		{
			fw = rnd.nextBool();
		}
		if (fw)
		{
			for (size_t i = 0; i < len; i++)
			{
				dst[i] = this->cs_[i + rndoff];
			}
		}
		else
		{
			for (size_t i = 0; i < len; i++)
			{
				dst[i] = maskcomp[(int)this->cs_[i + rndoff + (len - i - 1)]];
			}
		}
	}
	virtual const char *toZBuf() const { return this->toZBufXForm(iupacs); }
};
typedef SStringExpandable<char, 1024, 2> BTString;
typedef SDnaStringExpandable<1024, 2> BTDnaString;
typedef SDnaMaskString<32, 2> BTDnaMask;
#endif
