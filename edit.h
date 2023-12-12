#ifndef EDIT_H_
#define EDIT_H_
#include <iostream>
#include <stdint.h>
#include <limits>
#include "assert_helpers.h"
#include "filebuf.h"
#include "sstring.h"
#include "ds.h"
enum
{
	EDIT_TYPE_READ_GAP = 1,
	EDIT_TYPE_REF_GAP,
	EDIT_TYPE_MM,
	EDIT_TYPE_SNP
};
struct Edit
{
	Edit() { reset(); }
	Edit(
		uint32_t po,
		int ch,
		int qc,
		int ty,
		bool chrs = true)
	{
		init(po, ch, qc, ty, chrs);
	}
	void reset()
	{
		pos = pos2 = std::numeric_limits<uint32_t>::max();
		chr = qchr = type = 0;
	}
	bool inited() const
	{
		return pos != std::numeric_limits<uint32_t>::max();
	}
	void init(
		uint32_t po,
		int ch,
		int qc,
		int ty,
		bool chrs = true)
	{
		chr = ch;
		qchr = qc;
		type = ty;
		pos = po;
		if (qc == '-')
		{
			pos2 = std::numeric_limits<uint32_t>::max() >> 1;
		}
		else
		{
			pos2 = std::numeric_limits<uint32_t>::max();
		}
		if (!chrs)
		{
			assert_range(0, 4, (int)chr);
			assert_range(0, 4, (int)qchr);
			chr = "ACGTN"[chr];
			qchr = "ACGTN"[qchr];
		}
		assert_in(chr, "ACMGRSVTWYHKDBN-");
		assert_in(qchr, "ACGTN-");
		assert(chr != qchr || chr == 'N');
		assert(inited());
	}
	bool hasN() const
	{
		assert(inited());
		return chr == 'N' || qchr == 'N';
	}
	int operator<(const Edit &rhs) const
	{
		assert(inited());
		if (pos < rhs.pos)
			return 1;
		if (pos > rhs.pos)
			return 0;
		if (pos2 < rhs.pos2)
			return 1;
		if (pos2 > rhs.pos2)
			return 0;
		if (type < rhs.type)
			return 1;
		if (type > rhs.type)
			return 0;
		if (chr < rhs.chr)
			return 1;
		if (chr > rhs.chr)
			return 0;
		return (qchr < rhs.qchr) ? 1 : 0;
	}
	int operator==(const Edit &rhs) const
	{
		assert(inited());
		return (pos == rhs.pos &&
				pos2 == rhs.pos2 &&
				chr == rhs.chr &&
				qchr == rhs.qchr &&
				type == rhs.type);
	}
	bool isReadGap() const
	{
		assert(inited());
		return type == EDIT_TYPE_READ_GAP;
	}
	bool isRefGap() const
	{
		assert(inited());
		return type == EDIT_TYPE_REF_GAP;
	}
	bool isGap() const
	{
		assert(inited());
		return (type == EDIT_TYPE_REF_GAP || type == EDIT_TYPE_READ_GAP);
	}
	static size_t numGaps(const EList<Edit> &es)
	{
		size_t gaps = 0;
		for (size_t i = 0; i < es.size(); i++)
		{
			if (es[i].isGap())
				gaps++;
		}
		return gaps;
	}
	bool isMismatch() const
	{
		assert(inited());
		return type == EDIT_TYPE_MM;
	}
	static void sort(EList<Edit> &edits);
	static void invertPoss(
		EList<Edit> &edits,
		size_t sz,
		size_t ei,
		size_t en,
		bool sort = false);
	static void invertPoss(EList<Edit> &edits, size_t sz, bool sort = false)
	{
		invertPoss(edits, sz, 0, edits.size(), sort);
	}
	static void clipLo(EList<Edit> &edits, size_t len, size_t amt);
	static void clipHi(EList<Edit> &edits, size_t len, size_t amt);
	static void toRef(
		const BTDnaString &read,
		const EList<Edit> &edits,
		BTDnaString &ref,
		bool fw = true,
		size_t trim5 = 0,
		size_t trim3 = 0);
	static void printQAlign(
		std::ostream &os,
		const BTDnaString &read,
		const EList<Edit> &edits);
	static void printQAlign(
		std::ostream &os,
		const char *prefix,
		const BTDnaString &read,
		const EList<Edit> &edits);
	static void printQAlignNoCheck(
		std::ostream &os,
		const BTDnaString &read,
		const EList<Edit> &edits);
	static void printQAlignNoCheck(
		std::ostream &os,
		const char *prefix,
		const BTDnaString &read,
		const EList<Edit> &edits);
#ifndef NDEBUG
	bool repOk() const;
	static bool repOk(
		const EList<Edit> &edits,
		const BTDnaString &s,
		bool fw = true,
		size_t trim5 = 0,
		size_t trim3 = 0);
#endif
	uint8_t chr;
	uint8_t qchr;
	uint8_t type;
	uint32_t pos;
	uint32_t pos2;
	friend std::ostream &operator<<(std::ostream &os, const Edit &e);
	static void print(
		std::ostream &os,
		const EList<Edit> &edits,
		char delim = '\t');
	static void merge(EList<Edit> &dst, const EList<Edit> &src);
};
#endif
