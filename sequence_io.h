#ifndef SEQUENCE_IO_H_
#define SEQUENCE_IO_H_
#include <string>
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include "assert_helpers.h"
#include "ds.h"
#include "filebuf.h"
#include "sstring.h"
using namespace std;
template <typename TFnStr>
static void parseFastaLens(
	const TFnStr &infile,
	EList<size_t> &namelens, EList<size_t> &seqlens)
{
	FILE *in = fopen(sstr_to_cstr(infile), "r");
	if (in == NULL)
	{
		cerr << "Could not open sequence file" << endl;
		throw 1;
	}
	FileBuf fb(in);
	while (!fb.eof())
	{
		namelens.expand();
		namelens.back() = 0;
		seqlens.expand();
		seqlens.back() = 0;
		fb.parseFastaRecordLength(namelens.back(), seqlens.back());
		if (seqlens.back() == 0)
		{
			namelens.pop_back();
			seqlens.pop_back();
			continue;
		}
	}
	fb.close();
}
template <typename TFnStr, typename TNameStr, typename TSeqStr>
static void parseFasta(
	const TFnStr &infile,
	EList<TNameStr> &names, EList<size_t> &namelens, EList<TSeqStr> &seqs, EList<size_t> &seqlens)
{
	assert_eq(namelens.size(), seqlens.size());
	assert_eq(names.size(), namelens.size());
	assert_eq(seqs.size(), seqlens.size());
	size_t cur = namelens.size();
	parseFastaLens(infile, namelens, seqlens);
	FILE *in = fopen(sstr_to_cstr(infile), "r");
	if (in == NULL)
	{
		cerr << "Could not open sequence file" << endl;
		throw 1;
	}
	FileBuf fb(in);
	while (!fb.eof())
	{
		names.expand();
		seqs.expand();
		names.back() = new char[namelens[cur] + 1];
		seqs.back() = new char[seqlens[cur] + 1];
		fb.parseFastaRecord(names.back(), seqs.back());
		if (seqs.back().empty())
		{
			names.pop_back();
			seqs.pop_back();
			continue;
		}
	}
	fb.close();
}
template <typename TFnStr, typename TNameStr, typename TSeqStr>
static void parseFastas(
	const EList<TFnStr> &infiles,
	EList<TNameStr> &names, EList<size_t> &namelens, EList<TSeqStr> &seqs, EList<size_t> &seqlens)
{
	for (size_t i = 0; i < infiles.size(); i++)
	{
		parseFasta<TFnStr, TNameStr, TSeqStr>(
			infiles[i],
			names,
			namelens,
			seqs,
			seqlens);
	}
}
#endif
