#ifndef REFERENCE_H_
#define REFERENCE_H_
#include <stdexcept>
#include <fcntl.h>
#include <sys/stat.h>
#include <utility>
#ifdef _MM
#include <sys/mman.h>
#include <sys/shm.h>
#endif
#include "endian_swap.h"
#include "ref_read.h"
#include "sequence_io.h"
#include "mm.h"
#include "shmem.h"
#include "timer.h"
#include "sstring.h"
#include "btypes.h"
class BitPairReference
{
public:
	BitPairReference(
		const string &in,
		bool color,
		bool sanity = false,
		EList<string> *infiles = NULL,
		EList<SString<char>> *origs = NULL,
		bool infilesSeq = false,
		bool useMm = false,
		bool useShmem = false,
		bool mmSweep = false,
		bool verbose = false,
		bool startVerbose = false);
	~BitPairReference();
	int getBase(size_t tidx, size_t toff) const;
	int getStretchNaive(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count) const;
	int getStretch(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count
			ASSERT_ONLY(, SStringExpandable<uint32_t> &destU32_2)) const;
	TIndexOffU numRefs() const
	{
		return nrefs_;
	}
	TIndexOffU approxLen(TIndexOffU elt) const
	{
		assert_lt(elt, nrefs_);
		return refLens_[elt];
	}
	bool loaded() const
	{
		return loaded_;
	}
	TIndexOffU pastedOffset(TIndexOffU idx) const
	{
		return refOffs_[idx];
	}
	static std::pair<size_t, size_t>
	szsFromFasta(
		EList<FileBuf *> &is,
		const string &outfile,
		bool bigEndian,
		const RefReadInParams &refparams,
		EList<RefRecord> &szs,
		bool sanity);
protected:
	uint32_t byteToU32_[256];
	EList<RefRecord> recs_;
	EList<TIndexOffU> cumUnambig_;
	EList<TIndexOffU> cumRefOff_;
	EList<TIndexOffU> refLens_;
	EList<TIndexOffU> refOffs_;
	EList<TIndexOffU> refRecOffs_;
	uint8_t *buf_;
	uint8_t *sanityBuf_;
	TIndexOffU bufSz_;
	TIndexOffU bufAllocSz_;
	TIndexOffU nrefs_;
	bool loaded_;
	bool sanity_;
	bool useMm_;
	bool useShmem_;
	bool verbose_;
	ASSERT_ONLY(SStringExpandable<uint32_t> tmp_destU32_);
};
#endif
