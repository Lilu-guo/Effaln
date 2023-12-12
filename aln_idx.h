#ifndef EBWT_H_
#define EBWT_H_
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <fcntl.h>
#include <math.h>
#include <errno.h>
#include <stdexcept>
#include <sys/stat.h>
#ifdef _MM
#include <sys/mman.h>
#include <sys/shm.h>
#endif
#include "BitMap.h"
#include "InArray.h"
#include "loadkit.h"
#include "savekit.h"
#include "shmem.h"
#include "alphabet.h"
#include "assert_helpers.h"
#include "bitpack.h"
#include "blockwise_sa.h"
#include "endian_swap.h"
#include "word_io.h"
#include "random_source.h"
#include "ref_read.h"
#include "threading.h"
#include "str_util.h"
#include "mm.h"
#include "timer.h"
#include "reference.h"
#include "search_globals.h"
#include "ds.h"
#include "random_source.h"
#include "mem_ids.h"
#include "btypes.h"
#ifdef POPCNT_CAPABILITY
#include "processor_support.h"
#endif
#if __cplusplus <= 199711L
#define unique_ptr auto_ptr
#endif
using namespace std;
extern unsigned long long int countLF_Static;
extern unsigned long long int exact_Nums;
extern unsigned long long int onemm_Nums;
extern unsigned long long int allseed_Nums;
extern uint8_t cCntLUT_4[4][4][256];
static const uint64_t c_table[4] = {
	0xffffffffffffffff,
	0xaaaaaaaaaaaaaaaa,
	0x5555555555555555,
	0x0000000000000000};
#ifndef VMSG_NL
#define VMSG_NL(...)                \
	if (this->verbose())            \
	{                               \
		stringstream tmp;           \
		tmp << __VA_ARGS__ << endl; \
		this->verbose(tmp.str());   \
	}
#endif
#ifndef VMSG
#define VMSG(...)                 \
	if (this->verbose())          \
	{                             \
		stringstream tmp;         \
		tmp << __VA_ARGS__;       \
		this->verbose(tmp.str()); \
	}
#endif
enum EBWT_FLAGS
{
	EBWT_COLOR = 2,
	EBWT_ENTIRE_REV = 4
};
class EbwtParams
{
public:
	EbwtParams() {}
	EbwtParams(
		TIndexOffU len,
		int32_t lineRate,
		int32_t offRate,
		int32_t ftabChars,
		bool color,
		bool entireReverse)
	{
		init(len, lineRate, offRate, ftabChars, color, entireReverse);
	}
	EbwtParams(const EbwtParams &eh)
	{
		init(eh._len, eh._lineRate, eh._offRate,
			 eh._ftabChars, eh._color, eh._entireReverse);
	}
	void init(
		TIndexOffU len,
		int32_t lineRate,
		int32_t offRate,
		int32_t ftabChars,
		bool color,
		bool entireReverse)
	{
		_color = color;
		_entireReverse = entireReverse;
		_len = len;
		_bwtLen = _len + 1;
		_sz = (len + 3) / 4;
		_bwtSz = (len / 4 + 1);
		_lineRate = lineRate;
		_origOffRate = offRate;
		_offRate = offRate;
		_offMask = OFF_MASK << _offRate;
		_ftabChars = ftabChars;
		_eftabLen = _ftabChars * 2;
		_eftabSz = _eftabLen * OFF_SIZE;
		_ftabLen = (1 << (_ftabChars * 2)) + 1;
		_ftabSz = _ftabLen * OFF_SIZE;
		_offsLen = (_bwtLen + (1 << _offRate) - 1) >> _offRate;
		_offsSz = (uint64_t)_offsLen * OFF_SIZE;
		_lineSz = 1 << _lineRate;
		_sideSz = _lineSz * 1;
		_sideBwtSz = _sideSz - OFF_SIZE * 4;
		_sideBwtLen = _sideBwtSz * 4;
		_numSides = (_bwtSz + (_sideBwtSz)-1) / (_sideBwtSz);
		_numLines = _numSides * 1;
		_ebwtTotLen = _numSides * _sideSz;
		_ebwtTotSz = _ebwtTotLen;
		assert(repOk());
	}
	TIndexOffU len() const { return _len; }
	TIndexOffU lenNucs() const { return _len + (_color ? 1 : 0); }
	TIndexOffU bwtLen() const { return _bwtLen; }
	TIndexOffU sz() const { return _sz; }
	TIndexOffU bwtSz() const { return _bwtSz; }
	int32_t lineRate() const { return _lineRate; }
	int32_t origOffRate() const { return _origOffRate; }
	int32_t offRate() const { return _offRate; }
	TIndexOffU offMask() const { return _offMask; }
	int32_t ftabChars() const { return _ftabChars; }
	int32_t eftabLen() const { return _eftabLen; }
	int32_t eftabSz() const { return _eftabSz; }
	TIndexOffU ftabLen() const { return _ftabLen; }
	TIndexOffU ftabSz() const { return _ftabSz; }
	TIndexOffU offsLen() const { return _offsLen; }
	uint64_t offsSz() const { return _offsSz; }
	int32_t lineSz() const { return _lineSz; }
	int32_t sideSz() const { return _sideSz; }
	int32_t sideBwtSz() const { return _sideBwtSz; }
	int32_t sideBwtLen() const { return _sideBwtLen; }
	TIndexOffU numSides() const { return _numSides; }
	TIndexOffU numLines() const { return _numLines; }
	TIndexOffU ebwtTotLen() const { return _ebwtTotLen; }
	TIndexOffU ebwtTotSz() const { return _ebwtTotSz; }
	bool color() const { return _color; }
	bool entireReverse() const { return _entireReverse; }
	void setOffRate(int __offRate)
	{
		_offRate = __offRate;
		_offMask = OFF_MASK << _offRate;
		_offsLen = (_bwtLen + (1 << _offRate) - 1) >> _offRate;
		_offsSz = (uint64_t)_offsLen * OFF_SIZE;
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert_gt(_len, 0);
		assert_gt(_lineRate, 3);
		assert_geq(_offRate, 0);
		assert_leq(_ftabChars, 16);
		assert_geq(_ftabChars, 1);
		assert_lt(_lineRate, 32);
		assert_lt(_ftabChars, 32);
		assert_eq(0, _ebwtTotSz % _lineSz);
		return true;
	}
#endif
	void print(ostream &out) const
	{
		out << "Headers:" << endl
			<< "    len: " << _len << endl
			<< "    bwtLen: " << _bwtLen << endl
			<< "    sz: " << _sz << endl
			<< "    bwtSz: " << _bwtSz << endl
			<< "    lineRate: " << _lineRate << endl
			<< "    offRate: " << _offRate << endl
			<< "    offMask: 0x" << hex << _offMask << dec << endl
			<< "    ftabChars: " << _ftabChars << endl
			<< "    eftabLen: " << _eftabLen << endl
			<< "    eftabSz: " << _eftabSz << endl
			<< "    ftabLen: " << _ftabLen << endl
			<< "    ftabSz: " << _ftabSz << endl
			<< "    offsLen: " << _offsLen << endl
			<< "    offsSz: " << _offsSz << endl
			<< "    lineSz: " << _lineSz << endl
			<< "    sideSz: " << _sideSz << endl
			<< "    sideBwtSz: " << _sideBwtSz << endl
			<< "    sideBwtLen: " << _sideBwtLen << endl
			<< "    numSides: " << _numSides << endl
			<< "    numLines: " << _numLines << endl
			<< "    ebwtTotLen: " << _ebwtTotLen << endl
			<< "    ebwtTotSz: " << _ebwtTotSz << endl
			<< "    color: " << _color << endl
			<< "    reverse: " << _entireReverse << endl;
	}
	TIndexOffU _len;
	TIndexOffU _bwtLen;
	TIndexOffU _sz;
	TIndexOffU _bwtSz;
	int32_t _lineRate;
	int32_t _origOffRate;
	int32_t _offRate;
	TIndexOffU _offMask;
	int32_t _ftabChars;
	uint32_t _eftabLen;
	uint32_t _eftabSz;
	TIndexOffU _ftabLen;
	TIndexOffU _ftabSz;
	TIndexOffU _offsLen;
	uint64_t _offsSz;
	uint32_t _lineSz;
	uint32_t _sideSz;
	uint32_t _sideBwtSz;
	uint32_t _sideBwtLen;
	TIndexOffU _numSides;
	TIndexOffU _numLines;
	TIndexOffU _ebwtTotLen;
	TIndexOffU _ebwtTotSz;
	bool _color;
	bool _entireReverse;
};
class EbwtFileOpenException : public std::runtime_error
{
public:
	EbwtFileOpenException(const std::string &msg = "") : std::runtime_error(msg) {}
};
static inline int64_t fileSize(const char *name)
{
	std::ifstream f;
	f.open(name, std::ios_base::binary | std::ios_base::in);
	if (!f.good() || f.eof() || !f.is_open())
	{
		return 0;
	}
	f.seekg(0, std::ios_base::beg);
	std::ifstream::pos_type begin_pos = f.tellg();
	f.seekg(0, std::ios_base::end);
	return static_cast<int64_t>(f.tellg() - begin_pos);
}
struct SideLocus
{
	SideLocus() : _sideByteOff(0),
				  _sideNum(0),
				  _charOff(0),
				  _by(-1),
				  _bp(-1) {}
	SideLocus(TIndexOffU row, const EbwtParams &ep, const uint8_t *ebwt)
	{
		initFromRow(row, ep, ebwt);
	}
	static void initFromTopBot(
		TIndexOffU top,
		TIndexOffU bot,
		const EbwtParams &ep,
		const uint8_t *ebwt,
		SideLocus &ltop,
		SideLocus &lbot)
	{
		const TIndexOffU sideBwtLen = ep._sideBwtLen;
		assert_gt(bot, top);
		ltop.initFromRow(top, ep, ebwt);
		TIndexOffU spread = bot - top;
		if (ltop._charOff + spread < sideBwtLen)
		{
			lbot._charOff = (uint32_t)(ltop._charOff + spread);
			lbot._sideNum = ltop._sideNum;
			lbot._sideByteOff = ltop._sideByteOff;
			lbot._by = lbot._charOff >> 2;
			assert_lt(lbot._by, (int)ep._sideBwtSz);
			lbot._bp = lbot._charOff & 3;
		}
		else
		{
			lbot.initFromRow(bot, ep, ebwt);
		}
	}
	void initFromRow(TIndexOffU row, const EbwtParams &ep, const uint8_t *ebwt)
	{
		const int32_t sideSz = ep._sideSz;
		_sideNum = row / (112 * OFF_SIZE);
		assert_lt(_sideNum, ep._numSides);
		_charOff = row % (112 * OFF_SIZE);
		_sideByteOff = _sideNum * sideSz;
		assert_leq(row, ep._len);
		assert_leq(_sideByteOff + sideSz, ep._ebwtTotSz);
		_by = _charOff >> 2;
		assert_lt(_by, (int)ep._sideBwtSz);
		_bp = _charOff & 3;
	}
	void nextSide(const EbwtParams &ep)
	{
		assert(valid());
		_sideByteOff += ep.sideSz();
		_sideNum++;
		_by = _bp = _charOff = 0;
		assert(valid());
	}
	bool valid() const
	{
		if (_bp != -1)
		{
			return true;
		}
		return false;
	}
	TIndexOffU toBWRow() const
	{
		return _sideNum * 112 * OFF_SIZE + _charOff;
	}
#ifndef NDEBUG
	bool repOk(const EbwtParams &ep) const
	{
		ASSERT_ONLY(TIndexOffU row = _sideNum * 112 * OFF_SIZE + _charOff);
		assert_leq(row, ep._len);
		assert_range(-1, 3, _bp);
		assert_range(0, (int)ep._sideBwtSz, _by);
		return true;
	}
#endif
	void invalidate()
	{
		_bp = -1;
	}
	const uint8_t *side(const uint8_t *ebwt) const
	{
		return ebwt + _sideByteOff;
	}
	TIndexOffU _sideByteOff;
	TIndexOffU _sideNum;
	uint32_t _charOff;
	int32_t _by;
	int32_t _bp;
};
#ifdef POPCNT_CAPABILITY
struct USE_POPCNT_GENERIC
{
#endif
	inline static int pop64(uint64_t x)
	{
		x = x - ((x >> 1) & 0x5555555555555555llu);
		x = (x & 0x3333333333333333llu) + ((x >> 2) & 0x3333333333333333llu);
		x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Fllu;
		x = x + (x >> 8);
		x = x + (x >> 16);
		x = x + (x >> 32);
		return (int)(x & 0x3Fllu);
	}
#ifdef POPCNT_CAPABILITY
};
#endif
#ifdef POPCNT_CAPABILITY
struct USE_POPCNT_INSTRUCTION
{
	inline static int pop64(uint64_t x)
	{
		return __builtin_popcountll(x);
	}
};
#endif
#ifdef POPCNT_CAPABILITY
template <typename Operation>
#endif
inline static int countInU64(int c, uint64_t dw)
{
#ifdef POPCNT_CAPABILITY
	uint64_t tmp = Operation().pop64(dw);
#else
	uint64_t tmp = pop64(dw);
#endif
	return (int)tmp;
}
class EbwtSearchParams;
class Ebwt
{
public:
#define Ebwt_INITS                         \
	_toBigEndian(currentlyBigEndian()),    \
		_overrideOffRate(overrideOffRate), \
		_verbose(verbose),                 \
		_passMemExc(passMemExc),           \
		_sanity(sanityCheck),              \
		fw_(fw),                           \
		_in1(NULL),                        \
		_in2(NULL),                        \
		_zOff(OFF_MASK),                   \
		_zEbwtByteOff(OFF_MASK),           \
		_zEbwtBpOff(-1),                   \
		_nPat(0),                          \
		_nFrag(0),                         \
		_plen(EBWT_CAT),                   \
		_rstarts(EBWT_CAT),                \
		_fchr(EBWT_CAT),                   \
		_ftab(EBWT_CAT),                   \
		_eftab(EBWT_CAT),                  \
		_offs(EBWT_CAT),                   \
		_ebwt(EBWT_CAT),                   \
		_useMm(false),                     \
		useShmem_(false),                  \
		_refnames(EBWT_CAT),               \
		mmFile1_(NULL),                    \
		mmFile2_(NULL)
	Ebwt(const string &in,
		 int color,
		 int needEntireReverse,
		 bool fw,
		 int32_t overrideOffRate,
		 int32_t offRatePlus,
		 bool useMm,
		 bool useShmem,
		 bool mmSweep,
		 bool loadNames,
		 bool loadSASamp,
		 bool loadFtab,
		 bool loadRstarts,
		 bool verbose,
		 bool startVerbose,
		 bool passMemExc,
		 bool sanityCheck) : Ebwt_INITS
	{
		assert(!useMm || !useShmem);
#ifdef POPCNT_CAPABILITY
		ProcessorSupport ps;
		_usePOPCNTinstruction = ps.POPCNTenabled();
#endif
		packed_ = false;
		_useMm = useMm;
		useShmem_ = useShmem;
		_in1Str = in + ".1." + gEbwt_ext;
		_in2Str = in + ".2." + gEbwt_ext;
		readIntoMemory(
			color,
			fw ? -1 : needEntireReverse,
			loadSASamp,
			loadFtab,
			loadRstarts,
			true,
			&_eh,
			mmSweep,
			loadNames,
			startVerbose);
		if (offRatePlus > 0 && _overrideOffRate == -1)
		{
			_overrideOffRate = _eh._offRate + offRatePlus;
		}
		if (_overrideOffRate > _eh._offRate)
		{
			_eh.setOffRate(_overrideOffRate);
			assert_eq(_overrideOffRate, _eh._offRate);
		}
		assert(repOk());
	}
	template <typename TStr>
	Ebwt(
		TStr exampleStr,
		bool packed,
		int color,
		int needEntireReverse,
		int32_t lineRate,
		int32_t offRate,
		int32_t ftabChars,
		int nthreads,
		const string &file,
		bool fw,
		bool useBlockwise,
		TIndexOffU bmax,
		TIndexOffU bmaxSqrtMult,
		TIndexOffU bmaxDivN,
		int dcv,
		EList<FileBuf *> &is,
		EList<RefRecord> &szs,
		TIndexOffU sztot,
		const RefReadInParams &refparams,
		uint32_t seed,
		int32_t overrideOffRate = -1,
		bool doSaFile = false,
		bool doBwtFile = false,
		bool verbose = false,
		bool passMemExc = false,
		bool sanityCheck = false) : Ebwt_INITS,
									_eh(
										joinedLen(szs),
										lineRate,
										offRate,
										ftabChars,
										color,
										refparams.reverse == REF_READ_REVERSE)
	{
#ifdef POPCNT_CAPABILITY
		ProcessorSupport ps;
		_usePOPCNTinstruction = ps.POPCNTenabled();
#endif
		_in1Str = file + ".1." + gEbwt_ext;
		_in2Str = file + ".2." + gEbwt_ext;
		_in6Str = file + ".6." + gEbwt_ext;
		packed_ = packed;
		ofstream fout1(_in1Str.c_str(), ios::binary);
		if (!fout1.good())
		{
			cerr << "Could not open index file for writing: \"" << _in1Str.c_str() << "\"" << endl;
			throw 1;
		}
		ofstream fout2(_in2Str.c_str(), ios::binary);
		if (!fout2.good())
		{
			cerr << "Could not open index file for writing: \"" << _in2Str.c_str() << "\"" << endl;
			throw 1;
		}
		_inSaStr = file + ".sa";
		_inBwtStr = file + ".bwt";
		ofstream *saOut = NULL, *bwtOut = NULL;
		if (doSaFile)
		{
			saOut = new ofstream(_inSaStr.c_str(), ios::binary);
			if (!saOut->good())
			{
				cerr << "Could not open suffix-array file for writing: \"" << _inSaStr.c_str() << "\"" << endl;
				throw 1;
			}
		}
		if (doBwtFile)
		{
			bwtOut = new ofstream(_inBwtStr.c_str(), ios::binary);
			if (!bwtOut->good())
			{
				cerr << "Could not open suffix-array file for writing: \"" << _inBwtStr.c_str() << "\"" << endl;
				throw 1;
			}
		}
		initFromVector<TStr>(
			is,
			szs,
			sztot,
			refparams,
			fout1,
			fout2,
			_in6Str,
			file,
			saOut,
			bwtOut,
			nthreads,
			useBlockwise,
			bmax,
			bmaxSqrtMult,
			bmaxDivN,
			dcv,
			seed,
			verbose);
		fout1.flush();
		int64_t tellpSz1 = (int64_t)fout1.tellp();
		VMSG_NL("Wrote " << fout1.tellp() << " bytes to primary EBWT file: " << _in1Str.c_str());
		fout1.close();
		bool err = false;
		if (tellpSz1 > fileSize(_in1Str.c_str()))
		{
			err = true;
			cerr << "Index is corrupt: File size for " << _in1Str.c_str() << " should have been " << tellpSz1
				 << " but is actually " << fileSize(_in1Str.c_str()) << "." << endl;
		}
		fout2.flush();
		int64_t tellpSz2 = (int64_t)fout2.tellp();
		VMSG_NL("Wrote " << fout2.tellp() << " bytes to secondary EBWT file: " << _in2Str.c_str());
		fout2.close();
		if (tellpSz2 > fileSize(_in2Str.c_str()))
		{
			err = true;
			cerr << "Index is corrupt: File size for " << _in2Str.c_str() << " should have been " << tellpSz2
				 << " but is actually " << fileSize(_in2Str.c_str()) << "." << endl;
		}
		if (saOut != NULL)
		{
			int64_t tellpSzSa = (int64_t)saOut->tellp();
			VMSG_NL("Wrote " << tellpSzSa << " bytes to suffix-array file: " << _inSaStr.c_str());
			saOut->close();
			delete saOut;
			if (tellpSzSa > fileSize(_inSaStr.c_str()))
			{
				err = true;
				cerr << "Index is corrupt: File size for " << _inSaStr.c_str() << " should have been " << tellpSzSa
					 << " but is actually " << fileSize(_inSaStr.c_str()) << "." << endl;
			}
		}
		if (bwtOut != NULL)
		{
			int64_t tellpSzBwt = (int64_t)bwtOut->tellp();
			VMSG_NL("Wrote " << tellpSzBwt << " bytes to BWT file: " << _inBwtStr.c_str());
			bwtOut->close();
			delete bwtOut;
			if (tellpSzBwt > fileSize(_inBwtStr.c_str()))
			{
				err = true;
				cerr << "Index is corrupt: File size for " << _inBwtStr.c_str() << " should have been " << tellpSzBwt
					 << " but is actually " << fileSize(_inBwtStr.c_str()) << "." << endl;
			}
		}
		if (err)
		{
			cerr << "Please check if there is a problem with the disk or if disk is full." << endl;
			throw 1;
		}
		VMSG_NL("Re-opening _in1 and _in2 as input streams");
		if (_sanity)
		{
			VMSG_NL("Sanity-checking Bt2");
			assert(!isInMemory());
			readIntoMemory(
				color,
				fw ? -1 : needEntireReverse,
				true,
				true,
				true,
				false,
				NULL,
				false,
				true,
				false);
			sanityCheckAll(refparams.reverse);
			evictFromMemory();
			assert(!isInMemory());
		}
		VMSG_NL("Returning from Ebwt constructor");
	}
	template <typename TStr>
	static pair<Ebwt *, Ebwt *>
	fromString(
		const char *str,
		bool packed,
		int color,
		int reverse,
		bool bigEndian,
		int32_t lineRate,
		int32_t offRate,
		int32_t ftabChars,
		const string &file,
		bool useBlockwise,
		TIndexOffU bmax,
		TIndexOffU bmaxSqrtMult,
		TIndexOffU bmaxDivN,
		int dcv,
		uint32_t seed,
		bool verbose,
		bool autoMem,
		bool sanity)
	{
		EList<std::string> strs(EBWT_CAT);
		strs.push_back(std::string(str));
		return fromStrings<TStr>(
			strs,
			packed,
			color,
			reverse,
			bigEndian,
			lineRate,
			offRate,
			ftabChars,
			file,
			useBlockwise,
			bmax,
			bmaxSqrtMult,
			bmaxDivN,
			dcv,
			seed,
			verbose,
			autoMem,
			sanity);
	}
	template <typename TStr>
	static pair<Ebwt *, Ebwt *>
	fromStrings(
		const EList<std::string> &strs,
		bool packed,
		int color,
		int reverse,
		bool bigEndian,
		int32_t lineRate,
		int32_t offRate,
		int32_t ftabChars,
		const string &file,
		bool useBlockwise,
		TIndexOffU bmax,
		TIndexOffU bmaxSqrtMult,
		TIndexOffU bmaxDivN,
		int dcv,
		uint32_t seed,
		bool verbose,
		bool autoMem,
		bool sanity)
	{
		assert(!strs.empty());
		EList<FileBuf *> is(EBWT_CAT);
		RefReadInParams refparams(color, REF_READ_FORWARD, false, false);
		unique_ptr<stringstream> ss(new stringstream());
		for (TIndexOffU i = 0; i < strs.size(); i++)
		{
			(*ss) << ">" << i << endl
				  << strs[i] << endl;
		}
		unique_ptr<FileBuf> fb(new FileBuf(ss.get()));
		assert(!fb->eof());
		assert(fb->get() == '>');
		ASSERT_ONLY(fb->reset());
		assert(!fb->eof());
		is.push_back(fb.get());
		EList<RefRecord> szs(EBWT_CAT);
		std::pair<TIndexOffU, TIndexOffU> sztot;
		sztot = BitPairReference::szsFromFasta(is, file, bigEndian, refparams, szs, sanity);
		Ebwt *ebwtFw = new Ebwt(
			TStr(),
			packed,
			refparams.color ? 1 : 0,
			-1,
			lineRate,
			offRate,
			ftabChars,
			file,
			true,
			useBlockwise,
			bmax,
			bmaxSqrtMult,
			bmaxDivN,
			dcv,
			is,
			szs,
			sztot.first,
			refparams,
			seed,
			-1,
			verbose,
			autoMem,
			sanity);
		refparams.reverse = reverse;
		szs.clear();
		sztot = BitPairReference::szsFromFasta(is, file, bigEndian, refparams, szs, sanity);
		Ebwt *ebwtBw = new Ebwt(
			TStr(),
			packed,
			refparams.color ? 1 : 0,
			reverse == REF_READ_REVERSE,
			lineRate,
			offRate,
			ftabChars,
			file + ".rev",
			false,
			useBlockwise,
			bmax,
			bmaxSqrtMult,
			bmaxDivN,
			dcv,
			is,
			szs,
			sztot.first,
			refparams,
			seed,
			-1,
			verbose,
			autoMem,
			sanity);
		return make_pair(ebwtFw, ebwtBw);
	}
	bool isPacked() { return packed_; }
	void szsToDisk(const EList<RefRecord> &szs, ostream &os, int reverse);
	template <typename TStr>
	void initFromVector(EList<FileBuf *> &is,
						EList<RefRecord> &szs,
						TIndexOffU sztot,
						const RefReadInParams &refparams,
						ofstream &out1,
						ofstream &out2,
						string &_in6Str,
						const string &outfile,
						ofstream *saOut,
						ofstream *bwtOut,
						int nthreads,
						bool useBlockwise,
						TIndexOffU bmax,
						TIndexOffU bmaxSqrtMult,
						TIndexOffU bmaxDivN,
						int dcv,
						uint32_t seed,
						bool verbose)
	{
		VMSG_NL("Calculating joined length");
		TStr s;
		TIndexOffU jlen;
		jlen = joinedLen(szs);
		assert_geq(jlen, sztot);
		VMSG_NL("Writing header");
		writeFromMemory(true, out1, out2);
		try
		{
			VMSG_NL("Reserving space for joined string");
			s.resize(jlen);
			VMSG_NL("Joining reference sequences");
			if (refparams.reverse == REF_READ_REVERSE)
			{
				{
					Timer timer(cout, "  Time to join reference sequences: ", _verbose);
					joinToDisk(is, szs, sztot, refparams, s, out1, out2);
				}
				{
					Timer timer(cout, "  Time to reverse reference sequence: ", _verbose);
					EList<RefRecord> tmp(EBWT_CAT);
					s.reverse();
					reverseRefRecords(szs, tmp, false, verbose);
					szsToDisk(tmp, out1, refparams.reverse);
				}
			}
			else
			{
				Timer timer(cout, "  Time to join reference sequences: ", _verbose);
				joinToDisk(is, szs, sztot, refparams, s, out1, out2);
				szsToDisk(szs, out1, refparams.reverse);
			}
		}
		catch (bad_alloc &e)
		{
			cerr << "Could not allocate space for a joined string of " << jlen << " elements." << endl;
			if (!isPacked() && _passMemExc)
			{
				throw e;
			}
			if (isPacked())
			{
				cerr << "Please try running effaln-index on a computer with more memory." << endl;
			}
			else
			{
				cerr << "Please try running effaln-index in packed mode (-p/--packed) or in automatic" << endl
					 << "mode (-a/--auto), or try again on a computer with more memory." << endl;
			}
			if (sizeof(void *) == 4)
			{
				cerr << "If this computer has more than 4 GB of memory, try using a 64-bit executable;" << endl
					 << "this executable is 32-bit." << endl;
			}
			throw 1;
		}
		assert_geq(s.length(), jlen);
		if (bmax != OFF_MASK)
		{
			VMSG_NL("bmax according to bmax setting: " << bmax);
		}
		else if (bmaxSqrtMult != OFF_MASK)
		{
			bmax *= bmaxSqrtMult;
			VMSG_NL("bmax according to bmaxSqrtMult setting: " << bmax);
		}
		else if (bmaxDivN != OFF_MASK)
		{
			bmax = max<TIndexOffU>(jlen / bmaxDivN, 1);
			VMSG_NL("bmax according to bmaxDivN setting: " << bmax);
		}
		else
		{
			bmax = (TIndexOffU)sqrt(s.length());
			VMSG_NL("bmax defaulted to: " << bmax);
		}
		int iter = 0;
		bool first = true;
		streampos out1pos = out1.tellp();
		streampos out2pos = out2.tellp();
		while (true)
		{
			if (!first && bmax < 40 && _passMemExc)
			{
				cerr << "Could not find approrpiate bmax/dcv settings for building this index." << endl;
				if (!isPacked())
				{
					throw bad_alloc();
				}
				else
				{
					cerr << "Already tried a packed string representation." << endl;
				}
				cerr << "Please try indexing this reference on a computer with more memory." << endl;
				if (sizeof(void *) == 4)
				{
					cerr << "If this computer has more than 4 GB of memory, try using a 64-bit executable;" << endl
						 << "this executable is 32-bit." << endl;
				}
				throw 1;
			}
			if (!first)
			{
				out1.seekp(out1pos);
				out2.seekp(out2pos);
			}
			if (dcv > 4096)
				dcv = 4096;
			if ((iter % 6) == 5 && dcv < 4096 && dcv != 0)
			{
				dcv <<= 1;
			}
			else
			{
				bmax -= (bmax >> 2);
			}
			VMSG("Using parameters --bmax " << bmax);
			if (dcv == 0)
			{
				VMSG_NL(" and *no difference cover*");
			}
			else
			{
				VMSG_NL(" --dcv " << dcv);
			}
			iter++;
			try
			{
				{
					VMSG_NL("  Doing ahead-of-time memory usage test");
					dcv <<= 1;
					TIndexOffU sz = (TIndexOffU)DifferenceCoverSample<TStr>::simulateAllocs(s, dcv >> 1);
					if (nthreads > 1)
						sz *= (nthreads + 1);
					AutoArray<uint8_t> tmp(sz, EBWT_CAT);
					dcv >>= 1;
					sz = (TIndexOffU)KarkkainenBlockwiseSA<TStr>::simulateAllocs(s, bmax);
					AutoArray<uint8_t> tmp2(sz, EBWT_CAT);
					AutoArray<TIndexOffU> ftab(_eh._ftabLen * 2, EBWT_CAT);
					AutoArray<uint8_t> side(_eh._sideSz, EBWT_CAT);
					AutoArray<uint32_t> extra(20 * 1024 * 1024, EBWT_CAT);
					VMSG("  Passed!  Constructing with these parameters: --bmax " << bmax << " --dcv " << dcv);
					if (isPacked())
					{
						VMSG(" --packed");
					}
					VMSG_NL("");
				}
				VMSG_NL("Constructing suffix-array element generator");
				KarkkainenBlockwiseSA<TStr> bsa(s, bmax, nthreads, dcv, seed, _sanity, _passMemExc, _verbose, outfile);
				assert(bsa.suffixItrIsReset());
				assert_eq(bsa.size(), s.length() + 1);
				VMSG_NL("Converting suffix-array elements to index image");
				buildToDisk(bsa, s, out1, out2, _in6Str, saOut, bwtOut);
				out1.flush();
				out2.flush();
				bool failed = out1.fail() || out2.fail();
				if (saOut != NULL)
				{
					saOut->flush();
					failed = failed || saOut->fail();
				}
				if (bwtOut != NULL)
				{
					bwtOut->flush();
					failed = failed || bwtOut->fail();
				}
				if (failed)
				{
					cerr << "An error occurred writing the index to disk.  Please check if the disk is full." << endl;
					throw 1;
				}
				break;
			}
			catch (bad_alloc &e)
			{
				if (_passMemExc)
				{
					VMSG_NL("  Ran out of memory; automatically trying more memory-economical parameters.");
				}
				else
				{
					cerr << "Out of memory while constructing suffix array.  Please try using a smaller" << endl
						 << "number of blocks by specifying a smaller --bmax or a larger --bmaxdivn" << endl;
					throw 1;
				}
			}
			first = false;
		}
		assert(repOk());
		assert_eq(this->_refnames.size(), this->_nPat);
		for (TIndexOffU i = 0; i < this->_refnames.size(); i++)
		{
			out1 << this->_refnames[i].c_str() << endl;
		}
		out1 << '\0';
		out1.flush();
		out2.flush();
		if (out1.fail() || out2.fail())
		{
			cerr << "An error occurred writing the index to disk.  Please check if the disk is full." << endl;
			throw 1;
		}
		VMSG_NL("Returning from initFromVector");
	}
	TIndexOffU joinedLen(EList<RefRecord> &szs)
	{
		TIndexOffU ret = 0;
		for (unsigned int i = 0; i < szs.size(); i++)
		{
			ret += szs[i].len;
		}
		return ret;
	}
	~Ebwt()
	{
		_fchr.reset();
		_ftab.reset();
		_eftab.reset();
		_plen.reset();
		_rstarts.reset();
		_offs.reset();
		_ebwt.reset();
		if (offs() != NULL && useShmem_)
		{
			FREE_SHARED(offs());
		}
		if (ebwt() != NULL && useShmem_)
		{
			FREE_SHARED(ebwt());
		}
		if (_in1 != NULL)
			fclose(_in1);
		if (_in2 != NULL)
			fclose(_in2);
	}
	inline const EbwtParams &eh() const { return _eh; }
	TIndexOffU zOff() const { return _zOff; }
	TIndexOffU zEbwtByteOff() const { return _zEbwtByteOff; }
	TIndexOff zEbwtBpOff() const { return _zEbwtBpOff; }
	TIndexOffU nPat() const { return _nPat; }
	TIndexOffU nFrag() const { return _nFrag; }
	inline TIndexOffU *fchr() { return _fchr.get(); }
	inline TIndexOffU *ftab() { return _ftab.get(); }
	inline TIndexOffU *eftab() { return _eftab.get(); }
	inline TIndexOffU *offs() { return _offs.get(); }
	inline TIndexOffU *plen() { return _plen.get(); }
	inline TIndexOffU *rstarts() { return _rstarts.get(); }
	inline uint8_t *ebwt() { return _ebwt.get(); }
	inline const TIndexOffU *fchr() const { return _fchr.get(); }
	inline const TIndexOffU *ftab() const { return _ftab.get(); }
	inline const TIndexOffU *eftab() const { return _eftab.get(); }
	inline const TIndexOffU *offs() const { return _offs.get(); }
	inline const TIndexOffU *plen() const { return _plen.get(); }
	inline const TIndexOffU *rstarts() const { return _rstarts.get(); }
	inline const uint8_t *ebwt() const { return _ebwt.get(); }
	bool toBe() const { return _toBigEndian; }
	bool verbose() const { return _verbose; }
	bool sanityCheck() const { return _sanity; }
	EList<string> &refnames() { return _refnames; }
	bool fw() const { return fw_; }
#ifdef POPCNT_CAPABILITY
	bool _usePOPCNTinstruction;
#endif
	bool contains(
		const BTDnaString &str,
		TIndexOffU *top = NULL,
		TIndexOffU *bot = NULL) const;
	bool contains(
		const char *str,
		TIndexOffU *top = NULL,
		TIndexOffU *bot = NULL) const
	{
		return contains(BTDnaString(str, true), top, bot);
	}
	bool isInMemory() const
	{
		if (ebwt() != NULL)
		{
			assert(_eh.repOk());
			assert(fchr() != NULL);
			assert_neq(_zEbwtByteOff, OFF_MASK);
			assert_neq(_zEbwtBpOff, -1);
			return true;
		}
		else
		{
			assert(ftab() == NULL);
			assert(eftab() == NULL);
			assert(fchr() == NULL);
			assert(offs() == NULL);
			assert(rstarts() == NULL);
			assert_eq(_zEbwtByteOff, OFF_MASK);
			assert_eq(_zEbwtBpOff, -1);
			return false;
		}
	}
	bool isEvicted() const
	{
		return !isInMemory();
	}
	void loadIntoMemory(
		int color,
		int needEntireReverse,
		bool loadSASamp,
		bool loadFtab,
		bool loadRstarts,
		bool loadNames,
		bool verbose)
	{
		readIntoMemory(
			color,
			needEntireReverse,
			loadSASamp,
			loadFtab,
			loadRstarts,
			false,
			NULL,
			false,
			loadNames,
			verbose);
	}
	void evictFromMemory()
	{
		assert(isInMemory());
		_fchr.free();
		_ftab.free();
		_eftab.free();
		_rstarts.free();
		_offs.free();
		_ebwt.free();
		_zEbwtByteOff = OFF_MASK;
		_zEbwtBpOff = -1;
	}
	TIndexOffU ftabSeqToInt(
		const BTDnaString &seq,
		size_t off,
		bool rev) const
	{
		int fc = _eh._ftabChars;
		size_t lo = off, hi = lo + fc;
		assert_leq(hi, seq.length());
		TIndexOffU ftabOff = 0;
		for (int i = 0; i < fc; i++)
		{
			bool fwex = fw();
			if (rev)
				fwex = !fwex;
			int c = (fwex ? seq[lo + i] : seq[hi - i - 1]);
			if (c > 3)
			{
				return std::numeric_limits<TIndexOffU>::max();
			}
			assert_range(0, 3, c);
			ftabOff <<= 2;
			ftabOff |= c;
		}
		return ftabOff;
	}
	TIndexOffU ftabHi(TIndexOffU i) const
	{
		return Ebwt::ftabHi(
			ftab(),
			eftab(),
			_eh._len,
			_eh._ftabLen,
			_eh._eftabLen,
			i);
	}
	static TIndexOffU ftabHi(
		const TIndexOffU *ftab,
		const TIndexOffU *eftab,
		TIndexOffU len,
		TIndexOffU ftabLen,
		TIndexOffU eftabLen,
		TIndexOffU i)
	{
		assert_lt(i, ftabLen);
		if (ftab[i] <= len)
		{
			return ftab[i];
		}
		else
		{
			TIndexOffU efIdx = ftab[i] ^ OFF_MASK;
			assert_lt(efIdx * 2 + 1, eftabLen);
			return eftab[efIdx * 2 + 1];
		}
	}
	TIndexOffU ftabLo(TIndexOffU i) const
	{
		return Ebwt::ftabLo(
			ftab(),
			eftab(),
			_eh._len,
			_eh._ftabLen,
			_eh._eftabLen,
			i);
	}
	TIndexOffU ftabLo(const BTDnaString &seq, size_t off) const
	{
		return ftabLo(ftabSeqToInt(seq, off, false));
	}
	TIndexOffU ftabHi(const BTDnaString &seq, size_t off) const
	{
		return ftabHi(ftabSeqToInt(seq, off, false));
	}
	bool
	ftabLoHi(
		const BTDnaString &seq,
		size_t off,
		bool rev,
		TIndexOffU &top,
		TIndexOffU &bot) const
	{
		TIndexOffU fi = ftabSeqToInt(seq, off, rev);
		if (fi == std::numeric_limits<TIndexOffU>::max())
		{
			return false;
		}
		top = ftabHi(fi);
		bot = ftabLo(fi + 1);
		assert_geq(bot, top);
		return true;
	}
	static TIndexOffU ftabLo(
		const TIndexOffU *ftab,
		const TIndexOffU *eftab,
		TIndexOffU len,
		TIndexOffU ftabLen,
		TIndexOffU eftabLen,
		TIndexOffU i)
	{
		assert_lt(i, ftabLen);
		if (ftab[i] <= len)
		{
			return ftab[i];
		}
		else
		{
			TIndexOffU efIdx = ftab[i] ^ OFF_MASK;
			assert_lt(efIdx * 2 + 1, eftabLen);
			return eftab[efIdx * 2];
		}
	}
	TIndexOffU tryOffset(TIndexOffU elt, size_t i) const
	{
		assert(offs() != NULL);
		if (elt == _zOff)
			return 0;
		if ((elt & _eh._offMask) == elt)
		{
			TIndexOffU eltOff = elt >> _eh._offRate;
			assert_lt(eltOff, _eh._offsLen);
			TIndexOffU off = offs()[eltOff];
			assert_neq(OFF_MASK, off);
			return off;
		}
		else
		{
			return OFF_MASK;
		}
	}
	TIndexOffU tryOffset(
		TIndexOffU elt,
		bool fw,
		TIndexOffU hitlen) const
	{
		TIndexOffU off = tryOffset(elt, 0);
		if (off != OFF_MASK && !fw)
		{
			assert_lt(off, _eh._len);
			off = _eh._len - off - 1;
			assert_geq(off, hitlen - 1);
			off -= (hitlen - 1);
			assert_lt(off, _eh._len);
		}
		return off;
	}
	TIndexOffU walkLeft(TIndexOffU row, TIndexOffU steps) const;
	TIndexOffU getOffset(TIndexOffU row) const;
	TIndexOffU getOffset(
		TIndexOffU elt,
		bool fw,
		TIndexOffU hitlen) const;
	void postReadInit(EbwtParams &eh)
	{
		TIndexOffU sideNum = _zOff / eh._sideBwtLen;
		TIndexOffU sideCharOff = _zOff % eh._sideBwtLen;
		TIndexOffU sideByteOff = sideNum * eh._sideSz;
		_zEbwtByteOff = sideCharOff >> 2;
		assert_lt(_zEbwtByteOff, eh._sideBwtSz);
		_zEbwtBpOff = sideCharOff & 3;
		assert_lt(_zEbwtBpOff, 4);
		_zEbwtByteOff += sideByteOff;
		assert(repOk(eh));
	}
	static int32_t readFlags(const string &instr);
	void print(ostream &out) const
	{
		print(out, _eh);
	}
	void print(ostream &out, const EbwtParams &eh) const
	{
		eh.print(out);
		out << "Ebwt (" << (isInMemory() ? "memory" : "disk") << "):" << endl
			<< "    zOff: " << _zOff << endl
			<< "    zEbwtByteOff: " << _zEbwtByteOff << endl
			<< "    zEbwtBpOff: " << _zEbwtBpOff << endl
			<< "    nPat: " << _nPat << endl
			<< "    plen: ";
		if (plen() == NULL)
		{
			out << "NULL" << endl;
		}
		else
		{
			out << "non-NULL, [0] = " << plen()[0] << endl;
		}
		out << "    rstarts: ";
		if (rstarts() == NULL)
		{
			out << "NULL" << endl;
		}
		else
		{
			out << "non-NULL, [0] = " << rstarts()[0] << endl;
		}
		out << "    ebwt: ";
		if (ebwt() == NULL)
		{
			out << "NULL" << endl;
		}
		else
		{
			out << "non-NULL, [0] = " << ebwt()[0] << endl;
		}
		out << "    fchr: ";
		if (fchr() == NULL)
		{
			out << "NULL" << endl;
		}
		else
		{
			out << "non-NULL, [0] = " << fchr()[0] << endl;
		}
		out << "    ftab: ";
		if (ftab() == NULL)
		{
			out << "NULL" << endl;
		}
		else
		{
			out << "non-NULL, [0] = " << ftab()[0] << endl;
		}
		out << "    eftab: ";
		if (eftab() == NULL)
		{
			out << "NULL" << endl;
		}
		else
		{
			out << "non-NULL, [0] = " << eftab()[0] << endl;
		}
		out << "    offs: ";
		if (offs() == NULL)
		{
			out << "NULL" << endl;
		}
		else
		{
			out << "non-NULL, [0] = " << offs()[0] << endl;
		}
	}
	template <typename TStr>
	static TStr join(EList<TStr> &l, uint32_t seed);
	template <typename TStr>
	static TStr join(EList<FileBuf *> &l, EList<RefRecord> &szs, TIndexOffU sztot, const RefReadInParams &refparams, uint32_t seed);
	template <typename TStr>
	void joinToDisk(EList<FileBuf *> &l, EList<RefRecord> &szs, TIndexOffU sztot, const RefReadInParams &refparams, TStr &ret, ostream &out1, ostream &out2);
	template <typename TStr>
	void buildToDisk(InorderBlockwiseSA<TStr> &sa, const TStr &s, ostream &out1, ostream &out2, string &_in6Str, ostream *saOut, ostream *bwtOut);
	void readIntoMemory(int color, int needEntireRev, bool loadSASamp, bool loadFtab, bool loadRstarts, bool justHeader, EbwtParams *params, bool mmSweep, bool loadNames, bool startVerbose);
	void writeFromMemory(bool justHeader, ostream &out1, ostream &out2) const;
	void writeFromMemory(bool justHeader, const string &out1, const string &out2) const;
	void sanityCheckUpToSide(TIndexOff upToSide) const;
	void sanityCheckAll(int reverse) const;
	void restore(SString<char> &s) const;
	void checkOrigs(const EList<SString<char>> &os, bool color, bool mirror) const;
	void joinedToTextOff(TIndexOffU qlen, TIndexOffU off, TIndexOffU &tidx, TIndexOffU &textoff, TIndexOffU &tlen, bool rejectStraddle, bool &straddled) const;
#define WITHIN_BWT_LEN(x)                    \
	assert_leq(x[0], this->_eh._sideBwtLen); \
	assert_leq(x[1], this->_eh._sideBwtLen); \
	assert_leq(x[2], this->_eh._sideBwtLen); \
	assert_leq(x[3], this->_eh._sideBwtLen)
#define WITHIN_FCHR(x)                 \
	assert_leq(x[0], this->fchr()[1]); \
	assert_leq(x[1], this->fchr()[2]); \
	assert_leq(x[2], this->fchr()[3]); \
	assert_leq(x[3], this->fchr()[4])
#define WITHIN_FCHR_DOLLARA(x)             \
	assert_leq(x[0], this->fchr()[1] + 1); \
	assert_leq(x[1], this->fchr()[2]);     \
	assert_leq(x[2], this->fchr()[3]);     \
	assert_leq(x[3], this->fchr()[4])
	inline TIndexOffU countBt2Side(const SideLocus &l, int c) const
	{
#ifdef TIME_STATS
		countLF_Static++;
#endif
		assert_range(0, 3, c);
		assert_range(0, (int)this->_eh._sideBwtSz - 1, (int)l._by);
		assert_range(0, 3, (int)l._bp);
		const uint8_t *side = l.side(this->ebwt());
		TIndexOffU cCnt = countUpTo(l, c);
		assert_leq(cCnt, l.toBWRow());
		assert_leq(cCnt, this->_eh._sideBwtLen);
		if (c == 0 && l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff)
		{
			if ((l._sideByteOff + l._by > _zEbwtByteOff) ||
				(l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
			{
				cCnt--;
			}
		}
		TIndexOffU ret;
		const uint8_t *acgt8 = side + _eh._sideBwtSz;
		const TIndexOffU *acgt = reinterpret_cast<const TIndexOffU *>(acgt8);
		assert_leq(acgt[0], this->_eh._numSides * this->_eh._sideBwtLen);
		assert_leq(acgt[1], this->_eh._len);
		assert_leq(acgt[2], this->_eh._len);
		assert_leq(acgt[3], this->_eh._len);
		ret = acgt[c] + cCnt + this->fchr()[c];
#ifndef NDEBUG
		assert_leq(ret, this->fchr()[c + 1]);
		if (c == 0)
		{
			assert_leq(cCnt, this->_eh._sideBwtLen);
		}
		else
		{
			assert_leq(ret, this->_eh._bwtLen);
		}
#endif
		return ret;
	}
	inline void countBt2SideRange(
		SideLocus &l,
		TIndexOffU num,
		TIndexOffU *cntsUpto,
		TIndexOffU *cntsIn,
		EList<bool> *masks) const
	{
		assert_gt(num, 0);
		assert_range(0, (int)this->_eh._sideBwtSz - 1, (int)l._by);
		assert_range(0, 3, (int)l._bp);
		countUpToEx(l, cntsUpto);
		WITHIN_FCHR_DOLLARA(cntsUpto);
		WITHIN_BWT_LEN(cntsUpto);
		const uint8_t *side = l.side(this->ebwt());
		if (l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff)
		{
			if ((l._sideByteOff + l._by > _zEbwtByteOff) ||
				(l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
			{
				cntsUpto[0]--;
			}
		}
		const TIndexOffU *acgt = reinterpret_cast<const TIndexOffU *>(side + _eh._sideBwtSz);
		assert_leq(acgt[0], this->fchr()[1] + this->_eh.sideBwtLen());
		assert_leq(acgt[1], this->fchr()[2] - this->fchr()[1]);
		assert_leq(acgt[2], this->fchr()[3] - this->fchr()[2]);
		assert_leq(acgt[3], this->fchr()[4] - this->fchr()[3]);
		assert_leq(acgt[0], this->_eh._len + this->_eh.sideBwtLen());
		assert_leq(acgt[1], this->_eh._len);
		assert_leq(acgt[2], this->_eh._len);
		assert_leq(acgt[3], this->_eh._len);
		cntsUpto[0] += (acgt[0] + this->fchr()[0]);
		cntsUpto[1] += (acgt[1] + this->fchr()[1]);
		cntsUpto[2] += (acgt[2] + this->fchr()[2]);
		cntsUpto[3] += (acgt[3] + this->fchr()[3]);
		masks[0].resize(num);
		masks[1].resize(num);
		masks[2].resize(num);
		masks[3].resize(num);
		WITHIN_FCHR_DOLLARA(cntsUpto);
		WITHIN_FCHR_DOLLARA(cntsIn);
		TIndexOffU nm = 0;
		nm += countBt2SideRange2(l, true, num - nm, cntsIn, masks, nm);
		assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
		assert_leq(nm, num);
		SideLocus lcopy = l;
		while (nm < num)
		{
			lcopy.nextSide(this->_eh);
			nm += countBt2SideRange2(lcopy, false, num - nm, cntsIn, masks, nm);
			WITHIN_FCHR_DOLLARA(cntsIn);
			assert_leq(nm, num);
			assert_eq(nm, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
		}
		assert_eq(num, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
		WITHIN_FCHR_DOLLARA(cntsIn);
	}
	inline void countBt2SideEx(const SideLocus &l, TIndexOffU *arrs) const
	{
		assert_range(0, (int)this->_eh._sideBwtSz - 1, (int)l._by);
		assert_range(0, 3, (int)l._bp);
		countUpToEx(l, arrs);
		if (l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff)
		{
			if ((l._sideByteOff + l._by > _zEbwtByteOff) ||
				(l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
			{
				arrs[0]--;
			}
		}
		WITHIN_FCHR(arrs);
		WITHIN_BWT_LEN(arrs);
		const uint8_t *side = l.side(this->ebwt());
		const uint8_t *acgt16 = side + this->_eh._sideSz - OFF_SIZE * 4;
		const TIndexOffU *acgt = reinterpret_cast<const TIndexOffU *>(acgt16);
		assert_leq(acgt[0], this->fchr()[1] + this->_eh.sideBwtLen());
		assert_leq(acgt[1], this->fchr()[2] - this->fchr()[1]);
		assert_leq(acgt[2], this->fchr()[3] - this->fchr()[2]);
		assert_leq(acgt[3], this->fchr()[4] - this->fchr()[3]);
		assert_leq(acgt[0], this->_eh._len + this->_eh.sideBwtLen());
		assert_leq(acgt[1], this->_eh._len);
		assert_leq(acgt[2], this->_eh._len);
		assert_leq(acgt[3], this->_eh._len);
		arrs[0] += (acgt[0] + this->fchr()[0]);
		arrs[1] += (acgt[1] + this->fchr()[1]);
		arrs[2] += (acgt[2] + this->fchr()[2]);
		arrs[3] += (acgt[3] + this->fchr()[3]);
		WITHIN_FCHR(arrs);
	}
	inline TIndexOffU countUpTo(const SideLocus &l, int c) const
	{
		const uint8_t *side = l.side(this->ebwt());
		TIndexOffU i = 0, cCnt = 0, cCnt_add = 0;
		TIndexOffU p_by = l._charOff >> 3;
		for (; i + 7 < p_by; i += 8)
		{
			switch (c)
			{
			case 0:
				cCnt += countInU64<USE_POPCNT_INSTRUCTION>(c, ~((*(uint64_t *)&side[i]) | (*(uint64_t *)&side[i + 56])));
				break;
			case 1:
				cCnt += countInU64<USE_POPCNT_INSTRUCTION>(c, (~(*(uint64_t *)&side[i])) & (*(uint64_t *)&side[i + 56]));
				break;
			case 2:
				cCnt += countInU64<USE_POPCNT_INSTRUCTION>(c, (*(uint64_t *)&side[i]) & (~(*(uint64_t *)&side[i + 56])));
				break;
			case 3:
				cCnt += countInU64<USE_POPCNT_INSTRUCTION>(c, (*(uint64_t *)&side[i]) & (*(uint64_t *)&side[i + 56]));
				break;
			}
		}
		if ((l._charOff & 0xffffffc0) == l._charOff)
		{
			return cCnt;
		}
		TIndexOffU shift = 64 + (i << 3) - l._charOff;
		uint64_t side_i_hi = *(uint64_t *)&side[i];
		uint64_t side_i_lo = *(uint64_t *)&side[i + 56];
		side_i_hi = (_toBigEndian ? side_i_hi >> shift : side_i_hi << shift);
		side_i_lo = (_toBigEndian ? side_i_lo >> shift : side_i_lo << shift);
		switch (c)
		{
		case 0:
			cCnt_add = (countInU64<USE_POPCNT_INSTRUCTION>(0, ~(side_i_hi | side_i_lo)) - shift);
			break;
		case 1:
			cCnt_add = countInU64<USE_POPCNT_INSTRUCTION>(c, (~side_i_hi) & side_i_lo);
			break;
		case 2:
			cCnt_add = countInU64<USE_POPCNT_INSTRUCTION>(c, side_i_hi & (~side_i_lo));
			break;
		case 3:
			cCnt_add = countInU64<USE_POPCNT_INSTRUCTION>(c, side_i_hi & side_i_lo);
			break;
		}
		return cCnt + cCnt_add;
	}
#ifdef POPCNT_CAPABILITY
	template <typename Operation>
#endif
	inline static void countInU64Ex(uint64_t dw, TIndexOffU *arrs)
	{
		uint64_t c0 = c_table[0];
		uint64_t x0 = dw ^ c0;
		uint64_t x1 = (x0 >> 1);
		uint64_t x2 = x1 & (0x5555555555555555llu);
		uint64_t x3 = x0 & x2;
#ifdef POPCNT_CAPABILITY
		uint64_t tmp = Operation().pop64(x3);
#else
		uint64_t tmp = pop64(x3);
#endif
		arrs[0] += (uint32_t)tmp;
		c0 = c_table[1];
		x0 = dw ^ c0;
		x1 = (x0 >> 1);
		x2 = x1 & (0x5555555555555555llu);
		x3 = x0 & x2;
#ifdef POPCNT_CAPABILITY
		tmp = Operation().pop64(x3);
#else
		tmp = pop64(x3);
#endif
		arrs[1] += (uint32_t)tmp;
		c0 = c_table[2];
		x0 = dw ^ c0;
		x1 = (x0 >> 1);
		x2 = x1 & (0x5555555555555555llu);
		x3 = x0 & x2;
#ifdef POPCNT_CAPABILITY
		tmp = Operation().pop64(x3);
#else
		tmp = pop64(x3);
#endif
		arrs[2] += (uint32_t)tmp;
		c0 = c_table[3];
		x0 = dw ^ c0;
		x1 = (x0 >> 1);
		x2 = x1 & (0x5555555555555555llu);
		x3 = x0 & x2;
#ifdef POPCNT_CAPABILITY
		tmp = Operation().pop64(x3);
#else
		tmp = pop64(x3);
#endif
		arrs[3] += (uint32_t)tmp;
	}
	inline void countUpToEx(const SideLocus &l, TIndexOffU *arrs) const
	{
		const uint8_t *side = l.side(this->ebwt());
		TIndexOffU i = 0;
		TIndexOffU p_by = l._charOff >> 3;
		for (; i + 7 < p_by; i += 8)
		{
			arrs[0] += countInU64<USE_POPCNT_INSTRUCTION>(0, ~((*(uint64_t *)&side[i]) | (*(uint64_t *)&side[i + 56])));
			arrs[1] += countInU64<USE_POPCNT_INSTRUCTION>(1, (~(*(uint64_t *)&side[i])) & (*(uint64_t *)&side[i + 56]));
			arrs[2] += countInU64<USE_POPCNT_INSTRUCTION>(2, (*(uint64_t *)&side[i]) & (~(*(uint64_t *)&side[i + 56])));
			arrs[3] += countInU64<USE_POPCNT_INSTRUCTION>(3, (*(uint64_t *)&side[i]) & (*(uint64_t *)&side[i + 56]));
		}
		if ((l._charOff & 0xffffffc0) == l._charOff)
			return;
		TIndexOffU shift = 64 + (i << 3) - l._charOff;
		uint64_t side_i_hi = *(uint64_t *)&side[i];
		uint64_t side_i_lo = *(uint64_t *)&side[i + 56];
		side_i_hi = (_toBigEndian ? side_i_hi >> shift : side_i_hi << shift);
		side_i_lo = (_toBigEndian ? side_i_lo >> shift : side_i_lo << shift);
		arrs[0] += (countInU64<USE_POPCNT_INSTRUCTION>(0, ~(side_i_hi | side_i_lo)) - shift);
		arrs[1] += countInU64<USE_POPCNT_INSTRUCTION>(1, (~side_i_hi) & side_i_lo);
		arrs[2] += countInU64<USE_POPCNT_INSTRUCTION>(2, side_i_hi & (~side_i_lo));
		arrs[3] += countInU64<USE_POPCNT_INSTRUCTION>(3, side_i_hi & side_i_lo);
	}
#ifndef NDEBUG
	inline void mapLFEx(
		const SideLocus &l,
		TIndexOffU *arrs
			ASSERT_ONLY(, bool overrideSanity = false)) const
	{
		assert_eq(0, arrs[0]);
		assert_eq(0, arrs[1]);
		assert_eq(0, arrs[2]);
		assert_eq(0, arrs[3]);
		countBt2SideEx(l, arrs);
		if (_sanity && !overrideSanity)
		{
			assert_eq(mapLF(l, 0, true), arrs[0]);
			assert_eq(mapLF(l, 1, true), arrs[1]);
			assert_eq(mapLF(l, 2, true), arrs[2]);
			assert_eq(mapLF(l, 3, true), arrs[3]);
		}
	}
#endif
	inline void mapLFEx(
		TIndexOffU top,
		TIndexOffU bot,
		TIndexOffU *tops,
		TIndexOffU *bots
			ASSERT_ONLY(, bool overrideSanity = false)) const
	{
		SideLocus ltop, lbot;
		SideLocus::initFromTopBot(top, bot, _eh, ebwt(), ltop, lbot);
		mapLFEx(ltop, lbot, tops, bots ASSERT_ONLY(, overrideSanity));
	}
	inline void mapLFEx(
		const SideLocus &ltop,
		const SideLocus &lbot,
		TIndexOffU *tops,
		TIndexOffU *bots
			ASSERT_ONLY(, bool overrideSanity = false)) const
	{
		assert(ltop.repOk(this->eh()));
		assert(lbot.repOk(this->eh()));
		assert_eq(0, tops[0]);
		assert_eq(0, bots[0]);
		assert_eq(0, tops[1]);
		assert_eq(0, bots[1]);
		assert_eq(0, tops[2]);
		assert_eq(0, bots[2]);
		assert_eq(0, tops[3]);
		assert_eq(0, bots[3]);
		countBt2SideEx(ltop, tops);
		countBt2SideEx(lbot, bots);
#ifndef NDEBUG
		if (_sanity && !overrideSanity)
		{
			assert_eq(mapLF(ltop, 0, true), tops[0]);
			assert_eq(mapLF(ltop, 1, true), tops[1]);
			assert_eq(mapLF(ltop, 2, true), tops[2]);
			assert_eq(mapLF(ltop, 3, true), tops[3]);
			assert_eq(mapLF(lbot, 0, true), bots[0]);
			assert_eq(mapLF(lbot, 1, true), bots[1]);
			assert_eq(mapLF(lbot, 2, true), bots[2]);
			assert_eq(mapLF(lbot, 3, true), bots[3]);
		}
#endif
	}
	inline TIndexOffU countBt2SideRange2(
		const SideLocus &l,
		bool startAtLocus,
		TIndexOffU num,
		TIndexOffU *arrs,
		EList<bool> *masks,
		TIndexOffU maskOff) const
	{
		assert(!masks[0].empty());
		assert_eq(masks[0].size(), masks[1].size());
		assert_eq(masks[0].size(), masks[2].size());
		assert_eq(masks[0].size(), masks[3].size());
		ASSERT_ONLY(TIndexOffU myarrs[4] = {0, 0, 0, 0});
		TIndexOffU nm = 0;
		int iby = 0;
		int ibp = 0;
		if (startAtLocus)
		{
			iby = l._by;
			ibp = l._bp;
		}
		else
		{
		}
		int by = iby, bp = ibp;
		assert_lt(bp, 4);
		assert_lt(by, (int)this->_eh._sideBwtSz);
		const uint8_t *side = l.side(this->ebwt());
		int side_i_hi = 0, side_i_lo = 0, t_by = 0, t_bp = 0;
		while (nm < num)
		{
			t_by = by >> 1;
			t_bp = bp + ((by & 0x1) << 2);
			side_i_hi = unpack_1b_from_8b(side[t_by], t_bp);
			side_i_lo = unpack_1b_from_8b(side[t_by + 56], t_bp);
			int c = (side_i_hi << 1) | side_i_lo;
			assert_lt(maskOff + nm, masks[c].size());
			masks[0][maskOff + nm] = masks[1][maskOff + nm] =
				masks[2][maskOff + nm] = masks[3][maskOff + nm] = false;
			assert_range(0, 3, c);
			arrs[c]++;
			ASSERT_ONLY(myarrs[c]++);
			masks[c][maskOff + nm] = true;
			nm++;
			if (++bp == 4)
			{
				bp = 0;
				by++;
				assert_leq(by, (int)this->_eh._sideBwtSz);
				if (by == (int)this->_eh._sideBwtSz)
				{
					break;
				}
			}
		}
		WITHIN_FCHR_DOLLARA(arrs);
#ifndef NDEBUG
		if (_sanity)
		{
			TIndexOffU tops[4] = {0, 0, 0, 0};
			TIndexOffU bots[4] = {0, 0, 0, 0};
			TIndexOffU top = l.toBWRow();
			TIndexOffU bot = top + nm;
			mapLFEx(top, bot, tops, bots, false);
			assert(myarrs[0] == (bots[0] - tops[0]) || myarrs[0] == (bots[0] - tops[0]) + 1);
			assert_eq(myarrs[1], bots[1] - tops[1]);
			assert_eq(myarrs[2], bots[2] - tops[2]);
			assert_eq(myarrs[3], bots[3] - tops[3]);
		}
#endif
		return nm;
	}
	inline int rowL(const SideLocus &l) const
	{
		uint32_t p_by = l._charOff >> 3;
		uint32_t p_bp = l._charOff & 0x7;
		int side_i_hi = unpack_1b_from_8b(l.side(this->ebwt())[p_by], p_bp);
		int side_i_lo = unpack_1b_from_8b(l.side(this->ebwt())[p_by + 56], p_bp);
		return ((side_i_hi << 1) | side_i_lo);
	}
	inline int rowL(TIndexOffU i) const
	{
		SideLocus l;
		l.initFromRow(i, _eh, ebwt());
		return rowL(l);
	}
	inline void mapLFRange(
		SideLocus &ltop,
		SideLocus &lbot,
		TIndexOffU num,
		TIndexOffU *cntsUpto,
		TIndexOffU *cntsIn,
		EList<bool> *masks
			ASSERT_ONLY(, bool overrideSanity = false)) const
	{
		assert(ltop.repOk(this->eh()));
		assert(lbot.repOk(this->eh()));
		assert_eq(num, lbot.toBWRow() - ltop.toBWRow());
		assert_eq(0, cntsUpto[0]);
		assert_eq(0, cntsIn[0]);
		assert_eq(0, cntsUpto[1]);
		assert_eq(0, cntsIn[1]);
		assert_eq(0, cntsUpto[2]);
		assert_eq(0, cntsIn[2]);
		assert_eq(0, cntsUpto[3]);
		assert_eq(0, cntsIn[3]);
		countBt2SideRange(ltop, num, cntsUpto, cntsIn, masks);
		assert_eq(num, cntsIn[0] + cntsIn[1] + cntsIn[2] + cntsIn[3]);
#ifndef NDEBUG
		if (_sanity && !overrideSanity)
		{
			TIndexOffU tops[4] = {0, 0, 0, 0};
			TIndexOffU bots[4] = {0, 0, 0, 0};
			assert(ltop.repOk(this->eh()));
			assert(lbot.repOk(this->eh()));
			mapLFEx(ltop, lbot, tops, bots, false);
			for (int i = 0; i < 4; i++)
			{
				assert(cntsUpto[i] == tops[i] || tops[i] == bots[i]);
				if (i == 0)
				{
					assert(cntsIn[i] == bots[i] - tops[i] ||
						   cntsIn[i] == bots[i] - tops[i] + 1);
				}
				else
				{
					assert_eq(cntsIn[i], bots[i] - tops[i]);
				}
			}
		}
#endif
	}
	inline TIndexOffU mapLF(
		const SideLocus &l
			ASSERT_ONLY(, bool overrideSanity = false)) const
	{
		ASSERT_ONLY(TIndexOffU srcrow = l.toBWRow());
		TIndexOffU ret;
		assert(l.side(this->ebwt()) != NULL);
		int c = rowL(l);
		assert_lt(c, 4);
		assert_geq(c, 0);
		ret = countBt2Side(l, c);
		cerr << "2.mapLF -> " << c << "  " << ret << "  " << l._charOff << "  " << l._sideByteOff << "  " << l._sideNum << "  " << l._by << "  " << l._bp << endl
			 << endl;
		assert_lt(ret, this->_eh._bwtLen);
		assert_neq(srcrow, ret);
#ifndef NDEBUG
		if (_sanity && !overrideSanity)
		{
			TIndexOffU arrs[] = {0, 0, 0, 0};
			mapLFEx(l, arrs, true);
			assert_eq(arrs[c], ret);
		}
#endif
		return ret;
	}
	inline TIndexOffU mapLF(
		const SideLocus &l, int c
								ASSERT_ONLY(, bool overrideSanity = false)) const
	{
		TIndexOffU ret;
		assert_lt(c, 4);
		assert_geq(c, 0);
		ret = countBt2Side(l, c);
		assert_lt(ret, this->_eh._bwtLen);
#ifndef NDEBUG
		if (_sanity && !overrideSanity)
		{
			TIndexOffU arrs[] = {0, 0, 0, 0};
			mapLFEx(l, arrs, true);
			assert_eq(arrs[c], ret);
		}
#endif
		return ret;
	}
	inline void mapBiLFEx(
		const SideLocus &ltop,
		const SideLocus &lbot,
		TIndexOffU *tops,
		TIndexOffU *bots,
		TIndexOffU *topsP,
		TIndexOffU *botsP
			ASSERT_ONLY(, bool overrideSanity = false)) const
	{
#ifndef NDEBUG
		for (int i = 0; i < 4; i++)
		{
			assert_eq(0, tops[0]);
			assert_eq(0, bots[0]);
		}
#endif
		countBt2SideEx(ltop, tops);
		countBt2SideEx(lbot, bots);
#ifndef NDEBUG
		if (_sanity && !overrideSanity)
		{
			assert_eq(mapLF(ltop, 0, true), tops[0]);
			assert_eq(mapLF(ltop, 1, true), tops[1]);
			assert_eq(mapLF(ltop, 2, true), tops[2]);
			assert_eq(mapLF(ltop, 3, true), tops[3]);
			assert_eq(mapLF(lbot, 0, true), bots[0]);
			assert_eq(mapLF(lbot, 1, true), bots[1]);
			assert_eq(mapLF(lbot, 2, true), bots[2]);
			assert_eq(mapLF(lbot, 3, true), bots[3]);
		}
#endif
		botsP[0] = topsP[0] + (bots[0] - tops[0]);
		topsP[1] = botsP[0];
		botsP[1] = topsP[1] + (bots[1] - tops[1]);
		topsP[2] = botsP[1];
		botsP[2] = topsP[2] + (bots[2] - tops[2]);
		topsP[3] = botsP[2];
		botsP[3] = topsP[3] + (bots[3] - tops[3]);
	}
	inline TIndexOffU mapLF1(
		TIndexOffU row,
		const SideLocus &l,
		int c
			ASSERT_ONLY(, bool overrideSanity = false)) const
	{
		if (rowL(l) != c || row == _zOff)
			return OFF_MASK;
		TIndexOffU ret;
		assert_lt(c, 4);
		assert_geq(c, 0);
		ret = countBt2Side(l, c);
		assert_lt(ret, this->_eh._bwtLen);
#ifndef NDEBUG
		if (_sanity && !overrideSanity)
		{
			TIndexOffU arrs[] = {0, 0, 0, 0};
			mapLFEx(l, arrs, true);
			assert_eq(arrs[c], ret);
		}
#endif
		return ret;
	}
	inline int mapLF1(
		TIndexOffU &row,
		const SideLocus &l
			ASSERT_ONLY(, bool overrideSanity = false)) const
	{
		if (row == _zOff)
			return -1;
		int c = rowL(l);
		assert_range(0, 3, c);
		row = countBt2Side(l, c);
		assert_lt(row, this->_eh._bwtLen);
#ifndef NDEBUG
		if (_sanity && !overrideSanity)
		{
			TIndexOffU arrs[] = {0, 0, 0, 0};
			mapLFEx(l, arrs, true);
			assert_eq(arrs[c], row);
		}
#endif
		return c;
	}
#ifndef NDEBUG
	bool inMemoryRepOk(const EbwtParams &eh) const
	{
		assert_geq(_zEbwtBpOff, 0);
		assert_lt(_zEbwtBpOff, 4);
		assert_lt(_zEbwtByteOff, eh._ebwtTotSz);
		assert_lt(_zOff, eh._bwtLen);
		assert_geq(_nFrag, _nPat);
		return true;
	}
	bool inMemoryRepOk() const
	{
		return repOk(_eh);
	}
	bool repOk(const EbwtParams &eh) const
	{
		assert(_eh.repOk());
		if (isInMemory())
		{
			return inMemoryRepOk(eh);
		}
		return true;
	}
	bool repOk() const
	{
		return repOk(_eh);
	}
#endif
	bool _toBigEndian;
	int32_t _overrideOffRate;
	bool _verbose;
	bool _passMemExc;
	bool _sanity;
	bool fw_;
	FILE *_in1;
	FILE *_in2;
	string _in1Str;
	string _in2Str;
	string _in6Str;
	string _inSaStr;
	string _inBwtStr;
	TIndexOffU _zOff;
	TIndexOffU _zEbwtByteOff;
	TIndexOff _zEbwtBpOff;
	TIndexOffU _nPat;
	TIndexOffU _nFrag;
	APtrWrap<TIndexOffU> _plen;
	APtrWrap<TIndexOffU> _rstarts;
	APtrWrap<TIndexOffU> _fchr;
	APtrWrap<TIndexOffU> _ftab;
	APtrWrap<TIndexOffU> _eftab;
	APtrWrap<TIndexOffU> _offs;
	APtrWrap<uint8_t> _ebwt;
	bool _useMm;
	bool useShmem_;
	EList<string> _refnames;
	char *mmFile1_;
	char *mmFile2_;
	EbwtParams _eh;
	bool packed_;
	static const TIndexOffU default_bmax = OFF_MASK;
	static const TIndexOffU default_bmaxMultSqrt = OFF_MASK;
	static const TIndexOffU default_bmaxDivN = 4;
	static const int default_dcv = 1024;
	static const bool default_noDc = false;
	static const bool default_useBlockwise = true;
	static const uint32_t default_seed = 0;
#ifdef _64BIT_INDEX
	static const int default_lineRate = 7;
#else
	static const int default_lineRate = 7;
#endif
	static const int default_offRate = 5;
	static const int default_offRatePlus = 0;
	static const int default_ftabChars = 10;
	static const bool default_bigEndian = false;
private:
	ostream &log() const
	{
		return cout;
	}
	void verbose(const string &s) const
	{
		if (this->verbose())
		{
			this->log() << s.c_str();
			this->log().flush();
		}
	}
};
void readEbwtRefnames(FILE *fin, EList<string> &refnames);
void readEbwtRefnames(const string &instr, EList<string> &refnames);
bool readEbwtColor(const string &instr);
bool readEntireReverse(const string &instr);
template <typename TStr>
TStr Ebwt::join(EList<TStr> &l, uint32_t seed)
{
	RandomSource rand;
	rand.init(seed);
	TStr ret;
	TIndexOffU guessLen = 0;
	for (TIndexOffU i = 0; i < l.size(); i++)
	{
		guessLen += length(l[i]);
	}
	ret.resize(guessLen);
	TIndexOffU off = 0;
	for (size_t i = 0; i < l.size(); i++)
	{
		TStr &s = l[i];
		assert_gt(s.length(), 0);
		for (size_t j = 0; j < s.size(); j++)
		{
			ret.set(s[j], off++);
		}
	}
	return ret;
}
template <typename TStr>
TStr Ebwt::join(EList<FileBuf *> &l,
				EList<RefRecord> &szs,
				TIndexOffU sztot,
				const RefReadInParams &refparams,
				uint32_t seed)
{
	RandomSource rand;
	rand.init(seed);
	RefReadInParams rpcp = refparams;
	TStr ret;
	TIndexOffU guessLen = sztot;
	ret.resize(guessLen);
	ASSERT_ONLY(TIndexOffU szsi = 0);
	TIndexOffU dstoff = 0;
	for (TIndexOffU i = 0; i < l.size(); i++)
	{
		assert(!l[i]->eof());
		bool first = true;
		while (!l[i]->eof())
		{
			RefRecord rec = fastaRefReadAppend(*l[i], first, ret, dstoff, rpcp);
			first = false;
			if (rec.first && rec.len == 0)
			{
				continue;
			}
			TIndexOffU bases = rec.len;
			assert_eq(rec.off, szs[szsi].off);
			assert_eq(rec.len, szs[szsi].len);
			assert_eq(rec.first, szs[szsi].first);
			ASSERT_ONLY(szsi++);
			if (bases == 0)
				continue;
		}
	}
	return ret;
}
template <typename TStr>
void Ebwt::joinToDisk(
	EList<FileBuf *> &l,
	EList<RefRecord> &szs,
	TIndexOffU sztot,
	const RefReadInParams &refparams,
	TStr &ret,
	ostream &out1,
	ostream &out2)
{
	RefReadInParams rpcp = refparams;
	assert_gt(szs.size(), 0);
	assert_gt(l.size(), 0);
	assert_gt(sztot, 0);
	this->_nPat = 0;
	this->_nFrag = 0;
	for (TIndexOffU i = 0; i < szs.size(); i++)
	{
		if (szs[i].len > 0)
			this->_nFrag++;
		if (szs[i].first && szs[i].len > 0)
			this->_nPat++;
	}
	assert_gt(this->_nPat, 0);
	assert_geq(this->_nFrag, this->_nPat);
	_rstarts.reset();
	writeU<TIndexOffU>(out1, this->_nPat, this->toBe());
	try
	{
		this->_plen.init(new TIndexOffU[this->_nPat], this->_nPat);
	}
	catch (bad_alloc &e)
	{
		cerr << "Out of memory allocating plen[] in Ebwt::join()"
			 << " at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}
	TIndexOff npat = -1;
	for (TIndexOffU i = 0; i < szs.size(); i++)
	{
		if (szs[i].first && szs[i].len > 0)
		{
			if (npat >= 0)
			{
				writeU<TIndexOffU>(out1, this->plen()[npat], this->toBe());
			}
			this->plen()[++npat] = (szs[i].len + szs[i].off);
		}
		else if (!szs[i].first)
		{
			if (npat < 0)
				npat = 0;
			this->plen()[npat] += (szs[i].len + szs[i].off);
		}
	}
	assert_eq((TIndexOffU)npat, this->_nPat - 1);
	writeU<TIndexOffU>(out1, this->plen()[npat], this->toBe());
	writeU<TIndexOffU>(out1, this->_nFrag, this->toBe());
	TIndexOffU seqsRead = 0;
	ASSERT_ONLY(TIndexOffU szsi = 0);
	ASSERT_ONLY(TIndexOffU entsWritten = 0);
	TIndexOffU dstoff = 0;
	for (unsigned int i = 0; i < l.size(); i++)
	{
		assert(!l[i]->eof());
		bool first = true;
		TIndexOffU patoff = 0;
		while (!l[i]->eof())
		{
			string name;
			_refnames.push_back("");
			RefRecord rec = fastaRefReadAppend(
				*l[i], first, ret, dstoff, rpcp, &_refnames.back());
			first = false;
			TIndexOffU bases = rec.len;
			if (rec.first && rec.len > 0)
			{
				if (_refnames.back().length() == 0)
				{
					ostringstream stm;
					stm << seqsRead;
					_refnames.back() = stm.str();
				}
			}
			else
			{
				_refnames.pop_back();
			}
			if (rec.first && rec.len == 0)
			{
				continue;
			}
			assert_lt(szsi, szs.size());
			assert_eq(rec.off, szs[szsi].off);
			assert_eq(rec.len, szs[szsi].len);
			assert_eq(rec.first, szs[szsi].first);
			assert(rec.first || rec.off > 0);
			ASSERT_ONLY(szsi++);
			if (rec.first)
				seqsRead++;
			if (bases == 0)
				continue;
			assert_leq(bases, this->plen()[seqsRead - 1]);
			if (rec.first)
				patoff = 0;
			patoff += rec.off;
			ASSERT_ONLY(entsWritten++);
			patoff += bases;
		}
		assert_gt(szsi, 0);
		l[i]->reset();
		assert(!l[i]->eof());
#ifndef NDEBUG
		int c = l[i]->get();
		assert_eq('>', c);
		assert(!l[i]->eof());
		l[i]->reset();
		assert(!l[i]->eof());
#endif
	}
	assert_eq(entsWritten, this->_nFrag);
}
template <typename TStr>
void Ebwt::buildToDisk(
	InorderBlockwiseSA<TStr> &sa,
	const TStr &s,
	ostream &out1,
	ostream &out2,
	string &_in6Str,
	ostream *saOut,
	ostream *bwtOut)
{
	const EbwtParams &eh = this->_eh;
	assert(eh.repOk());
	assert_eq(s.length() + 1, sa.size());
	assert_eq(s.length(), eh._len);
	assert_gt(eh._lineRate, 3);
	assert(sa.suffixItrIsReset());
	TIndexOffU len = eh._len;
	TIndexOffU ftabLen = eh._ftabLen;
	TIndexOffU sideSz = eh._sideSz;
	TIndexOffU ebwtTotSz = eh._ebwtTotSz;
	TIndexOffU fchr[] = {0, 0, 0, 0, 0};
	EList<TIndexOffU> ftab(EBWT_CAT);
	TIndexOffU zOff = OFF_MASK;
	TIndexOffU occ[4] = {0, 0, 0, 0};
	TIndexOffU occSave[4] = {0, 0, 0, 0};
	uint8_t absorbCnt = 0;
	EList<uint8_t> absorbFtab(EBWT_CAT);
	try
	{
		VMSG_NL("Allocating ftab, absorbFtab");
		ftab.resize(ftabLen);
		ftab.fillZero();
		absorbFtab.resize(ftabLen);
		absorbFtab.fillZero();
	}
	catch (bad_alloc &e)
	{
		cerr << "Out of memory allocating ftab[] or absorbFtab[] "
			 << "in Ebwt::buildToDisk() at " << __FILE__ << ":"
			 << __LINE__ << endl;
		throw e;
	}
#ifdef SIXTY4_FORMAT
	EList<uint64_t> ebwtSide(EBWT_CAT);
#else
	EList<uint8_t> ebwtSide(EBWT_CAT);
#endif
	try
	{
#ifdef SIXTY4_FORMAT
		ebwtSide.resize(sideSz >> 3);
#else
		ebwtSide.resize(sideSz);
#endif
	}
	catch (bad_alloc &e)
	{
		cerr << "Out of memory allocating ebwtSide[] in "
			 << "Ebwt::buildToDisk() at " << __FILE__ << ":"
			 << __LINE__ << endl;
		throw e;
	}
	TIndexOffU side = 0;
	bool fw;
	TIndexOff sideCur = 0;
	fw = true;
	ASSERT_ONLY(bool dollarSkipped = false);
	TIndexOffU si = 0;
	ASSERT_ONLY(TIndexOffU lastSufInt = 0);
	ASSERT_ONLY(bool inSA = true);
	VMSG_NL("Entering Ebwt loop");
	ASSERT_ONLY(TIndexOffU beforeEbwtOff = (TIndexOffU)out1.tellp());
	if (saOut != NULL)
	{
		writeU<TIndexOffU>(*saOut, len + 1, this->toBe());
	}
	if (bwtOut != NULL)
	{
		writeU<TIndexOffU>(*bwtOut, len + 1, this->toBe());
	}
	while (side < ebwtTotSz)
	{
		assert_geq(sideCur, 0);
		assert_lt(sideCur, (int)eh._sideBwtSz);
		assert_eq(0, side % sideSz);
		ebwtSide[sideCur] = 0;
		ebwtSide[sideCur + 56] = 0;
		assert_lt(side + sideCur, ebwtTotSz);
#ifdef SIXTY4_FORMAT
		for (int bpi = 0; bpi < 32; bpi++, si++)
#else
		for (int bpi = 0; bpi < 8; bpi++, si++)
#endif
		{
			int bwtChar;
			bool count = true;
			if (si <= len)
			{
				TIndexOffU saElt = sa.nextSuffix();
				if (saOut != NULL)
				{
					writeU<TIndexOffU>(*saOut, saElt, this->toBe());
				}
				if (saElt == 0)
				{
					bwtChar = 0;
					count = false;
					ASSERT_ONLY(dollarSkipped = true);
					zOff = si;
				}
				else
				{
					bwtChar = (int)(s[saElt - 1]);
					assert_lt(bwtChar, 4);
					fchr[bwtChar]++;
				}
				if ((len - saElt) >= (TIndexOffU)eh._ftabChars)
				{
					TIndexOffU sufInt = 0;
					for (int i = 0; i < eh._ftabChars; i++)
					{
						sufInt <<= 2;
						assert_lt((TIndexOffU)i, len - saElt);
						sufInt |= (unsigned char)(s[saElt + i]);
					}
#ifndef NDEBUG
					if (lastSufInt > 0)
						assert_geq(sufInt, lastSufInt);
					lastSufInt = sufInt;
#endif
					assert_lt(sufInt + 1, ftabLen);
					ftab[sufInt + 1]++;
					if (absorbCnt > 0)
					{
						absorbFtab[sufInt] = absorbCnt;
						absorbCnt = 0;
					}
				}
				else
				{
					assert_lt(absorbCnt, 255);
					absorbCnt++;
				}
				if ((si & eh._offMask) == si)
				{
					assert_lt((si >> eh._offRate), eh._offsLen);
					writeU<TIndexOffU>(out2, saElt, this->toBe());
				}
			}
			else
			{
#ifndef NDEBUG
				if (inSA)
				{
					assert_eq(si, len + 1);
					inSA = false;
				}
#endif
				bwtChar = 0;
			}
			if (count)
				occ[bwtChar]++;
			if (fw)
			{
#ifdef SIXTY4_FORMAT
				ebwtSide[sideCur] |= ((uint64_t)bwtChar << (bpi << 1));
				if (bwtChar > 0)
					assert_gt(ebwtSide[sideCur], 0);
#else
				pack_1b_in_8b(bwtChar >> 1, ebwtSide[sideCur], bpi);
				pack_1b_in_8b(bwtChar & 0x1, ebwtSide[sideCur + 56], bpi);
				assert_eq((ebwtSide[sideCur] >> (bpi * 2)) & 3, bwtChar);
#endif
			}
			else
			{
#ifdef SIXTY4_FORMAT
				ebwtSide[sideCur] |= ((uint64_t)bwtChar << ((31 - bpi) << 1));
				if (bwtChar > 0)
					assert_gt(ebwtSide[sideCur], 0);
#else
				pack_1b_in_8b(bwtChar >> 1, ebwtSide[sideCur], 7 - bpi);
				pack_1b_in_8b(bwtChar & 0x1, ebwtSide[sideCur + 56], 7 - bpi);
				assert_eq((ebwtSide[sideCur] >> ((3 - bpi) * 2)) & 3, bwtChar);
#endif
			}
		}
		assert_eq(dollarSkipped ? 3 : 0, (occ[0] + occ[1] + occ[2] + occ[3]) & 3);
#ifdef SIXTY4_FORMAT
		assert_eq(0, si & 31);
#else
		assert_eq(0, si & 3);
#endif
		sideCur++;
		if (sideCur << 1 == (int)eh._sideBwtSz)
		{
			sideCur = 0;
			TIndexOffU *cpptr = reinterpret_cast<TIndexOffU *>(ebwtSide.ptr());
			side += sideSz;
			assert_leq(side, eh._ebwtTotSz);
#ifdef _64BIT_INDEX
			cpptr[(sideSz >> 3) - 4] = endianizeU<TIndexOffU>(occSave[0], this->toBe());
			cpptr[(sideSz >> 3) - 3] = endianizeU<TIndexOffU>(occSave[1], this->toBe());
			cpptr[(sideSz >> 3) - 2] = endianizeU<TIndexOffU>(occSave[2], this->toBe());
			cpptr[(sideSz >> 3) - 1] = endianizeU<TIndexOffU>(occSave[3], this->toBe());
#else
			cpptr[(sideSz >> 2) - 4] = endianizeU<TIndexOffU>(occSave[0], this->toBe());
			cpptr[(sideSz >> 2) - 3] = endianizeU<TIndexOffU>(occSave[1], this->toBe());
			cpptr[(sideSz >> 2) - 2] = endianizeU<TIndexOffU>(occSave[2], this->toBe());
			cpptr[(sideSz >> 2) - 1] = endianizeU<TIndexOffU>(occSave[3], this->toBe());
#endif
			occSave[0] = occ[0];
			occSave[1] = occ[1];
			occSave[2] = occ[2];
			occSave[3] = occ[3];
			out1.write((const char *)ebwtSide.ptr(), sideSz);
		}
	}
	const char *fn = _in6Str.c_str();
	savekit ss(fn);
	VMSG_NL("Exited Ebwt loop");
	assert_neq(zOff, OFF_MASK);
	if (absorbCnt > 0)
	{
		absorbFtab[ftabLen - 1] = absorbCnt;
	}
	assert_eq(side, eh._ebwtTotSz);
	assert_eq(((TIndexOffU)out1.tellp() - beforeEbwtOff), eh._ebwtTotSz);
	writeU<TIndexOffU>(out1, zOff, this->toBe());
	for (int i = 1; i < 4; i++)
	{
		fchr[i] += fchr[i - 1];
	}
	assert_eq(fchr[3], len);
	for (int i = 4; i >= 1; i--)
	{
		fchr[i] = fchr[i - 1];
	}
	fchr[0] = 0;
	if (_verbose)
	{
		for (int i = 0; i < 5; i++)
			cout << "fchr["
				 << "ACGT$"[i] << "]: " << fchr[i] << endl;
	}
	for (int i = 0; i < 5; i++)
	{
		writeU<TIndexOffU>(out1, fchr[i], this->toBe());
	}
	TIndexOffU eftabLen = 0;
	assert_eq(0, absorbFtab[0]);
	for (TIndexOffU i = 1; i < ftabLen; i++)
	{
		if (absorbFtab[i] > 0)
			eftabLen += 2;
	}
	assert_leq(eftabLen, (TIndexOffU)eh._ftabChars * 2);
	eftabLen = eh._ftabChars * 2;
	EList<TIndexOffU> eftab(EBWT_CAT);
	try
	{
		eftab.resize(eftabLen);
		eftab.fillZero();
	}
	catch (bad_alloc &e)
	{
		cerr << "Out of memory allocating eftab[] "
			 << "in Ebwt::buildToDisk() at " << __FILE__ << ":"
			 << __LINE__ << endl;
		throw e;
	}
	TIndexOffU eftabCur = 0;
	for (TIndexOffU i = 1; i < ftabLen; i++)
	{
		TIndexOffU lo = ftab[i] + Ebwt::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i - 1);
		if (absorbFtab[i] > 0)
		{
			TIndexOffU hi = lo + absorbFtab[i];
			assert_lt(eftabCur * 2 + 1, eftabLen);
			eftab[eftabCur * 2] = lo;
			eftab[eftabCur * 2 + 1] = hi;
			ftab[i] = (eftabCur++) ^ OFF_MASK;
			assert_eq(lo, Ebwt::ftabLo(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i));
			assert_eq(hi, Ebwt::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i));
		}
		else
		{
			ftab[i] = lo;
		}
	}
	assert_eq(Ebwt::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, ftabLen - 1), len + 1);
	for (TIndexOffU i = 0; i < ftabLen; i++)
	{
		writeU<TIndexOffU>(out1, ftab[i], this->toBe());
	}
	for (TIndexOffU i = 0; i < eftabLen; i++)
	{
		writeU<TIndexOffU>(out1, eftab[i], this->toBe());
	}
	assert(!isInMemory());
	VMSG_NL("Exiting Ebwt::buildToDisk()");
}
string adjustEbwtBase(const string &cmdline,
					  const string &ebwtFileBase,
					  bool verbose);
extern string gLastIOErrMsg;
inline bool is_read_err(int fdesc, ssize_t ret, size_t count)
{
	if (ret < 0)
	{
		std::stringstream sstm;
		sstm << "ERRNO: " << errno << " ERR Msg:" << strerror(errno) << std::endl;
		gLastIOErrMsg = sstm.str();
		return true;
	}
	return false;
}
inline bool is_fread_err(FILE *file_hd, size_t ret, size_t count)
{
	if (ferror(file_hd))
	{
		gLastIOErrMsg = "Error Reading File!";
		return true;
	}
	return false;
}
#endif
