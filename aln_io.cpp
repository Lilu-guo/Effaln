#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "aln_idx.h"
#include <iomanip>
using namespace std;
void Ebwt::readIntoMemory(
	int color,
	int needEntireRev,
	bool loadSASamp,
	bool loadFtab,
	bool loadRstarts,
	bool justHeader,
	EbwtParams *params,
	bool mmSweep,
	bool loadNames,
	bool startVerbose)
{
	bool switchEndian;
#ifdef _MM
	char *mmFile[] = {NULL, NULL};
#endif
	if (_in1Str.length() > 0)
	{
		if (_verbose || startVerbose)
		{
			cerr << "  About to open input files: ";
			logTime(cerr);
		}
		if (_in1 != NULL)
			fclose(_in1);
		if (_verbose || startVerbose)
			cerr << "Opening \"" << _in1Str.c_str() << "\"" << endl;
		if ((_in1 = fopen(_in1Str.c_str(), "rb")) == NULL)
		{
			cerr << "Could not open index file " << _in1Str.c_str() << endl;
		}
		if (loadSASamp)
		{
			if (_in2 != NULL)
				fclose(_in2);
			if (_verbose || startVerbose)
				cerr << "Opening \"" << _in2Str.c_str() << "\"" << endl;
			if ((_in2 = fopen(_in2Str.c_str(), "rb")) == NULL)
			{
				cerr << "Could not open index file " << _in2Str.c_str() << endl;
			}
		}
		if (_verbose || startVerbose)
		{
			cerr << "  Finished opening input files: ";
			logTime(cerr);
		}
#ifdef _MM
		if (_useMm)
		{
			const char *names[] = {_in1Str.c_str(), _in2Str.c_str()};
			int fds[] = {fileno(_in1), fileno(_in2)};
			for (int i = 0; i < (loadSASamp ? 2 : 1); i++)
			{
				if (_verbose || startVerbose)
				{
					cerr << "  Memory-mapping input file " << (i + 1) << ": ";
					logTime(cerr);
				}
				struct stat sbuf;
				if (stat(names[i], &sbuf) == -1)
				{
					perror("stat");
					cerr << "Error: Could not stat index file " << names[i] << " prior to memory-mapping" << endl;
					throw 1;
				}
				mmFile[i] = (char *)mmap((void *)0, (size_t)sbuf.st_size,
										 PROT_READ, MAP_SHARED, fds[(size_t)i], 0);
				if (mmFile[i] == (void *)(-1))
				{
					perror("mmap");
					cerr << "Error: Could not memory-map the index file " << names[i] << endl;
					throw 1;
				}
				if (mmSweep)
				{
					int sum = 0;
					for (off_t j = 0; j < sbuf.st_size; j += 1024)
					{
						sum += (int)mmFile[i][j];
					}
					if (startVerbose)
					{
						cerr << "  Swept the memory-mapped ebwt index file 1; checksum: " << sum << ": ";
						logTime(cerr);
					}
				}
			}
			mmFile1_ = mmFile[0];
			mmFile2_ = loadSASamp ? mmFile[1] : NULL;
		}
#endif
	}
#ifdef _MM
	else if (_useMm && !justHeader)
	{
		mmFile[0] = mmFile1_;
		mmFile[1] = mmFile2_;
	}
	if (_useMm && !justHeader)
	{
		assert(mmFile[0] == mmFile1_);
		assert(mmFile[1] == mmFile2_);
	}
#endif
	if (_verbose || startVerbose)
	{
		cerr << "  Reading header: ";
		logTime(cerr);
	}
	uint64_t bytesRead = 0;
	switchEndian = false;
	uint32_t one = readU<uint32_t>(_in1, switchEndian);
	bytesRead += 4;
	if (loadSASamp)
	{
#ifndef NDEBUG
		assert_eq(one, readU<uint32_t>(_in2, switchEndian));
#else
		readU<uint32_t>(_in2, switchEndian);
#endif
	}
	if (one != 1)
	{
		assert_eq((1u << 24), one);
		assert_eq(1, endianSwapU32(one));
		switchEndian = true;
	}
	if (switchEndian && _useMm)
	{
		cerr << "Error: Can't use memory-mapped files when the index is the opposite endianness" << endl;
		throw 1;
	}
	TIndexOffU len = readU<TIndexOffU>(_in1, switchEndian);
	bytesRead += OFF_SIZE;
	int32_t lineRate = readI<int32_t>(_in1, switchEndian);
	bytesRead += 4;
	readI<int32_t>(_in1, switchEndian);
	bytesRead += 4;
	int32_t offRate = readI<int32_t>(_in1, switchEndian);
	bytesRead += 4;
	int32_t ftabChars = readI<int32_t>(_in1, switchEndian);
	bytesRead += 4;
	int32_t flags = readI<int32_t>(_in1, switchEndian);
	bool entireRev = false;
	if (flags < 0 && (((-flags) & EBWT_COLOR) != 0))
	{
		if (color != -1 && !color)
		{
			cerr << "Error: -C was not specified." << endl;
			throw 1;
		}
		color = 1;
	}
	else if (flags < 0)
	{
		if (color != -1 && color)
		{
			cerr << "Error: -C was specified." << endl;
			throw 1;
		}
		color = 0;
	}
	if (flags < 0 && (((-flags) & EBWT_ENTIRE_REV) == 0))
	{
		if (needEntireRev != -1 && needEntireRev != 0)
		{
			cerr << "Error: This index is not compatible." << endl;
			throw 1;
		}
	}
	else
		entireRev = true;
	bytesRead += 4;
	EbwtParams *eh;
	bool deleteEh = false;
	if (params != NULL)
	{
		params->init(len, lineRate, offRate, ftabChars, color, entireRev);
		if (_verbose || startVerbose)
			params->print(cerr);
		eh = params;
	}
	else
	{
		eh = new EbwtParams(len, lineRate, offRate, ftabChars, color, entireRev);
		deleteEh = true;
	}
	TIndexOffU offsLen = eh->_offsLen;
	uint64_t offsSz = eh->_offsSz;
	TIndexOffU offRateDiff = 0;
	TIndexOffU offsLenSampled = offsLen;
	if (_overrideOffRate > offRate)
	{
		offRateDiff = _overrideOffRate - offRate;
	}
	if (offRateDiff > 0)
	{
		offsLenSampled >>= offRateDiff;
		if ((offsLen & ~(OFF_MASK << offRateDiff)) != 0)
		{
			offsLenSampled++;
		}
	}
	if (_useMm && (offRateDiff))
	{
		cerr << "Error: Can't use memory-mapped files when the offrate is overridden" << endl;
		throw 1;
	}
	this->_nPat = readI<TIndexOffU>(_in1, switchEndian);
	bytesRead += OFF_SIZE;
	_plen.reset();
	if (_useMm)
	{
#ifdef _MM
		_plen.init((TIndexOffU *)(mmFile[0] + bytesRead), _nPat, false);
		bytesRead += _nPat * OFF_SIZE;
		fseeko(_in1, _nPat * OFF_SIZE, SEEK_CUR);
#endif
	}
	else
	{
		try
		{
			if (_verbose || startVerbose)
			{
				cerr << "Reading plen (" << this->_nPat << "): ";
				logTime(cerr);
			}
			_plen.init(new TIndexOffU[_nPat], _nPat, true);
			if (switchEndian)
			{
				for (TIndexOffU i = 0; i < this->_nPat; i++)
				{
					plen()[i] = readU<TIndexOffU>(_in1, switchEndian);
				}
			}
			else
			{
				size_t r = MM_READ(_in1, (void *)(plen()), _nPat * OFF_SIZE);
				if (r != (size_t)(_nPat * OFF_SIZE))
				{
					cerr << "Error reading _plen[] array: " << r << ", " << _nPat * OFF_SIZE << endl;
					throw 1;
				}
			}
		}
		catch (bad_alloc &e)
		{
			cerr << "Out of memory allocating plen[] in Ebwt::read()"
				 << " at " << __FILE__ << ":" << __LINE__ << endl;
			throw e;
		}
	}
	bool shmemLeader;
	if (justHeader)
		goto done;
	this->_nFrag = readU<TIndexOffU>(_in1, switchEndian);
	bytesRead += OFF_SIZE;
	if (_verbose || startVerbose)
	{
		cerr << "Reading rstarts (" << this->_nFrag * 3 << "): ";
		logTime(cerr);
	}
	assert_geq(this->_nFrag, this->_nPat);
	_rstarts.reset();
	if (loadRstarts)
	{
		if (_useMm)
		{
#ifdef _MM
			_rstarts.init((TIndexOffU *)(mmFile[0] + bytesRead), _nFrag * 3, false);
			bytesRead += this->_nFrag * OFF_SIZE * 3;
			fseeko(_in1, this->_nFrag * OFF_SIZE * 3, SEEK_CUR);
#endif
		}
		else
		{
			_rstarts.init(new TIndexOffU[_nFrag * 3], _nFrag * 3, true);
			if (switchEndian)
			{
				for (TIndexOffU i = 0; i < this->_nFrag * 3; i += 3)
				{
					this->rstarts()[i] = readU<TIndexOffU>(_in1, switchEndian);
					this->rstarts()[i + 1] = readU<TIndexOffU>(_in1, switchEndian);
					this->rstarts()[i + 2] = readU<TIndexOffU>(_in1, switchEndian);
				}
			}
			else
			{
				size_t r = MM_READ(_in1, (void *)rstarts(), this->_nFrag * OFF_SIZE * 3);
				if (r != (size_t)(this->_nFrag * OFF_SIZE * 3))
				{
					cerr << "Error reading _rstarts[] array: " << r << ", " << (this->_nFrag * OFF_SIZE * 3) << endl;
					throw 1;
				}
			}
		}
	}
	else
	{
		assert(rstarts() == NULL);
		bytesRead += this->_nFrag * OFF_SIZE * 3;
		fseeko(_in1, this->_nFrag * OFF_SIZE * 3, SEEK_CUR);
	}
	_ebwt.reset();
	if (_useMm)
	{
#ifdef _MM
		_ebwt.init((uint8_t *)(mmFile[0] + bytesRead), eh->_ebwtTotLen, false);
		bytesRead += eh->_ebwtTotLen;
		fseek(_in1, eh->_ebwtTotLen, SEEK_CUR);
#endif
	}
	else
	{
		if (_verbose || startVerbose)
		{
			cerr << "Reading ebwt (" << eh->_ebwtTotLen << "): ";
			logTime(cerr);
		}
		bool shmemLeader = true;
		if (useShmem_)
		{
			uint8_t *tmp = NULL;
			shmemLeader = ALLOC_SHARED_U8(
				(_in1Str + "[ebwt]"), eh->_ebwtTotLen, &tmp,
				"ebwt[]", (_verbose || startVerbose));
			assert(tmp != NULL);
			_ebwt.init(tmp, eh->_ebwtTotLen, false);
			if (_verbose || startVerbose)
			{
				cerr << "  shared-mem " << (shmemLeader ? "leader" : "follower") << endl;
			}
		}
		else
		{
			try
			{
				_ebwt.init(new uint8_t[eh->_ebwtTotLen], eh->_ebwtTotLen, true);
			}
			catch (bad_alloc &e)
			{
				cerr << "Out of memory allocating the ebwt[] array for the index." << endl;
				throw 1;
			}
		}
		if (shmemLeader)
		{
			uint64_t bytesLeft = eh->_ebwtTotLen;
			char *pebwt = (char *)this->ebwt();
			while (bytesLeft > 0)
			{
				size_t r = MM_READ(_in1, (void *)pebwt, bytesLeft);
				if (MM_IS_IO_ERR(_in1, r, bytesLeft))
				{
					cerr << "Error reading _ebwt[] array: " << r << ", "
						 << bytesLeft << gLastIOErrMsg << endl;
					throw 1;
				}
				else if (r == 0)
				{
					cerr << "Error reading _ebwt[] array: no more data" << endl;
					throw 1;
				}
				pebwt += r;
				bytesLeft -= r;
			}
			if (switchEndian)
			{
				uint8_t *side = this->ebwt();
				for (size_t i = 0; i < eh->_numSides; i++)
				{
					TIndexOffU *cums = reinterpret_cast<TIndexOffU *>(side + eh->_sideSz - OFF_SIZE * 2);
					cums[0] = endianSwapU(cums[0]);
					cums[1] = endianSwapU(cums[1]);
					side += this->_eh._sideSz;
				}
			}
#ifdef _SHARED_MEM
			if (useShmem_)
				NOTIFY_SHARED(ebwt(), eh->_ebwtTotLen);
#endif
		}
		else
		{
			fseeko(_in1, eh->_ebwtTotLen, SEEK_CUR);
#ifdef _SHARED_MEM
			if (useShmem_)
				WAIT_SHARED(ebwt(), eh->_ebwtTotLen);
#endif
		}
	}
	_zOff = readU<TIndexOffU>(_in1, switchEndian);
	bytesRead += OFF_SIZE;
	assert_lt(_zOff, len);
	try
	{
		if (_verbose || startVerbose)
			cerr << "Reading fchr (5)" << endl;
		_fchr.reset();
		if (_useMm)
		{
#ifdef _MM
			_fchr.init((TIndexOffU *)(mmFile[0] + bytesRead), 5, false);
			bytesRead += 5 * OFF_SIZE;
			fseek(_in1, 5 * OFF_SIZE, SEEK_CUR);
#endif
		}
		else
		{
			_fchr.init(new TIndexOffU[5], 5, true);
			for (int i = 0; i < 5; i++)
			{
				this->fchr()[i] = readU<TIndexOffU>(_in1, switchEndian);
				assert_leq(this->fchr()[i], len);
				assert(i <= 0 || this->fchr()[i] >= this->fchr()[i - 1]);
			}
		}
		assert_gt(this->fchr()[4], this->fchr()[0]);
		if (_verbose || startVerbose)
		{
			if (loadFtab)
			{
				cerr << "Reading ftab (" << eh->_ftabLen << "): ";
				logTime(cerr);
			}
			else
			{
				cerr << "Skipping ftab (" << eh->_ftabLen << "): ";
			}
		}
		_ftab.reset();
		if (loadFtab)
		{
			if (_useMm)
			{
#ifdef _MM
				_ftab.init((TIndexOffU *)(mmFile[0] + bytesRead), eh->_ftabLen, false);
				bytesRead += eh->_ftabLen * OFF_SIZE;
				fseeko(_in1, eh->_ftabLen * OFF_SIZE, SEEK_CUR);
#endif
			}
			else
			{
				_ftab.init(new TIndexOffU[eh->_ftabLen], eh->_ftabLen, true);
				if (switchEndian)
				{
					for (TIndexOffU i = 0; i < eh->_ftabLen; i++)
						this->ftab()[i] = readU<TIndexOffU>(_in1, switchEndian);
				}
				else
				{
					size_t r = MM_READ(_in1, (void *)ftab(), eh->_ftabLen * OFF_SIZE);
					if (r != (size_t)(eh->_ftabLen * OFF_SIZE))
					{
						cerr << "Error reading _ftab[] array: " << r << ", " << (eh->_ftabLen * OFF_SIZE) << endl;
						throw 1;
					}
				}
			}
			if (_verbose || startVerbose)
			{
				if (loadFtab)
				{
					cerr << "Reading eftab (" << eh->_eftabLen << "): ";
					logTime(cerr);
				}
				else
				{
					cerr << "Skipping eftab (" << eh->_eftabLen << "): ";
				}
			}
			_eftab.reset();
			if (_useMm)
			{
#ifdef _MM
				_eftab.init((TIndexOffU *)(mmFile[0] + bytesRead), eh->_eftabLen, false);
				bytesRead += eh->_eftabLen * OFF_SIZE;
				fseeko(_in1, eh->_eftabLen * OFF_SIZE, SEEK_CUR);
#endif
			}
			else
			{
				_eftab.init(new TIndexOffU[eh->_eftabLen], eh->_eftabLen, true);
				if (switchEndian)
				{
					for (TIndexOffU i = 0; i < eh->_eftabLen; i++)
						this->eftab()[i] = readU<TIndexOffU>(_in1, switchEndian);
				}
				else
				{
					size_t r = MM_READ(_in1, (void *)this->eftab(), eh->_eftabLen * OFF_SIZE);
					if (r != (size_t)(eh->_eftabLen * OFF_SIZE))
					{
						cerr << "Error reading _eftab[] array: " << r << ", " << (eh->_eftabLen * OFF_SIZE) << endl;
						throw 1;
					}
				}
			}
			for (TIndexOffU i = 0; i < eh->_eftabLen; i++)
			{
				if (i > 0 && this->eftab()[i] > 0)
				{
					assert_geq(this->eftab()[i], this->eftab()[i - 1]);
				}
				else if (i > 0 && this->eftab()[i - 1] == 0)
				{
					assert_eq(0, this->eftab()[i]);
				}
			}
		}
		else
		{
			assert(ftab() == NULL);
			assert(eftab() == NULL);
			bytesRead += eh->_ftabLen * OFF_SIZE;
			fseeko(_in1, eh->_ftabLen * OFF_SIZE, SEEK_CUR);
			bytesRead += eh->_eftabLen * OFF_SIZE;
			fseeko(_in1, eh->_eftabLen * OFF_SIZE, SEEK_CUR);
		}
	}
	catch (bad_alloc &e)
	{
		cerr << "Out of memory allocating fchr[], ftab[] or eftab[] arrays for the index." << endl;
		throw 1;
	}
	if (loadNames)
	{
		while (true)
		{
			char c = '\0';
			if (MM_READ(_in1, (void *)(&c), (size_t)1) != (size_t)1)
				break;
			bytesRead++;
			if (c == '\0')
				break;
			else if (c == '\n')
			{
				this->_refnames.push_back("");
			}
			else
			{
				if (this->_refnames.size() == 0)
				{
					this->_refnames.push_back("");
				}
				this->_refnames.back().push_back(c);
			}
		}
	}
	_offs.reset();
	if (loadSASamp)
	{
		bytesRead = 4;
		shmemLeader = true;
		if (_verbose || startVerbose)
		{
			cerr << "Reading offs (" << offsLenSampled << std::setw(2) << OFF_SIZE * 8 << "-bit words): ";
			logTime(cerr);
		}
		if (!_useMm)
		{
			if (!useShmem_)
			{
				try
				{
					_offs.init(new TIndexOffU[offsLenSampled], offsLenSampled, true);
				}
				catch (bad_alloc &e)
				{
					cerr << "Out of memory allocating the offs[] array  for the index." << endl;
					throw 1;
				}
			}
			else
			{
				TIndexOffU *tmp = NULL;
				shmemLeader = ALLOC_SHARED_U(
					(_in2Str + "[offs]"), offsLenSampled * OFF_SIZE, &tmp,
					"offs", (_verbose || startVerbose));
				_offs.init((TIndexOffU *)tmp, offsLenSampled, false);
			}
		}
		if (_overrideOffRate < 32)
		{
			if (shmemLeader)
			{
				if (switchEndian || offRateDiff > 0)
				{
					assert(!_useMm);
					const TIndexOffU blockMaxSz = (2 * 1024 * 1024);
					const TIndexOffU blockMaxSzU = (blockMaxSz >> (OFF_SIZE / 4 + 1));
					char *buf;
					try
					{
						buf = new char[blockMaxSz];
					}
					catch (std::bad_alloc &e)
					{
						cerr << "Error: Out of memory allocating part of _offs array: '" << e.what() << "'" << endl;
						throw e;
					}
					for (TIndexOffU i = 0; i < offsLen; i += blockMaxSzU)
					{
						TIndexOffU block = min<TIndexOffU>(blockMaxSzU, offsLen - i);
						size_t r = MM_READ(_in2, (void *)buf, block << (OFF_SIZE / 4 + 1));
						if (r != (size_t)(block << (OFF_SIZE / 4 + 1)))
						{
							cerr << "Error reading block of _offs[] array: " << r << ", " << (block << (OFF_SIZE / 4 + 1)) << endl;
							throw 1;
						}
						TIndexOffU idx = i >> offRateDiff;
						for (TIndexOffU j = 0; j < block; j += (1 << offRateDiff))
						{
							assert_lt(idx, offsLenSampled);
							this->offs()[idx] = ((TIndexOffU *)buf)[j];
							if (switchEndian)
							{
								this->offs()[idx] = endianSwapU(this->offs()[idx]);
							}
							idx++;
						}
					}
					delete[] buf;
				}
				else
				{
					if (_useMm)
					{
#ifdef _MM
						_offs.init((TIndexOffU *)(mmFile[1] + bytesRead), offsLen, false);
						bytesRead += offsSz;
						fseeko(_in2, offsSz, SEEK_CUR);
#endif
					}
					else
					{
						uint64_t bytesLeft = offsSz;
						char *offs = (char *)this->offs();
						while (bytesLeft > 0)
						{
							size_t r = MM_READ(_in2, (void *)offs, bytesLeft);
							if (MM_IS_IO_ERR(_in2, r, bytesLeft))
							{
								cerr << "Error reading block of _offs[] array: "
									 << r << ", " << bytesLeft << gLastIOErrMsg << endl;
								throw 1;
							}
							else if (r == 0)
							{
								cerr << "Error reading block of _offs[] array: no more data" << endl;
								throw 1;
							}
							offs += r;
							bytesLeft -= r;
						}
					}
				}
#ifdef _SHARED_MEM
				if (useShmem_)
					NOTIFY_SHARED(offs(), offsLenSampled * OFF_SIZE);
#endif
			}
			else
			{
				fseeko(_in2, offsLenSampled * OFF_SIZE, SEEK_CUR);
#ifdef _SHARED_MEM
				if (useShmem_)
					WAIT_SHARED(offs(), offsLenSampled * OFF_SIZE);
#endif
			}
		}
	}
	this->postReadInit(*eh);
	if (_verbose || startVerbose)
		print(cerr, *eh);
done:
	if (deleteEh)
		delete eh;
	if (_in1 != NULL)
	{
		rewind(_in1);
	}
	if (_in2 != NULL)
	{
		rewind(_in2);
	}
}
void readEbwtRefnames(FILE *fin, EList<string> &refnames)
{
	assert(fin != NULL);
	assert_eq(ftello(fin), 0);
	bool switchEndian = false;
	uint32_t one = readU<uint32_t>(fin, switchEndian);
	if (one != 1)
	{
		assert_eq((1u << 24), one);
		switchEndian = true;
	}
	TIndexOffU len = readU<TIndexOffU>(fin, switchEndian);
	int32_t lineRate = readI<int32_t>(fin, switchEndian);
	readI<int32_t>(fin, switchEndian);
	int32_t offRate = readI<int32_t>(fin, switchEndian);
	int32_t ftabChars = readI<int32_t>(fin, switchEndian);
	int32_t flags = readI<int32_t>(fin, switchEndian);
	bool color = false;
	bool entireReverse = false;
	if (flags < 0)
	{
		color = (((-flags) & EBWT_COLOR) != 0);
		entireReverse = (((-flags) & EBWT_ENTIRE_REV) != 0);
	}
	EbwtParams eh(len, lineRate, offRate, ftabChars, color, entireReverse);
	TIndexOffU nPat = readI<TIndexOffU>(fin, switchEndian);
	fseeko(fin, nPat * OFF_SIZE, SEEK_CUR);
	TIndexOffU nFrag = readU<TIndexOffU>(fin, switchEndian);
	fseeko(fin, nFrag * OFF_SIZE * 3, SEEK_CUR);
	fseeko(fin, eh._ebwtTotLen, SEEK_CUR);
	readU<TIndexOffU>(fin, switchEndian);
	fseeko(fin, 5 * OFF_SIZE, SEEK_CUR);
	fseeko(fin, eh._ftabLen * OFF_SIZE, SEEK_CUR);
	fseeko(fin, eh._eftabLen * OFF_SIZE, SEEK_CUR);
	while (true)
	{
		char c = '\0';
		int read_value = 0;
		read_value = fgetc(fin);
		if (read_value == EOF)
			break;
		c = read_value;
		if (c == '\0')
			break;
		else if (c == '\n')
		{
			refnames.push_back("");
		}
		else
		{
			if (refnames.size() == 0)
			{
				refnames.push_back("");
			}
			refnames.back().push_back(c);
		}
	}
	if (refnames.back().empty())
	{
		refnames.pop_back();
	}
	fseeko(fin, 0, SEEK_SET);
	assert(ferror(fin) == 0);
}
void readEbwtRefnames(const string &instr, EList<string> &refnames)
{
	FILE *fin;
	fin = fopen((instr + ".1." + gEbwt_ext).c_str(), "rb");
	if (fin == NULL)
	{
		throw EbwtFileOpenException("Cannot open file " + instr);
	}
	assert_eq(ftello(fin), 0);
	readEbwtRefnames(fin, refnames);
	fclose(fin);
}
int32_t Ebwt::readFlags(const string &instr)
{
	ifstream in;
	in.open((instr + ".1." + gEbwt_ext).c_str(), ios_base::in | ios::binary);
	if (!in.is_open())
	{
		throw EbwtFileOpenException("Cannot open file " + instr);
	}
	assert(in.is_open());
	assert(in.good());
	bool switchEndian = false;
	uint32_t one = readU<uint32_t>(in, switchEndian);
	if (one != 1)
	{
		assert_eq((1u << 24), one);
		assert_eq(1, endianSwapU32(one));
		switchEndian = true;
	}
	readU<TIndexOffU>(in, switchEndian);
	readI<int32_t>(in, switchEndian);
	readI<int32_t>(in, switchEndian);
	readI<int32_t>(in, switchEndian);
	readI<int32_t>(in, switchEndian);
	int32_t flags = readI<int32_t>(in, switchEndian);
	return flags;
}
bool readEbwtColor(const string &instr)
{
	int32_t flags = Ebwt::readFlags(instr);
	if (flags < 0 && (((-flags) & EBWT_COLOR) != 0))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool readEntireReverse(const string &instr)
{
	int32_t flags = Ebwt::readFlags(instr);
	if (flags < 0 && (((-flags) & EBWT_ENTIRE_REV) != 0))
	{
		return true;
	}
	else
	{
		return false;
	}
}
void Ebwt::writeFromMemory(bool justHeader,
						   ostream &out1,
						   ostream &out2) const
{
	const EbwtParams &eh = this->_eh;
	assert(eh.repOk());
	uint32_t be = this->toBe();
	assert(out1.good());
	assert(out2.good());
	writeI<int32_t>(out1, 1, be);
	writeI<int32_t>(out2, 1, be);
	writeU<TIndexOffU>(out1, eh._len, be);
	writeI<int32_t>(out1, eh._lineRate, be);
	writeI<int32_t>(out1, 2, be);
	writeI<int32_t>(out1, eh._offRate, be);
	writeI<int32_t>(out1, eh._ftabChars, be);
	int32_t flags = 1;
	if (eh._color)
		flags |= EBWT_COLOR;
	if (eh._entireReverse)
		flags |= EBWT_ENTIRE_REV;
	writeI<int32_t>(out1, -flags, be);
	if (!justHeader)
	{
		assert(rstarts() != NULL);
		assert(offs() != NULL);
		assert(ftab() != NULL);
		assert(eftab() != NULL);
		assert(isInMemory());
		writeU<TIndexOffU>(out1, this->_nPat, be);
		for (TIndexOffU i = 0; i < this->_nPat; i++)
			writeU<TIndexOffU>(out1, this->plen()[i], be);
		assert_geq(this->_nFrag, this->_nPat);
		writeU<TIndexOffU>(out1, this->_nFrag, be);
		for (TIndexOffU i = 0; i < this->_nFrag * 3; i++)
			writeU<TIndexOffU>(out1, this->rstarts()[i], be);
		out1.write((const char *)this->ebwt(), eh._ebwtTotLen);
		writeU<TIndexOffU>(out1, this->zOff(), be);
		TIndexOffU offsLen = eh._offsLen;
		for (TIndexOffU i = 0; i < offsLen; i++)
			writeU<TIndexOffU>(out2, this->offs()[i], be);
		for (int i = 0; i < 5; i++)
			writeU<TIndexOffU>(out1, this->fchr()[i], be);
		for (TIndexOffU i = 0; i < eh._ftabLen; i++)
			writeU<TIndexOffU>(out1, this->ftab()[i], be);
		for (TIndexOffU i = 0; i < eh._eftabLen; i++)
			writeU<TIndexOffU>(out1, this->eftab()[i], be);
	}
}
void Ebwt::writeFromMemory(bool justHeader,
						   const string &out1,
						   const string &out2) const
{
	ASSERT_ONLY(const EbwtParams &eh = this->_eh);
	assert(isInMemory());
	assert(eh.repOk());
	ofstream fout1(out1.c_str(), ios::binary);
	ofstream fout2(out2.c_str(), ios::binary);
	writeFromMemory(justHeader, fout1, fout2);
	fout1.close();
	fout2.close();
	if (_sanity)
	{
#if 0
		if(_verbose)
			cout << "Re-reading \"" << out1 << "\"/\"" << out2 << "\" for sanity check" << endl;
		Ebwt copy(out1, out2, _verbose, _sanity);
		assert(!isInMemory());
		copy.loadIntoMemory(eh._color ? 1 : 0, true, false, false);
		assert(isInMemory());
	    assert_eq(eh._lineRate,     copy.eh()._lineRate);
	    assert_eq(eh._offRate,      copy.eh()._offRate);
	    assert_eq(eh._ftabChars,    copy.eh()._ftabChars);
	    assert_eq(eh._len,          copy.eh()._len);
	    assert_eq(_zOff,             copy.zOff());
	    assert_eq(_zEbwtBpOff,       copy.zEbwtBpOff());
	    assert_eq(_zEbwtByteOff,     copy.zEbwtByteOff());
		assert_eq(_nPat,             copy.nPat());
		for(TIndexOffU i = 0; i < _nPat; i++)
			assert_eq(this->_plen[i], copy.plen()[i]);
		assert_eq(this->_nFrag, copy.nFrag());
		for(TIndexOffU i = 0; i < this->nFrag*3; i++) {
			assert_eq(this->_rstarts[i], copy.rstarts()[i]);
		}
		for(int i = 0; i < 5; i++)
			assert_eq(this->_fchr[i], copy.fchr()[i]);
		for(TIndexOffU i = 0; i < eh._ftabLen; i++)
			assert_eq(this->ftab()[i], copy.ftab()[i]);
		for(TIndexOffU i = 0; i < eh._eftabLen; i++)
			assert_eq(this->eftab()[i], copy.eftab()[i]);
		for(TIndexOffU i = 0; i < eh._offsLen; i++)
			assert_eq(this->_offs[i], copy.offs()[i]);
		for(TIndexOffU i = 0; i < eh._ebwtTotLen; i++)
			assert_eq(this->ebwt()[i], copy.ebwt()[i]);
		copy.sanityCheckAll();
		if(_verbose)
			cout << "Read-in check passed for \"" << out1 << "\"/\"" << out2 << "\"" << endl;
#endif
	}
}
void Ebwt::szsToDisk(const EList<RefRecord> &szs, ostream &os, int reverse)
{
	TIndexOffU seq = 0;
	TIndexOffU off = 0;
	TIndexOffU totlen = 0;
	for (unsigned int i = 0; i < szs.size(); i++)
	{
		if (szs[i].len == 0)
			continue;
		if (szs[i].first)
			off = 0;
		off += szs[i].off;
		if (szs[i].first && szs[i].len > 0)
			seq++;
		TIndexOffU seqm1 = seq - 1;
		assert_lt(seqm1, _nPat);
		TIndexOffU fwoff = off;
		if (reverse == REF_READ_REVERSE)
		{
			seqm1 = _nPat - seqm1 - 1;
			assert_leq(off + szs[i].len, plen()[seqm1]);
			fwoff = plen()[seqm1] - (off + szs[i].len);
		}
		writeU<TIndexOffU>(os, totlen, this->toBe());
		writeU<TIndexOffU>(os, seqm1, this->toBe());
		writeU<TIndexOffU>(os, fwoff, this->toBe());
		totlen += szs[i].len;
		off += szs[i].len;
	}
}
