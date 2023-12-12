#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <getopt.h>
#include "assert_helpers.h"
#include "endian_swap.h"
#include "aln_idx.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "timer.h"
#include "ref_read.h"
#include "filebuf.h"
#include "reference.h"
#include "ds.h"
int verbose;
static int sanityCheck;
static int format;
static TIndexOffU bmax;
static TIndexOffU bmaxMultSqrt;
static uint32_t bmaxDivN;
static int dcv;
static int noDc;
static int entireSA;
static int seed;
static int showVersion;
static int32_t lineRate;
static int32_t linesPerSide;
static int32_t offRate;
static int32_t ftabChars;
static int bigEndian;
static bool nsToAs;
static bool doSaFile;
static bool doBwtFile;
static bool autoMem;
static bool packed;
static bool writeRef;
static bool justRef;
static bool reverseEach;
static int nthreads;
static string wrapper;
static void resetOptions()
{
	verbose = true;
	sanityCheck = 0;
	format = FASTA;
	bmax = OFF_MASK;
	bmaxMultSqrt = OFF_MASK;
	bmaxDivN = 4;
	dcv = 1024;
	noDc = 0;
	entireSA = 0;
	seed = 0;
	showVersion = 0;
	lineRate = Ebwt::default_lineRate;
	linesPerSide = 1;
	offRate = 5;
	ftabChars = 10;
	bigEndian = 0;
	nsToAs = false;
	doSaFile = false;
	doBwtFile = false;
	autoMem = true;
	packed = false;
	writeRef = true;
	justRef = false;
	reverseEach = false;
	nthreads = 1;
	wrapper.clear();
}
enum
{
	ARG_BMAX = 256,
	ARG_BMAX_MULT,
	ARG_BMAX_DIV,
	ARG_DCV,
	ARG_SEED,
	ARG_CUTOFF,
	ARG_PMAP,
	ARG_NTOA,
	ARG_USAGE,
	ARG_REVERSE_EACH,
	ARG_SA,
	ARG_THREADS,
	ARG_WRAPPER
};
static const char *short_options = "qraph?nscfl:i:o:t:h:3C";
static struct option long_options[] = {
	{(char *)"quiet", no_argument, 0, 'q'},
	{(char *)"sanity", no_argument, 0, 's'},
	{(char *)"packed", no_argument, 0, 'p'},
	{(char *)"little", no_argument, &bigEndian, 0},
	{(char *)"big", no_argument, &bigEndian, 1},
	{(char *)"bmax", required_argument, 0, ARG_BMAX},
	{(char *)"bmaxmultsqrt", required_argument, 0, ARG_BMAX_MULT},
	{(char *)"bmaxdivn", required_argument, 0, ARG_BMAX_DIV},
	{(char *)"dcv", required_argument, 0, ARG_DCV},
	{(char *)"nodc", no_argument, &noDc, 1},
	{(char *)"seed", required_argument, 0, ARG_SEED},
	{(char *)"entiresa", no_argument, &entireSA, 1},
	{(char *)"version", no_argument, &showVersion, 1},
	{(char *)"noauto", no_argument, 0, 'a'},
	{(char *)"noblocks", required_argument, 0, 'n'},
	{(char *)"linerate", required_argument, 0, 'l'},
	{(char *)"linesperside", required_argument, 0, 'i'},
	{(char *)"offrate", required_argument, 0, 'o'},
	{(char *)"ftabchars", required_argument, 0, 't'},
	{(char *)"help", no_argument, 0, 'h'},
	{(char *)"ntoa", no_argument, 0, ARG_NTOA},
	{(char *)"justref", no_argument, 0, '3'},
	{(char *)"noref", no_argument, 0, 'r'},
	{(char *)"sa", no_argument, 0, ARG_SA},
	{(char *)"reverse-each", no_argument, 0, ARG_REVERSE_EACH},
	{(char *)"threads", required_argument, 0, ARG_THREADS},
	{(char *)"usage", no_argument, 0, ARG_USAGE},
	{(char *)"wrapper", required_argument, 0, ARG_WRAPPER},
	{(char *)0, 0, 0, 0}};
template <typename T>
static T parseNumber(T lower, const char *errmsg)
{
	char *endPtr = NULL;
	T t = (T)strtoll(optarg, &endPtr, 10);
	if (endPtr != NULL)
	{
		if (t < lower)
		{
			cerr << errmsg << endl;
			throw 1;
		}
		return t;
	}
	cerr << errmsg << endl;
	throw 1;
	return -1;
}
static bool parseOptions(int argc, const char **argv)
{
	int option_index = 0;
	int next_option;
	bool bmaxDivNSet = false;
	bool abort = false;
	do
	{
		next_option = getopt_long(
			argc, const_cast<char **>(argv),
			short_options, long_options, &option_index);
		switch (next_option)
		{
		case ARG_WRAPPER:
			wrapper = optarg;
			break;
		case 'f':
			format = FASTA;
			break;
		case 'c':
			format = CMDLINE;
			break;
		case 'p':
			packed = true;
			break;
		case 'l':
			lineRate = parseNumber<int>(3, "-l/--lineRate arg must be at least 3");
			break;
		case 'i':
			linesPerSide = parseNumber<int>(1, "-i/--linesPerSide arg must be at least 1");
			break;
		case 'o':
			offRate = parseNumber<int>(0, "-o/--offRate arg must be at least 0");
			break;
		case '3':
			justRef = true;
			break;
		case 't':
			ftabChars = parseNumber<int>(1, "-t/--ftabChars arg must be at least 1");
			if (ftabChars > 16)
			{
				std::cerr << "-t/--ftabChars arg must not exceed 16" << std::endl;
				throw 1;
			}
			break;
		case 'n':
			bmax = 0xfffffffe;
			break;
		case 'h':
		case ARG_USAGE:
			abort = true;
			break;
		case ARG_BMAX:
			bmax = parseNumber<TIndexOffU>(1, "--bmax arg must be at least 1");
			bmaxMultSqrt = OFF_MASK;
			bmaxDivN = 0xffffffff;
			break;
		case ARG_BMAX_MULT:
			bmaxMultSqrt = parseNumber<TIndexOffU>(1, "--bmaxmultsqrt arg must be at least 1");
			bmax = OFF_MASK;
			bmaxDivN = 0xffffffff;
			break;
		case ARG_BMAX_DIV:
			bmaxDivNSet = true;
			bmaxDivN = parseNumber<uint32_t>(1, "--bmaxdivn arg must be at least 1");
			bmax = OFF_MASK;
			bmaxMultSqrt = OFF_MASK;
			break;
		case ARG_DCV:
			dcv = parseNumber<int>(3, "--dcv arg must be at least 3");
			break;
		case ARG_SEED:
			seed = parseNumber<int>(0, "--seed arg must be at least 0");
			break;
		case ARG_REVERSE_EACH:
			reverseEach = true;
			break;
		case ARG_SA:
			doSaFile = true;
			break;
		case ARG_NTOA:
			nsToAs = true;
			break;
		case ARG_THREADS:
			nthreads = parseNumber<int>(0, "--threads arg must be at least 1");
			break;
		case 'a':
			autoMem = false;
			break;
		case 'q':
			verbose = false;
			break;
		case 's':
			sanityCheck = true;
			break;
		case 'r':
			writeRef = false;
			break;
		case -1:
			break;
		case 0:
			if (long_options[option_index].flag != 0)
				break;
		default:
			throw 1;
		}
	} while (next_option != -1);
	if (bmax < 40)
	{
		cerr << "Warning: specified bmax is very small (" << bmax << ").  This can lead to" << endl
			 << "extremely slow performance and memory exhaustion.  Perhaps you meant to specify" << endl
			 << "a small --bmaxdivn?" << endl;
	}
	if (!bmaxDivNSet)
	{
		bmaxDivN *= nthreads;
	}
	return abort;
}
EList<string> filesWritten;
static void deleteIdxFiles(
	const string &outfile,
	bool doRef,
	bool justRef)
{
	for (size_t i = 0; i < filesWritten.size(); i++)
	{
		cerr << "Deleting \"" << filesWritten[i].c_str()
			 << "\" file written during aborted indexing attempt." << endl;
		remove(filesWritten[i].c_str());
	}
}
template <typename TStr>
static void driver(
	const string &infile,
	EList<string> &infiles,
	const string &outfile,
	bool packed,
	int reverse)
{
	EList<FileBuf *> is(MISC_CAT);
	bool bisulfite = false;
	RefReadInParams refparams(false, reverse, nsToAs, bisulfite);
	assert_gt(infiles.size(), 0);
	if (format == CMDLINE)
	{
		stringstream *ss = new stringstream();
		for (size_t i = 0; i < infiles.size(); i++)
		{
			(*ss) << ">" << i << endl
				  << infiles[i].c_str() << endl;
		}
		FileBuf *fb = new FileBuf(ss);
		assert(fb != NULL);
		assert(!fb->eof());
		assert(fb->get() == '>');
		ASSERT_ONLY(fb->reset());
		assert(!fb->eof());
		is.push_back(fb);
	}
	else
	{
		for (size_t i = 0; i < infiles.size(); i++)
		{
			FileBuf *fb;
			size_t idx = infiles[i].find_last_of(".");
			std::string ext = (idx == std::string::npos) ? "" : infiles[i].substr(idx + 1);
			if (ext == "" || ext == "gz" || ext == "Z")
			{
				gzFile zFp = gzopen(infiles[i].c_str(), "rb");
				if (zFp == NULL)
				{
					cerr << "Error: could not open " << infiles[i].c_str() << endl;
					throw 1;
				}
				fb = new FileBuf(zFp);
			}
			else
			{
				FILE *f = fopen(infiles[i].c_str(), "rb");
				if (f == NULL)
				{
					cerr << "Error: could not open " << infiles[i].c_str() << endl;
					throw 1;
				}
				fb = new FileBuf(f);
			}
			assert(fb != NULL);
			if (fb->peek() == -1 || fb->eof())
			{
				cerr << "Warning: Empty fasta file: '" << infile.c_str() << "'" << endl;
				continue;
			}
			assert(!fb->eof());
			assert(fb->get() == '>');
			ASSERT_ONLY(fb->reset());
			assert(!fb->eof());
			is.push_back(fb);
		}
	}
	if (is.empty())
	{
		cerr << "Warning: All fasta inputs were empty" << endl;
		throw 1;
	}
	if (!reverse)
	{
#ifdef _64BIT_INDEX
		if (verbose)
			cerr << "Building a LARGE index" << endl;
#else
		if (verbose)
			cerr << "Building a SMALL index" << endl;
#endif
	}
	EList<RefRecord> szs(MISC_CAT);
	std::pair<size_t, size_t> sztot;
	{
		if (verbose)
			cout << "Reading reference sizes" << endl;
		Timer _t(cout, "  Time reading reference sizes: ", verbose);
		if (!reverse && (writeRef || justRef))
		{
			filesWritten.push_back(outfile + ".3." + gEbwt_ext);
			filesWritten.push_back(outfile + ".4." + gEbwt_ext);
			sztot = BitPairReference::szsFromFasta(is, outfile, bigEndian, refparams, szs, sanityCheck);
		}
		else
		{
			sztot = BitPairReference::szsFromFasta(is, string(), bigEndian, refparams, szs, sanityCheck);
		}
	}
	if (justRef)
		return;
	assert_gt(sztot.first, 0);
	assert_gt(sztot.second, 0);
	assert_gt(szs.size(), 0);
	filesWritten.push_back(outfile + ".1." + gEbwt_ext);
	filesWritten.push_back(outfile + ".2." + gEbwt_ext);
	Ebwt ebwt(
		TStr(),
		packed,
		0,
		1,
		lineRate,
		offRate,
		ftabChars,
		nthreads,
		outfile,
		reverse == 0,
		!entireSA,
		bmax,
		bmaxMultSqrt,
		bmaxDivN,
		noDc ? 0 : dcv,
		is,
		szs,
		(TIndexOffU)sztot.first,
		refparams,
		seed,
		-1,
		doSaFile,
		doBwtFile,
		verbose,
		autoMem,
		sanityCheck);
	if (verbose)
	{
		ebwt.eh().print(cout);
	}
	if (sanityCheck)
	{
		ebwt.loadIntoMemory(
			0,
			reverse ? (refparams.reverse == REF_READ_REVERSE) : 0,
			true,
			true,
			true,
			false,
			false);
		SString<char> s2;
		ebwt.restore(s2);
		ebwt.evictFromMemory();
		{
			SString<char> joinedss = Ebwt::join<SString<char>>(
				is,
				szs,
				(TIndexOffU)sztot.first,
				refparams,
				seed);
			if (refparams.reverse == REF_READ_REVERSE)
			{
				joinedss.reverse();
			}
			assert_eq(joinedss.length(), s2.length());
			assert(sstr_eq(joinedss, s2));
		}
		if (verbose)
		{
			if (s2.length() < 1000)
			{
				cout << "Passed restore check: " << s2.toZBuf() << endl;
			}
			else
			{
				cout << "Passed restore check: (" << s2.length() << " chars)" << endl;
			}
		}
	}
	for (size_t i = 0; i < is.size(); ++i)
	{
		if (is[i] != NULL)
			delete is[i];
	}
}
static const char *argv0 = NULL;
extern "C"
{
	int build_index(int argc, const char **argv)
	{
		string outfile;
		try
		{
			opterr = optind = 1;
			resetOptions();
			string infile;
			EList<string> infiles(MISC_CAT);
			if (parseOptions(argc, argv))
			{
				return 0;
			}
			argv0 = argv[0];
			if (showVersion)
			{
				cout << argv0 << " version 1.0" << endl;
				if (sizeof(void *) == 4)
				{
					cout << "32-bit" << endl;
				}
				else if (sizeof(void *) == 8)
				{
					cout << "64-bit" << endl;
				}
				else
				{
					cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void *) << endl;
				}
				cout << "Built on " << BUILD_HOST << endl;
				cout << BUILD_TIME << endl;
				cout << "Compiler: " << COMPILER_VERSION << endl;
				cout << "Options: " << COMPILER_OPTIONS << endl;
				cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
					 << sizeof(int)
					 << ", " << sizeof(long) << ", " << sizeof(long long)
					 << ", " << sizeof(void *) << ", " << sizeof(size_t)
					 << ", " << sizeof(off_t) << "}" << endl;
				return 0;
			}
			if (optind >= argc)
			{
				cerr << "No input sequence or sequence file specified!" << endl;
				return 1;
			}
			infile = argv[optind++];
			if (optind >= argc)
			{
				cerr << "No output file specified!" << endl;
				return 1;
			}
			outfile = argv[optind++];
			tokenize(infile, ",", infiles);
			if (infiles.size() < 1)
			{
				cerr << "Tokenized input file list was empty!" << endl;
				return 1;
			}
			if (verbose)
			{
				cout << "Settings:" << endl
					 << "  Output files: \"" << outfile.c_str() << ".*." + gEbwt_ext + "\"" << endl
					 << "  Line rate: " << lineRate << " (line is " << (1 << lineRate) << " bytes)" << endl
					 << "  Lines per side: " << linesPerSide << " (side is " << ((1 << lineRate) * linesPerSide) << " bytes)" << endl
					 << "  Offset rate: " << offRate << " (one in " << (1 << offRate) << ")" << endl
					 << "  FTable chars: " << ftabChars << endl
					 << "  Strings: " << (packed ? "packed" : "unpacked") << endl;
				if (bmax == OFF_MASK)
				{
					cout << "  Max bucket size: default" << endl;
				}
				else
				{
					cout << "  Max bucket size: " << bmax << endl;
				}
				if (bmaxMultSqrt == OFF_MASK)
				{
					cout << "  Max bucket size, sqrt multiplier: default" << endl;
				}
				else
				{
					cout << "  Max bucket size, sqrt multiplier: " << bmaxMultSqrt << endl;
				}
				if (bmaxDivN == 0xffffffff)
				{
					cout << "  Max bucket size, len divisor: default" << endl;
				}
				else
				{
					cout << "  Max bucket size, len divisor: " << bmaxDivN << endl;
				}
				cout << "  Difference-cover sample period: " << dcv << endl;
				cout << "  Endianness: " << (bigEndian ? "big" : "little") << endl
					 << "  Actual local endianness: " << (currentlyBigEndian() ? "big" : "little") << endl
					 << "  Sanity checking: " << (sanityCheck ? "enabled" : "disabled") << endl;
#ifdef NDEBUG
				cout << "  Assertions: disabled" << endl;
#else
				cout << "  Assertions: enabled" << endl;
#endif
				cout << "  Random seed: " << seed << endl;
				cout << "  Sizeofs: void*:" << sizeof(void *) << ", int:" << sizeof(int) << ", long:" << sizeof(long) << ", size_t:" << sizeof(size_t) << endl;
				cout << "Input files DNA, " << file_format_names[format].c_str() << ":" << endl;
				for (size_t i = 0; i < infiles.size(); i++)
				{
					cout << "  " << infiles[i].c_str() << endl;
				}
			}
			srand(seed);
			{
				Timer timer(cout, "Total time for call to driver() for forward index: ", verbose);
				if (!packed)
				{
					try
					{
						driver<SString<char>>(infile, infiles, outfile, false, REF_READ_FORWARD);
					}
					catch (bad_alloc &e)
					{
						if (autoMem)
						{
							cerr << "Switching to a packed string representation." << endl;
							packed = true;
						}
						else
						{
							throw e;
						}
					}
				}
				if (packed)
				{
					driver<S2bDnaString>(infile, infiles, outfile, true, REF_READ_FORWARD);
				}
			}
			int reverseType = reverseEach ? REF_READ_REVERSE_EACH : REF_READ_REVERSE;
			srand(seed);
			Timer timer(cout, "Total time for backward call to driver() for mirror index: ", verbose);
			if (!packed)
			{
				try
				{
					driver<SString<char>>(infile, infiles, outfile + ".rev", false, reverseType);
				}
				catch (bad_alloc &e)
				{
					if (autoMem)
					{
						cerr << "Switching to a packed string representation." << endl;
						packed = true;
					}
					else
					{
						throw e;
					}
				}
			}
			if (packed)
			{
				driver<S2bDnaString>(infile, infiles, outfile + ".rev", true, reverseType);
			}
			return 0;
		}
		catch (std::exception &e)
		{
			cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
			cerr << "Command: ";
			for (int i = 0; i < argc; i++)
				cerr << argv[i] << " ";
			cerr << endl;
			deleteIdxFiles(outfile, writeRef || justRef, justRef);
			return 1;
		}
		catch (int e)
		{
			if (e != 0)
			{
				cerr << "Error: Encountered internal effaln exception (#" << e << ")" << endl;
				cerr << "Command: ";
				for (int i = 0; i < argc; i++)
					cerr << argv[i] << " ";
				cerr << endl;
			}
			deleteIdxFiles(outfile, writeRef || justRef, justRef);
			return e;
		}
	}
}
