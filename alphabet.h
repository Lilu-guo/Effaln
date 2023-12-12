#ifndef ALPHABETS_H_
#define ALPHABETS_H_
#include <stdexcept>
#include <string>
#include <sstream>
#include <stdint.h>
#include "assert_helpers.h"
using namespace std;
extern uint8_t asc2dnacat[];
extern char mask2dna[];
extern uint8_t asc2dnamask[];
extern int mask2popcnt[];
extern uint8_t asc2dna[];
extern uint8_t asc2col[];
extern uint8_t dinuc2color[5][5];
extern signed char mask2iupac[16];
extern int dnacomp[5];
extern const char *iupacs;
extern int maskcomp[16];
static inline bool isUnambigNuc(char c)
{
	return asc2dnacat[(int)c] == 1;
}
static inline char comp(char c)
{
	switch (c)
	{
	case 'a':
		return 't';
	case 'A':
		return 'T';
	case 'c':
		return 'g';
	case 'C':
		return 'G';
	case 'g':
		return 'c';
	case 'G':
		return 'C';
	case 't':
		return 'a';
	case 'T':
		return 'A';
	default:
		return c;
	}
}
static inline int compDna(int c)
{
	assert_leq(c, 4);
	return dnacomp[c];
}
extern uint8_t dinuc2color[5][5];
static inline void decodeNuc(char c, int &num, int *alts)
{
	switch (c)
	{
	case 'A':
		alts[0] = 0;
		num = 1;
		break;
	case 'C':
		alts[0] = 1;
		num = 1;
		break;
	case 'G':
		alts[0] = 2;
		num = 1;
		break;
	case 'T':
		alts[0] = 3;
		num = 1;
		break;
	case 'M':
		alts[0] = 0;
		alts[1] = 1;
		num = 2;
		break;
	case 'R':
		alts[0] = 0;
		alts[1] = 2;
		num = 2;
		break;
	case 'W':
		alts[0] = 0;
		alts[1] = 3;
		num = 2;
		break;
	case 'S':
		alts[0] = 1;
		alts[1] = 2;
		num = 2;
		break;
	case 'Y':
		alts[0] = 1;
		alts[1] = 3;
		num = 2;
		break;
	case 'K':
		alts[0] = 2;
		alts[1] = 3;
		num = 2;
		break;
	case 'V':
		alts[0] = 0;
		alts[1] = 1;
		alts[2] = 2;
		num = 3;
		break;
	case 'H':
		alts[0] = 0;
		alts[1] = 1;
		alts[2] = 3;
		num = 3;
		break;
	case 'D':
		alts[0] = 0;
		alts[1] = 2;
		alts[2] = 3;
		num = 3;
		break;
	case 'B':
		alts[0] = 1;
		alts[1] = 2;
		alts[2] = 3;
		num = 3;
		break;
	case 'N':
		alts[0] = 0;
		alts[1] = 1;
		alts[2] = 2;
		alts[3] = 3;
		num = 4;
		break;
	default:
	{
		std::cerr << "Bad IUPAC code: " << c << ", (int: " << (int)c << ")" << std::endl;
		throw std::runtime_error("");
	}
	}
}
extern void setIupacsCat(uint8_t cat);
#endif
