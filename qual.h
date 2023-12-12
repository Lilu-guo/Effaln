#ifndef QUAL_H_
#define QUAL_H_
#include <stdexcept>
#include "search_globals.h"
#include "sstring.h"
extern unsigned char qualRounds[];
extern unsigned char solToPhred[];
static inline uint8_t phredcToPhredq(char c)
{
	return ((uint8_t)c >= 33 ? ((uint8_t)c - 33) : 0);
}
static inline uint8_t solexaToPhred(int sol)
{
	assert_lt(sol, 256);
	if (sol < -10)
		return 0;
	return solToPhred[sol + 10];
}
class SimplePhredPenalty
{
public:
	static uint8_t mmPenalty(uint8_t qual)
	{
		return qual;
	}
	static uint8_t delPenalty(uint8_t qual)
	{
		return qual;
	}
	static uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right)
	{
		return std::max(qual_left, qual_right);
	}
};
class MaqPhredPenalty
{
public:
	static uint8_t mmPenalty(uint8_t qual)
	{
		return qualRounds[qual];
	}
	static uint8_t delPenalty(uint8_t qual)
	{
		return qualRounds[qual];
	}
	static uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right)
	{
		return qualRounds[std::max(qual_left, qual_right)];
	}
};
static inline uint8_t mmPenalty(bool maq, uint8_t qual)
{
	if (maq)
	{
		return MaqPhredPenalty::mmPenalty(qual);
	}
	else
	{
		return SimplePhredPenalty::mmPenalty(qual);
	}
}
static inline uint8_t delPenalty(bool maq, uint8_t qual)
{
	if (maq)
	{
		return MaqPhredPenalty::delPenalty(qual);
	}
	else
	{
		return SimplePhredPenalty::delPenalty(qual);
	}
}
static inline uint8_t insPenalty(bool maq, uint8_t qual_left, uint8_t qual_right)
{
	if (maq)
	{
		return MaqPhredPenalty::insPenalty(qual_left, qual_right);
	}
	else
	{
		return SimplePhredPenalty::insPenalty(qual_left, qual_right);
	}
}
inline static char charToPhred33(char c, bool solQuals, bool phred64Quals)
{
	using namespace std;
	if (c == ' ')
	{
		std::cerr << "Saw a space but expected an ASCII-encoded quality value." << endl
				  << "Are quality values formatted as integers?  If so, try --integer-quals." << endl;
		throw 1;
	}
	if (solQuals)
	{
		char cc = solexaToPhred((int)c - 64) + 33;
		if (cc < 33)
		{
			std::cerr << "Saw ASCII character "
					  << ((int)c)
					  << " but expected 64-based Solexa qual (converts to " << (int)cc << ")." << endl
					  << "Try not specifying --solexa-quals." << endl;
			throw 1;
		}
		c = cc;
	}
	else if (phred64Quals)
	{
		if (c < 64)
		{
			cerr << "Saw ASCII character "
				 << ((int)c)
				 << " but expected 64-based Phred qual." << endl
				 << "Try not specifying --solexa1.3-quals/--phred64-quals." << endl;
			throw 1;
		}
		c -= (64 - 33);
	}
	else
	{
		if (c < 33)
		{
			cerr << "Saw ASCII character "
				 << ((int)c)
				 << " but expected 33-based Phred qual." << endl;
			throw 1;
		}
	}
	return c;
}
inline static char intToPhred33(int iQ, bool solQuals)
{
	using namespace std;
	int pQ;
	if (solQuals)
	{
		pQ = solexaToPhred((int)iQ) + 33;
	}
	else
	{
		pQ = (iQ <= 93 ? iQ : 93) + 33;
	}
	if (pQ < 33)
	{
		cerr << "Saw negative Phred quality " << ((int)pQ - 33) << "." << endl;
		throw 1;
	}
	assert_geq(pQ, 0);
	return (int)pQ;
}
inline static uint8_t roundPenalty(uint8_t p)
{
	if (gNoMaqRound)
		return p;
	return qualRounds[p];
}
inline static uint8_t penaltiesAt(size_t off, uint8_t *q,
								  int alts,
								  const BTString &qual,
								  const BTDnaString *altQry,
								  const BTString *altQual)
{
	uint8_t primQ = qual[off];
	uint8_t bestPenalty = roundPenalty(phredcToPhredq(primQ));
	q[0] = q[1] = q[2] = q[3] = bestPenalty;
	for (int i = 0; i < alts; i++)
	{
		uint8_t altQ = altQual[i][off];
		if (altQ == 33)
			break;
		assert_leq(altQ, primQ);
		uint8_t pen = roundPenalty(primQ - altQ);
		if (pen < bestPenalty)
		{
			bestPenalty = pen;
		}
		int altC = (int)altQry[i][off];
		assert_lt(altC, 4);
		q[altC] = pen;
	}
	return bestPenalty;
}
inline static uint8_t loPenaltyAt(size_t off, int alts,
								  const BTString &qual,
								  const BTString *altQual)
{
	uint8_t primQ = qual[off];
	uint8_t bestPenalty = roundPenalty(phredcToPhredq(primQ));
	for (int i = 0; i < alts; i++)
	{
		uint8_t altQ = altQual[i][off];
		if (altQ == 33)
			break;
		assert_leq(altQ, primQ);
		uint8_t pen = roundPenalty(primQ - altQ);
		if (pen < bestPenalty)
		{
			bestPenalty = pen;
		}
	}
	return bestPenalty;
}
#endif
