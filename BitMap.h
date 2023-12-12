#ifndef WT_NODE_H
#define WT_NODE_H
#include <string.h>
#include "loadkit.h"
#include "savekit.h"
#include "InArray.h"
#include "BaseClass.h"
#include <math.h>
#include <iostream>
#include <bitset>
using namespace std;
class Dmap
{
public:
	Dmap();
	u64 GetMemorySize()
	{
		return SBd->GetMemorySize() + Bd->GetMemorySize() + datanum * datawidth / 8;
	};
	Dmap(u64 data_num, u32 data_width);
	Dmap(u64 len, u32 size, int Blength, int SBlength);
	void write(savekit &s)
	{
		s.writeu64(datanum);
		s.writeu32(datawidth);
		s.writei32(Blength);
		s.writei32(SBlength);
		SBd->write(s);
		Bd->write(s);
	};
	void load(loadkit &s)
	{
		s.loadu64(this->datanum);
		s.loadu32(this->datawidth);
		s.loadi32(this->Blength);
		s.loadi32(this->SBlength);
		this->SBd = new InArray();
		this->Bd = new InArray();
		this->SBd->load(s);
		this->Bd->load(s);
	};
	~Dmap(void);
	void SetValue(u64 index, u64 v);
	void constructrank(int Blength, int SBlength, int offrate = 5);
	u64 rank(u64 i)
	{
		u64 Mdnum = i >> Blength;
		u64 SBnum = Mdnum >> SBlength;
		u64 Sbdbefore = SBd->GetValue(SBnum);
		u64 Sdnowi = Sbdbefore + Bd->GetValue(Mdnum);
		u32 popc_n = (i - (i & 0xffffff00)) >> 5;
		u32 temp = i & 0xffffffe0;
		u32 left_s = i - temp;
		u32 d_s = Mdnum * (1 << (Blength - 5));
		for (u32 j = 0; j < popc_n; j++)
			Sdnowi += (u32)__builtin_popcount(data[d_s + j]);
		if (left_s != 0)
		{
			u32 left1 = data[d_s + popc_n] & (0xffffffff << (32 - left_s));
			Sdnowi += (u32)__builtin_popcount(left1);
		}
		return Sdnowi + 1;
	}
	int Save(savekit &s)
	{
		s.writei32(Blength);
		s.writei32(SBlength);
		SBd->write(s);
		Bd->write(s);
		s.writeu64(datanum);
		s.writeu32(datawidth);
		u64 len = (datanum * datawidth);
		if (len % 32 == 0)
			len = len / 32 + 1;
		else
			len = len / 32 + 2;
		s.writeu64(len);
		s.writeu32array(data, len);
		return 0;
	}
	int Load(loadkit &s)
	{
		s.loadi32(this->Blength);
		s.loadi32(this->SBlength);
		SBd = new InArray();
		SBd->load(s);
		Bd = new InArray();
		Bd->load(s);
		s.loadu64(datanum);
		s.loadu32(datawidth);
		u64 len = 0;
		s.loadu64(len);
		s.loadu32array(data, len);
		return 0;
	}
private:
	u64 mask;
	u32 *data;
	u64 datanum;
	u32 datawidth;
	InArray *SBd;
	InArray *Bd;
	int SBlength;
	int Blength;
};
class BitMap
{
public:
	BitMap(unsigned long long int *bitbuff, i64 bit_len, int level, int block_size = 1024, unsigned char label = '\0', uchar **tables = NULL);
	BitMap(){};
	BitMap(uchar **tables) : Z(tables[0]), R(tables[1]) {}
	~BitMap();
	i64 Rank(i64 pos);
	i64 Rank(i64 pos, int &bit);
	void Rank(i64 pos_left, i64 pos_right, i64 &rank_left, i64 &rank_right);
	void Left(BitMap *left);
	BitMap *Left() { return left; };
	void Right(BitMap *right);
	BitMap *Right() { return right; };
	void testselect();
	unsigned char Label();
	i64 GetBitLen();
	int Getlevel();
	int Load(loadkit &s);
	int Save(savekit &S);
	u64 SizeInByte();
	u64 All01Nums() { return all01; };
	u64 PlainNums() { return plain; };
	u64 RlgNums() { return rlg01; };
	u64 PlainSizeInByte() { return plainsize; };
	u64 RlgSizeInByte() { return rlgsize; };
private:
	uchar *Z;
	uchar *R;
	BitMap(const BitMap &);
	BitMap &operator=(const BitMap &right);
	void Coding();
	int GetBit(u64 *data, i64 index);
	u16 Zeros(u16 x) { return (Z[x >> 8] == 8) ? Z[x >> 8] + Z[(uchar)x] : Z[x >> 8]; }
	u64 GetBits(u64 *buff, i64 &index, int bits);
	int GetZerosRuns(u64 *buff, i64 &index);
	int FixedDecode(u64 *buff, i64 &index);
	int GammaDecode(u64 *buff, i64 &index);
	int GetRuns(u64 *data, i64 &index, int &bit);
	void Append_g(u64 *temp, i64 &index, u32 value);
	void BitCopy(u64 *temp, i64 &index, u64 value);
	void BitCopy(u64 *temp, i64 &index, u64 *data, i64 index1, int len);
	void RL_Rank(u64 *buff, i64 &index, int bits_left, int bits_right, i64 &rank_left, i64 &rank_right, int rl_type);
	i64 RL_Rank(u64 *buff, i64 &index, int bits_num, int rl_type);
	i64 RL_Rank(u64 *buff, i64 &index, int bits_num, int rl_type, int &bit);
	i64 RL0_Rank(u64 *buff, i64 &index, int bits_num);
	int RL0_Bit(u64 *buff, i64 &index, int bits);
	i64 RL0_Rank(u64 *buff, i64 &index, int bits, int &bit);
	i64 RL1_Rank(u64 *buff, i64 &index, int bits);
	int RL1_Bit(u64 *buff, i64 &index, int bits);
	i64 RL1_Rank(u64 *buff, i64 &index, int bits, int &bit);
	void Plain_Rank(u64 *buff, i64 &index, int bits_left, int bits_right, i64 &rank_left, i64 &rank_right);
	i64 Plain_Rank(u64 *buff, i64 &index, int bits);
	int Plain_Bit(u64 *buff, i64 &index, int bits);
	i64 Plain_Rank(u64 *buff, i64 &index, int bits, int &bit);
	int level;
	unsigned char label;
	unsigned long long int *data;
	i64 bitLen;
	u64 memorysize;
	int block_size;
	int block_width;
	int super_block_size;
	BitMap *left;
	BitMap *right;
	InArray *superblock;
	InArray *block;
	InArray *coding_style;
	u64 all01 = 0;
	u64 plain = 0;
	u64 rlg01 = 0;
	u64 plainsize = 0;
	u64 rlgsize = 0;
	unsigned long long int *buff;
};
#endif
