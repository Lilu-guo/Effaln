#include "InArray.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
using namespace std;
InArray::~InArray(void)
{
	delete[] data;
}
InArray::InArray()
{
	this->datanum = 0;
	this->datawidth = 0;
}
InArray::InArray(u64 data_num, u32 data_width)
{
	if (data_num <= 0 || data_width <= 0)
		cout << "InArray����������������" << endl;
	else
	{
		this->datanum = data_num;
		this->datawidth = data_width;
		u64 totlesize = datanum * datawidth;
		if (totlesize % 32 == 0)
			totlesize = totlesize / 32 + 1;
		else
			totlesize = (totlesize / 32) + 2;
		this->data = new u32[totlesize + 1];
		memset(data, 0, 4 * totlesize);
		mask = (u64)pow(2, datawidth) - 1;
	}
}
void InArray::SetValue(u64 index, u64 v)
{
	if (index > datanum - 1 || index < 0)
	{
		cerr << "InArray:index out of boundary:" << v << endl;
	}
	else if (v > mask)
	{
		cerr << "InArray:value:" << v << " is out of boundary" << endl;
	}
	else
	{
		u64 value = v;
		u64 anchor = (index * datawidth) >> 5;
		u64 temp1 = data[anchor];
		u32 temp2 = data[anchor + 1];
		temp1 = (temp1 << 32) + temp2;
		i32 overloop = ((anchor + 2) << 5) - (index + 1) * datawidth;
		if (overloop < 0)
		{
			value = (value >> (-overloop));
			temp1 = temp1 + value;
			data[anchor + 1] = (temp1 & (0xffffffff));
			data[anchor] = (temp1 >> 32) & (0xffffffff);
			data[anchor + 2] = (v << (32 + overloop)) & (0xffffffff);
		}
		else
		{
			value = (value << overloop);
			temp1 = temp1 + value;
			data[anchor + 1] = (temp1 & (0xffffffff));
			data[anchor] = (temp1 >> 32) & (0xffffffff);
		}
	}
}
u64 InArray::GetNum()
{
	return datanum;
}
u64 InArray::GetMemorySize()
{
	return (datanum * datawidth) / 8;
}
u32 InArray::GetDataWidth()
{
	return datawidth;
}
u64 InArray::GetValue(u64 index)
{
	if (index > datanum - 1 || index < 0)
	{
		cerr << datanum - 1 << "  " << index << " InArray:GetValue: index out of boundary" << endl;
		exit(0);
	}
	u64 anchor = (index * datawidth) >> 5;
	u64 temp1 = data[anchor];
	u32 temp2 = data[anchor + 1];
	u32 temp3 = data[anchor + 2];
	i32 overloop = ((anchor + 2) << 5) - (index + 1) * datawidth;
	temp1 = (temp1 << 32) + temp2;
	if (overloop < 0)
	{
		temp1 = (temp1 << (-overloop)) + (temp3 >> (32 + overloop));
		return temp1 & mask;
	}
	else
	{
		return (temp1 >> overloop) & mask;
	}
}
i64 InArray::write(savekit &s)
{
	s.writeu64(datanum);
	s.writeu32(datawidth);
	u64 len = (datanum * datawidth);
	if (len % 32 == 0)
		len = len / 32 + 1;
	else
		len = len / 32 + 2;
	s.writeu64(len);
	s.writeu32array(data, len);
	return 1;
}
i64 InArray::load(loadkit &s)
{
	s.loadu64(datanum);
	s.loadu32(datawidth);
	u64 len = 0;
	s.loadu64(len);
	data = new u32[len + 1];
	s.loadu32array(data, len);
	mask = (u64)pow(2, datawidth) - 1;
	return 1;
}
