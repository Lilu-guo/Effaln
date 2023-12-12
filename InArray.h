#ifndef _Inarray
#define _Inarray
#include "BaseClass.h"
#include "savekit.h"
#include "loadkit.h"
class InArray
{
public:
	InArray();
	InArray(u64 data_num, u32 data_width);
	~InArray(void);
	u64 GetValue(u64 index);
	void SetValue(u64 index, u64 value);
	u64 GetNum();
	u32 GetDataWidth();
	u64 GetMemorySize();
	u64 GetValue2(u64 index);
	i64 write(savekit &s);
	i64 load(loadkit &s);
private:
	u32 *data;
	u64 datanum;
	u32 datawidth;
	u64 mask;
};
#endif