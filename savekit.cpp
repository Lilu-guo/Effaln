#include "savekit.h"
void savekit::close()
{
	if (w != NULL)
		fclose(w);
	w = NULL;
}
savekit::~savekit()
{
	if (w != NULL)
		fclose(w);
}
savekit::savekit(const char *file)
{
	this->w = fopen(file, "wb");
	if (w == NULL)
	{
		cout << "Fopen error" << endl;
		exit(0);
	}
}
i32 savekit::writei64(i64 value)
{
	fwrite(&value, sizeof(i64), 1, w);
	return 0;
}
i32 savekit::writeu64(u64 value)
{
	fwrite(&value, sizeof(u64), 1, w);
	return 0;
}
i32 savekit::writei32(i32 value)
{
	fwrite(&value, sizeof(i32), 1, w);
	return 0;
}
i32 savekit::writeu32(u32 value)
{
	fwrite(&value, sizeof(u32), 1, w);
	return 0;
}
i32 savekit::writei16(i16 value)
{
	fwrite(&value, sizeof(i16), 1, w);
	return 0;
}
i32 savekit::writeu16(u16 value)
{
	fwrite(&value, sizeof(u16), 1, w);
	return 0;
}
i32 savekit::writeu8(u8 value)
{
	fwrite(&value, sizeof(u8), 1, w);
	return 0;
}
i32 savekit::writei64array(i64 *value, i32 len)
{
	fwrite(value, sizeof(i64), len, w);
	return 0;
}
i32 savekit::writeu64array(u64 *value, i32 len)
{
	fwrite(value, sizeof(u64), len, w);
	return 0;
}
i32 savekit::writei32array(i32 *value, i32 len)
{
	fwrite(value, sizeof(i32), len, w);
	return 0;
}
i32 savekit::writeu32array(u32 *value, i32 len)
{
	fwrite(value, sizeof(u32), len, w);
	return 0;
}
i32 savekit::writei16array(i16 *value, i32 len)
{
	fwrite(value, sizeof(i16), len, w);
	return 0;
}
i32 savekit::writeu16array(u16 *value, i32 len)
{
	fwrite(value, sizeof(u16), len, w);
	return 0;
}
i32 savekit::writeu8array(u8 *value, i32 len)
{
	fwrite(value, sizeof(u8), len, w);
	return 0;
}
