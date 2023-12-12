#ifndef STR_UTIL_H_
#define STR_UTIL_H_
#include <stdint.h>
#include <string>
static inline int
hash_string(const std::string &s)
{
	int ret = 0;
	int a = 63689;
	int b = 378551;
	for (size_t i = 0; i < s.length(); i++)
	{
		ret = (ret * a) + (int)s[i];
		if (a == 0)
		{
			a += b;
		}
		else
		{
			a *= b;
		}
		if (a == 0)
		{
			a += b;
		}
	}
	return ret;
}
static inline uint32_t hash_str(const char *str)
{
	int c;
	const uint32_t FNV_PRIME = 0x010000193;
	uint32_t hash = 0x811C9Dc5;
	while ((c = *str++) != '\0')
		hash = (hash ^ c) * FNV_PRIME;
	return hash;
}
#endif
