#ifndef UTIL_H_
#define UTIL_H_
#include <stdlib.h>
#include <limits>
template <typename T>
char *itoa10(const T &value, char *result)
{
	char *out = result;
	T quotient = value;
	if (std::numeric_limits<T>::is_signed)
	{
		if (quotient <= 0)
			quotient = -quotient;
	}
	do
	{
		*out = "0123456789"[quotient % 10];
		++out;
		quotient /= 10;
	} while (quotient > 0);
	if (std::numeric_limits<T>::is_signed)
	{
		if (value <= 0 && value != 0)
			*out++ = '-';
	}
	reverse(result, out);
	*out = 0;
	return out;
}
#endif
