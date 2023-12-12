#define bit_SSE3	(1 << 0)
#define bit_PCLMUL	(1 << 1)
#define bit_SSSE3	(1 << 9)
#define bit_FMA		(1 << 12)
#define bit_CMPXCHG16B	(1 << 13)
#define bit_SSE4_1	(1 << 19)
#define bit_SSE4_2	(1 << 20)
#define bit_MOVBE	(1 << 22)
#define bit_POPCNT	(1 << 23)
#define bit_AES		(1 << 25)
#define bit_XSAVE	(1 << 26)
#define bit_OSXSAVE	(1 << 27)
#define bit_AVX		(1 << 28)
#define bit_F16C	(1 << 29)
#define bit_RDRND	(1 << 30)
#define bit_CMPXCHG8B	(1 << 8)
#define bit_CMOV	(1 << 15)
#define bit_MMX		(1 << 23)
#define bit_FXSAVE	(1 << 24)
#define bit_SSE		(1 << 25)
#define bit_SSE2	(1 << 26)
#define bit_LAHF_LM	(1 << 0)
#define bit_ABM		(1 << 5)
#define bit_SSE4a	(1 << 6)
#define bit_XOP         (1 << 11)
#define bit_LWP 	(1 << 15)
#define bit_FMA4        (1 << 16)
#define bit_TBM         (1 << 21)
#define bit_LM		(1 << 29)
#define bit_3DNOWP	(1 << 30)
#define bit_3DNOW	(1 << 31)
#define bit_FSGSBASE	(1 << 0)
#define bit_BMI		(1 << 3)
#if defined(__i386__) && defined(__PIC__)
#if __GNUC__ >= 3
#define __cpuid(level, a, b, c, d)			\
  __asm__ ("xchg{l}\t{%%}ebx, %1\n\t"			\
	   "cpuid\n\t"					\
	   "xchg{l}\t{%%}ebx, %1\n\t"			\
	   : "=a" (a), "=r" (b), "=c" (c), "=d" (d)	\
	   : "0" (level))
#define __cpuid_count(level, count, a, b, c, d)		\
  __asm__ ("xchg{l}\t{%%}ebx, %1\n\t"			\
	   "cpuid\n\t"					\
	   "xchg{l}\t{%%}ebx, %1\n\t"			\
	   : "=a" (a), "=r" (b), "=c" (c), "=d" (d)	\
	   : "0" (level), "2" (count))
#else
#define __cpuid(level, a, b, c, d)			\
  __asm__ ("xchgl\t%%ebx, %1\n\t"			\
	   "cpuid\n\t"					\
	   "xchgl\t%%ebx, %1\n\t"			\
	   : "=a" (a), "=r" (b), "=c" (c), "=d" (d)	\
	   : "0" (level))
#define __cpuid_count(level, count, a, b, c, d)		\
  __asm__ ("xchgl\t%%ebx, %1\n\t"			\
	   "cpuid\n\t"					\
	   "xchgl\t%%ebx, %1\n\t"			\
	   : "=a" (a), "=r" (b), "=c" (c), "=d" (d)	\
	   : "0" (level), "2" (count))
#endif
#else
#define __cpuid(level, a, b, c, d)			\
  __asm__ ("cpuid\n\t"					\
	   : "=a" (a), "=b" (b), "=c" (c), "=d" (d)	\
	   : "0" (level))
#define __cpuid_count(level, count, a, b, c, d)		\
  __asm__ ("cpuid\n\t"					\
	   : "=a" (a), "=b" (b), "=c" (c), "=d" (d)	\
	   : "0" (level), "2" (count))
#endif
static __inline unsigned int
__get_cpuid_max (unsigned int __ext, unsigned int *__sig)
{
  unsigned int __eax, __ebx, __ecx, __edx;
#ifndef __x86_64__
#if __GNUC__ >= 3
  __asm__ ("pushf{l|d}\n\t"
	   "pushf{l|d}\n\t"
	   "pop{l}\t%0\n\t"
	   "mov{l}\t{%0, %1|%1, %0}\n\t"
	   "xor{l}\t{%2, %0|%0, %2}\n\t"
	   "push{l}\t%0\n\t"
	   "popf{l|d}\n\t"
	   "pushf{l|d}\n\t"
	   "pop{l}\t%0\n\t"
	   "popf{l|d}\n\t"
	   : "=&r" (__eax), "=&r" (__ebx)
	   : "i" (0x00200000));
#else
  __asm__ ("pushfl\n\t"
	   "pushfl\n\t"
	   "popl\t%0\n\t"
	   "movl\t%0, %1\n\t"
	   "xorl\t%2, %0\n\t"
	   "pushl\t%0\n\t"
	   "popfl\n\t"
	   "pushfl\n\t"
	   "popl\t%0\n\t"
	   "popfl\n\t"
	   : "=&r" (__eax), "=&r" (__ebx)
	   : "i" (0x00200000));
#endif
  if (!((__eax ^ __ebx) & 0x00200000))
    return 0;
#endif
  __cpuid (__ext, __eax, __ebx, __ecx, __edx);
  if (__sig)
    *__sig = __ebx;
  return __eax;
}
static __inline int
__get_cpuid (unsigned int __level,
	     unsigned int *__eax, unsigned int *__ebx,
	     unsigned int *__ecx, unsigned int *__edx)
{
  unsigned int __ext = __level & 0x80000000;
  if (__get_cpuid_max (__ext, 0) < __level)
    return 0;
  __cpuid (__level, *__eax, *__ebx, *__ecx, *__edx);
  return 1;
}
