#ifndef _FAST_MUTEX_H_
#define _FAST_MUTEX_H_
#if !defined(_TTHREAD_PLATFORM_DEFINED_)
#if defined(_WIN32) || defined(__WIN32__) || defined(__WINDOWS__)
#define _TTHREAD_WIN32_
#else
#define _TTHREAD_POSIX_
#endif
#define _TTHREAD_PLATFORM_DEFINED_
#endif
#if (defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))) || \
    (defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))) ||      \
    (defined(__GNUC__) && (defined(__ppc__)))
#define _FAST_MUTEX_ASM_
#else
#define _FAST_MUTEX_SYS_
#endif
#if defined(_TTHREAD_WIN32_)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#define __UNDEF_LEAN_AND_MEAN
#endif
#include <windows.h>
#ifdef __UNDEF_LEAN_AND_MEAN
#undef WIN32_LEAN_AND_MEAN
#undef __UNDEF_LEAN_AND_MEAN
#endif
#else
#ifdef _FAST_MUTEX_ASM_
#include <sched.h>
#else
#include <pthread.h>
#endif
#endif
namespace tthread
{
  class fast_mutex
  {
  public:
#if defined(_FAST_MUTEX_ASM_)
    fast_mutex() : mLock(0)
    {
    }
#else
    fast_mutex()
    {
#if defined(_TTHREAD_WIN32_)
      InitializeCriticalSection(&mHandle);
#elif defined(_TTHREAD_POSIX_)
      pthread_mutex_init(&mHandle, NULL);
#endif
    }
#endif
#if !defined(_FAST_MUTEX_ASM_)
    ~fast_mutex()
    {
#if defined(_TTHREAD_WIN32_)
      DeleteCriticalSection(&mHandle);
#elif defined(_TTHREAD_POSIX_)
      pthread_mutex_destroy(&mHandle);
#endif
    }
#endif
    inline void lock()
    {
#if defined(_FAST_MUTEX_ASM_)
      bool gotLock;
      do
      {
        gotLock = try_lock();
        if (!gotLock)
        {
#if defined(_TTHREAD_WIN32_)
          Sleep(0);
#elif defined(_TTHREAD_POSIX_)
          sched_yield();
#endif
        }
      } while (!gotLock);
#else
#if defined(_TTHREAD_WIN32_)
      EnterCriticalSection(&mHandle);
#elif defined(_TTHREAD_POSIX_)
      pthread_mutex_lock(&mHandle);
#endif
#endif
    }
    inline bool try_lock()
    {
#if defined(_FAST_MUTEX_ASM_)
      int oldLock;
#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
      asm volatile(
          "movl $1,%%eax\n\t"
          "xchg %%eax,%0\n\t"
          "movl %%eax,%1\n\t"
          : "=m"(mLock), "=m"(oldLock)
          :
          : "%eax", "memory");
#elif defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
      int *ptrLock = &mLock;
      __asm {
        mov eax,1
        mov ecx,ptrLock
        xchg eax,[ecx]
        mov oldLock,eax
      }
#elif defined(__GNUC__) && (defined(__ppc__))
      int newLock = 1;
      asm volatile(
          "\n1:\n\t"
          "lwarx  %0,0,%1\n\t"
          "cmpwi  0,%0,0\n\t"
          "bne-   2f\n\t"
          "stwcx. %2,0,%1\n\t"
          "bne-   1b\n\t"
          "isync\n"
          "2:\n\t"
          : "=&r"(oldLock)
          : "r"(&mLock), "r"(newLock)
          : "cr0", "memory");
#endif
      return (oldLock == 0);
#else
#if defined(_TTHREAD_WIN32_)
      return TryEnterCriticalSection(&mHandle) ? true : false;
#elif defined(_TTHREAD_POSIX_)
      return (pthread_mutex_trylock(&mHandle) == 0) ? true : false;
#endif
#endif
    }
    inline void unlock()
    {
#if defined(_FAST_MUTEX_ASM_)
#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
      asm volatile(
          "movl $0,%%eax\n\t"
          "xchg %%eax,%0\n\t"
          : "=m"(mLock)
          :
          : "%eax", "memory");
#elif defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
      int *ptrLock = &mLock;
      __asm {
        mov eax,0
        mov ecx,ptrLock
        xchg eax,[ecx]
      }
#elif defined(__GNUC__) && (defined(__ppc__))
      asm volatile(
          "sync\n\t"
          :
          :
          : "memory");
      mLock = 0;
#endif
#else
#if defined(_TTHREAD_WIN32_)
      LeaveCriticalSection(&mHandle);
#elif defined(_TTHREAD_POSIX_)
      pthread_mutex_unlock(&mHandle);
#endif
#endif
    }
  private:
#if defined(_FAST_MUTEX_ASM_)
    int mLock;
#else
#if defined(_TTHREAD_WIN32_)
    CRITICAL_SECTION mHandle;
#elif defined(_TTHREAD_POSIX_)
    pthread_mutex_t mHandle;
#endif
#endif
  };
}
#endif
