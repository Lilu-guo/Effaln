#ifndef _TINYTHREAD_H_
#define _TINYTHREAD_H_
#if !defined(_TTHREAD_PLATFORM_DEFINED_)
#if defined(_WIN32) || defined(__WIN32__) || defined(__WINDOWS__)
#define _TTHREAD_WIN32_
#else
#define _TTHREAD_POSIX_
#endif
#define _TTHREAD_PLATFORM_DEFINED_
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
#include <pthread.h>
#include <signal.h>
#include <sched.h>
#include <unistd.h>
#endif
#include <ostream>
#define TINYTHREAD_VERSION_MAJOR 1
#define TINYTHREAD_VERSION_MINOR 1
#define TINYTHREAD_VERSION (TINYTHREAD_VERSION_MAJOR * 100 + TINYTHREAD_VERSION_MINOR)
#if (__cplusplus > 199711L) || (defined(__STDCXX_VERSION__) && (__STDCXX_VERSION__ >= 201001L))
#define _TTHREAD_CPP11_
#endif
#if defined(_TTHREAD_CPP11_) || defined(__GXX_EXPERIMENTAL_CXX0X__) || defined(__GXX_EXPERIMENTAL_CPP0X__)
#define _TTHREAD_CPP11_PARTIAL_
#endif
#ifdef _TTHREAD_CPP11_PARTIAL_
#define _TTHREAD_DISABLE_ASSIGNMENT(name) \
  name(const name &) = delete;            \
  name &operator=(const name &) = delete;
#else
#define _TTHREAD_DISABLE_ASSIGNMENT(name) \
  name(const name &);                     \
  name &operator=(const name &);
#endif
#if !defined(_TTHREAD_CPP11_) && !defined(thread_local)
#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__SUNPRO_CC) || defined(__IBMCPP__)
#define thread_local __thread
#else
#define thread_local __declspec(thread)
#endif
#endif
namespace tthread
{
  class mutex
  {
  public:
    mutex()
#if defined(_TTHREAD_WIN32_)
        : mAlreadyLocked(false)
#endif
    {
#if defined(_TTHREAD_WIN32_)
      InitializeCriticalSection(&mHandle);
#else
      pthread_mutex_init(&mHandle, NULL);
#endif
    }
    ~mutex()
    {
#if defined(_TTHREAD_WIN32_)
      DeleteCriticalSection(&mHandle);
#else
      pthread_mutex_destroy(&mHandle);
#endif
    }
    inline void lock()
    {
#if defined(_TTHREAD_WIN32_)
      EnterCriticalSection(&mHandle);
      while (mAlreadyLocked)
        Sleep(1000);
      mAlreadyLocked = true;
#else
      pthread_mutex_lock(&mHandle);
#endif
    }
    inline bool try_lock()
    {
#if defined(_TTHREAD_WIN32_)
      bool ret = (TryEnterCriticalSection(&mHandle) ? true : false);
      if (ret && mAlreadyLocked)
      {
        LeaveCriticalSection(&mHandle);
        ret = false;
      }
      return ret;
#else
      return (pthread_mutex_trylock(&mHandle) == 0) ? true : false;
#endif
    }
    inline void unlock()
    {
#if defined(_TTHREAD_WIN32_)
      mAlreadyLocked = false;
      LeaveCriticalSection(&mHandle);
#else
      pthread_mutex_unlock(&mHandle);
#endif
    }
    _TTHREAD_DISABLE_ASSIGNMENT(mutex)
  private:
#if defined(_TTHREAD_WIN32_)
    CRITICAL_SECTION mHandle;
    bool mAlreadyLocked;
#else
    pthread_mutex_t mHandle;
#endif
    friend class condition_variable;
  };
  class recursive_mutex
  {
  public:
    recursive_mutex()
    {
#if defined(_TTHREAD_WIN32_)
      InitializeCriticalSection(&mHandle);
#else
      pthread_mutexattr_t attr;
      pthread_mutexattr_init(&attr);
      pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
      pthread_mutex_init(&mHandle, &attr);
#endif
    }
    ~recursive_mutex()
    {
#if defined(_TTHREAD_WIN32_)
      DeleteCriticalSection(&mHandle);
#else
      pthread_mutex_destroy(&mHandle);
#endif
    }
    inline void lock()
    {
#if defined(_TTHREAD_WIN32_)
      EnterCriticalSection(&mHandle);
#else
      pthread_mutex_lock(&mHandle);
#endif
    }
    inline bool try_lock()
    {
#if defined(_TTHREAD_WIN32_)
      return TryEnterCriticalSection(&mHandle) ? true : false;
#else
      return (pthread_mutex_trylock(&mHandle) == 0) ? true : false;
#endif
    }
    inline void unlock()
    {
#if defined(_TTHREAD_WIN32_)
      LeaveCriticalSection(&mHandle);
#else
      pthread_mutex_unlock(&mHandle);
#endif
    }
    _TTHREAD_DISABLE_ASSIGNMENT(recursive_mutex)
  private:
#if defined(_TTHREAD_WIN32_)
    CRITICAL_SECTION mHandle;
#else
    pthread_mutex_t mHandle;
#endif
    friend class condition_variable;
  };
  template <class T>
  class lock_guard
  {
  public:
    typedef T mutex_type;
    lock_guard() : mMutex(0) {}
    explicit lock_guard(mutex_type &aMutex)
    {
      mMutex = &aMutex;
      mMutex->lock();
    }
    ~lock_guard()
    {
      if (mMutex)
        mMutex->unlock();
    }
  private:
    mutex_type *mMutex;
  };
  class condition_variable
  {
  public:
#if defined(_TTHREAD_WIN32_)
    condition_variable();
#else
    condition_variable()
    {
      pthread_cond_init(&mHandle, NULL);
    }
#endif
#if defined(_TTHREAD_WIN32_)
    ~condition_variable();
#else
    ~condition_variable()
    {
      pthread_cond_destroy(&mHandle);
    }
#endif
    template <class _mutexT>
    inline void wait(_mutexT &aMutex)
    {
#if defined(_TTHREAD_WIN32_)
      EnterCriticalSection(&mWaitersCountLock);
      ++mWaitersCount;
      LeaveCriticalSection(&mWaitersCountLock);
      aMutex.unlock();
      _wait();
      aMutex.lock();
#else
      pthread_cond_wait(&mHandle, &aMutex.mHandle);
#endif
    }
#if defined(_TTHREAD_WIN32_)
    void notify_one();
#else
    inline void notify_one()
    {
      pthread_cond_signal(&mHandle);
    }
#endif
#if defined(_TTHREAD_WIN32_)
    void notify_all();
#else
    inline void notify_all()
    {
      pthread_cond_broadcast(&mHandle);
    }
#endif
    _TTHREAD_DISABLE_ASSIGNMENT(condition_variable)
  private:
#if defined(_TTHREAD_WIN32_)
    void _wait();
    HANDLE mEvents[2];
    unsigned int mWaitersCount;
    CRITICAL_SECTION mWaitersCountLock;
#else
    pthread_cond_t mHandle;
#endif
  };
  class thread
  {
  public:
#if defined(_TTHREAD_WIN32_)
    typedef HANDLE native_handle_type;
#else
    typedef pthread_t native_handle_type;
#endif
    class id;
    thread() : mHandle(0), mNotAThread(true)
#if defined(_TTHREAD_WIN32_)
               ,
               mWin32ThreadID(0)
#endif
    {
    }
    thread(void (*aFunction)(void *), void *aArg);
    ~thread();
    void join();
    bool joinable() const;
    void detach();
    id get_id() const;
    inline native_handle_type native_handle()
    {
      return mHandle;
    }
    static unsigned hardware_concurrency();
    _TTHREAD_DISABLE_ASSIGNMENT(thread)
  private:
    native_handle_type mHandle;
    mutable mutex mDataMutex;
    bool mNotAThread;
#if defined(_TTHREAD_WIN32_)
    unsigned int mWin32ThreadID;
#endif
#if defined(_TTHREAD_WIN32_)
    static unsigned WINAPI wrapper_function(void *aArg);
#else
    static void *wrapper_function(void *aArg);
#endif
  };
  class thread::id
  {
  public:
    id() : mId(0){};
    id(unsigned long int aId) : mId(aId){};
    id(const id &aId) : mId(aId.mId){};
    inline id &operator=(const id &aId)
    {
      mId = aId.mId;
      return *this;
    }
    inline friend bool operator==(const id &aId1, const id &aId2)
    {
      return (aId1.mId == aId2.mId);
    }
    inline friend bool operator!=(const id &aId1, const id &aId2)
    {
      return (aId1.mId != aId2.mId);
    }
    inline friend bool operator<=(const id &aId1, const id &aId2)
    {
      return (aId1.mId <= aId2.mId);
    }
    inline friend bool operator<(const id &aId1, const id &aId2)
    {
      return (aId1.mId < aId2.mId);
    }
    inline friend bool operator>=(const id &aId1, const id &aId2)
    {
      return (aId1.mId >= aId2.mId);
    }
    inline friend bool operator>(const id &aId1, const id &aId2)
    {
      return (aId1.mId > aId2.mId);
    }
    inline friend std::ostream &operator<<(std::ostream &os, const id &obj)
    {
      os << obj.mId;
      return os;
    }
  private:
    unsigned long int mId;
  };
  typedef long long __intmax_t;
  template <__intmax_t N, __intmax_t D = 1>
  class ratio
  {
  public:
    static double _as_double() { return double(N) / double(D); }
  };
  namespace chrono
  {
    template <class _Rep, class _Period = ratio<1>>
    class duration
    {
    private:
      _Rep rep_;
    public:
      typedef _Rep rep;
      typedef _Period period;
      template <class _Rep2>
      explicit duration(const _Rep2 &r) : rep_(r){};
      rep count() const
      {
        return rep_;
      }
    };
    typedef duration<__intmax_t, ratio<1, 1000000000>> nanoseconds;
    typedef duration<__intmax_t, ratio<1, 1000000>> microseconds;
    typedef duration<__intmax_t, ratio<1, 1000>> milliseconds;
    typedef duration<__intmax_t> seconds;
    typedef duration<__intmax_t, ratio<60>> minutes;
    typedef duration<__intmax_t, ratio<3600>> hours;
  }
  namespace this_thread
  {
    thread::id get_id();
    inline void yield()
    {
#if defined(_TTHREAD_WIN32_)
      Sleep(0);
#else
      sched_yield();
#endif
    }
    template <class _Rep, class _Period>
    void sleep_for(const chrono::duration<_Rep, _Period> &aTime)
    {
#if defined(_TTHREAD_WIN32_)
      Sleep(int(double(aTime.count()) * (1000.0 * _Period::_as_double()) + 0.5));
#else
      usleep(int(double(aTime.count()) * (1000000.0 * _Period::_as_double()) + 0.5));
#endif
    }
  }
}
#undef _TTHREAD_DISABLE_ASSIGNMENT
#endif
