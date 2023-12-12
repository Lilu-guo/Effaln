#include <exception>
#include "tinythread.h"
#if defined(_TTHREAD_POSIX_)
#include <unistd.h>
#include <map>
#elif defined(_TTHREAD_WIN32_)
#include <process.h>
#endif
namespace tthread
{
#if defined(_TTHREAD_WIN32_)
#define _CONDITION_EVENT_ONE 0
#define _CONDITION_EVENT_ALL 1
#endif
#if defined(_TTHREAD_WIN32_)
  condition_variable::condition_variable() : mWaitersCount(0)
  {
    mEvents[_CONDITION_EVENT_ONE] = CreateEvent(NULL, FALSE, FALSE, NULL);
    mEvents[_CONDITION_EVENT_ALL] = CreateEvent(NULL, TRUE, FALSE, NULL);
    InitializeCriticalSection(&mWaitersCountLock);
  }
#endif
#if defined(_TTHREAD_WIN32_)
  condition_variable::~condition_variable()
  {
    CloseHandle(mEvents[_CONDITION_EVENT_ONE]);
    CloseHandle(mEvents[_CONDITION_EVENT_ALL]);
    DeleteCriticalSection(&mWaitersCountLock);
  }
#endif
#if defined(_TTHREAD_WIN32_)
  void condition_variable::_wait()
  {
    int result = WaitForMultipleObjects(2, mEvents, FALSE, INFINITE);
    EnterCriticalSection(&mWaitersCountLock);
    --mWaitersCount;
    bool lastWaiter = (result == (WAIT_OBJECT_0 + _CONDITION_EVENT_ALL)) &&
                      (mWaitersCount == 0);
    LeaveCriticalSection(&mWaitersCountLock);
    if (lastWaiter)
      ResetEvent(mEvents[_CONDITION_EVENT_ALL]);
  }
#endif
#if defined(_TTHREAD_WIN32_)
  void condition_variable::notify_one()
  {
    EnterCriticalSection(&mWaitersCountLock);
    bool haveWaiters = (mWaitersCount > 0);
    LeaveCriticalSection(&mWaitersCountLock);
    if (haveWaiters)
      SetEvent(mEvents[_CONDITION_EVENT_ONE]);
  }
#endif
#if defined(_TTHREAD_WIN32_)
  void condition_variable::notify_all()
  {
    EnterCriticalSection(&mWaitersCountLock);
    bool haveWaiters = (mWaitersCount > 0);
    LeaveCriticalSection(&mWaitersCountLock);
    if (haveWaiters)
      SetEvent(mEvents[_CONDITION_EVENT_ALL]);
  }
#endif
#if defined(_TTHREAD_POSIX_)
  static thread::id _pthread_t_to_ID(const pthread_t &aHandle)
  {
    static mutex idMapLock;
    static std::map<pthread_t, unsigned long int> idMap;
    static unsigned long int idCount(1);
    lock_guard<mutex> guard(idMapLock);
    if (idMap.find(aHandle) == idMap.end())
      idMap[aHandle] = idCount++;
    return thread::id(idMap[aHandle]);
  }
#endif
  struct _thread_start_info
  {
    void (*mFunction)(void *);
    void *mArg;
    thread *mThread;
  };
#if defined(_TTHREAD_WIN32_)
  unsigned WINAPI thread::wrapper_function(void *aArg)
#elif defined(_TTHREAD_POSIX_)
  void *thread::wrapper_function(void *aArg)
#endif
  {
    _thread_start_info *ti = (_thread_start_info *)aArg;
    try
    {
      ti->mFunction(ti->mArg);
    }
    catch (...)
    {
      std::terminate();
    }
    lock_guard<mutex> guard(ti->mThread->mDataMutex);
    ti->mThread->mNotAThread = true;
    delete ti;
    return 0;
  }
  thread::thread(void (*aFunction)(void *), void *aArg)
  {
    lock_guard<mutex> guard(mDataMutex);
    _thread_start_info *ti = new _thread_start_info;
    ti->mFunction = aFunction;
    ti->mArg = aArg;
    ti->mThread = this;
    mNotAThread = false;
#if defined(_TTHREAD_WIN32_)
    mHandle = (HANDLE)_beginthreadex(0, 0, wrapper_function, (void *)ti, 0, &mWin32ThreadID);
#elif defined(_TTHREAD_POSIX_)
    if (pthread_create(&mHandle, NULL, wrapper_function, (void *)ti) != 0)
      mHandle = 0;
#endif
    if (!mHandle)
    {
      mNotAThread = true;
      delete ti;
    }
  }
  thread::~thread()
  {
    if (joinable())
      std::terminate();
  }
  void thread::join()
  {
    if (joinable())
    {
#if defined(_TTHREAD_WIN32_)
      WaitForSingleObject(mHandle, INFINITE);
      CloseHandle(mHandle);
#elif defined(_TTHREAD_POSIX_)
      pthread_join(mHandle, NULL);
#endif
    }
  }
  bool thread::joinable() const
  {
    mDataMutex.lock();
    bool result = !mNotAThread;
    mDataMutex.unlock();
    return result;
  }
  void thread::detach()
  {
    mDataMutex.lock();
    if (!mNotAThread)
    {
#if defined(_TTHREAD_WIN32_)
      CloseHandle(mHandle);
#elif defined(_TTHREAD_POSIX_)
      pthread_detach(mHandle);
#endif
      mNotAThread = true;
    }
    mDataMutex.unlock();
  }
  thread::id thread::get_id() const
  {
    if (!joinable())
      return id();
#if defined(_TTHREAD_WIN32_)
    return id((unsigned long int)mWin32ThreadID);
#elif defined(_TTHREAD_POSIX_)
    return _pthread_t_to_ID(mHandle);
#endif
  }
  unsigned thread::hardware_concurrency()
  {
#if defined(_TTHREAD_WIN32_)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return (int)si.dwNumberOfProcessors;
#elif defined(_SC_NPROCESSORS_ONLN)
    return (int)sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
    return (int)sysconf(_SC_NPROC_ONLN);
#else
    return 0;
#endif
  }
  thread::id this_thread::get_id()
  {
#if defined(_TTHREAD_WIN32_)
    return thread::id((unsigned long int)GetCurrentThreadId());
#elif defined(_TTHREAD_POSIX_)
    return _pthread_t_to_ID(pthread_self());
#endif
  }
}
