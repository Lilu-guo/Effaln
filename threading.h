#ifndef THREADING_H_
#define THREADING_H_
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef WITH_TBB
#include <mutex>
#include <tbb/spin_mutex.h>
#include <tbb/queuing_mutex.h>
#include <atomic>
#ifdef WITH_AFFINITY
#include <sched.h>
#include <tbb/task_group.h>
#include <tbb/task_scheduler_observer.h>
#include <tbb/task_scheduler_init.h>
#endif
#else
#include "tinythread.h"
#include "fast_mutex.h"
#endif
#ifdef NO_SPINLOCK
#ifdef WITH_TBB
#ifdef WITH_QUEUELOCK
#define MUTEX_T tbb::queuing_mutex
#else
#define MUTEX_T std::mutex
#endif
#else
#define MUTEX_T tthread::mutex
#endif
#else
#ifdef WITH_TBB
#define MUTEX_T tbb::spin_mutex
#else
#define MUTEX_T tthread::fast_mutex
#endif
#endif
#ifdef WITH_TBB
struct thread_tracking_pair
{
    int tid;
    std::atomic<int> *done;
};
#endif
class ThreadSafe
{
public:
    ThreadSafe(MUTEX_T &ptr_mutex) : mutex_(ptr_mutex)
    {
#if WITH_TBB && NO_SPINLOCK && WITH_QUEUELOCK
#else
        mutex_.lock();
#endif
    }
    ~ThreadSafe()
    {
#if WITH_TBB && NO_SPINLOCK && WITH_QUEUELOCK
#else
        mutex_.unlock();
#endif
    }
private:
#if WITH_TBB && NO_SPINLOCK && WITH_QUEUELOCK
    MUTEX_T::scoped_lock mutex_;
#else
    MUTEX_T &mutex_;
#endif
};
#if defined(_TTHREAD_WIN32_)
#define SLEEP(x) Sleep(x)
#else
#define SLEEP(x)                                          \
    do                                                    \
    {                                                     \
        const static timespec ts_tmp_ = {0, 1000000 * x}; \
        nanosleep(&ts_tmp_, NULL);                        \
    } while (false)
#endif
#ifdef WITH_TBB
#ifdef WITH_AFFINITY
class concurrency_tracker : public tbb::task_scheduler_observer
{
    std::atomic<int> num_threads;
public:
    concurrency_tracker() : num_threads() { observe(true); }
    void on_scheduler_entry(bool) { ++num_threads; }
    void on_scheduler_exit(bool) { --num_threads; }
    int get_concurrency() { return num_threads; }
};
class pinning_observer : public tbb::task_scheduler_observer
{
    cpu_set_t *mask;
    int ncpus;
    const int pinning_step;
    std::atomic<int> thread_index;
public:
    pinning_observer(int pinning_step = 1) : pinning_step(pinning_step), thread_index()
    {
        for (ncpus = sizeof(cpu_set_t) / CHAR_BIT; ncpus < 16 * 1024; ncpus <<= 1)
        {
            mask = CPU_ALLOC(ncpus);
            if (!mask)
                break;
            const size_t size = CPU_ALLOC_SIZE(ncpus);
            CPU_ZERO_S(size, mask);
            const int err = sched_getaffinity(0, size, mask);
            if (!err)
                break;
            CPU_FREE(mask);
            mask = NULL;
            if (errno != EINVAL)
                break;
        }
        if (!mask)
            std::cout << "Warning: Failed to obtain process affinity mask. Thread affinitization is disabled." << std::endl;
    }
    void on_scheduler_entry(bool)
    {
        if (!mask)
            return;
        const size_t size = CPU_ALLOC_SIZE(ncpus);
        const int num_cpus = CPU_COUNT_S(size, mask);
        int thr_idx =
#if USE_TASK_ARENA_CURRENT_SLOT
            tbb::task_arena::current_slot();
#else
            thread_index++;
#endif
#if __MIC__
        thr_idx += 1;
#endif
        thr_idx %= num_cpus;
        int cpu_idx = 0;
        for (int i = 0, offset = 0; i < thr_idx; ++i)
        {
            cpu_idx += pinning_step;
            if (cpu_idx >= num_cpus)
                cpu_idx = ++offset;
        }
        int mapped_idx = -1;
        while (cpu_idx >= 0)
        {
            if (CPU_ISSET_S(++mapped_idx, size, mask))
                --cpu_idx;
        }
        cpu_set_t *target_mask = CPU_ALLOC(ncpus);
        CPU_ZERO_S(size, target_mask);
        CPU_SET_S(mapped_idx, size, target_mask);
        const int err = sched_setaffinity(0, size, target_mask);
        if (err)
        {
            std::cout << "Failed to set thread affinity!n";
            exit(EXIT_FAILURE);
        }
#if LOG_PINNING
        else
        {
            std::stringstream ss;
            ss << "Set thread affinity: Thread " << thr_idx << ": CPU " << mapped_idx << std::endl;
            std::cerr << ss.str();
        }
#endif
        CPU_FREE(target_mask);
    }
    ~pinning_observer()
    {
        if (mask)
            CPU_FREE(mask);
    }
};
#endif
#endif
#endif
