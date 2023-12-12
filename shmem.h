#ifndef SHMEM_H_
#define SHMEM_H_
#ifdef _SHARED_MEM
#include <string>
#include <sys/shm.h>
#include <unistd.h>
#include <sys/shm.h>
#include <errno.h>
#include <stdint.h>
#include <stdexcept>
#include "str_util.h"
#include "btypes.h"
extern void notifySharedMem(void *mem, size_t len);
extern void waitSharedMem(void *mem, size_t len);
#define ALLOC_SHARED_U allocSharedMem<TIndexOffU>
#define ALLOC_SHARED_U8 allocSharedMem<uint8_t>
#define ALLOC_SHARED_U32 allocSharedMem<uint32_t>
#define FREE_SHARED shmdt
#define NOTIFY_SHARED notifySharedMem
#define WAIT_SHARED waitSharedMem
#define SHMEM_UNINIT 0xafba4242
#define SHMEM_INIT 0xffaa6161
template <typename T>
bool allocSharedMem(std::string fname,
					size_t len,
					T **dst,
					const char *memName,
					bool verbose)
{
	using namespace std;
	int shmid = -1;
	key_t key = (key_t)hash_string(fname);
	shmid_ds ds;
	int ret;
	size_t shmemLen = len + 4;
	if (verbose)
	{
		cerr << "Reading " << len << "+4 bytes into shared memory for " << memName << endl;
	}
	T *ptr = NULL;
	while (true)
	{
		if ((shmid = shmget(key, shmemLen, IPC_CREAT | 0666)) < 0)
		{
			if (errno == ENOMEM)
			{
				cerr << "Out of memory allocating shared area " << memName << endl;
			}
			else if (errno == EACCES)
			{
				cerr << "EACCES" << endl;
			}
			else if (errno == EEXIST)
			{
				cerr << "EEXIST" << endl;
			}
			else if (errno == EINVAL)
			{
				cerr << "Warning: shared-memory chunk's segment size doesn't match expected size (" << (shmemLen) << ")" << endl
					 << "Deleteing old shared memory block and trying again." << endl;
				shmid = shmget(key, 0, 0);
				if ((ret = shmctl(shmid, IPC_RMID, &ds)) < 0)
				{
					cerr << "shmctl returned " << ret
						 << " for IPC_RMID, errno is " << errno
						 << ", shmid is " << shmid << endl;
					throw 1;
				}
				else
				{
					cerr << "Deleted shared mem chunk with shmid " << shmid << endl;
				}
				continue;
			}
			else if (errno == ENOENT)
			{
				cerr << "ENOENT" << endl;
			}
			else if (errno == ENOSPC)
			{
				cerr << "ENOSPC" << endl;
			}
			else
			{
				cerr << "shmget returned " << shmid << " for and errno is " << errno << endl;
			}
			throw 1;
		}
		ptr = (T *)shmat(shmid, 0, 0);
		if (ptr == (void *)-1)
		{
			cerr << "Failed to attach " << memName << " to shared memory with shmat()." << endl;
			throw 1;
		}
		if (ptr == NULL)
		{
			cerr << memName << " pointer returned by shmat() was NULL." << endl;
			throw 1;
		}
		if ((ret = shmctl(shmid, IPC_STAT, &ds)) < 0)
		{
			cerr << "shmctl returned " << ret << " for IPC_STAT and errno is " << errno << endl;
			throw 1;
		}
		if (ds.shm_segsz != shmemLen)
		{
			cerr << "Warning: shared-memory chunk's segment size (" << ds.shm_segsz
				 << ") doesn't match expected size (" << shmemLen << ")" << endl
				 << "Deleteing old shared memory block and trying again." << endl;
			if ((ret = shmctl(shmid, IPC_RMID, &ds)) < 0)
			{
				cerr << "shmctl returned " << ret << " for IPC_RMID and errno is " << errno << endl;
				throw 1;
			}
		}
		else
		{
			break;
		}
	}
	*dst = ptr;
	bool initid = (((volatile uint32_t *)((char *)ptr + len))[0] == SHMEM_INIT);
	if (ds.shm_cpid == getpid() && !initid)
	{
		if (verbose)
		{
			cerr << "  I (pid = " << getpid() << ") created the "
				 << "shared memory for " << memName << endl;
		}
		((volatile uint32_t *)((char *)ptr + len))[0] = SHMEM_UNINIT;
		return true;
	}
	else
	{
		if (verbose)
		{
			cerr << "  I (pid = " << getpid()
				 << ") did not create the shared memory for "
				 << memName << ".  Pid " << ds.shm_cpid << " did." << endl;
		}
		return false;
	}
}
#else
#define ALLOC_SHARED_U(...) 0
#define ALLOC_SHARED_U8(...) 0
#define ALLOC_SHARED_U32(...) 0
#define FREE_SHARED(...)
#define NOTIFY_SHARED(...)
#define WAIT_SHARED(...)
#endif
#endif
