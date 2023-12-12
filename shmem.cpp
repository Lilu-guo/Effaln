#ifdef _SHARED_MEM
#include <iostream>
#include <string>
#include <unistd.h>
#include <sys/shm.h>
#include <errno.h>
#include "shmem.h"
using namespace std;
void notifySharedMem(void *mem, size_t len)
{
	((volatile uint32_t *)((char *)mem + len))[0] = SHMEM_INIT;
}
void waitSharedMem(void *mem, size_t len)
{
	while (((volatile uint32_t *)((char *)mem + len))[0] != SHMEM_INIT)
	{
		sleep(1);
	}
}
#endif
