#ifndef OUTQ_H_
#define OUTQ_H_
#include "assert_helpers.h"
#include "ds.h"
#include "sstring.h"
#include "read.h"
#include "threading.h"
#include "mem_ids.h"
#include <vector>
class OutputQueue
{
	static const size_t NFLUSH_THRESH = 8;
public:
	OutputQueue(
		OutFileBuf &obuf,
		bool reorder,
		size_t nthreads,
		bool threadSafe,
		int perThreadBufSize,
		TReadId rdid = 0) : obuf_(obuf),
							cur_(rdid),
							nfinished_(0),
							nflushed_(0),
							lines_(RES_CAT),
							started_(RES_CAT),
							finished_(RES_CAT),
							reorder_(reorder),
							threadSafe_(threadSafe),
							mutex_m(),
							nthreads_(nthreads),
							perThreadBuf(NULL),
							perThreadCounter(NULL),
							perThreadBufSize_(perThreadBufSize)
	{
		nstarted_ = 0;
		assert(nthreads_ <= 2 || threadSafe);
		if (!reorder)
		{
			perThreadBuf = new BTString *[nthreads_];
			perThreadCounter = new int[nthreads_];
			size_t i = 0;
			for (i = 0; i < nthreads_; i++)
			{
				perThreadBuf[i] = new BTString[perThreadBufSize_];
				perThreadCounter[i] = 0;
			}
		}
	}
	~OutputQueue()
	{
		if (perThreadBuf != NULL)
		{
			for (size_t i = 0; i < nthreads_; i++)
			{
				delete[] perThreadBuf[i];
			}
			delete[] perThreadBuf;
			delete[] perThreadCounter;
		}
	}
	void beginRead(TReadId rdid, size_t threadId);
	void finishRead(const BTString &rec, TReadId rdid, size_t threadId);
	size_t size() const
	{
		return lines_.size();
	}
	TReadId numFlushed() const
	{
		return nflushed_;
	}
	TReadId numStarted() const
	{
		return nstarted_;
	}
	TReadId numFinished() const
	{
		return nfinished_;
	}
	void flush(bool force = false, bool getLock = true);
protected:
	OutFileBuf &obuf_;
	TReadId cur_;
#ifdef WITH_TBB
	std::atomic<TReadId> nstarted_;
#else
	TReadId nstarted_;
#endif
	TReadId nfinished_;
	TReadId nflushed_;
	EList<BTString> lines_;
	EList<bool> started_;
	EList<bool> finished_;
	bool reorder_;
	bool threadSafe_;
	MUTEX_T mutex_m;
	size_t nthreads_;
	BTString **perThreadBuf;
	int *perThreadCounter;
	int perThreadBufSize_;
private:
	void flushImpl(bool force);
	void beginReadImpl(TReadId rdid, size_t threadId);
	void finishReadImpl(const BTString &rec, TReadId rdid, size_t threadId);
};
class OutputQueueMark
{
public:
	OutputQueueMark(
		OutputQueue &q,
		const BTString &rec,
		TReadId rdid,
		size_t threadId) : q_(q),
						   rec_(rec),
						   rdid_(rdid),
						   threadId_(threadId)
	{
		q_.beginRead(rdid, threadId);
	}
	~OutputQueueMark()
	{
		q_.finishRead(rec_, rdid_, threadId_);
	}
protected:
	OutputQueue &q_;
	const BTString &rec_;
	TReadId rdid_;
	size_t threadId_;
};
#endif
