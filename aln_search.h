#ifndef ALN_SEARCH_H_
#define ALN_SEARCH_H_
#ifdef WITH_TBB
#include <tbb/tbb.h>
#include <tbb/task_group.h>
class multiseedSearchWorker
{
	int tid;
public:
	multiseedSearchWorker(const multiseedSearchWorker &W) : tid(W.tid){};
	multiseedSearchWorker(int id) : tid(id){};
	void operator()() const;
};
class multiseedSearchWorker_2p5
{
	int tid;
public:
	multiseedSearchWorker_2p5(const multiseedSearchWorker_2p5 &W) : tid(W.tid){};
	multiseedSearchWorker_2p5(int id) : tid(id){};
	void operator()() const;
};
#endif
#endif
