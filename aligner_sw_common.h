#ifndef ALIGNER_SW_COMMON_H_
#define ALIGNER_SW_COMMON_H_
#include "aligner_result.h"
struct SwResult
{
	SwResult() : alres(),
				 sws(0),
				 swcups(0),
				 swrows(0),
				 swskiprows(0),
				 swskip(0),
				 swsucc(0),
				 swfail(0),
				 swbts(0)
	{
	}
	void reset()
	{
		sws = swcups = swrows = swskiprows = swskip = swsucc =
			swfail = swbts = 0;
		alres.reset();
	}
	void reverse()
	{
		alres.reverseEdits();
	}
	bool empty() const
	{
		return alres.empty();
	}
#ifndef NDEBUG
	bool repOk() const
	{
		assert(alres.repOk());
		return true;
	}
	bool repOk(const Read &rd) const
	{
		assert(alres.repOk(rd));
		return true;
	}
#endif
	AlnRes alres;
	uint64_t sws;
	uint64_t swcups;
	uint64_t swrows;
	uint64_t swskiprows;
	uint64_t swskip;
	uint64_t swsucc;
	uint64_t swfail;
	uint64_t swbts;
};
struct SwMetrics
{
	SwMetrics() : mutex_m()
	{
		reset();
	}
	void reset()
	{
		sws = swcups = swrows = swskiprows = swskip = swsucc = swfail = swbts =
			sws10 = sws5 = sws3 =
				rshit = ungapsucc = ungapfail = ungapnodec = 0;
		exatts = exranges = exrows = exsucc = exooms = 0;
		mm1atts = mm1ranges = mm1rows = mm1succ = mm1ooms = 0;
		sdatts = sdranges = sdrows = sdsucc = sdooms = 0;
	}
	void init(
		uint64_t sws_,
		uint64_t sws10_,
		uint64_t sws5_,
		uint64_t sws3_,
		uint64_t swcups_,
		uint64_t swrows_,
		uint64_t swskiprows_,
		uint64_t swskip_,
		uint64_t swsucc_,
		uint64_t swfail_,
		uint64_t swbts_,
		uint64_t rshit_,
		uint64_t ungapsucc_,
		uint64_t ungapfail_,
		uint64_t ungapnodec_,
		uint64_t exatts_,
		uint64_t exranges_,
		uint64_t exrows_,
		uint64_t exsucc_,
		uint64_t exooms_,
		uint64_t mm1atts_,
		uint64_t mm1ranges_,
		uint64_t mm1rows_,
		uint64_t mm1succ_,
		uint64_t mm1ooms_,
		uint64_t sdatts_,
		uint64_t sdranges_,
		uint64_t sdrows_,
		uint64_t sdsucc_,
		uint64_t sdooms_)
	{
		sws = sws_;
		sws10 = sws10_;
		sws5 = sws5_;
		sws3 = sws3_;
		swcups = swcups_;
		swrows = swrows_;
		swskiprows = swskiprows_;
		swskip = swskip_;
		swsucc = swsucc_;
		swfail = swfail_;
		swbts = swbts_;
		ungapsucc = ungapsucc_;
		ungapfail = ungapfail_;
		ungapnodec = ungapnodec_;
		exatts = exatts_;
		exranges = exranges_;
		exrows = exrows_;
		exsucc = exsucc_;
		exooms = exooms_;
		mm1atts = mm1atts_;
		mm1ranges = mm1ranges_;
		mm1rows = mm1rows_;
		mm1succ = mm1succ_;
		mm1ooms = mm1ooms_;
		sdatts = sdatts_;
		sdranges = sdranges_;
		sdrows = sdrows_;
		sdsucc = sdsucc_;
		sdooms = sdooms_;
	}
	void update(const SwResult &r)
	{
		sws += r.sws;
		swcups += r.swcups;
		swrows += r.swrows;
		swskiprows += r.swskiprows;
		swskip += r.swskip;
		swsucc += r.swsucc;
		swfail += r.swfail;
		swbts += r.swbts;
	}
	void merge(const SwMetrics &r)
	{
		ThreadSafe ts(mutex_m);
		sws += r.sws;
		sws10 += r.sws10;
		sws5 += r.sws5;
		sws3 += r.sws3;
		swcups += r.swcups;
		swrows += r.swrows;
		swskiprows += r.swskiprows;
		swskip += r.swskip;
		swsucc += r.swsucc;
		swfail += r.swfail;
		swbts += r.swbts;
		rshit += r.rshit;
		ungapsucc += r.ungapsucc;
		ungapfail += r.ungapfail;
		ungapnodec += r.ungapnodec;
		exatts += r.exatts;
		exranges += r.exranges;
		exrows += r.exrows;
		exsucc += r.exsucc;
		exooms += r.exooms;
		mm1atts += r.mm1atts;
		mm1ranges += r.mm1ranges;
		mm1rows += r.mm1rows;
		mm1succ += r.mm1succ;
		mm1ooms += r.mm1ooms;
		sdatts += r.sdatts;
		sdranges += r.sdranges;
		sdrows += r.sdrows;
		sdsucc += r.sdsucc;
		sdooms += r.sdooms;
	}
	void tallyGappedDp(size_t readGaps, size_t refGaps)
	{
		size_t mx = max(readGaps, refGaps);
		if (mx < 10)
			sws10++;
		if (mx < 5)
			sws5++;
		if (mx < 3)
			sws3++;
	}
	uint64_t sws;
	uint64_t sws10;
	uint64_t sws5;
	uint64_t sws3;
	uint64_t swcups;
	uint64_t swrows;
	uint64_t swskiprows;
	uint64_t swskip;
	uint64_t swsucc;
	uint64_t swfail;
	uint64_t swbts;
	uint64_t rshit;
	uint64_t ungapsucc;
	uint64_t ungapfail;
	uint64_t ungapnodec;
	uint64_t exatts;
	uint64_t exranges;
	uint64_t exrows;
	uint64_t exsucc;
	uint64_t exooms;
	uint64_t mm1atts;
	uint64_t mm1ranges;
	uint64_t mm1rows;
	uint64_t mm1succ;
	uint64_t mm1ooms;
	uint64_t sdatts;
	uint64_t sdranges;
	uint64_t sdrows;
	uint64_t sdsucc;
	uint64_t sdooms;
	MUTEX_T mutex_m;
};
enum
{
	SW_BT_OALL_DIAG,
	SW_BT_OALL_REF_OPEN,
	SW_BT_OALL_READ_OPEN,
	SW_BT_RDGAP_EXTEND,
	SW_BT_RFGAP_EXTEND
};
#endif
