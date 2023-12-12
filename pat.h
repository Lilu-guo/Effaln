#ifndef PAT_H_
#define PAT_H_
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/stat.h>
#include <zlib.h>
#include <cassert>
#include <string>
#include <ctype.h>
#include <vector>
#include "alphabet.h"
#include "assert_helpers.h"
#include "random_source.h"
#include "threading.h"
#include "qual.h"
#include "search_globals.h"
#include "sstring.h"
#include "ds.h"
#include "read.h"
#include "util.h"
#ifdef USE_SRA
#include <ncbi-vdb/NGS.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>
#endif
#ifdef _WIN32
#define getc_unlocked _fgetc_nolock
#endif
struct PatternParams
{
	PatternParams() {}
	PatternParams(
		int format_,
		bool interleaved_,
		bool fileParallel_,
		uint32_t seed_,
		size_t max_buf_,
		bool solexa64_,
		bool phred64_,
		bool intQuals_,
		int trim5_,
		int trim3_,
		pair<short, size_t> trimTo_,
		int sampleLen_,
		int sampleFreq_,
		size_t skip_,
		uint64_t upto_,
		int nthreads_,
		bool fixName_,
		bool preserve_tags_,
		bool align_paired_reads_) : format(format_),
									interleaved(interleaved_),
									fileParallel(fileParallel_),
									seed(seed_),
									max_buf(max_buf_),
									solexa64(solexa64_),
									phred64(phred64_),
									intQuals(intQuals_),
									trim5(trim5_),
									trim3(trim3_),
									trimTo(trimTo_),
									sampleLen(sampleLen_),
									sampleFreq(sampleFreq_),
									skip(skip_),
									upto(upto_),
									nthreads(nthreads_),
									fixName(fixName_),
									preserve_tags(preserve_tags_),
									align_paired_reads(align_paired_reads_) {}
	int format;
	bool interleaved;
	bool fileParallel;
	uint32_t seed;
	size_t max_buf;
	bool solexa64;
	bool phred64;
	bool intQuals;
	int trim5;
	int trim3;
	pair<short, size_t> trimTo;
	int sampleLen;
	int sampleFreq;
	size_t skip;
	uint64_t upto;
	int nthreads;
	bool fixName;
	bool preserve_tags;
	bool align_paired_reads;
};
struct PerThreadReadBuf
{
	PerThreadReadBuf(size_t max_buf, int tid) : max_buf_(max_buf),
												bufa_(max_buf),
												bufb_(max_buf),
												rdid_(),
												tid_(tid)
	{
		bufa_.resize(max_buf);
		bufb_.resize(max_buf);
		reset();
	}
	Read &read_a() { return bufa_[cur_buf_]; }
	Read &read_b() { return bufb_[cur_buf_]; }
	const Read &read_a() const { return bufa_[cur_buf_]; }
	const Read &read_b() const { return bufb_[cur_buf_]; }
	TReadId rdid() const
	{
		assert_neq(rdid_, std::numeric_limits<TReadId>::max());
		return rdid_ + cur_buf_;
	}
	void reset()
	{
		cur_buf_ = bufa_.size();
		for (size_t i = 0; i < max_buf_; i++)
		{
			bufa_[i].reset();
			bufb_[i].reset();
		}
		rdid_ = std::numeric_limits<TReadId>::max();
	}
	void next()
	{
		assert_lt(cur_buf_, bufa_.size());
		cur_buf_++;
	}
	bool exhausted()
	{
		assert_leq(cur_buf_, bufa_.size());
		return cur_buf_ >= bufa_.size() - 1 || bufa_[cur_buf_ + 1].readOrigBuf.empty();
	}
	void init()
	{
		cur_buf_ = 0;
	}
	void setReadId(TReadId rdid)
	{
		rdid_ = rdid;
	}
	const size_t max_buf_;
	EList<Read> bufa_;
	EList<Read> bufb_;
	size_t cur_buf_;
	TReadId rdid_;
	int tid_;
};
extern void wrongQualityFormat(const BTString &read_name);
extern void tooFewQualities(const BTString &read_name);
extern void tooManyQualities(const BTString &read_name);
class PatternSource
{
public:
	PatternSource(const PatternParams &p) : pp_(p),
											readCnt_(0),
											mutex() {}
	virtual ~PatternSource() {}
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf &pt,
		bool batch_a,
		bool lock = true) = 0;
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const = 0;
	virtual void reset() { readCnt_ = 0; }
	static PatternSource *patsrcFromStrings(
		const PatternParams &p,
		const EList<std::string> &qs);
	TReadId readCount() const { return readCnt_; }
protected:
	const PatternParams &pp_;
	volatile TReadId readCnt_;
	MUTEX_T mutex;
};
class VectorPatternSource : public PatternSource
{
public:
	VectorPatternSource(
		const EList<std::string> &v,
		const PatternParams &p);
	virtual ~VectorPatternSource() {}
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf &pt,
		bool batch_a,
		bool lock = true);
	virtual void reset()
	{
		PatternSource::reset();
		cur_ = skip_;
		paired_ = false;
	}
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
private:
	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf &pt,
		bool batch_a);
	size_t cur_;
	size_t skip_;
	bool paired_;
	EList<string> tokbuf_;
	EList<Read::TBuf> bufs_;
	char nametmp_[20];
};
class CFilePatternSource : public PatternSource
{
public:
	CFilePatternSource(
		const EList<std::string> &infiles,
		const PatternParams &p) : PatternSource(p),
								  infiles_(infiles),
								  filecur_(0),
								  fp_(NULL),
								  zfp_(NULL),
								  is_open_(false),
								  skip_(p.skip),
								  first_(true),
								  compressed_(false)
	{
		assert_gt(infiles.size(), 0);
		errs_.resize(infiles_.size());
		errs_.fill(0, infiles_.size(), false);
		open();
		filecur_++;
	}
	virtual ~CFilePatternSource()
	{
		if (is_open_)
		{
			if (compressed_)
			{
				assert(zfp_ != NULL);
				gzclose(zfp_);
			}
			else
			{
				assert(fp_ != NULL);
				fclose(fp_);
			}
		}
	}
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf &pt,
		bool batch_a,
		bool lock = true);
	virtual void reset()
	{
		PatternSource::reset();
		filecur_ = 0,
		open();
		filecur_++;
	}
protected:
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf &pt,
		bool batch_a,
		unsigned read_idx) = 0;
	virtual void resetForNextFile() {}
	void open();
	int getc_wrapper()
	{
		int c;
		do
		{
			c = compressed_ ? gzgetc(zfp_) : getc_unlocked(fp_);
		} while (c != EOF && c != '\t' && c != '\r' && c != '\n' && !isprint(c));
		return c;
	}
	int ungetc_wrapper(int c)
	{
		return compressed_ ? gzungetc(c, zfp_) : ungetc(c, fp_);
	}
	int zread(voidp buf, unsigned len)
	{
		int r = gzread(zfp_, buf, len);
		if (r < 0)
		{
			const char *err = gzerror(zfp_, NULL);
			if (err != NULL)
			{
				std::cerr << err << std::endl;
			}
		}
		return r;
	}
	bool is_gzipped_file(int fd)
	{
		if (fd == -1)
		{
			return false;
		}
		uint8_t byte1, byte2;
		ssize_t r1 = read(fd, &byte1, sizeof(uint8_t));
		ssize_t r2 = read(fd, &byte2, sizeof(uint8_t));
		lseek(fd, 0, SEEK_SET);
		if (r1 == 0 || r2 == 0)
		{
			std::cerr << "Unable to read file magic number" << std::endl;
			return false;
		}
		if (byte1 == 0x1f && byte2 == 0x8b)
		{
			return true;
		}
		return false;
	}
	EList<std::string> infiles_;
	EList<bool> errs_;
	size_t filecur_;
	FILE *fp_;
	gzFile zfp_;
	bool is_open_;
	TReadId skip_;
	bool first_;
	char buf_[64 * 1024];
	bool compressed_;
private:
	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf &pt,
		bool batch_a);
};
class FastaPatternSource : public CFilePatternSource
{
public:
	FastaPatternSource(
		const EList<std::string> &infiles,
		const PatternParams &p, bool interleaved) : CFilePatternSource(infiles, p),
													first_(true),
													interleaved_(interleaved) {}
	virtual void reset()
	{
		first_ = true;
		CFilePatternSource::reset();
	}
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
protected:
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf &pt,
		bool batch_a,
		unsigned read_idx);
	static int skipToNextFastaRecord(FileBuf &in)
	{
		int c;
		while ((c = in.get()) != '>')
		{
			if (in.eof())
				return -1;
		}
		return c;
	}
	virtual void resetForNextFile()
	{
		first_ = true;
	}
	bool first_;
	bool interleaved_;
};
class TabbedPatternSource : public CFilePatternSource
{
public:
	TabbedPatternSource(
		const EList<std::string> &infiles,
		const PatternParams &p,
		bool secondName) : CFilePatternSource(infiles, p),
						   secondName_(secondName) {}
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
protected:
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf &pt,
		bool batch_a,
		unsigned read_idx);
	bool secondName_;
};
class QseqPatternSource : public CFilePatternSource
{
public:
	QseqPatternSource(
		const EList<std::string> &infiles,
		const PatternParams &p) : CFilePatternSource(infiles, p) {}
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
protected:
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf &pt,
		bool batch_a,
		unsigned read_idx);
	EList<std::string> qualToks_;
};
class FastaContinuousPatternSource : public CFilePatternSource
{
public:
	FastaContinuousPatternSource(
		const EList<std::string> &infiles,
		const PatternParams &p) : CFilePatternSource(infiles, p),
								  length_(p.sampleLen),
								  freq_(p.sampleFreq),
								  eat_(length_ - 1),
								  beginning_(true),
								  bufCur_(0),
								  cur_(0llu),
								  last_(0llu)
	{
		assert_gt(freq_, 0);
		resetForNextFile();
		assert_leq(length_, 1024);
	}
	virtual void reset()
	{
		CFilePatternSource::reset();
		resetForNextFile();
	}
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
protected:
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf &pt,
		bool batch_a,
		unsigned read_idx);
	virtual void resetForNextFile()
	{
		eat_ = length_ - 1;
		name_prefix_buf_.clear();
		beginning_ = true;
		bufCur_ = 0;
		last_ = cur_;
	}
private:
	const size_t length_;
	const size_t freq_;
	size_t eat_;
	bool beginning_;
	char buf_[1024];
	Read::TBuf name_prefix_buf_;
	char name_int_buf_[20];
	size_t bufCur_;
	uint64_t cur_;
	uint64_t last_;
};
class FastqPatternSource : public CFilePatternSource
{
public:
	FastqPatternSource(
		const EList<std::string> &infiles,
		const PatternParams &p, bool interleaved) : CFilePatternSource(infiles, p),
													first_(true),
													interleaved_(interleaved) {}
	virtual void reset()
	{
		first_ = true;
		CFilePatternSource::reset();
	}
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
protected:
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf &pt,
		bool batch_a,
		unsigned read_idx);
	virtual void resetForNextFile()
	{
		first_ = true;
	}
	bool first_;
	bool interleaved_;
};
class BAMPatternSource : public CFilePatternSource
{
	struct hdr_t
	{
		uint8_t id1;
		uint8_t id2;
		uint8_t cm;
		uint8_t flg;
		uint32_t mtime;
		uint8_t xfl;
		uint8_t os;
		uint16_t xlen;
	};
	struct ftr_t
	{
		uint32_t crc32;
		uint32_t isize;
	};
	struct BGZF
	{
		hdr_t hdr;
		uint8_t cdata[1 << 16];
		ftr_t ftr;
	};
	struct orphan_mate_t
	{
		orphan_mate_t() : data(NULL),
						  size(0),
						  cap(0),
						  hash(0) {}
		void reset()
		{
			size = 0;
			hash = 0;
		}
		bool empty() const
		{
			return size == 0;
		}
		uint8_t *data;
		uint16_t size;
		uint16_t cap;
		uint32_t hash;
	};
	struct BAMField
	{
		enum aln_rec_field_name
		{
			refID,
			pos,
			l_read_name,
			mapq,
			bin,
			n_cigar_op,
			flag,
			l_seq,
			next_refID,
			next_pos,
			tlen,
			read_name,
		};
	};
public:
	BAMPatternSource(
		const EList<std::string> &infiles,
		const PatternParams &p) : CFilePatternSource(infiles, p),
								  first_(true),
								  bam_batches_(p.nthreads),
								  bam_batch_indexes_(p.nthreads),
								  orphan_mates(p.nthreads * 2),
								  orphan_mates_mutex_(),
								  pp_(p)
	{
		for (size_t i = 0; i < bam_batches_.size(); ++i)
		{
			bam_batches_[i].reserve(1 << 16);
		}
	}
	virtual void reset()
	{
		first_ = true;
		CFilePatternSource::reset();
	}
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
	~BAMPatternSource()
	{
		for (size_t i = 0; i < orphan_mates.size(); i++)
		{
			if (orphan_mates[i].data != NULL)
			{
				delete[] orphan_mates[i].data;
			}
		}
	}
protected:
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf &pt, bool batch_a, bool lock = true);
	uint16_t nextBGZFBlockFromFile(BGZF &block);
	virtual void resetForNextFile()
	{
		first_ = true;
	}
	bool first_;
private:
	virtual std::pair<bool, int> nextBatchFromFile(PerThreadReadBuf &, bool, unsigned)
	{
		return make_pair(true, 0);
	}
	int decompress_bgzf_block(uint8_t *dst, size_t dst_len, uint8_t *src, size_t src_len);
	std::pair<bool, int> get_alignments(PerThreadReadBuf &pt, bool batch_a, unsigned &readi, bool lock);
	void store_orphan_mate(const uint8_t *read, size_t read_len);
	void get_orphaned_pairs(EList<Read> &buf_a, EList<Read> &buf_b, const size_t max_buf, unsigned &readi);
	void get_or_store_orhaned_mate(EList<Read> &buf_a, EList<Read> &buf_b, unsigned &readi, const uint8_t *mate, size_t mate_len);
	size_t get_matching_read(const uint8_t *rec);
	static const int offset[];
	static const uint8_t EOF_MARKER[];
	std::vector<std::vector<uint8_t>> bam_batches_;
	std::vector<size_t> bam_batch_indexes_;
	std::vector<orphan_mate_t> orphan_mates;
	MUTEX_T orphan_mates_mutex_;
	PatternParams pp_;
};
class RawPatternSource : public CFilePatternSource
{
public:
	RawPatternSource(
		const EList<std::string> &infiles,
		const PatternParams &p) : CFilePatternSource(infiles, p), first_(true) {}
	virtual void reset()
	{
		first_ = true;
		CFilePatternSource::reset();
	}
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
protected:
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf &pt,
		bool batch_a,
		unsigned read_idx);
	virtual void resetForNextFile()
	{
		first_ = true;
	}
private:
	bool first_;
};
class PatternComposer
{
public:
	PatternComposer(const PatternParams &p) : mutex_m() {}
	virtual ~PatternComposer() {}
	virtual void reset() = 0;
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf &pt) = 0;
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) = 0;
	static PatternComposer *setupPatternComposer(
		const EList<std::string> &si,
		const EList<std::string> &m1,
		const EList<std::string> &m2,
		const EList<std::string> &m12,
		const EList<std::string> &q,
		const EList<std::string> &q1,
		const EList<std::string> &q2,
#ifdef USE_SRA
		const EList<string> &sra_accs,
#endif
		PatternParams &p,
		bool verbose);
protected:
	static void free_EList_pmembers(const EList<PatternSource *> &);
	MUTEX_T mutex_m;
};
class SoloPatternComposer : public PatternComposer
{
public:
	SoloPatternComposer(
		const EList<PatternSource *> *src,
		const PatternParams &p) : PatternComposer(p),
								  cur_(0),
								  src_(src)
	{
		assert(src_ != NULL);
		lock_ = p.nthreads > 1;
		for (size_t i = 0; i < src_->size(); i++)
		{
			assert((*src_)[i] != NULL);
		}
	}
	virtual ~SoloPatternComposer()
	{
		free_EList_pmembers(*src_);
		delete src_;
	}
	virtual void reset()
	{
		for (size_t i = 0; i < src_->size(); i++)
		{
			(*src_)[i]->reset();
		}
		cur_ = 0;
	}
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf &pt);
	virtual bool parse(Read &ra, Read &rb, TReadId rdid)
	{
		return (*src_)[0]->parse(ra, rb, rdid);
	}
protected:
	volatile bool lock_;
	volatile size_t cur_;
	const EList<PatternSource *> *src_;
};
class DualPatternComposer : public PatternComposer
{
public:
	DualPatternComposer(
		const EList<PatternSource *> *srca,
		const EList<PatternSource *> *srcb,
		const PatternParams &p) : PatternComposer(p), cur_(0), srca_(srca), srcb_(srcb)
	{
		assert(srca_ != NULL);
		assert(srcb_ != NULL);
		assert_eq(srca_->size(), srcb_->size());
		lock_ = p.nthreads > 1;
		for (size_t i = 0; i < srca_->size(); i++)
		{
			assert((*srca_)[i] != NULL);
			for (size_t j = 0; j < srcb_->size(); j++)
			{
				assert_neq((*srca_)[i], (*srcb_)[j]);
			}
		}
	}
	virtual ~DualPatternComposer()
	{
		free_EList_pmembers(*srca_);
		delete srca_;
		free_EList_pmembers(*srcb_);
		delete srcb_;
	}
	virtual void reset()
	{
		for (size_t i = 0; i < srca_->size(); i++)
		{
			(*srca_)[i]->reset();
			if ((*srcb_)[i] != NULL)
			{
				(*srcb_)[i]->reset();
			}
		}
		cur_ = 0;
	}
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf &pt);
	virtual bool parse(Read &ra, Read &rb, TReadId rdid)
	{
		return (*srca_)[0]->parse(ra, rb, rdid);
	}
protected:
	volatile bool lock_;
	volatile size_t cur_;
	const EList<PatternSource *> *srca_;
	const EList<PatternSource *> *srcb_;
};
class PatternSourcePerThread
{
public:
	PatternSourcePerThread(
		PatternComposer &composer,
		const PatternParams &pp, int tid) : composer_(composer),
											buf_(pp.max_buf, tid),
											pp_(pp),
											last_batch_(false),
											last_batch_size_(0) {}
	std::pair<bool, bool> nextReadPair();
	Read &read_a() { return buf_.read_a(); }
	Read &read_b() { return buf_.read_b(); }
	const Read &read_a() const { return buf_.read_a(); }
	const Read &read_b() const { return buf_.read_b(); }
private:
	std::pair<bool, int> nextBatch()
	{
		buf_.reset();
		std::pair<bool, int> res = composer_.nextBatch(buf_);
		buf_.init();
		return res;
	}
	void finalize(Read &ra);
	void finalizePair(Read &ra, Read &rb);
	bool parse(Read &ra, Read &rb)
	{
		return composer_.parse(ra, rb, buf_.rdid());
	}
	void trim(Read &r)
	{
		if (pp_.trimTo.second > 0)
		{
			switch (pp_.trimTo.first)
			{
			case 3:
				if (r.patFw.length() > pp_.trimTo.second)
				{
					r.trimmed5 = r.patFw.length() - pp_.trimTo.second;
					r.patFw.trimEnd(r.trimmed5);
					r.qual.trimEnd(r.trimmed5);
				}
				break;
			case 5:
				if (r.patFw.length() > pp_.trimTo.second)
				{
					r.trimmed3 = r.patFw.length() - pp_.trimTo.second;
					r.patFw.trimBegin(r.trimmed3);
					r.qual.trimBegin(r.trimmed3);
				}
				break;
			}
		}
	}
	PatternComposer &composer_;
	PerThreadReadBuf buf_;
	const PatternParams &pp_;
	bool last_batch_;
	int last_batch_size_;
};
class PatternSourcePerThreadFactory
{
public:
	PatternSourcePerThreadFactory(
		PatternComposer &composer,
		const PatternParams &pp, int tid) : composer_(composer),
											pp_(pp),
											tid_(tid) {}
	virtual PatternSourcePerThread *create() const
	{
		return new PatternSourcePerThread(composer_, pp_, tid_);
	}
	virtual EList<PatternSourcePerThread *> *create(uint32_t n) const
	{
		EList<PatternSourcePerThread *> *v = new EList<PatternSourcePerThread *>;
		for (size_t i = 0; i < n; i++)
		{
			v->push_back(new PatternSourcePerThread(composer_, pp_, tid_));
			assert(v->back() != NULL);
		}
		return v;
	}
	virtual ~PatternSourcePerThreadFactory() {}
private:
	PatternComposer &composer_;
	const PatternParams &pp_;
	int tid_;
};
#ifdef USE_SRA
namespace ngs
{
	class ReadCollection;
	class ReadIterator;
}
class SRAPatternSource : public PatternSource
{
public:
	SRAPatternSource(
		const EList<string> &sra_accs,
		const PatternParams &p) : PatternSource(p),
								  sra_accs_(sra_accs),
								  sra_acc_cur_(0),
								  cur_(0),
								  first_(true),
								  sra_its_(p.nthreads),
								  mutex_m(),
								  pp_(p)
	{
		assert_gt(sra_accs_.size(), 0);
		errs_.resize(sra_accs_.size());
		errs_.fill(0, sra_accs_.size(), false);
		open();
		sra_acc_cur_++;
	}
	virtual ~SRAPatternSource()
	{
		for (size_t i = 0; i < sra_its_.size(); i++)
		{
			if (sra_its_[i] != NULL)
			{
				delete sra_its_[i];
				sra_its_[i] = NULL;
			}
		}
	}
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf &pt,
		bool batch_a,
		bool lock);
	virtual bool parse(Read &ra, Read &rb, TReadId rdid) const;
	virtual void reset()
	{
		PatternSource::reset();
		sra_acc_cur_ = 0,
		open();
		sra_acc_cur_++;
	}
protected:
	std::pair<bool, int> nextBatchImpl(
		PerThreadReadBuf &pt,
		bool batch_a);
	void open();
	EList<string> sra_accs_;
	EList<bool> errs_;
	size_t sra_acc_cur_;
	size_t cur_;
	bool first_;
	std::vector<ngs::ReadIterator *> sra_its_;
	MUTEX_T mutex_m;
	PatternParams pp_;
};
#endif
#endif
