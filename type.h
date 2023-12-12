#pragma once
#include "const.h"
struct Alignment
{
  std::string cigar_string;
  int ref_begin;
  int mismatches;
};
struct Region
{
  uint32_t rs, re;
  uint32_t qs, qe;
  uint16_t cov;
  uint16_t embed_dist;
  int score;
  bool operator()(Region &X, Region &Y)
  {
    if (X.embed_dist == Y.embed_dist)
      return X.cov > Y.cov;
    else
      return X.embed_dist < Y.embed_dist;
  }
};
struct eRead
{
  char name[MAX_LEN], qua[MAX_LEN], seq[MAX_LEN], fwd[MAX_LEN], rev[MAX_LEN], rev_str[MAX_LEN], cigar[MAX_LEN];
  int tid, as, nm;
  uint32_t pos;
  char strand;
  short mapq;
  Region best_region;
  friend gzFile &operator>>(gzFile &in, eRead &r);
};
class Reference
{
public:
  void load_index(const char *F);
  Reference(const char *F);
  std::string ref;
  std::vector<std::string> name;
  std::vector<uint32_t> offset;
  uint32_t *keyv, *posv;
  uint32_t nposv, nkeyv;
  ~Reference();
};
typedef std::tuple<eRead *, eRead *, int> ReadCnt;
typedef std::tuple<eRead *, eRead *> ReadPair;
typedef struct
{
  uint32_t capacity;
  uint32_t n_cigar;
  int32_t dp_score;
  uint32_t cigar[];
} mm_extra_t;