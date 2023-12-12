#ifndef KSW2_H_
#define KSW2_H_
#include <stdint.h>
#define KSW_NEG_INF -0x40000000
#define KSW_EZ_SCORE_ONLY 0x01
#define KSW_EZ_RIGHT 0x02 #define KSW_EZ_GENERIC_SC 0x04 #define KSW_EZ_APPROX_MAX 0x08 #define KSW_EZ_APPROX_DROP 0x10 #define KSW_EZ_EXTZ_ONLY 0x40 #define KSW_EZ_REV_CIGAR 0x80 #define KSW_EZ_SPLICE_FOR 0x100
#define KSW_EZ_SPLICE_REV 0x200
#define KSW_EZ_SPLICE_FLANK 0x400
#ifdef __cplusplus
extern "C"
{
#endif
  typedef struct
  {
    uint32_t max : 31, zdropped : 1;
    int max_q, max_t;
    int mqe, mqe_t;
    int mte, mte_q;
    int score;
    int m_cigar, n_cigar;
    int reach_end;
    uint32_t *cigar;
  } ksw_extz_t;
  void ksw_extz(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);
  void ksw_extz2_sse(void *km,
                     int qlen,
                     const uint8_t *query,
                     int tlen,
                     const uint8_t *target,
                     int8_t m,
                     const int8_t *mat,
                     int8_t q,
                     int8_t e,
                     int w,
                     int zdrop,
                     int end_bonus,
                     int flag,
                     ksw_extz_t *ez);
  void ksw_extd(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int flag, ksw_extz_t *ez);
  void ksw_extd2_sse(void *km,
                     int qlen,
                     const uint8_t *query,
                     int tlen,
                     const uint8_t *target,
                     int8_t m,
                     const int8_t *mat,
                     int8_t gapo,
                     int8_t gape,
                     int8_t gapo2,
                     int8_t gape2,
                     int w,
                     int zdrop,
                     int end_bonus,
                     int flag,
                     ksw_extz_t *ez);
  void ksw_exts2_sse(void *km,
                     int qlen,
                     const uint8_t *query,
                     int tlen,
                     const uint8_t *target,
                     int8_t m,
                     const int8_t *mat,
                     int8_t gapo,
                     int8_t gape,
                     int8_t gapo2,
                     int8_t noncan,
                     int zdrop,
                     int flag,
                     ksw_extz_t *ez);
  void ksw_extf2_sse(void *km,
                     int qlen,
                     const uint8_t *query,
                     int tlen,
                     const uint8_t *target,
                     int8_t mch,
                     int8_t mis,
                     int8_t e,
                     int w,
                     int xdrop,
                     ksw_extz_t *ez);
  int ksw_gg(void *km,
             int qlen,
             const uint8_t *query,
             int tlen,
             const uint8_t *target,
             int8_t m,
             const int8_t *mat,
             int8_t gapo,
             int8_t gape,
             int w,
             int *m_cigar_,
             int *n_cigar_,
             uint32_t **cigar_);
  int ksw_gg2(void *km,
              int qlen,
              const uint8_t *query,
              int tlen,
              const uint8_t *target,
              int8_t m,
              const int8_t *mat,
              int8_t gapo,
              int8_t gape,
              int w,
              int *m_cigar_,
              int *n_cigar_,
              uint32_t **cigar_);
  int ksw_gg2_sse(void *km,
                  int qlen,
                  const uint8_t *query,
                  int tlen,
                  const uint8_t *target,
                  int8_t m,
                  const int8_t *mat,
                  int8_t gapo,
                  int8_t gape,
                  int w,
                  int *m_cigar_,
                  int *n_cigar_,
                  uint32_t **cigar_);
  void *ksw_ll_qinit(void *km, int size, int qlen, const uint8_t *query, int m, const int8_t *mat);
  int ksw_ll_i16(void *q, int tlen, const uint8_t *target, int gapo, int gape, int *qe, int *te);
#ifdef __cplusplus
}
#endif
#ifdef HAVE_KALLOC
#include "kalloc.h"
#else
#include <stdlib.h>
#define kmalloc(km, size) malloc((size))
#define kcalloc(km, count, size) calloc((count), (size))
#define krealloc(km, ptr, size) realloc((ptr), (size))
#define kfree(km, ptr) free((ptr))
#endif
static inline uint32_t *ksw_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len)
{
  if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1] & 0xf))
  {
    if (*n_cigar == *m_cigar)
    {
      *m_cigar = *m_cigar ? (*m_cigar) << 1 : 4;
      cigar = (uint32_t *)krealloc(km, cigar, (*m_cigar) << 2);
    }
    cigar[(*n_cigar)++] = len << 4 | op;
  }
  else
    cigar[(*n_cigar) - 1] += len << 4;
  return cigar;
}
static inline void ksw_backtrack(void *km,
                                 int is_rot,
                                 int is_rev,
                                 int min_intron_len,
                                 const uint8_t *p,
                                 const int *off,
                                 const int *off_end,
                                 int n_col,
                                 int i0,
                                 int j0,
                                 int *m_cigar_,
                                 int *n_cigar_,
                                 uint32_t **cigar_)
{
  int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
  uint32_t *cigar = *cigar_, tmp;
  while (i >= 0 && j >= 0)
  {
    int force_state = -1;
    if (is_rot)
    {
      r = i + j;
      if (i < off[r])
        force_state = 2;
      if (off_end && i > off_end[r])
        force_state = 1;
      tmp = force_state < 0 ? p[(size_t)r * n_col + i - off[r]] : 0;
    }
    else
    {
      if (j < off[i])
        force_state = 2;
      if (off_end && j > off_end[i])
        force_state = 1;
      tmp = force_state < 0 ? p[(size_t)i * n_col + j - off[i]] : 0;
    }
    if (state == 0)
      state = tmp & 7;
    else if (!(tmp >> (state + 2) & 1))
      state = 0;
    if (state == 0)
      state = tmp & 7;
    if (force_state >= 0)
      state = force_state;
    if (state == 0)
      cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j;
    else if (state == 1 || (state == 3 && min_intron_len <= 0))
      cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i;
    else if (state == 3 && min_intron_len > 0)
      cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 3, 1), --i;
    else
      cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j;
  }
  if (i >= 0)
    cigar = ksw_push_cigar(km,
                           &n_cigar,
                           &m_cigar,
                           cigar,
                           min_intron_len > 0 && i >= min_intron_len ? 3 : 2,
                           i + 1);
  if (j >= 0)
    cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1);
  if (!is_rev)
    for (i = 0; i < n_cigar >> 1; ++i)
      tmp = cigar[i], cigar[i] = cigar[n_cigar - 1 - i], cigar[n_cigar - 1 - i] = tmp;
  *m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}
static inline void ksw_reset_extz(ksw_extz_t *ez)
{
  ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
  ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
  ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
}
static inline int ksw_apply_zdrop(ksw_extz_t *ez, int is_rot, int32_t H, int a, int b, int zdrop, int8_t e)
{
  int r, t;
  if (is_rot)
    r = a, t = b;
  else
    r = a + b, t = a;
  if (H > (int32_t)ez->max)
  {
    ez->max = H, ez->max_t = t, ez->max_q = r - t;
  }
  else if (t >= ez->max_t && r - t >= ez->max_q)
  {
    int tl = t - ez->max_t, ql = (r - t) - ez->max_q, l;
    l = tl > ql ? tl - ql : ql - tl;
    if (zdrop >= 0 && ez->max - H > zdrop + l * e)
    {
      ez->zdropped = 1;
      return 1;
    }
  }
  return 0;
}
#endif
