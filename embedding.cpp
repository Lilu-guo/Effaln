#include "header.h"
#include <unistd.h>
#ifndef myNDEBUG
#define jj printf
#else
#define jj
#endif
#define EMBED_PAD 4
void Embedding::cgk2_embedQ(const char *oridata, unsigned rlen, int strid, char *embeddedQ)
{
  int j = 0;
  int elen = efactor * rlen;
  for (unsigned i = 0; i < rlen; i++)
  {
    uint8_t s = oridata[i];
    char bit = hash_eb[0][BITPOS(strid, j, s)];
    if (!bit)
    {
      embeddedQ[j] = s;
      j++;
    }
    else
    {
      embeddedQ[j + 1] = embeddedQ[j] = s;
      j += 2;
    }
  }
  for (; j < elen; j++)
  {
    embeddedQ[j] = EMBED_PAD;
  }
  assert(j <= elen);
}
int Embedding::cgk2_embed_nmismatch(const char *oridata, unsigned rlen, int threshold, int strid, char *embeddedQ)
{
  int nmismatch = 0;
  int j = 0;
  int elen = efactor * rlen;
  for (unsigned i = 0; i < rlen; i++)
  {
    uint8_t s = oridata[i];
    char bit = hash_eb[0][BITPOS(strid, j, s)];
    if (!bit)
    {
      nmismatch += (embeddedQ[j] == s ? 0 : 1);
      if (nmismatch > threshold)
      {
        nmismatch = elen;
        goto end;
      }
      j++;
    }
    else
    {
      nmismatch += (embeddedQ[j] == s ? 0 : 1);
      nmismatch += (embeddedQ[j + 1] == s ? 0 : 1);
      if (nmismatch > threshold)
      {
        goto end;
      }
      j += 2;
    }
  }
  for (; j < elen; j++)
    nmismatch += (embeddedQ[j] == EMBED_PAD ? 0 : 1);
  assert(j == elen);
end:
  return nmismatch;
}
int Embedding::cgk2_embed(const char **oridata, unsigned rlen, int threshold, int id,
                          int strid, char *embeddedQ)
{
  int nmismatch = 0;
  if (id == 0)
  {
    cgk2_embedQ(oridata[id], rlen, strid, embeddedQ);
  }
  else
  {
    nmismatch = cgk2_embed_nmismatch(oridata[id], rlen, threshold, strid, embeddedQ);
  }
  return nmismatch;
}
int Embedding::cgk_embed(const char **oridata, unsigned rlen, int threshold, int id, int strid, char *embeddedQ, int n)
{
  unsigned j = 0;
  int nmismatch = 0;
  int elen = 2 * rlen;
  assert(elen < MAX_ELEN);
  if (id != 2)
  {
    for (int i = 0; i < rlen; i++)
    {
      uint8_t s = oridata[id][i];
      char bit = hash_eb[n][BITPOS(strid, j, s)];
      if (!bit)
      {
        embeddedQ[j] = s;
        j = j + 1;
      }
      else
      {
        embeddedQ[j + 1] = embeddedQ[j] = s;
        j = j + 2;
      }
    }
    for (; j < elen; j++)
    {
      embeddedQ[j] = EMBED_PAD;
    }
  }
  else
  {
    for (int i = 0; i < rlen; i++)
    {
      uint8_t s = oridata[id][i];
      char bit = hash_eb[n][BITPOS(strid, j, s)];
      if (!bit)
      {
        nmismatch += (embeddedQ[j] == s ? 0 : 1);
        j = j + 1;
      }
      else
      {
        nmismatch += (embeddedQ[j] == s ? 0 : 1);
        nmismatch += (embeddedQ[j + 1] == s ? 0 : 1);
        j = j + 2;
      }
      if (nmismatch > threshold)
      {
        nmismatch = elen;
        goto end;
      }
    }
    for (; j < elen; j++)
    {
      nmismatch += (embeddedQ[j] == EMBED_PAD ? 0 : 1);
      if (nmismatch > threshold)
      {
        nmismatch = elen;
        goto end;
      }
    }
  }
end:
  return nmismatch;
}
int Embedding::embedstr(const char **oridata, unsigned rlen, int threshold, int id,
                        int strid, char *embeddedQ)
{
#ifdef CGK2_EMBED
  return cgk2_embed(oridata, rlen, threshold, id, strid, embeddedQ);
#else
  return cgk_embed(oridata, rlen, threshold, id, strid, embeddedQ, 0);
#endif
}
void Embedding::embed_two_pairs(vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                                unsigned rlen, int *thresholds, char embeddedQ_f1[], char embeddedQ_r2[],
                                const char **candidate_refs_f1, const char **candidate_refs_r2,
                                unsigned nregions_f1, unsigned nregions_r2, bool flag_f1[], bool flag_r2[],
                                unsigned best_f1, unsigned next_f1, unsigned best_r2, unsigned next_r2)
{
  if (best_f1 == nregions_f1 || best_r2 == nregions_r2)
    return;
  int elen = rlen * efactor;
  int nmismatch_f1 = embedstr(candidate_refs_f1, rlen, elen, best_f1 + 1, 0, embeddedQ_f1);
  candidate_regions_f1[best_f1].embed_dist = nmismatch_f1;
  int nmismatch_r2 = embedstr(candidate_refs_r2, rlen, elen, best_r2 + 1, 0, embeddedQ_r2);
  candidate_regions_r2[best_r2].embed_dist = nmismatch_r2;
  thresholds[0] = nmismatch_f1 + nmismatch_r2;
  flag_f1[best_f1] = 0;
  flag_r2[best_r2] = 0;
  if (next_f1 != nregions_f1 && next_r2 != nregions_r2)
  {
    int nmismatch;
    if (next_f1 == best_f1)
    {
      assert(next_r2 != best_r2);
      nmismatch = embedstr(candidate_refs_r2, rlen, elen, next_r2 + 1, 0, embeddedQ_r2);
      candidate_regions_r2[next_r2].embed_dist = nmismatch;
      thresholds[1] = nmismatch_f1 + nmismatch;
    }
    else if (next_r2 == best_r2)
    {
      assert(next_f1 != best_f1);
      nmismatch = embedstr(candidate_refs_f1, rlen, elen, next_f1 + 1, 0, embeddedQ_f1);
      candidate_regions_f1[next_f1].embed_dist = nmismatch;
      thresholds[1] = nmismatch_r2 + nmismatch;
    }
    else
    {
      nmismatch = embedstr(candidate_refs_f1, rlen, elen, next_f1 + 1, 0, embeddedQ_f1);
      candidate_regions_f1[next_f1].embed_dist = nmismatch;
      thresholds[1] = nmismatch;
      nmismatch = embedstr(candidate_refs_r2, rlen, elen, next_r2 + 1, 0, embeddedQ_r2);
      candidate_regions_r2[next_r2].embed_dist = nmismatch;
      thresholds[1] += nmismatch;
    }
    flag_f1[next_f1] = 0;
    flag_r2[next_r2] = 0;
  }
}
void Embedding::embeddata_pair(vector<Region> &candidate_regions, char embeddedQ[],
                               const char **candidate_refs, unsigned ncandidates, bool flag[],
                               unsigned rlen, int threshold)
{
  auto start = std::chrono::system_clock::now();
  int step = 1;
  for (unsigned i = 0; i < ncandidates; i += step)
  {
    if (!flag[i])
      continue;
    candidate_regions[i].embed_dist = embedstr(candidate_refs, rlen, threshold, i + 1, 0, embeddedQ);
  }
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  embed_time += elapsed.count();
}
void Embedding::embeddata_iterative_update(vector<Region> &candidate_regions,
                                           const char **input, unsigned ninput, unsigned rlen,
                                           int &best_threshold, int &next_threshold,
                                           bool max_rnd, unsigned &best_idx, unsigned &next_idx)
{
  int step = 1;
  int elen = rlen * efactor;
  char embeddedQ[elen];
  auto start = std::chrono::system_clock::now();
  int nmismatch = 0;
  if (best_threshold == 0 && next_threshold == 0)
  {
    nmismatch = (memcmp(input[1], input[0], rlen) == 0 ? 0 : elen);
  }
  else
  {
    embedstr(input, rlen, elen, 0, 0, embeddedQ);
    nmismatch = embedstr(input, rlen, next_threshold, 1, 0, embeddedQ);
  }
  candidate_regions[0].embed_dist = nmismatch;
  if (nmismatch < best_threshold)
  {
    next_threshold = best_threshold;
    best_threshold = nmismatch;
  }
  else if (nmismatch < next_threshold)
  {
    next_threshold = nmismatch;
  }
  int best_dist = nmismatch;
  int next_dist = numeric_limits<int>::max();
  best_idx = next_idx = 0;
  for (unsigned i = 2; i < ninput; i += step)
  {
    if (best_threshold == 0 && next_threshold == 0)
    {
      nmismatch = (memcmp(input[i], input[0], rlen) == 0 ? 0 : elen);
    }
    else
    {
      nmismatch = embedstr(input, rlen, next_threshold, i, 0, embeddedQ);
    }
    candidate_regions[i - 1].embed_dist = nmismatch;
    if (nmismatch < best_threshold)
    {
      next_threshold = best_threshold;
      best_threshold = nmismatch;
    }
    else if (nmismatch < next_threshold)
    {
      next_threshold = nmismatch;
    }
    if (nmismatch < best_dist ||
        (nmismatch == best_dist && best_idx != 0 && candidate_regions[i - 1].rs < candidate_regions[best_idx].rs))
    {
      next_dist = best_dist;
      next_idx = best_idx;
      best_dist = nmismatch;
      best_idx = i - 1;
    }
    else if (nmismatch < next_dist)
    {
      next_dist = nmismatch;
      next_idx = i - 1;
    }
  }
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  embed_time += elapsed.count();
}
Embedding::Embedding()
{
#ifdef CGK2_EMBED
  efactor = 2;
#else
  efactor = 3;
#endif
  time_t seed = time(NULL);
  srand(seed);
  embed_time = 0;
  int loop = CGK_LOOP;
  string cPath(getcwd(NULL, 0));
  string file = cPath + "/hasheb/hasheb" + to_string(loop);
  fstream fin(file);
  if (fin.fail())
  {
    cout << file << "  hasheb.txt does not exist!" << endl;
  }
  for (int i = 0; i < loop; i++)
  {
    string tmp;
    getline(fin, tmp);
    for (int j = 0; j < NUM_STR; j++)
      for (int t = 0; t < 1; t++) 
        for (int d = 0; d < MAX_ELEN; d++)
          hash_eb[i][BITPOS(j, d, t)] = tmp.at(t * MAX_ELEN + d) - '0';
  }
}
Embedding::Embedding(const char *fname)
{
#ifdef CGK2_EMBED
  efactor = 2;
#else
  efactor = 3;
#endif
  cerr << "Loading embedding from file " << fname << endl;
  ifstream input(fname);
  if (input)
  {
    char num_str, num_char;
    input.read((char *)&num_str, sizeof(int));
    input.read((char *)&num_char, sizeof(int));
    assert(num_str == NUM_STR && num_char == NUM_CHAR);
    cerr << "NUM_STR: " << NUM_STR << ", NUM_CHAR: " << NUM_CHAR << " ,MAX_ELEN: " << MAX_ELEN << endl;
    string str_hasheb;
    std::getline(input, str_hasheb);
    for (int i = 0; i < TOTAL_RBITS; i++)
      hash_eb[i] = str_hasheb[i];
  }
}
Embedding::~Embedding()
{
}
