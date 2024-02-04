#ifndef SUFFIX_ARRAY_H
#define SUFFIX_ARRAY_H

#include <stdbool.h>
#include <limits.h>
#include <stdint.h>
#include <math.h>

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }
#define NULL_ELEM 0
#define INDEX_SIZE 8
#define INDEX_SIZE_N 87380 // = 4^1 + 4^2 + ... + 4^(INDEX_SIZE-1) + 4^INDEX_SIZE

typedef struct {
  int* sa;
  int* s;
  int* lcp;
  int** rmq;
  int* sar;
  int*** index;
  int n;
  int K;
  //char* str;
} SuffixArray;

extern int NA_ENC[UCHAR_MAX];
extern int IUPAC[15][4];

extern int* get_na_enc();
extern void na_enc_init();
extern int encode_str(
  const char* str, int* encoding, int* s);
extern void suffix_array_init(
  const char* str, int* encoding, SuffixArray* sa);
extern void delete_suffix_array(
  SuffixArray* sa);
extern void print_sa_lcd(
  SuffixArray* sa, const char* str);
extern void radixPass(
  int* a, int* b, int* r, int n, int K);
extern void skew(
  int* s, int* SA, int n, int K);
extern void kasai(
  int* s, int* SA, int n, int* lcp, int* invertedSuffixArray);
extern void construct_rmq_array(
  int* lcp, int n, int** rmq);
extern int encoding_to_index(
  int* encoding, int n, int index_size);
extern void create_index(
  int* s, int* SA, int n, int* lcp, int*** index);
extern int find(
  int* q, int qlen, int* SA, int* lcp, int* s, int n, 
  int** rmq, int* res);
extern int find_all(
  int* q, int qlen, int* SA, int* lcp, int* s, int n, 
  int** rmq, int* res);
extern int find_all_indexed(
  int* q, int qlen, int* SA, int* lcp, int* s, int n, 
  int*** index, int* res);
extern int find_all_bipartite(
  int* q1, int q1len, int g, int* q2, int q2len, 
  int* SA, int* lcp, int* s, int n, int** rmq,
  int** idxs);
extern void get_indices(
  int k, int n, int* SA, int* indices);
extern void get_indices_from_bipartite_search(
  int* idxs, int N, int* indices);
extern int get_indices_bipartite(
  int k, int N, int* SA, int* s, int n, int q1len, 
  int g, int* q2, int q2len, 
  int* indices);
extern int find_motif(
  const char* motif, int mlen, int* SA, int* lcp, int* s, int n, int** rmq,
  int** res);
extern int find_motif_nonparallel_indexed(
  const char* motif, int mlen,
  int* SA, int* lcp, int* s, int n, int*** index, int* SAr,
  int** res);
//extern void motif_means(
//  const char* motifs_data, int motif_count, int max_mlen,
//  float* fwd, float* rev,
//  int* SA, int* lcp, int* s, int n, int** rmq,
//  float* mean_data, int* count_data);
extern void motif_means(
  const char* motifs_data, int motif_count, int max_mlen,
  float* fwd, float* rev,
  int* SA, int* lcp, int* s, int n, 
  int** rmq, int*** index, int* SAr,
  float* mean_data, int* count_data);
extern void motif_medians(
  const char* motifs_data, int motif_count, int max_mlen,
  float* fwd, float* rev,
  int* SA, int* lcp, int* s, int n, 
  int** rmq, int*** index, int* SAr,
  float* median_data, int* count_data);
extern float quick_select_median(
  float arr[], uint32_t n);
extern void write_sa(
  SuffixArray* sa, char* fn);
extern void read_sa(
  char* fn, int n, SuffixArray* sa);

#endif