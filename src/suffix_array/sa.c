#include "sa.h"
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <omp.h>
#include <stdio.h>

int comp (const void * elem1, const void * elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

// nucleic acid encoding
int NA_ENC[UCHAR_MAX];
int IUPAC[15][4];
void na_enc_init() {
  for (int i=0; i <= UCHAR_MAX; i++) {
    switch((char)i) {
      case 'A':
      case 'a': 
        NA_ENC[i] = 1; break;
      case 'C':
      case 'c': 
        NA_ENC[i] = 2; break;
      case 'G':
      case 'g': 
        NA_ENC[i] = 3; break;
      case 'T':
      case 't':
      case 'U':
      case 'u': 
        NA_ENC[i] = 4; break;
      default: NA_ENC[i] = NULL_ELEM; break;
    }
  }

  int iupac[15][4] = {
    {1,NULL_ELEM,NULL_ELEM,NULL_ELEM}, //  0 : A
    {2,NULL_ELEM,NULL_ELEM,NULL_ELEM}, //  1 : C
    {3,NULL_ELEM,NULL_ELEM,NULL_ELEM}, //  2 : G
    {4,NULL_ELEM,NULL_ELEM,NULL_ELEM}, //  3 : T,U
    {1,2,NULL_ELEM,NULL_ELEM}, //  4 : M  (A/C)
    {1,3,NULL_ELEM,NULL_ELEM}, //  5 : R  (A/G)
    {1,4,NULL_ELEM,NULL_ELEM}, //  6 : W  (A/T)
    {2,3,NULL_ELEM,NULL_ELEM}, //  7 : S  (C/G)
    {2,4,NULL_ELEM,NULL_ELEM}, //  8 : Y  (C/T)
    {3,4,NULL_ELEM,NULL_ELEM}, //  9 : K  (G/T)
    {1,2,3,NULL_ELEM}, // 10 : V  (A/C/G)
    {1,2,4,NULL_ELEM}, // 11 : H  (A/C/T)
    {1,3,4,NULL_ELEM}, // 12 : D  (A/G/T)
    {2,3,4,NULL_ELEM}, // 13 : B  (C/G/T)
    {1,2,3,4}  // 14 : N  (A/C/G/T)
  };
  int* p = &IUPAC[0][0];
  int* q = &iupac[0][0];
  memcpy(p, q, sizeof(*p) * 15 * 4);
}

void reverse_complement(char* motif, int mlen, char* motif_rc)
{
  // "ACGTMRWSYKVHDBN", "TGCAKYWSRMBDHVN"
  for (int i=0; i<mlen; i++) {
    switch (motif[i]) {
      case 'A' : motif_rc[mlen - i - 1] = 'T'; break;
      case 'C' : motif_rc[mlen - i - 1] = 'G'; break;
      case 'G' : motif_rc[mlen - i - 1] = 'C'; break;
      case 'T' : motif_rc[mlen - i - 1] = 'A'; break;
      case 'M' : motif_rc[mlen - i - 1] = 'K'; break;
      case 'R' : motif_rc[mlen - i - 1] = 'Y'; break;
      case 'W' : motif_rc[mlen - i - 1] = 'W'; break;
      case 'S' : motif_rc[mlen - i - 1] = 'S'; break;
      case 'Y' : motif_rc[mlen - i - 1] = 'R'; break;
      case 'K' : motif_rc[mlen - i - 1] = 'M'; break;
      case 'V' : motif_rc[mlen - i - 1] = 'B'; break;
      case 'H' : motif_rc[mlen - i - 1] = 'D'; break;
      case 'D' : motif_rc[mlen - i - 1] = 'H'; break;
      case 'B' : motif_rc[mlen - i - 1] = 'V'; break;
      case 'N' : motif_rc[mlen - i - 1] = 'N'; break;
    }
  }
  motif_rc[mlen] = '\0';
}

int* get_na_enc()
{
  return NA_ENC;
}

int encode_str(const char* str, int* encoding, int* s)
{
  int K = 0; // max encoded character value
  if (encoding == NULL) {
    for (uint32_t i=0; i < strlen(str); i++) {
      s[i] = (uint)str[i];
      if (K < s[i]) 
        K = s[i];
    }
  } else {
    for (int i=0; i <= UCHAR_MAX; i++) {
      if (K < encoding[i])
        K = encoding[i];
    }
#pragma omp parallel for if(strlen(str) > 10000)
    for (uint32_t i=0; i < strlen(str); i++) {
      s[i] = encoding[(uint)str[i]];
      if (s[i] == NULL_ELEM) {
        // character is not valid (non-ribonucleine character)
        // TODO: handle error
      }
    }
  }
  return K;
}

// initialize the SuffixArray given a string
void suffix_array_init(const char* str, int* encoding, SuffixArray* sa)
{
  int n = strlen(str) + 1;
  sa->sa = (int*)malloc(sizeof(int) * n);
  sa->s = (int*)malloc(sizeof(int) * (n + 3));
  sa->lcp = (int*)malloc(sizeof(int) * n);
  int h = floor(log((double)n) / log(2.)) + 1;
  sa->rmq = (int**)malloc(sizeof(int*) * (h+1) + sizeof(int) * n * h);
  sa->sar = (int*)malloc(sizeof(int) * n);
  sa->index = (int***)malloc(sizeof(int**) * INDEX_SIZE + sizeof(int*) * INDEX_SIZE_N + sizeof(int) * 2 * INDEX_SIZE_N);
  sa->n = n;
  sa->K = 0;

  sa->K = encode_str(str, encoding, sa->s);
  sa->s[strlen(str)] = sa->s[n] = sa->s[n+1] = sa->s[n+2] = 0;
}

void delete_suffix_array(SuffixArray* sa)
{
  free(sa->sa);
  free(sa->s);
  free(sa->lcp);
  free(sa->rmq);
  free(sa->sar);
  free(sa->index);
}

// lexicographic order for pairs
//inline
bool leq2(int a1, int a2, int b1, int b2) {
  return(a1 < b1 || (a1 == b1 && a2 <= b2));
}

// and triples
//inline
bool leq3(int a1, int a2, int a3, int b1, int b2, int b3) {
  return(a1 < b1 || (a1 == b1 && leq2(a2,a3, b2,b3)));
} // and triples

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
void radixPass(int* a, int* b, int* r, int n, int K) {// count occurrences
    //int c[K + 1]; // counter array
    int* c = (int*)malloc(sizeof(int) * (K + 1));
    for (int i = 0; i <= K; i++) c[i] = 0; // reset counters
    for (int i = 0; i < n; i++) c[r[a[i]]]++; // count occurrences
    for (int i = 0, sum = 0; i <= K; i++) { // exclusive prefix sums
      int t = c[i];
      c[i] = sum;
      sum += t;
    }
    for (int i = 0;  i < n; i++) b[c[r[a[i]]]++] = a[i]; // sort
    free(c);
}

// find the suffix array SA of s[0..n-1] in {1..K}ˆn
// require s[n]=s[n+1]=s[n+2]=0, n>=2
// taken from Kärkkäinen, Sanders https://www.cs.helsinki.fi/u/tpkarkka/publications/icalp03.pdf
void skew(int* s, int* SA, int n, int K) {
  int n0 = (n+2)/3, n1 = (n+1)/3, n2 = n/3, n02 = n0+n2;
  int* s12 = (int*)malloc(sizeof(int) * (n02+3)); // malloc necessary to store in heap instead of stack
  s12[n02] = s12[n02+1] = s12[n02+2] = NULL_ELEM;
  int* SA12;
  SA12 = (int*)malloc(sizeof(int) * (n02+3));
  SA12[n02] = SA12[n02+1] = SA12[n02+2] = NULL_ELEM;
  int* s0;
  s0 = (int*)malloc(sizeof(int) * (n0));
  int* SA0;
  SA0 = (int*)malloc(sizeof(int) * (n0));
  // generate positions of mod 1 and mod 2 suffixes
  // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
  for (int i=0, j=0; i < n + (n0-n1); i++)
    if (i%3 != 0) s12[j++] = i;
  // lsb radix sort the mod 1 and mod 2 triples
  radixPass(s12 , SA12, s+2, n02, K);
  radixPass(SA12, s12 , s+1, n02, K);
  radixPass(s12 , SA12, s  , n02, K);
  // find lexicographic names of triples
  int name = 0, c0 = -1, c1 = -1, c2 = -1;
  for (int i = 0; i < n02; i++) {
    if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) {
      name++;
      c0 = s[SA12[i]];
      c1 = s[SA12[i]+1];
      c2 = s[SA12[i]+2];
    }
    if (SA12[i]%3 == 1) s12[SA12[i]/3] = name; // left half
    else s12[SA12[i]/3 + n0] = name; // right half
  }
  // recurse if names are not yet unique
  if (name < n02) {
    skew(s12, SA12, n02, name);
    // store unique names in s12 using the suffix array
    for (int i = 0; i < n02; i++) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (int i = 0;  i < n02; i++) SA12[s12[i] - 1] = i;
  // stably sort the mod 0 suffixes from SA12 by their first character
  for (int i = 0, j = 0; i < n02; i++)
    if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  radixPass(s0, SA0, s, n0, K);
  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (int p = 0, t = n0-n1, k = 0; k < n; k++) {
    #define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    int i = GetI(); // pos of current offset 12 suffix
    int j = SA0[p]; // pos of current offset 0 suffix
    if (SA12[t] < n0 ? // different compares for mod 1 and mod 2 suffixes
      leq2(s[i], s12[SA12[t] + n0], s[j], s12[j/3]) :
      leq3(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
    {// suffix from SA12 is smaller
      SA[k] = i; t++;
      if (t == n02) // done --- only SA0 suffixes left
      for (k++; p < n0; p++, k++) SA[k] = SA0[p];
    } else {// suffix from SA0 is smaller
      SA[k] = j; p++;
      if (p == n0) // done --- only SA12 suffixes left
      for (k++; t < n02; t++, k++) SA[k] = GetI();
    }
  }
  free(SA12); free(SA0); free(s0); free(s12);
  //delete [] s12; delete [] SA12; delete [] SA0; delete [] s0;
}

// make lcp array using the Kasai algorithm
void kasai(int* s, int* SA, int n, int* lcp, int* invertedSuffixArray)
{
  //int* invertedSuffixArray;
  //invertedSuffixArray = (int*)malloc(sizeof(int) * n);
	//Construct Inverted Suffix Array
	for (int i = 0; i < n; i++)
		invertedSuffixArray[SA[i]] = i;

	//Kasai Algorithim
  lcp[0] = 0;
  int l = 0;
	for (int i=0; i < n - 1; i++)
	{
    int k = invertedSuffixArray[i];
    if (k != 0) {
      int j = SA[k - 1];
      while (s[i + l] == s[j + l]) //  && (i + l) < (n - 1) && (j + l) < (n - 1)
			  l++;
      lcp[k] = l;
      if (l > 0)
			  l--;
    }
	}
  //free(invertedSuffixArray);
}

void print_sa_lcd(SuffixArray* sa, const char* str)
{
  // display sa & lca in tabular form
  printf("%4s %4s %4s %s\n", "i", "sa", "lcd", "suffix");
  for (int i=0; i<(15+sa->n); i++) {printf("-");} printf("\n");
  printf("%4d %4d %4d %s -\n", 0, sa->sa[0], sa->lcp[0], &str[sa->sa[0]]);
  for (int i=1; i < sa->n; i++) {
    printf("%4d %4d %4d %s\n", i, sa->sa[i], sa->lcp[i], &str[sa->sa[i]]);
  }
}

// O(n+log(n)) construction & space complexity
// https://doi.org/10.1016/j.jalgor.2005.08.001
// O(n) possible, TODO: implement
void construct_rmq_array(int* lcp, int n, int** M)
{
  int h = floor(log((double)n) / log(2.)) + 1;
  // column major
  int* data = (int*)(M + h);
  for (int j=0; j<h; j++) {
    M[j] = data + (j * n);
  }

  // fill array
  for (int i=0; i<n; i++){
    int j = 0;
    M[j][i] = lcp[i];
  }
  for (int j=1; j<=h; j++){
    for (int i=0; i + floor(exp2(j)) < n; i++){ // (1 << j)
      int b = i + floor(exp2(j-1)); //(1 << (j - 1));
      //printf("j = %d, i = %d, b = %d\n", j,i,b);
      if (M[j-1][i] <= M[j-1][b]) {
        M[j][i] = M[j-1][i];
      } else {
        M[j][i] = M[j-1][b];
      }
    }
  }
}

int arr_lcp(int* a, int* b, int n) 
{
  int k = 0;
  for (int i=0; i<n ; i++) {
    if (a[i] != b[i]) break;
    k++;
  }
  return k;
}

int encoding_to_index(int* encoding, int n, int index_size)
{
  int res = encoding[0] - 1;
  for (int i=1; i<n; i++) {
    res = res << 2;
    res |= (encoding[i] - 1);
  }
  res = res << (index_size - n) * 2;
  return res;
}

void create_index(int* s, int* SA, int n, int* lcp, int*** index)
{
  /*
  int* data = (int*)(index + INDEX_SIZE_N);
  for (int i=0; i<INDEX_SIZE_N; i++) {
    index[i] = data + (i * 2);
  }
  */

  int* data = (int*)(index + INDEX_SIZE + INDEX_SIZE_N);
  for (int i=0; i<INDEX_SIZE; i++) {
    (index)[i] = (int**)(index) + INDEX_SIZE;
    int offset = 0;
    for (int j=0; j<i; j++) {
      offset += (int)pow(4, j+1);
    }
    (index)[i] += offset;
    //printf("%d: %d %d\n", i, offset, (int)pow(4, i+1));
    for (int j=0; j<(int)pow(4, i+1); j++) {
      (index)[i][j] = data + 2 * offset + (j * 2);
    }
  }
  // set defaults (for encoding not in genomic sequence)
  for (int i=0; i<INDEX_SIZE; i++) {
    //printf("%d: %d\n", i, (int)pow(4, i+1));
    for (int j=0; j<(int)pow(4, i+1); j++) {
      //printf("%d\n", j);
      index[i][j][0] = n;//(n-1);
      index[i][j][1] = n;//(n-1);
    }
  }
  
  //printf("go\n");
  // go through suffix array and assign start / end points for encodings
  int idx;
  for (int j=0; j<(n-1 < INDEX_SIZE ? n-1 : INDEX_SIZE); j++) {
    int l = j + 1;
    //printf("\nl = %d\n", l);
    int i=1; // first elem is empty suffix
    int suffix_len;
    for (; i<n; i++) {
      suffix_len = n - SA[i] - 1;
      if (suffix_len >= l)
        break;
    }
    //printf("%d\n", i);
    //printf("%d\n", SA[i]);
    //printf("%d\n", n);
    idx = encoding_to_index(&s[SA[i]], l, l);
    //printf("%#04X\n", idx);
    //int start = i;
    index[j][idx][0] = i;
    ////printf("%d\n", i);
    for (; i<n; i++) {
      //printf("i = %d\n", i);
      suffix_len = n - SA[i] - 1;
      //printf("suffix_len = %d\n", suffix_len);
      //printf("lcp[i] = %d\n", lcp[i]);
      if (suffix_len < l) {
        index[j][idx][1] = i;
        while ((suffix_len < l) && (i < n-1)) {
          i++;
          suffix_len = n - SA[i] - 1;
          //printf("*suffix_len = %d\n", suffix_len);
        }
        if (i == n)
          break;
        //printf("SA[i] = %d\n", SA[i]);
        idx = encoding_to_index(&s[SA[i]], l, l);
        //printf("%#04X\n", idx);
        index[j][idx][0] = i;
      } else if (lcp[i] < l) {
        index[j][idx][1] = i;
        //printf("SA[i] = %d\n", SA[i]);
        idx = encoding_to_index(&s[SA[i]], l, l);
        //printf("%#04X\n", idx);
        index[j][idx][0] = i;
      }
    }
    //index[idx][1] = n - 1;
  }
  /*/ print
  for (int idx=0; idx<INDEX_SIZE_N; idx++) {
    if (index[idx][0] != (n-1)) {
      printf("%#04X : (%d, %d)\n", idx, index[idx][0], index[idx][1]);
      for (int j=index[idx][0]; j<index[idx][1]; j++) {
        for (int c=0; c<INDEX_SIZE; c++) {
          printf("%d ", s[SA[j]+c]);
        }
        printf("\n");
      }
    }
  }
  
  //*/
  //printf("done\n");
}

int binary_search_r(int* q, int qlen, int* SA, int* lcp, int* s, int** rmq,
                     int l, int mp, int r, int m, int k,
                     int* res) 
{
  //printf("%d\n", mp);

  // find lcp between the previous and the current suffix
  // --> use precomputed RMQ for O(1) complexity
  int kp;
  if (mp == m) {
    kp = qlen;
  } else {
    int l_, r_;
    if (mp < m) {
      l_ = mp + 1;
      r_ = m;
    } else {
      l_ = m + 1;
      r_ = mp;
    }
    int h, a, b;
    h = floor(log2((double)r_ - l_ + 1));
    a = rmq[h][l_];
    b = rmq[h][r_ - (1 << h) + 1];
    kp = a < b ? a : b;
    if (kp > qlen)
      kp = qlen;
  }
  //printf("kp / rmq : %2d / %2d, k : %2d\n", kp, kp_rmq, k);

  // As in https://de.wikipedia.org/wiki/LCP-Array
  if (k == kp) { // k == kp
    int kmp = k + arr_lcp(&q[k], &s[SA[mp] + k], qlen - k); // will become the new k
    if (kmp == qlen) {
      // found
      *res = mp;
      return 1;
    } else if (q[kmp] < s[SA[mp] + kmp]) {
      // search left
      return binary_search_r(q, qlen, SA, lcp, s, rmq,
                      l, (l + (mp - l) / 2), mp, mp, kmp, res);
    } else {
      return binary_search_r(q, qlen, SA, lcp, s, rmq,
                      mp, (mp + (r - mp) / 2), r, mp, kmp, res);
    }
  } else if (k < kp) {
    // search is continued in the same direction
    if (mp == m) {
      // query not in s
      return 0;
    } else if (mp < m) {
      return binary_search_r(q, qlen, SA, lcp, s, rmq,
                      l, (l + (mp - l) / 2), mp, mp, k, res);
    } else {
      return binary_search_r(q, qlen, SA, lcp, s, rmq,
                      mp, (mp + (r - mp) / 2), r, mp, k, res);
    }
  } else { // k > kp
    // search is continued in the opposite direction
    if (mp == m) {
      // query not in s
      return 0;
    } else if (mp < m) {
      return binary_search_r(q, qlen, SA, lcp, s, rmq,
                      mp, (mp + (r - mp) / 2), r, mp, kp, res);
    } else {
      return binary_search_r(q, qlen, SA, lcp, s, rmq,
                      l, (l + (mp - l) / 2), mp, mp, kp, res);
    }
  }
}

int find(int* q, int qlen, int* SA, int* lcp, int* s, int n, int** rmq,
                   int* res)
{
  if (qlen > n) {
    // query cannot be longer than string that is searched
    return 0;
  } else if (qlen == 0) {
    // empty suffix is always at index 0 and always valid for qlen = 0
    *res = 0;
    return 1;
  }
  int l, m, r; 
  l = 0;
  r = n;
  m = l + (n - l) / 2;
  //printf("%d\n", m);
  int k = arr_lcp(q, &s[SA[m]], qlen);
  if (k == qlen) {
    // found
    *res = m;
    return 1;
  } else if (q[k] < s[SA[m] + k]) {
    return binary_search_r(q, qlen, SA, lcp, s, rmq,
                    l, (l + (m - l) / 2), m, m, k,
                    res);
  } else {
    return binary_search_r(q, qlen, SA, lcp, s, rmq,
                    m, (m + (r - m) / 2), r, m, k,
                    res);
  }
}

int find_all(int* q, int qlen, int* SA, int* lcp, int* s, int n, int** rmq,
              int* res)
{
  int any_hit = find(q, qlen, SA, lcp, s, n, rmq, res);
  if (any_hit == 0) {
    return 0;
  }
  // query was found, now find all suffixes that have q as their prefix
  int initial_hit = *res;
  int l, r;
  l = initial_hit; 
  r = initial_hit + 1;
  while (l > 0) {
    if (lcp[l] >= qlen) {
      l--;
    } else {
      break;
    }
  }
  while (r < n) {
    if (lcp[r] >= qlen) {
      r++;
    } else {
      break;
    }
  }
  *res = l;
  return r - l;
}

int find_all_indexed(
  int* q, int qlen, int* SA, int* lcp, int* s, int n, 
  int*** index, int* res)
{
  int idx = encoding_to_index(q, qlen, qlen);
  *res = index[qlen-1][idx][0];
  return index[qlen-1][idx][1] - index[qlen-1][idx][0];
}

int find_all_bipartite(int* q1, int q1len, int g, int* q2, int q2len, 
                       int* SA, int* lcp, int* s, int n, int** rmq,
                       int** idxs)
{
  int N1, i1, N2, i2;
  int* indices1;
  int* indices2;
#pragma omp parallel num_threads(2)
{
#pragma omp single
{
#pragma omp task
{
  i1 = 0;
  N1 = find_all(q1, q1len, SA, lcp, s, n, rmq, &i1);
}
#pragma omp task
{
  i2 = 0;
  N2 = find_all(q2, q2len, SA, lcp, s, n, rmq, &i2);
}
#pragma omp taskwait
  if (N1 > 0 && N2 > 0) {
#pragma omp task
{
    indices1 = (int*)malloc(sizeof(int) * N1);
    get_indices(i1, N1, SA, indices1);
    qsort(indices1, N1, sizeof(int), comp);
}
#pragma omp task
{
    indices2 = (int*)malloc(sizeof(int) * N2);
    get_indices(i2, N2, SA, indices2);
    qsort(indices2, N2, sizeof(int), comp);
}
  }
}
}
  if (N1 == 0 || N2 == 0) return 0;

  // identify bipartite hits
  int i, j, k;
  i = 0; j = 0; k = 0;
  int* indices;
  indices = (int*)malloc(sizeof(int) * ((N1 < N2) ? N1 : N2));
  while ((i < N1) && (j < N2)) {
    if (indices1[i] < indices2[j]) {
      if ((indices2[j] - indices1[i]) == (q1len + g)) {
        indices[k] = indices1[i];
        k++;
        i++;
        j++;
      } else if ((indices2[j] - indices1[i]) > (q1len + g)) {
        i++;
      } else {
        j++;
      }
    } else {
      j++;
    }
  }
  free(indices1); free(indices2);
  *idxs = indices;
  return k;
}

void get_indices(int k, int N, int* SA, int* indices)
{
  for (int i=0; i<N; i++){
    indices[i] = SA[k+i];
  }
}

void get_indices_from_bipartite_search(int* idxs, int N, int*indices)
{
  memcpy(indices, idxs, sizeof(int) * N);
  free(idxs);
}

int get_indices_bipartite(int k, int N, int* SA, int* s, int n, int q1len, int g, int* q2, int q2len, int* indices)
{
  // filter for indices with a mathing sub-query q2 that follows after a gap of length g
  int new_N = 0;
  if (g >= 0) { // search q2 DOWNSTREAM of q1
    int max_loc = n - (q1len + g + q2len) - 1;
#pragma omp parallel for if(N > 100)
    for (int j=k; j<(k+N); j++) {
      if (SA[j] <= max_loc) {
        int o = arr_lcp(q2, &(s[SA[j] + q1len + g]), q2len);
        if (o == q2len) {
#pragma omp critical
{
          indices[new_N] = SA[j];
          new_N++;
}
        }
      }
    }
  } else { // search q2 UPSTREAM of q1
    int max_loc = q2len - g;
#pragma omp parallel for if(N > 100)
    for (int j=k; j<(k+N); j++) {
      if (SA[j] >= max_loc) {
        int o = arr_lcp(q2, &(s[SA[j] + g - q2len]), q2len);
        if (o == q2len) {
#pragma omp critical
{
          indices[new_N] = SA[j];
          new_N++;
}
        }
      }
    }
  }
  return new_N;
}

void encode_IUPAC(const char* motif, int mlen, int n_enc, int* factors, int** encodings)
{
  int blocksize = n_enc;
  int subblocksize;
  //int* base_enc = NULL;
  int base_enc;

  for (int i=0; i<mlen; i++) {
    subblocksize = blocksize/factors[i];
    switch(motif[i]) {
      case 'A': base_enc = 0; break;
      case 'C': base_enc = 1; break;
      case 'G': base_enc = 2; break;
      case 'T':
      case 'U': base_enc = 3; break;
      case 'M': base_enc = 4; break;
      case 'R': base_enc = 5; break;
      case 'W': base_enc = 6; break;
      case 'S': base_enc = 7; break;
      case 'Y': base_enc = 8; break;
      case 'K': base_enc = 9; break;
      case 'V': base_enc = 10; break;
      case 'H': base_enc = 11; break;
      case 'D': base_enc = 12; break;
      case 'B': base_enc = 13; break;
      case 'N': base_enc = 14; break;
    }
    for (int j=0; j<n_enc; j+=blocksize) {
      for (int k=0; k<factors[i];k++) {
        for (int l=0; l<subblocksize; l++) {
          //printf("%d %d %d\n", j + k*subblocksize + l, i, IUPAC[base_enc][k]);
          encodings[j + k*subblocksize + l][i] = IUPAC[base_enc][k];
        }
      }
    }
    blocksize /= factors[i];
  }
  #ifdef DEBUG
  // print encodings
  for (int i=0; i<n_enc; i++) {
    printf("encoding #%3d : ", i);
    for (int j=0; j<mlen; j++) {
      printf("%d ", encodings[i][j]);
    }
    printf("\n");
  }
  #endif
}

int get_number_of_IUPAC_encodings(const char* motif, int mlen, int* factors)
{
  // determine number of encodings
  int n_enc = 1;
  for (int i=0; i<mlen; i++) {
    switch (motif[i]) {
      case 'R':
      case 'Y':
      case 'S':
      case 'W':
      case 'K':
      case 'M':
        n_enc *= 2;
        factors[i] = 2;
        break;
      case 'B':
      case 'D':
      case 'H':
      case 'V':
        n_enc *= 3;
        factors[i] = 3;
        break;
      case 'N':
        n_enc *= 4;
        factors[i] = 4;
        break;
      default:
        factors[i] = 1;
        break;
    }
  }
  return n_enc;
}

int find_IUPAC(const char* motif, int mlen, 
               int* SA, int* lcp, int* s, int n, int** rmq,
               int** Is, int** Ns)
{
  int factors[mlen];
  int n_enc = get_number_of_IUPAC_encodings(motif, mlen, factors);
  // define encodings array
  int** encodings = (int**)malloc(sizeof(int*) * n_enc + sizeof(int) * (n_enc * mlen));
  int* data = (int*)(encodings + n_enc);
  for (int i=0; i<n_enc; i++) {
    encodings[i] = data + (i * mlen);
  }
  // create encodings
  encode_IUPAC(motif, mlen, n_enc, factors, encodings);
  // search all encodings
  *Is = (int*)malloc(sizeof(int) * n_enc);
  *Ns = (int*)malloc(sizeof(int) * n_enc);
  for (int i=0; i<n_enc; i++) {
    (*Ns)[i] = find_all(encodings[i], mlen, SA, lcp, s, n, rmq, &((*Is)[i]));
  }
  // free memory
  free(encodings);
  return n_enc;
}

int find_IUPAC_indexed(
  const char* motif, int mlen, 
  int* SA, int* lcp, int* s, int n, int*** index,
  int** Is, int** Ns)
{
  int factors[mlen];
  int n_enc = get_number_of_IUPAC_encodings(motif, mlen, factors);
  // define encodings array
  int** encodings = (int**)malloc(sizeof(int*) * n_enc + sizeof(int) * (n_enc * mlen));
  int* data = (int*)(encodings + n_enc);
  for (int i=0; i<n_enc; i++) {
    encodings[i] = data + (i * mlen);
  }
  // create encodings
  encode_IUPAC(motif, mlen, n_enc, factors, encodings);
  // search all encodings
  *Is = (int*)malloc(sizeof(int) * n_enc);
  *Ns = (int*)malloc(sizeof(int) * n_enc);
  for (int i=0; i<n_enc; i++) {
    (*Ns)[i] = find_all_indexed(encodings[i], mlen, SA, lcp, s, n, index, &((*Is)[i]));
  }
  // free memory
  free(encodings);
  return n_enc;
}

bool IUPAC_match(const char* iupac, int* enc, int n)
{
  for (int i=0; i<n; i++) {
    switch (iupac[i]) {
      case 'A': if (enc[i] != 1) return false; break;
      case 'C': if (enc[i] != 2) return false; break;
      case 'G': if (enc[i] != 3) return false; break;
      case 'T': 
      case 'U': if (enc[i] != 4) return false; break;
      case 'M': if ((enc[i] != 1) && (enc[i] != 2)) return false; break;
      case 'R': if ((enc[i] != 1) && (enc[i] != 3)) return false; break;
      case 'W': if ((enc[i] != 1) && (enc[i] != 4)) return false; break;
      case 'S': if ((enc[i] != 2) && (enc[i] != 3)) return false; break;
      case 'Y': if ((enc[i] != 2) && (enc[i] != 4)) return false; break;
      case 'K': if ((enc[i] != 3) && (enc[i] != 4)) return false; break;
      case 'V': if ((enc[i] != 1) && (enc[i] != 2) && (enc[i] != 3)) return false; break;
      case 'H': if ((enc[i] != 1) && (enc[i] != 2) && (enc[i] != 4)) return false; break;
      case 'D': if ((enc[i] != 1) && (enc[i] != 3) && (enc[i] != 4)) return false; break;
      case 'B': if ((enc[i] != 2) && (enc[i] != 3) && (enc[i] != 4)) return false; break; 
      case 'N': break;
    }
  }
  return true;
}

// find indices of the given motif
// motif given in IUPAC must be uppercase and not contain terminal N characters
int find_motif(const char* motif, int mlen,
               int* SA, int* lcp, int* s, int n, int** rmq,
               int** res)
{
  int N_total = 0;
  // determine best motif-search strategy
  // find longest gap (stretch of N)
  int g, gloc, in_gap; // longest gap
  g = 0; gloc = 0; in_gap = 0;
  for (int i=1; i<mlen; i++) {
    if (motif[i] == 'N') {
      in_gap++;
    } else if (in_gap > 0) {
      if (in_gap >= g) {
        g = in_gap;
        gloc = i - g;
      }
      in_gap = 0;
    }
  }
  if (g <= 2) {
    // search as contiguous motif
    int* Is;
    int* Ns;
    int n_enc = find_IUPAC(motif, mlen, SA, lcp, s, n, rmq,
                           &Is, &Ns);
    for (int i=0; i<n_enc; i++) {
      N_total += Ns[i];
    }
    int* indices = (int*)malloc(sizeof(int) * N_total);
    int N_curr = 0;
    for (int i=0; i<n_enc; i++) {
      get_indices(Is[i], Ns[i], SA, &indices[N_curr]);
      //#ifdef DEBUG
      //printf("hits for encoding %d (n = %d):\n", i, Ns[i]);
      //for (int j=0; j<Ns[i]; j++) {
      //  printf("#%3d : SA[%3d] = %3d, substring: ", j, Is[i] + j, SA[Is[i] + j]);
      //  for (int k=0; k<mlen; k++) {
      //    printf("%d ", s[SA[Is[i] + j] + k]);
      //  }
      //  printf("\n");
      //}
      //#endif
      N_curr += Ns[i];
    }
    *res = indices;
    free(Is); free(Ns);
  } else {
    // search as bipartite motif
    // search the longer one of the two parts
    int* Is;
    int* Ns;
    const char* q1; const char* q2;
    int q1len, q2len;
    int offset = 0;
    if (gloc >= mlen - (gloc + g)) {
      // search for prefix
      q1 = motif;
      q1len = gloc;
      q2 = &(motif[gloc + g]);
      q2len = mlen - (gloc + g);
    } else {
      // search for suffix
      offset = -(gloc + g);
      q1 = &(motif[gloc + g]);
      q1len = mlen - (gloc + g);
      q2 = motif;
      q2len = gloc;
      g = -g;
    }
    int n_enc = find_IUPAC(q1, q1len, SA, lcp, s, n, rmq,
                           &Is, &Ns);
    for (int i=0; i<n_enc; i++) {
      N_total += Ns[i];
    }
    int* indices = (int*)malloc(sizeof(int) * N_total);
    
    int new_N = 0;
#pragma omp parallel
{
    for (int i=0; i<n_enc; i++) {
      if (g >= 0) { // search q2 DOWNSTREAM of q1
        int max_loc = n - (q1len + g + q2len) - 1;
#pragma omp for
        for (int j=Is[i]; j<(Is[i]+Ns[i]); j++) {
          if (SA[j] <= max_loc) {
            //int o = arr_lcp(q2, &(s[SA[j] + q1len + g]), q2len);
            if (IUPAC_match(q2, &(s[SA[j] + q1len + g]), q2len)) {
#pragma omp critical
{
              indices[new_N] = SA[j];
              new_N++;
}
            }       
          }
        }
      } else { // search q2 UPSTREAM of q1
        int max_loc = q2len - g;
#pragma omp for
        for (int j=Is[i]; j<(Is[i]+Ns[i]); j++) {
          if (SA[j] >= max_loc) {
            //int o = arr_lcp(q2, &(s[SA[j] + g - q2len]), q2len);
            if (IUPAC_match(q2, &(s[SA[j] + g - q2len]), q2len)) {
#pragma omp critical
{
              indices[new_N] = SA[j] + offset;
              new_N++;
}
            }
          }         
        }
      }
    }
}
    indices = (int*)realloc(indices, sizeof(int) * new_N);
    *res = indices;
    N_total = new_N;
    free(Is); free(Ns);
  }

  return N_total;
}

int find_motif_nonparallel(const char* motif, int mlen,
                           int* SA, int* lcp, int* s, int n, int** rmq,
                           int** res)
{
  int N_total = 0;
  // determine best motif-search strategy
  // find longest gap (stretch of N)
  int g, gloc, in_gap; // longest gap
  g = 0; gloc = 0; in_gap = 0;
  for (int i=1; i<mlen; i++) {
    if (motif[i] == 'N') {
      in_gap++;
    } else if (in_gap > 0) {
      if (in_gap >= g) {
        g = in_gap;
        gloc = i - g;
      }
      in_gap = 0;
    }
  }
  if (g <= 2) {
    // search as contiguous motif
    int* Is;
    int* Ns;
    int n_enc = find_IUPAC(motif, mlen, SA, lcp, s, n, rmq,
                           &Is, &Ns);
    for (int i=0; i<n_enc; i++) {
      N_total += Ns[i];
    }
    int* indices = (int*)malloc(sizeof(int) * N_total);
    int N_curr = 0;
    for (int i=0; i<n_enc; i++) {
      get_indices(Is[i], Ns[i], SA, &indices[N_curr]);
      //#ifdef DEBUG
      //printf("hits for encoding %d (n = %d):\n", i, Ns[i]);
      //for (int j=0; j<Ns[i]; j++) {
      //  printf("#%3d : SA[%3d] = %3d, substring: ", j, Is[i] + j, SA[Is[i] + j]);
      //  for (int k=0; k<mlen; k++) {
      //    printf("%d ", s[SA[Is[i] + j] + k]);
      //  }
      //  printf("\n");
      //}
      //#endif
      N_curr += Ns[i];
    }
    *res = indices;
    free(Is); free(Ns);
  } else {
    // search as bipartite motif
    // search the longer one of the two parts
    int* Is;
    int* Ns;
    const char* q1; const char* q2;
    int q1len, q2len;
    int offset = 0;
    if (gloc >= mlen - (gloc + g)) {
      // search for prefix
      q1 = motif;
      q1len = gloc;
      q2 = &(motif[gloc + g]);
      q2len = mlen - (gloc + g);
    } else {
      // search for suffix
      offset = -(gloc + g);
      q1 = &(motif[gloc + g]);
      q1len = mlen - (gloc + g);
      q2 = motif;
      q2len = gloc;
      g = -g;
    }
    int n_enc = find_IUPAC(q1, q1len, SA, lcp, s, n, rmq,
                           &Is, &Ns);
    for (int i=0; i<n_enc; i++) {
      N_total += Ns[i];
    }
    int* indices = (int*)malloc(sizeof(int) * N_total);
    
    int new_N = 0;
    for (int i=0; i<n_enc; i++) {
      if (g >= 0) { // search q2 DOWNSTREAM of q1
        int max_loc = n - (q1len + g + q2len) - 1;
        for (int j=Is[i]; j<(Is[i]+Ns[i]); j++) {
          if (SA[j] <= max_loc) {
            //int o = arr_lcp(q2, &(s[SA[j] + q1len + g]), q2len);
            if (IUPAC_match(q2, &(s[SA[j] + q1len + g]), q2len)) {
              indices[new_N] = SA[j];
              new_N++;
            }       
          }
        }
      } else { // search q2 UPSTREAM of q1
        int max_loc = q2len - g;
        for (int j=Is[i]; j<(Is[i]+Ns[i]); j++) {
          if (SA[j] >= max_loc) {
            //int o = arr_lcp(q2, &(s[SA[j] + g - q2len]), q2len);
            if (IUPAC_match(q2, &(s[SA[j] + g - q2len]), q2len)) {
              indices[new_N] = SA[j] + offset;
              new_N++;
            }
          }         
        }
      }
    }
    indices = (int*)realloc(indices, sizeof(int) * new_N);
    *res = indices;
    N_total = new_N;
    free(Is); free(Ns);
  }

  return N_total;
}


int find_motif_nonparallel_indexed(
  const char* motif, int mlen,
  int* SA, int* lcp, int* s, int n, int*** index, int* SAr,
  int** res)
{
  int N_total = 0;
  // determine best motif-search strategy
  // find longest gap (stretch of N)
  int g, gloc, in_gap; // longest gap
  g = 0; gloc = 0; in_gap = 0;
  for (int i=1; i<mlen; i++) {
    if (motif[i] == 'N') {
      in_gap++;
    } else if (in_gap > 0) {
      if (in_gap >= g) {
        g = in_gap;
        gloc = i - g;
      }
      in_gap = 0;
    }
  }
  if (g == 0) {
    // search as contiguous motif
    int* Is;
    int* Ns;
    int n_enc = find_IUPAC_indexed(motif, mlen, SA, lcp, s, n, index,
                                   &Is, &Ns);
    for (int i=0; i<n_enc; i++) {
      N_total += Ns[i];
    }
    int* indices = (int*)malloc(sizeof(int) * N_total);
    int N_curr = 0;
    for (int i=0; i<n_enc; i++) {
      get_indices(Is[i], Ns[i], SA, &indices[N_curr]);
      //#ifdef DEBUG
      //printf("hits for encoding %d (n = %d):\n", i, Ns[i]);
      //for (int j=0; j<Ns[i]; j++) {
      //  printf("#%3d : SA[%3d] = %3d, substring: ", j, Is[i] + j, SA[Is[i] + j]);
      //  for (int k=0; k<mlen; k++) {
      //    printf("%d ", s[SA[Is[i] + j] + k]);
      //  }
      //  printf("\n");
      //}
      //#endif
      N_curr += Ns[i];
    }
    *res = indices;
    free(Is); free(Ns);
  } else {
    // search as bipartite motif
    // search the longer one of the two parts
    int* Is;
    int* Ns;
    const char* q1; const char* q2;
    int q1len, q2len;
    int offset = 0;
    if (gloc >= mlen - (gloc + g)) {
      // search for prefix
      q1 = motif;
      q1len = gloc;
      q2 = &(motif[gloc + g]);
      q2len = mlen - (gloc + g);
    } else {
      // search for suffix
      offset = -(gloc + g);
      q1 = &(motif[gloc + g]);
      q1len = mlen - (gloc + g);
      q2 = motif;
      q2len = gloc;
      g = -g;
    }
    int n_enc = find_IUPAC_indexed(q1, q1len, SA, lcp, s, n, index,
                                   &Is, &Ns);
    for (int i=0; i<n_enc; i++) {
      N_total += Ns[i];
    }
    int* indices = (int*)malloc(sizeof(int) * N_total);

    int factors[q2len];
    int n_enc2 = get_number_of_IUPAC_encodings(motif, q2len, factors);
    // define encodings array
    int** encodings = (int**)malloc(sizeof(int*) * n_enc2 + sizeof(int) * (n_enc2 * q2len));
    int* data = (int*)(encodings + n_enc2);
    for (int i=0; i<n_enc2; i++) {
      encodings[i] = data + (i * q2len);
    }
    // create encodings
    encode_IUPAC(q2, q2len, n_enc2, factors, encodings);
    // get SA ranges
    int idx;
    int sa_ranges[n_enc2][2];
    for (int i=0; i<n_enc2; i++) {
      idx = encoding_to_index(encodings[i], q2len, q2len);
      sa_ranges[i][0] = index[q2len-1][idx][0];
      sa_ranges[i][1] = index[q2len-1][idx][1];
    }
    
    
    int new_N = 0;
    if (g >= 0) { // search q2 DOWNSTREAM of q1
      int max_loc = n - (q1len + g + q2len) - 1;
      for (int i=0; i<n_enc; i++) {
        int loc;
        for (int j=Is[i]; j<(Is[i]+Ns[i]); j++) {
          if (SA[j] <= max_loc) {
            for (int enc=0; enc<n_enc2; enc++) {
              loc = SAr[SA[j] + q1len + g];
              if (sa_ranges[enc][0] <= loc && loc < sa_ranges[enc][1]) {
                indices[new_N] = SA[j];
                new_N++;
                break;
              }
            }
          }
        }
      }
    } else { // search q2 UPSTREAM of q1
      int max_loc = q2len - g;
      for (int i=0; i<n_enc; i++) {
        int loc;
        for (int j=Is[i]; j<(Is[i]+Ns[i]); j++) {
          if (SA[j] >= max_loc) {
            for (int enc=0; enc<n_enc2; enc++) {
              loc = SAr[SA[j] + g - q2len];
              if (sa_ranges[enc][0] <= loc && loc < sa_ranges[enc][1]) {
                indices[new_N] = SA[j] + offset;
                new_N++;
                break;
              }
            }
          }
        }
      }
    }
    indices = (int*)realloc(indices, sizeof(int) * new_N);
    *res = indices;
    N_total = new_N;
    free(Is); free(Ns);
  }

  return N_total;
}

//void motif_means(const char* motifs_data, int motif_count, int max_mlen,
//                 float* fwd, float* rev,
//                 int* SA, int* lcp, int* s, int n, int** rmq,
//                 float* mean_data, int* count_data)
void motif_means(
  const char* motifs_data, int motif_count, int max_mlen, const char* bases,
  //float* fwd, float* rev,
  //int* SA, int* lcp, int* s, int n, 
  //int** rmq, int*** index, int* SAr,
  float** fwd, float** rev,
  SuffixArray** SAs, int n_SAs,
  float* mean_data, int* count_data)
{
  // create row pointers for numpy arrays
  char** motifs = (char**)malloc(sizeof(char*) * motif_count);
  float** mean = (float**)malloc(sizeof(float*) * motif_count);
  int** count = (int**)malloc(sizeof(int*) * motif_count);
  for (int i=0; i<motif_count; i++) {
    motifs[i] = ((char*)motifs_data) + i * max_mlen;
    mean[i] = mean_data + i * (max_mlen - 1);
    count[i] = count_data + i * (max_mlen - 1);
  }
  // process each motif
  int n_bases = strlen(bases);
#pragma omp parallel for schedule(dynamic, 1000) if(motif_count > 10000)
  for (int m=0; m<motif_count; m++) {
    char* motif = motifs[m];
    int mlen = strlen(motifs[m]);
    char motif_rc[mlen + 1];
    reverse_complement(motif, mlen, motif_rc);
    // get offset locations
    int n_offsets = 0;
    int offsets[mlen];
    for (int o=0; o<mlen; o++) {
      //if (motif[o] == 'A' || motif[o] == 'C') {// || motif[o] == 'M') {
      //  offsets[n_offsets] = o;
      //  n_offsets++;
      //}
      for (int b=0; b<n_bases; b++) {
        if (motif[o] == bases[b]) {
          offsets[n_offsets] = o;
          n_offsets++;
          break;
        }
      }
    }
    // 
    int *indices, *indices_rc;
    double total_sum[n_offsets];
    int total_n[n_offsets];
    for (int o=0; o<n_offsets; o++) {
      total_sum[o] = 0.0;
      total_n[o] = 0;
    }
    int n_indices, idx;
    for (int k=0; k<n_SAs; k++) {
      if (mlen > INDEX_SIZE) {
        n_indices = find_motif_nonparallel((const char*)motif, mlen,
                              SAs[k]->sa, SAs[k]->lcp, SAs[k]->s, SAs[k]->n, SAs[k]->rmq,
                              &indices);
      } else {
        n_indices = find_motif_nonparallel_indexed((const char*)motif, mlen,
                              SAs[k]->sa, SAs[k]->lcp, SAs[k]->s, SAs[k]->n, SAs[k]->index, SAs[k]->sar,
                              &indices);
      }
      for (int i=0; i<n_indices; i++) {
        for (int o=0; o<n_offsets; o++) {
          idx = indices[i] + offsets[o];
          if (fwd[k][idx] == fwd[k][idx]) { // non-nan
            total_sum[o] += fwd[k][idx];
            total_n[o]++;
          }
        }
      }
      if (mlen > INDEX_SIZE) {
        n_indices = find_motif_nonparallel((const char*)motif_rc, mlen,
                              SAs[k]->sa, SAs[k]->lcp, SAs[k]->s, SAs[k]->n, SAs[k]->rmq,
                              &indices_rc);
      } else {
        n_indices = find_motif_nonparallel_indexed((const char*)motif_rc, mlen,
                              SAs[k]->sa, SAs[k]->lcp, SAs[k]->s, SAs[k]->n, SAs[k]->index, SAs[k]->sar,
                              &indices_rc);
      }
      for (int i=0; i<n_indices; i++) {
        for (int o=0; o<n_offsets; o++) {
          idx = indices_rc[i] + mlen - 1 - offsets[o];
          if (rev[k][idx] == rev[k][idx]) { // non-nan
            total_sum[o] += rev[k][idx];
            total_n[o]++;
          }
        }
      }
    }
    for (int o=0; o<n_offsets; o++) {
      mean[m][offsets[o]] = (float)(total_sum[o] / total_n[o]);
      count[m][offsets[o]] = total_n[o];
    }
    free(indices); free(indices_rc);
  }
  // free memory
  free(motifs); free(mean); free(count);
}

// adapted from:
//Algorithm from Numerical recipes in C of 1992
//https://stackoverflow.com/questions/1961173/median-function-in-c-math-library
float quick_select_median(float arr[], uint32_t n)
{
  uint32_t low, high ;
  uint32_t median;
  uint32_t middle, ll, hh;
  low = 0 ; high = n-1 ; median = (low + high) / 2;
  for (;;) {
    if (high <= low) {/* One element only */
      if (n % 2) {
        return arr[median] ;
      } else {
        return (arr[median] + arr[median + 1]) / 2.0 ;
      }
    }
    if (high == low + 1) { /* Two elements only */
      if (arr[low] > arr[high])
      ELEM_SWAP(arr[low], arr[high]) ;
      if (n % 2) {
        return arr[median] ;
      } else {
        return (arr[median] + arr[median + 1]) / 2.0 ;
      }
    }
    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])
      ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])
      ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])
      ELEM_SWAP(arr[middle], arr[low]) ;
    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do ll++; while (arr[low] > arr[ll]) ;
      do hh--; while (arr[hh] > arr[low]) ;
      if (hh < ll)
        break;
      ELEM_SWAP(arr[ll], arr[hh]) ;
    }
    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;
    /* Re-set active partition */
    if (hh <= median)
      low = ll;
    if (hh >= median)
      high = hh - 1;
  }
  if (n % 2) {
    return arr[median] ;
  } else {
    return (arr[median] + arr[median + 1]) / 2.0 ;
  }
}

void motif_medians(
  const char* motifs_data, int motif_count, int max_mlen, const char* bases,
  //float* fwd, float* rev,
  //int* SA, int* lcp, int* s, int n, 
  //int** rmq, int*** index, int* SAr,
  float** fwd, float** rev,
  SuffixArray** SAs, int n_SAs,
  float* median_data, int* count_data)
{
  // create row pointers for numpy arrays
  char** motifs = (char**)malloc(sizeof(char*) * motif_count);
  float** median = (float**)malloc(sizeof(float*) * motif_count);
  int** count = (int**)malloc(sizeof(int*) * motif_count);
  for (int i=0; i<motif_count; i++) {
    motifs[i] = ((char*)motifs_data) + i * max_mlen;
    median[i] = median_data + i * (max_mlen - 1);
    count[i] = count_data + i * (max_mlen - 1);
  }
  // process each motif
  int n_bases = strlen(bases);
#pragma omp parallel for schedule(dynamic, 1000) if(motif_count > 10000)
  for (int m=0; m<motif_count; m++) {
    char* motif = motifs[m];
    int mlen = strlen(motifs[m]);
    char motif_rc[mlen + 1];
    reverse_complement(motif, mlen, motif_rc);
    //// get offset locations
    //int n_offsets = 0;
    //int offsets[mlen+4];
    //for (int o=-2; o<mlen+2; o++) {
    //  //if (motif[o] == 'A' || motif[o] == 'C' || motif[o] == 'M') {
    //  offsets[n_offsets] = o;
    //  n_offsets++;
    //  //}
    //}
    int n_offsets = 0;
    int offsets[mlen];
    for (int o=0; o<mlen; o++) {
      //if (motif[o] == 'A' || motif[o] == 'C') {// || motif[o] == 'M') {
      //  offsets[n_offsets] = o;
      //  n_offsets++;
      //}
      for (int b=0; b<n_bases; b++) {
        if (motif[o] == bases[b]) {
          offsets[n_offsets] = o;
          n_offsets++;
          break;
        }
      }
    }

    // search all contigs for all occurences
    int n_indices[n_SAs], n_indices_rc[n_SAs], idx;
    int *indices[n_SAs], *indices_rc[n_SAs];
    int n_indices_total = 0;
    //double total_sum[n_offsets];
    int total_n[n_offsets];
    for (int o=0; o<n_offsets; o++) {
      //total_sum[o] = 0.0;
      total_n[o] = 0;
    }
    for (int k=0; k<n_SAs; k++) {
      //printf("k: %d", k);
      //n_indices = find_motif_nonparallel((const char*)motif, mlen,
      //                       SA, lcp, s, n, rmq,
      //                       &indices);
      if (mlen > INDEX_SIZE) {
        n_indices[k] = find_motif_nonparallel((const char*)motif, mlen,
                              SAs[k]->sa, SAs[k]->lcp, SAs[k]->s, SAs[k]->n, SAs[k]->rmq,
                              &(indices[k]));
      } else {
        n_indices[k] = find_motif_nonparallel_indexed((const char*)motif, mlen,
                              SAs[k]->sa, SAs[k]->lcp, SAs[k]->s, SAs[k]->n, SAs[k]->index, SAs[k]->sar,
                              &(indices[k]));
      }
      //n_indices_rc = find_motif_nonparallel((const char*)motif_rc, mlen,
      //                       SA, lcp, s, n, rmq,
      //                       &indices_rc);
      if (mlen > INDEX_SIZE) {
        n_indices_rc[k] = find_motif_nonparallel((const char*)motif_rc, mlen,
                              SAs[k]->sa, SAs[k]->lcp, SAs[k]->s, SAs[k]->n, SAs[k]->rmq,
                              &(indices_rc[k]));
      } else {
        n_indices_rc[k] = find_motif_nonparallel_indexed((const char*)motif_rc, mlen,
                              SAs[k]->sa, SAs[k]->lcp, SAs[k]->s, SAs[k]->n, SAs[k]->index, SAs[k]->sar,
                              &(indices_rc[k]));
      }
      n_indices_total += n_indices[k] + n_indices_rc[k];
      //printf(" done \n");
    }
    //printf("n_indices_total: %d\n", n_indices_total);

    // define encodings array
    float** position_vals = (float**)malloc(sizeof(float*) * n_offsets + sizeof(float) * (n_offsets * n_indices_total));
    float* data = (float*)(position_vals + n_offsets);
    for (int i=0; i<n_offsets; i++) {
      position_vals[i] = data + (i * n_indices_total);
    }
    // gather non-nan values
    for (int k=0; k<n_SAs; k++) {
      for (int i=0; i<n_indices[k]; i++) {
        for (int o=0; o<n_offsets; o++) {
          idx = indices[k][i] + offsets[o];
          if (idx < 0 || idx >= (SAs[k]->n-1))
            continue;
          if (fwd[k][idx] == fwd[k][idx]) { // non-nan
            position_vals[o][total_n[o]] = fwd[k][idx]; //fabsf(fwd[idx]);
            total_n[o]++;
          }
        }
      }
      for (int i=0; i<n_indices_rc[k]; i++) {
        for (int o=0; o<n_offsets; o++) {
          idx = indices_rc[k][i] + mlen - 1 - offsets[o];
          if (idx < 0 || idx >= (SAs[k]->n-1))
            continue;
          if (rev[k][idx] == rev[k][idx]) { // non-nan
            position_vals[o][total_n[o]] = rev[k][idx]; //fabsf(rev[idx]);
            total_n[o]++;
          }
        }
      }
    }
    // calculate medians
    for (int o=0; o<n_offsets; o++) {
      if (total_n[o] > 0) {
        median[m][offsets[o]] = quick_select_median(position_vals[o], total_n[o]);
      } else {
        median[m][offsets[o]] = NAN;
      }
      count[m][offsets[o]] = total_n[o];
    }
    for (int k=0; k<n_SAs; k++) {
      free(indices[k]); free(indices_rc[k]);
    }
    free(position_vals);
  }
  // free memory
  free(motifs); free(median); free(count);
}

void all_motif_medians(
  const char* motifs_data, int motif_count, int max_mlen, int pad,
  float* fwd, float* rev,
  int* SA, int* lcp, int* s, int n, 
  int** rmq, int*** index, int* SAr,
  float* median_data, int* count_data)
{
  // create row pointers for numpy arrays
  char** motifs = (char**)malloc(sizeof(char*) * motif_count);
  float** median = (float**)malloc(sizeof(float*) * motif_count);
  int** count = (int**)malloc(sizeof(int*) * motif_count);
  for (int i=0; i<motif_count; i++) {
    motifs[i] = ((char*)motifs_data) + i * max_mlen;
    median[i] = median_data + i * (max_mlen - 1 + 2*pad);
    count[i] = count_data + i * (max_mlen - 1 + 2*pad);
  }
  // process each motif
#pragma omp parallel for schedule(dynamic, 1000) if(motif_count > 10000)
  for (int m=0; m<motif_count; m++) {
    char* motif = motifs[m];
    int mlen = strlen(motifs[m]);
    char motif_rc[mlen + 1];
    reverse_complement(motif, mlen, motif_rc);
    // get offset locations
    int n_offsets = 0;
    int offsets[mlen+2*pad];
    for (int o=-pad; o<mlen+pad; o++) {
      //if (motif[o] == 'A' || motif[o] == 'C' || motif[o] == 'M') {
      offsets[n_offsets] = o;
      n_offsets++;
      //}
    }
    // 
    int n_indices, n_indices_rc, idx;
    int *indices, *indices_rc;
    //double total_sum[n_offsets];
    int total_n[n_offsets];
    for (int o=0; o<n_offsets; o++) {
      //total_sum[o] = 0.0;
      total_n[o] = 0;
    }
    //n_indices = find_motif_nonparallel((const char*)motif, mlen,
    //                       SA, lcp, s, n, rmq,
    //                       &indices);
    if (mlen > INDEX_SIZE) {
      n_indices = find_motif_nonparallel((const char*)motif, mlen,
                             SA, lcp, s, n, rmq,
                             &indices);
    } else {
      n_indices = find_motif_nonparallel_indexed((const char*)motif, mlen,
                            SA, lcp, s, n, index, SAr,
                            &indices);
    }
    //n_indices_rc = find_motif_nonparallel((const char*)motif_rc, mlen,
    //                       SA, lcp, s, n, rmq,
    //                       &indices_rc);
    if (mlen > INDEX_SIZE) {
      n_indices_rc = find_motif_nonparallel((const char*)motif_rc, mlen,
                             SA, lcp, s, n, rmq,
                             &indices_rc);
    } else {
      n_indices_rc = find_motif_nonparallel_indexed((const char*)motif_rc, mlen,
                            SA, lcp, s, n, index, SAr,
                            &indices_rc);
    }
    int n_indices_total = n_indices + n_indices_rc;
    // define encodings array
    float** position_vals = (float**)malloc(sizeof(float*) * n_offsets + sizeof(float) * (n_offsets * n_indices_total));
    float* data = (float*)(position_vals + n_offsets);
    for (int i=0; i<n_offsets; i++) {
      position_vals[i] = data + (i * n_indices_total);
    }
    // gather non-nan values
    for (int i=0; i<n_indices; i++) {
      for (int o=0; o<n_offsets; o++) {
        idx = indices[i] + offsets[o];
        if (idx < 0 || idx >= (n-1))
          continue;
        if (fwd[idx] == fwd[idx]) { // non-nan
          position_vals[o][total_n[o]] = fwd[idx]; //fabsf(fwd[idx]);
          total_n[o]++;
        }
      }
    }
    for (int i=0; i<n_indices_rc; i++) {
      for (int o=0; o<n_offsets; o++) {
        idx = indices_rc[i] + mlen - 1 - offsets[o];
        if (idx < 0 || idx >= (n-1))
          continue;
        if (rev[idx] == rev[idx]) { // non-nan
          position_vals[o][total_n[o]] = rev[idx]; //fabsf(rev[idx]);
          total_n[o]++;
        }
      }
    }
    // calculate medians
    for (int o=0; o<n_offsets; o++) {
      if (total_n[o] > 0) {
        median[m][o] = quick_select_median(position_vals[o], total_n[o]);
      } else {
        median[m][o] = NAN;
      }
      count[m][o] = total_n[o];
    }
    free(indices); free(indices_rc); free(position_vals);
  }
  // free memory
  free(motifs); free(median); free(count);
}

void write_sa(SuffixArray* sa, char* fn)
{
  FILE *fh = fopen(fn, "wb");
  int h = floor(log((double)sa->n) / log(2.)) + 1;
  int* rmq_data = (int*)((sa->rmq) + h);
  int* index_data = (int*)((sa->index) + INDEX_SIZE + INDEX_SIZE_N);
  fwrite(sa->sa, sizeof(int), sa->n, fh);
  fwrite(sa->s, sizeof(int), sa->n + 3, fh);
  fwrite(sa->lcp, sizeof(int), sa->n, fh);
  fwrite(rmq_data, sizeof(int), sa->n * h, fh);
  fwrite(sa->sar, sizeof(int), sa->n, fh);
  fwrite(index_data, sizeof(int), INDEX_SIZE_N * 2, fh);
  fwrite(&(sa->n), sizeof(int), 1, fh);
  fwrite(&(sa->K), sizeof(int), 1, fh);
  fclose(fh);
}

void read_sa(char* fn, int n, SuffixArray* sa)
{
  int h = floor(log((double)n) / log(2.)) + 1;
  sa->sa = (int*)malloc(sizeof(int) * n);
  sa->s = (int*)malloc(sizeof(int) * (n + 3));
  sa->lcp = (int*)malloc(sizeof(int) * n);
  sa->rmq = (int**)malloc(sizeof(int*) * (h+1) + sizeof(int) * n * h);
  sa->sar = (int*)malloc(sizeof(int) * n);
  //sa->index = (int**)malloc(sizeof(int*) * INDEX_SIZE_N + sizeof(int) * 2 * INDEX_SIZE_N);//(int*)malloc(sizeof(int) * n);
  sa->index = (int***)malloc(sizeof(int**) * INDEX_SIZE + sizeof(int*) * INDEX_SIZE_N + sizeof(int) * 2 * INDEX_SIZE_N);

  int* rmq_data = (int*)((sa->rmq) + h);
  int* index_data = (int*)((sa->index) + INDEX_SIZE + INDEX_SIZE_N);
  FILE *fh = fopen(fn, "rb");
  fread(sa->sa, sizeof(int), n, fh);
  fread(sa->s, sizeof(int), n + 3, fh);
  fread(sa->lcp, sizeof(int), n, fh);
  fread(rmq_data, sizeof(int), n * h, fh);
  fread(sa->sar, sizeof(int), n, fh);
  fread(index_data, sizeof(int), INDEX_SIZE_N * 2, fh);
  fread(&(sa->n), sizeof(int), 1, fh);
  fread(&(sa->K), sizeof(int), 1, fh);
  fclose(fh);
  for (int j=0; j<h; j++) {
    (sa->rmq)[j] = rmq_data + (j * n);
  }
  for (int i=0; i<INDEX_SIZE; i++) {
    (sa->index)[i] = (int**)(sa->index) + INDEX_SIZE;
    int offset = 0;
    for (int j=0; j<i; j++) {
      offset += (int)pow(4, j+1);
    }
    (sa->index)[i] += offset;
    for (int j=0; j<(int)pow(4, i+1); j++) {
      (sa->index)[i][j] = index_data + 2 * offset + (j * 2);
    }
  }
}