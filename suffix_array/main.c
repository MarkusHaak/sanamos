#include "suffix_array.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char** argv) {
  char* str;
  char* query;
  char* query2;
  int g = 0;
  if (argc < 3) {
    return EXIT_FAILURE;
  } else {
    str = argv[1];
    query = argv[2];
    if (argc > 3) {
      g = atoi(argv[3]);
      query2 = argv[4];
    }
  }
  
  //printf("lets go\n");

  // initialize nucleic acid encoding array
  na_enc_init();
  
  // initialize suffix array
  //SuffixArray sa;
  //sa = suffix_array_init(str, NA_ENC);
  //printf("create sa\n");
  SuffixArray sa0;
  //printf("init sa\n");
  suffix_array_init(str, NA_ENC, &sa0);

  //printf("%s\n", str);
  //for (int i=0; i < strlen(str); i++) {
  //  printf("%d ", sa0.s[i]);
  //}
  //printf("\n");

  // construct suffix array
  //printf("run skew\n");
  skew(sa0.s, sa0.sa, sa0.n, sa0.K);

  // construct longest common prefix array
  //int* lcp = (int*)malloc(sizeof(int) * sa0.n);
  kasai(sa0.s, sa0.sa, sa0.n, sa0.lcp, sa0.sar);
  
  //// display sa & lca in tabular form
  //printf("%4s %4s %4s %s\n", "i", "sa", "lcd", "suffix");
  //for (int i=0; i<(15+sa0.n); i++) {printf("-");} printf("\n");
  //printf("%4d %4d %4d %s -\n", 0, sa0.sa[0], lcp[0], &sa0.str[sa0.sa[0]]);
  //for (int i=1; i < sa0.n; i++) {
  //  printf("%4d %4d %4d %s\n", i, sa0.sa[i], lcp[i], &sa0.str[sa0.sa[i]]);
  //}

  // create an RMQ array
  //int** rmq;
  construct_rmq_array(sa0.lcp, sa0.n, sa0.rmq);
  //for (int l=0; l<sa0.n; l++) {
  //  for (int r=l; r<sa0.n; r++) {
  //    int h = floor(log2((double)r - l + 1));
  //    int a = rmq[h][l];
  //    int b = rmq[h][r - (1 << h) + 1];
  //    int v = a < b ? a : b;
  //    printf("(%d, %d) : %d\n", l, r, v);
  //  }
  //}

  print_sa_lcd(&sa0, str);
  create_index(sa0.s, sa0.sa, sa0.n, sa0.lcp, sa0.index);

  // create and encode a search query
  int qlen = strlen(query);
  int q[qlen];
  encode_str(query, NA_ENC, q);

  // save and load sa
  int n = sa0.n;
  write_sa(&sa0, "test.index");
  SuffixArray sa;
  read_sa("test.index", n, &sa);


  
  //printf("%s (%d characters): ", query, qlen);
  //for (int i=0; i < qlen; i++) {
  //  printf(" %d", q[i]);
  //}
  //printf("\n");
  
  //// search the query
  //int n_found = 0;
  //int res = 0;
  //n_found = find(q, qlen, sa.sa, sa.lcp, sa.s, sa.n, sa.rmq,
  //     &res);
  //if (n_found == 0) {
  //  printf("not found\n");
  //} else {
  //  printf("found, idx: %d\n", res);
  //}
//
  //n_found = find_all(q, qlen, sa.sa, sa.lcp, sa.s, sa.n, sa.rmq,
  //                   &res);
  //int idx = res;
  //if (n_found == 0) {
  //  printf("not found\n");
  //} else {
  //  printf("found, idx: %d, n = %d\n", idx, n_found);
  //}
//
  //if (g != 0) {
  //  int indices[n_found];
  //  int q2len = strlen(query2);
  //  int q2[q2len];
  //  encode_str(query2, NA_ENC, q2);
  //  n_found = get_indices_bipartite(idx, n_found, sa.sa, sa.s, sa.n, qlen, g, q2, q2len, indices);
  //  printf("found %d bipartite motifs\n", n_found);
//
  //  //int** ret = (int**)malloc(sizeof(int*));
  //  int* ret = NULL;
  //  n_found = find_all_bipartite(q, qlen, g, q2, q2len, 
  //                    sa.sa, sa.lcp, sa.s, sa.n, sa.rmq,
  //                    &ret);
  //  printf("found %d bipartite motifs\n", n_found);
  //  free(ret);
  //}
  ////printf("huray\n");

  //char IUPAC_motif[5] = {'A','T','C','G','U'}; 
  char* IUPAC_motif = argv[5];
  int* indices3;
  int n_found = 0;
  n_found = find_motif(IUPAC_motif, strlen(IUPAC_motif), sa.sa, sa.lcp, sa.s, sa.n, sa.rmq,
             &indices3);
  printf("IUPAC %s found: %d\n", IUPAC_motif, n_found);
  n_found = find_motif_nonparallel_indexed(IUPAC_motif, strlen(IUPAC_motif), sa.sa, sa.lcp, sa.s, sa.n, sa.index, sa.sar,
             &indices3);
  printf("IUPAC %s found: %d\n", IUPAC_motif, n_found);

  //int encoding[8] = {4,4,4,4,4,4,4,4};
  //printf("%#04X\n", encoding_to_index(encoding, 8));
  int idx = encoding_to_index(q, qlen, qlen);
  printf("%d %d\n", sa0.index[qlen-1][idx][0], sa0.index[qlen-1][idx][1]);
  printf("%#04X\n", idx);
  for (int j=sa.index[qlen-1][idx][0]; j<sa.index[qlen-1][idx][1]; j++) {
    printf("SA[%d] = %d : ", j, sa.sa[j]);
    for (int c=0; c<qlen; c++) {
      printf("%c", str[sa.sa[j]+c]);
    }
    printf("\n");
  }

  return EXIT_SUCCESS;
}