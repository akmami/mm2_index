#ifndef INDEX_H
#define INDEX_H

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#define __STDC_LIMIT_MACROS
#include "kthread.h"
#include "bseq.h"
#include "minimap.h"
#include "sketch.h"
#include "kvec.h"
#include "khash.h"
#include "gfa.h"

#define idx_hash(a) ((a) >> 1)
#define idx_eq(a, b) ((a) >> 1 == (b) >> 1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x) |= (x) >> 1, (x) |= (x) >> 2, (x) |= (x) >> 4, (x) |= (x) >> 8, (x) |= (x) >> 16, (x) |= (x) >> 32, ++(x))

#define mm_seq4_set(s, i, c) ((s)[(i) >> 3] |= (uint32_t)(c) << (((i) & 7) << 2))
#define mm_seq4_get(s, i) ((s)[(i) >> 3] >> (((i) & 7) << 2) & 0xf)

typedef struct mm_idx_bucket_s
{
	mm128_v a;	 // (minimizer, position) array
	int32_t n;	 // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;	 // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct
{
	int32_t st, en, cnt;
	int32_t score : 30, strand : 2;
} mm_idx_intv1_t;

typedef struct mm_idx_intv_s
{
	int32_t n, m;
	mm_idx_intv1_t *a;
} mm_idx_intv_t;

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out);

void mm_idx_reader_close(mm_idx_reader_t *r);

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads);

void mm_idx_destroy(mm_idx_t *mi);

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f);

void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi);

void mm_idx_stat(const mm_idx_t *mi);

int64_t mm_idx_is_idx(const char *fn);

void mm_idxopt_init(mm_idxopt_t *opt);

mm_idx_t *mm_idx_gen_fa(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size);

mm_idx_reader_t *mm_idx_gfa_init(const mm_idxopt_t *opt);

mm_idx_t *mm_idx_gfa(graph_t *g, mm_idx_reader_t *r);

mm_idx_t *mm_idx_gen_gfa(graph_t *g, int w, int k, int b, int flag, int mini_batch_size, uint64_t batch_size);

/******************
 * Generate index *
 ******************/

typedef struct
{
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
} pipeline_t;

typedef struct
{
	int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

mm_idx_t *mm_idx_load(FILE *fp);

void mm_idx_dump(FILE *fp, const mm_idx_t *mi);

mm_idx_t *mm_idx_init(int w, int k, int b, int flag);

#endif