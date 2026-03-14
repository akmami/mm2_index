#include "index.h"

#include "ksort.h"
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8)

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)

KSORT_INIT_GENERIC(uint32_t)

/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int n, n_keys;
	size_t j, start_a, start_p;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t *)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0)
		return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j)
	{
		if (j == b->a.n || b->a.a[j].x >> 8 != b->a.a[j - 1].x >> 8)
		{
			++n_keys;
			if (n > 1)
				b->n += n;
			n = 1;
		}
		else
			++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t *)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j)
	{
		if (j == b->a.n || b->a.a[j].x >> 8 != b->a.a[j - 1].x >> 8)
		{
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j - 1];
			itr = kh_put(idx, h, p->x >> 8 >> mi->b << 1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1)
			{
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			}
			else
			{
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p << 32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		}
		else
			++n;
	}
	b->h = h;
	assert(b->n == (int32_t)start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}

static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1 << mi->b);
}

/******************
 * Generate index *
 ******************/

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1 << mi->b) - 1;
	for (i = 0; i < n; ++i)
	{
		mm128_v *p = &mi->B[a[i].x >> 8 & mask].a;
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
	pipeline_t *p = (pipeline_t *)shared;
	if (step == 0)
	{ // step 0: read sequences
		step_t *s;
		if (p->sum_len > p->batch_size)
			return 0;
		s = (step_t *)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq)
		{
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m);
			kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t *)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ))
			{
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i)
					sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len);
				kroundup64(max_len);
				if (old_max_len != max_len)
				{
					p->mi->S = (uint32_t *)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i)
			{
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME))
				{
					seq->name = (char *)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				}
				else
					seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				seq->is_alt = 0;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ))
				{
					for (j = 0; j < seq->len; ++j)
					{ // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		}
		else
			free(s);
	}
	else if (step == 1)
	{ // step 1: compute sketch
		step_t *s = (step_t *)in;
		for (i = 0; i < s->n_seq; ++i)
		{
			mm_bseq1_t *t = &s->seq[i];
			if (t->l_seq > 0)
				mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, p->mi->flag & MM_I_HPC, &s->a);
			else
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			free(t->seq);
			free(t->name);
		}
		free(s->seq);
		s->seq = 0;
		return s;
	}
	else if (step == 2)
	{ // dispatch sketch to buckets
		step_t *s = (step_t *)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		kfree(0, s->a.a);
		free(s);
	}
	return 0;
}

// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0)
		return 0; // failed to open the index
	r = (mm_idx_reader_t *)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt)
		r->opt = *opt;
	else
		mm_idxopt_init(&r->opt);
	if (r->is_idx)
	{
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	}
	else
		r->fp.seq = mm_bseq_open(fn);
	if (fn_out)
		r->fp_out = fopen(fn_out, "wb");
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r->is_idx)
		fclose(r->fp.idx);
	else
		mm_bseq_close(r->fp.seq);
	if (r->fp_out)
		fclose(r->fp_out);
	free(r);
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
	mm_idx_t *mi;
	if (r->is_idx)
	{
		mi = mm_idx_load(r->fp.idx);
		if (mi && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag & MM_I_HPC) != (r->opt.flag & MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	}
	else
		mi = mm_idx_gen_fa(r->fp.seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
	if (mi)
	{
		if (r->fp_out)
			mm_idx_dump(r->fp_out, mi);
		mi->index = r->n_parts++;
	}
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	uint32_t i;
	if (mi == 0)
		return;
	if (mi->h)
		kh_destroy(str, (khash_t(str) *)mi->h);
	if (mi->B)
	{
		for (i = 0; i < 1U << mi->b; ++i)
		{
			free(mi->B[i].p);
			free(mi->B[i].a.a);
			kh_destroy(idx, (idxhash_t *)mi->B[i].h);
		}
	}
	if (mi->spsc)
		free(mi->spsc);
	if (mi->I)
	{
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->I[i].a);
		free(mi->I);
	}
	if (!mi->km)
	{
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	}
	else
		km_destroy(mi->km);
	free(mi->B);
	free(mi->S);
	free(mi);
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.)
		return INT32_MAX;
	for (i = 0; i < 1 << mi->b; ++i)
		if (mi->B[i].h)
			n += kh_size((idxhash_t *)mi->B[i].h);
	if (n == 0)
		return INT32_MAX;
	a = (uint32_t *)malloc(n * 4);
	for (i = n = 0; i < 1 << mi->b; ++i)
	{
		idxhash_t *h = (idxhash_t *)mi->B[i].h;
		if (h == 0)
			continue;
		for (k = 0; k < kh_end(h); ++k)
		{
			if (!kh_exist(h, k))
				continue;
			a[n++] = kh_key(h, k) & 1 ? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)
{
	if ((opt->flag & MM_F_SPLICE_FOR) || (opt->flag & MM_F_SPLICE_REV))
		opt->flag |= MM_F_SPLICE;
	if (opt->mid_occ <= 0)
	{
		opt->mid_occ = mm_idx_cal_max_occ(mi, opt->mid_occ_frac);
		if (opt->mid_occ < opt->min_mid_occ)
			opt->mid_occ = opt->min_mid_occ;
		if (opt->max_mid_occ > opt->min_mid_occ && opt->mid_occ > opt->max_mid_occ)
			opt->mid_occ = opt->max_mid_occ;
	}
	if (opt->bw_long < opt->bw)
		opt->bw_long = opt->bw;

	fprintf(stderr, "[M::%s] mid_occ = %d\n", __func__, opt->mid_occ);
}

void mm_idx_stat(const mm_idx_t *mi)
{
	int64_t n = 0, n1 = 0;
	uint32_t i;
	uint64_t sum = 0, len = 0;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag & MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1U << mi->b; ++i)
		if (mi->B[i].h)
			n += kh_size((idxhash_t *)mi->B[i].h);
	for (i = 0; i < 1U << mi->b; ++i)
	{
		idxhash_t *h = (idxhash_t *)mi->B[i].h;
		khint_t k;
		if (h == 0)
			continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k))
			{
				sum += kh_key(h, k) & 1 ? 1 : (uint32_t)kh_val(h, k);
				if (kh_key(h, k) & 1)
					++n1;
			}
	}
	fprintf(stderr, "[M::%s] distinct minimizers: %ld (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf; total length: %ld\n",
			__func__, (long)n, 100.0 * n1 / n, (double)sum / n, (double)len / sum, (long)len);
}

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0)
		return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0)
		return -1; // error

	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4)
	{
		lseek(fd, 0, SEEK_SET);
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx ? off_end : 0;
}

void mm_idxopt_init(mm_idxopt_t *opt)
{
	memset(opt, 0, sizeof(mm_idxopt_t));
	opt->k = 15, opt->w = 10, opt->flag = 0;
	opt->bucket_bits = 14;
	opt->mini_batch_size = 50000000;
	opt->batch_size = 8000000000ULL;
}

mm_idx_t *mm_idx_gen_fa(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp))
		return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size ? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);

	kt_pipeline(n_threads < 3 ? n_threads : 3, worker_pipeline, &pl, 3);
	fprintf(stderr, "[M::%s] collected minimizers\n", __func__);

	mm_idx_post(pl.mi, n_threads);
	fprintf(stderr, "[M::%s] sorted minimizers\n", __func__);

	return pl.mi;
}

// for seeds
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}


mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	uint32_t x[5], i;
	uint64_t sum_len = 0;
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4)
		return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0)
		return 0;
	if (fread(x, 4, 5, fp) != 5)
		return 0;
	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t *)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i)
	{
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
		size_t _ = fread(&l, 1, 1, fp);
		(void)_;
		if (l)
		{
			s->name = (char *)kmalloc(mi->km, l + 1);
			_ = fread(s->name, 1, l, fp);
			(void)_;
			s->name[l] = 0;
		}
		_ = fread(&s->len, 4, 1, fp);
		(void)_;
		s->offset = sum_len;
		s->is_alt = 0;
		sum_len += s->len;
	}
	for (i = 0; i < 1 << mi->b; ++i)
	{
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		size_t _ = fread(&b->n, 4, 1, fp);
		(void)_;
		b->p = (uint64_t *)malloc(b->n * 8);
		_ = fread(b->p, 8, b->n, fp);
		(void)_;
		_ = fread(&size, 4, 1, fp);
		(void)_;
		if (size == 0)
			continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j)
		{
			uint64_t x[2];
			int absent;
			_ = fread(x, 8, 2, fp);
			(void)_;
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
	{
		mi->S = (uint32_t *)malloc((sum_len + 7) / 8 * 4);
		size_t _ = fread(mi->S, 4, (sum_len + 7) / 8, fp);
		(void)_;
	}
	return mi;
}

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[5], i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 5, fp);
	for (i = 0; i < mi->n_seq; ++i)
	{
		if (mi->seq[i].name)
		{
			uint8_t l = strlen(mi->seq[i].name);
			fwrite(&l, 1, 1, fp);
			fwrite(mi->seq[i].name, 1, l, fp);
		}
		else
		{
			uint8_t l = 0;
			fwrite(&l, 1, 1, fp);
		}
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1 << mi->b; ++i)
	{
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t *)b->h;
		uint32_t size = h ? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0)
			continue;
		for (k = 0; k < kh_end(h); ++k)
		{
			uint64_t x[2];
			if (!kh_exist(h, k))
				continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k * 2 < b)
		b = k * 2;
	if (w < 1)
		w = 1;
	mi = (mm_idx_t *)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t *)calloc(1 << b, sizeof(mm_idx_bucket_t));
	// if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

// -----------------------------------------------
// -----------------------------------------------
// -----------------------------------------------
// -----------------------------------------------

mm_idx_t *mm_idx_gfa(graph_t *g, mm_idx_reader_t *r)
{
	mm_idx_t *mi;
	mi = mm_idx_gen_gfa(g, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, r->opt.batch_size);
	return mi;
}

mm_idx_reader_t *mm_idx_gfa_init(const mm_idxopt_t *opt)
{
	mm_idx_reader_t *r;
	r = (mm_idx_reader_t *)calloc(1, sizeof(mm_idx_reader_t));
	if (opt)
		r->opt = *opt;
	else
		mm_idxopt_init(&r->opt);
	return r;
}

mm_idx_t *mm_idx_gen_gfa(graph_t *g, int w, int k, int b, int flag, int mini_batch_size, uint64_t batch_size)
{
	mm_idx_t *mi = mm_idx_init(w, k, b, flag);
	mm128_v a;
	int i;
	a.n = a.m = 0;
	a.a = NULL;

	// step 1: compute sketch
	for (i = 0; i <  g->n_nodes; ++i)
	{
		mm_sketch(0, g->char_table.data + g->nodes[i].char_offset, g->nodes[i].seq_len, w, k, i, flag & MM_I_HPC, &a);
	}

	fprintf(stderr, "[M::%s] collected minimizers\n", __func__);
	
	// dispatch sketch to buckets
	mm_idx_add(mi, a.n, a.a);
	kfree(0, a.a);
	
	// post processing
	mm_idx_post(mi, 1);
	fprintf(stderr, "[M::%s] sorted minimizers\n", __func__);

	return mi;
}