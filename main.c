#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "index.h"
#include "gfa.h"

int main(int argc, char *argv[])
{

    if (argc < 5) {
        printf("usage: ./segram <k> <w> <fasta> <gfa>\n");
        return -1;
    }

	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	int i, n_threads = 3;
	char *fnw = 0;
	mm_idx_reader_t *idx_rdr;
	mm_idx_t *mi;

    ipt.k = atoi(argv[1]);
    ipt.w = atoi(argv[2]);
    char *fasta = argv[3];
    char *gfa = argv[4];
	
    ipt.flag = 0, ipt.k = 21, ipt.w = 11;
    opt.flag |= MM_F_SR | MM_F_FRAG_MODE | MM_F_NO_PRINT_2ND | MM_F_2_IO_THREADS | MM_F_HEAP_SORT;
    opt.pe_ori = 0<<1|1; // FR
    opt.a = 2, opt.b = 8, opt.q = 12, opt.e = 2, opt.q2 = 24, opt.e2 = 1;
    opt.zdrop = opt.zdrop_inv = 100;
    opt.end_bonus = 10;
    opt.max_frag_len = 800;
    opt.max_gap = 100;
    opt.bw = opt.bw_long = 100;
    opt.pri_ratio = 0.5f;
    opt.min_cnt = 2;
    opt.min_chain_score = 25;
    opt.min_dp_max = 40;
    opt.best_n = 20;
    opt.mid_occ = 1000;
    opt.max_occ = 5000;
    opt.mini_batch_size = 50000000;

	idx_rdr = mm_idx_reader_open(fasta, &ipt, fnw);

	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", fasta, strerror(errno));
		return 1;
	}
	while ((mi = mm_idx_reader_read(idx_rdr, n_threads)) != 0) {
		if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
			fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
			mm_idx_destroy(mi);
			mm_idx_reader_close(idx_rdr);
			return 1;
		}
		fprintf(stderr, "[M::%s] loaded/built the index for %d target sequence(s)\n",
				__func__, mi->n_seq);
		mm_mapopt_update(&opt, mi);
		mm_idx_stat(mi);

		// read fastq and align to mi
		
		mm_idx_destroy(mi);
	}
	mm_idx_reader_close(idx_rdr);

	if (fflush(stdout) == EOF) {
		perror("[ERROR] failed to write the results");
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);

	printf("\n");
	graph_t *g = gfa_read(gfa);

	if (!g) {
		if (!g->nodes) free(g->nodes);
		if (!g->char_table.data) free(g->char_table.data);
		if (!g->edge_table.edges) free(g->edge_table.edges);
	}

	return 0;
}
