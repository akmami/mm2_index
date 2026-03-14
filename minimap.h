#ifndef MINIMAP2_H
#define MINIMAP2_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

#define MM_VERSION "2.30-r1290-dirty"

#define MM_F_NO_DIAG       (0x001LL) // no exact diagonal hit
#define MM_F_NO_DUAL       (0x002LL) // skip pairs where query name is lexicographically larger than target name
#define MM_F_CIGAR         (0x004LL)
#define MM_F_OUT_SAM       (0x008LL)
#define MM_F_NO_QUAL       (0x010LL)
#define MM_F_OUT_CG        (0x020LL)
#define MM_F_OUT_CS        (0x040LL)
#define MM_F_SPLICE        (0x080LL) // splice mode
#define MM_F_SPLICE_FOR    (0x100LL) // match GT-AG
#define MM_F_SPLICE_REV    (0x200LL) // match CT-AC, the reverse complement of GT-AG
#define MM_F_NO_LJOIN      (0x400LL)
#define MM_F_OUT_CS_LONG   (0x800LL)
#define MM_F_SR            (0x1000LL)
#define MM_F_FRAG_MODE     (0x2000LL)
#define MM_F_NO_PRINT_2ND  (0x4000LL)
#define MM_F_2_IO_THREADS  (0x8000LL)
#define MM_F_LONG_CIGAR    (0x10000LL)
#define MM_F_INDEPEND_SEG  (0x20000LL)
#define MM_F_SPLICE_FLANK  (0x40000LL)
#define MM_F_SOFTCLIP      (0x80000LL)
#define MM_F_FOR_ONLY      (0x100000LL)
#define MM_F_REV_ONLY      (0x200000LL)
#define MM_F_HEAP_SORT     (0x400000LL)
#define MM_F_ALL_CHAINS    (0x800000LL)
#define MM_F_OUT_MD        (0x1000000LL)
#define MM_F_COPY_COMMENT  (0x2000000LL)
#define MM_F_EQX           (0x4000000LL) // use =/X instead of M
#define MM_F_PAF_NO_HIT    (0x8000000LL) // output unmapped reads to PAF
#define MM_F_NO_END_FLT    (0x10000000LL)
#define MM_F_HARD_MLEVEL   (0x20000000LL)
#define MM_F_SAM_HIT_ONLY  (0x40000000LL)
#define MM_F_RMQ           (0x80000000LL)
#define MM_F_QSTRAND       (0x100000000LL)
#define MM_F_NO_INV        (0x200000000LL)
#define MM_F_NO_HASH_NAME  (0x400000000LL)
#define MM_F_SPLICE_OLD    (0x800000000LL)
#define MM_F_SECONDARY_SEQ (0x1000000000LL)	//output SEQ field for seqondary alignments using hard clipping
#define MM_F_OUT_DS        (0x2000000000LL)
#define MM_F_WEAK_PAIRING  (0x4000000000LL)
#define MM_F_SR_RNA        (0x8000000000LL)
#define MM_F_OUT_JUNC      (0x10000000000LL)

#define MM_I_HPC          0x1
#define MM_I_NO_SEQ       0x2
#define MM_I_NO_NAME      0x4

#define MM_IDX_MAGIC   "MMI\2"

#define MM_MAX_SEG       255

#define MM_CIGAR_MATCH      0
#define MM_CIGAR_INS        1
#define MM_CIGAR_DEL        2
#define MM_CIGAR_N_SKIP     3
#define MM_CIGAR_SOFTCLIP   4
#define MM_CIGAR_HARDCLIP   5
#define MM_CIGAR_PADDING    6
#define MM_CIGAR_EQ_MATCH   7
#define MM_CIGAR_X_MISMATCH 8

#define MM_CIGAR_STR  "MIDNSHP=XB"

#ifdef __cplusplus
extern "C" {
#endif

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;

// minimap2 index
typedef struct {
	char *name;      // name of the db sequence
	uint64_t offset; // offset in mm_idx_t::S
	uint32_t len;    // length
	uint32_t is_alt;
} mm_idx_seq_t;

typedef struct {
	int32_t b, w, k, flag;
	uint32_t n_seq;            // number of reference sequences
	int32_t index;
	int32_t n_alt;
	mm_idx_seq_t *seq;         // sequence name, length and offset
	uint32_t *S;               // 4-bit packed sequence
	struct mm_idx_bucket_s *B; // index (hidden)
	struct mm_idx_intv_s *I;   // intervals (hidden)
	struct mm_idx_spsc_s *spsc;// splice score (hidden)
	struct mm_idx_jjump_s *J;  // junctions to create jumps (hidden)
	void *km, *h;
} mm_idx_t;

// minimap2 alignment
typedef struct {
	uint32_t capacity;                  // the capacity of cigar[]
	int32_t dp_score, dp_max, dp_max2;  // DP score; score of the max-scoring segment; score of the best alternate mappings
	int32_t dp_max0;                    // DP score before mm_update_dp_max() adjustment
	uint32_t n_ambi:30, trans_strand:2; // number of ambiguous bases; transcript strand: 0 for unknown, 1 for +, 2 for -
	uint32_t n_cigar;                   // number of cigar operations in cigar[]
	uint32_t cigar[];
} mm_extra_t;

typedef struct {
	int32_t id;             // ID for internal uses (see also parent below)
	int32_t cnt;            // number of minimizers; if on the reverse strand
	int32_t rid;            // reference index; if this is an alignment from inversion rescue
	int32_t score;          // DP alignment score
	int32_t qs, qe, rs, re; // query start and end; reference start and end
	int32_t parent, subsc;  // parent==id if primary; best alternate mapping score
	int32_t as;             // offset in the a[] array (for internal uses only)
	int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
	int32_t n_sub;          // number of suboptimal mappings
	int32_t score0;         // initial chaining score (before chain merging/spliting)
	uint32_t mapq:8, split:2, rev:1, inv:1, sam_pri:1, proper_frag:1, pe_thru:1, seg_split:1, seg_id:8, split_inv:1, is_alt:1, strand_retained:1, is_spliced:1, dummy:4;
	uint32_t hash;
	float div;
	mm_extra_t *p;
} mm_reg1_t;

// indexing and mapping options
typedef struct {
	short k, w, flag, bucket_bits;
	int64_t mini_batch_size;
	uint64_t batch_size;
} mm_idxopt_t;

typedef struct {
	int64_t flag;    // see MM_F_* macros
	int seed;
	int sdust_thres; // score threshold for SDUST; 0 to disable

	int max_qlen;    // max query length

	int bw, bw_long; // bandwidth
	int max_gap, max_gap_ref; // break a chain if there are no minimizers in a max_gap window
	int max_frag_len;
	int max_chain_skip, max_chain_iter;
	int min_cnt;         // min number of minimizers on each chain
	int min_chain_score; // min chaining score
	float chain_gap_scale;
	float chain_skip_scale;
	int rmq_size_cap, rmq_inner_dist;
	int rmq_rescue_size;
	float rmq_rescue_ratio;

	float mask_level;
	int mask_len;
	float pri_ratio;
	int best_n;      // top best_n chains are subjected to DP alignment

	float alt_drop;

	int a, b, q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext penalties
	int transition; // transition mismatch score (A:G, C:T)
	int sc_ambi; // score when one or both bases are "N"
	int noncan;      // cost of non-canonical splicing sites
	int junc_bonus;  // bonus for a splice site in annotation
	int junc_pen;    // penalty for GT- or -AG not scored in --spsc
	int zdrop, zdrop_inv;   // break alignment if alignment score drops too fast along the diagonal
	int end_bonus;
	int min_dp_max;  // drop an alignment if the score of the max scoring segment is below this threshold
	int min_ksw_len;
	int anchor_ext_len, anchor_ext_shift;
	float max_clip_ratio; // drop an alignment if BOTH ends are clipped above this ratio

	int rank_min_len;
	float rank_frac;

	int pe_ori, pe_bonus;

	int32_t jump_min_match;

	float mid_occ_frac;  // only used by mm_mapopt_update(); see below
	float q_occ_frac;
	int32_t min_mid_occ, max_mid_occ;
	int32_t mid_occ;     // ignore seeds with occurrences above this threshold
	int32_t max_occ, max_max_occ, occ_dist;
	int64_t mini_batch_size; // size of a batch of query bases to process in parallel
	int64_t max_sw_mat;
	int64_t cap_kalloc;

	const char *split_prefix;
} mm_mapopt_t;

// index reader
typedef struct {
	int is_idx, n_parts;
	int64_t idx_size;
	mm_idxopt_t opt;
	FILE *fp_out;
	union {
		struct mm_bseq_file_s *seq;
		FILE *idx;
	} fp;
} mm_idx_reader_t;

// memory buffer for thread-local storage during mapping
struct mm_tbuf_s {
	void *km;
	int rep_len, frag_gap;
};

typedef struct mm_tbuf_s mm_tbuf_t;

// seeds
#define get_hash(minimizer) ((minimizer).x >> 8)
#define get_span(minimizer) ((minimizer).x & 0xFF)
#define get_rid(minimizer) ((minimizer).y >> 32)
#define get_index(minimizer) ((((minimizer).y & 0xFFFFFFFF) >> 1) - get_span(minimizer) + 1)
#define get_end(minimizer) (((minimizer).y & 0xFFFFFFFF) >> 1)
#define get_strand(minimizer) ((minimizer).y & 1)

#define get_seed_rid(minimizer) ((minimizer) >> 32)
#define get_seed_index(minimizer, span) ((((minimizer) & 0xFFFFFFFF) >> 1) - span + 1)
#define get_seed_end(minimizer) (((minimizer) & 0xFFFFFFFF) >> 1)
#define get_seed_strand(minimizer) ((minimizer) & 1)

#ifdef __cplusplus
}
#endif

#endif // MINIMAP2_H
