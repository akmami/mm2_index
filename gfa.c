#include "gfa.h"

#include "khashl.h"
KHASHL_MAP_INIT(KH_LOCAL, segment_idx_t, map64, uint32_t, char *, kh_hash_uint64, kh_eq_generic)

#include "kseq.h"
KSTREAM_INIT2(static, gzFile, gzread, 16384)

typedef struct
{
    uint32_t incoming;
    uint32_t outgoing;
} gfa_link;

int cmp_link(const void *a, const void *b)
{
    const gfa_link *x = (gfa_link *)a;
    const gfa_link *y = (gfa_link *)b;

    if (x->incoming < y->incoming)
        return -1;
    if (x->incoming > y->incoming)
        return 1;
    return 0;
}

graph_t *gfa_read(const char *fn)
{
    gzFile fp;
    kstring_t s = {0, 0, 0};
    kstream_t *ks;
    int dret;
    uint64_t lineno = 0;

    fp = fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(0, "r");
    if (fp == 0)
        return 0;
    ks = ks_init(fp);

    segment_idx_t *segments = map64_init();

    kvec_t(gfa_link) links;
    kv_init(links);

    uint32_t max_id = 0;
    uint64_t char_total = 0;

    while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0)
    {
        if (s.l < 3 || s.s[1] != '\t')
            continue; // empty line
        if (s.s[0] == 'S')
        {
            char *saveptr;
            char *token = strtok_r(s.s, "\t", &saveptr); // S

            token = strtok_r(NULL, "\t", &saveptr); // ID
            uint32_t id = strtoul(token, NULL, 10);

            token = strtok_r(NULL, "\t", &saveptr); // SEQ
            uint32_t len = strlen(token);
            char *seq = (char *)malloc(len + 1);
            memcpy(seq, token, len + 1);

            char_total += len;

            khint_t k;
            int absent;
            k = map64_put(segments, id, &absent);
            kh_val(segments, k) = seq;

            if (id > max_id)
                max_id = id;
        }
        else if (s.s[0] == 'L')
        {
            uint32_t id1, id2;
            char strand1, strand2;
            size_t overlap = 0;
            sscanf(s.s, "L\t%u\t%c\t%u\t%c\t%luM", &id1, &strand1, &id2, &strand2, &overlap);
            gfa_link l = (gfa_link){id1, id2};
            kv_push(gfa_link, NULL, links, l);
        }
        lineno++;
    }

    graph_t *g = (graph_t *)calloc(1, sizeof(graph_t));

    // edge table
    // edge order should be compacted with the same order that nodes appear
    qsort(links.a, links.n, sizeof(gfa_link), cmp_link);

    g->edge_table.n_edges = links.n;
    g->edge_table.edges = (uint32_t *)malloc(sizeof(uint32_t) * links.n);

    // nodes table
    g->n_nodes = max_id;
    g->nodes = (node_entry_t *)calloc(max_id + 1, sizeof(node_entry_t));

    // char table
    g->char_table.length = char_total;
    g->char_table.data = (char_data *)malloc(char_total);
    uint32_t char_offset = 0;
    uint32_t edge_offset = 0;

    for (uint32_t i = 1; i <= max_id; i++)
    {
        khint_t k = map64_get(segments, i);
        if (k == kh_end(segments))
            continue;

        char *seq = kh_val(segments, k);
        uint32_t len = strlen(seq);

        g->nodes[i].seq_len = len;
        g->nodes[i].char_offset = char_offset;

        memcpy(g->char_table.data + char_offset, seq, len);
        char_offset += len;

        // put edge offset and scan iteratively the outgoing edges (edge_offset++)
        g->nodes[i].edge_offset = edge_offset;
        g->nodes[i].edge_count = 0;

        while (edge_offset < links.n && links.a[edge_offset].incoming == i)
        {

            g->edge_table.edges[edge_offset] = links.a[edge_offset].outgoing;

            edge_offset++;
            g->nodes[i].edge_count++;
        }
    }

    fprintf(stderr, "[M::%s] gfa: S - %u, L - %u\n", __func__, g->n_nodes, g->edge_table.n_edges);

    for (khint_t k = 0; k < kh_end(segments); ++k)
        if (kh_exist(segments, k))
            free(kh_val(segments, k));

    ks_destroy(ks);

    map64_destroy(segments);
    kv_destroy(links);

    gzclose(fp);
    return g;
}

void free_graph(graph_t *g)
{
    if (!g)
    {
        if (!g->nodes)
            free(g->nodes);
        if (!g->char_table.data)
            free(g->char_table.data);
        if (!g->edge_table.edges)
            free(g->edge_table.edges);
        free(g);
        g = NULL;
    }
}