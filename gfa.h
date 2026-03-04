#ifndef GFA_H
#define GFA_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "minimap.h"
#include "sketch.h"
#include "kseq.h"
#include "khashl.h"

// Node table entry
typedef struct {
    uint32_t seq_len;
    uint32_t char_offset;
    uint32_t edge_count;
    uint32_t edge_offset;
} node_entry_t;

// Character table
typedef struct {
    uint8_t *data;   // 2-bit packed bases
    uint64_t length;
} char_table_t;

// Edge table
typedef struct {
    uint32_t *edges;
    uint64_t n_edges;
} edge_table_t;

// Full graph container
typedef struct {
    node_entry_t *nodes;
    uint32_t n_nodes;

    char_table_t char_table;
    edge_table_t edge_table;
} graph_t;


graph_t *gfa_read(const char *fn);


#endif
