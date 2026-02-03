#ifndef EXONBLOCKS_H
#define EXONBLOCKS_H

#include <R.h>
#include <Rinternals.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <string.h>
#include <stdlib.h>

// Convert BAM query sequence to char at qpos
static inline char base_at(const uint8_t *s, int qpos) {
    return seq_nt16_str[bam_seqi(s, qpos)];
}

// Build semicolon-joined blocks from one alignment
int write_blocks_one(const bam1_t *b, kstring_t *ks_start,
                    kstring_t *ks_end, kstring_t *ks_seq);

// Fetch Z (string) tag or return NULL
static inline const char *get_tagZ(const bam1_t *b, const char *tag) {
    uint8_t *p = bam_aux_get(b, tag);
    if (!p) return NULL;
    return bam_aux2Z(p);
}

// Fetch integer tag (i/I types)
static inline long get_tagi(const bam1_t *b, const char *tag, int *ok) {
    *ok = 0;
    uint8_t *p = bam_aux_get(b, tag);
    if (!p) return 0;
    *ok = 1;
    return bam_aux2i(p);
}

#endif // EXONBLOCKS_H
