#ifndef BAM2DB_H
#define BAM2DB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <htslib/sam.h>
#include <zlib.h>
#include <sqlite3.h>
#include "hashtable.h"

#define MAX_LINE_LENGTH 4096
#define BASE_BITS 2
#define BYTE_SIZE 8
#define BYTE_MASK 0xFF
#define BASE_MASK ((1 << BASE_BITS) - 1)

typedef enum { PROTO_10X, PROTO_ZUMIS } protocol_t;

typedef struct {
    const char *cb;
    const char *umi;
} read_identity;

typedef struct {
    char *id;
    char *symbol;
    char *chr;
    hts_pos_t start;
    hts_pos_t end;
    char *annotation1;
    char *annotation2;
} region_t;

typedef struct {
    region_t *items;
    size_t n;
    char **chr_names;
    size_t n_chr;
    size_t *chr_first;
    size_t *chr_count;
} region_registry;

uint8_t encode_base(char base);
uint8_t *encode_DNA(const char *DNA_seq);
char *decode_DNA(uint8_t *encoded_DNA, size_t n);

region_registry *region_registry_build(const char *file);
void region_registry_destroy(region_registry *reg);

int bam2db(
    char *bam_file,
    char *db_file,
    char *path_out,
    char *barcodes_file,
    char *regions_file,
    float rate_cell,
    float rate_depth,
    unsigned int seed,
    int umi_copies_flag,
    protocol_t protocol,
    const char *region_chr,
    hts_pos_t region_start,
    hts_pos_t region_end);

int table2gz(
    sqlite3 *db_handle,
    const char *table_name,
    gzFile gz_file_ptr,
    unsigned int header,
    const char *delim);

size_t nrow_sql_table(
    sqlite3 *db_handle,
    const char *table_name);

#endif
