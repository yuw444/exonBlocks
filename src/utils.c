#define _POSIX_C_SOURCE 200809L

#include "utils.h"
#include "bam2db.h"
#include <string.h>

#include <R.h>
#include <Rinternals.h>

size_t *GetSeqInt(size_t start, size_t end, size_t step)
{
    if (step == 0) {
        return NULL;
    }

    if ((end > start && step == 0) || (end < start && step == 0)) {
        return NULL;
    }

    size_t count = (size_t)(end + step - 1 - start) / step + 1;
    size_t *seq = (size_t *)calloc(count, sizeof(size_t));
    if (!seq) return NULL;

    size_t i = 0;
    for (size_t j = start; j <= end; j += step) {
        seq[i] = j;
        i++;
    }

    return seq;
}

size_t *SampleInt(size_t *arrayIn, size_t nTotal, size_t nSample, unsigned int replace, unsigned int seed)
{
    init_genrand(seed);

    size_t *sampleOut = calloc(nSample, sizeof(size_t));
    if (!sampleOut) return NULL;

    if (replace == 0) {
        if (nSample > nTotal) {
            free(sampleOut);
            return NULL;
        }

        size_t *arrayInCopy = calloc(nTotal, sizeof(size_t));
        if (!arrayInCopy) {
            free(sampleOut);
            return NULL;
        }
        memcpy(arrayInCopy, arrayIn, nTotal * sizeof(size_t));

        if (nTotal == nSample) {
            free(sampleOut);
            return arrayInCopy;
        }

        size_t currentTotal = nTotal;
        for (size_t i = 0; i < nSample; i++) {
            size_t index = genrand_int32() % currentTotal;
            sampleOut[i] = arrayInCopy[index];
            if (index != (currentTotal - 1)) {
                arrayInCopy[index] = arrayInCopy[currentTotal - 1];
            }
            currentTotal--;
        }

        free(arrayInCopy);
    } else {
        for (size_t i = 0; i < nSample; i++) {
            sampleOut[i] = arrayIn[genrand_int32() % nTotal];
        }
    }

    return sampleOut;
}

int vsI(const void *a, const void *b)
{
    size_t va = *(const size_t *)a;
    size_t vb = *(const size_t *)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

/* ── Region registry ─────────────────────────────────────────────────────── */

static int region_cmp_chr(const void *a, const void *b)
{
    return strcmp(*(const char **)a, *(const char **)b);
}

static int region_cmp_start(const void *a, const void *b)
{
    const region_t *ra = (const region_t *)a;
    const region_t *rb = (const region_t *)b;
    int c = strcmp(ra->chr, rb->chr);
    if (c != 0) return c;
    if (ra->start < rb->start) return -1;
    if (ra->start > rb->start) return 1;
    return 0;
}

static char *dup_trim(char *s)
{
    if (!s) return NULL;
    while (*s == ' ' || *s == '\t') s++;
    size_t len = strlen(s);
    while (len > 0 && (s[len-1] == '\n' || s[len-1] == '\r' || s[len-1] == ' ' || s[len-1] == '\t'))
        s[--len] = '\0';
    return strdup(s);
}

region_registry *region_registry_build(const char *file)
{
    FILE *fp = fopen(file, "r");
    if (!fp) {
        REprintf("ERROR: Cannot open regions file %s\n", file);
        return NULL;
    }

    size_t cap = 4096, n = 0;
    region_t *items = (region_t *)malloc(cap * sizeof(region_t));
    if (!items) { fclose(fp); return NULL; }

    char buf[MAX_LINE_LENGTH];
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        buf[strcspn(buf, "\n\r")] = '\0';
        if (strlen(buf) == 0) continue;
        if (buf[0] == '#') continue;

        if (n >= cap) {
            cap *= 2;
            items = (region_t *)realloc(items, cap * sizeof(region_t));
            if (!items) { fclose(fp); return NULL; }
        }

        char *tok[8];
        int nt = 0;
        char *p = buf;
        char *save;
        while ((tok[nt] = strtok_r(p, "\t", &save)) != NULL && nt < 8) {
            p = NULL;
            nt++;
        }

        if (nt < 7) continue;

        items[n].id = dup_trim(tok[0]);
        items[n].symbol = dup_trim(tok[1]);
        items[n].chr = dup_trim(tok[2]);
        items[n].start = (hts_pos_t)atol(tok[3]);
        items[n].end = (hts_pos_t)atol(tok[4]);
        items[n].annotation1 = dup_trim(tok[5]);
        items[n].annotation2 = (nt >= 8) ? dup_trim(tok[6]) : NULL;
        n++;
    }
    fclose(fp);

    if (n == 0) { free(items); return NULL; }

    qsort(items, n, sizeof(region_t), region_cmp_start);

    /* Build per-chr index */
    char **chr_names = (char **)malloc(n * sizeof(char *));
    size_t *chr_first = (size_t *)malloc(n * sizeof(size_t));
    size_t *chr_count = (size_t *)malloc(n * sizeof(size_t));
    size_t n_chr = 0;

    for (size_t i = 0; i < n; i++) {
        if (n_chr == 0 || strcmp(items[i].chr, chr_names[n_chr - 1]) != 0) {
            chr_names[n_chr] = strdup(items[i].chr);
            chr_first[n_chr] = i;
            chr_count[n_chr] = 1;
            n_chr++;
        } else {
            chr_count[n_chr - 1]++;
        }
    }

    region_registry *reg = (region_registry *)calloc(1, sizeof(region_registry));
    if (!reg) {
        for (size_t i = 0; i < n; i++) {
            free(items[i].id); free(items[i].symbol); free(items[i].chr);
            free(items[i].annotation1); free(items[i].annotation2);
        }
        free(items); free(chr_names); free(chr_first); free(chr_count);
        return NULL;
    }

    reg->items = items;
    reg->n = n;
    reg->chr_names = chr_names;
    reg->n_chr = n_chr;
    reg->chr_first = chr_first;
    reg->chr_count = chr_count;
    return reg;
}

void region_registry_destroy(region_registry *reg)
{
    if (!reg) return;
    for (size_t i = 0; i < reg->n; i++) {
        free(reg->items[i].id);
        free(reg->items[i].symbol);
        free(reg->items[i].chr);
        free(reg->items[i].annotation1);
        free(reg->items[i].annotation2);
    }
    free(reg->items);
    for (size_t i = 0; i < reg->n_chr; i++)
        free(reg->chr_names[i]);
    free(reg->chr_names);
    free(reg->chr_first);
    free(reg->chr_count);
    free(reg);
}
