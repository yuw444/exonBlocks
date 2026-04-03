#define _POSIX_C_SOURCE 200809L

#include "utils.h"
#include <string.h>

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
