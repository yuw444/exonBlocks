#ifndef UTILS_DS_H
#define UTILS_DS_H

#include <stdlib.h>
#include <stdio.h>
#include "mt19937ar.h"

/**
 * @brief Generate a sequence of integers from start to end with given step.
 */
size_t *GetSeqInt(size_t start, size_t end, size_t step);

/**
 * @brief Sample nSample integers from arrayIn (with/without replacement).
 */
size_t *SampleInt(size_t *arrayIn, size_t nTotal, size_t nSample, unsigned int replace, unsigned int seed);

/**
 * @brief Comparison function for qsort on size_t arrays (ascending).
 */
int vsI(const void *a, const void *b);

#endif
