#ifndef UTILS_H
#define UTILS_H

#include <GASPI.h>

#include <stdlib.h>
#include <stdio.h>

#include "common/streaming.h"

#define CHECK(f)                                                          \
    {                                                                     \
        const gaspi_return_t __r = f;                                     \
        if (__r != GASPI_SUCCESS) {                                       \
            printf("Error: '%s' [%s:%i]: %i\n",#f,__FILE__,__LINE__,__r); \
            exit (EXIT_FAILURE);                                          \
        }                                                                 \
    }

extern gaspi_rank_t rank;
extern gaspi_rank_t nranks;

// Segment identifier
extern gaspi_segment_id_t RECVSEG;
extern gaspi_segment_id_t SENDSEG;

typedef struct {
	gaspi_number_t nqueues;
} GASPInfo;

void setupGaspiInfo(const StreamingConfiguration *conf, uint64_t size, GASPInfo *info);
void freeGaspiInfo(GASPInfo *info);

#endif // UTILS_H

