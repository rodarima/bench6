#ifndef UTILS_H
#define UTILS_H

#include <GASPI.h>

#include "common/heat.h"

#if defined(NBUFFER) || defined(TAGASPI)
#define NUM_GASPI_QUEUES 64
#else
#define NUM_GASPI_QUEUES 1
#endif

extern gaspi_rank_t rank;
extern gaspi_rank_t nranks;

// Segment ids for upper and lower halos
extern gaspi_segment_id_t UPPER_HALO;
extern gaspi_segment_id_t LOWER_HALO;

typedef struct {
	gaspi_queue_id_t *queues;
	gaspi_number_t numQueues;
	gaspi_group_t upperGroup;
	gaspi_group_t lowerGroup;
} GASPInfo;

void broadcastConfiguration(HeatConfiguration *conf);
void setupGaspiInfo(const HeatConfiguration *conf, int64_t rows, int64_t cols, GASPInfo *info);
void freeGaspiInfo(GASPInfo *info);

#endif // UTILS_H

