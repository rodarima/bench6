#include <mpi.h>
#include <GASPI.h>

#ifdef TAGASPI
#include <TAGASPI.h>
#endif

#include <assert.h>
#include <math.h>

#include "common/heat.h"
#include "gaspi/macros.h"
#include "gaspi/utils.h"

gaspi_rank_t rank;
gaspi_rank_t nranks;

gaspi_segment_id_t UPPER_HALO;
gaspi_segment_id_t LOWER_HALO;

void broadcastConfiguration(HeatConfiguration *conf)
{
	// Switch to MPI to broadcast the configuration
	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	MPI_Bcast(conf, sizeof(HeatConfiguration), MPI_BYTE, 0, MPI_COMM_WORLD);

	const int heatSourcesSize = sizeof(HeatSource)*conf->numHeatSources;
	if (rank > 0) {
		// Received heat sources pointer is not valid
		conf->heatSources = (HeatSource *) malloc(heatSourcesSize);
		assert(conf->heatSources != NULL);
	}
	MPI_Bcast(conf->heatSources, heatSourcesSize, MPI_BYTE, 0, MPI_COMM_WORLD);

	// Switch back to GASPI
	MPI_Barrier(MPI_COMM_WORLD);
}

#ifdef NBUFFER
static void setupNBufferQueues(GASPInfo *info, int ncb)
{
	info->queues = (gaspi_queue_id_t *) malloc(ncb*sizeof(gaspi_queue_id_t));
	assert(info->queues != 0);

	const int blocks = ncb-2;
	const int numQueues = info->numQueues;
	const int blocksPerQueue = blocks/numQueues;
	const int remainingBlocks = blocks%numQueues;

	int b = 1;
	for (int queue = 0; queue < numQueues; ++queue) {
		const int blocksInQueue = blocksPerQueue + (queue < remainingBlocks);
		for (int i = 0; i < blocksInQueue; ++i) {
			info->queues[b] = queue;
			++b;
		}
		if (b == ncb-1) break;
	}
}
#endif

void setupGaspiInfo(const HeatConfiguration *conf, int64_t rows, int64_t cols, GASPInfo *info)
{
	UPPER_HALO = 0;
	LOWER_HALO = 1;

	// Use barrier between segment registrations due to some limitations of GPI-2
	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	// Register upper halo segment
	CHECK(gaspi_segment_use(UPPER_HALO, &conf->matrix[0],
			2*cols*sizeof(double), GASPI_GROUP_ALL,
			GASPI_BLOCK, 0));

	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	// Register lower halo segment
	CHECK(gaspi_segment_use(LOWER_HALO, &conf->matrix[(rows-2)*cols],
			2*cols*sizeof(double), GASPI_GROUP_ALL,
			GASPI_BLOCK, 0));

	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	CHECK(gaspi_queue_num(&info->numQueues));
	assert(info->numQueues > 0);

#ifdef NBUFFER
	setupNBufferQueues(info, (cols-2)/conf->cbs+2);
#endif
}

void freeGaspiInfo(GASPInfo *info)
{
#ifdef NBUFFER
	free(info->queues);
#endif

	CHECK(gaspi_segment_delete(UPPER_HALO));
	CHECK(gaspi_segment_delete(LOWER_HALO));
}
