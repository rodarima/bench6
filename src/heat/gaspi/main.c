#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "macros.h"
#include "common/heat.h"

#ifdef TAGASPI
#include <TAGASPI.h>
#endif

#ifdef _OMPSS_2
#include <nanos6/debug.h>
#endif


void generateImage(const HeatConfiguration *conf, int64_t rows, int64_t cols, int64_t rowsPerRank);

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided != MPI_THREAD_MULTIPLE) {
		fprintf(stderr, "Error: MPI threading level not supported!\n");
		return 1;
	}

	gaspi_config_t gaspiConfig;
	CHECK(gaspi_config_get(&gaspiConfig));
	gaspiConfig.build_infrastructure = GASPI_TOPOLOGY_DYNAMIC;
	gaspiConfig.queue_size_max = 2048;
	CHECK(gaspi_config_set(gaspiConfig));

#ifdef TAGASPI
	CHECK(tagaspi_proc_init(GASPI_BLOCK));
#else
	CHECK(gaspi_proc_init(GASPI_BLOCK));
#endif

	CHECK(gaspi_proc_rank(&rank));
	CHECK(gaspi_proc_num(&nranks));

	// Commit the default group
	CHECK(gaspi_group_commit(GASPI_GROUP_ALL, GASPI_BLOCK));

	HeatConfiguration conf;
	if (!rank) {
		readConfiguration(argc, argv, &conf);
		refineConfiguration(&conf, nranks*conf.rbs, conf.cbs);
		if (conf.verbose) printConfiguration(&conf);
	}
	broadcastConfiguration(&conf);

	int64_t rows = conf.rows+2;
	int64_t cols = conf.cols+2;
	int64_t rowsPerRank = conf.rows/nranks+2;

	int err = initialize(&conf, rowsPerRank, cols, (rowsPerRank-2)*rank);
	assert(!err);

	GASPInfo gaspiInfo;
	setupGaspiInfo(&conf, rowsPerRank, cols, &gaspiInfo);

	if (conf.warmup) {
		CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
		solve(&conf, rowsPerRank, cols, 1, &gaspiInfo);
	}

	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	// Solve the problem
	double start = getTime();
	solve(&conf, rowsPerRank, cols, conf.timesteps, &gaspiInfo);
	double end = getTime();

	if (!rank) {
		int64_t totalElements = conf.rows*conf.cols;
		double throughput = totalElements*conf.timesteps/(end-start);
		throughput = throughput/1000000.0;

#ifdef _OMPSS_2
		int threads = nanos6_get_num_cpus();
#else
		int threads = 1;
#endif

		fprintf(stdout, "rows, %ld, cols, %ld, rows/rank, %ld, total, %ld, total/rank, %ld, rbs, %d, "
				"cbs, %d, ranks, %d, threads, %d, timesteps, %d, time, %f, Mupdates/s, %f\n",
				conf.rows, conf.cols, conf.rows/nranks, totalElements, totalElements/nranks,
				conf.rbs, conf.cbs, nranks, threads, conf.timesteps, end-start, throughput);
	}

	if (conf.generateImage) {
		generateImage(&conf, rows, cols, rowsPerRank);
	}

	freeGaspiInfo(&gaspiInfo);

	err = finalize(&conf);
	assert(!err);

#ifdef TAGASPI
	CHECK(tagaspi_proc_term(GASPI_BLOCK));
#else
	CHECK(gaspi_proc_term(GASPI_BLOCK));
#endif

	MPI_Finalize();

	return 0;
}

void generateImage(const HeatConfiguration *conf, int64_t rows, int64_t cols, int64_t rowsPerRank)
{
	double *auxMatrix = NULL;

	if (!rank) {
		auxMatrix = (double *) malloc(rows*cols*sizeof(double));
		if (auxMatrix == NULL) {
			fprintf(stderr, "Memory cannot be allocated!\n");
			exit(1);
		}

		initializeMatrix(conf, auxMatrix, rows, cols, 0);
	}

	int count = (rowsPerRank-2)*cols;
	MPI_Gather(
		&conf->matrix[cols], count, MPI_DOUBLE,
		&auxMatrix[cols], count, MPI_DOUBLE,
		0, MPI_COMM_WORLD
	);

	if (!rank) {
		int err = writeImage(conf->imageFileName, auxMatrix, rows, cols);
		assert(!err);

		free(auxMatrix);
	}
}
