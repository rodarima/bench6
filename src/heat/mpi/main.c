#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "common/heat.h"

#ifdef TAMPI
#include <TAMPI.h>
#endif

#ifdef _OMPSS_2
#include <nanos6/debug.h>
#endif


void generateImage(const HeatConfiguration *conf, int64_t rows, int64_t cols, int64_t rowsPerRank);

int main(int argc, char **argv)
{
#if defined(TAMPI)
	// TAMPI+OmpSs-2 variants
	const int required = MPI_TASK_MULTIPLE;
#elif defined(_OMPSS_2)
	// MPI+OmpSs-2 variants
	const int required = MPI_THREAD_SERIALIZED;
#else
	// MPI-only variants
	const int required = MPI_THREAD_SINGLE;
#endif

	int provided;
	MPI_Init_thread(&argc, &argv, required, &provided);
	if (provided != required) {
		fprintf(stderr, "Error: MPI threading level not supported!\n");
		return 1;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

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

	initialize(&conf, rowsPerRank, cols, (rowsPerRank-2)*rank);

	void *infoptr = NULL;

#ifdef MPIRMA
	MPIRMAInfo info;
	initializeWindows(&conf, rowsPerRank, cols, &info);
	infoptr = &info;
#endif

	if (conf.warmup) {
		MPI_Barrier(MPI_COMM_WORLD);
		solve(&conf, rowsPerRank, cols, 1, infoptr);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Solve the problem
	double start = getTime();
	solve(&conf, rowsPerRank, cols, conf.timesteps, infoptr);
	double end = getTime();

	if (!rank) {
		int64_t totalElements = conf.rows*conf.cols;
		double throughput = (totalElements*conf.timesteps)/(end-start);
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

#ifdef MPIRMA
	finalizeWindows(&info);
#endif

	finalize(&conf);

	MPI_Finalize();

	return 0;
}

void generateImage(const HeatConfiguration *conf, int64_t rows, int64_t cols, int64_t rowsPerRank)
{
	int rank, nranks;
	double *auxMatrix = NULL;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

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
		writeImage(conf->imageFileName, auxMatrix, rows, cols);

		free(auxMatrix);
	}
}
