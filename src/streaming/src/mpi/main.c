#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "common/streaming.h"

#ifdef TAMPI
#include <TAMPI.h>
#endif

#ifdef _OMPSS_2
#include <nanos6/debug.h>
#endif

int rank;
int nranks;

int main(int argc, char **argv)
{
#ifdef TAMPI
#ifdef TAMPI_BLOCKING
	// TAMPI blocking
	const int required = MPI_TASK_MULTIPLE;
#else
	// TAMPI non-blocking
	const int required = MPI_THREAD_MULTIPLE;
#endif
#else
	// MPI-only variant
	const int required = MPI_THREAD_SINGLE;
#endif

	int provided;
	MPI_Init_thread(&argc, &argv, required, &provided);
	if (provided != required) {
		fprintf(stderr, "Error: MPI threading level not supported!\n");
		return 1;
	}

	// Optimizations for the MPI multi-threaded. Notify the MPI library that we
	// are not going to use the MPI_ANY_TAG and MPI_ANY_SOURCE features. Also
	// notify that messages and their corresponding receive buffers have the
	// exact length
	MPI_Info info;
	MPI_Info_create(&info);
	MPI_Info_set(info, "mpi_assert_no_any_tag", "true");
	MPI_Info_set(info, "mpi_assert_no_any_source", "true");
	MPI_Info_set(info, "mpi_assert_exact_length", "true");
	MPI_Comm_set_info(MPI_COMM_WORLD, info);
	MPI_Info_free(&info);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	StreamingConfiguration conf;
	if (!rank) {
		readConfiguration(argc, argv, &conf);
		checkConfiguration(&conf, nranks);
		if (conf.verbose) printConfiguration(&conf);
	}
	broadcastConfiguration(&conf);

	uint64_t sizePerRank = conf.size/nranks;

	initialize(&conf, sizePerRank);

	if (conf.warmup) {
		MPI_Barrier(MPI_COMM_WORLD);
		solve(&conf, sizePerRank, conf.bs, 1, NULL);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Solve the problem
	double start = getTime();
	solve(&conf, sizePerRank, conf.bs, conf.timesteps, NULL);
	double end = getTime();

	if (!rank) {
		double throughput = (conf.size*conf.timesteps)/(end-start);
		throughput = throughput/1000000.0;

		double bandwidth = (sizePerRank*conf.timesteps*sizeof(double))/(end-start);
		bandwidth = 8.0*bandwidth*1e-9;

//#ifdef _OMPSS_2
//		const int cpus = nanos6_get_num_cpus();
//#else
//		const int cpus = 1;
//#endif

		double t = end - start;
		//fprintf(stdout, "size, %ld, size/rank, %ld, bs, %ld, ranks, %d, cpus, %d, timesteps, %d, time, %f, Mupdates/s, %f, Gbits/s, %f\n",
		//		conf.size, conf.size/nranks, conf.bs, nranks, cpus, conf.timesteps, end-start, throughput, bandwidth);
		printf("%14e %14e %14e %14ld %14ld %14d\n",
				t, throughput, bandwidth, conf.size, conf.bs, conf.timesteps);
	}

	finalize(&conf);

	MPI_Finalize();

	return 0;
}
