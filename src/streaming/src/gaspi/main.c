#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "common/streaming.h"

#ifdef TAGASPI
#include <TAGASPI.h>
#endif

#ifdef _OMPSS_2
#include <nanos6/debug.h>
#endif

gaspi_rank_t rank;
gaspi_rank_t nranks;

int main(int argc, char **argv)
{
	const int required = MPI_THREAD_MULTIPLE;

	int provided;
	MPI_Init_thread(&argc, &argv, required, &provided);
	if (provided != required) {
		fprintf(stderr, "Error: MPI threading level not supported!\n");
		return 1;
	}

	gaspi_config_t gaspiConfig;
	CHECK(gaspi_config_get(&gaspiConfig));
	gaspiConfig.build_infrastructure = GASPI_TOPOLOGY_DYNAMIC;
	gaspiConfig.queue_size_max = 4096;
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

	StreamingConfiguration conf;
	if (!rank) {
		readConfiguration(argc, argv, &conf);
		checkConfiguration(&conf, nranks);
		if (conf.verbose) printConfiguration(&conf);
	}
	broadcastConfiguration(&conf);

	uint64_t sizePerRank = conf.size/nranks;

	initialize(&conf, sizePerRank);

	GASPInfo gaspiInfo;
	setupGaspiInfo(&conf, sizePerRank, &gaspiInfo);

	if (conf.warmup) {
		CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
		solve(&conf, sizePerRank, conf.bs, 1, &gaspiInfo);
	}

	CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	// Solve the problem
	double start = getTime();
	solve(&conf, sizePerRank, conf.bs, conf.timesteps, &gaspiInfo);
	double end = getTime();

	if (!rank) {
		double throughput = conf.size*conf.timesteps/(end-start);
		throughput = throughput/1000000.0;

		double bandwidth = (sizePerRank*conf.timesteps*sizeof(double))/(end-start);
		bandwidth = 8.0*bandwidth*1e-9;

#ifdef _OMPSS_2
		const int cpus = nanos6_get_num_cpus();
#else
		const int cpus = 1;
#endif

		fprintf(stdout, "size, %ld, size/rank, %ld, bs, %ld, ranks, %d, cpus, %d, timesteps, %d, time, %f, Mupdates/s, %f, Gbit/s %f\n",
				conf.size, conf.size/nranks, conf.bs, nranks, cpus, conf.timesteps, end-start, throughput, bandwidth);
	}

	freeGaspiInfo(&gaspiInfo);

	finalize(&conf);

#ifdef TAGASPI
	CHECK(tagaspi_proc_term(GASPI_BLOCK));
#else
	CHECK(gaspi_proc_term(GASPI_BLOCK));
#endif

	MPI_Finalize();

	return 0;
}
