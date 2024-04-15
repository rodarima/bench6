#include <mpi.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils_mpi.h"
#include "heat.h"

int rank;
int nranks;

void broadcastConfiguration(HeatConfiguration *conf)
{
	MPI_Bcast(conf, sizeof(HeatConfiguration), MPI_BYTE, 0, MPI_COMM_WORLD);

	const int heatSourcesSize = sizeof(HeatSource)*conf->numHeatSources;
	if (rank > 0) {
		// Received heat sources pointer is not valid
		conf->heatSources = (HeatSource *) malloc(heatSourcesSize);
		assert(conf->heatSources != NULL);
	}
	MPI_Bcast(conf->heatSources, heatSourcesSize, MPI_BYTE, 0, MPI_COMM_WORLD);
}

void initializeWindows(HeatConfiguration *conf, int64_t rows, int64_t cols, MPIRMAInfo *info)
{
	assert(info != NULL);

	double (*matrix)[cols] = (double (*)[cols]) conf->matrix;
	const int cbs = conf->cbs;
	const int ncb = (cols-2)/cbs;

	if (ncb % conf->nwins != 0) {
		fprintf(stdout, "Error: Column blocks are not divisible by the number of windows\n");
		exit(1);
	}

#ifdef _OMPSS_2
	const int nrealwins = conf->nwins;
#else
	// MPI-only variants may have multiple logical windows but they only have
	// a single physical window, since they cannot benefit from multiple physical
	// windows in each halo
	const int nrealwins = 1;
#endif

	const int ncbxw = ncb / nrealwins;

	info->upperWins = (MPI_Win *) malloc(nrealwins*sizeof(MPI_Win));
	info->lowerWins = (MPI_Win *) malloc(nrealwins*sizeof(MPI_Win));
	assert(info->upperWins && info->lowerWins);

	// Create a communicator for exchanging halos with the upper and lower rank. The first
	// and last ranks are special cases since they have a communicator with no other rank
	// (e.g. upper communicator for rank zero is empty)
	MPI_Comm_split(MPI_COMM_WORLD, rank/2, rank,
			(rank%2 == 0) ? &info->lowerComm : &info->upperComm);
	MPI_Comm_split(MPI_COMM_WORLD, (rank+1)/2, rank,
			(rank%2 == 0) ? &info->upperComm : &info->lowerComm);

	// Register the MPI windows for the upper halo (first row) and the lower halo (last row)
	// Each halo window can hold more than one BS block. Communicating though these windows
	// will be fine-grained because they are using the special communicators (two ranks)
	const size_t winSize = ncbxw*cbs*sizeof(double);
	for (int w = 0; w < nrealwins; ++w) {
		const int64_t colOffset = 1+w*ncbxw*cbs;
		MPI_Win_create(&matrix[0][colOffset], winSize, sizeof(double),
				MPI_INFO_NULL, info->upperComm, &info->upperWins[w]);
		MPI_Win_create(&matrix[rows-1][colOffset], winSize, sizeof(double),
				MPI_INFO_NULL, info->lowerComm, &info->lowerWins[w]);
	}

	info->nwins = conf->nwins;
}

void finalizeWindows(MPIRMAInfo *info)
{
	assert(info != NULL);

#ifdef _OMPSS_2
	const int nrealwins = info->nwins;
#else
	const int nrealwins = 1;
#endif

	for (int w = 0; w < nrealwins; ++w) {
		MPI_Win_free(&info->upperWins[w]);
		MPI_Win_free(&info->lowerWins[w]);
	}

	MPI_Comm_free(&info->upperComm);
	MPI_Comm_free(&info->lowerComm);

	free(info->upperWins);
	free(info->lowerWins);
}
