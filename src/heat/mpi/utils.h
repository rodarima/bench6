#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <mpi.h>

#include "common/heat.h"

extern int rank;
extern int nranks;

typedef struct {
	int nwins;
	MPI_Comm upperComm, lowerComm;
	MPI_Win *upperWins, *lowerWins;
} MPIRMAInfo;

void broadcastConfiguration(HeatConfiguration *configuration);
void initializeWindows(HeatConfiguration *configuration, int64_t rows, int64_t cols, MPIRMAInfo *info);
void finalizeWindows(MPIRMAInfo *info);

#endif // MPI_UTILS_H
