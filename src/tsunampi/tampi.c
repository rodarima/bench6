#include <TAMPI.h>
#include <math.h>
#include <mpi.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "api.h"
#include "common.h"

#define ENABLE_DEBUG 0
#include "log.h"

void
pre_init(int *argc, char **argv[])
{
    int provided;
    MPI_Init_thread(argc, argv, MPI_TASK_MULTIPLE, &provided);

    if (provided != MPI_TASK_MULTIPLE) {
        die("MPI_TASK_MULTIPLE not supported\n");
    }
}

void
init(void)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    /* Neighbors */
    if ((rank % 2) == 0) {
        dstrank[0] = (rank + 1) % nranks;
        dstrank[1] = (rank - 1 + nranks) % nranks;
    } else {
        dstrank[0] = (rank - 1 + nranks) % nranks;
        dstrank[1] = (rank + 1) % nranks;
    }
}

void
cleanup(void)
{
    MPI_Finalize();
}

static void
transmit(int itask, int iter, int wanted_nreq, double *dur, int dosend,
        atomic_int *nsamples)
{
    (void) itask;
    int buf[1] = { 1 };

    int tag = iter;

    double t0 = get_time();

    int sample_nreq = 1 + atomic_fetch_add(&cur_req, 1);

    if (dosend) {
        TAMPI_Isend(buf, 1, MPI_INT, dstrank[iter % 2], tag,
                MPI_COMM_WORLD);
    } else {
        TAMPI_Irecv(buf, 1, MPI_INT, dstrank[iter % 2], tag,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    atomic_fetch_sub(&cur_req, 1);

    double t1 = get_time();

    /* Reject the sample if nreq is incorrect */
    if (sample_nreq != wanted_nreq) {
        return;
    }
    
    int j = atomic_fetch_add(nsamples, 1);

    /* Too many samples already */
    if (j >= maxnsamples)
        return;

    /* Store the sample in microseconds */
    dur[j] = (t1 - t0) * 1e-3;
}

void
sendblock(int itask, int iter, int wanted_nreq, double *dur)
{
    transmit(itask, iter, wanted_nreq, dur, 1, &nsamples_send);
}

void
recvblock(int itask, int iter, int wanted_nreq, double *dur)
{
    transmit(itask, iter, wanted_nreq, dur, 0, &nsamples_recv);
}
