#include <math.h>
#include <mpi.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "api.h"

#define ENABLE_DEBUG 0
#include "log.h"

/* Command line options */
int maxnreq = 64;
int maxiter = 1000;
int maxnsamples = 500;
int onlysend = 0;
int nqueues = 0;

/* Global constants */
int rank;
int nranks;
int dstrank[2];

/* Atomic counters */
atomic_int cur_req = 0;
atomic_int nsamples_send = 0;
atomic_int nsamples_recv = 0;

static void
parse_args(int argc, char *argv[])
{
    int opt;
    while ((opt = getopt(argc, argv, "r:n:q:s")) != -1) {
        switch (opt) {
            case 'r':
                maxnreq = atoi(optarg);
                break;
            case 'n':
                maxnsamples = atoi(optarg);
                maxiter = maxnsamples * 2;
                break;
            case 's':
                onlysend = 1;
                break;
            case 'q':
                nqueues = atoi(optarg);
                break;
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-r maxnreq] [-n nsamples] [-q nqueues] [-s]\n",
                        argv[0]);
                exit(EXIT_FAILURE);
        }
    }
}

static void
run_nreq(int wanted_nreq, double *senddur, double *recvdur)
{
    nsamples_send = 0;
    nsamples_recv = 0;

    for (int i = 0; i < maxnsamples; i++) {
        senddur[i] = 0.0;
        recvdur[i] = 0.0;
    }

    if (rank == 0)
        dbg("running with nreq = %d\n", wanted_nreq);

    MPI_Barrier(MPI_COMM_WORLD);

    /* Continue running until we have enough samples */
    int cont = 1;
    for (int run = 0; cont; run++) {

        for (int i = 0; i < wanted_nreq; i++) {
            #pragma oss task label("send")
            for (int iter = 0; iter < maxiter; iter++) {
                sendblock(i, iter, wanted_nreq, senddur);
            }
        }

        #pragma oss taskwait

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < wanted_nreq; i++) {
            #pragma oss task label("recv")
            for (int iter = 0; iter < maxiter; iter++) {
                recvblock(i, iter, wanted_nreq, recvdur);
            }
        }

        #pragma oss taskwait

        MPI_Barrier(MPI_COMM_WORLD);

        /* Only the rank 0 determines if we want to continue */
        int nsend = atomic_load(&nsamples_send);
        int nrecv = atomic_load(&nsamples_recv);
        if (rank == 0) {
            if (nsend >= maxnsamples) {
                if (onlysend || nrecv >= maxnsamples) {
                    cont = 0;
                }
            } else {
                dbg("not enough samples send=%d recv=%d, repeating run %d\n",
                        nsend, nrecv, run);
            }
        }

        MPI_Bcast(&cont, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    #pragma oss taskwait

    MPI_Barrier(MPI_COMM_WORLD);

//    if (rank == 0) {
//        printf("duration for %d concurrent requests\n", wanted_nreq);
//        for (int i = 0; i < maxnsamples; i++) {
//            printf("%d %f\n", i, senddur[i]);
//        }
//    }
}

static void
print_table(int nreq, double *senddur, double *recvdur)
{
    double meansend = 0.0;
    double meanrecv = 0.0;
    for (int i = 0; i < maxnsamples; i++) {
        meansend += senddur[i];
        meanrecv += recvdur[i];
    }

    meansend /= maxnsamples;
    meanrecv /= maxnsamples;

    double stdsend = 0.0;
    double stdrecv = 0.0;

    for (int i = 0; i < maxnsamples; i++) {
        stdsend += (senddur[i] - meansend) * (senddur[i] - meansend);
        stdrecv += (recvdur[i] - meanrecv) * (recvdur[i] - meanrecv);
    }

    stdsend = sqrt(stdsend / (maxnsamples - 1));
    stdrecv = sqrt(stdrecv / (maxnsamples - 1));

    if (onlysend) {
        printf("%6d  %6d  %12e  %12e\n",
                nreq, maxnsamples,
                meansend, stdsend);
    } else {
        printf("%6d  %6d  %12e  %12e  %12e  %12e\n",
                nreq, maxnsamples,
                meansend, stdsend,
                meanrecv, stdrecv);
    }
    
}

static void
run(void)
{
    double *senddur = calloc(maxnsamples, sizeof(double));
    if (senddur == NULL)
        die("calloc failed\n");

    double *recvdur = calloc(maxnsamples, sizeof(double));
    if (recvdur == NULL)
        die("calloc failed\n");

    if (rank == 0) {
        if (onlysend) {
            printf("%6s  %6s  %12s  %12s\n",
                    "nreq", "nsamp",
                    "meansend", "stdsend");
        } else {
            printf("%6s  %6s  %12s  %12s  %12s  %12s\n",
                    "nreq", "nsamp",
                    "meansend", "stdsend",
                    "meanrecv", "stdrecv");
        }
    }

    for (int n = 1; n <= maxnreq; n++) {
        run_nreq(n, senddur, recvdur);
        if (rank == 0)
            print_table(n, senddur, recvdur);
    }

    free(senddur);
    free(recvdur);
}

int
main(int argc, char *argv[])
{
    pre_init(&argc, &argv);
    parse_args(argc, argv);
    init();
    run();
    cleanup();
}
