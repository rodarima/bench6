#include <TAGASPI.h>
#include <assert.h>
#include <mpi.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "api.h"
#include "common.h"

#define ENABLE_DEBUG 0
#include "log.h"

static void
check_gaspi(gaspi_return_t ret, const char *fn, const char *file, int line)
{
    if (ret == GASPI_SUCCESS)
        return;

    char *msg;
    gaspi_print_error(ret, &msg);
    die("%s (%s at %s:%d)\n", msg, fn, file, line);
}

#define CHECK(x) check_gaspi(x, #x, __FILE__, __LINE__)

/* Specific for TAGASPI */
enum { SENDSEG = 0, RECVSEG = 1 };

int *sendseg;
int *recvseg;

void
pre_init(int *argc, char **argv[])
{
    int provided;
    MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        die("MPI_THREAD_MULTIPLE not supported\n");
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

    /* Setup GASPI config */
	gaspi_config_t conf;
	CHECK(gaspi_config_get(&conf));
	conf.build_infrastructure = GASPI_TOPOLOGY_DYNAMIC;
	conf.queue_size_max = 4*1024;
    if (nqueues)
        conf.queue_num = nqueues;
	CHECK(gaspi_config_set(conf));

	CHECK(tagaspi_proc_init(GASPI_BLOCK));

    unsigned short g_rank, g_nranks;
	CHECK(gaspi_proc_rank(&g_rank));
	CHECK(gaspi_proc_num(&g_nranks));

    /* Should be the same as MPI */
    assert(g_rank == rank);
    assert(g_nranks == nranks);

    /* Setup send and recv segments */
    int maxseg = maxiter * maxnreq;

    if ((sendseg = calloc(maxseg, sizeof(int))) == NULL)
        die("calloc failed\n");

    CHECK(gaspi_segment_use(SENDSEG, sendseg, maxseg * sizeof(int),
                GASPI_GROUP_ALL, GASPI_BLOCK, 0));
    CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

    if ((recvseg = calloc(maxseg, sizeof(int))) == NULL)
        die("calloc failed\n");

    CHECK(gaspi_segment_use(RECVSEG, recvseg, maxseg * sizeof(int),
                GASPI_GROUP_ALL, GASPI_BLOCK, 0));
    CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

    /* Setup queues */
    gaspi_number_t g_nqueues;
	CHECK(gaspi_queue_num(&g_nqueues));
	CHECK(tagaspi_queue_group_create(0, 0, g_nqueues,
                GASPI_QUEUE_GROUP_POLICY_CPU_RR));

    if (nqueues != 0 && nqueues != (int) g_nqueues)
        die("requested nqueues mismatch\n");

    err("requested nqueues=%d allocated=%d\n", nqueues, (int) g_nqueues);

    nqueues = g_nqueues;
}

void
cleanup(void)
{
    if (rank == 0) {
        err("terminating gaspi, wait...\n");
        fflush(stderr);
    }

	CHECK(tagaspi_proc_term(GASPI_BLOCK));
    MPI_Finalize();
}


static void
transmit(int itask, int iter, int wanted_nreq, double *dur, int dosend,
        atomic_int *nsamples)
{
    int index = itask * maxiter + iter;
    int tag = index;

    /* Use the task index to select the queue */
    int queue = itask % nqueues;

    double t0 = get_time();
    int sample_nreq = 1 + atomic_fetch_add(&cur_req, 1);
    int dst = dstrank[iter % 2];

    if (dosend) {
        while (1) {
            gaspi_return_t ret = tagaspi_write_notify(
                    SENDSEG, index * sizeof(int), dst,
                    RECVSEG, index * sizeof(int),
                    sizeof(int), tag, 1, queue);

            if (ret == GASPI_SUCCESS)
                break;

            if (ret != GASPI_QUEUE_FULL) {
                check_gaspi(ret, "tagaspi_write_notify",
                        __FILE__, __LINE__);
            }
        }
    } else {
        while (1) {
            gaspi_return_t ret = tagaspi_notify_async_wait(RECVSEG, tag,
                    GASPI_NOTIFICATION_IGNORE);

            if (ret == GASPI_SUCCESS)
                break;

            if (ret != GASPI_QUEUE_FULL) {
                check_gaspi(ret, "tagaspi_notify_async_wait",
                        __FILE__, __LINE__);
            }
        }
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


