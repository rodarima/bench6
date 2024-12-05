#ifndef API_H
#define API_H

#include <stdatomic.h>

/* Command line options */
extern int maxnreq;
extern int maxiter;
extern int maxnsamples;
extern int nqueues;

/* Global constants */
extern int rank;
extern int nranks;
extern int dstrank[2];

/* Atomic counters */
extern atomic_int cur_req;
extern atomic_int nsamples_send;
extern atomic_int nsamples_recv;

/* Forward declarations */
void sendblock(int itask, int iter, int wanted_nreq, double *dur);
void recvblock(int itask, int iter, int wanted_nreq, double *dur);
void pre_init(int *argc, char **argv[]);
void init(void);
void cleanup(void);

#endif /* API_H */
