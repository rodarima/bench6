#include <cassert>
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <algorithm>

#include <mpi.h>

static void axpy_task(double *x, double *y, double alpha, long N)
{
	for (long i=0; i < N; ++i) {
		y[i] += alpha * x[i];
	}
}

static void axpy(double *x, double *y, double alpha, long N, long TS)
{
	for (long i=0; i < N; i+=TS) {
		#pragma oss task label("axpy_task") in(x[i]) inout(y[i])
		axpy_task(x+i, y+i, alpha, std::min(TS, N-i));
	}
}

static void multisaxpy(double *x, double *y, double alpha, long N, long TS, long its)
{
	(void) alpha;

	for (long iteration=0; iteration < its; iteration++) {
		#pragma oss task label("axpy") weakin({x[i*TS], i=0; N/TS}) weakinout({y[i*TS], i=0; N/TS})
		axpy(x, y, 1.0, N, TS);
	}
	#pragma oss taskwait
}


static void initialize(double *data, double value, long N, long)
{
	for (long i=0; i < N; i++) {
		data[i] = value;
	}
}

int main(int argc, char **argv)
{
	long n   = 2000000L;
	long ts  = 50000L;
	long its = 500L;

	if (argc > 4 || (argc >= 2 && strcmp(argv[1], "-h") == 0)) {
		std::cerr << "[USAGE] " << argv[0] << " [-h] elements chunksize iterations" << std::endl;
		return 1;
	}

	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	assert(provided == MPI_THREAD_MULTIPLE);

	int rank, rank_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	assert(rank_size > 0);

	if (argc >= 2)
		n = atol(argv[1]);
	if (argc >= 3)
		ts = atol(argv[2]);
	if (argc >= 4)
		its = atol(argv[3]);

	if (n % rank_size != 0) {
		std::cerr << "elements must be divisable by the number of MPI processes" << std::endl;
		return 1;
	}

	double *x = nullptr, *y = nullptr;
	if (rank == 0) {
		x = new double[n];
		y = new double[n];
		initialize(x, 1.0, n, ts);
		initialize(y, 0.0, n, ts);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	unsigned long long int local_size = n / rank_size;
	double *local_x = new double[local_size];
	double *local_y = new double[local_size];

	MPI_Scatter(x, local_size, MPI_DOUBLE, local_x, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(y, local_size, MPI_DOUBLE, local_y, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);

	multisaxpy(local_x, local_y, 1.0, local_size, ts, its);

	clock_gettime(CLOCK_MONOTONIC, &end);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(local_y, local_size, MPI_DOUBLE, y, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double duration = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
		duration /= 1000000;

		double performance = n;
		performance *= its;
		performance = performance / duration;
		performance /= 1000000000;

		printf("%14e %14e %14ld %14ld %14ld\n", duration, performance, n, ts, its);
	}

	MPI_Finalize();
	return 0;
}
