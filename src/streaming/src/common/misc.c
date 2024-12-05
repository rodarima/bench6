#define _POSIX_C_SOURCE 200809L

#include <mpi.h>

#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "streaming.h"

static void printUsage(int argc, char **argv)
{
	(void) argc;
	fprintf(stdout, "Usage: %s <-s size> <-b blocksize> <-t timesteps> [OPTION]...\n", argv[0]);
	fprintf(stdout, "Parameters:\n");
	fprintf(stdout, "  -s, --size=SIZE          use SIZE as the global size of the array\n");
	fprintf(stdout, "  -b, --bs=BS              use BS as the number of elements of each block\n");
	fprintf(stdout, "  -t, --timesteps=TS       use TS as the number of timesteps\n");
	fprintf(stdout, "Optional parameters:\n");
	fprintf(stdout, "  -o, --offset=OFFSET      use OFFSET as the rank offset for communication (default: 1)\n");
	fprintf(stdout, "  -W, --no-warmup          do not perform warmup timestep (warmup enabled by default)\n");
	fprintf(stdout, "  -c, --compfactor=FACTOR     use FACTOR as the factor to calculate the computation weight (default: 1.0)\n");
	fprintf(stdout, "  -v, --verbose            display additional information (disabled by default)\n");
	fprintf(stdout, "  -h, --help               display this help and exit\n\n");
}

static void setDefaultConfiguration(StreamingConfiguration *conf)
{
	conf->size = 0;
	conf->bs = 0;
	conf->timesteps = 0;
	conf->array = NULL;
	conf->source = NULL;
	conf->compfactor = 1.0;
	conf->offset = 1;
	conf->warmup = true;
	conf->verbose = false;
}

static void readParameters(int argc, char **argv, StreamingConfiguration *conf)
{
	static struct option long_options[] = {
		{"size",         required_argument,  0, 's'},
		{"bs",           required_argument,  0, 'b'},
		{"timesteps",    required_argument,  0, 't'},
		{"offset",       required_argument,  0, 'o'},
		{"compfactor",   required_argument,  0, 'c'},
		{"no-warmup",    no_argument,        0, 'W'},
		{"verbose",      no_argument,        0, 'v'},
		{"help",         no_argument,        0, 'h'},
		{0, 0, 0, 0}
	};

	int c, index;
	char *ptr;
	while ((c = getopt_long(argc, argv, "hs:t:b:o:c:Wv", long_options, &index)) != -1) {
		switch (c) {
			case 'h':
				printUsage(argc, argv);
				exit(0);
			case 'v':
				conf->verbose = true;
				break;
			case 'W':
				conf->warmup = false;
				break;
			case 'c':
				conf->compfactor = atof(optarg);
				assert(conf->compfactor >= 0);
				break;
			case 's':
				conf->size = strtoull(optarg, &ptr, 10);
				assert(*ptr == '\0');
				assert(conf->size > 0);
				break;
			case 'b':
				conf->bs = strtoull(optarg, &ptr, 10);
				assert(*ptr == '\0');
				assert(conf->bs > 0);
				break;
			case 't':
				conf->timesteps = atoi(optarg);
				assert(conf->timesteps > 0);
				break;
			case 'o':
				conf->offset = atoi(optarg);
				assert(conf->offset > 0);
				break;
			case '?':
				exit(1);
			default:
				abort();
		}
	}

	if (!conf->size || !conf->timesteps || !conf->bs) {
		printUsage(argc, argv);
		exit(1);
	}
}

static void initializeArray(double *array, uint64_t size, double value)
{
	for (uint64_t e = 0; e < size; ++e) {
		array[e] = value;
	}
}

void initialize(StreamingConfiguration *conf, uint64_t size)
{
	int pagesize = sysconf(_SC_PAGESIZE);
	assert(pagesize > 0);

	int err = posix_memalign((void **)&conf->array, pagesize, size*sizeof(double));
	err |= posix_memalign((void **)&conf->source, pagesize, size*sizeof(double));
	if (err || conf->array == NULL || conf->source == NULL) {
		fprintf(stderr, "Error: Memory cannot be allocated!\n");
		exit(1);
	}

	initializeArray(conf->array, size, 0.0);
	initializeArray(conf->source, size, 1.0);
}

void finalize(StreamingConfiguration *conf)
{
	assert(conf->array != NULL);
	assert(conf->source != NULL);
	free(conf->array);
	free(conf->source);
	conf->array = NULL;
	conf->source = NULL;
}

void readConfiguration(int argc, char **argv, StreamingConfiguration *conf)
{
	setDefaultConfiguration(conf);
	readParameters(argc, argv, conf);
}

void broadcastConfiguration(StreamingConfiguration *conf)
{
	MPI_Bcast(conf, sizeof(StreamingConfiguration), MPI_BYTE, 0, MPI_COMM_WORLD);
}

void checkConfiguration(const StreamingConfiguration *conf, int totalranks)
{
	assert(conf->size > 0);
	assert(conf->bs > 0);
	assert(conf->timesteps > 0);
	assert(totalranks > 0);

	if (conf->size%totalranks) {
		fprintf(stderr, "Error: The array size (%ld) is not divisible by the number of ranks %d\n", conf->size, totalranks);
		exit(1);
	}
}

void printConfiguration(const StreamingConfiguration *conf)
{
	fprintf(stdout, "Size              : %ld\n", conf->size);
	fprintf(stdout, "Block size        : %ld\n", conf->bs);
	fprintf(stdout, "Timesteps         : %d\n", conf->timesteps);
}

void checkResult(const StreamingConfiguration *conf, uint64_t size, int rank)
{
	if (conf->warmup) {
		if (!rank)
			fprintf(stderr, "Warning: Skipping result check because of the warmup\n");
		return;
	}

	const double expected = conf->timesteps+rank;
	double *array = conf->array;

	for (uint64_t e = 0; e < size; ++e) {
		if (array[e] != expected) {
			fprintf(stdout, "Error: Wrong result at rank %d, position %ld, value %f, expected %f\n", rank, e, array[e], expected);
			exit(1);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

double getTime(void)
{
	struct timespec tv;
	clock_gettime(CLOCK_MONOTONIC, &tv);
	return tv.tv_sec+1e-9*tv.tv_nsec;
}
