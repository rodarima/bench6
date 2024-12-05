#ifndef STREAMING_H
#define STREAMING_H

#include <stdbool.h>
#include <stdint.h>

#define MIN(a,b) (((a)<(b))?(a):(b))

typedef struct {
	uint64_t size;
	uint64_t bs;
	int timesteps;
	double *array;
	double *source;
	double compfactor;
	int offset;
	bool warmup;
	bool verbose;
} StreamingConfiguration;

void initialize(StreamingConfiguration *conf, uint64_t size);
void finalize(StreamingConfiguration *conf);
void readConfiguration(int argc, char **argv, StreamingConfiguration *conf);
void broadcastConfiguration(StreamingConfiguration *configuration);
void checkConfiguration(const StreamingConfiguration *conf, int totalranks);
void printConfiguration(const StreamingConfiguration *conf);
void checkResult(const StreamingConfiguration *conf, uint64_t size, int rank);
double getTime(void);

void solve(StreamingConfiguration *conf, uint64_t size, uint64_t bs, int timesteps, void *extraData);
void computeBlock(const uint64_t size, double block[size], double source[size], double compfactor);

#endif // STREAMING_H
