#include "streaming.h"

#include <stdint.h>

void computeBlock(const uint64_t size, double result[size], double source[size], double compfactor)
{
	(void) result;
	(void) source;
	uint64_t compsize = 1000*size*compfactor;
	for (uint64_t e = 0; e < compsize; ++e) {
		__asm__("");
	}
}
