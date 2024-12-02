//
// This file is part of NBody and is licensed under the terms contained
// in the LICENSE file.
//
// Copyright (C) 2021 Barcelona Supercomputing Center (BSC)
//

#ifndef NBODY_H
#define NBODY_H

#include "common/common.h"

#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>

//#define MIN_PARTICLES (4096 * BLOCK_SIZE / sizeof(particles_block_t))

// Solver structures
typedef struct {
	float *position_x; /* m   */
	float *position_y; /* m   */
	float *position_z; /* m   */
	float *velocity_x; /* m/s */
	float *velocity_y; /* m/s */
	float *velocity_z; /* m/s */
	float *mass;       /* kg  */
	float *weight;
} particles_block_t;

typedef struct {
	float *x; /* x   */
	float *y; /* y   */
	float *z; /* z   */
} forces_block_t;

// Forward declaration
typedef struct nbody_file_t nbody_file_t;
typedef struct nbody_t nbody_t;

// Solver function
void nbody_solve(nbody_t *nbody, const int blocksize, const int num_blocks, const int timesteps, const float time_interval);

// Auxiliary functions
nbody_t nbody_setup(const nbody_conf_t *conf);
void nbody_particle_init(const nbody_conf_t *conf, particles_block_t *part);
void nbody_save_particles(const nbody_t *nbody);
void nbody_free(nbody_t *nbody);
void nbody_check(const nbody_t *nbody);
int nbody_compare_particles(const particles_block_t *local, const particles_block_t *reference, int blocksize, int num_blocks);

#endif // NBODY_H

