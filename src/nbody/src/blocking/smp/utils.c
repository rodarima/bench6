//
// This file is part of NBody and is licensed under the terms contained
// in the LICENSE file.
//
// Copyright (C) 2021 Barcelona Supercomputing Center (BSC)
//

#include "blocking/smp/nbody.h"

#undef NDEBUG
#include <assert.h>
#include <ctype.h>
#include <fcntl.h>
#include <getopt.h>
#include <ieee754.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

/*
static void nbody_generate_particles(const nbody_conf_t *conf, const nbody_file_t *file)
{
	char fname[1024];
	sprintf(fname, "%s.in", file->name);
	
	if (!access(fname, F_OK)) {
		return;
	}
	
	struct stat st = {0};
	if (stat("data", &st) == -1) {
		mkdir("data", 0755);
	}
	
	const int fd = open(fname, O_RDWR|O_CREAT|O_TRUNC, S_IRUSR|S_IRGRP|S_IROTH);
	assert(fd >= 0);
	
	const int size = file->size;
	assert(size % PAGE_SIZE == 0);
	
	int err = ftruncate(fd, size);
	assert(!err);
	
	particles_block_t * const particles = mmap(NULL, size, PROT_WRITE|PROT_READ, MAP_SHARED, fd, 0);
	
	for(int i = 0; i < conf->num_blocks; i++) {
		nbody_particle_init(conf, particles+i);
	}
	
	err = munmap(particles, size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}
*/

void nbody_check(const nbody_t *nbody)
{
	char fname[1024];
	sprintf(fname, "%s.ref", nbody->file.name);
	if (access(fname, F_OK) != 0) {
		fprintf(stderr, "Warning: %s file does not exist. Skipping the check...\n", fname);
		return;
	}
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	particles_block_t *reference = mmap(NULL, nbody->file.size, PROT_READ, MAP_SHARED, fd, 0);
	assert(reference != MAP_FAILED);
	
	if (nbody_compare_particles(nbody->particles, reference, nbody->blocksize, nbody->num_blocks)) {
		printf("Result validation: OK\n");
	} else {
		printf("Result validation: ERROR\n");
	}
	
	int err = munmap(reference, nbody->file.size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

static nbody_file_t nbody_setup_file(const nbody_conf_t *conf)
{
	nbody_file_t file;
	file.size = conf->num_blocks * sizeof(particles_block_t);
	
	sprintf(file.name, "%s-%s-%d-%d-%d", conf->name, TOSTRING(BIGO), conf->num_blocks * conf->blocksize, conf->blocksize, conf->timesteps);
	return file;
}
/*
static particles_block_t *nbody_load_particles(const nbody_conf_t *conf, const nbody_file_t *file)
{
	(void) conf;

	char fname[1024];
	sprintf(fname, "%s.in", file->name);
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	void * const ptr = mmap(NULL, file->size, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, 0);
	assert(ptr != MAP_FAILED);
	
	int err = close(fd);
	assert(!err);
	
	return ptr;
}
*/

static void alloc_particles(const nbody_conf_t *conf, nbody_t *nbody)
{
	size_t bs = (size_t) conf->blocksize;
	size_t n = (size_t) conf->num_blocks;
	size_t block_floats = bs * (3 + 3 + 2);
	size_t block_bytes = block_floats * sizeof(float);
	size_t alloc_bytes = n * block_bytes;

	nbody->particles_map = nbody_alloc(alloc_bytes);
	nbody->particles = nbody_alloc(n * sizeof(particles_block_t));

	float *m = nbody->particles_map;

	/* Setup the particles_block_t pointers, so they point to the big chunk
	 * in particle_map:
	 *
	 *     [      particles[0]     ] [    particles[1]    ] ...                
	 *     [ px ][ py ][ pz ] [... ]
	 *       |     \_________________
	 *       |                       \
	 *       v                       v
	 *     [ px0 ] [ px1 ] [ ... ] [ py0 ] [ ... ] ...    
	 *     [             particles_map                    ]
	 */

	for (int i = 0; i < conf->num_blocks; i++) {
		particles_block_t *p = &nbody->particles[i];
		float *b = &m[block_floats * i];
		p->position_x = &b[bs * 0];
		p->position_y = &b[bs * 1];
		p->position_z = &b[bs * 2];
		p->velocity_x = &b[bs * 3];
		p->velocity_y = &b[bs * 4];
		p->velocity_z = &b[bs * 5];
		p->mass       = &b[bs * 6];
		p->weight     = &b[bs * 7];
	}
}

static void alloc_forces(const nbody_conf_t *conf, nbody_t *nbody)
{
	size_t bs = (size_t) conf->blocksize;
	size_t n = (size_t) conf->num_blocks;
	size_t block_floats = bs * 3;
	size_t block_bytes = block_floats * sizeof(float);
	size_t alloc_bytes = n * block_bytes;

	nbody->forces_map = nbody_alloc(alloc_bytes);
	nbody->forces = nbody_alloc(n * sizeof(forces_block_t));

	float *m = nbody->forces_map;

	/* Setup the forces_block_t pointers, so they point to the big chunk
	 * in forces_map:
	 *
	 *     [  forces[0]  ] [  forces[1]  ] ...                
	 *     [ x ][ y ][ z ]
	 *       |     \________________
	 *       |                      \
	 *       v                      v
	 *     [ x0 ] [ x1 ] [ ... ] [ y0 ] [ y1 ]  [ ... ]
	 *     [             forces_map                   ]
	 */

	for (int i = 0; i < conf->num_blocks; i++) {
		forces_block_t *f = &nbody->forces[i];
		float *b = &m[block_floats * i];
		f->x = &b[bs * 0];
		f->y = &b[bs * 1];
		f->z = &b[bs * 2];
	}
}

nbody_t nbody_setup(const nbody_conf_t *conf)
{
	nbody_t nbody;
	nbody.timesteps = conf->timesteps;
	nbody.blocksize = conf->blocksize;
	nbody.num_blocks = conf->num_blocks;
	
	nbody_file_t file = nbody_setup_file(conf);
	nbody.file = file;
	
	if (conf->force_generation) {
		alloc_particles(conf, &nbody);
		
		for (int i = 0; i < conf->num_blocks; i++)
			nbody_particle_init(conf, nbody.particles+i);
	} else {
		fprintf(stderr, "not implemented\n");
		exit(1);
		/*
		nbody_generate_particles(conf, &file);
		nbody.particles = nbody_load_particles(conf, &file);
		assert(nbody.particles != NULL);
		*/
	}

	alloc_forces(conf, &nbody);
	
	return nbody;
}

void nbody_save_particles(const nbody_t *nbody)
{
	char fname[1024];
	sprintf(fname, "%s.out", nbody->file.name);
	
	const int fd = open(fname, O_RDWR|O_CREAT|O_TRUNC, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH);
	assert(fd >= 0);
	
	const int size = nbody->file.size;
	assert(size % PAGE_SIZE == 0);
	
	int err = ftruncate(fd, size);
	assert(!err);
	
	particles_block_t * const particles = mmap(NULL, size, PROT_WRITE|PROT_READ, MAP_SHARED, fd, 0);
	assert(particles != MAP_FAILED);
	
	memcpy(particles, nbody->particles, size);
	
	err = munmap(particles, size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_free(nbody_t *nbody)
{
	int err = munmap(nbody->particles, nbody->num_blocks * sizeof(particles_block_t));
	err |= munmap(nbody->forces, nbody->num_blocks * sizeof(forces_block_t));
	assert(!err);
}

