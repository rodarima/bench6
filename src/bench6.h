/* Copyright (c) 2022 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#ifndef BENCH6_H
#define BENCH6_H

#define UNUSED(x) (void)(x)

double get_time(void);
int get_ncpus();

int bench6_creator(int argc, char *argv[]);
int bench6_sched_get(int argc, char *argv[]);
int bench6_sched_add(int argc, char *argv[]);
int bench6_register_deps(int argc, char *argv[]);

#endif /* BENCH6_H */
