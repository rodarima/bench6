/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National
   Laboratories ( http://www.mantevo.org ). The primary
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) , Paul Crozier
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later
   version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov).

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#ifndef ATOM_H
#define ATOM_H

#include "threadData.hpp"
#include "types.hpp"

#define NDIM 3
#define X 0
#define Y 1
#define Z 2

class Neighbor;
struct Box {
    double xprd, yprd, zprd;

    double xlo[NDIM];
    double xhi[NDIM];
};

class Atom {
public:
    int natoms;
    int nlocal, nghost;
    int nmax;

    double *x;
    double *v;
    double *f;

    int ntypes;
    int *type;

    double *xold;

    ThreadData *threads;
    double virial, mass;

    int comm_size, reverse_size, border_size;

    struct Box box;

    Atom(int ntypes_);
    ~Atom();
    void addatom(double, double, double, double, double, double);
    void pbc();
    void growarray();

    void copy(int, int);

    void pack_comm(int, int *, double *, int *);
    void unpack_comm(int, int, double *);
    void pack_reverse(int, int, double *);
    void unpack_reverse(int, int *, double *);

    int pack_border(int, double *, int *);
    int unpack_border(int, double *);
    int pack_exchange(int, double *);
    int unpack_exchange(int, double *);
    int skip_exchange(double *);

    double *realloc_2d_double_array(double *, int, int, int);
    double *create_2d_double_array(int, int);
    void destroy_2d_double_array(double *);

    int *realloc_1d_int_array(int *, int, int);
    int *create_1d_int_array(int);
    void destroy_1d_int_array(int *);

    void sort(Neighbor &neighbor);

private:
    int *binpos;
    int *bins;
    double *x_copy;
    double *v_copy;
    int *type_copy;
    int copy_size;
};

#endif
