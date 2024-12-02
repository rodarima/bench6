
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file CG.cpp

 HPCG routine
 */

#include <fstream>

#include <cmath>

#include "hpcg.hpp"

#include "CG.hpp"
#include "mytimer.hpp"
#include "ComputeSPMV.hpp"
#include "ComputeMG.hpp"
#include "ComputeDotProduct.hpp"
#include "ComputeWAXPBY.hpp"

#include "ExchangeHalo.hpp"

#include <mpi.h>


// Use TICK and TOCK to time a code section in MATLAB-like fashion
#define TICK()  t0 = mytimer() //!< record current time in 't0'
#define TOCK(t) t += mytimer() - t0 //!< store time difference in 't' using time in 't0'

/*!
  Routine to compute an approximate solution to Ax = b

  @param[in]    geom The description of the problem's geometry.
  @param[inout] A    The known system matrix
  @param[inout] data The data structure with all necessary CG vectors preallocated
  @param[in]    b    The known right hand side vector
  @param[inout] x    On entry: the initial guess; on exit: the new approximate solution
  @param[in]    max_iter  The maximum number of iterations to perform, even if tolerance is not met.
  @param[in]    tolerance The stopping criterion to assert convergence: if norm of residual is <= to tolerance.
  @param[out]   niters    The number of iterations actually performed.
  @param[out]   normr     The 2-norm of the residual vector after the last iteration.
  @param[out]   normr0    The 2-norm of the residual vector before the first iteration.
  @param[out]   times     The 7-element vector of the timing information accumulated during all of the iterations.
  @param[in]    doPreconditioning The flag to indicate whether the preconditioner should be invoked at each iteration.

  @return Returns zero on success and a non-zero value otherwise.

  @see CG_ref()
*/

int CG(const SparseMatrix & A, CGData & data, const Vector & b, Vector & x,
    const int max_iter, const double tolerance, int & niters, double & normr, double & normr0,
    double * times, int numNeighComms, bool doPreconditioning) {

  double t_begin = mytimer();  // Start timing right away
  niters = 0;
  normr = 0.0;
  double rtz[2];
  rtz[0] = rtz[1] = 0.0;

  double oldrtz[2];
  oldrtz[0] = oldrtz[1] = 0.0;

  double alpha[2];
  alpha[0] = alpha[1] = 0.0;

  double beta[2];
  beta[0] = beta[1] = 0.0;

  double pAp[2];
  pAp[0] = pAp[1] = 0.0;

  double normr_local[2];
  normr_local[0] = 0.0;
  normr_local[1] = 0.0;

  double rtz_local[2];
  rtz_local[0] = 0.0;
  rtz_local[1] = 0.0;

  double pAp_local[2];
  pAp_local[0] = 0.0;
  pAp_local[1] = 0.0;

  double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0;
  Vector & r = data.r; // Residual vector
  Vector & z = data.z; // Preconditioned residual vector
  Vector & p = data.p; // Direction vector (in MPI mode ncol>=nrow)
  Vector & Ap = data.Ap;

  Vector r_buff[2];
  InitializeVector(r_buff[0], r.localLength);
  InitializeVector(r_buff[1], r.localLength);
  CopyVector(r, r_buff[0]);
  CopyVector(r, r_buff[1]);

  Vector z_buff[2];
  InitializeVector(z_buff[0], z.localLength);
  InitializeVector(z_buff[1], z.localLength);
  CopyVector(z, z_buff[0]);
  CopyVector(z, z_buff[1]);

  Vector p_buff[2];
  InitializeVector(p_buff[0], p.localLength);
  InitializeVector(p_buff[1], p.localLength);
  CopyVector(p, p_buff[0]);
  CopyVector(p, p_buff[1]);

  Vector Ap_buff[2];
  InitializeVector(Ap_buff[0], Ap.localLength);
  InitializeVector(Ap_buff[1], Ap.localLength);
  CopyVector(Ap, Ap_buff[0]);
  CopyVector(Ap, Ap_buff[1]);

  Vector x_buff[2];
  InitializeVector(x_buff[0], x.localLength);
  InitializeVector(x_buff[1], x.localLength);
  CopyVector(x, x_buff[0]);
  CopyVector(x, x_buff[1]);

  double *rv_buff[2];
  rv_buff[0] = r_buff[0].values;
  rv_buff[1] = r_buff[1].values;

  double *zv_buff[2];
  zv_buff[0] = z_buff[0].values;
  zv_buff[1] = z_buff[1].values;

  double *pv_buff[2];
  pv_buff[0] = p_buff[0].values;
  pv_buff[1] = p_buff[1].values;

  double *Apv_buff[2];
  Apv_buff[0] = Ap_buff[0].values;
  Apv_buff[1] = Ap_buff[1].values;

  double *xv_buff[2];
  xv_buff[0] = x_buff[0].values;
  xv_buff[1] = x_buff[1].values;

  double *const &bv = b.values; 

  double normr_buff[2];
  normr_buff[0] = normr_buff[1] = 0.0;

  MPI_Comm world_dup[2];
  // Now we have 2 iteration being executed.
  // We can have iteration k-1 in MPI_Allreduce4
  // and iteration k in any other MPI_Allreduce
  MPI_Comm_dup(MPI_COMM_WORLD, &world_dup[0]);
  MPI_Comm_dup(MPI_COMM_WORLD, &world_dup[1]);

  const BlockLayout &blkLayout = A.blkLayout;
  const local_int_t numBlocks = blkLayout.numBlocks;
  const std::vector<BlockInfo> &blkInfo = blkLayout.blkInfo;
  const std::vector< std::vector<local_int_t> > &blkDepList = blkLayout.blkDepList;

  if (!doPreconditioning && A.geom->rank==0) HPCG_fout << "WARNING: PERFORMING UNPRECONDITIONED ITERATIONS" << std::endl;

#ifdef HPCG_DEBUG
  int print_freq = 1;
  if (print_freq>50) print_freq=50;
  if (print_freq<1)  print_freq=1;
#endif

  #pragma omp parallel
  #pragma omp single
  {

  // p is of length ncols, copy x to p for sparse MV operation
  for (int i = 0; i < numBlocks; i++) {
    int start = blkInfo[i].start;
    int end = blkInfo[i].end;
    #pragma omp task \
        default(shared) \
        firstprivate(start, end) \
        depend(in: xv_buff[0][start] ) \
        depend(out: pv_buff[0][start] )
    CopyVector(start, end, x_buff[0], p_buff[0]);
  }

  ExchangeHalo(A,p_buff[0],numNeighComms);

  for (int i = 0; i < numBlocks; i++) {
    int start = blkInfo[i].start;
    int end = blkInfo[i].end;
    #pragma omp task \
        default(shared) \
        shared(blkInfo, blkDepList) \
        firstprivate(start, end) \
        depend(iterator(j=0:blkDepList[i].size()), in: pv_buff[0][blkInfo[blkDepList[i][j]].start] ) \
        depend(out: Apv_buff[0][start] )
    ComputeSPMV(A, start, end, p_buff[0], Ap_buff[0]);  // Ap = A*p
  }

  for (int i = 0; i < numBlocks; i++) {
    int start = blkInfo[i].start;
    int end = blkInfo[i].end;
    #pragma omp task \
        default(shared) \
        firstprivate(start, end) \
        depend(in: bv[start], Apv_buff[0][start] ) \
        depend(out: rv_buff[0][start] )
    ComputeWAXPBY(start, end, 1.0, b, -1.0, Ap_buff[0], r_buff[0], A.isWaxpbyOptimized);   // r = b - Ax (x stored in p)
  }
  #pragma omp taskgroup task_reduction(+: normr_local[0])
  for (int i = 0; i < numBlocks; i++) {
    int start = blkInfo[i].start;
    int end = blkInfo[i].end;
    #pragma omp task \
        default(shared) \
        firstprivate(start, end) \
        depend(in: rv_buff[0][start] ) \
        in_reduction( + : normr_local[0] )
    ComputeDotProduct(start, end, r_buff[0], r_buff[0], normr_local[0], A.isDotProductOptimized);
  }

  #pragma omp task \
      default(shared) \
      depend(inout: *normr_local ) \
      depend(out: normr_buff[0], normr0 )
  {
    MPI_Allreduce(&normr_local[0], &normr_buff[0], 1, MPI_DOUBLE, MPI_SUM,
        world_dup[0]);
    normr_local[0] = 0.0;
    normr_buff[0] = sqrt(normr_buff[0]);

#ifdef HPCG_DEBUG
    if (A.geom->rank==0) HPCG_fout << "Initial Residual = "<< normr_buff[0] << std::endl;
#endif

    // Record initial residual for convergence testing
    normr0 = normr_buff[0];
  }

  #pragma omp taskwait

  // Start iterations

  if (max_iter >= 1 && normr_buff[0]/normr0 > tolerance) {
    int k = 1;
    do {
      double *cur_rtz_local = &rtz_local[k%2];
      double *cur_pAp_local = &pAp_local[k%2];
      double *cur_normr_local = &normr_local[k%2];
      if (doPreconditioning) {
          ComputeMG(A, r_buff[(k-1)%2], z_buff[k%2], numNeighComms); // Apply preconditioner
      } else {
        for (int i = 0; i < numBlocks; i++) {
          int start = blkInfo[i].start;
          int end = blkInfo[i].end;
          #pragma omp task \
              default(shared) \
              firstprivate(k) \
              firstprivate(start, end) \
              depend(in: rv_buff[(k-1)%2][start] ) \
              depend(out: zv_buff[k%2][start] )
          CopyVector (start, end, r_buff[(k-1)%2], z_buff[k%2]); // copy r to z (no preconditioning)
        }
      }

      if (k == 1) {
        for (int i = 0; i < numBlocks; i++) {
           int start = blkInfo[i].start;
           int end = blkInfo[i].end;
           #pragma omp task \
               default(shared) \
               firstprivate(k) \
               firstprivate(start, end) \
               depend(in: zv_buff[k%2][start] ) \
               depend(out: pv_buff[k%2][start] )
           ComputeWAXPBY(start, end, 1.0, z_buff[k%2], 0.0, z_buff[k%2], p_buff[k%2], A.isWaxpbyOptimized);  // Copy Mr to p
         }

         #pragma omp taskgroup task_reduction(+: cur_rtz_local[0])
         for (int i = 0; i < numBlocks; i++) {
           int start = blkInfo[i].start;
           int end = blkInfo[i].end;
           #pragma omp task \
               default(shared) \
               firstprivate(k) \
               firstprivate(start, end) \
               depend(in: rv_buff[(k-1)%2][start], zv_buff[k%2][start] ) \
               in_reduction( + : cur_rtz_local[0] )
           ComputeDotProduct(start, end, r_buff[(k-1)%2], z_buff[k%2], *cur_rtz_local, A.isDotProductOptimized);  // rtz = r'*z
         }

         #pragma omp task \
             default(shared) \
             firstprivate(k, cur_rtz_local) \
             depend(inout: cur_rtz_local[0] ) \
             depend(out: rtz[k%2] )
        {
          MPI_Allreduce(cur_rtz_local, &rtz[k%2], 1, MPI_DOUBLE, MPI_SUM,
              world_dup[k%2]);
          *cur_rtz_local = 0.0;
        }
      } else {
        #pragma omp taskgroup task_reduction(+: cur_rtz_local[0])
        for (int i = 0; i < numBlocks; i++) {
          int start = blkInfo[i].start;
          int end = blkInfo[i].end;
          #pragma omp task \
              default(shared) \
              firstprivate(k) \
              firstprivate(start, end) \
              depend(in: rv_buff[(k-1)%2][start], zv_buff[k%2][start] ) \
              in_reduction( + :  cur_rtz_local[0] )
          ComputeDotProduct(start, end, r_buff[(k-1)%2], z_buff[k%2], *cur_rtz_local, A.isDotProductOptimized);  // rtz = r'*z
        }
        #pragma omp task \
            default(shared) \
            firstprivate(k, cur_rtz_local) \
            depend(inout: cur_rtz_local[0] ) \
            depend(in: rtz[(k-1)%2] ) \
            depend(out: oldrtz[k%2], beta[k%2], rtz[k%2] )
        {
          oldrtz[k%2] = rtz[(k-1)%2];
          MPI_Allreduce(cur_rtz_local, &rtz[k%2], 1, MPI_DOUBLE, MPI_SUM,
              world_dup[k%2]);
          *cur_rtz_local = 0.0;
          beta[k%2] = rtz[k%2]/oldrtz[k%2];
        }
        for (int i = 0; i < numBlocks; i++) {
          int start = blkInfo[i].start;
          int end = blkInfo[i].end;
          #pragma omp task \
              default(shared) \
              firstprivate(k) \
              firstprivate(start, end) \
              depend(in: zv_buff[k%2][start], beta[k%2], pv_buff[(k-1)%2][start]) \
              depend(out: pv_buff[k%2][start] )
          ComputeWAXPBY(start, end, 1.0, z_buff[k%2], beta[k%2], p_buff[(k-1)%2], p_buff[k%2], A.isWaxpbyOptimized); // p = beta*p + z
        }
      }

      ExchangeHalo(A,p_buff[k%2],numNeighComms);

      for (int i = 0; i < numBlocks; i++) {
        int start = blkInfo[i].start;
        int end = blkInfo[i].end;
        #pragma omp task \
            default(shared) \
            firstprivate(k) \
            shared(blkInfo, blkDepList) \
            firstprivate(start, end) \
            depend(iterator(j=0:blkDepList[i].size()), in: pv_buff[k%2][blkInfo[blkDepList[i][j]].start] ) \
            depend(out: Apv_buff[k%2][start] )
        ComputeSPMV(A, start, end, p_buff[k%2], Ap_buff[k%2]);  // Ap = A*p
      }

      #pragma omp taskgroup task_reduction(+: cur_pAp_local[0])
      for (int i = 0; i < numBlocks; i++) {
        int start = blkInfo[i].start;
        int end = blkInfo[i].end;
        #pragma omp task \
            default(shared) \
            firstprivate(k) \
            firstprivate(start, end) \
            depend(in: pv_buff[k%2][start], Apv_buff[k%2][start] ) \
            in_reduction( + : cur_pAp_local[0] )
        ComputeDotProduct(start, end, p_buff[k%2], Ap_buff[k%2], *cur_pAp_local, A.isDotProductOptimized); // alpha = p'*Ap
      }

      #pragma omp task \
          default(shared) \
          firstprivate(k, cur_pAp_local) \
          depend(inout: cur_pAp_local[0] ) \
          depend(in:rtz[k%2]) \
          depend(out: alpha[k%2], pAp[k%2] )
      {
        MPI_Allreduce(cur_pAp_local, &pAp[k%2], 1, MPI_DOUBLE, MPI_SUM,
            world_dup[k%2]);
        *cur_pAp_local = 0.0;
        alpha[k%2] = rtz[k%2]/pAp[k%2];
      }
      for (int i = 0; i < numBlocks; i++) {
        int start = blkInfo[i].start;
        int end = blkInfo[i].end;
        #pragma omp task \
            default(shared) \
            firstprivate(k) \
            firstprivate(start, end) \
            depend(in: alpha[k%2], pv_buff[k%2][start], xv_buff[(k-1)%2][start]) \
            depend(out: xv_buff[k%2][start] )
        ComputeWAXPBY(start, end, 1.0, x_buff[(k-1)%2], alpha[k%2], p_buff[k%2], x_buff[k%2], A.isWaxpbyOptimized); // x = x + alpha*p
        #pragma omp task \
            default(shared) \
            firstprivate(k) \
            firstprivate(start, end) \
            depend(in: alpha[k%2], Apv_buff[k%2][start], rv_buff[(k-1)%2][start] ) \
            depend(out: rv_buff[k%2][start] )
        ComputeWAXPBY(start, end, 1.0, r_buff[(k-1)%2], -alpha[k%2], Ap_buff[k%2], r_buff[k%2], A.isWaxpbyOptimized); // r = r - alpha*Ap
      }

      #pragma omp taskgroup task_reduction(+: cur_normr_local[0])
      for (int i = 0; i < numBlocks; i++) {
        int start = blkInfo[i].start;
        int end = blkInfo[i].end;
        #pragma omp task \
            default(shared) \
            firstprivate(k) \
            firstprivate(start, end) \
            depend(in: rv_buff[k%2][start] ) \
            in_reduction( + : cur_normr_local[0] )
        ComputeDotProduct(start, end, r_buff[k%2], r_buff[k%2], *cur_normr_local, A.isDotProductOptimized);
      }
      #pragma omp task \
          default(shared) \
          firstprivate(k, cur_normr_local) \
          depend(inout: cur_normr_local[0] ) \
          depend(in: normr0 ) \
          depend(out:normr_buff[k%2])
      {
        MPI_Allreduce(cur_normr_local, &normr_buff[k%2], 1, MPI_DOUBLE, MPI_SUM,
            world_dup[k%2]);
        *cur_normr_local = 0.0;
        normr_buff[k%2] = sqrt(normr_buff[k%2]);
  #ifdef HPCG_DEBUG
        if (A.geom->rank==0 && (k%print_freq == 0 || k == max_iter))
          HPCG_fout << "Iteration = "<< k << "   Scaled Residual = "<< normr_buff[k%2]/normr0 << std::endl;
  #endif
        niters = k;
      }
      // TODO: replace task if(0) by taskwait deps, currently unsupported
      // #pragma omp taskwait depend(in:normr_buff[(k-1)%2])
      #pragma omp task if(0) firstprivate(k) depend(in:normr_buff[(k-1)%2])
      {}
      ++k;
    } while (k<=max_iter && normr_buff[(k-2)%2]/normr0 > tolerance);
  }

  }

  die_if(niters == 0, "niters is 0, which is not possible.");
  // Set the result value that satisfies early exit
  --niters;
  normr = normr_buff[niters%2];

  // Store times
  times[1] += t1; // dot-product time
  times[2] += t2; // WAXPBY time
  times[3] += t3; // SPMV time
  times[4] += t4; // AllReduce time
  times[5] += t5; // preconditioner apply time
  times[0] += mytimer() - t_begin;  // Total time. All done...
  return 0;
}
