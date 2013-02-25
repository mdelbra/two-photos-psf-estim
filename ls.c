/*----------------------------------------------------------------------------
 
 "Recovering the Subpixel PSF from Two Photographs at Different Distances"
 
 Copyright 2013 mauricio delbracio (mdelbra@gmail.com)
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU Affero General Public License for more details.
 
 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 
 ----------------------------------------------------------------------------*/


/**
 * @file ls.c
 * @brief library code with numerical algorithms for solving least squares
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 * @date Nov 24, 2011
 */


/*
 Double precision functions were added, there were some problems
 due to the float precision, now, all *d functions convert the
 input A and b matrix to double matrix and do all the computations with doubles
 finally the result is truncated into a float to be compatible with the other
 part of the program
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ls.h"
#include <cblas.h>
#include <float.h>

/** Buffer Size */
#define BUFFER_SIZE 113337

/** If the absolute value is less than EPS_ZERO consider it is zero */
#define EPS_ZERO 1e-6


/*Wrapper functions to use LAPACK*/
/*
 *SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
 $                   WORK, LWORK, IWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.2.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     June 2010
 *
 *     .. Scalar Arguments ..
 INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
 DOUBLE PRECISION   RCOND
 *     ..
 *     .. Array Arguments ..
 INTEGER            IWORK( * )
 DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DGELSD computes the minimum-norm solution to a real linear least
 *  squares problem:
 *      minimize 2-norm(| b - A*x |)
 *  using the singular value decomposition (SVD) of A. A is an M-by-N
 *  matrix which may be rank-deficient.
 *
 *  Several right hand side vectors b and solution vectors x can be
 *  handled in a single call; they are stored as the columns of the
 *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
 *  matrix X.
 *
 *  The problem is solved in three steps:
 *  (1) Reduce the coefficient matrix A to bidiagonal form with
 *      Householder transformations, reducing the original problem
 *      into a "bidiagonal least squares problem" (BLS)
 *  (2) Solve the BLS using a divide and conquer approach.
 *  (3) Apply back all the Householder tranformations to solve
 *      the original least squares problem.
 *
 *  The effective rank of A is determined by treating as zero those
 *  singular values which are less than RCOND times the largest singular
 *  value.
 *
 *  The divide and conquer algorithm makes very mild assumptions about
 *  floating point arithmetic. It will work on machines with a guard
 *  digit in add/subtract, or on those binary machines without guard
 *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
 *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
 *  without guard digits, but we know of none.
 *
 *  Arguments
 *  =========
 *
 *  M       (input) INTEGER
 *          The number of rows of A. M >= 0.
 *
 *  N       (input) INTEGER
 *          The number of columns of A. N >= 0.
 *
 *  NRHS    (input) INTEGER
 *          The number of right hand sides, i.e., the number of columns
 *          of the matrices B and X. NRHS >= 0.
 *
 *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
 *          On entry, the M-by-N matrix A.
 *          On exit, A has been destroyed.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *
 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
 *          On entry, the M-by-NRHS right hand side matrix B.
 *          On exit, B is overwritten by the N-by-NRHS solution
 *          matrix X.  If m >= n and RANK = n, the residual
 *          sum-of-squares for the solution in the i-th column is given
 *          by the sum of squares of elements n+1:m in that column.
 *
 *  LDB     (input) INTEGER
 *          The leading dimension of the array B. LDB >= max(1,max(M,N)).
 *
 *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
 *          The singular values of A in decreasing order.
 *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
 *
 *  RCOND   (input) DOUBLE PRECISION
 *          RCOND is used to determine the effective rank of A.
 *          Singular values S(i) <= RCOND*S(1) are treated as zero.
 *          If RCOND < 0, machine precision is used instead.
 *
 *  RANK    (output) INTEGER
 *          The effective rank of A, i.e., the number of singular values
 *          which are greater than RCOND*S(1).
 *
 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK. LWORK must be at least 1.
 *          The exact minimum amount of workspace needed depends on M,
 *          N and NRHS. As long as LWORK is at least
 *              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
 *          if M is greater than or equal to N or
 *              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
 *          if M is less than N, the code will execute correctly.
 *          SMLSIZ is returned by ILAENV and is equal to the maximum
 *          size of the subproblems at the bottom of the computation
 *          tree (usually about 25), and
 *             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
 *          For good performance, LWORK should generally be larger.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued by XERBLA.
 *
 *  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))
 *          LIWORK >= max(1, 3 * MINMN * NLVL + 11 * MINMN),
 *          where MINMN = MIN( M,N ).
 *          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          > 0:  the algorithm for computing the SVD failed to converge;
 *                if INFO = i, i off-diagonal elements of an intermediate
 *                bidiagonal form did not converge to zero.
 *
 *  Further Details
 *  ===============
 *
 *  Based on contributions by
 *     Ming Gu and Ren-Cang Li, Computer Science Division, University of
 *       California at Berkeley, USA
 *     Osni Marques, LBNL/NERSC, USA
 *
 *  =====================================================================
 */

static long dgelsd(int m, int n, int nrhs,
                   double *a, int lda,  double *b, int ldb,
                   double *s, double rcond, int *rank,
                   double *work, int lwork, int *iwork)
{
    extern void dgelsd_(const int *m, const int *n, const int *nrhs,
                double *a, const int *lda, double *b,
                const int *ldb, double *s, const double *rcond,
                int *rank, double *work, int *lwork, int *iwork,
            int *info);
    int info;
    dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank,
            work, &lwork, iwork, &info);
    return info;
}





/**
 * @brief Solve Least Squares problem (double precision) x such that Ax = b.
 * @param A  - Array cointaining matrix 'A' elements (column major)
 * @param b  - Array of observed values 'b'
 * @param x  - Array with the solution 'x' 
 * @param n  - number of columns of 'A'
 * @param m  - number of rows of 'A'
 */
void solve_lsd(float *Af, float *bf, float *x, int n, int m)
{
    double *work;
    double *s;
    int i ;

    double lwork;
    int iwork[BUFFER_SIZE];
    int rank, info;

    double *A, *b;

    A = (double *) malloc(m*n*sizeof(double));
    b = (double *) malloc(m*sizeof(double));


    for(i=0;i<m*n;i++)
        A[i] = (double) Af[i];

    for(i=0;i<m;i++)
        b[i] = (double) bf[i];


    /*Least Squares */
    s = (double *) malloc(n * sizeof(double));

    /*Do a query to know the optimum BufferSize */
    info =
    dgelsd(m, n, 1, A, m, b, m, s, EPS_ZERO, &rank, &lwork, -1, iwork);

    work = (double *) malloc(lwork * sizeof(double));
    info =
    dgelsd(m, n, 1, A, m, b, m, s, EPS_ZERO, &rank, work, (int) lwork,
           iwork);

    for(i=0;i<n;i++)   x[i] = (float) b[i];

    free((void *) s);
    free((void *) work);
    free((void *) A);
    free((void *) b);
    
}


