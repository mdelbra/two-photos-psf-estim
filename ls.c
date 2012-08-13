
/**
 * @file ls.c
 * @brief library code with numerical algorithms for solving least squares
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 * @date Jun 28, 2011 
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ls.h"
#include <cblas.h>
#include <float.h>

/** Buffer Size */
#define BUFFER_SIZE 113337

/** If the  absolute value is less than EPS_ZERO consider it is zero */
#define EPS_ZERO 0.0001 


/*Wrapper from LAPACK, CBLAS*/


/*SUBROUTINE SGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
*				  $                   WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
REAL               RCOND
*     ..
*     .. Array Arguments ..
INTEGER            JPVT( * )
REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SGELSY computes the minimum-norm solution to a real linear least
*  squares problem:
*      minimize || A * X - B ||
*  using a complete orthogonal factorization of A.  A is an M-by-N
*  matrix which may be rank-deficient.
*
*  Several right hand side vectors b and solution vectors x can be
*  handled in a single call; they are stored as the columns of the
*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*  matrix X.
*
*  The routine first computes a QR factorization with column pivoting:
*      A * P = Q * [ R11 R12 ]
*                  [  0  R22 ]
*  with R11 defined as the largest leading submatrix whose estimated
*  condition number is less than 1/RCOND.  The order of R11, RANK,
*  is the effective rank of A.
*
*  Then, R22 is considered to be negligible, and R12 is annihilated
*  by orthogonal transformations from the right, arriving at the
*  complete orthogonal factorization:
*     A * P = Q * [ T11 0 ] * Z
*                 [  0  0 ]
*  The minimum-norm solution is then
*     X = P * Z**T [ inv(T11)*Q1**T*B ]
*                  [        0         ]
*  where Q1 consists of the first RANK columns of Q.
*
*  This routine is basically identical to the original xGELSX except
*  three differences:
*    o The call to the subroutine xGEQPF has been substituted by the
*      the call to the subroutine xGEQP3. This subroutine is a Blas-3
*      version of the QR factorization with column pivoting.
*    o Matrix B (the right hand side) is updated with Blas-3.
*    o The permutation of matrix B (the right hand side) is faster and
*      more simple.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of
*          columns of matrices B and X. NRHS >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A has been overwritten by details of its
*          complete orthogonal factorization.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the M-by-NRHS right hand side matrix B.
*          On exit, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,M,N).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
*          to the front of AP, otherwise column i is a free column.
*          On exit, if JPVT(i) = k, then the i-th column of AP
*          was the k-th column of A.
*
*  RCOND   (input) REAL
*          RCOND is used to determine the effective rank of A, which
*          is defined as the order of the largest leading triangular
*          submatrix R11 in the QR factorization with pivoting of A,
*          whose estimated condition number < 1/RCOND.
*
*  RANK    (output) INTEGER
*          The effective rank of A, i.e., the order of the submatrix
*          R11.  This is the same as the order of the submatrix T11
*          in the complete orthogonal factorization of A.
*
*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          The unblocked strategy requires that:
*             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),
*          where MN = min( M, N ).
*          The block algorithm requires that:
*             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
*          where NB is an upper bound on the blocksize returned
*          by ILAENV for the routines SGEQP3, STZRZF, STZRQF, SORMQR,
*          and SORMRZ.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: If INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
*
*  =====================================================================
*/



static long sgelsy(int m, int n, int nrhs,
				   float *a, int lda,  float *b, int ldb, 
				   int *jpvt, float rcond, int *rank,
				   float *work, int lwork, int *iwork)
{
    extern void sgelsy_(int *m, int *n, int *nrhs,
						float *a, int *lda,  float *b, int *ldb, 
						int *jpvt, float *rcond, int *rank,
						float *work, int *lwork, int *iwork, int *info);
    int info;
	
	sgelsy_(&m, &n, &nrhs,
			a, &lda,  b, &ldb, 
			jpvt, &rcond, rank,
			work, &lwork, iwork, &info);

    return info;
}




/* SUBROUTINE SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
*				 $                  INFO )
*
*  -- LAPACK driver routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
CHARACTER          TRANS
INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SGELS solves overdetermined or underdetermined real linear systems
*  involving an M-by-N matrix A, or its transpose, using a QR or LQ
*  factorization of A.  It is assumed that A has full rank.
*
*  The following options are provided: 
*
*  1. If TRANS = 'N' and m >= n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || B - A*X ||.
*
*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*     an underdetermined system A * X = B.
*
*  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
*     an undetermined system A**T * X = B.
*
*  4. If TRANS = 'T' and m < n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || B - A**T * X ||.
*
*  Several right hand side vectors b and solution vectors x can be 
*  handled in a single call; they are stored as the columns of the
*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution 
*  matrix X.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          = 'N': the linear system involves A;
*          = 'T': the linear system involves A**T. 
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of
*          columns of the matrices B and X. NRHS >=0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*            if M >= N, A is overwritten by details of its QR
*                       factorization as returned by SGEQRF;
*            if M <  N, A is overwritten by details of its LQ
*                       factorization as returned by SGELQF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the matrix B of right hand side vectors, stored
*          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
*          if TRANS = 'T'.  
*          On exit, if INFO = 0, B is overwritten by the solution
*          vectors, stored columnwise:
*          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
*          squares solution vectors; the residual sum of squares for the
*          solution in each column is given by the sum of squares of
*          elements N+1 to M in that column;
*          if TRANS = 'N' and m < n, rows 1 to N of B contain the
*          minimum norm solution vectors;
*          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
*          minimum norm solution vectors;
*          if TRANS = 'T' and m < n, rows 1 to M of B contain the
*          least squares solution vectors; the residual sum of squares
*          for the solution in each column is given by the sum of
*          squares of elements M+1 to N in that column.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= MAX(1,M,N).
*
*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          LWORK >= max( 1, MN + max( MN, NRHS ) ).
*          For optimal performance,
*          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
*          where MN = min(M,N) and NB is the optimum block size.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO =  i, the i-th diagonal element of the
*                triangular factor of A is zero, so that A does not have
*                full rank; the least squares solution could not be
*                computed.
*
*  =====================================================================
*/


static long sgels(char *trans, int m, int n, int nrhs,
				   float *a, int lda,  float *b, int ldb,
				   float *work, int lwork, int *iwork)
{
    extern void sgels_(const char *trans, const int *m, const int *n, const int *nrhs,
						float *a, const int *lda, float *b,
						const int *ldb, float *work, int *lwork, int *iwork,
						int *info);
    int info;
    sgels_(trans, &m, &n, &nrhs, a, &lda, b, &ldb,
			work, &lwork, iwork, &info);
    return info;
}


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
static long sgelsd(int m, int n, int nrhs,
				   float *a, int lda,  float *b, int ldb,
				   float *s, float rcond, int *rank,
				   float *work, int lwork, int *iwork)
{
    extern void sgelsd_(const int *m, const int *n, const int *nrhs,
						float *a, const int *lda, float *b,
						const int *ldb, float *s, const float *rcond,
						int *rank, float *work, int *lwork, int *iwork,
						int *info);
    int info;
    sgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank,
			work, &lwork, iwork, &info);
    return info;
}

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
 * @brief Solve Least Squares problem x such that Ax = b. 
 * @param A  - Array cointaining matrix 'A' elements (column major) 
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A' 
 * @param n  - number of columns of 'A' 
 * @return Array 'x' with the solution
 */
float *solve_ls(float *A, float *b, int m, int n)
{
    float *work;
    float *s;
	float *x;
	int i ;
	
    float lwork;
    int iwork[BUFFER_SIZE];
    int rank, info;
	
    /*Least Squares */
    s = (float *) malloc(n * sizeof(float));
	
    /*Do a query to know the optimum BufferSize */
    info =
	sgelsd(m, n, 1, A, m, b, m, s, EPS_ZERO, &rank, &lwork, -1, iwork);
	
    work = (float *) malloc(lwork * sizeof(float));
    info =
	sgelsd(m, n, 1, A, m, b, m, s, EPS_ZERO, &rank, work, (int) lwork,
	       iwork);
	
    free((void *) s);
    free((void *) work);
	
	x  = (float *) malloc(n * sizeof(float));
	for(i=0;i<n;i++)
		x[i] = b[i];
	
    return x;
	
}

/**
 * @brief Solve Least Squares problem x such that Ax = b. 
 * @param A  - Array cointaining matrix 'A' elements (column major) 
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A' 
 * @param n  - number of columns of 'A' 
 * @return Array 'x' with the solution
 */
float *solve_lss(float *A, float *b, int m, int n)
{
    float *work;
  	float *x;
	int i ;
	
    float lwork;
    int iwork[BUFFER_SIZE];
    int info;
	
   
    /*Do a query to know the optimum BufferSize */
    info =
	sgels("N", m, n, 1, A, m, b, m,  &lwork, -1, iwork);
	
	
    work = (float *) malloc(lwork * sizeof(float));
    info =
	sgels("N", m, n, 1, A, m, b, m, work, (int) lwork,
	       iwork);
	
    free((void *) work);
	
	x  = (float *) malloc(n * sizeof(float));
	for(i=0;i<n;i++)
		x[i] = b[i];
	
    return x;
	
}


/**
* @brief Solve Least Squares problem x such that Ax = b. 
* @param A  - Array cointaining matrix 'A' elements (column major) 
* @param b  - Array of observed values 'b'
* @param m  - number of rows of 'A' 
* @param n  - number of columns of 'A' 
* @return Array 'x' with the solution
*/
float *solve_lsd(float *Af, float *bf, int m, int n)
{
    double *work;
    double *s;
	float *x;
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
	
    free((void *) s);
    free((void *) work);
	
	x  = (float *) malloc(n * sizeof(float));
	for(i=0;i<n;i++)
		x[i] = (float) b[i];
	
    return x;
	
}



/**
 * @brief Solve Least Squares problem x such that Ax = b and then
 *        thresholds the solution with 'th' such that x<=th is 0. 
 * @param A  - Array cointaining matrix 'A' elements (column major) 
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A' 
 * @param n  - number of columns of 'A' 
 * @param th - final threshold. 
 * @return Array 'x' with the solution.
 */
float *solve_ls_th(float *A, float *b, int m, int n, float th)
{
	
    float *x;
    int i;
	
    x = solve_ls(A, b, m, n);
	
    for (i = 0; i < n; i++)
		x[i] = (x[i] > th) ? x[i] : 0;
	
    return x;
	
}




/**
 * @brief Solve Least Squares problem x such that Ax = b. 
 * @param A  - Array cointaining matrix 'A' elements (column major) 
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A' 
 * @param n  - number of columns of 'A' 
 * @return Array 'x' with the solution
 */
float *solve_ls_pivot(float *A, float *b, int m, int n)
{
    float *work;
    int *jpvt;
	float *x;
	int i ;
	
    float lwork;
    int iwork[BUFFER_SIZE];
    int rank, info;
	
    /*Least Squares */
    jpvt = (int *) malloc(n * sizeof(int));
	
    /*Do a query to know the optimum BufferSize */
    info = sgelsy(m, n, 1, A, m, b, m, jpvt, EPS_ZERO, &rank, &lwork, -1, iwork);
	
    work = (float *) malloc(lwork * sizeof(float));
    
	info = sgelsy(m, n, 1, A, m, b, m, jpvt, EPS_ZERO, &rank, work, (int) lwork, iwork);
	
    free((void *) jpvt);
    free((void *) work);
	
	x  = (float *) malloc(n * sizeof(float));
	for(i=0;i<n;i++)
		x[i] = b[i];
	
    return x;
	
}



