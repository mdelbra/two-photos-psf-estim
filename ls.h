
/**
 * @file ls.h
 * @brief library header numerical algorithms for solving least squares
 *     
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 * @date Jun 28, 2011 
 */


#ifndef LS_H_
#define LS_H_


float *solve_ls_th(float *A, float *b, int m, int n, float th);
float *solve_ls(float *A, float *b, int m, int n);
float *solve_lsd(float *A, float *b, int m, int n);
float *solve_lss(float *A, float *b, int m, int n);
float *solve_ls_pivot(float *A, float *b, int m, int n);

#endif				/* NNLS_H_ */
