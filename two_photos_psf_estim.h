/*
 * two_photos_psf_estim.h
 *
 *  Created on: Oct 22, 2010
 *      Author: mdelbra
 */

#ifndef TWO_PHOTOS_PSF_ESTIM_H_
#define TWO_PHOTOS_PSF_ESTIM_H_


void write_ascii_matrix(float *M, int ncol, int nrow, char *name);

void two_photos_psf_estim(float *img, int nx, int ny,
		 float *pattern, int pat_nx, int pat_ny,
		 int s, int psf_nrow, int psf_ncol, int solver, float *h,
					 float *k, char* detected_pgm);



#endif				/* TWO_PHOTOS_PSF_ESTIM_H_ */
