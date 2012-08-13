/**
 * @file thin_plates.h
 * @brief library header to estimate/evaluate thin plates splines.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#ifndef HOMOGRAPHY_H
#define HOMOGRAPHY

int evaluate_homography(float* H, float *Pin, float *Pout, int np);

int inverse_homography(float* H, float* Hi);


#endif				/* THIN_PLATES_H_ */
