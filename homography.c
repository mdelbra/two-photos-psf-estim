/*----------------------------------------------------------------------------
 
 Point Spread Function Estimation from a random calibration pattern
 
 Copyright 2011 mauricio delbracio (mdelbra@gmail.com)
 
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
 * @file  thin_plates.c
 * @brief library code to estimate/evaluate thin plates splines.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */

/*In this version all the computations are done in 'doubles' in ordert to
  be more accurate in particular when inverting a matrix*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "homography.h"


/**
 * @brief Evaluate a ThinPlate in an array 'np' points 'Pin' and
 *        return the result in the array 'Pout' 
 * @param Pin - Array of input points 
 * @param Pou - Array of output points
 * @param np   - number of points 
 * @return EXIT_SUCCESS
 */
int evaluate_homography(float* H, float *Pin, float *Pout, int np)
{
    int i;
    float xPo, yPo, zPo, xP, yP;
	
	
    for (i = 0; i < np; i++) {
		/*Read the point */
		xP = Pin[2 * i];
		yP = Pin[2 * i + 1];
		
		/*Affine Part */
		xPo = H[0] * xP + H[1] * yP + H[2];
		yPo = H[3] * xP + H[4] * yP + H[5];
		zPo = H[6] * xP + H[7] * yP + H[8];
	
		
		Pout[2 * i] = xPo/zPo;
		Pout[2 * i + 1] = yPo/zPo;
    }
	
    return EXIT_SUCCESS;
}


int inverse_homography(float* H, float* Hi)
{
    double det;
	int i;
	
	det =  0;
	
		 
	/*
	 * | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
	 * | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
	 * | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |
	 *
	 * DET = a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)
	 */
	 
	Hi[0] =  H[8]*H[4] - H[7]*H[5];
	Hi[1] = -H[8]*H[1] + H[7]*H[2];
	Hi[2] =  H[5]*H[1] - H[4]*H[2];
	
	Hi[3] = -H[8]*H[3] + H[6]*H[5];
	Hi[4] =  H[8]*H[0] - H[6]*H[2];
	Hi[5] = -H[5]*H[0] + H[3]*H[2];
	
	Hi[6] =  H[7]*H[3] - H[6]*H[4];
	Hi[7] = -H[7]*H[0] + H[6]*H[1];
	Hi[8] =  H[4]*H[0] - H[4]*H[1];
	
	det = H[0]*Hi[0] - H[3]*Hi[1] + H[6]*Hi[2];
	
	for (i=0;i<9;i++)
		Hi[i] = Hi[i]/det;
	
	/*normalize to have an Homography*/
	/*for (i=0;i<9;i++)
		Hi[i] = Hi[i]/Hi[8];
	*/
	
    return EXIT_SUCCESS;
}

