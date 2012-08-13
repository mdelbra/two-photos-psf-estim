
/**
 * @file psf_estim_main.c
 * @brief main for psf estimation algorithm execution.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "io_pgm.h"



int orsa_homography_sift(float* img1, int nx1, int ny1,
									float* img2, int nx2, int ny2,
									double precision, float *hom,
									float **matchs1, float **matchs2, int* nm);

/**
 * @brief main function call
 */
int main(int argc, char *argv[])
{
	int nx_zin,ny_zin, nx_zout, ny_zout;
	
	float *img_zin, *img_zout;
	
	int nm;
	float *m1, *m2;
	
	float H[9];
	
	
	/* read the PGM image into data */
	if (NULL == (img_zout = read_pgm_float(argv[1], &nx_zout, &ny_zout)))
	{
		fprintf(stderr, "the image could not be properly read\n");
		return EXIT_FAILURE;
	}
	
	
	/* read the PGM image pattern*/
	if (NULL ==
		(img_zin = read_pgm_float(argv[2], &nx_zin, &ny_zin)))
	{
		fprintf(stderr, "the pattern image could not be properly read\n");
		return EXIT_FAILURE;
	}
	
	
		
	/* Call psf_estimation */
	orsa_homography_sift(img_zout, nx_zout, ny_zout,
						 img_zin,  nx_zin, ny_zin,
						 5, H, &m1, &m2, &nm );
	
	
	
	free(img_zin);
	free(img_zout);
	
	return EXIT_SUCCESS;
}

