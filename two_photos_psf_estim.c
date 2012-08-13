
/**
 * @file two_photos_psf_estim.c
 * @brief library code to two-photos psf estimation.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


/** @mainpage psf estimation code documentation
 * 
 * The following is an implementation of 
 * "Recovering the Subpixel PSF from Two Photographs at Different Distances"
 * 
 * A More detail description can be found on the portal IPOL www.ipol.im where
 * there is more information, including this code and an online demo version:
 * 
 * http://www.ipol.im/pub/algo/
 * 
 * 
 * HISTORY:
 * - version 0.4 - sep 2012: Second BETA Release Ansi C Language
 *	version.
 * 
 * @author mauricio delbracio (mdelbra@gmail.com)
 * @date feb 2012
 */


/*updated auxiliary images. Now the auxiliary images are in PGM format and are
 rescaled to imgC. Difference image is rescaled to 0.15 of the min and max of imgC*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "image.h"
#include "homography.h"
#include "ls.h"
#include "io_pgm.h"


int orsa_homography_sift(float* img1, int nx1, int ny1,
						float* img2, int nx2, int ny2,
						double precision, float *hom,
						 float **matchs1, float **matchs2, int* nm, char* reverse);


/*
 #DEFINE SAVE_INTERMEDIATE_IMGS
 */

/** If the  value of the PSF is less than POSITIVE_TOL consider it is zero */
#define POSITIVE_TOL 0.000001
#define BIG_NUMBER 1000000
/** @brief Error/Exit print a message and exit.
 *  @param msg 
 */
void error (char *msg)
{
	fprintf (stderr, "psf_estim Error: %s\n", msg);
	exit (EXIT_FAILURE);
}


/** 
 *  @brief Write an ImageFloat to an ASCII file
 *  @param image - image to write
 *  @param name  - filename
 */
void write_ascii_imageFloat (ImageFloat image, char *name)
{
	FILE *f;
	int x, y, n;
	
	
	/* open file */
	f = fopen (name, "w");
	if (f == NULL)
    {
		error ("Can't open output file.");
    }
	
	/* write header */
	
	/* write data */
	for (y = 0; y < image->nrow; y++)
    {
		for (x = 0; x < image->ncol; x++, n++){
			
			fprintf (f, "%f ", image->val[x + y * image->ncol]);
		}
		
		fprintf (f,"\n");
    }
	
	/* close file */
	fclose (f);
}

/* homemade round function: float to int*/
static int roundfi(float x) {
	if ((x <= INT_MIN-0.5) || (x >= INT_MAX+0.5))
		error("roundfi() Float to int conversion out of range");
	if (x >= 0)
		return (int) (x+0.5);
	return (int) (x-0.5);
}



/** 
 *  @brief Write an Array of floats (Matrix)  to an ASCII file
 *  @param M - array of floats to write
 *  @param ncol - number of columns of M
 *  @param nrow - number of rows of M
 *  @param name  - filename
 */
void write_ascii_matrix(float *M, int ncol, int nrow, char *name)
{
	FILE *f;
	int x, y, n;
	
	/* open file */
	f = fopen (name, "w");
	if (f == NULL)
		error ("Can't open output file.");
	
	/* write header */
	
	/* write data */
	for (y = 0; y < nrow; y++)
    {
		for (x = 0; x < ncol; x++, n++)
			fprintf (f, "%f ", M[y + x * nrow]);
		
		fprintf (f, "\n");
    }
	
	
	/* close file */
	fclose (f);
}



/** 
 *  @brief Low-Pass filter Image (with DCT Transform)
 */
ImageFloat lpf_image_dct (ImageFloat in, int fcx, int fcy)
{
	ImageFloat data_dct, out;
	int i, j;
	int nx = in->ncol;
	int ny = in->nrow;
	float k;
	
	data_dct = compute_dct_image (in);
	
	for (i = fcy; i < ny; i++)
		for (j = 0; j < nx; j++)
			data_dct->val[i * nx + j] = 0;
	
	for (i = 0; i < ny; i++)
		for (j = fcx; j < nx; j++)
			data_dct->val[i * nx + j] = 0;
	
	
	out = compute_idct_image (data_dct);
	
	/*NormaliHRg image because DCT introduces
	 * a constant factor 4*nx*ny
	 */
	k = 4 * nx * ny;
	for (i = 0; i < nx * ny; i++)
		out->val[i] = out->val[i] / k;
	
	
	/*
	 * cleanup
	 */
	free_imageFloat (data_dct);
	
	return out;
	
}

/** 
 *  @brief Generate the Linear System Ax = b where x is the kernel to find.
 */
int make_Ab (ImageFloat imgC, ImageFloat imgW, ImageFloat imgMask,
			 int q, int p, int s, float **A, float *b[], int *ncol, int *nrow)
{
	/*Generate a SsU from U and s. */
	int max_pq, u, v, i, j, mkMs, nkMs, mc, nc, r;
	
	float *kerode;
	ImageFloat mask;
	/* generate SsU matrix */
	/* p,q kernel size */
	
	mkMs = (imgC->nrow + p - 2) / s + 1;
	nkMs = (imgC->ncol + q - 2) / s + 1;
	
	*A = (float *) calloc (mkMs * nkMs * p * q, sizeof (float));	/*initialize with zeros */
	*b = (float *) calloc (mkMs * nkMs, sizeof (float));	/*initialize with zeros */
	
	/*CHECK THE VALUE OF r!!!!*/
	/*Erode the Mask in order take into account the boundary
	 * problems : square element of side r */
	max_pq = (p > q) ? p : q;
	r = (max_pq-1) / s; /* CHECK!!! d +2 because there were strange artifacts*/
	kerode = (float *) malloc ((2 * r + 1) * sizeof (float));	/*initialize with zeros */
	
	
	for (i = 0; i < 2 * r + 1; i++)
		kerode[i] = 1;
	
	mask = convol_sep2 (imgMask, kerode, 2 * r + 1, kerode, 2 * r + 1);
	
	for (i = 0; i < imgMask->nrow; i++)
		for (j = 0; j < imgMask->ncol; j++)
			imgMask->val[i * imgMask->ncol + j] =
			(mask->val[(i + r) * mask->ncol + j + r] >= 4 * r * r + 1) ? 1 : 0;
	
	
	/*Filling in (*A) */
	for (i = (p - 1) / 2; i < (imgC->nrow) + (p - 1) / 2; i = i + s)
    {
		for (j = (q - 1) / 2; j < (imgC->ncol) + (q - 1) / 2; j = j + s)
		{
			if (imgMask->val[(j - (q - 1) / 2) / s +
							 imgMask->ncol * (i - (p - 1) / 2) / s])
			{
				for (u = 0; u < p; u++)
				{
					for (v = 0; v < q; v++)
					{
						if ((i - u >= 0) && (j - v >= 0)
							&& (i - u < imgC->nrow) && (j - v < imgC->ncol))
							(*A)[mkMs * nkMs * (u + v * p) + (i / s + j / s * mkMs)] 
							= imgC->val[j - v + imgC->ncol * (i - u)]; /*Save by cols */
						
					}
				}
				
			}
		}
    }
	
	/*Put the image in he middle, just to have a center kenel */
	mc = (p - 1) / (2 * s) + 1;
	nc = (q - 1) / (2 * s) + 1;
	
	
	for (i = 0; i < imgW->nrow; i++)
		for (j = 0; j < imgW->ncol; j++)
		{
			/*Check if it is in the mask... */
			if (imgMask->val[i * imgMask->ncol + j])
				
				(*b)[(j + nc - 1) * mkMs + i + mc - 1] =
				imgW->val[i * imgW->ncol + j];
			
			
		}
	
	
	*ncol = p * q;
	*nrow = mkMs * nkMs;
	
	free_imageFloat (mask);
	
	free ((void *) kerode);
	
	return EXIT_SUCCESS;
	
	
}


/* We call:  HR: zoomed-in image (close to the object) - "Sharp image"
 *          LR: zoomed-out image (far from the object) - "Blur image"
 */


/** 
 *  @brief Extract only the working region of image, namely the smallest 
 *  rectangle that contains all 'num_points' points 'p'.
 */
ImageFloat extract_image_region (ImageFloat in, float *p, int num_points, 
								 float *offset)
{
	float xmin, xmax, ymin, ymax;
	int xmini, xmaxi, ymini, ymaxi;	
	int i;
	ImageFloat imgT;
	
	/*Get the minimum rectangular window that covers all detected points */
	
	xmin = p[0];
	xmax = p[0];
	ymin = p[1];
	ymax = p[1];
	
	for (i = 1; i < num_points; i++)
    {
		xmin = (p[2 * i] < xmin) ? p[2 * i] : xmin;
		xmax = (p[2 * i] > xmax) ? p[2 * i] : xmax;
		ymin = (p[2 * i + 1] < ymin) ? p[2 * i + 1] : ymin;
		ymax = (p[2 * i + 1] > ymax) ? p[2 * i + 1] : ymax;
    }
	
	
    /*round the values: extract a integer subimage (to avoid interpolation)*/
	xmini = (int) xmin;
	xmaxi = (int) (xmax+1);
	ymini = (int) ymin;
	ymaxi = (int) (ymax+1);
	
	/*Check values are in bound */
	xmini = xmini >= 0 ? xmini : 0;
	xmaxi = xmaxi < in->ncol ? xmaxi : in->ncol - 1;
	ymini = ymini >= 0 ? ymin : 0;
	ymaxi = ymaxi < in->nrow ? ymaxi : in->nrow - 1;
	
	
	
	imgT = extract_window (in,  xmini,  xmaxi, ymini, ymaxi);
	
	/*update checkpoints to the new reference (xmini,ymini) */
	for (i = 0; i < num_points; i++)
    {
		p[2 * i] -= (float) xmini;
		p[2 * i + 1] -= (float) ymini;
		
    }
	
	
	/*
	xmin = p[0];
	xmax = p[0];
	ymin = p[1];
	ymax = p[1];
	
	for (i = 1; i < num_points; i++)
    {
		xmin = (p[2 * i] < xmin) ? p[2 * i] : xmin;
		xmax = (p[2 * i] > xmax) ? p[2 * i] : xmax;
		ymin = (p[2 * i + 1] < ymin) ? p[2 * i + 1] : ymin;
		ymax = (p[2 * i + 1] > ymax) ? p[2 * i + 1] : ymax;
    }
	*/
	
	offset[0] =  (float) xmini;
	offset[1] =  (float) ymini;

	
	return imgT;
	
}




void inter_image_kernel_to_psf(float *H, float* xinter, float *h,int nx, int ny)
{
	
	int it_max = 3;
	int it=0;
	int i,j;
	float lambda_x, lambda_y;
	float scale_x, scale_y;
	int nxS,nyS;
	
	ImageFloat x, xn;
	ImageFloat xGrid,yGrid,scaled_x,aux;
	
	
	/*Create image with the kernel*/
	x = new_imageFloat(nx,ny);
	xn = new_imageFloat(nx,ny);
	
	memcpy (x->val, xinter, nx * ny * sizeof (float));
	memcpy (xn->val, xinter, nx * ny * sizeof (float));
	
	/*lambda from ThinPlate Affine Part*/
	lambda_x = H[0];
	lambda_y = H[4];
		
	/* Initialize scale=lambda^n, n=1*/
	scale_x = lambda_x;
	scale_y = lambda_y;
	
	nxS = (int)(nx*scale_x+1);
	nyS = (int)(ny*scale_y+1);
	
	while (scale_x < 50 && scale_y < 50 && it< it_max)
	{
		
		nxS = (int)(nx*scale_x+1);
		nyS = (int)(ny*scale_y+1);
		
		xGrid = new_imageFloat (nxS, nyS);
		yGrid = new_imageFloat (nxS, nyS);
		
		/*Start from the center so the (((float)nx)-1)/2; is needed*/
		for (i = 0; i < xGrid->nrow; i++)
			for (j = 0; j < xGrid->ncol; j++)
			{
				xGrid->val[i * xGrid->ncol + j] = (j-(((float)nxS)-1)/2)/scale_x + (((float)nx)-1)/2;
				yGrid->val[i * xGrid->ncol + j] = (i-(((float)nyS)-1)/2)/scale_y + (((float)ny)-1)/2;
			}
		
		
		scaled_x = bicubic (xGrid, yGrid, x, -0.5);

		free_imageFloat(xGrid);
		free_imageFloat(yGrid);
		
		printf("Iteration %d :: scale x: %f, y: %f\n",it, scale_x, scale_y);
		
		
		aux = convol(scaled_x,xn);
		free_imageFloat(xn);
		free_imageFloat(scaled_x);


		xn = new_imageFloat(aux->ncol,aux->nrow);
		memcpy (xn->val, aux->val, aux->ncol * aux->nrow * sizeof (float));
		free_imageFloat(aux);
		
		scale_x =  scale_x*lambda_x;
		scale_y =  scale_y*lambda_y;
		
		it = it+1;
	
	
	}
	
	/*I added one extra scale...*/
	scale_x =  scale_x/lambda_x;
	scale_y =  scale_y/lambda_y;
	
	xGrid = new_imageFloat (nx,ny);
	yGrid = new_imageFloat (nx,ny);
	
	/*I start from the center*/
	for (i = 0; i < ny; i++)
		for (j = 0; j < nx; j++)
		{
			xGrid->val[i * xGrid->ncol + j] = (j-(((float)nx)-1)/2)*scale_x + (((float)nxS)-1)/2;
			yGrid->val[i * xGrid->ncol + j] = (i-(((float)ny)-1)/2)*scale_y + (((float)nxS)-1)/2;
		}
	
	aux = bicubic (xGrid, yGrid, xn, -0.5);
	
	/*h = (float *) malloc (nx * ny * sizeof (float));*/

	memcpy (h, aux->val, nx * ny * sizeof (float));
	
	free_imageFloat(aux);
	free_imageFloat(xGrid);
	free_imageFloat(yGrid);
	free_imageFloat(xn);
		
}


/** 
 *  @brief PSF Estimation (Main Function)
 */
void two_photos_psf_estim (float *img_LR, int nx_LR, int ny_LR,
					  float *img_HR, int nx_HR, int ny_HR,
					  int s, int psf_nrow, int psf_ncol, int threshold, 
					  float *h, float *k, char *outprefix)
{
	ImageFloat z_HR, z_LR, imgW,  xGrid, yGrid, imgC, imgP, imgPf;
	int i, j, np;
	int ncs, nrs;
	
	time_t tstart, tend;	
	
	void *swap;
	
	float *p_HR, *p_LR;
	float xmin, xmax, ymin, ymax;
 	
	float fcx, fcy;
		
	float *cblur, *csharp;
	float ps;
	
	float *A, *b, *x;
	float acsum = 0;
	
	int ncol, nrow;
	float lambda_x,lambda_y;
	char file_name[80];
			
	float max_val1, min_val1, max_val2, min_val2, min_val ,max_val, v;
	
	ImageFloat imgMask, imgMaskC, imgMaskPf, imgCxs, imgCx, imgx;
		
	float  *H;
	
	float *offset_HR, *offset_LR;
	
	char reverse;
	
	offset_HR = (float*) malloc(2*sizeof(float));
	offset_LR = (float*) malloc(2*sizeof(float));

	
	H = (float*) malloc(9*sizeof(float));

	/* convert inputs image to ImageFloat */
	z_HR = new_imageFloat (nx_HR, ny_HR);
	memcpy (z_HR->val, img_HR, nx_HR * ny_HR * sizeof (float));
	
	z_LR = new_imageFloat (nx_LR, ny_LR);
	memcpy (z_LR->val, img_LR, nx_LR * ny_LR * sizeof (float));
	
	
	/*--------------------RUNNING SIFT to ALIGN the IMAGES-----------------------*/
	printf ("Running ORSA/Homography...\n");
	/*return 0 if images do not match n_points if fimages match*/
	/*np = match_SIFT (z_in, z_out, &p_HR, &p_LR);*/
	
	
	/*Normalize input images for SIFT (0,255)*/
	/* check min and max values */
	max_val1 = 0;
	min_val1 = BIG_NUMBER;
	for(i=0;i<nx_HR*ny_HR;i++)
		{
			v = img_HR[i];
			if( v > max_val1 ) max_val1 = v;
			if( v < min_val1 ) min_val1 = v;
		}
	
	
	/* check min and max values */
	max_val2 = 0;
	min_val2 = BIG_NUMBER;
	for(i=0;i<nx_LR*ny_LR;i++)
	{
		v = img_LR[i];
		if( v > max_val2 ) max_val2 = v;
		if( v < min_val2 ) min_val2 = v;
	}

	
	/*normalize both images with the same constants*/
	max_val = (max_val2<max_val1)?max_val2:max_val1;
	min_val = (min_val2>min_val1)?min_val2:min_val1;
	
	printf("Image 1: min: %f, max: %f\n", min_val1, max_val1);
	printf("Image 2: min: %f, max: %f\n", min_val2, max_val2);
	printf("Normal:  min: %f, max: %f\n", min_val, max_val);
	
	for(i=0;i<nx_HR*ny_HR;i++)
		img_HR[i] = 255*(img_HR[i] - min_val)/(max_val-min_val);
	
	for(i=0;i<nx_LR*ny_LR;i++)
		img_LR[i] = 255*(img_LR[i] - min_val)/(max_val-min_val);
		
	
	/*
	np  = match_asift(img_HR,  nx_HR, ny_HR, img_LR, nx_LR, ny_LR, 
					  &img_siftmatchs, &img_siftmatchs_nx, &img_siftmatchs_ny, 
					  &p_HR, &p_LR);
	*/
	
	
	orsa_homography_sift(img_LR,  nx_LR,   ny_LR,
						 img_HR, nx_HR,  ny_HR,
						 1.5, H, &p_LR, &p_HR, &np, &reverse);
	
	
	printf("Reverse order : %d\n\n", reverse);
	
	
	if(reverse)
	{
		swap = z_HR;
		z_HR = z_LR;
		z_LR = (ImageFloat) swap;
		
	}
	
	
	/*inverse_homography(H,Hi);*/
	
	
	
	if (np == 0)
	{
		printf("Images do not match.\n");
		exit(1);
	}
	
	/*
	imgSIFT = (ImageFloat) malloc(sizeof(struct imageFloatStruct));
	imgSIFT->ncol = img_siftmatchs_nx;
	imgSIFT->nrow = img_siftmatchs_ny;
	imgSIFT->val = img_siftmatchs;
	
	printf("Saving SIFT matches image...\n");
	
	if(outprefix)
	{
		strcpy(file_name,outprefix);
		strcat(file_name,"_SIFT.pgm");
		write_pgm_normalize_float(file_name,imgSIFT->val, imgSIFT->ncol, imgSIFT->nrow);
		
	}
	*/
	
		
	/*--------------------------EXTRACT SUB-IMAGES------------------------------*/
	/*Extract a rectangular sub image. The minimum rectangle that contains all 
	 *the 'p_HR' points. Update the checkpoints location relative to the
	 *extracted image.
	 */
	
	
	printf ("Extracting the pattern...\n");
	imgP = extract_image_region (z_HR, p_HR, np, offset_HR);
	imgW = extract_image_region (z_LR, p_LR, np, offset_LR);
	
	
	
	/*-----------------GEOMETRIC DISTORTION - THIN PLATES-----------------------*/
	tstart = time(0);
	
	/*chequear si np=3 da error, debería funcionar y encontrar una afinidad*/
	/*tp = calculate_thinPlate (p_HR, p_LR, np, -1); 
	  tpI = calculate_thinPlate(p_LR, p_HR, np, -1);	
	*/
    
	/*Check if the images are in the right order: scale should be > 1 */
	
	/*lambda from ThinPlate Affine Part*/
	
	
	lambda_x = H[0];
	lambda_y = H[4];
	
	
	printf("Homography: \n");
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			printf("%10.10f ", H[3*i+j]);
		}
		printf("\n");
	}
	
		
		
	printf("Lambda x: %f\n",lambda_x);
	printf("Lambda y: %f\n",lambda_y);
	
		
	printf("H\n");
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			printf("%f ", H[3*i+j]);
		}
		printf("\n");
	}
	
	
	xmin = p_LR[0];
	xmax = p_LR[0];
	ymin = p_LR[1];
	ymax = p_LR[1];
	
	
	for (i = 1; i < np; i++)
	{
		xmin = (p_LR[2 * i] < xmin) ? p_LR[2 * i] : xmin;
		xmax = (p_LR[2 * i] > xmax) ? p_LR[2 * i] : xmax;
		ymin = (p_LR[2 * i + 1] < ymin) ? p_LR[2 * i + 1] : ymin;
		ymax = (p_LR[2 * i + 1] > ymax) ? p_LR[2 * i + 1] : ymax;
	}
	
	
	
	/*-----------------------ZOOMED-IN IMAGE INTERPOLATION----------------------*/
	/*
	 * fcx = imgW->ncol*s/(2*UP_RES*512);
	 * fcy = imgW->nrow*s/(2*UP_RES*512);
	 */
	
	tend = time(0);
	printf("Elapsed time: %f seconds.\n", difftime(tend,tstart));
	
	printf ("Filtering and interpolating zoomed-in pattern image...\n");
	
	fcx = roundfi (imgW->ncol * s);
	fcy = roundfi (imgW->nrow * s);
	
	printf("fcx = %f\n",fcx);
	printf("fcy = %f\n",fcy);
	
	imgPf = lpf_image_dct (imgP, (int) fcx, (int) fcy);
	
	printf("End filtering...\n");
	
	/* Pattern sx Interpolation */
	ps = 1 / ((float) s);
	
	ncs = (imgW->ncol-1) * s + 1;
	nrs = (imgW->nrow-1) * s + 1;
	
	cblur = (float *) malloc (nrs * ncs * 2 * sizeof (float));
	csharp = (float *) malloc (nrs * ncs * 2 * sizeof (float));
	
	/*Generating the sx-sampling grid */
	for (i = 0; i < nrs; i++)
		for (j = 0; j < ncs; j++)
		{
			cblur[2 * i * ncs + 2 * j] =   j * ps + offset_LR[0];
			cblur[2 * i * ncs + 2 * j + 1] =  i * ps + offset_LR[1];
		}
	
	/*Applying TP to the sx-sampling grid and then interpolate the Sharp pattern */
	printf("Start Homography Evaluation...\n");
	evaluate_homography (H, cblur, csharp, nrs * ncs);
	printf("End Homography Evaluation...\n");

	xGrid = new_imageFloat (ncs, nrs);
	yGrid = new_imageFloat (ncs, nrs);
	
	for (i = 0; i < nrs; i++)
		for (j = 0; j < ncs; j++)
		{
			xGrid->val[i * ncs + j] = csharp[2 * i * ncs + 2 * j] -  offset_HR[0];
			yGrid->val[i * ncs + j] = csharp[2 * i * ncs + 2 * j + 1] -  offset_HR[1];
		}
	
	
	printf("Start bicubic Homography interpolation...\n");
    imgC = bicubic (xGrid, yGrid, imgPf, -0.5);

	
	/*Create a Mask where the interpolation is valid
	 of course recomputing the interpolation  is not optimal but it works
	 */
	
	imgMaskPf = new_imageFloat(imgPf->ncol,imgPf->nrow);
	
	for(i=0;i<imgMaskPf->ncol*imgMaskPf->nrow;i++)
		imgMaskPf->val[i] = 1;
	
	imgMaskC = bicubic (xGrid, yGrid, imgMaskPf, -0.5);

	
	
    imgMask = new_imageFloat(imgW->ncol,imgW->nrow);
	
	for(i=0;i<imgMask->nrow;i++)
		for(j=0;j<imgMask->ncol;j++)
			imgMask->val[i * imgMask->ncol +j] = imgMaskC->val[i*s*imgMaskC->ncol+j*s];
	printf("End bicubic Homography interpolation...\n");

	
	/*Generating A and b for Ax = b system */
	printf ("Generating A,b for Ax = b linear system...\n");
	
	make_Ab (imgC, imgW, imgMask, psf_nrow, psf_ncol, s, &A, &b, &ncol, &nrow);
	
	
	/*
	write_ascii_matrix(A, ncol, nrow, "A.txt");
	write_ascii_matrix(b, 1, nrow, "b.txt");
	*/
	
	/*Solve the inter-image kernel by using Least Squares algorithm*/
	printf ("Estimating the PSF: solving Ax = b...\n");
	x = solve_lsd(A, b, nrow, ncol);
	
	if(outprefix)
	{
		
		/*all output images are re-scaled to [0,255]
		 normalizing with the imgC values
		 Difference image is normalized to 15%of the max-min values of imgC.
		 
		 */
		
				
		/*The normalization is done with the max value of imgC*/
		max_val1 = 0;
		min_val1 = BIG_NUMBER;
		
		for(i=0; i< imgC->nrow; i++)
			for(j=0; j< imgC->ncol; j++)
			{
				v = imgC->val[ j + i * imgC->ncol];
				if( v > max_val1 ) max_val1 = v;
				if( v < min_val1 ) min_val1 = v;
			}
		
		max_val2 = 0;
		min_val2 = BIG_NUMBER;
		
		for(i=0; i< imgW->nrow; i++)
			for(j=0; j< imgW->ncol; j++)
			{
				v = imgW->val[ j + i * imgW->ncol];
				if( v > max_val2 ) max_val2 = v;
				if( v < min_val2 ) min_val2 = v;
			}
		
		
		/*normalize evey images with the same constants*/
		max_val = (max_val2>max_val1)?max_val2:max_val1;
		min_val = (min_val2<min_val1)?min_val2:min_val1;
		
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_imgC.pgm");
		write_pgm_normalize_given_minmax_float(file_name,imgC->val, imgC->ncol,
											   imgC->nrow, min_val, max_val);
		
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_imgC.txt");
		write_ascii_imageFloat (imgC, file_name);
		
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_imgW.pgm");
		write_pgm_normalize_given_minmax_float(file_name,imgW->val, imgW->ncol, imgW->nrow, 
								  min_val, max_val);
		
		printf("%f %f\n", min_val, max_val);
		strcpy(file_name,outprefix);
		strcat(file_name,"_imgW.txt");
		write_ascii_imageFloat (imgW, file_name);
		
		
		
		/*Write Mask image*/

		strcpy(file_name,outprefix);
		strcat(file_name,"_mask.txt");
		write_ascii_imageFloat (imgMask, file_name);
		
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_mask.pgm");
		write_pgm_normalize_float(file_name,imgMask->val, imgMask->ncol,
								  imgMask->nrow);
		
		
		
		/*Compute difference image */
		imgx = new_imageFloat(psf_ncol,psf_nrow);
		memcpy (imgx->val, x, psf_ncol * psf_nrow * sizeof (float));
		imgCx = convol(imgC, imgx);
		
		/*Subsampling & Difference*/
		imgCxs = new_imageFloat(imgW->ncol, imgW->nrow);
		for(i=0;i<imgW->nrow;i++)
			for(j=0;j<imgW->ncol;j++)
				if(imgMask->val[i*imgW->ncol + j])
					imgCxs->val[i*imgCxs->ncol +j] =
					imgCx->val[s*i*imgCx->ncol + s*j] - imgW->val[i*imgW->ncol + j]; 
				else 
					imgCxs->val[i*imgCxs->ncol +j] = 0;	
		
		
		/*The normalization is done with the 20% of min_val
		 and 10% of max_vale*/
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_diff.pgm");
		write_pgm_normalize_given_minmax_float(file_name,imgCxs->val, imgCxs->ncol, imgCxs->nrow, 
								  0.15*(-max_val + min_val) , 0.1*(max_val - min_val));
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_diff.txt");
		write_ascii_imageFloat (imgCxs, file_name);
	
	
		
		
	}
	
	
	
	inter_image_kernel_to_psf(H, x, h, psf_ncol, psf_nrow);
	
	/* threshold==1 :threshold the final psf estimation*/
	if (threshold)
	{
		for (i = 0; i < psf_ncol * psf_nrow; i++)
			h[i] = (h[i] > POSITIVE_TOL) ? h[i] : 0;
	}
	
	
	/*inter-image kernel should be sum_i x[i] = 1*/
	acsum = 0;
	for (i = 0; i < psf_ncol * psf_nrow; i++)
		acsum = acsum + x[i];
	
	for (i = 0; i < psf_ncol * psf_nrow; i++)
		k[i] = x[i]/acsum;
	
	
	
	/*PSF should be sum_i h[i] = 1*/
	acsum = 0;
	for (i = 0; i < psf_ncol * psf_nrow; i++)
		acsum = acsum + h[i];
	
	for (i = 0; i < psf_ncol * psf_nrow; i++)
		h[i] = h[i]/acsum;
	
	
	printf ("Cleaning the house...\n");
	
	free_imageFloat (imgW);
	free_imageFloat (xGrid);
	free_imageFloat (yGrid);
	free_imageFloat (imgC);
/*	free_imageFloat (imgP);
	free_imageFloat (imgPf);*/
	free_imageFloat (z_LR);
	free_imageFloat (z_HR);	
	
	free ((void *) cblur);
	free ((void *) csharp);
	free ((void *) p_HR);
	free ((void *) p_LR);	
	free ((void *) A);
	free((void*) b);
	
		
}


