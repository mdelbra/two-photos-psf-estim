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
 * @file two_photos_psf_estim.c
 * @brief library code to two-photos psf estimation.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */

/** @mainpage Recovering the Subpixel PSF from Two Photographs at Different 
 *   Distances
 *
 * The following is an implementation of 
 * \li  M. Delbracio, A. Almansa, J.-M. Morel and P. Muse.
 * "Subpixel Point Spread Function Estimation from Two Photographs 
 *  at Different Distances", SIAM Journal on Imaging Sciences (SIIMS)
 *  November 2012.
 *
 * A detail desription and an online demo can be accessed from:
 *
 * \li to_be_updated_final_url
 *
 * The source code consists of:
 *
 *    \li    image.c
 *    \li    image.h
 *    \li    io_pgm.c
 *    \li    io_pgm.h
 *    \li    ls.c
 *    \li    ls.h
 *    \li    two_photos_psf_estim.c
 *    \li    two_photos_psf_estim.h
 *    \li    two_photos_psf_estim_main.c
 *    \li    third_party: ORSA/Homography
 *
 * HISTORY:
 * - version 1.0 - feb 2013: First Release
 * - version 0.4 - sep 2012: Second BETA Release Ansi C Language
 *	version.
 *
 * @author mauricio delbracio (mdelbra@gmail.com)
 * @date feb 2011
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "image.h"
#include "ls.h"
#include "io_pgm.h"


/* (p-1)/2 / s = K
 2Ks+1 = p;
 K=1 - p = 2s+1 7x7
 K=2 - p = 4s+1 13x13
 K=3 - p = 6s+1 19x19
 K=4 - p = 8s+1 25x25
 */
 
/*External function FROM ORSA-Homography 
 "Automatic Homographic Registration of a Pair of Images, with A Contrario
  Elimination of Outliers" L. Moisan, P. Moulon, P. Monasse
  Image Processing On Line, 2012,  DOI: 10.5201/ipol.2012.mmm-oh
 */
int orsa_homography_sift(float* img1, 
                         int nx1, 
                         int ny1,
                         float* img2,
                         int nx2, 
                         int ny2,
                         double precision, 
                         float *hom,
						 float **matchs1,
                         float **matchs2,
                         int* nm, 
                         char* reverse);



/** If the  value of the PSF is less than POSITIVE_TOL consider it is zero */
#define POSITIVE_TOL 0.000001
#define BIG_NUMBER 1000000


/** 
 *  @brief Generate the Linear System Ax = b where x is the kernel to find.
 */
void make_Ab (ImageFloat imgC, ImageFloat imgW, ImageFloat imgMask,
             int q, int p, int s, float **A, float *b[], int *ncol, int *nrow)
{
    /*Generate a SsU from U and s. */
    int max_pq, u, v, i, j, mkMs, nkMs, mc, nc, r;
    
    float *kerode;
    ImageFloat mask;
    /* generate SsU matrix */
    /* p,q kernel size */
    
    /*size of the final matrix A: (mkMs*nKMs)x(p*q) */
    mkMs = (imgC->nrow + p - 2) / s + 1;
    nkMs = (imgC->ncol + q - 2) / s + 1;
    
    *A = (float *)calloc (mkMs * nkMs * p * q, sizeof (float));/*ini. to zero*/
    *b = (float *)calloc (mkMs * nkMs, sizeof (float));/*ini. to zero */
    
    
    /*Erode the Mask in order take into account the boundary
     * problems
     */
    /*square element of side r */
    max_pq = (p > q) ? p : q;
    r = (max_pq - 1) / (2 * s);
    kerode = (float *) malloc ((2 * r + 1) * sizeof (float));/*ini. to zero */
    
    for (i = 0; i < 2 * r + 1; i++)
        kerode[i] = 1;
    
    mask = convol_sep2 (imgMask, kerode, 2 * r + 1, kerode, 2 * r + 1);
    
    for (i = 0; i < imgMask->nrow; i++)
        for (j = 0; j < imgMask->ncol; j++)
            imgMask->val[i * imgMask->ncol + j] =
            (mask->val[(i + r) * mask->ncol + j + r] >= 4 * r * r + 1) ? 1 : 0;
    
    
    int p2 = (p-1)/2;
    int q2 = (q-1)/2;
    
    /*Filling in (*A) */
    for (i = p2; i < (imgC->nrow) + p2; i = i + s)
    {
        for (j = q2; j < (imgC->ncol) + q2; j = j + s)
        {
            if (imgMask->val[(j - q2) / s + imgMask->ncol * (i - p2) / s])
            {
                for (u = 0; u < p; u++)
                {
                    for (v = 0; v < q; v++)
                    {
                        if ((i - u >= 0) && (j - v >= 0)
                            && (i - u < imgC->nrow) && (j - v < imgC->ncol))
                            (*A)[mkMs * nkMs * (v + u * q) + (i/s + j/s*mkMs)]
                            = imgC->val[j - v + imgC->ncol * (i - u)];
                                /*Save by cols*/
                    }
                }
            }
        }
    }
    
    /* Put the image in he center (just to have a centered kenel) */
    mc = p2/s + 1;
    nc = q2/s + 1;
    
    
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
    
}



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
	xmini = floor(xmin);
	xmaxi = ceil(xmax);
	ymini = floor(ymin);
	ymaxi = ceil(ymax);
	
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
		
	offset[0] =  (float) xmini;
	offset[1] =  (float) ymini;

	
	return imgT;
	
}


/** 
 *  @brief Interpolate the  closest image \f$\textbf{v}_1\f$
 *  to generate \f$\textbf{H}_{\lambda/s} \tilde{\textbf{v}}_1\f$ 
 */
 void closest_image_interpolation(float* H, int s, ImageFloat imgW, 
                                  ImageFloat imgP, ImageFloat* imgC, 
                                  ImageFloat* imgMask, float* offset_LR,
                                  float* offset_HR)
{
 
    	
	
    ImageFloat imgPf, imgMaskPf, imgMaskC, xGrid, yGrid;
    float *coordImgLRsx, *coordImgHRsx;
    float ps;
    int ncs, nrs, i, j;

	/* Cut all frequencies above fcx or fcy:
	 * fcx = imgW->ncol*s/(2*UP_RES*512);
	 * fcy = imgW->nrow*s/(2*UP_RES*512);
	 */
	
    float fcx,fcy;
    
	fcx = roundfi (imgW->ncol * s);
	fcy = roundfi (imgW->nrow * s);
		
	/*Filter to avoid aliasing*/
    imgPf = lpf_image_dct (imgP, (int) fcx, (int) fcy);
		
	/* Image sx Interpolation H_lambda/s v1 */
	ps = 1 / ((float) s);
	
	ncs = (imgW->ncol-1) * s + 1;
	nrs = (imgW->nrow-1) * s + 1;
	
    /*Sampling Grids points*/
	coordImgLRsx = (float *) malloc (nrs * ncs * 2 * sizeof (float));
	coordImgHRsx = (float *) malloc (nrs * ncs * 2 * sizeof (float));
	
	/*Generating the sx-sampling grid */
	for (i = 0; i < nrs; i++)
		for (j = 0; j < ncs; j++)
		{
            /*Add the offset to reach the common part*/
			coordImgLRsx[2 * i * ncs + 2 * j] =   j * ps + offset_LR[0];
			coordImgLRsx[2 * i * ncs + 2 * j + 1] =  i * ps + offset_LR[1];
		}
	
	/*Apply Homography to the sx-sampling grid and then interpolate the 
     *closest image */
	evaluate_homography (H, coordImgLRsx, coordImgHRsx, nrs * ncs);
    
    /*Sampling Grids at the common region*/
	xGrid = new_imageFloat (ncs, nrs);
	yGrid = new_imageFloat (ncs, nrs);
	
	for (i = 0; i < nrs; i++)
		for (j = 0; j < ncs; j++)
		{
            /*Remove the offset to reach the common region*/
			xGrid->val[i * ncs + j] = coordImgHRsx[2 * i * ncs + 2 * j] 
            - offset_HR[0];
			yGrid->val[i * ncs + j] = coordImgHRsx[2 * i * ncs + 2 * j + 1]  
            - offset_HR[1];
		}
	
	
    /*Bicubic interpolation of the filtered closest image at the sampling grid
     *sx the  resolution of the farthest image
     */
    *imgC = bicubic (xGrid, yGrid, imgPf, -0.5);
    
	/*Create a Mask where the interpolation is valid
	 *of course recomputing the interpolation  is not optimal but it works.
	 */
	
	imgMaskPf = new_imageFloat(imgPf->ncol,imgPf->nrow);
	for(i=0;i<imgMaskPf->ncol*imgMaskPf->nrow;i++) imgMaskPf->val[i] = 1;
    
    /*Mask at the sx grid*/
	imgMaskC = bicubic (xGrid, yGrid, imgMaskPf, -0.5);
    
	
	/*Mask at the 1x grid*/
    *imgMask = new_imageFloat(imgW->ncol,imgW->nrow);
	for(i=0;i< (*imgMask)->nrow;i++)
		for(j=0;j< (*imgMask)->ncol;j++)
			(*imgMask)->val[i * (*imgMask)->ncol +j] = 
            imgMaskC->val[i*s*imgMaskC->ncol+j*s];
    
    
    /*Do the cleaning*/
    free_imageFloat (imgPf);
    free_imageFloat (imgMaskPf);
	free_imageFloat (imgMaskC);
	free_imageFloat (xGrid);
	free_imageFloat (yGrid);
	
	free ((void *) coordImgLRsx);
	free ((void *) coordImgHRsx);
    
    
}
    




void inter_image_kernel_to_psf(float *H, float* xinter, float *h, 
                               int nx, int ny)
{
	
	int it_max = 3;
	int it=0;
	int i,j;
	float lambda_x, lambda_y;
	float scale_x, scale_y, nx2, ny2, nxS2, nyS2, tx, ty;
	int nxS,nyS;
	
	ImageFloat x, xn, xGrid,yGrid,scaled_x,aux;
	
    nx2 = (((float)nx)-1)/2;
    ny2 = (((float)ny)-1)/2;
    
	
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
	
	/*nxS = (int)(nx*scale_x+1);
	nyS = (int)(ny*scale_y+1);
	*/
	while (scale_x < 50 && scale_y < 50 && it< it_max)
	{
		nxS = ceil(nx*scale_x);
		nyS = ceil(ny*scale_y);
     
        nxS2 = (((float)nxS)-1)/2;
        nyS2 = (((float)nyS)-1)/2;
        
		xGrid = new_imageFloat (nxS, nyS);
		yGrid = new_imageFloat (nxS, nyS);
		
        tx = nx2 - nxS2/scale_x;
        ty = ny2 - nyS2/scale_y;
                
		/*Start from the center so the (tx,ty); is needed*/
		for (i = 0; i < xGrid->nrow; i++)
			for (j = 0; j < xGrid->ncol; j++)
			{
				xGrid->val[i * xGrid->ncol + j] = ((float)j)/scale_x + tx; 
				yGrid->val[i * xGrid->ncol + j] = ((float)i)/scale_y + ty;
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
	
    tx = nxS2 - nx2*scale_x;
    ty = nyS2 - ny2*scale_y;
    
    
	/*I start from the center*/
	for (i = 0; i < ny; i++)
		for (j = 0; j < nx; j++)
		{
			xGrid->val[i * xGrid->ncol + j] = j*scale_x + tx;
			yGrid->val[i * xGrid->ncol + j] = i*scale_y + ty;
		}
	
	aux = bicubic (xGrid, yGrid, xn, -0.5);
	
	/*h = (float *) malloc (nx * ny * sizeof (float));*/

	memcpy (h, aux->val, nx * ny * sizeof (float));
	
	free_imageFloat(aux);
	free_imageFloat(xGrid);
	free_imageFloat(yGrid);
	free_imageFloat(xn);
	free_imageFloat(x);
		
    
}




/** 
 *  @brief PSF Estimation (Main Function)
 */
void two_photos_psf_estim (float *img1, int nx1, int ny1,
                           float *img2, int nx2, int ny2,
                           int s, int psf_ncol, int psf_nrow,
                           float *h, float *k, 
                           int threshold, char *outprefix)
{
	
    
    ImageFloat z_HR, z_LR, imgW,  imgC, imgP;
    ImageFloat imgMask,  imgCxs, imgCx, imgx;

    int i, j, np, ncol, nrow;
		
	float *p_HR, *p_LR;
    float max_val1, min_val1, max_val2, min_val2, min_val ,max_val, v;

	float *A, *b;

	char file_name[80];

	char reverse;
    float offset_HR[2], offset_LR[2], H[9];

	
    /*------------------------------------------------------------------------*/
	/*- STEP 1--- IMAGE ALIGNMENT, GEOMETRIC TRANSFORMATION ESTIOMATION-------*/
    /*------------------------------------------------------------------------*/	
    
    /* ORSA/Homography (SIFT) Estimate a Homography between the input images */

	printf ("Running ORSA/Homography...\n");	
	
	orsa_homography_sift(img1, nx1, ny1, img2, nx2, ny2, 1.5, H,
                         &p_LR, &p_HR, &np, &reverse);
	
	
    if (np == 0)
	{
		printf("Images do not match.\n");
		exit(1);
	}
    
	
    if(!reverse)
        
    {
        printf(" First input image: farthest image\n");
        printf("Second input image: closest image\n");
        
        /* convert input image 2 to ImageFloat */
        z_HR = new_imageFloat (nx2, ny2);
        memcpy (z_HR->val, img2, nx2 * ny2 * sizeof (float));
        
        /* convert input image 1 to ImageFloat */
        z_LR = new_imageFloat (nx1, ny1);
        memcpy (z_LR->val, img1, nx1 * ny1 * sizeof (float));
        
    } else
	{
        
        printf(" First input image: closest image\n");
        printf("Second input image: farthest image\n");

        /* convert input image 1 to ImageFloat */
        z_HR = new_imageFloat (nx1, ny1);
        memcpy (z_HR->val, img1, nx1 * ny1 * sizeof (float));
        
        /* convert input image 1 to ImageFloat */
        z_LR = new_imageFloat (nx2, ny2);
        memcpy (z_LR->val, img2, nx2 * ny2 * sizeof (float));
     	
	}
	
   
    /*------------------------------------------------------------------------*/
	/*- STEP 2-------EXTRACT COMMON REGION (SUBIMAGES)------------------------*/
    /*------------------------------------------------------------------------*/	
	
	/*Extract a rectangular sub image. The minimum rectangle that contains all 
	 *the 'p_HR' points in 'z_HR' and 'p_LR' in 'z_LR'. 
     *Update the checkpoints location relative to the extracted image.
	 */
	
	imgP = extract_image_region (z_HR, p_HR, np, offset_HR);
	imgW = extract_image_region (z_LR, p_LR, np, offset_LR);
	

    /*------------------------------------------------------------------------*/
	/*- STEP 3-------IMAGE INTERPOLATION: H_\frac{\lambda}{s} \tilde{v}_1-----*/
    /*------------------------------------------------------------------------*/	
	
    
    closest_image_interpolation(H, s, imgW, imgP, &imgC, &imgMask,
                                offset_LR, offset_HR);
	

    /*------------------------------------------------------------------------*/
	/*- STEP 4-------GENERATING LINEAR SYSTEM MS_sUk = M\tilde{v}_2-----------*/
    /*------------------------------------------------------------------------*/	
	
    /* A = MS_sUk,  b = M\tilde{v}_2*/
  	
	make_Ab (imgC, imgW, imgMask, psf_ncol, psf_nrow, s, &A, &b, &ncol, &nrow);
	

    /*------------------------------------------------------------------------*/
	/*- STEP 5-------SOLVING LINEAR SYSTEM k/ Ak = b -------------------------*/
    /*------------------------------------------------------------------------*/	
	
    solve_lsd(A, b, k, ncol, nrow);
	
	
    /*------------------------------------------------------------------------*/
	/*- STEP 6-------FROM INTER-IMAGE-KERNEL TO PSF---------------------------*/
    /*------------------------------------------------------------------------*/	
    
	inter_image_kernel_to_psf(H, k, h, psf_ncol, psf_nrow);

	/* threshold==1 :threshold the final psf estimation*/
	if (threshold)
	{
		for (i = 0; i < psf_ncol * psf_nrow; i++)
			h[i] = (h[i] > POSITIVE_TOL) ? h[i] : 0;
	}
	
	
	/*inter-image kernel should be sum_i x[i] = 1*/
    normalize_area(k,psf_ncol * psf_nrow);
    
    /*PSF should be sum_i h[i] = 1*/
    normalize_area(h,psf_ncol * psf_nrow);
	
    
    /*------------------------------------------------------------------------*/
	/*----------SAVE INTERMEDIATE IMAGES, RESIDUAL----------------------------*/
    /*------------------------------------------------------------------------*/	
    
    if(outprefix)
	{
		
		/*all output images are re-scaled to be in [0,255]
		 normalizing with the max and min values of [imgC,imgW] values
		 Difference image is normalized to [-0.05(max-min),0.05(max-min)]
         values of [imgC,imgW].
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
		write_pgm_normalize_given_minmax_float(file_name, imgC->val, 
                                               imgC->ncol, imgC->nrow, 
                                               min_val, max_val);
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_imgC.txt");
		write_ascii_imageFloat (imgC, file_name);
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_imgW.pgm");
		write_pgm_normalize_given_minmax_float(file_name,imgW->val, imgW->ncol, 
                                               imgW->nrow, min_val, max_val);
        strcpy(file_name,outprefix);
		strcat(file_name,"_imgW.txt");
		write_ascii_imageFloat (imgW, file_name);
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_mask.txt");
		write_ascii_imageFloat (imgMask, file_name);
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_mask.pgm");
		write_pgm_normalize_float(file_name,imgMask->val, imgMask->ncol,
								  imgMask->nrow);
		
		
		
		/*Compute difference image */
		imgx = new_imageFloat(psf_ncol,psf_nrow);
		memcpy (imgx->val, k, psf_ncol * psf_nrow * sizeof (float));
		imgCx = convol(imgC, imgx);
		
		/*Subsampling & Difference*/
		imgCxs = new_imageFloat(imgW->ncol, imgW->nrow);
		for(i=0;i<imgW->nrow;i++)
			for(j=0;j<imgW->ncol;j++)
				if(imgMask->val[i*imgW->ncol + j])
                    imgCxs->val[i*imgCxs->ncol +j] =
                    imgCx->val[s*i*imgCx->ncol + s*j] 
                        -imgW->val[i*imgW->ncol + j]; 
				else 
					imgCxs->val[i*imgCxs->ncol +j] = 0;	
		
		
		/*The normalization of the image difference is done with so that the 
         values are in [-0.05(max_val-min_val),0.05(max_val-min_val)] (i.e.
         the dynamical range is compressed to 10%)*/
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_diff.pgm");
        
		write_pgm_normalize_given_minmax_float(file_name,imgCxs->val, 
                                               imgCxs->ncol, imgCxs->nrow, 
                                               -0.05*(max_val - min_val),
                                               0.05*(max_val - min_val));
		
		strcpy(file_name,outprefix);
		strcat(file_name,"_diff.txt");
		write_ascii_imageFloat (imgCxs, file_name);
        
		
	}
    
        
    free_imageFloat (imgMask);	
	free_imageFloat (imgW);
	free_imageFloat (imgC);
	free_imageFloat (imgP);
	free_imageFloat (z_LR);
	free_imageFloat (z_HR);	
	
	free ((void *) p_HR);
	free ((void *) p_LR);	
	free ((void *) A);
	free ((void *) b);
	
		
}


