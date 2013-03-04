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
 * @file image.c
 * @brief library code for basic image processing.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "image.h"



#define min(a,b)        (((a)>(b)) ? (b) : (a))
#define max(a,b)        (((a)>(b)) ? (a) : (b))


/** @brief Error/Exit print a message and exit.
 *  @param msg
 */
 void error(char *msg)
{
    fprintf(stderr, "PSF_ESTIM Error: %s\n", msg);
    exit(EXIT_FAILURE);
}


/**
 * @brief Free memory used in ImageFloat 'i'
 * @param i
 */
void free_imageFloat(ImageFloat i)
{
    if (i == NULL || i->val == NULL)
    error("free_image_float: invalid input image.");
    free((void *) i->val);
    free((void *) i);
}

/**
 * @brief Create new ImageFloat of size 'nrow' x 'ncol'
 * @param ncol - number of columns
 * @param nrow - number of rows
 * @return created ImageFloat
 */
ImageFloat new_imageFloat(int ncol, int nrow)
{
    ImageFloat image;

    if (ncol == 0 || nrow == 0)
    error("new_image_float: invalid image size.");

    image = (ImageFloat) malloc(sizeof(struct imageFloatStruct));
    if (image == NULL)
    error("not enough memory.");
    image->val = (float *) calloc(ncol * nrow, sizeof(float));
    if (image->val == NULL)
    error("not enough memory.");

    image->ncol = ncol;
    image->nrow = nrow;

    return image;
}



/**
 * @brief Separable full convolution between an ImageFloat and a kernel
 * @param in - Input Image Float
 * @param xker - horizontal kernel (array of floats)
 * @param xsize - length of horizontal kernel
 * @param yker - vertical kernel (array of floats)
 * @param ysize - length of vertical kernel
 * @return convolved ImageFloat
 */
ImageFloat convol_sep2(ImageFloat in, float *xker, int xsize, float *yker,
                       int ysize)
{
    ImageFloat aux, out;
    int N, M, n, m, q;
    int h;
    float sum;

    /* get memory for images */
    N = in->ncol + xsize - 1;
    M = in->nrow + ysize - 1;
    aux = new_imageFloat(N, in->nrow);
    out = new_imageFloat(N, M);

    /*convolution boundry value: 0 */

    /* First sampling: x axis */
    for (n = 0; n < aux->ncol; n++)
    {
        for (m = 0; m < aux->nrow; m++)
        {
            sum = 0.0;
            for (q = 0; q < xsize; q++)
            {
                h = n - q;
                /* null boundary condition */
                if (h >= 0 && h < in->ncol)
                    sum += in->val[h + m * in->ncol] * xker[q];
            }
            aux->val[n + m * aux->ncol] = sum;
        }
    }

    /* Second sampling: y axis */
    for (m = 0; m < out->nrow; m++)
    {
        for (n = 0; n < out->ncol; n++)
        {
            sum = 0.0;
            for (q = 0; q < ysize; q++)
            {
                h = m - q;

                /* null boundary condition */
                if (h >= 0 && h < in->nrow)
                    sum += aux->val[n + h * aux->ncol] * yker[q];
            }
            out->val[n + m * out->ncol] = sum;
        }
    }

    /* free memory */ ;
    free_imageFloat(aux);

    return out;
}



/**
 * @brief Central part of the convolution between an ImageFloat and a kernel
 * @param in - Input Image Float
 * @param ker -  Input kernel
 * @return convolved ImageFloat (same size as in)
 */
ImageFloat convol(ImageFloat in, ImageFloat kernel)
{
    int n,m,k,l,dxS,dyS,kmax,kmin,lmax,lmin,K2,L2;
    double S;
    float *ptrO,*ptrI,*ptrF;
    ImageFloat out;
    
    /*create output image - central part of the full convolution*/
    out=new_imageFloat(in->ncol, in->nrow);
    if (!out) error("not enough memory\n"); 
    
    K2 = kernel->ncol / 2;
    L2 = kernel->nrow / 2;
    dxS = in->ncol;
    dyS = in->nrow;
    
    ptrO = out->val;
    ptrI = in->val;
    
    for(m=0;m<in->nrow;m++) 
        for (n=0;n<in->ncol;n++) 
        {
            S = 0.0;
            kmax = min(kernel->ncol-1,n+K2);
            kmin = max(0,1+n+K2-dxS);
            lmax = min(kernel->nrow-1,m+L2);
            lmin = max(0,1+m+L2-dyS);
            
            ptrF = kernel->val;
            for (l=lmin;l<=lmax;l++) 
                for (k=kmin;k<=kmax;k++) 
                    S += ptrI[dxS*(m-l+L2) + (n-k+K2)] * ptrF[kernel->ncol*l+k];
            *ptrO++ = (float) S;
        }
    
    return out;
}



/**
 * @brief Extract subimage from  integer pixels 'xmin' to 'xmax'
 *        and 'ymin' to 'ymax'
 * @param in - input ImageFloat
 * @param xmin
 * @param xmax
 * @param ymin
 * @param ymax
 * @return ImageFloat created from [xmin,xmax]X[ymin,ymax] region
 */
ImageFloat extract_window(ImageFloat in, int xmin,
                          int xmax, int ymin, int ymax)
{
    ImageFloat out;
    int i, j;

    out = new_imageFloat(xmax - xmin + 1, ymax - ymin + 1);

    for (i = 0; i < out->nrow; i++)
    for (j = 0; j < out->ncol; j++)
    {
        out->val[j + i * out->ncol] =
        in->val[j + xmin + (i + ymin) * in->ncol];
    }

    return out;
}




/**
 * @brief Bicubic image interpolation at coordinates given by 'X' and 'Y'.
 * @details It is assumed that the input image is sampled in a rectangular grid
 *          from (0,nc-1)X(0,nr-1). X and Y must be the same size and
 *          the output will be also the same size. Where there is not enough
 *          information to interpolate (i.e. in the boundary) the function
 *          returns 0
 *
 * @param X - input ImageFloat with the horizontal coordinates where the
 *            interpolation is demanded
 * @param Y - input ImageFloat with the vertical coordinates where the
 *        interpolation is demanded
 * @param in - input ImageFloat
 * @param a -  Bicubic inerpolator parameter (tipically -0.5)
 * @return ImageFloat with the interpolated values at (X(i,j),Y(i,j))
 *         positions
 */
ImageFloat bicubic(ImageFloat X, ImageFloat Y, ImageFloat in, float a)
{
    /* X and Y are the coordinates where the output image is going to be
     * interpolated. Both images should be of the same time and equal the
     * size of the output image. I assume that image in is sampled in a
     * rectangular regular grid from x=0:nc-1 and y=0:nr-1. Where there is
     * not enough information to interpolate (i.e. in the boundary) the
     * function returns 0.
     */

    ImageFloat out;

    int nco = X->ncol;
    int nro = X->nrow;
    int nci = in->ncol;
    int nri = in->nrow;

    float cx[4], cy[4];

    int xinf, yinf;
    float xp, yp, res, t, s, at, as, t2, s2;

    int i, j;

    /*Create Output image */
    out = new_imageFloat(nco, nro);

    /*----Bicubic interpolation - Central Loop----*/

    for (j = 0; j < nco; j++)
    {
        for (i = 0; i < nro; i++)
        {
            xp = X->val[i * nco + j];
            yp = Y->val[i * nco + j];

            xinf = floor(xp);
            yinf = floor(yp);

            t = xp - (float) xinf;
            t2 = t * t;
            at = a * t;
            cx[0] = a * t2 * (1.0 - t);
            cx[1] = (2.0 * a + 3.0 - (a + 2.0) * t) * t2 - at;
            cx[2] = ((a + 2.0) * t - a - 3.0) * t2 + 1.0;
            cx[3] = a * (t - 2.0) * t2 + at;
            s = yp - (float) yinf;
            s2 = s * s;
            as = a * s;
            cy[0] = a * s2 * (1.0 - s);
            cy[1] = (2.0 * a + 3.0 - (a + 2.0) * s) * s2 - as;
            cy[2] = ((a + 2.0) * s - a - 3.0) * s2 + 1.0;
            cy[3] = a * (s - 2.0) * s2 + as;

            /*Put 0 in the border.... */
            if (xinf > nci - 3 || xinf < 1 || yinf > nri - 3 || yinf < 1)
                res = 0;
            else
            /* Re-check if this is not separable...or something */
                res = cy[0] * (cx[0] * in->val[(yinf + 2) * nci + xinf + 2] +
                               cx[1] * in->val[(yinf + 2) * nci + xinf + 1] +
                               cx[2] * in->val[(yinf + 2) * nci + xinf] +
                               cx[3] * in->val[(yinf + 2) * nci + xinf - 1]
                               ) +
                    cy[1] * (cx[0] * in->val[(yinf + 1) * nci + xinf + 2] +
                             cx[1] * in->val[(yinf + 1) * nci + xinf + 1] +
                             cx[2] * in->val[(yinf + 1) * nci + xinf] +
                             cx[3] * in->val[(yinf + 1) * nci + xinf - 1]
                             ) +
                    cy[2] * (cx[0] * in->val[yinf * nci + xinf + 2] +
                             cx[1] * in->val[yinf * nci + xinf + 1] +
                             cx[2] * in->val[yinf * nci + xinf] +
                             cx[3] * in->val[yinf * nci + xinf - 1]
                             ) +
                    cy[3] * (cx[0] * in->val[(yinf - 1) * nci + xinf + 2] +
                             cx[1] * in->val[(yinf - 1) * nci + xinf + 1] +
                             cx[2] * in->val[(yinf - 1) * nci + xinf] +
                             cx[3] * in->val[(yinf - 1) * nci + xinf - 1]
                             );
            out->val[i * nco + j] = res;
        }
    }

    return out;
}





/**
 * @brief Compute the DCT Transform of ImageFloat 'in'
 * @param in - input ImageFloat
 * @return ImageFloat computed DCT image
 */
ImageFloat compute_dct_image(ImageFloat in)
{

    ImageFloat out;
    fftwf_plan fw_plan;
    int nx, ny;
    nx = in->ncol;
    ny = in->nrow;

    out = new_imageFloat(nx, ny);

    /*Be careful! the order of the parameters: ny then nx
     "The multi-dimensional arrays passed to fftw_plan_dft etcetera are
     expected to be stored as a single contiguous block in row-major order
     (sometimes called “C order”).
     Basically, this means that as you step through adjacent
     memory locations, the first dimension's index varies most slowly and the
     last dimension's index varies most quickly."*/
    fw_plan = fftwf_plan_r2r_2d(ny, nx, in->val, out->val, FFTW_REDFT10,
                                FFTW_REDFT10, FFTW_ESTIMATE);
    fftwf_execute(fw_plan);

    /* Do the cleaning */
    fftwf_destroy_plan(fw_plan);

    return out;
}

/**
 * @brief Compute  the DCT Inverse Transform of ImageFloat 'in'
 * @param in - input ImageFloat
 * @return ImageFloat computed DCT image
 */
ImageFloat compute_idct_image(ImageFloat in)
{

    ImageFloat out;
    fftwf_plan bw_plan;
    int nx, ny;
    nx = in->ncol;
    ny = in->nrow;

    out = new_imageFloat(nx, ny);

    /*Be careful! the order of the parameters: ny then nx
     "The multi-dimensional arrays passed to fftw_plan_dft etcetera are
     expected to be stored as a single contiguous block in row-major order
     (sometimes called “C order”).
     Basically, this means that as you step through adjacent
     memory locations, the first dimension's index varies most slowly and the
     last dimension's index varies most quickly."*/

    bw_plan = fftwf_plan_r2r_2d(ny, nx, in->val, out->val, FFTW_REDFT01,
                                FFTW_REDFT01, FFTW_ESTIMATE);

    fftwf_execute(bw_plan);

    /* Do the cleaning */
    fftwf_destroy_plan(bw_plan);

    return out;
}


/*----------------------------------------------------------------------------
 The followoing function for computing the samples of a Gaussian kernel is 
 based on those of LSD - Line Segment Detector on digital images 
 Copyright 2007-2010 rafael grompone von gioi (grompone@gmail.com)
 ----------------------------------------------------------------------------
 */

/**
 * @brief Compute a Gaussian kernel of length 'nsize' centered at zero.
 * @param sigma - standard deviation 
 * @param nsize -  lenght of the kernel 
 * @return array of length 'nsize' containing the gaussian samples
 */
float* gaussian_kernel(float sigma, int *nsize)
{
    
    /** Compute a Gaussian kernel of length 'n',
     standard deviation 'sigma', and centered at value zero.
     
     For example, if mean=0.5, the Gaussian will be centered
     in the middle point between values 'values[0]'
     and 'values[1]'.
     
     
     */
    
    float sum = 0.0;
    float *kernel;
    float val;
    int n, i, h;
    float prec;
    
    
    /* check parameters */
    if( sigma <= 0.0 ) error("gaussian_kernel: 'sigma' must be positive.");
    
    /*
     The size of the kernel is selected to guarantee that the
     the first discarded term is at least 10^prec times smaller
     than the central value. For that, h should be larger than x, with
     e^(-x^2/2sigma^2) = 1/10^prec.
     Then,
     x = sigma * sqrt( 2 * prec * ln(10) ).
     */
    prec = 3.0;
    h = (int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
    n = 1+2*h; /* kernel size */
    
    kernel = (float*) malloc( n * sizeof(float));
    
    /* compute Gaussian kernel */
    for(i=0;i<n;i++)
    {
        val = ( (float) i - h ) / sigma;
        kernel[i] = exp( -0.5 * val * val );
        sum += kernel[i];
    }
    
    /* normalization */
    if( sum >= 0.0 ) for(i=0;i<n;i++) kernel[i] /= sum;
    
    *nsize = n;
    return kernel;
}


/* homemade round function: float to int*/
/** 
 *  @brief Homemeda round function: float to closest int.
 *  @param x  - input float
 *  @return closest int [x] to x 
 */
 int roundfi(float x) {
    if ((x <= INT_MIN-0.5) || (x >= INT_MAX+0.5))
        error("roundfi() Float to int conversion out of range");
    if (x >= 0)
        return (int) (x+0.5);
    return (int) (x-0.5);
}




/** 
 *  @brief Low-Pass filter Image (with DCT Transform).
 *  Put to zero all frequencies higher than fcx or fcy 
 *  (Let c_ij = DCT(in) then c_ij = 0 for j>=fcx or i>=fcy)
 *  @param fcx   - Cut-off frequency in x
 *  @param fcy   - Cut-off frequency in y  
 *  @return Low-pass filtered image
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
    
    /*Normalize image because DCT introduces
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
 * @brief Normalize an array of floats to sum_i h[i] = 1 
 * @param h   - Array of floats 
 * @param n   - Size of the array 
 */
void normalize_area(float* h, int n)
{

    int i;
    float acsum = 0;
    
    for (i = 0; i < n; i++)
        acsum = acsum + h[i];
    
    for (i = 0; i < n; i++)
        h[i] = h[i]/acsum;
    
}


/**
 * @brief Evaluate a Homography in an array 'np' points 'Pin' and
 *        return the result in the array 'Pout' 
 * @param H   - Homography 
 * @param Pin - Array of input points 
 * @param Pou - Array of output points
 * @param np  - number of points 
 */
void evaluate_homography(float* H, float *Pin, float *Pout, int np)
{
    int i;
    float xPo, yPo, zPo, xP, yP;
    
    
    for (i = 0; i < np; i++) {
        /*Read the point */
        xP = Pin[2 * i];
        yP = Pin[2 * i + 1];
        
        /*Apply the homography to (xP,yP) */
        xPo = H[0] * xP + H[1] * yP + H[2];
        yPo = H[3] * xP + H[4] * yP + H[5];
        zPo = H[6] * xP + H[7] * yP + H[8];
        
        Pout[2 * i] = xPo/zPo;
        Pout[2 * i + 1] = yPo/zPo;
    }
    
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


