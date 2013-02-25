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
 * @file image.h
 * @brief library header for basic image processing.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */

#ifndef IMAGE_H_
#define IMAGE_H_



/*---------------------------------------------------------------------------*/
/** @brief double image val type
 
 The pixel value at (x,y) is accessed by:
 
 image->val[ x + y * image->ncol ]
 
 with x and y integer.
 */
typedef struct imageFloatStruct
{
    float *val;
    int ncol, nrow;
} *ImageFloat;

void free_imageFloat(ImageFloat i);
ImageFloat new_imageFloat(int ncol, int nrow);


/*---------------------------------------------------------------------------*/


void error(char *msg);


ImageFloat convol_sep2(ImageFloat in, float *xker, int xsize, float *yker,
                       int ysize);

ImageFloat convol(ImageFloat in, ImageFloat kernel);


ImageFloat extract_window(ImageFloat in, int xmin, int xmax,
                          int ymin, int ymax);

ImageFloat bicubic(ImageFloat X, ImageFloat Y, ImageFloat in, float a);


ImageFloat compute_dct_image(ImageFloat img);
ImageFloat compute_idct_image(ImageFloat img);

float* gaussian_kernel(float sigma, int *nsize);

ImageFloat lpf_image_dct (ImageFloat in, int fcx, int fcy);

int roundfi(float x);

void normalize_area(float* h, int n);

void evaluate_homography(float* H, float *Pin, float *Pout, int np);

void write_ascii_matrix(float *M, int ncol, int nrow, char *name);

void write_ascii_imageFloat (ImageFloat image, char *name);


#endif          /* IMAGE_H_ */

