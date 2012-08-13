
/**
 * @file image.h
 * @brief library header for basic image processing.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */

#ifndef IMAGE_H_
#define IMAGE_H_




/*----------------------------------------------------------------------------*/
/** @brief double image val type

    The pixel value at (x,y) is accessed by:

      image->val[ x + y * image->ncol ]

    with x and y integer.
 */
typedef struct imageFloatStruct {
    float *val;
    int ncol, nrow;
} *ImageFloat;

void free_imageFloat(ImageFloat i);
ImageFloat new_imageFloat(int ncol, int nrow);


/*----------------------------------------------------------------------------*/




float dist_l2(float x1, float y1, float x2, float y2);

ImageFloat convol_sep2(ImageFloat in, float *xker, int xsize, float *yker,
		       int ysize);

ImageFloat convol(ImageFloat in, ImageFloat kernel);


ImageFloat extract_window(ImageFloat in, int xmin, int xmax,
			  int ymin, int ymax);

ImageFloat extract_subpx_window(ImageFloat in, int wsize, float cx,
				float cy);

ImageFloat gradx(ImageFloat in);

ImageFloat grady(ImageFloat in);

ImageFloat bilinear(ImageFloat X, ImageFloat Y, ImageFloat in);

ImageFloat bicubic(ImageFloat X, ImageFloat Y, ImageFloat in, float a);

float mean_window(ImageFloat in, int xmin, int xmax, int ymin, int ymax);

float power_window(ImageFloat in, int xmin, int xmax, int ymin, int ymax);


float mean_subpx_window(ImageFloat in, int wsize, float cx, float cy);
float power_subpx_window(ImageFloat in, int wsize, float cx, float cy);

ImageFloat compute_dct_image(ImageFloat img);
ImageFloat compute_idct_image(ImageFloat img);

float* gaussian_kernel(float sigma, int *nsize);


#endif				/* IMAGE_H_ */
