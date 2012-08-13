

/**
 * @file io_pgm.h
 * @brief library header for basic image processing.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */

#ifndef IO_PGM_H_
#define IO_PGM_H_


float *read_pgm_float(const char * fname, int * ncol, int *nrow);

void write_pgm_float(const char *fname, const float *data, int ncol, int nrow);

void write_pgm_normalize_float(const char *fname, const float *data, int ncol, int nrow);

void write_pgm_normalize_given_minmax_float(const char *fname, 
									   const float *data, int ncol,
									   int nrow, float min,
											float max);

void write_pgm_lognormalize_float(const char *fname, const float *data, int ncol, int nrow);

void write_ppm_normalize_float(const char *fname, const float *rdata, const float *gdata,
							   const float *bdata, int ncol, int nrow);

	
/*int write_pgm_float(const char * fname, const float *data, int nx,
				  int ny);*/

#endif  