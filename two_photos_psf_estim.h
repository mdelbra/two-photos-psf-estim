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

/*
 * two_photos_psf_estim.h
 *
 *  Created on: Oct 22, 2010
 *      Author: mdelbra
 */

#ifndef TWO_PHOTOS_PSF_ESTIM_H_
#define TWO_PHOTOS_PSF_ESTIM_H_


void two_photos_psf_estim (float *img1, int nx1, int ny1,
                           float *img2, int nx2, int ny2,
                           int s, int psf_nrow, int psf_ncol,
                           float *h, float *k, 
                           int threshold, char *outprefix);


#endif				/* TWO_PHOTOS_PSF_ESTIM_H_ */
