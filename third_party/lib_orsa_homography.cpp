/*----------------------------------------------------------------------------
 *The following function for estimating an Homography between two images
 *is an adaptation of the <orsa_homography.cpp> file of Monasse/Moulon 
 *(ORSA-Homography) [1]. This funcion is only a wrapper to use the software 
 * published in IPOL as a llibrary.
 *
 * [1]
 * "Automatic Homographic Registration of a Pair of Images, with A Contrario 
 * Elimination of Outliers"
 * Moisan, Lionel and Moulon, Pierre and Monasse, Pascal
 * Image Processing On Line, 2012
 * DOI, 10.5201/ipol.2012.mmm-oh
 *
 *
 * by M.Delbracio 22 Feb 2013.
 */


/**
 * @file orsa_homography.cpp
 * @brief Homographic image registration
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2012 Pascal Monasse, Pierre Moulon
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <ctime>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "libImage/image_io.hpp"
#include "libImage/image_drawing.hpp"
#include "libImage/image_crop.hpp"

#include "libOrsa/homography_model.hpp"

#include "extras/libNumerics/numerics.h"
#include "extras/sift/library.h"

#include "demo/siftMatch.hpp"
#include "demo/warping.hpp"
//#include "demo/cmdLine.h"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

#define BIG_NUMBER 1000000


/// Display average/max error of inliers of homography H.
static void display_stats(const std::vector<Match>& vec_matchings,
                          const std::vector<int>& vec_inliers,
                          libNumerics::matrix<double>& H) {
  std::vector<int>::const_iterator it=vec_inliers.begin();
  double l2=0, linf=0;
  for(; it!=vec_inliers.end(); ++it) {
    const Match& m=vec_matchings[*it];
    double x1=m.x1, y1=m.y1;
    TransformH(H, x1, y1);
    double e = (m.x2-x1)*(m.x2-x1) + (m.y2-y1)*(m.y2-y1);
    l2 += e;
    if(linf < e)
      linf = e;
  }
  std::cout << "Average/max error: "
            << sqrt(l2/vec_inliers.size()) << "/"
            << sqrt(linf) <<std::endl;
}

/// ORSA homography estimation
bool ORSA(const std::vector<Match>& vec_matchings, int w1,int h1, int w2,int h2,
          double precision,
          libNumerics::matrix<double>& H, std::vector<int>& vec_inliers)
{
  const int n = static_cast<int>( vec_matchings.size() );
  if(n < 5)
  {
      std::cerr << "Error: ORSA needs 5 matches or more to proceed" <<std::endl;
      return false;
  }
  libNumerics::matrix<double> xA(2,n), xB(2,n);

  for (int i=0; i < n; ++i)
  {
    xA(0,i) = vec_matchings[i].x1;
    xA(1,i) = vec_matchings[i].y1;
    xB(0,i) = vec_matchings[i].x2;
    xB(1,i) = vec_matchings[i].y2;
  }

  orsa::HomographyModel model(xA, w1, h1, xB, w2, h2, true);
  //model.setConvergenceCheck(true);

  if(model.orsa(vec_inliers, ITER_ORSA, &precision, &H, true)>0.0)
    return false;
  std::cout << "Before refinement: ";
  display_stats(vec_matchings, vec_inliers, H);
  if( model.ComputeModel(vec_inliers,&H) ) // Re-estimate with all inliers
  {
    std::cout << "After  refinement: ";
    display_stats(vec_matchings, vec_inliers, H);
  } else
    std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
  return true;
}

/// Output inlier and oulier matches in image files.
void display_match(const std::vector<Match>& vec_matchings,
                   std::vector<int>& vec_inliers,
                   const libNumerics::matrix<double>* H,
                   const libNumerics::matrix<double>& H1,
                   const libNumerics::matrix<double>& H2,
                   Image<RGBColor>& in, Image<RGBColor>& out)
{
  std::sort(vec_inliers.begin(), vec_inliers.end());

  // For outliers, show vector (yellow) from prediction to observation
  const RGBColor col=YELLOW;
  std::vector<int>::const_iterator it = vec_inliers.begin();
  matchingslist::const_iterator m = vec_matchings.begin();
  if(H) // Otherwise, no prediction
    for(int i=0; m != vec_matchings.end(); ++m, ++i) {
      if(it != vec_inliers.end() && i==*it)
        ++it;
      else { //Outlier
          double x1=m->x1, y1=m->y1, x2=m->x2, y2=m->y2;
          TransformH(H2 * *H, x1, y1);
          TransformH(H2, x2, y2);
          libs::DrawLine((int)x1,(int)y1,(int)x2,(int)y2, col,&out);
      }
    }

  // Show link for inliers (green) and outliers (red)
  it = vec_inliers.begin();
  m = vec_matchings.begin();
  for(int i=0; m != vec_matchings.end(); ++m, ++i)
  {
    Image<RGBColor>* im=&out;
    RGBColor col=RED;
    if(it != vec_inliers.end() && i==*it) {
      ++it;
      im=&in;
      col=GREEN;
    }
    double x1=m->x1, y1=m->y1, x2=m->x2, y2=m->y2;
    TransformH(H1, x1, y1);
    TransformH(H2, x2, y2);
    libs::DrawLine((int)x1,(int)y1,(int)x2, (int)y2, col, im);
  }
}

/// Return 3x3 translation matrix
libNumerics::matrix<double> translation(double dx, double dy) {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}

/// Return 3x3 zoom-translation matrix
libNumerics::matrix<double> zoomtrans(double z, double dx, double dy) {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,0)=T(1,1) = z;
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}



extern "C" int orsa_homography_sift(float* img1, int nx1, int ny1,
									float* img2, int nx2, int ny2,
									double precision, float *hom,
									float **matchsHR, float **matchsLR, int* nm,
									char *reverse)
{
	
	
	/*BEGIN Normalization input images for SIFT (0,255)*/
    /* check min and max values img1*/
    float max_val1, max_val2, min_val1, min_val2, max_val, min_val; 
    float v;
    int i;
   
    max_val1 = max_val2 = 0;
    min_val1 = min_val2 = BIG_NUMBER;
    
    for(i=0;i<nx1*ny1;i++)
    {
        v = img1[i];
        if( v > max_val1 ) max_val1 = v;
        if( v < min_val1 ) min_val1 = v;
    }
    
    
    /* check min and max values img2*/
    for(i=0;i<nx2*ny2;i++)
    {
        v = img2[i];
        if( v > max_val2 ) max_val2 = v;
        if( v < min_val2 ) min_val2 = v;
    }
    
    
    /*normalize both images with the same constants*/
    max_val = (max_val2<max_val1)?max_val2:max_val1;
    min_val = (min_val2>min_val1)?min_val2:min_val1;
    
    printf("Image 1: min: %f, max: %f\n", min_val1, max_val1);
    printf("Image 2: min: %f, max: %f\n", min_val2, max_val2);
    printf("Normal:  min: %f, max: %f\n", min_val, max_val);
    
    for(i=0;i<nx1*ny1;i++)
        img1[i] = 255*(img1[i] - min_val)/(max_val-min_val);
    
    for(i=0;i<nx2*ny2;i++)
        img2[i] = 255*(img2[i] - min_val)/(max_val-min_val);
    
    
    /*END Normalization*/
    
	Image<float> image1(nx1, ny1, img1);
	Image<float> image2(nx2, ny2, img2);
	
	
	Image<unsigned char> image1Gray, image2Gray;
	libs::convertImage(image1, &image1Gray);
	libs::convertImage(image2, &image2Gray);
	
	const int w1=image1Gray.Width(), h1=image1Gray.Height();
	const int w2=image2Gray.Width(), h2=image2Gray.Height();
	
    
    
	// SIFT
	std::vector<Match> vec_matchings;
	SIFT(image1Gray, image2Gray, vec_matchings);
	
	// Remove duplicates (frequent with SIFT)
	rm_duplicates(vec_matchings);
	
	
	
	/*Matcher before ORSA consistency check!
     int np;
     float *matchings1;
     float *matchings2;
     
     np = vec_matchings.size();
     
     matchings1 = (float*) malloc( 2*np*sizeof(float));
     matchings2 = (float*) malloc( 2*np*sizeof(float));
     
     
     for (size_t i=0; i < np; ++i)
     {
     matchings1[2*i] = vec_matchings[i].x1;
     matchings1[2*i+1] = vec_matchings[i].y1;
     matchings2[2*i] = vec_matchings[i].x2;
     matchings2[2*i+1] = vec_matchings[i].y2;
     
     }
     
     
     *matchs1 = matchings1;
     *matchs2 = matchings2;
     *nm = np;
     */
    
	
	// Estimation of homography with ORSA
	/*Assume order ir correct*/
	*reverse = 0;
	
	libNumerics::matrix<double> H(3,3);
	std::vector<int> vec_inliers;
	bool ok = ORSA(vec_matchings, w1, h1, w2, h2, precision, H, vec_inliers);
	if(ok)
	{
		H /= H(2,2);
		//std::cout << "H=" << H <<std::endl;
	}
	
	/*Check if H(1,1) and H(2,2) are greater than 1*/
	if(H(0,0)<1 && H(1,1)<1)
	{
        
		*reverse = 1;
		
		/*change the order of the matches img1 --> img2 and viceverza*/
		for (size_t i=0; i < vec_matchings.size(); ++i)
		{
			float aux;
			
			aux = vec_matchings[i].x1;
			vec_matchings[i].x1 = vec_matchings[i].x2;
			vec_matchings[i].x2 = aux;
			
			aux = vec_matchings[i].y1;
			vec_matchings[i].y1 = vec_matchings[i].y2;
			vec_matchings[i].y2 = aux;
			
		}
		
		ok = ORSA(vec_matchings, w2, h2, w1, h1, precision, H, vec_inliers);
        
		if(ok)
		{
			H /= H(2,2);
			//std::cout << "H=" << H <<std::endl;
		}
		
	}
	
	
	std::vector<Match> good_match;
	
	std::vector<int>::const_iterator it = vec_inliers.begin();
	for(; it != vec_inliers.end(); it++)
		good_match.push_back(vec_matchings[*it]);
	
	
	/*Matcher before ORSA consistency check!*/
	if(matchsHR!=NULL && matchsLR!=NULL){
		
		int np;
		float *matchingsHR;
		float *matchingsLR;
		
		np = good_match.size();
		
		matchingsHR = (float*) malloc( 2*np*sizeof(float));
		matchingsLR = (float*) malloc( 2*np*sizeof(float));
		
		
		for (int i=0; i < np; ++i)
		{
			matchingsHR[2*i] = good_match[i].x1;
			matchingsHR[2*i+1] = good_match[i].y1;
			matchingsLR[2*i] = good_match[i].x2;
			matchingsLR[2*i+1] = good_match[i].y2;
			
		}
		
		*matchsHR = matchingsHR;
		*matchsLR = matchingsLR;
		*nm = np;
		
		
	}
	
	/*Copy homography*/
	if(hom!=NULL)
	{
		hom[0] = H(0,0);
		hom[1] = H(0,1);
		hom[2] = H(0,2);
		hom[3] = H(1,0);
		hom[4] = H(1,1);
		hom[5] = H(1,2);
		hom[6] = H(2,0);
		hom[7] = H(2,1);
		hom[8] = H(2,2);
	}
	else{
		std::cerr << "Failed to return the homography parameters" << std::endl;
		return 1;
	}
	
    
    
	if(!ok)
	{
		std::cerr << "Failed to estimate a model" << std::endl;
		return 1;
	}
	
	return 0;
}


