//Copyright (C) 2011 Pascal Monasse, Pierre Moulon
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "libImage/image_io.hpp"
#include "libImage/image_drawing.hpp"
#include "libImage/image_concat.hpp"

#include "libOrsa/homography_model.hpp"

#include "demo/warping.hpp"

#include "extras/libNumerics/numerics.h"
#include "extras/sift/library.h"
#include "extras/sift/demo_lib_sift.h"

#include <time.h>

/// SIFT matches
static void SIFT(const Image<unsigned char> &im1,
                 const Image<unsigned char> &im2,
                 std::vector<Match>& vec_matchings) {
	//Convert images to float
	Image<float> If1, If2;
	libs::convertImage(im1, &If1);
	libs::convertImage(im2, &If2);
	
	siftPar param;
	default_sift_parameters(param);
	param.DoubleImSize=0;
	
	keypointslist keyp1, keyp2;
	compute_sift_keypoints(If1.data(), keyp1, If1.Width(), If1.Height(), param);
	std::cout<< "sift:: 1st image: " << keyp1.size() << " keypoints"<<std::endl;
	compute_sift_keypoints(If2.data(), keyp2, If2.Width(), If2.Height(), param);
	std::cout<< "sift:: 2nd image: " << keyp2.size() << " keypoints"<<std::endl;
	
	// Find Putatives matches
	compute_sift_matches(keyp1, keyp2, vec_matchings, param);
	std::cout << "sift:: matches: " << vec_matchings.size() <<std::endl;
}

/// Remove multiple "same position" matches
static void rm_duplicates(std::vector<Match>& m) {
	std::sort(m.begin(), m.end());
	std::vector<Match>::iterator end = std::unique(m.begin(), m.end());
	if(end != m.end()) {
		std::cout << "Remove " << std::distance(end, m.end())
		<< "/" << m.size() << " duplicate matches, "
		<< " keeping " << std::distance(m.begin(), end) <<std::endl;
		m.erase(end, m.end());
	}
}

/// Display average/max error of inliers of homography H.
static void display_stats(const std::vector<Match>& vec_matchings,
                          const std::vector<size_t>& vec_inliers,
                          libNumerics::matrix<double>& H) {
	std::vector<size_t>::const_iterator it=vec_inliers.begin();
	double l2=0, linf=0;
	for(; it!=vec_inliers.end(); ++it) {
		const Match& m=vec_matchings[*it];
		libNumerics::vector<double> x(3);
		x(0)=m.x1;
		x(1)=m.y1;
		x(2)=1.0;
		x = H*x;
		x /= x(2);
		double e = (m.x2-x(0))*(m.x2-x(0)) + (m.y2-x(1))*(m.y2-x(1));
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
          libNumerics::matrix<double>& H, std::vector<size_t>& vec_inliers)
{
	const size_t n = vec_matchings.size();
	if(n < 5)
	{
		std::cerr << "Error: ORSA needs 5 matches or more to proceed" <<std::endl;
		return false;
	}
	libNumerics::matrix<double> xA(2,n), xB(2,n);
	
	for (size_t i=0; i < n; ++i)
	{
		xA(0,i) = vec_matchings[i].x1;
		xA(1,i) = vec_matchings[i].y1;
		xB(0,i) = vec_matchings[i].x2;
		xB(1,i) = vec_matchings[i].y2;
	}
	
	orsa::HomographyModel model(xA, w1, h1, xB, w2, h2);
	
	if(model.orsa(vec_inliers, 1000, &precision, &H, true)>0.0)
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







extern "C" int orsa_homography_sift(float* img1, int nx1, int ny1,
									float* img2, int nx2, int ny2,
									double precision, float *hom,
									float **matchsHR, float **matchsLR, int* nm,
									char *reverse)
{
	
	
	

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
	std::vector<size_t> vec_inliers;
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
	
	std::vector<size_t>::const_iterator it = vec_inliers.begin();
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

