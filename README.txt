Recovering the Subpixel PSF from Two Photographs at Different Distances
======================================================================
Version 1.0 - February 22, 2013

by    Mauricio Delbracio <mdelbra@gmail.com>
      Andres Almansa
      Pablo Muse


Introduction
-----------
In most typical digital cameras, even high-end digital single lens reflex 
cameras (DSLR), the acquired images are sampled at rates below the Nyquist
critical rate, causing aliasing effects. In this work we describe a new
algorithm presented in [1] for the estimation of the point spread function
(PSF) of a digital camera from aliased photographs, that achieves subpixel 
accuracy. The procedure is based on taking two parallel photographs of the 
same scene, from different distances leading to different geometric zooms,
and then estimating the kernel blur between them. 

[1] "Subpixel Point Spread Function Estimation from Two Photographs at 
 Different Distances"  M. Delbracio, A. Almansa, J.-M. Morel and P. Muse
 SIAM Journal on Imaging Sciences (SIIMS), November 2012.
 DOI: 10.1137/110848335


 Files
 -----
 Makefile
 doxygen.config
 VERSION
 COPYING
 README.txt
 ls.c
 ls.h
 image.c
 image.h
 io_pgm.c
 io_phm.h
 two_photos_psf_estim.c
 two_photos_psf_estim.h
 two_photos_estim_main.c
 third_party/lib_orsa_homography.cpp
 third_party/OrsaHomography_20120515 [IPOL published]

Requirements
------------
- The fftw3 header and libraries are required on the system for
compilation and execution. See http://www.fftw.org/

- The cblas header and libraries are required on the system for
compilation and execution.

- The lapack library is required on the system for
compilation and execution.

- HOMOGRAPHY + ORSA (an implementation form IPOL [2] is given with the proposed code)


[2] "Automatic Homographic Registration of a Pair of Images, with A 
 Contrario Elimination of Outliers" L. Moisan, P. Moulon, P. Monasse
 Image Processing On Line, 2012 
 DOI: 10.5201/ipol.2012.mmm-oh


Compilation
-----------
Simply use the provided makefile, with the command `make`. You need to set
the directory where the libraries: ffw3, cblas and lapack have the respective
header and libraries files.

Running
-------

Usage: ./two_photos_psf_estim [options] <input file 1> <input file 2> 
<outputPSF txt> <output inter-kernel txt>

Only  PGM 16/8 bits images are supported.

Options:
   -s <number>   The super-resolution factor, positive integer  (default 3)
   -k <number>   PSF support size (default 13)
   -o <file>     Estimated PSF save to a 8bits PGM image
   -i <file>     Estimated inter-image kernel save to a 8bits PGM image
   -d <file>     Save all the intermediate images in files with prefix <file>
   -t            1 - Threshold negative values to zero (0-default)

Parameter Explanation

-s 'number' : The superresolution factor, i.e. how many additional samples
per observed pixel will be computed. (default 3)

-k 'number' : The support size in the superresolved grid. For very sharp
images 4s+1 should be enough. (default 13)

-o 'filename' : Save the estimated psf into a PGM image file 'filename'. For
visualization purposes.

-i 'filename' : Save the estimated inter-image kernel into a PGM image file
'filename'. For visualization purposes.

-d 'filename' : Save all the intermediate images in files with prefix <file>

-t            : Threshold negative values to zero

Documentation
-------------
A detailed documentation is in (IPOL).


Please report bugs in two_scales_psf_estim to <mdelbra@gmail.com>.



Copyright and License
---------------------
 
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

