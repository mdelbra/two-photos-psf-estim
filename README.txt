Recovering the Subpixel PSF from Two Photographs at Different Distances
======================================================================
Version Beta 0.4 - August 13, 2012

by    Mauricio Delbracio <mdelbra@gmail.com>
      Andres Almansa
      Pablo Muse


Introduction
-----------
Extrinsic image blur can be observed when the camera's focal distance is
not correctly adjusted by the user, when the objects in the scene appear at
different depths, or when the relative motion between the camera and the scene
is faster than the shutter speed (motion blur). Besides these sources of blur,
even in ideal acquisition conditions there is a permanent intrinsic physical
camera blur due to light diffraction, sensor resolution, lens aberration,
and anti-aliasing filters. Our goal here is to accurately estimate the Point
Spread Function - PSF, that models the intrinsic camera blur. This function can
be locally interpreted as the response of the camera to a point light source.

In [1] we presented a theoretical study proving that the sub-pixel PSF
estimation problem is well-posed even with a single image capture, as long
as the captured scene is well chosen. Indeed, theoretical bounds show that
a near-optimal accuracy can be achieved by taking a single snapshot of a
calibration pattern mimicking a Bernoulli(0.5) white noise. We first use an
algorithm to accurately estimate the pattern position and its illumination
conditions. This allows for accurate geometric registration and radiometric
correction; Once these procedures have been applied, the local PSF can
be directly computed by inverting a linear system that is well-posed and
consequently its inversion does not require any regularization or prior model.


 Files
 -----
 Makefile
 README.txt
 ls.c
 ls.h
 homography.c
 homography.h
 image.c
 image.h
 io_pgm.c
 io_phm.h
 detect_pattern.c
 detect_pattern.h
 two_photos_psf_estim.c
 two_photos_psf_estim.h
 two_photos_estim_main.c

Requirements
------------
- The fftw3 header and libraries are required on the system for
compilation and execution. See http://www.fftw.org/

- The cblas header and libraries are required on the system for
compilation and execution.

- The lapack library is required on the system for
compilation and execution.

- HOMOGRAPHY + ORSA (an implementation form IPOL is given with the proposed
code)

Compilation
-----------
Simply use the provided makefile, with the command `make`. You need to set
the directory where the libraries: ffw3, cblas and lapack have the respective
header and libraries files.

Running
-------

Usage: ./two_photos_psf_estim [options] <input file 1> <input file 2> <output
PSF txt> <output inter-kernel txt>

Only  PGM 16/8 bits images are supported.

Options:
   -s <number>   The super-resolution factor, positive integer  (default 3)
   -k <number>   PSF support size (default 13)
   -o <file>     Estimated PSF save to a 8bits PGM image
   -i <file>     Estimated inter-image kernel save to a 8bits PGM image
   -d <file>     Save all the intermediate images in files with prefix <file>
   -t 0,1        LS: 0 (default), LS + thresholding: 1


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

-t <0,1> : Choose the numerical algorithm for solving Least Squares.
		     0 - Least Squares	(default)
		     1 - Least Squares + thresholding


Documentation
-------------
A detailed documentation is in (IPOL).


Please report bugs in psf_estim to <mdelbra@gmail.com>.


Copyright and License
---------------------
