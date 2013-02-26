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
 * @file two_photos_psf_estim_main.c
 * @brief main for psf estimation algorithm execution.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "io_pgm.h"
#include "image.h"
#include "two_photos_psf_estim.h"


/** @brief struct of program parameters */
typedef struct
{
	int s;
	int psf_nx;
	int psf_ny;
	char *img_zin;
	char *out_psf_pgm;
	char *out_intk_pgm;
	char *img_zout;
	char *output_psf;
	char *output_intk;
	char *outprefix;
	int threshold;
} program_argums;



static void usage(const char* name)
{
	printf("Recovering the Subpixel PSF from Two Photographs ");
    printf("at Different Distances\n");
	printf("Copyright M.Delbracio, A.Almansa, P.Muse. ");
    printf("Version 1.0 Feb 22, 2013\n\n");
	printf("Usage: %s [options] <input file 1> <input file 2> <output PSF txt> ", name);
    printf("<output inter-kernel txt>\n\n");
    printf("Only  PGM 16/8 bits images are supported.\n\n");
	printf("Options:\n");
	printf("   -s <number>   Super-resolution factor ");
    printf("(positive integer, default 3)\n");
	printf("   -k <number>   PSF support size (default 13)\n");
	printf("   -o <file>     Estimated PSF save to a 8bits PGM image \n");
	printf("   -i <file>     Estimated inter-image kernel save ");
    printf("to a 8bits PGM image \n");
	printf("   -d <file>     Save all the intermediate images in ");
    printf("files with prefix <file>\n");
	printf("   -t 0,1        1 - Threshold negative values to zero (0-default)\n");
}

static int parse_arguments(program_argums *param, int argc, char *argv[])
{
	char *OptionString;
	char OptionChar;
	int i;
	
	
	if(argc < 5)
	{
		usage(argv[0]);
		return 0;
	}
	
	/* Set parameter defaults */	
	param->s = 3;
	param->psf_nx = 13;
	param->psf_ny = 13;
	param->threshold = 0;
	param->img_zin = NULL;
	param->img_zout = NULL;
	param->output_psf = NULL;
	param->output_intk = NULL;
	param->outprefix = NULL;
	param->out_psf_pgm = NULL;
	param->out_intk_pgm = NULL;
	
	
	for(i = 1; i < argc;)
	{
		if(argv[i] && argv[i][0] == '-')
		{
            /*Read options that do not need extra arguments*/
			if((OptionChar = argv[i][1]) == 0)
			{
				printf("Invalid parameter format.\n");
				return 0;
			}
            
            
            if(argv[i][2])
                OptionString = &argv[i][2];
            else if(++i < argc)
                OptionString = argv[i];
            else
            {
                printf("Invalid parameter format.\n");
                return 0;
            }
            
            switch(OptionChar)
            {
                case 's':
                    param->s = atoi(OptionString);
                    if(param->s < 1)
                    {
                        printf("Invalid superresolution factor.\n");
                        return 0;
                    }
                    break;
                case 'k':
                    param->psf_nx = atoi(OptionString);
                    param->psf_ny = atoi(OptionString);
                    if(param->psf_nx <= 0)
                    {
                        printf("Invalid PSF support size.\n");
                        return 0;
                    }
                    break;
              
                case 't':
                    param->threshold = atoi(OptionString);
                    if((param->threshold != 0) && (param->threshold != 1))
                    {
                        printf("Invalid Threshold option.\n");
                        return 0;
                    }
                    break;
                    
                case 'o':
                    param->out_psf_pgm = OptionString;
                    break;
                    
                case 'i':
                    param->out_intk_pgm = OptionString;
                    break;
                    
                case 'd':
                    param->outprefix = OptionString;
                    break;
                    
                case '-':
                    usage(argv[0]);
                    return 0;
                default:
                    if(isprint(OptionChar))
                        printf("Unknown option \"-%c\".\n", OptionChar);
                    else
                        printf("Unknown option.\n");
                    
                    return 0;
            }
            
            
            i++;
        }
        else
        {
            if(!param->img_zin)
                param->img_zin = argv[i];
            else if(!param->img_zout)
                param->img_zout = argv[i];
            else if(!param->output_psf) 
                param->output_psf = argv[i];
            else
                param->output_intk = argv[i];
            i++;
        }
    }
    
    if(!param->img_zin || !param->img_zout ||
       !param->output_psf || !param->output_intk)
    {
        usage(argv[0]);
        return 0;
    }
    
    
    printf("\n");
    printf("   Loaded arguments: \n");
    printf("   ----------------- \n");
    printf("         Superresolution        s : %dx\n",
           param->s);
    printf("         PSF Support            k : %dx%d\n",
           param->psf_nx,param->psf_ny);
    printf("         PSF Output (PGM)       o : %s\n",
           param->out_psf_pgm?param->out_psf_pgm:"(no)");
    printf("         int-k Output (PGM)     i : %s\n",
           param->out_intk_pgm?param->out_intk_pgm:"(no)");
    printf("         Threshold Non-negative t : %s\n",
           param->threshold?"(yes)":"(no)");
    printf("         img1 (PGM)               : %s\n",
           param->img_zin);
    printf("         img2 (PGM)               : %s\n",
           param->img_zout);
    printf("         PSF Output (TXT)         : %s\n",
           param->output_psf);
    printf("         int-kernel Output (TXT)  : %s\n",
           param->output_intk);
    printf("         Debug Images (PGM)       : %s\n\n",
           param->outprefix?param->outprefix:"(no)");
    
    return 1;
}


/**
 * @brief main function call
 */
int main(int argc, char *argv[])
{
    int nx_zin,ny_zin, nx_zout, ny_zout;
    
    float *img_zin, *img_zout, *h, *k;
    program_argums param;
    
    
    /*Parse command-line arguments*/
    if(!parse_arguments(&param,argc,argv))
    {
        
        return EXIT_FAILURE;
    
    }
    /* read the first PGM image into data */
    if (NULL == (img_zout = read_pgm_float(param.img_zout, &nx_zout,
                                           &ny_zout))) {
        fprintf(stderr, 
                "first image could not be properly read\n");
        return EXIT_FAILURE;
    }
    
    
    /* read the second PGM image into data */
    if (NULL ==
        (img_zin = read_pgm_float(param.img_zin, &nx_zin, &ny_zin)))
    {
        fprintf(stderr, 
                "second image could not be properly read\n");
        return EXIT_FAILURE;
    }
    
    h = (float*)malloc(param.psf_nx * param.psf_ny * sizeof(float));
    k = (float*)malloc(param.psf_nx * param.psf_ny * sizeof(float));
    
    /* Call psf_estimation */
    two_photos_psf_estim(img_zout,  
                         nx_zout, 
                         ny_zout,
                         img_zin, 
                         nx_zin, 
                         ny_zin,
                         param.s, 
                         param.psf_nx, 
                         param.psf_ny, 
                         h,
                         k, 
                         param.threshold,
                         param.outprefix);
    
    
    
    /*Write the estimated PSF to a text file */
    write_ascii_matrix(h, param.psf_nx, 
                       param.psf_ny, param.output_psf);
    
    /*Write the estimated inter kernel to a text file */
    write_ascii_matrix(k, param.psf_nx,
                       param.psf_ny, param.output_intk);
    
	
    /* Write the estimated PSF to a 8bit-PGM image file if 
     required
     */
    if(param.out_psf_pgm)
    {
        /* Normalize and truncate negative values */
        write_pgm_normalize_float(param.out_psf_pgm,
                                  h, param.psf_nx, param.psf_ny);
    }
    
    /* Write the estimated inter image kernel to a 8bit-PGM 
     image file if required
     */
    if(param.out_intk_pgm)
    {
        /* Normalize and truncate negative values */
        write_pgm_normalize_float(param.out_intk_pgm,
                                  k, param.psf_nx, param.psf_ny);
    }
    
    /* do the cleanning */
    free(h);
    free(k);
    free(img_zin);
    free(img_zout);
    
    return EXIT_SUCCESS;
}

