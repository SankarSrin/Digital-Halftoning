/********************************************************************************/
/*                                                                              */
/* BAYER                                                                        */
/*	This function performs the dither array halftoning using Bayer's original   */
/*  8x8 threshold array.                                                        */
/*                                                                              */
/*  Synopsis:                                                                   */
/*  	Y = bayer(X)                                                            */
/*				X = continuous-tone image of type double                        */
/*				Y = binary halftone of type uint8                               */
/*                                                                              */
/* Daniel Leo Lau									                            */
/* Copyright April 23, 2004 								                    */
/*											                                    */
/********************************************************************************/ 
#include <math.h>
#include "mex.h"
#define mask_row 8
#define mask_col 8

void bayer(unsigned char *output_data, double *input_image, int image_row, int image_col)
{	
    int m_p, n_p, m, n;
	double bayer_mask[64]={0.0234, 0.5234, 0.1484, 0.6484, 0.0547, 0.5547, 0.1797, 0.6797, 
                           0.9297, 0.2734, 0.7734, 0.3984, 0.8984, 0.3047, 0.8047, 0.4297,
                           0.2422, 0.7422, 0.0859, 0.5859, 0.2109, 0.7109, 0.1172, 0.6172,
                           0.8672, 0.4922, 0.9922, 0.3359, 0.8359, 0.4609, 0.9609, 0.3672,
                           0.0391, 0.5391, 0.1641, 0.6641, 0.0078, 0.5078, 0.1328, 0.6328,
                           0.8828, 0.2891, 0.7891, 0.4141, 0.9141, 0.2578, 0.7578, 0.3828,
                           0.1953, 0.6953, 0.1016, 0.6016, 0.2266, 0.7266, 0.0703, 0.5703,
                           0.8203, 0.4453, 0.9453, 0.3516, 0.8516, 0.4766, 0.9766, 0.3203};

	for (m=0; m<image_row; m++){
		for (n=0; n<image_col; n++){
		    if (input_image[m+n*image_row] > bayer_mask[(m%8)+(n%8)*8])
		        output_data[m+n*image_row]=1;
		    else
		        output_data[m+n*image_row]=0;
		    }
		}
	return;
}

/*******************************************************************************/ 
/* mexFUNCTION                                                                 */
/* Gateway routine for use with MATLAB.                                        */
/*******************************************************************************/ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned char *output_data;
	double *input_image;
	int m, n, output_dims[2];
	int image_row, image_col;

    /****** Check for errors in user's call to function. ***********************/
    if (nrhs!=1)
            mexErrMsgTxt("Bayer requires exactly one input argument!");
    else if (nlhs!=1)
            mexErrMsgTxt("Bayer returns exactly one output argument!");
    else if (!mxIsNumeric(prhs[0]) ||
              mxIsComplex(prhs[0]) ||
              mxIsSparse(prhs[0])  ||
             !mxIsDouble(prhs[0]))
            mexErrMsgTxt("Input X must be a real matrix of type double!");

    /****** Get row and column sizes for input image and weight matrix. ********/
    image_row=mxGetM(prhs[0]);
    image_col=mxGetN(prhs[0]);

    /****** Extract pointers to input data. ************************************/
	input_image=mxGetPr(prhs[0]);

    /****** Create output variables and extract pointers to data. **************/
    output_dims[0]=image_row; output_dims[1]=image_col;
    plhs[0]=mxCreateLogicalArray(2, output_dims);
    output_data=mxGetLogicals(plhs[0]);

    /****** Call Bayer function. *************************************/
	bayer(output_data, input_image, image_row, image_col);

	return;
}
