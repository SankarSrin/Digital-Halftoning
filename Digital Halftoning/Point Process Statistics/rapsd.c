/*******************************************************************************/ 
/*
/*                            MATLAB FUNCTION: RAPSD
/*
/* Title:       Radially Averaged Power Spectral Density    
/* Description: From the Power Spectrum it estimates the radially averaged
/*              power spectral density by taking power samples over concentric
/*              annular rings of radius fr where for simplicity incrementsin
/*              fr are taken to be delta=1  
/*              Reference: Robert Ulichney, "Digital Halftoning"
/*
/* MATLAB Prompt Call: [P,fr,count] = RAPSD (X,[M],[FM],[DM DN])
/*
/* X     = Input Image Data (Matrix of type UINT8 (LOGICAL) and dimensions mXn
/* M     = Dimensions of the Sample Window in Spatial Domain
/* FM    = Dimensions of the Sample Window in Frequency Domain
/* DM,DN = Displacement parameters signifying the shift in sample windows
/*
/*******************************************************************************/ 

#include "mex.h"
#include <math.h>

/* Global Variable Declaration */

int    Image_M;
int    Image_N;
int    Spat_Win_Dim;
int    Spec_Win_Dim;
int    Shift_Row, Shift_Col;
int    fr_max;

void   absfft2 (double*,double*,int,int);
void   Compute_Power_Estimate (double*, double*, double*);
void   Rad_Avg_PS (double*, double*, double*, double*);

/*******************************************************************************/ 
/*                  C Gateway routine for use with MATLAB
/*
/*        Accepts the Image Data Matrix, Computes the Power Estimate over
/*        the given Spatial Window size  & returns the center radii fr,
/*        the average radial power contained in each annulus and the
/*                  sample count for Variance calculations
/*
/*     Annulus_Power = sum (Power of Samples in Annulus)/Number of Samples
/*
/*******************************************************************************/ 

void Rad_Avg_PS (double *Image_Data, double *Annulus_Power, double *fr,
double *count)

{
  int i,j;
  double  f1,f2,f3,f4,f_min;
  int     radial_index;
  double  *Power_Estimate;
  double  *Radial_Samples;
  double  g=0.0,variance;
  int     delta=1;

  Power_Estimate = mxCalloc(Spec_Win_Dim*Spec_Win_Dim,sizeof(double)); 
  Radial_Samples = mxCalloc(fr_max,sizeof(double));

  /* Computes the Power Estimate over the given Spectral Window */

  Compute_Power_Estimate(Image_Data,Power_Estimate,count);;

  /* Computes the gray level g of the halftoned image*/

  for (i=0;i<Image_M*Image_N;i++)
    g += Image_Data [i];

  g=g/(Image_M*Image_N);
  variance = g*(1-g);

  /* Creates the Annulii Radius Table */
  
  for (i=0;i<fr_max;i++)
    fr [i] = (0.5 + (double)i)/Spec_Win_Dim;

  /* Start actual computation */

  for (i=0;i<Spec_Win_Dim;i++)
    {
      for (j=0;j<Spec_Win_Dim;j++)
	{
	  f1 = sqrt (i*i+j*j);
	  f_min = f1;
	  f2 = sqrt (i*i+(j-Spec_Win_Dim)*(j-Spec_Win_Dim));
	  if (f2<f_min)
	    f_min = f2;
	  f3 = sqrt((i-Spec_Win_Dim)*(i-Spec_Win_Dim)+(j-Spec_Win_Dim)*(j-Spec_Win_Dim));
	  if (f3<f_min)
	    f_min = f3;
	  f4 = sqrt ((i-Spec_Win_Dim)*(i-Spec_Win_Dim)+j*j);
	  if (f4<f_min)
	    f_min = f4;
	  radial_index=floor (f_min/delta)-1;
	  
	  if (radial_index>=0){
	    Annulus_Power[radial_index] += Power_Estimate[i+j*Spec_Win_Dim];
	    Radial_Samples[radial_index] ++;
	  }
	}
    }
  
  for (i=0;i<fr_max;i++)
    Annulus_Power[i]=Annulus_Power[i]/(Radial_Samples[i]*variance);
  
  mxFree(Radial_Samples);
  return;
}



/*******************************************************************************/ 
/*
/*                   C Gateway routine for use with MATLAB
/*
/*        Accepts the Image Data Matrix & returns the Power_Estimate
/*          as well as the Sample Count for Variance calculations
/*
/*     Power_Estimate = sum (Sample_Power) / No. of Samples (Periodograms)
/*
/*******************************************************************************/ 

void Compute_Power_Estimate (double *Image_Data, double *Power_Estimate,
double *Num_Samples)

{
  int     m=0,n=0;
  int     i,j,count=(int)*(Num_Samples);
  double  *Sample_Window;
  double  *Frequency_Data;
 
  Sample_Window = mxCalloc (Spec_Win_Dim*Spec_Win_Dim,sizeof(double));
  Frequency_Data = mxCalloc (Spec_Win_Dim*Spec_Win_Dim,sizeof(double));

  while ((m+Spat_Win_Dim)<=Image_M)
    {
    n=0;
    while ((n+Spat_Win_Dim)<=Image_N)
      {
      count+=1;

      for (i=0;i<Spat_Win_Dim;i++)
	{
	for (j=0;j<Spat_Win_Dim;j++)
	  {
	  Sample_Window [i+j*Spec_Win_Dim] = Image_Data[(m+i)+(n+j)*Image_M];
	  }
	}
      absfft2 (Frequency_Data,Sample_Window,Spec_Win_Dim,Spec_Win_Dim);

      for (i=0;i<Spec_Win_Dim*Spec_Win_Dim;i++)
	Power_Estimate [i] += Frequency_Data [i];
      
      n=n+Shift_Col;
      }
    
    m=m+Shift_Row;
    } 

  for (i=0;i<Spec_Win_Dim*Spec_Win_Dim;i++)
    Power_Estimate[i] = Power_Estimate[i]/(double)count;

  mxFree (Frequency_Data);
  mxFree (Sample_Window);

  *(Num_Samples)= (double)count;
  return;
}

/*******************************************************************************/ 
/*                            C FUNCTION: absfft2 
/*
/* Title:       Magnitudes Squared of 2D FFT Coefficients    
/* Description: This Function accepts an spatial domain image window of size 
/*              FM X FN, calls the Matlab Routine "fft2" to find the 2-D FFT
/*              of the image, and computes the Magnitude Squared of the
/*              Complex FFT Coefficients. Returns the computed array.
/*
/* Call: absfft2 (X,Y,NR,NC);
/*
/* X  = Output Matrix containing the maginitudes (squared) of the FFT Coeff. 
/* Y  = Zero padded Input Image Matrix (Sample Window), whose FFT is desired
/* NR = Number of Rows in the Image MAtrix
/* NC = Number of Columns in the Image Matrix
/*
/*******************************************************************************/ 

void absfft2 (double *Image_PS, double *Input_Image, int Num_Row, int
Num_Col)
{
  int     i;
  mxArray *array_ptr;
  double  *Real_Data, *Imag_Data;
  mxArray *Input_Array[1], *Output_Array[1];

  array_ptr = mxCreateDoubleMatrix(Num_Row,Num_Col,mxREAL);
  Real_Data = mxGetPr(array_ptr);

  for (i=0; i<Num_Row*Num_Col; i++){
    Real_Data[i]=Input_Image[i];
  }

  Input_Array[0] = array_ptr;
  mexCallMATLAB (1, Output_Array, 1, Input_Array, "fft2");
  if (!mxGetM(Output_Array[0])==Num_Row ||
      !mxGetN(Output_Array[0])==Num_Col)
    mexErrMsgTxt ("Output Array Dimensions are not Correct");

  Real_Data = (double*) mxGetPr(Output_Array[0]);
  Imag_Data = (double*) mxGetPi(Output_Array[0]);

  if (mxIsComplex(Output_Array[0]))
      {
      for (i=0;i<Num_Row*Num_Col;i++)
      Image_PS[i] = (Real_Data[i]*Real_Data[i] +Imag_Data[i]*Imag_Data[i])/(Num_Row*Num_Col);
      }
  else
      {
      for (i=0;i<Num_Row*Num_Col;i++)
      Image_PS[i] = (Real_Data[i]*Real_Data[i])/(Num_Row*Num_Col);
      }
  mxDestroyArray (array_ptr);
  mxDestroyArray (Output_Array[0]);
  return;
}

/*******************************************************************************/ 
/*                                 mexFUNCTION 
/*                      Gateway routine for use with MATLAB.
/*******************************************************************************/ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int    Spat_Win_M, Spat_Win_N;
  int    Spec_Win_M, Spec_Win_N;
  int    Shift_M, Shift_N, i;
  
  double *Spat_Win_Data, *Spec_Win_Data, *Shift_Data;
  double *Image_Data, *Annulus_Power, *fr, *count;
  
  /* Check for errors in User Calls from Matlab pronpt. */
  
  if ((nrhs<2) || (nrhs>4))
    mexErrMsgTxt("This Function requires two to four Input Arguments!");
  
  else if (nlhs<2 || nlhs>3)
    mexErrMsgTxt("This Function returns 2 or 3 Output Argument!");
  
  else if (!mxIsNumeric(prhs[0]) ||
	   mxIsComplex(prhs[0]) ||
	   mxIsSparse(prhs[0])  ||
	   !mxIsDouble(prhs[0])   )
    mexErrMsgTxt("Input Argument X must be a REAL matrix of type DOUBLE!");
  
  else if (!mxIsNumeric(prhs[1]) ||
	   mxIsComplex(prhs[1]) ||
	   mxIsSparse(prhs[1])  ||
	   !mxIsDouble(prhs[1])   )
    mexErrMsgTxt("Input Argumnent [M] must be a REAL scalar of type DOUBLE!");
  
  else if ((nrhs>2) && (!mxIsNumeric(prhs[2]) ||
			mxIsComplex(prhs[2]) ||
			mxIsSparse(prhs[2])  ||
			!mxIsDouble(prhs[2])  ))
    mexErrMsgTxt("Input Argumnent [FM] must be a REAL scalar of type DOUBLE!");
  
  else if ((nrhs>3) && (!mxIsNumeric(prhs[3]) ||
			mxIsComplex(prhs[3]) ||
			mxIsSparse(prhs[3])  ||
			!mxIsDouble(prhs[3])  ))
    mexErrMsgTxt("Input Argumnent [DM DN] must be a REAL vector of type DOUBLE!");

  /* Assigns the pointer of type double to the Image_Data Matrix & */ 	  
  /* Get row & column sizes of the Image Matrix, Vectors [M],[FM] & [DM DN] */

  Image_Data = (double*)mxGetPr (prhs[0]);
  Image_M    = mxGetM  (prhs[0]);
  Image_N    = mxGetN  (prhs[0]);

  /* Checks on Argument 1 */
  
  Spat_Win_M = mxGetM( prhs[1] );
  Spat_Win_N = mxGetN( prhs[1] );
  
  if ((Spat_Win_M*Spat_Win_N)>1) mexErrMsgTxt ("Input Argument N must be scalar!!!");
  else if ((Spat_Win_M*Spat_Win_N)==1) Spat_Win_Dim = (int) mxGetScalar (prhs[1]);
  
  if (Spat_Win_Dim == 0) mexErrMsgTxt ("Spatial Window Dimensions cannot be 0");
  
  /* Set Default values for Arguments 3 & 4 */
  /* Assuming that the other arguments will be zero */          
  
  Spec_Win_Dim = Spat_Win_Dim;
  Shift_Row    = Spat_Win_Dim;
  Shift_Col    = Spat_Win_Dim;

  /* Checks on Argument 2 */
  
  if (nrhs>2) {
        Spec_Win_M = mxGetM  (prhs[2]);
        Spec_Win_N = mxGetN  (prhs[2]); 

        if ((Spec_Win_M*Spec_Win_N)>1) mexErrMsgTxt ("Input Argument 2 can have dimensions 1X1");         
        else if ((Spec_Win_M*Spec_Win_N)==1) Spec_Win_Dim = (int) mxGetScalar (prhs[2]); 

        if (Spec_Win_Dim == 0 || Spec_Win_Dim < 32) Spec_Win_Dim = Spat_Win_Dim;
  }
  
  /* Checks on Argument 3 */
  
  if (nrhs>3)
    {
      Shift_M = mxGetM  (prhs[3]);
      Shift_N = mxGetN  (prhs[3]); 
      
      if ((Shift_M*Shift_N)>2)
	mexErrMsgTxt ("Input Argument 3 can have dimensions 1X1");
      
      else if ((Shift_M*Shift_N)==1)
	{
	  Shift_Row = (int) mxGetScalar (prhs[3]);
	  Shift_Col = Shift_Row;
	} 
      
      if ((Shift_Row*Shift_Col) == 0)
	{
	  Shift_Row = Spat_Win_Dim;
	    Shift_Col = Spat_Win_Dim;
	}
    }

  /* Creates the  Output Variable & Assigns Pointer to the output */
  
  fr_max = floor(Spec_Win_Dim/sqrt(2));
  
  if (Spec_Win_Dim%2==1)
    fr_max = floor ((Spec_Win_Dim-1)/sqrt(2));
  
  plhs[0] = mxCreateDoubleMatrix (1,fr_max,mxREAL);
  Annulus_Power = mxGetPr(plhs[0]);
  
  plhs[1] = mxCreateDoubleMatrix (1,fr_max,mxREAL);
  fr = mxGetPr(plhs[1]);
  
  if (nlhs==2)
    count = mxCalloc (1,sizeof (double));
  else
    {
      plhs[2] = mxCreateDoubleMatrix (1,1,mxREAL);
      count = mxGetPr(plhs[2]);
    }

  /* C Computational Routine Call */
  
  Rad_Avg_PS (Image_Data, Annulus_Power, fr, count);
  
  return;
}



