/*******************************************************************************/ 
/*
/*                            MATLAB FUNCTION: ANIS
/*
/* Title:       ANISotropy    
/* Description: From the Power Spectrum, it estimates the anistropy metric
/*              by taking measuring the variance of power samples over concentric
/*              annular rings of radius fr where for simplicity increments in
/*              fr are taken to be delta=1  
/*              Reference: Robert Ulichney, "Digital Halftoning"
/*
/* MATLAB Prompt Call: [P,fr,count] = ANIS(X,[M],[FM],[DM DN])
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
void   Anisotropy_Metric (double*, double*, double*, double*);

/*******************************************************************************/ 
/*******************************************************************************/ 
/*******************************************************************************/ 
void Anisotropy_Metric (double *Image_Data, double *Anisotropy, double *fr, double *count)
{
  int i,j;
  double  f1,f2,f3,f4,f_min;
  int     radial_index;
  double  *Power_Estimate;
  double  *Radial_Samples, *Annulus_Power;
  int     delta=1;

  Power_Estimate = mxCalloc(Spec_Win_Dim*Spec_Win_Dim,sizeof(double)); 
  Radial_Samples = mxCalloc(fr_max,sizeof(double));
  Annulus_Power = mxCalloc(fr_max,sizeof(double));

  /* Computes the Power Estimate over the given Spectral Window */

  Compute_Power_Estimate(Image_Data,Power_Estimate,count);

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
	  
	  if (radial_index>=0)
	    {
	      Annulus_Power[radial_index] += Power_Estimate[i+j*Spec_Win_Dim];
	      Radial_Samples[radial_index] ++;
	    }
	}
    }

  for (i=0;i<fr_max;i++)
    Annulus_Power[i]=Annulus_Power[i]/(Radial_Samples[i]);

  for (i=1;i<Spec_Win_Dim;i++)
    {
      for (j=1;j<Spec_Win_Dim;j++)
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
	 
	  if (radial_index>=0)
	    {
	      Anisotropy[radial_index]+=(Power_Estimate[i+j*Spec_Win_Dim]-Annulus_Power[radial_index])
                                   *(Power_Estimate[i+j*Spec_Win_Dim]-Annulus_Power[radial_index]);
	
	    }
	}
    }

  for (i=0;i<fr_max;i++)
    {
    Anisotropy[i] = Anisotropy[i]/(Annulus_Power[i]*Annulus_Power[i]*(Radial_Samples[i]-1));
    Anisotropy[i] = 10*log10(Anisotropy[i]);
    }

  mxFree(Annulus_Power);  
  mxFree(Radial_Samples);
  mxFree(Power_Estimate);
  return;
}

/*******************************************************************************/ 
/*******************************************************************************/ 
/*******************************************************************************/ 
void Compute_Power_Estimate (double *Image_Data, double *Power_Estimate, double *Num_Samples)
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

  *(Num_Samples)= (double)count;

  mxFree (Frequency_Data);
  mxFree (Sample_Window);
  return;
}

/*******************************************************************************/ 
/*******************************************************************************/ 
/*******************************************************************************/ 
void absfft2 (double *Image_PS, double *Input_Image, int Num_Row, int Num_Col)
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
      Image_PS[i] = (Real_Data[i]*Real_Data[i]+Imag_Data[i]*Imag_Data[i])/(Num_Row*Num_Col);
      }
  else
      {
      for (i=0;i<Num_Row*Num_Col;i++)
      Image_PS[i] = (Real_Data[i]*Real_Data[i])/(Num_Row*Num_Col);
      }
      
  mxDestroyArray(Input_Array[0]);
  mxDestroyArray(Output_Array[0]);
  return;
}

/*******************************************************************************/ 
/*******************************************************************************/ 
/*******************************************************************************/ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int    Spat_Win_M, Spat_Win_N;
  int    Spec_Win_M, Spec_Win_N;
  int    Shift_M, Shift_N, i;
  
  double *Spat_Win_Data, *Spec_Win_Data, *Shift_Data;
  double *Image_Data, *Anisotropy, *fr, *count;

  /* Check for errors in User Calls from Matlab pronpt. */
  
  if ((nrhs<2) || (nrhs>4))
    mexErrMsgTxt("This Function requires two to four Input Arguments!");
  
  else if (nlhs<2 || nlhs>3)
    mexErrMsgTxt("This Function returns 2 or 3 Output Argument!");
  
  else if (!mxIsNumeric(prhs[0]) ||
	   mxIsComplex(prhs[0]) ||
	   mxIsSparse(prhs[0])  ||
	   !mxIsDouble(prhs[0])   )
    mexErrMsgTxt("Input Argument X must be a REAL matrix of type UINT8(LOGICAL)!");
  
  else if (!mxIsNumeric(prhs[1]) ||
	   mxIsComplex(prhs[1]) ||
	   mxIsSparse(prhs[1])  ||
	   !mxIsDouble(prhs[1])   )
    mexErrMsgTxt("Input Argumnent [M] must be a REAL scalar of typeDOUBLE!");
  
  else if ((nrhs>2) && (!mxIsNumeric(prhs[2]) ||
			mxIsComplex(prhs[2]) ||
			mxIsSparse(prhs[2])  ||
			!mxIsDouble(prhs[2])  ))
    mexErrMsgTxt("Input Argumnent [FM] must be a REAL scalar of typeDOUBLE!");
  
  else if ((nrhs>3) && (!mxIsNumeric(prhs[3]) ||
			mxIsComplex(prhs[3]) ||
			mxIsSparse(prhs[3])  ||
			!mxIsDouble(prhs[3])  ))
    mexErrMsgTxt("Input Argumnent [DM DN] must be a REAL vector of typeDOUBLE!");

  /* Assigns the pointer of type double to the Image_Data Matrix & */ 	  
  /* Get row & column sizes of the Image Matrix, Vectors [M],[FM] & [DM DN] */

  Image_Data = (double*)mxGetPr (prhs[0]);
  Image_M    = mxGetM  (prhs[0]);
  Image_N    = mxGetN  (prhs[0]);

  /* Checks on Argument 1 */
  
  Spat_Win_M = mxGetM( prhs[1] );
  Spat_Win_N = mxGetN( prhs[1] );
  
  if ((Spat_Win_M*Spat_Win_N)>1)
    mexErrMsgTxt ("Input Argument N must be scalar!!!");
  
  else if ((Spat_Win_M*Spat_Win_N)==1)
    Spat_Win_Dim = (int) mxGetScalar (prhs[1]);
  
  if (Spat_Win_Dim == 0)
    mexErrMsgTxt ("Spatial Window Dimensions cannot be 0");
  
  /* Set Default values for Arguments 3 & 4 */
  /* Assuming that the other arguments will be zero */          
  
  Spec_Win_Dim = Spat_Win_Dim;
  Shift_Row    = Spat_Win_Dim;
  Shift_Col    = Spat_Win_Dim;

  /* Checks on Argument 2 */
  
  if (nrhs>2)
    {
      Spec_Win_M = mxGetM  (prhs[2]);
      Spec_Win_N = mxGetN  (prhs[2]); 
      
      if ((Spec_Win_M*Spec_Win_N)>1)
	mexErrMsgTxt ("Input Argument 2 can have dimensions 1X1");         
      
      else if ((Spec_Win_M*Spec_Win_N)==1)
	Spec_Win_Dim = (int) mxGetScalar (prhs[2]); 
      
      if (Spec_Win_Dim == 0)
	Spec_Win_Dim = Spat_Win_Dim;
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
  Anisotropy = mxGetPr(plhs[0]);
  
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
  
  Anisotropy_Metric (Image_Data, Anisotropy, fr, count);
  
  return;
}




