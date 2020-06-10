/****************************************************************************************/
/*                                                                                      */
/* Directional Distribution Function							                        */
/*      This function calculates the directional distribution of points of a		    */
/*	point process returning the distributions in matrix T.				                */
/*                                                                                      */
/* Synopsis:                                                                            */
/*      T=ddf(H, r1, r2)								                                */
/*              T=directional distribution of points such that r1<=|x-y|<r2		        */
/*              H=Input image (~=0 a point)                                             */
/*      T=ddf(H, r1, r2, N)								                                */
/*              N=Number of angular bins (default=72)					                */
/*      [T, a]=ddf(H, r1, r2, N)							                            */
/*              a=vector returning center angle of each bin.				            */
/*      [T, a]=ddf(H, r1, r2, N, str)							                        */
/*              str=String specifying clustered point process when str='cluster'.	    */
/*                                                                                      */
/* Daniel Leo Lau                                                                       */
/* Copyright August 17, 1997                                                            */
/*                                                                                      */
/****************************************************************************************/ 
#include <math.h>
#include "mex.h"
#define pi 3.14159265358979

int image_row;
int image_col;
int *pixel_list_M;
int *pixel_list_N;
int maximum_rank;
int num_radial_bins;
int num_angular_bins;
double dR;
double dA;

typedef struct {  
        double row;
        double col;
        double cost;
} cost_array_node;

void hpsort_cost_array(int n, cost_array_node *ra)
{
        int i, ir, j, l;
        cost_array_node rra;

        if (n<2) return;
        l=(n >> 1)+1;
        ir=n;

        for( ; ; ){
                if (l>1){
                        --l;
                        rra.row=ra[l-1].row;
                        rra.col=ra[l-1].col;
                        rra.cost=ra[l-1].cost;
                        }
                else {
                        rra.row=ra[ir-1].row;
                        rra.col=ra[ir-1].col;
                        rra.cost=ra[ir-1].cost;

                        ra[ir-1].row=ra[0].row;
                        ra[ir-1].col=ra[0].col;
                        ra[ir-1].cost=ra[0].cost;
                        if (--ir==1){
                                ra[0].row=rra.row;
                                ra[0].col=rra.col;
                                ra[0].cost=rra.cost;
                                break;
                                }
                        }
                i=l;
                j=l+l;
                while(j<=ir){
                        if (j < ir && ra[j-1].cost < ra[j].cost)j++;
                        if (rra.cost < ra[j-1].cost){
                                ra[i-1].row=ra[j-1].row;
                                ra[i-1].col=ra[j-1].col;
                                ra[i-1].cost=ra[j-1].cost;
                                i=j;
                                j<<=1;
                                }
                        else j=ir+1;
                        }
                ra[i-1].row=rra.row;
                ra[i-1].col=rra.col;
                ra[i-1].cost=rra.cost;
                }
        return;
}

void create_cost_array(cost_array_node **cost_array)
{
        int i,j,m,n;

        m=-image_row/2+1;
        n=-image_col/2+1;
        for (j=0; j<image_col; j++){
                for (i=0; i<image_row; i++){
                        (*cost_array)[i+j*image_row].row=m+i;
                        (*cost_array)[i+j*image_row].col=n+j;
                        (*cost_array)[i+j*image_row].cost=sqrt((m+i)*(m+i)+(n+j)*(n+j));
                        }
                }
        hpsort_cost_array(image_row*image_col, *cost_array);
        return;
}

void directional_distribution_function(double output_angle[], double angular_data[],
                                       unsigned char input_image[], double inner_radii[],
                                       double outer_radii[])
{
        int r, m, n, j, md, nd, q;
        cost_array_node *cost_array;
        double angle, cost;
	double *area_bins, *area_rings, *point_rings;

        cost_array=(cost_array_node*)mxCalloc(image_row*image_col, sizeof(cost_array_node));
        create_cost_array(&cost_array);
	area_bins=(double*)mxCalloc(num_radial_bins*num_angular_bins, sizeof(double));
	area_rings=(double*)mxCalloc(num_radial_bins, sizeof(double));
	point_rings=(double*)mxCalloc(num_radial_bins, sizeof(double));

        for (m=0; m<num_angular_bins; m++) angular_data[m]=(double)m*dA;
	for (m=1; m<image_row*image_col; m++){
	        cost=cost_array[m].cost;
		angle=atan2(-1*(double)cost_array[m].row, (double)cost_array[m].col);
                if (angle<-dA/2.0) angle=angle+2*pi;
                angle+=dA/2.0;
		for (n=0; n<num_radial_bins; n++){
		        if (inner_radii[n] <= cost && cost < outer_radii[n]){
			        area_rings[n]++;
                                q=(int)floor(angle/dA);
			        area_bins[n+q*num_radial_bins]++;
			        }
		        }
	        }
        for (r=0; r<num_radial_bins; r++){
                for (n=0; n<image_col; n++){
                        for (m=0; m<image_row; m++){
                                if (input_image[m+n*image_row]!=0){
                                        if (m <= outer_radii[r] || 
                                            n <= outer_radii[r] || 
                                           (image_row-m) <= outer_radii[r] ||
                                           (image_col-n) <= outer_radii[r]) continue;
                                        j=1;
                                        while (cost_array[j].cost < outer_radii[r] && j < image_row*image_col){
                                                if (cost_array[j].cost >= inner_radii[r]){
                                                        md=m+cost_array[j].row;
                                                        nd=n+cost_array[j].col;
                                                        if (input_image[md+nd*image_row]!=0){
                                                                angle=atan2(-1*(double)cost_array[j].row, 
                                                                               (double)cost_array[j].col);
                                                                if (angle<-dA/2.0) angle=angle+2*pi;
                                                                angle+=dA/2.0;
                                                                q=(int)floor(angle/dA);
                                                                if (q<0 || q>=num_angular_bins){
                                                                        mexErrMsgTxt("Error in angular bin");
								        }
                                                                output_angle[r+q*num_radial_bins]++;
                                                                }
                                                        }
                                                j++;
                                                }
                                        }
                                }
                        }
                point_rings[r]=0.0;
                for (n=0; n<num_angular_bins; n++){
                        point_rings[r]+=output_angle[r+n*num_radial_bins];
                        }
                if (point_rings[r]>0.0){
                        for (n=0; n<num_angular_bins; n++){
			        if (area_bins[r+n*num_radial_bins]>0){
                                       output_angle[r+n*num_radial_bins]=(output_angle[r+n*num_radial_bins]/
                                                                  area_bins[r+n*num_radial_bins])/
				                                  (point_rings[r]/
				                                  area_rings[r]);				
                                       }
    			        else {
				       output_angle[r+n*num_radial_bins]=0.0;
				       }
                                }
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
        int m, n;
        unsigned char *input_data;
        double *output_cost, *output_angle, *angular_data;
	double *inner_radii, *outer_radii;
	char cluster_str[5];

        if (nrhs<3 || nrhs>4)
                mexErrMsgTxt("DDF accepts three or four input arguments!");
        else if (nlhs>2)
                mexErrMsgTxt("DDF outputs up to two output arguments!");
	else if (!mxIsLogical(prhs[0])) 
        mexErrMsgTxt("Input X must be a real matrix of type LOGICAL!");
	else if (!mxIsNumeric(prhs[1]) ||
		  mxIsComplex(prhs[1]) ||
                  mxIsSparse(prhs[1])  ||
		 (mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1))
		mexErrMsgTxt("Input Ri must be a real vector!");
	else if (!mxIsNumeric(prhs[2]) ||
		  mxIsComplex(prhs[2]) ||
                  mxIsSparse(prhs[2])  ||
		 (mxGetM(prhs[2])!=1 && mxGetN(prhs[2])!=1))
		mexErrMsgTxt("Input Ro must be a real vector!");
	else if (nrhs>3 && (!mxIsNumeric(prhs[3]) ||
			     mxIsComplex(prhs[3]) ||
                             mxIsSparse(prhs[3])  ||
			    (mxGetM(prhs[3])*mxGetN(prhs[3])!=1)))
		mexErrMsgTxt("Input N must be a real scalar!");

	if (nrhs>3) num_angular_bins=(int)mxGetScalar(prhs[3]); else num_angular_bins=72;
	dA=2.0*pi/(double)num_angular_bins;

        image_row=mxGetM(prhs[0]);
        image_col=mxGetN(prhs[0]);

	num_radial_bins=mxGetM(prhs[1])*mxGetN(prhs[1]);
	if (num_radial_bins!=mxGetM(prhs[2])*mxGetN(prhs[2]))
		mexErrMsgTxt("Vectors Ri and Ro must be same length!");
	inner_radii=mxGetPr(prhs[1]);
	outer_radii=mxGetPr(prhs[2]);
	
	input_data=(unsigned char*)mxGetLogicals(prhs[0]);

        plhs[0]=mxCreateDoubleMatrix(num_radial_bins, num_angular_bins, mxREAL);
        output_angle=mxGetPr(plhs[0]);
	if (nlhs>1){
        	plhs[1]=mxCreateDoubleMatrix(1, num_angular_bins, mxREAL);
        	angular_data=mxGetPr(plhs[1]);
		}
	else    angular_data=(double*)mxCalloc(num_angular_bins, sizeof(double));
        directional_distribution_function(output_angle, angular_data, input_data, inner_radii, outer_radii);
        return;
}
