#include <math.h>
#include "mex.h"
#define pi 3.14159265358979

typedef struct {  
        double row;
        double col;
        double cost;
        int num_points;
        int num_min_pixels;
        int radial_bin_index;
} cost_array_node;

int image_row;
int image_col;
int image_chn;
int chn_I;
int chn_J;
int num_radial_bins;
double max_distance;
double delta_r;

/********************************************************************************/
/*                                                                              */
/*                                                                              */
/********************************************************************************/
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
                        rra=ra[l-1];
                        }
                else {
                        rra=ra[ir-1];

                        ra[ir-1]=ra[0];
                        if (--ir==1){
                                ra[0]=rra;
                                break;
                                }
                        }
                i=l;
                j=l+l;
                while(j<=ir){
                        if (j < ir && ra[j-1].cost < ra[j].cost)j++;
                        if (rra.cost < ra[j-1].cost){
                                ra[i-1]=ra[j-1];
                                i=j;
                                j<<=1;
                                }
                        else j=ir+1;
                        }
                ra[i-1]=rra;
                }
        return;
}

void create_cost_array(cost_array_node **cost_array)
{
        int i,j,m,n,q;

        m=-image_row/2+1;
        n=-image_col/2+1;
        for (j=0; j<image_col; j++){
                for (i=0; i<image_row; i++){
                        (*cost_array)[i+j*image_row].row=m+i;
                        (*cost_array)[i+j*image_row].col=n+j;
                        (*cost_array)[i+j*image_row].cost=sqrt((m+i)*(m+i)+(n+j)*(n+j));
			(*cost_array)[i+j*image_row].num_points=0;
			(*cost_array)[i+j*image_row].num_min_pixels=0;

			q=(int)ceil((*cost_array)[i+j*image_row].cost/delta_r)-1; if (q==-1) q=0;
			(*cost_array)[i+j*image_row].radial_bin_index=q;
                        }
                }
        hpsort_cost_array(image_row*image_col, *cost_array);
        return;
}

void vpc(double pair_correlation[], unsigned char input_image[])
{
        int m, n, r, x_coor, y_coor;
	double gray_level=0.0;
	double max_cost;
        cost_array_node *cost_array;
	int *num_pixels, *num_bins;

        cost_array=(cost_array_node*)mxCalloc(image_row*image_col, sizeof(cost_array_node));
        create_cost_array(&cost_array);
        for (n=0; n<image_col; n++){
                for (m=0; m<image_row; m++){
                        if (input_image[m+(n+chn_J*image_col)*image_row]==(unsigned char)1){
 				for (r=0; r<image_row*image_col; r++){
				         if (cost_array[r].cost > max_distance) break;
				         if (cost_array[r].row >= m ||
                                             cost_array[r].row >= (image_row-m) ||
					     cost_array[r].col >= n ||
                                             cost_array[r].col >= (image_col-n)){
						 break;
					         }
					 if (cost_array[r].cost >= max_distance){
						 break;
					         }
				         x_coor=m+cost_array[r].row;
					 y_coor=n+cost_array[r].col;
				         cost_array[r].num_points++;
					 if (input_image[x_coor+(y_coor+chn_I*image_col)*image_row]==(unsigned char)1)
					        cost_array[r].num_min_pixels++;
				         }
			        }
		         }
                }
	num_pixels=(int*)mxCalloc(num_radial_bins, sizeof(int));
	num_bins=(int*)mxCalloc(num_radial_bins, sizeof(int));
	for (n=0; n<image_row*image_col; n++){
	        if (input_image[n+chn_I*image_col*image_row]==(unsigned char)1){
		         gray_level++;
	                 }
                }
	gray_level=gray_level/((double)image_row*(double)image_col);
	for (m=0; m<image_row*image_col; m++){
	        if (cost_array[m].radial_bin_index>=num_radial_bins) break;
	        pair_correlation[cost_array[m].radial_bin_index]+=
		         (double)cost_array[m].num_min_pixels/(gray_level*cost_array[m].num_points);
	        num_pixels[cost_array[m].radial_bin_index]++;
		}
	for (m=0; m<num_radial_bins; m++){
	        pair_correlation[m]=pair_correlation[m]/(double)num_pixels[m];
	        }
	mxFree(num_pixels);
	mxFree(num_bins);
	mxFree(cost_array);
	return;
}

/*******************************************************************************/ 
/* mexFUNCTION                                                                 */
/* Gateway routine for use with MATLAB.                                        */
/*******************************************************************************/ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        int m, n;
        unsigned char *input_image;
        double *pair_correlation, *radial_data, *channels;
	int number_of_dims;
        const int *dim_array;

        if (nrhs<2 | nrhs>4)
                mexErrMsgTxt("PC accepts one or two!");
        else if (nlhs>2)
                mexErrMsgTxt("PC outputs one or two arguments!");
        else if (!mxIsNumeric(prhs[0]) ||
                  mxIsComplex(prhs[0]) ||
                  mxIsSparse(prhs[0])  ||
                 !mxIsUint8(prhs[0])   ||
		 !mxIsLogical(prhs[0]))
                mexErrMsgTxt("Input X must be a real matrix of type uint8 (logical)!");
      	else if (!mxIsNumeric(prhs[1]) ||
                 mxIsComplex(prhs[1]) ||
                 mxIsSparse(prhs[1])  ||
		!mxIsDouble(prhs[1])  ||
		(mxGetM(prhs[1])*mxGetN(prhs[1])!=2))
	        mexErrMsgTxt("Input C must be a real vector of length 2!");
	else if (nrhs>2 && (!mxIsNumeric(prhs[2]) ||
		  	     mxIsComplex(prhs[2]) ||
                  	     mxIsSparse(prhs[2])  ||
		 	   !(mxGetM(prhs[2])==1 && mxGetN(prhs[2])==1)))
		mexErrMsgTxt("Input dr must be a real scalar!");
	else if (nrhs>3 && (!mxIsNumeric(prhs[3]) ||
		  	     mxIsComplex(prhs[3]) ||
                  	     mxIsSparse(prhs[3])  ||
		 	   !(mxGetM(prhs[3])==1 && mxGetN(prhs[3])==1)))
		mexErrMsgTxt("Input Rmx must be a real scalar!");

    /****** Get row and column sizes for input image **************************************/
	number_of_dims=mxGetNumberOfDimensions(prhs[0]);
	if (number_of_dims<2 || number_of_dims>3)
	        mexErrMsgTxt("Input image X must be 2 or 3 dimensional!");
	else if (number_of_dims==2){
                image_row=mxGetM(prhs[0]);
                image_col=mxGetN(prhs[0]);
                image_chn=1;
	        }
	else{
	        dim_array=mxGetDimensions(prhs[0]);
                image_row=dim_array[0];
                image_col=dim_array[1];
                image_chn=dim_array[2];
	        }
	input_image=(unsigned char*)mxGetPr(prhs[0]);

	channels=mxGetPr(prhs[1]);
	chn_I=(int)channels[0]-1;
	chn_J=(int)channels[1]-1;

	if (chn_I<0 || chn_I>=image_chn || chn_J<0 || chn_J>=image_chn)
	        mexErrMsgTxt("Non-existing channels specified by C!");

	if (nrhs>=3)
	        delta_r=mxGetScalar(prhs[2]);
	else
	        delta_r=0.5;

	if (nrhs==4)
	        max_distance=mxGetScalar(prhs[3]);
	else
	        max_distance=image_row/2;
	if (max_distance > image_col) max_distance=image_col/2;
	if (max_distance > image_row) max_distance=image_row/2;
	num_radial_bins=(int)ceil(max_distance/delta_r)-1;

        plhs[0]=mxCreateDoubleMatrix(1, num_radial_bins, mxREAL);
	pair_correlation=mxGetPr(plhs[0]);
	if (nlhs==2){
        	plhs[1]=mxCreateDoubleMatrix(1, num_radial_bins, mxREAL);
		radial_data=mxGetPr(plhs[1]);
		for (m=0; m<num_radial_bins; m++){
		        radial_data[m]=((double)m+0.5)*delta_r;
		        }
		}
	vpc(pair_correlation, input_image);
        return;
}






