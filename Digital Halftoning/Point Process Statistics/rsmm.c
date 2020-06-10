#include <math.h>
#include "mex.h"
#define pi 3.14159265358979

typedef struct {  
    double row;
    double col;
    double cost;
    int num_points;
    int num_min_pixels;
} cost_array_node;

int image_row;
int image_col;
double max_distance;

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

/********************************************************************************/
/*                                                                              */
/*                                                                              */
/********************************************************************************/
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
			(*cost_array)[i+j*image_row].num_points=0;
			(*cost_array)[i+j*image_row].num_min_pixels=0;
		}
	}
    hpsort_cost_array(image_row*image_col, *cost_array);
    return;
}

/********************************************************************************/
/*                                                                              */
/*                                                                              */
/********************************************************************************/
void srmm_function(double srmm[], unsigned char input_image[])
{
	int m, n, r, x_coor, y_coor;
	double gray_level=0.0;
	cost_array_node *cost_array;

    cost_array=(cost_array_node*)mxCalloc(image_row*image_col, sizeof(cost_array_node));
    create_cost_array(&cost_array);
    for (n=0; n<image_col; n++){
		for (m=0; m<image_row; m++){
			if (input_image[m+n*image_row]==(unsigned char)1){
				gray_level=gray_level+1.0;
				for (r=1; r<image_row*image_col; r++){
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
					if (input_image[x_coor+y_coor*image_row]==(unsigned char)1) cost_array[r].num_min_pixels++;
				}
			}
		}
    }
    
	gray_level=gray_level/((double)image_row*(double)image_col);
	for (m=0; m<image_row*image_col; m++){
	    x_coor=cost_array[m].row+image_row/2-1;
		y_coor=cost_array[m].col+image_col/2-1;
		if (x_coor<image_row & y_coor<image_col){
		    srmm[x_coor+y_coor*image_row]=(double)cost_array[m].num_min_pixels/((double)cost_array[m].num_points*gray_level);
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
	unsigned char *input_image;
	double *srmm;
	
	if (nrhs<1 | nrhs>2)
        mexErrMsgTxt("SRMM accepts one or two!");
	else if (nlhs>1)
        mexErrMsgTxt("SRMM outputs only one argument!");
	else if (!mxIsLogical(prhs[0])) 
        mexErrMsgTxt("Input X must be a real matrix of type LOGICAL!");
	else if (nrhs>1 && (!mxIsNumeric(prhs[1]) ||
      	     mxIsComplex(prhs[1]) ||
             mxIsSparse(prhs[1])  ||
     	   !(mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1)))
	    mexErrMsgTxt("Input Rmx must be a real scalar!");
	    
	image_row=mxGetM(prhs[0]);
	image_col=mxGetN(prhs[0]);
	input_image=(unsigned char*)mxGetLogicals(prhs[0]);

	if (nrhs==2)
	    max_distance=mxGetScalar(prhs[1]);
	else
	    max_distance=image_row;

    plhs[0]=mxCreateDoubleMatrix(image_row, image_col, mxREAL);
	srmm=mxGetPr(plhs[0]);
    srmm_function(srmm, input_image);
    return;
}
