#include <math.h>
#include "mex.h"
#define pi 3.14159265358979

typedef struct {  
        double row;
        double col;
        double cost;
} cost_array_node;

struct cluster{
        int cluster_number;
        int *M_coor; double M_center;
        int *N_coor; double N_center;
        int size;
        struct cluster *next;
};

struct cluster_image{
        unsigned char *image;
        int frequency;
        int image_number;
        struct cluster_image *next;
};

int num_rings;
int image_row;
int image_col;
int *pixel_list_M;
int *pixel_list_N;
int maximum_rank;
int num_radial_bins;
int num_clusters;
double dR;
struct cluster *cluster_list;

/********************************************************************************/
/*                                                                              */
/* MAGICWAND_MN: performs the search of all neighboring pixels to the top, left,*/
/*            bottom, and right of pixel(m,n) and returns the cooridinates.     */
/*                                                                              */
/********************************************************************************/
int magic_wand_mn(int *pixel_list_M, int *pixel_list_N, unsigned char *Y, unsigned char *X, int m, int n)
{
        int r, s, t, length_pixel_list;
        int first_previous_iteration, last_previous_iteration, next_available_slot;
        unsigned char fixed_level;

        length_pixel_list=image_row*image_col;
        fixed_level=X[m+n*image_row];
        Y[m+n*image_row]=1;
        
        pixel_list_M[0]=m;
        pixel_list_N[0]=n;
        first_previous_iteration=0;
        last_previous_iteration=0;
        next_available_slot=1;
        while(1){
                for (r=first_previous_iteration; r<=last_previous_iteration; r++){
                        s=pixel_list_M[r]-1; t=pixel_list_N[r];
                        if (s>=0 && Y[s+t*image_row]!=1 && (fixed_level==X[s+t*image_row])){
                                pixel_list_M[next_available_slot]=s;
                                pixel_list_N[next_available_slot]=t;
                                Y[s+t*image_row]=1;
                                next_available_slot++;
                                if (next_available_slot==length_pixel_list) break;
                                }
                        s=pixel_list_M[r]; t=pixel_list_N[r]-1;
                        if (t>=0 && Y[s+t*image_row]!=1 && (fixed_level==X[s+t*image_row])){
                                pixel_list_M[next_available_slot]=s;
                                pixel_list_N[next_available_slot]=t;
                                Y[s+t*image_row]=1;
                                next_available_slot++;
                                if (next_available_slot==length_pixel_list) break;
                                }
                        s=pixel_list_M[r]+1; t=pixel_list_N[r];
                        if (s<image_row && Y[s+t*image_row]!=1 && (fixed_level==X[s+t*image_row])){
                                pixel_list_M[next_available_slot]=s;
                                pixel_list_N[next_available_slot]=t;
                                Y[s+t*image_row]=1;
                                next_available_slot++;
                                if (next_available_slot==length_pixel_list) break;
                                }
                        s=pixel_list_M[r]; t=pixel_list_N[r]+1;
                        if (t<image_col && Y[s+t*image_row]!=1 && (fixed_level==X[s+t*image_row])){
                                pixel_list_M[next_available_slot]=s;
                                pixel_list_N[next_available_slot]=t;
                                Y[s+t*image_row]=1;
                                next_available_slot++;
                                if (next_available_slot==length_pixel_list) break;
                                }
                        }
                if (last_previous_iteration==next_available_slot-1) break;
                first_previous_iteration=last_previous_iteration+1;
                last_previous_iteration=next_available_slot-1;
                }
        return(next_available_slot);
}

/********************************************************************************/
/*                                                                              */
/*                                                                              */
/********************************************************************************/
void append_to_list(struct cluster *new_cluster)
{
        struct cluster *current_cluster;

        if (cluster_list==NULL){
                cluster_list=new_cluster;
                return;
                }
        current_cluster=cluster_list;
        while(current_cluster->next!=NULL){
                current_cluster=current_cluster->next;
                }
        current_cluster->next=new_cluster;
        new_cluster->cluster_number=current_cluster->cluster_number+1;
        return;
}

/********************************************************************************/
/*                                                                              */
/*                                                                              */
/********************************************************************************/
int get_all_clusters(unsigned char *input_image, unsigned char pixel_color)
{
        int *coor_m, *coor_n, number_of_clusters=0;
        int m, n, r, length_coor;
        unsigned char *Y;
        struct cluster *new_cluster;

        Y=(unsigned char*)mxCalloc(image_row*image_col, sizeof(unsigned char));
        coor_m=(int*)mxCalloc(image_row*image_col, sizeof(int));
        coor_n=(int*)mxCalloc(image_row*image_col, sizeof(int));
        
        for (n=0; n<image_col; n++){
                for (m=0; m<image_row; m++){
                        if (Y[m+n*image_row]!=1 && input_image[m+n*image_row]==pixel_color){
                                length_coor=magic_wand_mn(coor_m, coor_n, Y, input_image, m, n);
                                new_cluster=(struct cluster*)mxCalloc(1, sizeof(struct cluster));
                                new_cluster->M_coor=(int*)mxCalloc(length_coor, sizeof(int));
                                new_cluster->N_coor=(int*)mxCalloc(length_coor, sizeof(int));
                                new_cluster->M_center=0.0;
                                new_cluster->N_center=0.0;
                                new_cluster->size=length_coor;
                                for (r=0; r<length_coor; r++){
                                        new_cluster->M_coor[r]=coor_m[r];
                                        new_cluster->N_coor[r]=coor_n[r];
                                        new_cluster->M_center+=(double)coor_m[r];
                                        new_cluster->N_center+=(double)coor_n[r];
                                        }
                                new_cluster->M_center=new_cluster->M_center/length_coor;
                                new_cluster->N_center=new_cluster->N_center/length_coor;
                                new_cluster->next=NULL;
                                append_to_list(new_cluster);
                                number_of_clusters++;
                                }
                        }
                }
        mxFree(coor_m);
        mxFree(coor_n);
        return(number_of_clusters);
}

/********************************************************************************/
/*                                                                              */
/*                                                                              */
/********************************************************************************/
int get_all_pixels(unsigned char *input_image, unsigned char pixel_color)
{
	int m,n, number_of_clusters=0;
	struct cluster *new_cluster;

	for (n=0; n<image_col; n++){
		for (m=0; m<image_row; m++){
			if (input_image[m+n*image_row]==pixel_color){
				cluster_list=(struct cluster*)mxCalloc(1, sizeof(struct cluster));
                                cluster_list->M_center=(double)m;
                                cluster_list->N_center=(double)n;
                                cluster_list->size=0;
                                cluster_list->next=NULL;
				new_cluster=cluster_list;
				number_of_clusters=1;
				break;
				}
			}
		if (cluster_list!=NULL) break;
		}
	for (n=0; n<image_col; n++){
		for (m=0; m<image_row; m++){
			if (input_image[m+n*image_row]==pixel_color){
				new_cluster->next=(struct cluster*)mxCalloc(1, sizeof(struct cluster));
                                new_cluster->next->M_center=(double)m;
                                new_cluster->next->N_center=(double)n;
                                new_cluster->next->size=0;
                                new_cluster->next->next=NULL;
				new_cluster=new_cluster->next;
				number_of_clusters++;
				}
			}
		}
	return(number_of_clusters);
}

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

void pair_correlation_function(double pair_correlation[], double radial_data[], unsigned char input_image[])
{
    int r, m, n, q, pixel_list_length=0;
	double *num_points_data, *volume_data, max_distance, intensity;
    cost_array_node *cost_array, *pixel_list, *sort_list;
    struct cluster *current_cluster;

	volume_data=(double*)mxCalloc(num_radial_bins, sizeof(double));
	num_points_data=(double*)mxCalloc(num_radial_bins, sizeof(double));

    cost_array=(cost_array_node*)mxCalloc(image_row*image_col, sizeof(cost_array_node));
    create_cost_array(&cost_array);
	for (m=0; m<image_row*image_col; m++){
		q=(int)ceil(cost_array[m].cost/dR)-1; if (q==-1) q=0;
		if (q>=num_radial_bins) break;
		volume_data[q]++;
		}
    for (m=0; m<num_radial_bins; m++) radial_data[m]=((double)m+0.5)*dR;

    for (m=0; m<image_row*image_col; m++){
            pixel_list_length+=(input_image[m]==(unsigned char)1);
            }
    pixel_list=(cost_array_node*)mxCalloc(pixel_list_length, sizeof(cost_array_node));
    sort_list=(cost_array_node*)mxCalloc(pixel_list_length, sizeof(cost_array_node));
	intensity=(double)pixel_list_length/(double)image_row/(double)image_col;
    r=0;
    for (n=0; n<image_col; n++){
            for (m=0; m<image_row; m++){
                    if (input_image[m+n*image_row]!=(unsigned char)0){
                            pixel_list[r].row=(double)m;
                            pixel_list[r].col=(double)n;
                            r++;
                            }
                    }
            }

    current_cluster=cluster_list;
    while(current_cluster!=NULL){
		max_distance=(double)(num_radial_bins)*dR;
		if (max_distance>current_cluster->M_center)
        		max_distance=current_cluster->M_center;
        if (max_distance>current_cluster->N_center) 
                max_distance=current_cluster->N_center;
        if (max_distance>image_col-current_cluster->N_center)
                max_distance=image_col-current_cluster->N_center;
        if (max_distance>image_row-current_cluster->M_center)
                max_distance=image_row-current_cluster->M_center;
		q=(int)ceil(max_distance/dR)-1; if (q==-1) q=0;
		for (n=0; n<q; n++) num_points_data[n]++;
		max_distance=dR*q;
        r=0;
        for (m=0; m<pixel_list_length; m++){
            sort_list[r].row=pixel_list[m].row-current_cluster->M_center;
            sort_list[r].col=pixel_list[m].col-current_cluster->N_center;
            sort_list[r].cost=sqrt(sort_list[r].row*sort_list[r].row + sort_list[r].col*sort_list[r].col);
            if (sort_list[r].cost <= max_distance){
                    r++;
                    }
            }
        for (m=0; m<r; m++){
            q=(int)ceil(sort_list[m].cost/dR)-1; if (q==-1) q=0;
			pair_correlation[q]++;
            }
        current_cluster=current_cluster->next;
    }
	for (m=0; m<num_radial_bins; m++){
		/*pair_correlation[m]=num_points_data[m]; continue;*/
		/*pair_correlation[m]=volume_data[m]; continue;*/
		if (volume_data[m]==0.0 || num_points_data[m]==0.0){
			pair_correlation[m]=0.0;
			continue;
			}
		pair_correlation[m]=((pair_correlation[m]/volume_data[m])/num_points_data[m])/intensity;
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
    double *radial_data, *pair_correlation, max_distance;
	char cluster_str[5];

    if (nrhs<2 || nrhs>4)
        mexErrMsgTxt("PC accepts two to four input arguments!");
    else if (nlhs>2)
        mexErrMsgTxt("PC outputs up to two output arguments!");
	else if (!mxIsLogical(prhs[0])) 
        mexErrMsgTxt("Input X must be a real matrix of type LOGICAL!");
	else if (!mxIsNumeric(prhs[1]) ||
		      mxIsComplex(prhs[1]) ||
              mxIsSparse(prhs[1])  ||
		     (mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1))
		mexErrMsgTxt("Input dr must be a real scalar!");
	else if (nrhs>2 && (!mxIsNumeric(prhs[2]) ||
		  	             mxIsComplex(prhs[2]) ||
                  	     mxIsSparse(prhs[2])  ||
		 	            (mxGetM(prhs[2])!=1 && mxGetN(prhs[2])!=1)))
		mexErrMsgTxt("Input Rmx must be a real scalar!");

    image_row=mxGetM(prhs[0]);
    image_col=mxGetN(prhs[0]);
	input_image=(unsigned char*)mxGetLogicals(prhs[0]);

	dR=mxGetScalar(prhs[1]);
	if (nrhs>2){
		max_distance=mxGetScalar(prhs[2]);
		if (image_row<image_col && (image_row/2)<max_distance)
   			num_radial_bins=(int)ceil((double)image_row/(2.0*dR))-1;
		else if (image_col<image_col && (image_col/2)<max_distance)
   			num_radial_bins=(int)ceil((double)image_col/(2.0*dR))-1;
		else
   			num_radial_bins=(int)ceil((double)max_distance/dR)-1;
		}
	else{
       		if (image_row<image_col) num_radial_bins=(int)ceil((double)image_row/(2.0*dR))-1;
       		else                     num_radial_bins=(int)ceil((double)image_col/(2.0*dR))-1;
		}

    plhs[0]=mxCreateDoubleMatrix(1, num_radial_bins, mxREAL);
	pair_correlation=mxGetPr(plhs[0]);
	if (nlhs==2){
    	plhs[1]=mxCreateDoubleMatrix(1, num_radial_bins, mxREAL);
		radial_data=mxGetPr(plhs[1]);
		}
	else
		radial_data=(double*)mxCalloc(num_radial_bins, sizeof(double));

	cluster_list=NULL;
	num_clusters=0;
    if (nrhs>3){
        mxGetString(prhs[3], cluster_str, 5);
        if (cluster_str[0]==99 && cluster_str[1]==108 && cluster_str[2]==117 && cluster_str[3]==115)
    		num_clusters=get_all_clusters(input_image, 1);
		pair_correlation_function(pair_correlation, radial_data, input_image);
        return;
        }
   	num_clusters=get_all_pixels(input_image, 1);
    pair_correlation_function(pair_correlation, radial_data, input_image);
    return;
}
