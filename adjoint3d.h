#ifndef _adjoint3d_h

//set the token
#define _adjoint3d_h

#include "complex.h"
#include "matrix.h"
#include "mode.h"
#include "timefield3d.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

//define the structure and typedef it for ease of use
struct adjoint3d_
{
	int t; //length of time for storage
	int x_start, x_end, y_start, y_end, z_start, z_end; //limits of rectangular prism for storing fields
	double *forward_fields, *gradient; //matrix containing electromagnetic field components at each volume cell. Ordered as Ex[0,0], Ey[0,0]...
};
typedef struct adjoint3d_  adjoint3d;

adjoint3d ajoint3d_new(int x_start, int x_end, int y_start, int y_end, int z_start, int z_end, int t);
void adjoint3d_free(adjoint3d block);

//gradient calculation
void adjoint3d_store_fwd_fields(adjoint3d block, timefield3d main_fields, int timestep); //to be used at each forward timestep to store fields
void adjoint3d_update_gradient(adjoint3d block, timefield3d main_fields, int timestep); //to be used at each reverse (adjoint) timestep to accumulate gradient
void adjoint3d_reset_gradient(adjoint3d block); //to be used before reverse (ajoint) simulation to zero out gradient
void adjoint3d_reset_forward_fields(adjoint3d block); //to be use before running subsequent forward timestep

//structure adjustment
void adjoint3d_adjust_permittivity(adjoint3d block, timefield3d main_fields, double omega, double improvement, double ep_min, double ep_max);

//tapers
struct adjoint3d_taper_
{
	double ep; //real, isotropic permittivity
	double mu; //real, isotropic permeability
	int seg_num; //number of inidvidual linear taper segments
	int start; //start (z) of taper (in pixels)
	int end; //end (z) of taper (in pixels)
	int bottom; //bottom (y)
	int top; //top (y)
	int duration; //final time-step

	int *left; //seg_num+1 points defining the left edge
	int *right; //seg_num+1 points defining the right edge

	double *forward_fields; //forward-run stored fields along edges 
	double *adjoint_fields; //adjoint-run stored fields along edges 
	double *ep_gradient; //gradient with respect to ep/mu along edges

	double *left_grad; //gradient for adjustment
	double *right_grad; //gradient for adjustment
};
typedef struct adjoint3d_taper_  adjoint3d_taper;

adjoint3d_taper adjoint3d_taper_new(double ep, double mu, int seg_num, int start, int end, int *left, int *right, int bottom, int top, int duration);
void adjoint3d_taper_insert(timefield3d psi, adjoint3d_taper taper);
void adjoint3d_taper_store_forward_fields(adjoint3d_taper taper, timefield3d main_fields, int timestep);
void adjoint3d_taper_store_adjoint_fields(adjoint3d_taper taper, timefield3d main_fields, int timestep);
void adjoint3d_taper_calculate_ep_gradient(adjoint3d_taper taper);
void adjoint3d_taper_calculate_gradients_from_ep_gradient(adjoint3d_taper taper, double ep_diff);

void adjoint3d_taper_calculate_gradients(adjoint3d_taper taper, adjoint3d block, double ep_diff);
int adjoint3d_taper_update_symmetric(adjoint3d_taper taper, double improvement_goal, double disc, int min_width);
int adjoint3d_taper_update_asymmetric(adjoint3d_taper taper, double improvement_goal, double disc, int min_width);

#endif