#ifndef _timefield3d_h

//set the token
#define _timefield3d_h

#include "complex.h"
#include "matrix.h"
#include "mode.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

//define the structure and typedef it for ease of use
struct timefield3d_
{
	int i, j, k; //number of area cells in x, y, z directions
	double disc; //discretization in µm
	float *ep; //pointer to array of permittivity values at each volume cell. 
	float *mu; //pointer to array of permeability values at each volume cell.
	float *sig; //pointer to array of conductivity (and magnetic conductivity) values at each volume cell.
	double *fields; //matrix containing electromagnetic field components at each volume cell. Ordered as Ex[0,0], Ey[0,0]...
	double *old_fields; //matrix for storing previous timestep
	complex *fourier_fields; //matrix for storing fourier transform of fields
	float *sources; //matrix containing all sources (current density and magnetic current density)
};
typedef struct timefield3d_  timefield3d;

//initialing, setting and getting fields
timefield3d timefield3d_new(int i, int j, int k, double disc); //returns an timefield3d object with i*j*k area cells, discretization set to disc and all fields set to 0.
void timefield3d_free(timefield3d psi); //frees allocated memory associated with psi
timefield3d timefield3d_copy(timefield3d psi); //returns a copy of timefield3d structure psi
void timefield3d_zero_fields(timefield3d psi); //zeroes out all fields, leaves structure unchanged
double* timefield3d_get_fields(timefield3d psi); //returns a copy of the field array of psi
double* timefield3d_get_edge_fields(timefield3d psi); //returns an array with only the border (2 pixels) of the fields of psi copied
void timefield3d_update_old_edge_fields(timefield3d psi);

//plotting functions
void timefield3d_plot_Ex_x_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ex_y_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ex_z_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ey_x_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ey_y_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ey_z_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ez_x_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ez_y_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ez_z_slice(timefield3d psi, int slice_index);
void timefield3d_plot_Ex_x_slice_fourier(timefield3d psi, int slice_index);
void timefield3d_plot_Ex_y_slice_fourier(timefield3d psi, int slice_index);
void timefield3d_plot_Ex_z_slice_fourier(timefield3d psi, int slice_index);
void timefield3d_plot_permittivity_x_slice(timefield3d psi, int slice_index);
void timefield3d_plot_permittivity_y_slice(timefield3d psi, int slice_index);
void timefield3d_plot_permittivity_z_slice(timefield3d psi, int slice_index);

//setting dielectric structures
void timefield3d_set_background(timefield3d psi, float ep, float mu);
void timefield3d_set_top_half_plane(timefield3d psi, float ep, float mu, int bottom);
void timefield3d_set_background_permittivity(timefield3d psi, float ep);
void timefield3d_set_background_permeability(timefield3d psi, float mu);
void timefield3d_set_background_electric_conductivity(timefield3d psi, float sig);
void timefield3d_set_background_magnetic_conductivity(timefield3d psi, float sig);
void timefield3d_set_background_complex_permittivity(timefield3d psi, double omega, float real_ep, float imag_ep);

void timefield3d_set_prism_permittivity(timefield3d psi, float ep, int start_x, int end_x, int start_y, int end_y, int start_z, int end_z);
void timefield3d_set_prism_electric_conductivity(timefield3d psi, float sig, int start_x, int end_x, int start_y, int end_y, int start_z, int end_z);
void timefield3d_set_prism_complex_permittivity(timefield3d psi, double omega, float real_ep, float imag_ep, int start_x, int end_x, int start_y, int end_y, int start_z, int end_z);

//setting waveguide structures
void timefield3d_set_waveguide(timefield3d psi, float ep, int start_x, int end_x, int start_z, int end_z, int bottom, int top);
void timefield3d_set_taper(timefield3d psi, double ep, double mu, int start_left, int start_right, int end_left, int end_right, int start_z, int end_z, int bottom, int top);

//FDTD functions
void timefield3d_apply_H_timestep(timefield3d psi, double delta_t); //update H_fields using curl of E and magnetic current (source and induced)
void timefield3d_apply_E_timestep(timefield3d psi, double delta_t); //update E_fields using cul of H and elecric current (source and induced)
void timefield3d_apply_H_timestep_zsection(timefield3d psi, double delta_t, int start_index, int stop_index); //update H_fields only between z = start and z = stop
void timefield3d_apply_E_timestep_zsection(timefield3d psi, double delta_t, int start_index, int stop_index); //update E_fields only between z = start and z = stop
void timefield3d_apply_H_first_order_abc(timefield3d psi, double *old_field, double delta_t); //sets H on boundaries to absorb normal-incidence plane waves
void timefield3d_apply_E_first_order_abc(timefield3d psi, double *old_field, double delta_t); //sets E on boundaries to absorb normal-incidence plane waves
void timefield3d_apply_H_first_order_abc_zsection(timefield3d psi, double delta_t, int start_index, int stop_index); 
void timefield3d_apply_E_first_order_abc_zsection(timefield3d psi, double delta_t, int start_index, int stop_index); 
void timefield3d_apply_timestep(timefield3d psi, double delta_t); //propagates fields one full time step (uses first order a.b.c.)
void timefield3d_apply_timestep_omp(timefield3d psi, double delta_t); //timestep using openMP threading (default threads on given machine)
void timefield3d_apply_timestep_omp_full(timefield3d psi, double delta_t);
void timefield3d_update_fourier_fields(timefield3d psi, double phase);
void timefield3d_update_fourier_fields_y_slice(timefield3d psi, double phase, int slice_index);
void timefield3d_update_fourier_fields_zsection(timefield3d psi, double phase, int start_index, int stop_index);
void timefield3d_update_fourier_fields_y_slice_zsection(timefield3d psi, double phase, int slice_index, int start_index, int stop_index);
void timefield3d_update_fourier_fields_omp(timefield3d psi, double phase);
void timefield3d_update_fourier_fields_y_slice_omp(timefield3d psi, double phase, int slice_index);

void timefield3d_apply_E_timestep_use_polarization(timefield3d psi, double *old_field, double delta_t); //uses material polarization instead of permittivity
void timefield3d_apply_timestep_use_polarization(timefield3d psi, double *old_old_field, double delta_t); //uses material polarization instead of permittivity

//source functions
void timefield3d_mode_source(timefield3d psi, int center_x, int center_y, int center_z, mode source_mode, double phase, double amplitude, int reverse);
void timefield3d_mode_source_adjoint(timefield3d psi, int center_x, int center_y, int center_z, mode source_mode, double phase, int reverse);
void timefield3d_TE_gaussian_source(timefield3d psi, int center_x, int center_y, int center_z, double width, double omega, double phase, double amplitude);
void timefield3d_TM_gaussian_source(timefield3d psi, int center_x, int center_y, int center_z, double width, double omega, double phase, double amplitude);

//overlap functions
complex timefield3d_mode_overlap(timefield3d psi, int center_x, int center_y, int center_z, mode overlap_mode, double phase, int reverse, complex *cumulative_overlap);

//structure output functions
void timefield3d_file_write_permittivity_y_slice(timefield3d psi, char *filename, int left, int right, int start, int end, int y_slice, double ep_min, double ep_max);

#endif