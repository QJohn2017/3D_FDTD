#ifndef _mode_h

//set the token
#define _mode_h

#include "complex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//define the structure and typedef it for ease of use
struct mode_
{
	int i, j; //number of area cells in x, y
	double disc; //discretization in µm
	matrix pml; //matrix containing value of sx, sy for pml implementation.
	matrix ep; //pointer to array of permittivity values at each volume cell. may want to change this to complex in the future.
	matrix mu; //pointer to array of permeability values at each volume cell. may want to change this to complex in the future.
	matrix fields; //matrix containing electromagnetic field components at each volume cell. Ordered as Ex[0,0], Ey[0,0]...
};
typedef struct mode_  mode;

mode mode_new(int i, int j, double disc); //returns an mode object with i*j area cells, discretization set to disc and all fields set to 0.
void mode_free(mode psi); //frees allocated memory associated with psi
mode mode_copy(mode psi); //returns copy of psi stored in separate location in memory
void mode_copy_fields(mode psi, mode phi); //copies fields from phi into psi
void mode_smul(complex scalar, mode psi); //scalar multiplies all fields
mode mode_linear_combination(mode psi, mode phi, complex scalar); //returns new mode as complex linear combination psi + scalar*phi;
mode mode_linear_combination_3(mode psi, mode phi, mode alpha, complex scalar_phi, complex scalar_alpha);//returns new mode as complex linear combination psi + scalar_phi*phi + scalar_alpha*alpha;
mode mode_gaussian(int i, int j, double disc, complex ep, complex mu, double omega, double sigma, int polarization);//makes gaussian mode with sigma, centered in rectangle ixj, and with horizontal polarization 0, and vertical 1
void mode_chop_bottom_half(mode psi);

void mode_plot(mode psi); //uses gnuplot to plot image of index distribution and field distribution
void mode_plot_real(mode psi); //plots real part of field
void mode_plot_ex(mode psi); //uses gnuplot to plot Ex
void mode_plot_hy(mode psi); //uses gnuplot to plot Hy

void mode_file_print(FILE *mode_data, mode psi, int dimension); //prints all mode data to a file. first line is structure data (i, j and disc), subsequent lines are field data
void mode_angular_file_print(FILE *mode_data, mode psi, int dimension);
mode mode_file_load(FILE *mode_data); //returns mode with dimensions and field values specified in file
mode mode_angular_file_load(FILE *mode_data);

double mode_max_index(mode psi); //find the maximum real part of index in field psi

void mode_initialize_ones(mode psi); //initialize all field values to complex 1
void mode_initialize_rand(mode psi); //initialize fields to random complex values

complex mode_inner_product(mode psi, mode phi); //takes inner product as initegral of psi*x A x phi (A is as defined in problem set 5)
complex mode_H_inner_product(mode psi, mode phi); //takes inner product as integral of H_psi* x H_phi
complex mode_power_inner_product(mode psi, mode phi); //takes inner product as integral of z_component of H_psi* X E_phi + H_phi X E_psi* 
double mode_norm(mode psi); //returns norm of psi in inner product space
double mode_H_norm(mode psi); //returns norm of psi in H inner product space
double mode_power_norm(mode psi); //returns norm of psi in power inner product space
void mode_normalize(mode psi); //normalizes (for now set it so that the average value of mode component is 1)
void mode_H_normalize(mode psi); //normalizes transverse H-field
void mode_power_normalize(mode psi); //normalized fields so that 1 unit of power is carried

void mode_apply_zero_bc(mode psi); //sets all fields along the boundaries to zero
void mode_apply_zero_bc_vertical(mode psi); //sets all fields along the vertical (top and bottom) boundaries to zero
void mode_apply_zero_bc_horizontal(mode psi); //sets all fields along the horizontal (left and right) boundaries to zero
void mode_apply_zero_gradient_bc(mode psi); //sets gradients of all fields to zero along the boundary
void mode_apply_zero_gradient_bc_vertical(mode psi); //sets gradients of all fields to zero along the vertical (top and bottom) boundary
void mode_apply_zero_gradient_bc_horizontal(mode psi); //sets gradients of all fields to zero along the horizontal (left and right) boundary
void mode_apply_bc(mode psi, int confinement); //applies appropriate bc for 2d confinement or 1d confinement
void mode_apply_zero_bc_pml(mode psi);//applies zero field where imaginary part of s is > 1e-9

void mode_set_pml(mode psi, double max_ep_im, int order, int pml_thickness); //sets pml around outer thickness, with a polynomial of order 'order', starting from 0 ep_im and ending at max_ep_im
void mode_set_pml_bottom(mode psi, double max_ep_im, int order, int pml_thickness); //sets pml on bottom thickness, with a polynomial of order 'order', starting from 0 ep_im and ending at max_ep_im
void mode_set_pml_top(mode psi, double max_ep_im, int order, int pml_thickness); //sets pml on top thickness, with a polynomial of order 'order', starting from 0 ep_im and ending at max_ep_im
void mode_set_pml_left(mode psi, double max_ep_im, int order, int pml_thickness); //sets pml on left thickness, with a polynomial of order 'order', starting from 0 ep_im and ending at max_ep_im
void mode_set_pml_right(mode psi, double max_ep_im, int order, int pml_thickness); //sets pml on right thickness, with a polynomial of order 'order', starting from 0 ep_im and ending at max_ep_im
void mode_set_grade_to_constant_pml(mode psi, double max_ep_im, int grade_thickness, int pml_thickness); //sets pml around outer thickness, quadratic grade over grade thickness and constant max for the rest

void mode_set_background(mode psi, complex ep, complex mu); //sets entire mode area to material defined by ep, mu
void mode_set_bottom_half_plane(mode psi, complex ep, complex mu, int top); //sets half-plane below top to material defined by ep, mu
void mode_set_top_half_plane(mode psi, complex ep, complex mu, int bottom); //sets half-plane above bottom to material defined by ep, mu
void mode_set_left_half_plane(mode psi, complex ep, complex mu, int right); //sets half-plane left of right to material defined by ep, mu
void mode_set_right_half_plane(mode psi, complex ep, complex mu, int left); //sets half-plane right of left to material defined by ep, mu
void mode_make_rectangle(mode psi, complex ep, complex mu, int left, int right, int bottom, int top); //creates rectangle defined by sides with of material defined by ep, mu

void mode_apply_H_curl(mode psi, complex beta, double omega); //replaces E-field with (1/(i*omega*ep))*curl(H) in place, given that the z-dependence is (i*beta*z)
void mode_apply_E_curl(mode psi, complex beta, double omega); //replaces H-field with (1/(i*omega*mu))*curl(E) in place, given that the z-dependence is (i*beta*z)
void mode_apply_radial_H_curl(mode psi, complex beta, double omega, double radius_left); //solves for E fields using H-curl, using radial coordinates
void mode_apply_omega_squared_operator(mode psi); //applies omega-squared operator with zero boundary conditions 

mode mode_apply_beta_operator(mode psi, double omega); //returns the result of applying the beta operator to psi. assumes psi is linear guided mode psi=Psi(x,y)*exp(i(omega*t-beta*z)) 
mode mode_apply_radial_beta_squared_operator(mode psi, double omega, double radius_left);

void mode_apply_H_div(mode psi, complex beta);
void mode_apply_radial_H_div(mode psi, complex beta, double radius_left);

complex mode_solve_beta(mode *psi, double omega, double tol, mode *prev_modes, int prev_num, int size); //solve for complex beta and set field solutions (ignore spurious values of beta > n*omega)
complex mode_solve_beta_squared(mode *psi, double omega, double tol, mode *prev_modes, int prev_num, int size, int confinement, int lateral_symmetry); //solve for complex beta and set field solutions (ignore spurious values of beta > n*omega)
complex mode_solve_radial_beta_squared(mode *psi, double omega, double tol, mode *prev_modes, int prev_num, int size, int confinement, double radius_left); //solve for complex angular beta and set field solutions for radial mode with inner radius radius_left

mode mode_apply_beta_operator_cylindrical(mode psi, double omega); //returns the result of applying the beta operator to psi in cylindrical coordinates

double mode_arnoldi_eigensolver(mode psi); //returns the eigenvalue associated with the eigenmode psi

complex mode_true_beta(complex beta); //returns true beta in rad/m instead of rad/µm
mode mode_true_fields(mode psi); //returns true fields in V/m and A/m normalized to 1 Watt

double mode_n_susceptibility_rectangle(mode psi, int left, int right, int bottom, int top); //calculates susceptibility of mode to changes within rectangle
double mode_n_susceptibility_rectangle_exclude_rectangle(mode psi, int left, int right, int bottom, int top, int left_exclude, int right_exclude, int bottom_exclude, int top_exclude);
double mode_n_susceptibility_rectangle_top(mode psi, int left, int right, int bottom, int top); //calculates susceptibility of mode to changes on top, left and right edges of rectangle
double mode_n_susceptibility_rectangle_bottom(mode psi, int left, int right, int bottom, int top);
double mode_n_susceptibility_rectangle_sides(mode psi, int left, int right, int bottom, int top);//both inside and ouside layer
double mode_n_susceptibility_rectangle_sides_outside(mode psi, int left, int right, int bottom, int top);//only outside layer
double mode_n_susceptibility_rectangle_sides_normal(mode psi, int left, int right, int bottom, int top); //Ex^2 only
double mode_n_susceptibility_rectangle_sides_tangential(mode psi, int left, int right, int bottom, int top); //Ey^2 + Ez^2
double mode_n_susceptibility_rectangle_sides_longitudinal(mode psi, int left, int right, int bottom, int top);//Ez^2 only
double mode_n_susceptibility_rectangle_edge(mode psi, int left, int right, int bottom, int top); //calculates susceptibility of mode to changes on edges of rectangle
double mode_n_susceptibility_bottom_half_plane(mode psi, int top);

double mode_loss_susceptibility_rectangle(mode psi, int left, int right, int bottom, int top); //calulates ratio of effective loss over material loss for rectangle (2*n_material*index_susceptibility)
double mode_loss_susceptibility_rectangle_top(mode psi, int left, int right, int bottom, int top); //for top surface, assumes 1 pixel thick loss on either side of interface
double mode_loss_susceptibility_rectangle_bottom(mode psi, int left, int right, int bottom, int top); //for bottom surface, assumes 1 pixel thick loss on either side of interface
double mode_loss_susceptibility_rectangle_sides(mode psi, int left, int right, int bottom, int top); //for lateral surfaces, assumes 1 pixel thick loss on either side of interface
double mode_loss_susceptibility_bottom_half_plane(mode psi, int top); //for bottom half-planes

void mode_apply_lateral_symmetry(mode psi); //returns mode forced to be either symmetric or anti-symmetric (whichever is closer)
#endif