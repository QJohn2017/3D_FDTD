#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "matrix.h"
#include "complex.h"
#include "mode.h"
#include "timefield3d.h"
#include "adjoint3d.h"

#define PI 3.141592652

int main(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected, *unwanted;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode, unwanted_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 100;
	y_size = 100;
	z_size = 600;
	t_length = 3501;

	disc = 0.02;

	omega = 2*PI/1.310;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));
	unwanted = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/elara_3d_360ridge_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/elara_3d_720ridge_TE01.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);


	mode_data = fopen("mode_files/elara_3d_720ridge_TE00.txt", "r");
	unwanted_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(unwanted_mode);

	//build
	seg_num = 25;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));
	widths = malloc((seg_num+1)*sizeof(int));

	timefield3d_set_background(main_field, 1.447*1.447, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 64); //SiCOH

	timefield3d_set_taper(main_field, 3.51*3.51, 1.0, 50, 50, 41, 59, 0, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.51*3.51, 41, 59, 95, 150, 46, 54);

	//setup and insert taper
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 41 - i*9/seg_num;
		taper_right[i] = 59 + i*9/seg_num;
	}
	taper = adjoint3d_taper_new(3.51*3.51, 1.0, seg_num, 150, 400, taper_left, taper_right, 46, 54, t_length);
	adjoint3d_taper_insert(main_field, taper);

	timefield3d_set_waveguide(main_field, 3.51*3.51, 32, 68, 400, 455, 46, 54);
	timefield3d_set_taper(main_field, 3.51*3.51, 1.0, 32, 68, 50, 50, 455, 600, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 50);

	no_change_flag = 0;
	for(iteration = 0; iteration < 10; iteration++){

		mode_data = fopen("mode_files/elara_3d_720ridge_TE01.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		*unwanted = complex_zero();
		timefield3d_zero_fields(main_field);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 300)*(i - 300)/(2*150.0*150.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, x_size/2, y_size/2, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);
			timefield3d_update_fourier_fields_y_slice_omp(main_field, phase, y_size/2);

			instant_overlap = timefield3d_mode_overlap(main_field, x_size/2, y_size/2, 130, launch_mode, phase, 0, launched);
			instant_overlap = timefield3d_mode_overlap(main_field, x_size/2, y_size/2, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, x_size/2, y_size/2, 425, output_mode, phase, 0, transmitted);
			instant_overlap = timefield3d_mode_overlap(main_field, x_size/2, y_size/2, 425, unwanted_mode, phase, 0, unwanted);

			adjoint3d_taper_store_forward_fields(taper, main_field, i);

			if(i%100 == 0){

				//timefield3d_plot_Ex_x_slice(main_field, (y_size/2));
				//timefield3d_plot_Ex_y_slice_fourier(main_field, (y_size/2));

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\tunwanted field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected), complex_mag(*unwanted));
			}
		}

		timefield3d_plot_Ex_y_slice_fourier(main_field, (y_size/2));

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, x_size/2, y_size/2, 425, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_taper_store_adjoint_fields(taper, main_field, i);

			if(i%100 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//calculate gradients
		adjoint3d_taper_calculate_ep_gradient(taper);
		adjoint3d_taper_calculate_gradients_from_ep_gradient(taper, 3.51*3.51-1.447*1.447);

		//adjust taper
		improvement = 0.2*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		no_change_flag = adjoint3d_taper_update_asymmetric(taper, improvement, main_field.disc, 18);

		timefield3d_set_background(main_field, 1.447*1.447, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 64); //SiCOH
		timefield3d_set_taper(main_field, 3.51*3.51, 1.0, 50, 50, 41, 59, 0, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.51*3.51, 41, 59, 95, 150, 46, 54);
		adjoint3d_taper_insert(main_field, taper);
		timefield3d_set_waveguide(main_field, 3.51*3.51, 32, 68, 400, 455, 46, 54);
		timefield3d_set_taper(main_field, 3.51*3.51, 1.0, 32, 68, 50, 50, 455, 600, 46, 54);

		timefield3d_plot_permittivity_y_slice(main_field, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 50);

		printf("fwd taper edges:\n");	
		for(i = 0; i < seg_num+1; i++){

			printf("%i\t%i\t%i\n", i, (int)((taper.right[i]-x_size/2)*1000*disc), (int)((taper.left[i]-x_size/2)*1000*disc));
		}

		//free stuff
		mode_free(output_mode);

		if(no_change_flag == 1){

			break;
		}
	}

	printf("\nTaper Length = %i nm\n", 20*z_size);
	printf("Segment number = %i\n", seg_num);
	printf("Taper widths:\n");
	for(i = 0; i < seg_num+1; i++){

		printf("Width %i = %i nm\n", i, 20*(taper.right[i]-taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "elara_TE01_coupler_1.bin", 20, 80, 150, 400, 50, 1.447*1.447, 3.51*3.51);
}