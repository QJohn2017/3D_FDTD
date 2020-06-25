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

int main_test_taper(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 100;
	y_size = 50;
	z_size = 300;
	t_length = 1601;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/3d_ridge_400nm_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/3d_ridge_800nm_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(25, 75, 21, 29, 150, 175, t_length);

	//build
	seg_num = 10;
	taper_left = malloc(seg_num*sizeof(int));
	taper_right = malloc(seg_num*sizeof(int));

	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 50, 50, 40, 60, 25, 100, 21, 29);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 40, 60, 100, 150, 21, 29);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 40 - i*10/seg_num;
		taper_right[i] = 60 + i*10/seg_num;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 150, 175, taper_left, taper_right, 21, 29);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_waveguide(main_field, 3.476*3.476, 30, 70, 175, 300, 21, 29);
	//timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 30, 70, 50, 50, 250, 375, 21, 29);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);

	no_change_flag = 0;
	for(iteration = 0; iteration < 10; iteration++){

		mode_data = fopen("mode_files/3d_ridge_800nm_TE00.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 300)*(i - 300)/(2*150.0*150.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 50, 25, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 50, 25, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 50, 25, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 50, 25, 200, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%200 == 0){

				//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
				timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		timefield3d_zero_fields(main_field);
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 50, 25, 200, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%200 == 0){

				timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 20);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 50, 50, 40, 60, 25, 100, 21, 29);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 40, 60, 100, 150, 21, 29);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 30, 70, 175, 300, 21, 29);
		//timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 30, 70, 50, 50, 250, 375, 21, 29);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);

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

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}
}

int main_y_juntion(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 150;
	y_size = 100;
	z_size = 500;
	t_length = 3201;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/brain_y_junction_in_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(20, 130, 46, 54, 150, 300, t_length);

	//build
	seg_num = 10;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));
	widths = malloc((seg_num+1)*sizeof(int));

	widths[0] = 24;
	widths[1] = 28;
	widths[2] = 36;
	widths[3] = 44;
	widths[4] = 38;
	widths[5] = 32;
	widths[6] = 54;
	widths[7] = 64;
	widths[8] = 68;
	widths[9] = 70;
	widths[10] = 74;

	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 75 - widths[i]/2;
		taper_right[i] = 75 + widths[i]/2;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 150, 300, taper_left, taper_right, 46, 54);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 87+25, 63+25, 87+25, 300, 350, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 350, 405, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 405, 475, 46, 54);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 97-25, 63-25, 87-25, 300, 350, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 350, 405, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 405, 475, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

	no_change_flag = 0;
	for(iteration = 0; iteration < 10; iteration++){

		mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_forward_fields(opt_area);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 300)*(i - 300)/(2*150.0*150.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 375, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%100 == 0){

				//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 75, 50, 375, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%100 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.2*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 20);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 87+25, 63+25, 87+25, 300, 350, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 350, 405, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 405, 475, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 97-25, 63-25, 87-25, 300, 350, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 350, 405, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 405, 475, 46, 54);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

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

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "brain_y_junction_design_4.bin", 20, 130, 150, 350, 50, 1.444*1.444, 3.476*3.476);
}

int main_y_junction_design_4(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 150;
	y_size = 100;
	z_size = 600;
	t_length = 4001;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/brain_y_junction_in_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(20, 130, 46, 54, 150, 350, t_length);

	//build
	seg_num = 20;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));
	widths = malloc((seg_num+1)*sizeof(int));

	widths[0] = 24;
	widths[1] = 22;
	widths[2] = 24;
	widths[3] = 22;
	widths[4] = 26;
	widths[5] = 30;
	widths[6] = 34;
	widths[7] = 36;
	widths[8] = 40;
	widths[9] = 30;
	widths[10] = 34;
	widths[11] = 24;
	widths[12] = 50;
	widths[13] = 56;
	widths[14] = 72;
	widths[15] = 72;
	widths[16] = 72;
	widths[17] = 70;
	widths[18] = 72;
	widths[19] = 72;
	widths[20] = 74;

	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 75 - widths[i]/2;
		taper_right[i] = 75 + widths[i]/2;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 150, 350, taper_left, taper_right, 46, 54);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 87+25, 63+25, 87+25, 350, 450, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 450, 505, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 505, 575, 46, 54);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 97-25, 63-25, 87-25, 350, 450, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 450, 505, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 505, 575, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

	no_change_flag = 0;
	for(iteration = 0; iteration < 2; iteration++){

		mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_forward_fields(opt_area);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 600)*(i - 600)/(2*200.0*200.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 475, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 75, 50, 475, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.5*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 20);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 87+25, 63+25, 87+25, 350, 450, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 450, 505, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 505, 575, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 97-25, 63-25, 87-25, 350, 450, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 450, 505, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 505, 575, 46, 54);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

		//free stuff
		mode_free(output_mode);

		if(no_change_flag == 1){

			break;
		}
	}

	printf("\nTaper Length = %i nm\n", 20*200);
	printf("Segment number = %i\n", seg_num);
	printf("Taper widths:\n");
	for(i = 0; i < seg_num+1; i++){

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "brain_y_junction_design_4.bin", 20, 130, 150, 450, 50, 1.444*1.444, 3.476*3.476);
}

int main_y_junction_design_5(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 150;
	y_size = 100;
	z_size = 700;
	t_length = 4801;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/brain_y_junction_in_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(20, 130, 46, 54, 150, 400, t_length);

	//build
	seg_num = 20;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));
	widths = malloc((seg_num+1)*sizeof(int));

	widths[0] = 24;
	widths[1] = 20;
	widths[2] = 20;
	widths[3] = 20;
	widths[4] = 22;
	widths[5] = 40;
	widths[6] = 50;
	widths[7] = 52;
	widths[8] = 46;
	widths[9] = 44;
	widths[10] = 36;
	widths[11] = 48;
	widths[12] = 66;
	widths[13] = 74;
	widths[14] = 78;
	widths[15] = 80;
	widths[16] = 84;
	widths[17] = 88;
	widths[18] = 86;
	widths[19] = 82;
	widths[20] = 94;

	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 75 - widths[i]/2;
		taper_right[i] = 75 + widths[i]/2;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 150, 400, taper_left, taper_right, 46, 54);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 97+25, 63+25, 87+25, 400, 500, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 500, 605, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 605, 675, 46, 54);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53-25, 97-25, 63-25, 87-25, 400, 500, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 500, 605, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 605, 675, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

	no_change_flag = 0;
	for(iteration = 0; iteration < 10; iteration++){

		mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_forward_fields(opt_area);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 600)*(i - 600)/(2*200.0*200.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 575, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%400 == 0){

				if(iteration == 0){

					//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
					//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));
				}

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 75, 50, 575, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.1*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 20);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 97+25, 63+25, 87+25, 400, 500, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 500, 605, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 605, 675, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53-25, 97-25, 63-25, 87-25, 400, 500, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 500, 605, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 605, 675, 46, 54);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

		//free stuff
		mode_free(output_mode);

		if(no_change_flag == 1){

			break;
		}
	}

	printf("\nTaper Length = %i nm\n", 20*200);
	printf("Segment number = %i\n", seg_num);
	printf("Taper widths:\n");
	for(i = 0; i < seg_num+1; i++){

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "brain_y_junction_design_5.bin", 20, 130, 150, 500, 50, 1.444*1.444, 3.476*3.476);
}

int main_y_from_scratch_5um_full_width(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 150;
	y_size = 100;
	z_size = 700;
	t_length = 4801;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/brain_y_junction_in_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(20, 130, 46, 54, 150, 400, t_length);

	//build
	seg_num = 20;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));
	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 63 - i*35/seg_num;
		taper_right[i] = 87 + i*35/seg_num;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 150, 400, taper_left, taper_right, 46, 54);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 97+25, 63+25, 87+25, 400, 500, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 500, 605, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 605, 675, 46, 54);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53-25, 97-25, 63-25, 87-25, 400, 500, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 500, 605, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 605, 675, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

	no_change_flag = 0;
	for(iteration = 0; iteration < 30; iteration++){

		mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_forward_fields(opt_area);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 600)*(i - 600)/(2*200.0*200.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 575, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%400 == 0){

				if(iteration == 0){

					//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
					//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));
				}

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 75, 50, 575, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.1*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 20);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 97+25, 63+25, 87+25, 400, 500, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 500, 605, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 605, 675, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53-25, 97-25, 63-25, 87-25, 400, 500, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 500, 605, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 605, 675, 46, 54);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

		//free stuff
		mode_free(output_mode);

		if(no_change_flag == 1){

			break;
		}
	}

	printf("\nTaper Length = %i nm\n", 20*200);
	printf("Segment number = %i\n", seg_num);
	printf("Taper widths:\n");
	for(i = 0; i < seg_num+1; i++){

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "brain_y_junction_design_5.bin", 20, 130, 150, 500, 50, 1.444*1.444, 3.476*3.476);
}

int main_y_from_scratch_4um_partial_width(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 150;
	y_size = 100;
	z_size = 600;
	t_length = 4001;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/brain_y_junction_in_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(20, 130, 46, 54, 150, 350, t_length);

	//build
	seg_num = 20;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));

	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 63 - i*25/seg_num;
		taper_right[i] = 87 + i*25/seg_num;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 150, 350, taper_left, taper_right, 46, 54);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 87+25, 63+25, 87+25, 350, 450, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 450, 505, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 505, 575, 46, 54);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 97-25, 63-25, 87-25, 350, 450, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 450, 505, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 505, 575, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

	no_change_flag = 0;
	for(iteration = 0; iteration < 20; iteration++){

		printf("\nIteration = %i\n", iteration);

		mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_forward_fields(opt_area);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 600)*(i - 600)/(2*200.0*200.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 475, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 75, 50, 475, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.05*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 20);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 87+25, 63+25, 87+25, 350, 450, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 450, 505, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 505, 575, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 97-25, 63-25, 87-25, 350, 450, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 450, 505, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 505, 575, 46, 54);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

		//free stuff
		mode_free(output_mode);

		if(no_change_flag == 1){

			break;
		}
	}

	printf("\nTaper Length = %i nm\n", 20*200);
	printf("Segment number = %i\n", seg_num);
	printf("Taper widths:\n");
	for(i = 0; i < seg_num+1; i++){

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "brain_y_junction_design_4.bin", 20, 130, 150, 450, 50, 1.444*1.444, 3.476*3.476);
}

int main_y_from_scratch_3um_full_width_extra_space(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 150;
	y_size = 100;
	z_size = 600;
	t_length = 4001;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/brain_y_junction_in_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/brain_y_junction_out_TE00_1400nm_pitch.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(20, 130, 46, 54, 150, 300, t_length);

	//build
	seg_num = 20;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));

	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 63 - i*45/seg_num;
		taper_right[i] = 87 + i*45/seg_num;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 150, 300, taper_left, taper_right, 46, 54);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+35, 97+35, 63+35, 87+35, 300, 400, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63+35, 87+35, 400, 505, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+35, 87+35, 110, 110, 505, 575, 46, 54);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53-35, 97-35, 63-35, 87-35, 300, 400, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63-35, 87-35, 400, 505, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-35, 87-35, 40, 40, 505, 575, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

	no_change_flag = 0;
	for(iteration = 0; iteration < 20; iteration++){

		printf("\nIteration = %i\n", iteration);

		mode_data = fopen("mode_files/brain_y_junction_out_TE00_1400nm_pitch.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_forward_fields(opt_area);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 600)*(i - 600)/(2*200.0*200.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 475, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 75, 50, 475, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.15*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 20);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+35, 97+35, 63+35, 87+35, 300, 400, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63+35, 87+35, 400, 505, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+35, 87+35, 100, 100, 505, 575, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53-35, 97-35, 63-35, 87-35, 300, 400, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63-35, 87-35, 400, 505, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-35, 87-35, 50, 50, 505, 575, 46, 54);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

		//free stuff
		mode_free(output_mode);

		if(no_change_flag == 1){

			break;
		}
	}

	printf("\nTaper Length = %i nm\n", 20*200);
	printf("Segment number = %i\n", seg_num);
	printf("Taper widths:\n");
	for(i = 0; i < seg_num+1; i++){

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "brain_y_junction_design_4.bin", 20, 130, 150, 400, 50, 1.444*1.444, 3.476*3.476);
}

int main_y_from_scratch_wider_input(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 150;
	y_size = 100;
	z_size = 700;
	t_length = 4801;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/brain_y_junction_in_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(20, 130, 46, 54, 200, 400, t_length);

	//build
	seg_num = 20;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));
	widths = malloc((seg_num+1)*sizeof(int));

	widths[0] = 30;
	widths[1] = 32;
	widths[2] = 36;
	widths[3] = 38;
	widths[4] = 42;
	widths[5] = 44;
	widths[6] = 40;
	widths[7] = 34;
	widths[8] = 30;
	widths[9] = 58;
	widths[10] = 62;
	widths[11] = 64;
	widths[12] = 80;
	widths[13] = 76;
	widths[14] = 86;
	widths[15] = 84;
	widths[16] = 88;
	widths[17] = 88;
	widths[18] = 86;
	widths[19] = 90;
	widths[20] = 94;

	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63, 87, 60, 90, 150, 200, 46, 54);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 75 - widths[i]/2;
		taper_right[i] = 75 + widths[i]/2;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 200, 400, taper_left, taper_right, 46, 54);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 97+25, 63+25, 87+25, 400, 500, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 500, 605, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 605, 675, 46, 54);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53-25, 97-25, 63-25, 87-25, 400, 500, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 500, 605, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 605, 675, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

	no_change_flag = 0;
	for(iteration = 0; iteration < 10; iteration++){

		mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_forward_fields(opt_area);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 600)*(i - 600)/(2*200.0*200.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 575, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%400 == 0){

				if(iteration == 0){

					//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
					//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));
				}

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 75, 50, 575, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.1*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 30);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 25, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63, 87, 60, 90, 150, 200, 46, 54);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 97+25, 63+25, 87+25, 400, 500, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 500, 605, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 605, 675, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53-25, 97-25, 63-25, 87-25, 400, 500, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 500, 605, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 605, 675, 46, 54);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

		//free stuff
		mode_free(output_mode);

		if(no_change_flag == 1){

			break;
		}
	}

	//Final forward run
	mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);

	//forward run
	*launched = complex_zero();
	*transmitted = complex_zero();
	*reflected = complex_zero();
	timefield3d_zero_fields(main_field);
	adjoint3d_reset_forward_fields(opt_area);
	for(i = 0; i < t_length; i++){

		if(i < 1000){
		
			amplitude = exp(-(i - 600)*(i - 600)/(2*200.0*200.0));
		}
		else{

			amplitude = 0;
		}
		phase = i*delta_t*omega;

		timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
		timefield3d_apply_timestep_omp_full(main_field, delta_t);

		instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
		instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 575, output_mode, phase, 0, transmitted);


		if(i%400 == 0){

			timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));
			printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
		}
	}

	reflection_coeff = complex_div(*reflected,*launched);
	transmission_coeff = complex_div(*transmitted,*launched);
	printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
	printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

	printf("\nTaper Length = %i nm\n", 20*200);
	printf("Segment number = %i\n", seg_num);
	printf("Taper widths:\n");
	for(i = 0; i < seg_num+1; i++){

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "brain_y_junction_design_5.bin", 20, 130, 150, 500, 50, 1.444*1.444, 3.476*3.476);
}

int main_y_from_starting_point(){
	
	int i, x_size, y_size, z_size, t_length, iteration;
	double disc, omega, delta_t, amplitude, phase;
	timefield3d main_field;
	double transmission_phase, improvement;
	complex instant_overlap, *launched, *transmitted, *reflected;
	complex transmission_coeff, reflection_coeff;
	mode launch_mode, output_mode;
	FILE *mode_data;
	adjoint3d opt_area;
	int seg_num, *taper_left, *taper_right, *widths;
	adjoint3d_taper opt_taper;
	int no_change_flag;
	
	//simulation setup
	x_size = 150;
	y_size = 100;
	z_size = 650;
	t_length = 4801;

	disc = 0.02;

	omega = 2*PI/1.550;
	delta_t = 0.99*disc/sqrt(3);

	launched = malloc(sizeof(complex));
	transmitted = malloc(sizeof(complex));
	reflected = malloc(sizeof(complex));

	main_field = timefield3d_new(x_size, y_size, z_size, disc);
		
	mode_data = fopen("mode_files/brain_y_junction_in_TE00.txt", "r");
	launch_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(launch_mode);

	mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
	output_mode = mode_file_load(mode_data);
	fclose(mode_data);
	mode_plot(output_mode);

	//adjoint setup
	opt_area = ajoint3d_new(20, 130, 46, 54, 150, 400, t_length);

	//build
	seg_num = 20;
	taper_left = malloc((seg_num+1)*sizeof(int));
	taper_right = malloc((seg_num+1)*sizeof(int));
	widths = malloc((seg_num+1)*sizeof(int));

	widths[0] = 24;
	widths[1] = 20;
	widths[2] = 24;
	widths[3] = 28;
	widths[4] = 34;
	widths[5] = 36;
	widths[6] = 36;
	widths[7] = 36;
	widths[8] = 44;
	widths[9] = 54;
	widths[10] = 52;
	widths[11] = 48;
	widths[12] = 44;
	widths[13] = 40;
	widths[14] = 58;
	widths[15] = 60;
	widths[16] = 66;
	widths[17] = 68;
	widths[18] = 68;
	widths[19] = 70;
	widths[20] = 74;

	timefield3d_set_background(main_field, 1.444*1.444, 1.0);
	timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
	timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
	timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 5, 95, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
	for(i = 0; i < seg_num+1; i++){

		taper_left[i] = 75 - widths[i]/2;
		taper_right[i] = 75 + widths[i]/2;
	}
	opt_taper = adjoint3d_taper_new(3.476*3.476, 1.0, seg_num, 150, 400, taper_left, taper_right, 46, 54);
	adjoint3d_taper_insert(main_field, opt_taper);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 87+25, 63+25, 87+25, 400, 500, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 500, 555, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 555, 645, 46, 54);

	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 97-25, 63-25, 87-25, 400, 500, 46, 54);
	timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 500, 555, 46, 54);
	timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 555, 645, 46, 54);

	timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
	timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

	no_change_flag = 0;
	for(iteration = 0; iteration < 3; iteration++){

		mode_data = fopen("mode_files/brain_y_junction_out_TE00.txt", "r");
		output_mode = mode_file_load(mode_data);
		fclose(mode_data);

		//forward run
		*launched = complex_zero();
		*transmitted = complex_zero();
		*reflected = complex_zero();
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_forward_fields(opt_area);
		for(i = 0; i < t_length; i++){

			if(i < 1000){
			
				amplitude = exp(-(i - 600)*(i - 600)/(2*200.0*200.0));
			}
			else{

				amplitude = 0;
			}
			phase = i*delta_t*omega;

			timefield3d_mode_source(main_field, 75, 50, 125, launch_mode, phase, amplitude, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 0, launched);
			//instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 130, launch_mode, phase, 1, reflected);
			instant_overlap = timefield3d_mode_overlap(main_field, 75, 50, 525, output_mode, phase, 0, transmitted);

			adjoint3d_store_fwd_fields(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_x_slice(main_field, -1, 1, (y_size/2));
				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\tlaunched field = %f\ttransmitted field = %f\treflected field = %f\n", i, complex_mag(*launched), complex_mag(*transmitted), complex_mag(*reflected));
			}
		}

		reflection_coeff = complex_div(*reflected,*launched);
		transmission_coeff = complex_div(*transmitted,*launched);
		printf("r = %f at %fπ\n", complex_mag(reflection_coeff), complex_angle(reflection_coeff)/PI);
		printf("t = %f at %fπ\n", complex_mag(transmission_coeff), complex_angle(transmission_coeff)/PI);

		transmission_phase = complex_angle(*transmitted);
		printf("\ntransmission phase used for adjoint = %f\n", transmission_phase);

		//adjoint simulation
		mode_smul(complex_exp(complex_from_doubles(0,transmission_phase+PI/2)), output_mode);
		timefield3d_zero_fields(main_field);
		adjoint3d_reset_gradient(opt_area);

		for(i = t_length-1; i >= 0; i--){

			phase = (i+1)*delta_t*omega;

			timefield3d_mode_source_adjoint(main_field, 75, 50, 525, output_mode, phase, 0);
			timefield3d_apply_timestep_omp_full(main_field, delta_t);

			adjoint3d_update_gradient(opt_area, main_field, i);

			if(i%400 == 0){

				//timefield3d_plot_Ex_y_slice(main_field, -1, 1, (y_size/2));

				printf("i = %i\n", i);
			}
		}

		//adjust permittivity
		/*improvement = 0.3*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_adjust_permittivity(opt_area, main_field, omega, improvement, 1.444*1.444, 3.476*3.476);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 25);*/

		//adjust taper
		improvement = 0.1*(complex_mag(*launched)-complex_mag(*transmitted));
		printf("improvement goal = %f\n", improvement);
		adjoint3d_taper_calculate_gradients(opt_taper, opt_area, 3.476*3.476-1.444*1.444);
		no_change_flag = adjoint3d_taper_update_symmetric(opt_taper, improvement, main_field.disc, 20);

		timefield3d_set_background(main_field, 1.444*1.444, 1.0);
		timefield3d_set_top_half_plane(main_field, 1.93*1.93, 1, 55); //nitride
		timefield3d_set_top_half_plane(main_field, 1.47*1.47, 1, 59); //BPSG
		timefield3d_set_top_half_plane(main_field, 1.39*1.39, 1, 84); //SiCOH
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 75, 75, 63, 87, 5, 95, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63, 87, 95, 150, 46, 54);
		adjoint3d_taper_insert(main_field, opt_taper);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 53+25, 87+25, 63+25, 87+25, 400, 500, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63+25, 87+25, 500, 555, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63+25, 87+25, 100, 100, 555, 645, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 97-25, 63-25, 87-25, 400, 500, 46, 54);
		timefield3d_set_waveguide(main_field, 3.476*3.476, 63-25, 87-25, 500, 555, 46, 54);
		timefield3d_set_taper(main_field, 3.476*3.476, 1.0, 63-25, 87-25, 50, 50, 555, 645, 46, 54);
		timefield3d_plot_permittivity_y_slice(main_field, 2, 11, 50);
		timefield3d_plot_permittivity_x_slice(main_field, 2, 11, 75);

		//free stuff
		mode_free(output_mode);

		if(no_change_flag == 1){

			break;
		}
	}

	printf("\nTaper Length = %i nm\n", 20*200);
	printf("Segment number = %i\n", seg_num);
	printf("Taper widths:\n");
	for(i = 0; i < seg_num+1; i++){

		printf("Width %i = %i nm\n", i, 20*(opt_taper.right[i]-opt_taper.left[i]));
	}

	timefield3d_file_write_permittivity_y_slice(main_field, "brain_y_junction_design_9.bin", 20, 130, 150, 500, 50, 1.444*1.444, 3.476*3.476);
}