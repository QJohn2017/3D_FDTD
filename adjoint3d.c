#include "adjoint3d.h"

#define PI 3.141592652

//general functions
adjoint3d ajoint3d_new(int x_start, int x_end, int y_start, int y_end, int z_start, int z_end, int t){
	
	adjoint3d new;
	int x_size, y_size, z_size;
	long long n, index, size;

	new.t = t;

	new.x_start = x_start;
	new.x_end = x_end;
	new.y_start = y_start;
	new.y_end = y_end;
	new.z_start = z_start;
	new.z_end = z_end;

	x_size = x_end-x_start;
	y_size = y_end-y_start;
	z_size = z_end-z_start;

	n = x_size*y_size*z_size;
	
	new.forward_fields = malloc(3*n*t*sizeof(double));
	new.gradient = malloc(n*sizeof(double));
	
	for(index = 0; index < 3*n*t; index++){
		
		new.forward_fields[index] = 0;
	}
	for(index = 0; index < n; index++){
		
		new.gradient[index] = 0;
	}

	return(new);
}

void adjointd3d_free(adjoint3d block){
	
	free(block.forward_fields);
	free(block.gradient);
}

void adjoint3d_store_fwd_fields(adjoint3d block, timefield3d main_fields, int timestep){

	int i, j, k, m;
	long long x_size, y_size, z_size;
	int field_index;
	long long block_index;

	x_size = block.x_end-block.x_start;
	y_size = block.y_end-block.y_start;
	z_size = block.z_end-block.z_start;

	for(i = block.x_start; i < block.x_end; i++){
		for(j = block.y_start; j < block.y_end; j++){
			for(k = block.z_start; k < block.z_end; k++){

				field_index = 6*(i + main_fields.i*(j + main_fields.j*k));
				block_index = 3*((i-block.x_start) + x_size*((j-block.y_start) + y_size*((k-block.z_start) + z_size*timestep)));

				for(m = 0; m < 3; m++){

					block.forward_fields[block_index + m] = main_fields.fields[field_index + m];
				}
			}
		}
	}
}

void adjoint3d_update_gradient(adjoint3d block, timefield3d main_fields, int timestep){

	int i, j, k, m;
	long long x_size, y_size, z_size;
	int field_index, gradient_index;
	long long block_index;
	double sum;

	x_size = block.x_end-block.x_start;
	y_size = block.y_end-block.y_start;
	z_size = block.z_end-block.z_start;

	for(i = block.x_start; i < block.x_end; i++){
		for(j = block.y_start; j < block.y_end; j++){
			for(k = block.z_start; k < block.z_end; k++){

				field_index = 6*(i + main_fields.i*(j + main_fields.j*k));
				block_index = 3*((i-block.x_start) + x_size*((j-block.y_start) + y_size*((k-block.z_start) + z_size*timestep)));
				gradient_index = (i-block.x_start) + x_size*((j-block.y_start) + y_size*(k-block.z_start));

				sum = 0;
				for(m = 0; m < 3; m++){

					sum += block.forward_fields[block_index+m]*main_fields.fields[field_index+m];
				}
				block.gradient[gradient_index] += sum;

			}
		}
	}
}

void adjoint3d_reset_gradient(adjoint3d block){

	long long n, index;

	n = (block.x_end-block.x_start)*(block.y_end-block.y_start)*(block.z_end-block.z_start);

	for(index = 0; index < n; index++){
		
		block.gradient[index] = 0;
	}
}

void adjoint3d_reset_forward_fields(adjoint3d block){

	long long n, index;

	n = 3*((long long)block.t*(block.x_end-block.x_start)*(block.y_end-block.y_start)*(block.z_end-block.z_start));

	for(index = 0; index < n; index++){

		block.forward_fields[index] = 0;
	}
}

void adjoint3d_adjust_permittivity(adjoint3d block, timefield3d main_fields, double omega, double improvement, double ep_min, double ep_max){

	int i, j, k, m, n;
	int x_size, y_size, z_size;
	int gradient_index, ep_index;
	double grad_mag, *grad_dir;

	x_size = block.x_end-block.x_start;
	y_size = block.y_end-block.y_start;
	z_size = block.z_end-block.z_start;

	n = x_size*y_size*z_size;

	grad_mag = 0;
	for(m = 0; m < n; m++){

		grad_mag += block.gradient[m]*block.gradient[m];
	}
	grad_mag = sqrt(grad_mag)*main_fields.disc*main_fields.disc;

	grad_dir = malloc(n*sizeof(double));
	for(m = 0; m < n; m++){

		grad_dir[m] = block.gradient[m]*main_fields.disc*main_fields.disc/grad_mag;
	}
	grad_mag *= omega;

	for(i = block.x_start; i < block.x_end; i++){
		for(j = block.y_start; j < block.y_end; j++){
			for(k = block.z_start; k < block.z_end; k++){

				gradient_index = (i-block.x_start) + x_size*((j-block.y_start) + y_size*(k-block.z_start));
				ep_index = 3*(i + main_fields.i*(j + main_fields.j*k));

				for(m = 0; m < 3; m++){

					main_fields.ep[ep_index+m] += (improvement/grad_mag)*grad_dir[gradient_index];

					//ensure realistic permittivity values
					if(main_fields.ep[ep_index+m] < ep_min){

						main_fields.ep[ep_index+m] = ep_min;
					}

					if(main_fields.ep[ep_index+m] > ep_max){

						main_fields.ep[ep_index+m] = ep_max;
					}
				}
			}
		}
	}

	free(grad_dir);
}

adjoint3d_taper adjoint3d_taper_new(double ep, double mu, int seg_num, int start, int end, int *left, int *right, int bottom, int top, int duration){

	int i, n;
	int height, length;
	adjoint3d_taper new;

	new.ep = ep;
	new.mu = mu;
	new.seg_num = seg_num;
	new.start = start;
	new.end = end;
	new.bottom = bottom;
	new.top = top;
	new.duration = duration;

	//adjustable edges
	new.left = malloc((seg_num+1)*sizeof(int));
	new.right = malloc((seg_num+1)*sizeof(int));
	for(i = 0; i < seg_num+1; i++){

		new.left[i] = left[i];
		new.right[i] = right[i];
	}

	//fwd field storage along edges
	height = top-bottom;
	length = end-start;
	n = 6*2*height*length*duration;
	new.forward_fields = malloc(n*sizeof(double));
	new.adjoint_fields = malloc(n*sizeof(double));
	for(i = 0; i < n; i++){

		new.forward_fields[i] = 0;
		new.adjoint_fields[i] = 0;
	}

	//gradient w.r.t epsilon along edges
	n = 2*height*length;
	new.ep_gradient = malloc(n*sizeof(double));
	for(i = 0; i < n; i++){

		new.ep_gradient[i] = 0;
	}

	//edge gradients
	new.left_grad = malloc((seg_num+1)*sizeof(double));
	new.right_grad = malloc((seg_num+1)*sizeof(double));
	for(i = 0; i < seg_num+1; i++){

		new.left_grad[i] = 0;
		new.right_grad[i] = 0;
	}

	return(new);
}

void adjoint3d_taper_insert(timefield3d psi, adjoint3d_taper taper){

	int n;
	int seg_start, seg_end, taper_length;

	taper_length = taper.end - taper.start;

	for(n = 0; n < taper.seg_num; n++){

		seg_start = taper.start + (n*taper_length)/taper.seg_num;
		seg_end = taper.start + ((n+1)*taper_length)/taper.seg_num;

		timefield3d_set_taper(psi, taper.ep, taper.mu, taper.left[n], taper.right[n], taper.left[n+1], taper.right[n+1], seg_start, seg_end, taper.bottom, taper.top);
	}
}

void adjoint3d_taper_store_forward_fields(adjoint3d_taper taper, timefield3d main_fields, int timestep){

	int i, j, k, n, m;
	int instant_field_index, taper_field_index;
	int taper_length, taper_height;
	int seg_mid, seg_fwd;

	taper_length = taper.end-taper.start;
	taper_height = taper.top-taper.bottom;
	for(n = 0; n < taper.seg_num; n++){

		seg_mid = taper.start + ((n)*taper_length)/taper.seg_num;
		seg_fwd = taper.start + ((n+1)*taper_length)/taper.seg_num;
		for(k = seg_mid; k < seg_fwd; k++){
			for(j = taper.bottom; j < taper.top; j++){

				//right
				i = taper.right[n] + round(((taper.right[n+1]-taper.right[n])*(k-seg_mid))/(seg_fwd-seg_mid));
				instant_field_index = 6*(i + main_fields.i*(j + main_fields.j*k));
				taper_field_index = 6*(0 + 2*((j-taper.bottom) + taper_height*((k-taper.start)+ taper_length*timestep)));
				for(m = 0; m < 6; m++){

					taper.forward_fields[taper_field_index+m] = main_fields.fields[instant_field_index+m];
				}

				//left
				i = taper.left[n] + round(((taper.left[n+1]-taper.left[n])*(k-seg_mid))/(seg_fwd-seg_mid));
				instant_field_index = 6*(i + main_fields.i*(j + main_fields.j*k));
				taper_field_index = 6*(1 + 2*((j-taper.bottom) + taper_height*((k-taper.start)+ taper_length*timestep)));
				for(m = 0; m < 6; m++){

					taper.forward_fields[taper_field_index+m] = main_fields.fields[instant_field_index+m];
				}

			}
		}

	}
}

void adjoint3d_taper_store_adjoint_fields(adjoint3d_taper taper, timefield3d main_fields, int timestep){

	int i, j, k, n, m;
	int instant_field_index, taper_field_index;
	int taper_length, taper_height;
	int seg_mid, seg_fwd;

	taper_length = taper.end-taper.start;
	taper_height = taper.top-taper.bottom;
	for(n = 0; n < taper.seg_num; n++){

		seg_mid = taper.start + ((n)*taper_length)/taper.seg_num;
		seg_fwd = taper.start + ((n+1)*taper_length)/taper.seg_num;
		for(k = seg_mid; k < seg_fwd; k++){
			for(j = taper.bottom; j < taper.top; j++){

				//right
				i = taper.right[n] + round(((taper.right[n+1]-taper.right[n])*(k-seg_mid))/(seg_fwd-seg_mid));
				instant_field_index = 6*(i + main_fields.i*(j + main_fields.j*k));
				taper_field_index = 6*(0 + 2*((j-taper.bottom) + taper_height*((k-taper.start)+ taper_length*timestep)));
				for(m = 0; m < 6; m++){

					taper.adjoint_fields[taper_field_index+m] = main_fields.fields[instant_field_index+m];
				}

				//left
				i = taper.left[n] + round(((taper.left[n+1]-taper.left[n])*(k-seg_mid))/(seg_fwd-seg_mid));
				instant_field_index = 6*(i + main_fields.i*(j + main_fields.j*k));
				taper_field_index = 6*(1 + 2*((j-taper.bottom) + taper_height*((k-taper.start)+ taper_length*timestep)));
				for(m = 0; m < 6; m++){

					taper.adjoint_fields[taper_field_index+m] = main_fields.fields[instant_field_index+m];
				}

			}
		}

	}
}

void adjoint3d_taper_calculate_ep_gradient(adjoint3d_taper taper){

	int i, j, k, n, m, t;
	int field_index, ep_gradient_index;
	int taper_length, taper_height;
	double sum;

	taper_length = taper.end-taper.start;
	taper_height = taper.top-taper.bottom;

	//Update ep-gradient
	taper_length = taper.end-taper.start;
	taper_height = taper.top-taper.bottom;
	for(i = 0; i < 2; i++){
		for(j = 0 ; j < taper_height; j++){
			for(k = 0; k < taper_length; k++){

				ep_gradient_index = i + 2*(j + taper_height*k);

				sum = 0;
				for(t = 0; t < taper.duration; t++){

					field_index = 6*(i + 2*(j + taper_height*(k + taper_length*t)));
					
					for(m = 0; m < 3; m++){

						sum += taper.forward_fields[field_index+m]*taper.adjoint_fields[field_index+m];
					}
				}

				taper.ep_gradient[ep_gradient_index] = sum;
			}
		}
	}
}

void adjoint3d_taper_calculate_gradients_from_ep_gradient(adjoint3d_taper taper, double ep_diff){

	int i, j, k, n, grad_index;
	int taper_length, taper_height;
	int seg_bck, seg_mid, seg_fwd;

	taper_length = taper.end-taper.start;
	taper_height = taper.top-taper.bottom;

	taper.left_grad[0] = 0;
	taper.right_grad[0] = 0;
	for(n = 1; n < taper.seg_num; n++){

		taper.left_grad[n] = 0;
		taper.right_grad[n] = 0;

		seg_bck = ((n-1)*taper_length)/taper.seg_num;
		seg_mid = ((n)*taper_length)/taper.seg_num;
		seg_fwd = ((n+1)*taper_length)/taper.seg_num;

		for(k = seg_bck; k < seg_mid; k++){
			for(j = 0; j < taper_height; j++){

				grad_index = 0 + 2*(j + taper_height*k);
				taper.right_grad[n] += ((double)(k-seg_bck)/(double)(seg_mid-seg_bck))*taper.ep_gradient[grad_index]*ep_diff;

				grad_index = 1 + 2*(j + taper_height*k);
				taper.left_grad[n] -= ((double)(k-seg_bck)/(double)(seg_mid-seg_bck))*taper.ep_gradient[grad_index]*ep_diff;
			}
		}
		for(k = seg_mid; k < seg_fwd; k++){
			for(j = 0; j < taper_height; j++){

				grad_index = 0 + 2*(j + taper_height*k);
				taper.right_grad[n] += ((double)(seg_fwd-k)/(double)(seg_fwd-seg_mid))*taper.ep_gradient[grad_index]*ep_diff;

				grad_index = 1 + 2*(j + taper_height*k);
				taper.left_grad[n] -= ((double)(seg_fwd-k)/(double)(seg_fwd-seg_mid))*taper.ep_gradient[grad_index]*ep_diff;
			}
		}
	}
	taper.left_grad[taper.seg_num] = 0;
	taper.right_grad[taper.seg_num] = 0;
}

void adjoint3d_taper_calculate_gradients(adjoint3d_taper taper, adjoint3d block, double ep_diff){

	int i, j, k, n, grad_index;
	int x_size, y_size, z_size, taper_length;
	int seg_bck, seg_mid, seg_fwd;

	x_size = block.x_end-block.x_start;
	y_size = block.y_end-block.y_start;
	z_size = block.z_end-block.z_start;
	taper_length = taper.end-taper.start;

	taper.left_grad[0] = 0;
	taper.right_grad[0] = 0;
	for(n = 1; n < taper.seg_num; n++){

		taper.left_grad[n] = 0;
		taper.right_grad[n] = 0;

		seg_bck = taper.start + ((n-1)*taper_length)/taper.seg_num;
		seg_mid = taper.start + ((n)*taper_length)/taper.seg_num;
		seg_fwd = taper.start + ((n+1)*taper_length)/taper.seg_num;

		for(k = seg_bck; k < seg_mid; k++){
			for(j = taper.bottom; j < taper.top; j++){

				i = taper.right[n-1] + ((taper.right[n]-taper.right[n-1])*(k-seg_bck))/(seg_mid-seg_bck);
				grad_index = (i-block.x_start) + x_size*((j-block.y_start) + y_size*(k-block.z_start));
				taper.right_grad[n] += ((double)(k-seg_bck)/(double)(seg_mid-seg_bck))*block.gradient[grad_index]*ep_diff;

				i = taper.left[n-1] + ((taper.left[n]-taper.left[n-1])*(k-seg_bck))/(seg_mid-seg_bck);
				grad_index = (i-block.x_start) + x_size*((j-block.y_start) + y_size*(k-block.z_start));
				taper.left_grad[n] -= ((double)(k-seg_bck)/(double)(seg_mid-seg_bck))*block.gradient[grad_index]*ep_diff;
			}
		}
		for(k = seg_mid; k < seg_fwd; k++){
			for(j = taper.bottom; j < taper.top; j++){

				i = taper.right[n] + ((taper.right[n+1]-taper.right[n])*(k-seg_mid))/(seg_fwd-seg_mid);
				grad_index = (i-block.x_start) + x_size*((j-block.y_start) + y_size*(k-block.z_start));
				taper.right_grad[n] += ((double)(seg_fwd-k)/(double)(seg_fwd-seg_mid))*block.gradient[grad_index]*ep_diff;

				i = taper.left[n] + ((taper.left[n+1]-taper.left[n])*(k-seg_mid))/(seg_fwd-seg_mid);
				grad_index = (i-block.x_start) + x_size*((j-block.y_start) + y_size*(k-block.z_start));
				taper.left_grad[n] -= ((double)(seg_fwd-k)/(double)(seg_fwd-seg_mid))*block.gradient[grad_index]*ep_diff;
			}
		}
	}
	taper.left_grad[taper.seg_num] = 0;
	taper.right_grad[taper.seg_num] = 0;
}

int adjoint3d_taper_update_symmetric(adjoint3d_taper taper, double improvement_goal, double disc, int min_width){

	int i;
	int min_index, no_change_flag, zero_counter;
	double grad_mag, *grad_dir;
	double min_component;

	//get magnitude and direction of gradient
	grad_dir = malloc((taper.seg_num+1)*sizeof(double));
	grad_mag = 0;
	for(i = 0; i < taper.seg_num+1; i++){

		grad_mag += (taper.right_grad[i]-taper.left_grad[i])*(taper.right_grad[i]-taper.left_grad[i]);
	}
	grad_mag = sqrt(grad_mag)*0.5*disc*disc;

	for(i = 0; i < taper.seg_num+1; i++){

		grad_dir[i] = (taper.right_grad[i]-taper.left_grad[i])*0.5*disc*disc/grad_mag;
	}

	//if all improvements are smaller than 1 pixel, zero out smallest gradient components and recalculate
	no_change_flag = 1;
	for(i = 1; i < taper.seg_num; i++){

		if(fabs((improvement_goal/grad_mag)*grad_dir[i]) >= 0.5){

			no_change_flag = 0;
		}
	}

	zero_counter = 0;
	while((no_change_flag == 1)&&(zero_counter < taper.seg_num-1)){

		min_component = 1000;
		min_index = 0;
		for(i = 1; i < taper.seg_num; i++){

			grad_mag = fabs(grad_dir[i])*2*disc;
			if((grad_mag < min_component)&&(grad_mag > 0)){

				min_component = grad_mag;
				min_index = i;
			}
		}
	
		taper.right_grad[min_index] = 0;
		taper.left_grad[min_index] = 0;

		grad_mag = 0;
		for(i = 0; i < taper.seg_num+1; i++){

			grad_mag += (taper.right_grad[i]-taper.left_grad[i])*(taper.right_grad[i]-taper.left_grad[i]);
		}
		grad_mag = sqrt(grad_mag)*0.5*disc*disc;

		for(i = 0; i < taper.seg_num+1; i++){

			grad_dir[i] = (taper.right_grad[i]-taper.left_grad[i])*0.5*disc*disc/grad_mag;
		}

		//if all improvements are smaller than 1 pixel, zero out smallest gradient components and recalculate
		no_change_flag = 1;
		for(i = 1; i < taper.seg_num; i++){

			if(fabs((improvement_goal/grad_mag)*grad_dir[i]) >= 0.5){

				no_change_flag = 0;
			}
		}

		printf("zeroed out index %i\n", min_index);
		zero_counter++;
	}

	//update
	for(i = 1; i < taper.seg_num; i++){

		taper.left[i] -= (int)(round((improvement_goal/grad_mag)*grad_dir[i]));
		taper.right[i] += (int)(round((improvement_goal/grad_mag)*grad_dir[i]));

		if(taper.right[i]-taper.left[i] < min_width){

			taper.right[i] = (taper.right[i]+taper.left[i]+min_width)/2;
			taper.left[i] = (taper.right[i]+taper.left[i]-min_width)/2;
		}
	}

	free(grad_dir);

	return(no_change_flag);
}

int adjoint3d_taper_update_asymmetric(adjoint3d_taper taper, double improvement_goal, double disc, int min_width){

	int i;
	int min_index, no_change_flag, zero_counter;
	double grad_mag, *grad_dir;
	double min_component;

	//get magnitude and direction of gradient
	grad_dir = malloc(2*(taper.seg_num+1)*sizeof(double));
	grad_mag = 0;
	for(i = 0; i < taper.seg_num+1; i++){

		grad_mag += taper.right_grad[i]*taper.right_grad[i];
		grad_mag += taper.left_grad[i]*taper.left_grad[i];
	}
	grad_mag = sqrt(grad_mag)*0.5*disc*disc;

	for(i = 0; i < taper.seg_num+1; i++){

		grad_dir[2*i] = taper.right_grad[i]*0.5*disc*disc/grad_mag;
		grad_dir[2*i+1] = taper.left_grad[i]*0.5*disc*disc/grad_mag;
	}

	//if all improvements are smaller than 1 pixel, zero out smallest gradient components and recalculate
	no_change_flag = 1;
	for(i = 2; i < 2*taper.seg_num; i++){

		if(round(fabs((improvement_goal/grad_mag)*grad_dir[i])) >= 1){

			no_change_flag = 0;
		}
	}

	zero_counter = 0;
	while((no_change_flag == 1)&&(zero_counter < 2*(taper.seg_num-2))){

		min_component = 1e10;
		min_index = 0;
		for(i = 2; i < 2*taper.seg_num; i++){

			grad_mag = fabs(grad_dir[i]);
			if((grad_mag < min_component)&&(grad_mag > 0)){

				min_component = grad_mag;
				min_index = i;
			}
		}

		if(min_index%2 == 0){

			taper.right_grad[min_index/2] = 0;
		}
		else{

			taper.left_grad[min_index/2] = 0;
		}

		grad_mag = 0;
		for(i = 0; i < taper.seg_num+1; i++){

			grad_mag += taper.right_grad[i]*taper.right_grad[i];
			grad_mag += taper.left_grad[i]*taper.left_grad[i];
		}
		grad_mag = sqrt(grad_mag)*0.5*disc*disc;

		for(i = 0; i < taper.seg_num+1; i++){

			grad_dir[2*i] = taper.right_grad[i]*0.5*disc*disc/grad_mag;
			grad_dir[2*i+1] = taper.left_grad[i]*0.5*disc*disc/grad_mag;
		}

		//if all improvements are smaller than 1 pixel, zero out smallest gradient components and recalculate
		no_change_flag = 1;
		for(i = 2; i < 2*taper.seg_num; i++){

			if(round(fabs((improvement_goal/grad_mag)*grad_dir[i])) >= 1){

				no_change_flag = 0;
			}
		}

		printf("zeroed out index %i\n", min_index);
		zero_counter++;
	}

	//update
	for(i = 1; i < taper.seg_num; i++){

		taper.right[i] += (int)(round((improvement_goal/grad_mag)*grad_dir[2*i]));
		taper.left[i] += (int)(round((improvement_goal/grad_mag)*grad_dir[2*i+1]));

		if(taper.right[i]-taper.left[i] < min_width){

			taper.right[i] = (taper.right[i]+taper.left[i]+min_width)/2;
			taper.left[i] = (taper.right[i]+taper.left[i]-min_width)/2;
		}
	}

	free(grad_dir);

	return(no_change_flag);
}