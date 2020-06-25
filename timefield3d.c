#include "timefield3d.h"

#define PI 3.141592652

//general functions
timefield3d timefield3d_new(int i, int j, int k, double disc){
	
	timefield3d new;
	int n, index;

	new.i = i;
	new.j = j;
	new.k = k;
	new.disc = disc;

	n = i*j*k;
	
	new.ep = malloc(3*n*sizeof(double));
	new.mu = malloc(3*n*sizeof(double));
	new.sig = malloc(6*n*sizeof(double));
	new.fields = malloc(6*n*sizeof(double));
	new.old_fields = malloc(6*n*sizeof(double));
	new.fourier_fields = malloc(6*n*sizeof(complex));
	new.sources = malloc(6*n*sizeof(double));
	
	for(index = 0; index < 3*n; index++){
		
		new.ep[index] = 1;
		new.mu[index] = 1;
	}
	for(index = 0; index < 6*n; index++){
		
		new.sig[index] = 0;
		new.fields[index] = 0;
		new.old_fields[index] = 0;
		new.fourier_fields[index] = complex_zero();
		new.sources[index] = 0;
	}

	return(new);
}

void timefield3d_free(timefield3d psi){
	
	free(psi.ep);
	free(psi.mu);
	free(psi.sig);
	free(psi.fields);
	free(psi.old_fields);
	free(psi.fourier_fields);
	free(psi.sources);
}

timefield3d timefield3d_copy(timefield3d psi){
	
	int i, field_index, ep_index;
	timefield3d copy;
	
	copy = timefield3d_new(psi.i, psi.j, psi.k, psi.disc);
	
	//copy contents into new matrices
	for(i = 0; i < psi.i*psi.j*psi.k; i++){
		
		field_index = 6*i;
		ep_index = 3*i;

		copy.ep[ep_index] = psi.ep[ep_index];
		copy.ep[ep_index + 1] = psi.ep[ep_index + 1];
		copy.ep[ep_index + 2] = psi.ep[ep_index + 2];
		copy.mu[ep_index] = psi.mu[ep_index];
		copy.mu[ep_index + 1] = psi.mu[ep_index + 1];
		copy.mu[ep_index + 2] = psi.mu[ep_index + 2];
		copy.sig[field_index] = psi.fields[field_index];
		copy.sig[field_index + 1] = psi.sig[field_index + 1];
		copy.sig[field_index + 2] = psi.sig[field_index + 2];
		copy.sig[field_index + 3] = psi.sig[field_index + 3];
		copy.sig[field_index + 4] = psi.sig[field_index + 4];
		copy.sig[field_index + 5] = psi.sig[field_index + 5];
		copy.fields[field_index] = psi.fields[field_index];
		copy.fields[field_index + 1] = psi.fields[field_index + 1];
		copy.fields[field_index + 2] = psi.fields[field_index + 2];
		copy.fields[field_index + 3] = psi.fields[field_index + 3];
		copy.fields[field_index + 4] = psi.fields[field_index + 4];
		copy.fields[field_index + 5] = psi.fields[field_index + 5];
		copy.old_fields[field_index] = psi.old_fields[field_index];
		copy.old_fields[field_index + 1] = psi.old_fields[field_index + 1];
		copy.old_fields[field_index + 2] = psi.old_fields[field_index + 2];
		copy.old_fields[field_index + 3] = psi.old_fields[field_index + 3];
		copy.old_fields[field_index + 4] = psi.old_fields[field_index + 4];
		copy.old_fields[field_index + 5] = psi.old_fields[field_index + 5];
		copy.fourier_fields[field_index] = psi.fourier_fields[field_index];
		copy.fourier_fields[field_index + 1] = psi.fourier_fields[field_index + 1];
		copy.fourier_fields[field_index + 2] = psi.fourier_fields[field_index + 2];
		copy.fourier_fields[field_index + 3] = psi.fourier_fields[field_index + 3];
		copy.fourier_fields[field_index + 4] = psi.fourier_fields[field_index + 4];
		copy.fourier_fields[field_index + 5] = psi.fourier_fields[field_index + 5];
		copy.sources[field_index] = psi.sources[field_index];
		copy.sources[field_index + 1] = psi.sources[field_index + 1];
		copy.sources[field_index + 2] = psi.sources[field_index + 2];
		copy.sources[field_index + 3] = psi.sources[field_index + 3];
		copy.sources[field_index + 4] = psi.sources[field_index + 4];
		copy.sources[field_index + 5] = psi.sources[field_index + 5];
	}
	
	return(copy);
}

void timefield3d_zero_fields(timefield3d psi){

	int n, index;

	n = psi.i*psi.j*psi.k;

	for(index = 0; index < 6*n; index++){

		psi.fields[index] = 0.0;
		psi.old_fields[index] = 0.0;
		psi.fourier_fields[index] = complex_zero();
		psi.sources[index] = 0.0;
	}
}

double* timefield3d_get_fields(timefield3d psi){
	
	double *fields;
	int i, n;
	
	n = 6*psi.i*psi.j*psi.k;
	fields = malloc(n*sizeof(double));
	for(i = 0; i < n; i++){
		
		fields[i] = psi.fields[i];
	}
	
	return(fields);
}

double* timefield3d_get_edge_fields(timefield3d psi){
	
	double *fields;
	int i, j, k, m, n;
	int field_index, field_index_inner;
	
	n = 6*psi.i*psi.j*psi.k;
	fields = malloc(n*sizeof(double));
	
	//x = 0 boundary
	i = 0;
	for(j = 0; j < psi.j; j++){
		for(k = 0; k < psi.k; k++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i+1) + psi.i*(j + psi.j*k));
			
			for(m = 0; m < 6; m++){

				fields[field_index + m] = psi.fields[field_index + m];
				fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//y = 0 boundary
	j = 0;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j+1) + psi.j*k));
			
			for(m = 0; m < 6; m++){

				fields[field_index + m] = psi.fields[field_index + m];
				fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//z = 0 boundary
	k = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*(j + psi.j*(k+1)));
			
			for(m = 0; m < 6; m++){

				fields[field_index + m] = psi.fields[field_index + m];
				fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//x = end boundary
	i = psi.i - 1;
	for(j = 0; j < psi.j; j++){
		for(k = 0; k < psi.k; k++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i-1) + psi.i*(j + psi.j*k));
			
			for(m = 0; m < 6; m++){

				fields[field_index + m] = psi.fields[field_index + m];
				fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//y = end boundary
	j = psi.j - 1;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j-1) + psi.j*k));
			
			for(m = 0; m < 6; m++){

				fields[field_index + m] = psi.fields[field_index + m];
				fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//z = end boundary
	k = psi.k - 1;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*(j + psi.j*(k-1)));
			
			for(m = 0; m < 6; m++){

				fields[field_index + m] = psi.fields[field_index + m];
				fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	return(fields);
}

void timefield3d_update_old_edge_fields(timefield3d psi){

	int i, j, k, m, n;
	int field_index, field_index_inner;
		
	//x = 0 boundary
	i = 0;
	for(j = 0; j < psi.j; j++){
		for(k = 0; k < psi.k; k++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i+1) + psi.i*(j + psi.j*k));
			
			for(m = 0; m < 6; m++){

				psi.old_fields[field_index + m] = psi.fields[field_index + m];
				psi.old_fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//y = 0 boundary
	j = 0;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j+1) + psi.j*k));
			
			for(m = 0; m < 6; m++){

				psi.old_fields[field_index + m] = psi.fields[field_index + m];
				psi.old_fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//z = 0 boundary
	k = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*(j + psi.j*(k+1)));
			
			for(m = 0; m < 6; m++){

				psi.old_fields[field_index + m] = psi.fields[field_index + m];
				psi.old_fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//x = end boundary
	i = psi.i - 1;
	for(j = 0; j < psi.j; j++){
		for(k = 0; k < psi.k; k++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i-1) + psi.i*(j + psi.j*k));
			
			for(m = 0; m < 6; m++){

				psi.old_fields[field_index + m] = psi.fields[field_index + m];
				psi.old_fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//y = end boundary
	j = psi.j - 1;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j-1) + psi.j*k));
			
			for(m = 0; m < 6; m++){

				psi.old_fields[field_index + m] = psi.fields[field_index + m];
				psi.old_fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
	
	//z = end boundary
	k = psi.k - 1;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*(j + psi.j*(k-1)));
			
			for(m = 0; m < 6; m++){

				psi.old_fields[field_index + m] = psi.fields[field_index + m];
				psi.old_fields[field_index_inner + m] = psi.fields[field_index_inner + m];
			}
		}
	}
}

void timefield3d_plot_Ex_x_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.k; i++){
			
			index = slice_index + psi.i*(j + psi.j*i);
			fprintf(Ex, "%lf\t", psi.fields[6*index]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Ex_y_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *X, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and X axis values
	X = malloc(psi.i*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.i; j++){
		
		fprintf(Ex, "%lf\t", X[j]);
		for(i = 0; i < psi.k; i++){
			
			index = j + psi.i*(slice_index + psi.j*i);
			fprintf(Ex, "%lf\t", psi.fields[6*index]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.i+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Ex_z_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *X, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	X = malloc(psi.i*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.i; j++){
		
		X[j] = (j + 0.5 - psi.i/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + psi.i*(j + psi.j*slice_index);
			fprintf(Ex, "%lf\t", psi.fields[6*index]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.i+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Ey_x_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.k; i++){
			
			index = slice_index + psi.i*(j + psi.j*i);
			fprintf(Ex, "%lf\t", psi.fields[6*index+1]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Ey_y_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *X, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and X axis values
	X = malloc(psi.i*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.i; j++){
		
		fprintf(Ex, "%lf\t", X[j]);
		for(i = 0; i < psi.k; i++){
			
			index = j + psi.i*(slice_index + psi.j*i);
			fprintf(Ex, "%lf\t", psi.fields[6*index+1]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.i+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Ey_z_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *X, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	X = malloc(psi.i*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.i; j++){
		
		X[j] = (j + 0.5 - psi.i/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + psi.i*(j + psi.j*slice_index);
			fprintf(Ex, "%lf\t", psi.fields[6*index+1]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.i+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Hx_x_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.k; i++){
			
			index = slice_index + psi.i*(j + psi.j*i);
			fprintf(Ex, "%lf\t", psi.fields[6*index+3]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Hx_y_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *X, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and X axis values
	X = malloc(psi.i*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.i; j++){
		
		fprintf(Ex, "%lf\t", X[j]);
		for(i = 0; i < psi.k; i++){
			
			index = j + psi.i*(slice_index + psi.j*i);
			fprintf(Ex, "%lf\t", psi.fields[6*index+3]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.i+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Hx_z_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *X, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	X = malloc(psi.i*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.i; j++){
		
		X[j] = (j + 0.5 - psi.i/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + psi.i*(j + psi.j*slice_index);
			fprintf(Ex, "%lf\t", psi.fields[6*index+3]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.i+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Hy_x_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.k; i++){
			
			index = slice_index + psi.i*(j + psi.j*i);
			fprintf(Ex, "%lf\t", psi.fields[6*index+4]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Hy_y_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *X, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and X axis values
	X = malloc(psi.i*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.i; j++){
		
		fprintf(Ex, "%lf\t", X[j]);
		for(i = 0; i < psi.k; i++){
			
			index = j + psi.i*(slice_index + psi.j*i);
			fprintf(Ex, "%lf\t", psi.fields[6*index+4]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.i+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Hy_z_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *X, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	X = malloc(psi.i*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.i; j++){
		
		X[j] = (j + 0.5 - psi.i/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + psi.i*(j + psi.j*slice_index);
			fprintf(Ex, "%lf\t", psi.fields[6*index+4]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.i+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Ex_x_slice_fourier(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.k; i++){
			
			index = slice_index + psi.i*(j + psi.j*i);
			fprintf(Ex, "%lf\t", complex_mag(psi.fourier_fields[6*index]));
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Ex_y_slice_fourier(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *X, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and X axis values
	X = malloc(psi.i*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.i; j++){
		
		fprintf(Ex, "%lf\t", X[j]);
		for(i = 0; i < psi.k; i++){
			
			index = j + psi.i*(slice_index + psi.j*i);
			fprintf(Ex, "%lf\t", complex_mag(psi.fourier_fields[6*index]));
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.i+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_Ex_z_slice_fourier(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *X, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	X = malloc(psi.i*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.i; j++){
		
		X[j] = (j + 0.5 - psi.i/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + psi.i*(j + psi.j*slice_index);
			fprintf(Ex, "%lf\t", complex_mag(psi.fourier_fields[6*index]));
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.i+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_permittivity_x_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.k; i++){
			
			index = slice_index + psi.i*(j + psi.j*i);
			fprintf(Ex, "%lf\t", psi.ep[3*index]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_permittivity_y_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *X, *Z, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and X axis values
	X = malloc(psi.i*sizeof(double));
	Z = malloc(psi.k*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.k; j++){
		
		Z[j] = (j + 0.5 - psi.k/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.k; i++){
		
		fprintf(Ex, "%lf\t", Z[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.i; j++){
		
		fprintf(Ex, "%lf\t", X[j]);
		for(i = 0; i < psi.k; i++){
			
			index = j + psi.i*(slice_index + psi.j*i);
			fprintf(Ex, "%lf\t", psi.ep[3*index]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.k+100, 2*psi.i+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_plot_permittivity_z_slice(timefield3d psi, int slice_index){
	
	int i, j, index;
	double *Y, *X, y_height;
	FILE *gnuplot;
	FILE *Ex;
	
	//y_height = (2*0.33*psi.j)/psi.k;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	Y = malloc(psi.j*sizeof(double));
	X = malloc(psi.i*sizeof(double));
	
	for(i = 0; i < psi.j; i++){
		
		Y[i] = (i + 0.5 - psi.j/2.0)*psi.disc;
	}
	for(j = 0; j < psi.i; j++){
		
		X[j] = (j + 0.5 - psi.i/2.0)*psi.disc;
	}
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.k);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + psi.i*(j + psi.j*slice_index);
			fprintf(Ex, "%lf\t", psi.ep[3*index]);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);

	
	//plot the index and field distributions
	fprintf(gnuplot, "unset colorbox\n");
	fprintf(gnuplot, "set term x11 size %i,%i background \"gray\";\n", 2*psi.i+100, 2*psi.j+100);
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	//fprintf(gnuplot, "set cbrange[%f:%f];\n", color_min, color_max);
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void timefield3d_set_background(timefield3d psi, float ep, float mu){
	
	int i, n;
	
	n = 3*psi.i*psi.j*psi.k;
	for(i = 0; i < n; i++){
		
		psi.ep[i] = ep;
		psi.mu[i] = mu;
	}
}

void timefield3d_set_top_half_plane(timefield3d psi, float ep, float mu, int bottom){

	int i, j, k, index;

	for(i = 0; i < psi.i; i++){
		for(j = bottom; j < psi.j; j++){
			for(k = 0; k < psi.k; k++){

				index = 3*(i + psi.i*(j + psi.j*k));
				if(j >= bottom){

					psi.ep[index] = ep;
					psi.mu[index] = mu;
					psi.ep[index + 1] = ep;
					psi.mu[index + 1] = mu;
					psi.ep[index + 2] = ep;
					psi.mu[index + 2] = mu;

				}
			}
		}
	}
}

void timefield3d_set_background_permittivity(timefield3d psi, float ep){
	
	int i, n;
	
	n = 3*psi.i*psi.j*psi.k;
	for(i = 0; i < n; i++){
		
		psi.ep[i] = ep;
	}
}

void timefield3d_set_background_permeability(timefield3d psi, float mu){
	
	int i, n;
	
	n = 3*psi.i*psi.j*psi.k;
	for(i = 0; i < n; i++){
		
		psi.mu[i] = mu;
	}
}

void timefield3d_set_background_electric_conductivity(timefield3d psi, float sig){
	
	int i, n;
	int sig_index;
	
	n = psi.i*psi.j*psi.k;
	for(i = 0; i < n; i++){
		
		sig_index = 6*i;
		
		psi.sig[sig_index] = sig;
		psi.sig[sig_index + 1] = sig;
		psi.sig[sig_index + 2] = sig;
	}
}

void timefield3d_set_background_magnetic_conductivity(timefield3d psi, float sig){
	
	int i, n;
	int sig_index;
	
	n = psi.i*psi.j*psi.k;
	for(i = 0; i < n; i++){
		
		sig_index = 6*i;
		
		psi.sig[sig_index + 3] = sig;
		psi.sig[sig_index + 4] = sig;
		psi.sig[sig_index + 5] = sig;
	}
}

void timefield3d_set_background_complex_permittivity(timefield3d psi, double omega, float real_ep, float imag_ep){
	
	//set real part of epsilon
	timefield3d_set_background_permittivity(psi, real_ep);
	
	//set imaginary part as electric conductivity
	timefield3d_set_background_electric_conductivity(psi, omega*imag_ep);
}

void timefield3d_set_prism_permittivity(timefield3d psi, float ep, int start_x, int end_x, int start_y, int end_y, int start_z, int end_z){
	
	int i, j, k, ep_index;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			for(k = 0; k < psi.k; k++){
				
				if((i >= start_x)&&(i < end_x)&&(j >= start_y)&&(j < end_y)&&(k >= start_z)&&(k < end_z)){
					
					ep_index = 3*(i + psi.i*(j + psi.j*k));
					
					psi.ep[ep_index] = ep;
					psi.ep[ep_index + 1] = ep;
					psi.ep[ep_index + 2] = ep;
				}
			}
		}
	}
}

void timefield3d_set_prism_electric_conductivity(timefield3d psi, float sig, int start_x, int end_x, int start_y, int end_y, int start_z, int end_z){
	
	int i, j, k, sig_index;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			for(k = 0; k < psi.k; k++){
				
				if((i >= start_x)&&(i < end_x)&&(j >= start_y)&&(j < end_y)&&(k >= start_z)&&(k < end_z)){
					
					sig_index = 6*(i + psi.i*(j + psi.j*k));
					
					psi.sig[sig_index] = sig;
					psi.sig[sig_index + 1] = sig;
					psi.sig[sig_index + 2] = sig;
				}
			}
		}
	}
}

void timefield3d_set_prism_complex_permittivity(timefield3d psi, double omega, float real_ep, float imag_ep, int start_x, int end_x, int start_y, int end_y, int start_z, int end_z){
	
	//set real permittivity
	timefield3d_set_prism_permittivity(psi, real_ep, start_x, end_x, start_y, end_y, start_z, end_z);
	
	//set imag permittivity
	timefield3d_set_prism_electric_conductivity(psi, omega*imag_ep, start_x, end_x, start_y, end_y, start_z, end_z);
}

void timefield3d_set_waveguide(timefield3d psi, float ep, int start_x, int end_x, int start_z, int end_z, int bottom, int top){

	timefield3d_set_prism_permittivity(psi, ep, start_x, end_x, bottom, top, start_z, end_z);
}

void timefield3d_set_taper(timefield3d psi, double ep, double mu, int start_left, int start_right, int end_left, int end_right, int start_z, int end_z, int bottom, int top){

	int i, j, k, index;
	double right_slope, left_slope;

	right_slope = (double)(end_right-start_right)/((double)(end_z-start_z));
	left_slope = (double)(end_left-start_left)/((double)(end_z-start_z));

	for(j = bottom; j < top; j++){
		for(k = start_z; k < end_z; k++){
			for(i = (start_left + round(left_slope*(k-start_z))); i < (start_right + round(right_slope*(k - start_z))); i++){

				index = 3*(i + psi.i*(j + psi.j*k));

				psi.ep[index] = ep;
				psi.mu[index] = mu;
				psi.ep[index + 1] = ep;
				psi.mu[index + 1] = mu;
				psi.ep[index + 2] = ep;
				psi.mu[index + 2] = mu;

			}
		}
	}
}

void timefield3d_apply_H_timestep(timefield3d psi, double delta_t){
	
	int i, j, k, ep_index, field_index, field_index_100, field_index_010, field_index_001;
	double divisor;
	double Hx, Hy, Hz, Mx, My, Mz;
	double Ex_000, Ex_010, Ex_001, Ey_000, Ey_100, Ey_001, Ez_000, Ez_100, Ez_010;
	double curl_x, curl_y, curl_z;
	double mu_x, mu_y, mu_z;
	double dEz_dy, dEy_dz, dEx_dz, dEz_dx, dEy_dx, dEx_dy;
	
	divisor = 1/(psi.disc);
	for(i = 1; i < psi.i - 1; i++){
		for(j = 1; j < psi.j - 1; j++){
			for(k = 1; k < psi.k - 1; k++){
		
				ep_index = 3*(i + psi.i*(j + psi.j*k));
				field_index = 6*(i + psi.i*(j + psi.j*k));
				field_index_100 = 6*((i+1) + psi.i*(j + psi.j*k));
				field_index_010 = 6*(i + psi.i*((j+1) + psi.j*k));
				field_index_001 = 6*(i + psi.i*(j + psi.j*(k+1)));
	
				Hx = psi.fields[field_index + 3];
				Hy = psi.fields[field_index + 4];
				Hz = psi.fields[field_index + 5];
		
				Mx = psi.sources[field_index + 3] + Hx*psi.sig[field_index + 3];
				My = psi.sources[field_index + 4] + Hy*psi.sig[field_index + 4];
				Mz = psi.sources[field_index + 5] + Hz*psi.sig[field_index + 5];
		
				Ex_000 = psi.fields[field_index];
				Ex_010 = psi.fields[field_index_010];
				Ex_001 = psi.fields[field_index_001];
		
				Ey_000 = psi.fields[field_index + 1];
				Ey_100 = psi.fields[field_index_100 + 1];
				Ey_001 = psi.fields[field_index_001 + 1];
		
				Ez_000 = psi.fields[field_index + 2];
				Ez_100 = psi.fields[field_index_100 + 2];
				Ez_010 = psi.fields[field_index_010 + 2];
				
				mu_x = psi.mu[ep_index];
				mu_y = psi.mu[ep_index + 1];
				mu_z = psi.mu[ep_index + 2];
		
				//take partials
				dEz_dy = (Ez_010 - Ez_000);
				dEy_dz = (Ey_001 - Ey_000);
				dEx_dz = (Ex_001 - Ex_000);
				dEz_dx = (Ez_100 - Ez_000);
				dEy_dx = (Ey_100 - Ey_000);
				dEx_dy = (Ex_010 - Ex_000);
				
				//take curls
				curl_x = divisor*(dEz_dy - dEy_dz);
				curl_y = divisor*(dEx_dz - dEz_dx);
				curl_z = divisor*(dEy_dx - dEx_dy);
	
				//Hx
				psi.fields[field_index + 3] = Hx - delta_t*(curl_x + Mx)/mu_x;
		
				//Hy
				psi.fields[field_index + 4] = Hy - delta_t*(curl_y + My)/mu_y;
		
				//Hz
				psi.fields[field_index + 5] = Hz - delta_t*(curl_z + Mz)/mu_z;			
			}
		}
	}
}

void timefield3d_apply_E_timestep(timefield3d psi, double delta_t){
	
	int i, j, k, ep_index, field_index, field_index_100, field_index_010, field_index_001;
	double divisor;
	double Ex, Ey, Ez, Jx, Jy, Jz;
	double Hx_000, Hx_010, Hx_001, Hy_000, Hy_100, Hy_001, Hz_000, Hz_100, Hz_010;
	double curl_x, curl_y, curl_z;
	double ep_x, ep_y, ep_z;
	double dHz_dy, dHy_dz, dHx_dz, dHz_dx, dHy_dx, dHx_dy;
	
	divisor = 1/(psi.disc);
	for(i = 1; i < psi.i - 1; i++){
		for(j = 1; j < psi.j - 1; j++){
			for(k = 1; k < psi.k - 1; k++){
		
				ep_index = 3*(i + psi.i*(j + psi.j*k));
				field_index = 6*(i + psi.i*(j + psi.j*k));
				field_index_100 = 6*((i-1) + psi.i*(j + psi.j*k));
				field_index_010 = 6*(i + psi.i*((j-1) + psi.j*k));
				field_index_001 = 6*(i + psi.i*(j + psi.j*(k-1)));
	
				Ex = psi.fields[field_index];
				Ey = psi.fields[field_index + 1];
				Ez = psi.fields[field_index + 2];
		
				Jx = psi.sources[field_index] + Ex*psi.sig[field_index];
				Jy = psi.sources[field_index + 1] + Ey*psi.sig[field_index + 1];
				Jz = psi.sources[field_index + 2] + Ez*psi.sig[field_index + 2];
		
				Hx_000 = psi.fields[field_index + 3];
				Hx_010 = psi.fields[field_index_010 + 3];
				Hx_001 = psi.fields[field_index_001 + 3];
		
				Hy_000 = psi.fields[field_index + 4];
				Hy_100 = psi.fields[field_index_100 + 4];
				Hy_001 = psi.fields[field_index_001 + 4];
		
				Hz_000 = psi.fields[field_index + 5];
				Hz_100 = psi.fields[field_index_100 + 5];
				Hz_010 = psi.fields[field_index_010 + 5];
				
				ep_x = psi.ep[ep_index];
				ep_y = psi.ep[ep_index + 1];
				ep_z = psi.ep[ep_index + 2];
		
				//take partials
				dHz_dy = (Hz_000 - Hz_010);
				dHy_dz = (Hy_000 - Hy_001);
				dHx_dz = (Hx_000 - Hx_001);
				dHz_dx = (Hz_000 - Hz_100);
				dHy_dx = (Hy_000 - Hy_100);
				dHx_dy = (Hx_000 - Hx_010);
				
				//take curls
				curl_x = divisor*(dHz_dy - dHy_dz);
				curl_y = divisor*(dHx_dz - dHz_dx);
				curl_z = divisor*(dHy_dx - dHx_dy);
	
				//Ex
				psi.fields[field_index] = Ex + delta_t*(curl_x - Jx)/ep_x;
		
				//Ey
				psi.fields[field_index + 1] = Ey + delta_t*(curl_y - Jy)/ep_y;
		
				//Ez
				psi.fields[field_index + 2] = Ez + delta_t*(curl_z - Jz)/ep_z;			
			}
		}
	}
}

void timefield3d_apply_H_timestep_zsection(timefield3d psi, double delta_t, int start_index, int stop_index){
	
	int i, j, k, ep_index, field_index, field_index_100, field_index_010, field_index_001;
	double divisor;
	double Hx, Hy, Hz, Mx, My, Mz;
	double Ex_000, Ex_010, Ex_001, Ey_000, Ey_100, Ey_001, Ez_000, Ez_100, Ez_010;
	double curl_x, curl_y, curl_z;
	double mu_x, mu_y, mu_z, sig_x, sig_y, sig_z;
	double dEz_dy, dEy_dz, dEx_dz, dEz_dx, dEy_dx, dEx_dy;
	double main_coeff, H_coeff;
	
	divisor = 1/(psi.disc);
	for(i = 1; i < psi.i - 1; i++){
		for(j = 1; j < psi.j - 1; j++){
			for(k = start_index + 1; k < stop_index - 1; k++){
		
				ep_index = 3*(i + psi.i*(j + psi.j*k));
				field_index = 6*(i + psi.i*(j + psi.j*k));
				field_index_100 = 6*((i+1) + psi.i*(j + psi.j*k));
				field_index_010 = 6*(i + psi.i*((j+1) + psi.j*k));
				field_index_001 = 6*(i + psi.i*(j + psi.j*(k+1)));
	
				Hx = psi.fields[field_index + 3];
				Hy = psi.fields[field_index + 4];
				Hz = psi.fields[field_index + 5];
		
				Mx = psi.sources[field_index + 3] + Hx*psi.sig[field_index + 3];
				My = psi.sources[field_index + 4] + Hy*psi.sig[field_index + 4];
				Mz = psi.sources[field_index + 5] + Hz*psi.sig[field_index + 5];
		
				Ex_000 = psi.fields[field_index];
				Ex_010 = psi.fields[field_index_010];
				Ex_001 = psi.fields[field_index_001];
		
				Ey_000 = psi.fields[field_index + 1];
				Ey_100 = psi.fields[field_index_100 + 1];
				Ey_001 = psi.fields[field_index_001 + 1];
		
				Ez_000 = psi.fields[field_index + 2];
				Ez_100 = psi.fields[field_index_100 + 2];
				Ez_010 = psi.fields[field_index_010 + 2];
				
				mu_x = psi.mu[ep_index];
				mu_y = psi.mu[ep_index + 1];
				mu_z = psi.mu[ep_index + 2];
				
				sig_x = psi.sig[field_index + 3];
				sig_y = psi.sig[field_index + 4];
				sig_z = psi.sig[field_index + 5];
		
				//take partials
				dEz_dy = (Ez_010 - Ez_000);
				dEy_dz = (Ey_001 - Ey_000);
				dEx_dz = (Ex_001 - Ex_000);
				dEz_dx = (Ez_100 - Ez_000);
				dEy_dx = (Ey_100 - Ey_000);
				dEx_dy = (Ex_010 - Ex_000);
				
				//take curls
				curl_x = divisor*(dEz_dy - dEy_dz);
				curl_y = divisor*(dEx_dz - dEz_dx);
				curl_z = divisor*(dEy_dx - dEx_dy);
				
				//update fields
				main_coeff = 1.0/(mu_x/delta_t + 0.5*sig_x);
				H_coeff = mu_x/delta_t - 0.5*sig_x;
				psi.fields[field_index + 3] = main_coeff*(H_coeff*Hx - curl_x - Mx);
				
				main_coeff = 1.0/(mu_y/delta_t + 0.5*sig_y);
				H_coeff = mu_y/delta_t - 0.5*sig_y;
				psi.fields[field_index + 4] = main_coeff*(H_coeff*Hy - curl_y - My);
				
				main_coeff = 1.0/(mu_z/delta_t + 0.5*sig_z);
				H_coeff = mu_z/delta_t - 0.5*sig_z;
				psi.fields[field_index + 5] = main_coeff*(H_coeff*Hz - curl_z - Mz);			
			}
		}
	}
}

void timefield3d_apply_E_timestep_zsection(timefield3d psi, double delta_t, int start_index, int stop_index){
	
	int i, j, k, ep_index, field_index, field_index_100, field_index_010, field_index_001;
	double divisor;
	double Ex, Ey, Ez, Jx, Jy, Jz;
	double Hx_000, Hx_010, Hx_001, Hy_000, Hy_100, Hy_001, Hz_000, Hz_100, Hz_010;
	double curl_x, curl_y, curl_z;
	double ep_x, ep_y, ep_z, sig_x, sig_y, sig_z;
	double dHz_dy, dHy_dz, dHx_dz, dHz_dx, dHy_dx, dHx_dy;
	double main_coeff, E_coeff;
	
	divisor = 1/(psi.disc);
	for(i = 1; i < psi.i - 1; i++){
		for(j = 1; j < psi.j - 1; j++){
			for(k = start_index + 1; k < stop_index - 1; k++){
		
				ep_index = 3*(i + psi.i*(j + psi.j*k));
				field_index = 6*(i + psi.i*(j + psi.j*k));
				field_index_100 = 6*((i-1) + psi.i*(j + psi.j*k));
				field_index_010 = 6*(i + psi.i*((j-1) + psi.j*k));
				field_index_001 = 6*(i + psi.i*(j + psi.j*(k-1)));
	
				Ex = psi.fields[field_index];
				Ey = psi.fields[field_index + 1];
				Ez = psi.fields[field_index + 2];
		
				Jx = psi.sources[field_index];
				Jy = psi.sources[field_index + 1];
				Jz = psi.sources[field_index + 2];
		
				Hx_000 = psi.fields[field_index + 3];
				Hx_010 = psi.fields[field_index_010 + 3];
				Hx_001 = psi.fields[field_index_001 + 3];
		
				Hy_000 = psi.fields[field_index + 4];
				Hy_100 = psi.fields[field_index_100 + 4];
				Hy_001 = psi.fields[field_index_001 + 4];
		
				Hz_000 = psi.fields[field_index + 5];
				Hz_100 = psi.fields[field_index_100 + 5];
				Hz_010 = psi.fields[field_index_010 + 5];
				
				ep_x = psi.ep[ep_index];
				ep_y = psi.ep[ep_index + 1];
				ep_z = psi.ep[ep_index + 2];
				
				sig_x = psi.sig[field_index];
				sig_y = psi.sig[field_index + 1];
				sig_z = psi.sig[field_index + 2];
		
				//take partials
				dHz_dy = (Hz_000 - Hz_010);
				dHy_dz = (Hy_000 - Hy_001);
				dHx_dz = (Hx_000 - Hx_001);
				dHz_dx = (Hz_000 - Hz_100);
				dHy_dx = (Hy_000 - Hy_100);
				dHx_dy = (Hx_000 - Hx_010);
				
				//take curls
				curl_x = divisor*(dHz_dy - dHy_dz);
				curl_y = divisor*(dHx_dz - dHz_dx);
				curl_z = divisor*(dHy_dx - dHx_dy);
	
				//update fields
				main_coeff = 1.0/(ep_x/delta_t + 0.5*sig_x);
				E_coeff = ep_x/delta_t - 0.5*sig_x;
				psi.fields[field_index] = main_coeff*(E_coeff*Ex + curl_x - Jx);
				
				main_coeff = 1.0/(ep_y/delta_t + 0.5*sig_y);
				E_coeff = ep_y/delta_t - 0.5*sig_y;
				psi.fields[field_index + 1] = main_coeff*(E_coeff*Ey + curl_y - Jy);
				
				main_coeff = 1.0/(ep_z/delta_t + 0.5*sig_z);
				E_coeff = ep_z/delta_t - 0.5*sig_z;
				psi.fields[field_index + 2] = main_coeff*(E_coeff*Ez + curl_z - Jz);
			}
		}
	}
}

void timefield3d_apply_H_first_order_abc(timefield3d psi, double *old_field, double delta_t){
	
	int i, j, k, ep_index, field_index, field_index_inner;
	double Hx_old, Hx_inner, Hx_inner_old, Hy_old, Hy_inner, Hy_inner_old, Hz_old, Hz_inner, Hz_inner_old;
	double ep_x, ep_y, ep_z, mu_x, mu_y, mu_z;
	double speed, scalar;
	
	//x = 0 boundary
	i = 0;
	for(j = 0; j < psi.j; j++){
		for(k = 0; k < psi.k; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i+1) + psi.i*(j + psi.j*k));
			
			Hy_old = old_field[field_index + 4];
			Hy_inner = psi.fields[field_index_inner + 4];
			Hy_inner_old = old_field[field_index_inner + 4];
			
			Hz_old = old_field[field_index + 5];
			Hz_inner = psi.fields[field_index_inner + 5];
			Hz_inner_old = old_field[field_index_inner + 5];
			
			ep_y = psi.ep[ep_index + 1];
			ep_z = psi.ep[ep_index + 2];
			
			mu_y = psi.mu[ep_index + 1];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Hy-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 4] = Hy_inner_old + scalar*(Hy_inner - Hy_old);
			
			//Absorb Hz-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 5] = Hz_inner_old + scalar*(Hz_inner - Hz_old);
		}
	}
	
	//y = 0 boundary
	j = 0;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j+1) + psi.j*k));
			
			Hx_old = old_field[field_index + 3];
			Hx_inner = psi.fields[field_index_inner + 3];
			Hx_inner_old = old_field[field_index_inner + 3];
			
			Hz_old = old_field[field_index + 5];
			Hz_inner = psi.fields[field_index_inner + 5];
			Hz_inner_old = old_field[field_index_inner + 5];
			
			ep_x = psi.ep[ep_index];
			ep_z = psi.ep[ep_index + 2];
			
			mu_x = psi.mu[ep_index];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Hx-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 3] = Hx_inner_old + scalar*(Hx_inner - Hx_old);
			
			//Absorb Hz-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 5] = Hz_inner_old + scalar*(Hz_inner - Hz_old);
		}
	}
	
	//z = 0 boundary
	k = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*(j + psi.j*(k+1)));
			
			Hx_old = old_field[field_index + 3];
			Hx_inner = psi.fields[field_index_inner + 3];
			Hx_inner_old = old_field[field_index_inner + 3];
			
			Hy_old = old_field[field_index + 4];
			Hy_inner = psi.fields[field_index_inner + 4];
			Hy_inner_old = old_field[field_index_inner + 4];
			
			ep_x = psi.ep[ep_index];
			ep_y = psi.ep[ep_index + 1];
			
			mu_x = psi.mu[ep_index];
			mu_y = psi.mu[ep_index + 1];
			
			//Absorb Hx-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 3] = Hx_inner_old + scalar*(Hx_inner - Hx_old);
			
			//Absorb Hy-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 4] = Hy_inner_old + scalar*(Hy_inner - Hy_old);
		}
	}
	
	//x = end boundary
	i = psi.i - 1;
	for(j = 0; j < psi.j; j++){
		for(k = 0; k < psi.k; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i-1) + psi.i*(j + psi.j*k));
			
			Hy_old = old_field[field_index + 4];
			Hy_inner = psi.fields[field_index_inner + 4];
			Hy_inner_old = old_field[field_index_inner + 4];
			
			Hz_old = old_field[field_index + 5];
			Hz_inner = psi.fields[field_index_inner + 5];
			Hz_inner_old = old_field[field_index_inner + 5];
			
			ep_y = psi.ep[ep_index + 1];
			ep_z = psi.ep[ep_index + 2];
			
			mu_y = psi.mu[ep_index + 1];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Hy-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 4] = Hy_inner_old + scalar*(Hy_inner - Hy_old);
			
			//Absorb Hz-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 5] = Hz_inner_old + scalar*(Hz_inner - Hz_old);
		}
	}
	
	//y = end boundary
	j = psi.j - 1;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j-1) + psi.j*k));
			
			Hx_old = old_field[field_index + 3];
			Hx_inner = psi.fields[field_index_inner + 3];
			Hx_inner_old = old_field[field_index_inner + 3];
			
			Hz_old = old_field[field_index + 5];
			Hz_inner = psi.fields[field_index_inner + 5];
			Hz_inner_old = old_field[field_index_inner + 5];
			
			ep_x = psi.ep[ep_index];
			ep_z = psi.ep[ep_index + 2];
			
			mu_x = psi.mu[ep_index];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Hx-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 3] = Hx_inner_old + scalar*(Hx_inner - Hx_old);
			
			//Absorb Hz-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 5] = Hz_inner_old + scalar*(Hz_inner - Hz_old);
		}
	}
	
	//z = end boundary
	k = psi.k - 1;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*(j + psi.j*(k-1)));
			
			Hx_old = old_field[field_index + 3];
			Hx_inner = psi.fields[field_index_inner + 3];
			Hx_inner_old = old_field[field_index_inner + 3];
			
			Hy_old = old_field[field_index + 4];
			Hy_inner = psi.fields[field_index_inner + 4];
			Hy_inner_old = old_field[field_index_inner + 4];
			
			ep_x = psi.ep[ep_index];
			ep_y = psi.ep[ep_index + 1];
			
			mu_x = psi.mu[ep_index];
			mu_y = psi.mu[ep_index + 1];
			
			//Absorb Hx-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 3] = Hx_inner_old + scalar*(Hx_inner - Hx_old);
			
			//Absorb Hy-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 4] = Hy_inner_old + scalar*(Hy_inner - Hy_old);
		}
	}
}

void timefield3d_apply_E_first_order_abc(timefield3d psi, double *old_field, double delta_t){
	
	int i, j, k, ep_index, field_index, field_index_inner;
	double Ex_old, Ex_inner, Ex_inner_old, Ey_old, Ey_inner, Ey_inner_old, Ez_old, Ez_inner, Ez_inner_old;
	double ep_x, ep_y, ep_z, mu_x, mu_y, mu_z;
	double speed, scalar;
	
	//x = 0 boundary
	i = 0;
	for(j = 0; j < psi.j; j++){
		for(k = 0; k < psi.k; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i+1) + psi.i*(j + psi.j*k));
			
			Ey_old = old_field[field_index + 1];
			Ey_inner = psi.fields[field_index_inner + 1];
			Ey_inner_old = old_field[field_index_inner + 1];
			
			Ez_old = old_field[field_index + 2];
			Ez_inner = psi.fields[field_index_inner + 2];
			Ez_inner_old = old_field[field_index_inner + 2];
			
			ep_y = psi.ep[ep_index + 1];
			ep_z = psi.ep[ep_index + 2];

			mu_y = psi.mu[ep_index + 1];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Ey-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 1] = Ey_inner_old + scalar*(Ey_inner - Ey_old);
			
			//Absorb Ez-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 2] = Ez_inner_old + scalar*(Ez_inner - Ez_old);
		}
	}
	
	//y = 0 boundary
	j = 0;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j+1) + psi.j*k));
			
			Ex_old = old_field[field_index];
			Ex_inner = psi.fields[field_index_inner];
			Ex_inner_old = old_field[field_index_inner];
			
			Ez_old = old_field[field_index + 2];
			Ez_inner = psi.fields[field_index_inner + 2];
			Ez_inner_old = old_field[field_index_inner + 2];
			
			ep_x = psi.ep[ep_index];
			ep_z = psi.ep[ep_index + 2];
			
			mu_x = psi.mu[ep_index];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Ex-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index] = Ex_inner_old + scalar*(Ex_inner - Ex_old);
			
			//Absorb Ez-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 2] = Ez_inner_old + scalar*(Ez_inner - Ez_old);
		}
	}
	
	//z = 0 boundary
	k = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*(j + psi.j*(k+1)));
			
			Ex_old = old_field[field_index];
			Ex_inner = psi.fields[field_index_inner];
			Ex_inner_old = old_field[field_index_inner];
			
			Ey_old = old_field[field_index + 1];
			Ey_inner = psi.fields[field_index_inner + 1];
			Ey_inner_old = old_field[field_index_inner + 1];
			
			ep_x = psi.ep[ep_index];
			ep_y = psi.ep[ep_index + 1];
			
			mu_x = psi.mu[ep_index];
			mu_y = psi.mu[ep_index + 1];
			
			//Absorb Ex-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index] = Ex_inner_old + scalar*(Ex_inner - Ex_old);
			
			//Absorb Ey-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 1] = Ey_inner_old + scalar*(Ey_inner - Ey_old);
		}
	}
	
	//x = end boundary
	i = psi.i - 1;
	for(j = 0; j < psi.j; j++){
		for(k = 0; k < psi.k; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i-1) + psi.i*(j + psi.j*k));
			
			Ey_old = old_field[field_index + 1];
			Ey_inner = psi.fields[field_index_inner + 1];
			Ey_inner_old = old_field[field_index_inner + 1];
			
			Ez_old = old_field[field_index + 2];
			Ez_inner = psi.fields[field_index_inner + 2];
			Ez_inner_old = old_field[field_index_inner + 2];
			
			ep_y = psi.ep[ep_index + 1];
			ep_z = psi.ep[ep_index + 2];

			mu_y = psi.mu[ep_index + 1];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Ey-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 1] = Ey_inner_old + scalar*(Ey_inner - Ey_old);
			
			//Absorb Ez-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 2] = Ez_inner_old + scalar*(Ez_inner - Ez_old);
		}
	}
	
	//y = end boundary
	j = psi.j - 1;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j-1) + psi.j*k));
			
			Ex_old = old_field[field_index];
			Ex_inner = psi.fields[field_index_inner];
			Ex_inner_old = old_field[field_index_inner];
			
			Ez_old = old_field[field_index + 2];
			Ez_inner = psi.fields[field_index_inner + 2];
			Ez_inner_old = old_field[field_index_inner + 2];
			
			ep_x = psi.ep[ep_index];
			ep_z = psi.ep[ep_index + 2];
			
			mu_x = psi.mu[ep_index];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Ex-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index] = Ex_inner_old + scalar*(Ex_inner - Ex_old);
			
			//Absorb Ez-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 2] = Ez_inner_old + scalar*(Ez_inner - Ez_old);
		}
	}
	
	//z = end boundary
	k = psi.k - 1;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*(j + psi.j*(k-1)));
			
			Ex_old = old_field[field_index];
			Ex_inner = psi.fields[field_index_inner];
			Ex_inner_old = old_field[field_index_inner];
			
			Ey_old = old_field[field_index + 1];
			Ey_inner = psi.fields[field_index_inner + 1];
			Ey_inner_old = old_field[field_index_inner + 1];
			
			ep_x = psi.ep[ep_index];
			ep_y = psi.ep[ep_index + 1];
			
			mu_x = psi.mu[ep_index];
			mu_y = psi.mu[ep_index + 1];
			
			//Absorb Ex-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index] = Ex_inner_old + scalar*(Ex_inner - Ex_old);
			
			//Absorb Ey-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 1] = Ey_inner_old + scalar*(Ey_inner - Ey_old);
		}
	}
}

void timefield3d_apply_H_first_order_abc_zsection(timefield3d psi, double delta_t, int start_index, int stop_index){
	
	int i, j, k, ep_index, field_index, field_index_inner;
	double Hx_old, Hx_inner, Hx_inner_old, Hy_old, Hy_inner, Hy_inner_old, Hz_old, Hz_inner, Hz_inner_old;
	double ep_x, ep_y, ep_z, mu_x, mu_y, mu_z;
	double speed, scalar;
	
	//x = 0 boundary
	i = 0;
	for(j = 0; j < psi.j; j++){
		for(k = start_index; k < stop_index; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i+1) + psi.i*(j + psi.j*k));
			
			Hy_old = psi.old_fields[field_index + 4];
			Hy_inner = psi.fields[field_index_inner + 4];
			Hy_inner_old = psi.old_fields[field_index_inner + 4];
			
			Hz_old = psi.old_fields[field_index + 5];
			Hz_inner = psi.fields[field_index_inner + 5];
			Hz_inner_old = psi.old_fields[field_index_inner + 5];
			
			ep_y = psi.ep[ep_index + 1];
			ep_z = psi.ep[ep_index + 2];
			
			mu_y = psi.mu[ep_index + 1];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Hy-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 4] = Hy_inner_old + scalar*(Hy_inner - Hy_old);
			
			//Absorb Hz-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 5] = Hz_inner_old + scalar*(Hz_inner - Hz_old);
		}
	}
	
	//y = 0 boundary
	j = 0;
	for(i = 0; i < psi.i; i++){
		for(k = start_index; k < stop_index; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j+1) + psi.j*k));
			
			Hx_old = psi.old_fields[field_index + 3];
			Hx_inner = psi.fields[field_index_inner + 3];
			Hx_inner_old = psi.old_fields[field_index_inner + 3];
			
			Hz_old = psi.old_fields[field_index + 5];
			Hz_inner = psi.fields[field_index_inner + 5];
			Hz_inner_old = psi.old_fields[field_index_inner + 5];
			
			ep_x = psi.ep[ep_index];
			ep_z = psi.ep[ep_index + 2];
			
			mu_x = psi.mu[ep_index];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Hx-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 3] = Hx_inner_old + scalar*(Hx_inner - Hx_old);
			
			//Absorb Hz-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 5] = Hz_inner_old + scalar*(Hz_inner - Hz_old);
		}
	}
	
	//z = 0 boundary
	if(start_index == 0){

		k = 0;
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){
				
				ep_index = 3*(i + psi.i*(j + psi.j*k));
				field_index = 6*(i + psi.i*(j + psi.j*k));
				field_index_inner = 6*(i + psi.i*(j + psi.j*(k+1)));
				
				Hx_old = psi.old_fields[field_index + 3];
				Hx_inner = psi.fields[field_index_inner + 3];
				Hx_inner_old = psi.old_fields[field_index_inner + 3];
				
				Hy_old = psi.old_fields[field_index + 4];
				Hy_inner = psi.fields[field_index_inner + 4];
				Hy_inner_old = psi.old_fields[field_index_inner + 4];
				
				ep_x = psi.ep[ep_index];
				ep_y = psi.ep[ep_index + 1];
				
				mu_x = psi.mu[ep_index];
				mu_y = psi.mu[ep_index + 1];
				
				//Absorb Hx-polarized, normal-incidence plane waves
				speed = 1/sqrt(ep_y*mu_x);
				scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
				psi.fields[field_index + 3] = Hx_inner_old + scalar*(Hx_inner - Hx_old);
				
				//Absorb Hy-polarized, normal-incidence plane waves
				speed = 1/sqrt(ep_x*mu_y);
				scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
				psi.fields[field_index + 4] = Hy_inner_old + scalar*(Hy_inner - Hy_old);
			}
		}
	}
	
	//x = end boundary
	i = psi.i - 1;
	for(j = 0; j < psi.j; j++){
		for(k = start_index; k < stop_index; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i-1) + psi.i*(j + psi.j*k));
			
			Hy_old = psi.old_fields[field_index + 4];
			Hy_inner = psi.fields[field_index_inner + 4];
			Hy_inner_old = psi.old_fields[field_index_inner + 4];
			
			Hz_old = psi.old_fields[field_index + 5];
			Hz_inner = psi.fields[field_index_inner + 5];
			Hz_inner_old = psi.old_fields[field_index_inner + 5];
			
			ep_y = psi.ep[ep_index + 1];
			ep_z = psi.ep[ep_index + 2];
			
			mu_y = psi.mu[ep_index + 1];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Hy-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 4] = Hy_inner_old + scalar*(Hy_inner - Hy_old);
			
			//Absorb Hz-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 5] = Hz_inner_old + scalar*(Hz_inner - Hz_old);
		}
	}
	
	//y = end boundary
	j = psi.j - 1;
	for(i = 0; i < psi.i; i++){
		for(k = start_index; k < stop_index; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j-1) + psi.j*k));
			
			Hx_old = psi.old_fields[field_index + 3];
			Hx_inner = psi.fields[field_index_inner + 3];
			Hx_inner_old = psi.old_fields[field_index_inner + 3];
			
			Hz_old = psi.old_fields[field_index + 5];
			Hz_inner = psi.fields[field_index_inner + 5];
			Hz_inner_old = psi.old_fields[field_index_inner + 5];
			
			ep_x = psi.ep[ep_index];
			ep_z = psi.ep[ep_index + 2];
			
			mu_x = psi.mu[ep_index];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Hx-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 3] = Hx_inner_old + scalar*(Hx_inner - Hx_old);
			
			//Absorb Hz-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 5] = Hz_inner_old + scalar*(Hz_inner - Hz_old);
		}
	}
	
	//z = end boundary
	if(stop_index == psi.k){

		k = psi.k - 1;
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){
				
				ep_index = 3*(i + psi.i*(j + psi.j*k));
				field_index = 6*(i + psi.i*(j + psi.j*k));
				field_index_inner = 6*(i + psi.i*(j + psi.j*(k-1)));
				
				Hx_old = psi.old_fields[field_index + 3];
				Hx_inner = psi.fields[field_index_inner + 3];
				Hx_inner_old = psi.old_fields[field_index_inner + 3];
				
				Hy_old = psi.old_fields[field_index + 4];
				Hy_inner = psi.fields[field_index_inner + 4];
				Hy_inner_old = psi.old_fields[field_index_inner + 4];
				
				ep_x = psi.ep[ep_index];
				ep_y = psi.ep[ep_index + 1];
				
				mu_x = psi.mu[ep_index];
				mu_y = psi.mu[ep_index + 1];
				
				//Absorb Hx-polarized, normal-incidence plane waves
				speed = 1/sqrt(ep_y*mu_x);
				scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
				psi.fields[field_index + 3] = Hx_inner_old + scalar*(Hx_inner - Hx_old);
				
				//Absorb Hy-polarized, normal-incidence plane waves
				speed = 1/sqrt(ep_x*mu_y);
				scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
				psi.fields[field_index + 4] = Hy_inner_old + scalar*(Hy_inner - Hy_old);
			}
		}
	}
}

void timefield3d_apply_E_first_order_abc_zsection(timefield3d psi, double delta_t, int start_index, int stop_index){
	
	int i, j, k, ep_index, field_index, field_index_inner;
	double Ex_old, Ex_inner, Ex_inner_old, Ey_old, Ey_inner, Ey_inner_old, Ez_old, Ez_inner, Ez_inner_old;
	double ep_x, ep_y, ep_z, mu_x, mu_y, mu_z;
	double speed, scalar;
	
	//x = 0 boundary
	i = 0;
	for(j = 0; j < psi.j; j++){
		for(k = start_index; k < stop_index; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i+1) + psi.i*(j + psi.j*k));
			
			Ey_old = psi.old_fields[field_index + 1];
			Ey_inner = psi.fields[field_index_inner + 1];
			Ey_inner_old = psi.old_fields[field_index_inner + 1];
			
			Ez_old = psi.old_fields[field_index + 2];
			Ez_inner = psi.fields[field_index_inner + 2];
			Ez_inner_old = psi.old_fields[field_index_inner + 2];
			
			ep_y = psi.ep[ep_index + 1];
			ep_z = psi.ep[ep_index + 2];

			mu_y = psi.mu[ep_index + 1];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Ey-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 1] = Ey_inner_old + scalar*(Ey_inner - Ey_old);
			
			//Absorb Ez-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 2] = Ez_inner_old + scalar*(Ez_inner - Ez_old);
		}
	}
	
	//y = 0 boundary
	j = 0;
	for(i = 0; i < psi.i; i++){
		for(k = start_index; k < stop_index; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j+1) + psi.j*k));
			
			Ex_old = psi.old_fields[field_index];
			Ex_inner = psi.fields[field_index_inner];
			Ex_inner_old = psi.old_fields[field_index_inner];
			
			Ez_old = psi.old_fields[field_index + 2];
			Ez_inner = psi.fields[field_index_inner + 2];
			Ez_inner_old = psi.old_fields[field_index_inner + 2];
			
			ep_x = psi.ep[ep_index];
			ep_z = psi.ep[ep_index + 2];
			
			mu_x = psi.mu[ep_index];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Ex-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index] = Ex_inner_old + scalar*(Ex_inner - Ex_old);
			
			//Absorb Ez-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 2] = Ez_inner_old + scalar*(Ez_inner - Ez_old);
		}
	}
	
	//z = 0 boundary
	if(start_index == 0){

		k = 0;
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){
				
				ep_index = 3*(i + psi.i*(j + psi.j*k));
				field_index = 6*(i + psi.i*(j + psi.j*k));
				field_index_inner = 6*(i + psi.i*(j + psi.j*(k+1)));
				
				Ex_old = psi.old_fields[field_index];
				Ex_inner = psi.fields[field_index_inner];
				Ex_inner_old = psi.old_fields[field_index_inner];
				
				Ey_old = psi.old_fields[field_index + 1];
				Ey_inner = psi.fields[field_index_inner + 1];
				Ey_inner_old = psi.old_fields[field_index_inner + 1];
				
				ep_x = psi.ep[ep_index];
				ep_y = psi.ep[ep_index + 1];
				
				mu_x = psi.mu[ep_index];
				mu_y = psi.mu[ep_index + 1];
				
				//Absorb Ex-polarized, normal-incidence plane waves
				speed = 1/sqrt(ep_x*mu_y);
				scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
				psi.fields[field_index] = Ex_inner_old + scalar*(Ex_inner - Ex_old);
				
				//Absorb Ey-polarized, normal-incidence plane waves
				speed = 1/sqrt(ep_y*mu_x);
				scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
				psi.fields[field_index + 1] = Ey_inner_old + scalar*(Ey_inner - Ey_old);
			}
		}
	}
	
	//x = end boundary
	i = psi.i - 1;
	for(j = 0; j < psi.j; j++){
		for(k = start_index; k < stop_index; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*((i-1) + psi.i*(j + psi.j*k));
			
			Ey_old = psi.old_fields[field_index + 1];
			Ey_inner = psi.fields[field_index_inner + 1];
			Ey_inner_old = psi.old_fields[field_index_inner + 1];
			
			Ez_old = psi.old_fields[field_index + 2];
			Ez_inner = psi.fields[field_index_inner + 2];
			Ez_inner_old = psi.old_fields[field_index_inner + 2];
			
			ep_y = psi.ep[ep_index + 1];
			ep_z = psi.ep[ep_index + 2];

			mu_y = psi.mu[ep_index + 1];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Ey-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_y*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 1] = Ey_inner_old + scalar*(Ey_inner - Ey_old);
			
			//Absorb Ez-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_y);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 2] = Ez_inner_old + scalar*(Ez_inner - Ez_old);
		}
	}
	
	//y = end boundary
	j = psi.j - 1;
	for(i = 0; i < psi.i; i++){
		for(k = start_index; k < stop_index; k++){
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			field_index = 6*(i + psi.i*(j + psi.j*k));
			field_index_inner = 6*(i + psi.i*((j-1) + psi.j*k));
			
			Ex_old = psi.old_fields[field_index];
			Ex_inner = psi.fields[field_index_inner];
			Ex_inner_old = psi.old_fields[field_index_inner];
			
			Ez_old = psi.old_fields[field_index + 2];
			Ez_inner = psi.fields[field_index_inner + 2];
			Ez_inner_old = psi.old_fields[field_index_inner + 2];
			
			ep_x = psi.ep[ep_index];
			ep_z = psi.ep[ep_index + 2];
			
			mu_x = psi.mu[ep_index];
			mu_z = psi.mu[ep_index + 2];
			
			//Absorb Ex-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_x*mu_z);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index] = Ex_inner_old + scalar*(Ex_inner - Ex_old);
			
			//Absorb Ez-polarized, normal-incidence plane waves
			speed = 1/sqrt(ep_z*mu_x);
			scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
			psi.fields[field_index + 2] = Ez_inner_old + scalar*(Ez_inner - Ez_old);
		}
	}
	
	//z = end boundary
	if(stop_index == psi.k){

		k = psi.k - 1;
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){
				
				ep_index = 3*(i + psi.i*(j + psi.j*k));
				field_index = 6*(i + psi.i*(j + psi.j*k));
				field_index_inner = 6*(i + psi.i*(j + psi.j*(k-1)));
				
				Ex_old = psi.old_fields[field_index];
				Ex_inner = psi.fields[field_index_inner];
				Ex_inner_old = psi.old_fields[field_index_inner];
				
				Ey_old = psi.old_fields[field_index + 1];
				Ey_inner = psi.fields[field_index_inner + 1];
				Ey_inner_old = psi.old_fields[field_index_inner + 1];
				
				ep_x = psi.ep[ep_index];
				ep_y = psi.ep[ep_index + 1];
				
				mu_x = psi.mu[ep_index];
				mu_y = psi.mu[ep_index + 1];
				
				//Absorb Ex-polarized, normal-incidence plane waves
				speed = 1/sqrt(ep_x*mu_y);
				scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
				psi.fields[field_index] = Ex_inner_old + scalar*(Ex_inner - Ex_old);
				
				//Absorb Ey-polarized, normal-incidence plane waves
				speed = 1/sqrt(ep_y*mu_x);
				scalar = (speed*delta_t - psi.disc)/(speed*delta_t + psi.disc);
				psi.fields[field_index + 1] = Ey_inner_old + scalar*(Ey_inner - Ey_old);
			}
		}
	}
}

void timefield3d_apply_timestep(timefield3d psi, double delta_t){
	
	double *old_field;
	
	old_field = timefield3d_get_edge_fields(psi);

	timefield3d_apply_E_timestep(psi, delta_t);
	timefield3d_apply_E_first_order_abc(psi, old_field, delta_t);
	timefield3d_apply_H_timestep(psi, delta_t);
	timefield3d_apply_H_first_order_abc(psi, old_field, delta_t);
	
	free(old_field);
}

void timefield3d_apply_timestep_omp(timefield3d psi, double delta_t){
	
	double *old_field;
	int tid, nthreads;
	int start_index, stop_index;
	
	old_field = timefield3d_get_edge_fields(psi);
	
	#pragma omp parallel private(nthreads, tid)
	{
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
			
		start_index = tid*psi.k/nthreads;
		stop_index = ((tid+1)*psi.k)/nthreads;
		if(tid != 0){
				
			start_index--;
		}
		if(tid != (nthreads - 1)){
				
			stop_index++;
		}
		timefield3d_apply_E_timestep_zsection(psi, delta_t, start_index, stop_index);
	}
				
	timefield3d_apply_E_first_order_abc(psi, old_field, delta_t);
	
	#pragma omp parallel private(nthreads, tid)
	{
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
	
		start_index = tid*psi.k/nthreads;
		stop_index = ((tid+1)*psi.k)/nthreads;
		if(tid != 0){
				
			start_index--;
		}
		if(tid != (nthreads - 1)){
				
			stop_index++;
		}
		timefield3d_apply_H_timestep_zsection(psi, delta_t, start_index, stop_index);
	}
				
	timefield3d_apply_H_first_order_abc(psi, old_field, delta_t);
	free(old_field);
}

void timefield3d_apply_timestep_omp_full(timefield3d psi, double delta_t){
	
	int tid, nthreads;
	int start_index, stop_index;
	
	timefield3d_update_old_edge_fields(psi);
	
	#pragma omp parallel private(nthreads, tid)
	{
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
			
		start_index = tid*psi.k/nthreads;
		stop_index = ((tid+1)*psi.k)/nthreads;
		if(tid != 0){
				
			start_index--;
		}
		if(tid != (nthreads - 1)){
				
			stop_index++;
		}
		timefield3d_apply_E_timestep_zsection(psi, delta_t, start_index, stop_index);
	}

	#pragma omp parallel private(nthreads, tid)
	{
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
			
		start_index = tid*psi.k/nthreads;
		stop_index = ((tid+1)*psi.k)/nthreads;
		timefield3d_apply_E_first_order_abc_zsection(psi, delta_t, start_index, stop_index);
	}
	
	#pragma omp parallel private(nthreads, tid)
	{
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
	
		start_index = tid*psi.k/nthreads;
		stop_index = ((tid+1)*psi.k)/nthreads;
		if(tid != 0){
				
			start_index--;
		}
		if(tid != (nthreads - 1)){
				
			stop_index++;
		}
		timefield3d_apply_H_timestep_zsection(psi, delta_t, start_index, stop_index);
	}
				
	#pragma omp parallel private(nthreads, tid)
	{
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
			
		start_index = tid*psi.k/nthreads;
		stop_index = ((tid+1)*psi.k)/nthreads;
		timefield3d_apply_H_first_order_abc_zsection(psi, delta_t, start_index, stop_index);
	}
}

void timefield3d_update_fourier_fields(timefield3d psi, double phase){

	int i, j, k, m, field_index;
	double reduced_phase;
	complex twiddle;

	reduced_phase = phase - ((int)(phase/(2*PI)))*2*PI;
	twiddle = complex_exp(complex_from_doubles(0,reduced_phase));

	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			for(k = 0; k < psi.k; k++){

				field_index = 6*(i + psi.i*(j + psi.j*k));

				for(m = 0; m < 6; m++){

					psi.fourier_fields[field_index + m] = complex_add(psi.fourier_fields[field_index + m], complex_smul(psi.fields[field_index + m], twiddle));
				}
			}
		}
	}
}

void timefield3d_update_fourier_fields_y_slice(timefield3d psi, double phase, int slice_index){

	int i, j, k, m, field_index;
	double reduced_phase;
	complex twiddle;

	reduced_phase = phase - ((int)(phase/(2*PI)))*2*PI;
	twiddle = complex_exp(complex_from_doubles(0,reduced_phase));

	j = slice_index;
	for(i = 0; i < psi.i; i++){
		for(k = 0; k < psi.k; k++){

			field_index = 6*(i + psi.i*(j + psi.j*k));

			for(m = 0; m < 6; m++){

				psi.fourier_fields[field_index + m] = complex_add(psi.fourier_fields[field_index + m], complex_smul(psi.fields[field_index + m], twiddle));
			}
		}
	}
}

void timefield3d_update_fourier_fields_zsection(timefield3d psi, double phase, int start_index, int stop_index){

	int i, j, k, m, field_index;
	double reduced_phase;
	complex twiddle;

	reduced_phase = phase - ((int)(phase/(2*PI)))*2*PI;
	twiddle = complex_exp(complex_from_doubles(0,reduced_phase));

	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			for(k = start_index; k < stop_index; k++){

				field_index = 6*(i + psi.i*(j + psi.j*k));

				for(m = 0; m < 6; m++){

					psi.fourier_fields[field_index + m] = complex_add(psi.fourier_fields[field_index + m], complex_smul(psi.fields[field_index + m], twiddle));
				}
			}
		}
	}
}

void timefield3d_update_fourier_fields_y_slice_zsection(timefield3d psi, double phase, int slice_index, int start_index, int stop_index){

	int i, j, k, m, field_index;
	double reduced_phase;
	complex twiddle;

	reduced_phase = phase - ((int)(phase/(2*PI)))*2*PI;
	twiddle = complex_exp(complex_from_doubles(0,reduced_phase));

	j = slice_index;
	for(i = 0; i < psi.i; i++){
		for(k = start_index; k < stop_index; k++){

			field_index = 6*(i + psi.i*(j + psi.j*k));

			for(m = 0; m < 6; m++){

				psi.fourier_fields[field_index + m] = complex_add(psi.fourier_fields[field_index + m], complex_smul(psi.fields[field_index + m], twiddle));
			}
		}
	}
}

void timefield3d_update_fourier_fields_omp(timefield3d psi, double phase){
	
	int tid, nthreads;
	int start_index, stop_index;
	
	#pragma omp parallel private(nthreads, tid)
	{
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
			
		start_index = tid*psi.k/nthreads;
		stop_index = ((tid+1)*psi.k)/nthreads;
		timefield3d_update_fourier_fields_zsection(psi, phase, start_index, stop_index);
	}
}

void timefield3d_update_fourier_fields_y_slice_omp(timefield3d psi, double phase, int slice_index){
	
	int tid, nthreads;
	int start_index, stop_index;
	
	#pragma omp parallel private(nthreads, tid)
	{
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
			
		start_index = tid*psi.k/nthreads;
		stop_index = ((tid+1)*psi.k)/nthreads;
		timefield3d_update_fourier_fields_y_slice_zsection(psi, phase, slice_index, start_index, stop_index);
	}
}

void timefield3d_mode_source(timefield3d psi, int center_x, int center_y, int center_z, mode source_mode, double phase, double amplitude, int reverse){

	int i, j, k, m, ep_index, source_index, mode_x_index, mode_y_index, mode_index, M_sign;
	double x_distance, y_distance;
	double slope, inv_slope;
	double norm, inv_norm;
	complex mode_Ex, mode_Ey, mode_Hx, mode_Hy;
	double permittivity, permeability;

	M_sign = -1 + 2*reverse;

	//power-normalize
	mode_power_normalize(source_mode);

	//set sources
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			k = center_z;

			x_distance = (i - center_x)*(psi.disc/source_mode.disc);
			y_distance = (j - center_y)*(psi.disc/source_mode.disc);

			for(m = 0; m < 6; m++){

				source_index = 6*(i + psi.i*(j + k*(psi.j))) + m;
				mode_x_index = (int)x_distance + source_mode.i/2;
				mode_y_index = (int)y_distance + source_mode.j/2;
				ep_index = 3*(i + psi.i*(j + psi.j*k)) + m%3;

				if((mode_x_index >= 0)&&(mode_x_index < source_mode.i)&&(mode_y_index >= 0)&&(mode_y_index < source_mode.j)){

					permittivity = psi.ep[ep_index];
					permeability = psi.mu[ep_index];

					mode_index = 6*(mode_x_index + source_mode.i*mode_y_index) + m;

					if(m < 3){

						psi.sources[source_index] = amplitude*permittivity*(source_mode.fields.m[mode_index].re*cos(phase) + source_mode.fields.m[mode_index].im*sin(phase));
					}
					else{

						psi.sources[source_index] = -M_sign*amplitude*permeability*(source_mode.fields.m[mode_index].re*cos(phase) + source_mode.fields.m[mode_index].im*sin(phase));
					}
				}
			}
		}
	}
}

void timefield3d_mode_source_adjoint(timefield3d psi, int center_x, int center_y, int center_z, mode source_mode, double phase, int reverse){

	int i, j, k, m, source_index, mode_x_index, mode_y_index, mode_index, M_sign;
	double x_distance, y_distance;
	double slope, inv_slope;
	double norm, inv_norm;
	complex mode_Ex, mode_Ey, mode_Hx, mode_Hy;

	M_sign = -1 + 2*reverse;

	//power-normalize
	mode_power_normalize(source_mode);

	//set sources
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			k = center_z;

			x_distance = (i - center_x)*(psi.disc/source_mode.disc);
			y_distance = (j - center_y)*(psi.disc/source_mode.disc);

			source_index = 6*(i + psi.i*(j + k*(psi.j)));
			mode_x_index = (int)x_distance + source_mode.i/2;
			mode_y_index = (int)y_distance + source_mode.j/2;

			if((mode_x_index >= 0)&&(mode_x_index < source_mode.i)&&(mode_y_index >= 0)&&(mode_y_index < source_mode.j)){

				mode_index = 6*(mode_x_index + source_mode.i*mode_y_index);

				psi.sources[source_index] = 1*(source_mode.fields.m[mode_index+4].re*cos(phase) + source_mode.fields.m[mode_index+4].im*sin(phase));
				psi.sources[source_index+1] = -1*(source_mode.fields.m[mode_index+3].re*cos(phase) + source_mode.fields.m[mode_index+3].im*sin(phase));
				psi.sources[source_index+2] = 0;

				psi.sources[source_index+3] = -1*M_sign*(source_mode.fields.m[mode_index+1].re*cos(phase) + source_mode.fields.m[mode_index+1].im*sin(phase));
				psi.sources[source_index+4] = 1*M_sign*(source_mode.fields.m[mode_index].re*cos(phase) + source_mode.fields.m[mode_index].im*sin(phase));
				psi.sources[source_index+5] = 0;
			}
		}
	}
}

void timefield3d_TE_gaussian_source(timefield3d psi, int center_x, int center_y, int center_z, double width, double omega, double phase, double amplitude){
	
	int i, j, k, m, mode_size, source_index, mode_index, ep_index, sig_index;
	mode source_mode;
	complex ep, mu;
	double mag_ep, mag_mu;
	double source_re, source_im;

	//calculate source mode (use complex ep and mu and center of mode)
	ep_index = 3*(center_x + psi.i*(center_y + psi.j*center_z));
	sig_index = 6*(center_x + psi.i*(center_y + psi.j*center_z));
	ep = complex_from_doubles(psi.ep[ep_index], psi.sig[sig_index]/omega);
	mu = complex_from_doubles(psi.mu[ep_index], psi.sig[sig_index + 3]/omega);
	
	mode_size = (int)(3.0*width/psi.disc);
	source_mode = mode_gaussian(mode_size, mode_size, psi.disc, ep, mu, omega, (width/(2.35482*psi.disc)), 0);
	//mode_plot(source_mode);

	//set sources
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			k = center_z;
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			sig_index = 6*(i + psi.i*(j + psi.j*k));
			source_index = 6*(i + psi.i*(j + psi.j*k));
			mode_index = 6*((i - center_x + mode_size/2) + source_mode.i*(j - center_y + mode_size/2));

			if((mode_index >= 0)&&(mode_index < 6*source_mode.i*source_mode.j)){
				for(m = 0; m < 3; m++){
					
					ep = complex_from_doubles(psi.ep[ep_index + m], psi.sig[sig_index + m]/omega);
					mu = complex_from_doubles(psi.mu[ep_index + m], psi.sig[sig_index + m + 3]/omega);
					
					//E
					source_re = amplitude*(ep.re*source_mode.fields.m[mode_index + m].re - ep.im*source_mode.fields.m[mode_index + m].im);
					source_im = amplitude*(ep.re*source_mode.fields.m[mode_index + m].im + ep.im*source_mode.fields.m[mode_index + m].re);
					psi.sources[source_index + m] = source_re*cos(phase) + source_im*sin(phase);
					
					//H
					source_re = amplitude*(mu.re*source_mode.fields.m[mode_index + m + 3].re - mu.im*source_mode.fields.m[mode_index + m + 3].im);
					source_im = amplitude*(mu.re*source_mode.fields.m[mode_index + m + 3].im + mu.im*source_mode.fields.m[mode_index + m + 3].re);
					psi.sources[source_index + m + 3] = source_re*cos(phase) + source_im*sin(phase);
				
				}
			}
		}
	}
	
	mode_free(source_mode);
}

void timefield3d_TM_gaussian_source(timefield3d psi, int center_x, int center_y, int center_z, double width, double omega, double phase, double amplitude){
	
	int i, j, k, m, mode_size, source_index, mode_index, ep_index, sig_index;
	mode source_mode;
	complex ep, mu;
	double mag_ep, mag_mu;
	double source_re, source_im;

	//calculate source mode (use complex ep and mu and center of mode)
	ep_index = 3*(center_x + psi.i*(center_y + psi.j*center_z));
	sig_index = 6*(center_x + psi.i*(center_y + psi.j*center_z));
	ep = complex_from_doubles(psi.ep[ep_index], psi.sig[sig_index]/omega);
	mu = complex_from_doubles(psi.mu[ep_index], psi.sig[sig_index + 3]/omega);
	
	mode_size = (int)(3.0*width/psi.disc);
	source_mode = mode_gaussian(mode_size, mode_size, psi.disc, ep, mu, omega, (width/(2.35482*psi.disc)), 1);
	//mode_plot(source_mode);

	//set sources
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			k = center_z;
			
			ep_index = 3*(i + psi.i*(j + psi.j*k));
			sig_index = 6*(i + psi.i*(j + psi.j*k));
			source_index = 6*(i + psi.i*(j + psi.j*k));
			mode_index = 6*((i - center_x + mode_size/2) + source_mode.i*(j - center_y + mode_size/2));

			if((mode_index >= 0)&&(mode_index < 6*source_mode.i*source_mode.j)){
				for(m = 0; m < 3; m++){
					
					ep = complex_from_doubles(psi.ep[ep_index + m], psi.sig[sig_index + m]/omega);
					mu = complex_from_doubles(psi.mu[ep_index + m], psi.sig[sig_index + m + 3]/omega);
					
					//E
					source_re = amplitude*(ep.re*source_mode.fields.m[mode_index + m].re - ep.im*source_mode.fields.m[mode_index + m].im);
					source_im = amplitude*(ep.re*source_mode.fields.m[mode_index + m].im + ep.im*source_mode.fields.m[mode_index + m].re);
					psi.sources[source_index + m] = source_re*cos(phase) + source_im*sin(phase);
					
					//H
					source_re = amplitude*(mu.re*source_mode.fields.m[mode_index + m + 3].re - mu.im*source_mode.fields.m[mode_index + m + 3].im);
					source_im = amplitude*(mu.re*source_mode.fields.m[mode_index + m + 3].im + mu.im*source_mode.fields.m[mode_index + m + 3].re);
					psi.sources[source_index + m + 3] = source_re*cos(phase) + source_im*sin(phase);
				
				}
			}
		}
	}
	
	mode_free(source_mode);
}

complex timefield3d_mode_overlap(timefield3d psi, int center_x, int center_y, int center_z, mode overlap_mode, double phase, int reverse, complex *cumulative_overlap){

	int i, j, k, m, ep_index, field_index;
	int mode_x_index, mode_y_index, mode_index;
	double x_distance, y_distance, reduced_phase;
	double slope, inv_slope;
	double norm, inv_norm;
	complex mode_Ex, mode_Ey, mode_Hx, mode_Hy;
	double Ex, Ey, Hx, Hy;
	complex product1, product2, product3, product4;
	complex inner_product, instant_overlap;

	mode_power_normalize(overlap_mode);

	//now take overlap with timefield3d structure (1D inner product)
	k = center_z;
	instant_overlap = complex_zero();
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			field_index = 6*(i + psi.i*(j + k*(psi.j)));

			x_distance = (i - center_x)*(psi.disc/overlap_mode.disc);
			y_distance = (j - center_y)*(psi.disc/overlap_mode.disc);
			mode_x_index = (int)x_distance + overlap_mode.i/2;
			mode_y_index = (int)y_distance + overlap_mode.j/2;
			mode_index = 6*(mode_x_index + overlap_mode.i*mode_y_index);

			if((mode_x_index >= 0)&&(mode_x_index < overlap_mode.i)&&(mode_y_index >= 0)&&(mode_y_index < overlap_mode.j)){

				//pick out field elements
				mode_Ex = overlap_mode.fields.m[mode_index];
				mode_Ey = overlap_mode.fields.m[mode_index + 1];
				mode_Hx = overlap_mode.fields.m[mode_index + 3];
				mode_Hy = overlap_mode.fields.m[mode_index + 4];

				if(reverse == 1){

					mode_Hx = complex_sub(complex_zero(), complex_conj(mode_Hx));
					mode_Hy = complex_sub(complex_zero(), complex_conj(mode_Hy));
				}

				Ex = psi.fields[field_index];
				Ey = psi.fields[field_index + 1];
				Hx = psi.fields[field_index + 3];
				Hy = psi.fields[field_index + 4];

				product1 = complex_smul(Hx, mode_Ey);
				product2 = complex_sub(complex_zero(), complex_smul(Hy, mode_Ex));
				product3 = complex_smul(Ey, mode_Hx);
				product4 = complex_sub(complex_zero(), complex_smul(Ex, mode_Hy));

				instant_overlap = complex_add(instant_overlap, complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			}
		}
	}
	instant_overlap = complex_smul((-0.5*psi.disc*psi.disc), instant_overlap);

	//add to cumulative overlap with phase
	reduced_phase = phase - ((int)(phase/(2*PI)))*2*PI;
	*cumulative_overlap = complex_add(*cumulative_overlap, complex_mul(complex_exp(complex_from_doubles(0,reduced_phase)), instant_overlap));

	return(instant_overlap);
}

void timefield3d_file_write_permittivity_y_slice(timefield3d psi, char *filename, int left, int right, int start, int end, int y_slice, double ep_min, double ep_max){

	int i, j;
	int ep_index, map_index;
	int *y_size, *z_size, *permittivity_map;
	double *discretization;
	FILE *output_file;

	y_size = malloc(sizeof(int));
	z_size = malloc(sizeof(int));
	discretization = malloc(sizeof(double));

	*y_size = right-left;
	*z_size = end-start;
	*discretization = psi.disc;

	permittivity_map = malloc((*y_size)*(*z_size)*sizeof(double));
	for(i = 0; i < *y_size; i++){
		for(j = 0; j < *z_size; j++){

			ep_index = 3*((i+left) + psi.i*(y_slice + psi.j*(j+start)));
			map_index = i + j*(*y_size);

			if(psi.ep[ep_index] > (ep_min+ep_max)/2){

				permittivity_map[map_index] = 1;
			}
			else{

				permittivity_map[map_index] = 0;
			}

		}
	}

	output_file = fopen(filename, "wb");

	fwrite(y_size, sizeof(int), 1, output_file);
	fwrite(z_size, sizeof(int), 1, output_file);
	fwrite(discretization, sizeof(double), 1, output_file);
	fwrite(permittivity_map, sizeof(int), (*y_size)*(*z_size), output_file);

	fclose(output_file);
}
