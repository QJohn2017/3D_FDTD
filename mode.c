#include "mode.h"

//general functions
mode mode_new(int i, int j, double disc){
	
	mode new;
	int n, index;

	new.i = i;
	new.j = j;
	new.disc = disc;

	n = i*j;
	
	new.pml = matrix_new(2*n, 1);
	new.ep = matrix_new(3*n, 1);
	new.mu = matrix_new(3*n, 1);
	new.fields = matrix_new(6*n, 1);
	
	for(index = 0; index < 2*n; index++){
		
		new.pml.m[index] = complex_from_doubles(1, 0);
	}
	for(index = 0; index < 3*n; index++){
		
		new.ep.m[index] = complex_zero();
		new.mu.m[index] = complex_zero();
	}
	for(index = 0; index < 6*n; index++){
		
		new.fields.m[index] = complex_zero();
	}

	return(new);
}

void mode_free(mode psi){
	
	matrix_free(psi.pml);
	matrix_free(psi.ep);
	matrix_free(psi.mu);
	matrix_free(psi.fields);
}

mode mode_copy(mode psi){
	
	int i, field_index, ep_index, pml_index;
	mode copy;
	
	copy = mode_new(psi.i, psi.j, psi.disc);
	
	//copy contents into new matrices
	for(i = 0; i < psi.i*psi.j; i++){
		
		field_index = 6*i;
		ep_index = 3*i;
		pml_index = 2*i;

		copy.pml.m[pml_index] = psi.pml.m[pml_index];
		copy.pml.m[pml_index + 1] = psi.pml.m[pml_index + 1];
		copy.ep.m[ep_index] = psi.ep.m[ep_index];
		copy.ep.m[ep_index + 1] = psi.ep.m[ep_index + 1];
		copy.ep.m[ep_index + 2] = psi.ep.m[ep_index + 2];
		copy.mu.m[ep_index] = psi.mu.m[ep_index];
		copy.mu.m[ep_index + 1] = psi.mu.m[ep_index + 1];
		copy.mu.m[ep_index + 2] = psi.mu.m[ep_index + 2];
		copy.fields.m[field_index] = psi.fields.m[field_index];
		copy.fields.m[field_index + 1] = psi.fields.m[field_index + 1];
		copy.fields.m[field_index + 2] = psi.fields.m[field_index + 2];
		copy.fields.m[field_index + 3] = psi.fields.m[field_index + 3];
		copy.fields.m[field_index + 4] = psi.fields.m[field_index + 4];
		copy.fields.m[field_index + 5] = psi.fields.m[field_index + 5];
	}
	
	return(copy);
}

void mode_smul(complex scalar, mode psi){

	int i, field_index;

	//copy contents into new matrices
	for(i = 0; i < psi.i*psi.j; i++){

		field_index = 6*i;

		psi.fields.m[field_index] = complex_mul(scalar, psi.fields.m[field_index]);
		psi.fields.m[field_index + 1] = complex_mul(scalar, psi.fields.m[field_index+1]);
		psi.fields.m[field_index + 2] = complex_mul(scalar, psi.fields.m[field_index+2]);
		psi.fields.m[field_index + 3] = complex_mul(scalar, psi.fields.m[field_index+3]);
		psi.fields.m[field_index + 4] = complex_mul(scalar, psi.fields.m[field_index+4]);
		psi.fields.m[field_index + 5] = complex_mul(scalar, psi.fields.m[field_index+5]);
	}
}

mode mode_linear_combination(mode psi, mode phi, complex scalar){
	
	int i, field_index, ep_index, pml_index;
	mode linear_combination;
	
	if((psi.i != phi.i)||(psi.j != phi.j)){
		
		printf("Mode dimensions do not match for linear combination!\n");
		abort();
	}
	
	linear_combination = mode_new(psi.i, psi.j, psi.disc);
	
	//copy contents into new matrices
	for(i = 0; i < psi.i*psi.j; i++){
		
		field_index = 6*i;
		ep_index = 3*i;
		pml_index = 2*i;

		linear_combination.pml.m[pml_index] = psi.pml.m[pml_index];
		linear_combination.pml.m[pml_index + 1] = psi.pml.m[pml_index + 1];
		linear_combination.ep.m[ep_index] = psi.ep.m[ep_index];
		linear_combination.ep.m[ep_index + 1] = psi.ep.m[ep_index + 1];
		linear_combination.ep.m[ep_index + 2] = psi.ep.m[ep_index + 2];
		linear_combination.mu.m[ep_index] = psi.mu.m[ep_index];
		linear_combination.mu.m[ep_index + 1] = psi.mu.m[ep_index + 1];
		linear_combination.mu.m[ep_index + 2] = psi.mu.m[ep_index + 2];
		linear_combination.fields.m[field_index] = complex_add(psi.fields.m[field_index], complex_mul(scalar, phi.fields.m[field_index]));
		linear_combination.fields.m[field_index + 1] = complex_add(psi.fields.m[field_index + 1], complex_mul(scalar, phi.fields.m[field_index + 1]));
		linear_combination.fields.m[field_index + 2] = complex_add(psi.fields.m[field_index + 2], complex_mul(scalar, phi.fields.m[field_index + 2]));
		linear_combination.fields.m[field_index + 3] = complex_add(psi.fields.m[field_index + 3], complex_mul(scalar, phi.fields.m[field_index + 3]));
		linear_combination.fields.m[field_index + 4] = complex_add(psi.fields.m[field_index + 4], complex_mul(scalar, phi.fields.m[field_index + 4]));
		linear_combination.fields.m[field_index + 5] = complex_add(psi.fields.m[field_index + 5], complex_mul(scalar, phi.fields.m[field_index + 5]));
	}
	
	return(linear_combination);
}

mode mode_linear_combination_3(mode psi, mode phi, mode alpha, complex scalar_phi, complex scalar_alpha){
	
	int i, field_index, ep_index, pml_index;
	mode linear_combination;
	
	if((psi.i != phi.i)||(psi.j != phi.j)){
		
		printf("Mode dimensions do not match for linear combination!\n");
		abort();
	}
	
	linear_combination = mode_new(psi.i, psi.j, psi.disc);
	
	//copy contents into new matrices
	for(i = 0; i < psi.i*psi.j; i++){
		
		field_index = 6*i;
		ep_index = 3*i;
		pml_index = 2*i;

		linear_combination.pml.m[pml_index] = psi.pml.m[pml_index];
		linear_combination.pml.m[pml_index + 1] = psi.pml.m[pml_index + 1];
		linear_combination.ep.m[ep_index] = psi.ep.m[ep_index];
		linear_combination.ep.m[ep_index + 1] = psi.ep.m[ep_index + 1];
		linear_combination.ep.m[ep_index + 2] = psi.ep.m[ep_index + 2];
		linear_combination.mu.m[ep_index] = psi.mu.m[ep_index];
		linear_combination.mu.m[ep_index + 1] = psi.mu.m[ep_index + 1];
		linear_combination.mu.m[ep_index + 2] = psi.mu.m[ep_index + 2];
		linear_combination.fields.m[field_index] = complex_add(complex_add(psi.fields.m[field_index], complex_mul(scalar_phi, phi.fields.m[field_index])), complex_mul(scalar_alpha, alpha.fields.m[field_index]));
		linear_combination.fields.m[field_index + 1] = complex_add(complex_add(psi.fields.m[field_index + 1], complex_mul(scalar_phi, phi.fields.m[field_index + 1])), complex_mul(scalar_alpha, alpha.fields.m[field_index + 1]));
		linear_combination.fields.m[field_index + 2] = complex_add(complex_add(psi.fields.m[field_index + 2], complex_mul(scalar_phi, phi.fields.m[field_index + 2])), complex_mul(scalar_alpha, alpha.fields.m[field_index + 2]));
		linear_combination.fields.m[field_index + 3] = complex_add(complex_add(psi.fields.m[field_index + 3], complex_mul(scalar_phi, phi.fields.m[field_index + 3])), complex_mul(scalar_alpha, alpha.fields.m[field_index + 3]));
		linear_combination.fields.m[field_index + 4] = complex_add(complex_add(psi.fields.m[field_index + 4], complex_mul(scalar_phi, phi.fields.m[field_index + 4])), complex_mul(scalar_alpha, alpha.fields.m[field_index + 4]));
		linear_combination.fields.m[field_index + 5] = complex_add(complex_add(psi.fields.m[field_index + 5], complex_mul(scalar_phi, phi.fields.m[field_index + 5])), complex_mul(scalar_alpha, alpha.fields.m[field_index + 5]));
	}
	
	return(linear_combination);
}

mode mode_gaussian(int i, int j, double disc, complex ep, complex mu, double omega, double sigma, int polarization){
	
	mode gaussian;
	int n, index, field_index;
	double x, y, x_center, y_center, del_x, del_y, r_squared;

	gaussian.i = i;
	gaussian.j = j;
	gaussian.disc = disc;

	n = i*j;
	x_center = (double)(i)/2;
	y_center = (double)(j)/2;
	
	gaussian.pml = matrix_new(2*n, 1);
	gaussian.ep = matrix_new(3*n, 1);
	gaussian.mu = matrix_new(3*n, 1);
	gaussian.fields = matrix_new(6*n, 1);
	
	for(index = 0; index < 2*n; index++){
		
		gaussian.pml.m[index] = complex_from_doubles(1, 0);
	}
	for(index = 0; index < 3*n; index++){
		
		gaussian.ep.m[index] = ep;
		gaussian.mu.m[index] = mu;
	}
	for(index = 0; index < n; index++){
		
		x = (double)(index%i);
		y = (double)(index)/i;
		del_x = x - x_center;
		del_y = y - y_center;
		r_squared = del_x*del_x + del_y*del_y;
		
		field_index = 6*index;
		
		gaussian.fields.m[field_index] = complex_zero();
		gaussian.fields.m[field_index + 1] = complex_zero();
		gaussian.fields.m[field_index + 2] = complex_zero();
		gaussian.fields.m[field_index + 3] = complex_zero();
		gaussian.fields.m[field_index + 4] = complex_zero();
		gaussian.fields.m[field_index + 5] = complex_zero();

		if(polarization){

			gaussian.fields.m[field_index + 3] = complex_from_doubles(exp(-r_squared/(2*sigma*sigma)),0);	
		}
		else{
			
			gaussian.fields.m[field_index + 4] = complex_from_doubles(exp(-r_squared/(2*sigma*sigma)),0);
		}
	}

	mode_apply_H_curl(gaussian, complex_from_doubles(omega, 0), omega);
	mode_apply_H_div(gaussian, complex_from_doubles(omega, 0));
	mode_power_normalize(gaussian);
		
	return(gaussian);
}

void mode_chop_bottom_half(mode psi){
	
	int i, j, field_index;
	
	//zero out bottom half
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j/2; j++){
			
			field_index = 6*(i + j*psi.i);

			psi.fields.m[field_index] = complex_zero();
			psi.fields.m[field_index + 1] = complex_zero();
			psi.fields.m[field_index + 2] = complex_zero();
			psi.fields.m[field_index + 3] = complex_zero();
			psi.fields.m[field_index + 4] = complex_zero();
			psi.fields.m[field_index + 5] = complex_zero();
		}
	}
}

void mode_copy_fields(mode psi, mode phi){
	
	int i;
	
	//ensure that psi and phi are the same size
	if((psi.i != phi.i)||(psi.j != phi.j)||(psi.disc != phi.disc)){
		
		printf("Mode dimensions do not match for addition!\n");
		abort();
	}
	
	for(i = 0; i < 6*psi.i*psi.j; i++){
		
		psi.fields.m[i] = phi.fields.m[i];
	}
}

void mode_plot(mode psi){
	

	int i, j, index;
	double *X, *Y;
	FILE *gnuplot;
	FILE *ep, *mu, *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;

	gnuplot = popen("/usr/bin/gnuplot5-qt -persist", "w");
	//set arrays corresponding to the X and Y axis values
	X = malloc(psi.i*sizeof(double));
	Y = malloc(psi.j*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.j; j++){
		
		Y[j] = (j + 0.5 - psi.j/2.0)*psi.disc;
	}
	
	//print matrices to files
	ep = fopen("ep.txt", "w");
	fprintf(ep, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(ep, "%lf\t", X[i]);
	}
	fprintf(ep, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(ep, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(ep, "%lf\t", psi.ep.m[3*index].re);
		}	
		fprintf(ep, "\n");
	}
	fclose(ep);
	
	//print matrices to files
	mu = fopen("mu.txt", "w");
	fprintf(mu, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(mu, "%lf\t", X[i]);
	}
	fprintf(mu, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(mu, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(mu, "%lf\t", psi.mu.m[3*index].re);
		}	
		fprintf(mu, "\n");
	}
	fclose(mu);
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Ex, "%lf\t", complex_mag(psi.fields.m[6*index]));
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);
	
	Ey = fopen("ey.txt", "w");
	fprintf(Ey, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ey, "%lf\t", X[i]);
	}
	fprintf(Ey, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ey, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Ey, "%lf\t", complex_mag(psi.fields.m[6*index + 1]));
		}	
		fprintf(Ey, "\n");
	}
	fclose(Ey);

	Ez = fopen("ez.txt", "w");
	fprintf(Ez, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ez, "%lf\t", X[i]);
	}
	fprintf(Ez, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ez, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Ez, "%lf\t", complex_mag(psi.fields.m[6*index + 2]));
		}	
		fprintf(Ez, "\n");
	}
	fclose(Ez);
	
	Hx = fopen("hx.txt", "w");
	fprintf(Hx, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Hx, "%lf\t", X[i]);
	}
	fprintf(Hx, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Hx, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Hx, "%lf\t", complex_mag(psi.fields.m[6*index + 3]));
		}	
		fprintf(Hx, "\n");
	}
	fclose(Hx);
	
	Hy = fopen("hy.txt", "w");
	fprintf(Hy, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Hy, "%lf\t", X[i]);
	}
	fprintf(Hy, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Hy, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Hy, "%lf\t", complex_mag(psi.fields.m[6*index + 4]));
		}	
		fprintf(Hy, "\n");
	}
	fclose(Hy);

	Hz = fopen("hz.txt", "w");
	fprintf(Hz, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Hz, "%lf\t", X[i]);
	}
	fprintf(Hz, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Hz, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Hz, "%lf\t", complex_mag(psi.fields.m[6*index + 5]));
		}	
		fprintf(Hz, "\n");
	}
	fclose(Hz);
	
	//plot the index and field distributions
	fprintf(gnuplot, "set term x11 size 900,700 background \"gray\";\n");
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set multiplot; set size 0.33,0.33;\n");
	
	fprintf(gnuplot, "set origin 0.167,0.66; set title \"Permittivity\";\n");
	fprintf(gnuplot, "set cbrange[0:11];\n");
	fprintf(gnuplot, "plot \"ep.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.5,0.66; set title \"Permeability\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"mu.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0,0.33; set title \"Ex\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.33,0.33; set title \"Ey\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"ey.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.66,0.33; set title \"Ez\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"ez.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0,0; set title \"Hx\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"hx.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.33,0; set title \"Hy\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"hy.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.66,0; set title \"Hz\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"hz.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ep.txt");	
	remove("ex.txt");
	remove("ey.txt");
	remove("ez.txt");
	remove("hx.txt");
	remove("hy.txt");
	remove("hz.txt");
}

void mode_plot_real(mode psi){
	

	int i, j, index;
	double *X, *Y;
	FILE *gnuplot;
	FILE *ep, *mu, *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;

	gnuplot = popen("/usr/local/bin/gnuplot -persist", "w");
	//set arrays corresponding to the X and Y axis values
	X = malloc(psi.i*sizeof(double));
	Y = malloc(psi.j*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.j; j++){
		
		Y[j] = (j + 0.5 - psi.j/2.0)*psi.disc;
	}
	
	//print matrices to files
	ep = fopen("ep.txt", "w");
	fprintf(ep, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(ep, "%lf\t", X[i]);
	}
	fprintf(ep, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(ep, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(ep, "%lf\t", psi.ep.m[3*index].re);
		}	
		fprintf(ep, "\n");
	}
	fclose(ep);
	
	//print matrices to files
	mu = fopen("mu.txt", "w");
	fprintf(mu, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(mu, "%lf\t", X[i]);
	}
	fprintf(mu, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(mu, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(mu, "%lf\t", psi.mu.m[3*index].re);
		}	
		fprintf(mu, "\n");
	}
	fclose(mu);
	
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Ex, "%lf\t", psi.fields.m[6*index].re);
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);
	
	Ey = fopen("ey.txt", "w");
	fprintf(Ey, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ey, "%lf\t", X[i]);
	}
	fprintf(Ey, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ey, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Ey, "%lf\t", psi.fields.m[6*index + 1].re);
		}	
		fprintf(Ey, "\n");
	}
	fclose(Ey);

	Ez = fopen("ez.txt", "w");
	fprintf(Ez, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ez, "%lf\t", X[i]);
	}
	fprintf(Ez, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ez, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Ez, "%lf\t", psi.fields.m[6*index + 2].re);
		}	
		fprintf(Ez, "\n");
	}
	fclose(Ez);
	
	Hx = fopen("hx.txt", "w");
	fprintf(Hx, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Hx, "%lf\t", X[i]);
	}
	fprintf(Hx, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Hx, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Hx, "%lf\t", psi.fields.m[6*index + 3].re);
		}	
		fprintf(Hx, "\n");
	}
	fclose(Hx);
	
	Hy = fopen("hy.txt", "w");
	fprintf(Hy, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Hy, "%lf\t", X[i]);
	}
	fprintf(Hy, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Hy, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Hy, "%lf\t", psi.fields.m[6*index + 4].re);
		}	
		fprintf(Hy, "\n");
	}
	fclose(Hy);

	Hz = fopen("hz.txt", "w");
	fprintf(Hz, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Hz, "%lf\t", X[i]);
	}
	fprintf(Hz, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Hz, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Hz, "%lf\t", psi.fields.m[6*index + 5].re);
		}	
		fprintf(Hz, "\n");
	}
	fclose(Hz);
	
	//plot the index and field distributions
	fprintf(gnuplot, "set term x11 size 900,700 background \"gray\";\n");
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set multiplot; set size 0.33,0.33;\n");
	
	fprintf(gnuplot, "set origin 0.167,0.66; set title \"Permittivity\";\n");
	fprintf(gnuplot, "set cbrange[0:11];\n");
	fprintf(gnuplot, "plot \"ep.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.5,0.66; set title \"Permeability\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"mu.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0,0.33; set title \"Ex\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.33,0.33; set title \"Ey\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"ey.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.66,0.33; set title \"Ez\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"ez.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0,0; set title \"Hx\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"hx.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.33,0; set title \"Hy\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"hy.txt\" nonuniform matrix with image notitle\n");
	
	fprintf(gnuplot, "set origin 0.66,0; set title \"Hz\";\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"hz.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ep.txt");	
	remove("ex.txt");
	remove("ey.txt");
	remove("ez.txt");
	remove("hx.txt");
	remove("hy.txt");
	remove("hz.txt");
}

void mode_plot_ex(mode psi){
	
	int i, j, index;
	double *X, *Y;
	FILE *gnuplot;
	FILE *ep, *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;

	gnuplot = popen("/usr/local/bin/gnuplot -persist", "w");
	//set arrays corresponding to the X and Y axis values
	X = malloc(psi.i*sizeof(double));
	Y = malloc(psi.j*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.j; j++){
		
		Y[j] = (j + 0.5 - psi.j/2.0)*psi.disc;
	}
	
	//print matrices to files
	Ex = fopen("ex.txt", "w");
	fprintf(Ex, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Ex, "%lf\t", X[i]);
	}
	fprintf(Ex, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Ex, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Ex, "%lf\t", complex_mag(psi.fields.m[6*index]));
		}	
		fprintf(Ex, "\n");
	}
	fclose(Ex);
	

	
	//plot the index and field distributions
	fprintf(gnuplot, "set term x11 size 900,700 background \"gray\";\n");
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"ex.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("ex.txt");
}

void mode_plot_hy(mode psi){
	
	int i, j, index;
	double *X, *Y;
	FILE *gnuplot;
	FILE *ep, *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;

	gnuplot = popen("/usr/local/bin/gnuplot -persist", "w");
	//set arrays corresponding to the X and Y axis values
	X = malloc(psi.i*sizeof(double));
	Y = malloc(psi.j*sizeof(double));
	
	for(i = 0; i < psi.i; i++){
		
		X[i] = (i + 0.5 - psi.i/2.0)*psi.disc;
	}
	for(j = 0; j < psi.j; j++){
		
		Y[j] = (j + 0.5 - psi.j/2.0)*psi.disc;
	}
	
	//print matrices to files
	Hy = fopen("hy.txt", "w");
	fprintf(Hy, "%lf\t", (double)psi.j);
	for(i = 0; i < psi.i; i++){
		
		fprintf(Hy, "%lf\t", X[i]);
	}
	fprintf(Hy, "\n");
	
	for(j = 0; j < psi.j; j++){
		
		fprintf(Hy, "%lf\t", Y[j]);
		for(i = 0; i < psi.i; i++){
			
			index = i + j*psi.i;
			fprintf(Hy, "%lf\t", complex_mag(psi.fields.m[6*index + 4]));
		}	
		fprintf(Hy, "\n");
	}
	fclose(Hy);
	

	
	//plot the index and field distributions
	fprintf(gnuplot, "set term x11 size 900,700 background \"gray\";\n");
	fprintf(gnuplot, "set palette defined (-8 \"black\", 1 '#00008f', 8 '#0000ff', 22 '#00ffff', 42 '#ffff00', 56 '#ff0000', 64 '#800000');\n");
	fprintf(gnuplot, "set autoscale;\n");
	fprintf(gnuplot, "plot \"hy.txt\" nonuniform matrix with image notitle\n");
	
	pclose(gnuplot);
	remove("hy.txt");
}

void mode_file_print(FILE *mode_data, mode psi, int dimension){
	
	int i, j, field_index, ep_index;
	
	fprintf(mode_data, "%d\t%d\t%f\t%d\n", psi.i, psi.j, psi.disc, dimension);
	
	if(dimension == 2){
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				field_index = 6*(i + j*psi.i);
				ep_index = 3*(i + j*psi.i);

				fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index].re, psi.ep.m[ep_index].im);
				fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index + 1].re, psi.ep.m[ep_index + 1].im);
				fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index + 2].re, psi.ep.m[ep_index + 2].im);
				fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index].re, psi.mu.m[ep_index].im);
				fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index + 1].re, psi.mu.m[ep_index + 1].im);
				fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index + 2].re, psi.mu.m[ep_index + 2].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index].re, psi.fields.m[field_index].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 1].re, psi.fields.m[field_index + 1].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 2].re, psi.fields.m[field_index + 2].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 3].re, psi.fields.m[field_index + 3].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 4].re, psi.fields.m[field_index + 4].im);
				fprintf(mode_data, "%f\t%f\n", psi.fields.m[field_index + 5].re, psi.fields.m[field_index + 5].im);
			}
		}
	}
	else if(dimension == 1){
		i = psi.i/2;
		for(j = 0; j < psi.j; j++){

			field_index = 6*(i + j*psi.i);

			field_index = 6*(i + j*psi.i);
			ep_index = 3*(i + j*psi.i);

			fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index].re, psi.ep.m[ep_index].im);
			fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index + 1].re, psi.ep.m[ep_index + 1].im);
			fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index + 2].re, psi.ep.m[ep_index + 2].im);
			fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index].re, psi.mu.m[ep_index].im);
			fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index + 1].re, psi.mu.m[ep_index + 1].im);
			fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index + 2].re, psi.mu.m[ep_index + 2].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index].re, psi.fields.m[field_index].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 1].re, psi.fields.m[field_index + 1].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 2].re, psi.fields.m[field_index + 2].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 3].re, psi.fields.m[field_index + 3].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 4].re, psi.fields.m[field_index + 4].im);
			fprintf(mode_data, "%f\t%f\n", psi.fields.m[field_index + 5].re, psi.fields.m[field_index + 5].im);
		}
	}
}

void mode_angular_file_print(FILE *mode_data, mode psi, int dimension){
	
	int i, j, field_index, ep_index;
	
	fprintf(mode_data, "%d\t%d\t%f\t%d\n", psi.i, psi.j, psi.disc, dimension);
	
	if(dimension == 2){
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				field_index = 6*(i + j*psi.i);
				ep_index = 3*(i + j*psi.i);

				fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index].re, psi.ep.m[ep_index].im);
				fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index + 1].re, psi.ep.m[ep_index + 1].im);
				fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index + 2].re, psi.ep.m[ep_index + 2].im);
				fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index].re, psi.mu.m[ep_index].im);
				fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index + 1].re, psi.mu.m[ep_index + 1].im);
				fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index + 2].re, psi.mu.m[ep_index + 2].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index].re, psi.fields.m[field_index].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 1].re, psi.fields.m[field_index + 1].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 2].re, psi.fields.m[field_index + 2].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 3].re, psi.fields.m[field_index + 3].im);
				fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 4].re, psi.fields.m[field_index + 4].im);
				fprintf(mode_data, "%f\t%f\n", psi.fields.m[field_index + 5].re, psi.fields.m[field_index + 5].im);
			}
		}
	}
	else if(dimension == 1){
		j = psi.j/2;
		for(i = 0; i < psi.i; i++){

			field_index = 6*(i + j*psi.i);

			field_index = 6*(i + j*psi.i);
			ep_index = 3*(i + j*psi.i);

			fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index].re, psi.ep.m[ep_index].im);
			fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index + 1].re, psi.ep.m[ep_index + 1].im);
			fprintf(mode_data, "%f\t%f\t", psi.ep.m[ep_index + 2].re, psi.ep.m[ep_index + 2].im);
			fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index].re, psi.mu.m[ep_index].im);
			fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index + 1].re, psi.mu.m[ep_index + 1].im);
			fprintf(mode_data, "%f\t%f\t", psi.mu.m[ep_index + 2].re, psi.mu.m[ep_index + 2].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index].re, psi.fields.m[field_index].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 1].re, psi.fields.m[field_index + 1].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 2].re, psi.fields.m[field_index + 2].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 3].re, psi.fields.m[field_index + 3].im);
			fprintf(mode_data, "%f\t%f\t", psi.fields.m[field_index + 4].re, psi.fields.m[field_index + 4].im);
			fprintf(mode_data, "%f\t%f\n", psi.fields.m[field_index + 5].re, psi.fields.m[field_index + 5].im);
		}
	}
}

mode mode_file_load(FILE *mode_data){
	
	int i, j, k, field_index, ep_index;
	int dimension, x, y;
	double disc;
	double real_val, imag_val;
	mode file_mode;
	
	fscanf(mode_data, "%d", &x);
	fscanf(mode_data, "%d", &y);
	fscanf(mode_data, "%lf", &disc);
	fscanf(mode_data, "%d", &dimension);
	
	file_mode = mode_new(x, y, disc);
	
	if(dimension == 1){
		
		for(j = 0; j < file_mode.j; j++){

			for(k = 0; k < 3; k++){

				fscanf(mode_data, "%lf", &real_val);
				fscanf(mode_data, "%lf", &imag_val); 

					for(i = 0; i < file_mode.i; i++){

						ep_index = 3*(i + j*file_mode.i) + k;
						file_mode.ep.m[ep_index] = complex_from_doubles(real_val, imag_val);
					}
			}
			for(k = 0; k < 3; k++){

				fscanf(mode_data, "%lf", &real_val);
				fscanf(mode_data, "%lf", &imag_val); 

					for(i = 0; i < file_mode.i; i++){

						ep_index = 3*(i + j*file_mode.i) + k;
						file_mode.mu.m[ep_index] = complex_from_doubles(real_val, imag_val);
					}
			}
			for(k = 0; k < 6; k++){

				fscanf(mode_data, "%lf", &real_val);
				fscanf(mode_data, "%lf", &imag_val); 

					for(i = 0; i < file_mode.i; i++){

						field_index = 6*(i + j*file_mode.i) + k;
						file_mode.fields.m[field_index] = complex_from_doubles(real_val, imag_val);
					}
			}

		}
	}
	if(dimension == 2){

		for(i = 0; i < file_mode.i; i++){		
			for(j = 0; j < file_mode.j; j++){
				for(k = 0; k < 3; k++){

					fscanf(mode_data, "%lf", &real_val);
					fscanf(mode_data, "%lf", &imag_val); 

					ep_index = 3*(i + j*file_mode.i) + k;
					file_mode.ep.m[ep_index] = complex_from_doubles(real_val, imag_val);
				}
				for(k = 0; k < 3; k++){

					fscanf(mode_data, "%lf", &real_val);
					fscanf(mode_data, "%lf", &imag_val); 

					ep_index = 3*(i + j*file_mode.i) + k;
					file_mode.mu.m[ep_index] = complex_from_doubles(real_val, imag_val);
				}
				for(k = 0; k < 6; k++){

					fscanf(mode_data, "%lf", &real_val);
					fscanf(mode_data, "%lf", &imag_val); 

					field_index = 6*(i + j*file_mode.i) + k;
					file_mode.fields.m[field_index] = complex_from_doubles(real_val, imag_val);
				}
			}

		}
	}
	
	return(file_mode);
}

mode mode_angular_file_load(FILE *mode_data){
	
	int i, j, k, field_index, ep_index;
	int dimension, x, y;
	double disc;
	double real_val, imag_val;
	mode file_mode;
	
	fscanf(mode_data, "%d", &x);
	fscanf(mode_data, "%d", &y);
	fscanf(mode_data, "%lf", &disc);
	fscanf(mode_data, "%d", &dimension);
	
	file_mode = mode_new(x, y, disc);
	
	if(dimension == 1){
		
		for(i = 0; i < file_mode.i; i++){

			for(k = 0; k < 3; k++){

				fscanf(mode_data, "%lf", &real_val);
				fscanf(mode_data, "%lf", &imag_val); 

					for(j = 0; j < file_mode.j; j++){

						ep_index = 3*(i + j*file_mode.i) + k;
						file_mode.ep.m[ep_index] = complex_from_doubles(real_val, imag_val);
					}
			}
			for(k = 0; k < 3; k++){

				fscanf(mode_data, "%lf", &real_val);
				fscanf(mode_data, "%lf", &imag_val); 

					for(j = 0; j < file_mode.j; j++){

						ep_index = 3*(i + j*file_mode.i) + k;
						file_mode.mu.m[ep_index] = complex_from_doubles(real_val, imag_val);
					}
			}
			for(k = 0; k < 6; k++){

				fscanf(mode_data, "%lf", &real_val);
				fscanf(mode_data, "%lf", &imag_val); 

					for(j = 0; j < file_mode.j; j++){

						field_index = 6*(i + j*file_mode.i) + k;
						file_mode.fields.m[field_index] = complex_from_doubles(real_val, imag_val);
					}
			}

		}
	}
	if(dimension == 2){

		for(i = 0; i < file_mode.i; i++){		
			for(j = 0; j < file_mode.j; j++){
				for(k = 0; k < 3; k++){

					fscanf(mode_data, "%lf", &real_val);
					fscanf(mode_data, "%lf", &imag_val); 

					ep_index = 3*(i + j*file_mode.i) + k;
					file_mode.ep.m[ep_index] = complex_from_doubles(real_val, imag_val);
				}
				for(k = 0; k < 3; k++){

					fscanf(mode_data, "%lf", &real_val);
					fscanf(mode_data, "%lf", &imag_val); 

					ep_index = 3*(i + j*file_mode.i) + k;
					file_mode.mu.m[ep_index] = complex_from_doubles(real_val, imag_val);
				}
				for(k = 0; k < 6; k++){

					fscanf(mode_data, "%lf", &real_val);
					fscanf(mode_data, "%lf", &imag_val); 

					field_index = 6*(i + j*file_mode.i) + k;
					file_mode.fields.m[field_index] = complex_from_doubles(real_val, imag_val);
				}
			}

		}
	}
	
	return(file_mode);
}

double mode_max_index(mode psi){
	
	int i;
	complex n;
	double max;
	
	max = 0;
	for(i = 0; i < psi.i*psi.j; i++){
		
		n = complex_sqrt(psi.ep.m[i]);
		if(n.re > max){
			
			max = n.re;
		}
	}
	
	return(max);
}

void mode_initialize_ones(mode psi){
	
	int i, index;
	
	for(i = 0; i < psi.i*psi.j; i++){
		
		index = 6*i;
		psi.fields.m[index] = complex_from_doubles(1, 0);
		psi.fields.m[index + 1] = complex_from_doubles(1, 0);
		psi.fields.m[index + 3] = complex_from_doubles(1, 0);
		psi.fields.m[index + 4] = complex_from_doubles(1, 0);
	}
}

void mode_initialize_general(mode psi){
	
	int i, index, x_center, y_center;
	double x, y;
	double value;
	
	x_center = psi.i/2;
	y_center = psi.j/2;
	
	//initialize to 1 + x + y centered at (x_center, y_center)
	
	for(i = 0; i < psi.i*psi.j; i++){
		
		index = 6*i;
		x = (double)(i%psi.i);
		y = (double)(i)/psi.i;
		
		value = (1 + ((x - x_center)*3)/psi.i + ((y - y_center)*1)/psi.j);
		
		//psi.fields.m[index] = complex_from_doubles(value, 0);
		//psi.fields.m[index + 1] = complex_from_doubles(value, 0);
		//psi.fields.m[index + 2] = complex_zero();
		psi.fields.m[index + 3] = complex_from_doubles(value, 0);
		psi.fields.m[index + 4] = complex_from_doubles(value, 0);
		//psi.fields.m[index + 5] = complex_zero();
	}
}

void mode_initialize_rand(mode psi){
	
	int i, index;
	
	for(i = 0; i < psi.i*psi.j; i++){
		
		index = 6*i;
		//psi.fields.m[index] = complex_from_doubles(rand(), 0);
		//psi.fields.m[index + 1] = complex_from_doubles(rand(), 0);
		//psi.fields.m[index + 2] = complex_from_doubles(rand(), 0);
		psi.fields.m[index + 3] = complex_from_doubles(rand(), 0);
		psi.fields.m[index + 4] = complex_from_doubles(rand(), 0);
		//psi.fields.m[index + 5] = complex_from_doubles(rand(), 0);
	}
}

complex mode_inner_product(mode psi, mode phi){
	
	int i, j, index;
	complex sum;
	
	//ensure that psi and phi are the same size
	if((psi.i != phi.i)||(psi.j != phi.j)||(psi.disc != phi.disc)){
		
		printf("Mode dimensions do not match for addition!\n");
		abort();
	}

	sum = complex_zero();
	for(i = 0; i < psi.i*psi.j; i++){
		
		index = 6*i;
		
		sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 0]), phi.fields.m[index + 0]));
		sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 1]), phi.fields.m[index + 1]));
		sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 2]), phi.fields.m[index + 2]));
		sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 3]), phi.fields.m[index + 3]));
		sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 4]), phi.fields.m[index + 4]));
		sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 5]), phi.fields.m[index + 5]));
	}
	
	sum = complex_smul((psi.disc*psi.disc), sum);
	return(sum);	
}

complex mode_H_inner_product(mode psi, mode phi){
	
	int i, index;
	int field_index;
	complex sum;
	
	//ensure that psi and phi are the same size
	if((psi.i != phi.i)||(psi.j != phi.j)||(psi.disc != phi.disc)){
		
		printf("Mode dimensions do not match for addition!\n");
		abort();
	}
	
	
	//Sum H_fields only
	sum = complex_zero();
	for(i = 0; i < psi.i*psi.j; i++){
		
		index = 6*i;
		
		sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 3]), phi.fields.m[index + 3]));
		sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 4]), phi.fields.m[index + 4]));
		//sum = complex_add(sum, complex_mul(complex_conj(psi.fields.m[index + 5]), phi.fields.m[index + 5]));
	}
	
	
	sum = complex_smul((psi.disc*psi.disc), sum);
	return(sum);
}

complex mode_power_inner_product(mode psi, mode phi){
	
	int i, j;
	int field_index;
	complex Ex_phi, Ey_phi, Hx_phi, Hy_phi;
	complex Ex_psi, Ey_psi, Hx_psi, Hy_psi;
	complex product1, product2, product3, product4, sum;
	
	//ensure that psi and phi are the same size
	if((psi.i != phi.i)||(psi.j != phi.j)||(psi.disc != phi.disc)){
		
		printf("Mode dimensions do not match for addition!\n");
		abort();
	}
	
	sum = complex_zero();
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
				
			field_index = 6*(i + j*psi.i);
			
			//pick out field elements
			Ex_psi = psi.fields.m[field_index];
			Ey_psi = psi.fields.m[field_index + 1];
			Hx_psi = psi.fields.m[field_index + 3];
			Hy_psi = psi.fields.m[field_index + 4];
			
			Ex_phi = phi.fields.m[field_index];
			Ey_phi = phi.fields.m[field_index + 1];
			Hx_phi = phi.fields.m[field_index + 3];
			Hy_phi = phi.fields.m[field_index + 4];
			
			product1 = complex_mul(complex_conj(Hx_psi), Ey_phi);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Hy_psi), Ex_phi));
			product3 = complex_mul(Hx_phi, complex_conj(Ey_psi));
			product4 = complex_sub(complex_zero(), complex_mul(Hy_phi, complex_conj(Ex_psi)));
			
			sum = complex_add(sum, complex_add(complex_add(product1, product2), complex_add(product3, product4)));
		}
	}
	
	sum = complex_smul((-0.5*psi.disc*psi.disc), sum);
	return(sum);
}

double mode_norm(mode psi){
	
	double norm;
	
	norm = sqrt(complex_mag(mode_inner_product(psi, psi)));
	
	return(norm);
}

double mode_H_norm(mode psi){
	
	double norm;
	
	norm = sqrt(complex_mag(mode_H_inner_product(psi, psi)));
	
	return(norm);
}

double mode_power_norm(mode psi){
	
	double norm;
	
	norm = sqrt(complex_mag(mode_power_inner_product(psi, psi)));
	
	return(norm);
}

void mode_normalize(mode psi){
	
	int i;
	double norm, divisor;
	
	norm = mode_norm(psi);
	divisor = 1.0/norm;
	
	for(i = 0; i < 6*psi.i*psi.j; i++){
		
		psi.fields.m[i] = complex_smul(divisor, psi.fields.m[i]);
	}
}

void mode_H_normalize(mode psi){
	
	int i, index;
	double norm, divisor;
	
	norm = mode_H_norm(psi);
	divisor = 1.0/norm;
	
	for(i = 0; i < 6*psi.i*psi.j; i++){
		
		psi.fields.m[i] = complex_smul(divisor, psi.fields.m[i]);
	}
}

void mode_power_normalize(mode psi){
	
	int i;
	double norm, divisor;
	
	norm = mode_power_norm(psi);
	divisor = 1.0/norm;
	
	for(i = 0; i < 6*psi.i*psi.j; i++){
		
		psi.fields.m[i] = complex_smul(divisor, psi.fields.m[i]);
	}
}

void mode_apply_zero_bc(mode psi){
	
	int i, j;
	int right_index, left_index, top_index, bottom_index;
	
	//left and right edges
	for(j = 0; j < psi.j; j++){
		
		right_index =  6*((j + 1)*psi.i - 1);
		left_index = 6*j*psi.i;
		
		psi.fields.m[right_index] = complex_zero();
		psi.fields.m[right_index + 1] = complex_zero();
		psi.fields.m[right_index + 2] = complex_zero();
		psi.fields.m[right_index + 3] = complex_zero();
		psi.fields.m[right_index + 4] = complex_zero();
		psi.fields.m[right_index + 5] = complex_zero();
		
		psi.fields.m[left_index] = complex_zero();
		psi.fields.m[left_index + 1] = complex_zero();
		psi.fields.m[left_index + 2] = complex_zero();
		psi.fields.m[left_index + 3] = complex_zero();
		psi.fields.m[left_index + 4] = complex_zero();
		psi.fields.m[left_index + 5] = complex_zero();
	}
	
	//top and bottom edges
	for(i = 1; i < (psi.i - 1); i++){
		
		top_index =  6*(i + (psi.j - 1)*psi.i);
		bottom_index = 6*i;
		
		psi.fields.m[top_index] = complex_zero();
		psi.fields.m[top_index + 1] = complex_zero();
		psi.fields.m[top_index + 2] = complex_zero();
		psi.fields.m[top_index + 3] = complex_zero();
		psi.fields.m[top_index + 4] = complex_zero();
		psi.fields.m[top_index + 5] = complex_zero();
		
		psi.fields.m[bottom_index] = complex_zero();
		psi.fields.m[bottom_index + 1] = complex_zero();
		psi.fields.m[bottom_index + 2] = complex_zero();
		psi.fields.m[bottom_index + 3] = complex_zero();
		psi.fields.m[bottom_index + 4] = complex_zero();
		psi.fields.m[bottom_index + 5] = complex_zero();
	}
}

void mode_apply_zero_bc_vertical(mode psi){
	
	int i, j;
	int top_index, bottom_index;
	
	//top and bottom edges
	for(i = 1; i < (psi.i - 1); i++){
		
		top_index =  6*(i + (psi.j - 1)*psi.i);
		bottom_index = 6*i;
		
		psi.fields.m[top_index] = complex_zero();
		psi.fields.m[top_index + 1] = complex_zero();
		psi.fields.m[top_index + 2] = complex_zero();
		psi.fields.m[top_index + 3] = complex_zero();
		psi.fields.m[top_index + 4] = complex_zero();
		psi.fields.m[top_index + 5] = complex_zero();
		
		psi.fields.m[bottom_index] = complex_zero();
		psi.fields.m[bottom_index + 1] = complex_zero();
		psi.fields.m[bottom_index + 2] = complex_zero();
		psi.fields.m[bottom_index + 3] = complex_zero();
		psi.fields.m[bottom_index + 4] = complex_zero();
		psi.fields.m[bottom_index + 5] = complex_zero();
	}
	
}

void mode_apply_zero_bc_horizontal(mode psi){
	
	int i, j;
	int right_index, left_index;
	
	//left and right edges
	for(j = 0; j < psi.j; j++){
		
		right_index =  6*((j + 1)*psi.i - 1);
		left_index = 6*j*psi.i;
		
		psi.fields.m[right_index] = complex_zero();
		psi.fields.m[right_index + 1] = complex_zero();
		psi.fields.m[right_index + 2] = complex_zero();
		psi.fields.m[right_index + 3] = complex_zero();
		psi.fields.m[right_index + 4] = complex_zero();
		psi.fields.m[right_index + 5] = complex_zero();
		
		psi.fields.m[left_index] = complex_zero();
		psi.fields.m[left_index + 1] = complex_zero();
		psi.fields.m[left_index + 2] = complex_zero();
		psi.fields.m[left_index + 3] = complex_zero();
		psi.fields.m[left_index + 4] = complex_zero();
		psi.fields.m[left_index + 5] = complex_zero();
	}
}

void mode_apply_zero_gradient_bc(mode psi){
	
	int i, j;
	int right_index, left_index, top_index, bottom_index;
	
	//left and right edges
	for(j = 0; j < psi.j; j++){
		
		right_index =  6*((j + 1)*psi.i - 1);
		left_index = 6*j*psi.i;
		
		psi.fields.m[right_index] = psi.fields.m[right_index - 6];
		psi.fields.m[right_index + 1] = psi.fields.m[right_index + 1 - 6];
		psi.fields.m[right_index + 2] = psi.fields.m[right_index + 2 - 6];
		psi.fields.m[right_index + 3] = psi.fields.m[right_index + 3 - 6];
		psi.fields.m[right_index + 4] = psi.fields.m[right_index + 4 - 6];
		psi.fields.m[right_index + 5] = psi.fields.m[right_index + 5 - 6];
		
		psi.fields.m[left_index] = psi.fields.m[left_index + 6];
		psi.fields.m[left_index + 1] = psi.fields.m[left_index + 2 + 6];
		psi.fields.m[left_index + 2] = psi.fields.m[left_index + 3 + 6];
		psi.fields.m[left_index + 3] = psi.fields.m[left_index + 4 + 6];
		psi.fields.m[left_index + 4] = psi.fields.m[left_index + 5 + 6];
		psi.fields.m[left_index + 5] = psi.fields.m[left_index + 6 + 6];
	}
	
	//top and bottom edges
	for(i = 1; i < (psi.i - 1); i++){
		
		top_index =  6*(i + (psi.j - 1)*psi.i);
		bottom_index = 6*i;
		
		psi.fields.m[top_index] = psi.fields.m[top_index - 6*psi.i];
		psi.fields.m[top_index + 1] = psi.fields.m[top_index + 1 - 6*psi.i];
		psi.fields.m[top_index + 2] = psi.fields.m[top_index + 2 - 6*psi.i];
		psi.fields.m[top_index + 3] = psi.fields.m[top_index + 3 - 6*psi.i];
		psi.fields.m[top_index + 4] = psi.fields.m[top_index + 4 - 6*psi.i];
		psi.fields.m[top_index + 5] = psi.fields.m[top_index + 5 - 6*psi.i];
		
		psi.fields.m[bottom_index] = psi.fields.m[bottom_index + 6*psi.i];
		psi.fields.m[bottom_index + 1] = psi.fields.m[bottom_index + 1 + 6*psi.i];
		psi.fields.m[bottom_index + 2] = psi.fields.m[bottom_index + 2 + 6*psi.i];
		psi.fields.m[bottom_index + 3] = psi.fields.m[bottom_index + 3 + 6*psi.i];
		psi.fields.m[bottom_index + 4] = psi.fields.m[bottom_index + 4 + 6*psi.i];
		psi.fields.m[bottom_index + 5] = psi.fields.m[bottom_index + 5 + 6*psi.i];
	}
}

void mode_apply_zero_gradient_bc_vertical(mode psi){

	int i, j;
	int top_index, bottom_index;
	
	//top and bottom edges
	for(i = 1; i < (psi.i - 1); i++){
		
		top_index =  6*(i + (psi.j - 1)*psi.i);
		bottom_index = 6*i;
		
		psi.fields.m[top_index] = psi.fields.m[top_index - 6*psi.i];
		psi.fields.m[top_index + 1] = psi.fields.m[top_index + 1 - 6*psi.i];
		psi.fields.m[top_index + 2] = psi.fields.m[top_index + 2 - 6*psi.i];
		psi.fields.m[top_index + 3] = psi.fields.m[top_index + 3 - 6*psi.i];
		psi.fields.m[top_index + 4] = psi.fields.m[top_index + 4 - 6*psi.i];
		psi.fields.m[top_index + 5] = psi.fields.m[top_index + 5 - 6*psi.i];
		
		psi.fields.m[bottom_index] = psi.fields.m[bottom_index + 6*psi.i];
		psi.fields.m[bottom_index + 1] = psi.fields.m[bottom_index + 1 + 6*psi.i];
		psi.fields.m[bottom_index + 2] = psi.fields.m[bottom_index + 2 + 6*psi.i];
		psi.fields.m[bottom_index + 3] = psi.fields.m[bottom_index + 3 + 6*psi.i];
		psi.fields.m[bottom_index + 4] = psi.fields.m[bottom_index + 4 + 6*psi.i];
		psi.fields.m[bottom_index + 5] = psi.fields.m[bottom_index + 5 + 6*psi.i];
	}
	
}

void mode_apply_zero_gradient_bc_horizontal(mode psi){
	
	int i, j;
	int right_index, left_index;
	
	//left and right edges
	for(j = 0; j < psi.j; j++){
		
		right_index =  6*((j + 1)*psi.i - 1);
		left_index = 6*j*psi.i;
		
		psi.fields.m[right_index] = psi.fields.m[right_index - 6];
		psi.fields.m[right_index + 1] = psi.fields.m[right_index + 1 - 6];
		psi.fields.m[right_index + 2] = psi.fields.m[right_index + 2 - 6];
		psi.fields.m[right_index + 3] = psi.fields.m[right_index + 3 - 6];
		psi.fields.m[right_index + 4] = psi.fields.m[right_index + 4 - 6];
		psi.fields.m[right_index + 5] = psi.fields.m[right_index + 5 - 6];
		
		psi.fields.m[left_index] = psi.fields.m[left_index + 6];
		psi.fields.m[left_index + 1] = psi.fields.m[left_index + 1 + 6];
		psi.fields.m[left_index + 2] = psi.fields.m[left_index + 2 + 6];
		psi.fields.m[left_index + 3] = psi.fields.m[left_index + 3 + 6];
		psi.fields.m[left_index + 4] = psi.fields.m[left_index + 4 + 6];
		psi.fields.m[left_index + 5] = psi.fields.m[left_index + 5 + 6];
	}
}

void mode_apply_zero_bc_pml(mode psi){
	
	int i, j;
	int index, field_index, pml_index;
	
	//top and bottom edges
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = i + j*psi.i;
			pml_index = 2*index;
			field_index = 6*index;

			if((psi.pml.m[pml_index].im > 1e-9)||(psi.pml.m[pml_index + 1].im > 1e-9)){

				psi.fields.m[field_index] = complex_zero();
				psi.fields.m[field_index + 1] = complex_zero();
				psi.fields.m[field_index + 2] = complex_zero();
				psi.fields.m[field_index + 3] = complex_zero();
				psi.fields.m[field_index + 4] = complex_zero();
				psi.fields.m[field_index + 5] = complex_zero();
			}
		}
	}
	
}

void mode_apply_bc(mode psi, int confinement){
	
	//ensure that confinement is 1d or 2d
	if((confinement != 1)&&(confinement != 2)){
		
		printf("Mode confinement dimension is not 1 or 2!\n");
		abort();
	}
	
	if(confinement == 1){
		
		mode_apply_zero_bc_vertical(psi);
		mode_apply_zero_gradient_bc_horizontal(psi);
	}
	else{
		
		mode_apply_zero_bc(psi);
	}
}

void mode_set_pml(mode psi, double max_ep_im, int order, int pml_thickness){
	
	int i, j, pml_index;
	int x_size, y_size;
	double pml_step;
	
	//ensure that order is 1, 2 or 3;
	if((order != 0)&&(order != 1)&&(order != 2)&&(order != 3)&&(order != 4)){
		
		printf("Order not valid for set pml!\n");
		abort();
	}
	
	x_size = psi.i;
	y_size = psi.j;
	
	if(order == 0){
		
		pml_step = max_ep_im;
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, pml_step);
				}
				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, pml_step);
				}
				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, pml_step);
				} 
				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, pml_step);
				}
			}
		}
	}
	else if(order == 1){
		
		pml_step = max_ep_im/pml_thickness;
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*pml_step);
				}
				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*pml_step);
				}
				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*pml_step);
				} 
				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 2){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*(pml_thickness - i)*pml_step);
				}
				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*pml_step);
				}
				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*(pml_thickness - j)*pml_step);
				} 
				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 3){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*(pml_thickness - i)*(pml_thickness - i)*pml_step);
				}
				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*pml_step);
				}
				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*(pml_thickness - j)*(pml_thickness - j)*pml_step);
				} 
				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 4){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*(pml_thickness - i)*(pml_thickness - i)*(pml_thickness - i)*pml_step);
				}
				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*pml_step);
				}
				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*(pml_thickness - j)*(pml_thickness - j)*(pml_thickness - j)*pml_step);
				} 
				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
}

void mode_set_pml_bottom(mode psi, double max_ep_im, int order, int pml_thickness){
	
	int i, j, pml_index;
	int x_size, y_size;
	double pml_step;

		//ensure that order is 1, 2 or 3;
	if((order != 1)&&(order != 2)&&(order != 3)&&(order != 4)){

		printf("Order not valid for set pml!\n");
		abort();
	}

	x_size = psi.i;
	y_size = psi.j;

	if(order == 1){

		pml_step = max_ep_im/pml_thickness;

		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*pml_step);
				} 
			}
		}
	}
	else if(order == 2){

		pml_step = max_ep_im/(pml_thickness*pml_thickness);

		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*(pml_thickness - j)*pml_step);
				}
			}
		}
	}
	else if(order == 3){

		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness);

		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*(pml_thickness - j)*(pml_thickness - j)*pml_step);
				}
			}
		}
	}
	else if(order == 4){

		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness*pml_thickness);

		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(j < pml_thickness){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*(pml_thickness - j)*(pml_thickness - j)*(pml_thickness - j)*pml_step);
				}
			}
		}
	}
}

void mode_set_pml_top(mode psi, double max_ep_im, int order, int pml_thickness){
	
	int i, j, pml_index;
	int x_size, y_size;
	double pml_step;
	
	//ensure that order is 1, 2 or 3;
	if((order != 1)&&(order != 2)&&(order != 3)&&(order != 4)){
		
		printf("Order not valid for set pml!\n");
		abort();
	}
	
	x_size = psi.i;
	y_size = psi.j;
	
	if(order == 1){
		
		pml_step = max_ep_im/pml_thickness;
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 2){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 3){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 4){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(j >= (y_size - pml_thickness)){

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
}

void mode_set_pml_left(mode psi, double max_ep_im, int order, int pml_thickness){
	
	int i, j, pml_index;
	int x_size, y_size;
	double pml_step;
	
	//ensure that order is 1, 2 or 3;
	if((order != 1)&&(order != 2)&&(order != 3)&&(order != 4)){
		
		printf("Order not valid for set pml!\n");
		abort();
	}
	
	x_size = psi.i;
	y_size = psi.j;
	
	if(order == 1){
		
		pml_step = max_ep_im/pml_thickness;
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*pml_step);
				}
			}
		}
	}
	else if(order == 2){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*(pml_thickness - i)*pml_step);
				}
			}
		}
	}
	else if(order == 3){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*(pml_thickness - i)*(pml_thickness - i)*pml_step);
				}
			}
		}
	}
	else if(order == 4){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i < pml_thickness){

					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*(pml_thickness - i)*(pml_thickness - i)*(pml_thickness - i)*pml_step);
				}
			}
		}
	}
}

void mode_set_pml_right(mode psi, double max_ep_im, int order, int pml_thickness){
	
	int i, j, pml_index;
	int x_size, y_size;
	double pml_step;
	
	//ensure that order is 1, 2 or 3;
	if((order != 1)&&(order != 2)&&(order != 3)&&(order != 4)){
		
		printf("Order not valid for set pml!\n");
		abort();
	}
	
	x_size = psi.i;
	y_size = psi.j;
	
	if(order == 1){
		
		pml_step = max_ep_im/pml_thickness;
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 2){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 3){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
	else if(order == 4){
		
		pml_step = max_ep_im/(pml_thickness*pml_thickness*pml_thickness*pml_thickness);
		
		for(i = 0; i < psi.i; i++){
			for(j = 0; j < psi.j; j++){

				pml_index = 2*(i + j*(psi.i));

				if(i >= (x_size - pml_thickness)){

					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
}

void mode_set_grade_to_constant_pml(mode psi, double max_ep_im, int grade_thickness, int pml_thickness){
	
	int i, j, pml_index;
	int x_size, y_size;
	int constant_thickness;
	double pml_step;
	
	//ensure that order is 1, 2 or 3;
	if(pml_thickness < grade_thickness){
		
		printf("PML thickness is less than grade thickness!\n");
		abort();
	}
	
	x_size = psi.i;
	y_size = psi.j;
	
	pml_step = max_ep_im/(grade_thickness*grade_thickness);
	//pml_step = max_ep_im/grade_thickness;
	
	constant_thickness = pml_thickness - grade_thickness;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			pml_index = 2*(i + j*(psi.i));

			if(i < pml_thickness){
				
				if(i < constant_thickness){
					
					psi.pml.m[pml_index] = complex_from_doubles(1, max_ep_im);
				}
				else{
					
					psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*(pml_thickness - i)*pml_step);
					//psi.pml.m[pml_index] = complex_from_doubles(1, (pml_thickness - i)*pml_step);
				}
			}
			if(i >= (x_size - pml_thickness)){

				if(i >= (x_size - constant_thickness)){
					
					psi.pml.m[pml_index] = complex_from_doubles(1, max_ep_im);
				}
				else{
					
					psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*(i - (x_size - 1 - pml_thickness))*pml_step);
					//psi.pml.m[pml_index] = complex_from_doubles(1, (i - (x_size - 1 - pml_thickness))*pml_step);
				}
			}
			if(j < pml_thickness){

				if(j < constant_thickness){
					
					psi.pml.m[pml_index + 1] = complex_from_doubles(1, max_ep_im);
				}
				else{

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*(pml_thickness - j)*pml_step);
					//psi.pml.m[pml_index + 1] = complex_from_doubles(1, (pml_thickness - j)*pml_step);
				}
			} 
			if(j >= (y_size - pml_thickness)){

				if(j >= (y_size - constant_thickness)){
					
					psi.pml.m[pml_index + 1] = complex_from_doubles(1, max_ep_im);
				}
				else{

					psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*(j - (y_size - 1 - pml_thickness))*pml_step);
					//psi.pml.m[pml_index + 1] = complex_from_doubles(1, (j - (y_size - 1 - pml_thickness))*pml_step);
				}
			}
		}
	}
}

void mode_set_background(mode psi, complex ep, complex mu){
	
	int i, j, index;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			index = 3*(i + j*(psi.i));
		
			psi.ep.m[index] = ep;
			psi.mu.m[index] = mu;
			psi.ep.m[index + 1] = ep;
			psi.mu.m[index + 1] = mu;
			psi.ep.m[index + 2] = ep;
			psi.mu.m[index + 2] = mu;
	

		}
	}
}

void mode_set_bottom_half_plane(mode psi, complex ep, complex mu, int top){
	
	int i, j, index;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			index = 3*(i + j*(psi.i));
			if(j < top){
					
				psi.ep.m[index] = ep;
				psi.mu.m[index] = mu;
				psi.ep.m[index + 1] = ep;
				psi.mu.m[index + 1] = mu;
				psi.ep.m[index + 2] = ep;
				psi.mu.m[index + 2] = mu;
	
			}
		}
	}
}

void mode_set_top_half_plane(mode psi, complex ep, complex mu, int bottom){
	
	int i, j, index;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			index = 3*(i + j*(psi.i));
			if(j >= bottom){
					
				psi.ep.m[index] = ep;
				psi.mu.m[index] = mu;
				psi.ep.m[index + 1] = ep;
				psi.mu.m[index + 1] = mu;
				psi.ep.m[index + 2] = ep;
				psi.mu.m[index + 2] = mu;
	
			}
		}
	}
}

void mode_set_left_half_plane(mode psi, complex ep, complex mu, int right){
	
	int i, j, index;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			index = 3*(i + j*(psi.i));
			if(i < right){
					
				psi.ep.m[index] = ep;
				psi.mu.m[index] = mu;
				psi.ep.m[index + 1] = ep;
				psi.mu.m[index + 1] = mu;
				psi.ep.m[index + 2] = ep;
				psi.mu.m[index + 2] = mu;
	
			}
		}
	}
}

void mode_set_right_half_plane(mode psi, complex ep, complex mu, int left){
	
	int i, j, index;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			index = 3*(i + j*(psi.i));
			if(i >= left){
					
				psi.ep.m[index] = ep;
				psi.mu.m[index] = mu;
				psi.ep.m[index + 1] = ep;
				psi.mu.m[index + 1] = mu;
				psi.ep.m[index + 2] = ep;
				psi.mu.m[index + 2] = mu;
	
			}
		}
	}
}

void mode_make_rectangle(mode psi, complex ep, complex mu, int left, int right, int bottom, int top){
	
	int i, j, index;
	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			index = 3*(i + j*(psi.i));
			if((i >= left)&&(i < right)&&(j >= bottom)&&(j < top)){
					
				psi.ep.m[index] = ep;
				psi.mu.m[index] = mu;
				psi.ep.m[index + 1] = ep;
				psi.mu.m[index + 1] = mu;
				psi.ep.m[index + 2] = ep;
				psi.mu.m[index + 2] = mu;
	
			}
		}
	}
}

void mode_apply_H_curl(mode psi, complex beta, double omega){
	
	int i, j;
	int ep_index, pml_index;
	int field_index, field_index_right, field_index_left, field_index_up, field_index_down;
	double divisor;
	complex Hx, Hy, Hz;
	complex Hz_up, Hz_down, Hz_right, Hz_left;
	complex Hx_up, Hx_down;
	complex Hy_right, Hy_left;
	complex dx_Hy, dx_Hz, dy_Hx, dy_Hz;
	complex inv_ep_x, inv_ep_y, inv_ep_z;
	complex sx, sy;
	
	divisor = 1/(psi.disc);
	
	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
		
			ep_index = 3*(i + j*psi.i);
			
			field_index = 2*ep_index;
			field_index_right = field_index + 6;
			field_index_left = field_index - 6;
			field_index_up = field_index + 6*psi.i;
			field_index_down = field_index - 6*psi.i;
			
			//pick out field elements
			Hx = psi.fields.m[field_index + 3];
			Hy = psi.fields.m[field_index + 4];
			Hz = psi.fields.m[field_index + 5];
			
			Hy_right = psi.fields.m[field_index_right + 4];
			Hz_right = psi.fields.m[field_index_right + 5];	
			
			Hy_left = psi.fields.m[field_index_left + 4];
			Hz_left = psi.fields.m[field_index_left + 5];
			
			Hx_up = psi.fields.m[field_index_up + 3];
			Hz_up = psi.fields.m[field_index_up + 5];

			Hx_down = psi.fields.m[field_index_down + 3];
			Hz_down = psi.fields.m[field_index_down + 5];	
			
			sx = psi.pml.m[pml_index];
			sy = psi.pml.m[pml_index + 1];	
			
			//pick out derivatives as central differences
			dx_Hy = complex_div(complex_smul(divisor, complex_sub(Hy, Hy_left)), sx);
			dx_Hz = complex_div(complex_smul(divisor, complex_sub(Hz, Hz_left)), sx);
			
			dy_Hx = complex_div(complex_smul(divisor, complex_sub(Hx, Hx_down)), sy);
			dy_Hz = complex_div(complex_smul(divisor, complex_sub(Hz, Hz_down)), sy);
			
			//insert psi E-field values i*omega*E = 1/ep*curl(H)
			inv_ep_x = complex_reciprocal(psi.ep.m[ep_index]);
			inv_ep_y = complex_reciprocal(psi.ep.m[ep_index + 1]);
			inv_ep_z = complex_reciprocal(psi.ep.m[ep_index + 2]);
			
			psi.fields.m[field_index] = complex_mul(complex_smul(1/omega, complex_from_doubles(0, -1)), complex_mul(inv_ep_x, complex_add(dy_Hz, complex_mul(complex_mul(complex_from_doubles(0,1), beta), Hy))));
			psi.fields.m[field_index + 1] = complex_mul(complex_smul(1/omega, complex_from_doubles(0, -1)), complex_smul(-1, complex_mul(inv_ep_y, complex_add(dx_Hz, complex_mul(complex_mul(complex_from_doubles(0,1), beta), Hx)))));
			psi.fields.m[field_index + 2] = complex_mul(complex_smul(1/omega, complex_from_doubles(0, -1)), complex_mul(inv_ep_z, complex_sub(dx_Hy, dy_Hx)));
		}
	}
}

void mode_apply_E_curl(mode psi, complex beta, double omega){
	
	int i, j;
	int ep_index, pml_index;
	int field_index, field_index_right, field_index_left, field_index_up, field_index_down;
	double divisor;
	complex Ex, Ey, Ez;
	complex Ez_up, Ez_down, Ez_right, Ez_left;
	complex Ex_up, Ex_down;
	complex Ey_right, Ey_left;
	complex dx_Ey, dx_Ez, dy_Ex, dy_Ez;
	complex inv_mu_x, inv_mu_y, inv_mu_z;
	complex sx, sy;
	
	divisor = 1/(psi.disc);
	
	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
		
			ep_index = 3*(i + j*psi.i);
			
			field_index = 2*ep_index;
			field_index_right = field_index + 6;
			field_index_left = field_index - 6;
			field_index_up = field_index + 6*psi.i;
			field_index_down = field_index - 6*psi.i;
			
			//pick out field elements
			Ex = psi.fields.m[field_index];
			Ey = psi.fields.m[field_index + 1];
			Ez = psi.fields.m[field_index + 2];
			
			Ey_right = psi.fields.m[field_index_right + 1];
			Ez_right = psi.fields.m[field_index_right + 2];	
			
			Ey_left = psi.fields.m[field_index_left + 1];
			Ez_left = psi.fields.m[field_index_left + 2];
			
			Ex_up = psi.fields.m[field_index_up];
			Ez_up = psi.fields.m[field_index_up + 2];

			Ex_down = psi.fields.m[field_index_down];
			Ez_down = psi.fields.m[field_index_down + 2];	
			
			sx = psi.pml.m[pml_index];
			sy = psi.pml.m[pml_index + 1];	
			
			//pick out derivatives as central differences
			dx_Ey = complex_div(complex_smul(divisor, complex_sub(Ey, Ey_left)), sx);
			dx_Ez = complex_div(complex_smul(divisor, complex_sub(Ez, Ez_left)), sx);
			
			dy_Ex = complex_div(complex_smul(divisor, complex_sub(Ex, Ex_down)), sy);
			dy_Ez = complex_div(complex_smul(divisor, complex_sub(Ez, Ez_down)), sy);
			
			//insert psi E-field values i*omega*H = -1/mu*curl(E)
			inv_mu_x = complex_reciprocal(psi.mu.m[ep_index]);
			inv_mu_y = complex_reciprocal(psi.mu.m[ep_index + 1]);
			inv_mu_z = complex_reciprocal(psi.mu.m[ep_index + 2]);
			
			psi.fields.m[field_index + 3] = complex_mul(complex_smul(-1/omega, complex_from_doubles(0, -1)), complex_mul(inv_mu_x, complex_add(dy_Ez, complex_mul(complex_mul(complex_from_doubles(0,1), beta), Ey))));
			psi.fields.m[field_index + 4] = complex_mul(complex_smul(-1/omega, complex_from_doubles(0, -1)), complex_smul(-1, complex_mul(inv_mu_y, complex_add(dx_Ez, complex_mul(complex_mul(complex_from_doubles(0,1), beta), Ex)))));
			psi.fields.m[field_index + 5] = complex_mul(complex_smul(-1/omega, complex_from_doubles(0, -1)), complex_mul(inv_mu_z, complex_sub(dx_Ey, dy_Ex)));
		}
	}
}

void mode_apply_radial_H_curl(mode psi, complex beta, double omega, double radius_left){
	
	int i, j;
	int ep_index, pml_index;
	int field_index, field_index_right, field_index_left, field_index_up, field_index_down;
	double radius, divisor;
	complex Hx, Hy, Hz;
	complex Hz_up, Hz_down, Hz_right, Hz_left;
	complex Hx_up, Hx_down;
	complex Hy_right, Hy_left;
	complex dx_Hy, dx_Hz, dy_Hx, dy_Hz;
	complex inv_ep_x, inv_ep_y, inv_ep_z;
	complex sx, sy;
	
	divisor = 1/(psi.disc);
	
	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		
		radius = radius_left + i*psi.disc;
		
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
		
			ep_index = 3*(i + j*psi.i);
			
			field_index = 2*ep_index;
			field_index_right = field_index + 6;
			field_index_left = field_index - 6;
			field_index_up = field_index + 6*psi.i;
			field_index_down = field_index - 6*psi.i;
			
			//pick out field elements
			Hx = psi.fields.m[field_index + 3];
			Hy = psi.fields.m[field_index + 4];
			Hz = psi.fields.m[field_index + 5];
			
			Hy_right = psi.fields.m[field_index_right + 4];
			Hz_right = psi.fields.m[field_index_right + 5];	
			
			Hy_left = psi.fields.m[field_index_left + 4];
			Hz_left = psi.fields.m[field_index_left + 5];
			
			Hx_up = psi.fields.m[field_index_up + 3];
			Hz_up = psi.fields.m[field_index_up + 5];

			Hx_down = psi.fields.m[field_index_down + 3];
			Hz_down = psi.fields.m[field_index_down + 5];	
			
			sx = psi.pml.m[pml_index];
			sy = psi.pml.m[pml_index + 1];	
			
			//pick out derivatives as central differences
			dx_Hy = complex_div(complex_smul(divisor, complex_sub(Hy, Hy_left)), sx);
			dx_Hz = complex_div(complex_smul(divisor, complex_sub(Hz, Hz_left)), sx);
			
			dy_Hx = complex_div(complex_smul(divisor, complex_sub(Hx, Hx_down)), sy);
			dy_Hz = complex_div(complex_smul(divisor, complex_sub(Hz, Hz_down)), sy);
			
			//insert psi E-field values i*omega*E = 1/ep*curl(H)
			inv_ep_x = complex_reciprocal(psi.ep.m[ep_index]);
			inv_ep_y = complex_reciprocal(psi.ep.m[ep_index + 1]);
			inv_ep_z = complex_reciprocal(psi.ep.m[ep_index + 2]);
			
			psi.fields.m[field_index] = complex_mul(complex_smul(1/omega, complex_from_doubles(0, -1)), complex_mul(inv_ep_x, complex_add(dy_Hz, complex_mul(complex_mul(complex_from_doubles(0,1), complex_smul(1/radius, beta)), Hy))));
			psi.fields.m[field_index + 1] = complex_mul(complex_smul(1/omega, complex_from_doubles(0, -1)), complex_smul(-1, complex_mul(inv_ep_y, complex_add(dx_Hz, complex_mul(complex_mul(complex_from_doubles(0,1), complex_smul(1/radius, beta)), Hx)))));
			psi.fields.m[field_index + 2] = complex_mul(complex_smul(1/omega, complex_from_doubles(0, -1)), complex_mul(inv_ep_z, complex_sub(dx_Hy, dy_Hx)));
		}
	}
}

mode mode_apply_beta_operator(mode psi, double omega){

	int i, j;
	int pml_index;
	int ep_index, ep_index_right, ep_index_left, ep_index_up, ep_index_down;
	int field_index, field_index_right, field_index_left, field_index_up, field_index_down;
	double divisor;	
	complex dx_Dx, dx_Ez, dy_Dy, dy_Ez;
	complex dx_Bx, dx_Hz, dy_By, dy_Hz;
	complex Ex, Ey, Ez, Hx, Hy, Hz;
	complex Dx_right, Ez_right, Bx_right, Hz_right;
	complex Dx_left, Ez_left, Bx_left, Hz_left;
	complex Dy_up, Ez_up, By_up, Hz_up;
	complex Dy_down, Ez_down, By_down, Hz_down;
	complex sx, sy;
	mode beta_psi;
	
	beta_psi = mode_new(psi.i, psi.j, psi.disc);
	for(i = 0; i < psi.ep.i; i++){
		
		beta_psi.ep.m[i] = psi.ep.m[i];
		beta_psi.mu.m[i] = psi.mu.m[i];
	}
	for(i = 0; i < 2*psi.i*psi.j; i++){
		
		beta_psi.pml.m[i] = psi.pml.m[i];
	}
	
	divisor = 1/(2*psi.disc);

	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
		
			ep_index = 3*(i + j*psi.i);
			ep_index_right = ep_index + 3;
			ep_index_left = ep_index - 3;
			ep_index_up = ep_index + 3*psi.i;
			ep_index_down = ep_index - 3*psi.i;
			
			field_index = 6*(i + j*psi.i);
			field_index_right = field_index + 6;
			field_index_left = field_index - 6;
			field_index_up = field_index + 6*psi.i;
			field_index_down = field_index - 6*psi.i;
			
			//pick out field elements
			Ex = psi.fields.m[field_index];
			Ey = psi.fields.m[field_index + 1];
			Ez = psi.fields.m[field_index + 2];
			Hx = psi.fields.m[field_index + 3];
			Hy = psi.fields.m[field_index + 4];
			Hz = psi.fields.m[field_index + 5];
			
			Dx_right = complex_mul(psi.ep.m[ep_index_right], psi.fields.m[field_index_right]);
			Ez_right = psi.fields.m[field_index_right + 2];
			Bx_right = complex_mul(psi.mu.m[ep_index_right], psi.fields.m[field_index_right + 3]);
			Hz_right = psi.fields.m[field_index_right + 5];	
			
			Dx_left = complex_mul(psi.ep.m[ep_index_left], psi.fields.m[field_index_left]);
			Ez_left = psi.fields.m[field_index_left + 2];
			Bx_left = complex_mul(psi.mu.m[ep_index_left], psi.fields.m[field_index_left + 3]);
			Hz_left = psi.fields.m[field_index_left + 5];
			
			Dy_up = complex_mul(psi.ep.m[ep_index_up + 1], psi.fields.m[field_index_up + 1]);
			Ez_up = psi.fields.m[field_index_up + 2];
			By_up = complex_mul(psi.mu.m[ep_index_up + 1], psi.fields.m[field_index_up + 4]);
			Hz_up = psi.fields.m[field_index_up + 5];
			
			Dy_down = complex_mul(psi.ep.m[ep_index_down + 1], psi.fields.m[field_index_down + 1]);
			Ez_down = psi.fields.m[field_index_down + 2];
			By_down = complex_mul(psi.mu.m[ep_index_down + 1], psi.fields.m[field_index_down + 4]);
			Hz_down = psi.fields.m[field_index_down + 5];		
			
			sx = psi.pml.m[pml_index];
			sy = psi.pml.m[pml_index + 1];
			
			//pick out derivatives as central differences
			dx_Dx = complex_div(complex_smul(divisor, complex_sub(Dx_right, Dx_left)), sx);
			dx_Ez = complex_div(complex_smul(divisor, complex_sub(Ez_right, Ez_left)), sx);
			dx_Bx = complex_div(complex_smul(divisor, complex_sub(Bx_right, Bx_left)), sx);
			dx_Hz = complex_div(complex_smul(divisor, complex_sub(Hz_right, Hz_left)), sx);
			
			dy_Dy = complex_div(complex_smul(divisor, complex_sub(Dy_up, Dy_down)), sy);
			dy_Ez = complex_div(complex_smul(divisor, complex_sub(Ez_up, Ez_down)), sy);
			dy_By = complex_div(complex_smul(divisor, complex_sub(By_up, By_down)), sy);
			dy_Hz = complex_div(complex_smul(divisor, complex_sub(Hz_up, Hz_down)), sy);
			
			//insert beta_psi values beta_psi = BETA*psi
			beta_psi.fields.m[field_index] = complex_add(complex_mul(complex_smul(omega, psi.mu.m[ep_index + 1]), Hy), complex_mul(complex_from_doubles(0, 1), dx_Ez));
			beta_psi.fields.m[field_index + 1] = complex_add(complex_mul(complex_smul(-omega, psi.mu.m[ep_index]), Hx), complex_mul(complex_from_doubles(0, 1), dy_Ez));
			beta_psi.fields.m[field_index + 2] = complex_div(complex_mul(complex_from_doubles(0, -1), complex_add(dx_Dx, dy_Dy)), psi.ep.m[ep_index + 2]);

			beta_psi.fields.m[field_index + 3] = complex_add(complex_mul(complex_smul(-omega, psi.ep.m[ep_index + 1]), Ey), complex_mul(complex_from_doubles(0, 1), dx_Hz));
			beta_psi.fields.m[field_index + 4] = complex_add(complex_mul(complex_smul(omega, psi.ep.m[ep_index]), Ex), complex_mul(complex_from_doubles(0, 1), dy_Hz));
			beta_psi.fields.m[field_index + 5] = complex_div(complex_mul(complex_from_doubles(0, -1), complex_add(dx_Bx, dy_By)), psi.mu.m[ep_index + 2]);
		}
	}
	
	return(beta_psi);
}

complex mode_solve_beta(mode *psi, double omega, double tol, mode *prev_modes, int prev_num, int size){
	
	int i, j, k, m, spurs, H_index, mode_count, orthogonalizations;
	int max_index;
	double beta_mag, n_mag, n_old, n_diff;
	double max_real_n, overlap, norm, max_2;
	double max;
	complex *eigenvalues;
	complex beta, n_eff, inner_prod;
	mode temp1;
	mode *krylov;
	matrix H, D;
	matrix reduced_eigenvector;
	matrix TEMP;
	
	orthogonalizations = 5;
	
	eigenvalues = malloc(size*sizeof(complex));
	krylov = malloc((size + 1)*sizeof(mode));

	H = matrix_new(size, size);
	
	mode_initialize_general(*psi);
	mode_normalize(*psi);
	
	n_diff = 100;
	n_mag = 100;
	
	m = 0;
	mode_count = 1;
	beta = complex_from_doubles(1, 0);
	while(fabs(n_diff) > tol){
		

		//orthogonalize psi
		/*for(j = 0; j < orthogonalizations; j++){

			mode_H_normalize(*psi);
			for(i = prev_num - 1; i >= 0; i--){

				inner_prod = mode_H_inner_product(prev_modes[i], *psi);
				TEMP = matrix_complex_smul(inner_prod, prev_modes[i].fields);
				matrix_sub_in_place(psi->fields, TEMP);
				mode_H_normalize(*psi);
				matrix_free(TEMP);
			}

			mode_apply_H_curl(*psi, beta, omega);
		}*/
		
		//build Arnoldi matrix
		krylov[0] = mode_copy(*psi);
		mode_normalize(krylov[0]);
		
		for(i = 1; i <= size; i++){

			temp1 = mode_apply_beta_operator(krylov[i - 1], omega);
			mode_apply_zero_bc(temp1);

			krylov[i] = mode_apply_beta_operator(temp1, omega);
			mode_apply_zero_bc(krylov[i]);
					
			mode_free(temp1);

			for(j = 0; j < i; j++){

				H_index = (i - 1) + j*size;

				inner_prod = mode_inner_product(krylov[j], krylov[i]);
				H.m[H_index] = inner_prod;
				TEMP = matrix_complex_smul(inner_prod, krylov[j].fields);
				matrix_sub_in_place(krylov[i].fields, TEMP);
				matrix_free(TEMP);

			}

			H_index = (i - 1) + i*size;

			norm = mode_norm(krylov[i]);
			H.m[H_index] = complex_from_doubles(norm, 0);
			mode_normalize(krylov[i]);
		}

		D = matrix_diagonalize(H, 1e-3, 100000);
		
		//count positive eigenvalues and obtain ordered list
		mode_count = 0;
		max = 0;
		max_index = 0;
		for(i = 0; i < size; i++){

			eigenvalues[i] = D.m[i + i*size];
			
			if(eigenvalues[i].re > max){
						
				max = eigenvalues[i].re;
				max_index = i;

			}
		}

		reduced_eigenvector = matrix_eigenvector(H, eigenvalues[max_index]);


		beta = complex_sqrt(eigenvalues[max_index]);
		n_eff = complex_smul(1.0/omega, beta);


		//zero out fields for psi
		for(j = 0; j < 6*psi->i*psi->j; j++){

			psi->fields.m[j] = complex_zero();

		}

		//set fields equal to aproximated eigenvector
		for(j = 0; j < size; j++){

			TEMP = matrix_complex_smul(reduced_eigenvector.m[j], krylov[j].fields);			
			matrix_add_in_place(psi->fields,TEMP);
			matrix_free(TEMP);
		}
		matrix_free(reduced_eigenvector);

		printf("\n#%i: beta = %f + i%f\nn_eff = %f + i%f\n\n", m, beta.re, beta.im, n_eff.re, n_eff.im);


		for(i = 0; i < size; i++){
			
			mode_free(krylov[i]);
		}
		matrix_free(D);

		mode_normalize(*psi);
		//mode_apply_H_curl(*psi, beta, omega);
		//mode_normalize(*psi);
		//mode_plot(*psi);

		n_old = n_mag;
		n_mag = complex_mag(n_eff);
		n_diff = n_mag - n_old;
		m++;	
	}

	mode_normalize(*psi);
	
	matrix_free(H);
	free(eigenvalues);
	free(krylov);
	return(beta);		
}

mode mode_apply_beta_squared_operator(mode psi, double omega){

	int i, j;
	int pml_index, pml_index_up, pml_index_down, pml_index_right, pml_index_left;
	int ep_index, ep_index_right, ep_index_left, ep_index_up, ep_index_down, ep_index_down_right, ep_index_up_left;
	int field_index, field_index_right, field_index_left, field_index_up, field_index_down, field_index_down_right, field_index_up_left;
	double divisor;
	complex Hx, Hx_right, Hx_down, Hx_left, Hx_up, Hx_down_right;
	complex Hy, Hy_left, Hy_up, Hy_right, Hy_down, Hy_up_left;
	complex Bx, Bx_right, Bx_left, Bx_down, Bx_down_right;
	complex By, By_up, By_down, By_left, By_up_left;
	complex Bx_dx, Bx_dx_left, Bx_dx_down;
	complex By_dy, By_dy_left, By_dy_down;
	complex Hx_dy, Hx_dy_up, Hx_dy_right;
	complex Hy_dx, Hy_dx_up, Hy_dx_right;
	complex inv_mu_z, inv_mu_z_left, inv_mu_z_down;
	complex inv_ep_z, inv_ep_z_right, inv_ep_z_up;
	complex sx, sx_up, sx_down, sx_right, sx_left;
	complex sy, sy_up, sy_down, sy_right, sy_left;
	complex x_term1, x_term2, x_term3, x_term4, x_term5;
	complex y_term1, y_term2, y_term3, y_term4, y_term5;
	mode beta_psi;
	
	
	beta_psi = mode_new(psi.i, psi.j, psi.disc);
	for(i = 0; i < psi.ep.i; i++){
		
		beta_psi.ep.m[i] = psi.ep.m[i];
		beta_psi.mu.m[i] = psi.mu.m[i];
	}
	for(i = 0; i < 2*psi.i*psi.j; i++){
		
		beta_psi.pml.m[i] = psi.pml.m[i];
	}
	
	divisor = 1/(psi.disc);

	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
			pml_index_up = pml_index + 2*psi.i;
			pml_index_down = pml_index - 2*psi.i;
			pml_index_right = pml_index + 2;
			pml_index_left = pml_index - 2;
		
			ep_index = 3*(i + j*psi.i);
			ep_index_right = ep_index + 3;
			ep_index_left = ep_index - 3;
			ep_index_up = ep_index + 3*psi.i;
			ep_index_down = ep_index - 3*psi.i;
			ep_index_down_right = ep_index - 3*psi.i + 3;
			ep_index_up_left = ep_index + 3*psi.i - 3;
			
			field_index = 6*(i + j*psi.i);
			field_index_right = field_index + 6;
			field_index_left = field_index - 6;
			field_index_up = field_index + 6*psi.i;
			field_index_down = field_index - 6*psi.i;
			field_index_down_right = field_index - 6*psi.i + 6;
			field_index_up_left = field_index + 6*psi.i - 6;
			
			//pick out field elements
			Hx = psi.fields.m[field_index + 3];
			Hx_right = psi.fields.m[field_index_right + 3];
			Hx_left = psi.fields.m[field_index_left + 3];
			Hx_up = psi.fields.m[field_index_up + 3];
			Hx_down = psi.fields.m[field_index_down + 3];
			Hx_down_right = psi.fields.m[field_index_down_right + 3];
			
			Hy = psi.fields.m[field_index + 4];
			Hy_right = psi.fields.m[field_index_right + 4];
			Hy_left = psi.fields.m[field_index_left + 4];
			Hy_up = psi.fields.m[field_index_up + 4];
			Hy_down = psi.fields.m[field_index_down + 4];
			Hy_up_left = psi.fields.m[field_index_up_left + 4];
			
			Bx = complex_mul(psi.mu.m[ep_index], Hx);
			Bx_right = complex_mul(psi.mu.m[ep_index_right], Hx_right);
			Bx_left = complex_mul(psi.mu.m[ep_index_left], Hx_left);
			Bx_down = complex_mul(psi.mu.m[ep_index_down], Hx_down);
			Bx_down_right = complex_mul(psi.mu.m[ep_index_down_right], Hx_down_right);
			
			By = complex_mul(psi.mu.m[ep_index + 1], Hy);
			By_up = complex_mul(psi.mu.m[ep_index_up + 1], Hy_up);
			By_down = complex_mul(psi.mu.m[ep_index_down + 1], Hy_down);
			By_left = complex_mul(psi.mu.m[ep_index_left + 1], Hy_left);
			By_up_left = complex_mul(psi.mu.m[ep_index_up_left + 1], Hy_up_left);

			inv_mu_z = complex_reciprocal(psi.mu.m[ep_index + 2]);
			inv_mu_z_left = complex_reciprocal(psi.mu.m[ep_index_left + 2]);	
			inv_mu_z_down = complex_reciprocal(psi.mu.m[ep_index_down + 2]);
			
			inv_ep_z = complex_reciprocal(psi.ep.m[ep_index + 2]);
			inv_ep_z_right = complex_reciprocal(psi.ep.m[ep_index_right + 2]);	
			inv_ep_z_up = complex_reciprocal(psi.ep.m[ep_index_up + 2]);			
					
			sx = psi.pml.m[pml_index];
			sx_up = psi.pml.m[pml_index_up];
			sx_down = psi.pml.m[pml_index_down];
			sx_right = psi.pml.m[pml_index_right];
			sx_left = psi.pml.m[pml_index_left];
			
			sy = psi.pml.m[pml_index + 1];
			sy_up = psi.pml.m[pml_index_up + 1];
			sy_down = psi.pml.m[pml_index_down + 1];
			sy_right = psi.pml.m[pml_index_right + 1];
			sy_left = psi.pml.m[pml_index_left + 1];
			
			//pick out derivatives as central differences
			Bx_dx = complex_div(complex_smul(divisor, complex_sub(Bx_right, Bx)), sx);
			Bx_dx_left = complex_div(complex_smul(divisor, complex_sub(Bx, Bx_left)), sx_left);
			Bx_dx_down = complex_div(complex_smul(divisor, complex_sub(Bx_down_right, Bx_down)), sx_down);
			
			By_dy = complex_div(complex_smul(divisor, complex_sub(By_up, By)), sy);
			By_dy_down = complex_div(complex_smul(divisor, complex_sub(By, By_down)), sy_down);
			By_dy_left = complex_div(complex_smul(divisor, complex_sub(By_up_left, By_left)), sy_left);
			
			Hx_dy = complex_div(complex_smul(divisor, complex_sub(Hx, Hx_down)), sy);
			Hx_dy_up = complex_div(complex_smul(divisor, complex_sub(Hx_up, Hx)), sy_up);
			Hx_dy_right = complex_div(complex_smul(divisor, complex_sub(Hx_right, Hx_down_right)), sy_right);

			Hy_dx = complex_div(complex_smul(divisor, complex_sub(Hy, Hy_left)), sx);
			Hy_dx_right = complex_div(complex_smul(divisor, complex_sub(Hy_right, Hy)), sx_right);
			Hy_dx_up = complex_div(complex_smul(divisor, complex_sub(Hy_up, Hy_up_left)), sx_up);
			
			//take appropriate second derivatives to construct the terms (again as central differences)
			x_term1 = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Bx_dx), complex_mul(inv_mu_z_left, Bx_dx_left))), sx);	
			x_term2 = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, By_dy), complex_mul(inv_mu_z_left, By_dy_left))), sx);	
			x_term3 = complex_div(complex_mul(complex_smul(-divisor, psi.ep.m[ep_index + 1]), complex_sub(complex_mul(inv_ep_z_up, Hy_dx_up), complex_mul(inv_ep_z, Hy_dx))), sy);
			x_term4 = complex_div(complex_mul(complex_smul(divisor, psi.ep.m[ep_index + 1]), complex_sub(complex_mul(inv_ep_z_up, Hx_dy_up), complex_mul(inv_ep_z, Hx_dy))), sy);
			x_term5 = complex_smul(omega*omega, complex_mul(complex_mul(psi.ep.m[ep_index + 1], psi.mu.m[ep_index]), Hx));
			
			y_term1 = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Bx_dx), complex_mul(inv_mu_z_down, Bx_dx_down))), sy);	
			y_term2 = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, By_dy), complex_mul(inv_mu_z_down, By_dy_down))), sy);	
			y_term3 = complex_div(complex_mul(complex_smul(divisor, psi.ep.m[ep_index]), complex_sub(complex_mul(inv_ep_z_right, Hy_dx_right), complex_mul(inv_ep_z, Hy_dx))), sx);
			y_term4 = complex_div(complex_mul(complex_smul(-divisor, psi.ep.m[ep_index]), complex_sub(complex_mul(inv_ep_z_right, Hx_dy_right), complex_mul(inv_ep_z, Hx_dy))), sx);
			y_term5 = complex_smul(omega*omega, complex_mul(complex_mul(psi.ep.m[ep_index], psi.mu.m[ep_index + 1]), Hy));	

			//insert beta_psi values beta_psi = BETA*psi
			beta_psi.fields.m[field_index + 3] = complex_add(complex_add(complex_add(x_term1, x_term2), complex_add(x_term3, x_term4)), x_term5);
			beta_psi.fields.m[field_index + 4] = complex_add(complex_add(complex_add(y_term1, y_term2), complex_add(y_term3, y_term4)), y_term5);
		}
	}
	
	return(beta_psi);
}

/*mode mode_apply_radial_beta_squared_operator(mode psi, double omega, double radius_left){

	int i, j;
	int pml_index;
	int ep_index, ep_index_right, ep_index_left, ep_index_up, ep_index_down, ep_index_down_right, ep_index_up_left;
	int field_index, field_index_right, field_index_left, field_index_up, field_index_down, field_index_down_right, field_index_up_left;
	int field_index_up_left_left, field_index_left_left, field_index_down_left, field_index_down_down;
	double radius, divisor;
	complex Hx, Hx_right, Hx_down, Hx_left, Hx_up, Hx_down_right, Hx_up_left, Hx_down_left, Hx_down_down;
	complex Hy, Hy_left, Hy_up, Hy_right, Hy_down, Hy_up_left, Hy_up_left_left, Hy_left_left, Hy_down_left;
	complex Bx, Bx_right, Bx_left, Bx_down, Bx_down_right;
	complex By, By_up, By_down, By_left, By_up_left;
	complex Bx_dx, Bx_dx_left, Bx_dx_down;
	complex By_dy, By_dy_left, By_dy_down;
	complex Hx_dy, Hx_dy_up, Hx_dy_right, Hx_dy_left, Hx_dy_up_left, Hx_dy_down;
	complex Hy_dx, Hy_dx_up, Hy_dx_right, Hy_dx_up_left, Hy_dx_left, Hy_dx_down;
	complex Hy_dx_dy, Hy_dx_dy_left, Hy_dx_dy_down;
	complex Hx_dy_dy, Hx_dy_dy_left, Hx_dy_dy_down;
	complex inv_mu_z, inv_mu_z_left, inv_mu_z_down;
	complex inv_ep_z, inv_ep_z_right, inv_ep_z_up, inv_ep_z_left, inv_ep_z_up_left, inv_ep_z_down;
	complex sx, sy;
	complex x_term1, x_term2, x_term3, x_term4, x_term5, x_term6, x_term7, x_term8, x_term9, x_term10;
	complex y_term1, y_term2, y_term3, y_term4, y_term5, y_term6, y_term7, y_term8;
	mode beta_psi;
	
	
	beta_psi = mode_new(psi.i, psi.j, psi.disc);
	for(i = 0; i < psi.ep.i; i++){
		
		beta_psi.ep.m[i] = psi.ep.m[i];
		beta_psi.mu.m[i] = psi.mu.m[i];
	}
	for(i = 0; i < 2*psi.i*psi.j; i++){
		
		beta_psi.pml.m[i] = psi.pml.m[i];
	}
	
	divisor = 1/(psi.disc);

	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		
		radius = radius_left + i*psi.disc;
		
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
		
			ep_index = 3*(i + j*psi.i);
			ep_index_right = ep_index + 3;
			ep_index_left = ep_index - 3;
			ep_index_up = ep_index + 3*psi.i;
			ep_index_down = ep_index - 3*psi.i;
			ep_index_down_right = ep_index - 3*psi.i + 3;
			ep_index_up_left = ep_index + 3*psi.i - 3;
			
			field_index = 6*(i + j*psi.i);
			field_index_right = field_index + 6;
			field_index_left = field_index - 6;
			field_index_up = field_index + 6*psi.i;
			field_index_down = field_index - 6*psi.i;
			field_index_down_right = field_index - 6*psi.i + 6;
			field_index_up_left = field_index + 6*psi.i - 6;
			field_index_down_left = field_index - 6*psi.i - 6;
			field_index_up_left_left = field_index + 6*psi.i - 12;
			field_index_left_left = field_index - 12;
			field_index_down_down = field_index - 12*psi.i;
			
			//pick out field elements
			Hx = psi.fields.m[field_index + 3];
			Hx_right = psi.fields.m[field_index_right + 3];
			Hx_left = psi.fields.m[field_index_left + 3];
			Hx_up = psi.fields.m[field_index_up + 3];
			Hx_down = psi.fields.m[field_index_down + 3];
			Hx_down_right = psi.fields.m[field_index_down_right + 3];
			Hx_up_left = psi.fields.m[field_index_up_left + 3];
			Hx_down_left = psi.fields.m[field_index_down_left + 3];
			
			if(j >= 2){
				
				Hx_down_down = psi.fields.m[field_index_down_down + 3];
			}
			else{
			
				Hx_down_down = complex_zero();
			}
			
			Hy = psi.fields.m[field_index + 4];
			Hy_right = psi.fields.m[field_index_right + 4];
			Hy_left = psi.fields.m[field_index_left + 4];
			Hy_up = psi.fields.m[field_index_up + 4];
			Hy_down = psi.fields.m[field_index_down + 4];
			Hy_down_left = psi.fields.m[field_index_down_left + 4];
			Hy_up_left = psi.fields.m[field_index_up_left + 4];
			
			if(i >= 2){
			
				Hy_up_left_left = psi.fields.m[field_index_up_left_left + 4];
				Hy_left_left = psi.fields.m[field_index_left_left + 4];
			}
			else{
				
				Hy_up_left_left = complex_zero();
				Hy_left_left = complex_zero();
			}

			Bx = complex_mul(psi.mu.m[ep_index], Hx);
			Bx_right = complex_mul(psi.mu.m[ep_index_right], Hx_right);
			Bx_left = complex_mul(psi.mu.m[ep_index_left], Hx_left);
			Bx_down = complex_mul(psi.mu.m[ep_index_down], Hx_down);
			Bx_down_right = complex_mul(psi.mu.m[ep_index_down_right], Hx_down_right);
			
			By = complex_mul(psi.mu.m[ep_index + 1], Hy);
			By_up = complex_mul(psi.mu.m[ep_index_up + 1], Hy_up);
			By_down = complex_mul(psi.mu.m[ep_index_down + 1], Hy_down);
			By_left = complex_mul(psi.mu.m[ep_index_left + 1], Hy_left);
			By_up_left = complex_mul(psi.mu.m[ep_index_up_left + 1], Hy_up_left);

			inv_mu_z = complex_reciprocal(psi.mu.m[ep_index + 2]);
			inv_mu_z_left = complex_reciprocal(psi.mu.m[ep_index_left + 2]);	
			inv_mu_z_down = complex_reciprocal(psi.mu.m[ep_index_down + 2]);
			
			inv_ep_z = complex_reciprocal(psi.ep.m[ep_index + 2]);
			inv_ep_z_right = complex_reciprocal(psi.ep.m[ep_index_right + 2]);	
			inv_ep_z_up = complex_reciprocal(psi.ep.m[ep_index_up + 2]);
			inv_ep_z_left = complex_reciprocal(psi.ep.m[ep_index_left + 2]);
			inv_ep_z_down = complex_reciprocal(psi.ep.m[ep_index_down + 2]);
			inv_ep_z_up_left = complex_reciprocal(psi.ep.m[ep_index_up_left + 2]);
							
			sx = psi.pml.m[pml_index];
			sy = psi.pml.m[pml_index + 1];
			
			//pick out derivatives as central differences
			Bx_dx = complex_div(complex_smul(divisor, complex_sub(Bx_right, Bx)), sx);
			Bx_dx_left = complex_div(complex_smul(divisor, complex_sub(Bx, Bx_left)), sx);
			Bx_dx_down = complex_div(complex_smul(divisor, complex_sub(Bx_down_right, Bx_down)), sx);
			
			By_dy = complex_div(complex_smul(divisor, complex_sub(By_up, By)), sy);
			By_dy_down = complex_div(complex_smul(divisor, complex_sub(By, By_down)), sy);
			By_dy_left = complex_div(complex_smul(divisor, complex_sub(By_up_left, By_left)), sy);
			
			Hx_dy = complex_div(complex_smul(divisor, complex_sub(Hx, Hx_down)), sy);
			Hx_dy_up = complex_div(complex_smul(divisor, complex_sub(Hx_up, Hx)), sy);
			Hx_dy_down = complex_div(complex_smul(divisor, complex_sub(Hx_down, Hx_down_down)), sy);
			Hx_dy_left = complex_div(complex_smul(divisor, complex_sub(Hx_left, Hx_down_left)), sy);
			Hx_dy_up_left = complex_div(complex_smul(divisor, complex_sub(Hx_up_left, Hx_left)), sy);
			Hx_dy_right = complex_div(complex_smul(divisor, complex_sub(Hx_right, Hx_down_right)), sy);

			Hy_dx = complex_div(complex_smul(divisor, complex_sub(Hy, Hy_left)), sx);
			Hy_dx_right = complex_div(complex_smul(divisor, complex_sub(Hy_right, Hy)), sx);
			Hy_dx_up = complex_div(complex_smul(divisor, complex_sub(Hy_up, Hy_up_left)), sx);
			Hy_dx_down = complex_div(complex_smul(divisor, complex_sub(Hy_down, Hy_down_left)), sx);
			Hy_dx_up_left = complex_div(complex_smul(divisor, complex_sub(Hy_up_left, Hy_up_left_left)), sx);
			Hy_dx_left = complex_div(complex_smul(divisor, complex_sub(Hy_left, Hy_left_left)), sx);
			
			//second_derivatives
			Hy_dx_dy = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z_up, Hy_dx_up), complex_mul(inv_ep_z, Hy_dx))), sy);
			Hy_dx_dy_left = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z_up_left, Hy_dx_up_left), complex_mul(inv_ep_z_left, Hy_dx_left))), sy);
			Hy_dx_dy_down = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z, Hy_dx), complex_mul(inv_ep_z_down, Hy_dx_down))), sy);
			
			Hx_dy_dy = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z_up, Hx_dy_up), complex_mul(inv_ep_z, Hx_dy))), sy);
			Hx_dy_dy_left = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z_up_left, Hx_dy_up_left), complex_mul(inv_ep_z_left, Hx_dy_left))), sy);
			Hx_dy_dy_down = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z, Hx_dy), complex_mul(inv_ep_z_down, Hx_dy_down))), sy);
			
			//take appropriate second derivatives to construct the terms (again as central differences)
			//r_squared terms
			x_term1 = complex_smul(radius*radius, complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Bx_dx), complex_mul(inv_mu_z_left, Bx_dx_left))), sx));	
			x_term2 = complex_smul(radius*radius, complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, By_dy), complex_mul(inv_mu_z_left, By_dy_left))), sx));	
			x_term3 = complex_smul(radius*radius, complex_div(complex_mul(complex_smul(-divisor, psi.ep.m[ep_index + 1]), complex_sub(complex_mul(inv_ep_z_up, Hy_dx_up), complex_mul(inv_ep_z, Hy_dx))), sy));
			x_term4 = complex_smul(radius*radius, complex_div(complex_mul(complex_smul(divisor, psi.ep.m[ep_index + 1]), complex_sub(complex_mul(inv_ep_z_up, Hx_dy_up), complex_mul(inv_ep_z, Hx_dy))), sy));
			x_term5 = complex_smul(radius*radius*omega*omega, complex_mul(complex_mul(psi.ep.m[ep_index + 1], psi.mu.m[ep_index]), Hx));
			
			//r terms
			x_term6 = complex_smul(radius, complex_mul(inv_mu_z, Bx_dx));
			x_term7 = complex_smul(radius, complex_div(complex_smul(divisor, complex_sub(Hx_right, Hx)), sx));
			x_term8 = complex_smul(radius, complex_mul(inv_mu_z, By_dy));
			x_term9 = complex_smul(radius/(omega*omega), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z_left, Hy_dx_dy_left), complex_mul(inv_mu_z, Hy_dx_dy))), sx));
			x_term10 = complex_smul(radius/(omega*omega), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Hx_dy_dy), complex_mul(inv_mu_z_left, Hx_dy_dy_left))), sx));
			
			//r_squared terms
			y_term1 = complex_smul(radius*radius, complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Bx_dx), complex_mul(inv_mu_z_down, Bx_dx_down))), sy));	
			y_term2 = complex_smul(radius*radius, complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, By_dy), complex_mul(inv_mu_z_down, By_dy_down))), sy));	
			y_term3 = complex_smul(radius*radius, complex_div(complex_mul(complex_smul(divisor, psi.ep.m[ep_index]), complex_sub(complex_mul(inv_ep_z_right, Hy_dx_right), complex_mul(inv_ep_z, Hy_dx))), sx));
			y_term4 = complex_smul(radius*radius, complex_div(complex_mul(complex_smul(-divisor, psi.ep.m[ep_index]), complex_sub(complex_mul(inv_ep_z_right, Hx_dy_right), complex_mul(inv_ep_z, Hx_dy))), sx));
			y_term5 = complex_smul(radius*radius*omega*omega, complex_mul(complex_mul(psi.ep.m[ep_index], psi.mu.m[ep_index + 1]), Hy));	
			
			//r_terms
			y_term6 = complex_smul(radius, Hx_dy);
			y_term7 = complex_smul(radius/(omega*omega), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Hx_dy_dy), complex_mul(inv_mu_z_down, Hx_dy_dy_down))), sy));
			y_term8 = complex_smul(radius/(omega*omega), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z_down, Hy_dx_dy_down), complex_mul(inv_mu_z, Hy_dx_dy))), sy));
			
			//insert beta_psi values beta_psi = BETA*psi
			beta_psi.fields.m[field_index + 3] = complex_add(complex_add(complex_add(complex_add(x_term1, x_term2), complex_add(x_term3, x_term4)), complex_add(complex_add(x_term5, x_term6), complex_add(x_term7, x_term8))), complex_add(x_term9, x_term10));
			beta_psi.fields.m[field_index + 4] = complex_add(complex_add(complex_add(y_term1, y_term2), complex_add(y_term3, y_term4)), complex_add(complex_add(y_term5, y_term6), complex_add(y_term7, y_term8)));
		}
	}
	
	return(beta_psi);
}
*/

mode mode_apply_radial_beta_squared_operator(mode psi, double omega, double radius_left){

	int i, j, k;
	int pml_index;
	int ep_index, ep_index_right, ep_index_left, ep_index_up, ep_index_down, ep_index_down_right, ep_index_up_left;
	int field_index, field_index_right, field_index_left, field_index_up, field_index_down, field_index_down_right, field_index_up_left;
	int field_index_up_left_left, field_index_left_left, field_index_down_left, field_index_down_down;
	double divisor;
	complex Hx, Hx_right, Hx_down, Hx_left, Hx_up, Hx_down_right, Hx_up_left, Hx_down_left, Hx_down_down;
	complex Hy, Hy_left, Hy_up, Hy_right, Hy_down, Hy_up_left, Hy_up_left_left, Hy_left_left, Hy_down_left;
	complex Bx, Bx_right, Bx_left, Bx_down, Bx_down_right;
	complex By, By_up, By_down, By_left, By_up_left;
	complex Bx_dx, Bx_dx_left, Bx_dx_down;
	complex By_dy, By_dy_left, By_dy_down;
	complex Hx_dy, Hx_dy_up, Hx_dy_right, Hx_dy_left, Hx_dy_up_left, Hx_dy_down;
	complex Hy_dx, Hy_dx_up, Hy_dx_right, Hy_dx_up_left, Hy_dx_left, Hy_dx_down;
	complex Hy_dx_dy, Hy_dx_dy_left, Hy_dx_dy_down;
	complex Hx_dy_dy, Hx_dy_dy_left, Hx_dy_dy_down;
	complex inv_mu_z, inv_mu_z_left, inv_mu_z_down;
	complex inv_ep_z, inv_ep_z_right, inv_ep_z_up, inv_ep_z_left, inv_ep_z_up_left, inv_ep_z_down;
	complex rad, sx, sy;
	complex x_term1, x_term2, x_term3, x_term4, x_term5, x_term6, x_term7, x_term8, x_term9, x_term10;
	complex y_term1, y_term2, y_term3, y_term4, y_term5, y_term6, y_term7, y_term8;
	mode beta_psi;
	
	
	beta_psi = mode_new(psi.i, psi.j, psi.disc);
	for(i = 0; i < psi.ep.i; i++){
		
		beta_psi.ep.m[i] = psi.ep.m[i];
		beta_psi.mu.m[i] = psi.mu.m[i];
	}
	for(i = 0; i < 2*psi.i*psi.j; i++){
		
		beta_psi.pml.m[i] = psi.pml.m[i];
	}
	
	divisor = 1/(psi.disc);

	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
		
			ep_index = 3*(i + j*psi.i);
			ep_index_right = ep_index + 3;
			ep_index_left = ep_index - 3;
			ep_index_up = ep_index + 3*psi.i;
			ep_index_down = ep_index - 3*psi.i;
			ep_index_down_right = ep_index - 3*psi.i + 3;
			ep_index_up_left = ep_index + 3*psi.i - 3;
			
			field_index = 6*(i + j*psi.i);
			field_index_right = field_index + 6;
			field_index_left = field_index - 6;
			field_index_up = field_index + 6*psi.i;
			field_index_down = field_index - 6*psi.i;
			field_index_down_right = field_index - 6*psi.i + 6;
			field_index_up_left = field_index + 6*psi.i - 6;
			field_index_down_left = field_index - 6*psi.i - 6;
			field_index_up_left_left = field_index + 6*psi.i - 12;
			field_index_left_left = field_index - 12;
			field_index_down_down = field_index - 12*psi.i;
			
			//pick out field elements
			Hx = psi.fields.m[field_index + 3];
			Hx_right = psi.fields.m[field_index_right + 3];
			Hx_left = psi.fields.m[field_index_left + 3];
			Hx_up = psi.fields.m[field_index_up + 3];
			Hx_down = psi.fields.m[field_index_down + 3];
			Hx_down_right = psi.fields.m[field_index_down_right + 3];
			Hx_up_left = psi.fields.m[field_index_up_left + 3];
			Hx_down_left = psi.fields.m[field_index_down_left + 3];
			
			if(j >= 2){
				
				Hx_down_down = psi.fields.m[field_index_down_down + 3];
			}
			else{
			
				Hx_down_down = complex_zero();
			}
			
			Hy = psi.fields.m[field_index + 4];
			Hy_right = psi.fields.m[field_index_right + 4];
			Hy_left = psi.fields.m[field_index_left + 4];
			Hy_up = psi.fields.m[field_index_up + 4];
			Hy_down = psi.fields.m[field_index_down + 4];
			Hy_down_left = psi.fields.m[field_index_down_left + 4];
			Hy_up_left = psi.fields.m[field_index_up_left + 4];
			
			if(i >= 2){
			
				Hy_up_left_left = psi.fields.m[field_index_up_left_left + 4];
				Hy_left_left = psi.fields.m[field_index_left_left + 4];
			}
			else{
				
				Hy_up_left_left = complex_zero();
				Hy_left_left = complex_zero();
			}

			Bx = complex_mul(psi.mu.m[ep_index], Hx);
			Bx_right = complex_mul(psi.mu.m[ep_index_right], Hx_right);
			Bx_left = complex_mul(psi.mu.m[ep_index_left], Hx_left);
			Bx_down = complex_mul(psi.mu.m[ep_index_down], Hx_down);
			Bx_down_right = complex_mul(psi.mu.m[ep_index_down_right], Hx_down_right);
			
			By = complex_mul(psi.mu.m[ep_index + 1], Hy);
			By_up = complex_mul(psi.mu.m[ep_index_up + 1], Hy_up);
			By_down = complex_mul(psi.mu.m[ep_index_down + 1], Hy_down);
			By_left = complex_mul(psi.mu.m[ep_index_left + 1], Hy_left);
			By_up_left = complex_mul(psi.mu.m[ep_index_up_left + 1], Hy_up_left);

			inv_mu_z = complex_reciprocal(psi.mu.m[ep_index + 2]);
			inv_mu_z_left = complex_reciprocal(psi.mu.m[ep_index_left + 2]);	
			inv_mu_z_down = complex_reciprocal(psi.mu.m[ep_index_down + 2]);
			
			inv_ep_z = complex_reciprocal(psi.ep.m[ep_index + 2]);
			inv_ep_z_right = complex_reciprocal(psi.ep.m[ep_index_right + 2]);	
			inv_ep_z_up = complex_reciprocal(psi.ep.m[ep_index_up + 2]);
			inv_ep_z_left = complex_reciprocal(psi.ep.m[ep_index_left + 2]);
			inv_ep_z_down = complex_reciprocal(psi.ep.m[ep_index_down + 2]);
			inv_ep_z_up_left = complex_reciprocal(psi.ep.m[ep_index_up_left + 2]);
							
			sx = psi.pml.m[pml_index];
			sy = psi.pml.m[pml_index + 1];
			
			//radius in stretched coordinates
			rad = complex_from_doubles(radius_left, 0);
			for(k = 1; k <= i; k++){
				
				rad = complex_add(rad, complex_smul(psi.disc, sx));
			}
			
			//pick out derivatives as central differences
			Bx_dx = complex_div(complex_smul(divisor, complex_sub(Bx_right, Bx)), sx);
			Bx_dx_left = complex_div(complex_smul(divisor, complex_sub(Bx, Bx_left)), sx);
			Bx_dx_down = complex_div(complex_smul(divisor, complex_sub(Bx_down_right, Bx_down)), sx);
			
			By_dy = complex_div(complex_smul(divisor, complex_sub(By_up, By)), sy);
			By_dy_down = complex_div(complex_smul(divisor, complex_sub(By, By_down)), sy);
			By_dy_left = complex_div(complex_smul(divisor, complex_sub(By_up_left, By_left)), sy);
			
			Hx_dy = complex_div(complex_smul(divisor, complex_sub(Hx, Hx_down)), sy);
			Hx_dy_up = complex_div(complex_smul(divisor, complex_sub(Hx_up, Hx)), sy);
			Hx_dy_down = complex_div(complex_smul(divisor, complex_sub(Hx_down, Hx_down_down)), sy);
			Hx_dy_left = complex_div(complex_smul(divisor, complex_sub(Hx_left, Hx_down_left)), sy);
			Hx_dy_up_left = complex_div(complex_smul(divisor, complex_sub(Hx_up_left, Hx_left)), sy);
			Hx_dy_right = complex_div(complex_smul(divisor, complex_sub(Hx_right, Hx_down_right)), sy);

			Hy_dx = complex_div(complex_smul(divisor, complex_sub(Hy, Hy_left)), sx);
			Hy_dx_right = complex_div(complex_smul(divisor, complex_sub(Hy_right, Hy)), sx);
			Hy_dx_up = complex_div(complex_smul(divisor, complex_sub(Hy_up, Hy_up_left)), sx);
			Hy_dx_down = complex_div(complex_smul(divisor, complex_sub(Hy_down, Hy_down_left)), sx);
			Hy_dx_up_left = complex_div(complex_smul(divisor, complex_sub(Hy_up_left, Hy_up_left_left)), sx);
			Hy_dx_left = complex_div(complex_smul(divisor, complex_sub(Hy_left, Hy_left_left)), sx);
			
			//second_derivatives
			Hy_dx_dy = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z_up, Hy_dx_up), complex_mul(inv_ep_z, Hy_dx))), sy);
			Hy_dx_dy_left = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z_up_left, Hy_dx_up_left), complex_mul(inv_ep_z_left, Hy_dx_left))), sy);
			Hy_dx_dy_down = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z, Hy_dx), complex_mul(inv_ep_z_down, Hy_dx_down))), sy);
			
			Hx_dy_dy = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z_up, Hx_dy_up), complex_mul(inv_ep_z, Hx_dy))), sy);
			Hx_dy_dy_left = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z_up_left, Hx_dy_up_left), complex_mul(inv_ep_z_left, Hx_dy_left))), sy);
			Hx_dy_dy_down = complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_ep_z, Hx_dy), complex_mul(inv_ep_z_down, Hx_dy_down))), sy);
			
			//take appropriate second derivatives to construct the terms (again as central differences)
			//r_squared terms
			x_term1 = complex_mul(complex_mul(rad, rad), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Bx_dx), complex_mul(inv_mu_z_left, Bx_dx_left))), sx));	
			x_term2 = complex_mul(complex_mul(rad, rad), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, By_dy), complex_mul(inv_mu_z_left, By_dy_left))), sx));	
			x_term3 = complex_mul(complex_mul(rad, rad), complex_div(complex_mul(complex_smul(-divisor, psi.ep.m[ep_index + 1]), complex_sub(complex_mul(inv_ep_z_up, Hy_dx_up), complex_mul(inv_ep_z, Hy_dx))), sy));
			x_term4 = complex_mul(complex_mul(rad, rad), complex_div(complex_mul(complex_smul(divisor, psi.ep.m[ep_index + 1]), complex_sub(complex_mul(inv_ep_z_up, Hx_dy_up), complex_mul(inv_ep_z, Hx_dy))), sy));
			x_term5 = complex_mul(complex_smul(omega*omega, complex_mul(rad, rad)), complex_mul(complex_mul(psi.ep.m[ep_index + 1], psi.mu.m[ep_index]), Hx));
			
			//r terms
			x_term6 = complex_mul(rad, complex_mul(inv_mu_z, Bx_dx));
			x_term7 = complex_mul(rad, complex_div(complex_smul(divisor, complex_sub(Hx_right, Hx)), sx));
			x_term8 = complex_mul(rad, complex_mul(inv_mu_z, By_dy));
			x_term9 = complex_mul(complex_smul(1/(omega*omega),rad), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z_left, Hy_dx_dy_left), complex_mul(inv_mu_z, Hy_dx_dy))), sx));
			x_term10 = complex_mul(complex_smul(1/(omega*omega),rad), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Hx_dy_dy), complex_mul(inv_mu_z_left, Hx_dy_dy_left))), sx));
			
			//r_squared terms
			y_term1 = complex_mul(complex_mul(rad, rad), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Bx_dx), complex_mul(inv_mu_z_down, Bx_dx_down))), sy));	
			y_term2 = complex_mul(complex_mul(rad, rad), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, By_dy), complex_mul(inv_mu_z_down, By_dy_down))), sy));	
			y_term3 = complex_mul(complex_mul(rad, rad), complex_div(complex_mul(complex_smul(divisor, psi.ep.m[ep_index]), complex_sub(complex_mul(inv_ep_z_right, Hy_dx_right), complex_mul(inv_ep_z, Hy_dx))), sx));
			y_term4 = complex_mul(complex_mul(rad, rad), complex_div(complex_mul(complex_smul(-divisor, psi.ep.m[ep_index]), complex_sub(complex_mul(inv_ep_z_right, Hx_dy_right), complex_mul(inv_ep_z, Hx_dy))), sx));
			y_term5 = complex_mul(complex_smul(omega*omega, complex_mul(rad, rad)), complex_mul(complex_mul(psi.ep.m[ep_index], psi.mu.m[ep_index + 1]), Hy));	
			
			//r_terms
			y_term6 = complex_mul(rad, Hx_dy);
			y_term7 = complex_mul(complex_smul(1/(omega*omega),rad), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z, Hx_dy_dy), complex_mul(inv_mu_z_down, Hx_dy_dy_down))), sy));
			y_term8 = complex_mul(complex_smul(1/(omega*omega),rad), complex_div(complex_smul(divisor, complex_sub(complex_mul(inv_mu_z_down, Hy_dx_dy_down), complex_mul(inv_mu_z, Hy_dx_dy))), sy));
			
			//insert beta_psi values beta_psi = BETA*psi
			beta_psi.fields.m[field_index + 3] = complex_add(complex_add(complex_add(complex_add(x_term1, x_term2), complex_add(x_term3, x_term4)), complex_add(complex_add(x_term5, x_term6), complex_add(x_term7, x_term8))), complex_add(x_term9, x_term10));
			beta_psi.fields.m[field_index + 4] = complex_add(complex_add(complex_add(y_term1, y_term2), complex_add(y_term3, y_term4)), complex_add(complex_add(y_term5, y_term6), complex_add(y_term7, y_term8)));
		}
	}
	
	return(beta_psi);
}

void mode_apply_H_div(mode psi, complex beta){
	
	int i, j;
	int ep_index, ep_index_up, ep_index_right;
	int field_index, field_index_right, field_index_up;
	int pml_index;
	double divisor;
	complex Hx, Hy, Hx_right, Hy_up;
	complex Bx, By, Bx_right, By_up;
	complex Bx_dx, By_dy;
	complex sx, sy;
	
	divisor = 1/(psi.disc);
	
	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
		
			ep_index = 3*(i + j*psi.i);
			ep_index_right = ep_index + 3;
			ep_index_up = ep_index + 3*psi.i;
			
			field_index = 2*ep_index;
			field_index_right = field_index + 6;
			field_index_up = field_index + 6*psi.i;
			
			//pick out field elements
			Hx = psi.fields.m[field_index + 3];
			Hy = psi.fields.m[field_index + 4];
			
			Hx_right = psi.fields.m[field_index_right + 3];
			Hy_up = psi.fields.m[field_index_up + 4];
			
			Bx = complex_mul(psi.mu.m[ep_index], Hx);
			By = complex_mul(psi.mu.m[ep_index + 1], Hy);
			
			Bx_right = complex_mul(psi.mu.m[ep_index_right], Hx_right);
			By_up = complex_mul(psi.mu.m[ep_index_up + 1], Hy_up);
			
			sx = psi.pml.m[pml_index];
			sy = psi.pml.m[pml_index + 1];
			
			//pick out derivatives as central differences
			Bx_dx = complex_div(complex_smul(divisor, complex_sub(Bx_right, Bx)), sx);
			By_dy = complex_div(complex_smul(divisor, complex_sub(By_up, By)), sy);
			
			//insert psi Hz value  as Hz = -i/beta*(Bx_dx + By_dy)
			psi.fields.m[field_index + 5] = complex_mul(complex_div(complex_from_doubles(0, -1), beta), complex_add(Bx_dx, By_dy));
		}
	}
}

void mode_apply_radial_H_div(mode psi, complex beta, double radius_left){
	
	int i, j;
	int ep_index, ep_index_up, ep_index_right;
	int field_index, field_index_right, field_index_up;
	int pml_index;
	double radius, divisor;
	complex Hx, Hy, Hx_right, Hy_up;
	complex Bx, By, Bx_right, By_up;
	complex Bx_dx, By_dy;
	complex sx, sy;
	
	divisor = 1/(psi.disc);
	
	//move through entire mode vector
	for(i = 1; i < (psi.i - 1); i++){
		
		radius = radius_left + i*psi.disc;
		
		for(j = 1; j < (psi.j - 1); j++){
			
			pml_index = 2*(i + j*psi.i);
		
			ep_index = 3*(i + j*psi.i);
			ep_index_right = ep_index + 3;
			ep_index_up = ep_index + 3*psi.i;
			
			field_index = 2*ep_index;
			field_index_right = field_index + 6;
			field_index_up = field_index + 6*psi.i;
			
			//pick out field elements
			Hx = psi.fields.m[field_index + 3];
			Hy = psi.fields.m[field_index + 4];
			
			Hx_right = psi.fields.m[field_index_right + 3];
			Hy_up = psi.fields.m[field_index_up + 4];
			
			Bx = complex_mul(psi.mu.m[ep_index], Hx);
			By = complex_mul(psi.mu.m[ep_index + 1], Hy);
			
			Bx_right = complex_mul(psi.mu.m[ep_index_right], Hx_right);
			By_up = complex_mul(psi.mu.m[ep_index_up + 1], Hy_up);
			
			sx = psi.pml.m[pml_index];
			sy = psi.pml.m[pml_index + 1];
			
			//pick out derivatives as central differences
			Bx_dx = complex_div(complex_smul(divisor, complex_sub(Bx_right, Bx)), sx);
			By_dy = complex_div(complex_smul(divisor, complex_sub(By_up, By)), sy);
			
			//insert psi Hz value  as Hz = -i/beta*(Bx_dx + By_dy)
			psi.fields.m[field_index + 5] = complex_mul(complex_div(complex_from_doubles(0, -radius), beta), complex_add(Bx_dx, By_dy));
		}
	}
}

complex mode_solve_beta_squared(mode *psi, double omega, double tol, mode *prev_modes, int prev_num, int size, int confinement, int lateral_symmetry){
	
	int i, j, k, m, spurs, H_index, mode_count, orthogonalizations;
	int max_index;
	double beta_mag, n_mag, n_old, n_diff;
	double max_real_n, overlap, norm, max_2;
	double max;
	complex *eigenvalues;
	complex beta, n_eff, inner_prod;
	mode temp1;
	mode *krylov;
	matrix H, D;
	matrix reduced_eigenvector;
	matrix TEMP;
	
	orthogonalizations = 1;
	
	eigenvalues = malloc(size*sizeof(complex));
	krylov = malloc((size + 1)*sizeof(mode));

	H = matrix_new(size, size);
	
	//mode_initialize_general(*psi);
	mode_initialize_rand(*psi);
	mode_H_normalize(*psi);
	
	n_diff = 100;
	n_mag = 100;
	
	m = 0;
	mode_count = 1;
	while(fabs(n_diff) > tol){

		//orthogonalize psi
		for(j = 0; j < orthogonalizations; j++){

			mode_H_normalize(*psi);
			for(i = prev_num - 1; i >= 0; i--){

				inner_prod = mode_H_inner_product(prev_modes[i], *psi);
				TEMP = matrix_complex_smul(inner_prod, prev_modes[i].fields);
				matrix_sub_in_place(psi->fields, TEMP);
				mode_H_normalize(*psi);
				matrix_free(TEMP);
			}

			mode_apply_H_curl(*psi, beta, omega);
			mode_apply_bc(*psi, confinement);
		}
		
		//symmetrize
		if(lateral_symmetry){

			mode_apply_lateral_symmetry(*psi);
		}
		
		//zero out fields in PML region
		mode_apply_zero_bc_pml(*psi);
		
		//build Arnoldi matrix
		krylov[0] = mode_copy(*psi);
		mode_H_normalize(krylov[0]);
		mode_apply_bc(krylov[0], confinement);
		
		for(i = 1; i <= size; i++){

			krylov[i] = mode_apply_beta_squared_operator(krylov[i - 1], omega);
			mode_apply_bc(krylov[i], confinement);

			for(j = 0; j < i; j++){

				H_index = (i - 1) + j*size;

				inner_prod = mode_H_inner_product(krylov[j], krylov[i]);
				H.m[H_index] = inner_prod;
				TEMP = matrix_complex_smul(inner_prod, krylov[j].fields);
				matrix_sub_in_place(krylov[i].fields, TEMP);
				matrix_free(TEMP);

			}

			H_index = (i - 1) + i*size;

			norm = mode_H_norm(krylov[i]);
			H.m[H_index] = complex_from_doubles(norm, 0);
			mode_H_normalize(krylov[i]);
		}

		D = matrix_diagonalize(H, 1e-3, 1000);
		
		//count positive eigenvalues and obtain ordered list
		mode_count = 0;
		max = 0;
		max_index = 0;
		for(i = 0; i < size; i++){

			eigenvalues[i] = D.m[i + i*size];
			
			if(eigenvalues[i].re > max){
						
				max = eigenvalues[i].re;
				max_index = i;

			}
		}

		reduced_eigenvector = matrix_eigenvector(H, eigenvalues[max_index]);


		beta = complex_sqrt(eigenvalues[max_index]);
		n_eff = complex_smul(1.0/omega, beta);


		//zero out fields for psi
		for(j = 0; j < 6*psi->i*psi->j; j++){

			psi->fields.m[j] = complex_zero();

		}

		//set fields equal to aproximated eigenvector
		for(j = 0; j < size; j++){

			TEMP = matrix_complex_smul(reduced_eigenvector.m[j], krylov[j].fields);			
			matrix_add_in_place(psi->fields,TEMP);
			matrix_free(TEMP);
		}
		matrix_free(reduced_eigenvector);

		printf("\n#%i: beta = %f + i%f\nn_eff = %f + i%f\n\n", m, beta.re, beta.im, n_eff.re, n_eff.im);


		for(i = 0; i <= size; i++){
			
			mode_free(krylov[i]);
		}
		matrix_free(D);

		n_old = n_mag;
		n_mag = complex_mag(n_eff);
		n_diff = n_mag - n_old;
		m++;	
	}
	
	//mode_apply_zero_bc(*psi);
	mode_apply_H_div(*psi, beta);
	mode_apply_bc(*psi, confinement);
	mode_apply_H_curl(*psi, beta, omega);
	mode_apply_bc(*psi, confinement);
	mode_H_normalize(*psi);
	
	matrix_free(H);
	free(eigenvalues);
	free(krylov);
	return(beta);		
}

complex mode_solve_radial_beta_squared(mode *psi, double omega, double tol, mode *prev_modes, int prev_num, int size, int confinement, double radius_left){
	
	int i, j, k, m, spurs, H_index, mode_count, orthogonalizations;
	int max_index;
	double beta_mag, n_mag, n_old, n_diff;
	double max_real_n, overlap, norm, max_2;
	double max;
	complex *eigenvalues;
	complex beta, n_eff, inner_prod;
	mode temp1;
	mode *krylov;
	matrix H, D;
	matrix reduced_eigenvector;
	matrix TEMP;
	
	orthogonalizations = 1;
	
	eigenvalues = malloc(size*sizeof(complex));
	krylov = malloc((size + 1)*sizeof(mode));

	H = matrix_new(size, size);
	
	mode_initialize_general(*psi);
	mode_H_normalize(*psi);
	
	n_diff = 100;
	n_mag = 100;
	
	m = 0;
	mode_count = 1;
	while(fabs(n_diff) > tol){
		

		//orthogonalize psi
		for(j = 0; j < orthogonalizations; j++){

			mode_H_normalize(*psi);
			for(i = prev_num - 1; i >= 0; i--){

				inner_prod = mode_H_inner_product(prev_modes[i], *psi);
				TEMP = matrix_complex_smul(inner_prod, prev_modes[i].fields);
				matrix_sub_in_place(psi->fields, TEMP);
				mode_H_normalize(*psi);
				matrix_free(TEMP);
			}

			mode_apply_H_curl(*psi, beta, omega);
			mode_apply_bc(*psi, confinement);
		}
		
		//zero out fields in PML region
		mode_apply_zero_bc_pml(*psi);
		
		//build Arnoldi matrix
		krylov[0] = mode_copy(*psi);
		mode_H_normalize(krylov[0]);
		mode_apply_bc(krylov[0], confinement);
		
		for(i = 1; i <= size; i++){

			krylov[i] = mode_apply_radial_beta_squared_operator(krylov[i - 1], omega, radius_left);
			mode_apply_bc(krylov[i], confinement);

			for(j = 0; j < i; j++){

				H_index = (i - 1) + j*size;

				inner_prod = mode_H_inner_product(krylov[j], krylov[i]);
				H.m[H_index] = inner_prod;
				TEMP = matrix_complex_smul(inner_prod, krylov[j].fields);
				matrix_sub_in_place(krylov[i].fields, TEMP);
				matrix_free(TEMP);

			}

			H_index = (i - 1) + i*size;

			norm = mode_H_norm(krylov[i]);
			H.m[H_index] = complex_from_doubles(norm, 0);
			mode_H_normalize(krylov[i]);
		}

		/*printf("\nH = \n");
		matrix_print(H);
		printf("\n\n");*/

		//D = matrix_diagonalize_non_hermitian(H, omega*omega/1000, 1000);
		D = matrix_diagonalize(H, 1e-3, 1000);
		/*printf("\nD = \n");
		matrix_print(D);
		printf("\n\n");*/
		
		//count positive eigenvalues and obtain ordered list
		mode_count = 0;
		max = 1000;
		max_index = 0;
		for(i = 0; i < size; i++){
  
			eigenvalues[i] = D.m[i + i*size];
			
			if(fabs(eigenvalues[i].im) < max){
						
				max = fabs(eigenvalues[i].im);
				max_index = i;
			}
		}

		/*max = 0;
		max_index = 0;
		for(i = 0; i < size; i++){

			eigenvalues[i] = D.m[i + i*size];

			if(eigenvalues[i].re > max){

				max = eigenvalues[i].re;
				max_index = i;
			}
		}*/

		reduced_eigenvector = matrix_eigenvector(H, eigenvalues[max_index]);


		beta = complex_sqrt(eigenvalues[max_index]);
		n_eff = complex_smul(1.0/(omega*(radius_left + psi->i*psi->disc/2)), beta);


		//zero out fields for psi
		for(j = 0; j < 6*psi->i*psi->j; j++){

			psi->fields.m[j] = complex_zero();

		}

		//set fields equal to aproximated eigenvector
		for(j = 0; j < size; j++){

			TEMP = matrix_complex_smul(reduced_eigenvector.m[j], krylov[j].fields);			
			matrix_add_in_place(psi->fields,TEMP);
			matrix_free(TEMP);
		}
		matrix_free(reduced_eigenvector);

		printf("\n#%i: beta = %f + i%e\nn_eff = %f + i%e\n\n", m, beta.re, beta.im, n_eff.re, n_eff.im);


		for(i = 0; i <= size; i++){
			
			mode_free(krylov[i]);
		}
		matrix_free(D);

		n_old = n_mag;
		n_mag = complex_mag(n_eff);
		n_diff = n_mag - n_old;
		m++;	
	}
	
	//mode_apply_zero_bc(*psi);
	mode_apply_radial_H_div(*psi, beta, radius_left);
	mode_apply_bc(*psi, confinement);
	mode_apply_radial_H_curl(*psi, beta, omega, radius_left);
	mode_apply_bc(*psi, confinement);
	mode_H_normalize(*psi);
	
	matrix_free(H);
	free(eigenvalues);
	free(krylov);
	return(beta);		
}

complex mode_true_beta(complex beta){
	
	complex true_beta;
	
	true_beta = complex_smul(1e6, beta);
	
	return(true_beta);
}

mode mode_true_fields(mode psi){
	
	int i, index;
	double ep0, mu0;
	mode true_fields;
	double norm;
	
	ep0 = 8.85418782e-12;
	mu0 = 1.25663706e-6;
	
	true_fields = mode_copy(psi);
	
	for(i = 0; i < psi.i*psi.j; i++){
		
		index = 6*i;
		
		true_fields.fields.m[index] = complex_smul(sqrt(mu0/ep0), psi.fields.m[index]);
		true_fields.fields.m[index + 1] = complex_smul(sqrt(mu0/ep0), psi.fields.m[index + 1]);
		true_fields.fields.m[index + 2] = complex_smul(sqrt(mu0/ep0), psi.fields.m[index + 2]);
		true_fields.fields.m[index + 3] = psi.fields.m[index + 3];
		true_fields.fields.m[index + 4] = psi.fields.m[index + 4];
		true_fields.fields.m[index + 5] = psi.fields.m[index + 5];
	}

	mode_power_normalize(true_fields);
	
	return(true_fields);
}

double mode_n_susceptibility_rectangle(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if((i >= left)&&(i < right)&&(j >= bottom)&&(j < top)){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_exclude_rectangle(mode psi, int left, int right, int bottom, int top, int left_exclude, int right_exclude, int bottom_exclude, int top_exclude){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((i >= left)&&(i < right)&&(j >= bottom)&&(j < top))&&((i < left_exclude)||(i >= right_exclude)||(j < bottom_exclude)||(j >= top_exclude))){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_top(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((j == top)||(j == top - 1))&&(i >= left)&&(i < right)){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_bottom(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((j == bottom)||(j == bottom - 1))&&(i >= left)&&(i < right)){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_sides(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((i == left)||(i == left - 1)||(i == right)||(i == right - 1))&&(j >= bottom)&&(j < top)){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_sides_outside(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((i == left - 1)||(i == right))&&(j >= bottom)&&(j < top)){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_sides_normal(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((i == left)||(i == left - 1)||(i == right)||(i == right - 1))&&(j >= bottom)&&(j < top)){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_sides_tangential(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((i == left)||(i == left - 1)||(i == right)||(i == right - 1))&&(j >= bottom)&&(j < top)){
				
				E_squared += complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_sides_longitudinal(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((i == left)||(i == left - 1)||(i == right)||(i == right - 1))&&(j >= bottom)&&(j < top)){
				
				E_squared += complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_rectangle_edge(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((i == left)||(i == left - 1)||(i == right)||(i == right - 1))&&(j >= bottom - 1)&&(j <= top)){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
			if(((j == bottom)||(j == bottom - 1)||(j == top)||(j == top - 1))&&(i >= left - 1)&&(i <= right)){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_n_susceptibility_bottom_half_plane(mode psi, int top){
	
	int i, j, index;
	double norm, divisor;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(j < top){
				
				E_squared += complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez);
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_loss_susceptibility_rectangle(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index, ep_index;
	double norm, divisor, real_n;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			ep_index = 3*(i + j*psi.i);
			
			real_n = complex_real(complex_sqrt(psi.ep.m[ep_index]));
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if((i >= left)&&(i < right)&&(j >= bottom)&&(j < top)){
				
				E_squared += 2*real_n*(complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez));
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_loss_susceptibility_rectangle_top(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index, ep_index;
	double norm, divisor, real_n;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			ep_index = 3*(i + j*psi.i);
			
			real_n = complex_real(complex_sqrt(psi.ep.m[ep_index]));
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((j == top)||(j == top - 1))&&(i >= left)&&(i < right)){
				
				E_squared += 2*real_n*(complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez));
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_loss_susceptibility_rectangle_bottom(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index, ep_index;
	double norm, divisor, real_n;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			ep_index = 3*(i + j*psi.i);
			
			real_n = complex_real(complex_sqrt(psi.ep.m[ep_index]));
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((j == bottom)||(j == bottom - 1))&&(i >= left)&&(i < right)){
				
				E_squared += 2*real_n*(complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez));
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_loss_susceptibility_rectangle_sides(mode psi, int left, int right, int bottom, int top){
	
	int i, j, index, ep_index;
	double norm, divisor, real_n;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			ep_index = 3*(i + j*psi.i);
			
			real_n = complex_real(complex_sqrt(psi.ep.m[ep_index]));
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(((i == left)||(i == left - 1)||(i == right)||(i == right - 1))&&(j >= bottom)&&(j < top)){
				
				E_squared += 2*real_n*(complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez));
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

double mode_loss_susceptibility_bottom_half_plane(mode psi, int top){
	
	int i, j, index, ep_index;
	double norm, divisor, real_n;
	complex Ex, Ey, Ez, Hx, Hy;
	complex product1, product2, product3, product4;
	double E_squared, real_cross_prod, index_susceptibility;

	E_squared = 0;
	real_cross_prod = 0;
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){

			index = 6*(i + j*psi.i);
			ep_index = 3*(i + j*psi.i);
			
			real_n = complex_real(complex_sqrt(psi.ep.m[ep_index]));
			
			Ex = psi.fields.m[index];
			Ey = psi.fields.m[index + 1];
			Ez = psi.fields.m[index + 2];
			Hx = psi.fields.m[index + 3];
			Hy = psi.fields.m[index + 4];
			
			product1 = complex_mul(complex_conj(Ex), Hy);
			product2 = complex_sub(complex_zero(), complex_mul(complex_conj(Ey), Hx));
			product3 = complex_sub(complex_zero(), complex_mul(complex_conj(Hx), Ey));
			product4 = complex_mul(complex_conj(Hy), Ex);
			
			real_cross_prod += complex_mag(complex_add(complex_add(product1, product2), complex_add(product3, product4)));
			
			if(j < top){
				
				E_squared += 2*real_n*(complex_mag(Ex)*complex_mag(Ex) + complex_mag(Ey)*complex_mag(Ey) + complex_mag(Ez)*complex_mag(Ez));
			}
		}
	}
	
	index_susceptibility = E_squared/real_cross_prod;
	
	return(index_susceptibility);
}

void mode_apply_lateral_symmetry(mode psi){
	
	int i, j, field_index, symmetry_index;
	mode reflection, symmetric_mode;
	complex inner_product;
	
	//create reflection
	reflection = mode_new(psi.i, psi.j, psi.disc);	
	for(i = 0; i < psi.i; i++){
		for(j = 0; j < psi.j; j++){
			
			field_index = 6*(i + j*psi.i);
			symmetry_index = 6*((psi.i - 1 - i) + j*psi.i);
			
			reflection.fields.m[symmetry_index + 3] = complex_sub(complex_zero(), psi.fields.m[field_index + 3]);
			reflection.fields.m[symmetry_index + 4] = psi.fields.m[field_index + 4];
		}
	}
	
	//perform symmetrization
	inner_product = mode_H_inner_product(psi, reflection);
	if(inner_product.re > 0){
	
		symmetric_mode = mode_linear_combination(psi, reflection, complex_from_doubles(1,0));
	}
	else{
		
		symmetric_mode = mode_linear_combination(psi, reflection, complex_from_doubles(-1,0));	
	}
	mode_free(reflection);
	
	//copy fields and release
	mode_copy_fields(psi, symmetric_mode);
	mode_free(symmetric_mode);	
}