#include "complex.h"

#define PI 3.14159265359
#define SERIES_TOLERANCE 1e-30

double complex_mag(complex x){
	
	double magnitude;
	
	magnitude = sqrt(x.re*x.re + x.im*x.im);
	
	return(magnitude);
}

double complex_angle(complex x){
	
	double ratio, principal_angle, angle;
	
	//calculate principal angle on {-Pi/2, Pi/2}
	ratio = (x.im)/(x.re);
	principal_angle = atan(ratio);
	
	//map to actual angle on {0, 2Pi}
	angle = principal_angle;
	if(x.re < 0){

		angle += PI;		
	}
	else if(x.im < 0){
		
		angle += 2*PI;
	}
	
	return(angle);
}

double complex_real(complex x){
	
	double real;
	
	real = x.re;
	
	return(real);
}

double complex_imag(complex x){
	
	double imag;
	
	imag = x.im;
	
	return(imag);
}

complex complex_zero(){
	
	complex zero;
	
	zero.re = 0;
	zero.im = 0;
	
	return(zero);
}

complex complex_from_doubles(double real, double imag){
	
	complex val;
	
	val.re = real;
	val.im = imag;
	
	return(val);
}

complex complex_conj(complex x){
	
	complex conjugate;
	
	conjugate.re = x.re;
	conjugate.im = -x.im;
	
	return(conjugate);
}

complex complex_add(complex x, complex y){
	
	complex sum;
	
	sum.re = x.re + y.re;
	sum.im = x.im + y.im;
	
	return(sum);
}

complex complex_sub(complex x, complex y){
	
	complex diff;
	
	diff.re = x.re - y.re;
	diff.im = x.im - y.im;
	
	return(diff);
}

complex complex_smul(double scalar, complex x){
	
	complex prod;
	
	prod.re = scalar*x.re;
	prod.im = scalar*x.im;
	
	return(prod);
}

complex complex_mul(complex x, complex y){
	
	complex prod;
	
	prod.re = x.re*y.re - x.im*y.im;
	prod.im = x.re*y.im + x.im*y.re;
	
	return(prod);
}

complex complex_div(complex x, complex y){
	
	double denom, scalar;
	complex num, quotient;
	
	//multiply top and bottom by conjugate of y
	num = complex_mul(x, complex_conj(y));
	denom = complex_mag(y)*complex_mag(y); 
	
	//scalar multiply numerator by 1/denom
	scalar = 1.0/denom;
	quotient = complex_smul(scalar, num);
	
	return(quotient);	
}

complex complex_reciprocal(complex x){
	
	complex reciprocal, one;
	
	one = complex_from_doubles(1, 0);
	reciprocal = complex_div(one, x);
	
	return(reciprocal);
}

complex complex_sqrt(complex x){
	
	double mag, angle;
	complex square_root;
	
	mag = sqrt(complex_mag(x));
	
	angle = complex_angle(x);
	if(angle > PI){
		
		angle -= 2*PI;
	}
	angle /= 2;
	
	square_root = complex_smul(mag, complex_exp(complex_from_doubles(0, angle)));
	return(square_root);
}

complex complex_exp(complex x){
	
	int i;
	double factorial, scalar;
	complex power, exponential, diff;
	
	//compute as power series, starting with the zero term
	diff = complex_from_doubles(1, 0);
	exponential = complex_zero();
	i = 0;
	factorial = 1;
	power = complex_from_doubles(1, 0);
	
	while(complex_mag(diff) >= SERIES_TOLERANCE){
		
		scalar = 1.0/factorial;
		diff = complex_smul(scalar, power);
		exponential = complex_add(exponential, complex_smul(scalar, power));		
		i++;
		factorial *= i;
		power = complex_mul(power, x);
	}
	
	return(exponential);
}

complex complex_cos(complex x){
	
	int i;
	double factorial, scalar;
	complex power, cosine, diff;
	
	//compute cos(x) as power series, starting with the zero term
	diff = complex_from_doubles(1, 0);
	cosine = complex_zero();
	i = 0;
	factorial = 1;
	power = complex_from_doubles(1, 0);
	
	while(complex_mag(diff) >= SERIES_TOLERANCE){
		
		if(i%2 == 0){
			
			scalar = 1.0/factorial;
			diff = complex_smul(scalar, power);
			cosine = complex_add(cosine, complex_smul(scalar, power));
			power = complex_smul(-1, power);		
		}
		
		i++;
		factorial *= i;
		power = complex_mul(power, x);
	}
	
	return(cosine);
}

complex complex_sin(complex x){
	
	int i;
	double factorial, scalar;
	complex power, sine, diff;
	
	//compute sin(x) as power series, starting with the 1st order term
	diff = complex_from_doubles(1, 0);
	sine = complex_zero();
	i = 0;
	factorial = 1;
	power = complex_from_doubles(1, 0);
	
	while(complex_mag(diff) >= SERIES_TOLERANCE){
		
		if(i%2 == 1){
			
			scalar = 1.0/factorial;
			diff = complex_smul(scalar, power);
			sine = complex_add(sine, complex_smul(scalar, power));
			power = complex_smul(-1, power);		
		}
		
		i++;
		factorial *= i;
		power = complex_mul(power, x);
	}
	
	return(sine);
}