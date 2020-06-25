#ifndef _complex_h

//set the token
#define _complex_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//define the structure and typedef it for ease of use
struct complex_
{
  double re, im; //real part and imaginary parts of complex number
};
typedef struct complex_  complex;

//list the functions written in the source file

//self-contained functions
double complex_mag(complex x); //returns the absolute value of complex argument
double complex_angle(complex x); //returns the phase of complex argument (0 <= angle < 2Pi)
double complex_real(complex x); //returns the real part of complex argument
double complex_imag(complex x); //returns the imaginary part of complex argument
complex complex_zero(); //returns a complex zero
complex complex_from_doubles(double real, double imag); //returns c = real + i*imag
complex complex_conj(complex x); //returns the complex conjugate of x
complex complex_add(complex x, complex y); //returns c = x + y
complex complex_sub(complex x, complex y); //returns c = x - y
complex complex_smul(double scalar, complex x); //returns c = scalar*x
complex complex_mul(complex x, complex y); //returns c = x*y

//dependent functions
complex complex_div(complex x, complex y); //returns c = x/y
complex complex_reciprocal(complex x); //returns c = 1/x
complex complex_sqrt(complex x); //returns principal value of c = sqrt(x) (value in right-half of complex plane)
complex complex_exp(complex x); //returns c = exp(x);
complex complex_cos(complex x); //returns c = cos(x);
complex complex_sin(complex x); //returns c = sin(x);
 
#endif