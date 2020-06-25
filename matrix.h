#ifndef _matrix_h

//set the token
#define _matrix_h

#include "complex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//define the structure and typedef it for ease of use
struct matrix_
{
	int i, j; //number of rows, number of columns
	complex *m; //pointer to matrix values
};
typedef struct matrix_  matrix;

//list the functions written in the source file

//general matrix functions
matrix matrix_new(int i, int j); //returns a matrix with memory allocated and all entries equal to 0
void matrix_free(matrix A); //frees the array associated with A
void matrix_print(matrix A); //prints our each entry
matrix matrix_copy(matrix A); //returns a copy of matrix A
void matrix_copy_in_place(matrix A, matrix B); //copies the contents of B into A in place
matrix matrix_get_column(matrix A, int j); //returns the j'th column of A
matrix matrix_get_row(matrix A, int i); //returns the i'th row of A
void matrix_set_column(matrix A, matrix col, int j); //sets j'th column of A to col; requires that col be a column vector
void matrix_set_row(matrix A, matrix row, int i); //sets the i'th row of A to row; requires that row be a row vector
matrix matrix_get_submatrix(matrix A, int i, int j); //returns the submatrix obtained by removing the i'th row, and j'th column
complex matrix_determinant(matrix A); //returns the complex determinant of A; requires that A be a square matrix
complex matrix_trace(matrix A); //returns the complex trace of A; requires that A be a square matrix
matrix matrix_hermitian(matrix A); //returns the hermitian of A
matrix matrix_smul(double scalar, matrix A); //returns C = scalar*A
matrix matrix_complex_smul(complex scalar, matrix A); //returns C = scalar*A, where scalar is complex
matrix matrix_add(matrix A, matrix B); //returns C = A + B
matrix matrix_sub(matrix A, matrix B); //returns C = A - B
void matrix_smul_in_place(double scalar, matrix A); //multiplies A by scalar in place
void matrix_complex_smul_in_place(complex scalar, matrix A); //multiplies A by a complex scalar in place
void matrix_add_in_place(matrix A, matrix B); //add B to A in place
void matrix_sub_in_place(matrix A, matrix B); //subtracts B from A in place
matrix matrix_mul(matrix A, matrix B); //returns C = AB
matrix matrix_linear_solve(matrix A, matrix B); //returns x that solves Ax = b (uses LU decomposition)
matrix matrix_inverse(matrix A); //returns the inverse of A (uses LU decomposition)
void matrix_gram_schmidt(matrix A); //orthonormalizes columns of A to form an orthogonal basis for col(A), destroys A
matrix matrix_diagonalize(matrix A, double tol, int max_iterations); //returns similar diagonal matrix D =  Q_inv*A*Q using QL method (less round-off errors for larger numbers in bottom right); requires that A be square
matrix matrix_diagonalize_non_hermitian(matrix A, double tol, int max_iterations); //uses QL decomp and 2x2's for non-symmetric matrices
matrix matrix_diagonalize_2_by_2(matrix A); //solves quadratic for eigenvalues of any 2x2 matrix
matrix matrix_eigenvector(matrix A, complex eigenvalue); //returns normalized eigenvector corresponding to eigenvalue

//matrix matrix_choleski_linear_solve(matrix A, matrix b); //returns x that solves Ax = b (uses choleski, only works for symmetric positive-def)
//matrix matrix_choleski_inverse(matrix A, matrix b); //returns the inverse of A (uses choleski, only works for symmetric positive-def)

//vector functions
complex matrix_vector_inner_product(matrix A, matrix B); //returns the inner product of a A_herm*B; requires that A and B are column vectors.
double matrix_vector_norm(matrix A); //returns the L_2 norm of A; requires that A be a column vector
matrix matrix_vector_normalize(matrix A); //returns normalized vector A
void matrix_vector_normalize_in_place(matrix A); //normalizes vector A in place

#endif