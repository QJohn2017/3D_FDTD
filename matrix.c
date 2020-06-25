#include "matrix.h"

//general matrix functions
matrix matrix_new(int i, int j){
	
	int counter;
	matrix A;
	complex *m;
	
	m = malloc(i*j*sizeof(complex));
	for(counter = 0; counter < i*j; counter++){
		
		m[counter] = complex_zero();
	}
	
	A.i = i;
	A.j = j;
	A.m = m;
	
	return(A);
}

void matrix_free(matrix A){
	
	free(A.m);
}

void matrix_print(matrix A){
	
	int i, j, index;
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			//printf("A(%i, %i) = %f + i%f\n", (i + 1), (j + 1), A.m[index].re, A.m[index].im);
			printf("%f + i%f\t\t", A.m[index].re, A.m[index].im);
		}
		
		printf("\n");
	}
}

matrix matrix_copy(matrix A){
	
	int i, j, index;
	matrix copy;
	
	copy = matrix_new(A.i, A.j);
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			copy.m[index] = A.m[index];
		}
	}
	
	return(copy);
}

void matrix_copy_in_place(matrix A, matrix B){

	int i, j, index;
	
	//ensure that A and B have the same dimensions
	if((A.i != B.i)||(A.j != B.j)){
		
		printf("Matrix dimensions do not match for addition!\n");
		abort();
	}
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			A.m[index] = B.m[index];
		}
	}	
}

matrix matrix_get_column(matrix A, int j){
	
	int i, index;
	matrix col;
	
	col = matrix_new(A.i, 1);
	for(i = 0; i < A.i; i++){
		
		index = j + i*A.j;
		col.m[i] = A.m[index];
	}
	
	return(col);
}

matrix matrix_get_row(matrix A, int i){
	
	int j, index;
	matrix row;
	
	row = matrix_new(1, A.j);
	for(j = 0; j < A.j; j++){
		
		index = j + i*A.j;
		row.m[j] = A.m[index];
	}
	
	return(row);
}

void matrix_set_column(matrix A, matrix col, int j){
	
	int i, n, index;
	
	//ensure that col is a column vector
	if(col.j != 1){
		
		printf("Column is not a column vector for set_column!\n");
		abort();
	}
	//ensure that col is the same length as A
	if(col.i != A.i){
		
		printf("Column is the wrong length for set_column!");
		abort();
	}
	
	n = A.i;
	for(i = 0; i < n; i++){
		
		index = j + i*A.j;
		A.m[index] = col.m[i];
	}
}

void matrix_set_row(matrix A, matrix row, int i){
	
	int j, n, index;
	
	//ensure that row is a row vector
	if(row.i != 1){
		
		printf("Row is not a row vector for set_row!\n");
		abort();
	}
	//ensure that row is the same length as A
	if(row.j != A.j){
		
		printf("Row is the wrong length for set_row!");
		abort();
	}
	
	n = A.j;
	for(j = 0; j < n; j++){
		
		index = j + i*A.j;
		A.m[index] = row.m[j];
	}
}

matrix matrix_get_submatrix(matrix A, int i, int j){
	
	int k, l, src_index, dst_index;
	matrix M;
	
	//ensure that i and j are in-bounds
	if((i >= A.i)||(i >= A.j)){
		
		printf("Row or column is out of bounds for getting minor matrix!\n");
		abort();
	}
	
	M = matrix_new((A.i - 1), (A.j - 1));
	for(k = 0; k < A.i; k++){
		for(l = 0; l < A.j; l++){
			
			src_index = l + k*A.j;
			dst_index = l + k*M.j;
			
			if((k != i)&&(l != j)){
				
				if(k > i){
				
					dst_index -= M.j;
				}
				if(l > j){

					dst_index--;
				}

				M.m[dst_index] = A.m[src_index];
			}
		}
	}
	
	return(M);
}

complex matrix_determinant(matrix A){
	
	int j, n, index;
	double sign;
	complex det, addition_term, sub_det;
	matrix M;
	
	//make sure A is square
	if(A.i != A.j){
		
		printf("A is not square for finding determinant!\n");
		abort();
	}
	
	n = A.i;
	//going along first row, the determinant is A(1,j)*det(M(1,j)), using recursion
	if(n == 1){
		
		det = A.m[0];
	}
	else{
		
		det = complex_zero();
		sign = 1;
		for(j = 0; j < n; j++){
			
			index = j;
			
			M = matrix_get_submatrix(A, 0, j);
			sub_det = matrix_determinant(M);
			sub_det = complex_mul(sub_det, A.m[index]);
			sub_det = complex_smul(sign, sub_det);
			det = complex_add(det, sub_det);
			
			sign *= -1;
		}
	}
	
	return(det);
}

complex matrix_trace(matrix A){
	
	int i, n, index;
	complex trace;
	
	//ensure that A is square
	//make sure A is square
	if(A.i != A.j){
		
		printf("A is not square for finding trace!\n");
		abort();
	}
	
	n = A.i;
	trace = complex_zero();
	for(i = 0; i < n; i++){
		
		index = i + i*n;
		trace = complex_add(trace, A.m[index]);
	}
	
	return(trace);
}

matrix matrix_hermitian(matrix A){
	
	int i, j, A_index, B_index;
	matrix B;
	
	B = matrix_new(A.j, A.i);
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			A_index = j + i*A.j;
			B_index = i + j*B.j;
			
			B.m[B_index] = complex_conj(A.m[A_index]);
		}
	}
	
	return(B);
}

matrix matrix_smul(double scalar, matrix A){
	
	int i, j, index;
	matrix B;
	
	B = matrix_new(A.i, A.j);
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			B.m[index] = complex_smul(scalar, A.m[index]);
		}
	}
	
	return(B);	
}

matrix matrix_complex_smul(complex scalar, matrix A){
	
	int i, j, index;
	matrix B;
	
	B = matrix_new(A.i, A.j);
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			B.m[index] = complex_mul(scalar, A.m[index]);
		}
	}	
	
	return(B);
}

matrix matrix_add(matrix A, matrix B){
	
	int i, j, index;
	matrix C;
	
	//ensure that A and B have the same dimensions
	if((A.i != B.i)||(A.j != B.j)){
		
		printf("Matrix dimensions do not match for addition!\n");
		abort();
	}
	
	C = matrix_new(A.i, A.j);
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			C.m[index] = complex_add(A.m[index], B.m[index]);
		}
	}
	
	return(C);
}

matrix matrix_sub(matrix A, matrix B){
	
	int i, j, index;
	matrix C;
	
	//ensure that A and B have the same dimensions
	if((A.i != B.i)||(A.j != B.j)){
		
		printf("Matrix dimensions do not match for subtraction!\n");
		abort();
	}
	
	C = matrix_new(A.i, A.j);
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			C.m[index] = complex_sub(A.m[index], B.m[index]);
		}
	}
	
	return(C);
}

void matrix_smul_in_place(double scalar, matrix A){
	
	int i, j, index;
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			A.m[index] = complex_smul(scalar, A.m[index]);
		}
	}	
}

void matrix_complex_smul_in_place(complex scalar, matrix A){
	
	int i, j, index;

	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			A.m[index] = complex_mul(scalar, A.m[index]);
		}
	}	
}

void matrix_add_in_place(matrix A, matrix B){
	
	int i, j, index;
	
	//ensure that A and B have the same dimensions
	if((A.i != B.i)||(A.j != B.j)){
		
		printf("Matrix dimensions do not match for addition!\n");
		abort();
	}
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			A.m[index] = complex_add(A.m[index], B.m[index]);
		}
	}
}

void matrix_sub_in_place(matrix A, matrix B){
	
	int i, j, index;
	
	//ensure that A and B have the same dimensions
	if((A.i != B.i)||(A.j != B.j)){
		
		printf("Matrix dimensions do not match for subtraction!\n");
		abort();
	}
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < A.j; j++){
			
			index = j + i*A.j;
			A.m[index] = complex_sub(A.m[index], B.m[index]);
		}
	}
}

matrix matrix_mul(matrix A, matrix B){
	
	matrix C;
	int i, j, k, A_index, B_index, C_index;
	
	//ensure that A and B have the right dimensions
	if(A.j != B.i){
		
		printf("matrix inner dimensions do not match for multiplication!\n");
		abort();
	}
	
	//do multiplication
	C = matrix_new(A.i, B.j);
	
	for(i = 0; i < A.i; i++){
		for(j = 0; j < B.j; j++){

			C_index = j + i*B.j;
			C.m[C_index] = complex_zero();
			for(k = 0; k < A.j; k++){
			
				A_index = k + i*A.j;
				B_index = j + k*B.j;
				C.m[C_index] = complex_add(C.m[C_index], complex_mul(A.m[A_index], B.m[B_index]));
			}
		}
	}
	
	return(C);
}

matrix matrix_linear_solve(matrix A, matrix B){
	
	int i, j, k;
	int index, L_index, U_index, B_index, y_index, x_index;
	int pivot_index, max_index, temp;
	int *row_map;
	double element_mag, max;
	matrix LU, x, y_hat;

	//ensure that A is square
	if(A.i != A.j){
	
		printf("A matrix is not square for LU linear solve!\n");
		abort();
	}
	//ensure B is the right dimension
	if(A.i != B.i){
		
		printf("matrix dimensions do not match LU linear solve!\n");
		abort();
	}
	
	//allocate matrices and row_map array
	LU = matrix_copy(A);	
	x = matrix_new(LU.i, B.j);
	y_hat = matrix_new(LU.i, B.j);
	
	row_map = malloc(LU.i*sizeof(int));
	for(i = 0; i < LU.i; i++){
		
		row_map[i] = i;
	}

	//perform LU decomp column-by-column
	for(j = 0; j < LU.j; j++){
		
		//start above the diagonal
		for(i = 0; i < j; i++){
			
			//subtract L(i,k)*U(k,j) for k < i
			index = j + row_map[i]*LU.j;
			for(k = 0; k < i; k++){
				
				L_index = k + row_map[i]*LU.j;
				U_index = j + row_map[k]*LU.j;
				
				LU.m[index] = complex_sub(LU.m[index], complex_mul(LU.m[L_index], LU.m[U_index]));	
			}
		}
		
		//diagonal and below, make the subtractions and find pivot
		max = 0;
		for(i = j; i < LU.i; i++){
			
			//subtract L(i,k)*U(k,j) for k < j 
			index = j + row_map[i]*LU.j;
			for(k = 0; k < j; k++){
				
				L_index = k + row_map[i]*LU.j;
				U_index = j + row_map[k]*LU.j;

				LU.m[index] = complex_sub(LU.m[index], complex_mul(LU.m[L_index], LU.m[U_index]));
			}
			//pivot is chosen as largest-magnitude element in column and that row becomes the j'th row in row_map
			element_mag = complex_mag(LU.m[index]);
			if(element_mag > max){
					
				max = element_mag;
				max_index = i;
			}
		}
		//change rows in row map to match pivot choice
		temp = row_map[j];
		row_map[j] = row_map[max_index];
		row_map[max_index] = temp;		
		
		//below diagonal, divide by pivot
		for(i = j + 1; i < LU.i; i++){
	
			index = j + row_map[i]*LU.j;
			pivot_index = j + row_map[j]*LU.j;
			LU.m[index] = complex_div(LU.m[index], LU.m[pivot_index]);
		}
	}
	
	//solve for y_hat (L*y_hat = B_hat); fill in columns from top to bottom	
	for(i = 0; i < LU.i; i++){
		for(j = 0; j < B.j; j++){
		
			index = j + i*B.j;
			B_index = j + row_map[i]*B.j;
			
			//subtract L(i,k)*y_hat(k,j) for k < i 
			y_hat.m[index] = B.m[B_index];
			for(k = 0; k < i; k++){
			
				L_index = k + row_map[i]*LU.i;
				y_index = j + k*y_hat.j;
				y_hat.m[index] = complex_sub(y_hat.m[index], complex_mul(LU.m[L_index], y_hat.m[y_index]));
			}
		}
	}
	
	//solve for x (U*x = y_hat); fill in columns from bottom to top
	for(i = LU.i - 1; i >= 0; i--){
		for(j = 0; j < B.j; j++){
		
			index = j + i*B.j;
			
			//subtract U(i,k)*x(k,j) for k > i
			x.m[index] = y_hat.m[index];
			for(k = LU.i - 1; k > i; k--){
			
				U_index = k + row_map[i]*LU.i;
				x_index = j + k*x.j;
				x.m[index] = complex_sub(x.m[index], complex_mul(LU.m[U_index], x.m[x_index]));
			}
			//divide by U(i,i)
			U_index = i + row_map[i]*LU.i;
			x.m[index] = complex_div(x.m[index], LU.m[U_index]);

		}
	}

	free(row_map);
	matrix_free(LU);
	matrix_free(y_hat);
	return(x);
}

matrix matrix_inverse(matrix A){
	
	int i, n, index;
	matrix I, A_inv;
	
	//ensure that A is square
	if(A.i != A.j){
	
		printf("A matrix is not square for LU linear solve!\n");
		abort();
	}
	
	//build identity (same size as A)
	n = A.i;
	I = matrix_new(n, n);
	for(i = 0; i < n; i++){
		
		index = i + i*n;
		I.m[index] = complex_from_doubles(1, 0);
	}
	
	//solve for right-inverse of A --> A*A_inv = I
	A_inv = matrix_linear_solve(A, I);
	
	matrix_free(I);
	return(A_inv);
}

void matrix_gram_schmidt(matrix A){
	
	int j, k;
	double norm;
	complex inner_product;
	matrix current, prev, proj;
	
	//ensure that there are fewer columns than or equal number of columns to number of rows --> dim(row(A)) >= dim(col(A))
	if(A.j > A.i){
		
		printf("Number of columns exceeds dimension of rowspace for Gram-Schmidt!\n");
		abort();
	}
	
	//go through column by column
	for(j = 0; j < A.j; j++){
		
		//normalize, (subtract projection onto a previous column, normalize) --> until all previous columns are done
		current = matrix_get_column(A, j);
		norm = matrix_vector_norm(current);
		matrix_smul_in_place(1.0/norm, current);
		
		for(k = 0; k < j; k++){
			
			//projection
			prev = matrix_get_column(A, k);
			inner_product = matrix_vector_inner_product(prev, current);
			proj = matrix_complex_smul(inner_product, prev);
			matrix_free(prev);
			
			//subtraction
			matrix_sub_in_place(current, proj);
			matrix_free(proj);
			
			//normalization
			norm = matrix_vector_norm(current);
			matrix_smul_in_place(1.0/norm, current);
		}
				
		//replace column with orthonormalized column
		matrix_set_column(A, current, j);		
		matrix_free(current);
	}		
}

matrix matrix_diagonalize(matrix A, double tol, int max_iterations){
	
	int i, j, k, n, index, major_index, minor_index, counter;
	double max_upper, element_mag, diag_mag, max_diag, ratio;
	complex inner_product, unit_phasor;
	matrix current_column, minor_column, median, proj, v, v_herm;
	matrix D, I, V_minor, V, H, H_herm, Q, L, TEMP;
	
	//ensure that A is square
	if(A.i != A.j){
		
		printf("Matrix is not square for QL diagonalization!\n");
		abort();
	}
	
	n = A.i;
	D = matrix_copy(A);
	
	//going to use QL method (smaller roundoff error that QR for larger values in lower right-hand corner) while large off-diagonals exist
	max_upper = 0;
	max_diag = 0;
	for(i = 0; i < n; i++){
		
		index = i + i*n;
		diag_mag = complex_mag(D.m[index]);
		if(diag_mag > max_diag){
		
			max_diag = diag_mag;
		}
		
		for(j = 0; j < n; j++){
			
			index = j + i*n;
			element_mag = complex_mag(D.m[index]);
			if((i < j)&&(element_mag > max_upper)){
				
				max_upper = element_mag;
			}
		}
	}
	index = n*n - 1;
	ratio = max_upper/max_diag;
		
	counter = 0;
	while(ratio > tol){

		//set L to D, and Q to I initially
		L = matrix_copy(D);
		Q = matrix_new(n, n);
		for(i = 0; i < n; i++){
		
			index = i + i*n;
			Q.m[index] = complex_from_doubles(1, 0);
		}
	
		//going from right-most column to left-most column
		for(j = (n - 1); j > 0; j--){
			
			//get current column, then only the part above and including diagonal element
			current_column = matrix_get_column(L, j);
			minor_column = matrix_new((j + 1), 1);
			for(i = 0; i <= j; i++){
			
				minor_column.m[i] = current_column.m[i];
			}
			matrix_free(current_column);
			
			//from the minor column, obtain the householder vector to map to transpose(..., 0, 0, 0, norm1)
			//obtain normalized median vector (unit vector pointing halfway between minor column and desired (..., 0, 0, 0, norm1))
			median = matrix_vector_normalize(minor_column);
			element_mag = complex_mag(L.m[j + j*L.j]);
			
			if(element_mag > 0){
				
				unit_phasor = complex_smul(1.0/element_mag, L.m[j + j*L.j]);
			}
			else{
				
				unit_phasor = complex_from_doubles(1, 0);
			}
			//printf("unit_phasor = %f + i%f\n", unit_phasor.re, unit_phasor.im);
			
			median.m[j] = complex_add(median.m[j], unit_phasor);
			matrix_vector_normalize_in_place(median);

			//obtain householder vector v as the component of minor column vector orthogonal to the median, then normalized
			inner_product = matrix_vector_inner_product(median, minor_column);
			proj = matrix_complex_smul(inner_product, median);
			matrix_free(median);
		
			v = matrix_sub(minor_column, proj);
			matrix_free(minor_column);
			matrix_free(proj);

			if(matrix_vector_norm(v) > 0){
				
				matrix_vector_normalize_in_place(v);	
			}

			v_herm = matrix_hermitian(v);
			//build appropriate householder matrix H = I - 2*v*v_herm
			I = matrix_new(n, n);
			V = matrix_new(n, n);
			for(i = 0; i < n; i++){

				index = i + i*n;
				I.m[index] = complex_from_doubles(1, 0);
			}

			V_minor = matrix_mul(v, v_herm);
			matrix_free(v);
			matrix_free(v_herm);
			
			matrix_smul_in_place(2, V_minor);
			for(i = 0; i < (j + 1); i++){
				for(k = 0; k < (j + 1); k++){

					minor_index = k + i*(j+1);
					major_index = k + i*n;
					V.m[major_index] = V_minor.m[minor_index];
				}
			}

			H = matrix_sub(I, V);

			matrix_free(I);
			matrix_free(V_minor);
			matrix_free(V);

			//update L by L = H*L
			TEMP = matrix_mul(H, L);
			matrix_copy_in_place(L, TEMP);
			matrix_free(TEMP);

			//update Q by Q = Q*H_herm
			H_herm = matrix_hermitian(H);
			TEMP = matrix_mul(Q, H_herm);
			matrix_copy_in_place(Q, TEMP);
			matrix_free(TEMP);

			matrix_free(H);
			matrix_free(H_herm);
		}
	
		TEMP = matrix_mul(L, Q);
		matrix_copy_in_place(D, TEMP);
		matrix_free(TEMP);
		
		matrix_free(L);
		matrix_free(Q);
		
		max_upper = 0;
		max_diag = 0;
		for(i = 0; i < n; i++){
			
			index = i + i*n;
			diag_mag = complex_mag(D.m[index]);
			if(diag_mag > max_diag){

				max_diag = diag_mag;
			}
			
			for(j = 0; j < n; j++){

				index = j + i*n;
				element_mag = complex_mag(D.m[index]);
				
				if((i < j)&&(element_mag > max_upper)){

					max_upper = element_mag;
				}
			}
		}
		index = n*n - 1;
		ratio = max_upper/max_diag;
		
		counter++;
		
		//printf("QL method off_diagonals = %f\n", max_upper);
		if(counter > max_iterations){
			
			printf("QL algorithm reached maximum iterations!\n");
			//abort();
			break;
		}
		
	}
	
	//zero all but the diagonals
	TEMP = matrix_new(n, n);
	for(i = 0; i < n; i++){
		
		index = i + i*n;
		TEMP.m[index] = D.m[index];
	}
	matrix_copy_in_place(D, TEMP);
	matrix_free(TEMP);
	
	return(D);
}

matrix matrix_diagonalize_non_hermitian(matrix A, double tol, int max_iterations){
	
	int i, j, k, n, index, major_index, minor_index, counter, first_upper_switch;
	double max_upper, element_mag;
	complex inner_product, unit_phasor, a, b, c, d, lamda1, lamda2;
	matrix current_column, minor_column, median, proj, v, v_herm;
	matrix D, I, V_minor, V, H, H_herm, Q, Q_herm, L, TEMP;
	
	//ensure that A is square
	if(A.i != A.j){
		
		printf("Matrix is not square for QL diagonalization!\n");
		abort();
	}
	
	n = A.i;
	D = matrix_copy(A);
	
	//going to use QL method (smaller roundoff error that QR for larger values in lower right-hand corner) while large off-diagonals exist
	max_upper = 0;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			
			index = j + i*n;
			element_mag = complex_mag(D.m[index]);
			if((i < j)&&(element_mag > max_upper)){
				
				max_upper = element_mag;
			}
		}
	}
	
	counter = 0;
	while(max_upper > tol){

		//set L to D, and Q to I initially
		L = matrix_copy(D);
		Q = matrix_new(n, n);
		for(i = 0; i < n; i++){
		
			index = i + i*n;
			Q.m[index] = complex_from_doubles(1, 0);
		}
	
		//going from right-most column to left-most column
		for(j = (n - 1); j > 0; j--){
		
			//get current column, then only the part above and including diagonal element
			current_column = matrix_get_column(L, j);
			minor_column = matrix_new((j + 1), 1);
			for(i = 0; i <= j; i++){
			
				minor_column.m[i] = current_column.m[i];
			}
			matrix_free(current_column);
		
			//from the minor column, obtain the householder vector to map to transpose(..., 0, 0, 0, norm1)
			//obtain normalized median vector (unit vector pointing halfway between minor column and desired (..., 0, 0, 0, norm1))
			median = matrix_vector_normalize(minor_column);
			element_mag = complex_mag(L.m[j + j*L.j]);
			
			if(element_mag > 0){
				
				unit_phasor = complex_smul(1.0/element_mag, L.m[j + j*L.j]);
			}
			else{
				
				unit_phasor = complex_from_doubles(1, 0);
			}
			//printf("unit_phasor = %f + i%f\n", unit_phasor.re, unit_phasor.im);
			
			median.m[j] = complex_add(median.m[j], unit_phasor);
			matrix_vector_normalize_in_place(median);

			//obtain householder vector v as the component of minor column vector orthogonal to the median, then normalized
			inner_product = matrix_vector_inner_product(median, minor_column);
			proj = matrix_complex_smul(inner_product, median);
			matrix_free(median);
		
			v = matrix_sub(minor_column, proj);
			matrix_free(minor_column);
			matrix_free(proj);

			if(matrix_vector_norm(v) > 0){
				
				matrix_vector_normalize_in_place(v);	
			}

			v_herm = matrix_hermitian(v);
			//build appropriate householder matrix H = I - 2*v*v_herm
			I = matrix_new(n, n);
			V = matrix_new(n, n);
			for(i = 0; i < n; i++){

				index = i + i*n;
				I.m[index] = complex_from_doubles(1, 0);
			}

			V_minor = matrix_mul(v, v_herm);
			matrix_smul_in_place(2, V_minor);
			for(i = 0; i < (j + 1); i++){
				for(k = 0; k < (j + 1); k++){

					minor_index = k + i*(j+1);
					major_index = k + i*n;
					V.m[major_index] = V_minor.m[minor_index];
				}
			}

			H = matrix_sub(I, V);

			matrix_free(I);
			matrix_free(V_minor);
			matrix_free(V);

			//update L by L = H*L
			TEMP = matrix_mul(H, L);
			matrix_copy_in_place(L, TEMP);
			matrix_free(TEMP);

			//update Q by Q = Q*H_herm
			H_herm = matrix_hermitian(H);
			TEMP = matrix_mul(Q, H_herm);
			matrix_copy_in_place(Q, TEMP);
			matrix_free(TEMP);

			matrix_free(H);
			matrix_free(H_herm);
		}
	
		//Q_herm = matrix_hermitian(Q);
		TEMP = matrix_mul(L, Q);
		//matrix_free(Q_herm);
		matrix_copy_in_place(D, TEMP);
		matrix_free(TEMP);
		matrix_free(L);
		matrix_free(Q);
		
		//take note of largest upper. don't count isolated non-zeroes on first superdiagonal
		max_upper = 0;
		first_upper_switch = 0;
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){

				index = j + i*n;
				element_mag = complex_mag(D.m[index]);
				if((i < j)&&(element_mag > max_upper)){

					if(j == i + 1){
						
						if((first_upper_switch == 1)&&(element_mag > tol)){
							
							max_upper = element_mag;
							first_upper_switch = 0;
						}
						else{
							
							first_upper_switch = 1;
						}
					}
					else{

						max_upper = element_mag;
					}
				}
			}
		}
		
		counter++;
		if(counter > max_iterations){
			
			printf("QL algorithm reached maximum iterations!\n");
			//abort();
			break;
		}
	}
	
	//take note of 2-by-2s
	for(i = 0; i < (n - 1); i++){
		
		index = (i + 1) + i*n;
		/*if(complex_mag(D.m[index]) > 1.01*tol){*/
		if(i == 0){
			
			a = D.m[index - 1];
			b = D.m[index];
			c = D.m[index + n - 1];
			d = D.m[index + n];

			lamda1 = complex_smul(0.5, complex_add(complex_add(a,d), complex_sqrt(complex_add(complex_mul(complex_sub(a, d), complex_sub(a, d)), complex_smul(4, complex_mul(b, c))))));
			lamda2 = complex_smul(0.5, complex_sub(complex_add(a,d), complex_sqrt(complex_add(complex_mul(complex_sub(a, d), complex_sub(a, d)), complex_smul(4, complex_mul(b, c))))));

			D.m[index - 1] = lamda1;
			D.m[index] = complex_zero();
			D.m[index + n - 1] = complex_zero();
			D.m[index + n] = lamda2;
		}
	}
	
	//zero all but the diagonals
	/*TEMP = matrix_new(n, n);
	for(i = 0; i < n; i++){
		
		index = i + i*n;
		TEMP.m[index] = D.m[index];
	}
	matrix_copy_in_place(D, TEMP);*/
	
	return(D);
}

matrix matrix_diagonalize_2_by_2(matrix A){
	
	matrix D;
	complex a, b, c, d;
	complex lamda1, lamda2;
	
	if((A.i != 2)||(A.j != 2)){
		
		printf("Matrix is not 2x2 for 2x2 diagonalization!\n");
		abort();
	}
	
	D = matrix_new(2, 2);

	a = A.m[0];
	b = A.m[1];
	c = A.m[2];
	d = A.m[3];

	lamda1 = complex_smul(0.5, complex_add(complex_add(a,d), complex_sqrt(complex_add(complex_mul(complex_sub(a, d), complex_sub(a, d)), complex_smul(4, complex_mul(b, c))))));
	lamda2 = complex_smul(0.5, complex_sub(complex_add(a,d), complex_sqrt(complex_add(complex_mul(complex_sub(a, d), complex_sub(a, d)), complex_smul(4, complex_mul(b, c))))));

	D.m[0] = lamda1;
	D.m[1] = complex_zero();
	D.m[2] = complex_zero();
	D.m[3] = lamda2;

	return(D);
}

matrix matrix_eigenvector(matrix A, complex eigenvalue){
	
	int i, j, n, index;
	matrix V, B, bottom_row;
	matrix b, eigenvector;
	
	//ensure that A is square
	if(A.i != A.j){
	
		printf("A matrix is not square for eigenvector solve!\n");
		abort();
	}
	
	n = A.i;
	V = matrix_new(n, n);
	bottom_row = matrix_new(1, n);
	b = matrix_new(n, 1);
	eigenvector = matrix_new(n, 1);
	
	//V is I*lambda
	for(i = 0; i < n; i++){
		
		index = i + i*n;
		V.m[index] = eigenvalue;
		bottom_row.m[i] = complex_from_doubles(1, 0);
	}
	B = matrix_sub(A, V);
	
	//replace bottom row of B with [1, 1, 1, ...]
	matrix_set_row(B, bottom_row, (n - 1));
	
	//set bottom value of b to n
	b.m[n - 1] = complex_from_doubles((double)n, 0);
	
	eigenvector = matrix_linear_solve(B, b);
	matrix_vector_normalize_in_place(eigenvector);
	
	matrix_free(V);
	matrix_free(bottom_row);
	matrix_free(b);
	matrix_free(B);
	
	return(eigenvector);
}

//vector functions
complex matrix_vector_inner_product(matrix A, matrix B){
	
	int i;
	complex inner_product;
	matrix A_herm, matrix_product;	
	
	//ensure that A and B are both column vectors
	if((A.j != 1)||(B.j != 1)){
	
		printf("Matrices are not column vectors for inner product!\n");
		abort();
	}
	//ensure that A and B are the same lenght
	if(A.i != B.i){
		
		printf("Vectors are not the same length for inner product!\n");
		abort();
	}
	
	//calculate inner product as hermitian(A)*B
	A_herm = matrix_hermitian(A);
	matrix_product = matrix_mul(A_herm, B);
	inner_product = matrix_product.m[0];
	
	matrix_free(A_herm);
	matrix_free(matrix_product);
	return(inner_product);
}

double matrix_vector_norm(matrix A){
	
	int i;
	double norm;
	
	//ensure that A is a column vector
	if(A.j != 1){
	
		printf("Matrix is not a vector for norm!\n");
		abort();
	}
	
	norm = 0;
	for(i = 0; i < A.i; i++){
		
		norm += A.m[i].re*A.m[i].re + A.m[i].im*A.m[i].im;
	}
	
	norm = sqrt(norm);
	
	return(norm);
}

matrix matrix_vector_normalize(matrix A){
	
	matrix normalized_A;
	double norm;
	
	norm = matrix_vector_norm(A);
	normalized_A = matrix_smul(1.0/norm, A);
	
	return(normalized_A);
}

void matrix_vector_normalize_in_place(matrix A){
	
	double norm;
	
	norm = matrix_vector_norm(A);
	matrix_smul_in_place(1.0/norm, A);
}