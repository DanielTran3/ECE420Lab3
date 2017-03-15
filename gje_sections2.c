#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab3IO.h"
#include "timer.h"

void fill_in_matrix(int* row_index, int size);
void gaus_elim(double** initial_mat, int* row_index, int size);
void find_max_column(double** initial_mat, int* row_index, int size);
void elim_under_diag(double** initial_mat, int* row_index, int size);
void jord_elim(double** initial_mat, int* row_index, int size);
void store_result(double **initial_mat, double* result_mat, int* row_index, int size);

int main(int argc, char *argv[]) {
	
	int thread_count = atoi(argv[1]);
	double **initial_mat;
	double* result_mat;
	int* row_index;
	int size_check;
	int size;
	double start_time, end_time;
	FILE* fp;
	
	Lab3LoadInput(&initial_mat, &size);
	if ((fp = fopen("data_input","r")) == NULL){
		printf("Fail to open the result data file!\n");
		return 2;
	}
	fscanf(fp, "%d\n\n", &size_check);
	if (size_check != size){
		printf("The problem size of the input file and result file does not match!\n");
		return -1;
	}

	result_mat = CreateVec(size);

	row_index = malloc(size * sizeof(int));
	fill_in_matrix(row_index, size);

	if (size == 1) {
        result_mat[0] = initial_mat[0][1] / initial_mat[0][0];
	}

	else {
		GET_TIME(start_time);
		// Gaussian elimination
		# pragma omp parallel num_threads(thread_count)
		{
			//# pragma omp section
			gaus_elim(initial_mat, row_index, size);
			
			//# pragma omp section
			// Jordan Elimination
			jord_elim(initial_mat, row_index, size);		

			//# pragma omp section
			// Store Result
			store_result(initial_mat, result_mat, row_index, size);			
		}
		GET_TIME(end_time);
	}

	Lab3SaveOutput(result_mat, size, end_time - start_time);
	return 0;
}

void gaus_elim(double** initial_mat, int* row_index, int size) {
	find_max_column(initial_mat, row_index, size);
	elim_under_diag(initial_mat, row_index, size);
}

void jord_elim(double** initial_mat, int* row_index, int size) {
	int diagonal_index, row_elim_index;
	double factorial_elim;
	for (diagonal_index = size - 1; diagonal_index > 0; diagonal_index--) {
		# pragma omp sections private(row_elim_index, factorial_elim)
		{		
			# pragma omp section
			{		
				for (row_elim_index = diagonal_index - 1; row_elim_index >= 0; row_elim_index--) {
					factorial_elim = initial_mat[row_index[row_elim_index]][diagonal_index] / 
									 initial_mat[row_index[diagonal_index]][diagonal_index];
					initial_mat[row_index[row_elim_index]][diagonal_index] = 0;
					initial_mat[row_index[row_elim_index]][size] -= factorial_elim * initial_mat[row_index[diagonal_index]][size];		
				}
			}
		}
	}
}

void store_result(double **initial_mat, double* result_mat, int* row_index, int size) {
	int i;
	# pragma omp sections private(i)
	{
		# pragma omp section
		{	
			for (i = 0; i < size; i++) {
				result_mat[i] = initial_mat[row_index[i]][size] / initial_mat[row_index[i]][i];
			}
		}
	}
}

void fill_in_matrix(int* row_index, int size) {
	int i;
	for (i = 0; i < size; ++i) {
        	row_index[i] = i;
	}
}

void find_max_column(double** initial_mat, int* row_index, int size) {
	int diagonal_index, row, temp, temp_row_index_max;
	# pragma omp sections private(diagonal_index, temp_row_index_max, temp)
	{
		# pragma omp section 
		{
			for (diagonal_index = 0; diagonal_index < size - 1; diagonal_index++){
				// Pivoting
				temp_row_index_max = diagonal_index;
				for (row = diagonal_index; row < size; row++) {
					if (fabs(initial_mat[row][temp_row_index_max]) < fabs(initial_mat[row][diagonal_index])) {
						temp_row_index_max = row;
					}
				}	
				// Swap
				if (diagonal_index != temp_row_index_max) {
					temp = row_index[diagonal_index];
					row_index[diagonal_index] = row_index[temp_row_index_max];
					row_index[temp_row_index_max] = temp;
				}
			}
		}
	}
}

void elim_under_diag(double** initial_mat, int* row_index, int size) {
	int diagonal_index, row_elim_index, column_elim_index;
	double factorial_elim;

	for (diagonal_index = 0; diagonal_index < size - 1; diagonal_index++){
		for (row_elim_index = diagonal_index + 1; row_elim_index < size; row_elim_index++) {
			factorial_elim = initial_mat[row_index[row_elim_index]][diagonal_index] / 
							 initial_mat[row_index[diagonal_index]][diagonal_index];
				
					for (column_elim_index = diagonal_index; column_elim_index < size + 1; column_elim_index++) {
						initial_mat[row_index[row_elim_index]][column_elim_index] -= factorial_elim * 
						initial_mat[row_index[diagonal_index]][column_elim_index];
					}
				}
			}
		}
	}

}
