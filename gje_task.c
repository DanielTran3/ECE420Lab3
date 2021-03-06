#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "Lab3IO.h"
#include "timer.h"


int main(int argc, char *argv[]) {
	
	//# pragma omp_set_num_threads(atoi(argv[1]));
	int thread_count = atoi(argv[1]);
	double **initial_mat;
	double* result_mat;
	int* row_index;
	int size_check;
	int i, row, temp, temp_row_index_max, temp_row_value_max, diagonal_index, row_elim_index, column_elim_index, size;
	double start_time, end_time, factorial_elim, largest;
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
	
	// See about making this part parallel (Starting)
	# pragma omp parallel for 
	for (i = 0; i < size; ++i) {
		row_index[i] = i;
	}


	if (size == 1) {
        result_mat[0] = initial_mat[0][1] / initial_mat[0][0];
	}

	else {
		GET_TIME(start_time);
		/*Gaussian elimination*/
		# pragma omp parallel private(temp_row_index_max, temp_row_value_max) num_threads(thread_count)
		{
			# pragma omp single
			{
		    	for (diagonal_index = 0; diagonal_index < size - 1; diagonal_index++){
					// Pivoting
					temp_row_index_max = diagonal_index;
					temp_row_value_max = 0;
					# pragma omp task
					for (row = diagonal_index; row < size; row++) {
						if (temp_row_value_max < fabs(initial_mat[row][diagonal_index])) {
							temp_row_index_max = row;
							temp_row_value_max = fabs(initial_mat[row][temp_row_index_max]);
						}
					}
					if (temp_row_value_max > largest) {
						#pragma critical
						{
							if (temp_row_value_max > largest) {
								largest = temp_row_value_max;					
							}
						}
					}
					// Swap
					if (diagonal_index != temp_row_index_max) {
						temp = row_index[diagonal_index];
						row_index[diagonal_index] = row_index[temp_row_index_max];
						row_index[temp_row_index_max] = temp;
					}

					# pragma omp task private(factorial_elim)
					for (row_elim_index = diagonal_index + 1; row_elim_index < size; row_elim_index++) {				
						factorial_elim = initial_mat[row_index[row_elim_index]][diagonal_index] / 
										 initial_mat[row_index[diagonal_index]][diagonal_index];
						# pragma omp task
						for (column_elim_index = diagonal_index; column_elim_index < size + 1; column_elim_index++) {
							initial_mat[row_index[row_elim_index]][column_elim_index] -= factorial_elim *
							initial_mat[row_index[diagonal_index]][column_elim_index];
						}
					}
				}

				// Jordan Elimination
				for (diagonal_index = size - 1; diagonal_index > 0; diagonal_index--) {
					# pragma omp task private(factorial_elim)
					for (row_elim_index = diagonal_index - 1; row_elim_index >= 0; row_elim_index--) {
						factorial_elim = initial_mat[row_index[row_elim_index]][diagonal_index] / 
										 initial_mat[row_index[diagonal_index]][diagonal_index];
						initial_mat[row_index[row_elim_index]][diagonal_index] = 0;
						initial_mat[row_index[row_elim_index]][size] -= factorial_elim * initial_mat[row_index[diagonal_index]][size];		
					}
				}
		
				# pragma omp task
				for (i = 0; i < size; i++) {
					result_mat[i] = initial_mat[row_index[i]][size] / initial_mat[row_index[i]][i];
				}
			}
		}
		GET_TIME(end_time);
	}
	Lab3SaveOutput(result_mat, size, end_time - start_time);
	return 0;
}
