#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "Lab3IO.h"
#include "timer.h"

int main(int argc, char *argv[]) {
	
	// Get thread number from user
	int thread_count = atoi(argv[1]);

	// Matrix variables
	double **initial_mat;
	double* result_mat;

	// Row pointers variables
	int* row_index;

	// Size obtained from file
	int size_check;

	// Arbitary loop variables and size
	int i, row, temp, temp_row_index_max, temp_row_value_max, diagonal_index, row_elim_index, column_elim_index, size;

	// Start times, end times, factorial multipler for row subtraction, and temp variable to find largest value
	double start_time, end_time, factorial_elim, largest;

	// File pointer
	FILE* fp;

	// Load the matrix from file and check if loaded correctly
	Lab3LoadInput(&initial_mat, &size);
	if ((fp = fopen("data_input","r")) == NULL){
		printf("Fail to open the result data file!\n");
		return 2;
	}
	
	// Check if the size is correct
	fscanf(fp, "%d\n\n", &size_check);
	if (size_check != size){
		printf("The problem size of the input file and result file does not match!\n");
		return -1;
	}

	// Create result matrix
	result_mat = CreateVec(size);

	// Create row pointer indexes
	row_index = malloc(size * sizeof(int));
	
	// Set inital row indexes in parallel
	# pragma omp parallel for num_threads(thread_count)
	for (i = 0; i < size; ++i) {
		row_index[i] = i;
	}

	// Check for case if size is equal to 1
	if (size == 1) {
        result_mat[0] = initial_mat[0][1] / initial_mat[0][0];
	}

	else {
		// Get start time
		GET_TIME(start_time);

		/* Gaussian Elimination */
    	for (diagonal_index = 0; diagonal_index < size - 1; diagonal_index++){

			// Keep track of value and index of highest value
			temp_row_index_max = diagonal_index;
			temp_row_value_max = 0;

			// Start team for parallel execution
			# pragma omp parallel num_threads(thread_count) reduction(max:largest)
			{	
	
				// Cycle through loop in parallel to find max value and row index		
				# pragma omp for
				for (row = diagonal_index; row < size; row++) {
					if (temp_row_value_max < fabs(initial_mat[row][diagonal_index])) {
						temp_row_index_max = row;
						temp_row_value_max = fabs(initial_mat[row][temp_row_index_max]);
					}
				}

				// Finds the largest value from team
				if (temp_row_value_max > largest) {
					if (temp_row_value_max > largest) {
						largest = temp_row_value_max;					
					}
				}
			}

			// Swap the row with the highest value to correct spot
			if (diagonal_index != temp_row_index_max) {
				temp = row_index[diagonal_index];
				row_index[diagonal_index] = row_index[temp_row_index_max];
				row_index[temp_row_index_max] = temp;
			}

			// Run team in parallel to create an upper triangle by manipulating rows to add/subtract
			# pragma omp parallel for private(factorial_elim, column_elim_index) num_threads(thread_count) shared(initial_mat, row_index, diagonal_index)
			for (row_elim_index = diagonal_index + 1; row_elim_index < size; row_elim_index++) {	
				factorial_elim = initial_mat[row_index[row_elim_index]][diagonal_index] / 
								 initial_mat[row_index[diagonal_index]][diagonal_index];			
				for (column_elim_index = diagonal_index; column_elim_index < size + 1; column_elim_index++) {
					initial_mat[row_index[row_elim_index]][column_elim_index] -= factorial_elim *
					initial_mat[row_index[diagonal_index]][column_elim_index];
				}
			}
		}

		/* Jordan Elimination */
		// Start from bottom-up and remove portions in upper triangle so that only one variable is left. (Diagonal of values)
		for (diagonal_index = size - 1; diagonal_index > 0; diagonal_index--) {
			
			// Run team in parallel to remove corresponding variables in the rows above the current one
			# pragma omp parallel for private(factorial_elim) num_threads(thread_count)			
			for (row_elim_index = diagonal_index - 1; row_elim_index >= 0; row_elim_index--) {
				factorial_elim = initial_mat[row_index[row_elim_index]][diagonal_index] / 
								 initial_mat[row_index[diagonal_index]][diagonal_index];
				initial_mat[row_index[row_elim_index]][diagonal_index] = 0;
				initial_mat[row_index[row_elim_index]][size] -= factorial_elim * initial_mat[row_index[diagonal_index]][size];		
			}
		}
	
		// Run team in parallel to place the solutions into the result matrix
		# pragma omp parallel for num_threads(thread_count)
		for (i = 0; i < size; i++) {
			result_mat[i] = initial_mat[row_index[i]][size] / initial_mat[row_index[i]][i];
		}

		// Get endtime
		GET_TIME(end_time);
	}
	// Save result to file with the correct result and time elapsed
	Lab3SaveOutput(result_mat, size, end_time - start_time);
	return 0;
}
