#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab3IO.h"
#include "timer.h"

// Function definitions
void fill_in_matrix(int* row_index, int size);
void gaus_elim(double** initial_mat, int* row_index, int size);
void jord_elim(double** initial_mat, int* row_index, int size);
void store_result(double **initial_mat, double* result_mat, int* row_index, int size);

int main(int argc, char *argv[]) {
	
	// Matrix variables
	double **initial_mat;
	double* result_mat;

	// Row pointers variables
	int* row_index;

	// Size obtained from file
	int size_check;

	// Size obtained from array
	int size;

	// Start times, end times
	double start_time, end_time;

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
	fill_in_matrix(row_index, size);

	// Check for case if size is equal to 1
	if (size == 1) {
        result_mat[0] = initial_mat[0][1] / initial_mat[0][0];
	}

	else {

		// Get start time
		GET_TIME(start_time);

		/* Gaussian Elimination */
		gaus_elim(initial_mat, row_index, size);
		
		/* Jordan Elimination */
		jord_elim(initial_mat, row_index, size);
		
		// Place the solutions into the result matrix
		store_result(initial_mat, result_mat, row_index, size);
	
		// Get endtime
		GET_TIME(end_time);
	}

	// Save result to file with the correct result and time elapsed
	Lab3SaveOutput(result_mat, size, end_time - start_time);

	return 0;
}

/* Gaussian Elimination */
void gaus_elim(double** initial_mat, int* row_index, int size) {

	// Arbitary loop variables
	int diagonal_index, row, temp, row_elim_index, column_elim_index, temp_row_index_max;

	// Factorial multipler for row subtraction
	double factorial_elim;

	// Cycle through diagonal indexes
	for (diagonal_index = 0; diagonal_index < size - 1; diagonal_index++){
		
		// Keep track of the index with the highest value
		temp_row_index_max = diagonal_index;
		
		// Cycle through loop to find row index	with max value
		for (row = diagonal_index; row < size; row++) {
			if (fabs(initial_mat[row][temp_row_index_max]) < fabs(initial_mat[row][diagonal_index])) {
				temp_row_index_max = row;
			}
		}
	
		// Swap the row with the highest value to correct spot
		if (diagonal_index != temp_row_index_max) {
			temp = row_index[diagonal_index];
			row_index[diagonal_index] = row_index[temp_row_index_max];
			row_index[temp_row_index_max] = temp;
		}
		
		// Create an upper triangle by manipulating rows to add/subtract
		for (row_elim_index = diagonal_index + 1; row_elim_index < size; row_elim_index++) {
			factorial_elim = initial_mat[row_index[row_elim_index]][diagonal_index] / 
							 initial_mat[row_index[diagonal_index]][diagonal_index];
			for (column_elim_index = diagonal_index; column_elim_index < size + 1; column_elim_index++) {
				initial_mat[row_index[row_elim_index]][column_elim_index] -= factorial_elim * initial_mat[row_index[diagonal_index]][column_elim_index];
			}
		}
	}
}

/* Jordan Elimination */
void jord_elim(double** initial_mat, int* row_index, int size) {

	// Arbitary loop variables
	int diagonal_index, row_elim_index;

	// Factorial multipler for row subtraction
	double factorial_elim;

	// Start from bottom-up and remove portions in upper triangle so that only one variable is left. (Diagonal of values)
	for (diagonal_index = size - 1; diagonal_index > 0; diagonal_index--) {

		// Cycle in for loop to remove corresponding variables in the rows above the current one
		for (row_elim_index = diagonal_index - 1; row_elim_index >= 0; row_elim_index--) {
			factorial_elim = initial_mat[row_index[row_elim_index]][diagonal_index] / 
							 initial_mat[row_index[diagonal_index]][diagonal_index];
			initial_mat[row_index[row_elim_index]][diagonal_index] = 0;
			initial_mat[row_index[row_elim_index]][size] -= factorial_elim * initial_mat[row_index[diagonal_index]][size];		
		}
	}
}

// Place the solutions into the result matrix
void store_result(double **initial_mat, double* result_mat, int* row_index, int size) {

	// Arbitary loop variable
	int i;

	// Store solution in result_mat
	for (i = 0; i < size; i++) {
		result_mat[i] = initial_mat[row_index[i]][size] / initial_mat[row_index[i]][i];
	}
}

// Set inital row indexes
void fill_in_matrix(int* row_index, int size) {

	// Arbitary loop variable
	int i;

	// Set inital row indexes
	for (i = 0; i < size; ++i) {
        	row_index[i] = i;
	}
}
