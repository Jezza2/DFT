/*
DISCRETE FOURIER TRANSFORM

Author:  Jeremy Stanger
Date:		 22/09/16 

This file is basically the housekeeping; no interesting computing
going on here.
*/

#include <stdio.h>
#include <stdlib.h>
// C s native support for complex numbers is ideal
#include <complex.h>
#include "header.h"

/* Write a gnuplot script to produce the required plot
	 Plot depends on programme execution mode.
	 Gnuplot can then be called externally

	 name[]: string identifier for the plot
	 */
void plot (char name[]) {

	// File pointer
	FILE *fp;
	// Buffer to store filename
	char filename[35] = { };

	// Set filename according to execution mode, string
	// identifier and N.
	sprintf(filename, "plots/plot_%s_m%d_N%d.p", name, mode, N);
	// Open file for writing if possible.  If not, quit
	// with error message
	if ((fp = fopen(filename, "w")) == NULL) {
		printf("Unable to open file for plotting\n");
		_exit(3);
	}

	fprintf(fp,
		"set terminal jpeg size 1000,750\n"
		
		"set autoscale\n"
		"set zeroaxis\n"
		"set xlabel \"u\"\n"
		"set xtic auto\n"
		"set ytic auto\n"
		"set ztic auto\n"
		"set key off\n"
		"set grid\n"
	);

	if (mode == 8) {
		fprintf(fp, "set ylabel \"Re(h(X))\"\n");
		fprintf(fp, "set xlabel \"X\"\n");
		fprintf(fp, "set output \"real_%s_m%d_N%d.jpg\"\n", name, mode, N);
		fprintf(fp, "plot \"../data/data_%s_m%d_N%d.dat\" u 2:3 with lines\n", name, mode, N);
		fprintf(fp, "set ylabel \"|h(X)|\"\n");
		fprintf(fp, "set xlabel \"X\"\n");
		fprintf(fp, "set output \"abs_%s_m%d_N%d.jpg\"\n", name, mode, N);
		fprintf(fp, "plot \"../data/data_%s_m%d_N%d.dat\" u 2:5 with lines\n", name, mode, N);
	} else {
		fprintf(fp, "set ylabel \"Re(F(u))\"\n");
		fprintf(fp, "set output \"real_%s_m%d_N%d.jpg\"\n", name, mode, N);
		fprintf(fp, "plot \"../data/data_%s_m%d_N%d.dat\" u 2:3 with lines\n", name, mode, N);
		fprintf(fp, "set ylabel \"|F(u)|\"\n");
		fprintf(fp, "set output \"abs_%s_m%d_N%d.jpg\"\n", name, mode, N);
		fprintf(fp, "plot \"../data/data_%s_m%d_N%d.dat\" u 2:5 with lines\n", name, mode, N);
	}

	// Close file
	fclose(fp);
}

/* Writes data in array to a file for plotting purposes
	 Space delimited file:
	 0 index real_part imag_part magnitude

	 The 0 column is to facilitate a vector plot, though this
	 is not done from this programme

	 *array: data to be written to the file
	 name[]: string identifier for data file.
	 */
void write_datafile(double complex *array, char name[]) {
	// File pointer
	FILE *fp;
	// index variable
	int i;
	// buffer for the filename
	char filename[50] = { };

	// Set filename according to identifier, N and execution mode
	sprintf(filename, "data/data_%s_m%d_N%d.dat", name, mode, N);
	// Open file if possible, otherwise quit with error message
	if ((fp = fopen(filename, "w")) == NULL) {
		printf("Unable to open file to write data\n");
		_exit(4);
	}

	// Write data to file
	for (i = 0; i < 2*N; i++) {
		fprintf(fp, "0 %d %.9g %.9g %.9g\n", i-N,
			creal(*(array+i)), cimag(*(array+i)), cabs(*(array+i)));
	}
	// Close file
	fclose(fp);
}

/* Get user parameters and set global variables accordingly.
	 Also allocate memory required for the data storage.

	 count:  Number of parameters
	 *argvec[]: pointer to string array of the parameters
	 */
void set_params(int count, char *argvec[]) {
	// If _exit() is called before these are assigned then a free() will be
	// attempted with garbage pointers.  If they are NULL, free() will
	// do nothing.
	real_space = NULL;
	freq_space = NULL;
	convolved = NULL;

	// Check the number of inputs is correct.
	// If not, display help and exit cleanly.
	if (count != 2) {
		help();
		_exit(2);
	}

	// Assign input parameters.  Display help and exit cleanly if error.
	// Mode is the only parameter that can sensibly be zero, otherwise
	// return value 0 indicates error.
	if ((mode = atoi(*(++argvec))) == 0 && (*argvec)[0] != '0') {
		help();
		_exit(2);
	}

	// Allocate memory for data arrays. Exit cleanly if there is an error.
	if ( ( real_space = malloc(2*N * sizeof(double complex)) ) == NULL
		|| ( freq_space = malloc(2*N * sizeof(double complex)) ) == NULL ) {
		printf("Unable to allocate memory for data array(s)");
		_exit(1);
	}

	// When in convolution mode we need a third array for data storage
	if ( mode == 8 && (convolved = malloc(2*N * sizeof(double complex))) == NULL ) {
		printf("Unable to allocate memory for data array(s)");
		_exit(1);
	}
}

/* Exit cleanly, freeing any dynamically allocated memory first.
	 Good practice.  Display exit status.

	 code: Exit status
	 */
void _exit(int code) {
	free(real_space);
	free(freq_space);
	free(convolved);
	printf("Exit status: %d\n", code);
	exit(code);
}

/* Display parameter list and modes
	 */
void help(void) {
	printf("3 input parameters are required as follows in the following order:\n\n"
				 "int mode            The mode number\n\n\n");

	printf("EXIT STATUSES:\n\n"
				 "0 -                 Successfully executed\n"
				 "1 -                 Unable to allocate memory for data arrays\n"
				 "2 -                 Input error\n"
				 "3 -                 Unable to open file for plotting\n"
				 "4 -                 Unable to open file to write data\n\n\n");

	printf("MODES:\n\n"
				 "0 -                 Single slit width 10, height 1.0, centre 0\n"
				 "1 -                 Single slit width 10, height 1.0, centre -10\n"
				 "2 -                 Single slit width 20, height 1.0, centre 0\n"
				 "3 -                 Single slit width 20, height 0.5, centre 0\n"
				 "4 -                 Double slit width 20, height 1.0, centres +/- 15\n"
				 "5 -                 Double slit width 20, height 1.0, centres +/- 25\n"
				 "6 -                 Double slit width 40, height 1.0, centres +/- 25\n"
				 "7 -                 Double slit width 40, height 0.5, centres +/- 25\n"
				 "8 -                 Single slit width 20, height 1.0, centre 0\n"
				 "                    Convolve with itself, square FT\n");
}
