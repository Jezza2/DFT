/*
DISCRETE FOURIER TRANSFORM

Author:  Jeremy Stanger
Date:  

This file is basically the housekeeping; no interesting computing
going on here.
*/

#include <stdio.h>
#include <stdlib.h>
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
	// Buffer for file name
	char filename[35] = { };

	// Set filename according to identifier, mode number and N
	sprintf(filename, "plots/plot_%s_m%d_N%d.p", name, mode, N);
	// Open file for writing if possible, otherwise exit with
	// error message
	if ((fp = fopen(filename, "w")) == NULL) {
		printf("Unable to open file for plotting\n");
		_exit(3);
	}

	fprintf(fp,
		"set terminal jpeg size 1000,750\n"
		
		"set autoscale\n"
		"set zeroaxis\n"
		"set xlabel \"u\"\n"
		"set ylabel \"v\"\n"
		"set xtic auto\n"
		"set ytic auto\n"
		"set ztic auto\n"
		"set key off\n"
		"set contour\n"
		"set grid\n"
		"set dgrid3d %d,%d\n"
		"set pm3d corners2color mean\n", 2*N, 2*M
	);

	if (name == "freq") {
		fprintf(fp, "set zlabel \"Re(F(u, v))\"\n");
		fprintf(fp, "set output \"real_%s_m%d_N%d_M%d.jpg\"\n", name, mode, N, M);
		fprintf(fp, "splot \"../data/data_%s_m%d_N%d_M%d.dat\" u 1:2:3 with pm3d\n", name, mode, N, M);
		fprintf(fp, "set zlabel \"|F(u, v)|\"\n");
		fprintf(fp, "set output \"abs_%s_m%d_N%d_M%d.jpg\"\n", name, mode, N, M);
		fprintf(fp, "splot \"../data/data_%s_m%d_N%d_M%d.dat\" u 1:2:5 with pm3d\n", name, mode, N, M);
	} else {
		fprintf(fp, "set zlabel \"Re(f(x,y))\"\n");
		fprintf(fp, "set xlabel \"x\"\n");
		fprintf(fp, "set ylabel \"y\"\n");
		fprintf(fp, "set output \"real_%s_m%d_N%d_M%d.jpg\"\n", name, mode, N, M);
		fprintf(fp, "splot \"../data/data_%s_m%d_N%d_M%d.dat\" u 1:2:3 with pm3d\n", name, mode, N, M);
		fprintf(fp, "set zlabel \"|f(x, y)|\"\n");
		fprintf(fp, "set output \"abs_%s_m%d_N%d_M%d.jpg\"\n", name, mode, N, M);
		fprintf(fp, "splot \"../data/data_%s_m%d_N%d_M%d.dat\" u 1:2:5 with pm3d\n", name, mode, N, M);
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
	// index variables
	int i, j;
	// Buffer for file name
	char filename[50] = { };

	// Set filename according to identifier, execution mode, N and M
	sprintf(filename, "data/data_%s_m%d_N%d_M%d.dat", name, mode, N, M);
	// Open file for writing if possible, otherwise quit with
	// error message.
	if ((fp = fopen(filename, "w")) == NULL) {
		printf("Unable to open file to write data\n");
		_exit(4);
	}

	// Loop through arrays and write data to file	
	for (i = 0; i < 2*N; i++) {
		for (j = 0; j < 2*M; j++) {
			fprintf(fp, "%d %d %.9g %.9g %.9g\n", i-N, j-M,
				creal(*(array+i*2*N+j)), cimag(*(array+i*2*N+j)), cabs(*(array+i*2*N+j)));
		}
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
	// If this are not NULL before an _exit() is called, the programme
	// will attempt to free memory at a garbage pointer.
	real_space = NULL;
	freq_space = NULL;

	// Check the number of inputs is correct.
	// If not, display help and exit cleanly.
	if (count != 2) {
		help();
		_exit(2);
	}

	// Allocate memory for data arrays if possible, otherwise quit
	// error message
	if ( (real_space = malloc( 2*N * 2*M * sizeof(double complex) )) == NULL
		|| (freq_space = malloc( 2*N * 2*M * sizeof(double complex) )) == NULL ) {
			printf("Unable to allocate memory for data storage");
			_exit(1);
		}

	// Assign input parameters.  Display help and exit cleanly if error.
	// Mode is the only parameter that can sensibly be zero, otherwise
	// return value 0 indicates error.
	if ((mode = atoi(*(++argvec))) == 0 && (*argvec)[0] != '0') {
		help();
		_exit(2);
	}
}

/* Exit cleanly, freeing any dynamically allocated memory first.
	 Good practice.  Display exit status.

	 code: Exit status
	 */
void _exit(int code) {
	free(real_space);
	free(freq_space);
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
				 "0 -                 Single cross source at 0, strength 1.0\n"
				 "1 -                 Two square sources, arbitrary position, strength 1.0\n"
				 "2 -                 Single slit, centre 0, width 3, length 20, strength 1.0\n"
				 "3 -                 Double slit, centre 0, width 3, length 20,\n"
				 "                    centre distance 10, centre distance 10\n\n");
}
