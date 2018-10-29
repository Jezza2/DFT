/*
2D DISCRETE FOURIER TRANSFORM

Author:  Jeremy Stanger
Date:		 22/09/16 

Main file; where the business happens.

The programme finds the fourier transform of a function
and plots it.

The function f(x) depends on the mode in which the 
programme executes.
*/


#include <stdio.h>
#include <stdlib.h>
// C s native support for complex numbers is ideal
#include <complex.h>
#include <string.h>
#include <math.h>
#include "header.h"


/* Global variable definitions */
// Execution mode number
int mode;
// Defines f(x,y)
double complex *real_space;
// Holds fourier transform, F(u,v)
double complex *freq_space;


int main(int argc, char *argv[]) {

	set_params(argc, argv);

	/* Set f(x,y) according to the execution mode */
	switch (mode) {
		case 0:
			construct_single(real_space, 0, 0, 1.0);	
			break;
		case 1:
			construct_squares(real_space, 0, 0, 3, 3, 1.0);
			break;
		case 2:
			construct_slit(real_space, 0, 0, 3, 20, 1.0);
			break;
		case 3:
			construct_doubleslit(real_space, 0, 0, 3, 30, 10, 1.0);
			break;
	}	

	dft(real_space, freq_space);

	write_datafile(real_space, "real");
	write_datafile(freq_space, "freq");
	plot("freq");
	plot("real");

	printf("Successfully executed!\n");
	_exit(0);
}

/* This function finds the FT of the discrete function
	 defined in the input array, and stores it in the output
	 array. 
	 
	 *input: pointer to the start of the array containing the 
	 				 function to be transformed.
	 *output: pointer to the start of the array where the FT
	 					will be stored.
	 */
void dft(double complex *input, double complex *output) {
	// index variables
	int u, v, i, j;

	// Initialize to 0.  Memset is fast
	memset(output, 0, 4*N*M*sizeof(double complex));

	// Loop through all values of input array
	for (i = -N; i < N; i++) {
		for (j = -M; j < M; j++) {
	
			// If the function value is 0 then  there's no
			// point in calculating sum.  This optimisation
			// improves execution speed by orders of magnitude
			if (*(input+indexof(i,j)) != 0.0) {
				for (u = -N; u < N; u++) {
					for (v = -N; v < M; v++) {
						
						// Implement sum
						*(output + indexof(u,v)) += ((*(input + indexof(i,j))) 
							* cexp(-M_PI * I * (((double)(i*u))/((double)N) 
							+ ((double)(j*v))/((double)M))))
								/ ((double)(4.0 * M * N));
					}
				}
			}
		}
	}

}

/* This function populates output with a simulated
	 single cross light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 c_x,c_y: The position of the centre of the slit
	 strength: intensity of light from slit
	 */
void construct_single(double complex *output, int c_x, int c_y, double strength) {
	memset(output, 0, 4*M*N*sizeof(double complex));
	*(output+indexof(c_x  ,c_y  ))   = strength;
	*(output+indexof(c_x+1,c_y  ))   = strength;
	*(output+indexof(c_x  ,c_y+1))   = strength;
	*(output+indexof(c_x-1,c_y  ))   = strength;
	*(output+indexof(c_x  ,c_y-1))   = strength;
}

/* This function populates output with a simulated
	 single square light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 c _x,c _y: The positions of the centres of each of the squares
	 strength: intensity of light from slit
	 */
void construct_squares(double complex *output, int c1_x, int c1_y, 
												int c2_x, int c2_y, double strength) {

	memset(output, 0, 4*M*N*sizeof(double complex));
	*(output+indexof(c1_x  , c1_y  )) = strength;
	*(output+indexof(c1_x+1, c1_y  )) = strength;
	*(output+indexof(c1_x  , c1_y+1)) = strength;
	*(output+indexof(c1_x+1, c1_y+1)) = strength;

	*(output+indexof(c2_x  , c2_y  )) = strength;
	*(output+indexof(c2_x+1, c2_y  )) = strength;
	*(output+indexof(c2_x  , c2_y+1)) = strength;
	*(output+indexof(c2_x+1, c2_y+1)) = strength;
}

/* This function populates output with a simulated
	 single slit light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 c_x,c_y: coords of centre of the slit
	 width: The width of the slit
	 length: The length of the slit
	 strength: The intensity of the light from the slit
	 */
void construct_slit(double complex *output, int c_x, int c_y, 
										int width, int length, double strength) {

	// Most of the function is 0.  Memset is fast.
	memset(output, 0, 4*M*N*sizeof(double complex));

	// index variables
	int i, j;
	// flag - 1 if width odd
	int width_odd;
	// flag - 1 if length odd
	int length_odd;
	// set half_width and half_length.  If odd, then make the respective
	// quantity the floor of half
	int half_width = (width_odd = width%2) ? (width-1)/2 : width/2;
	int half_length = (length_odd = length%2) ? (length-1)/2 : length/2;

	// Set half_width and half_length elements either side of the
	// respective centres to strength.  If either is odd then we want
	// 1 element extra on the positive side of centre.
	for (i = -half_width; i < half_width + width_odd; i++) {
		for (j = -half_length; j < half_length + length_odd; j++) {
			*(output+indexof(c_x+i, c_y+j)) = strength;	
		}
	}
}

/* This function populates output with a simulated
	 double slit light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 c_x,c_y: coords of centre of the slits
	 width: The width of the slit
	 length: The length of the slit
	 centres: Distance of each slit from centre in x-direction
	 strength: The intensity of the light from the slit
	 */
void construct_doubleslit(double complex *output, int c_x, int c_y,
												int width, int length, int centres, double strength) {

	// Most of the function is 0.  Memset is fast
	memset(output, 0, 4*M*N*sizeof(double complex));

	// index variables
	int i, j;
	// flag - 1 if width odd
	int width_odd;
	// flag - 1 if length odd
	int length_odd;
	// set half_width and half_length.  If odd, then make the respective
	// quantity the floor of half
	int half_width = (width_odd = width%2) ? (width-1)/2 : width/2;
	int half_length = (length_odd = length%2) ? (length-1)/2 : length/2;

	// Set half_width and half_length elements either side of the
	// respective centres to strength.  If either is odd then we want
	// 1 element extra on the positive side of centre.
	for (i = -half_width; i < half_width + width_odd; i++) {
		for (j = -half_length; j < half_length + length_odd; j++) {
			*(output+indexof(c_x+i+centres, c_y+j)) = strength;	
		}
	}
	for (i = -half_width; i < half_width + width_odd; i++) {
		for (j = -half_length; j < half_length + length_odd; j++) {
			*(output+indexof(c_x+i-centres, c_y+j)) = strength;	
		}
	}
}
