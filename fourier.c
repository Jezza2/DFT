/*
DISCRETE FOURIER TRANSFORM

Author:  Jeremy Stanger
Date:		 22/09/16 

Main file; where the business happens.

The programme finds the fourier transform of a function
and plots it.  It also finds the convolution of two
functions.

The function f(x) depends on the mode in which the 
programme executes.  In one mode, the convolution of
a function with itself is found, and then the Fourier
transform of it is found.  By the convolution theorem,
this is compared to the square of the fourier transform
of the original function.
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
// Defines f(x)
double complex *real_space;
// Holds fourier transform, F(u)
double complex *freq_space;
// Holds convulution of two functions
double complex *convolved;

int main(int argc, char *argv[]) {

	set_params(argc, argv);

	/* Set f(x) according to the execution mode */
	switch (mode) {
		case 0:
			construct_slit(real_space, 10, 1.0, 0);	
			break;
		case 1:
			construct_slit(real_space, 10, 1.0, -10);
			break;
		case 2:
			construct_slit(real_space, 20, 1.0, 0);
			break;
		case 3:
			construct_slit(real_space, 20, 0.5, 0);
			break;
		case 4:
			construct_double_slit(real_space, 20, 1.0, 15);
			break;
		case 5:
			construct_double_slit(real_space, 20, 1.0, 25);
			break;
		case 6:
			construct_double_slit(real_space, 40, 1.0, 25);
			break;
		case 7:
			construct_double_slit(real_space, 40, 0.5, 25);
			break;
		case 8:
			/* This mode finds the convolution of a single 
				 slit with itself, FTs and plots it, before 
				 proceeding as normal. */
			construct_slit(real_space, 20, 1.0, 0);
			convolve(convolved, real_space, real_space);
			write_datafile(convolved, "conv");
			plot("conv");
			dft(convolved, freq_space);
			write_datafile(freq_space, "convfreq");
			plot("convfreq");
			break;
	}	

	dft(real_space, freq_space);
	// If in convolution mode, square the FT.
	if (mode == 8) {
		multiply(freq_space, freq_space, freq_space);
	}

	write_datafile(freq_space, "freq");
	plot("freq");

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
	int i, j;

	memset(output, 0, 2*N*sizeof(double complex));

	/* We loop through all points of the FT, and then
		 sum over the values of the function*exp(...) */
	for (j = -N; j < N; j++) {
		// There is no point in doing the summation if f(x)
		// is 0.  Since f(x) is often largely 0, this offers
		// huge performance improvements, allowing for much
		// higher resolution in our plots.
		if (*(input + indexof(j)) != 0) {
			for (i = -N; i < N; i++) {
				// Implement DFT formula
				*(output + indexof(i)) += (*(input + indexof(j))
														* cexp(M_PI * I * ((double)j) * ((double)i) 
														/ (((double)N)))) / (2.0*((double)N));
			}
		}
	}
}

/* This function populates output with a simulated
	 single slit light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 width: The width of the slit
	 height: The intensity of the light from the slit
	 centre: The position of the centre of the slit
	 */
void construct_slit(double complex *output, int width, double height, int centre) {
	// index variable
	int i;
	// flag - 1 for width odd
	int width_odd;
	// Set width_odd.  If odd, take half_width to be floor(width/2)
	int half_width = (width_odd = width%2) ? (width-1)/2 : width/2;

	// Most of the array is 0.  Memset is faster tan looping
	memset(output, 0, 2*N * sizeof(double complex));

	// Set output array to height for half_wudth either side
	// of centre.  If width is odd, then we need one extra
	// on the positive side.
	for (i = -half_width; i < half_width + width_odd; i++) {
		*(output + indexof(centre + i)) = height;
	}
}

/* This function populates output with a simulated
	 double slit light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 width: The width of the slits
	 height: The intensity of the light from the slits
	 centre_distance: The distance of the centres of 
	 									the slits from x=0
	 */
void construct_double_slit(double complex *output, int width, double height, int centre_distance) {
	// index variable
	int i;
	// flag - 1 for width odd.
	int width_odd;
	// Set width_odd.  If odd, take half_width to be floor(width/2)
	int half_width = (width_odd = width%2) ? (width-1)/2 : width/2;

	// Most of the array is 0.  Memset is faster tan looping
	memset(output, 0, 2*N * sizeof(double complex));

	// Set output array to height for half_wudth either side
	// of centres.  If width is odd, then we need one extra
	// on the positive side.
	for (i = -half_width; i < half_width + width_odd; i++) {
		*(output + indexof(-centre_distance+i)) = height;
	}
	for (i = -half_width; i < half_width + width_odd; i++) {
		*(output + indexof(centre_distance+i)) = height;
	}
}

/* Produce the convolution of the functions defined in input1
	 and input2.
	 
	 *output: pointer to start of the array where convolution is
	 					to be stored
	 *input : pointer to start of the arrays to be convolved
	 */
void convolve(double complex *output, double complex *input1, double complex *input2) {
	// index variables
	int X, i;

	// Most of the array is 0.  Memset is faster tan looping
	memset(output, 0, 2*N * sizeof(double complex));

	// Loop over X and perform integral for each value of X
	for (X = -N; X < N; X++) {
		for (i = -N; i < N; i++) {
			// Can't access arrays out of bounds
			if (X-i >= -N && X-i < N) {
				*(output+indexof(X)) += (*(input1+indexof(i))) * (*(input2 + indexof(X-i))); 	
			}
		}
	}
}

/* Multiply, element-by-element, input1 and input2.
	 Store the result in output.

	 *output: pointer to start of array in which result
	 					is to be stored.
	 *input : pointers to starts of arrays to multiplied.
	 */
void multiply(double complex *output, double complex *input1, double complex *input2) {
	// index variable
	int i;

	for (i = 0; i < 2*N; i++) {
		*(output+i) = (*(input1+i)) * (*(input2+i));
	}
}
