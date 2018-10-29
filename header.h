/*
DISCRETE FOURIER TRANSFORM

Author:  Jeremy Stanger
Date:    22/09/2016

Header file.
Function prototypes, macros and external variable declarations

*/


// Some definitions and macros
// x ranges between -N to N-1 whereas the array is zero-indexed.
// This macro maps x to array index
#define indexof(x) ((x)+N)
#define N 1000

// Function prototypes
// Functions in io.c
/* Display parameter list and modes
	 */
void help(void);   

/* Exit cleanly, freeing any dynamically allocated memory first.
	 Good practice.  Display exit status.

	 code: Exit status
	 */
void _exit(int code);

/* Get user parameters and set global variables accordingly.
	 Also allocate memory required for the data storage.

	 count:  Number of parameters
	 *argvec[]: pointer to string array of the parameters
	 */
void set_params(int count, char *argvec[]);

/* Writes data in array to a file for plotting purposes
	 Space delimited file:
	 0 index real_part imag_part magnitude

	 The 0 column is to facilitate a vector plot, though this
	 is not done from this programme

	 *array: data to be written to the file
	 name[]: string identifier for data file.
	 */
void write_datafile(double complex *array, char name[]);

/* Write a gnuplot script to produce the required plot
	 Plot depends on programme execution mode.
	 Gnuplot can then be called externally

	 name[]: string identifier for the plot
	 */
void plot(char name[]);


// Functions in schrodinger.c
/* This function finds the FT of the discrete function
	 defined in the input array, and stores it in the output
	 array. 
	 
	 *input: pointer to the start of the array containing the 
	 				 function to be transformed.
	 *output: pointer to the start of the array where the FT
	 					will be stored.
	 */
void dft(double complex *input, double complex *output);

/* This function populates output with a simulated
	 single slit light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 width: The width of the slit
	 height: The intensity of the light from the slit
	 centre: The position of the centre of the slit
	 */
void construct_slit(double complex *output, int width, double height, int centre);

/* This function populates output with a simulated
	 double slit light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 width: The width of the slits
	 height: The intensity of the light from the slits
	 centre_distance: The distance of the centres of 
	 									the slits from x=0
	 */
void construct_double_slit(double complex *output, int width, double height, int centre_distance);

/* Produce the convolution of the functions defined in input1
	 and input2.
	 
	 *output: pointer to start of the array where convolution is
	 					to be stored
	 *input : pointer to start of the arrays to be convolved
	 */
void convolve(double complex *output, double complex *input1, double complex *input2);

/* Multiply, element-by-element, input1 and input2.
	 Store the result in output.

	 *output: pointer to start of array in which result
	 					is to be stored.
	 *input : pointers to starts of arrays to multiplied.
	 */
void multiply(double complex *output, double complex *input1, double complex *input2);


/* Declare global variables */
// Execution mode number
extern int mode;

// Data storage arrays
// The native complex data type is ideal for this application
// Defines f(x)
extern double complex *real_space;
// Holds fourier transform, F(u)
extern double complex *freq_space;
// Holds convulution of two functions
extern double complex *convolved;
