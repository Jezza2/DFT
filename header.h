/*
2D DISCRETE FOURIER TRANSFORM

Author:  Jeremy Stanger
Date:    22/09/2016

Header file.
Function prototypes, macros and external variable declarations

*/

// Some definitions and macros
// This macro maps cordinates to array index
#define indexof(i, j) ( ((i)+N)*2*N + ((j)+M) )
#define N 100 // x dimension
#define M 100 // y dimension

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
	 single cross light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 c_x,c_y: The position of the centre of the slit
	 strength: intensity of light from slit
	 */
void construct_single(double complex *output, int c_x, int c_y, double strength);

/* This function populates output with a simulated
	 single square light source.
	 
	 *output: pointer to start of the array where the
	 					function is to be stored
	 c _x,c _y: The positions of the centres of each of the squares
	 strength: intensity of light from slit
	 */
void construct_squares(double complex *output, int c1_x, int c1_y, 
											int c2_x, int c2_y, double strength);

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
										int width, int height, double strength);

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
										int width, int height, int centres, double strength);


/* Declare global variables */
// Execution mode number
extern int mode;

// Data storage arrays
// The native complex data type is ideal for this application
// Defines f(x)
extern double complex *real_space;
// Holds fourier transform, F(u)
extern double complex *freq_space;
