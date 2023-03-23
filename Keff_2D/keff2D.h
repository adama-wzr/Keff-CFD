#include<omp.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include <vector>
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include<stdbool.h>



int readImage(unsigned char* imageAddress, int* Width, int* Height, int* NumOfChannels, char* ImageName){
	/*
		readImage Function:
		Inputs:
			- imageAddress: unsigned char pointer to address to read the image to.
			- Width: pointer to variable to store image width
			- Height: pointer to variable to store image height
			- NumofChannels: pointer to variable to store number of channels in the image.
					-> NumofChannels has to be 1, otherwise code is terminated. Please enter grayscale
						images with NumofChannels = 1.
		Outputs: None

		Function reads the image into the pointer to the array to store it.
	*/

	imageAddress = stbi_load(ImageName, Width, Height, NumOfChannels, 1);

	return 0;
}


double WeightedHarmonicMean(double w1, double w2, double x1, double x2){
	/*
		WeightedHarmonicMean Function:
		Inputs:
			-w1: weight of the first number
			-w2: weight of the second number
			-x1: first number to be averaged
			-x2: second number to be averaged
		Output:
			- returns H, the weighted harmonic mean between x1 and x2, using weights w1 and w2.
	*/
	double H = (w1 + w2)/(w1/x1 + w2/x2);
	return H;
}


int DiscretizeMatrixCD2D(double* K, int numRows, int numCols, double* A, double* b, double tL, double tR, double* xCenter, double* yCenter){
	/*
	DiscretizeMatrixCD_2D Function:
	Inputs:
		- K: Pointer to dynamically allocated array K. Each point in K gives the thermal conudctivity 
			of the point at given index.
		- numRows: Number of rows
		- numCols: Number of columns
		- A: pointer to array to store coefficient matrix
		- b: pointer to array to store the right-hand side of the system Ax = b
		- tR: temperature at the right side of the domain
		- tL: temperature at the left side of the domain
		- xCenter: array containing the x-coordinate of the center 

	Outputs: None

	Function will populate the dynamically allocated vectors A and b on the Ax = b system, all from
	the provided information, using central-differencing approximation for the variables, 2D domain.

	Coefficient Matrix in cardinal indexing notation:
	0: P
	1: W
	2: E
	3: S
	4: N

	*/

	int index;
	double dxw, dxe, dy;
	double kw, ke, ks, kn;

	for(int i = 0; i<numRows; i++){
		for(int j = 0; j<numCols; j++){
			// initialize everything to zeroes
			index = (i*numCols + j); 
			b[index] = 0;
			for(int k = 0; k<5; k++){
				A[index + k] = 0;
			}

			// left boundary, only P and E
			if (j == 0){
				dxe = xCenter[j+1] - xCenter[j];
				ke = WeightedHarmonicMean(dxe/2,dxe/2, K[index], K[index+1]);
				dxw = xCenter[j];
				kw = K[index];
				A[index*5 + 2] = ke/dxe;
				A[index*5 + 0] += -(ke/dxe + kw/dxw);
				b[index] += -tL*kw/dxw;
			} else if(j == numCols - 1){		// Right boundary, only P and W
				dxw = xCenter[index] - xCenter[index - 1];
				kw = WeightedHarmonicMean(dxw/2,dxw/2, K[index], K[index-1]);
				dxe = 1 - xCenter[index];
				ke = K[index];
				A[index*5 + 1] = kw/dxw;
				A[index*5 + 0] += -(ke/dxe + kw/dwx);
				b[index] += -tR*ke/dxe;
			} else{								// P, W, and E
				dxw = xCenter[index] - xCenter[index - 1];
				kw = WeightedHarmonicMean(dxw/2,dxw/2, K[index], K[index-1]);
				dxe = xCenter[j+1] - xCenter[j];
				ke = WeightedHarmonicMean(dxe/2,dxe/2, K[index], K[index+1]);
				A[index*5 + 1] = kw/dxw;
				A[index*5 + 2] = ke/dxe;
				A[index*5 + 0] = -(ke/dxe + kw/dxw);
			}
		}
	}

}