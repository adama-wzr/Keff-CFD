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



int readImage(unsigned char** imageAddress, int* Width, int* Height, int* NumOfChannels, char* ImageName){
	/*
		readImage Function:
		Inputs:
			- imageAddress: unsigned char reference to the pointer in which the image will be read to.
			- Width: pointer to variable to store image width
			- Height: pointer to variable to store image height
			- NumofChannels: pointer to variable to store number of channels in the image.
					-> NumofChannels has to be 1, otherwise code is terminated. Please enter grayscale
						images with NumofChannels = 1.
		Outputs: None

		Function reads the image into the pointer to the array to store it.
	*/

	*imageAddress = stbi_load(ImageName, Width, Height, NumOfChannels, 1);

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
	double dxw, dxe, dys, dyn;
	double kw, ke, ks, kn;

	// Remove following lines of code if not structured steady mesh:

	double dx, dy;
	dx = (double)1.0/numCols;
	dy = (double)1.0/numRows;

	for(int i = 0; i<numRows; i++){
		for(int j = 0; j<numCols; j++){
			// initialize everything to zeroes
			index = (i*numCols + j); 
			b[index] = 0;
			for(int k = 0; k<5; k++){
				A[index*5 + k] = 0;
			}
			// left boundary, only P and E
			if (j == 0){
				dxe = xCenter[j+1] - xCenter[j];
				ke = WeightedHarmonicMean(dxe/2,dxe/2, K[index], K[index+1]);
				dxw = xCenter[j];
				kw = K[index];
				A[index*5 + 2] = -ke*dy/dxe;
				A[index*5 + 0] += (ke*dy/dxe + kw*dy/dxw);
				b[index] += tL*kw*dy/dxw;
			} else if(j == numCols - 1){		// Right boundary, only P and W
				dxw = xCenter[j] - xCenter[j - 1];
				kw = WeightedHarmonicMean(dxw/2,dxw/2, K[index], K[index-1]);
				dxe = 1 - xCenter[j];
				ke = K[index];
				A[index*5 + 1] = -kw*dy/dxw;
				A[index*5 + 0] += (ke*dy/dxe + kw*dy/dxw);
				b[index] += tR*ke*dy/dxe;
			} else{								// P, W, and E
				dxw = xCenter[j] - xCenter[j - 1];
				kw = WeightedHarmonicMean(dxw/2,dxw/2, K[index], K[index-1]);
				dxe = xCenter[j+1] - xCenter[j];
				ke = WeightedHarmonicMean(dxe/2,dxe/2, K[index], K[index+1]);
				A[index*5 + 1] = -kw*dy/dxw;
				A[index*5 + 2] = -ke*dy/dxe;
				A[index*5 + 0] += (ke*dy/dxe + kw*dy/dxw);
			}
			// top boundary, only S and P
			if (i == 0){
				dyn = yCenter[i];
				kn = K[index];
				dys = yCenter[i + 1] - yCenter[i];
				ks = WeightedHarmonicMean(dys/2, dys/2, K[index + numCols], K[index]);
				A[index*5 + 3] = -ks*dx/dys;
				A[index*5 + 0] += (ks*dx/dys);
			}else if(i == numRows - 1){
				dyn = yCenter[i] - yCenter[i-1];
				kn = WeightedHarmonicMean(dyn/2, dyn/2, K[index], K[index - numCols]);
				dys = 1 - yCenter[i];
				ks = K[index];
				A[index*5 + 4] = -kn*dx/dyn;
				A[index*5 + 0] += kn*dx/dyn;
			} else{
				dyn = yCenter[i] - yCenter[i-1];
				kn = WeightedHarmonicMean(dyn/2, dyn/2, K[index], K[index - numCols]);
				dys = yCenter[i + 1] - yCenter[i];
				ks = WeightedHarmonicMean(dys/2, dys/2, K[index + numCols], K[index]);
				A[index*5 + 3] = -ks*dx/dys;
				A[index*5 + 4] = -kn*dx/dyn;
				A[index*5 + 0] += (kn*dx/dyn + ks*dx/dys);
			}
		}
	}

	return 0;
}


int JacobiIteration(double *arr, double *sol, double *x_vec, int iter_limit, double tolerance, int size){
	int iteration_count = 0;
	double *temp_x_vec = (double *)malloc(sizeof(double)*size*size);
	double conv_stat = 1;
	double sigma = 0;
	double norm_diff = 1;
	printf("size = %d\n", size);
	for (int i = 0; i<size*size; i++){
		x_vec[i] = 0.5;
		temp_x_vec[i] = x_vec[i];
	}

	while(norm_diff > tolerance && iteration_count < iter_limit){
		for(int i = 0; i<size*size; i++){
			sigma = 0;
			for(int j = 1; j<5; j++){
				if(arr[i*5 + j] != 0){
					if(j == 1){
						sigma += arr[i*5 + j]*temp_x_vec[i - 1];
					} else if(j == 2){
						sigma += arr[i*5 + j]*temp_x_vec[i + 1];
					} else if(j == 3){
						sigma += arr[i*5 + j]*temp_x_vec[i + size];
					} else if(j == 4){
						sigma += arr[i*5 + j]*temp_x_vec[i - size];
					}
				}
			}
			x_vec[i] = 1/arr[i*5 + 0]*(sol[i] - sigma);
			// printf("x_vec[%d] = %f, sigma = %f\n",i,x_vec[i],sigma);
		}

		norm_diff = 0;

		for (int i = 0; i < size*size; i++){
			norm_diff += sqrt((x_vec[i] - temp_x_vec[i])*(x_vec[i] - temp_x_vec[i]));
			temp_x_vec[i] = x_vec[i];
		}
		iteration_count++;
	}
	printf("Norm diff = %f\n",norm_diff);
	free(temp_x_vec);
	return iteration_count;
}