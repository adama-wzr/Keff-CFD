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

int TDMA(int n, double *A, double *b, double *x, int rowNumber){
	/*
	Function TDMA:
	Inputs:
		- n, integer: size of matrix A (nxn)
		- *A: pointer to the A matrix on a Ax = b system.
		- *b: pointer to the RHS of the system, b, in Ax=b
		- *x: pointer the values to be modifyed in the Ax = b system.
		- rowNumber: number of the current row being solver.

	Outputs: None

	Function uses the Thomas algorithm or TriDiagonal matrix algorithm
	to solve for x in the Ax = b system
	*/
	double m;

	// with this method, solution is direct.
	// Analogous to Gaussian Elimination, but only with 3 diagonals.

	for(int i = 1; i<n; i++){
		m = A[i*n + i - 1]/A[(i-1)*n + i - 1];
		A[i*n + i] = A[i*n + i] - m*A[(i-1)*n + i];
		b[i] -= m*b[i-1];
	}

	// The solution to the last value of x is simply b[n]/A[n,n]
	// but with the indexing starting at 0
	x[rowNumber*n + n-1] = b[n-1]/A[n*n - 1];

	// now work the way back through the indexes solving it.
	// Remember A is a upper triangular matrix at this stage.

	for(int i = n-2; i >= 0; i--){
		x[rowNumber*n + i] = (b[i] - A[i*n + i + 1]*x[rowNumber*n + i + 1])/A[i*n + i];
	}

	return 0;
}


int TDMAscanner(double *A, double *b, double *X, int maxIter, double Threshold, int numCols, int numRows, double *qR, double *qL,
	double tR, double tL, double *K, double *XC, double *YC){

	/*
	Function TDMAscanner:
	Inputs:
		- *A: pointer to the discretization matrix, size (numCols*numRows, 5)
		- *b: pointert to the RHS, size (numCols*numRows)
		- *X: pointer to the solution vector, size (numCols*numRows)
		- maxIter: integer with the maximum number of allowed iterations
		- Threshold: value of convergence criteria
		- numCols: numCols of the domain
		- numRows: numRows of the domain
		- *qR: pointer a vector of size (numRows) containing heat flux from each cell
			at the right side of the domain.
		- *qL: pointer a vector of size (numRows) containing heat flux from each cell
			at the left side of the domain.
		- tR: temperature in the right boundary.
		- tL: temperature in the left boundary.
		- *K: vector containing the thermal conductivity of each cell in the domain.
		- *XC: pointer to vector containing the x-coordinate of the center of each cell (length numCols)
		- *YC: pointer to vector containing the y-coordinate of the center of each cell (length numRows)

	Outputs: None

	Function extends the functionality of the regular TDMA to 5 diagonal systems, using 
	the iterative TDMA algorithm. Basically it scans one row at a time, taken the others as knowns.
	Modifies the vector X.

	*/

	// Create TDMA_A and TDMA_b matrices to be modified by the TDMA function, subsets of the
	// entire A and b system.

	double *TDMA_A = (double *)malloc(sizeof(double)*numCols*numCols);
	double *TDMA_b = (double *)malloc(sizeof(double)*numCols);

	// initialize to 0

	for(int i = 0; i<numCols; i++){
		TDMA_b[i] = 0;
		for(int j = 0; j<numCols; j++){
			TDMA_A[i*numCols + j] = 0;
		}
	}

	// Create 3 vectors to store and modify the values that will become part of TDMA_A and TDMA_b

	double *upperDiag = (double *)malloc(sizeof(double)*(numCols-1));
	double *lowerDiag = (double *)malloc(sizeof(double)*(numCols-1));
	double *mainDiag = (double *)malloc(sizeof(double)*(numCols));

	// more helper variables

	int iterCount = 0;
	double keffOld = 1;
	double keffNew = 1;
	int iterToCheck = 500;
	double qAvg = 0;
	double Q1 = 0;
	double Q2 = 0;

	double percentChange = 100;

	while(iterCount < maxIter && Threshold < percentChange){
		for(int i = 0; i<numRows; i++){
			for(int j = 0; j<numCols; j++){
				if(j == 0){ // no lower diagonal
					upperDiag[j] = A[(i)*numCols*5 + j*5 + 2];
					mainDiag[j] = A[(i)*numCols*5 + j*5 + 0];
					TDMA_A[j*numCols + j] = mainDiag[j];
					TDMA_A[j*numCols + j + 1] = upperDiag[j];
					TDMA_b[j] = b[i*numCols + j];
				}else if(j == numCols-1){ 	// no upper diagonal
					lowerDiag[j-1] = A[i*numCols*5 + j*5 + 1];
					mainDiag[j] = A[i*numCols*5 + j*5 + 0];
					TDMA_A[j*numCols + j] = mainDiag[j];
					TDMA_A[j*numCols + j - 1] = lowerDiag[j-1];
					TDMA_b[j] = b[i*numCols + j]; 
				} else{
					lowerDiag[j-1] = A[i*numCols*5 + j*5 + 1];
					mainDiag[j] = A[i*numCols*5 + j*5 + 0];
					upperDiag[j] = A[(i)*numCols*5 + j*5 + 2];
					TDMA_A[j*numCols + j] = mainDiag[j];
					TDMA_A[j*numCols + j + 1] = upperDiag[j];
					TDMA_A[j*numCols + j - 1] = lowerDiag[j-1];
					TDMA_b[j] = b[i*numCols + j];
				}
			}

			for(int j = 0; j<numCols; j++){
				if(i == 0){
					TDMA_b[j] -= A[i*numCols*5 + j*5 + 3]*X[(i+1)*numCols + j];
				} else if(i == numRows - 1){
					TDMA_b[j] -= A[i*numCols*5 + j*5 + 4]*X[(i-1)*numCols + j];
				} else{
					TDMA_b[j] -= A[i*numCols*5 + j*5 + 3]*X[(i+1)*numCols + j] + A[i*numCols*5 + j*5 + 4]*X[(i-1)*numCols + j];
				}

			}
			TDMA(numCols, TDMA_A, TDMA_b, X, i);
		}
		iterCount++;

		if (iterCount % iterToCheck == 0){
			for (int j = 0; j<numRows; j++){
				double dy = YC[0]*2;
				qL[j] = K[j*numRows]*dy*(X[j*numCols] - tL)/(XC[0]);
				qR[j] = K[(j + 1)*numRows - 1]*dy*(tR - X[(j+1)*numCols -1])/(1 - XC[numCols - 1]);
				Q1 += qL[j];
				Q2 += qR[j];
			}
			Q1 = Q1;
			Q2 = Q2;
			qAvg = (Q1 + Q2)/2;
			keffNew = qAvg/((tR - tL)/numRows);
			percentChange = fabs((keffNew - keffOld)/keffOld);
			keffOld = keffNew;

			if (percentChange < 0.01){
				iterToCheck = 100;
			} else if(percentChange < 0.001){
				iterToCheck = 10;
			}
		}
	}

	free(mainDiag);
	free(upperDiag);
	free(lowerDiag);

	return iterCount;
}