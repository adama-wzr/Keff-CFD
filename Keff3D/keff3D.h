#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <stdbool.h>
#include <fstream>



int readImage3D(unsigned char* structureAddress, int width, int height, int depth, char* filename){
	/*
		readImage3D Function:
		Inputs:
			- structureAddress: pointer to unsigned char where the coordinates of solid pixels will be saved.
			- width: integer, width of the original domain (x-direction).
			- height: integer, height of the original domain (y-direction).
			- depth: integer, depth of the original domain (z-direction).
			- filename: pointer to char containing name of the file to read.
		Outputs:
			- None

		Function will read the structure with the given number of voxels in each direction and store it in the 
		allocated array.
	*/

	int *x = (int *)malloc(sizeof(int)*width*height*depth);
    int *y = (int *)malloc(sizeof(int)*width*height*depth);
    int *z = (int *)malloc(sizeof(int)*width*height*depth);

    FILE *target_data;

	target_data = fopen(filename, "r");

	size_t count = 0;

	if (target_data == NULL) {
        fprintf(stderr, "Error reading file\n");
        return 1;
    }

    while(fscanf(target_data, " %d,%d,%d", &x[count], &y[count], &z[count]) == 3) {
		count++;
    }

	fclose(target_data);

	for(int i = 0; i< count; i++){
		structureAddress[z[i]*width*height + y[i]*width + x[i]] = 1;
	}

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


int DiscretizeMatrixCD3D(double* K, int numRows, int numCols, int numSlices, double* A, double* b, double tL, double tR, double* xCenter, double* yCenter, double* zCenter){

	/*
	
	Coefficient Matrix notation:
	0: P
	1: W
	2: E
	3: S
	4: N
	5: F
	6: B
	*/

	int index;
	double kw, ke, ks, kn, kb, kf;
	double dx, dy, dz;

	dx = (double)1.0/numCols;
	dy = (double)1.0/numRows;
	dz = (double)1.0/numSlices;

	for(int i = 0; i<numSlices; i++){
		for(int j = 0; j<numRows; j++){
			for(int k = 0; k<numCols; k++){
				index = i*numCols*numRows + j*numCols + k;
			}
		}
	}

	return 0;
}