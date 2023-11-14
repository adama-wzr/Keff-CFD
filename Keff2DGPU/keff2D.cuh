#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <stdbool.h>
#include <fstream>
#include <cfloat>
#include <set>
#include <string>


typedef struct{
  double TCsolid;
  double TCfluid;
  int MeshIncreaseX;
  int MeshIncreaseY;
  double TempLeft;
  double TempRight;
  long int MAX_ITER;
  double ConvergeCriteria;
  char* inputFilename;
  char* outputFilename;
  int printTmap;
  char* TMapName;
  int printQmap;
  char* QMapName;
  int verbose;
  int numCores;
  int MeshFlag;
  int MaxMesh;
  int BatchFlag;
  int NumImg;
}options;


typedef struct{
    int Width;
	int Height;
	int nChannels;
	float porosity;
	float gpuTime;
	unsigned char *target_data;
	float keff;
	bool PathFlag;
	float conv;
    int numCellsX;
    int numCellsY;
    int nElements;
    float dx;
    float dy;
}simulationInfo;


int printOptions(options* opts){
	if(opts->MeshFlag == 0 && opts->BatchFlag == 0){
		printf("--------------------------------------\n\n");
		printf("Current selected options:\n\n");
		printf("--------------------------------------\n");
		printf("Number of Threads Allocated: %d\n", opts->numCores);
		printf("TC Fluid = %.2lf\n", opts->TCfluid);
		printf("TC Solid = %.2lf\n", opts->TCsolid);
		printf("Temperature Left = %.2lf\n", opts->TempLeft);
		printf("Temperature Right = %.2lf\n", opts->TempRight);
		printf("Mesh Amp. X = %d\n", opts->MeshIncreaseX);
		printf("Mesh Amp. Y = %d\n", opts->MeshIncreaseY);
		printf("Maximum Iterations = %ld\n", opts->MAX_ITER);
		printf("Convergence = %.10lf\n", opts->ConvergeCriteria);
		printf("Name of input image: %s\n", opts->inputFilename);
		printf("Name of output file: %s\n", opts->outputFilename);

		if(opts->printTmap == 0){
			printf("Print Temperature Map = False\n");
		} else{
			printf("Temperature Map Name = %s\n", opts->TMapName);
		}
		if(opts->printQmap == 0){
			printf("Print Heat Flux Map = False\n");
		} else{
			printf("Heat Flux Map Name = %s\n", opts->QMapName);
		}
		printf("--------------------------------------\n\n");
	} else if(opts->MeshFlag == 1 && opts->BatchFlag == 0){
		printf("--------------------------------------\n\n");
		printf("Mesh Convergence Test:\n\n");
		printf("Output of Mesh convergence Test: %s\n", opts->outputFilename);
		printf("Number of Threads Allocated: %d\n", opts->numCores);
		printf("TC Fluid = %.2lf\n", opts->TCfluid);
		printf("TC Solid = %.2lf\n", opts->TCsolid);
		printf("Temperature Left = %.2lf\n", opts->TempLeft);
		printf("Temperature Right = %.2lf\n", opts->TempRight);
		printf("Maximum Iterations = %ld\n", opts->MAX_ITER);
		printf("Convergence = %.10lf\n", opts->ConvergeCriteria);
		printf("Maximum Mesh increase: %d\n", opts->MaxMesh);
		printf("--------------------------------------\n\n");
	} else if(opts->MeshFlag == 0 && opts->BatchFlag == 1){
		printf("--------------------------------------\n\n");
		printf("Running Image Batch:\n\n");
		printf("Number of Threads Allocated: %d\n", opts->numCores);
		printf("TC Fluid = %.2lf\n", opts->TCfluid);
		printf("TC Solid = %.2lf\n", opts->TCsolid);
		printf("Temperature Left = %.2lf\n", opts->TempLeft);
		printf("Temperature Right = %.2lf\n", opts->TempRight);
		printf("Mesh Amp. X = %d\n", opts->MeshIncreaseX);
		printf("Mesh Amp. Y = %d\n", opts->MeshIncreaseY);
		printf("Maximum Iterations = %ld\n", opts->MAX_ITER);
		printf("Convergence = %.10lf\n", opts->ConvergeCriteria);
		printf("Name of output file: %s\n", opts->outputFilename);
		printf("Number of files to run: %d\n", opts->NumImg);
		if (opts->printTmap == 1){
			printf("Printing Temperature Distribution for all images.\n");
		} else{
			printf("No temperature maps will be printed.\n");
		}
		if (opts->printQmap == 1){
			printf("Printing Heat Flux Distribution for all images.\n");
		} else{
			printf("No Heat Flux maps will be printed.\n");
		}
		printf("--------------------------------------\n\n");
	} else{
		printf("Options entered are not valid, code will exit.\n");
		return 1;
	}	

	return 0;
}

int readInputFile(char* FileName, options* opts){

	/*
		readInputFile Function:
		Inputs:
			- FileName: pointer to where the input file name is stored.
			- struct options: pass a struct with the options.
		Outputs: None

		Function reads the input file and stores the options in the opts struct.
	*/

	std::string myText;

	char tempC[1000];
	double tempD;
	char tempFilenames[1000];
	std::ifstream InputFile(FileName);

	// initialize the pointers so they are not random

	opts->inputFilename=(char*)malloc(1000*sizeof(char));
	opts->outputFilename=(char*)malloc(1000*sizeof(char));
	opts->TMapName=(char*)malloc(1000*sizeof(char));
	opts->QMapName=(char*)malloc(1000*sizeof(char));
	while(std::getline(InputFile, myText)){

	 	sscanf(myText.c_str(), "%s %lf", tempC, &tempD);
	 	if (strcmp(tempC, "ks:") == 0){
	 		opts->TCsolid = tempD;
	 	}else if(strcmp(tempC, "kf:") == 0){
	 		opts->TCfluid = tempD;

	 	}else if(strcmp(tempC, "MeshAmpX:") == 0){
	 		opts->MeshIncreaseX = (int)tempD;

	 	}else if(strcmp(tempC, "MeshAmpY:") == 0){
	 		opts->MeshIncreaseY = (int)tempD;

	 	}else if(strcmp(tempC, "InputName:") == 0){
	 		sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
	 		strcpy(opts->inputFilename, tempFilenames);

	 	}else if(strcmp(tempC, "TR:") == 0){
	 		opts->TempRight = tempD;

	 	}else if(strcmp(tempC, "TL:") == 0){
	 		opts->TempLeft = tempD;

	 	}else if(strcmp(tempC, "OutputName:") == 0){
	 		sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
	 		strcpy(opts->outputFilename, tempFilenames);

	 	}else if(strcmp(tempC, "printTMap:") == 0){
	 		opts->printTmap = (int)tempD;

	 	}else if(strcmp(tempC, "TMapName:") == 0){
	 		sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
	 		strcpy(opts->TMapName, tempFilenames);

	 	}else if(strcmp(tempC, "printQMap:") == 0){
	 		opts->printQmap = (int)tempD;

	 	}else if(strcmp(tempC, "QMapName:") == 0){
	 		sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
	 		strcpy(opts->QMapName, tempFilenames);

	 	}else if(strcmp(tempC, "Convergence:") == 0){
	 		opts->ConvergeCriteria = tempD;

	 	}else if(strcmp(tempC, "MaxIter:") == 0){
	 		opts->MAX_ITER = (long int)tempD;

	 	}else if(strcmp(tempC, "Verbose:") == 0){
	 		opts->verbose = (int)tempD;
	 	}else if(strcmp(tempC, "NumCores:") == 0){
	 		if(tempD<1){
	 			printf("Entered number of cores not valid.\n");
	 			printf("Default = 1.");
	 			opts->numCores = 1;
	 		}else{
	 			opts->numCores = (int)tempD;
	 		}
	 	} else if(strcmp(tempC, "MeshFlag:") == 0){
	 		opts->MeshFlag = (int)tempD;
	 	} else if(strcmp(tempC, "MaxMesh:") == 0){
	 		opts->MaxMesh = (int)tempD;
	 	} else if(strcmp(tempC, "RunBatch:") == 0){
	 		opts->BatchFlag = (int)tempD;
	 	} else if(strcmp(tempC, "NumImages:") == 0){
	 		opts->NumImg = (int)tempD;
	 	}
	}
	
	InputFile.close();



	if(opts->verbose == 1){
		printOptions(opts);
	} else if(opts->verbose != 0){
		printf("Please enter a value of 0 or 1 for 'verbose'. Default = 0.\n");
	}
	return 0;
}


float WeightedHarmonicMean(float w1, float w2, float x1, float x2){
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
	float H = (w1 + w2)/(w1/x1 + w2/x2);
	return H;
}


int DiscretizeMatrix2D(float* K, float* A, float* b, simulationInfo simInfo, options opts){
	/*
		DiscretizeMatrix2D

		Inputs:
			- pointer to float array K, where local diffusion coefficients are stored
			- pointer to empty coefficient matrix
			- pointer to RHS of the system of equations
			- datastructure containing simulation information
			- datastructure with the user-entered options
		Output:
			- None

			Function creates the CoeffMatrix and RHS of the system of equations and stores them
			on the appropriate memory spaces.
	*/

	int index;
	float dxw, dxe, dys, dyn;
	float kw, ke, ks, kn;

	float dx, dy;
	dx = simInfo.dx;
	dy = simInfo.dy;

	for(int i = 0; i<simInfo.numCellsY; i++){
		for(int j = 0; j<simInfo.numCellsX; j++){
			// initialize everything to zeroes
			index = (i*simInfo.numCellsX + j); 
			b[index] = 0;
			for(int k = 0; k<5; k++){
				A[index*5 + k] = 0;
			}
			// left boundary, only P and E
			if (j == 0){
				dxe = dx;
				ke = WeightedHarmonicMean(dxe/2,dxe/2, K[index], K[index+1]);
				dxw = dx/2;
				kw = K[index];
				A[index*5 + 2] = -ke*dy/dxe;
				A[index*5 + 0] += (ke*dy/dxe + kw*dy/dxw);
				b[index] += opts.TempLeft*kw*dy/dxw;
			} else if(j == simInfo.numCellsX - 1){		// Right boundary, only P and W
				dxw = dx;
				kw = WeightedHarmonicMean(dxw/2,dxw/2, K[index], K[index-1]);
				dxe = dx/2;
				ke = K[index];
				A[index*5 + 1] = -kw*dy/dxw;
				A[index*5 + 0] += (ke*dy/dxe + kw*dy/dxw);
				b[index] += opts.TempRight*ke*dy/dxe;
			} else{								// P, W, and E
				dxw = dx;
				kw = WeightedHarmonicMean(dxw/2,dxw/2, K[index], K[index-1]);
				dxe = dx;
				ke = WeightedHarmonicMean(dxe/2,dxe/2, K[index], K[index+1]);
				A[index*5 + 1] = -kw*dy/dxw;
				A[index*5 + 2] = -ke*dy/dxe;
				A[index*5 + 0] += (ke*dy/dxe + kw*dy/dxw);
			}
			// top boundary, only S and P
			if (i == 0){
				dyn = dy/2;
				kn = K[index];
				dys = dy;
				ks = WeightedHarmonicMean(dys/2, dys/2, K[index + simInfo.numCellsX], K[index]);
				A[index*5 + 3] = -ks*dx/dys;
				A[index*5 + 0] += (ks*dx/dys);
			}else if(i == simInfo.numCellsY - 1){
				dyn = dy;
				kn = WeightedHarmonicMean(dyn/2, dyn/2, K[index], K[index - simInfo.numCellsX]);
				dys = dy/2;
				ks = K[index];
				A[index*5 + 4] = -kn*dx/dyn;
				A[index*5 + 0] += kn*dx/dyn;
			} else{
				dyn = dy;
				kn = WeightedHarmonicMean(dyn/2, dyn/2, K[index], K[index - simInfo.numCellsX]);
				dys = dy;
				ks = WeightedHarmonicMean(dys/2, dys/2, K[index + simInfo.numCellsX], K[index]);
				A[index*5 + 3] = -ks*dx/dys;
				A[index*5 + 4] = -kn*dx/dyn;
				A[index*5 + 0] += (kn*dx/dyn + ks*dx/dys);
			}
		}
	}
	return 0;
}


int initializeGPU(float **d_x_vec, float **d_temp_x_vec, float **d_RHS, float **d_Coeff, simulationInfo simInfo){

	// Set device, when cudaStatus is called give status of assigned device.
	// This is important to know if we are running out of GPU space
	cudaError_t cudaStatus = cudaSetDevice(0);

	// Start by allocating space in GPU memory

	if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		getchar();
        return 0;
    }

    cudaStatus = cudaMalloc((void**)&(*d_x_vec), simInfo.nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		getchar();
        return 0;
    }

    cudaStatus = cudaMalloc((void**)&(*d_temp_x_vec), simInfo.nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		getchar();
        return 0;
    }

    cudaStatus = cudaMalloc((void**)&(*d_RHS), simInfo.nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		getchar();
        return 0;
    }

    cudaStatus = cudaMalloc((void**)&(*d_Coeff), simInfo.nElements*sizeof(float)*5);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		getchar();
        return 0;
    }

    // Set GPU buffers (initializing matrices to 0)

     // Memset GPU buffers
    cudaStatus = cudaMemset((*d_x_vec),0, simInfo.nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
		getchar();
        return 0;
    }

	// Memset GPU buffers
    cudaStatus = cudaMemset((*d_temp_x_vec),0, simInfo.nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
		getchar();
        return 0;
    }

     // Memset GPU buffers
    cudaStatus = cudaMemset((*d_RHS),0, simInfo.nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
		getchar();
        return 0;
    }

	// Memset GPU buffers
    cudaStatus = cudaMemset((*d_Coeff),0, 5*simInfo.nElements*sizeof(float));		// coefficient matrix has the 5 main diagonals for all elements
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
		getchar();
        return 0;
    }

    return 1;
}

void unInitializeGPU(float **d_x_vec, float **d_temp_x_vec, float **d_RHS, float **d_Coeff)
{
	cudaError_t cudaStatus;

	if((*d_x_vec)!=NULL)
    cudaStatus = cudaFree((*d_x_vec));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFree failed!");
        return;
    }

	if((*d_temp_x_vec)!=NULL)
    cudaStatus = cudaFree((*d_temp_x_vec));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFree failed!");
        return;
    }

	if((*d_Coeff)!=NULL)
    cudaStatus = cudaFree((*d_Coeff));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFree failed!");
        return;
    }

	if((*d_RHS)!=NULL)
    cudaStatus = cudaFree((*d_RHS));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFree failed!");
        return;
    }    

	cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
		getchar();
        return;
    }
}