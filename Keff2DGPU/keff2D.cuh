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
  float TCsolid;
  float TCfluid;
  int MeshIncreaseX;
  int MeshIncreaseY;
  float TempLeft;
  float TempRight;
  long int MAX_ITER;
  float ConvergeCriteria;
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
	float conv;
    int numCellsX;
    int numCellsY;
    int nElements;
    float dx;
    float dy;
}simulationInfo;


__global__ void updateX_V1(float* A, float* x, float* b, float* xNew, int nElements, int numCellsX)
{
	unsigned int myRow = blockIdx.x * blockDim.x + threadIdx.x;

	if (myRow < nElements){
		float sigma = 0;
		for(int j = 1; j<5; j++){
			if(A[myRow*5 + j] !=0){
				if(j == 1){
					sigma += A[myRow*5 + j]*x[myRow - 1];
				} else if(j == 2){
					sigma += A[myRow*5 + j]*x[myRow + 1];
				} else if(j == 3){
					sigma += A[myRow*5 + j]*x[myRow + numCellsX];
				} else if(j == 4){
					sigma += A[myRow*5 + j]*x[myRow - numCellsX];
				}
			}
		}
		xNew[myRow] = 1/A[myRow*5 + 0] * (b[myRow] - sigma);
	}
		
}


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
	float tempD;
	char tempFilenames[1000];
	std::ifstream InputFile(FileName);

	// initialize the pointers so they are not random

	opts->inputFilename=(char*)malloc(1000*sizeof(char));
	opts->outputFilename=(char*)malloc(1000*sizeof(char));
	opts->TMapName=(char*)malloc(1000*sizeof(char));
	opts->QMapName=(char*)malloc(1000*sizeof(char));
	while(std::getline(InputFile, myText)){

	 	sscanf(myText.c_str(), "%s %f", tempC, &tempD);
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

int readImage(options opts, simulationInfo* myImg){
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

	myImg->target_data = stbi_load(opts.inputFilename, &myImg->Width, &myImg->Height, &myImg->nChannels, 1);

	return 0;
}


float calcPorosity(unsigned char* imageAddress, int Width, int Height){
	/*
		calcPorosity
		Inputs:
			- imageAddress: pointer to the read image.
			- Width: original width from std_image
			- Height: original height from std_image

		Output:
			- porosity: float containing porosity.

		Function calculates porosity by counting pixels.
	*/

	float totalCells = (float)Height*Width;
	float porosity = 0;
	for(int i = 0; i<Height; i++){
		for(int j = 0; j<Width; j++){
			if(imageAddress[i*Width + j] < 150){
				porosity += 1.0/totalCells;
			}
		}
	}

	return porosity;
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

int outputSingle(options opts, simulationInfo simInfo){
	/*
		outputSingle:
		Inputs:
			- struct opts: data structure containing user entered options
			- struct simInfo: data structure containing simulation domain information
				and results information
		Outputs:
			- none
		
		Function will create a file or append an existing file save information related to a single simulation.
	*/
	FILE *OUTPUT;

	// imgNum, porosity,keff,Time,nElements,converge,ks,kf

	OUTPUT = fopen(opts.outputFilename, "a+");
	fprintf(OUTPUT,"imgNum,porosity,keff,Time,nElements,converge,ks,kf\n");
	fprintf(OUTPUT, "%s,%f,%2.3f,%f,%d,%f,%f,%f\n", opts.inputFilename, simInfo.porosity, simInfo.keff, simInfo.gpuTime/1000, simInfo.nElements, simInfo.conv,
		opts.TCsolid, opts.TCfluid);
	fclose(OUTPUT);
	printf("Final Keff = %2.3f\n", simInfo.keff);
	return 0;
}


int printTMAP(options* o, float* x, int numRows, int numCols, int imgNum){
	/*
		printTMAP:
		Inputs:
			- o -> options datastructure
			- x -> pointer to temperature distribution map.
			- numRows -> number of rows
			- numCols -> number of columns
			- imgNum: this option is only useful when running a batch, this is the image number.
		Outputs:
			- none
		Function creates and saves a temperature map onto a .csv file
	*/

	FILE *T_OUT;
	char filename[100];
	if(o->BatchFlag == 0){
		strcpy(filename, o->TMapName);
	} else{
		sprintf(filename, "TMAP_%05d.csv\n", imgNum);
	}
  	T_OUT = fopen(filename, "w+");
  	fprintf(T_OUT,"x,y,T\n");

	for(int i = 0; i<numRows; i++){
		for(int j = 0; j<numCols; j++){
			fprintf(T_OUT,"%d,%d,%f\n",j,i,x[i*numCols + j]);
		}
	}

	fclose(T_OUT);
	return 0;
}


int printQMAP(options* o, float* x, float* K, int numRows, int numCols, float dx, float dy, float* qR, float* qL, int imgNum){

	/*
		printQMAP:
		Inputs:
			- o -> options datastructure
			- x -> pointer to temperature distribution map.
			- numRows -> number of rows
			- numCols -> number of columns
			- float dx -> size of one grid element (for steady regular grids)
			- float dy -> size of one grid element (for steady regular grids)
			- qL = pointer to array heat transfer across the left boundary
			- qR = pointer to array with heat transfer across the right boundary
			- imgNum: this option is only useful when running a batch, this is the image number.
		Outputs:
			- none
		Function creates and saves a heat flux map onto a .csv file.
	*/


	FILE *Q_OUT;

	char filename[100];
	if(o->BatchFlag == 0){
		strcpy(filename, o->QMapName);
	} else{
		sprintf(filename, "QMAP_%05d.csv\n",imgNum);
	}
	Q_OUT = fopen(filename, "w+");
	fprintf(Q_OUT,"x,y,Q\n");

	float localQ;
	float localK;

	for(int i = 0; i<numRows; i++){
		for(int j = 0; j<numCols + 1; j++){
			if(j == 0){
				fprintf(Q_OUT,"%d,%d,%f\n",j,i, qL[i]);
			} else if(j == numCols){
				fprintf(Q_OUT, "%d,%d,%f\n",j,i, qR[i]);
			} else{
				localK = WeightedHarmonicMean(dx/2, dx/2, K[i*numCols + j], K[i*numCols + j - 1]);
				localQ = dy/dx*localK*(x[i*numCols + j] - x[i*numCols + j - 1]);
				fprintf(Q_OUT, "%d,%d,%f\n",j,i, localQ);
			}
		}
	}

	fclose(Q_OUT);
	return 0;
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


int JacobiGPU(float *arr, float *sol, float *x_vec, float *temp_x_vec, options opts,
	float *d_x_vec, float *d_temp_x_vec, float *d_Coeff, float *d_RHS, float *QL, float *QR, float *K, simulationInfo* simInfo)
{

	int iterCount = 0;
	float percentChange = 1;
	int threads_per_block = 160;
	int numBlocks = simInfo->nElements/threads_per_block + 1;
	float keffOld = 1;
	float keffNew = 1;
	int iterToCheck = 1000;
	float Q1,Q2;
	float qAvg = 0;
	float dx,dy;
	int numRows = simInfo->numCellsY;
	int numCols = simInfo->numCellsX;
	const char *str = (char*) malloc(1024); // To store error string

	dx = simInfo->dx;
	dy = simInfo->dy;

	int nRows = simInfo->nElements;	// number of rows in the coefficient matrix
	int nCols = 5;							// number of cols in the coefficient matrix

	// Initialize temp_x_vec

	for(int i = 0; i<nRows; i++){
		temp_x_vec[i] = x_vec[i];
	}

	//Copy arrays into GPU memory

	cudaError_t cudaStatus = cudaMemcpy(d_temp_x_vec, temp_x_vec, sizeof(float) * nRows, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "temp_x_vec cudaMemcpy failed!");
		str = cudaGetErrorString(cudaStatus);
		fprintf(stderr, "CUDA Error!:: %s\n", str);
	}
	cudaStatus = cudaMemcpy(d_RHS, sol, sizeof(float)*nRows, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "d_RHS cudaMemcpy failed!");
		str = cudaGetErrorString(cudaStatus);
		fprintf(stderr, "CUDA Error!:: %s\n", str);
	}
	cudaStatus = cudaMemcpy(d_Coeff, arr, sizeof(float)*nRows*nCols, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "d_Coeff cudaMemcpy failed!");
		str = cudaGetErrorString(cudaStatus);
		fprintf(stderr, "CUDA Error!:: %s\n", str);
	}

	// Declare event to get time
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start, 0);

	while(iterCount < opts.MAX_ITER && opts.ConvergeCriteria < percentChange)
	{
		// Call Kernel to Calculate new x-vector
		
		updateX_V1<<<numBlocks, threads_per_block>>>(d_Coeff, d_temp_x_vec, d_RHS, d_x_vec, simInfo->nElements, simInfo->numCellsX);

		// update x vector

		d_temp_x_vec = d_x_vec;

		// Convergence related material

		if (iterCount % iterToCheck == 0){
			cudaStatus = cudaMemcpy(x_vec, d_x_vec, sizeof(float) * nRows, cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "x_vec cudaMemcpy failed!");
				str = cudaGetErrorString(cudaStatus);
				fprintf(stderr, "CUDA Error!:: %s\n", str);
			}
			Q1 = 0;
			Q2 = 0;
			for (int j = 0; j<numRows; j++){
				QL[j] = K[j*numRows]*dy*(x_vec[j*numCols] - opts.TempLeft)/(dx/2);
				QR[j] = K[(j + 1)*numRows - 1]*dy*(opts.TempRight - x_vec[(j+1)*numCols -1])/(dx/2);
				// printf("T(0,%d) = %2.3f\n", j, x_vec[j*numCols]);
				Q1 += QL[j];
				Q2 += QR[j];
			}
			Q1 = Q1;
			Q2 = Q2;
			qAvg = (Q1 + Q2)/2;
			keffNew = qAvg/((opts.TempRight - opts.TempLeft));
			percentChange = fabs((keffNew - keffOld)/keffOld);
			keffOld = keffNew;

			printf("Iteration = %d, Keff = %2.3f\n", iterCount, keffNew);

			if (percentChange < 0.001){
				iterToCheck = 100;
			} else if(percentChange < 0.0001){
				iterToCheck = 10;
			}
			simInfo->conv = percentChange;
		}

		// Update iteration count
		iterCount++;
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);

	cudaStatus = cudaMemcpy(x_vec, d_x_vec, sizeof(float)*nRows, cudaMemcpyDeviceToHost);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "x_vec cudaMemcpy failed!");
		str = cudaGetErrorString(cudaStatus);
		fprintf(stderr, "CUDA Error!:: %s\n", str);
	}

	simInfo->keff = keffNew;

	simInfo->gpuTime += elapsedTime;

	return iterCount;
}


int SingleSim(options opts){
	/*
		Function to read a single image and simulate the effective diffusivity. Results
		 are stored on the output file

		Inputs:
			Datastructure with user-defined simulation options
		Outputs:
			none
	*/

	// Define data structures

	simulationInfo simInfo;

	// first step is to read the image properly and calculate the porosity

	readImage(opts, &simInfo);

	simInfo.porosity = calcPorosity(simInfo.target_data, simInfo.Width, simInfo.Height);

	// right now the program only deals with grayscale binary images, so we need to make sure to return that to the user

	if(opts.verbose == 1){
		std::cout << "Width = " << simInfo.Width << " Height = " << simInfo.Height << " Channel = " << simInfo.nChannels << std::endl;
		std::cout << "Porosity = " << simInfo.porosity << std::endl;
	}

	if (simInfo.nChannels != 1){
		printf("Error: please enter a grascale image with 1 channel.\n Current number of channels = %d\n", simInfo.nChannels);
		return 1;
	}

	// Sort out the current mesh

	if(opts.MeshIncreaseX < 1 || opts.MeshIncreaseY < 1){						// Return error if mesh refinement is smaller than 1
		printf("MeshIncrease has to be an integer greater than 1.\n");
		return 1;
	}

	// Define number of cells in each direction

	simInfo.numCellsX = simInfo.Width*opts.MeshIncreaseX;
	simInfo.numCellsY = simInfo.Height*opts.MeshIncreaseY;
	simInfo.nElements = simInfo.numCellsX*simInfo.numCellsY;
	simInfo.dx = 1.0/simInfo.numCellsX;
	simInfo.dy = 1.0/simInfo.numCellsY;

	// Diffusion coefficients

	float ks = opts.TCsolid;
	float kf = opts.TCfluid;

	// We will use an artificial scaling of the diffusion coefficient to converge to the correct solution

	// Declare useful arrays
	float *K = (float*)malloc(sizeof(float)*simInfo.numCellsX*simInfo.numCellsY); 			// Grid matrix containing the diffusion coefficient of each cell with appropriate mesh
	float *QL = (float*)malloc(sizeof(float)*simInfo.numCellsY);										// mass flux in the left boundary
	float *QR = (float*)malloc(sizeof(float)*simInfo.numCellsY);										// mass flux in the right boundary

	float *CoeffMatrix = (float *)malloc(sizeof(float)*simInfo.nElements*5);					// array will be used to store our coefficient matrix
	float *RHS = (float *)malloc(sizeof(float)*simInfo.nElements);										// array used to store RHS of the system of equations
	float *TemperatureMap = (float *)malloc(sizeof(float)*simInfo.nElements);			// array used to store the solution to the system of equations
	float *temp_TMap = (float *)malloc(sizeof(float)*simInfo.nElements);			// array used to store the solution to the system of equations

	// Initialize the concentration map with a linear gradient between the two boundaries
	for(int i = 0; i<simInfo.numCellsY; i++){
		for(int j = 0; j<simInfo.numCellsX; j++){
			TemperatureMap[i*simInfo.numCellsX + j] = (float)j/simInfo.numCellsX*(opts.TempRight - opts.TempLeft) + opts.TempLeft;
		}
	}

	// Zero the time

	simInfo.gpuTime = 0;

	// Declare GPU arrays

	float *d_x_vec = NULL;
	float *d_temp_x_vec = NULL;
	
	float *d_Coeff = NULL;
	float *d_RHS = NULL;

	// Initialize the GPU arrays

	if(!initializeGPU(&d_x_vec, &d_temp_x_vec, &d_RHS, &d_Coeff, simInfo))
	{
		printf("\n Error when allocating space in GPU");
		unInitializeGPU(&d_x_vec, &d_temp_x_vec, &d_RHS, &d_Coeff);
		return 0;
	}

	// Populate K array

	for(int i = 0; i<simInfo.numCellsY; i++){
		QL[i] = 0;
		QR[i] = 0;
		for(int j = 0; j<simInfo.numCellsX; j++){
			int targetIndexRow = i/opts.MeshIncreaseY;
			int targetIndexCol = j/opts.MeshIncreaseX;
			if(simInfo.target_data[targetIndexRow*simInfo.Width + targetIndexCol] < 150){
				K[i*simInfo.numCellsX + j] = kf;
			} else{
				K[i*simInfo.numCellsX + j] = ks;
			}
		}
	}

	// Initialize arrays for discretization

	memset(RHS, 0, sizeof(RHS));
	memset(CoeffMatrix, 0, sizeof(CoeffMatrix));

	// Discretize

	DiscretizeMatrix2D(K, CoeffMatrix, RHS, simInfo, opts);

	// Solve

	int iter_taken = 0;
	iter_taken = JacobiGPU(CoeffMatrix, RHS, TemperatureMap, temp_TMap, opts, 
		d_x_vec, d_temp_x_vec, d_Coeff, d_RHS, QL, QR, K, &simInfo);

	// create output file 
	printf("Final Keff = %2.3f\n", simInfo.keff);
	outputSingle(opts, simInfo);

	// create tmap

	if(opts.printTmap == 1){
		printTMAP(&opts, TemperatureMap, simInfo.numCellsY, simInfo.numCellsY, 0);
	}

	if(opts.printQmap == 1){
		printQMAP(&opts, TemperatureMap, K, simInfo.numCellsY, simInfo.numCellsX, simInfo.dx, simInfo.dy, QR, QL, 0);
	}
	// Free everything

	unInitializeGPU(&d_x_vec, &d_temp_x_vec, &d_RHS, &d_Coeff);
	free(QL);
	free(QR);
	free(CoeffMatrix);
	free(RHS);
	free(TemperatureMap);
	free(temp_TMap);
	free(K);

	return 0;
}