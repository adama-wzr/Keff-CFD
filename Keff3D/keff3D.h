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


typedef struct{
  double TCsolid;
  double TCfluid;
  int Width;
  int Height;
  int Depth;
  int MeshIncreaseX;
  int MeshIncreaseY;
  int MeshIncreaseZ;
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
		printf("Domain Width = %d, Height = %d, Depth = %d\n",opts->Width, opts->Height, opts->Depth);
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

	 	}else if(strcmp(tempC, "Width:") == 0){
	 		opts->Width = (int)tempD;

	 	}else if(strcmp(tempC, "Height:") == 0){
	 		opts->Height = (int)tempD;

	 	}else if(strcmp(tempC, "Depth:") == 0){
	 		opts->Depth = (int)tempD;

	 	}else if(strcmp(tempC, "MeshAmpX:") == 0){
	 		opts->MeshIncreaseX = (int)tempD;

	 	}else if(strcmp(tempC, "MeshAmpY:") == 0){
	 		opts->MeshIncreaseY = (int)tempD;

	 	}else if(strcmp(tempC, "MeshAmpZ:") == 0){
	 		opts->MeshIncreaseZ = (int)tempD;

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


int createOutput(options* o, double keff, double Q1, double Q2, long int iter, int nElements, double porosity, double time){
	/*
		createOutput Function:
		Inputs:
			- *o -> pointer to structure containing the options array.
			- keff -> double with calculated keff
			- Q1 -> heat flux through the left boundary
			- Q2 -> heat flux through the right boundary
			- iter -> number of iterations
			- nElements -> total number of elements
			- porosity -> calculated pixel based porosity
			- time -> time in seconds

		Outputs:
			-None.

		Function creates an output file, stored in the user entered address. The line of code below contains the correct order in which they are stored:

		keff,QL,QR,Iter,ConvergeCriteria,inputName,nElements,MeshIncreaseX,MeshIncreaseY,porosity,ks,kf,tl,tr,time,nCores

	*/

	FILE *OUTPUT;

  OUTPUT = fopen(o->outputFilename, "a+");
  fprintf(OUTPUT, "%.10lf,%f,%f,%ld,%.10lf,%s,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%d\n",keff, Q1, Q2, iter, o->ConvergeCriteria, o->inputFilename,
  	nElements, o->MeshIncreaseX, o->MeshIncreaseY, o->MeshIncreaseZ, porosity, o->TCsolid, o->TCfluid, o->TempLeft, o->TempRight, time, o->numCores);

  fclose(OUTPUT);
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


int createOutputBatch(options* o, double* Stats){
	/*
		createOutputBatch Function:
		Inputs:
			- *o -> pointer to structure containing the options array.
			- *Stats -> pointer to array containing statistics about each individual run using parallel computing.

		Outputs:
			-None.

		Function creates an output file, stored in the user entered address. The line of code below contains the correct order in which they are stored:

		keff,QL,QR,Iter,ConvergeCriteria,inputName,nElements,MeshIncreaseX,MeshIncreaseY,porosity,ks,kf,tl,tr,time,nCores

	*/

	FILE *OUTPUT;

  OUTPUT = fopen(o->outputFilename, "a+");
  char filename[100];
  for(int i = 0; i< o->NumImg; i++){
  	sprintf(filename,"%05d.csv",i);
  	fprintf(OUTPUT, "%.10lf,%f,%f,%ld,%.10lf,%s,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%d\n",Stats[i*7 + 0], Stats[i*7 + 1], Stats[i*7 + 2], (long int)Stats[i*7 + 3], o->ConvergeCriteria, filename,
  		(int)Stats[i*7 + 4], o->MeshIncreaseX, o->MeshIncreaseY, o->MeshIncreaseZ, Stats[i*7 + 6], o->TCsolid, o->TCfluid, o->TempLeft, o->TempRight, Stats[i*7 + 5], 1);
  }
  

  fclose(OUTPUT);
  return 0;
}


double calcPorosity3D(unsigned char* imageAddress, int Width, int Height, int Depth, options* opts){
	/*
		calcPorosity3D
		Inputs:
			- imageAddress: pointer to the read image.
			- Width: Number of cols
			- Height: number of rows
			- Depth: number of slices
			- opts: pointer to options struc

		Output:
			- porosity: double containing porosity.

		Function calculates porosity by counting pixels.
	*/

	double totalCells = (double)Height*Width*Depth;
	double porosity = 0;

	for(int i = 0; i<Depth; i++){
		for(int j = 0; j<Height; j++){
			for(int k = 0; k<Width; k++){
				if(imageAddress[i*Height*Width + j*Width + k] == opts->TCfluid){
					porosity += 1.0/totalCells;
				}
			}
		}
	}

	return porosity;
}


int printQMAP(options* o, double* x, double* K, int numRows, int numCols, int numSlices, double* xCenter, double* yCenter, double* zCenter, double* qR, double* qL, int imgNum){

	/*
		printQMAP:
		Inputs:
			- o -> options datastructure
			- x -> pointer to temperature distribution map.
			- numRows -> number of rows
			- numCols -> number of columns
			- numSlices -> number of slices
			- xCenter = pointer to array with center of each cell (x-coordinate)
			- yCenter = pointer to array with center of each cell (y-coordinate)
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

	double localQ;
	double localK;
	double dy = yCenter[2] - yCenter[1];
	double dx = xCenter[2] - xCenter[1];
	double dz = zCenter[2] - zCenter[1];

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


int printTMAP(options* o, double* x, int numRows, int numCols, int numSlices, int imgNum){
	/*
		printTMAP:
		Inputs:
			- o -> options datastructure
			- x -> pointer to temperature distribution map.
			- numRows -> number of rows
			- numCols -> number of columns
			- numSlices -> number of slices
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
  fprintf(T_OUT,"x,y,z,T\n");

  for(int i = 0; i<numSlices; i++){
		for(int j = 0; j<numRows; j++){
			for(int k = 0; k<numCols; k++){
				fprintf(T_OUT,"%d,%d,%d,%f\n",k,j,i,x[i*numCols*numRows + j*numCols + k]);
			}
		}
	}

  fclose(T_OUT);
  return 0;
}


int ParallelGS3D(double *arr, double *sol, double *x_vec, double *qL, double *qR, double *K, long int iterLimit, double tolerance, int numCols, int numRows,
	int numSlices, double tL, double tR, double *XC, double *YC, double *ZC, int nCores, options* opts){
	/*
	Function ParallelGS3D:
	Inputs:
		- *arr: pointer to the discretization matrix, size (numCols*numRows, 5)
		- *sol: pointert to the RHS, size (numCols*numRows)
		- *x_vec: pointer to the solution vector, size (numCols*numRows)
		- *qL: pointer a vector of size (numRows) containing heat flux from each cell
			at the left side of the domain.
		- *qR: pointer a vector of size (numRows) containing heat flux from each cell
			at the right side of the domain.
		- *K: vector containing the thermal conductivity of each cell in the domain.
		- iterLimit: integer with the maximum number of allowed iterations
		- tolerance: value of convergence criteria
		- numCols: numCols of the domain
		- numRows: numRows of the domain
		- numSlices: numSlices of the domain
		- tL: temperature in the left boundary.
		- tR: temperature in the right boundary.
		- *XC: pointer to vector containing the x-coordinate of the center of each cell (length numCols)
		- *YC: pointer to vector containing the y-coordinate of the center of each cell (length numRows)
		- *ZC: pointer to vector containing the z-coordinate of the center of each cell (length numSlices)
		- nCores: integer, number of cores.

	Outputs:
		- int IterationCount

	Gauss-Seidel iteration implemented in parallel using OpenMP. This function is specifically
	meant for the 7 diagonal discretization. 

	*/
	int iterCount = 0;
	double sigma = 0;
	double norm_diff = 1;
	double percentChange = 1;
	int i;
	double keffOld = 1;
	double keffNew = 1;
	int iterToCheck = 1000;
	double Q1,Q2;
	double qAvg = 0;
	#pragma omp parallel private(i, sigma)
	

	while(percentChange > tolerance && iterCount < iterLimit){
		#pragma omp parallel private(i, sigma)
		#pragma omp for schedule(dynamic, numRows*numCols*numSlices/nCores)
		for(i = 0; i<numRows*numCols*numSlices; i++){
			sigma = 0;
			for(int j = 1; j<7; j++){
				if(arr[i*7 + j] != 0){
					if(j == 1){
						sigma += arr[i*7 + j]*x_vec[i - 1];
					} else if(j == 2){
						sigma += arr[i*7 + j]*x_vec[i + 1];
					} else if(j == 3){
						sigma += arr[i*7 + j]*x_vec[i + numCols];
					} else if(j == 4){
						sigma += arr[i*7 + j]*x_vec[i - numCols];
					} else if(j==5){
						sigma += arr[i*7 + j]*x_vec[i - numCols*numRows];
					}else if(j==6){
						sigma += arr[i*7 + j]*x_vec[i + numCols*numRows];
					}
				}
			}
			x_vec[i] = 1/arr[i*7 + 0]*(sol[i] - sigma);
		}

		iterCount++;

		if (iterCount %iterToCheck == 0){
			Q1 = 0;
			Q2 = 0;
			double dx = XC[0];
			double dy = YC[0]*2;
			double dz = ZC[0]*2;
			for(int j = 0; j<numSlices; j++){
				for(int k = 0; k<numRows; k++){
					qL[j*numRows + k] = K[j*numRows*numCols + k*numCols]*dy*dz*(x_vec[j*numRows*numCols + k*numCols] - tL)/dx;
					qR[j*numRows + k] = K[j*numRows*numCols + (k+1)*numCols - 1]*dy*dz*(tR - x_vec[j*numRows*numCols + (k+1)*numCols - 1])/dx;

					Q1+=qL[j*numRows + k];
					Q2+=qR[j*numRows + k];
				}
			}

			Q1 = Q1;
			Q2 = Q2;
			qAvg = (Q1 + Q2)/2;
			keffNew = qAvg/(tR - tL);
			percentChange = fabs((keffNew - keffOld)/keffOld);
			keffOld = keffNew;

			if (percentChange < 0.001){
				iterToCheck = 100;
			} else if(percentChange < 0.0001){
				iterToCheck = 10;
			}
			if(opts->verbose == 1){
				printf("Iteration = %d, Keff = %f, Convergence = %f\n", iterCount, keffNew, percentChange);
			}
		}
	}
	return iterCount;
}


int SolInitLinear(double *xVec, double tL, double tR, int numCols, int numRows, int numSlices){
	/*
	Function SolInitLinear:
	Inputs:
		- *xVec: pointer to the temperature map allocated vector
		- double tL: temperature at the left boundary
		- double tR: temperature at the right boundary
		- int numCols: number of cells in the x-direction
		- int numRows: number of cells in the y-direction
		- int numSlices: number of cells in the z-direction

	Outputs: none

	SolInitLinear initializes the temperature assuming a linear distribution between
	the right and left boundaries.
	*/
	for(int i = 0; i<numSlices; i++){
		for(int j = 0; j<numRows; j++){
			for(int k = 0; k<numCols; k++){
				if(tR > tL){
					xVec[i*numCols*numSlices + j*numCols + k] = 1/numCols*k*(tR - tL) + tL;
				} else if(tR < tL){
					xVec[i*numCols*numSlices + j*numCols + k] = 1/numCols*(numCols - (k+1))*(tL - tR) + tR;
				}
			}
		}
	}

	return 0;
}


int SolExpandSimple(double *OGTemperature, double *ExpandedT, int NativeWidth, int NativeHeight, int NativeDepth, int AmpFactor){
	/*
	Function SolExpandSimple:
	Inputs:
		- *OGTemperature: pointer to the solution temperature map of the original mesh size
		- *ExpandedT: pointer to array that will receive the new temperature distribution
				and use it as initial guess at the solution.
		- int NativeWidth: width of the original solution
		- int NativeHeight: height of the original solution
		- int NativeDepth: depth of original solution
		- int AmpFactor: increase in size from the original

	Outputs: none

	SolExpandSimple grabs the final solution temperature distribution and simply expands
	it onto a larger mesh. No calculations or interpolation is done, this is just to roughly
	initialize it close to a solution.
	*/

	for(int i = 0; i<NativeDepth*AmpFactor; i++){
		for(int j = 0; j< NativeHeight*AmpFactor; j++){
			for(int k = 0; k<NativeWidth*AmpFactor; k++){
				int targetIndexSlice = i/AmpFactor;
				int targetIndexRow = j/AmpFactor;
				int targetIndexCol = k/AmpFactor;

				ExpandedT[i*(NativeWidth*NativeHeight*AmpFactor*AmpFactor) + j*(NativeWidth*AmpFactor) + k] = 
					OGTemperature[targetIndexSlice*NativeWidth*NativeHeight + targetIndexRow*NativeWidth + targetIndexCol];
			}
		}
	}
	
	return 0;
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
				b[index] = 0;
				for(int temp = 0; temp<7; temp++){
					A[index*7 + temp] = 0;
				}

				// i == 0 is the slice with the least depth, only P and B
				if (i == 0){
					kb = WeightedHarmonicMean(dz/2,dz/2, K[index], K[index+numRows*numCols]);
					A[index*7 + 6] = -kb*dy*dx/dz;
					A[index*7 + 0] += (kb*dy*dx/dz);
				} else if(i == numSlices - 1){
					// only P and F
					kf = WeightedHarmonicMean(dz/2,dz/2, K[index], K[index-numRows*numCols]);
					A[index*7 + 5] = -kf*dy*dx/dz;
					A[index*7 + 0] += (kf*dy*dx/dz);
				} else{
					kb = WeightedHarmonicMean(dz/2,dz/2, K[index], K[index+numRows*numCols]);
					kf = WeightedHarmonicMean(dz/2,dz/2, K[index], K[index-numRows*numCols]);

					A[index*7 + 5] = -kf*dy*dx/dz;
					A[index*7 + 6] = -kb*dy*dx/dz;
					A[index*7 + 0] += (kf*dy*dx/dz + kb*dy*dx/dz);
				}

				// j == 0 is the top boundary, P and S
				if(j == 0){
					ks = WeightedHarmonicMean(dy/2, dy/2, K[index], K[index + numCols]);
					A[index*7 + 3] = -ks*dx*dz/dy;
					A[index*7 + 0] += (ks*dx*dz/dy);
				} else if(j == numRows - 1){
				// bottom boundary, P and N
					kn = WeightedHarmonicMean(dy/2, dy/2, K[index], K[index - numCols]);
					A[index*7 + 4] = -kn*dx*dz/dy;
					A[index*7 + 0] += (kn*dx*dz/dy);
				} else{
					ks = WeightedHarmonicMean(dy/2, dy/2, K[index], K[index + numCols]);
					kn = WeightedHarmonicMean(dy/2, dy/2, K[index], K[index - numCols]);

					A[index*7 + 3] = -ks*dx*dz/dy;
					A[index*7 + 4] = -kn*dx*dz/dy;
					A[index*7 + 0] += (kn*dx*dz/dy + ks*dx*dz/dy);
				}

				// k == 0 is the left boundary, P and E
				if(k == 0){
					ke = WeightedHarmonicMean(dx/2, dx/2, K[index], K[index + 1]);
					kw = K[index];

					A[index*7 + 2] = -ke*dy*dz/dx;
					A[index*7 + 0] += (ke*dy*dz/dx + kw*dy*dz/(dx/2));
					b[index] += tL*kw*dy*dz/(dx/2);
				} else if(k == numCols - 1){
					// only P and W
					ke = K[index];
					kw = WeightedHarmonicMean(dx/2, dx/2, K[index], K[index - 1]);

					A[index*7 + 1] = -kw*dy*dz/dx;
					A[index*7 + 0] += (ke*dy*dz/(dx/2) + kw*dy*dz/dx);
					b[index] += tR*ke*dy*dz/(dx/2);
				} else{
					ke = WeightedHarmonicMean(dx/2, dx/2, K[index], K[index + 1]);
					kw = WeightedHarmonicMean(dx/2, dx/2, K[index], K[index - 1]);
					A[index*7 + 1] = -kw*dy*dz/dx;
					A[index*7 + 2] = -ke*dy*dz/dx;
					A[index*7 + 0] += (ke*dy*dz/dx + kw*dy*dz/dx);
				}
			}
		}
	}

	return 0;
}