# include "keff3D.h"

int main(int argc, char const *argv[])
{
	// Declare initialization variables

	int numThreads;
	double startTime, runTime;

	startTime = omp_get_wtime();

	fflush(stdout);

	if(argc == 1){
		printf("Number of threads set to default = 1.\n");
		numThreads = 1;
		omp_set_num_threads(numThreads);
	} else{
		numThreads = strtol(argv[1], NULL, 10);
		printf("Number of Threads = %d\n", numThreads);
		omp_set_num_threads(numThreads);
	}

	double k_Fluid = 1;
	double k_Solid = 10.0;

	long int iterLimit = 100000;
	double thresh = 0.00001;

	double TR = 1;
	double TL = 0;

	int Width = 128;
	int Height = 128;
	int Depth = 128;

	char filename[100];

	sprintf(filename, "Schwarz163.csv");

	// Declare strructure, set all to fluid

	unsigned char *myStructure = (unsigned char*)malloc(sizeof(unsigned char)*Height*Width*Depth);

	for(int i = 0; i<Depth; i++){
		for(int j = 0; j<Height; j++){
			for(int k = 0; k<Width; k++){
				myStructure[i*Height*Width + j*Width + k] = 0;
			}
		}
	}

	// Read input file

	readImage3D(myStructure, Width, Height, Depth, filename);

	// Mesh Options and size

	int AmpFactorX = 1;
	int AmpFactorY = 1;
	int AmpFactorZ = 1;

	if (AmpFactorX < 1 || AmpFactorY < 1 || AmpFactorZ < 1){
		printf("Feature not currently available, AmpFactor has to be integer larger than 1.\n");
		return 1;
	}		

	// calculate number of cells

	int numCellsX = Width*AmpFactorX;
	int numCellsY = Height*AmpFactorY;
	int numCellsZ = Depth*AmpFactorZ;

	int nElements = numCellsX*numCellsY*numCellsZ;

	// get centers

	double *xCenters = (double *)malloc(sizeof(double)*numCellsX);
	double *yCenters = (double *)malloc(sizeof(double)*numCellsY);
	double *zCenters = (double *)malloc(sizeof(double)*numCellsZ);

	for(int i = 0; i<numCellsX; i++){
		xCenters[i] = (double)i*(1.0/numCellsX) + 0.5*(1.0/numCellsX);
	}

	for(int i = 0; i<numCellsY; i++){
		yCenters[i] = (double)i*(1.0/numCellsY) + 0.5*(1.0/numCellsY);
	}

	for(int i = 0; i<numCellsZ; i++){
		zCenters[i] = (double)i*(1.0/numCellsZ) + 0.5*(1.0/numCellsZ);
	}

	// Define more useful arrays

	double *k_Matrix = (double *)malloc(sizeof(double)*numCellsX*numCellsY*numCellsZ);
	double *QL = (double *)malloc(sizeof(double)*numCellsY*numCellsZ);
	double *QR = (double *)malloc(sizeof(double)*numCellsY*numCellsZ);

	// Populate arrays

	for (int i = 0; i<numCellsZ; i++){
		for(int j = 0; j<numCellsY; j++){
			QL[i*Height + j] = 0;
			QR[i*Height + j] = 0;
			for(int k = 0; k<numCellsZ; k++){
				int targetIndCol = k/AmpFactorX;
				int targetIndRow = j/AmpFactorY;
				int targetIndSlice = i/AmpFactorZ;
				if(myStructure[targetIndSlice*Width*Height + targetIndRow*Width + targetIndCol] == 1){
					k_Matrix[i*numCellsY*numCellsX + j*numCellsX + k] = k_Fluid;
				} else{
					k_Matrix[i*numCellsY*numCellsX + j*numCellsX + k] = k_Solid;
				}
			}
		}
	}

	// Define arrays for discretized matrix and solution

	double *CoeffMatrix = (double *)malloc(sizeof(double)*numCellsY*numCellsX*numCellsZ*7);
	double *TemperatureDist = (double *)malloc(sizeof(double)*numCellsY*numCellsX*numCellsZ);
	double *RHS = (double *)malloc(sizeof(double)*numCellsX*numCellsY*numCellsZ);

	// Discretize

	DiscretizeMatrixCD3D(k_Matrix, numCellsY, numCellsX, numCellsZ, CoeffMatrix, RHS, TL, TR, xCenters, yCenters, zCenters);

	printf("Done discretizing, solving now.\n");

	SolInitLinear(TemperatureDist, TL, TR, numCellsX, numCellsY, numCellsZ);

	ParallelGS3D(CoeffMatrix, RHS, TemperatureDist, QL, QR, k_Matrix, iterLimit, thresh, numCellsX, numCellsY,
	numCellsZ, TL, TR, xCenters, yCenters, zCenters, numThreads);
	

	printf("Done solving.\n");




	return 0;
}