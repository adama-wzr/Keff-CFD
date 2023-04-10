# include "keff2D.h"

int main(int argc, char const *argv[])
{
	// initialization variables for openmp

	int num_threads;

	// start timer

	double start_time, run_time;

	start_time = omp_get_wtime();

	fflush(stdout);

	//
	options opts;
	// user input number of threads and default

	char inputFilename[30];

	sprintf(inputFilename, "input.txt");

	readInputFile(inputFilename, &opts);

	num_threads = opts.numCores;

	if(opts.verbose == 1){
		printf("Number of threads = %d\n", num_threads);
	}

	// Initialize domain variables

	double kFluid = opts.TCfluid;
	double kSolid = opts.TCsolid;

	double TL = opts.TempLeft;
	double TR = opts.TempRight;

	// Read image

	unsigned char *target_data;

	int width, height, channel;

	// Call function to read the image:

	readImage(&target_data, &width, &height, &channel, opts.inputFilename);

	if(opts.verbose == 1){
		std::cout << "width = " << width << " height = " << height << " channel = " << channel << std::endl;
	}

	if (channel != 1){
		printf("Error: please enter a grascale image with 1 channel.\n Current number of channels = %d\n", channel);
		return 1;
	}

	// Define mesh related variables:

	bool structuredMesh = true;

	// Define amplification factors for mesh resolution
	int AmpFactorX = opts.MeshIncreaseX;
	int AmpFactorY = opts.MeshIncreaseY;

	if (AmpFactorX < 1 || AmpFactorY < 1){
		printf("Feature not currently available, AmpFactor has to be integer larger than 1.\n");
		return 1;
	}

	// calculate number of cells
	int numCellsX = width*AmpFactorX;
	int numCellsY = height*AmpFactorY;

	// Now get the center of each cell:
	double *xCenters = (double *)malloc(sizeof(double)*numCellsX);
	double *yCenters = (double *)malloc(sizeof(double)*numCellsY);

	// #pragma omp parallel for schedule(auto)
	for(int i = 0; i<numCellsX; i++){
		xCenters[i] = (double)i*(1.0/numCellsX) + 0.5*(1.0/numCellsX);
	}

	#pragma omp parallel for schedule(auto)
	for(int i = 0; i<numCellsY; i++){
		yCenters[i] = (double)i*(1.0/numCellsY) + 0.5*(1.0/numCellsY);
	}

	// Define variables that are useful for the solution

	double *kMatrix = (double *)malloc(sizeof(double)*numCellsX*numCellsY);
	double *QL = (double *)malloc(sizeof(double)*numCellsY);
	double *QR = (double *)malloc(sizeof(double)*numCellsY);

	// Populate the arrays:

	#pragma omp parallel for schedule(auto)
	for(int i = 0; i<numCellsY; i++){
		QL[i] = 0;
		QR[i] = 0;
		for (int j = 0; j<numCellsX; j++){
			int targetIndex_Row = i/AmpFactorY;
			int targetIndex_Col = j/AmpFactorX;
			if(target_data[targetIndex_Row*width + targetIndex_Col] < 150){
				kMatrix[i*numCellsX + j] = kFluid; 			// black => fluid => 0 => void
			} else{
				kMatrix[i*numCellsX + j] = kSolid;			// white => solid => 1 => material
			}
		}
	}
	
	// Get number of elements for coefficient matrix, create more useful arrays for solution

	double *CoeffMatrix = (double *)malloc(sizeof(double)*numCellsY*numCellsX*5); 		// Only store 5 diagonals
	double *TemperatureDist = (double *)malloc(sizeof(double)*numCellsY*numCellsX);
	double *RHS = (double *)malloc(sizeof(double)*numCellsY*numCellsX);

	// Discretize

	DiscretizeMatrixCD2D(kMatrix, numCellsY, numCellsX, CoeffMatrix, RHS, TL, TR, xCenters, yCenters);

	// Initialize x-vector

	for (int i = 0; i< numCellsY; i++){
		for(int j = 0; j<numCellsX; j++){
			double a = (double)1.0/numCellsX;
			TemperatureDist[i*numCellsX + j] = (double)a*j;
		}
	}

	// Solve

	int iterTaken = 0;

	iterTaken = TDMAscanner(CoeffMatrix, RHS, TemperatureDist, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY, QR, QL,
	TR, TL, kMatrix, xCenters, yCenters);

	//
	double dy = yCenters[0]*2;
	for (int j = 0; j<numCellsY; j++){
		QL[j] = kMatrix[j*numCellsX]*dy*(TemperatureDist[j*numCellsX]-TL)/(xCenters[0]);
		QR[j] = kMatrix[(j+1)*numCellsX - 1]*dy*(TR - TemperatureDist[(j+1)*numCellsX-1])/(1-xCenters[numCellsX-1]);
	}

	double Q1 = 0;
	double Q2 = 0;

	for(int i =0; i<numCellsY; i++){
		Q1 += QL[i];
		Q2 += QR[i];
	}

	Q1 = Q1;
	Q2 = Q2;

	double Q_avg = (Q2 + Q1)/2;

	double k_eff = Q_avg/(TR-TL);

	printf("Keff = %f, iterTaken = %d\n", k_eff, iterTaken);
	run_time = omp_get_wtime() - start_time;

    printf("Time spent = %lf\n", run_time);


	return 0;
}