# include "keff2D.h"

#define IMG_SIZE 128
#define ITER_LIMIT 500000
#define CONVERGE  0.000001

int main(int argc, char const *argv[])
{
	// initialization variables for openmp

	int num_threads;

	// start timer

	double start_time, run_time;

	start_time = omp_get_wtime();

	fflush(stdout);

	// user input number of threads and default

	if(argc == 1){
		printf("Number of threads set to default = 1.\n");
		num_threads = 1;
		omp_set_num_threads(num_threads);
	} else{
		num_threads = strtol(argv[1], NULL, 10);
		printf("Number of Threads = %d\n", num_threads);
		omp_set_num_threads(num_threads);
	}

	// Initialize domain variables

	double kFluid = 1;
	double kSolid = 10.0;

	double TL = 0;
	double TR = 1;

	std::cout << "kFluid = " << kFluid << std::endl;

	// Read image

	unsigned char *target_data;

	int width, height, channel;

	char filename[30];

	sprintf(filename, "00001.jpg");

	// Call function to read the image:

	readImage(target_data, &width, &height, &channel, filename);

	// Print if successful

	if (channel != 1){
		printf("Error: please enter a grascale image with 1 channel.\n Current image Number of Channels = %d", channel);
		return 1;
	}

	std::cout << "width = " << width << " height = " << height << " channel = " << channel << std::endl;

	// Define mesh related variables:

	bool structuredMesh = true;

	// Define amplification factors for mesh resolution
	int AmpFactorX = 1;
	int AmpFactorY = 1;

	if (AmpFactorX < 1 || AmpFactorY < 1){
		printf("Feature not currently available, variable has to be integer larger than 1.\n");
		return 1;
	}

	// calculate number of cells
	int numCellsX = width*AmpFactorX;
	int numCellsY = height*AmpFactorY;

	// Now get the center of each cell:
	double *xCenters = (double *)malloc(sizeof(double)*numCellsX);
	double *yCenters = (double *)malloc(sizeof(double)*numCellsY);


	#pragma omp parallel for schedule(auto)
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
				k[i*numCellsX + j] = kFluid; 			// black => fluid => 0 => void
			} else{
				k[i*numCellsX + j] = kSolid;			// white => solid => 1 => material
			}
		}
	}
	
	// Get number of elements for coefficient matrix, create more useful arrays for solution

	double *CoeffMatrix = (double *)malloc(sizeof(double)*numCellsY*numCellsX*5); 		// Only store 5 diagonals
	double *TemperatureDist = (double *)malloc(sizeof(double)*numCellsY*numCellsX);
	double *RHS = (double *)malloc(sizeof(double)*numCellsY*numCellsX);
	

	return 0;
}