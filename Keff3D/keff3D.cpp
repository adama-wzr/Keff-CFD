# include "keff3D.h"

int main(int argc, char const *argv[])
{
	// Declare initialization variables

	int numThreads;
	double startTime, runTime;

	options opts;
	char inputName[100];
	sprintf(inputName,"input.txt");
	readInputFile(inputName, &opts);

	startTime = omp_get_wtime();

	fflush(stdout);

	numThreads = opts.numCores;

	double k_Fluid = opts.TCfluid;
	double k_Solid = opts.TCsolid;

	long int iterLimit = opts.MAX_ITER;
	double thresh = opts.ConvergeCriteria;

	double TR = opts.TempRight;
	double TL = opts.TempLeft;

	int Width = opts.Width;
	int Height = opts.Height;
	int Depth = opts.Depth;

	char filename[100];

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

	readImage3D(myStructure, Width, Height, Depth, opts.inputFilename);

	double porosity = calcPorosity3D(myStructure, Width, Height, Depth, &opts);

	if(opts.verbose == 1){
		std::cout << "Widht: " << Width << " Height: " << Height << " Depth: " << Depth << std::endl;
		std::cout << "Porosity" << porosity << std::endl;
	}

	if(opts.MeshFlag == 0 && opts.BatchFlag == 0){
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
			xCenters[i] = (double)i*(1.0) + 0.5*(1.0/numCellsX);
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
				for(int k = 0; k<numCellsX; k++){
					int targetIndCol = k/AmpFactorX;
					int targetIndRow = j/AmpFactorY;
					int targetIndSlice = i/AmpFactorZ;
					if(myStructure[targetIndSlice*Width*Height + targetIndRow*Width + targetIndCol] == 1){
						k_Matrix[i*numCellsY*numCellsX + j*numCellsX + k] = k_Solid;
					} else{
						k_Matrix[i*numCellsY*numCellsX + j*numCellsX + k] = k_Fluid;
					}
				}
			}
		}

		// Define arrays for discretized matrix and solution

		double *CoeffMatrix = (double *)malloc(sizeof(double)*numCellsY*numCellsX*numCellsZ*7);
		double *TemperatureDist = (double *)malloc(sizeof(double)*numCellsY*numCellsX*numCellsZ);
		double *RHS = (double *)malloc(sizeof(double)*numCellsX*numCellsY*numCellsZ);

		// Discretize

		if(opts.verbose == 1){
			std::cout << "Discretizing" << std::endl;
		}

		DiscretizeMatrixCD3D(k_Matrix, numCellsY, numCellsX, numCellsZ, CoeffMatrix, RHS, TL, TR, xCenters, yCenters, zCenters);

		SolInitLinear(TemperatureDist, TL, TR, numCellsX, numCellsY, numCellsZ);

		if(opts.verbose == 1){
			std::cout << "Solving" << std::endl;
		}

		int iterTaken = ParallelGS3D(CoeffMatrix, RHS, TemperatureDist, QL, QR, k_Matrix, iterLimit, thresh, numCellsX, numCellsY,
		numCellsZ, TL, TR, xCenters, yCenters, zCenters, numThreads, &opts);

		double Q1 = 0;
		double Q2 = 0;
		double dx = xCenters[0];
		double dy = yCenters[0]*2;
		double dz = zCenters[0]*2;
		for(int j = 0; j<numCellsZ; j++){
			for(int k = 0; k<numCellsY; k++){
				QL[j*numCellsY + k] = k_Matrix[j*numCellsY*numCellsX + k*numCellsX]*dy*dz*(TemperatureDist[j*numCellsY*numCellsX + k*numCellsX] - TL)/dx;
				QR[j*numCellsY + k] = k_Matrix[j*numCellsY*numCellsX + (k+1)*numCellsX - 1]*dy*dz*(TR - TemperatureDist[j*numCellsY*numCellsX + (k+1)*numCellsX - 1])/dx;

				Q1+=QL[j*numCellsY + k];
				Q2+=QR[j*numCellsY + k];
			}
		}

		Q1 = Q1;
		Q2 = Q2;
		double qAvg = (Q1 + Q2)/2;
		double k_eff = qAvg/(TR-TL);

		if(opts.verbose == 1){
			std::cout << "Done." << std::endl;
		}

		printf("--------------------------------------\n\n");

		if(opts.verbose == 1){
			std::cout << "Effective Thermal Conductivity = " << k_eff << std::endl;
		}

		if(opts.printTmap == 1){
			printTMAP(&opts, TemperatureDist, numCellsY, numCellsX, numCellsZ, 0);
		}

		if(opts.printQmap == 1){
			printf("Heat Flux Map: Feature not currently available.\n");
		}

		runTime = omp_get_wtime() - startTime;

		createOutput(&opts, k_eff, Q1, Q2, iterTaken, nElements, porosity, runTime);

	} else if(opts.MeshFlag == 1 && opts.BatchFlag == 0){

		double *NativeMeshTemperature = (double *)malloc(sizeof(double)*Width*Depth*Height);

		for(int m = 1; m<=opts.MaxMesh; m++){
			// Mesh convergence mode
			int AmpFactorX = m;
			int AmpFactorY = m;
			int AmpFactorZ = m;

			int numCellsX = Width*AmpFactorX;
			int numCellsY = Height*AmpFactorY;
			int numCellsZ = Depth*AmpFactorZ;

			int nElements = numCellsX*numCellsY*numCellsZ;

			// get centers

			double *xCenters = (double *)malloc(sizeof(double)*numCellsX);
			double *yCenters = (double *)malloc(sizeof(double)*numCellsY);
			double *zCenters = (double *)malloc(sizeof(double)*numCellsZ);

			for(int i = 0; i<numCellsX; i++){
				xCenters[i] = (double)i*(1.0) + 0.5*(1.0/numCellsX);
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
					for(int k = 0; k<numCellsX; k++){
						int targetIndCol = k/AmpFactorX;
						int targetIndRow = j/AmpFactorY;
						int targetIndSlice = i/AmpFactorZ;
						if(myStructure[targetIndSlice*Width*Height + targetIndRow*Width + targetIndCol] == 1){
							k_Matrix[i*numCellsY*numCellsX + j*numCellsX + k] = k_Solid;
						} else{
							k_Matrix[i*numCellsY*numCellsX + j*numCellsX + k] = k_Fluid;
						}
					}
				}
			}

			// Define arrays for discretized matrix and solution

			double *CoeffMatrix = (double *)malloc(sizeof(double)*numCellsY*numCellsX*numCellsZ*7);
			double *TemperatureDist = (double *)malloc(sizeof(double)*numCellsY*numCellsX*numCellsZ);
			double *RHS = (double *)malloc(sizeof(double)*numCellsX*numCellsY*numCellsZ);

			// Discretize

			if(opts.verbose == 1){
				std::cout << "Discretizing" << std::endl;
			}

			DiscretizeMatrixCD3D(k_Matrix, numCellsY, numCellsX, numCellsZ, CoeffMatrix, RHS, TL, TR, xCenters, yCenters, zCenters);

			if(m == 1){
				SolInitLinear(TemperatureDist, TL, TR, numCellsX, numCellsY, numCellsZ);
			} else{
				SolExpandSimple(NativeMeshTemperature, TemperatureDist, Width, Height, Depth, m);
			}

			SolInitLinear(TemperatureDist, TL, TR, numCellsX, numCellsY, numCellsZ);

			if(opts.verbose == 1){
				std::cout << "Solving" << std::endl;
			}

			int iterTaken = ParallelGS3D(CoeffMatrix, RHS, TemperatureDist, QL, QR, k_Matrix, iterLimit, thresh, numCellsX, numCellsY,
			numCellsZ, TL, TR, xCenters, yCenters, zCenters, numThreads, &opts);

			double Q1 = 0;
			double Q2 = 0;
			double dx = xCenters[0];
			double dy = yCenters[0]*2;
			double dz = zCenters[0]*2;
			for(int j = 0; j<numCellsZ; j++){
				for(int k = 0; k<numCellsY; k++){
					QL[j*numCellsY + k] = k_Matrix[j*numCellsY*numCellsX + k*numCellsX]*dy*dz*(TemperatureDist[j*numCellsY*numCellsX + k*numCellsX] - TL)/dx;
					QR[j*numCellsY + k] = k_Matrix[j*numCellsY*numCellsX + (k+1)*numCellsX - 1]*dy*dz*(TR - TemperatureDist[j*numCellsY*numCellsX + (k+1)*numCellsX - 1])/dx;

					Q1+=QL[j*numCellsY + k];
					Q2+=QR[j*numCellsY + k];
				}
			}

			Q1 = Q1;
			Q2 = Q2;
			double qAvg = (Q1 + Q2)/2;
			double k_eff = qAvg/(TR-TL);

			if(opts.verbose == 1){
				std::cout << "Done." << std::endl;
			}

			printf("--------------------------------------\n\n");

			if(opts.verbose == 1){
				std::cout << "Effective Thermal Conductivity = " << k_eff << std::endl;
				std::cout << "Time (s) = " << runTime << ", Number of Elements = " << nElements << std::endl;
			}

			runTime = omp_get_wtime() - startTime;

			createOutput(&opts, k_eff, Q1, Q2, iterTaken, nElements, porosity, runTime);

			if(m == 1){
				for(int i = 0; i<numCellsZ; i++){
					for(int j = 0; j<numCellsY; j++){
						for(int k = 0; k<numCellsX; k++){
							NativeMeshTemperature[i*numCellsX*numCellsY + j*numCellsX + k] = TemperatureDist[i*numCellsX*numCellsY + j*numCellsX + k];
						}
					}
				}
			}


			free(TemperatureDist);
		    free(k_Matrix);
		    free(QL);
		    free(QR);
		    free(CoeffMatrix);
		    free(RHS);
		    free(xCenters);
		    free(yCenters);
		    free(zCenters);

		}
	} else if(opts.MeshFlag == 0 && opts.BatchFlag == 1){

		// keff,q1,q2,iter,nElements,time,porosity
		double* Results = (double *)malloc(sizeof(double)*7*opts.NumImg);

		int AmpFactorX = 1;
		int AmpFactorY = 1;
		int AmpFactorZ = 1;

		#pragma omp parallel for schedule(auto)
		for(int img = 0; img <  opts.NumImg; img++){
			if(opts.verbose == 1){
				printf("Image number %05d.csv", img);
			}

			double privateStartTime;

			privateStartTime = omp_get_wtime();

			int Width = opts.Width;
			int Height = opts.Height;
			int Depth = opts.Depth;

			char filename[100];

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

			readImage3D(myStructure, Width, Height, Depth, opts.inputFilename);

			double porosity = calcPorosity3D(myStructure, Width, Height, Depth, &opts);

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
				xCenters[i] = (double)i*(1.0) + 0.5*(1.0/numCellsX);
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
					for(int k = 0; k<numCellsX; k++){
						int targetIndCol = k/AmpFactorX;
						int targetIndRow = j/AmpFactorY;
						int targetIndSlice = i/AmpFactorZ;
						if(myStructure[targetIndSlice*Width*Height + targetIndRow*Width + targetIndCol] == 1){
							k_Matrix[i*numCellsY*numCellsX + j*numCellsX + k] = k_Solid;
						} else{
							k_Matrix[i*numCellsY*numCellsX + j*numCellsX + k] = k_Fluid;
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

			// initialize solution

			SolInitLinear(TemperatureDist, TL, TR, numCellsX, numCellsY, numCellsZ);

			// solve

			int iterTaken = ParallelGS3D(CoeffMatrix, RHS, TemperatureDist, QL, QR, k_Matrix, iterLimit, thresh, numCellsX, numCellsY,
			numCellsZ, TL, TR, xCenters, yCenters, zCenters, numThreads, &opts);

			double Q1 = 0;
			double Q2 = 0;
			double dx = xCenters[0];
			double dy = yCenters[0]*2;
			double dz = zCenters[0]*2;
			for(int j = 0; j<numCellsZ; j++){
				for(int k = 0; k<numCellsY; k++){
					QL[j*numCellsY + k] = k_Matrix[j*numCellsY*numCellsX + k*numCellsX]*dy*dz*(TemperatureDist[j*numCellsY*numCellsX + k*numCellsX] - TL)/dx;
					QR[j*numCellsY + k] = k_Matrix[j*numCellsY*numCellsX + (k+1)*numCellsX - 1]*dy*dz*(TR - TemperatureDist[j*numCellsY*numCellsX + (k+1)*numCellsX - 1])/dx;

					Q1+=QL[j*numCellsY + k];
					Q2+=QR[j*numCellsY + k];
				}
			}

			Q1 = Q1;
			Q2 = Q2;
			double qAvg = (Q1 + Q2)/2;
			double k_eff = qAvg/(TR-TL);

			if(opts.printTmap == 1){
				printTMAP(&opts, TemperatureDist, numCellsY, numCellsX, numCellsZ, 0);
			}

			if(opts.printQmap == 1){
				// printf("Heat Flux Map: Feature not currently available.\n");
			}

			Results[img*7 + 0] = qAvg/(TR-TL);
			Results[img*7 + 1] = Q1;
			Results[img*7 + 2] = Q2;
			Results[img*7 + 3] = iterTaken;
			Results[img*7 + 4] = numCellsX*numCellsY;
			Results[img*7 + 5] = omp_get_wtime() - privateStartTime;
			Results[img*7 + 6] = porosity;

			free(TemperatureDist);
		    free(k_Matrix);
		    free(QL);
		    free(QR);
		    free(CoeffMatrix);
		    free(RHS);
		    free(xCenters);
		    free(yCenters);
		    free(zCenters);
		}

		printf("Done, saving output\n");
		createOutputBatch(&opts, Results);
		free(Results);
	} else{
		printf("Invalid Input options, exiting.\n");

		return 1;
	}

	return 0;
}