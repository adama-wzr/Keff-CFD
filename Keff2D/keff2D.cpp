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

	// Initialize domain variables

	omp_set_num_threads(opts.numCores);

	double kFluid = opts.TCfluid;
	double kSolid = opts.TCsolid;

	double TL = opts.TempLeft;
	double TR = opts.TempRight;

	double porosity = 0;


	// Read image

	int width, height, channel;
	unsigned char *target_data;

	if(opts.BatchFlag == 0){

		// Call function to read the image:

		readImage(&target_data, &width, &height, &channel, opts.inputFilename);

		porosity = calcPorosity(target_data, width, height);

		if(opts.verbose == 1){
			std::cout << "width = " << width << " height = " << height << " channel = " << channel << std::endl;
			std::cout << "Porosity = " << porosity << std::endl;
		}

		if (channel != 1){
			printf("Error: please enter a grascale image with 1 channel.\n Current number of channels = %d\n", channel);
			return 1;
		}
	}

	// Now we split the code into its different modes:


	if(opts.MeshFlag == 0 && opts.BatchFlag == 0){
		// run normal code

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
		int nElements = numCellsX*numCellsY;

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

		SolInitLinear(TemperatureDist, TR, TL, numCellsX, numCellsY);

		// Solve

		int iterTaken = 0;

		// iterTaken = ParallelJacobi(CoeffMatrix, RHS, TemperatureDist, QL, QR, kMatrix, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY,
		// 	TL, TR, xCenters, yCenters);

		// iterTaken = TDMAscanner(CoeffMatrix, RHS, TemperatureDist, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY, QR, QL,
		// 	TR, TL, kMatrix, xCenters, yCenters);

		iterTaken = ParallelGS(CoeffMatrix, RHS, TemperatureDist, QL, QR, kMatrix, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY,
			TL, TR, xCenters, yCenters, num_threads);

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

		printf("Keff = %f, Total Iterations = %d\n", k_eff, iterTaken);

		if(opts.printTmap == 1){
			printTMAP(&opts, TemperatureDist, numCellsY, numCellsX, 0);
		}

		if(opts.printQmap == 1){
			printQMAP(&opts, TemperatureDist, kMatrix, numCellsY, numCellsX, xCenters, yCenters, QL, QR, 0);
		}

		

		run_time = omp_get_wtime() - start_time;

		createOutput(&opts, k_eff, Q1, Q2, iterTaken, nElements, porosity, run_time);

	    printf("Time spent = %lf\n", run_time);
	} else if(opts.MeshFlag == 1 && opts.BatchFlag == 0){
		// run mesh convergence test

		// Allocate an array to store the temperature dist with original mesh.
		// This is used to estimate the initial temp. dist at
		// other mesh sizes.

		double *NativeMeshTemperature = (double *)malloc(sizeof(double)*width*height);

		for(int m = 1; m <= opts.MaxMesh; m++){
			start_time = omp_get_wtime();
			int AmpFactorX = m;
			int AmpFactorY = m;
		

			// calculate number of cells
			int numCellsX = width*AmpFactorX;
			int numCellsY = height*AmpFactorY;
			int nElements = numCellsX*numCellsY;

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

			if(m == 1){
				SolInitLinear(TemperatureDist, TR, TL, numCellsX, numCellsY);
			} else{
				SolExpandSimple(NativeMeshTemperature, TemperatureDist, width, height, m);
			}

			

			// Solve

			int iterTaken = 0;

			// iterTaken = ParallelJacobi(CoeffMatrix, RHS, TemperatureDist, QL, QR, kMatrix, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY,
			// 	TL, TR, xCenters, yCenters);

			// iterTaken = TDMAscanner(CoeffMatrix, RHS, TemperatureDist, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY, QR, QL,
			// 	TR, TL, kMatrix, xCenters, yCenters);

			iterTaken = ParallelGS(CoeffMatrix, RHS, TemperatureDist, QL, QR, kMatrix, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY,
				TL, TR, xCenters, yCenters, num_threads);

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

			printf("-----------------------------------------\n");

			printf("Keff = %f, Total Iterations = %d\n", k_eff, iterTaken);
			

			run_time = omp_get_wtime() - start_time;

			createOutput(&opts, k_eff, Q1, Q2, iterTaken, nElements, porosity, run_time);

		    printf("Time (s) = %lf, Number of Elements = %d\n\n", run_time, nElements);

		    if(m == 1){
		    	for(int i = 0; i<numCellsY; i++){
		    		for(int j = 0; j<numCellsX; j++){
		    			NativeMeshTemperature[i*numCellsX + j] = TemperatureDist[i*numCellsX + j];
		    		}
		    	}
		    }

		    free(TemperatureDist);
		    free(kMatrix);
		    free(QL);
		    free(QR);
		    free(CoeffMatrix);
		    free(RHS);
		    free(xCenters);
		    free(yCenters);
		}
	} else if(opts.MeshFlag == 0 && opts.BatchFlag == 1){

		// keff,q1,q2,iter,nElements,time,porosity
		double* Results = (double *)malloc(sizeof(double)*7*opts.NumImg);

		int AmpFactorX = opts.MeshIncreaseX;
		int AmpFactorY = opts.MeshIncreaseY;

		#pragma omp parallel for schedule(auto)
		for(int img = 0; img < opts.NumImg; img++){

			printf("Image %05d.jpg\n",img);

			double privateStartTime;

			privateStartTime = omp_get_wtime();

			unsigned char *target_data;

			int width, height, channel;

			char filename[100];

			sprintf(filename, "%05d.jpg", img);

			// Call function to read the image:

			readImage(&target_data, &width, &height, &channel, filename);

			double porosity = calcPorosity(target_data, width, height);

			if (channel != 1){
				printf("Error: please enter a grascale image with 1 channel.\n Current number of channels = %d\n", channel);
				printf("Consider stopping the code manually, the RGB image will likely result in error.\n");
			}

			int numCellsX = width*AmpFactorX;
			int numCellsY = height*AmpFactorY;
			int nElements = numCellsX*numCellsY;

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

			SolInitLinear(TemperatureDist, TR, TL, numCellsX, numCellsY);

			// Solve

			int iterTaken = 0;

			// iterTaken = ParallelJacobi(CoeffMatrix, RHS, TemperatureDist, QL, QR, kMatrix, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY,
			// 	TL, TR, xCenters, yCenters);

			// iterTaken = TDMAscanner(CoeffMatrix, RHS, TemperatureDist, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY, QR, QL,
			// 	TR, TL, kMatrix, xCenters, yCenters);

			iterTaken = SerialGS(CoeffMatrix, RHS, TemperatureDist, QL, QR, kMatrix, opts.MAX_ITER, opts.ConvergeCriteria, numCellsX, numCellsY,
				TL, TR, xCenters, yCenters);

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

			if(opts.printTmap == 1){
				printTMAP(&opts, TemperatureDist, numCellsY, numCellsX, img);
			}

			if(opts.printQmap == 1){
				printQMAP(&opts, TemperatureDist, kMatrix, numCellsY, numCellsX, xCenters, yCenters, QL, QR, img);
			}

			Results[img*7 + 0] = Q_avg/(TR-TL);
			Results[img*7 + 1] = Q1;
			Results[img*7 + 2] = Q2;
			Results[img*7 + 3] = iterTaken;
			Results[img*7 + 4] = numCellsX*numCellsY;
			Results[img*7 + 5] = omp_get_wtime() - privateStartTime;
			Results[img*7 + 6] = porosity;

			free(TemperatureDist);
		    free(kMatrix);
		    free(QL);
		    free(QR);
		    free(CoeffMatrix);
		    free(RHS);
		    free(xCenters);
		    free(yCenters);
		}
		printf("Done, saving output\n");
		createOutputBatch(&opts, Results);
		free(Results);
	}



		


	return 0;
}