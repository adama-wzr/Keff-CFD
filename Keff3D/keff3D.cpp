# include "keff3D.h"

int main(int argc, char const *argv[])
{
	// Declare initialization variables

	int numThreads;
	double startTime, runTime;

	start_time = omp_get_wtime();

	fflush(stdout);

	if(argc == 1){
		printf("Number of threads set to default = 1.\n");
		num_threads = 1;
		omp_set_num_threads(num_threads);
	} else{
		num_threads = strtol(argv[1], NULL, 10);
		printf("Number of Threads = %d\n", num_threads);
		omp_set_num_threads(num_threads);
	}

	double k_fluid = 1;
	double k_solid = 10.0;

	double TR = 1;
	double TL = 0;

	int img_size = 128;

	

	return 0;
}