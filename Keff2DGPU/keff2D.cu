#include "keff2D.cuh"

int main(void){

	// Important call for efficiency on Linux
	fflush(stdout);

	//	Declare data structure
	options opts;

	char inputFilename[30];

	sprintf(inputFilename, "input.txt");

	readInputFile(inputFilename, &opts);

	SingleSim(opts);

	// if(opts.BatchFlag == 0){
	// 	SingleSim(opts);
	// } else if(opts.BatchFlag == 1){
	// 	BatchSim(opts);
	// } else if(opts.MeshFlag == 1){
    //     MeshStudy(opts);
    // } 
    // else{
	// 	std::cout << "Error: no valid BatchFlag option, check input file." << std::endl;
	// }

	return 0;

}