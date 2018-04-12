//#include "parVectorMap.cc"
//#include "../parMatrix/parMatrixSparse.cc"
#include "specGen.h"
#include <math.h>
#include <complex.h>

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;

    double start, end;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    //MPI_Barrier(MPI_COMM_WORLD);
    
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);


    int probSize = 20;
    int span, lower_b, upper_b;

    span = int(floor(double(probSize)/double(world_size)));

    printf("span = %d\n", span);
    if(world_rank == world_size - 1){
	lower_b = world_rank * span;
	upper_b = probSize - 1 + 1;
    }
    else{
	lower_b = world_rank * span;
	upper_b = (world_rank + 1) * span - 1 + 1;
    }

    printf("Proc. %d   Lower bound = %d   Upper bound = %d \n",world_rank, lower_b , upper_b ) ;

    double a = 1.0, c = 2.0; float b =2.0;
    int lbandwidth = 3;
  
    parVector<double,int> *vec = new parVector<double,int>(MPI_COMM_WORLD, lower_b, upper_b);
    parVector<double,int> *prod = new parVector<double,int>(MPI_COMM_WORLD, lower_b, upper_b);


    MPI_Barrier(MPI_COMM_WORLD);

    //generate vec containing the given spectra
    specGen<double,int>(vec);

    //vec->VecView();

    //Matrix Initialization

    parMatrixSparse<double,int> *Am = new parMatrixSparse<double,int>(vec,prod);
    parMatrixSparse<double,int> *Bm = new parMatrixSparse<double,int>(vec,prod);

    if(world_rank == 0){printf("Matrix Initialized\n");}

    MPI_Barrier(MPI_COMM_WORLD);

    //setup the lower part of initial matrix

    start = MPI_Wtime();

    for(int i = 0; i < probSize; i++){
        for(int j = i - lbandwidth; j < i; j++){
            if(j >= 0){
                Am->SetValue(i,j,1);   
            }
        }
    }

    Am->SetDiagonal(vec);

    end = MPI_Wtime();

    double t2 = end - start;

    if(world_rank == 0) {printf("Initial matrix generation time = %1.6f\n", t2);}

    MPI_Barrier(MPI_COMM_WORLD);

    start = MPI_Wtime();

    Am->glocPlusLloc();

    MPI_Barrier(MPI_COMM_WORLD);

    Am->LOC_MatView();

    end = MPI_Wtime();

    t2 = end - start;
    
    if(world_rank == 0) {printf("Combination time = %1.6f\n", t2);}

    //insert the diagonal of initial matrix with given spectra.


    MPI_Barrier(MPI_COMM_WORLD);

//    Am->MatView();

    MPI_Finalize();

    return 0;
}

