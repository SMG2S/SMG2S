//#include "parVectorMap.cc"
//#include "../parMatrix/parMatrixSparse.cc"
#include "specGen.h"
#include <math.h>
#include <complex.h>
#include "../utils/logo.h"

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

    if(world_rank == 0) {logo(1.0);}

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


    std::complex<double> a = 1.0, c = 2.0; float b =2.0;
    int lbandwidth = 3;
  
    parVector<std::complex<double>,int> *vec = new parVector<std::complex<double>,int>(MPI_COMM_WORLD, lower_b, upper_b);
    parVector<std::complex<double>,int> *prod = new parVector<std::complex<double>,int>(MPI_COMM_WORLD, lower_b, upper_b);


    MPI_Barrier(MPI_COMM_WORLD);

    //generate vec containing the given spectra
    specGen<std::complex<double>,int>(vec);

    //vec->VecView();

    //Matrix Initialization

    parMatrixSparse<std::complex<double>,int> *Am = new parMatrixSparse<std::complex<double>,int>(vec,prod);
    parMatrixSparse<std::complex<double>,int> *MA = new parMatrixSparse<std::complex<double>,int>(vec,prod);
    parMatrixSparse<std::complex<double>,int> *MB = new parMatrixSparse<std::complex<double>,int>(vec,prod);
    parMatrixSparse<std::complex<double>,int> *MC = new parMatrixSparse<std::complex<double>,int>(vec,prod);
    parMatrixSparse<std::complex<double>,int> *AM = new parMatrixSparse<std::complex<double>,int>(vec,prod);

    if(world_rank == 0){printf("Matrix Initialized\n");}

    MPI_Barrier(MPI_COMM_WORLD);

    //setup the lower part of initial matrix

    start = MPI_Wtime();

    for(int i = 0; i < probSize; i++){
        for(int j = i - lbandwidth; j < i; j++){
            if(j >= 0){
                Am->Loc_SetValue(i,j,1);   
            }
        }
    }

    //insert the diagonal of initial matrix with given spectra.

    Am->Loc_SetDiagonal(vec);

    end = MPI_Wtime();

    double t2 = end - start;

    if(world_rank == 0) {printf("Initial matrix generation time = %1.6f\n", t2);}

    MPI_Barrier(MPI_COMM_WORLD);

    start = MPI_Wtime();

//    Am->glocPlusLloc();

    MPI_Barrier(MPI_COMM_WORLD);

    //Am->LOC_MatView();

    end = MPI_Wtime();

    t2 = end - start;
    
    if(world_rank == 0) {printf("Combination time = %1.6f\n", t2);}

    Nilpotency<int> nilp, nps;
    
    nilp.NilpType1(2,probSize);

    start = MPI_Wtime();

    Am->MA(nilp, MA);

    end = MPI_Wtime();

    t2 = end - start;

    if(world_rank == 0) {printf("MA time = %1.6f\n", t2);}

    MPI_Barrier(MPI_COMM_WORLD);

    start = MPI_Wtime();

    Am->AM(nilp, MA);

    end = MPI_Wtime();

    t2 = end - start;

    if(world_rank == 0) {printf("AM time = %1.6f\n", t2);}


//    MA->LOC_MatView();

/*
    MPI_Barrier(MPI_COMM_WORLD);

    MA->llocToGlocLoc();

    MPI_Barrier(MPI_COMM_WORLD);

    //MA->MatView();

    Am->Loc_ZeroEntries();

    MPI_Barrier(MPI_COMM_WORLD);

    MA->MA(nilp, MB);

    MPI_Barrier(MPI_COMM_WORLD);

    //MB->LOC_MatView();

    MPI_Barrier(MPI_COMM_WORLD);

    MB->MA(nilp, MC);

    MPI_Barrier(MPI_COMM_WORLD);

    //MC->LOC_MatView();
//    Am->MatView();
*/

    MPI_Finalize();

    return 0;
}

