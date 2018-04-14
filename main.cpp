//#include "parVectorMap.cc"
//#include "../parMatrix/parMatrixSparse.cc"
#include "src/matgen/smg2s.h"
#include <math.h>
#include <complex.h>
#include "../utils/logo.h"
#include "confg/config.h"

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int size;

    double start, end;

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) {logo(1.0);}

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    //MPI_Barrier(MPI_COMM_WORLD);
    
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, rank, size);

    
    int probSize, lbandwidth;

    Nilpotency<int> nilp;
    
    nilp.NilpType1(2,probSize);

    smg2s(probSize, nilp, lbandwidth);

    MPI_Finalize();

    return 0;
}

