/*

MIT License

Copyright (c) 2019 Xinzhe WU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

//#include "parVectorMap.cc"
//#include "../parMatrix/parMatrixSparse.cc"
#include <math.h>
#include "petsc_interface.h"

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc,&argv,(char *)0,NULL);
  //  PetscInitialize(&argc,&args,(char*)0,help);

    // Get the number of processes
    int world_size;

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


    int probSize = 5;
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
  
    parVector<std::complex<double>,int> *vec = new parVector<std::complex<double>,int>(MPI_COMM_WORLD, lower_b, upper_b);
    parVector<std::complex<double>,int> *prod = new parVector<std::complex<double>,int>(MPI_COMM_WORLD, lower_b, upper_b);


    MPI_Barrier(MPI_COMM_WORLD);

    parMatrixSparse<std::complex<double>,int> *Am = new parMatrixSparse<std::complex<double>,int>(vec,prod);


    if(world_rank == 0){printf("Matrix Initialized\n");}

    MPI_Barrier(MPI_COMM_WORLD);

    //setup the lower part of initial matrix


    for(int i = 0; i < probSize; i++){
        for(int j = i - 3; j < i; j++){
            if(j >= 0){
                Am->Loc_SetValue(i,j,0.5);   
            }
        }
    }



    for(int i = 0; i < probSize; i++){
        Am->Loc_SetValue(i,i,1.0*i+1);
    }  


    MPI_Barrier(MPI_COMM_WORLD);

    Am->LOC_MatView();

    Am->Loc_ConvertToCSR();

    Mat M;
    MatCreate(PETSC_COMM_WORLD,&M);

    M = ConvertToPETSCMat(Am);

    MatView(M, PETSC_VIEWER_STDOUT_WORLD);
    PetscFree(M);
    PetscFinalize();

    MPI_Finalize();

    return 0;
}

