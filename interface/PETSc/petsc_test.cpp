/*
   This file is part of SMG2S.
   Author(s): Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr or xinzhe.wu1990@gmail.com>
        Date: 2018-04-20
   Copyright (C) 2018-     Xinzhe WU
   
   SMG2S is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   SMG2S is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with SMG2S.  If not, see <http://www.gnu.org/licenses/>.
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

