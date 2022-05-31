#include <C/c-smg2s.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

    int world_size;
    int world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    nilp_t *nilp = newNilp_1(2, 8);
    nilp_show(nilp);


    int probSize = 7;
    int span, lower_b, upper_b;

    span = (int)ceil((double)probSize/(double)world_size);


    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    parVecMap_t *p = newParVecMap(MPI_COMM_WORLD, lower_b, upper_b);

	ds_parVec_t *vec = new_ds_ParVec_2(p);

	MPI_Finalize();
}