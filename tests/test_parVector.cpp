#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include <string>
#include <typeinfo>  
#include <mpi.h>

#include <smg2s-interface.hpp>


int main(int argc, char** argv) 
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    int world_size;
    int world_rank;
    int probSize = 7;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int span, lower_b, upper_b;

    span = int(ceil(double(probSize)/double(world_size)));

    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    auto parVecMap = parVectorMap<int>(MPI_COMM_WORLD, lower_b, upper_b);

    auto spec3 = specNonSymmCplex<double, int>(parVecMap, "v2.txt");
    spec3.VecView();

    //for(auto i = 0 ; i < spec3.size(); i++){
    //	std::cout << spec3[i] << std::endl;
    //}

	MPI_Finalize();

	return 0;
}