#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include <string>
#include <typeinfo>  
#include <mpi.h>

#include "../parMatrix/MatrixCSR.h"
#include "../parMatrix/parMatrixSparse.h"

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
    parVector<double, int> vec = parVector<double, int>(parVecMap);

    auto Matrix = parMatrixSparse<double, int>(vec);
   	//initMat(S diag_l, S diag_u, Base<T> scale, Base<T> shift, Base<T> sparsity );

    auto spec2 = specNonSymm<double, int>(parVecMap, "v1.txt");
    //spec2.VecView();

    //auto spec3 = specNonSymmCplex<double, int>("v3.txt");
    //spec3.VecView();
    //Matrix.setSpecNonSymm(spec3);
    //Matrix.show();
    //Matrix.MatView();

    
    Matrix.initMat(-4, -2, 2.0, 1.0, 0.5);
    //Matrix.show();
    //Matrix.MatView();

    auto nilp = Nilpotent<int>(2, 2, probSize);
	auto iz = nilp.getIndOfZeros();

	if(world_rank == 0){
		for(auto i = 0; i < iz.size(); i++){
			std::cout << "getIndOfZeros ]>" << i << ": " << iz[i] << std::endl;
		}
    }


/*
    auto prod = Matrix.MA(nilp);
	//prod.show();
    //prod.MatView();
    auto prod2 = Matrix.AM(nilp);
    prod2.MatView();
*/
	MPI_Finalize();

	return 0;
}