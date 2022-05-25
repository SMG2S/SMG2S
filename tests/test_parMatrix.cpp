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
    int probSize = 11;

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
    vec.SetToValue(3.3);
    vec.VecView();

    auto Matrix = parMatrixSparse<double, int>(vec);

    Matrix.initMat(-5, -3, 1.0, 0.0, 0.0);
    //Matrix.MatView();

    Matrix.SetValue(3,0,2.4);
    //Matrix.MatView();

    Matrix.SetDiagonal(vec);
    Matrix.MatScale(2.0);

    auto X = parMatrixSparse<double, int>(vec);
    auto Y = parMatrixSparse<double, int>(vec);
    X.initMat(-6,-4, 1.0, 0.0, 0.0);
    Y.initMat(-7,-4, 1.0, 0.0, 0.0);
    //X.MatView();
    //Y.MatView();
    //Matrix.MatView();

    Matrix.MatAXPY(X, 2);
    //Matrix.MatView();

    Matrix.MatAYPX(Y, 0.0001);
    //Matrix.MatView();

    //Matrix.writeToMatrixMarket("matrix.txt");

    Nilpotent<int> nilp = Nilpotent<int>(2,2, probSize);
    nilp.show();
    auto iz = nilp.getIndOfZeros();

    if(world_rank == 0){
        for(auto i = 0; i < iz.size(); i++){
            std::cout << iz[i] << std::endl;
        }
    }

    //auto MA = Matrix.MA(nilp);
    auto AM = Matrix.AM(nilp);
    //Matrix.copy(Y);
    Matrix.MatView();
    //MA.MatView();
    AM.MatView();

	MPI_Finalize();

	return 0;
}