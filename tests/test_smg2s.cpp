#include <mpi.h>
#include <smg2s-interface.hpp>

int main(int argc, char** argv) 
{
    MPI_Init(&argc, &argv);
    
    int world_size;
    int world_rank;
    int probSize = 7;
    int l_diag = -7;
    int u_diag = -3;
    int nbOne = 2;
    int offset = 1;
    double sparsity = 0.5;

    /* construct a nilpotent object for generation */
    Nilpotent<int> nilp = Nilpotent<int>(nbOne, offset, probSize);
    
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int span, lower_b, upper_b;
    span = int(floor(double(probSize)/double(world_size)));

    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    /* construct a parVecMap object which determines the distribution scheme of vectors and matrices*/
    auto parVecMap = parVectorMap<int>(MPI_COMM_WORLD, lower_b, upper_b);
  	
  	/* example 1, generation of a non-Hermtian matrix */
  	// 1. generate the spectrum on the fly
    parVector<std::complex<double>, int> spec1 = parVector<std::complex<double>, int>(parVecMap);
    for(int i = lower_b; i < upper_b; i++){
    	std::complex<double> v(i+1, i+2);
        spec1.SetValueGlobal(i, v);
    }
 	// 2. generation 
    auto mat = nonherm<std::complex<double>, int>(probSize, nilp, initMat<int>(l_diag, u_diag, sparsity), spec1);

  	/* example 2, generation of a non-Symmetric matrix with real eigenvalues */
  	// 1. generate the spectrum on the fly
	parVector<double, int> spec2 = parVector<double, int>(parVecMap);
    for(int i = lower_b; i < upper_b; i++){
        spec2.SetValueGlobal(i, i+1);
    }
	// 2. generation 
    auto mat2 = nonsymm<double , int>(probSize, nilp, initMat<int>(l_diag, u_diag, sparsity), spec2);

  	/* example 3, generation of a non-Symmetric matrix with conjugated eigenvalues */
  	// 1. generate the spectrum on the fly
    parVector<std::complex<double>, int> spec3 = parVector<std::complex<double>, int>(parVecMap);

    for(int i = lower_b; i < upper_b; i++){
        if(i % 2 == 0){
            std::complex<double> v(i/2 + 1, i/2 + 2);
            spec3.SetValueGlobal(i, v);
        }else{
            std::complex<double> v(i/2 + 1, -i/2 - 2);
            spec3.SetValueGlobal(i, v);
        }
        if(i == probSize - 1){
            std::complex<double> v(i + 1, 0);
            spec3.SetValueGlobal(i, v);
        }
    }
	// 2. generation 
    auto mat3 = nonsymmconj<double , int>(probSize, nilp, initMat<int>(l_diag, u_diag, sparsity), spec3);

	MPI_Finalize();
}
