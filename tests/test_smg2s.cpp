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
    MPI_Init(&argc, &argv);
    int probSize = 7;

    auto nilp = Nilpotent<int>(2, 2, probSize);

 /*   auto mat = nonherm<std::complex<double>, int>(probSize, nilp, initMat<int>(-4, -2), "v2.txt", MPI_COMM_WORLD);
    mat.MatView();
    mat.show();

    mat.writeToMatrixMarketCmplx("mat_v2_nonherm.mtx");
*/
    auto mat2 = nonsymm<double, int>(probSize, nilp, initMat<int>(-4, -2), "v1.txt", MPI_COMM_WORLD);
    mat2.MatView();

 
/*    auto mat3 = nonsymm<double, int>(probSize, nilp, initMat<int>(-4, -2), "v2.txt", MPI_COMM_WORLD);
    mat3.MatView();

    mat2.writeToMatrixMarket("mat_v1_nonsym.mtx");
*/
   // mat3.writeToMatrixMarket("mat_v2_nonsym.mtx");

	MPI_Finalize();
    

	return 0;
}
