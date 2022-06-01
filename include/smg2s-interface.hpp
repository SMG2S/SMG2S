/*
MIT License
Copyright (c) 2019 Xinzhe WU @ Maison de la Simulation, France
Copyright (c) 2019-2022, Xinzhe Wu @ Simulation and Data Laboratory Quantum 
									 Materials,  Forschungszentrum Juelich GmbH.
									 
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

#ifndef __SMG2S_INTERFACE_HPP_
#define __SMG2S_INTERFACE_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>
#include <vector>
#include <mpi.h>
#include <map>
#include <exception>
#include <complex>

#include <parMatrix/MatrixCSR.hpp>
#include <parVector/parVectorMap.hpp>
#include <parVector/parVector.hpp>
#include <smg2s/spectrum.hpp>
#include <parMatrix/parMatrixSparse.hpp>
#include <smg2s/nilpotent.hpp>
#include <smg2s/initMat.hpp>
#include <smg2s/smg2s.hpp>



/** @mainpage Overview
*
* @authors Xinzhe Wu @ Simulation and Data Laboratory Quantum 
									 Materials,  Forschungszentrum Juelich GmbH.
*
* @section intro Introduction
* 
* SMG2S (Sparse Matrix Generator with Given Spectrum) is a package which provides to generator the non-Hermitian and non-Symmetric sparse Matrices with user-customized eigenvalues. The generated sparse matrices contain still 
* the user-provided spectrum. SMG2S is implemented in parallel based on MPI (Message Passing Interface) and C++14 to support efficiently the generation of test matrices on distributed memory platforms.
* 
* The initial motivation of the creation of SMG2S comes from the studay of iterative methods
* for solving linear systems and eigenvalue problems. Iterative linear algebra methods are important for the applications in various fields. The analysis of the iterative method behaviors is complex, and it is necessary to evaluate their convergence to solve extremely large non-Hermitian/non-symmtric eigenvalue and linear problems on parallel and/or distributed machines. This convergence depends on the properties of spectra. Thus, we propse SMG2S to generate large matrices with known spectra to benchmark these methods. The generated matrices are non-trivial with very high dimension.
* 
* As a matrix generator, SMG2S provides:
* 
* - generating of both Non-Hermitian and Non-Symmetric sparse matrix
* 
* - generated matrices are naturally Sparse with non-trivial sparsity pattern
* 
* - Given Spectrum: the spectrum of generated matrix is the same as the one specified by the users
* 
* - Sparsity Patterns are diverse and controllable
* 
* 
* As a software, SMG2S provides:
* 
* - a collection of c++ header files
* 
* - parallel implementation based on MPI which is able to efficiently generate very large sparse matrices in parallel on supercomputers
* 
* - a verification module based on Python for the plotting and verification of small size of generated matrix.
* 
* - an easy-to-use C interface 
*
* 
* 
* \image html ../docs/figure/verification_5.png width=1000px
* \image html ./docs/logo.png width=0px
* <hr>
* @section notes Release Notes
* This is the release note
* <hr>
* @section Requirements
* <hr> 
* Some requirements
*
*/

/** @page started User guide
* <hr>
* @section quick Quick start
\code{.sh}
  git clone https://github.com/SMG2S/SMG2S.git
  cd SMG2S
  mkdir build
  cd build
  cmake ..
  make
\endcode

\code{.sh}
  mpirun -np ${PROCS} ./smg2s.exe -D ${dim} -L ${diag_l} -U ${diag_u} -O ${offset} -C ${nbOne} -S ${sparsity} -M {no-herm or non-symm}
\endcode

\code{.sh}
usage: ./smg2s.exe [options] ...
options:
  -D, --dim           Dimension of matrix to be generated (int [=1000])
  -L, --diagL         offset of lower diagonal of initial matrix (int [=-10])
  -U, --diagU         offset of upper diagonal of initial matrix (int [=-5])
  -O, --nilpOffset    offset of diagonal of a nilpotent (int [=5])
  -C, --continous     Continuous length in Nilpotent matrix (int [=2])
  -S, --sparsity      sparsity of initial matrix (NOT THE FINAL GENERATED ONES) (double [=0.95])
  -M, --mattype       Matrix type to be generated: non-symmetric or non-Hermitian (string [=non-herm])
  -?, --help          print this message
\endcode
* <hr>
* @section install Installation
\code{.sh}
cmake .. -DCMAKE_INSTALL_PREFIX=${PATH_TO_INSTALL}
make -j install
\endcode
* <hr>
* @section use Usage
* @subsection cmake CMake
\code{.sh}
cmake_minimum_required(VERSION 3.8)
project(xxx)
#find installation of SMG2S
find_package( smg2s REQUIRED CONFIG)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
     message([FATAL_ERROR] "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()
# for C++ code
add_executable(smg2s-app test_parMatrix.cpp)
target_link_libraries(smg2s-app PUBLIC SMG2S::smg2s)
# for C code
add_executable(test_c.exe test_c.c)
target_link_libraries(test_c.exe PRIVATE SMG2S::smg2s2c)
\endcode

*@subsection link Direct linking

* @section example Mini example
\code{.cpp}
#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include <string>
#include <typeinfo>  
#include <mpi.h>
#include <random>
#include <cmath>

#include <smg2s-interface.hpp>


int main(int argc, char** argv) 
{
    MPI_Init(&argc, &argv);
    
    int world_size;
    int world_rank;
    int probSize = 100;
    int l_diag = -60;
    int u_diag = -20;
    int nbOne = 5;
    int offset = 50;
    double sparsity = 0.5;

    Nilpotent<int> nilp = Nilpotent<int>(nbOne, offset, probSize);

    nilp.show();

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

    auto parVecMap = parVectorMap<int>(MPI_COMM_WORLD, lower_b, upper_b);
    parVector<std::complex<double>, int> spec1 = parVector<std::complex<double>, int>(parVecMap);
    specGenNonHerm2(&spec1);
    spec1.writeToTxtCmplx("specNonHerm1.txt");

    auto mat = nonherm<std::complex<double>, int>(probSize, nilp, initMat<int>(-l_diag, -u_diag, 0.1, sparsity), spec1);
    mat.writeToMatrixMarketCmplx("matGenNonHerm1.mtx");

    parVector<double, int> spec2 = parVector<double, int>(parVecMap);
    specGenNonSymm(&spec2, 10.0, 0.0);
    spec2.writeToTxt("specNonSymm1.txt");
    auto mat2 = nonsymm<double , int>(probSize, nilp, initMat<int>(-l_diag, -u_diag, 1, sparsity), spec2);
    mat2.writeToMatrixMarket("matGenNonSymm1.mtx");

    parVector<std::complex<double>, int> spec3 = parVector<std::complex<double>, int>(parVecMap);
    specGenNonSymmConj(&spec3);
    spec3.writeToTxtCmplx("specNonSymmConj1.txt");    
//    spec3.VecView();
    auto mat3 = nonsymmconj<double , int>(probSize, nilp, initMat<int>(-l_diag, -u_diag, 1, sparsity), spec3);
    mat3.writeToMatrixMarket("matGenNonSymmConj1.mtx");

    parVector<std::complex<double>, int> spec4 = parVector<std::complex<double>, int>(parVecMap);
    specGenNonSymm2(&spec4, 100000.0, 0.0);
    spec4.writeToTxtCmplx("specNonSymm2.txt");    
//    spec3.VecView();
    auto mat4 = nonsymmconj<double , int>(probSize, nilp, initMat<int>(-l_diag, -u_diag, 1, sparsity), spec4);
    mat4.writeToMatrixMarket("matGenNonSymm2.mtx");


	MPI_Finalize();
    

	return 0;
}
\endcode

* <hr>
* @section custom Customization
* @subsection nilp Nilpotent matrix
* @subsection initmat Initialization of matrix
* @subsection spec Set spectrum
* @subsection asm Assembling the customizations
* <hr>
* @section verify Plotting and verification
* * \image html ../docs/figure/verification_1png width=1000px
* \image html ../docs/figure/verification_2.png width=1000px
* \image html ../docs/figure/verification_3.png width=1000px
* \image html ../docs/figure/verification_4.png width=1000px
* \image html ../docs/figure/verification_5.png width=1000px
* <hr>
* @section c-interface C interface
\code{.c}
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

    int probSize = 7;
    int span, lower_b, upper_b;

    int diagl = -5;
    int diagu = -3;
    double sparsity = 0.5;

    nilp_t *nilp = newNilp_1(2, 7);
    nilp_show(nilp);

    span = (int)ceil((double)probSize/(double)world_size);


    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    initMatrix_t *initMat = newInitMatrix_3(diagl, diagu, sparsity);  

    parVecMap_t *p = newParVecMap(MPI_COMM_WORLD, lower_b, upper_b);

    //non-symm case 1
	ds_parVec_t *spec = new_ds_ParVec_2(p);

    for(int i = lower_b; i < upper_b; i++){
        ds_parVecSetVal(spec, i, i + 1);
    }

    ds_parVecView(spec);

    //initMatrix_show (initMat);
    ds_parMatSparse_t *mat = ds_nonsymm_2(probSize, nilp, initMat, spec);    

    ds_parMatSparse_View(mat);  
    ds_parMatSparse_destory(mat);

    //non-symm case 2 with conjugate eigenvalues
    zs_parVec_t *spec2 = new_zs_ParVec_2(p);
    for(int i = lower_b; i < upper_b; i++){
        if(i % 2 == 0){
            dcomplex_t v = {i/2 + 1, i/2 + 2};
            zs_parVecSetVal(spec2, i, v);
        }else{
            dcomplex_t v = {i/2 + 1, -i/2 - 2};
            zs_parVecSetVal(spec2, i, v);
        }
        if(i == probSize - 1){
            dcomplex_t v = {i + 1, 0};
            zs_parVecSetVal(spec2, i, v);
        }
    }

    zs_parVecView(spec2);
    ds_parMatSparse_t *mat2 = ds_nonsymmconj_2(probSize, nilp, initMat, spec2);    
    ds_parMatSparse_View(mat2);  
    ds_parMatSparse_destory(mat2);

    //non-herm case 
    zs_parVec_t *spec3 = new_zs_ParVec_2(p);
    for(int i = lower_b; i < upper_b; i++){
        dcomplex_t v = {i+1, i+2};
        zs_parVecSetVal(spec3, i, v);
    }

    zs_parVecView(spec3);

    zs_parMatSparse_t *mat3 = zs_nonherm_2(probSize, nilp, initMat, spec3);    
    zs_parMatSparse_View(mat3);  
    zs_parMatSparse_destory(mat3);
    initMatrix_destory(initMat);

    ds_parVec_destory(spec);
    zs_parVec_destory(spec2);
    zs_parVec_destory(spec3);

    nilp_destory (nilp);

	MPI_Finalize();
}
\endcode
* <hr>
* @section moreexample More examples
* <hr>
*/

/** @page Gallery Gallery: Sparsity Patterns
* <hr>
n=100, diag_l=20, diag_u=5, nbOne=10, offset=20, sparsity=0.9
 * \image{inline} html ../docs/figure/100_20_5_10_20_0.9.png
* <hr>
n=100, diag_l=60, diag_u=2, nbOne=5, offset=5, sparsity=0.9
* \image{inline} html ../docs/figure/100_60_2_5_5_0.9.png
* <hr>
n=100, diag_l=60, diag_u=20, nbOne=5, offset=20, sparsity=0.9
* \image{inline} html ../docs/figure/100_60_20_5_20_0.9.png
* <hr>
n=100, diag_l=60, diag_u=20, nbOne=5, offset=50, sparsity=0.9
* \image{inline} html ../docs/figure/100_60_20_5_50_0.9.png
* <hr>
n=100, diag_l=80, diag_u=5, nbOne=10, offset=20, sparsity=0.9
* \image{inline} html ../docs/figure/100_80_5_10_20_0.9.png
* <hr>
n=100, diag_l=80, diag_u=20, nbOne=10, offset=50, sparsity=0.5
* \image{inline} html ../docs/figure/100_80_20_10_50_0.5.png
* <hr>
n=100, diag_l=80, diag_u=20, nbOne=10, offset=50, sparsity=0.9
* \image{inline} html ../docs/figure/100_80_20_10_50_0.9.png
* <hr>
n=100, diag_l=80, diag_u=70, nbOne=20, offset=10, sparsity=0.5
* \image{inline} html ../docs/figure/100_80_70_20_10_0.5.png
*/

/** @page contact Contact and Contributation
* <hr>
* @section developpers Developpers

The following people are involved in the development of SMG2S:

	- Xinzhe Wu (main development and algorithms)
	- Serge G. Petiton (maths and algorithms)
	- Hervé Gachlier (maths)
...

* <hr>
* <hr>
* @section contact Contact and Contributation
* Feel free to contact us: xin DOT wu AT fz BAR juelich DOT de
* <hr>
*
*/

/** @page Citing Citing SMG2S
 * If you find SMG2S useful in your project, we kindly request that you cite the following paper:
 * 
* - WU, Xinzhe, PETITON, Serge G., et LU, Yutong. 
* A parallel generator of non‐Hermitian matrices computed from given spectra. 
* Concurrency and Computation: Practice and Experience, 2020, vol. 32, no 20, p. e571.
*  <a href="https://doi.org/10.1002/cpe.5710">[doi]</a> <a href="https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/cpe.5710?casa_token=UUntHdbHvo4AAAAA:CHJa3O1_B-15_eHKY09LuWdh5TNs_trh_IXa_qDuNZLeTKcxa4CQt9WzrNsU1XSWxunknU8GeXP9Ihv9">[pdf]</a> 
*
*/

/** @page License
 * MIT License
 * 
Copyright (c) 2019 Xinzhe WU @ Maison de la Simulation, France

Copyright (c) 2019-2022, Xinzhe Wu @ Simulation and Data Laboratory Quantum 
									 Materials,  Forschungszentrum Juelich GmbH.
									 
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
*
*/

#endif