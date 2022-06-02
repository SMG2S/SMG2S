# Sparse Matrix Generator with Given Spectrum

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2694506.svg)](https://doi.org/10.5281/zenodo.2694506)

-------------------------------------------------------------------------------

## Overview

Author [Xinzhe Wu](https://brunowu.github.io) @ [Maison de la Simulation](http://www.maisondelasimulation.fr), France (2016-2019).

                                  @ [SDL Quantum Materials](https://www.fz-juelich.de/en/ias/jsc/about-us/structure/simulation-and-data-labs/sdl-quantum-materials), Forschungszentrum Juelich GmbH, Germany (2019-present).

****

### What is SMG2S?

**SMG2S** is able to generate large-scale non-Hermitian and non-Symmetric matrices in parallel with the spectral distribution functions or eigenvalues given by users. SMG2S can be used to benchmark the iterative solvers for both linear systems and eigenvalue problems on supercomputers using the generated very large test matrices with customized spectral properties.

As a matrix generator, SMG2S provides:

- generating of both Non-Hermitian and Non-Symmetric sparse matrix

- generated matrices are naturally Sparse with non-trivial sparsity pattern

- Given Spectrum: the spectrum of generated matrix is the same as the one specified by the users

- Sparsity Patterns are diverse and controllable


As a software, SMG2S provides:

* a collection of C++ header only files
* C++ templated implementation for different data type
* parallel implementation based on [[MPI]](https://en.wikipedia.org/wiki/Message_Passing_Interface) which is able to efficiently generate very large sparse matrices in parallel on supercomputers
* an easy-to-use C interface
* a verification module based on Python for the sparsity pattern plotting and spectrum verification of small size of generated matrix.
* Efficient parallel IO to store the generated matrix into [MatrixMarket format](https://math.nist.gov/MatrixMarket/formats.html)


![Matrix Generation Pattern](https://github.com/SMG2S/SMG2S/tree/devel/docs/figure/matgen.png)

### Cite SMG2S

If you find SMG2S useful in your project, we kindly request that you cite the following paper:

Wu, Xinzhe, Serge G. Petiton, and Yutong Lu. "A Parallel Generator of Non-Hermitian Matrices computed from Given Spectra." Concurrency and Computation: Practice and Experience, 32(20), e5710, 2020. [[DOI]](https://doi.org/10.1002/cpe.5710) [[PDF]](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/cpe.5710?casa_token=UUntHdbHvo4AAAAA:CHJa3O1_B-15_eHKY09LuWdh5TNs_trh_IXa_qDuNZLeTKcxa4CQt9WzrNsU1XSWxunknU8GeXP9Ihv9)

### Gallery: Sparsity Patterns

Please refer to [docs/gallery](https://github.com/SMG2S/SMG2S/tree/devel/docs/gallery) for more examples.

### Contact and Contributation

Feel free to contact by email address: **xin DOT wu AT fz BAR juelich DOT de***

## Documentation

### Getting SMG2S

SMG2S is able to available on the Github. The most updated version of SMG2S can be gotten either by the following `git` command:

```bash
git clone https://github.com/SMG2S/SMG2S.git
```

Moreover a released version can be downloaded [here](http)

### Dependencies

SMG2S is developed in C++14 and MPI, and it is compiled with CMake. So the following software and compiler should be available before the installation of SMG2S.

1. a C++ compiler with C++14 support

2. MPI: message passing interface

3. CMake: version >= 3.8

### Quick start

SMG2S provides an executable `smg2s.exe` that the users can compile and start to play with SMG2S without installation as follows. 

```bash
cd SMG2S
mkdir build & cd build
cmake .. 
make -j
```

Then the executable `smg2s.exe`is available, and it can be run as follows:

```bash
  mpirun -np ${PROCS} ./smg2s.exe -D ${dim} -L ${diag_l} -U ${diag_u} -O ${offset} -C ${nbOne} -S ${sparsity} -M {no-herm or non-symm}
```

in which the command line parsers provides the customization of following parameters:

```bash
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
```

### Installation

SMG2S relies on CMake for compiling and installation. A CMake flag `CMAKE_INSTALL_PREFIX` should be provided for the path of installation.

```bash
cd SMG2S
mkdir build & cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PATH_TO_INSTALL}
make -j install
```

### Use SMG2S with own project

#### header-only

SMG2S is a collection of C++ header files. If users want to use SMG2S with C++, they can just copy SMG2S headers into their project.

#### CMake

SMG2S is installed as a CMake package, and it can be detected by the classic CMake `find_package` command. If the installation path is not in the default searching path of CMake, a CMake flag `CMAKE_PREFIX_PATH` should be provided which links to the installation path of SMG2S.

So in your own project which want to use SMG2S:

```bash
mkdir build & cd build
cmake .. -DCMAKE_PREFIX_PATH=${INSTALLED_PATH_OF_SMG2S}
make -j
```

and in the `CMakeLists.txt` of own project, it should provide some content as follows:

```cmake
cmake_minimum_required(VERSION 3.8)
project(YOUR-OWN-PROJECT)
#find installation of SMG2S
find_package( smg2s REQUIRED CONFIG)
# for C++ code
add_executable(smg2s-app test_parMatrix.cpp)
target_link_libraries(smg2s-app PUBLIC SMG2S::smg2s)
# for C-interface code
add_executable(test_c.exe test_c.c)
target_link_libraries(test_c.exe PRIVATE SMG2S::smg2s2c)
```

In case that the support of `C++14` is disabled by some `C++` compiler, please insert also the following lines into your `CMakeLists.txt` before the usage of SMG2S.

```cmake
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
     message([FATAL_ERROR] "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()
```

### Usage

#### Parallel vector and sparse matrix

#### Building blocks SMG2S

##### 1. Distribution of sparse matrix and vector

##### 2. User-provided spectrum

##### 3. Nilpotent matrix

##### 4. Initial matrix

SMG2S provides the generation of matrices in three different categories:

1. non-Hermtian matrices with complex eigenvalues

2. non-Symmetric matrices with real eigenvalues

3. non-Symmetric matrices with conjugated eigenvalues

For each categories, SMG2S provides three functions which allows the users having different levels of controlling and customizing the properties of generated matrices.

1. the first level, users need to provide multiple simple parameters and a local text file containing the eigenvalues. For the format of local file, please check the next section.

2. the second level, 

3. the third level

### Format of Given Spectrum Files

SMG2S is able to load user-provided spectrum in parallel from local text files. However, the provided files should conform into a specific format.

1. The first line is the comment part which includes the scalar types of given spectrum. This line should be: `%%SMG2S vector in complex scalar` and `%%SMG2S vector in real scalar` for the eigenvalues in complex or real scalar type, explicitly. **Attention**, for this line, the key word `complex` or `real` should always be there and conform with the type of user-provided spectrum. The parallel IO of SMG2S queries at first this line to check if the provided eigenvalues are complex or real.

2. The second line indicates the number of given eigenvaues in the files. For the ones with `3` complex values, it is `3 3 3`, and for the ones with `3` real eigenvalues, it should be `3 3`.

3. Starting from the `3rd` line, it is the main content of this file. It can have either `2` or `3` columns, which depends on the scalar types of eigenvalues. For the case with complex values, the first column indicates the coordinates for each eigenvalue, the second column contains the real part of eigenvalues, and the third column is for the imaginary part of eigenvalues. For the case with real values, the two columns contain the indexing and values of eigenvalues, respectively. **Attention**, the indexing is `1`-based, rather than `0`-based.  

##### Real eigenvalues for non-Symmetric matrices

For the case with real eigenvalues for non-Symmetric matrices, the given spectrum file format should be as follows:

```
%%SMG2S vector in real scalar
3 3 
1 10
2 3.4790
3 5.0540
```

##### Complex eigenvalues for non-Hermtian matrices

For the complex values for non-Hermitian matrices which are not supposed to be conjugated, the given spectrum is stored in three columns, the first column is the coordinates, the second column is the real part of complex values, and the third column is the imaginary part of complex values. Here is an example with `3` eigenvalues:

    %%SMG2S vector in complex scalar
    3 3 3
    1 10 6.5154
    2 10.6288 3.4790
    3 10.7621 5.0540

##### Conjugate eigenvalues for non-Symmetric matrices

For the non-Symmetric matrices whose entries are all in real scalar, they can have conjugate eigenvalues which are in complex scalar. So in order to generate non-Symmetric test matrices with given conjugated eigenvalues, the give spectrum are always stored in complex form, with three columns.

##### Attention:

For the non-Symmetric matrices, if one eigenvalue is complex, there is another value that they two are symmetric to the real axis in the real-imaginary plain, this is their conjugated eigenvalue. So when setting up the spectral file, one eigenvalue `a+bi` with `b != 0` should be closely followed by another eigenvalue `a-bi`. For the eigenvalues with their imaginary part to be `0`, they are stored with their imaginary part being 0. Here is an example

    %%SMG2S vector in complex scalar
    9 9 9
    1 10.6288 -3.4790
    2 10.6288 3.4790
    3 2.332 0
    4 10.7621 5.0540
    5 10.7621 -5.0540
    6 -2.332 0
    7 -11.02 0
    8 21.21 4.4
    9 21.21 -4.4

### Interface

#### Interface to C

A basic example of usge:

```c
#include "interface/C/c_wrapper.h"
#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
}
```

### 

![](https://github.com/SMG2S/SMG2S/tree/devel/docs/figure/verification_5.png)

![](https://github.com/SMG2S/SMG2S/tree/devel/docs/figure/verification_4.png)

![](https://github.com/SMG2S/SMG2S/tree/devel/docs/figure/verification_1.png)

![](https://github.com/SMG2S/SMG2S/tree/devel/docs/figure/verification_2.png)

![](https://github.com/SMG2S/SMG2S/tree/devel/docs/figure/verification_3.png)

