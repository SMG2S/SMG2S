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
* \image html ../docs/figure/matgen.png
* <hr>
* @section notes Release Notes
* This is the release note
* <hr>
* @section Requirements
* @verbinclude requirements
* <hr> 
* Some requirements
*
*/

/** @page started User guide
* <hr>
* @section quick Quick start
* <hr>
* @section install Installation
* <hr>
* @section quickexample Quick example
* <hr>
* @section custom Customization
* @subsection nilp Nilpotent matrix
* @subsection initmat Initialization of matrix
* @subsection spec Set spectrum
* @subsection asm Assembling the customizations
* <hr>
* @section verify Plotting and verification
* <hr>
* @section c-interface C interface
* <hr>
* @section moreexample More examples
* <hr>
*/

/** @page Gallery Gallery: Sparsity Patterns
*
*/

/** @page contact Contact and Contributation
* <hr>
* @section developpers Developpers
* <hr>
* <hr>
* @section contact Contact and Contributation
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