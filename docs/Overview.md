# Overview

SMG2S = Sparse Matrix Generator with Given Spectrum

It is a software which provides to generator the non-Hermitian and non-Symmetric sparse Matrices with user-customized eigenvalues. The generated sparse matrices contain still 
the user-provided spectrum. SMG2S is implemented in parallel based on MPI (Message Passing Interface) and C++14 to support efficiently the generation of test matrices on distributed memory platforms.

The initial motivation of the creation of SMG2S comes from the studay of iterative methods
for solving linear systems and eigenvalue problems. Iterative linear algebra methods are important for the applications in various fields. The analysis of the iterative method behaviors is complex, and it is necessary to evaluate their convergence to solve extremely large non-Hermitian/non-symmtric eigenvalue and linear problems on parallel and/or distributed machines. This convergence depends on the properties of spectra. Thus, we propse SMG2S to generate large matrices with known spectra to benchmark these methods. The generated matrices are non-trivial with very high dimension.

As a matrix generator, SMG2S provides:

- generating of both Non-Hermitian and Non-Symmetric sparse matrix

- generated matrices are naturally Sparse with non-trivial sparsity pattern

- Given Spectrum: the spectrum of generated matrix is the same as the one specified by the users

- Sparsity Patterns are diverse and controllable


As a software, SMG2S provides:

- a collection of c++ header files

- parallel implementation based on MPI which is able to efficiently generate very large sparse matrices in parallel on supercomputers

- a verification module based on Python for the plotting and verification of small size of generated matrix.

- an easy-to-use C interface 


