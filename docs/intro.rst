Introduction
===============

**SMG2S** is an implementation of the Sparse Matrix Generator with Given Spectrum (SMG2S). It is implemented based on MPI (Message Passing Interface) and C++, which is able to generate very large-scale non-Hermitian/Symmetric Sparse matrices in parallel on modern supercomputers. 

The idea of creating a sparse matrix generator came from the fact that the spectrum of matrix have large impacts on the convergence behaviour of iterative linear solvers, such as the Krylov subspace method. Generating very large sparse with given spectrum would be beneficial both for the study/research on the numerical methods and benchmarking of the parallel performance of existing iterative solvers on supercomputers.
