SMG2S

Iterative linear algebra methods are the important parts of the overall computing
time  of  applications  in  various  fields  since  decades.   Recent  research  related
to social networking, big data, machine learning and artificial intelligence has
increased  the  necessity  for  non-hermitian  solvers  associated  with  much  larger
sparse matrices and graphs.  The analysis of the iterative method behaviors for
such problems is complex, and it is necessary to evaluate their convergence to
solve extremely large non-Hermitian eigenvalue and linear problems on parallel
and/or  distributed  machines.   This  convergence  depends  on  the  properties  of
spectra.  Then,  it is necessary to generate large matrices with known spectra
to benchmark the methods.  These matrices should be non-Hermitian and non-
trivial, with very high dimension.  A scalable parallel matrix generator SMG2S
that uses the user-defined spectrum to construct large-scale sparse matrices and
ensures their eigenvalues as the given ones with high accuracy is implemented
based on MPI and C++11.  This is the Python interface of SMG2S.


