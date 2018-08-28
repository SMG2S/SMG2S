%module smg2s
%{
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <sstream> 
#include <string>
#include <map>
#include "include/utils/utils.h"
#include "include/parVector/parVectorMap.h"
#include <mpi.h>
#include <fstream>
#include "include/config/config.h"
#include <vector>
#include "include/parMatrix/MatrixCSR.h"
#include "include/parVector/parVectorMap.h"
#include "include/parVector/parVector.h"
#include "include/config/config.h"
#include "include/utils/MPI_DataType.h"
#include "include/parMatrix/parMatrixSparse.h"
#include "include/smg2s/specGen.h"
#include "include/smg2s/smg2s.h"
%}

%include <stl.i>
%include "include/utils/utils.h"
%include "include/parVector/parVectorMap.h"
%include "include/parMatrix/MatrixCSR.h"
%include "include/config/config.h"
%include "include/parVector/parVector.h"
%include "include/utils/MPI_DataType.h"
%include "include/parMatrix/parMatrixSparse.h"
%include "include/smg2s/specGen.h"
%include "include/smg2s/smg2s.h"

%include "mpi4py/mpi4py.i"
%include "mpi4py/mpi4py.h"

%mpi4py_typemap(Comm, MPI_Comm);

%template(NilpotencyInt) Nilpotency<int>;

%template(parMatrixSparseRealDoubleInt) parMatrixSparse<double,int>;
%template(smg2sRealDoubleInt) smg2s<double,int>;

%template(parMatrixSparseRealDoubleLongInt) parMatrixSparse<double,__int64_t>;
%template(smg2sRealDoubleLongInt) smg2s<double,__int64_t>;

%template(parMatrixSparseRealSingleInt) parMatrixSparse<float,int>;
%template(smg2sRealSingleInt) smg2s<float,int>;

%template(parMatrixSparseRealSingleLongInt) parMatrixSparse<float,__int64_t>;
%template(smg2sRealSingleLongInt) smg2s<float,__int64_t>;

%template(parMatrixSparseComplexDoubleInt) parMatrixSparse<std::complex<double>,int>;
%template(smg2sComplexDoubleInt) smg2s<std::complex<double>,int>;

%template(parMatrixSparseComplexDoubleLongInt) parMatrixSparse<std::complex<double>,__int64_t>;
%template(smg2sComplexDoubleLongInt) smg2s<std::complex<double>,__int64_t>;

%template(parMatrixSparseComplexSingleInt) parMatrixSparse<std::complex<float>,int>;
%template(smg2sComplexSingleInt) smg2s<std::complex<float>,int>;

%template(parMatrixSparseComplexSingleLongInt) parMatrixSparse<std::complex<float>,__int64_t>;
%template(smg2sComplexSingleLongInt) smg2s<std::complex<float>,__int64_t>;
