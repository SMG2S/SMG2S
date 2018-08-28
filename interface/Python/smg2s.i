%module smg2s
%{
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <sstream> 
#include <string>
#include <map>
#include "../../utils/utils.h"
#include "../../parVector/parVectorMap.h"
#include <mpi.h>
#include <fstream>
#include "../../config/config.h"
#include <vector>
#include "../../parMatrix/MatrixCSR.h"
#include "../../parVector/parVectorMap.h"
#include "../../parVector/parVector.h"
#include "../../config/config.h"
#include "../../utils/MPI_DataType.h"
#include "../../parMatrix/parMatrixSparse.h"
#include "../../smg2s/specGen.h"
#include "../../smg2s/smg2s.h"
%}

%include <stl.i>
%include "../../utils/utils.h"
%include "../../parVector/parVectorMap.h"
%include "../../parMatrix/MatrixCSR.h"
%include "../../config/config.h"
%include "../../parVector/parVector.h"
%include "../../utils/MPI_DataType.h"
%include "../../parMatrix/parMatrixSparse.h"
%include "../../smg2s/specGen.h"
%include "../../smg2s/smg2s.h"

%include "mpi4py/mpi4py.i"
%include "mpi4py/mpi4py.h"

%mpi4py_typemap(Comm, MPI_Comm);
%template(NilpotencyInt) Nilpotency<int>;
%template(MatrixCSRDoubleInt) MatrixCSR<double,int>;
%template(parVectorMapInt) parVectorMap<int>;
%template(parVectorDoubleInt) parVector<double,int>;
%template(parMatrixSparseDoubleInt) parMatrixSparse<double,int>;

%template(smg2sDoubleInt) smg2s<double,int>;


