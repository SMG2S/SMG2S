#include "../../utils/csrMatrix.h"
#include <cusparse_v2.h>

int csr2d_csr(csrMatrix<double,int> A, d_csrMatrix<double,int> *B);

int d_csr2csr(d_csrMatrix<double,int> B, csrMatrix<double,int> *A);

int d_spgemm_d(d_csrMatrix<double,int> A, d_csrMatrix<double,int> B, d_csrMatrix<double,int> *C, cusparseHandle_t hndl);

