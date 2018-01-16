#include "d_spgemm_d.h"

// error check macros
#define CUSPARSE_CHECK(x) {cusparseStatus_t _c=x; if (_c != CUSPARSE_STATUS_SUCCESS) {printf("cusparse fail: %d, line: %d\n", (int)_c, __LINE__); exit(-1);}}

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

int csr2d_csr(csrMatrix<double,int> A, d_csrMatrix<double,int> *B){

	B->d_row = A.row;
	B->d_column = A.column;
	B->d_nnz = A.nnz;

	cudaMalloc(&B->d_csrRowPtr, (B->d_row+1)*sizeof(int));
        cudaMalloc(&B->d_csrColInd, (B->d_nnz)*sizeof(int));
        cudaMalloc(&B->d_csrVal, (B->d_nnz)*sizeof(double));

	cudaMemcpy(B->d_csrRowPtr, A.csrRowPtr, (B->d_row+1)*sizeof(int), cudaMemcpyHostToDevice);
  	cudaMemcpy(B->d_csrColInd, A.csrColInd, B->d_nnz*sizeof(int), cudaMemcpyHostToDevice);
  	cudaMemcpy(B->d_csrVal, A.csrVal, B->d_nnz*sizeof(double), cudaMemcpyHostToDevice);

	return 0;
}

int d_csr2csr(d_csrMatrix<double,int> B, csrMatrix<double,int> *A){

	A->row = B.d_row;
	A->column = B.d_column;
	A->nnz = B.d_nnz;

	A->csrRowPtr = (int *)malloc((B.d_row+1)*sizeof(int));
	A->csrColInd = (int *)malloc(B.d_nnz *sizeof(int));
	A->csrVal  = (double *)malloc(B.d_nnz *sizeof(double));

	cudaMemcpy(A->csrRowPtr, B.d_csrRowPtr, (B.d_row+1)*sizeof(int), cudaMemcpyDeviceToHost);
  	cudaMemcpy(A->csrColInd, B.d_csrColInd,  B.d_nnz*sizeof(int), cudaMemcpyDeviceToHost);
  	cudaMemcpy(A->csrVal, B.d_csrVal, B.d_nnz*sizeof(double), cudaMemcpyDeviceToHost);
  	cudaCheckErrors("cudaMemcpy fail");

	return 0;
}

int d_spgemm_d(d_csrMatrix<double,int> A, d_csrMatrix<double,int> B, d_csrMatrix<double,int> *C, cusparseHandle_t hndl)
{
	int baseC;
  	int *nnzTotalDevHostPtr = &C->d_nnz;
  	cusparseMatDescr_t descrA, descrB, descrC;
  	cusparseStatus_t stat;
  	CUSPARSE_CHECK(cusparseCreate(&hndl));
  	stat = cusparseCreateMatDescr(&descrA);
  	CUSPARSE_CHECK(stat);
  	stat = cusparseCreateMatDescr(&descrB);
  	CUSPARSE_CHECK(stat);
  	stat = cusparseCreateMatDescr(&descrC);
  	CUSPARSE_CHECK(stat);
  	stat = cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
  	CUSPARSE_CHECK(stat);
  	stat = cusparseSetMatType(descrB, CUSPARSE_MATRIX_TYPE_GENERAL);
  	CUSPARSE_CHECK(stat);
  	stat = cusparseSetMatType(descrC, CUSPARSE_MATRIX_TYPE_GENERAL);
  	CUSPARSE_CHECK(stat);
  	stat = cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
  	CUSPARSE_CHECK(stat);
  	stat = cusparseSetMatIndexBase(descrB, CUSPARSE_INDEX_BASE_ZERO);
  	CUSPARSE_CHECK(stat);
  	stat = cusparseSetMatIndexBase(descrC, CUSPARSE_INDEX_BASE_ZERO);
  	CUSPARSE_CHECK(stat);

  	cusparseOperation_t transA = CUSPARSE_OPERATION_NON_TRANSPOSE;
  	cusparseOperation_t transB = CUSPARSE_OPERATION_NON_TRANSPOSE;

  	stat = cusparseSetPointerMode(hndl, CUSPARSE_POINTER_MODE_HOST);

  	CUSPARSE_CHECK(stat);
  	cudaMalloc((void**)&C->d_csrRowPtr, sizeof(int)*(A.d_row+1));

  	cudaCheckErrors("cudaMalloc fail");

  	C->d_row = A.d_row;
  	C->d_column = B.d_column;
  	stat = cusparseXcsrgemmNnz(hndl, transA, transB, A.d_row, B.d_column, A.d_column,
        	descrA, A.d_nnz, A.d_csrRowPtr, A.d_csrColInd,
        	descrB, B.d_nnz, B.d_csrRowPtr, B.d_csrColInd,
        	descrC, C->d_csrRowPtr, nnzTotalDevHostPtr );

  	CUSPARSE_CHECK(stat);
	if (NULL != nnzTotalDevHostPtr){
    	C->d_nnz = *nnzTotalDevHostPtr;}
  	else{
    		cudaMemcpy(&C->d_nnz, C->d_csrRowPtr+A.d_row, sizeof(int), cudaMemcpyDeviceToHost);
    		cudaMemcpy(&baseC, C->d_csrRowPtr, sizeof(int), cudaMemcpyDeviceToHost);
    		cudaCheckErrors("cudaMemcpy fail");
    		C->d_nnz -= baseC;
	}

	//  cudaCheckErrors("cudaMalloc fail");
  	cudaMalloc((void**)&C->d_csrColInd, sizeof(int)*C->d_nnz);
  	cudaCheckErrors("cudaMalloc fail");
  	cudaMalloc((void**)&C->d_csrVal, sizeof(double)*C->d_nnz);
  	cudaCheckErrors("cudaMalloc fail");
	// perform multiplication C = A*B
  	stat = cusparseDcsrgemm(hndl, transA, transB, A.d_row, B.d_column, A.d_column,
        	descrA, A.d_nnz, A.d_csrVal, A.d_csrRowPtr, A.d_csrColInd,
        	descrB, B.d_nnz,B.d_csrVal, B.d_csrRowPtr, B.d_csrColInd,
        	descrC, C->d_csrVal, C->d_csrRowPtr, C->d_csrColInd);

 	CUSPARSE_CHECK(stat);

	return 0;
}




