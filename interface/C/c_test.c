#include "c_wrapper.h"
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

	int size,rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	struct NilpotencyInt *n;
	n = newNilpotencyInt();
	NilpType1(n, 2, 10);
	if(rank == 0){
		showNilpotencyInt(n);
	}

  int rs, cs;
	int *rows, *cols;
	int size_row, size_col;
	double *real, *imag;

	#if defined (__USE_COMPLEX__)

	struct parMatrixSparseComplexDoubleInt *m;
	m = newParMatrixSparseComplexDoubleInt();
	smg2sComplexDoubleInt(m, 10, n, 3 ," ", MPI_COMM_WORLD);
	LOC_MatViewComplexDoubleInt(m);
	GetLocalSizeComplexDoubleInt(m,&rs, &cs);
	Loc_ConvertToCSRComplexDoubleInt(m);
	MPI_Barrier(MPI_COMM_WORLD);
	Loc_CSRGetRowsArraySizes(m, &size_row, &size_col);
	printf("size 1 = %d, SIZE 2 = %d\n", size_row,size_col);

	#else

	struct parMatrixSparseDoubleInt *m;
	m = newParMatrixSparseDoubleInt();
	smg2sDoubleInt(m, 10, n, 3 ," ");
	LOC_MatViewDoubleInt(m);
	GetLocalSizeDoubleInt(m,&rs, &cs);
	Loc_ConvertToCSRDoubleInt(m);
	ReleaseParMatrixSparseDoubleInt(&m);

	#endif
	ReleaseNilpotencyInt(&n);

	MPI_Finalize();
	return 0;
}
