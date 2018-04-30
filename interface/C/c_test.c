#include "c_wrapper.h"
#include <stdio.h>
#include <mpi.h>

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
	#if defined (__USE_COMPLEX__)

	struct parMatrixSparseComplexDoubleInt *m;
	m = newParMatrixSparseComplexDoubleInt();
	smg2sComplexDoubleInt(m, 10, n, 3 ," ");
	LOC_MatViewComplexDoubleInt(m);

	ReleaseParMatrixSparseComplexDoubleInt(&m);

	#else

	struct parMatrixSparseDoubleInt *m;
	m = newParMatrixSparseDoubleInt();
	smg2sDoubleInt(m, 10, n, 3 ," ");
	LOC_MatViewDoubleInt(m);

	ReleaseParMatrixSparseDoubleInt(&m);

	#endif



	ReleaseNilpotencyInt(&n);

	MPI_Finalize();
	return 0;
}
