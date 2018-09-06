#include "interface/C/c_wrapper.h"
#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

	int size,rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*create Nilpotency object*/
	struct NilpotencyInt *n;
  /*create Instance*/
	n = newNilpotencyInt();
  /*setup Nilpotency object*/
	NilpType1(n, 2, 10);
	if(rank == 0){
		showNilpotencyInt(n);
	}

  /*Create parMatrixSparse Object*/
	struct parMatrixSparseRealDoubleInt *m;
  /*create Instance*/
	m = newParMatrixSparseRealDoubleInt();
  /*Generate by SMG2S*/
	smg2sRealDoubleInt(m, 10, n, 3 ," ",MPI_COMM_WORLD);
  /*Matrix View*/
	LOC_MatViewRealDoubleInt(m);

  /*Release Nilpotency Object
  Release parMatrixSparse Object*/
	ReleaseNilpotencyInt(&n);
	ReleaseParMatrixSparseRealDoubleInt(&m);

	MPI_Finalize();
	return 0;
}


