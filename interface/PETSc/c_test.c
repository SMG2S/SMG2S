#include "../C/c_wrapper.h"
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "petscmat.h"

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

        PetscInitialize(&argc,&argv,(char *)0,NULL);

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
	smg2sComplexDoubleInt(m, 10, n, 3 ," ");
	LOC_MatViewComplexDoubleInt(m);
	GetLocalSizeComplexDoubleInt(m,&rs, &cs);
	Loc_ConvertToCSRComplexDoubleInt(m);
	MPI_Barrier(MPI_COMM_WORLD);
	Loc_CSRGetRowsArraySizes(m, &size_row, &size_col);
	printf("size 1 = %d, SIZE 2 = %d\n", size_row,size_col);
	rows = (int *) malloc(size_row*sizeof(int));
        cols = (int *) malloc(size_col*sizeof(int));

        real = (double *) malloc(size_col*sizeof(double));
        imag = (double *) malloc(size_col*sizeof(double));

        Loc_CSRGetRowsArrays(m, size_row, &rows, size_col, &cols, &real, &imag);
	printf("row = %d, cols = %d real = %f, imag = %f\n", rows[0], cols[0], real[0],imag[0]);
	//Loc_CSRGetColsArray(m, &cols, &size_col);
	//ReleaseParMatrixSparseComplexDoubleInt(&m);

	PetscInt *i, *j, q;
        PetscScalar *a;

        Mat A;

        PetscMalloc1((PetscInt)size_row, &i);
        PetscMalloc1((PetscInt)size_col, &j);
        PetscMalloc1((PetscScalar)size_col, &a);


        for(q =0; q < size_row; q++){
                i[q] = rows[q];
//		printf("proc %d, i[%d] = %d\n", rank, q, i[q]);
        }

        for(q = 0; q < size_col; q++){
                j[q] = cols[q];
                a[q] = real[q] + PETSC_i * imag[q];
//		printf("proc %d, cols[%d] = %d\n", rank, q, cols[q]);
    //            printf("proc %d, j[%d] = %d, a[%d] = %f+%f \n", rank, q, j[q], q, real[q], imag[q]);
        }

	MatCreate(PETSC_COMM_WORLD, &A);
	MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, size_row-1, size_row-1, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

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

	MatView(A, PETSC_VIEWER_STDOUT_WORLD);
        PetscFree(A);
        PetscFinalize();

	MPI_Finalize();

	return 0;
}
