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
#if defined(__USE_COMPLEX__) && defined(__USE_DOUBLE__) && defined (__USE_64BIT__)
//complex double long int

	struct parMatrixSparseComplexDoubleLongInt *m;
	m = newParMatrixSparseComplexDoubleLongInt();
	smg2sComplexDoubleLongInt(m, 10, n, 3 ," ", MPI_COMM_WORLD);
	LOC_MatViewComplexDoubleLongInt(m);
	GetLocalSizeComplexDoubleLongInt(m,&rs, &cs);
	Loc_ConvertToCSRComplexDoubleLongInt(m);
	MPI_Barrier(MPI_COMM_WORLD);
	Loc_CSRGetRowsArraySizes(m, &size_row, &size_col);
	printf("size 1 = %d, SIZE 2 = %d\n", size_row,size_col);
void ReleaseParMatrixSparseComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt **ppInstance);
void LOC_MatViewComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m);
void GetLocalSizeComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, int *rs, int *cs);
void Loc_ConvertToCSRComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m);

void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexDoubleLongInt *m, int *size,int *size2);
void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexDoubleLongInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag);

void smg2sComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

#elif defined (__USE_COMPLEX__) && defined(__USE_DOUBLE__)
//complex double int

	struct parMatrixSparseComplexDoubleInt *m;
	m = newParMatrixSparseComplexDoubleInt();
	smg2sComplexDoubleInt(m, 10, n, 3 ," ", MPI_COMM_WORLD);
	LOC_MatViewComplexDoubleInt(m);
	GetLocalSizeComplexDoubleInt(m,&rs, &cs);
	Loc_ConvertToCSRComplexDoubleInt(m);
	MPI_Barrier(MPI_COMM_WORLD);
	Loc_CSRGetRowsArraySizes(m, &size_row, &size_col);
	printf("size 1 = %d, SIZE 2 = %d\n", size_row,size_col);

#elif defined (__USE_COMPLEX__) && defined(__USE_64BIT__)
//complex  single long int

struct parMatrixSparseComplexLongInt *newParMatrixSparseComplexLongInt(void);
void ReleaseParMatrixSparseComplexLongInt(struct parMatrixSparseComplexLongInt **ppInstance);
void LOC_MatViewComplexLongInt(struct parMatrixSparseComplexLongInt *m);
void GetLocalSizeComplexLongInt(struct parMatrixSparseComplexLongInt *m, int *rs, int *cs);
void Loc_ConvertToCSRComplexLongInt(struct parMatrixSparseComplexLongInt *m);

void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexLongInt *m, int *size,int *size2);
void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexLongInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag);

void smg2sComplexLongInt(struct parMatrixSparseComplexLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

#elif defined (__USE_DOUBLE__) && defined(__USE_64BIT__)
//real double long int

/*parMatrixSparse complex<double> int C wrapper*/

struct parMatrixSparseDoubleLongInt *newParMatrixSparseDoubleLongInt(void);
void ReleaseParMatrixSparseDoubleLongInt(struct parMatrixSparseDoubleLongInt **ppInstance);
void LOC_MatViewDoubleLongInt(struct parMatrixSparseDoubleLongInt *m);
void GetLocalSizeDoubleLongInt(struct parMatrixSparseDoubleLongInt *m, int *rs, int *cs);
void Loc_ConvertToCSRDoubleLongInt(struct parMatrixSparseDoubleLongInt *m);
void Loc_CSRGetRowsArraySizes(struct parMatrixSparseDoubleLongInt *m, int *size,int *size2);
void Loc_CSRGetRowsArrays(struct parMatrixSparseDoubleLongInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag);
void smg2sDoubleLongInt(struct parMatrixSparseDoubleLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);


#elif defined (__USE_COMPLEX__)
//complex single int

/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseComplexInt *newParMatrixSparseComplexInt(void);
void ReleaseParMatrixSparseComplexInt(struct parMatrixSparseComplexInt **ppInstance);
void LOC_MatViewComplexInt(struct parMatrixSparseComplexInt *m);

void GetLocalSizeComplexInt(struct parMatrixSparseComplexInt *m, int *rs, int *cs);
void Loc_ConvertToCSRComplexInt(struct parMatrixSparseComplexInt *m);
void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexInt *m, int *size,int *size2);
void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag);
void smg2sComplexInt(struct parMatrixSparseComplexInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

#elif defined (__USE_DOUBLE__)
//real double int

struct parMatrixSparseDoubleInt *newParMatrixSparseDoubleInt(void);
void ReleaseParMatrixSparseDoubleInt(struct parMatrixSparseDoubleInt **ppInstance);
void LOC_MatViewDoubleInt(struct parMatrixSparseDoubleInt *m);
void GetLocalSizeDoubleInt(struct parMatrixSparseDoubleInt *m, int *rs, int *cs);
void Loc_ConvertToCSRDoubleInt(struct parMatrixSparseDoubleInt *m);

void Loc_RealCSRGetRowsArraySizes(struct parMatrixSparseDoubleInt *m, int *size, int *size2);
void Loc_RealCSRGetRowsArrays(struct parMatrixSparseDoubleInt *m, int size, int **rows, int size2, int **cols, double **vals);

/*SMG2S C wrapper*/
void smg2sDoubleInt(struct parMatrixSparseDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

#elif defined (__USE_64BIT__)
//real single long int

struct parMatrixSparseLongInt *newParMatrixSparseLongInt(void);
void ReleaseParMatrixSparseLongInt(struct parMatrixSparseLongInt **ppInstance);
void LOC_MatViewLongInt(struct parMatrixSparseLongInt *m);
void GetLocalSizeLongInt(struct parMatrixSparseLongInt *m, int *rs, int *cs);
void Loc_ConvertToCSRLongInt(struct parMatrixSparseLongInt *m);

void Loc_LongCSRGetRowsArraySizes(struct parMatrixSparseLongInt *m, int *size, int *size2);
void Loc_RealCSRGetRowsArrays(struct parMatrixSparseLongInt *m, int size, int **rows, int size2, int **cols, double **vals);

void smg2sLongInt(struct parMatrixSparseLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

#else
/*parMatrixSparse int C wrapper*/
struct parMatrixSparseInt *newParMatrixSparseInt(void);
void ReleaseParMatrixSparseInt(struct parMatrixSparseInt **ppInstance);
void LOC_MatViewInt(struct parMatrixSparseInt *m);
void GetLocalSizeInt(struct parMatrixSparseInt *m, int *rs, int *cs);
void Loc_ConvertToCSRInt(struct parMatrixSparseInt *m);

void Loc_RealCSRGetRowsArraySizes(struct parMatrixSparseInt *m, int *size, int *size2);
void Loc_RealCSRGetRowsArrays(struct parMatrixSparseInt *m, int size, int **rows, int size2, int **cols, double **vals);

void smg2sInt(struct parMatrixSparseInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);












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
