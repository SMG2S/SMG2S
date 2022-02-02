/*

MIT License

Copyright (c) 2019 Xinzhe WU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __C_WRAPPER_H__
#define __C_WRAPPER_H__

#include <mpi.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"{
#endif

/*Nilpotency Matrix C Wrapper*/
// int
struct NilpInt;
typedef struct NilpInt NilpInt_t;

NilpInt_t *newNilpInt();
void nilpInt_destory(NilpInt_t *nilp);
void nilpIntType1(NilpInt_t *nilp, int num, int size);
void nilpIntShow(struct NilpInt *nilp);
//long int
struct NilpLong;
typedef struct NilpLong NilpLong_t;

NilpLong_t *newNilpLong();
void nilpLong_destory(NilpLong_t *nilp);
void NilpLongType1(NilpLong_t *nilp, long num, long size);
void NilpLongShow(struct NilpLong *nilp);


/*parMatrixSparse C wrapper*/
// double int
struct parMatrixSparseDoubleInt;
typedef struct parMatrixSparseDoubleInt parMatrixSparseDoubleInt_t;

parMatrixSparseDoubleInt_t *newparMatrixSparseDoubleInt();
void parMatrixSparseDoubleInt_destory(parMatrixSparseDoubleInt_t *mat);
void parMatrixSparseDoubleInt_LocMatView(parMatrixSparseDoubleInt_t *mat);
void GetLocalSizeDoubleInt(parMatrixSparseDoubleInt_t *mat, int *rs, int *cs);
void parMatrixSparseDoubleInt_LocToCSR(parMatrixSparseDoubleInt_t *mat);
void parMatrixSparseDoubleInt_LocGetCSRSize(parMatrixSparseDoubleInt_t *mat, int *size, int *size2);
void parMatrixSparseDoubleInt_LocGetCSRArrays(parMatrixSparseDoubleInt_t *mat, int size, int size2, int **rows, int **cols, double **vals);
void parMatrixSparseDoubleInt_smg2s(parMatrixSparseDoubleInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);
void parMatrixSparseDoubleInt_nonsym_smg2s(parMatrixSparseDoubleInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);


// double long
struct parMatrixSparseDoubleLong;
typedef struct parMatrixSparseDoubleLong parMatrixSparseDoubleLong_t;

parMatrixSparseDoubleLong_t *newparMatrixSparseDoubleLong();
void parMatrixSparseDoubleLong_destory(parMatrixSparseDoubleLong_t *mat);
void parMatrixSparseDoubleLong_LocMatView(parMatrixSparseDoubleLong_t *mat);
void GetLocalSizeDoubleLong(parMatrixSparseDoubleLong_t *mat, __int64_t *rs, __int64_t *cs);
void parMatrixSparseDoubleLong_LocToCSR(parMatrixSparseDoubleLong_t *mat);
void parMatrixSparseDoubleLong_LocGetCSRSize(parMatrixSparseDoubleLong_t *mat, __int64_t *size, __int64_t *size2);
void parMatrixSparseDoubleLong_LocGetCSRArrays(parMatrixSparseDoubleLong_t *mat, __int64_t size, __int64_t size2, __int64_t **rows, __int64_t **cols, double **vals);
void parMatrixSparseDoubleLong_smg2s(parMatrixSparseDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm);
void parMatrixSparseDoubleLong_nonsym_smg2s(parMatrixSparseDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm);


// float int
struct parMatrixSparseFloatInt;
typedef struct parMatrixSparseFloatInt parMatrixSparseFloatInt_t;

parMatrixSparseFloatInt_t *newparMatrixSparseFloatInt();
void parMatrixSparseFloatInt_destory(parMatrixSparseFloatInt_t *mat);
void parMatrixSparseFloatInt_LocMatView(parMatrixSparseFloatInt_t *mat);
void GetLocalSizeFloatInt(parMatrixSparseFloatInt_t *mat, int *rs, int *cs);
void parMatrixSparseFloatInt_LocToCSR(parMatrixSparseFloatInt_t *mat);
void parMatrixSparseFloatInt_LocGetCSRSize(parMatrixSparseFloatInt_t *mat, int *size, int *size2);
void parMatrixSparseFloatInt_LocGetCSRArrays(parMatrixSparseFloatInt_t *mat, int size, int size2, int **rows, int **cols, double **vals);
void parMatrixSparseFloatInt_smg2s(parMatrixSparseFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);
void parMatrixSparseFloatInt_nonsym_smg2s(parMatrixSparseFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

// float long
struct parMatrixSparseFloatLong;
typedef struct parMatrixSparseFloatLong parMatrixSparseFloatLong_t;

parMatrixSparseFloatLong_t *newparMatrixSparseFloatLong();
void parMatrixSparseFloatLong_destory(parMatrixSparseFloatLong_t *mat);
void parMatrixSparseFloatLong_LocMatView(parMatrixSparseFloatLong_t *mat);
void GetLocalSizeFloatLong(parMatrixSparseFloatLong_t *mat, __int64_t *rs, __int64_t *cs);
void parMatrixSparseFloatLong_LocToCSR(parMatrixSparseFloatLong_t *mat);
void parMatrixSparseFloatLong_LocGetCSRSize(parMatrixSparseFloatLong_t *mat, __int64_t *size, __int64_t *size2);
void parMatrixSparseFloatLong_LocGetCSRArrays(parMatrixSparseFloatLong_t *mat, __int64_t size, __int64_t size2, __int64_t **rows, __int64_t **cols, double **vals);
void parMatrixSparseFloatLong_smg2s(parMatrixSparseFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm);
void parMatrixSparseFloatLong_nonsym_smg2s(parMatrixSparseFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm);



///
// complex double int
struct parMatrixSparseCmplxDoubleInt;
typedef struct parMatrixSparseCmplxDoubleInt parMatrixSparseCmplxDoubleInt_t;

parMatrixSparseCmplxDoubleInt_t *newparMatrixSparseCmplxDoubleInt();
void parMatrixSparseCmplxDoubleInt_destory(parMatrixSparseCmplxDoubleInt_t *mat);
void parMatrixSparseCmplxDoubleInt_LocMatView(parMatrixSparseCmplxDoubleInt_t *mat);
void GetLocalSizeCmplxDoubleInt(parMatrixSparseCmplxDoubleInt_t *mat, int *rs, int *cs);
void parMatrixSparseCmplxDoubleInt_LocToCSR(parMatrixSparseCmplxDoubleInt_t *mat);
void parMatrixSparseCmplxDoubleInt_LocGetCSRSize(parMatrixSparseCmplxDoubleInt_t *mat, int *size, int *size2);
void parMatrixSparseCmplxDoubleInt_LocGetCSRArrays(parMatrixSparseCmplxDoubleInt_t *mat, int size, int size2, int **rows, int **cols, double **vals);
void parMatrixSparseCmplxDoubleInt_smg2s(parMatrixSparseCmplxDoubleInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);


//complex double long
struct parMatrixSparseCmplxDoubleLong;
typedef struct parMatrixSparseCmplxDoubleLong parMatrixSparseCmplxDoubleLong_t;

parMatrixSparseCmplxDoubleLong_t *newparMatrixSparseCmplxDoubleLong();
void parMatrixSparseCmplxDoubleLong_destory(parMatrixSparseCmplxDoubleLong_t *mat);
void parMatrixSparseCmplxDoubleLong_LocMatView(parMatrixSparseCmplxDoubleLong_t *mat);
void GetLocalSizeCmplxDoubleLong(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t *rs, __int64_t *cs);
void parMatrixSparseCmplxDoubleLong_LocToCSR(parMatrixSparseCmplxDoubleLong_t *mat);
void parMatrixSparseCmplxDoubleLong_LocGetCSRSize(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t *size, __int64_t *size2);
void parMatrixSparseCmplxDoubleLong_LocGetCSRArrays(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t size, __int64_t size2, __int64_t **rows, __int64_t **cols, double **vals);
void parMatrixSparseCmplxDoubleLong_smg2s(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm);

// complex float int
struct parMatrixSparseCmplxFloatInt;
typedef struct parMatrixSparseCmplxFloatInt parMatrixSparseCmplxFloatInt_t;

parMatrixSparseCmplxFloatInt_t *newparMatrixSparseCmplxFloatInt();
void parMatrixSparseCmplxFloatInt_destory(parMatrixSparseCmplxFloatInt_t *mat);
void parMatrixSparseCmplxFloatInt_LocMatView(parMatrixSparseCmplxFloatInt_t *mat);
void GetLocalSizeCmplxFloatInt(parMatrixSparseCmplxFloatInt_t *mat, int *rs, int *cs);
void parMatrixSparseCmplxFloatInt_LocToCSR(parMatrixSparseCmplxFloatInt_t *mat);
void parMatrixSparseCmplxFloatInt_LocGetCSRSize(parMatrixSparseCmplxFloatInt_t *mat, int *size, int *size2);
void parMatrixSparseCmplxFloatInt_LocGetCSRArrays(parMatrixSparseCmplxFloatInt_t *mat, int size, int size2, int **rows, int **cols, double **vals);
void parMatrixSparseCmplxFloatInt_smg2s(parMatrixSparseCmplxFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);


// complex float long
struct parMatrixSparseCmplxFloatLong;
typedef struct parMatrixSparseCmplxFloatLong parMatrixSparseCmplxFloatLong_t;

parMatrixSparseCmplxFloatLong_t *newparMatrixSparseCmplxFloatLong();
void parMatrixSparseCmplxFloatLong_destory(parMatrixSparseCmplxFloatLong_t *mat);
void parMatrixSparseCmplxFloatLong_LocMatView(parMatrixSparseCmplxFloatLong_t *mat);
void GetLocalSizeCmplxFloatLong(parMatrixSparseCmplxFloatLong_t *mat, __int64_t *rs, __int64_t *cs);
void parMatrixSparseCmplxFloatLong_LocToCSR(parMatrixSparseCmplxFloatLong_t *mat);
void parMatrixSparseCmplxFloatLong_LocGetCSRSize(parMatrixSparseCmplxFloatLong_t *mat, __int64_t *size, __int64_t *size2);
void parMatrixSparseCmplxFloatLong_LocGetCSRArrays(parMatrixSparseCmplxFloatLong_t *mat, __int64_t size, __int64_t size2, __int64_t **rows, __int64_t **cols, double **vals);
void parMatrixSparseCmplxFloatLong_smg2s(parMatrixSparseCmplxFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm);


#ifdef __cplusplus
}
#endif




#endif
