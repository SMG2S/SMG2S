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
void GetLocalSizeDoubleLong(parMatrixSparseDoubleLong_t *mat, long *rs, long *cs);
void parMatrixSparseDoubleLong_LocToCSR(parMatrixSparseDoubleLong_t *mat);
void parMatrixSparseDoubleLong_LocGetCSRSize(parMatrixSparseDoubleLong_t *mat, long *size, long *size2);
void parMatrixSparseDoubleLong_LocGetCSRArrays(parMatrixSparseDoubleLong_t *mat, long size, long size2, long **rows, long **cols, double **vals);
void parMatrixSparseDoubleLong_smg2s(parMatrixSparseDoubleLong_t *mat, long probSize, struct NilpInt *nilp, long lbandwidth, char *spectrum, MPI_Comm comm);
void parMatrixSparseDoubleLong_nonsym_smg2s(parMatrixSparseDoubleLong_t *mat, long probSize, struct NilpInt *nilp, long lbandwidth, char *spectrum, MPI_Comm comm);



#ifdef __cplusplus
}
#endif




#endif
