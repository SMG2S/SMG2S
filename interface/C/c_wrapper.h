/*
   This file is part of SMG2S.
   Author(s): Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr or xinzhe.wu1990@gmail.com>
        Date: 2018-04-20
   Copyright (C) 2018-     Xinzhe WU

   SMG2S is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   SMG2S is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with SMG2S.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __C_WRAPPER_H__
#define __C_WRAPPER_H__

#include "../../parMatrix/parMatrixSparse.h"
#include "../../config/config.h"
#include <mpi.h>
#include <stdint.h>


struct NilpotencyInt;

struct parVectorMapInt;

//complex
struct parMatrixSparseComplexInt;

//double
struct parMatrixSparseDoubleInt;

//double + complex + int64
struct parMatrixSparseComplexDoubleLongInt;

//double+int64
struct parMatrixSparseDoubleLongInt;
//double + complex
struct parMatrixSparseComplexDoubleInt;

//int64 + complex
struct parMatrixSparseComplexLongInt;

//int64
struct parMatrixSparseLongInt;

// *void*
struct parMatrixSparseInt;


#ifdef __cplusplus
extern "C" {
#endif

/*Nilpotency Matrix C Wrapper*/
struct NilpotencyInt *newNilpotencyInt(void);
void ReleaseNilpotencyInt(struct NilpotencyInt **ppInstance);
//setup NilpotencyInt Type 1 and Type 2
extern void NilpType1(struct NilpotencyInt *n, int num, int size);
extern void NilpType2(struct NilpotencyInt *n, int num, int size);
void showNilpotencyInt(struct NilpotencyInt *n);

/*parVectorMap C wrapper*/
struct parVectorMapInt *newparVectorMapInt(void);

#if defined(__USE_COMPLEX__) && defined(__USE_DOUBLE__) && defined (__USE_64BIT__)
//complex double long int

struct parMatrixSparseComplexDoubleLongInt *newParMatrixSparseComplexDoubleLongInt(void);
void ReleaseParMatrixSparseComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt **ppInstance);
void LOC_MatViewComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m);
void GetLocalSizeComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, int *rs, int *cs);
void Loc_ConvertToCSRComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m);

void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexDoubleLongInt *m, int *size,int *size2);
void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexDoubleLongInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag);

void smg2sComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

#elif defined (__USE_COMPLEX__) && defined(__USE_DOUBLE__)
//complex double int

struct parMatrixSparseComplexDoubleInt *newParMatrixSparseComplexDoubleInt(void);
void ReleaseParMatrixSparseComplexDoubleInt(struct parMatrixSparseComplexDoubleInt **ppInstance);
void LOC_MatViewComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m);
void GetLocalSizeComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int *rs, int *cs);
void Loc_ConvertToCSRComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m);

void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexDoubleInt *m, int *size, int *size2);
void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexDoubleInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag);

/*SMG2S C wrapper*/
void smg2sComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);


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
void GetLocalSizeDoubleLongInt(struct parMatrixSparseDoubleLongInt *m, __int64_t *rs, __int64_t *cs);
void Loc_ConvertToCSRDoubleLongInt(struct parMatrixSparseDoubleLongInt *m);
void Loc_CSRGetRowsArraySizes(struct parMatrixSparseDoubleLongInt *m, __int64_t *size,__int64_t *size2);
void Loc_CSRGetRowsArrays(struct parMatrixSparseDoubleLongInt *m, int size, int **rows, int size2, int **cols, double **vals);
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


#endif

#ifdef __cplusplus
};
#endif

#endif
