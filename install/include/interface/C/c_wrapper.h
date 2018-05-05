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

#include "../../config/config.h"
#include <mpi.h>

struct NilpotencyInt;

struct parVectorMapInt;

struct parMatrixSparseDoubleInt;

struct parMatrixSparseComplexDoubleInt;


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

#if defined (__USE_COMPLEX__)
/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseComplexDoubleInt *newParMatrixSparseComplexDoubleInt(void);
void ReleaseParMatrixSparseComplexDoubleInt(struct parMatrixSparseComplexDoubleInt **ppInstance);
void LOC_MatViewComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m);
void GetLocalSizeComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int *rs, int *cs);
void Loc_ConvertToCSRComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m);
void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexDoubleInt *m, int *size, int *size2);
void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexDoubleInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag);

/*SMG2S C wrapper*/
void smg2sComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

#else
/*parMatrixSparse double int C wrapper*/
struct parMatrixSparseDoubleInt *newParMatrixSparseDoubleInt(void);
void ReleaseParMatrixSparseDoubleInt(struct parMatrixSparseDoubleInt **ppInstance);
void LOC_MatViewDoubleInt(struct parMatrixSparseDoubleInt *m);
void GetLocalSizeDoubleInt(struct parMatrixSparseDoubleInt *m, int *rs, int *cs);
void Loc_ConvertToCSRDoubleInt(struct parMatrixSparseDoubleInt *m);

void Loc_RealCSRGetRowsArraySizes(struct parMatrixSparseDoubleInt *m, int *size, int *size2);
void Loc_RealCSRGetRowsArrays(struct parMatrixSparseDoubleInt *m, int size, int **rows, int size2, int **cols, double **vals);

/*SMG2S C wrapper*/
void smg2sDoubleInt(struct parMatrixSparseDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm);

#endif

#ifdef __cplusplus
};
#endif

#endif
