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

struct dcomplex{
	double real;
	double imag;
};

typedef struct dcomplex dcomplex_t;

struct fcomplex{
	float real;
	float imag;
};

typedef struct fcomplex fcomplex_t;

//interface of classes and structs
/*Nilpotency Matrix C Wrapper*/
// int
struct nilp;
typedef struct nilp nilp_t;

nilp_t *newNilp_1(int nbOne, int size);
nilp_t *newNilp_2(int nbOne, int diag, int size);
nilp_t *newNilp_3(int *nilpvec, int size);
nilp_t *newNilp_4(int *nilpvec, int diag, int size);
void nilp_destory(nilp_t *nilp);
int nilp_getDegree(nilp_t *nilp);
int* nilp_getIndOfZeros(nilp_t *nilp);
void nilp_show(nilp_t *nilp);


// long int
struct nilpL;
typedef struct nilpL nilpL_t;

nilpL_t *newNilpL_1(long nbOne, long size);
nilpL_t *newNilpL_2(long nbOne, long diag, long size);
nilpL_t *newNilpL_3(long *nilpvec, long size);
nilpL_t *newNilpL_4(long *nilpvec, long diag, long size);
void nilpL_destory(nilpL_t *nilp);
long nilpL_getDegree(nilpL_t *nilp);
long* nilpL_getIndOfZeros(nilpL_t *nilp);
void nilpL_show(nilpL_t *nilp);



/*init matrix struct*/
//int
struct initMatrix;
typedef struct initMatrix initMatrix_t;

initMatrix_t *newInitMatrix_1();
initMatrix_t *newInitMatrix_2(int diagl, int diagu);
initMatrix_t *newInitMatrix_3(int diagl, int diagu, double Sparsity);
initMatrix_t *newInitMatrix_4(int diagl, int diagu, double Scale, double Sparsity);
void initMatrix_show(initMatrix_t *init);
void initMatrix_destory(initMatrix_t *init);

//long
struct initMatrixL;
typedef struct initMatrixL initMatrixL_t;

initMatrixL_t *newInitMatrixL_1();
initMatrixL_t *newInitMatrixL_2(long diagl, long diagu);
initMatrixL_t *newInitMatrixL_3(long diagl, long diagu, double Sparsity);
initMatrixL_t *newInitMatrixL_4(long diagl, long diagu, double Scale, double Sparsity);
void initMatrixL_show(initMatrixL_t *init);
void initMatrixL_destory(initMatrixL_t *init);


/*parVectorMap*/
struct parVecMap;
typedef struct parVecMap parVecMap_t;
parVecMap_t *newParVecMap(MPI_Comm ncomm, int lbound, int ubound);
MPI_Comm parVecMapGetComm(parVecMap_t *pv);
int parVecMapL2G(parVecMap_t *pv, int local_index);
int parVecMapG2L(parVecMap_t *pv, int global_index);
int parVecMapGetLocSize(parVecMap_t *pv);
int parVecMapGetGlobSize(parVecMap_t *pv);
void parVecMap_destory(parVecMap_t *pv);

struct parVecMapL;
typedef struct parVecMapL parVecMapL_t;
parVecMapL_t *newParVecMapL(MPI_Comm ncomm, long lbound, long ubound);
MPI_Comm parVecMapLGetComm(parVecMapL_t *pv);
long parVecMapLL2G(parVecMapL_t *pv, long local_index);
long parVecMapLG2L(parVecMapL_t *pv, long global_index);
long parVecMapLGetLocSize(parVecMapL_t *pv);
long parVecMapLGetGlobSize(parVecMapL_t *pv);
void parVecMapL_destory(parVecMapL_t *pv);


/*parVector*/
//double int
struct ds_parVec;
typedef struct ds_parVec ds_parVec_t;
ds_parVec_t *new_ds_ParVec_1(MPI_Comm ncomm, int lbound, int ubound);
ds_parVec_t *new_ds_ParVec_2(parVecMap_t *map);
void ds_parVec_destory(ds_parVec_t *pv);
int ds_parVecGetLowerBound(ds_parVec_t *pv);
int ds_parVecGetUpperBound(ds_parVec_t *pv);
int ds_parVecGetLocSize(ds_parVec_t *pv);
int ds_parVecGetGlobSize(ds_parVec_t *pv);
double ds_parVecGetVal(ds_parVec_t *pv, int index);
double ds_parVecGetValLoc(ds_parVec_t *pv, int lindex);
double *ds_parVecGetArray(ds_parVec_t *pv);
MPI_Comm ds_parVecGetComm(ds_parVec_t *pv);
int ds_parVecL2G(ds_parVec_t *pv, int local_index);
int ds_parVecG2L(ds_parVec_t *pv, int global_index);
void ds_parVecSetToVal(ds_parVec_t *pv, double val);
void ds_parVecView(ds_parVec_t *pv);
void ds_parVecSetVal(ds_parVec_t *pv, int index, double val);
void ds_parVecSetValLoc(ds_parVec_t *pv, int lindex, double val);
void ds_parVecAdd(ds_parVec_t *pv, ds_parVec_t *pv2);
void ds_parVecScale(ds_parVec_t *pv, double scale);
double ds_parVecDot(ds_parVec_t *pv, ds_parVec_t *pv2);
void ds_parVecReadExtVec(ds_parVec_t *pv, char* spectrum);

//double long
struct dl_parVec;
typedef struct dl_parVec dl_parVec_t;
dl_parVec_t *new_dl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound);
dl_parVec_t *new_dl_ParVec_2(parVecMap_t *map);
void dl_parVec_destory(dl_parVec_t *pv);
long dl_parVecGetLowerBound(dl_parVec_t *pv);
long dl_parVecGetUpperBound(dl_parVec_t *pv);
long dl_parVecGetLocSize(dl_parVec_t *pv);
long dl_parVecGetGlobSize(dl_parVec_t *pv);
double dl_parVecGetVal(dl_parVec_t *pv, long index);
double dl_parVecGetValLoc(dl_parVec_t *pv, long lindex);
double *dl_parVecGetArray(dl_parVec_t *pv);
MPI_Comm dl_parVecGetComm(dl_parVec_t *pv);
long dl_parVecL2G(dl_parVec_t *pv, long local_index);
long dl_parVecG2L(dl_parVec_t *pv, long global_index);
void dl_parVecSetToVal(dl_parVec_t *pv, double val);
void dl_parVecView(dl_parVec_t *pv);
void dl_parVecSetVal(dl_parVec_t *pv, long index, double val);
void dl_parVecSetValLoc(dl_parVec_t *pv, long lindex, double val);
void dl_parVecAdd(dl_parVec_t *pv, dl_parVec_t *pv2);
void dl_parVecScale(dl_parVec_t *pv, double scale);
double dl_parVecDot(dl_parVec_t *pv, dl_parVec_t *pv2);
void dl_parVecReadExtVec(dl_parVec_t *pv, char* spectrum);

//float int
struct ss_parVec;
typedef struct ss_parVec ss_parVec_t;
ss_parVec_t *new_ss_ParVec_1(MPI_Comm ncomm, int lbound, int ubound);
ss_parVec_t *new_ss_ParVec_2(parVecMap_t *map);
void ss_parVec_destory(ss_parVec_t *pv);
int ss_parVecGetLowerBound(ss_parVec_t *pv);
int ss_parVecGetUpperBound(ss_parVec_t *pv);
int ss_parVecGetLocSize(ss_parVec_t *pv);
int ss_parVecGetGlobSize(ss_parVec_t *pv);
float ss_parVecGetVal(ss_parVec_t *pv, int index);
float ss_parVecGetValLoc(ss_parVec_t *pv, int lindex);
float *ss_parVecGetArray(ss_parVec_t *pv);
MPI_Comm ss_parVecGetComm(ss_parVec_t *pv);
int ss_parVecL2G(ss_parVec_t *pv, int local_index);
int ss_parVecG2L(ss_parVec_t *pv, int global_index);
void ss_parVecSetToVal(ss_parVec_t *pv, float val);
void ss_parVecView(ss_parVec_t *pv);
void ss_parVecSetVal(ss_parVec_t *pv, int index, float val);
void ss_parVecSetValLoc(ss_parVec_t *pv, int lindex, float val);
void ss_parVecAdd(ss_parVec_t *pv, ss_parVec_t *pv2);
void ss_parVecScale(ss_parVec_t *pv, float scale);
float ss_parVecDot(ss_parVec_t *pv, ss_parVec_t *pv2);
void ss_parVecReadExtVec(ss_parVec_t *pv, char* spectrum);


//float long
struct sl_parVec;
typedef struct sl_parVec sl_parVec_t;
sl_parVec_t *new_sl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound);
sl_parVec_t *new_sl_ParVec_2(parVecMap_t *map);
void sl_parVec_destory(sl_parVec_t *pv);
long sl_parVecGetLowerBound(sl_parVec_t *pv);
long sl_parVecGetUpperBound(sl_parVec_t *pv);
long sl_parVecGetLocSize(sl_parVec_t *pv);
long sl_parVecGetGlobSize(sl_parVec_t *pv);
float sl_parVecGetVal(sl_parVec_t *pv, long index);
float sl_parVecGetValLoc(sl_parVec_t *pv, long lindex);
float *sl_parVecGetArray(sl_parVec_t *pv);
MPI_Comm sl_parVecGetComm(sl_parVec_t *pv);
long sl_parVecL2G(sl_parVec_t *pv, long local_index);
long sl_parVecG2L(sl_parVec_t *pv, long global_index);
void sl_parVecSetToVal(sl_parVec_t *pv, float val);
void sl_parVecView(sl_parVec_t *pv);
void sl_parVecSetVal(sl_parVec_t *pv, long index, float val);
void sl_parVecSetValLoc(sl_parVec_t *pv, long lindex, float val);
void sl_parVecAdd(sl_parVec_t *pv, sl_parVec_t *pv2);
void sl_parVecScale(sl_parVec_t *pv, float scale);
float sl_parVecDot(sl_parVec_t *pv, sl_parVec_t *pv2);
void sl_parVecReadExtVec(sl_parVec_t *pv, char* spectrum);


//double complex int
struct zs_parVec;
typedef struct zs_parVec zs_parVec_t;
zs_parVec_t *new_zs_ParVec_1(MPI_Comm ncomm, int lbound, int ubound);
zs_parVec_t *new_zs_ParVec_2(parVecMap_t *map);
void zs_parVec_destory(zs_parVec_t *pv);
int zs_parVecGetLowerBound(zs_parVec_t *pv);
int zs_parVecGetUpperBound(zs_parVec_t *pv);
int zs_parVecGetLocSize(zs_parVec_t *pv);
int zs_parVecGetGlobSize(zs_parVec_t *pv);
dcomplex_t zs_parVecGetVal(zs_parVec_t *pv, int index);
dcomplex_t zs_parVecGetValLoc(zs_parVec_t *pv, int lindex);
dcomplex_t *zs_parVecGetArray(zs_parVec_t *pv);
MPI_Comm zs_parVecGetComm(zs_parVec_t *pv);
int zs_parVecL2G(zs_parVec_t *pv, int local_index);
int zs_parVecG2L(zs_parVec_t *pv, int global_index);
void zs_parVecSetToVal(zs_parVec_t *pv, dcomplex_t val);
void zs_parVecView(zs_parVec_t *pv);
void zs_parVecSetVal(zs_parVec_t *pv, int index, dcomplex_t val);
void zs_parVecSetValLoc(zs_parVec_t *pv, int lindex, dcomplex_t val);
void zs_parVecAdd(zs_parVec_t *pv, zs_parVec_t *pv2);
void zs_parVecScale(zs_parVec_t *pv, dcomplex_t scale);
dcomplex_t zs_parVecDot(zs_parVec_t *pv, zs_parVec_t *pv2);
void zs_parVecReadExtVec(zs_parVec_t *pv, char* spectrum);



//interface of selected functions

#ifdef __cplusplus
}
#endif




#endif
