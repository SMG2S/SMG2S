/*
MIT License
Copyright (c) 2019 Xinzhe WU @ Maison de la Simulation, France
Copyright (c) 2019-2022, Xinzhe Wu @ Simulation and Data Laboratory Quantum 
									 Materials,  Forschungszentrum Juelich GmbH.
									 
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

#ifndef __C_SMG2S_H__
#define __C_SMG2S_H__

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


//double complex long
struct zl_parVec;
typedef struct zl_parVec zl_parVec_t;
zl_parVec_t *new_zl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound);
zl_parVec_t *new_zl_ParVec_2(parVecMap_t *map);
void zl_parVec_destory(zl_parVec_t *pv);
long zl_parVecGetLowerBound(zl_parVec_t *pv);
long zl_parVecGetUpperBound(zl_parVec_t *pv);
long zl_parVecGetLocSize(zl_parVec_t *pv);
long zl_parVecGetGlobSize(zl_parVec_t *pv);
dcomplex_t zl_parVecGetVal(zl_parVec_t *pv, long index);
dcomplex_t zl_parVecGetValLoc(zl_parVec_t *pv, long lindex);
dcomplex_t *zl_parVecGetArray(zl_parVec_t *pv);
MPI_Comm zl_parVecGetComm(zl_parVec_t *pv);
long zl_parVecL2G(zl_parVec_t *pv, long local_index);
long zl_parVecG2L(zl_parVec_t *pv, long global_index);
void zl_parVecSetToVal(zl_parVec_t *pv, dcomplex_t val);
void zl_parVecView(zl_parVec_t *pv);
void zl_parVecSetVal(zl_parVec_t *pv, long index, dcomplex_t val);
void zl_parVecSetValLoc(zl_parVec_t *pv, long lindex, dcomplex_t val);
void zl_parVecAdd(zl_parVec_t *pv, zl_parVec_t *pv2);
void zl_parVecScale(zl_parVec_t *pv, dcomplex_t scale);
dcomplex_t zl_parVecDot(zl_parVec_t *pv, zl_parVec_t *pv2);
void zl_parVecReadExtVec(zl_parVec_t *pv, char* spectrum);


//float complex int
struct cs_parVec;
typedef struct cs_parVec cs_parVec_t;
cs_parVec_t *new_cs_ParVec_1(MPI_Comm ncomm, int lbound, int ubound);
cs_parVec_t *new_cs_ParVec_2(parVecMap_t *map);
void cs_parVec_destory(cs_parVec_t *pv);
int cs_parVecGetLowerBound(cs_parVec_t *pv);
int cs_parVecGetUpperBound(cs_parVec_t *pv);
int cs_parVecGetLocSize(cs_parVec_t *pv);
int cs_parVecGetGlobSize(cs_parVec_t *pv);
fcomplex_t cs_parVecGetVal(cs_parVec_t *pv, int index);
fcomplex_t cs_parVecGetValLoc(cs_parVec_t *pv, int lindex);
fcomplex_t *cs_parVecGetArray(cs_parVec_t *pv);
MPI_Comm cs_parVecGetComm(cs_parVec_t *pv);
int cs_parVecL2G(cs_parVec_t *pv, int local_index);
int cs_parVecG2L(cs_parVec_t *pv, int global_index);
void cs_parVecSetToVal(cs_parVec_t *pv, fcomplex_t val);
void cs_parVecView(cs_parVec_t *pv);
void cs_parVecSetVal(cs_parVec_t *pv, int index, fcomplex_t val);
void cs_parVecSetValLoc(cs_parVec_t *pv, int lindex, fcomplex_t val);
void cs_parVecAdd(cs_parVec_t *pv, cs_parVec_t *pv2);
void cs_parVecScale(cs_parVec_t *pv, fcomplex_t scale);
fcomplex_t cs_parVecDot(cs_parVec_t *pv, cs_parVec_t *pv2);
void cs_parVecReadExtVec(cs_parVec_t *pv, char* spectrum);


//float complex long
struct cl_parVec;
typedef struct cl_parVec cl_parVec_t;
cl_parVec_t *new_cl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound);
cl_parVec_t *new_cl_ParVec_2(parVecMap_t *map);
void cl_parVec_destory(cl_parVec_t *pv);
long cl_parVecGetLowerBound(cl_parVec_t *pv);
long cl_parVecGetUpperBound(cl_parVec_t *pv);
long cl_parVecGetLocSize(cl_parVec_t *pv);
long cl_parVecGetGlobSize(cl_parVec_t *pv);
fcomplex_t cl_parVecGetVal(cl_parVec_t *pv, long index);
fcomplex_t cl_parVecGetValLoc(cl_parVec_t *pv, long lindex);
fcomplex_t *cl_parVecGetArray(cl_parVec_t *pv);
MPI_Comm cl_parVecGetComm(cl_parVec_t *pv);
long cl_parVecL2G(cl_parVec_t *pv, long local_index);
long cl_parVecG2L(cl_parVec_t *pv, long global_index);
void cl_parVecSetToVal(cl_parVec_t *pv, fcomplex_t val);
void cl_parVecView(cl_parVec_t *pv);
void cl_parVecSetVal(cl_parVec_t *pv, long index, fcomplex_t val);
void cl_parVecSetValLoc(cl_parVec_t *pv, long lindex, fcomplex_t val);
void cl_parVecAdd(cl_parVec_t *pv, cl_parVec_t *pv2);
void cl_parVecScale(cl_parVec_t *pv, fcomplex_t scale);
fcomplex_t cl_parVecDot(cl_parVec_t *pv, cl_parVec_t *pv2);
void cl_parVecReadExtVec(cl_parVec_t *pv, char* spectrum);

/*parMatrixSparse*/

//double int
struct ds_parMatSparse;
typedef struct ds_parMatSparse ds_parMatSparse_t;
ds_parMatSparse_t *new_ds_ParMatSparse_1(ds_parVec_t *pv);
ds_parMatSparse_t *new_ds_ParMatSparse_2(parVecMap_t *map);
void ds_parMatSparse_destory(ds_parMatSparse_t *pm);
int ds_parMatSparse_GetNRows(ds_parMatSparse_t *pm);
int ds_parMatSparse_GetNCols(ds_parMatSparse_t *pm);
int ds_parMatSparse_GetNNzLoc(ds_parMatSparse_t *pm);
int ds_parMatSparse_GetLowerBound(ds_parMatSparse_t *pm);
int ds_parMatSparse_GetUpperBound(ds_parMatSparse_t *pm);
MPI_Comm ds_parMatSparse_GetComm(ds_parMatSparse_t *pm);
void ds_parMatSparse_SetValLocal(ds_parMatSparse_t *pm, int row, int col, double val);
void ds_parMatSparse_SetVal(ds_parMatSparse_t *pm, int row, int col, double val);
void ds_parMatSparse_AddValLocal(ds_parMatSparse_t *pm, int row, int col, double val);
void ds_parMatSparse_AddVal(ds_parMatSparse_t *pm, int row, int col, double val);
double ds_parMatSparse_GetValLocal(ds_parMatSparse_t *pm, int row, int col);
double ds_parMatSparse_GetVal(ds_parMatSparse_t *pm, int row, int col);
void ds_parMatSparse_toCSRArray(ds_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, double **val);
void ds_parMatSparse_Show(ds_parMatSparse_t *pm);
void ds_parMatSparse_View(ds_parMatSparse_t *pm);
void ds_parMatSparse_Wrt2MM(ds_parMatSparse_t *pm, char *file_name);

//double long
struct dl_parMatSparse;
typedef struct dl_parMatSparse dl_parMatSparse_t;
dl_parMatSparse_t *new_dl_ParMatSparse_1(dl_parVec_t *pv);
dl_parMatSparse_t *new_dl_ParMatSparse_2(parVecMap_t *map);
void dl_parMatSparse_destory(dl_parMatSparse_t *pm);
long dl_parMatSparse_GetNRows(dl_parMatSparse_t *pm);
long dl_parMatSparse_GetNCols(dl_parMatSparse_t *pm);
long dl_parMatSparse_GetNNzLoc(dl_parMatSparse_t *pm);
long dl_parMatSparse_GetLowerBound(dl_parMatSparse_t *pm);
long dl_parMatSparse_GetUpperBound(dl_parMatSparse_t *pm);
MPI_Comm dl_parMatSparse_GetComm(dl_parMatSparse_t *pm);
void dl_parMatSparse_SetValLocal(dl_parMatSparse_t *pm, long row, long col, double val);
void dl_parMatSparse_SetVal(dl_parMatSparse_t *pm, long row, long col, double val);
void dl_parMatSparse_AddValLocal(dl_parMatSparse_t *pm, long row, long col, double val);
void dl_parMatSparse_AddVal(dl_parMatSparse_t *pm, long row, long col, double val);
double dl_parMatSparse_GetValLocal(dl_parMatSparse_t *pm, long row, long col);
double dl_parMatSparse_GetVal(dl_parMatSparse_t *pm, long row, long col);
void dl_parMatSparse_toCSRArray(dl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, double **val);
void dl_parMatSparse_Show(dl_parMatSparse_t *pm);
void dl_parMatSparse_View(dl_parMatSparse_t *pm);
void dl_parMatSparse_Wrt2MM(dl_parMatSparse_t *pm, char *file_name);

//float int
struct ss_parMatSparse;
typedef struct ss_parMatSparse ss_parMatSparse_t;
ss_parMatSparse_t *new_ss_ParMatSparse_1(ss_parVec_t *pv);
ss_parMatSparse_t *new_ss_ParMatSparse_2(parVecMap_t *map);
void ss_parMatSparse_destory(ss_parMatSparse_t *pm);
int ss_parMatSparse_GetNRows(ss_parMatSparse_t *pm);
int ss_parMatSparse_GetNCols(ss_parMatSparse_t *pm);
int ss_parMatSparse_GetNNzLoc(ss_parMatSparse_t *pm);
int ss_parMatSparse_GetLowerBound(ss_parMatSparse_t *pm);
int ss_parMatSparse_GetUpperBound(ss_parMatSparse_t *pm);
MPI_Comm ss_parMatSparse_GetComm(ss_parMatSparse_t *pm);
void ss_parMatSparse_SetValLocal(ss_parMatSparse_t *pm, int row, int col, float val);
void ss_parMatSparse_SetVal(ss_parMatSparse_t *pm, int row, int col, float val);
void ss_parMatSparse_AddValLocal(ss_parMatSparse_t *pm, int row, int col, float val);
void ss_parMatSparse_AddVal(ss_parMatSparse_t *pm, int row, int col, float val);
float ss_parMatSparse_GetValLocal(ss_parMatSparse_t *pm, int row, int col);
float ss_parMatSparse_GetVal(ss_parMatSparse_t *pm, int row, int col);
void ss_parMatSparse_toCSRArray(ss_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, float **val);
void ss_parMatSparse_Show(ss_parMatSparse_t *pm);
void ss_parMatSparse_View(ss_parMatSparse_t *pm);
void ss_parMatSparse_Wrt2MM(ss_parMatSparse_t *pm, char *file_name);


//float long
struct sl_parMatSparse;
typedef struct sl_parMatSparse sl_parMatSparse_t;
sl_parMatSparse_t *new_sl_ParMatSparse_1(sl_parVec_t *pv);
sl_parMatSparse_t *new_sl_ParMatSparse_2(parVecMap_t *map);
void sl_parMatSparse_destory(sl_parMatSparse_t *pm);
long sl_parMatSparse_GetNRows(sl_parMatSparse_t *pm);
long sl_parMatSparse_GetNCols(sl_parMatSparse_t *pm);
long sl_parMatSparse_GetNNzLoc(sl_parMatSparse_t *pm);
long sl_parMatSparse_GetLowerBound(sl_parMatSparse_t *pm);
long sl_parMatSparse_GetUpperBound(sl_parMatSparse_t *pm);
MPI_Comm sl_parMatSparse_GetComm(sl_parMatSparse_t *pm);
void sl_parMatSparse_SetValLocal(sl_parMatSparse_t *pm, long row, long col, float val);
void sl_parMatSparse_SetVal(sl_parMatSparse_t *pm, long row, long col, float val);
void sl_parMatSparse_AddValLocal(sl_parMatSparse_t *pm, long row, long col, float val);
void sl_parMatSparse_AddVal(sl_parMatSparse_t *pm, long row, long col, float val);
float sl_parMatSparse_GetValLocal(sl_parMatSparse_t *pm, long row, long col);
float sl_parMatSparse_GetVal(sl_parMatSparse_t *pm, long row, long col);
void sl_parMatSparse_toCSRArray(sl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, float **val);
void sl_parMatSparse_Show(sl_parMatSparse_t *pm);
void sl_parMatSparse_View(sl_parMatSparse_t *pm);
void sl_parMatSparse_Wrt2MM(sl_parMatSparse_t *pm, char *file_name);

//dcomplex_t int
struct zs_parMatSparse;
typedef struct zs_parMatSparse zs_parMatSparse_t;
zs_parMatSparse_t *new_zs_ParMatSparse_1(zs_parVec_t *pv);
zs_parMatSparse_t *new_zs_ParMatSparse_2(parVecMap_t *map);
void zs_parMatSparse_destory(zs_parMatSparse_t *pm);
int zs_parMatSparse_GetNRows(zs_parMatSparse_t *pm);
int zs_parMatSparse_GetNCols(zs_parMatSparse_t *pm);
int zs_parMatSparse_GetNNzLoc(zs_parMatSparse_t *pm);
int zs_parMatSparse_GetLowerBound(zs_parMatSparse_t *pm);
int zs_parMatSparse_GetUpperBound(zs_parMatSparse_t *pm);
MPI_Comm zs_parMatSparse_GetComm(zs_parMatSparse_t *pm);
void zs_parMatSparse_SetValLocal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val);
void zs_parMatSparse_SetVal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val);
void zs_parMatSparse_AddValLocal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val);
void zs_parMatSparse_AddVal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val);
dcomplex_t zs_parMatSparse_GetValLocal(zs_parMatSparse_t *pm, int row, int col);
dcomplex_t zs_parMatSparse_GetVal(zs_parMatSparse_t *pm, int row, int col);
void zs_parMatSparse_toCSRArray(zs_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, dcomplex_t **val);
void zs_parMatSparse_Show(zs_parMatSparse_t *pm);
void zs_parMatSparse_View(zs_parMatSparse_t *pm);
void zs_parMatSparse_Wrt2MMCmplx(zs_parMatSparse_t *pm, char *file_name);

//dcomplex_t long
struct zl_parMatSparse;
typedef struct zl_parMatSparse zl_parMatSparse_t;
zl_parMatSparse_t *new_zl_ParMatSparse_1(zl_parVec_t *pv);
zl_parMatSparse_t *new_zl_ParMatSparse_2(parVecMap_t *map);
void zl_parMatSparse_destory(zl_parMatSparse_t *pm);
long zl_parMatSparse_GetNRows(zl_parMatSparse_t *pm);
long zl_parMatSparse_GetNCols(zl_parMatSparse_t *pm);
long zl_parMatSparse_GetNNzLoc(zl_parMatSparse_t *pm);
long zl_parMatSparse_GetLowerBound(zl_parMatSparse_t *pm);
long zl_parMatSparse_GetUpperBound(zl_parMatSparse_t *pm);
MPI_Comm zl_parMatSparse_GetComm(zl_parMatSparse_t *pm);
void zl_parMatSparse_SetValLocal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val);
void zl_parMatSparse_SetVal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val);
void zl_parMatSparse_AddValLocal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val);
void zl_parMatSparse_AddVal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val);
dcomplex_t zl_parMatSparse_GetValLocal(zl_parMatSparse_t *pm, long row, long col);
dcomplex_t zl_parMatSparse_GetVal(zl_parMatSparse_t *pm, long row, long col);
void zl_parMatSparse_toCSRArray(zl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, dcomplex_t **val);
void zl_parMatSparse_Show(zl_parMatSparse_t *pm);
void zl_parMatSparse_View(zl_parMatSparse_t *pm);
void zl_parMatSparse_Wrt2MMCmplx(zl_parMatSparse_t *pm, char *file_name);

//fcomplex_t int
struct cs_parMatSparse;
typedef struct cs_parMatSparse cs_parMatSparse_t;
cs_parMatSparse_t *new_cs_ParMatSparse_1(cs_parVec_t *pv);
cs_parMatSparse_t *new_cs_ParMatSparse_2(parVecMap_t *map);
void cs_parMatSparse_destory(cs_parMatSparse_t *pm);
int cs_parMatSparse_GetNRows(cs_parMatSparse_t *pm);
int cs_parMatSparse_GetNCols(cs_parMatSparse_t *pm);
int cs_parMatSparse_GetNNzLoc(cs_parMatSparse_t *pm);
int cs_parMatSparse_GetLowerBound(cs_parMatSparse_t *pm);
int cs_parMatSparse_GetUpperBound(cs_parMatSparse_t *pm);
MPI_Comm cs_parMatSparse_GetComm(cs_parMatSparse_t *pm);
void cs_parMatSparse_SetValLocal(cs_parMatSparse_t *pm, int row, int col, fcomplex_t val);
void cs_parMatSparse_SetVal(cs_parMatSparse_t *pm, int row, int col, fcomplex_t val);
void cs_parMatSparse_AddValLocal(cs_parMatSparse_t *pm, int row, int col, fcomplex_t val);
void cs_parMatSparse_AddVal(cs_parMatSparse_t *pm, int row, int col, fcomplex_t val);
fcomplex_t cs_parMatSparse_GetValLocal(cs_parMatSparse_t *pm, int row, int col);
fcomplex_t cs_parMatSparse_GetVal(cs_parMatSparse_t *pm, int row, int col);
void cs_parMatSparse_toCSRArray(cs_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, fcomplex_t **val);
void cs_parMatSparse_Show(cs_parMatSparse_t *pm);
void cs_parMatSparse_View(cs_parMatSparse_t *pm);
void cs_parMatSparse_Wrt2MMCmplx(cs_parMatSparse_t *pm, char *file_name);


//fcomplex_t long
struct cl_parMatSparse;
typedef struct cl_parMatSparse cl_parMatSparse_t;
cl_parMatSparse_t *new_cl_ParMatSparse_1(cl_parVec_t *pv);
cl_parMatSparse_t *new_cl_ParMatSparse_2(parVecMap_t *map);
void cl_parMatSparse_destory(cl_parMatSparse_t *pm);
long cl_parMatSparse_GetNRows(cl_parMatSparse_t *pm);
long cl_parMatSparse_GetNCols(cl_parMatSparse_t *pm);
long cl_parMatSparse_GetNNzLoc(cl_parMatSparse_t *pm);
long cl_parMatSparse_GetLowerBound(cl_parMatSparse_t *pm);
long cl_parMatSparse_GetUpperBound(cl_parMatSparse_t *pm);
MPI_Comm cl_parMatSparse_GetComm(cl_parMatSparse_t *pm);
void cl_parMatSparse_SetValLocal(cl_parMatSparse_t *pm, long row, long col, fcomplex_t val);
void cl_parMatSparse_SetVal(cl_parMatSparse_t *pm, long row, long col, fcomplex_t val);
void cl_parMatSparse_AddValLocal(cl_parMatSparse_t *pm, long row, long col, fcomplex_t val);
void cl_parMatSparse_AddVal(cl_parMatSparse_t *pm, long row, long col, fcomplex_t val);
fcomplex_t cl_parMatSparse_GetValLocal(cl_parMatSparse_t *pm, long row, long col);
fcomplex_t cl_parMatSparse_GetVal(cl_parMatSparse_t *pm, long row, long col);
void cl_parMatSparse_toCSRArray(cl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, fcomplex_t **val);
void cl_parMatSparse_Show(cl_parMatSparse_t *pm);
void cl_parMatSparse_View(cl_parMatSparse_t *pm);
void cl_parMatSparse_Wrt2MMCmplx(cl_parMatSparse_t *pm, char *file_name);

//interface of selected functions
cs_parVec_t *cs_specNonHermLoad(parVecMap_t *map, char *spectrum);
cl_parVec_t *cl_specNonHermLoad(parVecMap_t *map, char *spectrum);
zs_parVec_t *zs_specNonHermLoad(parVecMap_t *map, char *spectrum);
zl_parVec_t *zl_specNonHermLoad(parVecMap_t *map, char *spectrum);

cs_parVec_t *cs_specNonSymmLoad(parVecMap_t *map, char *spectrum);
cl_parVec_t *cl_specNonSymmLoad(parVecMap_t *map, char *spectrum);
zs_parVec_t *zs_specNonSymmLoad(parVecMap_t *map, char *spectrum);
zl_parVec_t *zl_specNonSymmLoad(parVecMap_t *map, char *spectrum);

ds_parVec_t *ds_specNonSymmLoad(parVecMap_t *map, char *spectrum);
dl_parVec_t *dl_specNonSymmLoad(parVecMap_t *map, char *spectrum);
ss_parVec_t *ss_specNonSymmLoad(parVecMap_t *map, char *spectrum);
sl_parVec_t *sl_specNonSymmLoad(parVecMap_t *map, char *spectrum);

cs_parMatSparse_t *cs_nonherm(int probSize, nilp_t *nilp, initMatrix_t *init, char *spectrum, MPI_Comm comm);
cl_parMatSparse_t *cl_nonherm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm);
zs_parMatSparse_t *zs_nonherm(int probSize, nilp_t *nilp, initMatrix_t *init,  char *spectrum, MPI_Comm comm);
zl_parMatSparse_t *zl_nonherm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm);

cs_parMatSparse_t *cs_nonherm_2(int probSize, nilp_t *nilp, initMatrix_t *init, cs_parVec_t *spec);
cl_parMatSparse_t *cl_nonherm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, cl_parVec_t *spec);
zs_parMatSparse_t *zs_nonherm_2(int probSize, nilp_t *nilp, initMatrix_t *init, zs_parVec_t *spec);
zl_parMatSparse_t *zl_nonherm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, zl_parVec_t *spec);

void cs_nonherm_3(int probSize, nilp_t *nilp, cs_parMatSparse_t *Am, cs_parVec_t *spec);
void cl_nonherm_3(long probSize, nilpL_t *nilp, cl_parMatSparse_t *Am, cl_parVec_t *spec);
void zs_nonherm_3(int probSize, nilp_t *nilp, zs_parMatSparse_t *Am, zs_parVec_t *spec);
void zl_nonherm_3(long probSize, nilpL_t *nilp, zl_parMatSparse_t *Am, zl_parVec_t *spec);

//
ss_parMatSparse_t *ss_nonsymm(int probSize, nilp_t *nilp, initMatrix_t *init, char *spectrum, MPI_Comm comm);
sl_parMatSparse_t *sl_nonsymm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm);
ds_parMatSparse_t *ds_nonsymm(int probSize, nilp_t *nilp, initMatrix_t *init,  char *spectrum, MPI_Comm comm);
dl_parMatSparse_t *dl_nonsymm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm);

//
ss_parMatSparse_t *ss_nonsymm_2(int probSize, nilp_t *nilp, initMatrix_t *init, ss_parVec_t *spec);
sl_parMatSparse_t *sl_nonsymm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, sl_parVec_t *spec);
ds_parMatSparse_t *ds_nonsymm_2(int probSize, nilp_t *nilp, initMatrix_t *init, ds_parVec_t *spec);
dl_parMatSparse_t *dl_nonsymm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, dl_parVec_t *spec);

void ss_nonsymm_3(int probSize, nilp_t *nilp, ss_parMatSparse_t *Am, ss_parVec_t *spec);
void sl_nonsymm_3(long probSize, nilpL_t *nilp, sl_parMatSparse_t *Am, sl_parVec_t *spec);
void ds_nonsymm_3(int probSize, nilp_t *nilp, ds_parMatSparse_t *Am, ds_parVec_t *spec);
void dl_nonsymm_3(long probSize, nilpL_t *nilp, dl_parMatSparse_t *Am, dl_parVec_t *spec);

//
ss_parMatSparse_t *ss_nonsymmconj_2(int probSize, nilp_t *nilp, initMatrix_t *init, cs_parVec_t *spec);
sl_parMatSparse_t *sl_nonsymmconj_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, cl_parVec_t *spec);
ds_parMatSparse_t *ds_nonsymmconj_2(int probSize, nilp_t *nilp, initMatrix_t *init, zs_parVec_t *spec);
dl_parMatSparse_t *dl_nonsymmconj_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, zl_parVec_t *spec);

void ss_nonsymmconj_3(int probSize, nilp_t *nilp, ss_parMatSparse_t *Am, cs_parVec_t *spec);
void sl_nonsymmconj_3(long probSize, nilpL_t *nilp, sl_parMatSparse_t *Am, cl_parVec_t *spec);
void ds_nonsymmconj_3(int probSize, nilp_t *nilp, ds_parMatSparse_t *Am, zs_parVec_t *spec);
void dl_nonsymmconj_3(long probSize, nilpL_t *nilp, dl_parMatSparse_t *Am, zl_parVec_t *spec);



#ifdef __cplusplus
}
#endif




#endif