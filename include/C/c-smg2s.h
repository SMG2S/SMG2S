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

/** @defgroup group5 C interface
 *  This module provides the C interface of SMG2S.
 *  - For the scalars, only `double precision` are supported.
 *  - For the integer, both `int` and `long` are supported.
 *  - Namings for C-intefaces:
 * 		- for Nilpotent, initMat and parVectorMap, the interfaces for the case with `long` integer has a letter `L` in some part of related function names.
 * 		- for parVector and parMatrixSparse and other independant functions which don't belong to any classes and structs, their function names starts with 
 * 			one of the four prefix `ds_`, `dl_`, `zs_` and `zl_`. Here, `d` refers to `double`, `s` refers to `int`, `z` refers to `dcomplex_t` and `l` 
 * 			refers to `long`.
 * 		- for the meaning of each function, please refer to the coresponding C++ functions for more details.
 * 		- for the usage, please check the provided examples. 
 *  @{
 */
//!  @brief A struct determines a complex number with double precision.
/*!
	This struct is only userful for the C interface of SMG2S. 
*/
struct dcomplex{
	/*! The real part of a complex scalar */
	double real;
	/*! The imaginary part of a complex scalar */	
	double imag;
};

typedef struct dcomplex dcomplex_t;

//interface of classes and structs
/*Nilpotency Matrix C Wrapper*/
// int
//! A struct which provides a C-interface for Nilpotent<S> with `S=int`
struct nilp;
//! Creates a type name for nilp which is a C-interface for Nilpotent with `S=int` 
typedef struct nilp nilp_t;

//! C-interface of Nilpotent#Nilpotent(S nbOne, S size)
nilp_t *newNilp_1(int nbOne, int size);
//! C-interface of Nilpotent#Nilpotent(S nbOne, S diag, S size)
nilp_t *newNilp_2(int nbOne, int diag, int size);
//! C-interface of Nilpotent#Nilpotent(std::vector<S> nilpvec, S size)
nilp_t *newNilp_3(int *nilpvec, int size);
//! C-interface of Nilpotent#Nilpotent(std::vector<S> nilpvec, S diag, S size)
nilp_t *newNilp_4(int *nilpvec, int diag, int size);
//! Destory a nilp struct 
void nilp_destory(nilp_t *nilp);
//! C-interface of Nilpotent#getDegree 
int nilp_getDegree(nilp_t *nilp);
//! C-interface of Nilpotent#getIndOfZeros 
int* nilp_getIndOfZeros(nilp_t *nilp);
//! C-interface of Nilpotent#show 
void nilp_show(nilp_t *nilp);


// long int
struct nilpL;
//! Creates a type name for nilp which is a C-interface for Nilpotent with `S=long` 
typedef struct nilpL nilpL_t;
//! C-interface of Nilpotent#Nilpotent(S nbOne, S size)
nilpL_t *newNilpL_1(long nbOne, long size);
//! C-interface of Nilpotent#Nilpotent(S nbOne, S diag, S size)
nilpL_t *newNilpL_2(long nbOne, long diag, long size);
//! C-interface of Nilpotent#Nilpotent(std::vector<S> nilpvec, S size)
nilpL_t *newNilpL_3(long *nilpvec, long size);
//! C-interface of Nilpotent#Nilpotent(std::vector<S> nilpvec, S diag, S size)
nilpL_t *newNilpL_4(long *nilpvec, long diag, long size);
//! Destory a nilp struct 
void nilpL_destory(nilpL_t *nilp);
//! C-interface of Nilpotent#getDegree 
long nilpL_getDegree(nilpL_t *nilp);
//! C-interface of Nilpotent#getIndOfZeros 
long* nilpL_getIndOfZeros(nilpL_t *nilp);
//! C-interface of Nilpotent#show 
void nilpL_show(nilpL_t *nilp);


/*init matrix struct*/
//int
struct initMatrix;
//! Creates a type name for initMatrix which is a C-interface for initMat with `S=int` 
typedef struct initMatrix initMatrix_t;
//! C-interface of initMat#initMat()
initMatrix_t *newInitMatrix_1();
//! C-interface of initMat#initMat(S diagl, S diagu)
initMatrix_t *newInitMatrix_2(int diagl, int diagu);
//! C-interface of initMat#initMat(S diagl, S diagu, double Sparsity)
initMatrix_t *newInitMatrix_3(int diagl, int diagu, double Sparsity);
//! C-interface of initMat#initMat(S diagl, S diagu, double Scale, double Sparsity)
initMatrix_t *newInitMatrix_4(int diagl, int diagu, double Scale, double Sparsity);
//! C-interface of initMat#show()
void initMatrix_show(initMatrix_t *init);
//! Destory a initMatrix struct 
void initMatrix_destory(initMatrix_t *init);

//long
struct initMatrixL;
//! Creates a type name for initMatrix which is a C-interface for initMat with `S=long` 
typedef struct initMatrixL initMatrixL_t;
//! C-interface of initMat#initMat()
initMatrixL_t *newInitMatrixL_1();
//! C-interface of initMat#initMat(S diagl, S diagu)
initMatrixL_t *newInitMatrixL_2(long diagl, long diagu);
//! C-interface of initMat#initMat(S diagl, S diagu, double Sparsity)
initMatrixL_t *newInitMatrixL_3(long diagl, long diagu, double Sparsity);
//! C-interface of initMat#initMat(S diagl, S diagu, double Scale, double Sparsity)
initMatrixL_t *newInitMatrixL_4(long diagl, long diagu, double Scale, double Sparsity);
//! C-interface of initMat#show()
void initMatrixL_show(initMatrixL_t *init);
//! Destory a initMatrixL struct 
void initMatrixL_destory(initMatrixL_t *init);


/*parVectorMap*/
struct parVecMap;
//! Creates a type name for parVecMap which is a C-interface for parVectorMap with `S=int` 
typedef struct parVecMap parVecMap_t;
//! C-interfance of parVectorMap#parVectorMap()
parVecMap_t *newParVecMap_empty();
//! C-interfance of parVectorMap#parVectorMap(MPI_Comm ncomm, S lbound, S ubound)
parVecMap_t *newParVecMap(MPI_Comm ncomm, int lbound, int ubound);
//! C-interfance of parVectorMap#GetCurrentComm()
MPI_Comm parVecMapGetComm(parVecMap_t *pv);
//! C-interfance of parVectorMap#Loc2Glob(S local_index)
int parVecMapL2G(parVecMap_t *pv, int local_index);
//! C-interfance of parVectorMap#Glob2Loc(S global_index)
int parVecMapG2L(parVecMap_t *pv, int global_index);
//! C-interfance of parVectorMap#GetLocalSize()
int parVecMapGetLocSize(parVecMap_t *pv);
//! C-interfance of parVectorMap#GetGlobalSize()
int parVecMapGetGlobSize(parVecMap_t *pv);
//! Destory a parVecMap struct 
void parVecMap_destory(parVecMap_t *pv);

struct parVecMapL;
//! Creates a type name for parVecMap which is a C-interface for parVectorMap with `S=long` typedef struct parVecMapL parVecMapL_t;
typedef struct parVecMapL parVecMapL_t;
//! C-interfance of parVectorMap#parVectorMap()
parVecMapL_t *newParVecMapL_empty();
//! C-interfance of parVectorMap#parVectorMap(MPI_Comm ncomm, S lbound, S ubound)
parVecMapL_t *newParVecMapL(MPI_Comm ncomm, long lbound, long ubound);
//! C-interfance of parVectorMap#GetCurrentComm()
MPI_Comm parVecMapLGetComm(parVecMapL_t *pv);
//! C-interfance of parVectorMap#Loc2Glob(S local_index)
long parVecMapLL2G(parVecMapL_t *pv, long local_index);
//! C-interfance of parVectorMap#Glob2Loc(S global_index)
long parVecMapLG2L(parVecMapL_t *pv, long global_index);
//! C-interfance of parVectorMap#GetLocalSize()
long parVecMapLGetLocSize(parVecMapL_t *pv);
//! C-interfance of parVectorMap#GetGlobalSize()
long parVecMapLGetGlobSize(parVecMapL_t *pv);
//! Destory a parVecMapL struct 
void parVecMapL_destory(parVecMapL_t *pv);


/*parVector*/
//double int
struct ds_parVec;
//! Creates a type name for ds_parVec which is a C-interface for parVector with `T=double` and `S=int` 
typedef struct ds_parVec ds_parVec_t;
//! C-interface of parVector#parVector(MPI_Comm ncomm, S lbound, S ubound)
ds_parVec_t *new_ds_ParVec_1(MPI_Comm ncomm, int lbound, int ubound);
//! C-interface of parVector#parVector(parVectorMap<S> map)
ds_parVec_t *new_ds_ParVec_2(parVecMap_t *map);
//! Destory the struct
void ds_parVec_destory(ds_parVec_t *pv);
//! C-interface of parVector#GetLowerBound()
int ds_parVecGetLowerBound(ds_parVec_t *pv);
//! C-interface of parVector#GetUpperBound()
int ds_parVecGetUpperBound(ds_parVec_t *pv); 
//! C-interface of parVector#GetLocalSize()
int ds_parVecGetLocSize(ds_parVec_t *pv);
//! C-interface of parVector#GetGlobalSize()
int ds_parVecGetGlobSize(ds_parVec_t *pv);
//! C-interface of parVector#GetValue(S index)
double ds_parVecGetVal(ds_parVec_t *pv, int index);
//! C-interface of parVector#GetValueLocal(S lindex)
double ds_parVecGetValLoc(ds_parVec_t *pv, int lindex);
//! C-interface of parVector#GetArray()
double *ds_parVecGetArray(ds_parVec_t *pv);
//! C-interface of parVector#GetVecMap()
parVecMap_t *ds_parVecGetMap(ds_parVec_t *pv);
//! C-interface of parVector#GetComm()
MPI_Comm ds_parVecGetComm(ds_parVec_t *pv);
//! C-interface of parVector#Loc2Glob(S local_index)
int ds_parVecL2G(ds_parVec_t *pv, int local_index);
//! C-interface of parVector#Glob2Loc(S global_index)
int ds_parVecG2L(ds_parVec_t *pv, int global_index);
//! C-interface of parVector#SetToValue(T value)
void ds_parVecSetToVal(ds_parVec_t *pv, double val);
//! C-interface of parVector#VecView()
void ds_parVecView(ds_parVec_t *pv);
//! C-interface of parVector#SetValueGlobal(S index, T value)
void ds_parVecSetVal(ds_parVec_t *pv, int index, double val);
//! C-interface of parVector#SetValueLocal(S row, T value)
void ds_parVecSetValLoc(ds_parVec_t *pv, int lindex, double val);
//! C-interface of parVector#VecAdd(parVector v)
void ds_parVecAdd(ds_parVec_t *pv, ds_parVec_t *pv2);
//! C-interface of parVector#VecScale(T scale)
void ds_parVecScale(ds_parVec_t *pv, double scale);
//! C-interface of parVector#VecDot(parVector v)
double ds_parVecDot(ds_parVec_t *pv, ds_parVec_t *pv2);
//! C-interface of parVector#ReadExtVec(std::string spectrum)
void ds_parVecReadExtVec(ds_parVec_t *pv, char* spectrum);
//! C-interface of parVector#writeToTxt (std::string file_name)
void ds_writeToTxt(ds_parVec_t *pv, char* spectrum);

//double long
struct dl_parVec;
//! Creates a type name for dl_parVec which is a C-interface for parVector with `T=double` and `S=long` 
typedef struct dl_parVec dl_parVec_t;
//! C-interface of parVector#parVector(MPI_Comm ncomm, S lbound, S ubound)
dl_parVec_t *new_dl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound);
//! C-interface of parVector#parVector(parVectorMap<S> map)
dl_parVec_t *new_dl_ParVec_2(parVecMapL_t *map);
//! Destory the struct
void dl_parVec_destory(dl_parVec_t *pv);
//! C-interface of parVector#GetLowerBound()
long dl_parVecGetLowerBound(dl_parVec_t *pv);
//! C-interface of parVector#GetUpperBound()
long dl_parVecGetUpperBound(dl_parVec_t *pv) ;
//! C-interface of parVector#GetLocalSize()
long dl_parVecGetLocSize(dl_parVec_t *pv);
//! C-interface of parVector#GetGlobalSize()
long dl_parVecGetGlobSize(dl_parVec_t *pv);
//! C-interface of parVector#GetValue(S index)
double dl_parVecGetVal(dl_parVec_t *pv, long index);
//! C-interface of parVector#GetValueLocal(S lindex)
double dl_parVecGetValLoc(dl_parVec_t *pv, long lindex);
//! C-interface of parVector#GetArray()
double *dl_parVecGetArray(dl_parVec_t *pv);
//! C-interface of parVector#GetVecMap()
parVecMapL_t *dl_parVecGetMap(dl_parVec_t *pv);
//! C-interface of parVector#GetComm()
MPI_Comm dl_parVecGetComm(dl_parVec_t *pv);
//! C-interface of parVector#Loc2Glob(S local_index)
long dl_parVecL2G(dl_parVec_t *pv, long local_index);
//! C-interface of parVector#Glob2Loc(S global_index)
long dl_parVecG2L(dl_parVec_t *pv, long global_index);
//! C-interface of parVector#SetToValue(T value)
void dl_parVecSetToVal(dl_parVec_t *pv, double val);
//! C-interface of parVector#VecView()
void dl_parVecView(dl_parVec_t *pv);
//! C-interface of parVector#SetValueGlobal(S index, T value)
void dl_parVecSetVal(dl_parVec_t *pv, long index, double val);
//! C-interface of parVector#SetValueLocal(S row, T value)
void dl_parVecSetValLoc(dl_parVec_t *pv, long lindex, double val);
//! C-interface of parVector#VecAdd(parVector v)
void dl_parVecAdd(dl_parVec_t *pv, dl_parVec_t *pv2);
//! C-interface of parVector#VecScale(T scale)
void dl_parVecScale(dl_parVec_t *pv, double scale);
//! C-interface of parVector#VecDot(parVector v)
double dl_parVecDot(dl_parVec_t *pv, dl_parVec_t *pv2);
//! C-interface of parVector#ReadExtVec(std::string spectrum)
void dl_parVecReadExtVec(dl_parVec_t *pv, char* spectrum);
//! C-interface of parVector#writeToTxt (std::string file_name)
void dl_writeToTxt(dl_parVec_t *pv, char* spectrum);


//double complex int
struct zs_parVec;
//! Creates a type name for zs_parVec which is a C-interface for parVector with `T=std::complex<double>` and `S=int` 
typedef struct zs_parVec zs_parVec_t;
//! C-interface of parVector#parVector(MPI_Comm ncomm, S lbound, S ubound)
zs_parVec_t *new_zs_ParVec_1(MPI_Comm ncomm, int lbound, int ubound);
//! C-interface of parVector#parVector(parVectorMap<S> map)
zs_parVec_t *new_zs_ParVec_2(parVecMap_t *map);
//! Destory the struct
void zs_parVec_destory(zs_parVec_t *pv);
//! C-interface of parVector#GetLowerBound()
int zs_parVecGetLowerBound(zs_parVec_t *pv);
//! C-interface of parVector#GetUpperBound()
int zs_parVecGetUpperBound(zs_parVec_t *pv); 
//! C-interface of parVector#GetLocalSize()
int zs_parVecGetLocSize(zs_parVec_t *pv);
//! C-interface of parVector#GetGlobalSize()
int zs_parVecGetGlobSize(zs_parVec_t *pv);
//! C-interface of parVector#GetValue(S index)
dcomplex_t zs_parVecGetVal(zs_parVec_t *pv, int index);
//! C-interface of parVector#GetValueLocal(S lindex)
dcomplex_t zs_parVecGetValLoc(zs_parVec_t *pv, int lindex);
//! C-interface of parVector#GetArray()
dcomplex_t *zs_parVecGetArray(zs_parVec_t *pv);
//! C-interface of parVector#GetVecMap()
parVecMap_t *zs_parVecGetMap(zs_parVec_t *pv);
//! C-interface of parVector#GetComm()
MPI_Comm zs_parVecGetComm(zs_parVec_t *pv);
//! C-interface of parVector#Loc2Glob(S local_index)
int zs_parVecL2G(zs_parVec_t *pv, int local_index);
//! C-interface of parVector#Glob2Loc(S global_index)
int zs_parVecG2L(zs_parVec_t *pv, int global_index);
//! C-interface of parVector#SetToValue(T value)
void zs_parVecSetToVal(zs_parVec_t *pv, dcomplex_t val);
//! C-interface of parVector#VecView()
void zs_parVecView(zs_parVec_t *pv);
//! C-interface of parVector#SetValueGlobal(S index, T value)
void zs_parVecSetVal(zs_parVec_t *pv, int index, dcomplex_t val);
//! C-interface of parVector#SetValueLocal(S row, T value)
void zs_parVecSetValLoc(zs_parVec_t *pv, int lindex, dcomplex_t val);
//! C-interface of parVector#VecAdd(parVector v)
void zs_parVecAdd(zs_parVec_t *pv, zs_parVec_t *pv2);
//! C-interface of parVector#VecScale(T scale)
void zs_parVecScale(zs_parVec_t *pv, dcomplex_t scale);
//! C-interface of parVector#VecDot(parVector v)
dcomplex_t zs_parVecDot(zs_parVec_t *pv, zs_parVec_t *pv2);
//! C-interface of parVector#ReadExtVec(std::string spectrum)
void zs_parVecReadExtVec(zs_parVec_t *pv, char* spectrum);
//! C-interface of parVector#writeToCmplx(std::string file_name)
void zs_writeToTxtCmplx(zs_parVec_t *pv, char* spectrum);


//double complex long
struct zl_parVec;
//! Creates a type name for zl_parVec which is a C-interface for parVector with `T=std::complex<double>` and `S=long` 
typedef struct zl_parVec zl_parVec_t;
//! C-interface of parVector#parVector(MPI_Comm ncomm, S lbound, S ubound)
zl_parVec_t *new_zl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound);
//! C-interface of parVector#parVector(parVectorMap<S> map)
zl_parVec_t *new_zl_ParVec_2(parVecMapL_t *map);
//! Destory the struct
void zl_parVec_destory(zl_parVec_t *pv);
//! C-interface of parVector#GetLowerBound()
long zl_parVecGetLowerBound(zl_parVec_t *pv);
//! C-interface of parVector#GetUpperBound()
long zl_parVecGetUpperBound(zl_parVec_t *pv) ;
//! C-interface of parVector#GetLocalSize()
long zl_parVecGetLocSize(zl_parVec_t *pv);
//! C-interface of parVector#GetGlobalSize()
long zl_parVecGetGlobSize(zl_parVec_t *pv);
//! C-interface of parVector#GetValue(S index)
dcomplex_t zl_parVecGetVal(zl_parVec_t *pv, long index);
//! C-interface of parVector#GetValueLocal(S lindex)
dcomplex_t zl_parVecGetValLoc(zl_parVec_t *pv, long lindex);
//! C-interface of parVector#GetArray()
dcomplex_t *zl_parVecGetArray(zl_parVec_t *pv);
//! C-interface of parVector#GetVecMap()
parVecMapL_t *zl_parVecGetMap(zl_parVec_t *pv);
//! C-interface of parVector#GetComm()
MPI_Comm zl_parVecGetComm(zl_parVec_t *pv);
//! C-interface of parVector#Loc2Glob(S local_index)
long zl_parVecL2G(zl_parVec_t *pv, long local_index);
//! C-interface of parVector#Glob2Loc(S global_index)
long zl_parVecG2L(zl_parVec_t *pv, long global_index);
//! C-interface of parVector#SetToValue(T value)
void zl_parVecSetToVal(zl_parVec_t *pv, dcomplex_t val);
//! C-interface of parVector#VecView()
void zl_parVecView(zl_parVec_t *pv);
//! C-interface of parVector#SetValueGlobal(S index, T value)
void zl_parVecSetVal(zl_parVec_t *pv, long index, dcomplex_t val);
//! C-interface of parVector#SetValueLocal(S row, T value)
void zl_parVecSetValLoc(zl_parVec_t *pv, long lindex, dcomplex_t val);
//! C-interface of parVector#VecAdd(parVector v)
void zl_parVecAdd(zl_parVec_t *pv, zl_parVec_t *pv2);
//! C-interface of parVector#VecScale(T scale)
void zl_parVecScale(zl_parVec_t *pv, dcomplex_t scale);

//! C-interface of parVector#VecDot(parVector v)
dcomplex_t zl_parVecDot(zl_parVec_t *pv, zl_parVec_t *pv2);
//! C-interface of parVector#ReadExtVec(std::string spectrum)
void zl_parVecReadExtVec(zl_parVec_t *pv, char* spectrum);
//! C-interface of parVector#writeToCmplx(std::string file_name)
void zl_writeToTxtCmplx(zl_parVec_t *pv, char* spectrum);




/*parMatrixSparse*/

//double int
struct ds_parMatSparse;
//! Creates a type name for ds_parMatSparse which is a C-interface for parMatrixSparse with `T=double` and `S=int` 
typedef struct ds_parMatSparse ds_parMatSparse_t;
//! C-interface of parMatrixSparse#parMatrixSparse(parVector< T, S > vec)
ds_parMatSparse_t *new_ds_ParMatSparse_1(ds_parVec_t *pv);
//! C-interface of parMatrixSparse#parMatrixSparse(parVectorMap< S > map)
ds_parMatSparse_t *new_ds_ParMatSparse_2(parVecMap_t *map);
//! Destory the struct
void ds_parMatSparse_destory(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNRows()
int ds_parMatSparse_GetNRows(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNCols()
int ds_parMatSparse_GetNCols(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNnzLoc()
int ds_parMatSparse_GetNNzLoc(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetLowerBound()
int ds_parMatSparse_GetLowerBound(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetUpperBound()
int ds_parMatSparse_GetUpperBound(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetComm()
MPI_Comm ds_parMatSparse_GetComm(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#SetValueLocal(S row, S col, T value)
void ds_parMatSparse_SetValLocal(ds_parMatSparse_t *pm, int row, int col, double val);
//! C-interface of parMatrixSparse#SetValue (S row, S col, T value)
void ds_parMatSparse_SetVal(ds_parMatSparse_t *pm, int row, int col, double val);
//! C-interface of parMatrixSparse#AddValueLocal(S row, S col, T value)
void ds_parMatSparse_AddValLocal(ds_parMatSparse_t *pm, int row, int col, double val);
//! C-interface of parMatrixSparse#AddValue(S row, S col, T value)
void ds_parMatSparse_AddVal(ds_parMatSparse_t *pm, int row, int col, double val);
//! C-interface of parMatrixSparse#GetValueLocal (S row, S col)
double ds_parMatSparse_GetValLocal(ds_parMatSparse_t *pm, int row, int col);
//! C-interface of parMatrixSparse#GetValue (S row, S col)
double ds_parMatSparse_GetVal(ds_parMatSparse_t *pm, int row, int col);
//! convert a parMatrixSparse matrix into CSR format
/*!
  * @param[in] pm the struct of parMatrixSparse to be converted
  * @param[out] nrows number of rows of the converted matrix in CSR format on each MPI proc
  * @param[out] nnz number of nnz of the converted matrix in CSR format on each MPI proc
  * @param[out] roffs an array containing the offsets of each row of CSR format
  * @param[out] cols an array containing the column index of each nnz of CSR format
  * @param[out] val an array containing all the nnz of CSR format
*/	
void ds_parMatSparse_toCSRArray(ds_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, double **val);
//! C-interface of parMatrixSparse#show()
void ds_parMatSparse_Show(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#MatView()
void ds_parMatSparse_View(ds_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#writeToMatrixMarket(std::string file_name)
void ds_parMatSparse_Wrt2MM(ds_parMatSparse_t *pm, char *file_name);
//! C-interface of parMatrixSparse#initMat(S diag_l, S diag_u, Base< T > scale, T shift, Base< T > sparsity)
void ds_parMatSparse_initMat(ds_parMatSparse_t *pm, int diagl, int diagu, double Scale, double Sparsity);


//double long
struct dl_parMatSparse;
//! Creates a type name for dl_parMatSparse which is a C-interface for parMatrixSparse with `T=double` and `S=long` 
typedef struct dl_parMatSparse dl_parMatSparse_t;
//! C-interface of parMatrixSparse#parMatrixSparse(parVector< T, S > vec)
dl_parMatSparse_t *new_dl_ParMatSparse_1(dl_parVec_t *pv);
//! C-interface of parMatrixSparse#parMatrixSparse(parVectorMap< S > map)
dl_parMatSparse_t *new_dl_ParMatSparse_2(parVecMapL_t *map);
//! Destory the struct
void dl_parMatSparse_destory(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNRows()
long dl_parMatSparse_GetNRows(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNCols()
long dl_parMatSparse_GetNCols(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNnzLoc()
long dl_parMatSparse_GetNNzLoc(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetLowerBound()
long dl_parMatSparse_GetLowerBound(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetUpperBound()
long dl_parMatSparse_GetUpperBound(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetComm()
MPI_Comm dl_parMatSparse_GetComm(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#SetValueLocal(S row, S col, T value)
void dl_parMatSparse_SetValLocal(dl_parMatSparse_t *pm, long row, long col, double val);
//! C-interface of parMatrixSparse#SetValue (S row, S col, T value)
void dl_parMatSparse_SetVal(dl_parMatSparse_t *pm, long row, long col, double val);
//! C-interface of parMatrixSparse#AddValueLocal(S row, S col, T value)
void dl_parMatSparse_AddValLocal(dl_parMatSparse_t *pm, long row, long col, double val);
//! C-interface of parMatrixSparse#AddValue(S row, S col, T value)
void dl_parMatSparse_AddVal(dl_parMatSparse_t *pm, long row, long col, double val);
//! C-interface of parMatrixSparse#GetValueLocal (S row, S col)
double dl_parMatSparse_GetValLocal(dl_parMatSparse_t *pm, long row, long col);
//! C-interface of parMatrixSparse#GetValue (S row, S col)
double dl_parMatSparse_GetVal(dl_parMatSparse_t *pm, long row, long col);
//! convert a parMatrixSparse matrix into CSR format
/*!
  * @param[in] pm the struct of parMatrixSparse to be converted
  * @param[out] nrows number of rows of the converted matrix in CSR format on each MPI proc
  * @param[out] nnz number of nnz of the converted matrix in CSR format on each MPI proc
  * @param[out] roffs an array containing the offsets of each row of CSR format
  * @param[out] cols an array containing the column index of each nnz of CSR format
  * @param[out] val an array containing all the nnz of CSR format
*/	
void dl_parMatSparse_toCSRArray(dl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, double **val);
//! C-interface of parMatrixSparse#show()
void dl_parMatSparse_Show(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#MatView()
void dl_parMatSparse_View(dl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#writeToMatrixMarket(std::string file_name)
void dl_parMatSparse_Wrt2MM(dl_parMatSparse_t *pm, char *file_name);
//! C-interface of parMatrixSparse#initMat(S diag_l, S diag_u, Base< T > scale, T shift, Base< T > sparsity)
void dl_parMatSparse_initMat(dl_parMatSparse_t *pm, long diagl, long diagu, double Scale, double Sparsity);

//dcomplex_t int
struct zs_parMatSparse;
//! Creates a type name for zs_parMatSparse which is a C-interface for parMatrixSparse with `T=std::complex<double>` and `S=int` 
typedef struct zs_parMatSparse zs_parMatSparse_t;
//! C-interface of parMatrixSparse#parMatrixSparse(parVector< T, S > vec)
zs_parMatSparse_t *new_zs_ParMatSparse_1(zs_parVec_t *pv);
//! C-interface of parMatrixSparse#parMatrixSparse(parVectorMap< S > map)
zs_parMatSparse_t *new_zs_ParMatSparse_2(parVecMapL_t *map);
//! Destory the struct
void zs_parMatSparse_destory(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNRows()
int zs_parMatSparse_GetNRows(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNCols()
int zs_parMatSparse_GetNCols(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNnzLoc()
int zs_parMatSparse_GetNNzLoc(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetLowerBound()
int zs_parMatSparse_GetLowerBound(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetUpperBound()
int zs_parMatSparse_GetUpperBound(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetComm()
MPI_Comm zs_parMatSparse_GetComm(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#SetValueLocal(S row, S col, T value)
void zs_parMatSparse_SetValLocal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val);
//! C-interface of parMatrixSparse#SetValue (S row, S col, T value)
void zs_parMatSparse_SetVal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val);
//! C-interface of parMatrixSparse#AddValueLocal(S row, S col, T value)
void zs_parMatSparse_AddValLocal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val);
//! C-interface of parMatrixSparse#AddValue(S row, S col, T value)
void zs_parMatSparse_AddVal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val);
//! C-interface of parMatrixSparse#GetValueLocal (S row, S col)
dcomplex_t zs_parMatSparse_GetValLocal(zs_parMatSparse_t *pm, int row, int col);
//! C-interface of parMatrixSparse#GetValue (S row, S col)
dcomplex_t zs_parMatSparse_GetVal(zs_parMatSparse_t *pm, int row, int col);
//! convert a parMatrixSparse matrix into CSR format
/*!
  * @param[in] pm the struct of parMatrixSparse to be converted
  * @param[out] nrows number of rows of the converted matrix in CSR format on each MPI proc
  * @param[out] nnz number of nnz of the converted matrix in CSR format on each MPI proc
  * @param[out] roffs an array containing the offsets of each row of CSR format
  * @param[out] cols an array containing the column index of each nnz of CSR format
  * @param[out] val an array containing all the nnz of CSR format
*/	
void zs_parMatSparse_toCSRArray(zs_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, dcomplex_t **val);
//! C-interface of parMatrixSparse#show()
void zs_parMatSparse_Show(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#MatView()
void zs_parMatSparse_View(zs_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#writeToMatrixMarketCmplx(std::string file_name)
void zs_parMatSparse_Wrt2MMCmplx(zs_parMatSparse_t *pm, char *file_name);
//! C-interface of parMatrixSparse#initMat(S diag_l, S diag_u, Base< T > scale, T shift, Base< T > sparsity)
void zs_parMatSparse_initMat(zs_parMatSparse_t *pm, int diagl, int diagu, double Scale, double Sparsity);

//dcomplex_t long
struct zl_parMatSparse;
//! Creates a type name for zl_parMatSparse which is a C-interface for parMatrixSparse with `T=std::complex<double>` and `S=long` 
typedef struct zl_parMatSparse zl_parMatSparse_t;
//! C-interface of parMatrixSparse#parMatrixSparse(parVector< T, S > vec)
zl_parMatSparse_t *new_zl_ParMatSparse_1(zl_parVec_t *pv);
//! C-interface of parMatrixSparse#parMatrixSparse(parVectorMap< S > map)
zl_parMatSparse_t *new_zl_ParMatSparse_2(parVecMap_t *map);
//! Destory the struct
void zl_parMatSparse_destory(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNRows()
long zl_parMatSparse_GetNRows(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNCols()
long zl_parMatSparse_GetNCols(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetNnzLoc()
long zl_parMatSparse_GetNNzLoc(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetLowerBound()
long zl_parMatSparse_GetLowerBound(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetUpperBound()
long zl_parMatSparse_GetUpperBound(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#GetComm()
MPI_Comm zl_parMatSparse_GetComm(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#SetValueLocal(S row, S col, T value)
void zl_parMatSparse_SetValLocal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val);
//! C-interface of parMatrixSparse#SetValue (S row, S col, T value)
void zl_parMatSparse_SetVal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val);
//! C-interface of parMatrixSparse#AddValueLocal(S row, S col, T value)
void zl_parMatSparse_AddValLocal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val);
//! C-interface of parMatrixSparse#AddValue(S row, S col, T value)
void zl_parMatSparse_AddVal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val);
//! C-interface of parMatrixSparse#GetValueLocal (S row, S col)
dcomplex_t zl_parMatSparse_GetValLocal(zl_parMatSparse_t *pm, long row, long col);
//! C-interface of parMatrixSparse#GetValue (S row, S col)
dcomplex_t zl_parMatSparse_GetVal(zl_parMatSparse_t *pm, long row, long col);
//! convert a parMatrixSparse matrix into CSR format
/*!
  * @param[in] pm the struct of parMatrixSparse to be converted
  * @param[out] nrows number of rows of the converted matrix in CSR format on each MPI proc
  * @param[out] nnz number of nnz of the converted matrix in CSR format on each MPI proc
  * @param[out] roffs an array containing the offsets of each row of CSR format
  * @param[out] cols an array containing the column index of each nnz of CSR format
  * @param[out] val an array containing all the nnz of CSR format
*/	
void zl_parMatSparse_toCSRArray(zl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, dcomplex_t **val);
//! C-interface of parMatrixSparse#show()
void zl_parMatSparse_Show(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#MatView()
void zl_parMatSparse_View(zl_parMatSparse_t *pm);
//! C-interface of parMatrixSparse#writeToMatrixMarketCmplx(std::string file_name)
void zl_parMatSparse_Wrt2MMCmplx(zl_parMatSparse_t *pm, char *file_name);
//! C-interface of parMatrixSparse#initMat(S diag_l, S diag_u, Base< T > scale, T shift, Base< T > sparsity)
void zl_parMatSparse_initMat(zl_parMatSparse_t *pm, long diagl, long diagu, double Scale, double Sparsity);

//interface of selected functions
//! C-interface of function specNonHerm(parVectorMap<S> index_map, std::string spectrum)
zs_parVec_t *zs_specNonHermLoad(parVecMap_t *map, char *spectrum);
//! C-interface of function specNonHerm(parVectorMap<S> index_map, std::string spectrum)
zl_parVec_t *zl_specNonHermLoad(parVecMapL_t *map, char *spectrum);
//! C-interface of function specNonSymmCplex(parVectorMap<S> index_map, std::string spectrum)
zs_parVec_t *zs_specNonSymmLoad(parVecMap_t *map, char *spectrum);
//! C-interface of function specNonSymmCplex(parVectorMap<S> index_map, std::string spectrum)
zl_parVec_t *zl_specNonSymmLoad(parVecMapL_t *map, char *spectrum);
//! C-interface of function specNonSymm(parVectorMap<S> index_map, std::string spectrum)
ds_parVec_t *ds_specNonSymmLoad(parVecMap_t *map, char *spectrum);
//! C-interface of function specNonSymm(parVectorMap<S> index_map, std::string spectrum)
dl_parVec_t *dl_specNonSymmLoad(parVecMapL_t *map, char *spectrum);

//! C-interface of function nonherm (S probSize, Nilpotent< S > nilp, initMat< S > init, std::string spectrum, MPI_Comm comm)
zs_parMatSparse_t *zs_nonherm(int probSize, nilp_t *nilp, initMatrix_t *init,  char *spectrum, MPI_Comm comm);
//! C-interface of function nonherm (S probSize, Nilpotent< S > nilp, initMat< S > init, std::string spectrum, MPI_Comm comm)
zl_parMatSparse_t *zl_nonherm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm);

//! C-interface of function nonherm (S probSize, Nilpotent< S > nilp, initMat< S > init, parVector< T, S > spec)
zs_parMatSparse_t *zs_nonherm_2(int probSize, nilp_t *nilp, initMatrix_t *init, zs_parVec_t *spec);
//! C-interface of function nonherm (S probSize, Nilpotent< S > nilp, initMat< S > init, parVector< T, S > spec)
zl_parMatSparse_t *zl_nonherm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, zl_parVec_t *spec);

//! C-interface of function nonherm (S probSize, Nilpotent< S > nilp, parMatrixSparse< T, S > *Am, parVector< T, S > spec)
void zs_nonherm_3(int probSize, nilp_t *nilp, zs_parMatSparse_t *Am, zs_parVec_t *spec);
//! C-interface of function nonherm (S probSize, Nilpotent< S > nilp, parMatrixSparse< T, S > *Am, parVector< T, S > spec)
void zl_nonherm_3(long probSize, nilpL_t *nilp, zl_parMatSparse_t *Am, zl_parVec_t *spec);

//
//! C-interface of function nonsymm (S probSize, Nilpotent< S > nilp, initMat< S > init, std::string spectrum, MPI_Comm comm) with real eigenvalues
ds_parMatSparse_t *ds_nonsymm(int probSize, nilp_t *nilp, initMatrix_t *init,  char *spectrum, MPI_Comm comm);
//! C-interface of function nonsymm (S probSize, Nilpotent< S > nilp, initMat< S > init, std::string spectrum, MPI_Comm comm) with real eigenvalues
dl_parMatSparse_t *dl_nonsymm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm);

//
//! C-interface of function nonsymm (S probSize, Nilpotent< S > nilp, initMat< S > init, parVector< T, S > spec)
ds_parMatSparse_t *ds_nonsymm_2(int probSize, nilp_t *nilp, initMatrix_t *init, ds_parVec_t *spec);
//! C-interface of function nonsymm (S probSize, Nilpotent< S > nilp, initMat< S > init, parVector< T, S > spec)
dl_parMatSparse_t *dl_nonsymm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, dl_parVec_t *spec);

//! C-interface of function nonsymm (S probSize, Nilpotent< S > nilp, parMatrixSparse< T, S > *Am, parVector< T, S > spec)
void ds_nonsymm_3(int probSize, nilp_t *nilp, ds_parMatSparse_t *Am, ds_parVec_t *spec);
//! C-interface of function nonsymm (S probSize, Nilpotent< S > nilp, parMatrixSparse< T, S > *Am, parVector< T, S > spec)
void dl_nonsymm_3(long probSize, nilpL_t *nilp, dl_parMatSparse_t *Am, dl_parVec_t *spec);
//
//! C-interface of function nonsymm (S probSize, Nilpotent< S > nilp, initMat< S > init, std::string spectrum, MPI_Comm comm) with conjugate eigenvalues
ds_parMatSparse_t *ds_nonsymmconj(int probSize, nilp_t *nilp, initMatrix_t *init,  char *spectrum, MPI_Comm comm);
//! C-interface of function nonsymm (S probSize, Nilpotent< S > nilp, initMat< S > init, std::string spectrum, MPI_Comm comm) with conjugate eigenvalues
dl_parMatSparse_t *dl_nonsymmconj(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm);

//
//! C-interface of function nonsymmconj (S probSize, Nilpotent< S > nilp, initMat< S > init, parVector< std::complex< T >, S > spec)
ds_parMatSparse_t *ds_nonsymmconj_2(int probSize, nilp_t *nilp, initMatrix_t *init, zs_parVec_t *spec);
//! C-interface of function nonsymmconj (S probSize, Nilpotent< S > nilp, initMat< S > init, parVector< std::complex< T >, S > spec)
dl_parMatSparse_t *dl_nonsymmconj_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, zl_parVec_t *spec);
//! C-interface of function nonsymmconj (S probSize, Nilpotent< S > nilp, parMatrixSparse< T, S > *Am, parVector< std::complex< T >, S > spec)
void ds_nonsymmconj_3(int probSize, nilp_t *nilp, ds_parMatSparse_t *Am, zs_parVec_t *spec);
//! C-interface of function nonsymmconj (S probSize, Nilpotent< S > nilp, parMatrixSparse< T, S > *Am, parVector< std::complex< T >, S > spec)
void dl_nonsymmconj_3(long probSize, nilpL_t *nilp, dl_parMatSparse_t *Am, zl_parVec_t *spec);

/** @} */ // end of group5

#ifdef __cplusplus
}
#endif




#endif