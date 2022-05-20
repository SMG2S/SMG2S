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

#include<cmath>
#include "c_wrapper.h"
#include "../../nilpotent/nilpotent.h"
#include "../../parVector/parVector.h"
#include "../../parMatrix/parMatrixSparse.h"
#include "../../parMatrix/initMat.h"
#include "../../smg2s/smg2s.h"


nilp_t *newNilp_1(int nbOne, int size){
	return reinterpret_cast<nilp_t *>( new Nilpotent<int>( nbOne, size ) );
}
nilp_t *newNilp_2(int nbOne, int diag, int size){
	return reinterpret_cast<nilp_t *>( new Nilpotent<int>( nbOne, diag, size ) );
}
nilp_t *newNilp_3(int *nilpvec, int size){
	std::vector<int> nv(nilpvec, nilpvec+size-1);
	return reinterpret_cast<nilp_t *>( new Nilpotent<int>( nv, size ) );
}
nilp_t *newNilp_4(int *nilpvec, int diag, int size){
	std::vector<int> nv(nilpvec, nilpvec+size-abs(diag));
	return reinterpret_cast<nilp_t *>( new Nilpotent<int>( nv, diag, size ) );
}
void nilp_destory(nilp_t *nilp){
	delete reinterpret_cast<Nilpotent<int> *>(nilp);
}
void nilp_show(nilp_t *nilp){
	return reinterpret_cast<Nilpotent<int> *>(nilp)->show();	
}
int nilp_getDegree(nilp_t *nilp){
	return reinterpret_cast<Nilpotent<int> *>(nilp)->getDegree();
}
int* nilp_getIndOfZeros(nilp_t *nilp){
	return reinterpret_cast<Nilpotent<int> *>(nilp)->getIndOfZeros().data();	
}

// long int
nilpL_t *newNilpL_1(long nbOne, long size){
	return reinterpret_cast<nilpL_t *>( new Nilpotent<long>( nbOne, size ) );	
}
nilpL_t *newNilpL_2(long nbOne, long diag, long size){
	return reinterpret_cast<nilpL_t *>( new Nilpotent<long>( nbOne, diag, size ) );
}
nilpL_t *newNilpL_3(long *nilpvec, long size){
	std::vector<long> nv(nilpvec, nilpvec+size-1);
	return reinterpret_cast<nilpL_t *>( new Nilpotent<long>( nv, size ) );	
}
nilpL_t *newNilpL_4(long *nilpvec, long diag, long size){
	std::vector<long> nv(nilpvec, nilpvec+size-abs(diag));
	return reinterpret_cast<nilpL_t *>( new Nilpotent<long>( nv, diag, size ) );
}
void nilpL_destory(nilpL_t *nilp){
	delete reinterpret_cast<Nilpotent<long> *>(nilp);
}
long nilpL_getDegree(nilpL_t *nilp){
	return reinterpret_cast<Nilpotent<long> *>(nilp)->getDegree();
}
long* nilpL_getIndOfZeros(nilpL_t *nilp){
	return reinterpret_cast<Nilpotent<long> *>(nilp)->getIndOfZeros().data();	
}
void nilpL_show(nilpL_t *nilp){
	return reinterpret_cast<Nilpotent<long> *>(nilp)->show();	
}


/*init matrix struct*/
//int
initMatrix_t *newInitMatrix_1(){
	return reinterpret_cast<initMatrix_t *>( new initMat<int>() );
}
initMatrix_t *newInitMatrix_2(int diagl, int diagu){
	return reinterpret_cast<initMatrix_t *>( new initMat<int>(diagl, diagu) );	
}
initMatrix_t *newInitMatrix_3(int diagl, int diagu, double Sparsity){
	return reinterpret_cast<initMatrix_t *>( new initMat<int>(diagl, diagu, Sparsity) );		
}
initMatrix_t *newInitMatrix_4(int diagl, int diagu, double Scale, double Sparsity){
	return reinterpret_cast<initMatrix_t *>( new initMat<int>(diagl, diagu, Scale, Sparsity) );			
}
void initMatrix_show(initMatrix_t *init){
	return reinterpret_cast<initMat<int> *>(init)->show();		
}
void initMatrix_destory(initMatrix_t *init){
	delete reinterpret_cast<initMat<int> *>(init);
}

//long
initMatrixL_t *newInitMatrixL_1(){
	return reinterpret_cast<initMatrixL_t *>( new initMat<long>() );
}
initMatrixL_t *newInitMatrixL_2(long diagl, long diagu){
	return reinterpret_cast<initMatrixL_t *>( new initMat<long>(diagl, diagu) );		
}	
initMatrixL_t *newInitMatrixL_3(long diagl, long diagu, double Sparsity){
	return reinterpret_cast<initMatrixL_t *>( new initMat<long>(diagl, diagu, Sparsity) );			
}
initMatrixL_t *newInitMatrixL_4(long diagl, long diagu, double Scale, double Sparsity){
	return reinterpret_cast<initMatrixL_t *>( new initMat<long>(diagl, diagu, Scale, Sparsity) );			
}
void initMatrixL_show(initMatrixL_t *init){
	return reinterpret_cast<initMat<long> *>(init)->show();			
}
void initMatrixL_destory(initMatrixL_t *init){
	delete reinterpret_cast<initMat<long> *>(init);
}

/*parVectorMap*/
parVecMap_t *newParVecMap(MPI_Comm ncomm, int lbound, int ubound){
	return reinterpret_cast<parVecMap_t *>( new parVectorMap<int>(ncomm, lbound, ubound) );
}
MPI_Comm parVecMapGetComm(parVecMap_t *pv){
	return reinterpret_cast<parVectorMap<int> *>(pv)->GetCurrentComm();				
}
int parVecMapL2G(parVecMap_t *pv, int local_index){
	return reinterpret_cast<parVectorMap<int> *>(pv)->Loc2Glob(local_index);					
}
int parVecMapG2L(parVecMap_t *pv, int global_index){
	return reinterpret_cast<parVectorMap<int> *>(pv)->Glob2Loc(global_index);						
}
int parVecMapGetLocSize(parVecMap_t *pv){
	return reinterpret_cast<parVectorMap<int> *>(pv)->GetLocalSize();							
}
int parVecMapGetGlobSize(parVecMap_t *pv){
	return reinterpret_cast<parVectorMap<int> *>(pv)->GetGlobalSize();								
}
void parVecMap_destory(parVecMap_t *pv){
	delete reinterpret_cast<parVectorMap<int> *>(pv);
}

parVecMapL_t *newParVecMapL(MPI_Comm ncomm, long lbound, long ubound){
	return reinterpret_cast<parVecMapL_t *>( new parVectorMap<long>(ncomm, lbound, ubound) );	
}
MPI_Comm parVecMapLGetComm(parVecMapL_t *pv){
	return reinterpret_cast<parVectorMap<long> *>(pv)->GetCurrentComm();					
}
long parVecMapLL2G(parVecMapL_t *pv, long local_index){
	return reinterpret_cast<parVectorMap<long> *>(pv)->Loc2Glob(local_index);						
}
long parVecMapLG2L(parVecMapL_t *pv, long global_index){
	return reinterpret_cast<parVectorMap<long> *>(pv)->Glob2Loc(global_index);							
}
long parVecMapLGetLocSize(parVecMapL_t *pv){
	return reinterpret_cast<parVectorMap<long> *>(pv)->GetLocalSize();								
}
long parVecMapLGetGlobSize(parVecMapL_t *pv){
	return reinterpret_cast<parVectorMap<long> *>(pv)->GetGlobalSize();									
}
void parVecMapL_destory(parVecMapL_t *pv){
	delete reinterpret_cast<parVectorMap<long> *>(pv);	
}

/*parVector*/
ds_parVec_t *new_ds_ParVec_1(MPI_Comm ncomm, int lbound, int ubound){
	return reinterpret_cast<ds_parVec_t *>( new parVector<double, int>(ncomm, lbound, ubound) );		
}
ds_parVec_t *new_ds_ParVec_2(parVecMap_t *map){
	return reinterpret_cast<ds_parVec_t *>( new parVector<double,int>(*reinterpret_cast<parVectorMap<int> *>(map)) );
}
void ds_parVec_destory(ds_parVec_t *pv){
	delete reinterpret_cast<parVector<double,int> *>(pv);
}
int ds_parVecGetLowerBound(ds_parVec_t *pv){
	return reinterpret_cast<parVector<double,int> *>(pv)->GetLowerBound();
}
int ds_parVecGetUpperBound(ds_parVec_t *pv){
	return reinterpret_cast<parVector<double, int> *>(pv)->GetUpperBound();	
}

int ds_parVecGetLocSize(ds_parVec_t *pv){
	return reinterpret_cast<parVector<double, int> *>(pv)->GetLocalSize();		
}
int ds_parVecGetGlobSize(ds_parVec_t *pv){
	return reinterpret_cast<parVector<double, int> *>(pv)->GetGlobalSize();			
}
double ds_parVecGetVal(ds_parVec_t *pv, int index){
	return reinterpret_cast<parVector<double, int> *>(pv)->GetValue(index);			
}
double ds_parVecGetValLoc(ds_parVec_t *pv, int lindex){
	return reinterpret_cast<parVector<double, int> *>(pv)->GetValueLocal(lindex);				
}
double *ds_parVecGetArray(ds_parVec_t *pv){
	return reinterpret_cast<parVector<double, int> *>(pv)->GetArray();				
}
MPI_Comm ds_parVecGetComm(ds_parVec_t *pv){
	return reinterpret_cast<parVector<double, int> *>(pv)->GetComm();				
}
int ds_parVecL2G(ds_parVec_t *pv, int local_index){
	return reinterpret_cast<parVector<double, int> *>(pv)->Loc2Glob(local_index);				
}
int ds_parVecG2L(ds_parVec_t *pv, int global_index){
	return reinterpret_cast<parVector<double, int> *>(pv)->Glob2Loc(global_index);					
}
void ds_parVecSetToVal(ds_parVec_t *pv, double val){
	return reinterpret_cast<parVector<double, int> *>(pv)->SetToValue(val);					
}
void ds_parVecView(ds_parVec_t *pv){
	return reinterpret_cast<parVector<double, int> *>(pv)->VecView();						
}
void ds_parVecSetVal(ds_parVec_t *pv, int index, double val){
	return reinterpret_cast<parVector<double, int> *>(pv)->SetValueGlobal(index, val);						
}
void ds_parVecSetValLoc(ds_parVec_t *pv, int lindex, double val){
	return reinterpret_cast<parVector<double, int> *>(pv)->SetValueLocal(lindex, val);						
}
void ds_parVecAdd(ds_parVec_t *pv, ds_parVec_t *pv2){
	return reinterpret_cast<parVector<double, int> *>(pv)->VecAdd(* reinterpret_cast<parVector<double, int> *>(pv2));						
}
void ds_parVecScale(ds_parVec_t *pv, double scale){
	reinterpret_cast<parVector<double, int> *>(pv)->VecScale(scale);	
}
double ds_parVecDot(ds_parVec_t *pv, ds_parVec_t *pv2){
	return reinterpret_cast<parVector<double, int> *>(pv)->VecDot(* reinterpret_cast<parVector<double, int> *>(pv2));							
}
void ds_parVecReadExtVec(ds_parVec_t *pv, char* spectrum){
	std::string spec(spectrum);
	reinterpret_cast<parVector<double, int> *>(pv)->ReadExtVec(spec);		
}


//
dl_parVec_t *new_dl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound){
	return reinterpret_cast<dl_parVec_t *>( new parVector<double, long>(ncomm, lbound, ubound) );		
}
dl_parVec_t *new_dl_ParVec_2(parVecMap_t *map){
	return reinterpret_cast<dl_parVec_t *>( new parVector<double,long>(*reinterpret_cast<parVectorMap<long> *>(map)) );
}
void dl_parVec_destory(dl_parVec_t *pv){
	delete reinterpret_cast<parVector<double,long> *>(pv);
}
long dl_parVecGetLowerBound(dl_parVec_t *pv){
	return reinterpret_cast<parVector<double,long> *>(pv)->GetLowerBound();
}
long dl_parVecGetUpperBound(dl_parVec_t *pv){
	return reinterpret_cast<parVector<double, long> *>(pv)->GetUpperBound();	
}

long dl_parVecGetLocSize(dl_parVec_t *pv){
	return reinterpret_cast<parVector<double, long> *>(pv)->GetLocalSize();		
}
long dl_parVecGetGlobSize(dl_parVec_t *pv){
	return reinterpret_cast<parVector<double, long> *>(pv)->GetGlobalSize();			
}
double dl_parVecGetVal(dl_parVec_t *pv, long index){
	return reinterpret_cast<parVector<double, long> *>(pv)->GetValue(index);			
}
double dl_parVecGetValLoc(dl_parVec_t *pv, long lindex){
	return reinterpret_cast<parVector<double, long> *>(pv)->GetValueLocal(lindex);				
}
double *dl_parVecGetArray(dl_parVec_t *pv){
	return reinterpret_cast<parVector<double, long> *>(pv)->GetArray();				
}
MPI_Comm dl_parVecGetComm(dl_parVec_t *pv){
	return reinterpret_cast<parVector<double, long> *>(pv)->GetComm();				
}
long dl_parVecL2G(dl_parVec_t *pv, long local_index){
	return reinterpret_cast<parVector<double, long> *>(pv)->Loc2Glob(local_index);				
}
long dl_parVecG2L(dl_parVec_t *pv, long global_index){
	return reinterpret_cast<parVector<double, long> *>(pv)->Glob2Loc(global_index);					
}
void dl_parVecSetToVal(dl_parVec_t *pv, double val){
	return reinterpret_cast<parVector<double, long> *>(pv)->SetToValue(val);					
}
void dl_parVecView(dl_parVec_t *pv){
	return reinterpret_cast<parVector<double, long> *>(pv)->VecView();						
}
void dl_parVecSetVal(dl_parVec_t *pv, long index, double val){
	return reinterpret_cast<parVector<double, long> *>(pv)->SetValueGlobal(index, val);						
}
void dl_parVecSetValLoc(dl_parVec_t *pv, long lindex, double val){
	return reinterpret_cast<parVector<double, long> *>(pv)->SetValueLocal(lindex, val);						
}
void dl_parVecAdd(dl_parVec_t *pv, dl_parVec_t *pv2){
	return reinterpret_cast<parVector<double, long> *>(pv)->VecAdd(* reinterpret_cast<parVector<double, long> *>(pv2));						
}
void dl_parVecScale(dl_parVec_t *pv, double scale){
	reinterpret_cast<parVector<double, long> *>(pv)->VecScale(scale);	
}
double dl_parVecDot(dl_parVec_t *pv, dl_parVec_t *pv2){
	return reinterpret_cast<parVector<double, long> *>(pv)->VecDot(* reinterpret_cast<parVector<double, long> *>(pv2));							
}
void dl_parVecReadExtVec(dl_parVec_t *pv, char* spectrum){
	std::string spec(spectrum);
	reinterpret_cast<parVector<double, long> *>(pv)->ReadExtVec(spec);		
}


//
//
ss_parVec_t *new_ss_ParVec_1(MPI_Comm ncomm, int lbound, int ubound){
	return reinterpret_cast<ss_parVec_t *>( new parVector<float, int>(ncomm, lbound, ubound) );		
}
ss_parVec_t *new_ss_ParVec_2(parVecMap_t *map){
	return reinterpret_cast<ss_parVec_t *>( new parVector<float,int>(*reinterpret_cast<parVectorMap<int> *>(map)) );
}
void ss_parVec_destory(ss_parVec_t *pv){
	delete reinterpret_cast<parVector<float,int> *>(pv);
}
int ss_parVecGetLowerBound(ss_parVec_t *pv){
	return reinterpret_cast<parVector<float,int> *>(pv)->GetLowerBound();
}
int ss_parVecGetUpperBound(ss_parVec_t *pv){
	return reinterpret_cast<parVector<float, int> *>(pv)->GetUpperBound();	
}

int ss_parVecGetLocSize(ss_parVec_t *pv){
	return reinterpret_cast<parVector<float, int> *>(pv)->GetLocalSize();		
}
int ss_parVecGetGlobSize(ss_parVec_t *pv){
	return reinterpret_cast<parVector<float, int> *>(pv)->GetGlobalSize();			
}
float ss_parVecGetVal(ss_parVec_t *pv, int index){
	return reinterpret_cast<parVector<float, int> *>(pv)->GetValue(index);			
}
float ss_parVecGetValLoc(ss_parVec_t *pv, int lindex){
	return reinterpret_cast<parVector<float, int> *>(pv)->GetValueLocal(lindex);				
}
float *ss_parVecGetArray(ss_parVec_t *pv){
	return reinterpret_cast<parVector<float, int> *>(pv)->GetArray();				
}
MPI_Comm ss_parVecGetComm(ss_parVec_t *pv){
	return reinterpret_cast<parVector<float, int> *>(pv)->GetComm();				
}
int ss_parVecL2G(ss_parVec_t *pv, int local_index){
	return reinterpret_cast<parVector<float, int> *>(pv)->Loc2Glob(local_index);				
}
int ss_parVecG2L(ss_parVec_t *pv, int global_index){
	return reinterpret_cast<parVector<float, int> *>(pv)->Glob2Loc(global_index);					
}
void ss_parVecSetToVal(ss_parVec_t *pv, float val){
	return reinterpret_cast<parVector<float, int> *>(pv)->SetToValue(val);					
}
void ss_parVecView(ss_parVec_t *pv){
	return reinterpret_cast<parVector<float, int> *>(pv)->VecView();						
}
void ss_parVecSetVal(ss_parVec_t *pv, int index, float val){
	return reinterpret_cast<parVector<float, int> *>(pv)->SetValueGlobal(index, val);						
}
void ss_parVecSetValLoc(ss_parVec_t *pv, int lindex, float val){
	return reinterpret_cast<parVector<float, int> *>(pv)->SetValueLocal(lindex, val);						
}
void ss_parVecAdd(ss_parVec_t *pv, ss_parVec_t *pv2){
	return reinterpret_cast<parVector<float, int> *>(pv)->VecAdd(* reinterpret_cast<parVector<float, int> *>(pv2));						
}
void ss_parVecScale(ss_parVec_t *pv, float scale){
	reinterpret_cast<parVector<float, int> *>(pv)->VecScale(scale);	
}
float ss_parVecDot(ss_parVec_t *pv, ss_parVec_t *pv2){
	return reinterpret_cast<parVector<float, int> *>(pv)->VecDot(* reinterpret_cast<parVector<float, int> *>(pv2));							
}
void ss_parVecReadExtVec(ss_parVec_t *pv, char* spectrum){
	std::string spec(spectrum);
	reinterpret_cast<parVector<float, int> *>(pv)->ReadExtVec(spec);		
}

//
sl_parVec_t *new_sl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound){
	return reinterpret_cast<sl_parVec_t *>( new parVector<float, long>(ncomm, lbound, ubound) );		
}
sl_parVec_t *new_sl_ParVec_2(parVecMap_t *map){
	return reinterpret_cast<sl_parVec_t *>( new parVector<float,long>(*reinterpret_cast<parVectorMap<long> *>(map)) );
}
void sl_parVec_destory(sl_parVec_t *pv){
	delete reinterpret_cast<parVector<float,long> *>(pv);
}
long sl_parVecGetLowerBound(sl_parVec_t *pv){
	return reinterpret_cast<parVector<float,long> *>(pv)->GetLowerBound();
}
long sl_parVecGetUpperBound(sl_parVec_t *pv){
	return reinterpret_cast<parVector<float, long> *>(pv)->GetUpperBound();	
}

long sl_parVecGetLocSize(sl_parVec_t *pv){
	return reinterpret_cast<parVector<float, long> *>(pv)->GetLocalSize();		
}
long sl_parVecGetGlobSize(sl_parVec_t *pv){
	return reinterpret_cast<parVector<float, long> *>(pv)->GetGlobalSize();			
}
float sl_parVecGetVal(sl_parVec_t *pv, long index){
	return reinterpret_cast<parVector<float, long> *>(pv)->GetValue(index);			
}
float sl_parVecGetValLoc(sl_parVec_t *pv, long lindex){
	return reinterpret_cast<parVector<float, long> *>(pv)->GetValueLocal(lindex);				
}
float *sl_parVecGetArray(sl_parVec_t *pv){
	return reinterpret_cast<parVector<float, long> *>(pv)->GetArray();				
}
MPI_Comm sl_parVecGetComm(sl_parVec_t *pv){
	return reinterpret_cast<parVector<float, long> *>(pv)->GetComm();				
}
long sl_parVecL2G(sl_parVec_t *pv, long local_index){
	return reinterpret_cast<parVector<float, long> *>(pv)->Loc2Glob(local_index);				
}
long sl_parVecG2L(sl_parVec_t *pv, long global_index){
	return reinterpret_cast<parVector<float, long> *>(pv)->Glob2Loc(global_index);					
}
void sl_parVecSetToVal(sl_parVec_t *pv, float val){
	return reinterpret_cast<parVector<float, long> *>(pv)->SetToValue(val);					
}
void sl_parVecView(sl_parVec_t *pv){
	return reinterpret_cast<parVector<float, long> *>(pv)->VecView();						
}
void sl_parVecSetVal(sl_parVec_t *pv, long index, float val){
	return reinterpret_cast<parVector<float, long> *>(pv)->SetValueGlobal(index, val);						
}
void sl_parVecSetValLoc(sl_parVec_t *pv, long lindex, float val){
	return reinterpret_cast<parVector<float, long> *>(pv)->SetValueLocal(lindex, val);						
}
void sl_parVecAdd(sl_parVec_t *pv, sl_parVec_t *pv2){
	return reinterpret_cast<parVector<float, long> *>(pv)->VecAdd(* reinterpret_cast<parVector<float, long> *>(pv2));						
}
void sl_parVecScale(sl_parVec_t *pv, float scale){
	reinterpret_cast<parVector<float, long> *>(pv)->VecScale(scale);	
}
float sl_parVecDot(sl_parVec_t *pv, sl_parVec_t *pv2){
	return reinterpret_cast<parVector<float, long> *>(pv)->VecDot(* reinterpret_cast<parVector<float, long> *>(pv2));							
}
void sl_parVecReadExtVec(sl_parVec_t *pv, char* spectrum){
	std::string spec(spectrum);
	reinterpret_cast<parVector<float, long> *>(pv)->ReadExtVec(spec);		
}


//
zs_parVec_t *new_zs_ParVec_1(MPI_Comm ncomm, int lbound, int ubound){
	return reinterpret_cast<zs_parVec_t *>( new parVector<std::complex<double>, int>(ncomm, lbound, ubound) );		
}
zs_parVec_t *new_zs_ParVec_2(parVecMap_t *map){
	return reinterpret_cast<zs_parVec_t *>( new parVector<std::complex<double>,int>(*reinterpret_cast<parVectorMap<int> *>(map)) );
}
void zs_parVec_destory(zs_parVec_t *pv){
	delete reinterpret_cast<parVector<std::complex<double>,int> *>(pv);
}
int zs_parVecGetLowerBound(zs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>,int> *>(pv)->GetLowerBound();
}
int zs_parVecGetUpperBound(zs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->GetUpperBound();	
}

int zs_parVecGetLocSize(zs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->GetLocalSize();		
}
int zs_parVecGetGlobSize(zs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->GetGlobalSize();			
}
dcomplex_t zs_parVecGetVal(zs_parVec_t *pv, int index){
	std::complex<double> v = reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->GetValue(index);
	dcomplex_t val;
	val.real = v.real();
	val.imag = v.imag();
	return val;			
}
dcomplex_t zs_parVecGetValLoc(zs_parVec_t *pv, int lindex){
	std::complex<double> v = reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->GetValueLocal(lindex);
	dcomplex_t val;
	val.real = v.real();
	val.imag = v.imag();
	return val;	
}
dcomplex_t *zs_parVecGetArray(zs_parVec_t *pv){
	int nnz = reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->GetLocalSize();
	std::complex<double> *a = reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->GetArray();
	dcomplex_t *array = new dcomplex_t [nnz];
	for(int i = 0; i < nnz; i++){
		array[i].real = a[i].real();
		array[i].imag = a[i].imag();
	}
	return array;				
}

MPI_Comm zs_parVecGetComm(zs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->GetComm();				
}
int zs_parVecL2G(zs_parVec_t *pv, int local_index){
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->Loc2Glob(local_index);				
}
int zs_parVecG2L(zs_parVec_t *pv, int global_index){
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->Glob2Loc(global_index);					
}
void zs_parVecSetToVal(zs_parVec_t *pv, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->SetToValue(v);					
}
void zs_parVecView(zs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->VecView();						
}
void zs_parVecSetVal(zs_parVec_t *pv, int index, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->SetValueGlobal(index, v);						
}
void zs_parVecSetValLoc(zs_parVec_t *pv, int lindex, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->SetValueLocal(lindex, v);						
}
void zs_parVecAdd(zs_parVec_t *pv, zs_parVec_t *pv2){
	return reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->VecAdd(* reinterpret_cast<parVector<std::complex<double>, int> *>(pv2));						
}
void zs_parVecScale(zs_parVec_t *pv, dcomplex_t scale){
	std::complex<double> s(scale.real, scale.imag);
	reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->VecScale(s);	
}
dcomplex_t zs_parVecDot(zs_parVec_t *pv, zs_parVec_t *pv2){
	std::complex<double> d = reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->VecDot(* reinterpret_cast<parVector<std::complex<double>, int> *>(pv2));	
	dcomplex_t pp;
	pp.real = d.real();
	pp.imag = d.imag();
	return pp;							
}
void zs_parVecReadExtVec(zs_parVec_t *pv, char* spectrum){
	std::string spec(spectrum);
	reinterpret_cast<parVector<std::complex<double>, int> *>(pv)->ReadExtVec(spec);		
}


