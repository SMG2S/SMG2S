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

#include <smg2s-interface.hpp>
#include<cmath>
#include <C/c-smg2s.h>

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


//
zl_parVec_t *new_zl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound){
	return reinterpret_cast<zl_parVec_t *>( new parVector<std::complex<double>, long>(ncomm, lbound, ubound) );		
}
zl_parVec_t *new_zl_ParVec_2(parVecMap_t *map){
	return reinterpret_cast<zl_parVec_t *>( new parVector<std::complex<double>,long>(*reinterpret_cast<parVectorMap<long> *>(map)) );
}
void zl_parVec_destory(zl_parVec_t *pv){
	delete reinterpret_cast<parVector<std::complex<double>,long> *>(pv);
}
long zl_parVecGetLowerBound(zl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>,long> *>(pv)->GetLowerBound();
}
long zl_parVecGetUpperBound(zl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->GetUpperBound();	
}

long zl_parVecGetLocSize(zl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->GetLocalSize();		
}
long zl_parVecGetGlobSize(zl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->GetGlobalSize();			
}
dcomplex_t zl_parVecGetVal(zl_parVec_t *pv, long index){
	std::complex<double> v = reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->GetValue(index);
	dcomplex_t val;
	val.real = v.real();
	val.imag = v.imag();
	return val;			
}
dcomplex_t zl_parVecGetValLoc(zl_parVec_t *pv, long lindex){
	std::complex<double> v = reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->GetValueLocal(lindex);
	dcomplex_t val;
	val.real = v.real();
	val.imag = v.imag();
	return val;	
}
dcomplex_t *zl_parVecGetArray(zl_parVec_t *pv){
	long nnz = reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->GetLocalSize();
	std::complex<double> *a = reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->GetArray();
	dcomplex_t *array = new dcomplex_t [nnz];
	for(long i = 0; i < nnz; i++){
		array[i].real = a[i].real();
		array[i].imag = a[i].imag();
	}
	return array;				
}

MPI_Comm zl_parVecGetComm(zl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->GetComm();				
}
long zl_parVecL2G(zl_parVec_t *pv, long local_index){
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->Loc2Glob(local_index);				
}
long zl_parVecG2L(zl_parVec_t *pv, long global_index){
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->Glob2Loc(global_index);					
}
void zl_parVecSetToVal(zl_parVec_t *pv, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->SetToValue(v);					
}
void zl_parVecView(zl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->VecView();						
}
void zl_parVecSetVal(zl_parVec_t *pv, long index, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->SetValueGlobal(index, v);						
}
void zl_parVecSetValLoc(zl_parVec_t *pv, long lindex, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->SetValueLocal(lindex, v);						
}
void zl_parVecAdd(zl_parVec_t *pv, zl_parVec_t *pv2){
	return reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->VecAdd(* reinterpret_cast<parVector<std::complex<double>, long> *>(pv2));						
}
void zl_parVecScale(zl_parVec_t *pv, dcomplex_t scale){
	std::complex<double> s(scale.real, scale.imag);
	reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->VecScale(s);	
}
dcomplex_t zl_parVecDot(zl_parVec_t *pv, zl_parVec_t *pv2){
	std::complex<double> d = reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->VecDot(* reinterpret_cast<parVector<std::complex<double>, long> *>(pv2));	
	dcomplex_t pp;
	pp.real = d.real();
	pp.imag = d.imag();
	return pp;							
}
void zl_parVecReadExtVec(zl_parVec_t *pv, char* spectrum){
	std::string spec(spectrum);
	reinterpret_cast<parVector<std::complex<double>, long> *>(pv)->ReadExtVec(spec);		
}

//
//
cs_parVec_t *new_cs_ParVec_1(MPI_Comm ncomm, int lbound, int ubound){
	return reinterpret_cast<cs_parVec_t *>( new parVector<std::complex<float>, int>(ncomm, lbound, ubound) );		
}
cs_parVec_t *new_cs_ParVec_2(parVecMap_t *map){
	return reinterpret_cast<cs_parVec_t *>( new parVector<std::complex<float>,int>(*reinterpret_cast<parVectorMap<int> *>(map)) );
}
void cs_parVec_destory(cs_parVec_t *pv){
	delete reinterpret_cast<parVector<std::complex<float>,int> *>(pv);
}
int cs_parVecGetLowerBound(cs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>,int> *>(pv)->GetLowerBound();
}
int cs_parVecGetUpperBound(cs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->GetUpperBound();	
}

int cs_parVecGetLocSize(cs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->GetLocalSize();		
}
int cs_parVecGetGlobSize(cs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->GetGlobalSize();			
}
fcomplex_t cs_parVecGetVal(cs_parVec_t *pv, int index){
	std::complex<float> v = reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->GetValue(index);
	fcomplex_t val;
	val.real = v.real();
	val.imag = v.imag();
	return val;			
}
fcomplex_t cs_parVecGetValLoc(cs_parVec_t *pv, int lindex){
	std::complex<float> v = reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->GetValueLocal(lindex);
	fcomplex_t val;
	val.real = v.real();
	val.imag = v.imag();
	return val;	
}
fcomplex_t *cs_parVecGetArray(cs_parVec_t *pv){
	int nnz = reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->GetLocalSize();
	std::complex<float> *a = reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->GetArray();
	fcomplex_t *array = new fcomplex_t [nnz];
	for(int i = 0; i < nnz; i++){
		array[i].real = a[i].real();
		array[i].imag = a[i].imag();
	}
	return array;				
}

MPI_Comm cs_parVecGetComm(cs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->GetComm();				
}
int cs_parVecL2G(cs_parVec_t *pv, int local_index){
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->Loc2Glob(local_index);				
}
int cs_parVecG2L(cs_parVec_t *pv, int global_index){
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->Glob2Loc(global_index);					
}
void cs_parVecSetToVal(cs_parVec_t *pv, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->SetToValue(v);					
}
void cs_parVecView(cs_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->VecView();						
}
void cs_parVecSetVal(cs_parVec_t *pv, int index, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->SetValueGlobal(index, v);						
}
void cs_parVecSetValLoc(cs_parVec_t *pv, int lindex, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->SetValueLocal(lindex, v);						
}
void cs_parVecAdd(cs_parVec_t *pv, cs_parVec_t *pv2){
	return reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->VecAdd(* reinterpret_cast<parVector<std::complex<float>, int> *>(pv2));						
}
void cs_parVecScale(cs_parVec_t *pv, fcomplex_t scale){
	std::complex<float> s(scale.real, scale.imag);
	reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->VecScale(s);	
}
fcomplex_t cs_parVecDot(cs_parVec_t *pv, cs_parVec_t *pv2){
	std::complex<float> d = reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->VecDot(* reinterpret_cast<parVector<std::complex<float>, int> *>(pv2));	
	fcomplex_t pp;
	pp.real = d.real();
	pp.imag = d.imag();
	return pp;							
}
void cs_parVecReadExtVec(cs_parVec_t *pv, char* spectrum){
	std::string spec(spectrum);
	reinterpret_cast<parVector<std::complex<float>, int> *>(pv)->ReadExtVec(spec);		
}


//
cl_parVec_t *new_cl_ParVec_1(MPI_Comm ncomm, long lbound, long ubound){
	return reinterpret_cast<cl_parVec_t *>( new parVector<std::complex<float>, long>(ncomm, lbound, ubound) );		
}
cl_parVec_t *new_cl_ParVec_2(parVecMap_t *map){
	return reinterpret_cast<cl_parVec_t *>( new parVector<std::complex<float>,long>(*reinterpret_cast<parVectorMap<long> *>(map)) );
}
void cl_parVec_destory(cl_parVec_t *pv){
	delete reinterpret_cast<parVector<std::complex<float>,long> *>(pv);
}
long cl_parVecGetLowerBound(cl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>,long> *>(pv)->GetLowerBound();
}
long cl_parVecGetUpperBound(cl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->GetUpperBound();	
}

long cl_parVecGetLocSize(cl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->GetLocalSize();		
}
long cl_parVecGetGlobSize(cl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->GetGlobalSize();			
}
fcomplex_t cl_parVecGetVal(cl_parVec_t *pv, long index){
	std::complex<float> v = reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->GetValue(index);
	fcomplex_t val;
	val.real = v.real();
	val.imag = v.imag();
	return val;			
}
fcomplex_t cl_parVecGetValLoc(cl_parVec_t *pv, long lindex){
	std::complex<float> v = reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->GetValueLocal(lindex);
	fcomplex_t val;
	val.real = v.real();
	val.imag = v.imag();
	return val;	
}
fcomplex_t *cl_parVecGetArray(cl_parVec_t *pv){
	long nnz = reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->GetLocalSize();
	std::complex<float> *a = reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->GetArray();
	fcomplex_t *array = new fcomplex_t [nnz];
	for(long i = 0; i < nnz; i++){
		array[i].real = a[i].real();
		array[i].imag = a[i].imag();
	}
	return array;				
}

MPI_Comm cl_parVecGetComm(cl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->GetComm();				
}
long cl_parVecL2G(cl_parVec_t *pv, long local_index){
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->Loc2Glob(local_index);				
}
long cl_parVecG2L(cl_parVec_t *pv, long global_index){
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->Glob2Loc(global_index);					
}
void cl_parVecSetToVal(cl_parVec_t *pv, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->SetToValue(v);					
}
void cl_parVecView(cl_parVec_t *pv){
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->VecView();						
}
void cl_parVecSetVal(cl_parVec_t *pv, long index, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->SetValueGlobal(index, v);						
}
void cl_parVecSetValLoc(cl_parVec_t *pv, long lindex, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->SetValueLocal(lindex, v);						
}
void cl_parVecAdd(cl_parVec_t *pv, cl_parVec_t *pv2){
	return reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->VecAdd(* reinterpret_cast<parVector<std::complex<float>, long> *>(pv2));						
}
void cl_parVecScale(cl_parVec_t *pv, fcomplex_t scale){
	std::complex<float> s(scale.real, scale.imag);
	reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->VecScale(s);	
}
fcomplex_t cl_parVecDot(cl_parVec_t *pv, cl_parVec_t *pv2){
	std::complex<float> d = reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->VecDot(* reinterpret_cast<parVector<std::complex<float>, long> *>(pv2));	
	fcomplex_t pp;
	pp.real = d.real();
	pp.imag = d.imag();
	return pp;							
}
void cl_parVecReadExtVec(cl_parVec_t *pv, char* spectrum){
	std::string spec(spectrum);
	reinterpret_cast<parVector<std::complex<float>, long> *>(pv)->ReadExtVec(spec);		
}


/*parMatrixSparse*/
//double int
ds_parMatSparse_t *new_ds_ParMatSparse_1(ds_parVec_t *pv){
	return reinterpret_cast<ds_parMatSparse_t *>( new parMatrixSparse<double,int >(*reinterpret_cast<parVector<double, int > *>(pv)) );
}
ds_parMatSparse_t *new_ds_ParMatSparse_2(parVecMap_t *map){
	return reinterpret_cast<ds_parMatSparse_t *>( new parMatrixSparse<double,int >(*reinterpret_cast<parVectorMap<int > *>(map)) );
}
void ds_parMatSparse_destory(ds_parMatSparse_t *pm){
	delete reinterpret_cast<parMatrixSparse<double, int > *>(pm);
}
int ds_parMatSparse_GetNRows(ds_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetNRows();
}
int ds_parMatSparse_GetNCols(ds_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetNCols();	
}
int ds_parMatSparse_GetNNzLoc(ds_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetNNzLoc();	
}
int ds_parMatSparse_GetLowerBound(ds_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetLowerBound();		
}
int ds_parMatSparse_GetUpperBound(ds_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetUpperBound();			
}
MPI_Comm ds_parMatSparse_GetComm(ds_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetComm();				
}
void ds_parMatSparse_SetValLocal(ds_parMatSparse_t *pm, int row, int col, double val){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->SetValueLocal(row, col, val);					
}
void ds_parMatSparse_SetVal(ds_parMatSparse_t *pm, int row, int col, double val){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->SetValue(row, col, val);						
}
void ds_parMatSparse_AddValLocal(ds_parMatSparse_t *pm, int row, int col, double val){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->AddValueLocal(row, col, val);						
}
void ds_parMatSparse_AddVal(ds_parMatSparse_t *pm, int row, int col, double val){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->AddValue(row, col, val);						
}
double ds_parMatSparse_GetValLocal(ds_parMatSparse_t *pm, int row, int col){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetValueLocal(row, col);							
}
double ds_parMatSparse_GetVal(ds_parMatSparse_t *pm, int row, int col){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetValue(row, col);								
}
void ds_parMatSparse_toCSRArray(ds_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, double **val){
	MatrixCSR<double, int > csr = reinterpret_cast<parMatrixSparse<double, int > *>(pm)->ConvertToCSR();
	*nrows = csr.nrows;
	*nnz = csr.nnz;
	*roffs = (int *)malloc( (csr.nrows + 1) * sizeof(int ) );
	*cols = (int *)malloc( csr.nnz * sizeof(int ) );
	*val = (double *)malloc(csr.nnz * sizeof(double));

	memcpy(*roffs, csr.rows.data(), (csr.nrows + 1) * sizeof(int ));
	memcpy(*cols, csr.cols.data(), csr.nnz * sizeof(int ));
	memcpy(*val, csr.vals.data(), csr.nnz * sizeof(double));

}
void ds_parMatSparse_Show(ds_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->show();						
}
void ds_parMatSparse_View(ds_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->MatView();						
}
void ds_parMatSparse_Wrt2MM(ds_parMatSparse_t *pm, char *file_name){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->writeToMatrixMarket(file_name);						
}

//
//double long
dl_parMatSparse_t *new_dl_ParMatSparse_1(dl_parVec_t *pv){
	return reinterpret_cast<dl_parMatSparse_t *>( new parMatrixSparse<double,long >(*reinterpret_cast<parVector<double, long > *>(pv)) );
}
dl_parMatSparse_t *new_dl_ParMatSparse_2(parVecMap_t *map){
	return reinterpret_cast<dl_parMatSparse_t *>( new parMatrixSparse<double,long >(*reinterpret_cast<parVectorMap<long > *>(map)) );
}
void dl_parMatSparse_destory(dl_parMatSparse_t *pm){
	delete reinterpret_cast<parMatrixSparse<double, long > *>(pm);
}
long dl_parMatSparse_GetNRows(dl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->GetNRows();
}
long dl_parMatSparse_GetNCols(dl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->GetNCols();	
}
long dl_parMatSparse_GetNNzLoc(dl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->GetNNzLoc();	
}
long dl_parMatSparse_GetLowerBound(dl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->GetLowerBound();		
}
long dl_parMatSparse_GetUpperBound(dl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->GetUpperBound();			
}
MPI_Comm dl_parMatSparse_GetComm(dl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->GetComm();				
}
void dl_parMatSparse_SetValLocal(dl_parMatSparse_t *pm, long row, long col, double val){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->SetValueLocal(row, col, val);					
}
void dl_parMatSparse_SetVal(dl_parMatSparse_t *pm, long row, long col, double val){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->SetValue(row, col, val);						
}
void dl_parMatSparse_AddValLocal(dl_parMatSparse_t *pm, long row, long col, double val){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->AddValueLocal(row, col, val);						
}
void dl_parMatSparse_AddVal(dl_parMatSparse_t *pm, long row, long col, double val){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->AddValue(row, col, val);						
}
double dl_parMatSparse_GetValLocal(dl_parMatSparse_t *pm, long row, long col){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetValueLocal(row, col);								
}
double dl_parMatSparse_GetVal(dl_parMatSparse_t *pm, long row, long col){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetValue(row, col);								
}
void dl_parMatSparse_toCSRArray(dl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, double **val){
	MatrixCSR<double, long > csr = reinterpret_cast<parMatrixSparse<double, long > *>(pm)->ConvertToCSR();
	*nrows = csr.nrows;
	*nnz = csr.nnz;
	*roffs = (long *)malloc( (csr.nrows + 1) * sizeof(long ) );
	*cols = (long *)malloc( csr.nnz * sizeof(long ) );
	*val = (double *)malloc(csr.nnz * sizeof(double));

	memcpy(*roffs, csr.rows.data(), (csr.nrows + 1) * sizeof(long ));
	memcpy(*cols, csr.cols.data(), csr.nnz * sizeof(long));
	memcpy(*val, csr.vals.data(), csr.nnz * sizeof(double));

}
void dl_parMatSparse_Show(dl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->show();						
}
void dl_parMatSparse_View(dl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->MatView();						
}
void dl_parMatSparse_Wrt2MM(dl_parMatSparse_t *pm, char *file_name){
	return reinterpret_cast<parMatrixSparse<double, long > *>(pm)->writeToMatrixMarket(file_name);						
}

//
//float int
ss_parMatSparse_t *new_ss_ParMatSparse_1(ss_parVec_t *pv){
	return reinterpret_cast<ss_parMatSparse_t *>( new parMatrixSparse<float,int >(*reinterpret_cast<parVector<float, int > *>(pv)) );
}
ss_parMatSparse_t *new_ss_ParMatSparse_2(parVecMap_t *map){
	return reinterpret_cast<ss_parMatSparse_t *>( new parMatrixSparse<float,int >(*reinterpret_cast<parVectorMap<int > *>(map)) );
}
void ss_parMatSparse_destory(ss_parMatSparse_t *pm){
	delete reinterpret_cast<parMatrixSparse<float, int > *>(pm);
}
int ss_parMatSparse_GetNRows(ss_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->GetNRows();
}
int ss_parMatSparse_GetNCols(ss_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->GetNCols();	
}
int ss_parMatSparse_GetNNzLoc(ss_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->GetNNzLoc();	
}
int ss_parMatSparse_GetLowerBound(ss_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->GetLowerBound();		
}
int ss_parMatSparse_GetUpperBound(ss_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->GetUpperBound();			
}
MPI_Comm ss_parMatSparse_GetComm(ss_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->GetComm();				
}
void ss_parMatSparse_SetValLocal(ss_parMatSparse_t *pm, int row, int col, float val){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->SetValueLocal(row, col, val);					
}
void ss_parMatSparse_SetVal(ss_parMatSparse_t *pm, int row, int col, float val){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->SetValue(row, col, val);						
}
void ss_parMatSparse_AddValLocal(ss_parMatSparse_t *pm, int row, int col, float val){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->AddValueLocal(row, col, val);						
}
void ss_parMatSparse_AddVal(ss_parMatSparse_t *pm, int row, int col, float val){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->AddValue(row, col, val);						
}
float ss_parMatSparse_GetValLocal(ss_parMatSparse_t *pm, int row, int col){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetValueLocal(row, col);	
}
float ss_parMatSparse_GetVal(ss_parMatSparse_t *pm, int row, int col){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetValue(row, col);								
}
void ss_parMatSparse_toCSRArray(ss_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, float **val){
	MatrixCSR<float, int > csr = reinterpret_cast<parMatrixSparse<float, int > *>(pm)->ConvertToCSR();
	*nrows = csr.nrows;
	*nnz = csr.nnz;
	*roffs = (int *)malloc( (csr.nrows + 1) * sizeof(int ) );
	*cols = (int *)malloc( csr.nnz * sizeof(int ) );
	*val = (float *)malloc(csr.nnz * sizeof(float));

	memcpy(*roffs, csr.rows.data(), (csr.nrows + 1) * sizeof(int ));
	memcpy(*cols, csr.cols.data(), csr.nnz * sizeof(int ));
	memcpy(*val, csr.vals.data(), csr.nnz * sizeof(float));

}
void ss_parMatSparse_Show(ss_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->show();						
}
void ss_parMatSparse_View(ss_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->MatView();						
}
void ss_parMatSparse_Wrt2MM(ss_parMatSparse_t *pm, char *file_name){
	return reinterpret_cast<parMatrixSparse<float, int > *>(pm)->writeToMatrixMarket(file_name);						
}

//float long
sl_parMatSparse_t *new_sl_ParMatSparse_1(sl_parVec_t *pv){
	return reinterpret_cast<sl_parMatSparse_t *>( new parMatrixSparse<float,long >(*reinterpret_cast<parVector<float, long > *>(pv)) );
}
sl_parMatSparse_t *new_sl_ParMatSparse_2(parVecMap_t *map){
	return reinterpret_cast<sl_parMatSparse_t *>( new parMatrixSparse<float,long >(*reinterpret_cast<parVectorMap<long > *>(map)) );
}
void sl_parMatSparse_destory(sl_parMatSparse_t *pm){
	delete reinterpret_cast<parMatrixSparse<float, long > *>(pm);
}
long sl_parMatSparse_GetNRows(sl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->GetNRows();
}
long sl_parMatSparse_GetNCols(sl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->GetNCols();	
}
long sl_parMatSparse_GetNNzLoc(sl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->GetNNzLoc();	
}
long sl_parMatSparse_GetLowerBound(sl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->GetLowerBound();		
}
long sl_parMatSparse_GetUpperBound(sl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->GetUpperBound();			
}
MPI_Comm sl_parMatSparse_GetComm(sl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->GetComm();				
}
void sl_parMatSparse_SetValLocal(sl_parMatSparse_t *pm, long row, long col, float val){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->SetValueLocal(row, col, val);					
}
void sl_parMatSparse_SetVal(sl_parMatSparse_t *pm, long row, long col, float val){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->SetValue(row, col, val);						
}
void sl_parMatSparse_AddValLocal(sl_parMatSparse_t *pm, long row, long col, float val){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->AddValueLocal(row, col, val);						
}
void sl_parMatSparse_AddVal(sl_parMatSparse_t *pm, long row, long col, float val){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->AddValue(row, col, val);						
}
float sl_parMatSparse_GetValLocal(sl_parMatSparse_t *pm, long row, long col){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetValueLocal(row, col);								
}
float sl_parMatSparse_GetVal(sl_parMatSparse_t *pm, long row, long col){
	return reinterpret_cast<parMatrixSparse<double, int > *>(pm)->GetValue(row, col);								
}
void sl_parMatSparse_toCSRArray(sl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, float **val){
	MatrixCSR<float, long > csr = reinterpret_cast<parMatrixSparse<float, long > *>(pm)->ConvertToCSR();
	*nrows = csr.nrows;
	*nnz = csr.nnz;
	*roffs = (long *)malloc( (csr.nrows + 1) * sizeof(long ) );
	*cols = (long *)malloc( csr.nnz * sizeof(long ) );
	*val = (float *)malloc(csr.nnz * sizeof(float));

	memcpy(*roffs, csr.rows.data(), (csr.nrows + 1) * sizeof(long ));
	memcpy(*cols, csr.cols.data(), csr.nnz * sizeof(long));
	memcpy(*val, csr.vals.data(), csr.nnz * sizeof(float));

}
void sl_parMatSparse_Show(sl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->show();						
}
void sl_parMatSparse_View(sl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->MatView();						
}
void sl_parMatSparse_Wrt2MM(sl_parMatSparse_t *pm, char *file_name){
	return reinterpret_cast<parMatrixSparse<float, long > *>(pm)->writeToMatrixMarket(file_name);						
}


//dcomplex_t int
zs_parMatSparse_t *new_zs_ParMatSparse_1(zs_parVec_t *pv){
	return reinterpret_cast<zs_parMatSparse_t *>( new parMatrixSparse<std::complex<double>,int >(*reinterpret_cast<parVector<std::complex<double>, int > *>(pv)) );
}
zs_parMatSparse_t *new_zs_ParMatSparse_2(parVecMap_t *map){
	return reinterpret_cast<zs_parMatSparse_t *>( new parMatrixSparse<std::complex<double>,int >(*reinterpret_cast<parVectorMap<int > *>(map)) );
}
void zs_parMatSparse_destory(zs_parMatSparse_t *pm){
	delete reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm);
}
int zs_parMatSparse_GetNRows(zs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->GetNRows();
}
int zs_parMatSparse_GetNCols(zs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->GetNCols();	
}
int zs_parMatSparse_GetNNzLoc(zs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->GetNNzLoc();	
}
int zs_parMatSparse_GetLowerBound(zs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->GetLowerBound();		
}
int zs_parMatSparse_GetUpperBound(zs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->GetUpperBound();			
}
MPI_Comm zs_parMatSparse_GetComm(zs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->GetComm();				
}
void zs_parMatSparse_SetValLocal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->SetValueLocal(row, col, v);					
}
void zs_parMatSparse_SetVal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->SetValue(row, col, v);						
}
void zs_parMatSparse_AddValLocal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->AddValueLocal(row, col, v);						
}
void zs_parMatSparse_AddVal(zs_parMatSparse_t *pm, int row, int col, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->AddValue(row, col, v);						
}
dcomplex_t zs_parMatSparse_GetValLocal(zs_parMatSparse_t *pm, int row, int col){
	auto val = reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->GetValueLocal(row, col);	
	return {val.real(), val.imag()};								
}
dcomplex_t zs_parMatSparse_GetVal(zs_parMatSparse_t *pm, int row, int col){
	auto val = reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->GetValue(row, col);	
	return {val.real(), val.imag()};																
}
void zs_parMatSparse_toCSRArray(zs_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, dcomplex_t **val){
	MatrixCSR<std::complex<double>, int > csr = reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->ConvertToCSR();
	*nrows = csr.nrows;
	*nnz = csr.nnz;
	*roffs = (int *)malloc( (csr.nrows + 1) * sizeof(int ) );
	*cols = (int *)malloc( csr.nnz * sizeof(int ) );
	*val = (dcomplex_t *)malloc(csr.nnz * sizeof(dcomplex_t));

	memcpy(*roffs, csr.rows.data(), (csr.nrows + 1) * sizeof(int ));
	memcpy(*cols, csr.cols.data(), csr.nnz * sizeof(int ));
	for(auto i = 0; i < csr.nnz; i++){
		*(*val+i) = {csr.vals[i].real(), csr.vals[i].imag()};
	}	
}
void zs_parMatSparse_Show(zs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->show();						
}
void zs_parMatSparse_View(zs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->MatView();						
}
void zs_parMatSparse_Wrt2MMCmplx(zs_parMatSparse_t *pm, char *file_name){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, int > *>(pm)->writeToMatrixMarketCmplx(file_name);						
}

//dcomplex_t long
zl_parMatSparse_t *new_zl_ParMatSparse_1(zl_parVec_t *pv){
	return reinterpret_cast<zl_parMatSparse_t *>( new parMatrixSparse<std::complex<double>,long >(*reinterpret_cast<parVector<std::complex<double>, long > *>(pv)) );
}
zl_parMatSparse_t *new_zl_ParMatSparse_2(parVecMap_t *map){
	return reinterpret_cast<zl_parMatSparse_t *>( new parMatrixSparse<std::complex<double>,long >(*reinterpret_cast<parVectorMap<long > *>(map)) );
}
void zl_parMatSparse_destory(zl_parMatSparse_t *pm){
	delete reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm);
}
long zl_parMatSparse_GetNRows(zl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->GetNRows();
}
long zl_parMatSparse_GetNCols(zl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->GetNCols();	
}
long zl_parMatSparse_GetNNzLoc(zl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->GetNNzLoc();	
}
long zl_parMatSparse_GetLowerBound(zl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->GetLowerBound();		
}
long zl_parMatSparse_GetUpperBound(zl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->GetUpperBound();			
}
MPI_Comm zl_parMatSparse_GetComm(zl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->GetComm();				
}
void zl_parMatSparse_SetValLocal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->SetValueLocal(row, col, v);					
}
void zl_parMatSparse_SetVal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->SetValue(row, col, v);						
}
void zl_parMatSparse_AddValLocal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->AddValueLocal(row, col, v);						
}
void zl_parMatSparse_AddVal(zl_parMatSparse_t *pm, long row, long col, dcomplex_t val){
	std::complex<double> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->AddValue(row, col, v);						
}
dcomplex_t zl_parMatSparse_GetValLocal(zl_parMatSparse_t *pm, long row, long col){
	auto val = reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->GetValueLocal(row, col);	
	return {val.real(), val.imag()};								
}
dcomplex_t zl_parMatSparse_GetVal(zl_parMatSparse_t *pm, long row, long col){
	auto val = reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->GetValue(row, col);	
	return {val.real(), val.imag()};																
}
void zl_parMatSparse_toCSRArray(zl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, dcomplex_t **val){
	MatrixCSR<std::complex<double>, long > csr = reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->ConvertToCSR();
	*nrows = csr.nrows;
	*nnz = csr.nnz;
	*roffs = (long *)malloc( (csr.nrows + 1) * sizeof(long ) );
	*cols = (long *)malloc( csr.nnz * sizeof(long ) );
	*val = (dcomplex_t *)malloc(csr.nnz * sizeof(dcomplex_t));

	memcpy(*roffs, csr.rows.data(), (csr.nrows + 1) * sizeof(long ));
	memcpy(*cols, csr.cols.data(), csr.nnz * sizeof(long ));
	for(auto i = 0; i < csr.nnz; i++){
		*(*val+i) = {csr.vals[i].real(), csr.vals[i].imag()};
	}	
}
void zl_parMatSparse_Show(zl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->show();						
}
void zl_parMatSparse_View(zl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->MatView();						
}
void zl_parMatSparse_Wrt2MMCmplx(zl_parMatSparse_t *pm, char *file_name){
	return reinterpret_cast<parMatrixSparse<std::complex<double>, long > *>(pm)->writeToMatrixMarketCmplx(file_name);						
}


//fcomplex_t int
cs_parMatSparse_t *new_cs_ParMatSparse_1(cs_parVec_t *pv){
	return reinterpret_cast<cs_parMatSparse_t *>( new parMatrixSparse<std::complex<float>,int >(*reinterpret_cast<parVector<std::complex<float>, int > *>(pv)) );
}
cs_parMatSparse_t *new_cs_ParMatSparse_2(parVecMap_t *map){
	return reinterpret_cast<cs_parMatSparse_t *>( new parMatrixSparse<std::complex<float>,int >(*reinterpret_cast<parVectorMap<int > *>(map)) );
}
void cs_parMatSparse_destory(cs_parMatSparse_t *pm){
	delete reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm);
}
int cs_parMatSparse_GetNRows(cs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->GetNRows();
}
int cs_parMatSparse_GetNCols(cs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->GetNCols();	
}
int cs_parMatSparse_GetNNzLoc(cs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->GetNNzLoc();	
}
int cs_parMatSparse_GetLowerBound(cs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->GetLowerBound();		
}
int cs_parMatSparse_GetUpperBound(cs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->GetUpperBound();			
}
MPI_Comm cs_parMatSparse_GetComm(cs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->GetComm();				
}
void cs_parMatSparse_SetValLocal(cs_parMatSparse_t *pm, int row, int col, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->SetValueLocal(row, col, v);					
}
void cs_parMatSparse_SetVal(cs_parMatSparse_t *pm, int row, int col, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->SetValue(row, col, v);						
}
void cs_parMatSparse_AddValLocal(cs_parMatSparse_t *pm, int row, int col, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->AddValueLocal(row, col, v);						
}
void cs_parMatSparse_AddVal(cs_parMatSparse_t *pm, int row, int col, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->AddValue(row, col, v);						
}
fcomplex_t cs_parMatSparse_GetValLocal(cs_parMatSparse_t *pm, int row, int col){
	auto val = reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->GetValueLocal(row, col);	
	return {val.real(), val.imag()};								
}
fcomplex_t cs_parMatSparse_GetVal(cs_parMatSparse_t *pm, int row, int col){
	auto val = reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->GetValue(row, col);	
	return {val.real(), val.imag()};																
}
void cs_parMatSparse_toCSRArray(cs_parMatSparse_t *pm, int *nrows, int *nnz, int **roffs, int **cols, fcomplex_t **val){
	MatrixCSR<std::complex<float>, int > csr = reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->ConvertToCSR();
	*nrows = csr.nrows;
	*nnz = csr.nnz;
	*roffs = (int *)malloc( (csr.nrows + 1) * sizeof(int ) );
	*cols = (int *)malloc( csr.nnz * sizeof(int ) );
	*val = (fcomplex_t *)malloc(csr.nnz * sizeof(fcomplex_t));

	memcpy(*roffs, csr.rows.data(), (csr.nrows + 1) * sizeof(int ));
	memcpy(*cols, csr.cols.data(), csr.nnz * sizeof(int ));
	for(auto i = 0; i < csr.nnz; i++){
		*(*val+i) = {csr.vals[i].real(), csr.vals[i].imag()};
	}	
}
void cs_parMatSparse_Show(cs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->show();						
}
void cs_parMatSparse_View(cs_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->MatView();						
}
void cs_parMatSparse_Wrt2MMCmplx(cs_parMatSparse_t *pm, char *file_name){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, int > *>(pm)->writeToMatrixMarketCmplx(file_name);						
}

//fcomplex_t long
cl_parMatSparse_t *new_cl_ParMatSparse_1(cl_parVec_t *pv){
	return reinterpret_cast<cl_parMatSparse_t *>( new parMatrixSparse<std::complex<float>,long >(*reinterpret_cast<parVector<std::complex<float>, long > *>(pv)) );
}
cl_parMatSparse_t *new_cl_ParMatSparse_2(parVecMap_t *map){
	return reinterpret_cast<cl_parMatSparse_t *>( new parMatrixSparse<std::complex<float>,long >(*reinterpret_cast<parVectorMap<long > *>(map)) );
}
void cl_parMatSparse_destory(cl_parMatSparse_t *pm){
	delete reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm);
}
long cl_parMatSparse_GetNRows(cl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->GetNRows();
}
long cl_parMatSparse_GetNCols(cl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->GetNCols();	
}
long cl_parMatSparse_GetNNzLoc(cl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->GetNNzLoc();	
}
long cl_parMatSparse_GetLowerBound(cl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->GetLowerBound();		
}
long cl_parMatSparse_GetUpperBound(cl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->GetUpperBound();			
}
MPI_Comm cl_parMatSparse_GetComm(cl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->GetComm();				
}
void cl_parMatSparse_SetValLocal(cl_parMatSparse_t *pm, long row, long col, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->SetValueLocal(row, col, v);					
}
void cl_parMatSparse_SetVal(cl_parMatSparse_t *pm, long row, long col, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->SetValue(row, col, v);						
}
void cl_parMatSparse_AddValLocal(cl_parMatSparse_t *pm, long row, long col, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->AddValueLocal(row, col, v);						
}
void cl_parMatSparse_AddVal(cl_parMatSparse_t *pm, long row, long col, fcomplex_t val){
	std::complex<float> v(val.real, val.imag);	
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->AddValue(row, col, v);						
}
fcomplex_t cl_parMatSparse_GetValLocal(cl_parMatSparse_t *pm, long row, long col){
	auto val = reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->GetValueLocal(row, col);	
	return {val.real(), val.imag()};								
}
fcomplex_t cl_parMatSparse_GetVal(cl_parMatSparse_t *pm, long row, long col){
	auto val = reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->GetValue(row, col);	
	return {val.real(), val.imag()};																
}
void cl_parMatSparse_toCSRArray(cl_parMatSparse_t *pm, long *nrows, long *nnz, long **roffs, long **cols, fcomplex_t **val){
	MatrixCSR<std::complex<float>, long > csr = reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->ConvertToCSR();
	*nrows = csr.nrows;
	*nnz = csr.nnz;
	*roffs = (long *)malloc( (csr.nrows + 1) * sizeof(long ) );
	*cols = (long *)malloc( csr.nnz * sizeof(long ) );
	*val = (fcomplex_t *)malloc(csr.nnz * sizeof(fcomplex_t));

	memcpy(*roffs, csr.rows.data(), (csr.nrows + 1) * sizeof(long ));
	memcpy(*cols, csr.cols.data(), csr.nnz * sizeof(long ));
	for(auto i = 0; i < csr.nnz; i++){
		*(*val+i) = {csr.vals[i].real(), csr.vals[i].imag()};
	}	
}
void cl_parMatSparse_Show(cl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->show();						
}
void cl_parMatSparse_View(cl_parMatSparse_t *pm){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->MatView();						
}
void cl_parMatSparse_Wrt2MMCmplx(cl_parMatSparse_t *pm, char *file_name){
	return reinterpret_cast<parMatrixSparse<std::complex<float>, long > *>(pm)->writeToMatrixMarketCmplx(file_name);						
}

//
cs_parVec_t *cs_specNonHermLoad(parVecMap_t *map, char *spectrum){
	cs_parVec_t *pv = new_cs_ParVec_2(map);
	cs_parVecReadExtVec(pv, spectrum);
	return pv;
}
cl_parVec_t *cl_specNonHermLoad(parVecMap_t *map, char *spectrum){
	cl_parVec_t *pv = new_cl_ParVec_2(map);
	cl_parVecReadExtVec(pv, spectrum);
	return pv;	
}
zs_parVec_t *zs_specNonHermLoad(parVecMap_t *map, char *spectrum){
	zs_parVec_t *pv = new_zs_ParVec_2(map);
	zs_parVecReadExtVec(pv, spectrum);
	return pv;	
}
zl_parVec_t *zl_specNonHermLoad(parVecMap_t *map, char *spectrum){
	zl_parVec_t *pv = new_zl_ParVec_2(map);
	zl_parVecReadExtVec(pv, spectrum);
	return pv;	
}

cs_parVec_t *cs_specNonSymmLoad(parVecMap_t *map, char *spectrum){
	cs_parVec_t *pv = new_cs_ParVec_2(map);
	cs_parVecReadExtVec(pv, spectrum);
	return pv;	
}
cl_parVec_t *cl_specNonSymmLoad(parVecMap_t *map, char *spectrum){
	cl_parVec_t *pv = new_cl_ParVec_2(map);
	cl_parVecReadExtVec(pv, spectrum);
	return pv;	
}
zs_parVec_t *zs_specNonSymmLoad(parVecMap_t *map, char *spectrum){
	zs_parVec_t *pv = new_zs_ParVec_2(map);
	zs_parVecReadExtVec(pv, spectrum);
	return pv;	
}
zl_parVec_t *zl_specNonSymmLoad(parVecMap_t *map, char *spectrum){
	zl_parVec_t *pv = new_zl_ParVec_2(map);
	zl_parVecReadExtVec(pv, spectrum);
	return pv;		
}

ds_parVec_t *ds_specNonSymmLoad(parVecMap_t *map, char *spectrum){
	ds_parVec_t *pv = new_ds_ParVec_2(map);
	ds_parVecReadExtVec(pv, spectrum);
	return pv;	
}
dl_parVec_t *dl_specNonSymmLoad(parVecMap_t *map, char *spectrum){
	dl_parVec_t *pv = new_dl_ParVec_2(map);
	dl_parVecReadExtVec(pv, spectrum);
	return pv;	
}
ss_parVec_t *ss_specNonSymmLoad(parVecMap_t *map, char *spectrum){
	ss_parVec_t *pv = new_ss_ParVec_2(map);
	ss_parVecReadExtVec(pv, spectrum);
	return pv;	
}
sl_parVec_t *sl_specNonSymmLoad(parVecMap_t *map, char *spectrum){
	sl_parVec_t *pv = new_sl_ParVec_2(map);
	sl_parVecReadExtVec(pv, spectrum);
	return pv;	
}

//
cs_parMatSparse_t *cs_nonherm(int probSize, nilp_t *nilp, initMatrix_t *init, char *spectrum, MPI_Comm comm){
	parMatrixSparse<std::complex<float>,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	std::string spec(spectrum);

	*matobj = nonherm<std::complex<float>, int>(probSize, *nilpobj, *initMobj, spec, comm);
	return reinterpret_cast<cs_parMatSparse_t *>(matobj);
}

cl_parMatSparse_t *cl_nonherm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm){
	parMatrixSparse<std::complex<float>,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	std::string spec(spectrum);

	*matobj = nonherm<std::complex<float>, long>(probSize, *nilpobj, *initMobj, spec, comm);
	return reinterpret_cast<cl_parMatSparse_t *>(matobj);	
}

zs_parMatSparse_t *zs_nonherm(int probSize, nilp_t *nilp, initMatrix_t *init, char *spectrum, MPI_Comm comm){
	parMatrixSparse<std::complex<double>,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	std::string spec(spectrum);

	*matobj = nonherm<std::complex<double>, int>(probSize, *nilpobj, *initMobj, spec, comm);
	return reinterpret_cast<zs_parMatSparse_t *>(matobj);
}

zl_parMatSparse_t *zl_nonherm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm){
	parMatrixSparse<std::complex<double>,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	std::string spec(spectrum);

	*matobj = nonherm<std::complex<double>, long>(probSize, *nilpobj, *initMobj, spec, comm);
	return reinterpret_cast<zl_parMatSparse_t *>(matobj);		
}


cs_parMatSparse_t *cs_nonherm_2(int probSize, nilp_t *nilp, initMatrix_t *init, cs_parVec_t *spec){
	parMatrixSparse<std::complex<float>,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	parVector<std::complex<float>, int> *spectrum = reinterpret_cast<parVector<std::complex<float>, int> *>(spec);
	*matobj = nonherm<std::complex<float>, int>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<cs_parMatSparse_t *>(matobj);	
}
cl_parMatSparse_t *cl_nonherm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, cl_parVec_t *spec){
	parMatrixSparse<std::complex<float>,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	parVector<std::complex<float>, long> *spectrum = reinterpret_cast<parVector<std::complex<float>, long> *>(spec);
	*matobj = nonherm<std::complex<float>, long>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<cl_parMatSparse_t *>(matobj);		
}
zs_parMatSparse_t *zs_nonherm_2(int probSize, nilp_t *nilp, initMatrix_t *init, zs_parVec_t *spec){
	parMatrixSparse<std::complex<double>,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	parVector<std::complex<double>, int> *spectrum = reinterpret_cast<parVector<std::complex<double>, int> *>(spec);
	*matobj = nonherm<std::complex<double>, int>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<zs_parMatSparse_t *>(matobj);	
}
zl_parMatSparse_t *zl_nonherm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, zl_parVec_t *spec){
	parMatrixSparse<std::complex<double>,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	parVector<std::complex<double>, long> *spectrum = reinterpret_cast<parVector<std::complex<double>, long> *>(spec);
	*matobj = nonherm<std::complex<double>, long>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<zl_parMatSparse_t *>(matobj);		
}

void cs_nonherm_3(int probSize, nilp_t *nilp, cs_parMatSparse_t *Am, cs_parVec_t *spec){
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	parVector<std::complex<float>, int> *spectrum = reinterpret_cast<parVector<std::complex<float>, int> *>(spec);
	parMatrixSparse<std::complex<float>, int> *Mat = reinterpret_cast<parMatrixSparse<std::complex<float>, int> *>(Am);
	nonherm<std::complex<float>, int>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<cs_parMatSparse_t *>(Mat);
}

void cl_nonherm_3(long probSize, nilpL_t *nilp, cl_parMatSparse_t *Am, cl_parVec_t *spec){
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	parVector<std::complex<float>, long> *spectrum = reinterpret_cast<parVector<std::complex<float>, long> *>(spec);
	parMatrixSparse<std::complex<float>, long> *Mat = reinterpret_cast<parMatrixSparse<std::complex<float>, long> *>(Am);
	nonherm<std::complex<float>, long>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<cl_parMatSparse_t *>(Mat);	
}
void zs_nonherm_3(int probSize, nilp_t *nilp, zs_parMatSparse_t *Am, zs_parVec_t *spec){
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	parVector<std::complex<double>, int> *spectrum = reinterpret_cast<parVector<std::complex<double>, int> *>(spec);
	parMatrixSparse<std::complex<double>, int> *Mat = reinterpret_cast<parMatrixSparse<std::complex<double>, int> *>(Am);
	nonherm<std::complex<double>, int>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<zs_parMatSparse_t *>(Mat);	
}
void zl_nonherm_3(long probSize, nilpL_t *nilp, zl_parMatSparse_t *Am, zl_parVec_t *spec){
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	parVector<std::complex<double>, long> *spectrum = reinterpret_cast<parVector<std::complex<double>, long> *>(spec);
	parMatrixSparse<std::complex<double>, long> *Mat = reinterpret_cast<parMatrixSparse<std::complex<double>, long> *>(Am);
	nonherm<std::complex<double>, long>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<zl_parMatSparse_t *>(Mat);		
}

//
ss_parMatSparse_t *ss_nonsymm(int probSize, nilp_t *nilp, initMatrix_t *init, char *spectrum, MPI_Comm comm){
	parMatrixSparse<float,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	std::string spec(spectrum);

	*matobj = nonsymm<float, int>(probSize, *nilpobj, *initMobj, spec, comm);
	return reinterpret_cast<ss_parMatSparse_t *>(matobj);	
}
sl_parMatSparse_t *sl_nonsymm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm){
	parMatrixSparse<float,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	std::string spec(spectrum);

	*matobj = nonsymm<float, long>(probSize, *nilpobj, *initMobj, spec, comm);
	return reinterpret_cast<sl_parMatSparse_t *>(matobj);		
}
ds_parMatSparse_t *ds_nonsymm(int probSize, nilp_t *nilp, initMatrix_t *init,  char *spectrum, MPI_Comm comm){
	parMatrixSparse<double,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	std::string spec(spectrum);

	*matobj = nonsymm<double, int>(probSize, *nilpobj, *initMobj, spec, comm);
	return reinterpret_cast<ds_parMatSparse_t *>(matobj);	
}
dl_parMatSparse_t *dl_nonsymm(long probSize, nilpL_t *nilp, initMatrixL_t *init, char *spectrum, MPI_Comm comm){
	parMatrixSparse<double,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	std::string spec(spectrum);

	*matobj = nonsymm<double, long>(probSize, *nilpobj, *initMobj, spec, comm);
	return reinterpret_cast<dl_parMatSparse_t *>(matobj);		
}

//
ss_parMatSparse_t *ss_nonsymm_2(int probSize, nilp_t *nilp, initMatrix_t *init, ss_parVec_t *spec){
	parMatrixSparse<float,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	parVector<float, int> *spectrum = reinterpret_cast<parVector<float, int> *>(spec);
	*matobj = nonsymm<float, int>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<ss_parMatSparse_t *>(matobj);		
}
sl_parMatSparse_t *sl_nonsymm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, sl_parVec_t *spec){
	parMatrixSparse<float,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	parVector<float, long> *spectrum = reinterpret_cast<parVector<float, long> *>(spec);
	*matobj = nonsymm<float, long>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<sl_parMatSparse_t *>(matobj);		
}
ds_parMatSparse_t *ds_nonsymm_2(int probSize, nilp_t *nilp, initMatrix_t *init, ds_parVec_t *spec){
{
	parMatrixSparse<double,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	parVector<double, int> *spectrum = reinterpret_cast<parVector<double, int> *>(spec);
	*matobj = nonsymm<double, int>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<ds_parMatSparse_t *>(matobj);		
}	
}
dl_parMatSparse_t *dl_nonsymm_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, dl_parVec_t *spec){
	parMatrixSparse<double,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	parVector<double, long> *spectrum = reinterpret_cast<parVector<double, long> *>(spec);
	*matobj = nonsymm<double, long>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<dl_parMatSparse_t *>(matobj);		
}

void ss_nonsymm_3(int probSize, nilp_t *nilp, ss_parMatSparse_t *Am, ss_parVec_t *spec){
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	parVector<float, int> *spectrum = reinterpret_cast<parVector<float, int> *>(spec);
	parMatrixSparse<float, int> *Mat = reinterpret_cast<parMatrixSparse<float, int> *>(Am);
	nonsymm<float, int>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<ss_parMatSparse_t *>(Mat);
}
void sl_nonsymm_3(int probSize, nilpL_t *nilp, sl_parMatSparse_t *Am, sl_parVec_t *spec){
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	parVector<float, long> *spectrum = reinterpret_cast<parVector<float, long> *>(spec);
	parMatrixSparse<float, long> *Mat = reinterpret_cast<parMatrixSparse<float, long> *>(Am);
	nonsymm<float, long>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<sl_parMatSparse_t *>(Mat);
}
void ds_nonsymm_3(int probSize, nilp_t *nilp, ds_parMatSparse_t *Am, ds_parVec_t *spec){
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	parVector<double, int> *spectrum = reinterpret_cast<parVector<double, int> *>(spec);
	parMatrixSparse<double, int> *Mat = reinterpret_cast<parMatrixSparse<double, int> *>(Am);
	nonsymm<double, int>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<ds_parMatSparse_t *>(Mat);	
}

void dl_nonsymm_3(int probSize, nilpL_t *nilp, dl_parMatSparse_t *Am, dl_parVec_t *spec){
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	parVector<double, long> *spectrum = reinterpret_cast<parVector<double, long> *>(spec);
	parMatrixSparse<double, long> *Mat = reinterpret_cast<parMatrixSparse<double, long> *>(Am);
	nonsymm<double, long>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<dl_parMatSparse_t *>(Mat);
}

//
ss_parMatSparse_t *ss_nonsymmconj_2(int probSize, nilp_t *nilp, initMatrix_t *init, cs_parVec_t *spec){
	parMatrixSparse<float,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	parVector<std::complex<float>, int> *spectrum = reinterpret_cast<parVector<std::complex<float>, int> *>(spec);
	*matobj = nonsymmconj<float, int>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<ss_parMatSparse_t *>(matobj);		
}
sl_parMatSparse_t *sl_nonsymmconj_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, cl_parVec_t *spec){
	parMatrixSparse<float,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	parVector<std::complex<float>, long> *spectrum = reinterpret_cast<parVector<std::complex<float>, long> *>(spec);
	*matobj = nonsymmconj<float, long>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<sl_parMatSparse_t *>(matobj);	
}
ds_parMatSparse_t *ds_nonsymmconj_2(int probSize, nilp_t *nilp, initMatrix_t *init, zs_parVec_t *spec){
	parMatrixSparse<double,int> *matobj;
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	initMat<int> *initMobj = reinterpret_cast<initMat<int> *>(init);
	parVector<std::complex<double>, int> *spectrum = reinterpret_cast<parVector<std::complex<double>, int> *>(spec);
	*matobj = nonsymmconj<double, int>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<ds_parMatSparse_t *>(matobj);		
}
dl_parMatSparse_t *dl_nonsymmconj_2(long probSize, nilpL_t *nilp, initMatrixL_t *init, zl_parVec_t *spec){
	parMatrixSparse<double,long> *matobj;
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	initMat<long> *initMobj = reinterpret_cast<initMat<long> *>(init);
	parVector<std::complex<double>, long> *spectrum = reinterpret_cast<parVector<std::complex<double>, long> *>(spec);
	*matobj = nonsymmconj<double, long>(probSize, *nilpobj, *initMobj, *spectrum);
	return reinterpret_cast<dl_parMatSparse_t *>(matobj);	
}

void ss_nonsymmconj_3(int probSize, nilp_t *nilp, ss_parMatSparse_t *Am, cs_parVec_t *spec){
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	parVector<std::complex<float>, int> *spectrum = reinterpret_cast<parVector<std::complex<float>, int> *>(spec);
	parMatrixSparse<float, int> *Mat = reinterpret_cast<parMatrixSparse<float, int> *>(Am);
	nonsymmconj<float, int>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<ss_parMatSparse_t *>(Mat);
}
void sl_nonsymmconj_3(long probSize, nilpL_t *nilp, sl_parMatSparse_t *Am, cl_parVec_t *spec){
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	parVector<std::complex<float>, long> *spectrum = reinterpret_cast<parVector<std::complex<float>, long> *>(spec);
	parMatrixSparse<float, long> *Mat = reinterpret_cast<parMatrixSparse<float, long> *>(Am);
	nonsymmconj<float, long>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<sl_parMatSparse_t *>(Mat);
}
void ds_nonsymmconj_3(int probSize, nilp_t *nilp, ds_parMatSparse_t *Am, zs_parVec_t *spec){
	Nilpotent<int> *nilpobj = reinterpret_cast<Nilpotent<int> *>(nilp);
	parVector<std::complex<double>, int> *spectrum = reinterpret_cast<parVector<std::complex<double>, int> *>(spec);
	parMatrixSparse<double, int> *Mat = reinterpret_cast<parMatrixSparse<double, int> *>(Am);
	nonsymmconj<double, int>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<ds_parMatSparse_t *>(Mat);	
}
void dl_nonsymmconj_3(long probSize, nilpL_t *nilp, dl_parMatSparse_t *Am, zl_parVec_t *spec){
	Nilpotent<long> *nilpobj = reinterpret_cast<Nilpotent<long> *>(nilp);
	parVector<std::complex<double>, long> *spectrum = reinterpret_cast<parVector<std::complex<double>, long> *>(spec);
	parMatrixSparse<double, long> *Mat = reinterpret_cast<parMatrixSparse<double, long> *>(Am);
	nonsymmconj<double, long>(probSize, *nilpobj, Mat ,*spectrum);
	Am = reinterpret_cast<dl_parMatSparse_t *>(Mat);
}


