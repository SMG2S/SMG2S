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


#ifndef __PAR_MATRIX_SPARSE_H__
#define __PAR_MATRIX_SPARSE_H__

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include "../parVector/parVector.h"
#include "../nilpotent/nilpotent.h"
#include "../utils/MPI_DataType.h"
#include "../utils/utils.h"

#include "MatrixCSR.h"

#ifdef __USE_COMPLEX__
#include <complex>
#endif

template<typename T, typename S>
class parMatrixSparse
{
    private:
	std::map<S, T> *dynmat_loc;

	//size of local matrix
	S	ncols, nrows;
	//number of non-zeros in local
	S	nnz_loc;

	parVectorMap<S>	index_map;

	S	lower_b, upper_b;

	MPI_Comm comm;

	// mpi size and rank
	int ProcID, nProcs;

    public:
    	parMatrixSparse();
    	parMatrixSparse(parVector<T,S> vec);
    	parMatrixSparse(parVectorMap<S> map);

    	//get
    	S GetNRows(){return nrows;};
    	S GetNCols(){return ncols;};
    	S GetNNzLoc(){return nnz_loc;};
    	parVectorMap<S> GetMap(){return index_map;};
    	S GetLowerBound(){return lower_b;};
    	S GetUpperBound(){return upper_b;};
    	MPI_Comm GetComm(){return comm;};
    	int GetProcId(){return ProcID;};
    	int GetNProcs(){return nProcs;};

	std::map<S, T> *GetDynMatLoc(){return dynmat_loc;};

	//set value
	void	SetValueLocal(S row, S col, T value);
	void	SetValuesLocal(S nindex, S *rows, S *cols, T *values);
	void 	SetValue(S row, S col, T value);
	void    SetDiagonal(parVector<T,S> diag);

	//add value
	void    AddValueLocal(S row, S col, T value);
	void 	AddValue(S row, S col, T value);

	//get value
	T 	GetValueLocal(S row, S col);
	T 	GetValue(S row, S col);

	void	MatScale(T scale);
	void    MatAXPY(parMatrixSparse<T,S> X, T scale);
	void    MatAYPX(parMatrixSparse<T,S> X, T scale);

	//init matrix lower part of matrix 
	void	initMat(S diag_l, S diag_u, Base<T> scale, T shift, Base<T> sparsity );
	void	initMat(S diag_l, S diag_u);
	
	void setSpecNonHerm(parVector<T,S> spectrum);
	void setSpecNonSymm(parVector<T,S> spectrum);
	void setSpecNonSymmCmplx(parVector<std::complex<Base<T>>,S> spectrum);

	void    updateNnz();

	void    copy(parMatrixSparse<T,S> X);

	void rmNNz();

   	//matrix multiplies a nilpotent matrix
	parMatrixSparse<T,S>	MA(Nilpotent<S> nilp);

   	//a nilpotent matrix multiplies a matrix
	parMatrixSparse<T,S>	AM(Nilpotent<S> nilp);


	void	ZeroEntries();
	MatrixCSR<T,S>	ConvertToCSR();

    	void show();
    	void MatView();

    	void writeToMatrixMarket(std::string file_name);
    	void writeToMatrixMarketCmplx(std::string file_name);

};

template<typename T, typename S>
parMatrixSparse<T,S>::parMatrixSparse(parVector<T,S> vec)
{
    index_map = vec.GetVecMap();
    lower_b = index_map.GetLowerBound();
    upper_b = index_map.GetUpperBound();
    comm = index_map.GetCurrentComm();
    MPI_Comm_rank(comm, &ProcID);
    MPI_Comm_size(comm, &nProcs);
    nrows = index_map.GetLocalSize();
    ncols = index_map.GetGlobalSize();

    //at the beginning, empty matrix, and elements will be dynamically added
    nnz_loc = 0;
    dynmat_loc = NULL;
}

template<typename T, typename S>
parMatrixSparse<T,S>::parMatrixSparse(parVectorMap<S> map)
{
    index_map = map;
    lower_b = index_map.GetLowerBound();
    upper_b = index_map.GetUpperBound();
    comm = index_map.GetCurrentComm();
    MPI_Comm_rank(comm, &ProcID);
    MPI_Comm_size(comm, &nProcs);
    nrows = index_map.GetLocalSize();
    ncols = index_map.GetGlobalSize();

    //at the beginning, empty matrix, and elements will be dynamically added
    nnz_loc = 0;
    dynmat_loc = NULL;
}

template<typename T, typename S>
void parMatrixSparse<T,S>::show()
{
    std::cout << "Id " << ProcID << ": local matrix size = " << nrows << "x" << ncols << ", nnz = " << nnz_loc << ", lower_b = "  << lower_b << ", upper_b = " << upper_b << std::endl;
}

template<typename T, typename S>
void parMatrixSparse<T,S>::MatView()
{
    typename std::map<S,T>::iterator it;
    if(ProcID == 0) {std::cout << "Parallel MatView: " << std::endl;}

    for (S i = 0; i < nrows; i++){
    	if(dynmat_loc != NULL){
    	    std::cout << "row " << index_map.Loc2Glob(i) << ": ";
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
		std::cout <<"("<<it->first << "," << it->second << "); ";	
	    }
	std::cout << std::endl;    
	}
    }

}

template<typename T,typename S>
void parMatrixSparse<T,S>::SetValueLocal(S row, S col, T value)
{
    typename std::map<S,T>::iterator it;

    if(dynmat_loc == NULL){
    	dynmat_loc = new std::map<S,T> [nrows];
    }
    it = dynmat_loc[row].find(col);

    if(value != T(0)){
    	if(it == dynmat_loc[row].end()){
    	    dynmat_loc[row][col] = value;
    	    nnz_loc++;
    	}
    	else{
    	    it->second = value;
        }    
    }
}


template<typename T,typename S>
void parMatrixSparse<T,S>::SetValuesLocal(S nindex, S *rows, S *cols, T *values)
{
    for( S i = 0; i < nindex; i++){
    	SetValueLocal(rows[i],cols[i],values[i]);
    }
}

template<typename T,typename S>
void parMatrixSparse<T,S>::SetValue(S row, S col, T value)
{
    S local_row = index_map.Glob2Loc(row);

    if(local_row >= 0 && local_row < nrows){
    	SetValueLocal(local_row, col, value);
    }
}


template<typename T,typename S>
void parMatrixSparse<T,S>::SetDiagonal(parVector<T,S> diag)
{
    T *a = diag.GetArray();
    S local_size = diag.GetLocalSize();

    for(auto i = 0; i < local_size; i++){
    	SetValueLocal(i, diag.Loc2Glob(i), a[i]);
    }
}


template<typename T,typename S>
void parMatrixSparse<T,S>::AddValueLocal(S row, S col, T value)
{
    typename std::map<S,T>::iterator it;

    if(dynmat_loc == NULL){
    	dynmat_loc = new std::map<S,T> [nrows];
    }
    it = dynmat_loc[row].find(col);

    if(it == dynmat_loc[row].end()){
    	dynmat_loc[row][col] = value;
    	nnz_loc++;
    }
    else{
    	it->second += value;
    }    
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AddValue(S row, S col, T value)
{
    S local_row = index_map.Glob2Loc(row);

    if(local_row >= 0 && local_row < nrows){
    	AddValueLocal(local_row, col, value);
    }
}


template<typename T,typename S>
T parMatrixSparse<T,S>::GetValueLocal(S row, S col)
{
    typename std::map<S,T>::iterator it;

    if(dynmat_loc != NULL){
    	it = dynmat_loc[row].find(col);
    	if(it == dynmat_loc[row].end()){
    	    return T(0);
    	}else{
    	    return dynmat_loc[row][col];
    	}
    }else{
    	return T(0);
    }    
}

template<typename T,typename S>
T parMatrixSparse<T,S>::GetValue(S row, S col)
{
    S local_row = index_map.Glob2Loc(row);
    if(local_row >= 0 && local_row < nrows){
    	return GetLocalValue(local_row, col);
    }
    else{
    	return 0.0;
    }
}

template<typename T,typename S>
void parMatrixSparse<T,S>::MatScale(T scale){
    typename std::map<S,T>::iterator it;
    S i;
    if(dynmat_loc != NULL){
    	for(i = 0; i < nrows; i++){
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
		it->second = it->second*scale;
	    }
	}
    }
}


template<typename T,typename S>
void parMatrixSparse<T,S>::MatAXPY(parMatrixSparse<T,S> X, T scale){
    typename std::map<S,T>::iterator it;
    typename std::map<S,T>::iterator it1;

    S i, k;

    if(index_map != X.GetMap()){
    	try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for AXPY, matrix should be in the same vectorMap" << std::endl;
    	}
    }

    if(X.GetDynMatLoc() != NULL && scale != T(1)){
    	for(i = 0; i < nrows; i++){
    	    for(it1 = X.GetDynMatLoc()[i].begin(); it1 != X.GetDynMatLoc()[i].end();++it1){
    	    	X.GetDynMatLoc()[i][it1->first] = scale * it1->second;
    	    }   
    	    	
    	}
    }
    if(dynmat_loc != NULL && X.GetDynMatLoc() != NULL){
    	for(i = 0; i < nrows; i++){

    	    dynmat_loc[i] = std::accumulate( X.GetDynMatLoc()[i].begin(), X.GetDynMatLoc()[i].end(), dynmat_loc[i],
        		[]( std::map<S, T> &m, const std::pair<const S, T> &p)
        			{
            			    return ( m[p.first] += p.second, m );
        			} );

    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end();){
    	    	k = it->first;
    	    	if(it->second == T(0)){
    	    	     it=dynmat_loc[i].erase(it);
    	    	}else{
    	    	     it++;
    	    	}
	    }

    	}

    	updateNnz();

    } else if(dynmat_loc == NULL && X.GetDynMatLoc() != NULL){
    	dynmat_loc = new std::map<S,T> [nrows];
    	for(i = 0; i < nrows; i++){
    	    dynmat_loc[i].insert(X.GetDynMatLoc()[i].begin(), X.GetDynMatLoc()[i].end());		
	}
    	nnz_loc = X.nnz_loc;
    }
}


template<typename T,typename S>
void parMatrixSparse<T,S>::MatAYPX(parMatrixSparse<T,S> X, T scale){
    typename std::map<S,T>::iterator it;

    S i, k;
    if(dynmat_loc != NULL && scale != T(1)){
    	for(i = 0; i < nrows; i++){
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
    	    	k = it->first;
    	    	dynmat_loc[i][k] = dynmat_loc[i][k]*scale;
	    }
	}
    }

    if(dynmat_loc != NULL && X.GetDynMatLoc() != NULL){
    	for(i = 0; i < nrows; i++){

    	    dynmat_loc[i] = std::accumulate( X.GetDynMatLoc()[i].begin(), X.GetDynMatLoc()[i].end(), dynmat_loc[i],
        		[]( std::map<S, T> &m, const std::pair<const S, T> &p)
        			{
            			    return ( m[p.first] += p.second, m );
        			} );

    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end();){
    	    	k = it->first;
    	    	if(it->second == T(0)){
    	    	     it=dynmat_loc[i].erase(it);
    	    	}else{
    	    	     it++;
    	    	}
	    }
    	}

    	updateNnz();
    }

    if(dynmat_loc == NULL && X.GetDynMatLoc() != NULL){
    	dynmat_loc = new std::map<S,T> [nrows];
    	for(i = 0; i < nrows; i++){
    	    dynmat_loc[i].insert(X.GetDynMatLoc()[i].begin(), X.GetDynMatLoc()[i].end() );		
	}
	nnz_loc = X.nnz_loc;
    }
}

template<typename T, typename S>
void parMatrixSparse<T,S>::copy(parMatrixSparse<T,S> X){
    
    dynmat_loc = new std::map<S,T> [nrows];
    
    typename std::map<S,T>::iterator it;

    S i, k;

    if(X.GetDynMatLoc() != NULL){
    	for(i = 0; i < nrows; i++){
    	    dynmat_loc[i].insert(X.GetDynMatLoc()[i].begin(),X.GetDynMatLoc()[i].end());
	}
    }
}


template<typename T, typename S>
void parMatrixSparse<T,S>::updateNnz(){
        
    typename std::map<S,T>::iterator it;

    S i, k;

    nnz_loc = 0;

    if(dynmat_loc != NULL){
    	for(i = 0; i < nrows; i++){
    	    nnz_loc += dynmat_loc[i].size();
	}
    }
}


template<typename T, typename S>
void parMatrixSparse<T,S>::rmNNz(){
    typename std::map<S,T>::iterator it;
    S i, k;
    if(dynmat_loc != NULL){
    	for(i = 0; i < nrows; i++){
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end();){
    	    	k = it->first;
    	    	if(it->second == T(0)){
    	    	     it=dynmat_loc[i].erase(it);
    	    	}else{
    	    	     it++;
    	    	}
	    }
	}
    
    	updateNnz();
    }

}

template<typename T,typename S>
void parMatrixSparse<T,S>::ZeroEntries()
{
    typename std::map<S,T>::iterator it;
    S i;
    if(dynmat_loc != NULL){
	for(i = 0; i < nrows; i++){
	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
		it->second = T(0);
	    }
	}
    }
}

template<typename T, typename S>
void parMatrixSparse<T,S>::initMat(S diag_l, S diag_u, Base<T> scale, T shift, Base<T> sparsity){
   
    bool isnonsymm = false;

    if (sizeof(scale)/sizeof(sparsity) == 1){
        //std::cout << "using real scalar" << std::endl;
    	isnonsymm = true;
    }

    S size = ncols;

    diag_l = abs(diag_l);
    diag_u = abs(diag_u);


    if (isnonsymm && diag_u < 2){
    	try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for initialisationg of non-symmetric matrix, please ensure abs(diag_u) >= 2" << std::endl;
    	}
    }

    if (diag_l < diag_u){
    	try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for initialisationg of matrix, please ensure abs(diag_l) < abs(diag_u)" << std::endl;
    	}    	
    }
                  
    if (diag_u >= size){
    	try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for initialisationg of matrix, please ensure abs(diag_u) < size" << std::endl;
    	}     	
    }

    //std::random_device rd;
    std::mt19937_64 rd(ProcID);
    std::uniform_real_distribution<> d(0, 1);
            
    for (auto i = std::max(diag_u, lower_b); i < upper_b; i++){
    	for (auto j = std::max(0, i - diag_l); j <= i - diag_u; j++){
    	    auto sampling = d(rd);
    	    if (sampling < 1.0 - sparsity){
    	    	SetValue(i, j, scale * d(rd) + shift);
    	    }
    	}
    }

}

template<typename T, typename S>
void parMatrixSparse<T,S>::initMat(S diag_l, S diag_u){
    initMat(diag_l, diag_u, 1.0, 0.0, 0.9);
}


template<typename T, typename S>
void parMatrixSparse<T,S>::setSpecNonHerm(parVector<T,S> spectrum){
    int size1 = sizeof(T) / sizeof(Base<T>);

    if(size1 == 1){
        printf("Info ]> For generating non-Hermtian matrices, the matrices should be real\n");
        return;
    }

    SetDiagonal(spectrum);
}

template<typename T, typename S>
void parMatrixSparse<T,S>::setSpecNonSymm(parVector<T,S> spectrum){
    
    int size1 = sizeof(T) / sizeof(Base<T>);

    if(size1 == 2){
        printf("Info ]> For generating non-Symmetric matrices, the matrices should be real\n");
        return;
    }

    SetDiagonal(spectrum);
}

template<typename T, typename S>
void parMatrixSparse<T,S>::setSpecNonSymmCmplx(parVector<std::complex<Base<T>>,S> spectrum){
    
    int size1 = sizeof(T) / sizeof(Base<T>);

    if(size1 == 2){
        printf("Info ]> For generating non-Symmetric matrices, the matrices should be real\n");
        return;
    }

    auto overlap = checkNonSymmSpec<T,S>(spectrum);

    auto *array = spectrum.GetArray();

    auto local_size = spectrum.GetLocalSize();
    auto nu = spectrum.GetUpperNeighbor(1);
    auto nl = spectrum.GetLowerNeighbor(1);

    S idx = 0;
    S step;
    S nzeros;
    if(overlap == 1){
    	idx = 1;
    	SetValueLocal(0, lower_b, array[0].real());
    	SetValueLocal(0, lower_b-1, -abs(array[0].imag()));
    }else if(overlap == 0){
    	idx = 0;
    }

    nzeros = 0;

    while(idx < (local_size - 1) ){
    	if(array[idx].imag() == 0){
    	    SetValueLocal(idx, idx + lower_b, array[idx].real());
    	    step = 1;
    	    nzeros++;
    	}else{
    	    SetValueLocal(idx, idx + lower_b, array[idx].real());
    	    SetValueLocal(idx, idx + lower_b + 1, abs(array[idx].imag()) );
    	    SetValueLocal(idx + 1, idx + lower_b + 1, array[idx+1].real());
    	    SetValueLocal(idx + 1, idx + lower_b , -abs(array[idx+1].imag()) );    	    
    	    step = 2;
    	}
    	idx += step;
    }

    if( (local_size - nzeros - overlap) % 2 != 0 ){
    	if(array[local_size - 1].imag() != 0){
    	    SetValueLocal(local_size - 1, lower_b+local_size-1, array[local_size - 1].real());
    	    SetValueLocal(local_size - 1, lower_b+local_size, abs(array[local_size - 1].imag()));
    	}else{
    	    SetValueLocal(local_size - 1, local_size - 1 + lower_b, array[local_size - 1].real());
    	}
    }
}

template<typename T,typename S>
MatrixCSR<T,S> parMatrixSparse<T,S>::ConvertToCSR()
{
    S 	count, i, j;
    T	v;
    typename std::map<S, T>::iterator it;

    MatrixCSR<T,S> csr = MatrixCSR<T,S>(nnz_loc, nrows);

    if(dynmat_loc != NULL){
	count = 0;
	csr.rows.push_back(count);
	for(i = 0; i < nrows; i++){	    
	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
	    	if(it->second != T(0)){
	    	    j = it->first;
	    	    v = it->second;
	    	    csr.vals.push_back(v);
	    	    csr.cols.push_back(j);
	    	    count++;
		}
	    }
	    csr.rows.push_back(count);
	}
    }
    return csr;
}

//matrix multiplies a nilpotent matrix
template<typename T,typename S>
parMatrixSparse<T,S> parMatrixSparse<T,S>::MA(Nilpotent<S> nilp){

    auto offset = nilp.getOffset();
    auto IndOfZeros = nilp.getIndOfZeros();
    auto prod = parMatrixSparse<T, S>(index_map);

    typename std::map<S,T>::iterator it;

    if(dynmat_loc != NULL){
    	for(auto i = 0; i < nrows; i++){
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
    	    	S j =  it->first + offset;
    	    	if(j < ncols){
    	    	    //if not the index of zeros in nilpotent
    	    	    if (std::find(IndOfZeros.begin(), IndOfZeros.end(), j-offset) == IndOfZeros.end() && it->second != T(0) ){
    	    	    	prod.SetValueLocal(i, j, it->second);
    	    	    }
    	    	}
    	    }	
    	}
    }
    return prod;
}

//a nilpotent matrix multiplies a matrix
template<typename T,typename S>
parMatrixSparse<T,S> parMatrixSparse<T,S>::AM(Nilpotent<S> nilp){
    auto offset = nilp.getOffset();
    auto IndOfZeros = nilp.getIndOfZeros();
    auto prod = parMatrixSparse<T, S>(index_map);

    auto lprocbound_map = index_map.GetLBoundMap();
    auto uprocbound_map = index_map.GetUBoundMap();

    std::vector<S> targetProcs(ncols);
    std::vector<S> originProcs(ncols);
    std::vector<S> procsDiff(ncols);
    
    typename std::map<S,T>::iterator it;

    if(prod.dynmat_loc == NULL){
    	prod.dynmat_loc = new std::map<S,T> [nrows];
    }
    

    //shift the in-memory part
    if(dynmat_loc != NULL){
    	for(auto row = 0; row < nrows; row++){
    	    for(it = dynmat_loc[row].begin(); it != dynmat_loc[row].end(); ++it){
    	    	auto i = row - offset;
    	    	auto j =  it->first;
    	    	if(i >= 0){
    	    	    //if not the index of zeros in nilpotent
    	    	    if (std::find(IndOfZeros.begin(), IndOfZeros.end(), index_map.Loc2Glob(row)-offset) == IndOfZeros.end()&& it->second != T(0) ){
    	    	    	prod.SetValueLocal(i, j, it->second);
    	    	    }
    	    	}
    	    }	
    	}
    }


    //communication part
    std::vector<std::vector<T>> sBufs;
    std::vector<std::vector<S>> sIndices; //    sRoffsets + sColIndx
    std::vector<std::vector<S>> sSize;


    std::vector<std::vector<T>> rBufs;
    std::vector<std::vector<S>> rIndices; //    rRoffsets + rColIndx
    std::vector<std::vector<S>> rSize;   

    for(S i = 0; i < ncols; i++){
    	for(S j = 0; j < lprocbound_map.size(); j++){
    	    if( (i - offset >= lprocbound_map[j]) && (i - offset < uprocbound_map[j]) ){
    	    	targetProcs[i] = j;
    	    	break;
    	    }
    	}
    }

    for(S i = 0; i < ncols; i++){
    	for(S j = 0; j < lprocbound_map.size(); j++){
    	    if( (i  >= lprocbound_map[j]) && (i < uprocbound_map[j]) ){
    	    	originProcs[i] = j;
    	    	break;
    	    }
    	}
    }

    for(S i = 0; i < ncols; i++){
    	procsDiff[i] = originProcs[i] - targetProcs[i];
    }

    auto sendrecv_cnt = distinct<S,S>(procsDiff.data(), procsDiff.size(), S(0));

    sBufs.resize(sendrecv_cnt);
    sIndices.resize(sendrecv_cnt);
    sSize.resize(sendrecv_cnt);

    rBufs.resize(sendrecv_cnt);
    rIndices.resize(sendrecv_cnt);
    rSize.resize(sendrecv_cnt);

    for(auto i = 0; i < sendrecv_cnt; i++){
    	sSize[i].push_back(0);
	sSize[i].push_back(0);
	rSize[i].push_back(0);
	rSize[i].push_back(0);
    }

    MPI_Request	*stypereq = new MPI_Request[sendrecv_cnt];
    MPI_Request *rtypereq = new MPI_Request[sendrecv_cnt];
    MPI_Status	*typestat = new MPI_Status[sendrecv_cnt];

    int *up = new int[sendrecv_cnt];
    int *down = new int[sendrecv_cnt];

    for(auto i = 0; i < sendrecv_cnt; i++){
    	up[i] = ProcID - i - 1;
    	down[i] = ProcID + i + 1;
    }


    if(ProcID != 0 && dynmat_loc != NULL){
    	S count = lprocbound_map[ProcID];
    	for(auto k = 0; k < sendrecv_cnt; k++){
    	    sIndices[k].insert(sIndices[k].begin(), count);
    	    for( S i = lprocbound_map[ProcID]; i < uprocbound_map[ProcID]; i++){
    	    	if( (originProcs[i] - targetProcs[i] - 1) == k ){
    	    	    S j = index_map.Glob2Loc(i);
    	    	    for(it = dynmat_loc[j].begin(); it != dynmat_loc[j].end(); ++it){
    	    	    	sIndices[k].push_back(it->first);
    	    	    	sBufs[k].push_back(it->second);
    	    	    	count ++;
    	    	    }
    	    	    sIndices[k].insert(sIndices[k].begin(), count);
    	    	}
    	    }

    	    sSize[k][0] = sBufs[k].size();
    	    sSize[k][1] = sIndices[k].size() - sBufs[k].size() - 1;
    	}

    }


    for(auto k = 0; k < sendrecv_cnt; k++){
    	if(ProcID != 0 && dynmat_loc != NULL){
    	    MPI_Isend(sSize[k].data(), 2, getMPI_Type<S>(), up[k], k, comm, &stypereq[k]);
    	}
    	if(ProcID != nProcs - 1){
    	    MPI_Irecv(rSize[k].data(), 2, getMPI_Type<S>(), down[k], k, comm, &rtypereq[k]);
    	    MPI_Wait(&rtypereq[k], &typestat[k]);

    	    rBufs[k].resize(rSize[k][0]);
    	    rIndices[k].resize(rSize[k][0]+rSize[k][1]+1);
    	}
    }


    for(auto k = 0; k < sendrecv_cnt; k++){
    	if(ProcID != 0 && dynmat_loc != NULL){
    	    MPI_Isend(sBufs[k].data(), sSize[k][0], getMPI_Type<T>(), up[k], k, comm, &stypereq[k]);
    	}
    	if(ProcID != nProcs - 1){
    	    MPI_Irecv(rBufs[k].data(), rSize[k][0], getMPI_Type<T>(), down[k], k, comm, &rtypereq[k]);	
    	    MPI_Wait(&rtypereq[k], &typestat[k]);
    	}
    }


    for(auto k = 0; k < sendrecv_cnt; k++){
    	if(ProcID != 0 && dynmat_loc != NULL){
    	    MPI_Isend(sIndices[k].data(), sSize[k][0]+sSize[k][1]+1, getMPI_Type<S>(), up[k], k, comm, &stypereq[k]);
    	}
    	if(ProcID != nProcs - 1){
    	    MPI_Irecv(rIndices[k].data(), rSize[k][0]+rSize[k][1]+1, getMPI_Type<S>(), down[k], k, comm, &rtypereq[k]);	
    	    MPI_Wait(&rtypereq[k], &typestat[k]);
    	}
    }

    for(auto k = 0; k < sendrecv_cnt; k++){
    	if(ProcID != nProcs - 1){
    	    for(auto row = 0; row < rSize[k][1] ; row++){
    	        for(auto cnt = rIndices[k][rSize[k][1] - row] - rIndices[k][rSize[k][1]] ; cnt < rIndices[k][rSize[k][1] - row-1] - rIndices[k][rSize[k][1]]; cnt++){
    		    auto col = rIndices[k][cnt + rSize[k][1] + 1];
    		    auto val = rBufs[k][cnt];
    		    auto i = index_map.Glob2Loc(row + rIndices[k][rSize[k][1]] - offset);
    		    auto j = rIndices[k][cnt + rSize[k][1] + 1];
    		    if (std::find(IndOfZeros.begin(), IndOfZeros.end(), row + rIndices[k][rSize[k][1]] - offset) == IndOfZeros.end()&& rBufs[k][cnt] != T(0) ){
    		        prod.SetValueLocal(i, j, rBufs[k][cnt]);
    		    }
    	        }
    	    }
    	}
    }

    return prod;
}

template<typename T,typename S>
void parMatrixSparse<T,S>::writeToMatrixMarket(std::string file_name){
	// 1. generate header
        // 2. generate output data
        // 3. perform write
    std::string header;
    std::string data;

    if( (sizeof(T)/sizeof(Base<T>) != 1) ){
        try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for complex matrix, please use writeToMatrixMarketCmplx" << std::endl;
    	}       	   	
    }

    //generate MatrixMarket header
    header.append("%%MatrixMarket matrix coordinate ");

    header.append("real general\n");
    
    header.append("%=================================================================================\n%\n");
    header.append("% This matrix is generated by SMG2S: Scalable Matrix with Given Spectrum\n");
    header.append("% Opensource codes is available: https://github.com/SMG2S/SMG2S\n");
    header.append("% @article{wu2020parallel,\n");
    header.append("%	title={A parallel generator of non-Hermitian matrices computed from given spectra},\n");
    header.append("%	author={Wu, Xinzhe and Petiton, Serge G and Lu, Yutong},\n");
    header.append("%	journal={Concurrency and Computation: Practice and Experience},\n");
    header.append("%	volume={32},\n");
    header.append("%	number={20},\n");
    header.append("%	pages={e5710},\n");
    header.append("%	year={2020},\n");
    header.append("%	publisher={Wiley Online Library},\n");
    header.append("%	doi={https://doi.org/10.1002/cpe.5710}\n");
    header.append("% }\n");
    header.append("%\n%=================================================================================\n");


    S gsize = index_map.GetGlobalSize();
    S gnnz;

    MPI_Allreduce(&nnz_loc, &gnnz, 1, getMPI_Type<S>(), MPI_SUM, comm);

    std::string dim_info = std::to_string(gsize) + " " + std::to_string(gsize) + " " + std::to_string(gnnz) + "\n";
    header.append(dim_info);

    typename std::map<S,T>::iterator it;

    if(dynmat_loc != NULL){
    	for(auto i = 0; i < nrows; i++){
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
    	    	auto k = it->first;
    	    	auto v = it->second;
    	    	data += std::to_string(index_map.Loc2Glob(i)+1) + " " + std::to_string(k+1) + " " + std::to_string(v) + "\n";
	    }
	}
    }

    MPI_File fh;
    MPI_Offset write_offset;
    MPI_Offset text_size;
    MPI_Offset *write_size_per_proc;
    MPI_Status sts;

    write_size_per_proc = (MPI_Offset *)malloc(sizeof(MPI_Offset) * nProcs);
    
    text_size = data.size();

    MPI_Allgather(&text_size, 1, MPI_OFFSET, write_size_per_proc, 1, MPI_OFFSET, comm);


    write_offset = header.size();

    for (auto i = 0; i < ProcID; ++i) {
    	write_offset += write_size_per_proc[i];
    }


    MPI_File_delete (file_name.c_str(), MPI_INFO_NULL);
    MPI_File_open(comm, file_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);


    if (ProcID == 0) {
    	MPI_File_write_at(fh, 0, header.c_str(), header.size(), MPI_CHAR, &sts);
    }

    MPI_File_write_at(fh, write_offset, data.c_str(), data.size(), MPI_CHAR, &sts);

    MPI_File_close(&fh);
    free(write_size_per_proc);
}



template<typename T,typename S>
void parMatrixSparse<T,S>::writeToMatrixMarketCmplx(std::string file_name){
	// 1. generate header
        // 2. generate output data
        // 3. perform write
    std::string header;
    std::string data;

    if( (sizeof(T)/sizeof(Base<T>) != 2) ){
        try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for real matrix, please use writeToMatrixMarket" << std::endl;
    	}  
    }

    //generate MatrixMarket header
    header.append("%%MatrixMarket matrix coordinate ");

    header.append("complex general");
    

    header.append("%=================================================================================\n%\n");
    header.append("% This matrix is generated by SMG2S: Scalable Matrix with Given Spectrum\n");
    header.append("% Opensource codes is available: https://github.com/SMG2S/SMG2S\n");
    header.append("% @article{wu2020parallel,\n");
    header.append("%	title={A parallel generator of non-Hermitian matrices computed from given spectra},\n");
    header.append("%	author={Wu, Xinzhe and Petiton, Serge G and Lu, Yutong},\n");
    header.append("%	journal={Concurrency and Computation: Practice and Experience},\n");
    header.append("%	volume={32},\n");
    header.append("%	number={20},\n");
    header.append("%	pages={e5710},\n");
    header.append("%	year={2020},\n");
    header.append("%	publisher={Wiley Online Library},\n");
    header.append("%	doi={https://doi.org/10.1002/cpe.5710}\n");
    header.append("% }\n");
    header.append("%\n%=================================================================================\n");


    S gsize = index_map.GetGlobalSize();
    S gnnz;

    MPI_Allreduce(&nnz_loc, &gnnz, 1, getMPI_Type<S>(), MPI_SUM, comm);

    std::string dim_info = std::to_string(gsize) + " " + std::to_string(gsize) + " " + std::to_string(gnnz) + "\n";
    header.append(dim_info);

    typename std::map<S,T>::iterator it;

    if(dynmat_loc != NULL){
    	for(auto i = 0; i < nrows; i++){
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
    	    	auto k = it->first;
    	    	auto v = it->second;
    	    	data += std::to_string(index_map.Loc2Glob(i)+1) + " " + std::to_string(k+1) + " " + std::to_string(v.real()) + " " + std::to_string(v.imag()) + "\n";
    	    }	
	}
    }

    MPI_File fh;
    MPI_Offset write_offset;
    MPI_Offset text_size;
    MPI_Offset *write_size_per_proc;
    MPI_Status sts;

    write_size_per_proc = (MPI_Offset *)malloc(sizeof(MPI_Offset) * nProcs);
    
    text_size = data.size();

    MPI_Allgather(&text_size, 1, MPI_OFFSET, write_size_per_proc, 1, MPI_OFFSET, comm);


    write_offset = header.size();

    for (auto i = 0; i < ProcID; ++i) {
    	write_offset += write_size_per_proc[i];
    }


    MPI_File_delete (file_name.c_str(), MPI_INFO_NULL);
    MPI_File_open(comm, file_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);


    if (ProcID == 0) {
    	MPI_File_write_at(fh, 0, header.c_str(), header.size(), MPI_CHAR, &sts);
    }

    MPI_File_write_at(fh, write_offset, data.c_str(), data.size(), MPI_CHAR, &sts);

    MPI_File_close(&fh);
    free(write_size_per_proc);
}

#endif
