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
#include "../utils/MPI_DataType.h"
#include "../parVector/parVector.h"
#include "../nilpotent/nilpotent.h"

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
	void setSpecNonSymm(std::vector<std::complex<Base<T>>> spectrum);

   	//matrix multiplies a nilpotent matrix
	parMatrixSparse<T,S>	MA(Nilpotent<S> nilp);

   	//a nilpotent matrix multiplies a matrix
	parMatrixSparse<T,S>	AM(Nilpotent<S> nilp);


	void	ZeroEntries();
	MatrixCSR<T,S>	ConvertToCSR();

    	void show();
    	void MatView();
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
		if(it->second != T(0)){
		    std::cout <<"("<<it->first << "," << it->second << "); ";
		}	
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

    if(it == dynmat_loc[row].end()){
    	dynmat_loc[row][col] = value;
    	nnz_loc++;
    }
    else{
    	it->second = value;
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
    S i, k;
    if(index_map != X.GetMap()){
    	throw "for AXPY, matrix should be in the same vectorMap";
    }

    if(dynmat_loc != NULL && X.GetDynMatLoc() != NULL){
    	for(i = 0; i < nrows; i++){
    	    std::map<S,T> merge;
    	    merge.insert(dynmat_loc[i].begin(),dynmat_loc[i].end());
    	    merge.insert(X.GetDynMatLoc()[i].begin(),X.GetDynMatLoc()[i].end());
    	    for(it = merge.begin(); it != merge.end(); ++it){
    	    	k = it->first;
    	    	dynmat_loc[i][k] = dynmat_loc[i][k]+X.GetDynMatLoc()[i][k]*scale;
    	    }
    	    merge.clear();	
    	}
    } else if(dynmat_loc == NULL && X.GetDynMatLoc() != NULL){
    	for(i = 0; i < nrows; i++){
    	    std::map<S,T> merge;
    	    merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
    	    for(it = merge.begin(); it != merge.end(); ++it){
    	    	k = it->first;
    	    	dynmat_loc[i][k] = X.GetDynMatLoc()[i][k]*scale;
    	    }
    	    merge.clear();
    	}
    }
}


template<typename T,typename S>
void parMatrixSparse<T,S>::MatAYPX(parMatrixSparse<T,S> X, T scale){
    typename std::map<S,T>::iterator it;

    S i, k;
    if(dynmat_loc != NULL){
    	for(i = 0; i < nrows; i++){
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
    	    	k = it->first;
    	    	dynmat_loc[i][k] = dynmat_loc[i][k]*scale;
	    }
	}
    }

    if(dynmat_loc != NULL && X.GetDynMatLoc() != NULL){
    	for(i = 0; i < nrows; i++){
    	    std::map<S,T> merge;
    	    merge.insert(dynmat_loc[i].begin(),dynmat_loc[i].end());
    	    merge.insert(X.GetDynMatLoc().begin(),X.GetDynMatLoc().end());
    	    for(it = merge.begin(); it != merge.end(); ++it){
    	    	k = it->first;
    	    	dynmat_loc[i][k] = dynmat_loc[i][k]+X.GetDynMatLoc()[i][k];
    	    }
	    merge.clear();
	}
    }

    if(dynmat_loc == NULL && X.GetDynMatLoc() != NULL){
    	for(i = 0; i < nrows; i++){
    	    std::map<S,T> merge;
    	    merge.insert(X.GetDynMatLoc().begin(),X.GetDynMatLoc().end());
    	    for(it = merge.begin(); it != merge.end(); ++it){
    	    	k = it->first;
		dynmat_loc[i][k] = X.GetDynMatLoc()[i][k];
	    }
	    merge.clear();
	}
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
		it->second = 0;
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


    if (isnonsymm && diag_u < 2)
        throw  "for initialisationg of non-symmetric matrix, please ensure abs(diag_u) >= 2 ";


    if (diag_l < diag_u)
        throw "for initialisationg of matrix, please ensure abs(diag_l) < abs(diag_u) ";
            
    

    if (diag_u >= size)
        throw "for initialisationg of matrix, please ensure abs(diag_u) < size ";

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
void parMatrixSparse<T,S>::setSpecNonSymm(std::vector<std::complex<Base<T>>> spectrum){
    
    int size1 = sizeof(T) / sizeof(Base<T>);

    if(size1 == 2){
        printf("Info ]> For generating non-Symmetric matrices, the matrices should be real\n");
        return;
    }

    S ind = 0;
    for(S i = 0; i < spectrum.size();i++){

    	if(spectrum[i].imag() == 0){
    		if(ind >= lower_b && ind < upper_b){
    			SetValue(ind, ind, spectrum[i].real());
    		}
    		ind = ind + 1;
    	}else{
    		if(ind >= lower_b && ind < upper_b){
    			SetValue(ind, ind, spectrum[i].real());
    			SetValue(ind, ind + 1, spectrum[i].imag());
    		}

    		if(ind+1 >= lower_b && ind+1 < upper_b){
    			SetValue(ind+1, ind+1, spectrum[i].real());
    			SetValue(ind+1, ind, -spectrum[i].imag());
    		}
    		ind = ind + 2;
    	}    	
    }

    S size;
    if (spectrum[spectrum.size()-1].imag() == 0){
    	size = ind;
    }else{
    	size = ind - 1;
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
    	    	    if (std::find(IndOfZeros.begin(), IndOfZeros.end(), j-offset) == IndOfZeros.end()){
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

    std::vector<int> targetProcs(ncols);
    std::vector<int> originProcs(ncols);
    
    typename std::map<S,T>::iterator it;

    if(prod.dynmat_loc == NULL){
    	prod.dynmat_loc = new std::map<S,T> [nrows];
    }
    

    if(dynmat_loc != NULL){
    	for(auto row = 0; row < nrows; row++){
    	    for(it = dynmat_loc[row].begin(); it != dynmat_loc[row].end(); ++it){
    	    	auto i = row - offset;
    	    	auto j =  it->first;
    	    	if(i >= 0){
    	    	    //if not the index of zeros in nilpotent
    	    	    if (std::find(IndOfZeros.begin(), IndOfZeros.end(), index_map.Loc2Glob(row)-offset) == IndOfZeros.end()){
    	    	    	prod.SetValueLocal(i, j, it->second);
    	    	    }
    	    	}
    	    }	
    	}
    }


    //communication part
    std::vector<T> sBuf;
    std::vector<S> sIndex; //    sRoffsets + sColIndx
    S sSize[2]={0, 0};

    std::vector<T> rBuf;
    std::vector<S> rIndex; //    rRoffsets + rColIndx

    S rSize[2]={0, 0};

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
/*
    for(S i = 0; i < ncols; i++){
    	if(ProcID != 0){
	    	if(i < uprocbound_map[ProcID] && i >= lprocbound_map[ProcID]){
	    	    if(originProcs[i] != targetProcs[i]){
	    	    	std::cout << "Sending: ID " << ProcID << ": origin row "<< i << ": " << originProcs[i] << " -> " << targetProcs[i]  << " :" << i - offset << " target row"<< std::endl;
	    	    }
	    	}
    	}
    	if(ProcID != nProcs - 1){
	    	if(i < uprocbound_map[ProcID+1] && i >= lprocbound_map[ProcID+1]){
	    	    if(originProcs[i] != targetProcs[i]){
	    	    	std::cout << "Receiving: ID " << ProcID << ": origin row "<< i << ": " << originProcs[i] << " -> " << targetProcs[i]  << " :" << i - offset << " target row"<< std::endl;
	    	    }
	    	}
    	}    	
  	
    }
*/  
    /*
    S count = 0;
    sRoffsets.push_back(count);
    for(S i = 0; i < ncols; i++){
    	if(ProcID != 0){
    	    if(i < uprocbound_map[ProcID] && i >= lprocbound_map[ProcID]){
    	    	
    	        if(originProcs[i] != targetProcs[i]){
    	            for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){	
    	            	sColIndx.push_back(it->first);   
    	            	sBuf.push_back(it->second);
    	            	count++;
    	    	    }
    	    	}
    	    }
    	}
    }
    */

    MPI_Request	rtypereq, stypereq;
    MPI_Status	typestat;
    
    int up = ProcID - 1;
    int down = ProcID + 1;

    if(ProcID != 0 && dynmat_loc != NULL){
    	S count = lprocbound_map[ProcID];
    	sIndex.insert(sIndex.begin(), count);
    	for( S i = lprocbound_map[ProcID]; i < uprocbound_map[ProcID]; i++){
    	    if(originProcs[i] != targetProcs[i]){
    	    	S j = index_map.Glob2Loc(i);
    	    	for(it = dynmat_loc[j].begin(); it != dynmat_loc[j].end(); ++it){	
    	    	    sIndex.push_back(it->first);
    	    	    sBuf.push_back(it->second);
    	    	    count++;
    	    	}
    	    	sIndex.insert(sIndex.begin(), count);
    	    }
    	}
    	
    	sSize[0] = sBuf.size();
    	sSize[1] = sIndex.size() - sBuf.size() - 1;
    	
    	for(auto row = 0; row < sSize[1] ; row++){
    	    for(auto cnt = sIndex[sSize[1] - row] - sIndex[sSize[1]] ; cnt < sIndex[sSize[1] - row-1] - sIndex[sSize[1]]; cnt++){
    		auto col = sIndex[cnt + sSize[1] + 1];
    		auto val = sBuf[cnt];
    		//std::cout << "Send: locind (" << row << "," << col << ") and "  << "globind (" << row + sIndex[sSize[1]] << "," << col << "): " << val << std::endl;
    	    }
    	}
    	
    }

    if(ProcID != 0 && dynmat_loc != NULL){
    	MPI_Isend(sSize, 2, getMPI_Type<S>(), up, 1, comm, &stypereq);
    }

    if(ProcID != nProcs - 1){
    	MPI_Irecv(rSize, 2, getMPI_Type<S>(), down, 1, comm, &rtypereq);
    	MPI_Wait(&rtypereq,&typestat);
    	rBuf.resize(rSize[0]);
    	rIndex.resize(rSize[0]+rSize[1]+1);
    }

    if(ProcID != 0 && dynmat_loc != NULL){
    	MPI_Isend(sBuf.data(), sSize[0], getMPI_Type<T>(), up, 1, comm, &stypereq);
    }

    if(ProcID != nProcs - 1){
    	MPI_Irecv(rBuf.data(), rSize[0], getMPI_Type<T>(), down, 1, comm, &rtypereq);
    	MPI_Wait(&rtypereq,&typestat);
    }

    if(ProcID != 0 && dynmat_loc != NULL){
    	MPI_Isend(sIndex.data(), sSize[0]+sSize[1]+1, getMPI_Type<S>(), up, 1, comm, &stypereq);
    }

    if(ProcID != nProcs - 1){
    	MPI_Irecv(rIndex.data(), rSize[0]+rSize[1]+1, getMPI_Type<S>(), down, 1, comm, &rtypereq);
    	MPI_Wait(&rtypereq,&typestat);

    	for(auto row = 0; row < rSize[1] ; row++){
    	    for(auto cnt = rIndex[rSize[1] - row] - rIndex[rSize[1]] ; cnt < rIndex[rSize[1] - row-1] - rIndex[rSize[1]]; cnt++){
    		auto col = rIndex[cnt + rSize[1] + 1];
    		auto val = rBuf[cnt];
    		//std::cout << "Recv: locind (" << row << "," << col << ") and "  << "globind (" << index_map.Glob2Loc(row + rIndex[rSize[1]] - offset) << "," << col << "): " << val << std::endl;
    	    }
    	}

    	for(auto row = 0; row < rSize[1] ; row++){
    	    for(auto cnt = rIndex[rSize[1] - row] - rIndex[rSize[1]] ; cnt < rIndex[rSize[1] - row-1] - rIndex[rSize[1]]; cnt++){
    		auto col = rIndex[cnt + rSize[1] + 1];
    		auto val = rBuf[cnt];
    		auto i = index_map.Glob2Loc(row + rIndex[rSize[1]] - offset);
    		auto j = rIndex[cnt + rSize[1] + 1];
    		if (std::find(IndOfZeros.begin(), IndOfZeros.end(), row + rIndex[rSize[1]] - offset) == IndOfZeros.end()){
    		    prod.SetValueLocal(i, j, rBuf[cnt]);
    		}
    	    }
    	}
    }


    return prod;
}

#endif
