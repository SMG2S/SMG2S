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

   Part of basic data structures' implementation of this file refers to this technical report
   (http://orbit.dtu.dk/files/51272329/tr12_10_Alexandersen_Lazarov_Dammann_1.pdf)
*/


#ifndef __PAR_MATRIX_SPARSE_H__
#define __PAR_MATRIX_SPARSE_H__

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "../utils/MPI_DataType.h"
#include "../parVector/parVector.h"
//#include "../utils/utils.h"
#include "MatrixCSR.h"

#ifdef __USE_COMPLEX__
#include <complex>
#endif

template<typename T, typename S>
class parMatrixSparse
{

	private:
		std::map<S, T> *dynmat_lloc, *dynmat_gloc;

		//size of local matrix
		S	ncols, nrows;

		S	nnz_lloc, nnz_gloc, nnz_loc;

		parVectorMap<S>	*x_index_map;
		parVectorMap<S>	*y_index_map;

		S	njloc;
		S	lower_x, lower_y, upper_x, upper_y;

		MPI_Comm comm;

		// mpi size and rank
		int ProcID, nProcs;

		// pointers and buffers for setting up communications
		S	*VNumRecv, *VNumSend;
		S	maxRecv, maxSend;
		S	**Rbuffer, **Sbuffer;

		MPI_Datatype *DTypeRecv , *DTypeSend ;

	public:

		MatrixCSR<T,S> *CSR_lloc, *CSR_gloc, *CSR_loc;

		std::map<S, T> *dynmat_loc;

		//constructor
		parMatrixSparse();

		parMatrixSparse(parVector<T,S> *XVec, parVector<T,S> *YVec);
		//deconstructor
		~parMatrixSparse();

		//get
		parVectorMap<S> *GetXMap(){return x_index_map;};
		parVectorMap<S> *GetYMap(){return y_index_map;};

		MPI_Comm GetComm(){
			return x_index_map->GetCurrentComm();
		}
		S	GetXLowerBound();
		S	GetYLowerBound();

	    S   GetXUpperBound();
        S   GetYUpperBound();

		void	GetTrueLocalSize(S &rs, S &cs){
			rs = nrows;
			cs = njloc;
		};

		void	GetLocalSize(S &rs, S &cs){
			rs = nrows;
			cs = ncols;
		};

		std::map<S,T>	*GetDynMatGLobLoc(){return dynmat_lloc;};
		std::map<S,T>	*GetDynMatGlobLoc(){return dynmat_gloc;};

		std::map<S,T>	*GetDynMatLoc(){return dynmat_loc;};

		MatrixCSR<T,S>	*GetCSRLocLoc(){return CSR_lloc;};
		MatrixCSR<T,S>	*GetCSRGlobLoc(){return CSR_gloc;};

		//add
		void	AddValueLocal( S row, S col, T value);
  		void 	AddValuesLocal( S nindex, S *rows, S *cols, T *values);

		//global add
		void	AddValue(S row, S col, T value);

  		//set
		void	SetValueLocal( S row, S col, T value);
  		void 	SetValuesLocal( S nindex, S *rows, S *cols, T *values);

		//global set
		void	SetValue(S row, S col, T value);

		//get
		T		GetLocalValue(S row, S col);
		T		GetValue(S row, S col);


		//combine gloc + lloc -> loc together
		void	glocPlusLloc();

		void	llocToGlocLoc();

		//show
		void	MatView();

		void	LOC_MatView();


  		//LOC set
		void	Loc_SetValueLocal( S row, S col, T value);
  	void 	Loc_SetValuesLocal( S nindex, S *rows, S *cols, T *values);

		//LOC global set
		void	Loc_SetValue(S row, S col, T value);

		//LOC get
		T		Loc_GetLocalValue(S row, S col);
		T		Loc_GetValue(S row, S col);


		//set mat diagonal by vector given
		void	SetDiagonal(parVector<T,S> *diag);

		//Loc set mat diagonal by vector given
		void	Loc_SetDiagonal(parVector<T,S> *diag);

		//Mat Scale
		void	MatScale(T scale);

		//Loc Mat Scale
		void	Loc_MatScale(T scale);

		//Loc AXPY
		void	Loc_MatAXPY(parMatrixSparse<T,S> *X, T scale);

		//Loc AYPX
		void    Loc_MatAYPX(parMatrixSparse<T,S> *X, T scale);

		//Reader
		void	ReadExtMat();

		//Writer

		// convert from dyn to csr
		void	ConvertToCSR();

		// convert from dyn to csr
		void	Loc_ConvertToCSR();

		// Zeros all entries with keeping the previous matrix pattern
		void	ZeroEntries();

		// Loc: Zeros all entries with keeping the previous matrix pattern
		void	Loc_ZeroEntries();

   	//matrix multiple a special nilpotent matrix
		void	MA(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod);

		//special nilpotent matrix multiple another matrix
		void	AM(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod);


};


template<typename T, typename S>
parMatrixSparse<T,S>::parMatrixSparse()
{
	dynmat_lloc = NULL;
	dynmat_gloc = NULL;
	dynmat_loc  = NULL;

	CSR_lloc = NULL;
	CSR_gloc = NULL;
	CSR_loc = NULL;

	nnz_lloc = 0;
	nnz_gloc = 0;
	nnz_loc = 0;

	ncols = 0;
	nrows = 0;
	njloc = 0;

	lower_x = 0;
	lower_y = 0;
	upper_x = 0;
	upper_y = 0;

	x_index_map = NULL;
	y_index_map = NULL;

//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
//	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	VNumRecv = NULL;
	VNumSend = NULL;
	Rbuffer = NULL;
	Sbuffer = NULL;
	DTypeRecv = NULL;
	DTypeSend = NULL;
}

template<typename T, typename S>
parMatrixSparse<T,S>::parMatrixSparse(parVector<T,S> *XVec, parVector<T,S> *YVec)
{
	dynmat_lloc = NULL;
	dynmat_gloc = NULL;
	dynmat_loc  = NULL;

	CSR_lloc = NULL;
	CSR_gloc = NULL;
    CSR_loc = NULL;

	nnz_lloc = 0;
	nnz_gloc = 0;
	nnz_loc  = 0;

	x_index_map = NULL;
	y_index_map = NULL;

	ncols = 0;
	nrows = 0;
	njloc = 0;

	lower_x = 0;
	lower_y = 0;

	upper_x = 0;
	upper_y = 0;

	VNumRecv = NULL;
	VNumSend = NULL;
	Rbuffer = NULL;
	Sbuffer = NULL;

	DTypeRecv = NULL;
	DTypeSend = NULL;

	//get vector map for x and y direction
	x_index_map = XVec->GetVecMap();
	x_index_map->AddUser();
	y_index_map = YVec->GetVecMap();
	y_index_map->AddUser();

	if(x_index_map != NULL && y_index_map != NULL){
		//get num of rows and cols in this mpi procs
		ncols = x_index_map->GetGlobalSize();
		nrows = y_index_map->GetLocalSize();
		njloc = x_index_map->GetLocalSize();
		//get upper and lower bounds
		lower_x = x_index_map->GetLowerBound();
		lower_y = y_index_map->GetLowerBound();
		upper_x = x_index_map->GetUpperBound();
		upper_y = y_index_map->GetUpperBound();
	}

	comm = GetComm();
	MPI_Comm_rank(comm, &ProcID);
	MPI_Comm_size(comm, &nProcs);

//	printf("ProId = %d, ProcsSize = %d\n",ProcID,nProcs);

}

template<typename T, typename S>
parMatrixSparse<T,S>::~parMatrixSparse()
{
	//if index map is defined
	if(x_index_map != NULL){
		x_index_map->DeleteUser();

//		if(x_index_map->GetUser() == 0){
			delete x_index_map;
//		}
	}


	if(y_index_map != NULL){
		y_index_map->DeleteUser();

//		if(x_index_map->GetUser() == 0){
			delete y_index_map;
//		}
	}

	//if dynmat has been defined
	if(dynmat_lloc != NULL){
		delete [] dynmat_lloc;
	}
	if(dynmat_gloc != NULL){
		delete [] dynmat_gloc;
	}
	if(CSR_lloc != NULL){
		delete [] CSR_lloc;
	}
	if(CSR_gloc != NULL){
		delete [] CSR_gloc;
	}
	if(CSR_gloc != NULL){
		delete [] CSR_loc;
	}
	if(VNumRecv != NULL){
		delete [] VNumRecv;
	}
	if(VNumSend != NULL){
		delete [] VNumSend;
	}
	if(Rbuffer != NULL){
		int i;
		for(i = 0; i < nProcs; i++){
			if (Rbuffer[i] != NULL){
				delete [] Rbuffer[i];
			}
		}
		delete [] Rbuffer;
	}

	if(Sbuffer != NULL){
		int i;
	        for(i = 0; i < nProcs; i++){
                        if (Sbuffer[i] != NULL){
                                delete [] Sbuffer[i];
                        }
                }
                delete [] Sbuffer;
	}


	if(DTypeRecv != NULL){
		int i;
		for(i = 0; i < nProcs; i++){
			if(DTypeSend[i] != MPI_DATATYPE_NULL){
				MPI_Type_free(&DTypeSend[i]);
			}
		}
		delete [] DTypeSend;
	}
}


template<typename T, typename S>
S parMatrixSparse<T,S>::GetXLowerBound(){
	if(x_index_map != NULL){
		return x_index_map->GetLowerBound();
	}
	else{
		return 0;
	}
}

template<typename T, typename S>
S parMatrixSparse<T,S>::GetXUpperBound(){
        if(x_index_map != NULL){
                return x_index_map->GetUpperBound();
        }
        else{
                return 0;
        }
}

template<typename T,typename S>
S parMatrixSparse<T,S>::GetYLowerBound(){
        if(y_index_map != NULL){
                return y_index_map->GetLowerBound();
        }
        else{
                return 0;
        }
}

template<typename T, typename S>
S parMatrixSparse<T,S>::GetYUpperBound(){
        if(y_index_map != NULL){
                return y_index_map->GetUpperBound();
        }
        else{
                return 0;
        }
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AddValueLocal(S row, S col, T value)
{
	typename std::map<S,T>::iterator it;
	//if location is inside of local area then add to local dynamic map
	if((row < nrows && row >= 0) && (col < upper_x && col >= lower_x && col >= 0)){
		if(dynmat_lloc == NULL){
			dynmat_lloc = new std::map<S,T> [nrows];
		}
		it = dynmat_lloc[row].find(col);
		if(it == dynmat_lloc[row].end()){
			dynmat_lloc[row][col] = value;
			nnz_lloc++;
		}
		else{
			it->second = it->second + value;
		}
	//if location is inside of local-global area
	}
	else if ((row < nrows && row >= 0) && (col >= upper_x || col < lower_x) && (col >= 0)){
		if(dynmat_gloc == NULL){
			dynmat_gloc = new std::map<S,T> [nrows];
		}
		it = dynmat_gloc[row].find(col);
		if(it == dynmat_gloc[row].end()){
			dynmat_gloc[row][col] = value;
			nnz_gloc++;
		}
		else{
			it->second = it->second + value;
		}
	}
}

template<typename T,typename S>
T parMatrixSparse<T,S>::GetLocalValue(S row, S col)
{
	if((row < nrows && row >= 0) && (col < upper_x && col >= lower_x && col >= 0)){
		return dynmat_lloc[row][col];
	}
	else if ((row < nrows && row >= 0) && (col >= upper_x || col < lower_x) && (col >= 0)){
		return dynmat_gloc[row][col];
	}
	else return 0.0;
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AddValuesLocal(S nindex, S *rows, S *cols, T *values)
{
	typename std::map<S,T>::iterator it;

	for( S i = 0; i < nindex; i++){
		AddValueLocal(rows[i],cols[i],values[i]);
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AddValue(S row, S col, T value)
{

	if((row >= lower_y) && (row < upper_y) && (col < ncols)){
		AddValueLocal(y_index_map->Glob2Loc(row),col,value);
	}
}

//set
template<typename T,typename S>
void parMatrixSparse<T,S>::SetValueLocal( S row, S col, T value)
{
	typename std::map<S,T>::iterator it;
	//if location is inside of local area then add to local dynamic map
	if((row < nrows && row >= 0) && (col < upper_x && col >= lower_x && col >= 0)){
		if(dynmat_lloc == NULL){
			dynmat_lloc = new std::map<S,T> [nrows];
		}
		it = dynmat_lloc[row].find(col);
		if(it == dynmat_lloc[row].end()){
			dynmat_lloc[row][col] = value;
			nnz_lloc++;
		}
		else{
	//		it->second = it->second + value;
			it->second = value;
		}
	//if location is inside of local-global area
	}
	else if ((row < nrows && row >= 0) && (col >= upper_x || col < lower_x) && (col >= 0)){
		if(dynmat_gloc == NULL){
			dynmat_gloc = new std::map<S,T> [nrows];
		}
		it = dynmat_gloc[row].find(col);
		if(it == dynmat_gloc[row].end()){
			dynmat_gloc[row][col] = value;
			nnz_gloc++;
		}
		else{
//			it->second = it->second + value;
			it->second = value;
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::SetValuesLocal( S nindex, S *rows, S *cols, T *values)
{
	typename std::map<S,T>::iterator it;

	for( S i = 0; i < nindex; i++){
		SetValueLocal(rows[i],cols[i],values[i]);
	}
}

//global set
template<typename T,typename S>
void parMatrixSparse<T,S>::SetValue(S row, S col, T value)
{
	S local_row = y_index_map->Glob2Loc(row);

	if((row >= lower_y) && (row < upper_y) && (col < ncols)){
		SetValueLocal(local_row,col,value);
	}
}


template<typename T,typename S>
T parMatrixSparse<T,S>::GetValue(S row, S col)
{
	S local_row = y_index_map->Glob2Loc(row);
	if(local_row >= 0 && local_row < nrows){
		return GetLocalValue(local_row, col);
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::glocPlusLloc(){

	S i,j;
	T v;
	typename std::map<S,T>::iterator it;

	if(ProcID == 0) {std::cout << "Combine the block-diagonal part and non block-diagonal part of parallel matrix together " << std::endl;}

	if(dynmat_loc == NULL){
		dynmat_loc = new std::map<S, T> [nrows];
	}
	for(i = 0; i < nrows; i++){
		if((dynmat_gloc != NULL) && (dynmat_lloc != NULL)){
			dynmat_loc[i].insert(dynmat_lloc[i].begin(),dynmat_lloc[i].end());
			dynmat_loc[i].insert(dynmat_gloc[i].begin(),dynmat_gloc[i].end());
		}
		else if ((dynmat_lloc != NULL)){
			dynmat_loc[i].insert(dynmat_lloc[i].begin(),dynmat_lloc[i].end());
		}
		else if (dynmat_gloc != NULL){
			dynmat_loc[i].insert(dynmat_gloc[i].begin(),dynmat_gloc[i].end());
		}
		else{
			if(ProcID == 0) {std::cout << "Cannot execute the function glocPlusLloc since the given matrix is NULL " << std::endl;}
		}
	}
	nnz_loc = nnz_gloc + nnz_lloc;

}

template<typename T,typename S>
void parMatrixSparse<T,S>::llocToGlocLoc()
{
	typename std::map<S,T>::iterator it;

	S col;

	//if location is inside of local area then add to local dynamic map
	if(dynmat_loc != NULL){
		if(dynmat_lloc == NULL){
			dynmat_lloc = new std::map<S,T> [nrows];
		}
		if(dynmat_gloc == NULL){
			dynmat_gloc = new std::map<S,T> [nrows];
		}

		for(S i = 0; i < nrows; i++){
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
				col = it->first;
				if(col < upper_x && col >= lower_x && col >= 0){
					dynmat_lloc[i][col] = it->second;
				}
				else if((col >= upper_x || col < lower_x) && (col >= 0)){
					dynmat_gloc[i][col] = it->second;
				}
			}
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::MatView(){

	S i, j;
	T v;
	typename std::map<S,T>::iterator it;

	if(ProcID == 0) {std::cout << "Parallel MatView: " << std::endl;}

	for (i = 0; i < nrows; i++){
		std::map<S,T> merge;
		std::cout << "row " << y_index_map->Loc2Glob(i) << ": ";

		if((dynmat_gloc != NULL) && (dynmat_lloc != NULL)){
			merge.insert(dynmat_lloc[i].begin(),dynmat_lloc[i].end());
			merge.insert(dynmat_gloc[i].begin(),dynmat_gloc[i].end());

			for(it = merge.begin(); it != merge.end(); ++it){
				std::cout <<"("<<it->first << "," << it->second << "); ";
			}
			merge.clear();
		}
		else if ((dynmat_lloc != NULL)){
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); ++it){
				std::cout <<"("<<it->first << "," << it->second << "); ";
			}
		}
		else if (dynmat_gloc != NULL){
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); ++it){
				std::cout <<"("<<it->first << "," << it->second << "); ";
			}
		}

		std::cout << std::endl;
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::LOC_MatView(){

	S i;
	T v;
	typename std::map<S,T>::iterator it;

	if(ProcID == 0) {std::cout << "LOC MODE Parallel MatView: " << std::endl;}

	for (i = 0; i < nrows; i++){
		if(dynmat_loc != NULL){
			std::cout << "row " << y_index_map->Loc2Glob(i) << ": ";
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
#if defined (__USE_COMPLEX__)
				if(it->second.real() != 0 || it->second.imag() != 0){
					std::cout <<"("<<it->first << "," << it->second << "); ";
				}
#else
				if(it->second != 0){
					std::cout <<"("<<it->first << "," << it->second << "); ";
				}
#endif
			}
		}
	std::cout << std::endl;
	}
}


//Loc set
template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetValueLocal( S row, S col, T value)
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
void parMatrixSparse<T,S>::Loc_SetValuesLocal( S nindex, S *rows, S *cols, T *values)
{
	typename std::map<S,T>::iterator it;

	for( S i = 0; i < nindex; i++){
		Loc_SetValueLocal(rows[i],cols[i],values[i]);
	}
}

//Loc global set
template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetValue(S row, S col, T value)
{
	S local_row = y_index_map->Glob2Loc(row);

	if (local_row >= 0 && local_row < nrows){
		Loc_SetValueLocal(local_row,col,value);
	}

}

template<typename T,typename S>
T parMatrixSparse<T,S>::Loc_GetLocalValue(S row, S col)
{
	if(dynmat_loc != NULL){
		return dynmat_loc[row][col];
	}
	else return 0.0;
}


template<typename T,typename S>
T parMatrixSparse<T,S>::Loc_GetValue(S row, S col)
{
	S local_row = y_index_map->Glob2Loc(row);
	if(local_row >= 0 && local_row < nrows){
		return Loc_GetLocalValue(local_row, col);
	}
}



template<typename T,typename S>
void parMatrixSparse<T,S>::SetDiagonal(parVector<T,S> *diag)
{
	S    lbound;
	S 	 ubound;

	if (nrows != njloc ){
		if(ProcID == 0){
			printf("ERROR: cannot set diagonal for non-square matrix.");
		}
	}
	else{
		T *a = diag->GetArray();
		S local_size = diag->GetArraySize();

		lbound = diag->GetLowerBound();
		ubound = diag->GetUpperBound();

		//std::cout << "local size  = " << local_size << "  lower bound = " << lbound << "  upper bound = " << ubound << std::endl;
		for(S i = 0; i < local_size; i++){
			SetValueLocal(i,diag->Loc2Glob(i),a[i]);
		}

	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetDiagonal(parVector<T,S> *diag)
{
	S    lbound;
	S 	 ubound;

	if (nrows != njloc ){
		if(ProcID == 0){
			printf("ERROR: cannot set diagonal for non-square matrix.");
		}
	}
	else{
		T *a = diag->GetArray();
		S local_size = diag->GetArraySize();

		lbound = diag->GetLowerBound();
		ubound = diag->GetUpperBound();

		//std::cout << "local size  = " << local_size << "  lower bound = " << lbound << "  upper bound = " << ubound << std::endl;
		for(S i = 0; i < local_size; i++){
			Loc_SetValueLocal(i,diag->Loc2Glob(i),a[i]);
		}

	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::ConvertToCSR()
{
	S 	count, i, j;
	T	v;
	typename std::map<S,T>::iterator it;

	if(dynmat_lloc != NULL){
		//allocate csr matrix


		CSR_lloc = new MatrixCSR<T,S>(nnz_lloc, nrows);

//		CSR_lloc = new MatrixCSR<T,S>(nrows);

		//convert local local to csr
		count = 0;
		//CSR_lloc->rows.push_back(0);

		for(i = 0; i < nrows; i++){
			CSR_lloc->rows.push_back(count);
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); it++){
				j = it->first;
				v = it->second;
				CSR_lloc->vals.push_back(v);
				CSR_lloc->cols.push_back(j);
				count++;
			}
		}
		CSR_lloc->rows.push_back(nnz_lloc);
	}

	if(dynmat_gloc != NULL){


		CSR_gloc = new MatrixCSR<T,S>(nnz_gloc, nrows);
//		CSR_gloc = new MatrixCSR<T,S>(nrows);

		//convert global-local to CSR
		count = 0;
		//CSR_gloc->rows.push_back(0);
		for(i = 0; i < nrows; i++){
			CSR_gloc->rows.push_back(count);
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); it++){
				j = it->first;
				v = it->second;
				CSR_gloc->vals.push_back(v);
				CSR_gloc->cols.push_back(j);
				count++;
			}
		}
		CSR_gloc->rows.push_back(nnz_gloc);
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_ConvertToCSR()
{
	S 	count, i, j;
	T	v;
	typename std::map<S,T>::iterator it;

	if(dynmat_loc != NULL){
		//allocate csr matrix


		CSR_loc = new MatrixCSR<T,S>(nnz_loc, nrows);

//		CSR_lloc = new MatrixCSR<T,S>(nrows);

		count = 0;
		//CSR_lloc->rows.push_back(0);

		for(i = 0; i < nrows; i++){
			CSR_loc->rows.push_back(count);
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
#if defined (__USE_COMPLEX__)
				if(it->second.imag() != 0 || it->second.real() != 0){
					j = it->first;
					v = it->second;
					CSR_loc->vals.push_back(v);
					CSR_loc->cols.push_back(j);
					count++;
				}
#else
				if(it->second != 0){
					j = it->first;
					v = it->second;
					CSR_loc->vals.push_back(v);
					CSR_loc->cols.push_back(j);
					count++;
				}
#endif

			}
		}
		CSR_loc->rows.push_back(count);
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::ReadExtMat()
{
	//Reader
	std::ifstream file("matrix.mat");
	std::string line;
	S	row, col;
	T	value;
	row = 0;
	col = 0;
	value = 0.0;
	S quit = 0;

	if(ProcID == 0) printf("start reading matrix !!!!\n");
	//Past no number content in the file
	while(std::getline(file,line)){
		row = 0;
		col = 0;
		value = 0.0;
		std::stringstream linestream(line);
		linestream >> row >> col >> value;
		if(row != 0 && col != 0 && value != 0.0){
			break;
		}
	}

	//read values and add them to matrix
	while(std::getline(file,line)){
		std::stringstream linestream(line);
		linestream >> row >> col >> value;
		row = row - 1;
		col = col - 1;
		if((row >= lower_y && row < upper_y) && (col < ncols)){
			AddValueLocal(y_index_map->Glob2Loc(row),col,value);
		}
	}
	if(ProcID == 0) printf("The matrix is readed !!!!\n");
}





template<typename T,typename S>
void parMatrixSparse<T,S>::MatScale(T scale){
	typename std::map<S,T>::iterator it;

	S i;
	T v;

	if(dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}

	if(dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_MatAXPY(parMatrixSparse<T,S> *X, T scale){

	typename std::map<S,T>::iterator it, itv, itvv;

	S i, k;
	T v;

	if(dynmat_loc != NULL && X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(dynmat_loc[i].begin(),dynmat_loc[i].end());
			merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_loc[i][k] = dynmat_loc[i][k]+X->dynmat_loc[i][k]*scale;
			}
			merge.clear();
		}
	}

	if(dynmat_loc == NULL && X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_loc[i][k] = X->dynmat_loc[i][k]*scale;
			}
			merge.clear();
		}
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_MatScale(T scale){
	typename std::map<S,T>::iterator it;

	S i;
	T v;

	if(dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_MatAYPX(parMatrixSparse<T,S> *X, T scale){

	typename std::map<S,T>::iterator it, itv, itvv;

	S i, k;
	T v;
	if(dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
				k = it->first;
				dynmat_loc[i][k] = dynmat_loc[i][k]*scale;
			}
		}
	}

	if(dynmat_loc != NULL && X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(dynmat_loc[i].begin(),dynmat_loc[i].end());
			merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_loc[i][k] = dynmat_loc[i][k]+X->dynmat_loc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_loc == NULL && X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_loc[i][k] = X->dynmat_loc[i][k];
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
	T v;

	if(dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); it++){
				it->second = 0;
			}
		}
	}

	if(dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); it++){
				it->second = 0;
			}
		}
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_ZeroEntries()
{
	typename std::map<S,T>::iterator it;

	S i;
	T v;

	if(dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
				it->second = 0;
			}
		}
	}
}


//matrix multiple a special nilpotent matrix
template<typename T,typename S>
void parMatrixSparse<T,S>::MA(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod)
{
	S i, j, k;
	T v;

	typename std::map<S,T>::iterator it;

	//use the given nilpotency matrix, MA operation will make elements of matrix right move diaPosition-1 offset.
	//And the positions of 0: pos = nbOne*integer - 1

#ifndef _OPENMP
	for(i = 0; i < nrows; i++){
		if(dynmat_loc == NULL) {
			return;
		}
		if(dynmat_loc != NULL && prod->dynmat_loc == NULL){
			prod->dynmat_loc = new std::map<S,T> [nrows];
		}

		for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
			j = it->first + nilp.diagPosition - 1;
			k = (j+1)%(nilp.nbOne + 1);
			if(j < ncols && k != 0){
				prod->dynmat_loc[i][j] = it->second;
			}
		}
	}

#else

	if(dynmat_loc == NULL) {
		return;
	}
	else if (dynmat_loc != NULL && prod->dynmat_loc == NULL){
		prod->dynmat_loc = new std::map<S,T> [nrows];
	}
#pragma omp parallel for
	for(i = 0; i < nrows; i++){
		for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
			j = it->first + nilp.diagPosition - 1;
			k = (j+1)%(nilp.nbOne + 1);
			if(j < ncols && k != 0){
				prod->dynmat_loc[i][j] = it->second;
			}
		}
	}
#endif

}

//special nilpotent matrix multiple another matrix
template<typename T,typename S>
void parMatrixSparse<T,S>::AM(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod)
{
	S i, j, k, p, q, loc;
	T v;
	bool sendflg = true, recvflag = true;

	typename std::map<S,T>::iterator it;

	MPI_Request	*Rreqs, *Sreqs, rtypereq, stypereq, *idxRreqs, *idxSreqs;

	MPI_Status	*Rstat, *Sstat, typestat, *idxRstat, *idxSstat;

	int up, down;

	int			tag1[nilp.diagPosition - 1], tag2[nilp.diagPosition - 1], tagtype = nilp.diagPosition - 1 + nilp.diagPosition ;
	int			Rtag, Stag;
	S           count, count1;

	count = 0;
	count1 = 0;

	up = ProcID - 1;
	down =ProcID + 1;
	S size[nilp.diagPosition - 1], rsize[nilp.diagPosition - 1];

	for(S a = 0; a < nilp.diagPosition - 1; a++){
		tag1[a] = a;
		tag2[a] = nilp.diagPosition - 1 + a;
	}

	if(ProcID != 0){
		for(S a = 0; a < nilp.diagPosition - 1; a++){
			size[a] = dynmat_loc[a].size();
		}
	}

	if(ProcID != nProcs - 1){
		for(S a = 0; a < nilp.diagPosition - 1; a++){
			rsize[a] = 0;
		}
	}



	if(ProcID != 0){
		MPI_Isend(size, nilp.diagPosition - 1, MPI_INT, up, tagtype, comm, &stypereq);
	}

	if(ProcID != nProcs - 1){
		MPI_Irecv(rsize,nilp.diagPosition - 1, MPI_INT, down, tagtype, comm, &rtypereq);
	}

	if(ProcID != nProcs - 1){
		MPI_Wait(&rtypereq,&typestat);
	}

	Rtag = 0; Stag = 1;

	//MPI non-blocking requests and status

	Rreqs = new MPI_Request [nilp.diagPosition - 1];
	Sreqs = new MPI_Request [nilp.diagPosition - 1];
	Rstat = new MPI_Status [nilp.diagPosition - 1];
	Sstat = new MPI_Status [nilp.diagPosition - 1];

	idxRreqs = new MPI_Request [nilp.diagPosition - 1];
	idxSreqs = new MPI_Request [nilp.diagPosition - 1];
	idxRstat = new MPI_Status [nilp.diagPosition - 1];
	idxSstat = new MPI_Status [nilp.diagPosition - 1];

	if(dynmat_loc == NULL){
		return;
	}

	for(p = nilp.diagPosition - 1; p < nrows; p++){
		if(prod->dynmat_loc == NULL){
			prod->dynmat_loc = new std::map<S,T> [nrows];
		}

		i = p - nilp.diagPosition + 1;
		q = y_index_map->Loc2Glob(i);

		k = (q + 1)%(nilp.nbOne + 1);

		for(it = dynmat_loc[p].begin(); it != dynmat_loc[p].end(); ++it){

			j = it->first;

			if(i >= 0 && i < nrows && k != 0){
				prod->dynmat_loc[i][j] = it->second;
			}
		}

	}

	MPI_Datatype MPI_SCALAR = MPI_Scalar<T>();

	MPI_Barrier(comm);

#if defined (__USE_COMPLEX__) && defined(__USE_DOUBLE__)
//complex double
	double *sBuf, *rBuf;

	S *sIndx, *rIndx;
	S cnt = 0, cnt2 = 0;

	for(S b = 0; b < nilp.diagPosition -1; b++){

		if(ProcID != 0){
			loc = y_index_map->Loc2Glob(b);
			q = loc - nilp.diagPosition + 1 + b;
			if((q + 1)%(nilp.nbOne + 1) == 0){
				sendflg = false;
			}
		}

		if(ProcID != nProcs - 1){
			loc = y_index_map->Loc2Glob(nrows - b - 1);
			if((loc + 1)%(nilp.nbOne + 1) == 0){
				recvflag = false;
			}
		}

		if(ProcID != 0){
			if(sendflg == true){
				sBuf  = new double [2*size[b]];
				sIndx  = new S [size[b]];
				for(it = dynmat_loc[b].begin(); it != dynmat_loc[b].end(); ++it){
					sBuf[cnt] = (it->second).real();
					sBuf[cnt+1] = (it->second).imag();
					sIndx[cnt2] = it->first;
					cnt = cnt + 2;
					cnt2++;
				}

				MPI_Isend(sBuf, size[b], MPI_SCALAR, up, tag1[b], comm, &Sreqs[b]);
				MPI_Isend(sIndx, size[b], MPI_INT, up, tag2[b], comm, &idxSreqs[b]);

			}
		}

		if(ProcID != nProcs - 1){
			if(recvflag == true){
				rBuf = new double [2*rsize[b]];
				rIndx = new S [2*rsize[b]];
				for(S tt =0; tt < 2*rsize[b]; tt++){
					rBuf[tt] = 0.0;
				}

				for(S tt =0; tt < rsize[b]; tt++){
					rIndx[tt] = 0;
				}
				MPI_Irecv(rBuf,rsize[b], MPI_SCALAR, down, tag1[b], comm, &Rreqs[b]);
				MPI_Irecv(rIndx,rsize[b], MPI_INT, down, tag2[b], comm, &idxRreqs[b]);
				MPI_Wait(&Rreqs[b],&Rstat[b]);
				MPI_Wait(&idxRreqs[b],&idxRstat[b]);
			}
		}

		if(ProcID != nProcs - 1){
			if(recvflag == true){
				for(S tt = 0; tt < rsize[b]; tt++){
						(prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]]).real(rBuf[2*tt]);
						(prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]]).imag(rBuf[2*tt + 1]);
				}
			}
		}

	}

#elif defined (__USE_COMPLEX__)
//complex single

	float *sBuf, *rBuf;

	S *sIndx, *rIndx;
	S cnt = 0, cnt2 = 0;

	for(S b = 0; b < nilp.diagPosition -1; b++){

		if(ProcID != 0){
			loc = y_index_map->Loc2Glob(b);
			q = loc - nilp.diagPosition + 1 + b;
			if((q + 1)%(nilp.nbOne + 1) == 0){
				sendflg = false;
			}
		}

		if(ProcID != nProcs - 1){
			loc = y_index_map->Loc2Glob(nrows - b - 1);
			if((loc + 1)%(nilp.nbOne + 1) == 0){
				recvflag = false;
			}
		}

		if(ProcID != 0){
			if(sendflg == true){
				sBuf  = new float [2*size[b]];
				sIndx  = new S [size[b]];
				for(it = dynmat_loc[b].begin(); it != dynmat_loc[b].end(); ++it){
					sBuf[cnt] = (it->second).real();
					sBuf[cnt+1] = (it->second).imag();
					sIndx[cnt2] = it->first;
					cnt = cnt + 2;
					cnt2++;
				}

				MPI_Isend(sBuf, size[b], MPI_SCALAR, up, tag1[b], comm, &Sreqs[b]);
				MPI_Isend(sIndx, size[b], MPI_INT, up, tag2[b], comm, &idxSreqs[b]);
			}
		}

		if(ProcID != nProcs - 1){
			if(recvflag == true){
				rBuf = new float [2*rsize[b]];
				rIndx = new S [2*rsize[b]];
				for(S tt =0; tt < 2*rsize[b]; tt++){
					rBuf[tt] = 0.0;
				}

				for(S tt =0; tt < rsize[b]; tt++){
					rIndx[tt] = 0;
				}

				MPI_Irecv(rBuf,rsize[b], MPI_SCALAR, down, tag1[b], comm, &Rreqs[b]);
				MPI_Irecv(rIndx,rsize[b], MPI_INT, down, tag2[b], comm, &idxRreqs[b]);
				MPI_Wait(&Rreqs[b],&Rstat[b]);
				MPI_Wait(&idxRreqs[b],&idxRstat[b]);
			}

		}

		if(ProcID != nProcs - 1){
			if(recvflag == true){
				for(S tt = 0; tt < rsize[b]; tt++){
						(prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]]).real(rBuf[2*tt]);
						(prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]]).imag(rBuf[2*tt + 1]);
				}
			}

		}

	}
}

#else
//real

	T *sBuf, *rBuf;
	S *sIndx, *rIndx;
	S cnt = 0;

	for(S b = 0; b < nilp.diagPosition -1; b++){

		if(ProcID != 0){
			loc = y_index_map->Loc2Glob(b);
			q = loc - nilp.diagPosition + 1 + b;
			if((q + 1)%(nilp.nbOne + 1) == 0){
				sendflg = false;
			}
		}

		if(ProcID != nProcs - 1){
			loc = y_index_map->Loc2Glob(nrows - b - 1);
			if((loc + 1)%(nilp.nbOne + 1) == 0){
				recvflag = false;
			}
		}

		if(ProcID != 0){
			if(sendflg == true){
				sBuf  = new T [size[b]];
				sIndx  = new S [size[b]];
				for(it = dynmat_loc[b].begin(); it != dynmat_loc[b].end(); ++it){
					sBuf[cnt] = it->second;
					sIndx[cnt] = it->first;
					cnt++;
			    }
			    MPI_Isend(sBuf, size[b], MPI_SCALAR, up, tag1[b], comm, &Sreqs[b]);
			    MPI_Isend(sIndx, size[b], MPI_INT, up, tag2[b], comm, &idxSreqs[b]);
			}
		}

		if(ProcID != nProcs - 1){
			if(recvflag == true){
				rBuf = new T [rsize[b]];
				rIndx = new S [rsize[b]];
				for(S tt =0; tt < rsize[b]; tt++){
					rBuf[tt] = 0.0;
					rIndx[tt] = 0;
				}
				MPI_Irecv(rBuf,rsize[b], MPI_SCALAR, down, tag1[b], comm, &Rreqs[b]);
				MPI_Irecv(rIndx,rsize[b], MPI_INT, down, tag2[b], comm, &idxRreqs[b]);
				MPI_Wait(&Rreqs[b],&Rstat[b]);
				MPI_Wait(&idxRreqs[b],&idxRstat[b]);
			}
		}

		if(ProcID != nProcs - 1){
			if(recvflag == true){
				for(S tt = 0; tt < rsize[b]; tt++){
						prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]] = rBuf[tt];
				}
			}

		}

	}

#endif

}


#endif
