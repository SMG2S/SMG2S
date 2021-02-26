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

		//Loc set mat diagonal by vector given
		void	Loc_SetDiagonal_index(parVector<T,S> *diag, S index);

		//Mat Scale
		void	MatScale(T scale);

		//Loc Mat Scale
		void	Loc_MatScale(T scale);

		//Loc AXPY
		void	Loc_MatAXPY(parMatrixSparse<T,S> *X, T scale);

		//Loc AYPX
		void    Loc_MatAYPX(parMatrixSparse<T,S> *X, T scale);


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

}

template<typename T, typename S>
parMatrixSparse<T,S>::~parMatrixSparse()
{
	//if index map is defined
	if(x_index_map != NULL){
	        x_index_map->DeleteUser();
   	        delete x_index_map;
	}


	if(y_index_map != NULL){
		y_index_map->DeleteUser();
	    	delete y_index_map;
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
	else{
		return 0.0;
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::glocPlusLloc(){

	S i;
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

	S i;
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
	typename std::map<S,T>::iterator it;

	if(ProcID == 0) {std::cout << "LOC MODE Parallel MatView: " << std::endl;}

	for (i = 0; i < nrows; i++){
		if(dynmat_loc != NULL){
			std::cout << "row " << y_index_map->Loc2Glob(i) << ": ";
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
				if(it->second != T(0)){
					std::cout <<"("<<it->first << "," << it->second << "); ";
				}	
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
	else{
		return 0.0;
	}
}



template<typename T,typename S>
void parMatrixSparse<T,S>::SetDiagonal(parVector<T,S> *diag)
{

	if (nrows != njloc ){
		if(ProcID == 0){
			printf("ERROR: cannot set diagonal for non-square matrix.");
		}
	}
	else{
		T *a = diag->GetArray();
		S local_size = diag->GetArraySize();


		for(S i = 0; i < local_size; i++){
			SetValueLocal(i,diag->Loc2Glob(i),a[i]);
		}

	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetDiagonal(parVector<T,S> *diag)
{

	if (nrows != njloc ){
		if(ProcID == 0){
			printf("ERROR: cannot set diagonal for non-square matrix.");
		}
	}
	else{
		T *a = diag->GetArray();
		S local_size = diag->GetArraySize();


		for(S i = 0; i < local_size; i++){
			Loc_SetValueLocal(i,diag->Loc2Glob(i),a[i]);
		}

	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetDiagonal_index(parVector<T,S> *diag, S index)
{

	if (nrows != njloc ){
		if(ProcID == 0){
			printf("ERROR: cannot set diagonal for non-square matrix.");
		}
	}
	else{
		T *a = diag->GetArray();
		S local_size = diag->GetArraySize();	
		S global_index;

		for(S i = 0; i < local_size; i++){
			global_index = diag->Loc2Glob(i);
			if(index >= 0){
				if(global_index < diag->GetGlobalSize() - index){
					Loc_SetValueLocal(i,global_index - index, a[i]);
				}
			}
			else if(index <= 0){
				if(global_index > -index){
					Loc_SetValueLocal(i, global_index + index, a[i]);
				}
			}
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
void parMatrixSparse<T,S>::Loc_ConvertToCSR(){
	S 	count, i, j;
	T	v;

	typename std::map<S, T>::iterator it;

	if(dynmat_loc != NULL){
		//allocate csr matrix

		CSR_loc = new MatrixCSR<T, S>(nnz_loc, nrows);

		count = 0;

		for(i = 0; i < nrows; i++){
			CSR_loc->rows.push_back(count);
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
				if(it->second != T(0)){
					j = it->first;
					v = it->second;
					CSR_loc->vals.push_back(v);
					CSR_loc->cols.push_back(j);
					count++;
				}
			}
		}
		CSR_loc->rows.push_back(count);
	}
}

///////////////////////////////

template<typename T,typename S>
void parMatrixSparse<T,S>::MatScale(T scale){
	typename std::map<S,T>::iterator it;

	S i;

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

	typename std::map<S,T>::iterator it;

	//use the given nilpotency matrix, MA operation will make elements of matrix right move diaPosition-1 offset.
	//And the positions of 0: pos = nbOne*integer - 1

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
}


template<typename T,typename S>
void parMatrixSparse<T,S>::AM(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod)
{
	S i, j, k, p, q, loc;

	typename std::map<S,T>::iterator it;

	MPI_Request	rtypereq, stypereq;

	MPI_Status	typestat;

	int up, down;
	int	tagtype = nilp.diagPosition - 1 + nilp.diagPosition ;


	S *sIndx, *rIndx;
	S gSize = 0, gRsize = 0;
	S cnt = 0, cnt2 = 0;

	MPI_Request	indSReqs, indRReqs, valSReqs, valRReqs;
	MPI_Status indRStats, valRStats;

	int tagv = 0;
	int tagi = 1;

	up = ProcID - 1;
	down = ProcID + 1;

	S size[nilp.diagPosition - 1], rsize[nilp.diagPosition - 1];

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
		MPI_Isend(size, nilp.diagPosition - 1, getMPI_Type<S>(), up, tagtype, comm, &stypereq);
	}

	if(ProcID != nProcs - 1){
		MPI_Irecv(rsize,nilp.diagPosition - 1, getMPI_Type<S>(), down, tagtype, comm, &rtypereq);
		MPI_Wait(&rtypereq,&typestat);
	}


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

	MPI_Barrier(comm);

	// sending and receving

	if(ProcID != 0){
		for(S a = 0; a < nilp.diagPosition - 1; a++){
			size[a] = dynmat_loc[a].size();
		}
	}

	if(ProcID != 0){
		for(S b = 0; b < nilp.diagPosition - 1; b++){
			gSize += size[b];
		}
	}

	if(ProcID != nProcs - 1){
		for(S b = 0; b < nilp.diagPosition - 1; b++){
			gRsize += rsize[b];
		}
	}


	T *sBuf, *rBuf;

	if(ProcID != 0){
		sBuf  = new T [gSize];
		sIndx  = new S [gSize];
	}

	if(ProcID != nProcs - 1){
		rBuf = new T [gRsize];
		rIndx = new S [gRsize];
	}

	if(ProcID != 0){
		for(S b = 0; b < nilp.diagPosition - 1; b++){
			for(it = dynmat_loc[b].begin(); it != dynmat_loc[b].end(); ++it){
				sBuf[cnt] = it->second;
				sIndx[cnt] = it->first;
				cnt++;
			}
		}

		MPI_Isend(sIndx, gSize, getMPI_Type<S>(), up, tagi, comm, &indSReqs);
		MPI_Isend(sBuf, gSize, getMPI_Type<T>(), up, tagv, comm, &valSReqs);
	}

	if(ProcID != nProcs - 1){
		for(S tt = 0; tt < 2*gRsize; tt++){
			rBuf[tt] = 0.0;
		}

		for(S tt = 0; tt < gRsize; tt++){
			rIndx[tt] = 0;
		}

		MPI_Irecv(rIndx, gRsize, getMPI_Type<S>(), down, tagi, comm, &indRReqs);
		MPI_Irecv(rBuf, gRsize, getMPI_Type<T>(), down, tagv, comm, &valRReqs);
		
		MPI_Wait(&indRReqs,&indRStats);
		MPI_Wait(&valRReqs,&valRStats);

	}

	if(ProcID != nProcs - 1){
		for(S b = 0; b < nilp.diagPosition - 1; b++){
			loc = y_index_map->Loc2Glob(nrows - nilp.diagPosition + 1 + b);
			if((loc + 1)%(nilp.nbOne + 1) != 0){
				if(b == 0){
					for(S tt = 0; tt < rsize[b]; tt++){
						prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]] = rBuf[tt];
					}
				} else{
					for(S tt = rsize[b - 1]; tt < rsize[b] + rsize[b - 1]; tt++){
						prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]] = rBuf[tt];
					}
				}
			}
		}
	}
}

#endif
