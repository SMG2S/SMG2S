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
#include <parVector/parVector.hpp>
#include <smg2s/nilpotent.hpp>
#include <smg2s/spectrum.hpp>
#include <utils/MPI_DataType.hpp>
#include <utils/utils.hpp>
#include <parMatrix/MatrixCSR.hpp>

#ifdef __USE_COMPLEX__
#include <complex>
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//!  A class which defines a sparse matrix distributed across 1D MPI grid.
/*!
 * @ingroup group6 
  - This class can be constructed with a given parVector with the same distribution scheme.
  - This class can be constructed with a distribution scheme by a given parVectorMap object.

  @tparam T describes the scalar types of the entries of a sparse matrix.   
  @tparam S type of integer to describes the dimension of vector to be generated. 
*/
template<typename T, typename S>
class parMatrixSparse
{
    private:
    //! An array of `std::map` in which each map stores the column index and related non-zero entry value of each row.
	std::map<S, T> *dynmat_loc;

	//size of local matrix
	//! number of columns of local matrix on each MPI proc
	S	ncols;
	//! number of rows of local matrix on each MPI proc
	S   nrows;
	//number of non-zeros in local
	//! number of non-zero entries of local matrix on each MPI proc
	S	nnz_loc;
	//! the parVectorMap object used to distribute the global matrix across 1D MPI grid
	parVectorMap<S>	index_map;
	//! the smallest index of row of a distributed matrix on each MPI proc
	S	lower_b;
	//! upper_b-1 = the largest index of row of a distributed matrix on each MPI proc
	S   upper_b;
	//! the working MPI communicator
	MPI_Comm comm;

	// mpi size and rank
	//! rank of each MPI procs within the working MPI communicator	
	int ProcID;
	//! number of MPI procs within the working MPI communicator
	int nProcs;

    public:
    	parMatrixSparse();
	    //! A constructor of `parMatrixSparse`. 
	    /*!
	      * @param[in] vec a given parVector, a parMatrixSparse object is constructed with the same distribution scheme of `vec`

	      - parMatrixSparse::dynmat_loc is not allocated. 
	    */	    	
    	parMatrixSparse(parVector<T,S> vec);
	    //! A constructor of `parMatrixSparse`. 
	    /*!
	      * @param[in] map the distribution scheme determined by this object of type parVectorMap

	      - parMatrixSparse::dynmat_loc is not allocated. 
	    */	    	
    	parMatrixSparse(parVectorMap<S> map);

    	//get
    	//! Return the number of rows of local matrix on each MPI proc
    	S GetNRows(){return nrows;};
    	//! Return the number of columns of local matrix on each MPI proc    	
    	S GetNCols(){return ncols;};
    	//! Return the number of non-zeros entries of local matrix on each MPI proc    	    	
    	S GetNNzLoc(){return nnz_loc;};
    	//! Return the parVectorMap object used to distribute sparse Matrix 
    	parVectorMap<S> GetMap(){return index_map;};
    	//! Return parMatrixSparse#lower_b
    	S GetLowerBound(){return lower_b;};
    	//! Return parMatrixSparse#upper_b    	
    	S GetUpperBound(){return upper_b;};
    	//! Return parMatrixSparse#comm    	
    	MPI_Comm GetComm(){return comm;};
    	//! Return parMatrixSparse#ProcID    	    	
    	int GetProcId(){return ProcID;};
    	//! Return parMatrixSparse#nProcs    	    	    	
    	int GetNProcs(){return nProcs;};

    	//! Return parMatrixSparse#dynmat_loc 
		std::map<S, T> *GetDynMatLoc(){return dynmat_loc;};

		//set value
		//! Set a value to a distributed sparse matrix with a given local index on each MPI proc
	    /*!
	      * @param[in] row the local index of row
	      * @param[in] col the local index of column
	      * @param[in] value the scalar to be set
	    */			
		void	SetValueLocal(S row, S col, T value);
		//set value
		//! Set multiple values to a distributed sparse matrix with given multiple local indices on each MPI proc
	    /*!
	      * @param[in] nindex number of values to be set
	      * @param[in] rows an array storing multiple local indices of row
	      * @param[in] cols an array storing multiple local indices of column
	      * @param[in] values an array storing multiple scalars to be set
	    */			
		void	SetValuesLocal(S nindex, S *rows, S *cols, T *values);
		//set value
		//! Set a value to a distributed sparse matrix with a given global index on each MPI proc
	    /*!
	      * @param[in] row the global index of row
	      * @param[in] col the global index of column
	      * @param[in] value the scalar to be set
	    */			
		void 	SetValue(S row, S col, T value);
		//! Set the diagonal of a distributed sparse matrix with a given parVector object
	    /*!
	      * @param[in] diag a given parVector object to be set on the diagonal
	    */			
		void    SetDiagonal(parVector<T,S> diag);

		//add value
		//! Add a value to an entry of a distributed sparse matrix with a given local index on each MPI proc
	    /*!
	      * @param[in] row the local index of row
	      * @param[in] col the local index of column
	      * @param[in] value the scalar to be added
	    */		
		void    AddValueLocal(S row, S col, T value);
		//! Add a value to an entry of a distributed sparse matrix with a given global index on each MPI proc
	    /*!
	      * @param[in] row the global index of row
	      * @param[in] col the global index of column
	      * @param[in] value the scalar to be added
	    */			
		void 	AddValue(S row, S col, T value);

		//get value
		//! Return a value of an entry of a distributed sparse matrix with a given local index on each MPI proc
	    /*!
	      * @param[in] row the local index of row
	      * @param[in] col the local index of column
	    */			
		T 	GetValueLocal(S row, S col);
		//! Return a value to an entry of a distributed sparse matrix with a given global index on each MPI proc
	    /*!
	      * @param[in] row the global index of row
	      * @param[in] col the global index of column
	    */				
		T 	GetValue(S row, S col);

		//! Perform an operation `A = scale * A`, in which `scale` is a scalar
	    /*!
	      * @param[in] scale the scalar to be multiplied on the distributed sparse matrix
	    */			
		void	MatScale(T scale);
		//! Perform `A = A + scale * X`
	    /*!
	      * @param[in] X another parMatrixSparse object with a same distribution scheme
	      * @param[in] scale the scalar to be multiplied on `X`
	    */			
		void    MatAXPY(parMatrixSparse<T,S> X, T scale);
		//! Perform `A = scale * A + X`
	    /*!
	      * @param[in] X another parMatrixSparse object with a same distribution scheme
	      * @param[in] scale the scalar to be multiplied on `A`
	    */			
		void    MatAYPX(parMatrixSparse<T,S> X, T scale);

		//init matrix lower part of matrix 
		//! filling the lower part of matrix between diagonal of offset `diag_l` and diagonal of offset `diag_u` with random values `scale * rnd + shift` where `rnd` is a randomly generated value between `0` and `1`
	    /*!
	      * @param[in] diag_l the offset of lower diagonal 
	      * @param[in] diag_u the offset of lower diagonal 
	      * @param[in] scale a scalar to be multiplied on the randomly generated value 
	      * @param[in] shift a scalar to be added on the randomly generated value 
	      * @param[in] sparsity the probability that a entry in the range set to be zero
	    */			
		void	initMat(S diag_l, S diag_u, Base<T> scale, T shift, Base<T> sparsity );
		//! filling the lower part of matrix between diagonal of offset `diag_l` and diagonal of offset `diag_u` with random values `rnd`, which is a randomly generated value between `0` and `1`
	    /*!
	      * @param[in] diag_l the offset of lower diagonal 
	      * @param[in] diag_u the offset of lower diagonal 
	    */			
		void	initMat(S diag_l, S diag_u);
	   	//! Set a given spectrum onto a parMatrixSparse object for constructing a non-Hermitian matrix
		/*!
	      * @param[in] spectrum a parVector object which stores the given spectrum in complex scalars

		  - This operation is naturally in parallel since the spectrum is stored in a parVector object

	      - Attention, this member function works only with complex scalar
		*/		
		void setSpecNonHerm(parVector<T,S> spectrum);
	   	//! Set a given spectrum (all eigenvalues are in real scalars) onto a parMatrixSparse object for constructing a non-Symmetric matrix
		/*!
	      * @param[in] spectrum a parVector object which stores the given spectrum in real scalars

		  - This operation is naturally in parallel since the spectrum is stored in a parVector object
		  
	      - Attention, this member function works only with real scalar
		*/			
		void setSpecNonSymm(parVector<T,S> spectrum);
	   	//! Set a given spectrum (all eigenvalues can be in real scalars or appear as pairs of conjugate complex scalars) onto a parMatrixSparse object for constructing a non-Symmetric matrix
		/*!
	      * @param[in] spectrum a parVector object which stores the given spectrum with conjugate eigenvalues

		  - This operation is naturally in parallel since the spectrum is stored in a parVector object
		  
	      - Attention, for each conjugate pairs of eigenvalues, they should be placed one after another, they cannot be placed in a random order.
	        Before setting the spectrum to the matrix, the spectrum will be checked if it satisfies this condition
		*/			
		void setSpecNonSymmCmplx(parVector<std::complex<Base<T>>,S> spectrum);

		//! Explicitly re-compute the number of nnz on each MPI proc
		/*!
		  - This member function is introduced in case any matrix operation would failed to update the number of nnz 
		*/
		void    updateNnz();

		//! Duplicate from another parMatrixSparse object
		/*!
	      * @param[in] X another parMatrixSparse object to be duplicated with a same distribution scheme
		*/
		void    copy(parMatrixSparse<T,S> X);
		//! Explicitly remove all the zeros of parMatrixSparse object
		/*!
		  - This member function is introduced in case any matrix operation would introduce some (explicit) zeros.
		*/
		void rmZeros();

	   	//matrix multiplies a nilpotent matrix
	   	//! Perform `M*A`, in which `A` is a nilpotent matrix
		/*!
	      * @param[in] nilp a Nilpotent<S> object which determines a nilpotent matrix
		*/	   	
		parMatrixSparse<T,S>	MA(Nilpotent<S> nilp);

	   	//a nilpotent matrix multiplies a matrix
	   	//! Perform `A*M`, in which `A` is a nilpotent matrix
		/*!
	      * @param[in] nilp a Nilpotent<S> object which determines a nilpotent matrix
		*/		   	
		parMatrixSparse<T,S>	AM(Nilpotent<S> nilp);

		//! Make all the entries a sparse matrix to be zeros, but keep the sparsity structure as it is
		void	ZeroEntries();

		//! Convert a parMatrixSparse with dynamic memory into a distributed CSR matrix
		/*!
		  - Attention, after the converting, any operations on the matrix can be implicity updated on the CSR matrix
		*/		
		MatrixCSR<T,S>	ConvertToCSR();

		//! Display multiple information of a parMatrixSparse object
		/*!
		  - number of rows of local matrix
		  - number of columns of local matrix
		  - number of nnz of local matrix
		  - parMatrixSparse#lower_b
		  - parMatrixSparse#upper_b
		*/		
    	void show();
		//! Print a parMatrixSparse object in a distributed COO format
		/*!
		  - This is a distributed function that each MPI proc can only display the piece of local matrix on itself.
		*/	    	
    	void MatView();
                //! Print a parMatrixSparse object in a distributed COO format with matrix name indicated at the starting point of each line
                /*!
                  - This is a distributed function that each MPI proc can only display the piece of local matrix on itself.
                */
	void MatView(std::string matName);
    	//! A parallel IO to write a parMatrixSparse object into a file of MatrixMarket format
   		/*!
	      * @param[in] file_name the path and file name to write into

	      - Attention, this method works only of sparse matrix with real scalar (`double`, `float`...)
		*/ 	
    	void writeToMatrixMarket(std::string file_name);
    	//! A parallel IO to write a parMatrixSparse object with complex scalar into a file of MatrixMarket format
   		/*!
	      * @param[in] file_name the path and file name to write into

	      - Attention, this method works only of sparse matrix with scalar scalar (std::complex<double>, std::complex<float>...)
		*/ 	    	
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

template<typename T, typename S>
void parMatrixSparse<T,S>::MatView(std::string matName)
{
    typename std::map<S,T>::iterator it;
    if(ProcID == 0) {std::cout << "Parallel MatView: " << std::endl;}

    for (S i = 0; i < nrows; i++){
        if(dynmat_loc != NULL){
            std::cout <<  matName + " ]> row " << index_map.Loc2Glob(i) << ": ";
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
    	return GetValueLocal(local_row, col);
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
void parMatrixSparse<T,S>::rmZeros(){
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
            
    for (auto i = MAX(diag_u, lower_b); i < upper_b; i++){
    	for (auto j = MAX(0, i - diag_l); j <= i - diag_u; j++){
    	    auto sampling = d(rd);
    	    if (sampling < 1.0 - sparsity){
    	    	SetValue(i, j, T(scale * d(rd)) + shift);
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

    std::vector<S> targetProcs;
    std::vector<S> originProcs;
    std::vector<S> procsDiff;
    
    typename std::map<S,T>::iterator it;

    if(prod.dynmat_loc == NULL){
    	prod.dynmat_loc = new std::map<S,T> [nrows];
    }
    

    //shift the in-memory part
    if(dynmat_loc != NULL){
    	for(auto row = lprocbound_map[ProcID] + offset; row < uprocbound_map[ProcID]; row++){
	    auto i = index_map.Glob2Loc(row);
    	    for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
    	    	auto j =  it->first;
    	    	    //if not the index of zeros in nilpotent
    	    	if (std::find(IndOfZeros.begin(), IndOfZeros.end(), row-offset-offset ) == IndOfZeros.end()&& it->second != T(0) ){
    	    	    prod.SetValueLocal(i-offset, j, it->second);
    	    	}
    	    }	
    	}
    }

    std::map<std::pair<int, int>, std::vector<S>> rowSendToProc;
    typename std::map<std::pair<int, int>, std::vector<S>>::iterator itm;

    for(auto id = 0; id < nProcs; id++){
	for(auto rid = 0; rid < id; rid++){   
             std::vector<S> rcollect;	       	
             for(auto r = lprocbound_map[id]; r < lprocbound_map[id] + MIN(offset, uprocbound_map[id] - lprocbound_map[id]); r++){
                if( (r - offset >= lprocbound_map[rid]) && (r - offset < uprocbound_map[rid]) ){
                    if (std::find(IndOfZeros.begin(), IndOfZeros.end(), r - offset - offset) == IndOfZeros.end() ){
		        rcollect.push_back(r);
                    }
                }
            }
            if(rcollect.size() != 0){
                rowSendToProc.insert(std::pair<std::pair<int, int>, std::vector<S>>( std::pair<int, S>(id, rid), rcollect ) );
            }	     
        }
    }

    std::vector<std::map<std::pair<int, int>,std::vector<S>>> rowSendToProc_Path;
    while(rowSendToProc.size() != 0){
        std::map<std::pair<int, int>,std::vector<S>> single_path;
        std::vector<int> first_collected;
	std::vector<int> second_collected;
	for(itm = rowSendToProc.begin(); itm != rowSendToProc.end(); ++itm){
	    if(std::find(first_collected.begin(), first_collected.end(), itm->first.first) == first_collected.end()){
		if(std::find(second_collected.begin(), second_collected.end(), itm->first.second) == second_collected.end()){
	            first_collected.push_back(itm->first.first);
		    second_collected.push_back(itm->first.second);
                    single_path.insert(std::pair<std::pair<int, int>,std::vector<S>>(itm->first, itm->second) );

		}
	    }
	}

	for(auto i = 0; i < first_collected.size(); i++){
	    rowSendToProc.erase(std::pair<int, int>(first_collected[i], second_collected[i]) );
	}

	rowSendToProc_Path.push_back(single_path);
    }
  
    if(ProcID == 0){
	for(auto i = 0; i < rowSendToProc_Path.size();i++){
		std::cout << "Sending path #" << i << ": ";
		for(itm = rowSendToProc_Path[i].begin(); itm != rowSendToProc_Path[i].end(); ++itm){
		    std::cout << itm->first.first << " ===> " << itm->first.second << "     ";
		}
		std::cout << std::endl;
    	} 
    }
    
    int comm_path_nb = rowSendToProc_Path.size();
    //package for sending and receving

    for(auto path = 0; path < comm_path_nb; path++){
    	
	std::vector<S> sSize;
    	std::vector<S> sIndices;
    	std::vector<T> sBufs;

        std::vector<S> rSize;
        std::vector<S> rIndices;
        std::vector<T> rBufs;

        for(itm = rowSendToProc_Path[path].begin(); itm != rowSendToProc_Path[path].end(); ++itm){
	     MPI_Request stypereq;
	     MPI_Request rtypereq;
	     MPI_Status  typestat;
	     if(ProcID == itm->first.first){
	         S count = 0;
		 sIndices.insert(sIndices.begin(), count);
		 for(auto i = 0; i < itm->second.size(); i++){
		     S j = index_map.Glob2Loc(itm->second[i]);
		     for(it = dynmat_loc[j].begin(); it != dynmat_loc[j].end(); ++it){
		          sIndices.push_back(it->first);
			  sBufs.push_back(it->second);
			  count++;
		     }
		     sIndices.insert(sIndices.begin(), count);
		 }
		 sSize.resize(2);
		 sSize[0] = sBufs.size();
		 sSize[1] = sIndices.size() - sBufs.size() - 1;

		 MPI_Isend(sSize.data(), 2, getMPI_Type<S>(), itm->first.second, itm->first.second, comm, &stypereq);

	     }

	     if(ProcID == itm->first.second){
	         rSize.resize(2);
                 MPI_Irecv(rSize.data(), 2, getMPI_Type<S>(), itm->first.first, ProcID, comm, &rtypereq);
                 MPI_Wait(&rtypereq, &typestat);

		 rBufs.resize(rSize[0]);
		 rIndices.resize(rSize[1]+rSize[0]+1);
	     }

	     if(ProcID == itm->first.first){
	         MPI_Isend(sBufs.data(), sSize[0], getMPI_Type<T>(), itm->first.second, itm->first.second, comm, &stypereq);
	     }

	     if(ProcID == itm->first.second){
     	    	MPI_Irecv(rBufs.data(), rSize[0], getMPI_Type<T>(), itm->first.first, ProcID, comm, &rtypereq);	
    	    	MPI_Wait(&rtypereq, &typestat);
	     }

             if(ProcID == itm->first.first){
                 MPI_Isend(sIndices.data(), sSize[0]+sSize[1]+1, getMPI_Type<S>(), itm->first.second, itm->first.second, comm, &stypereq);
             }

             if(ProcID == itm->first.second){
                MPI_Irecv(rIndices.data(), rSize[0]+rSize[1]+1, getMPI_Type<S>(), itm->first.first, ProcID, comm, &rtypereq);
                MPI_Wait(&rtypereq, &typestat);
     
                for(auto i = 0; i < itm->second.size(); i++){
		    auto row = index_map.Glob2Loc(itm->second[i] - offset);
		    for(auto cnt = rIndices[rSize[1] - i]  ; cnt < rIndices[rSize[1] - i - 1]; cnt++){
		        auto col = rIndices[rSize[1]+1+cnt];
			prod.SetValueLocal(row, col, rBufs[cnt]);
		    }
		}


	     }
	}
    }
    
    return prod;
}

template<typename T,typename S>
void parMatrixSparse<T,S>::writeToMatrixMarket(std::string file_name){

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
