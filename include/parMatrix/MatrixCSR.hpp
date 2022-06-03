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

#ifndef __MATRIXCSR_H__
#define __MATRIXCSR_H__

#include <vector>

/** @defgroup group6 Parallel Vector and Matrix
 *  Implementation of paralel vector and matrix distributed across 1D MPI grids.
 *  @{
 */
//!  @brief A struct which stores a sparse matrix in CSR format.
/*!
  This struct defines a CSR format to store a sparse matrix. This format stores a sparse matrix
  with three arrays
  
  - `rows` points to row starts in indices and data
  - `cols` is array of column indices
  - `val` is array of corresponding nonzero values

  @tparam T describes the scalar types of the entries of a sparse matrix.   
  @tparam S type of integer to describes the dimension of matrices to be generated. 
*/
template<typename T, typename S>
struct MatrixCSR
{
 	//! number of rows of a sparse matrix
    /*!
      number of rows of a sparse matrix
    */	
	S	nrows;
	//! number of non-zero entries of a sparse matrix
    /*!
      number of non-zero entries of a sparse matrix
    */		
	S	nnz;
	//! dimension of the array MatrixCSR#cols
	/*!
	  dimension of the array MatrixCSR#cols, which stores the column indices of all non-zero entries.
	  
	  - It equals to MatrixCSR#nnz
	*/
	S   ncols;

	//! An array points to row starts in indices and data
    /*!
      Its size is `nrows+1`
      - non-zero entry of the `i`-th row are `vals[rows[i]:rows[i+1]]`, with column indices `cols[rows[i]:rows[i+1]]`
      - item `(i,j)` can be accessed as `vals[rows[i]+k]`, where `k` is the position of `j` in `cols[rows[i]:rows[i+1]]`
    */		
	std::vector<S> rows;
	//! An array stores the column indices of non-zero elements
    /*!
      Its size is `nnz`
      - non-zero entry of the `i`-th row are `vals[rows[i]:rows[i+1]]`, with column indices `cols[rows[i]:rows[i+1]]`
      - item `(i,j)` can be accessed as `vals[rows[i]+k]`, where `k` is the position of `j` in `cols[rows[i]:rows[i+1]]`
    */		
	std::vector<S> cols;
	//! An array stores the values of non-zero elements
    /*!
      Its size is `nnz`
      - non-zero entry of the `i`-th row are `vals[rows[i]:rows[i+1]]`, with column indices `cols[rows[i]:rows[i+1]]`
      - item `(i,j)` can be accessed as `vals[rows[i]+k]`, where `k` is the position of `j` in `cols[rows[i]:rows[i+1]]`
    */		
	std::vector<T> vals;

	MatrixCSR(){};

    //! A constructor of `MatrixCSR`. The three arrays MatrixCSR#rows, MatrixCSR#cols and MatrixCSR#vals are allocated, without further filling in.
    /*!
      * @param[in] nnz_in number of non-zero elements in the sparse matrix to be allocated
      * @param[in] nrows_in number of rows of the sparse matrix to be allocated
    */	
	MatrixCSR(S nnz_in, S nrows_in)
	{
	
		nnz = nnz_in;
		ncols = nrows_in;
		nrows = nrows_in;

		rows.reserve(nrows_in+1);
		cols.reserve(nnz);
		vals.reserve(nnz);	
	
	};
    //! A constructor of `MatrixCSR`. The three arrays MatrixCSR#rows, MatrixCSR#cols and MatrixCSR#vals are filled with user-provided arrays.
    /*!
      * @param[in] rowoffs user-provided array stored in a `std::vector` for MatrixCSR#rows
      * @param[in] colidxs user-provided array stored in a `std::vector` for MatrixCSR#cols
      * @param[in] values user-provided array stored in a `std::vector` for MatrixCSR#vals
    */	
	MatrixCSR(std::vector<S> rowoffs, std::vector<S> colidxs, std::vector<T> values)
	{
	
		nnz = colidxs.size();
		ncols = rowoffs.size() - 1;
		nrows = rowoffs.size() - 1;

		rows = rowoffs;
		cols = colidxs;
		vals = values;	
	};

	//! This a destructor of MatrixCSR
	~MatrixCSR()
	{
		if(nnz != 0){
			rows.empty();
			cols.empty();
			vals.empty();
		}
	};

	//! Get the value of position `(row, col)` from sparse matrix
    /*!
      * @param[in] row index of row of the poisition in sparse matrix
      * @param[in] col index of column of the poisition in sparse matrix
    */		
	T GetValue(S row, S col)
	{
		T val;
		S currCol;

		for(S pos = rows[row - 1] - 1; pos < rows[row] - 1; ++pos){
			currCol = cols[pos];
			if (currCol == col){
				val = vals[pos];
			} else{
				val = 0.0;
			}
		}
		return val;
	};


	//! Overwrite the value of entry of position  `(row, col)` in sparse matrix by `val`
    /*!
      * @param[in] row index of row of the poisition in sparse matrix
      * @param[in] col index of column of the poisition in sparse matrix
      * @param[in] val the value should be set on the position `(row, col)` of sparse matrix
    */	
	void SetValue(S row, S col, T val)
	{
		S currCol;

		for(S pos = rows[row - 1] - 1; pos < rows[row] - 1; ++pos){
			currCol = cols[pos];
			if (currCol == col){
				vals[pos] = val;
			} 
		}
	};

	//! Insert a value `val` in the position `(row, col)` in sparse matrix
    /*!
      - unlike MatrixCSR<T>#SetValue, in this function, the position `(row, col)` can be anyone which is still not available in the matrix.
      * @param[in] row index of row of the poisition in sparse matrix
      * @param[in] col index of column of the poisition in sparse matrix
      * @param[in] val the value should be set on the position `(row, col)` of sparse matrix
    */		
	void InsertValue(S row, S col, T val)
	{
		S currCol, prevCol;
		typename std::vector<S>::iterator it=cols.begin();
		typename std::vector<T>::iterator itv=vals.begin();

		for(S pos = rows[row - 1] - 1; pos < rows[row] - 1; ++pos){
			currCol = cols[pos];
			prevCol = cols[pos-1];
			if (currCol == col){
				SetValue(row, col, val);
			} else if((currCol > col) && (prevCol < col)){	
				std::advance(it,currCol);
				std::advance(itv,currCol);
				cols.insert(it,col);
				vals.insert(itv, val);
				nnz++;
				return;
			}
		}
	};

	//! Print the sparse matrix in COO format
	/*!
		Print the sparse matrix in COO format
	 */
	void show(){
		std::cout << "MatrixCSR display:" << std::endl; 
		for(auto i = 0; i < rows.size() - 1; i++){
			for(auto j = rows[i]; j < rows[i+1]; j++){
				std::cout << "row " << i << ": (" << cols[j] << "," << vals[j]<< ")";
			}
			std::cout << std::endl;
		}
	}

};
/** @} */ // end of group6
#endif