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

   Part of basic data structures implementation of this file refers to this technical report 
   (http://orbit.dtu.dk/files/51272329/tr12_10_Alexandersen_Lazarov_Dammann_1.pdf)
*/

#ifndef __MATRIXCSR_H__
#define __MATRIXCSR_H__

#include <vector>

template<typename T, typename S>
struct MatrixCSR
{
	S	nrows;
	S	nnz;
	S   ncols;

	std::vector<S> rows;
	std::vector<S> cols;
	std::vector<T> vals;

	MatrixCSR()
	{
		nnz = 0;
//		rows = NULL;
//		cols = NULL; 
//		vals = NULL;
	};

	MatrixCSR(S nnz_in, S nrows_in)
	{
	
		nnz = nnz_in;
		ncols = nrows_in;
		nrows = nrows_in;

		rows.reserve(nrows_in+1);
		cols.reserve(nnz);
		vals.reserve(nnz);	
	
	};

	~MatrixCSR()
	{
		if(nnz != 0){
			rows.empty();
			cols.empty();
			vals.empty();
		}
	};

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
//				it = currCol;
//				itv = currCol;		
				std::advance(it,currCol);
				std::advance(itv,currCol);
				cols.insert(it,col);
				vals.insert(itv, val);
				nnz++;
				return;
			}
		}
	};

	void Add(MatrixCSR<T,S> m){
		S row, nnz;
		T v1, v2, v;
		
		if (nrows != m.nrows || ncols != m.ncols) {
			printf("Cannot add: matrices dimensions don't match.\n");
			return;
		}

		for (S i = 1; i <= nrows; i++){
			for (S j = 1; j <= ncols; j++) {
				v1 = GetValue(i,j);
				v2 = m.GetValue(i,j);
				v = v1 + v2;
				if(v != 0.0){
					InsertValue(i,j,v1+v2);
				}
			}
		}

	};


	void Free()
	{
		if(nnz != 0){
			std::vector<S>().swap(rows);
			std::vector<S>().swap(cols);
			std::vector<T>().swap(vals);
			
			nnz = 0;
//			rows = NULL;
//			cols = NULL;
//			vals = NULL;
		}
	};

};

#endif
