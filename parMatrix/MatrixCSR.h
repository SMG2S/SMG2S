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
