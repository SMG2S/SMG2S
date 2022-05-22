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

	MatrixCSR(){};

	MatrixCSR(S nnz_in, S nrows_in)
	{
	
		nnz = nnz_in;
		ncols = nrows_in;
		nrows = nrows_in;

		rows.reserve(nrows_in+1);
		cols.reserve(nnz);
		vals.reserve(nnz);	
	
	};

	MatrixCSR(std::vector<S> rowoffs, std::vector<S> colidxs, std::vector<T> values)
	{
	
		nnz = colidxs.size();
		ncols = rowoffs.size() - 1;
		nrows = rowoffs.size() - 1;

		rows = rowoffs;
		cols = colidxs;
		vals = values;	
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

	void show(){
		std::cout << "MatrixCSR display:" << std::endl; 
		for(auto i = 0; i < rows.size() - 1; i++){
			for(auto j = rows[i]; j < rows[i+1]; j++){
				std::cout << "row " << i << ": (" << cols[j] << "," << vals[j]<< ")";
			}
			std::cout << std::endl;
		}
	}

	void Free()
	{
		if(nnz != 0){
			std::vector<S>().swap(rows);
			std::vector<S>().swap(cols);
			std::vector<T>().swap(vals);
			
			nnz = 0;
		}
	};

};

#endif