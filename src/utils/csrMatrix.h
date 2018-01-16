#ifndef __CSRMATRIX_H__

#define __CSRMATRIX_H__

#include <iostream>

using namespace std;

template <typename T, typename S>
class csrMatrix
{
	public:
		S row, column, nnz;
		S *csrRowPtr, *csrColInd;
		T *csrVal;
		
};

template <typename T, typename S>
class d_csrMatrix
{
        public:
                S d_row, d_column, d_nnz;
                S *d_csrRowPtr, *d_csrColInd;
                T *d_csrVal;

};

#endif
