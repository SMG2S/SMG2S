#include "csrMatrix.h"
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;

int main()
{
	csrMatrix<double,int> A;
	A.row = 4;
	A.column = 4;
	A.nnz = 9;

	A.csrRowPtr = (int *)malloc((A.row+1)*sizeof(int));
        A.csrColInd = (int *)malloc(A.nnz*sizeof(int));
        A.csrVal = (double *)malloc(A.nnz*sizeof(double));
        
	std::vector<int> rp = {0,2,4,7,9};
	std::vector<int> ci = {0,1,1,2,0,2,3,1,3};
	std::vector<double> val = {1,7,2,8,5,3,9,6,4};


	int i;

	for(i = 0; i < A.row+1; i++){
		A.csrRowPtr[i] = rp[i];
		cout << "rP[" << i << "]=" << A.csrRowPtr[i] << ";\n";
	}


	for(i = 0; i < A.nnz; i++){
		A.csrColInd[i] = ci[i];
                cout << "cI[" << i << "]=" << A.csrColInd[i] << ";\n";
	}

        for(i = 0; i < A.nnz; i++){

		A.csrVal[i] = val[i];
		cout << "val[" << i << "]=" << A.csrVal[i] << ";\n";
        }

return 0;

}


