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

#include "interface/C/c_wrapper.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int size,rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	NilpInt_t *nilp = newNilpInt();
	nilpIntType1(nilp, 5, 100);
	nilpIntShow(nilp);
	nilpInt_destory(nilp);

	parMatrixSparseDoubleInt_t *mat = newparMatrixSparseDoubleInt();

	char spectrum[] = " ";

	parMatrixSparseDoubleInt_nonsym_smg2s(mat, 100, nilp, 2, spectrum, MPI_COMM_WORLD);
//	parMatrixSparseDoubleInt_LocMatView(mat);
/*
	int rs, cs;
	GetLocalSizeDoubleInt(mat, &rs, &cs);

	printf("rank %d: rs = %d, cs = %d\n", rank, rs, cs);

	int size1, size2;

	parMatrixSparseDoubleInt_LocToCSR(mat);

	parMatrixSparseDoubleInt_LocGetCSRSize(mat, &size1, &size2);

	printf("rank %d: size1 = %d, size2 = %d\n", rank, size1, size2);

	int *rowoffsets = (int *)malloc(size1 * sizeof(int));
	int *colinds = (int *)malloc(size2 * sizeof(int));
	double *vals = (double *)malloc(size2 * sizeof(double));

	parMatrixSparseDoubleInt_LocGetCSRArrays(mat, size1, size2, &rowoffsets, &colinds, &vals);			
*/
	parMatrixSparseDoubleInt_destory(mat);

	MPI_Finalize();
}


