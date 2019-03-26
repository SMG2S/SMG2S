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
#include <mpi.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

	int size,rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*create Nilpotency object*/
	struct NilpotencyInt *n;
  /*create Instance*/
	n = newNilpotencyInt();
  /*setup Nilpotency object*/
	NilpType1(n, 2, 10);
	if(rank == 0){
		showNilpotencyInt(n);
	}

  /*Create parMatrixSparse Object*/
	struct parMatrixSparseRealDoubleInt *m;
  /*create Instance*/
	m = newParMatrixSparseRealDoubleInt();
  /*Generate by SMG2S*/
	smg2sRealDoubleInt(m, 10, n, 3 ," ",MPI_COMM_WORLD);
  /*Matrix View*/
	LOC_MatViewRealDoubleInt(m);

  /*Release Nilpotency Object
  Release parMatrixSparse Object*/
	ReleaseNilpotencyInt(&n);
	ReleaseParMatrixSparseRealDoubleInt(&m);

	MPI_Finalize();
	return 0;
}


