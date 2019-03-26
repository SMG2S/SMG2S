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

#ifndef __PETSC_INTERFACE_H__
#define __PETSC_INTERFACE_H__

#include "petscmat.h"
#include "../../smg2s/smg2s.h"
#include <complex>
#include <map>


//#if defined (PETSC_USE_COMPLEX) && defined (PETSC_USE_64BIT_INDICES)
//complex long int

#if defined (PETSC_USE_COMPLEX)
//complex int

	Mat ConvertToPETSCMat(parMatrixSparse<std::complex<double>, int > *M){

		int rank, size;
		std::vector<int>::iterator it;
		std::vector<std::complex<double> >::iterator itm;

		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		int g, b;
		PetscInt m, n;
		PetscInt *i, *j;
		PetscScalar *a;

		Mat A;

		M->GetLocalSize(g, b);

		m = (PetscInt)g;
		n = (PetscInt)m;


		int count = 0;
		int count2 = 0;

		for(it = M->CSR_loc->rows.begin(); it != M->CSR_loc->rows.end(); ++it){
				count++;
		}

		for(itm = M->CSR_loc->vals.begin(); itm != M->CSR_loc->vals.end(); ++itm){
			if(itm->imag() !=0 || itm->real() != 0){
				count2++;
			}
		}

		PetscMalloc1((PetscInt)count, &i);
		PetscMalloc1((PetscInt)count2, &j);
		PetscMalloc1((PetscInt)count2, &a);

		for(PetscInt q =0; q < count; q++){
			i[q] = M->CSR_loc->rows[q];
		}

		for(PetscInt q = 0; q < count2; q++){
			j[q] = M->CSR_loc->cols[q];
			a[q] = M->CSR_loc->vals[q].real() + PETSC_i*M->CSR_loc->vals[q].imag();
			}

		MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, count-1, count-1, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

		return A;
	}

#elif defined (PETSC_USE_REAL_DOUBLE) && defined (PETSC_USE_64BIT_INDICES)
//real long int double
#elif defined (PETSC_USE_REAL_SINGLE) && defined (PETSC_USE_64BIT_INDICES)
//real long int single

#elif defined (PETSC_USE_REAL_DOUBLE)

	Mat ConvertToPETSCMat(parMatrixSparse<double, int > *M){

	int rank, size;
		typename std::vector<int>::iterator it;
		typename std::vector<std::complex<double> >::iterator itm;

		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		int g, b, c, d;
		PetscInt m, n;
		PetscInt *i, *j, ii;
		PetscScalar *a;

		Mat A;

		M->GetLocalSize(g, b);

		m = (PetscInt)g;
		n = (PetscInt)m;

		int count = 0;
		int count2 = 0;

		for(it = M->CSR_loc->rows.begin(); it != M->CSR_loc->rows.end(); ++it){
				count++;
		}

		for(itm = M->CSR_loc->vals.begin(); itm != M->CSR_loc->vals.end(); ++itm){
			if(itm != 0){
				count2++;
			}
		}

		PetscMalloc1((PetscInt)count, &i);
		PetscMalloc1((PetscInt)count2, &j);
		PetscMalloc1((PetscInt)count2, &a);

		for(PetscInt q =0; q < count; q++){
			i[q] = M->CSR_loc->rows[q];
		}

		for(PetscInt q = 0; q < count2; q++){
			j[q] = M->CSR_loc->cols[q];
			a[q] = M->CSR_loc->vals[q];
			}

		MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, count-1, count-1, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

		return A;

	}

//real double int
#elif defined (PETSC_USE_REAL_SINGLE)
//real single int

#endif



#endif
