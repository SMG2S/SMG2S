#include "petscmat.h"
#include "../../parMatrix/parMatrixSparse.cc"
#include <complex>


//#if defined (PETSC_USE_COMPLEX) && defined (PETSC_USE_64BIT_INDICES)
//complex long int

#if defined (PETSC_USE_COMPLEX)
//complex int

	//complex long int
	Mat ConvertToPETSCMat(parMatrixSparse<std::complex<double>, int > *M){

		int rank, size;

		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		int g, b, c, d, e, f;
		PetscInt m, n;
		PetscInt *i, *j, *oi, *oj;
		PetscScalar *a, *oa;

		Mat A;

		M->GetLocalSize(g, b);

		m = (PetscInt)g;
		n = (PetscInt)b;

		if(M->CSR_lloc != NULL){
			c = M->CSR_lloc->rows.size();
			d = M->CSR_lloc->cols.size();
			PetscMalloc1((PetscInt)c, &i);
			PetscMalloc1((PetscInt)d, &j);
			PetscMalloc1((PetscInt)d, &a);

			//printf("c = %d, d = %d\n",c, d);
			for(PetscInt q =0; q < c; q++){
				i[q] = M->CSR_lloc->rows[q];
				//printf("Proc: %d, i[%d] = %d\n", rank, q, i[q]);
			}

			for(PetscInt q = 0; q < d; q++){
				j[q] = M->CSR_lloc->cols[q];
				a[q] = M->CSR_lloc->vals[q].real() + PETSC_i * M->CSR_lloc->vals[q].imag();
				//printf("Proc: %d, j[%d] = %d, a[%d] = %f+%fi\n", rank, q, j[q],q, M->CSR_lloc->vals[q].real(), M->CSR_lloc->vals[q].imag());
			}
		}

		if(M->CSR_gloc != NULL){
			e = M->CSR_gloc->nrows;
			f = M->CSR_gloc->ncols;
			PetscMalloc1((PetscInt)e, &oi);
			PetscMalloc1((PetscInt)f, &oj);
			PetscMalloc1((PetscInt)f, &oa);

			for(PetscInt q = 0; q < e; q++){
				oi[q] = M->CSR_gloc->rows[q];
			}

			for(PetscInt q = 0; q < f; q++){
				oj[q] = M->CSR_gloc->cols[q];
				oa[q] = M->CSR_gloc->vals[q].real() + PETSC_i * M->CSR_gloc->vals[q].imag();
			}
		}

		//MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE,  M->CSR_lloc->rows.data(), M->CSR_lloc->cols.data(), M->CSR_lloc->vals.data(), M->CSR_gloc->rows.data(), M->CSR_gloc->cols.data(), M->CSR_gloc->vals.data(), &A);

		if(M->CSR_gloc == NULL){
			MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

		}
		else if(M->CSR_lloc == NULL){
			MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE,  oi, oj, oa, &A);
		}
		else if (M->CSR_gloc != NULL && M->CSR_lloc != NULL){
			MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, oi, oj, oa, &A);
		}
		return A;

	}

#elif defined (PETSC_USE_REAL_DOUBLE) && defined (PETSC_USE_64BIT_INDICES)
//real long int double
#elif defined (PETSC_USE_REAL_SINGLE) && defined (PETSC_USE_64BIT_INDICES)
//real long int single

#elif defined (PETSC_USE_REAL_DOUBLE) 

	Mat ConvertToPETSCMat(parMatrixSparse<double, int > *M){

		int rank, size;

		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		int g, b, c, d, e, f;
		PetscInt m, n;
		PetscInt *i, *j, *oi, *oj;
		PetscScalar *a, *oa;

		Mat A;

		M->GetLocalSize(g, b);

		m = (PetscInt)g;
		n = (PetscInt)b;

		if(M->CSR_lloc != NULL){
			c = M->CSR_lloc->rows.size();
			d = M->CSR_lloc->cols.size();
			PetscMalloc1((PetscInt)c, &i);
			PetscMalloc1((PetscInt)d, &j);
			PetscMalloc1((PetscInt)d, &a);

			for(PetscInt q =0; q < c; q++){
				i[q] = M->CSR_lloc->rows[q];
			}

			for(PetscInt q = 0; q < d; q++){
				j[q] = M->CSR_lloc->cols[q];
				a[q] = M->CSR_lloc->vals[q].real();
			}
		}

		if(M->CSR_gloc != NULL){
			e = M->CSR_gloc->nrows;
			f = M->CSR_gloc->ncols;
			PetscMalloc1((PetscInt)e, &oi);
			PetscMalloc1((PetscInt)f, &oj);
			PetscMalloc1((PetscInt)f, &oa);

			for(PetscInt q = 0; q < e; q++){
				oi[q] = M->CSR_gloc->rows[q];
			}

			for(PetscInt q = 0; q < f; q++){
				oj[q] = M->CSR_gloc->cols[q];
				oa[q] = M->CSR_gloc->vals[q];
			}
		}

		//MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE,  M->CSR_lloc->rows.data(), M->CSR_lloc->cols.data(), M->CSR_lloc->vals.data(), M->CSR_gloc->rows.data(), M->CSR_gloc->cols.data(), M->CSR_gloc->vals.data(), &A);

		if(M->CSR_gloc == NULL){
			MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

		}
		else if(M->CSR_lloc == NULL){
			MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE,  oi, oj, oa, &A);
		}
		else if (M->CSR_gloc != NULL && M->CSR_lloc != NULL){
			MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, oi, oj, oa, &A);
		}
		return A;

	}

//real double int
#elif defined (PETSC_USE_REAL_SINGLE)
//real single int

#endif
