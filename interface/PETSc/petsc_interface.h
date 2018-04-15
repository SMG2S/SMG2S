#include "petscmat.h"
#include "../../src/matgen/smg2s.cc"
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

		int g, b, c, d;
		PetscInt m, n;
		PetscInt *i, *j;
		PetscScalar *a;

		Mat A;

		MatCreate(PETSC_COMM_WORLD,&A);

		M->GetLocalSize(g, b);

		m = (PetscInt)g;
		n = (PetscInt)m;

//		printf("Proc: %d, m = %d, n = %d\n", rank, m, n);

		if(M->CSR_loc != NULL){
			c = M->CSR_loc->rows.size();
			d = M->CSR_loc->cols.size();

//			printf("Proc: %d, c = %d, d = %d\n",rank, c, d);


			PetscMalloc1((PetscInt)c, &i);
			PetscMalloc1((PetscInt)d, &j);
			PetscMalloc1((PetscInt)d, &a);

			for(PetscInt q =0; q < c; q++){
				i[q] = M->CSR_loc->rows[q];
//				printf("Proc: %d, i[%d] = %d\n", rank, q, i[q]);
			}

			for(PetscInt q = 0; q < d; q++){
				j[q] = M->CSR_loc->cols[q];
				a[q] = M->CSR_loc->vals[q].real() + PETSC_i * M->CSR_loc->vals[q].imag();
//				printf("Proc: %d, j[%d] = %d, a[%d] = %f+%fi\n", rank, q, j[q],q, M->CSR_loc->vals[q].real(), M->CSR_loc->vals[q].imag());
			}

			MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, m, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

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

		int g, b, c, d;
		PetscInt m, n;
		PetscInt *i, *j;
		PetscScalar *a;

		Mat A;

		M->GetLocalSize(g, b);

		m = (PetscInt)g;
		n = (PetscInt)m;

		//printf("Proc: %d, m = %d, n = %d\n", rank, m, n);

		if(M->CSR_loc != NULL){
			c = M->CSR_loc->rows.size();
			d = M->CSR_loc->cols.size();

//			printf("Proc: %d, c = %d, d = %d\n",rank, c, d);


			PetscMalloc1((PetscInt)c, &i);
			PetscMalloc1((PetscInt)d, &j);
			PetscMalloc1((PetscInt)d, &a);

			for(PetscInt q =0; q < c; q++){
				i[q] = M->CSR_loc->rows[q];
//				printf("Proc: %d, i[%d] = %d\n", rank, q, i[q]);
			}

			for(PetscInt q = 0; q < d; q++){
				j[q] = M->CSR_loc->cols[q];
				a[q] = M->CSR_loc->vals[q];
//				printf("Proc: %d, j[%d] = %d, a[%d] = %f\n", rank, q, j[q],q, M->CSR_loc->vals[q]);
			}

			MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, m, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

		}

		return A;

	}

//real double int
#elif defined (PETSC_USE_REAL_SINGLE)
//real single int

#endif
