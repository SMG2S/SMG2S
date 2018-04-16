#include "petscmat.h"
#include "../../src/matgen/smg2s.cc"
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

		int g, b, c, d;
		PetscInt m, n;
		PetscInt *i, *j, ii;
		PetscScalar *a;

		Mat A;

		M->GetLocalSize(g, b);

		m = (PetscInt)g;
		n = (PetscInt)m;

		//printf("Proc: %d, m = %d, n = %d\n", rank, m, n);

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

//			printf("Proc: %d, c = %d, d = %d\n",rank, c, d);


		PetscMalloc1((PetscInt)count, &i);
		PetscMalloc1((PetscInt)count2, &j);
		PetscMalloc1((PetscInt)count2, &a);

		for(PetscInt q =0; q < count; q++){
			i[q] = M->CSR_loc->rows[q];
		}

		for(PetscInt q = 0; q < count2; q++){
			j[q] = M->CSR_loc->cols[q];
			a[q] = M->CSR_loc->vals[q];
//			printf("Proc: %d, j[%d] = %d, a[%d] = %f\n", rank, q, j[q],q, M->CSR_loc->vals[q]);
			}

		MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, count-1, count-1, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

		return A;

/*
		int rank, size;

		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		int g, b;
		PetscInt m, n;
		PetscInt *j;
		PetscScalar *a;
		PetscInt rindx;

		Mat A;

		typename std::map<int,std::complex<double> >::iterator it;
		parVectorMap<int>	*y_index_map = M->GetYMap();


		M->GetLocalSize(g, b);

		m = (PetscInt)g;
		n = (PetscInt)b;

		MatCreate(PETSC_COMM_WORLD,&A);
		MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
		MatSetType(A,MATMPIAIJ);
		MatSetFromOptions(A);
		MatSetUp(A);


		PetscInt count, count2;


		for(PetscInt ii = 0; ii < g; ii++){
			count = 0;
			count2 = 0;
			rindx = (PetscInt)y_index_map->Loc2Glob(ii);

			for(it = M->dynmat_loc[ii].begin(); it != M->dynmat_loc[ii].end(); ++it){
				if(it->second.real() != 0 || it->second.imag() != 0){
					count2++;
				}
		    }

			PetscMalloc1(count2, &j);
		    PetscMalloc1(count2, &a);
		    for(it = M->dynmat_loc[ii].begin(); it != M->dynmat_loc[ii].end(); ++it){
		    	if(it->second.real() != 0 || it->second.imag() != 0){
		    		j[count] = (PetscInt)it->first;
			    	a[count] = (PetscReal)it->second.real()+PETSC_i*(PetscReal)it->second.imag();
			    	count++;
		    	}
		    }
			//std::cout << rank << "   " << rindx << "    " << s << std::endl;
			MatSetValues(A,1, &rindx, count2, j, a, INSERT_VALUES);

		}

		MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

		return A;
		*/

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

		//printf("Proc: %d, m = %d, n = %d\n", rank, m, n);

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

//			printf("Proc: %d, c = %d, d = %d\n",rank, c, d);


		PetscMalloc1((PetscInt)count, &i);
		PetscMalloc1((PetscInt)count2, &j);
		PetscMalloc1((PetscInt)count2, &a);

		for(PetscInt q =0; q < count; q++){
			i[q] = M->CSR_loc->rows[q];
		}

		for(PetscInt q = 0; q < count2; q++){
			j[q] = M->CSR_loc->cols[q];
			a[q] = M->CSR_loc->vals[q];
//			printf("Proc: %d, j[%d] = %d, a[%d] = %f\n", rank, q, j[q],q, M->CSR_loc->vals[q]);
			}

		MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, count-1, count-1, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, &A);

		return A;

	}

//real double int
#elif defined (PETSC_USE_REAL_SINGLE)
//real single int

#endif
