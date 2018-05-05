/*
   This file is part of SMG2S.
   Author(s): Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr or xinzhe.wu1990@gmail.com>
        Date: 2018-04-20
   Copyright (C) 2018-     Xinzhe WU

   SMG2S is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   SMG2S is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with SMG2S.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "petscmat.h"
#include "petscerror.h"
#include <string.h>
#include "petsc.h"
#include "petscvec.h"
#include "petscksp.h"
#include "mpi.h"
#include <time.h>
#include <unistd.h>
#include <interface/C/c_wrapper.h>

static const char help[] = "Solve";


PetscReal uniform_distribution(PetscReal rangeLow, PetscReal rangeHigh) {
  PetscReal myRand = rand()/(1.0 + RAND_MAX);
  PetscReal range = rangeHigh - rangeLow;
  PetscReal myRand_scaled = (myRand * range) + rangeLow;
  return myRand_scaled;
}

PetscErrorCode generate_random_seed_vector(PetscInt size, PetscReal low, PetscReal up, PetscReal seed, Vec * v){
	PetscErrorCode ierr;
	PetscInt   i;
	PetscScalar  tmp;

	VecCreate(PETSC_COMM_WORLD, v);
	VecSetSizes(*v, PETSC_DECIDE, size);
	VecSetFromOptions(*v);

	srand(seed);

	for(i = 0; i < size; i++){
		tmp = uniform_distribution(low,up);
		VecSetValues(*v, 1, &i, &tmp, INSERT_VALUES);
	}

  VecAssemblyBegin(*v);
  VecAssemblyEnd(*v);
  return ierr;
}

PetscErrorCode mat_generate(Mat * A, MPI_Comm comm){

	int *rows, *cols;
  int size_row, size_col;
  double *real, *imag;

	char spectra[PETSC_MAX_PATH_LEN];
	PetscBool flgsmg2s, flglb, flgsize, flgno, flgsptr;
	PetscErrorCode ierr;

	PetscInt lbandwidth, size, nbOne;
	ierr= PetscOptionsHasName(NULL,NULL,"-smg2s",&flgsmg2s);CHKERRQ(ierr);
	ierr= PetscOptionsGetInt(NULL,NULL,"-lbandwidth",&lbandwidth,&flglb);CHKERRQ(ierr);
	ierr= PetscOptionsGetInt(NULL,NULL,"-size",&size,&flgsize);CHKERRQ(ierr);
  	ierr= PetscOptionsGetInt(NULL,NULL,"-continous1",&nbOne,&flgno);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-sptr",spectra,PETSC_MAX_PATH_LEN-1,&flgsptr);CHKERRQ(ierr);
	if(flgsmg2s) PetscPrintf(comm, "Info ]> Using SMG2S to generate test matrices\n");
	if(!flglb){
		lbandwidth = 3;
		PetscPrintf(comm, "Remainder ]> low bandwidth in SMG2S is not set, use the default value = 3\n");
	}

        if(!flgsize){
                size = 100;
                PetscPrintf(comm, "Remainder ]> matrix size in SMG2S is not set, use the default value = 100\n");
        }

        if(!flgno){
                nbOne = 2;
                PetscPrintf(comm, "Remainder ]> Continus 1 in SMG2S is not set, use the default value = 2\n");
        }

	if(!flgsptr){
                strncpy(spectra, " ", sizeof(spectra));
		PetscPrintf(comm, "Remainder ]> given spectra file in SMG2S is not set, use the default method to generate\n");
        }

	struct NilpotencyInt *n;
        n = newNilpotencyInt();
        NilpType1(n, nbOne, size);
//        showNilpotencyInt(n);

        struct parMatrixSparseComplexDoubleInt *m;
        m = newParMatrixSparseComplexDoubleInt();

	double t1 = MPI_Wtime();
	smg2sComplexDoubleInt(m, size, n, lbandwidth ,spectra, comm);
        double t2 = MPI_Wtime();
	PetscPrintf(comm, "SMG2S time = %f\n", t2-t1);
        Loc_ConvertToCSRComplexDoubleInt(m);
        Loc_CSRGetRowsArraySizes(m, &size_row, &size_col);

        rows = (int *) malloc(size_row*sizeof(int));
        cols = (int *) malloc(size_col*sizeof(int));

        real = (double *) malloc(size_col*sizeof(double));
        imag = (double *) malloc(size_col*sizeof(double));

        Loc_CSRGetRowsArrays(m, size_row, &rows, size_col, &cols, &real, &imag);

        PetscInt *i, *j, q, p;
        PetscScalar *a;

        PetscMalloc1((PetscInt)size_row, &i);
        PetscMalloc1((PetscInt)size_col, &j);
        PetscMalloc1((PetscInt)size_col, &a);


	for(q =0; q < size_row; q++){
                i[q] = rows[q];
        }

        for(p = 0; p < size_col; p++){
                j[p] = cols[p];
                a[p] = (PetscReal) real[p] + PETSC_i * (PetscReal) imag[p];
	}


        MatCreate(PETSC_COMM_WORLD, A);
        MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, size_row-1, size_row-1, PETSC_DETERMINE, PETSC_DETERMINE,  i, j, a, A);


//	free(rows); free(cols); free(real); free(imag);
//	PetscFree(i); PetscFree(j); PetscFree(a);

//	ReleaseParMatrixSparseComplexDoubleInt(&m);
	ReleaseNilpotencyInt(&n);

	return ierr;

}

PetscErrorCode read_matrix_vector(Mat * A, Vec * v){
	char filea[PETSC_MAX_PATH_LEN];
	char fileb[PETSC_MAX_PATH_LEN];
	char err[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flaga,flagb, flgsmg2s;
	PetscViewer fd;
	PetscInt size,sizea;
	PetscScalar scal;
	PetscReal vnorm;

	ierr= PetscOptionsHasName(NULL,NULL,"-smg2s",&flgsmg2s);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-mfile",filea,PETSC_MAX_PATH_LEN-1,&flaga);CHKERRQ(ierr);
	if (!flaga) {
		if(!flgsmg2s){
		        sprintf(err,"Error : mfile is not properly set -> %s\n",filea);
		}
		else{
			PetscPrintf(PETSC_COMM_WORLD,"Remainder ]> USING SMG2S to generate test matrix\n");
			mat_generate(A, MPI_COMM_WORLD);
		}
	}
	else{
		/* read matrix file */
	        PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",filea);
       		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,filea,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        	ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
        	ierr=MatLoad(*A,fd);CHKERRQ(ierr);
        	ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);

	}

	/* read matrix file */
	ierr=MatGetSize(*A,&size,&sizea);CHKERRQ(ierr);

	if (!flagb) {
		/* the user did not provide a vector, so generate it*/
		generate_random_seed_vector(size, -10.0, 10.0, 0, v);
		VecNorm(*v, NORM_2, &vnorm);
		VecScale(*v, 0.01/vnorm);
		PetscPrintf(PETSC_COMM_WORLD,"Generated right hand side matrix b\n");
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Loading Vector : %s\n",fileb);
		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileb,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		VecCreate(PETSC_COMM_WORLD,v);
		ierr=VecLoad(*v,fd);CHKERRQ(ierr);
		ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Loaded Vector of size : %d\n",size);
	}

	return 0;
}

PetscErrorCode classicalGMRES(Vec * b, Mat * A){
	PetscErrorCode ierr;
	KSP ksp, kspp;
	Vec x, c;
	KSPConvergedReason reason;
	PetscInt its, nols, ntimes;
	int i, size;
	PetscBool flagls, flagtimes;
	double cost_time;
	clock_t start, end;

	PetscReal norm;
	VecGetSize(*b, &size);
	VecDuplicate(*b, &c);
        PetscPrintf(PETSC_COMM_WORLD,"#} GMRES Creating and setting vector x\n");
	ierr = VecDuplicate(*b, &x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
	ierr = VecSet(x, 0.0); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,PETSC_NULL,"-ksp_ls_nopc",&nols,&flagls);CHKERRQ(ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
//	ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);
//	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	ierr = KSPSetOperators(ksp, *A, *A);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

	PetscOptionsGetInt(NULL,NULL,"-ntimes",&ntimes,&flagtimes);

	if(!flagtimes){
		ntimes=1;
	}

	for (i=1; i<=ntimes; i++) {
		if(i==1){
			VecCopy(*b,c);
			VecNorm(c, NORM_2,&norm);
			VecScale(c, 0.01/norm);
		}
		else{
			generate_random_seed_vector(size, -10,10, i,&c);
			VecNorm(c, NORM_2,&norm);
			VecScale(c, 0.01/norm);
		}
		start=clock();
		ierr = KSPSolve(ksp, c, x); CHKERRQ(ierr);
		end=clock();
  		cost_time = (double)(end - start)/CLOCKS_PER_SEC;
		KSPGetConvergedReason(ksp,&reason);
		if (reason<0) {
			PetscPrintf(PETSC_COMM_WORLD,"\nResolution %d: Divergence in acceptable iteration steps.\n", i);
		}
		else {
			KSPGetIterationNumber(ksp,&its);
			PetscPrintf(PETSC_COMM_WORLD,"\nResolution %d: Convergence in %d iterations. \n", i, its);
		}
	}
//	int exit_type=666;
//	mpi_lsa_com_type_send(com,&exit_type);
	PetscPrintf(PETSC_COMM_WORLD,"\n\n#}Finish the resolution\n");
  	return ierr;
}


int main(int argc, char ** argv){
	Mat A;
	Vec v;
	PetscErrorCode ierr;
	PetscBool flag, gft_flg;
	int	size, rank;
	PetscLogDouble st, ed;

	/* init of MPI and MPI Utils */
	MPI_Init(&argc,&argv);

	ierr=PetscInitialize(&argc,&argv,(char *)0,help);       CHKERRQ(ierr);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  PetscPrintf(PETSC_COMM_WORLD,"MPI WORLD SIZE = %d \n", size);
  PetscTime(&st);
	read_matrix_vector(&A, &v);
  PetscTime(&ed);
  PetscPrintf(PETSC_COMM_WORLD, "The test matrix generation time is %f \n", ed-st);
  PetscTime(&st);
	classicalGMRES(&v, &A);
  PetscTime(&ed);
  PetscPrintf(PETSC_COMM_WORLD, "The Krylov method resolving time is %f \n", ed-st);
  MPI_Barrier(MPI_COMM_WORLD);
  VecDestroy(&v);
  MatDestroy(&A);
  MPI_Barrier(MPI_COMM_WORLD); //wait for all
  PetscFinalize(); //finalize petsc



	MPI_Finalize(); //finalize

	return 0;
}
