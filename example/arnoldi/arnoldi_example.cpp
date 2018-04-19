/*
   This file is part of SMG2S.
   Author(s): Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr or xinzhe.wu1990@gmail.com>
        Date: 2018-04-20
   Copyright (C) 2018-     Xinzhe WU
   
   SMG2S is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   HPDDM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with HPDDM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "arnoldi_example.h"


static const char help[] = "Solve Non-Hermitian eigenvalues problem by Arnoldi, options array_in_received_buffer\n\
\t-mfile matrix_file (matrix file in PETSc bin format, this is mandatory)\n\
\t-xfile initial_guess_file (in PETSc bin format)\n";



int main(int argc, char **argv){

	MPI_Init(&argc, &argv);
	
	PetscErrorCode ierr;
	Mat A;
	PetscInt degree, n, lbandwidth;
	PetscBool degree_flg, n_flg, l_flg;	
	PetscInt its, nev, nconv;
	EPS eps;
	EPSType type;

	ierr=SlepcInitialize(&argc,&argv,PETSC_NULL,help);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Initializing PETSc/SLEPc\n");
	
	int world_size;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    //MPI_Barrier(MPI_COMM_WORLD);
    
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);


    if(world_rank == 0){printf("Matrix Initialized\n");}

    MPI_Barrier(MPI_COMM_WORLD);

	ierr = PetscOptionsGetInt(NULL, PETSC_NULL, "-n", &n, &n_flg);
    ierr = PetscOptionsGetInt(NULL, PETSC_NULL, "-degree", &degree, &degree_flg);
    ierr = PetscOptionsGetInt(NULL, PETSC_NULL, "-l", &lbandwidth, &l_flg);

    if(!n_flg){
	  	PetscPrintf(PETSC_COMM_WORLD, "ERROR: Set Matrix dimension to generate... \n");
        PetscPrintf(PETSC_COMM_WORLD, "ERROR: Exit with errors ... \n\n");
	  return 0;
    } 

	if(!degree_flg){
	  	PetscPrintf(PETSC_COMM_WORLD, "ERROR: Please set the parameter degree to test... \n");
      	PetscPrintf(PETSC_COMM_WORLD, "ERROR: Exit with errors ... \n\n");
	  return 0;
    } 

	if(!l_flg){
	  	PetscPrintf(PETSC_COMM_WORLD, "ERROR: Please set the parameter lbandwidth to test... \n");
      	PetscPrintf(PETSC_COMM_WORLD, "ERROR: Exit with errors ... \n\n");	
    }

	Nilpotency<int> nilp;
 	nilp.NilpType1(degree,n);

#if defined (PETSC_USE_COMPLEX)

 	parMatrixSparse<std::complex<double>,int> *Mt;

 	Mt =  smg2s<std::complex<double>,int> (n, nilp,lbandwidth);

#else

 	parMatrixSparse<double,int> *Mt;

 	Mt =  smg2s<double,int> (n, nilp,lbandwidth);

#endif
  
    Mt->LOC_MatView();

    Mt->Loc_ConvertToCSR();


    MatCreate(PETSC_COMM_WORLD,&A);

    A = ConvertToPETSCMat(Mt);

    MatView(A, PETSC_VIEWER_STDOUT_WORLD);

	PetscPrintf(PETSC_COMM_WORLD,"]> Data Generated\n");

	/*Create the KSP context and setup*/

	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	ierr = EPSSetOperators(eps, A, NULL);CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps, EPS_NHEP);CHKERRQ(ierr);
	ierr = EPSSetType(eps,EPSARNOLDI);CHKERRQ(ierr);
	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver settings done\n");

	/*Solve the system*/
	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver Launching solving process\n");
	ierr = EPSSolve(eps);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver System solved\n");

/*Get some informations of resolution*/

	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);
	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);
	ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);
	/*Display the solution*/
	EPSGetConverged(eps,&nconv);
	PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);
	/*Clean*/
	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Cleaned structures, finalizing\n");

	/*Finalize PETSc*/
	SlepcFinalize(); 

	delete Mt;

	MPI_Finalize();

	return 0;
}
