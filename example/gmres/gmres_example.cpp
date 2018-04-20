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


#include "gmres_example.h"
#include "string"

static const char help[] = "Solve Ax=b using GMRES, options array_in_received_buffer\n\
\t-mfile matrix_file (matrix file in PETSc bin format, this is mandatory)\n\
\t-bfile right_and_side_file (also in PETSc bin format)\n\
\t-xfile initial_solution_file (in PETSc bin format)\n";



int main(int argc, char **argv){

	MPI_Init(&argc, &argv);
	
	PetscErrorCode ierr;
	Vec x,b;
	Mat A;
	KSP ksp;
	PetscInt degree, n, lbandwidth;
	PetscBool degree_flg, n_flg, l_flg;	
    
    std::string spec;
    PetscBool sptr_flg;

	char	spectrum[PETSC_MAX_PATH_LEN]=" ";

	ierr=PetscInitialize(&argc,&argv,PETSC_NULL,help);CHKERRQ(ierr);
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
    ierr = PetscOptionsGetString(NULL,NULL,"-sptr",spectrum,sizeof(spectrum),&sptr_flg);


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

    spec.assign(spectrum);

	Nilpotency<int> nilp;
 	nilp.NilpType1(degree,n);

#if defined (PETSC_USE_COMPLEX)

 	parMatrixSparse<std::complex<double>,int> *Mt;

 	Mt =  smg2s<std::complex<double>,int> (n, nilp,lbandwidth,spec);

#else

 	parMatrixSparse<double,int> *Mt;

 	Mt =  smg2s<double,int> (n, nilp,lbandwidth,spec);

#endif
  
    Mt->LOC_MatView();

    Mt->Loc_ConvertToCSR();


    MatCreate(PETSC_COMM_WORLD,&A);

    A = ConvertToPETSCMat(Mt);

    MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    PetscInt sizex, sizey;

    ierr=MatGetSize(A,&sizex,&sizey);CHKERRQ(ierr);

	ierr=generateVectorRandom(sizex,&b);CHKERRQ(ierr);
	ierr=generateVectorRandom(sizex,&x);CHKERRQ(ierr);


	PetscPrintf(PETSC_COMM_WORLD,"]> Data Generated\n");

	/*Create the KSP context and setup*/

	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);	
	ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);	
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);	
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);	
	ierr = KSPSetUp(ksp);CHKERRQ(ierr);	
	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver settings done\n");

	/*Solve the system*/
	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver Launching solving process\n");
	ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver System solved\n");

	/*Clean*/
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Cleaned structures, finalizing\n");

	/*Finalize PETSc*/
	PetscFinalize(); 

	delete Mt;

	MPI_Finalize();

	return 0;
}

PetscErrorCode generateVectorRandom(int size, Vec * v){
	PetscErrorCode ierr;

	ierr=PetscPrintf(PETSC_COMM_WORLD,"Generating Vector \n");CHKERRQ(ierr);
	ierr=generateVector(size,v);CHKERRQ(ierr);
	ierr=VecSetRandom(*v,PETSC_NULL);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Generated Random Vector of size : %d\n",size);	

	return 0;
}


PetscErrorCode generateVectorNorm(int size, Vec * v){
	PetscScalar scal;
	PetscErrorCode ierr;

	ierr=PetscPrintf(PETSC_COMM_WORLD,"Generating Vector \n");CHKERRQ(ierr);
	ierr=generateVector(size,v);CHKERRQ(ierr);
	scal=1.0/size;
	ierr=VecSet(*v,scal);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Generated Norm Vector of size : %d\n",size);	

	return 0;
}


PetscErrorCode generateVector(int size, Vec * v){
	PetscErrorCode ierr;

	ierr=VecCreate(PETSC_COMM_WORLD,v);CHKERRQ(ierr);
	ierr=VecSetSizes(*v,PETSC_DECIDE,size);CHKERRQ(ierr);
	ierr=VecSetFromOptions(*v);CHKERRQ(ierr);
	/*initiate the vector to its norm*/

	return 0;
}


