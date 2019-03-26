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

 	Mt =  smg2s<std::complex<double>,int> (n, nilp,lbandwidth,spec,MPI_COMM_WORLD);

#else

 	parMatrixSparse<double,int> *Mt;

 	Mt =  smg2s<double,int> (n, nilp,lbandwidth,spec,MPI_COMM_WORLD);

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
	//ierr = KSPSetType(ksp, KSPTFQMR);CHKERRQ(ierr);
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
