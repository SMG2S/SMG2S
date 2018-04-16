#include "gmres_example.h"


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


PetscErrorCode loadInputs(Mat * A, Vec * b, Vec * x){
	PetscErrorCode ierr;
	PetscInt sizex,sizey;
	char bfile[]="-bfile";
	char xfile[]="-xfile";
	
	//load data files
	ierr=loadMatrix(A);CHKERRQ(ierr);
	ierr=loadVector(bfile,b);CHKERRQ(ierr);
	if(*b==NULL) {
		PetscPrintf(PETSC_COMM_WORLD,"]> Creating vector b\n");
		ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
		ierr=generateVectorRandom(sizex,b);CHKERRQ(ierr);
	}
	ierr=loadVector(xfile,x);CHKERRQ(ierr);
	if(*x==NULL) {
		PetscPrintf(PETSC_COMM_WORLD,"]> Creating vector x\n");
		ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
		ierr=generateVectorRandom(sizex,x);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode loadMatrix(Mat * A){
	char file[PETSC_MAX_PATH_LEN];
	char err[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flag;
	PetscViewer fd;
	PetscInt sizex,sizey;

	/*check args, if no matrix then no work... matrix file is mandatory*/
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-mfile",file,PETSC_MAX_PATH_LEN-1,&flag);CHKERRQ(ierr);
	if (!flag) {		
		sprintf(err,"Error : mfile is not properly set -> %s\n",file);
		SETERRQ(PETSC_COMM_WORLD,(PetscErrorCode)83,err);
	}

	/* read matrix file */
	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",file);

	ierr=MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
	ierr=MatSetType(*A,MATAIJ);CHKERRQ(ierr);

	ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr=MatLoad(*A,fd);CHKERRQ(ierr);
	ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Loaded Matrix of size : %d %d\n",sizex,sizey);

	return 0;
}


PetscErrorCode loadVector(char * type_v,Vec * b){
	char file[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flag;
	PetscViewer fd;
	PetscInt size;

	// check if there is a vec file, vector is not mandatory
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,type_v,file,PETSC_MAX_PATH_LEN-1,&flag);CHKERRQ(ierr);
	if (!flag) {		
		PetscPrintf(PETSC_COMM_WORLD,"Error : %s is not properly set\n",type_v);
		*b = NULL;
	}else{
		PetscPrintf(PETSC_COMM_WORLD,"Loading Vector : %s\n",file);
		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr=VecLoad(*b,fd);CHKERRQ(ierr);
		ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
		ierr=VecGetSize(*b,&size);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Loaded Vector of size : %d\n",size);
	}

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


