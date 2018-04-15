#include "powerInverse.h"


static const char help[] = "Solve Non-Hermitian eigenvalues problem by Power Iterative method, options array_in_received_buffer\n\
\t-mfile matrix_file (matrix file in PETSc bin format, this is mandatory)\n\
\t-xfile initial_guess_file (in PETSc bin format)\n";

int main(int argc, char **argv){
	PetscErrorCode ierr;
	Vec x;
	Mat A;
	EPS eps;
	PetscInt its, nev, nconv;
	EPSType type;
	PetscScalar vtest, target, eigr, eigi;
	PetscBool   vtest_flg, tol_flg, degree_flg;	
	Vec Ax;
	Vec xr, xi, eigenvector;
	PetscReal norm, vnorm, residual;
	PetscReal test_tol;
	PetscInt  degree, fac;

	ierr=SlepcInitialize(&argc,&argv,PETSC_NULL,help);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\n]> Initializing SLEPc\n");

	ierr = PetscOptionsGetScalar(NULL, PETSC_NULL, "-exact_value", &vtest, &vtest_flg);
        ierr = PetscOptionsGetScalar(NULL, PETSC_NULL, "-test_tol", &test_tol, &tol_flg);
        ierr = PetscOptionsGetInt(NULL, PETSC_NULL, "-degree", &degree, &degree_flg);
/*        
        if(!degree_flg){
          PetscPrintf(PETSC_COMM_WORLD, "ERROR: Please set the degree used by matGen... \n");
          PetscPrintf(PETSC_COMM_WORLD, "ERROR: Exit with errors ... \n\n");
          return 0;
       } else{
          PetscPrintf(PETSC_COMM_WORLD, " @>Remainder: degree = %d ...\n", degree);
        }

        fac = factorial(1, 2*degree-2);
        PetscPrintf(PETSC_COMM_WORLD, " @>Remainder: fac =  %d ...\n", fac);
*/
	if(!vtest_flg){
	  PetscPrintf(PETSC_COMM_WORLD, "ERROR: Please set the exact expacted eigenvalues to test... \n");
          PetscPrintf(PETSC_COMM_WORLD, "ERROR: Exit with errors ... \n\n");
	  return 0;
       } 
//	else{
//	  vtest = vtest * fac;
//      	  vtest = vtest;
//          PetscPrintf(PETSC_COMM_WORLD, " @>Remainder: test value = %f + %fi...\n", PetscRealPart(vtest), PetscImaginaryPart(vtest));
//	}
	if(!tol_flg){
		test_tol = 1e-8;
          	PetscPrintf(PETSC_COMM_WORLD, " @>Remainder: Not set the tolerance for validation, use the defaut value tol =  %.2e ...\n",test_tol);
	}
	/*Load data*/
	ierr=loadInputs(&A,&x);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Data loaded\n");

	/*Create the EPS context and setup*/
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	ierr = EPSSetOperators(eps, A, NULL);CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps, EPS_NHEP);CHKERRQ(ierr);
	ierr = EPSSetType(eps,EPSPOWER);CHKERRQ(ierr);
	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
	EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
	#ifdef PETSC_USE_COMPLEX
		target = PetscRealPart(vtest)+0.001+PetscImaginaryPart(vtest)*PETSC_i;
	#else
                target = vtest+0.001;
	#endif
	EPSSetTarget(eps,target);
	#ifdef PETSC_USE_COMPLEX
        	PetscPrintf(PETSC_COMM_WORLD, " @>Remainder: The input target value is %f + %fi ...\n", PetscRealPart(target), PetscImaginaryPart(target));
	#else
		PetscPrintf(PETSC_COMM_WORLD, " @>Remainder: The input target value is %f ...\n", target);
	#endif
	ierr=EPSSetInitialSpace(eps,1,&x);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver settings done\n");

	/*Solve the problem*/
	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver Launching solving process\n");
        ierr = EPSSolve(eps);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Krylov Solver System solved\n");

	/*Get some informations of resolution*/

	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," @> Number of iterations of the method: %D\n",its);
	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," @> Solution method: %s\n",type);
	ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD," @> Number of requested eigenvalues: %D\n",nev);
	/*Display the solution*/
	EPSGetConverged(eps,&nconv);
	if(nconv > 0){
		PetscPrintf(PETSC_COMM_WORLD," @> Number of converged eigenpairs: %D\n",nconv);
	}else {
                PetscPrintf(PETSC_COMM_WORLD," @> No converged eigenpairs, exit without validatint input value\n");
		return 0;
	}
        /*Validation*/
	PetscPrintf(PETSC_COMM_WORLD,"]> Validation the input test value\n");
	ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
	ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);
        ierr = MatCreateVecs(A,NULL,&eigenvector);CHKERRQ(ierr);
	ierr = MatCreateVecs(A,NULL,&Ax);CHKERRQ(ierr);
	ierr = EPSGetEigenvector(eps,0,xr,xi);CHKERRQ(ierr);
	ierr = EPSGetEigenvalue(eps,0,&eigr,&eigi);CHKERRQ(ierr);
//	VecView(xr,PETSC_VIEWER_STDOUT_SELF);
//	VecView(xi,PETSC_VIEWER_STDOUT_SELF);
//	PetscScalarView(1,&eigr,PETSC_VIEWER_STDOUT_SELF);
	PetscReal re, im;
	#ifdef PETSC_USE_COMPLEX
		re = PetscRealPart(eigr);
		im = PetscImaginaryPart(eigr);
		PetscPrintf(PETSC_COMM_WORLD,"@> The eigenvalue is %f + %fi \n", re, im);
	#else
		re = eigr;
		im = eigi;
		if(im == 0)
			PetscPrintf(PETSC_COMM_WORLD,"@> The eigenvalue is %f \n", eigr);
	        else
			PetscPrintf(PETSC_COMM_WORLD,"@> The eigenvalue is %f + %fi \n", re, im);
	#endif
	ierr = VecWAXPY(eigenvector,1,xr,xi);CHKERRQ(ierr);
	ierr = MatMult(A, eigenvector, Ax);CHKERRQ(ierr);
	ierr = VecScale(eigenvector, vtest);CHKERRQ(ierr);
	ierr = VecNorm(Ax, NORM_2, &vnorm);
	ierr = VecAYPX(Ax, -1, eigenvector);CHKERRQ(ierr);
	ierr = VecNorm(Ax, NORM_2, &norm);CHKERRQ(ierr);
//	ierr = VecNorm(eigenvector, NORM_2, &vnorm);
	residual = norm / vnorm;
//        residual = norm;
	PetscPrintf(PETSC_COMM_WORLD," \n     Residual:||Av-kv||/||Av||  \n");
        PetscPrintf(PETSC_COMM_WORLD,"     ------------------------- \n");
        PetscPrintf(PETSC_COMM_WORLD," \n           %.5e  \n\n    ", residual);
	if(residual <= test_tol){
        PetscPrintf(PETSC_COMM_WORLD," @> Residual:= %.5e < test_tol %.2e, validated!!!\n\n", residual, test_tol);
	}else{
        PetscPrintf(PETSC_COMM_WORLD," @> Residual:= %.5e > test_tol %.2e, not validated!!!\n\n", residual, test_tol);
	}
	/*Clean*/
	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = VecDestroy(&Ax);CHKERRQ(ierr);
        ierr = VecDestroy(&xi);CHKERRQ(ierr);
        ierr = VecDestroy(&xr);CHKERRQ(ierr);
        ierr = VecDestroy(&eigenvector);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"]> Cleaned structures, finalized ... \n\n");

	/*Finalize SLEPc*/
	SlepcFinalize(); 

	return 0;
}


PetscErrorCode loadInputs(Mat * A, Vec * x){
	PetscErrorCode ierr;
	PetscInt sizex,sizey;
	char xfile[]="-xfile";
	
	//load data files
	ierr=loadMatrix(A);CHKERRQ(ierr);
	ierr=loadVector(xfile,x);CHKERRQ(ierr);
	if(*x==NULL) {
		PetscPrintf(PETSC_COMM_WORLD,"]> Creating initial guessed vector x\n");
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

PetscInt factorial(PetscInt low, PetscInt up){

	PetscInt       fac = 1, size, k;

	size = up -low + 1;
	
	for(k = 0; k < size; k++){
		fac = fac * (low + k);
	}
	
	return fac;

}	
