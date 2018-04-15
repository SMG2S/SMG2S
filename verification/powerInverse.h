#include "petsc.h"
#include <slepceps.h>
#include "../interface/PETSc/petsc_interface.h"

//#include "../src/matgen/smg2s.cc"
#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include "../config/config.h"
#include "../utils/utils.h"


PetscErrorCode loadInputs(Mat * A, Vec * x);
PetscErrorCode loadMatrix(Mat * A);
PetscErrorCode loadVector(char * type_v,Vec * b);
PetscErrorCode generateVector(int size, Vec * v);
PetscErrorCode generateVectorRandom(int size, Vec * v);
PetscErrorCode generateVectorNorm(int size, Vec * v);
PetscInt factorial(PetscInt low, PetscInt up);


