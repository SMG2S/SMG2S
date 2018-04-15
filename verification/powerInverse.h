#include "petsc.h"
#include <slepceps.h>
#include "../interface/PETSc/petsc_interface.h"

PetscErrorCode loadInputs(Mat * A, Vec * x);
PetscErrorCode loadMatrix(Mat * A);
PetscErrorCode loadVector(char * type_v,Vec * b);
PetscErrorCode generateVector(int size, Vec * v);
PetscErrorCode generateVectorRandom(int size, Vec * v);
PetscErrorCode generateVectorNorm(int size, Vec * v);
PetscInt factorial(PetscInt low, PetscInt up);


