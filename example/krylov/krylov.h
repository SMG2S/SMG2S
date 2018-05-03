#ifndef _KRYLOV_H
#define _KRYLOV_H

#include "petscmat.h"
#include "petscerror.h"
#include <string.h>
#include "petsc.h"
#include "petscvec.h"
#include "petscksp.h"
#include "mpi.h"
#include <time.h>
#include <unistd.h>

#endif

PetscErrorCode mat_generate(Mat * A, MPI_Comm comm);
PetscErrorCode generate_random_seed_vector(PetscInt size, PetscReal low, PetscReal up, PetscReal seed, Vec * v);
PetscErrorCode read_matrix_vector(Mat * A, Vec * v);
PetscErrorCode classicalGMRES(Vec * b, Mat * A);
