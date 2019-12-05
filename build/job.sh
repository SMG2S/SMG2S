#!/bin/bash -x
#SBATCH --account=slai
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi-err.%j
#SBATCH --time=00:15:00
#SBATCH --partition=devel

#export SCOREP_EXPERIMENT_DIRECTORY=scorep_sum
#export SCOREP_ENABLE_TRACING=true
#export SCOREP_ENABLE_PROFILING=false
#export SCOREP_TOTAL_MEMORY=140M

PROCS=4
MAT_SIZE=1000
LOW_BANDWIDTH=10
CONTINUOUS_ONES=4
FLOATTYPE=CPLX_DOUBLE
INTEGERTYPE=INT

srun -n ${PROCS} ./smg2s.exe -SIZE ${MAT_SIZE} -L ${LOW_BANDWIDTH} -C ${CONTINUOUS_ONES} -floattype ${FLOATTYPE} -integertype ${INTEGERTYPE}

