#!/bin/bash
#SBATCH -n 16

export OMP_NUM_THREADS=2
 	 
srun -n 16 ./a.out

