#!/bin/bash
#SBATCH -n 1

export OMP_NUM_THREADS=4
 	 
./a.out

