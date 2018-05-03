#!/bin/bash
#SBATCH --comment "Hello ROMEO!"
#SBATCH -J "TEST 1"

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=02:30:00

#SBATCH -N 1
#SBATCH -n 2


#make runexec
make runsmg2s

