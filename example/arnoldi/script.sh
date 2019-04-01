#!/bin/bash
#SBATCH --comment "Hello ROMEO!"
#SBATCH -J "TEST 1"

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=02:30:00

#SBATCH -n 2




make runa
#srun -n 1 ./powerInverse.exe  -n 1000 -l 10 -eps_monitor_conv -eps_power_shift_type constant \
#	-st_type sinvert -exact_value 230+230i -test_tol 0.00001 -degree 4


