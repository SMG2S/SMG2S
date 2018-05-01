#!/bin/bash
#SBATCH --comment "Hello ROMEO!"
#SBATCH -J "TEST 1"

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=02:30:00

#SBATCH -n 2


#./ex12 -ksp_monitor_true_residual -pc_type none -m 200 -n 200 -log_summary
#./ex7 -ksp_monitor_true_residual -log_summary

#make runa
#srun -n 1 ./powerInverse.exe  -n 100 -l 10 -eps_monitor_conv -eps_power_shift_type constant \
#	-st_type sinvert -exact_value 14.4167+5.3509i -test_tol 0.00001 -degree 8


srun -n 2 ./test.exe


