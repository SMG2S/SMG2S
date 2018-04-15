#!/bin/bash
#SBATCH --comment "Hello ROMEO!"
#SBATCH -J "TEST 1"

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=02:30:00

#SBATCH -n 4
#SBATCH -N 1
#SBATCH --gres=gpu:2


#./ex12 -ksp_monitor_true_residual -pc_type none -m 200 -n 200 -log_summary
#./ex7 -ksp_monitor_true_residual -log_summary

make runb
