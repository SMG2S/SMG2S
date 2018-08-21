#!/bin/bash
#SBATCH --comment "Hello ROMEO!"
#SBATCH -J "VALIDATION"

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=02:30:00

#SBATCH -n 1

EXEC=./powerInverse.exe
N=100
L=10
TEST_TOL=0.00001
DEGREE=4

LENGTH=$(awk 'NR==2{print $1}' vector.txt)
echo "Test Eigenvalues number = "${LENGTH}

for((i=3;i<=${LENGTH}+2;i++))
do
	real=$(awk 'NR=='$i'{print $2}' vector.txt)
	imag=$(awk 'NR=='$i'{print $3}' vector.txt)
	srun -n 1 ${EXEC} -n ${N} -l ${L} -eps_monitor_conv -eps_power_shift_type constant -st_type sinvert -exact_value ${real}+${imag}i -test_tol ${TEST_TOL} -degree ${DEGREE}
done
