#!/bin/bash
#SBATCH --comment "Hello ROMEO!"
#SBATCH -J "TEST 1"

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=100:30:00

#SBATCH -n 1
#SBATCH -N 1

MPI_NODES=1
MDIR=./data
#MAT=Rectangular_nb_50_50x50_1278_nnz
#MAT=Rectangular_nb_50_50x50_1306_nnz
#MAT=Rectangular_nb_50_50x50_1288_nnz
#MAT=Rectangular_nb_50_50x50_1280_nnz
#MAT=Rectangular_nb_50_50x50_1281_nnz
#MAT=Rectangular_nb_50_50x50_1305_nnz
#MAT=Rectangular_nb_50_50x50_1288_nnz
#MAT=Rectangular_nb_50_50x50_1279_nnz
#MAT=Rectangular_nb_50_50x50_1288_nnz
#MAT=Rectangular_nb_50_50x50_1299_nnz
#MAT=Rectangular_nb_50_50x50_1296_nnz
#MAT=Rectangular_nb_100_100x100_4134_nnz
#MAT=Rectangular_nb_100_100x100_4155_nnz
MAT=Rectangular_nb_100_100x100_4122_nnz
#MAT=Rectangular_nb_10_10x10_59_nnz
#LOG_VIEW=-log_view
SHIFT_TYPE=constant
ST_TYPE=sinvert
TEST_TOL=0.0001
EXACT_VALUE=5+5i
DEGREE=4
EPS_MONITOR=-eps_monitor_conv
COMPLEX=1
for I in {5..500..5}
do
  if [ ${COMPLEX} -eq 1 ]
  then
  	EXACT_VALUE=$I+${I}i
  else
	EXACT_VALUE=$I
  fi
  echo EXACT_VALUE = $EXACT_VALUE
  srun  -n ${MPI_NODES} ./powerInverse.exe -mfile ${MDIR}/${MAT} ${EPS_MONITOR} ${LOG_VIEW} -eps_power_shift_type ${SHIFT_TYPE} -st_type ${ST_TYPE} -exact_value ${EXACT_VALUE} -test_tol ${TEST_TOL} -degree ${DEGREE}
done

echo All done!


#srun  -n ${MPI_NODES} ./powerInverse.exe -mfile ${MDIR}/${MAT} ${LOG_VIEW} -eps_power_shift_type ${SHIFT_TYPE} -st_type ${ST_TYPE} -exact_value ${EXACT_VALUE} -test_tol ${TEST_TOL} -degree ${DEGREE}
