ALL: blib exec

#compilation and various flags
EXEC    = powerInverse.exe
CLEANFILES  = ${EXEC}
OFILES= ${wildcard ./*.o}
CXXFLAGS = -O3 --std=c++0x

MDIR=./data
###Tuning Parameters###

MPI_NODES=1
ARNOLDI_PRECISION=1e-10
#MAT=utm300.mtx_300x300_3155nnz
#MAT=matBlock_nb_300_90000x90000_1.88984e+06_nnz
#MAT = matLine_nb_600_180000x180000_1.91215e+06_nnz
MAT = EBMG_matrix_nb_10_10x10_58_nnz
ARNOLDI_NBEIGEN= 10
EPS_MONITOR = -eps_monitor_conv
ARNOLDI_NCV = 11
SHIFT_TYPE = constant
ST_TYPE = sinvert
TEST_VALUE =2
TEST_TOL = 0.00001
DEGREE = 2
SIZE = 10
LBANDWIDTH = 3
#LOG_VIEW = -log_view

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

blib :
	-@echo "BEGINNING TO COMPILE LIBRARIES "
	-@echo "========================================="
	-@${OMAKE}  PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ACTION=libfast tree
	-@echo "Completed building libraries"
	-@echo "========================================="

distclean :
	-@echo "Cleaning application and libraries"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH}  PETSC_DIR=${PETSC_DIR} clean
	-${RM} ${OFILES}
	-@echo "Finised cleaning application and libraries"
	-@echo "========================================="	

exec: powerInverse.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@${CLINKER} -o ${EXEC} powerInverse.o ${SLEPC_LIB} -L${SLEPC_DIR}/${PETSC_ARCH}/lib
	-@echo "Completed building application"
	-@echo "========================================="

runa:
	-@${MPIEXEC} -np ${MPI_NODES} ./powerInverse.exe ${LOG_VIEW} ${EPS_MONITOR} -eps_power_shift_type ${SHIFT_TYPE} -st_type ${ST_TYPE} -exact_value ${TEST_VALUE} -test_tol ${TEST_TOL} -degree ${DEGREE} -n ${SIZE} -l ${LBANDWIDTH}



