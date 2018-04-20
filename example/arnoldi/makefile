ALL: blib exec 

#compilation and various flags
EXEC    = ./slepcbench.exe
CLEANFILES  = ${EXEC}
OFILES= ${wildcard ./*.o}
CXXFLAGS = -O3 -std=c++0x

###Tuning Parameters###

MPI_NODES=2
ARNOLDI_PRECISION=1e-1
ARNOLDI_NBEIGEN= 18
ARNOLDI_MONITOR = -eps_monitor_conv
ARNOLDI_NCV = 100
LOG_VIEW = -log_view
DEGREE = 2
N = 1000
L = 3


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

exec: arnoldi_example.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@${CLINKER} -o ${EXEC} arnoldi_example.o ${PETSC_LIB} ${SLEPC_LIB} -L${SLEPC_DIR}/${PETSC_ARCH}/lib
	-@echo "Completed building application"
	-@echo "========================================="

#-ksp_monitor_true_residual -eps_monitor

runa:
	-@${MPIEXEC} -np ${MPI_NODES} ./slepcbench.exe -eps_ncv ${ARNOLDI_NCV} -eps_type arnoldi -eps_true_residual -eps_largest_imaginary -eps_nev ${ARNOLDI_NBEIGEN} -eps_tol ${ARNOLDI_PRECISION} ${ARNOLDI_MONITOR} -eps_max_it 50 -degree ${DEGREE} -n ${N} -l ${L} -log_view 

