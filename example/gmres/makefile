ALL: blib exec 

#compilation and various flags
EXEC    = gmres.exe
CLEANFILES  = ${EXEC}
OFILES= ${wildcard ./*.o}
CXXFLAGS = -O3 -std=c++0x

###Tuning Parameters###

MPI_NODES=2
GMRES_RESTART=200
KSP_MONITOR=-ksp_monitor_true_residual
RTOL=1e-100
DIVTOL=1e1000
MAX_ITE=20000
PC=none
ATOL=1e-10
MAT=utm300_300x300.dat
DEGREE = 2
N = 1000
L = 3


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

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

exec: gmres_example.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@${CLINKER} -o ${EXEC} gmres_example.o ${PETSC_LIB}
	-@echo "Completed building application"
	-@echo "========================================="

#-ksp_monitor_true_residual -eps_monitor

runa:
	-@${MPIEXEC} -np ${MPI_NODES} ./gmres.exe -ksp_gmres_restart ${GMRES_RESTART} ${KSP_MONITOR} -ksp_rtol ${RTOL} \
		-ksp_divtol ${DIVTOL} -ksp_max_it ${MAX_ITE} -pc_type ${PC} -ksp_atol ${ATOL} -degree ${DEGREE} -n ${N} -l ${L} -log_view 

