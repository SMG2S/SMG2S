ALL: blib exec2

#compilation and various flags

#SMG2S = SMG2S INSTALL PATH, which is already set to .bashrc

DIRS    = 
EXEC    = view
CFLAGS	= -L${SMG2S}/lib -I${SMG2S}/include -lsmg2s2c
FFLAGS	= 
CXXFLAGS	=
FPPFLAGS	=
CLEANFILES  = ${EXEC}
OFILES= ${wildcard ./*.o}

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

exec2: c_test.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@${CLINKER} -o ${EXEC} c_test.o ${PETSC_LIB}
	-@echo "Completed building application"
	-@echo "========================================="
#-ksp_monitor_true_residual -eps_monitor



