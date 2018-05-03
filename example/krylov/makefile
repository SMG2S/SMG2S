##
ALL: blib exec
EXEC  = krylov
CFLAGS	= -O0 -O3 -L/home/xinzhewu/workstations/SMG2Spy/SMG2S/install/lib -I/home/xinzhewu/workstations/SMG2Spy/SMG2S/install/include -lsmg2s2c
LDFLAGS = -L${SMG2S}/lib
CLEANFILES  = ${EXEC}
OFILES = krylov.o

NTIMES = 1

#################      SMG2S Flags       ########################
SMG2S = -smg2s
LBANDWIDTH = -lbandwidth 4
SIZE=-size 200
NBONE=-continous1 3
#SPECTRA= -sptr vector.txt
SMG2S_FLAGS = ${SMG2S} ${LBANDWIDTH} ${SIZE} ${NBONE} ${SPECTRA}

#################       MPI Flags        ########################
MPI_NODES = 2

#################       GMRES Flags       ########################

RESTART_MAX =40
GMRES_PRECISION= 1e-10
GMRES_RESTART= ${RESTART_MAX}
GMRES_MONITOR= -ksp_monitor_true_residual
KSP_MAX_ITS = 20000
PC_TYPE = none
#PC_TYPE=jacobi -pc_jacobi_type rowmax -pc_jacobi_abs
#PC_TYPE = asm -pc_asm_block 2 -pc_asm_type none
#PC_TYPE = sor -pc_sor_its 5 -pc_sor_omega 0.179
#CUDA_TYPE = -mat_type aijcusparse -vec_type cuda
#GMRES_FT = -GMRES_FT
GMRES_FLAGS= -ksp_rtol 1e-100 -ksp_divtol 1e1000 -ksp_max_it ${KSP_MAX_ITS} -pc_type ${PC_TYPE} -ksp_atol ${GMRES_PRECISION} -ksp_gmres_restart ${GMRES_RESTART}\
		${GMRES_MONITOR} ${LSA_GMRES} ${GMRES_NB_NODES} -ntimes ${NTIMES} ${CUDA_TYPE}\
		${GMRES_FT} ${INITIAL_GUESS_SEQ_RHS} -initial_guess_power ${INITIAL_GUESS_POWER}

GLSA_FLAGS = ${GMRES_FLAGS} ${SMG2S_FLAGS}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


blib:
	@for a in $(DIRS); do \
        if [ -d $$a ]; then \
        echo "processing folder $$a"; \
	$(MAKE) -C $$a;\
        fi; \
        done;
	@echo "$$a Done!"
	-@echo "Completed building libraries"
	-@echo "========================================="


distclean :
	-@echo "Cleaning application and libraries"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH}  PETSC_DIR=${PETSC_DIR} clean
	-${RM} ${OFILES}
	-@echo "Finised cleaning application and libraries"
	-@echo "========================================="

exec: krylov.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	@${CLINKER} -I${PETSC_DIR}/include ${PETSC_LIB} -L${PETSC_DIR}/${PETSC_ARCH}/lib -g -v -o ${EXEC} krylov.o
	-@echo "Completed building application"
	-@echo "========================================="

runsmg2s:
		-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./krylov  ${GLSA_FLAGS}
