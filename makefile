FC=mpif90
SUBDIR=./utils
OBJDIR=./obj
MODDIR=./modules

UTIL_LIB=${SUBDIR}/lib/libmylib.a
MYLIB=-I ${SUBDIR}/modules -L${SUBDIR}/lib -lmylib
FLAGS=-O2 -J $(MODDIR) ${MYLIB} # ${LAP}
# setenv LAP "-I/opt/intel/composer_xe_2011_sp1.7.256/mkl/include/intel64/ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5"

OBJ=${OBJDIR}/random_number_mod.o \
    ${OBJDIR}/global_constant.o \
    ${OBJDIR}/input_data.o \
    ${OBJDIR}/potentials.o \
    ${OBJDIR}/relative_potential.o \
    ${OBJDIR}/level_density.o \
    ${OBJDIR}/coupling_matrix.o \
    ${OBJDIR}/coupled_channels.o \
    ${OBJDIR}/calc_profile.o \
    ${OBJDIR}/scat_noncoll.o

a.out: ${UTIL_LIB} ${OBJ}
	${FC} -o a.out ${OBJ} ${FLAGS}
${UTIL_LIB}: utils
	cd utils && make
${OBJDIR}/random_number_mod.o : random_number_mod.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/random_number_mod.o random_number_mod.f90
${OBJDIR}/global_constant.o : global_constant.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/global_constant.o global_constant.f90
${OBJDIR}/input_data.o : input_data.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/input_data.o input_data.f90
${OBJDIR}/potentials.o : potentials.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/potentials.o potentials.f90
${OBJDIR}/relative_potential.o : relative_potential.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/relative_potential.o relative_potential.f90
${OBJDIR}/level_density.o : level_density.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/level_density.o level_density.f90
${OBJDIR}/coupling_matrix.o : coupling_matrix.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/coupling_matrix.o coupling_matrix.f90
${OBJDIR}/coupled_channels.o : coupled_channels.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/coupled_channels.o coupled_channels.f90
${OBJDIR}/calc_profile.o : calc_profile.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/calc_profile.o calc_profile.f90
${OBJDIR}/scat_noncoll.o : scat_noncoll.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/scat_noncoll.o scat_noncoll.f90

.PHONY : clean
clean :
	$(RM) ${OBJDIR}/*.o ${MODDIR}/*.mod *~
