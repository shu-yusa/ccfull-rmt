IFC=`which mpif90`
FLAGS=${LAP} ${MYLIB}#-O2 -check all -traceback -fpe0
# 環境変数LAP, MYLIBの定義
# setenv LAP "-I/opt/intel/composer_xe_2011_sp1.7.256/mkl/include/intel64/ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5"
# setenv MYLIB "-I ${HOME}/usr/include/module -L${HOME}/usr/lib -lmylib"
# libmylib.soはディレクトリlib以下のファイルを自分で作成したライブラリ化したもの.

# できた実行ファイルは
# > exe.sh (並列プロセス数) 実行ファイル名
# により実行していた.
# 使用するマシンはmachinefileにより指定.
a.out : random_number_mod.o global_constant.o input_data.o potentials.o \
        relative_potential.o level_density.o coupling_matrix.o coupled_channels.o calc_profile.o \
	  scat_exact_noncoll_ver6.o
	${IFC} -o a.out random_number_mod.o global_constant.o input_data.o potentials.o \
        relative_potential.o level_density.o coupling_matrix.o coupled_channels.o calc_profile.o \
	  scat_exact_noncoll_ver6.o ${FLAGS}
random_number_mod.o : random_number_mod.f90
	${IFC} ${FLAGS} -c random_number_mod.f90
global_constant.o : global_constant.f90
	${IFC} ${FLAGS} -c global_constant.f90
input_data.o : input_data.f90
	${IFC} ${FLAGS} -c input_data.f90
potentials.o : potentials.f90
	${IFC} ${FLAGS} -c potentials.f90
relative_potential.o : relative_potential.f90
	${IFC} ${FLAGS} -c relative_potential.f90
level_density.o : level_density.f90
	${IFC} ${FLAGS} -c level_density.f90
coupling_matrix.o : coupling_matrix.f90
	${IFC} ${FLAGS} -c coupling_matrix.f90
coupled_channels.o : coupled_channels.f90
	${IFC} ${FLAGS} -c coupled_channels.f90
calc_profile.o : calc_profile.f90
	${IFC} ${FLAGS} -c calc_profile.f90
scat_exact_noncoll_ver6.o : scat_exact_noncoll_ver6.f90
	${IFC} ${FLAGS} -c scat_exact_noncoll_ver6.f90

.PHONY : clean
clean :
	$(RM) *.o *.mod


