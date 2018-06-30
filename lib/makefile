SHELL=/bin/sh
AR=/usr/bin/ar

SRCDIR=$(HOME)/usr/src/fortran
INCDIR=$(HOME)/usr/include
LIBDIR=$(HOME)/usr/lib64
GFMODDIR=$(INCDIR)/module_gf64
OBJDIRGF=$(LIBDIR)/object_files/gfortran

GFCMPL=/usr/bin/gfortran
LAP=-I/opt/intel/composer_xe_2011_sp1.7.256/mkl/include/intel64/ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
GFFLAG=-O2 -J $(GFMODDIR) 

SHARED_LIB=libmylib64.so
STATIC_LIB=libmylib64gf.a

OBJIC=$(OBJDIRIC)/Coulomb_mod2.o \
	$(OBJDIRIC)/Coulomb_mod3.o \
	$(OBJDIRIC)/Eigen_mod.o \
	$(OBJDIRIC)/angular_coupling_mod.o \
	$(OBJDIRIC)/derive_mod.o \
	$(OBJDIRIC)/gamma_mod.o \
	$(OBJDIRIC)/gaus_int.o \
	$(OBJDIRIC)/gaus_weights.o \
	$(OBJDIRIC)/lin_coup.o \
	$(OBJDIRIC)/print_array_mod.o \
	$(OBJDIRIC)/special_fct_mod.o \
	$(OBJDIRIC)/fitting_mod.o \

OBJGF=$(OBJDIRGF)/Coulomb_mod2.o \
	$(OBJDIRGF)/Coulomb_mod3.o \
	$(OBJDIRGF)/Eigen_mod.o \
	$(OBJDIRGF)/angular_coupling_mod.o \
	$(OBJDIRGF)/derive_mod.o \
	$(OBJDIRGF)/gamma_mod.o \
	$(OBJDIRGF)/gaus_int.o \
	$(OBJDIRGF)/gaus_weights.o \
	$(OBJDIRGF)/lin_coup.o \
	$(OBJDIRGF)/print_array_mod.o \
	$(OBJDIRGF)/special_fct_mod.o \

lib : $(STATIC_LIB)

$(STATIC_LIB) : $(OBJGF)
	$(AR) rcs $(LIBDIR)/$(STATIC_LIB) $(OBJGF)


$(OBJDIRGF)/Coulomb_mod2.o : $(SRCDIR)/Coulomb_mod2.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/Coulomb_mod2.o \
	$(SRCDIR)/Coulomb_mod2.f90
$(OBJDIRGF)/Coulomb_mod3.o : $(SRCDIR)/Coulomb_mod3.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/Coulomb_mod3.o \
	$(SRCDIR)/Coulomb_mod3.f90
$(OBJDIRGF)/Eigen_mod.o : $(SRCDIR)/Eigen_mod.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/Eigen_mod.o \
	$(SRCDIR)/Eigen_mod.f90
$(OBJDIRGF)/angular_coupling_mod.o : $(SRCDIR)/angular_coupling_mod.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/angular_coupling_mod.o \
	$(SRCDIR)/angular_coupling_mod.f90
$(OBJDIRGF)/derive_mod.o : $(SRCDIR)/derive_mod.f90
	$(GFCMPL) -c $(GFFLAG)  \
	-o $(OBJDIRGF)/derive_mod.o \
	$(SRCDIR)/derive_mod.f90
$(OBJDIRGF)/gamma_mod.o : $(SRCDIR)/gamma_mod.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/gamma_mod.o \
	$(SRCDIR)/gamma_mod.f90
$(OBJDIRGF)/gaus_int.o : $(SRCDIR)/gaus_int.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/gaus_int.o \
	$(SRCDIR)/gaus_int.f90
$(OBJDIRGF)/gaus_weights.o : $(SRCDIR)/gaus_weights.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/gaus_weights.o \
	$(SRCDIR)/gaus_weights.f90
$(OBJDIRGF)/lin_coup.o : $(SRCDIR)/lin_coup.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/lin_coup.o \
	$(SRCDIR)/lin_coup.f90
$(OBJDIRGF)/print_array_mod.o : $(SRCDIR)/print_array_mod.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/print_array_mod.o \
	$(SRCDIR)/print_array_mod.f90
$(OBJDIRGF)/special_fct_mod.o : $(SRCDIR)/special_fct_mod.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIRGF)/special_fct_mod.o \
	$(SRCDIR)/special_fct_mod.f90

.PHONY : clean
clean :
	rm $(OBJGF)
