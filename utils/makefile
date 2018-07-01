AR=/usr/bin/ar

SRCDIR=./src
LIBDIR=./lib
MODDIR=./modules
OBJDIR=./obj

GFCMPL=gfortran
# LAP=-I/opt/intel/composer_xe_2011_sp1.7.256/mkl/include/intel64/ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
GFFLAG=-O2 -J $(MODDIR) 

STATIC_LIB=libmylib.a

OBJGF=$(OBJDIR)/coulomb.o \
	$(OBJDIR)/eigen.o \
	$(OBJDIR)/angular_coupling.o \
	$(OBJDIR)/diff.o \
	$(OBJDIR)/gammaln_fct.o \
	$(OBJDIR)/gaus_int.o \
	$(OBJDIR)/gaus_weights.o \
	$(OBJDIR)/lin_coup.o \
	$(OBJDIR)/print_array.o \
	$(OBJDIR)/special_fct.o \

lib : $(STATIC_LIB)

$(STATIC_LIB) : $(OBJGF)
	$(AR) rcs $(LIBDIR)/$(STATIC_LIB) $(OBJGF)


$(OBJDIR)/coulomb.o : $(SRCDIR)/coulomb.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/coulomb.o \
	$(SRCDIR)/coulomb.f90
$(OBJDIR)/eigen.o : $(SRCDIR)/eigen.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/eigen.o \
	$(SRCDIR)/eigen.f90
$(OBJDIR)/angular_coupling.o : $(SRCDIR)/angular_coupling.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/angular_coupling.o \
	$(SRCDIR)/angular_coupling.f90
$(OBJDIR)/diff.o : $(SRCDIR)/diff.f90
	$(GFCMPL) -c $(GFFLAG)  \
	-o $(OBJDIR)/diff.o \
	$(SRCDIR)/diff.f90
$(OBJDIR)/gammaln_fct.o : $(SRCDIR)/gammaln_fct.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/gammaln_fct.o \
	$(SRCDIR)/gammaln_fct.f90
$(OBJDIR)/gaus_int.o : $(SRCDIR)/gaus_int.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/gaus_int.o \
	$(SRCDIR)/gaus_int.f90
$(OBJDIR)/gaus_weights.o : $(SRCDIR)/gaus_weights.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/gaus_weights.o \
	$(SRCDIR)/gaus_weights.f90
$(OBJDIR)/lin_coup.o : $(SRCDIR)/lin_coup.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/lin_coup.o \
	$(SRCDIR)/lin_coup.f90
$(OBJDIR)/print_array.o : $(SRCDIR)/print_array.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/print_array.o \
	$(SRCDIR)/print_array.f90
$(OBJDIR)/special_fct.o : $(SRCDIR)/special_fct.f90
	$(GFCMPL) -c $(GFFLAG) \
	-o $(OBJDIR)/special_fct.o \
	$(SRCDIR)/special_fct.f90

.PHONY : clean
clean :
	rm $(OBJDIR)/*.o *~ $(MODDIR)/*.mod
