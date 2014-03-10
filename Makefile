ifeq ($(COMP),intel)
	FC=ifort
	FFLAGS=
	OFLAGS=-openmp
else ifeq ($(COMP),pgi)
	FC=pgf90
	FFLAGS=
	OFLAGS=-mp
else ifeq ($(COMP),ibm)
	FC=xlf90
	FFLAGS=
	OFLAGS=-qsmp=omp
else
	FC=gfortran
	FFLAGS=
	OFLAGS=-fopenmp
endif

OBJS = md.o

all : v0 v1 v2 v3 v4 v5 omp

v0 : 
	$(FC) $(FFLAGS) -o md0 md-orig.f90

v1 : 
	$(FC) $(FFLAGS) -o md1 md-v1.f90

v2 : 
	$(FC) $(FFLAGS) -o md2 md-v2.f90

v3 : 
	$(FC) $(FFLAGS) -o md3 md-v3.f90

v4 : 
	$(FC) $(FFLAGS) -o md4 md-v4.f90

v5 : 
	$(FC) $(FFLAGS) -o md5 md-v5.f90

omp :
	$(FC) $(FFLAGS) $(OFLAGS) -o mdo md-omp.f90

%.o :%.f90
	$(FC) $(FFLAGS) -c $<

