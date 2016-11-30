include ../../make.sys

F90 = mpif90
F77 = mpif90 
F90FLAGS = -I../modules

OBJS = \
../modules/modules.o \
lisaqmc.o \
cputim.o \
nfourier2.o \
sym.o \
program_end.o \
read_umatrix.o


default: qmc.x

qmc.x:	$(OBJS) qmc.o
	$(LD) -o $@ $(OBJS) qmc.o  $(F90FLAGS) $(LAPACK_LIBS)
	cd ../../bin ; ln -fs ../DMFT/qmc/$@

clean: 
	rm -f *.o *.mod

	
