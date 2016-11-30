include ../../make.sys

F90FLAGS = -I../modules

OBJS = \
../modules/modules.o \
../dmft/program_end.o \
maxent_calc.o

default: maxent.x

maxent.x:	$(OBJS) maxent.o 
	$(LD) -o $@ $(OBJS) maxent.o $(F90FLAGS) $(LAPACK_LIBS) ../../Modules/libqemod.a ../../flib/flib.a
	cd ../../bin ; ln -fs ../DMFT/maxent/$@

clean: 
	rm -f *.o *.mod
