include ../../make.sys

F90FLAGS = -I../modules 

OBJS = \
../modules/modules.o \
program_end.o \
read_hamilt.o \
read_sigma.o \
asymptotics.o \
nelectrons.o \
chem_pot.o \
print_input.o \
gloc.o \
w2t.o

default: dmft.x mix.x

dmft.x:	$(OBJS) dmft.o 
	$(LD) -o $@ $(OBJS) dmft.o $(F90FLAGS) $(LAPACK_LIBS)
	cd ../../bin ; ln -fs ../DMFT/dmft/$@

mix.x: mix.o
	$(LD) -o $@ mix.o ../modules/modules.o $(F90FLAGS) $(LAPACK_LIBS)
	cd ../../bin ; ln -fs ../DMFT/dmft/$@

clean: 
	rm -f *.o *.mod


#dmft.o : modules.o
	
