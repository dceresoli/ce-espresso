#!/bin/sh

QE_PATH=$HOME/Codes/qe-espresso/trunk/bin #path to QE CVS version bin directory
PROC_NUMBER=1 #Number of processes for parallel execution
PARA="mpirun" #prefix for parallel execution


#date
#echo 'Starting scf calculation'
#$PARA -np $PROC_NUMBER $QE_PATH/pw.x < scf.in > scf.out
#
#echo 'Starting nscf calculation'
#$PARA -np $PROC_NUMBER $QE_PATH/pw.x < nscf.in > nscf.out
#
#echo 'Starting Hamiltonian generation'
#$PARA -np 1 $QE_PATH/wannier_ham.x < hamilt.in > hamilt.out

echo 'Starting DMFT loop'
maxiter=20 # number of DNFT iterations
for times in $(seq 1 $maxiter); do

	echo DMFT
	$PARA -np 1 $QE_PATH/dmft.x
	
	cp gtau gtau.$times
	cp gw gw.$times

	echo QMC
	$PARA -np $PROC_NUMBER $QE_PATH/qmc.x

	cp gtau.out gtau.out.$times
	cp gw.out gw.out.$times

	mv sigma sigma.old
	echo MIX

	$PARA -np 1 $QE_PATH/mix.x

	cp sigma sg$times

	echo $times ITERATION OF $maxiter FINISHED
done

echo 'Starting spectral function calculation'

$PARA -np 1 $QE_PATH/maxent.x

echo 'Resulting spectra is in dos file'
date

