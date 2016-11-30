! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine print_input
	use kinds
	use io
	use system
	use constants
	implicit none

	integer :: i,j

	write(stdout,*)
	write(stdout,'(5x,a16,i2)') 'Number of spins:', nspin
	write(stdout,'(5x,a23,2x,i3,a,i3)') 'Hamiltonian dimentions:', hdim,'x',hdim
	write(stdout,'(5x,a19,i5)') 'Number of k-points:', nkp
	write(stdout,'(5x,a28,i2)') 'Number of iteraction blocks:', nsigma
	
	do i=1,nsigma
		write(stdout,'(7x,a4,i3,a3,i3)') 'from', block_start(i), ' to', block_start(i)+block_size(i)-1
	end do
	
	write(stdout,*)
	write(stdout,'(5x,a22,i4)') 'Number of time-slices:', time_slices
	write(stdout,'(5x,a5,f6.2)') 'Beta:', beta
	write(stdout,'(5x,a26,i4)') 'Number of Matsubara freq.:', nomega
	
	if(fix_mu) then
		write(stdout,'(5x,a22,f6.3)') 'Fixed chem. potential:', mu
	else 
		write(stdout,'(5x,a20,f6.3)') 'Number of electrons:', nel
	end if
	
	write(stdout,*)
	
end subroutine print_input

