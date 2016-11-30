! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.


! This program uses few .f90 procedures to read the input
! and then for QMC calculation calls lisaqmc subroutine.
! See lisaqmc.f

program qmc
	use kinds
	use io
	use system
	
	implicit none

	real(DP) :: U_dmft(5,20,20)
	integer :: i,j,k	
		
	call read_input()	

	!U matrix
	U_dmft = 0.d0
	call read_umatrix()
	do i=1,nsigma
		do j=1, block_size(i)*(2/nspin)
			do k=1, block_size(i)*(2/nspin)
				U_dmft(i,j,k) = DBLE(U(i,j,k))
			end do
		end do
	end do
	
	call lisaqmc(beta,ecut,time_slices,nsweeps,nbin,nwarmup,aaa,nspin,nsigma,block_size,iqmc,U_dmft)
	
	call program_end()
end program qmc

