! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine read_sigma
! Reads sigma from 'sigma' file
	use kinds
	use io
	use system
	use constants
	implicit none

	integer :: ispin, isigma, iomega, i, j
	real(DP) :: SRe, SIm, dummy
	
	allocate( sigma(nspin,nsigma,0:nomega,MAXVAL(block_size),MAXVAL(block_size)) )
	sigma = CZERO
	
	open (unit = file_sigma, file = 'sigma', status = 'old', form = &
       				'formatted', err = 102, iostat = ios)
       				
	do ispin = 1, nspin
		do isigma = 1, nsigma
			do i=1, block_size(isigma)
				do j=1, block_size(isigma)
					do iomega = 0, nomega
						read(file_sigma,*) dummy, SRe, SIm
						sigma(ispin,isigma,iomega,i,j) = SRe + CI*SIm
					end do
				end do
			end do
		end do 
	end do
	
	if(chatter) write(stdout,*) 'Sigma successfully readed'

	close(file_sigma)
	
	return

! Zero sigma	
102 	write(stdout,*) 'Zero sigma assumed'
	sigma = CZERO
	return
	
end subroutine read_sigma

