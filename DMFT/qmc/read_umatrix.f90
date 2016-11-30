! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine read_umatrix
! Reads interaction matrix from umatrix file.
	use kinds
	use constants
	use system
	use io
	
	implicit none
	
	integer :: i,j,k
	
	
	allocate(U(nsigma,MAXVAL(block_size(:))*(2/nspin),MAXVAL(block_size(:))*(2/nspin)))

	U = 0.d0

	open (unit = file_umatrix, file = 'umatrix', status = 'old', form = &
					'formatted', err = 100, iostat = ios)
					
	do i=1,nsigma
		do j=1,block_size(i)*(2/nspin)
				read(file_umatrix,*) (U(i,j,k),k=1,block_size(i)*(2/nspin))
		end do
	end do
				
	close(file_umatrix)
		
	return

100 call error('Error opening umatrix', ios)

end subroutine read_umatrix
