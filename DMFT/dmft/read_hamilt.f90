! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine read_hamilt
! Reads hamiltonian from 'hamilt' file
	use kinds
	use io
	use system, only: nkp, wk, H, hdim, nspin
	use constants
	implicit none

	integer :: ispin, ikp, i, j, dummy
	real(DP) :: HRe, HIm, summ
	
	open (unit = file_hamilt, file = 'hamilt', status = 'old', form = &
       				'formatted', err = 101, iostat = ios)

	read(file_hamilt,*) nkp, dummy
	if(dummy /= hdim) call error("Check hdim in dmft.in and hamilt",13)
	
	allocate(wk(nkp))
	allocate(H(nspin,nkp,hdim,hdim))
	
	H = CZERO
	wk = 0.d0
	
	do ispin = 1, nspin
		do ikp=1, nkp
			read(file_hamilt,*), wk(ikp)
			do i = 1, hdim
				do j=1, hdim
					read(file_hamilt,*) HRe, HIm
					H(ispin,ikp,i,j) = HRe + CI*HIm
				end do
			end do 
		end do
	end do
	
	if(chatter) then
		write(stdout,*) 'HMLT successfully readed'
		write(stdout,*) 'Sum of k-weights:', SUM(wk)
	end if
		

	! normalization of k-points weight
	wk(:) = SUM(wk)/(1.d0*nkp*(2.d0/nspin))
	!wk(:) = wk(:)/SUM(wk)
	
	close(file_hamilt)
	
	return
	
	101 call error('Error opening hamilt', ios)

	
end subroutine read_hamilt

