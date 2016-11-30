! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine asymptotics
	use kinds
	use io
	use system
	use constants
	implicit none

	integer :: ispin, isigma, iomega, i, j
	
	integer :: nfreq ! we will calculate asymptotics from nfreq largest frequencies

	if(chatter) write(stdout,*) 'Sigma asymptotics calculation start'
	
	nfreq = 10 ! This is a some start
	
	
	allocate( sigma_one(nspin,nsigma,MAXVAL(block_size),MAXVAL(block_size)) )
	allocate( sigma_inf(nspin,nsigma,MAXVAL(block_size),MAXVAL(block_size)) )
	
	sigma_one = CZERO
	sigma_inf = CZERO
	
	do ispin = 1, nspin
		do isigma = 1, nsigma
			do iomega = nomega-nfreq+1,nomega
				do i=1, block_size(isigma)
					do j=1, block_size(isigma)
						sigma_one(ispin,isigma,i,j) = sigma_one(ispin,isigma,i,j) + &
							0.5d0*(sigma(ispin,isigma,iomega,i,j)-dconjg(sigma(ispin,isigma,iomega,j,i)))*DBLE(2*iomega-1)*PI*CI/DBLE(beta)
						sigma_inf(ispin,isigma,i,j) = sigma_inf(ispin,isigma,i,j) + &
							0.5d0*(sigma(ispin,isigma,iomega,i,j)+dconjg(sigma(ispin,isigma,iomega,j,i)))
					end do	
				end do
				
!				sigma_one(ispin,isigma,:,:) = sigma_one(ispin,isigma,:,:) + 0.5d0*(sigma(ispin,isigma,iomega,:,:)+dconjg(sigma(ispin,isigma,iomega,:,:)))
!				sigma_one(ispin,isigma,:,:) = sigma_one(ispin,isigma,:,:) + dreal(sigma(ispin,isigma,iomega,:,:))
			end do
			sigma_one(ispin,isigma,:,:) = sigma_one(ispin,isigma,:,:)/DBLE(nfreq)
			sigma_inf(ispin,isigma,:,:) = sigma_inf(ispin,isigma,:,:)/DBLE(nfreq)
		end do
	end do
	
	write(stdout,*) 'Sigma(inf)',dreal(sigma_inf)
	write(stdout,*) 'Sigma(one)',dreal(sigma_one)
	

	if(chatter) write(stdout,*) 'Sigma asymptotics calculation end'

end subroutine asymptotics

