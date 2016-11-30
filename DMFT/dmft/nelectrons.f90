! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine nelectrons(ispin,a0,mu_in)
	use kinds
	use io
	use system
	use constants
	implicit none

	integer, intent(in) :: ispin
	real(DP), intent(in) :: mu_in
	real(DP), intent(out) :: a0
	
	integer :: i, j, ik, iomega, info, first, last, isigma
	complex(DP), allocatable :: T(:,:),work(:)
    complex(DP) :: ham_trace, sigma_trace, omega, suma
    real(DP) :: om,tmp,tail
    integer, allocatable :: ipiv(:)
	
	allocate(T(hdim,hdim))
	allocate(ipiv(hdim))
	allocate(work(hdim))

	om=(2*nomega+1)*PI/beta

!.....asymptotics for sigma
!..... Tr [S ReSigma]
	sigma_trace=CZERO
    ham_trace=CZERO
    
	do ik=1,nkp
		do i=1,hdim
	 		ham_trace=ham_trace+dreal(H(ispin,ik,i,i))*wk(ik)
	 	end do
		
		do j=1, nsigma
			do i=1,block_size(j)
	  			sigma_trace=sigma_trace+dreal(sigma_inf(ispin,j,i,i))*wk(ik)
	  		end do
        end do
	end do
	
    tmp=0.d0
    do iomega=0,nomega
    	tmp=tmp+1.d0/(2.d0*iomega+1.d0)**2
    enddo
    
    tail=beta*(0.25d0-2.d0/(PI*PI)*tmp)
!.....end asymptotics

	a0=0.d0
    do ik=1, nkp
    	suma=CZERO
        do iomega=0,nomega
	    	omega=(2*iomega+1)*PI*CI/beta
	    	
	    	T(:,:) = -H(ispin,ik,:,:)
	      	do i=1,hdim
	      		T(i,i)=T(i,i)+mu_in+omega
	    	enddo
	    	
	    	do isigma=1,nsigma
	    		first = block_start(isigma)
	    		last = first+block_size(isigma)-1
	    		T(first:last,first:last) = &
	    			 T(first:last,first:last) - sigma(ispin,isigma,iomega,:,:)
	    	end do
	    	
       		call ZGETRF(hdim,hdim,T,hdim,ipiv,info)
	   		call ZGETRI(hdim,T,hdim,ipiv,work,hdim,info)
	    	do i=1,hdim
	     		suma=suma+T(i,i)
 			end do
 		end do
 		a0=a0+2.d0*suma*wk(ik)
 	end do
	
	a0=a0/beta+real(hdim)/2.d0+(hdim*mu_in-ham_trace-sigma_trace)*tail
	a0=(2.d0/(1.d0*nspin))*a0
!	spin factor
!	a0=2.d0*a0
!	write(*,*) 'mu_in', mu_in
!	write(*,*) 'a0', a0
!       write(6,*)'tail contribution:', 2.d0*(hdim*mu_in-ham_trace-sigma_trace)*tail
!       write(6,*)

	deallocate(T)
	deallocate(ipiv)
	deallocate(work)

	return
end subroutine nelectrons

