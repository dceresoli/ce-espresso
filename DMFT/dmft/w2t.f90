! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine w2t(Gw,isize,xm1,xm2)

	use kinds
	use io
	use system
	use constants
	implicit none
	
	integer :: isigma, isize, i, j, k,l, it, iomega
	complex(DP), intent(in) :: Gw(0:nomega,isize,isize), xm1(isize,isize), xm2(isize,isize)
	complex(DP), allocatable :: Gt(:,:,:), xm0(:,:)
	real(DP) :: beta_pi, tau, omega,tail1,tail2
	complex(DP) :: sum,cdummy,dummy
	
	
	allocate(Gt(time_slices,isize,isize))
	allocate(xm0(isize,isize))
	
	xm0(:,:) = CZERO
	do i=1,isize
		xm0(i,i) = dcmplx(1.d0,0.d0)
	end do

    beta_pi=PI/beta
    do it=1, time_slices
        tau=(it-1)*beta/dble(time_slices)
!.........tail functions
        tail1=0.d0
        tail2=0.d0
        
        do iomega=0,nomega
            omega=(2.d0*iomega+1)*beta_pi
            tail1=tail1+dcos(omega*tau)/(omega*omega)
            tail2=tail2+dsin(omega*tau)/(omega*omega*omega)
        enddo
        
        tail1=beta*(1.d0-2.d0*(tau/beta))/4.d0-2.d0*tail1/beta
        tail2=beta**2*(tau/beta)*(1.d0-(tau/beta))/4.d0-2.d0*tail2/beta
        
!..........numerical summation to the cut-off          
        do k=1,isize
            do l=1,isize
                sum=CZERO
                do j=0,nomega
                    omega=(2.d0*j+1)*beta_pi
                    cdummy=CI*dmod(omega*tau,2.d0*PI)
                    dummy=(Gw(j,k,l)+xm0(k,l)*CI/omega)*zexp(-cdummy)+dconjg(Gw(j,l,k)+xm0(k,l)*CI/omega)*zexp(cdummy)
                    sum=sum+dummy
                end do
!........add assymptotic sums
	            Gt(it,k,l)=sum/beta-xm1(k,l)*tail1+xm2(k,l)*tail2-xm0(k,l)*0.5d0
            end do
        end do
    end do
    
    do i=1,isize
        do j=1,isize
            do k=1,time_slices
!.....sign change for QMC definition of G(t)
                write(file_Gtau,*)(k-1),-dreal(Gt(k,i,j)), -dimag(Gt(k,i,j))
            end do
        end do
    end do
    
    deallocate(Gt)
    deallocate(xm0)
 
end subroutine
