! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine gloc()
	use kinds
	use io
	use system
	use constants
	implicit none
	
	complex(DP), allocatable :: T(:,:), T2(:,:), Hk(:,:), Hk2(:,:)
	complex(DP), allocatable :: GLT(:,:), GL(:,:),Gw(:,:,:,:)
	complex(DP), allocatable :: xm1(:,:,:), xm2(:,:,:), ym1(:,:), ym2(:,:)
	complex(DP), allocatable :: work(:), work2(:)
	complex(DP) :: sum, omega
	real(DP), allocatable :: shift(:,:)
	
	integer,allocatable :: ipiv(:)
	integer :: i,j,k,l,iomega,isigma,ispin,ik,isize, first,last,info
	
	real(DP) :: omax,beta_pi
	
	isize = MAXVAL(block_size(:))
	
	allocate(T(hdim,hdim))
	allocate(T2(hdim,hdim))
	allocate(Hk(hdim,hdim))
	allocate(Hk2(hdim,hdim))
	allocate(GLT(hdim,hdim))
	allocate(GL(isize,isize))
	allocate(Gw(nsigma,0:nomega,isize,isize))
	allocate(xm1(nsigma,isize,isize))
	allocate(xm2(nsigma,isize,isize))
	allocate(ym1(isize,isize))
	allocate(ym2(isize,isize))
	
	allocate(shift(nsigma,isize))
	
	allocate(work(hdim))
	allocate(work2(isize))
	allocate(ipiv(hdim))
	
	open(file_Gomega,file='gw',status='unknown',form='formatted')
	open(file_Gtau,file='gtau',status='unknown',form='formatted')
	open(file_moments,file='moments',status='unknown',form='formatted')
	
	write(6,*)'***********************************************'
    write(6,*)'GREEN FUNCTION CALCULATION'
    write(6,*)'***********************************************'
    write(6,*)
    
	beta_pi=pi/beta
    omax=(2*nomega+1)*beta_pi


	do ispin = 1,nspin
	
		Hk(:,:)=CZERO
		Hk2(:,:)=CZERO
		T = CZERO
		T2 = CZERO
		Gw = CZERO
		xm1 = CZERO
		xm2 = CZERO
		ym1 = CZERO
		ym2 = CZERO
		
		shift = 0.d0

!.........preparation for moment calculation        
		do ik=1,nkp
			do i=1,hdim
				do j=1,hdim
					T(i,j)=H(ispin,ik,i,j)
				enddo
				T(i,i)=T(i,i)-mu
			enddo
			
			do isigma=1,nsigma
				first = block_start(isigma)
				last = first+block_size(isigma)-1
				T(first:last,first:last) =  T(first:last,first:last) + sigma_inf(ispin,isigma,:,:)
			end do
				
			do i=1,hdim
				do j=1,hdim
					sum=CZERO
					do l=1,hdim
						sum=sum+T(i,l)*T(l,j)
					enddo
					T2(i,j)=sum
				enddo
			enddo
	!........Hk=sum_k (H(k)-mu+eta(inf))
	!........Hk2=sum_k (H(k)-mu+eta(inf))^2
			do i=1,hdim
				do j=1,hdim
					Hk(i,j)=Hk(i,j)+T(i,j)*wk(ik)
					Hk2(i,j)=Hk2(i,j)+T2(i,j)*wk(ik)
				enddo
			enddo
		end do
		
	!.......moments of local & bath green functions
	!.......ym1, ym2 1st and 2nd moment of the local spectal function
	!.......xm1, xm2 1st and 2nd moment of the local bath function
		do isigma = 1,nsigma

			do i=1,block_size(isigma)
				do j=1,block_size(isigma)
					ym1(i,j)=Hk(block_start(isigma)-1+i,block_start(isigma)-1+j)
					ym2(i,j)=sigma_one(ispin,isigma,i,j)+Hk2(block_start(isigma)-1+i,block_start(isigma)-1+j)
					xm1(isigma,i,j)=ym1(i,j)-sigma_inf(ispin,isigma,i,j)
				end do
				shift(isigma,i)=dreal(xm1(isigma,i,i))
			end do
			
			do i=1,block_size(isigma)
				do j=1,block_size(isigma)
					sum = CZERO
					do k=1,block_size(isigma)
						sum=sum+xm1(isigma,i,k)*xm1(isigma,k,j)-ym1(i,k)*ym1(k,j)			
					end do
					xm2(isigma,i,j)=sum+ym2(i,j)-sigma_one(ispin,isigma,i,j)
				end do
			end do
			
		end do

	!-------GREEN FUNCTION - OMEGA LOOP
		do iomega = 0, nomega

			omega=(2.d0*iomega+1)*beta_pi*CI

			GLT(:,:)=CZERO
	!********k-point summation

			do ik=1,nkp
			
				T(:,:) = -H(ispin,ik,:,:)
				do i=1,hdim
					T(i,i) = T(i,i)+omega+mu
				end do
				
	!.........add -self-energy 
				do isigma=1,nsigma
					first = block_start(isigma)
					last = first+block_size(isigma)-1
					T(first:last,first:last) = &
						 T(first:last,first:last) - sigma(ispin,isigma,iomega,:,:)
				end do
	!.........invertion
				call ZGETRF(hdim,hdim,T,hdim,ipiv,info)
				call ZGETRI(hdim,T,hdim,ipiv,work,hdim,info)
				
				do i=1,hdim
					do j=1,hdim
						GLT(i,j)=GLT(i,j)+T(i,j)*wk(ik)
					enddo
				enddo
			end do

	!........local (small size) GF for each interacting site
			do isigma = 1, nsigma
			
				GL = CZERO
				
				do i=1, block_size(isigma)
					do j=1, block_size(isigma)
						GL(i,j) = GLT(block_start(isigma)-1+i,block_start(isigma)-1+j)
					end do
				end do
	!.........inversion
				call ZGETRF(isize,isize,GL,isize,ipiv,info)
				call ZGETRI(isize,GL,isize,ipiv,work2,isize,info)	 
	!.........add self-energy to G^-1 -> G_bath^-1 
				do i=1, block_size(isigma)
					do j=1, block_size(isigma)
						GL(i,j) = GL(i,j)+sigma(ispin,isigma,iomega,i,j)
						if(i.eq.j) GL(i,j) = GL(i,j)+shift(isigma,i)
					end do
				end do
	!.........inversion
				call ZGETRF(isize,isize,GL,isize,ipiv,info)
				call ZGETRI(isize,GL,isize,ipiv,work2,isize,info)
				
				do i=1, block_size(isigma)
					do j=1, block_size(isigma)
						Gw(isigma,iomega,i,j) = GL(i,j)
					end do
				end do
			  
				
			end do

		end do
	!---------
		do isigma=1,nsigma
			do i=1, block_size(isigma)
				do j=1, block_size(isigma)
					sum = CZERO
					do k=1, block_size(isigma)
						sum=sum+xm1(isigma,i,k)*xm1(isigma,k,j)
					end do
					xm2(isigma,i,j)=xm2(isigma,i,j)-sum
				end do
			end do
			
			do i=1,block_size(isigma)
				xm1(isigma,i,i)=xm1(isigma,i,i)-shift(isigma,i)
			end do
			
			do i=1, block_size(isigma)
				do j=1, block_size(isigma)
					sum = CZERO
					do k=1, block_size(isigma)
						sum=sum+xm1(isigma,i,k)*xm1(isigma,k,j)
					end do
					xm2(isigma,i,j)=xm2(isigma,i,j)-sum
				end do
			end do
			
			do i=1, block_size(isigma)
				do j=1, block_size(isigma)					
					do iomega=0,nomega
						omega = (2*iomega+1)*beta_pi
						write(file_Gomega,*) dreal(omega), dreal(Gw(isigma,iomega,i,j)), imag(Gw(isigma,iomega,i,j))
					end do
				end do
			end do
			
			do i=1,block_size(isigma)
				do j=1,block_size(isigma)
					write(file_moments,*) dreal(xm1(isigma,i,j)),imag(xm1(isigma,i,j))
				end do
			end do
			
			do i=1,block_size(isigma)
				do j=1,block_size(isigma)
					write(file_moments,*) dreal(xm2(isigma,i,j)),imag(xm2(isigma,i,j))
				end do
			end do
			
			do i=1,block_size(isigma)
				write(file_moments,*) shift(isigma,i)
			end do
				
			
	!.........Fourier transform w->tau		
			call w2t(Gw(isigma,:,:,:),isize,xm1(isigma,:,:),xm2(isigma,:,:))
			
		end do
	end do	
	
	deallocate(T)
	deallocate(T2)
	deallocate(Hk)
	deallocate(Hk2)
	deallocate(GLT)
	deallocate(GL)
	deallocate(Gw)
	deallocate(xm1)
	deallocate(xm2)
	deallocate(ym1)
	deallocate(ym2)
	
	deallocate(shift)
	
	deallocate(work)
	deallocate(work2)
	deallocate(ipiv)
	
	close(file_Gomega)
	close(file_Gtau)
	close(file_moments)
	
end subroutine gloc

