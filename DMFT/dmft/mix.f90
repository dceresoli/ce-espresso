! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

program mix
	use kinds
	use io
	use system
	use constants
	
	implicit none
	
	integer :: i,j,k,isigma,iomega, isize, info, ispin
	integer, allocatable :: ipiv(:)
	complex(DP), allocatable :: Gw(:,:,:,:), Gbath(:,:,:,:), sigma_old(:,:,:,:,:), work(:,:)
	complex(DP) :: s
	logical :: old, need_smooth, file_exists
	real(DP) :: DR,DI, DUM, omega, d, stot
	real(DP), allocatable :: shift(:,:), dmax(:,:), sum(:,:)
	
	call read_input()
	
	isize = MAXVAL(block_size(:))
	
	allocate(shift(nsigma,isize))
	allocate(ipiv(isize))
	allocate(work(isize,isize))
	
	allocate(Gw(nsigma,0:nomega,isize,isize))
	Gw = CZERO
	
	allocate(Gbath(nsigma,0:nomega,isize,isize))
	Gbath = CZERO
	
	allocate( sigma(nspin,nsigma,0:nomega,isize,isize) )
	sigma = CZERO
	
	allocate( sigma_old(nspin,nsigma,0:nomega,isize,isize) )
	sigma_old = CZERO
	
	
	allocate(dmax(isize,isize), sum(isize,isize))
	dmax = 0.d0
	sum = 0.d0
	
	open (unit = file_Gomega, file = 'gw', status = 'old', form = 'formatted', err = 100, iostat = ios)
	open (unit = file_Gomegaout, file = 'gw.out', status = 'old', form = 'formatted', err = 101, iostat = ios)
	open (unit = file_sigma, file = 'sigma', status = 'unknown', form = 'formatted', err = 102, iostat = ios)
	open (unit = file_moments, file = 'moments', status = 'old', form = 'formatted', err = 103, iostat = ios)
	
	file_exists = .false.
	inquire (file='sigma.old',exist=file_exists)
	if(file_exists) then
		open (unit = file_sigma_old, file = 'sigma.old', status = 'old', form = 'formatted', err = 103, iostat = ios)
		old = .true.
	else
		old = .false.
	end if
	
!	some checks
    if (nomega > 3*time_slices) then
        need_smooth=.true.
    else
        need_smooth=.false.
    endif
	
	do ispin = 1, nspin

		do isigma=1,nsigma
		
			do i=1,block_size(isigma)
				do j=1,block_size(isigma)
					do iomega = 0, nomega
						read(file_Gomega,*) DUM, DR, DI
						Gbath(isigma,iomega,i,j) = dcmplx(DR,DI)
						read(file_Gomegaout,*) DUM, DR, DI
						Gw(isigma,iomega,i,j) = dcmplx(DR,DI)
						if(old) then
							read(file_sigma_old,*) DUM, DR, DI
							sigma_old(ispin,isigma,iomega,i,j) = dcmplx(DR,DI)
						end if
					end do
				end do
			end do
			
			!skip few lines
			do i=1,block_size(isigma)
				do j=1,block_size(isigma)
					read(file_moments,*) 
				end do
			end do
			
			do i=1,block_size(isigma)
				do j=1,block_size(isigma)
					read(file_moments,*) 
				end do
			end do
			
			!now read
			do i=1,block_size(isigma)
				read(file_moments,*) shift(isigma,i)
			end do
					
			do iomega = 0,nomega

			  call ZGETRF(isize,isize,Gbath(isigma,iomega,:,:),isize,ipiv,info)
			  call ZGETRI(isize,Gbath(isigma,iomega,:,:),isize,ipiv,work,isize,info)	
			  call ZGETRF(isize,isize,Gw(isigma,iomega,:,:),isize,ipiv,info)
			  call ZGETRI(isize,Gw(isigma,iomega,:,:),isize,ipiv,work,isize,info)
!........calculate sigma
				
			  sum = 0.d0
			  
			  do i=1, block_size(isigma)
				do j=1, block_size(isigma)
					s = Gbath(isigma,iomega,i,j)-Gw(isigma,iomega,i,j)
					
					if (i.eq.j) s = s - DCMPLX(shift(isigma,i),0.d0)
					!write(stdout,*) s
					d = cdabs(sigma_old(ispin,isigma,iomega,i,j)-s)
					sum(i,j) = sum(i,j)+d*d
					if (d.gt.dmax(i,j)) dmax(i,j)=d
					
					sigma(ispin,isigma,iomega,i,j)=mixing*s+(1.d0-mixing)*sigma_old(ispin,isigma,iomega,i,j)
			  	enddo
			  enddo
      
			end do
			
			if (need_smooth) then
				write(stdout,*) 'need smooth'
				do i=1,block_size(isigma)
					do j=1,block_size(isigma)
						call smooth(sigma(ispin,isigma,:,i,j),nomega,time_slices)
			   		enddo
			  	enddo
			end if
				 
			write(stdout,*)'L_inf norm:'
			do i=1,block_size(isigma)
				write(stdout,*)(dmax(i,j),j=1,block_size(isigma))
			end do
			write(6,*)

			stot = 0.d0
			write(6,*)'L_2 norm:'
			do i=1, block_size(isigma)
				do j=1, block_size(isigma)
			   		sum(i,j)=sum(i,j)/DBLE(nomega+1)
			   		stot = stot+sum(i,j)
			   		sum(i,j)=sqrt(sum(i,j))
			  	enddo
				write(stdout,*)(sum(i,j),j=1,block_size(isigma))
			enddo
			write(6,*)
			stot=sqrt(stot/DBLE(block_size(isigma)**2))

			write(6,*)':DEV_max:',MAXVAL(dmax(:,:))
			write(6,*)':DEV_nor:',stot

		
			!and output
			do i=1,block_size(isigma)
				do j=1,block_size(isigma)
					do iomega = 0, nomega
						omega=DBLE(2*iomega+1)*PI/beta
						write(file_sigma,*) omega , dreal(sigma(ispin,isigma,iomega,i,j)),imag(sigma(ispin,isigma,iomega,i,j))
					enddo
				enddo
			enddo
			
		end do

	end do
	
	
	
	close(file_Gomega)
	close(file_Gomegaout)
	close(file_sigma)
	if(file_exists) close(file_sigma_old)
	close(file_moments)
	
	deallocate(Gw)
	deallocate(Gbath)
	deallocate(sigma)
	deallocate(sigma_old)
	deallocate(shift)
	
	deallocate(dmax)
	deallocate(sum)
	deallocate(work)
	deallocate(ipiv)
	
100 call error('Error opening gw', ios)
101 call error('Error opening gw.out', ios)
102 call error('Error opening sigma', ios)
103 call error('Error opening sigma.old', ios)
104 call error('Error opening moments', ios)
end program mix

subroutine smooth(ss,Im,LL)
	use kinds
	implicit none
	
	complex(DP) :: ss(0:Im),tmp(0:Im),sum,ti, tr
	integer :: i,j,n,l2, Im, LL

	l2=LL/2

	do n=0,l2
	 tmp(n)=ss(n)
	enddo

	do i=1,l2
	 n=l2+i
	 sum=(0.d0,0.d0)
	 do j=-l2,l2
	  sum=sum+ss(n+j)
	 enddo
	 tmp(n)=sum/(2*l2+1)
	 tmp(n)=((l2-i)*ss(n)+i*tmp(n))/DBLE(l2)
	enddo

	do n=2*l2+1,Im-l2
	 sum=(0.d0,0.d0)
	 do j=-l2,l2
	  sum=sum+ss(n+j)
	 enddo
	 tmp(n)=sum/DBLE(2*l2+1)
	enddo

	do n=Im-l2+1,Im
	 tr=DBLE(tmp(Im-l2))
	 ti=imag(tmp(Im-l2))*DBLE(Im-l2)/DBLE(n)
	 tmp(n)=tr+(0.,1.)*ti
	enddo

	do n=0,Im
	 ss(n)=tmp(n)
	enddo
	return
    end
