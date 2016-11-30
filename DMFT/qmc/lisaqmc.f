C=======+=========+=========+=========+=========+=========+=========+=$
C
C This subroutine was created by Jan Kunes (jan.kunes@physik.uni-augsburg.de)
C on the base of original lisaqmc code by W. Krauth:
C http://www.physics.rutgers.edu/~udo/qmc/readme_lisaqmc.html
C 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine lisaqmc(i_beta,ecut,i_LL,maxt,nbin,n_warmup,aaa,
     & i_spinfac,nsigma,block_size,iqmc,i_U)
      include 'lisaqmc.dat' 
      include 'mpif.h'
      external cputim
      real ranw
      logical logbin,det_only
      integer i,iii,j,k,kk,kkk,l,l1,l2,n,nfast,nran,isi,isj,iw,ix
      integer imove,isig,ibin,isize,iorb,icpu,ierr
      integer nleft,nl,iposition,isum,index
      integer iunit,ios,natom_int,jatom,iat,in,iout
      integer naccept,nswap,nfastup,iter,irec,iatom
      integer maxt,nbin,n_warmup
      integer ifnumeric
      double precision greenm,greensum,gnew,gbin,gtmp
      double precision greensum0
      double precision xx,yy,xn,xnn,xnn0,xm0,xm1,xm2,sum,sumi,sumij
      double precision xm0_n,xm1_n,xm2_n
      double precision det,det0,detrat,dd,fac,dummy,dtmp,dtmp0
      double precision t_sweep,t_swap,t_measure,t_total,t0,t1,t2
      double precision tol,diff,ddt
      double precision aaa,xpi,xlten
      double precision xs
      double precision ecut,eu,deltau,dn0,facbin
      double precision facnorm,facnorm2,detsum,sigsum
      double precision xmm,xm,xxu,xxd
      double precision facnorm_tot,facnorm2_tot,detsum_tot
      double precision sigsum_tot,xnn_tot,xm_tot,facbin_tot
      double precision greensum_tot,gbin_tot
      double precision input_data
       
      dimension imove(2,12),isig(NORBSQ)
      dimension greenm(0:MTAU-1),greensum_tot(NORB,0:MTAU-1),
     &          greensum(NORB,0:MTAU-1),greensum0(NORB,0:MTAU-1),
     &          gnew(MTAU,MTAU,NORB),gbin(NORB,0:MTAU-1),
     &          gtmp(0:MTAU-1),gbin_tot(NORB,0:MTAU-1)
      dimension xs(0:MTAU,NORBSQ)
      dimension xx(NORB,NORB),yy(NORB,NORB)
      dimension xn(NORB),xnn(NORB,NORB),xnn0(NORB,NORB)
      dimension xnn_tot(NORB,NORB),xm_tot(0:MTAU-1)
      dimension det(2),det0(2),dtmp(2),dtmp0(2)
      dimension xmm(0:mtau-1),xm(0:mtau-1),xxu(mtau),xxd(mtau)
      dimension input_data(5000)
      CHARACTER*80 DEFFN,FNAME,STATUS,FORM

C Array contains Nint for each interaction block
      integer block_size, i_LL, nsigma, i_spinfac, iqmc
      double precision i_beta, i_U
      dimension block_size(10), iqmc(10), i_U(5,20,20)

      nran(i)=mod(int(i*ranw(Idum)),i) + 1 

C========================================================================
C      input and output files linked
C========================================================================
      CALL MPI_Init( ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)  
      write(*,*)myid,' process from ',nprocs,' is alive'
      if (myid.eq.0) then    
       write(*,*)'# of cpus:',nprocs
       open(3,file='gtau',status='old',form='formatted')
       open(8,file='gw.out',status='unknown',form='formatted')
       open(6,file='qmc.out',status='unknown',form='formatted')
       open(7,file='qmc.scf',status='unknown',form='formatted')
       open(15,file='iqmc',status='unknown',form='formatted')
       open(17,file='gtau.out',status='unknown',form='formatted')
       open(18,file='gtau.nosym',status='unknown',form='formatted')
       open(31,file='moments',status='old',form='formatted')
       open(50,file='gbin',status='unknown',form='formatted')
       open(60,file='xnbin',status='unknown',form='formatted')
       open(70,file='umatrix',status='unknown',form='formatted')
!       open(82,file='ns',status='unknown',form='formatted')
       binname='gbin'
       write(*,*)'end oppening files' 
      endif
c======================================================================
c     initial set-up
c======================================================================
        xpi=acos(-One)
        xlten=log(Ten)
c......calculation for int. block iatom        
      if (myid.eq.0) then
       read(15,*,iostat=ios)iatom
       if (ios.eq.0) then
        read(15,*)Idum
        write(*,*)'IDUM=',Idum,ranw(Idum)
       else
	iatom=1
        Idum=iqmc(1)
       endif
       rewind(15)
       write(15,*)iatom
       write(15,'(I10)')int(-ranw(Idum)*100000)
18     format(a9,5i9)
       Beta = DBLE(i_beta)
       LL = i_LL

        write(6,*)'aaa=',aaa

       if (nbin.gt.0) then
        logbin=.true.
       else
        logbin=.false.
       endif

       if (LL.gt.MTAU) stop 'too many time slices'
       Spinfac = 2/i_spinfac
       Im=int(Ecut*Beta/(Two*acos(-1.d0)))
       write(6,*)'Im=',Im
       write(6,*)'Beta=',Beta
       write(6,*)'# of dtau intervals L=',LL
       write(6,*)'length of the bin:',nbin
       write(6,*)'# of warmup runs:',n_warmup
       natom_int = nsigma
       if (iatom.gt.natom_int) stop 'iatom > natom_int'
       do jatom=1,iatom
	
         Nint = block_size(jatom)
         if (jatom.lt.iatom) then
          do i=1,Nint
           do j=1,Nint
            do l=1,LL
             read(3,*)
            enddo
           enddo
          enddo
         endif

         if (Nint*Spinfac.gt.norb) stop 'too many orbitals'

	  do i=1,Nint*Spinfac
	    do j=1,Nint*Spinfac
             U(i,j) = i_U(jatom,i,j)
	    end do
	  enddo
          do i=1,Nint
           do j=1,Nint
            read(31,*)xx(i,j)
           enddo
          enddo
          do i=1,Nint
           do j=1,Nint
            read(31,*)yy(i,j)
           enddo
          enddo
          do i=1,Nint
           read(31,*)shift(i)
          enddo

        enddo
        if (Spinfac.eq.2) then
         do i=1,Nint
          shift(i+Nint)=shift(i)
          do j=1,Nint
           xx(i+Nint,j+Nint)=xx(i,j)
           xx(i+Nint,j)=Zero
           xx(i,j+Nint)=Zero
           yy(i+Nint,j+Nint)=yy(i,j)
           yy(i+Nint,j)=Zero
           yy(i,j+Nint)=Zero 
          enddo
         enddo
        endif
        goto 1002
 1001   write(6,*)'error reading input',iat
 1002	write(6,*)'Interaction:',Nint
	write(6,*)'# of sweeps maxt=',maxt
        write(6,*)
	write(6,*)'Size of the interacting block:',nint
c       endif
c------end myid.eq.0
cc
cc    
cc----distribution of matrices among CPUs for inversion
      isize=10
      call MPI_PACK_SIZE(1,MPI_LOGICAL,MPI_COMM_WORLD,
     &iout,ierr)
      isize=isize+iout
      call MPI_PACK_SIZE(8,MPI_INTEGER,MPI_COMM_WORLD,
     &iout,ierr)
      isize=isize+iout
      in=1+Norb*(Norb+1)
      call MPI_PACK_SIZE(in,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,
     &iout,ierr)
      isize=isize+iout
      if (isize.gt.5000) then
        write(*,*)'isize = ',isize,' > 5000. 
     .    Increase dimension of input_data. Stop'
        stop
      endif
      endif
c end if (myid.eq.0)      
      call MPI_BCAST(isize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (myid.eq.0) then
        iposition = 0
        write(6,*)'isize=',isize
       call MPI_PACK(logbin,1,MPI_LOGICAL,input_data,isize,iposition,
     & MPI_COMM_WORLD,ierr)
       call MPI_PACK(LL,1,MPI_INTEGER,input_data,isize,iposition,
     & MPI_COMM_WORLD,ierr)
       call MPI_PACK(Spinfac,1,MPI_INTEGER,input_data,isize,iposition,
     & MPI_COMM_WORLD,ierr)
       call MPI_PACK(Idum,1,MPI_INTEGER,input_data,isize,iposition,
     & MPI_COMM_WORLD,ierr)
       call MPI_PACK(Nint,1,MPI_INTEGER,input_data,isize,iposition,
     & MPI_COMM_WORLD,ierr)
       call MPI_PACK(maxt,1,MPI_INTEGER,input_data,isize,iposition,
     & MPI_COMM_WORLD,ierr)
       call MPI_PACK(n_warmup,1,MPI_INTEGER,input_data,isize,iposition,
     & MPI_COMM_WORLD,ierr)
       call MPI_PACK(nbin,1,MPI_INTEGER,input_data,isize,iposition,
     & MPI_COMM_WORLD,ierr)
       
       call MPI_PACK(Beta,1,MPI_DOUBLE_PRECISION,input_data,isize,
     & iposition,MPI_COMM_WORLD,ierr)
       n=Norb*Norb
       call MPI_PACK(U,n,MPI_DOUBLE_PRECISION,input_data,isize,
     & iposition,MPI_COMM_WORLD,ierr)
       call MPI_PACK(Shift,NORB,MPI_DOUBLE_PRECISION,input_data,isize,
     & iposition,MPI_COMM_WORLD,ierr)
       call MPI_PACK(aaa,1,MPI_DOUBLE_PRECISION,input_data,isize,
     & iposition,MPI_COMM_WORLD,ierr)
      endif
!       write(6,*)'label 0'
      
!      if (myid.ne.0) write(6,*)'on CPU ',myid,' isize = ',isize       
       call MPI_BCAST(input_data,isize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
!       write(6,*)'label 1'
!       stop
!      else
!       call MPI_BCAST(input_data,isize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
      if (myid.ne.0) then
!       write(6,*)'label 2 '
!       stop
       iposition = 0
       call MPI_UNPACK(input_data,isize,iposition,logbin,1,MPI_LOGICAL,
     & MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,LL,1,MPI_INTEGER,
     & MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,Spinfac,1,MPI_INTEGER,
     & MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,Idum,1,MPI_INTEGER,
     & MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,Nint,1,MPI_INTEGER,
     & MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,maxt,1,MPI_INTEGER,
     & MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,n_warmup,1,
     & MPI_INTEGER,MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,nbin,1,MPI_INTEGER,
     & MPI_COMM_WORLD,ierr)

       call MPI_UNPACK(input_data,isize,iposition,Beta,1,
     & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
       n=Norb*Norb
       call MPI_UNPACK(input_data,isize,iposition,U,n,
     & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,Shift,NORB,
     & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
       call MPI_UNPACK(input_data,isize,iposition,aaa,1,
     & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      endif
      deltau=Beta/real(LL)
      n=0
      do i=2,Nint*Spinfac
       do j=1,i-1
        if (abs(U(i,j)).gt.1.D-5) then
         n=n+1
         dummy=exp(deltau*U(i,j)/Two)
         Xlambda(n)=log(dummy+sqrt(dummy**2-One))
        endif
      enddo
      enddo
      LSQ=n
      if (lsq.gt.norbsq) stop 'too many interactions'
      call fill_index
      if (myid.eq.0) then
       write(6,*)'# of interactions (pairs) LSQ=',LSQ
       write(6,*)'H-S field structure:'
       write(6,*)'  HSindex   orbital 1  orbital 2   Lambda(n)'
       do n=1,LSQ 
        i=ind(1,n)
        j=ind(2,n)
        write(6,1)n,i,j,Xlambda(n)
 1      format(3(i5,7x),f8.5)
       enddo
      endif
      if (myid.eq.0) call seedgen(Idum,nprocs,seed)
      call MPI_BCAST(seed,5001,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Idum=seed(myid+1)
      if (myid.eq.0) then
        write(6,*)'seeds:'
        do i=1,nprocs
         write(6,*)i,seed(i)
        enddo   
       endif
cc    INITIALIZATION
      call initial(det)
      do n=1,lsq
       isig(n)=0
       xs(0,n)=Zero
       do l=1,LL
        xs(l,n)=Zero
        isig(n)=isig(n)+Is(n,l)
       enddo
      enddo

      dn0=det(2)
      do l=0,LL
       gtmp(l)=zero
      enddo
      do 810 k=0,LL-1
      do 810 i=1,Nint
       gbin(i,k)=Zero
       greensum(i,k)=Zero
810   continue

      do j=1,Nint
       do i=1,Nint
        xnn(i,j)=Zero
       enddo
      enddo
	
      sigsum=zero
      detsum=zero
      facnorm=zero
      facnorm2=zero
      facbin=zero
      iter=0
      irec=0
      ibin=0
      naccept=0
      nswap=0
      nfastup=100
      if (logbin) then
       if (.not.(mod(maxt-n_warmup,nbin).eq.0)) then
        write(*,*)'maxt=',maxt,'  nbin=',nbin
        stop 'maxt not a multiple of nbin'
       endif
      endif

      t_sweep=Zero
      t_swap=Zero
      t_measure=Zero
c***********************************************************************
c		start simulation  
c***********************************************************************
         call cputim(t0)
         do 1000 iii=1,1000000000
          if (myid.eq.0) write(6,18)':RUN_nfs:',nfastup
c.........fast update loop          
          do 2000 kkk=1,nfastup
           if (iter.gt.maxt) goto 3120
           iter=iter+1
!           if (myid.eq.0) then
!            rewind(14)
!            write(14,*)iter
!           endif
         call cputim(t1)
c
c........try orbital swap
c
c         goto 1234
c         n=nran(12)
c         call swap(imove(1,n))
c         call swap(imove(2,n))
c         det0(1)=det(1)
c         det0(2)=det(2)
c         det_only=.true.
c         call update(det_only,gnew,det)
c         if ((det(2)-det0(2)).gt.20) then
c           det_only=.false.
c           call update(det_only,Green,det)
c           nswap=nswap+1
c         elseif ((det(2)-det0(2)).lt.-20) then
c           det(1)=det0(1)
c           det(2)=det0(2)
c           call swap(imove(2,n))
c           call swap(imove(1,n))     
c         else     
c          dd=det(1)/det0(1)*Ten**(det(2)-det0(2))    
c          dummy=abs(dd**(One-aaa))
c          if (ranw(Idum).lt.dummy/(One+dummy)) then
c           do i=1,LSQ
c            isig(i)=0
c            do l=1,LL
c             isig(i)=isig(i)+Is(i,l)
c            enddo
c           enddo
c           det_only=.false.
c           call update(det_only,Green,det)
c           nswap=nswap+1
c          else
c           det(1)=det0(1)
c           det(2)=det0(2)
c           call swap(imove(2,n))
c           call swap(imove(1,n))
c          endif 
c         endif
c 1234    continue
c
c........end swap  
c      
         call cputim(t2)
         t_swap=t_swap+t2-t1
         t1=t2
c
cc       sweep loop   
c      
        do 1900 kk=1,LL*LSQ
         k=nran(LL*LSQ)
         l=(k-1)/LSQ+1
         n=k-(l-1)*LSQ
c.........fc-weighted transition rate
         dd=detrat(l,n)
         dummy=abs(dd**(One-aaa))
c........ heat bath
         if (ranw(Idum).lt.dummy/(One+dummy)) then
          det(1)=det(1)*dd
c.........adjust detetminat
          if (det(1).eq.Zero) goto 160
  110     if (abs(det(1)).ge.One) goto 120
          det(1)=Ten*det(1)
          det(2)=det(2)-One
          goto 110
  120     continue
  130     if (abs(det(1)).lt.Ten) goto 140
          det(1)=det(1)/Ten
          det(2)=det(2)+One
          goto 130
  140     continue  
  150     continue
  160     continue
c........end adjust
c.......fast update
          call record(l,n)
          naccept=naccept+1
          isig(n)=isig(n)-2*Is(n,l)
         endif
1900    continue
c
c.........end sweep
          call cputim(t2)
	  t_sweep=t_sweep+t2-t1
          t1=t2
c.........measurement, skipped for first n_warmup runs   
        if (iter.eq.n_warmup) then
         dtmp(1)=log(det(1))/xlten+det(2)
         dtmp(2)=myid
         call MPI_reduce(dtmp,dtmp0,1,MPI_2DOUBLE_PRECISION,
     &   MPI_MAXLOC,0,MPI_COMM_WORLD,ierr)
         icpu=int(dtmp0(2)+0.1)
         if (myid.eq.0) 
     &   write(6,*)'largest dlog(det)=',dtmp0(1),' on CPU #',icpu
         call MPI_BCAST(icpu,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

         isize=MTAU*NORBSQ
         call MPI_BCAST(Is,isize,MPI_INTEGER,icpu,
     &   MPI_COMM_WORLD,ierr)

         det_only=.false.
         call update(det_only,green,det)
         dn0=det(2)
         if (myid.eq.0) write(*,*)'det0:',det(1),det(2)
        elseif (iter.gt.n_warmup) then	  
         irec=irec+1
         fac=(det(1)**aaa)*ten**(aaa*(det(2)-dn0))
         do n=1,LSQ
          xs((isig(n)+LL)/2,n)=xs((isig(n)+LL)/2,n)+fac
         enddo
         sigsum=sigsum+det(1)/abs(det(1))*fac
         detsum=detsum+(log(det(1))/xlten+det(2))*fac
         facnorm=facnorm+fac
         facnorm2=facnorm2+fac*fac
         facbin=facbin+fac

         do i=1,Nint
          do l=0,LL-1
           greenm(l)=Zero
          enddo
	  do l=0,LL-1
	   do l1=1,LL-l
	    greenm(l)=greenm(l)+
     &      green(l1+l,l1,i)
           enddo
	   do l1=1,l
            greenm(l)=greenm(l)-
     &      green(l1,LL-l+l1,i)
           enddo
           greenm(l)=greenm(l)/real(LL)
          enddo
          do l=0,LL-1
           greensum(i,l)=greensum(i,l)+greenm(l)*fac
           gbin(i,l)=gbin(i,l)+greenm(l)*fac
          enddo
         enddo
	     
c........<n1*n2>
         do j=1,Nint
          do i=1,Nint
           sum=Zero
           do l=1,LL
            sum=sum+(One-green(l,l,i))*(One-green(l,l,j))
           enddo
	   xnn(i,j)=xnn(i,j)+sum/real(LL)*fac
          enddo
         enddo

c..... dynamic local susceptibility <Mz(tau)Mz(0)>
         do l=1,LL
          xxu(l)=Zero
          xxd(l)=Zero
          xmm(l-1)=Zero
          do i=1,Nint/2
           xxu(l)=xxu(l)+One-green(l,l,i)
           xxd(l)=xxd(l)+One-green(l,l,i+Nint/2)
          enddo
         enddo

         do l=0,LL-1 
          do l1=1,LL-l
           xmm(l)=xmm(l)+xxu(l+l1)*xxu(l1)+xxd(l+l1)*xxd(l1)-
     &     xxu(l+l1)*xxd(l1)-xxd(l+l1)*xxu(l1)
          enddo 
          do l1=LL-l+1,LL
           xmm(l)=xmm(l)+xxu(l+l1-LL)*xxu(l1)+xxd(l+l1-LL)*
     &     xxd(l1)-xxu(l+l1-LL)*xxd(l1)-xxd(l+l1-LL)*xxu(l1)
          enddo
         enddo 

         do l=0,LL-1
          do i=1,Nint
           do l1=1,LL-l
            xmm(l)=xmm(l)-green(l1+l,l1,i)*
     &      green(l1,l1+l,i)
           enddo
           do l1=1,l
            xmm(l)=xmm(l)-green(l1,LL-l+l1,i)*
     &      green(LL-l+l1,l1,i)
           enddo
          enddo
          xm(l)=xm(l)+xmm(l)/real(LL)*fac
         enddo
         if (logbin.and.(mod(irec,nbin).eq.0)) then
          ibin=ibin+1
          call MPI_reduce(facbin,facbin_tot,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          isize=NORB*MTAU
          call MPI_reduce(gbin,gbin_tot,isize,MPI_DOUBLE_PRECISION,
     &    MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if (myid.eq.0) then
           write(50,*)'bin',ibin,facbin_tot,dn0
           do i=1,Nint
            write(50,567)(gbin_tot(i,l)/facbin_tot,l=0,LL-1)
           enddo
          endif
          do i=1,Nint
           do l=0,LL-1
            gbin(i,l)=Zero
           enddo
          enddo
          facbin=Zero
         endif
        endif
       call cputim(t2)
       t_measure=t_measure+t2-t1
2000   continue
 567   format(200f12.8)
	      
cc
cc       update greens functions from scratch
cc
         det0(1)=det(1)
         det0(2)=det(2)
         det_only=.false.
         call update(det_only,gnew,det)
         ddt=abs(det(1)*Ten**(det(2)-det0(2))-det0(1))/abs(det(1))
         diff=ddt
         
         do i=1,Nint
          do l=1,LL
           do k=1,LL
            ddt=abs(Green(k,l,i)-gnew(k,l,i)) 
     &      /abs(Green(k,l,i)+gnew(k,l,i))
            diff=max(diff,ddt)
            Green(k,l,i)=gnew(k,l,i)
           enddo
          enddo
         enddo
         tol=1.d-6
         if (myid.eq.0) write(6,*)'fast update error:',diff, tol
         if (diff.gt.tol) then
          nfastup=max(1,nfastup/2)
         endif
         if (diff.lt.tol) nfastup=nfastup*2
1000     continue
3120  continue
        call cputim(t1)
        t_total=t1-t0
        if (myid.eq.0) then
         rewind(60)
         do k=1,LL
          write(60,*)(Is(n,k),n=1,LSQ)
         enddo
        endif
        
c====================================================================
c                      END OF SIMULLATION
c====================================================================
      call MPI_reduce(facnorm,facnorm_tot,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_reduce(facnorm2,facnorm2_tot,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_reduce(sigsum,sigsum_tot,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_reduce(detsum,detsum_tot,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      isize=MTAU*NORB
      call MPI_reduce(greensum,greensum_tot,isize,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      isize=NORB*NORB
      call MPI_reduce(xnn,xnn_tot,isize,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      isize=MTAU
      call MPI_reduce(xm,xm_tot,isize,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if (myid.eq.0) then
        write(*,*)'test:',BINNAME
       write(6,*)'***   REWEIGHTING STATISTICS:'
       write(6,*)'<f>=',facnorm_tot/(irec*nprocs)
       write(6,*)'f-r.m.s.=',sqrt((facnorm2_tot/irec/nprocs-
     & (facnorm_tot/irec/nprocs)**2)/(nprocs*irec-1))
       facnorm_tot=One/facnorm_tot
       write(6,*)'<sign>=',sigsum_tot*facnorm_tot 
       write(6,*)'<LOG_10(det)>=',detsum_tot*facnorm_tot
c.....calculation of <G(S)f(S)>
       write(17,*)irec*nprocs
       do i=1,Nint
        do l=0,LL-1
 	 greensum_tot(i,l)=-greensum_tot(i,l)*facnorm_tot 
         write(18,'(2f20.10)')l*beta/LL,greensum_tot(i,l) 
        enddo       
        write(18,'(2f20.10)')l*beta/LL,-One-greensum_tot(i,0)
        write(18,*)
        xn(i)=One+greensum_tot(i,0)
       enddo
c..........print occupations
       write(6,*)'orbital occupation'
       write(6,*)' i        <n_i> '     
       do i=1,Nint
        write(6,56)i,xn(i)
 56     format(i3,3x,f12.7)
 57     format(i3,3x,f10.4)
       enddo
       do i=1,Nint
        do j=1,Nint
         if (i.eq.j) then
          xnn_tot(i,i)=xn(i)
         else
	  xnn_tot(i,j)=xnn_tot(i,j)*facnorm_tot
         endif
        enddo
       enddo
c......unsymmetrized output
       write(6,*)'unsymmetrized output:'
       write(6,*)'double occupation'
       write(6,'(14i12)')(i,i=1,Nint)
       do i=1,Nint
        write(6,66)i,(xnn_tot(i,j),j=1,Nint)
       enddo
       write(7,*)
 66    format(i3,14f12.7)
 67    format(i3,14f10.4)
       write(*,*)'binname=',binname
       call symmetry(greensum_tot,xnn_tot,logbin,binname)
c........symmetrized output      
       do i=1,Nint/Spinfac
        do l=0,LL-1
         write(17,'(2f20.10)')l*beta/LL,greensum_tot(i,l)
        enddo
        write(17,'(2f20.10)')l*beta/LL,-One-greensum_tot(i,0)
       enddo
       write(7,*)'symmetrized output:'
       write(7,*)
       write(7,*)'orbital occupation:'
       write(7,*)' i        <n_i> '
       do i=1,Nint
        write(7,56)i,xnn_tot(i,i)
       enddo
       write(7,*)
       write(7,*)'double occupation'
       write(7,'(14i12)')(i,i=1,Nint)
       do i=1,Nint
        write(7,66)i,(xnn_tot(i,j),j=1,Nint)
       enddo
       write(17,*)
       do i=1,Nint
        write(17,*)(xnn_tot(i,j),j=1,Nint)
       enddo
       EU=0.d0
       do i=1,Nint
        do j=1,Nint
         EU=EU+U(i,j)*xnn_tot(i,j)
        enddo
       enddo
       EU=EU/2
c..... Fourier transformation
       write(6,*)'Fourier transform - moment spline'
       write(6,*)'Moments:'
c.....moment calculation
       do i=1,Nint/Spinfac
        do j=1,Nint/Spinfac
         if (i.eq.j) then
          do L=0,LL-1
           gtmp(L)=greensum_tot(i,L)
          enddo
          gtmp(LL)=-One-greensum_tot(i,0)
          xm0=-1.d0
          sumi=0.d0
          sumij=0.d0
          do k=1,Nint
           sumi=sumi+U(i,k)*xnn_tot(k,k)
           do l=1,Nint
            sumij=sumij+U(i,k)*U(i,l)*xnn_tot(k,l)
           enddo
          enddo
c..........xm1=M1 <= dG/dt=M1
          xm1=xx(i,i)
          xm1=xm1+sumi+shift(i)
c.........xm2=M2 <= d2G/dt2=-M2
          xm2=-yy(i,i)-shift(i)**2-Two*(xx(i,i)+shift(i))*sumi
     &    -sumij
          ifnumeric=0
          write(6,*)
          write(6,*)'orbital #',i
          write(6,*)'Derivatives(moments) 0 1 2 :    ',xm0,xm1,xm2
          write(6,*)'Derivatives(numeric) 0 1 2 :    ',xm0_n,xm1_n,xm2_n
          if(ifnumeric.eq.1)then 
           write(6,*)'use numeric derivatives'
          else
           write(6,*)'use moment derivatives'
          endif
          write(6,*)'sumi,sumij:',sumi,sumij
          write(6,*)'M10,M20: ',xx(i,i),yy(i,i)
          write(6,*)'M1 ,M2 : ',xm1,-xm2
          if(ifnumeric.eq.1)then
           call nfourier2(xm0_n,xm1_n,xm2_n,gtmp)            
          else
           call nfourier2(xm0,xm1,xm2,gtmp)
          endif 
         else
          do iw=0,IM
           write(8,*)(Two*iw+One)*xpi/Beta,Zero,Zero
          enddo
         endif
        enddo
       enddo
c.....calculation of M_z^2
       sum=Zero
       do i=1,Nint
        isi=1
        if (i.gt.(Nint/2)) isi=-1
        do j=1,Nint
         isj=1
         if (j.gt.(Nint/2)) isj=-1       
         sum=sum+xnn_tot(i,j)*isi*isj
        enddo
       enddo
!       write(95,*)Zero,sum
!       do l=1,LL-1
!        write(95,*)beta/LL*l,xm_tot(l)*facnorm_tot
!       enddo
!       write(95,*)beta,sum
       write(6,*)'<U-qmc>=',EU
       write(7,*)'<U-qmc>=',EU
       write(6,*)'sqrt(<m_z^2>)=',sqrt(sum)
       write(7,*)'sqrt(<m_z^2>)=',sqrt(sum)
c..........end Fourier transform
       write(6,'(a60)')'========================================'
       write(6,'(a20)')'         Summary'
       write(6,*)'        acc. steps     sweeps'
       write(6,18)':RUN_itr:',naccept,iter-1
       write(6,*)'acceptance probability'
       write(6,*)':RUN_acc:',naccept/real(LL*LSQ*irec)
       write(6,*)':RUN_swap:',nswap,nswap/real(irec)
       write(6,*)
       write(6,*)':CPU_sweep:',t_sweep,t_sweep/irec
       write(6,*)':CPU_swap:',t_swap,t_swap/irec
       write(6,*)':CPU_measurement:',t_measure
       write(6,'(a60)')'========================================'
      endif
      call MPI_FINALIZE(ierr)
	close(3)
	close(8)
	close(6)
	close(7)
	close(15)	
	close(17)
	close(18)
	close(31)
	close(50)
	close(60)
	close(70)
      end subroutine
C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: detrat.f  
C     TYPE   : function
C     PURPOSE: calculate the ratio of the new and
C              old determinants (cf. eq. (\ref{detrat})
C========+=========+=========+=========+=========+=========+=========+=$
      double precision function detrat(k,n)
      include 'lisaqmc.dat'
      integer i1,i2,n,k
      double precision fac1,fac2
        i1=ind(1,n)
        i2=ind(2,n)

        fac1=exp(-Two*Xlambda(n)*real(Is(n,k)))-One
        fac2=exp(Two*Xlambda(n)*real(Is(n,k)))-One

      detrat=(One+(One-Green(k,k,i1))*fac1)*
     &       (One+(One-Green(k,k,i2))*fac2)
      end
C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: swap.f
C     TYPE   : subroutine
C     PURPOSE: calculate the ratio of the new and
C              old determinants (cf. eq. (\ref{detrat})
C     I/O    :
C     VERSION: 30-Sep-95
C     COMMENT:
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine swap(n)
cc      exchange orbitals i1 <----> i2 by flipping and swapping the 
cc      corresponding fields
        include 'lisaqmc.dat'
        integer n,i1,i2,m,n1,n2,itmp
        integer l

        i1=ind(1,n)
        i2=ind(2,n)

        do m=1,i2-1
         n1=((i1-1)*(i1-2))/2+m
         n2=((i2-1)*(i2-2))/2+m
         do l=1,LL
          itmp=Is(n1,l)
          Is(n1,l)=Is(n2,l)
          Is(n2,l)=itmp
         enddo
        enddo

        do m=i2+1,i1-1
         n1=((i1-1)*(i1-2))/2+m
         n2=((m-1)*(m-2))/2+i2
         do l=1,LL
          itmp=-Is(n1,l)
          Is(n1,l)=-Is(n2,l)
          Is(n2,l)=itmp
         enddo
        enddo

        do m=i1+1,Nint
         n1=((m-1)*(m-2))/2+i1
         n2=((m-1)*(m-2))/2+i2
         do l=1,LL
          itmp=Is(n1,l)
          Is(n1,l)=Is(n2,l)
          Is(n2,l)=itmp
         enddo
        enddo

        do l=1,LL
         Is(n,l)=-Is(n,l)
        enddo

        end
C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: initial.f
C     TYPE   : subroutine
C     PURPOSE: read in initial configuration of bath Green's function
C              and of Ising spins, expand G(i-j) into matrix G(i,j).
C              invoke subroutine Update to calculate
C              Green's function for the initial choice of
C              Ising spins.
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine initial(det)
      include 'lisaqmc.dat'
      include 'mpif.h'
      real ranw
      logical det_only,restart
      integer i,j,k,l,n,ios,ix,isize,ierr
      double precision gtemp,gnew,dummy,det
      dimension gtemp(NORB,-MTAU+1:MTAU-1),gnew(MTAU,MTAU,NORB),
     &          det(2)

cc    read in G(sigma=0)
cc    
      if (myid.eq.0) then
       do i=1,Nint
        do j=1,Nint
         do k=0,LL-1
          if (i.eq.j) then
           read(3,*)dummy,gtemp(i,k)
          else
            read(3,*)dummy,dummy
          endif
         enddo
        enddo
       enddo
       rewind(3)
      endif
      isize=NORB*(2*MTAU-1)
      call MPI_BCAST(gtemp,isize,MPI_DOUBLE_PRECISION,0,
     & MPI_COMM_WORLD,ierr)
      
      if (Spinfac.eq.2) then
       do i=1,Nint
        do l=0,LL-1
         gtemp(i+Nint,l)=gtemp(i,l)
        enddo
       enddo
       Nint=Nint*Spinfac
      endif 

      do k=1,LL-1
       do i=1,Nint
        gtemp(i,-k)=-gtemp(i,LL-k)
       enddo
      enddo
cc
cc    initialize the hs fields
cc
      if (myid.eq.0) then
       read(60,*,iostat=ios)
       if (ios.eq.0) then
         rewind(60)
        do k=1,LL
         read(60,*,err=213,end=213)(Is(n,k),n=1,LSQ)
        enddo
        write(6,*)
        write(6,*)'***HS FIELD RESTART***'
        write(6,*)
        restart=.true.
        goto 214
       endif
213    continue
       restart=.false.
       write(6,*)
       write(6,*)'***RANDOM START****'
       write(6,*)
214    continue	
      endif
      call MPI_BCAST(restart,1,MPI_LOGICAL,0,
     & MPI_COMM_WORLD,ierr)

      if (restart) then
       isize=MTAU*NORBSQ
       call MPI_BCAST(Is,isize,MPI_INTEGER,0,
     & MPI_COMM_WORLD,ierr)

      else
       do k=1,LL
        do n=1,Lsq
         ix=mod(int(2*ranw(Idum)),2)
         Is(n,k)=(-1)**ix
        enddo
       enddo
      endif
cc
      do i=1,Nint
       do l=1,LL
        do k=1,LL
         Green0(k,l,i)=gtemp(i,k-l)
        enddo
       enddo
      enddo
      call update_ini(gnew)
      do i=1,Nint
       do l=1,LL
        do k=1,LL
         Green0(k,l,i)=Gnew(k,l,i)
        enddo
       enddo
      enddo
      det_only=.false.
      call update(det_only,green,det)   
      end
C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: ranw.f
C     TYPE   : real function
C     PURPOSE: produce uniformly distributed random numbers
C              following the algorithm of Mitchell and Moore
C========+=========+=========+=========+=========+=========+=========+=$
      real function ranw(idum)
      Parameter (Mbig=2**30-2, Xinvers=1./Mbig)
      data ibit/ 1/
      Integer IX(55)
      save
      if (ibit.ne.0) then
         ibit=0
cc
cc       fill up the vector ix with some random integers, which are
cc       not all even
cc       
         if (idum.eq.0) pause 'use nonzero value of idum'
         idum=abs(mod(idum,Mbig))
         ibit=0
         Ix(1)=871871
         Do i=2,55
            Ix(i)=mod(Ix(i-1)+idum,Ix(i-1))
            Ix(i)=max(mod(Ix(i),Mbig),idum)
         enddo
         j=24
         k=55
cc
cc       warm up the generator     
cc
         do i=1,1258
            Ix(k)=mod(Ix(k)+Ix(j),Mbig)
            j=j-1
            if (j.eq.0) j=55 
            k=k-1
            if (k.eq.0) k=55 
         enddo
      endif
cc
cc    this is where execution usually starts:
cc
      Ix(k)=mod(Ix(k)+Ix(j),Mbig)
      j=j-1
      if (j.eq.0) j=55 
      k=k-1
      if (k.eq.0) k=55 
      ranw=Ix(k)*Xinvers 
      end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: record.f  
C       TYPE   : subroutine 
C       PURPOSE: record changes of accepted move on the Green's function
C                (cf  eq. (\ref{fastupdate}))
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine record(k,n)
        include 'lisaqmc.dat'
        integer ir,i,l,l1,l2,k,n
        double precision  V,gl,gr
        dimension gl(MTAU),gr(MTAU)
cc
        do ir=1,2
         i=ind(ir,n)
         V=exp((-One)**ir*Two*Xlambda(n)*Is(n,k))-One
         V=V/(One+(One-Green(k,k,i))*V)
         do l=1,LL
          gl(l)=Green(l,k,i)
          gr(l)=V*Green(k,l,i)
         enddo
         gl(k)=gl(k)-One
         do l2=1,LL
          do l1=1,LL
           Green(l1,l2,i)=Green(l1,l2,i)+gl(l1)*gr(l2)
          enddo
         enddo
        enddo
cc
cc   update HS-spin
cc
        Is(n,k)=-Is(n,k)
        end


C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: update.f  
C     TYPE   : subroutine 
C     PURPOSE: calculate the Green's function 
C              for a given configuration of spins 
C              (in vector Is) from the Green's function
C              for spins set equal to zero  (eq. (\ref{inversion}))
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine update(det_only,gnew,det)
      include 'lisaqmc.dat'
      integer job,i,j,n,k,l,kk
      double precision a,ainv,v,adet
      double precision gnew,det
      double precision sum
      logical det_only
      dimension a(MTAU,MTAU),ainv(MTAU,MTAU),
     &          v(MTAU,NORB),adet(2)
      dimension gnew(MTAU,MTAU,NORB),det(2)
cc
cc    calculate the matrix a=1-(g-1)(exp(v')-1)
cc
        job=11
        do  n=1,Nint
	 do i=1,LL
	  v(i,n)=Zero
         enddo
        enddo
        det(1)=One
        det(2)=Zero

	do  l=1,LL
         do n=1,Lsq
	  i=ind(1,n)
	  j=ind(2,n)
	  v(l,i)=v(l,i)+real(Is(n,l))*Xlambda(n)
          v(l,j)=v(l,j)-real(Is(n,l))*Xlambda(n)
	 enddo
	enddo

      do i=1,Nint  
       do k=1,LL
         do l=1,LL
          a(k,l)=-Green0(k,l,i)*(exp(v(l,i))-One)
         enddo
         a(k,k)=One-(Green0(k,k,i)-One)*(exp(v(k,i))-One)
       enddo
       if (det_only) job=10
c
c------calculate inverse job=11 (inverse+det), job=10 (det only)
       call inverse(a,ainv,adet,job)
c
       det(1)=det(1)*adet(1)
       det(2)=det(2)+adet(2)
       if (.not.det_only) then 
        do k=1,LL
         do l=1,LL
          sum=Zero
          do kk=1,LL
           sum=sum+ainv(k,kk)*Green0(kk,l,i)
          enddo
          Gnew(k,l,i)=sum
         enddo
        enddo
       endif
      enddo
      end
C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: update_ini.f
C     TYPE   : subroutine
C     PURPOSE: adds the single particle correction to the non-interacting 
C              Green's function
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine update_ini(gnew)
      include 'lisaqmc.dat'
      integer job,i,j,k,l,kk     
      double precision a,ainv,w,adet,sum,gnew
      dimension a(MTAU,MTAU),ainv(MTAU,MTAU),
     &          w(NORB),adet(2)
      dimension gnew(MTAU,MTAU,NORB)

c     uses update scheme to renormalize the Green0 -
c     account for the exp(-dtau*U*N/2) term
c
      do i=1,Nint
       sum=Zero
       do j=1,i-1
        sum=sum+U(i,j)
       enddo
       do j=i+1,Nint
        sum=sum+U(j,i)
       enddo
       w(i)=(-shift(i)-sum/Two)*Beta/LL
      enddo
      job=11

      do i=1,Nint
       do l=1,LL
        do k=1,LL
         a(k,l)=-Green0(k,l,i)*(exp(w(i))-One)
        enddo
        a(l,l)=One-(Green0(l,l,i)-One)*(exp(w(i))-One)
       enddo
       call inverse(a,ainv,adet,job)
       do k=1,LL
        do l=1,LL
         sum=Zero
         do kk=1,LL
          sum=sum+ainv(k,kk)*Green0(kk,l,i)
         enddo
         Gnew(k,l,i)=sum
        enddo
       enddo
      enddo
      end

c======================================================================
C PROGRAM: inverse.f
C TYPE   : subroutine
C PURPOSE: calculate inverse of matrix
C========+=========+=========+=========+=========+=========+=========+=$
      Subroutine inverse(a,y,det,job)
      include 'lisaqmc.dat'
      integer nn,i,j,nsize,job,ipvt
      double precision a,y,z,det,rcond
      dimension a(MTAU,MTAU),y(MTAU,MTAU)
      dimension z(MTAU),ipvt(MTAU),det(2)
      nsize=MTAU
      nn=LL
      do i=1,nn
       do j=1,nn
        y(i,j)=a(i,j)
       enddo 
      enddo
      call dgeco(y,nsize,nn,ipvt,rcond,z)
      call dgedi(y,nsize,nn,ipvt,det,z,job)
      end
c======================================================================
C========+=========+=========+=========+=========+=========+=========+=$
C	PROGRAM: fill_index
C 	TYPE   : subroutine
C       PURPOSE: conversion for linear to 2D index
C========+=========+=========+=========+=========+=========+=========+=$
	subroutine fill_index
        include 'lisaqmc.dat'
	integer i,j,r
	r=0
	do i=1,Nint*Spinfac
	 do j=1,i-1
	  if (abs(U(i,j)).gt.1.D-5) then
	  r=r+1
	  ind(1,r)=i
	  ind(2,r)=j
	  endif
	 enddo
	enddo
	end
C========+=========+=========+=========+=========+=========+=========+=$
C    Seed generator for parallel running
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine seedgen(Idum,ncpu,seed)
        integer seed(0:5000)
        nran(i)=mod(int(i*ranw(Idum)),i) + 1

        seed(0)=-1
        seed(1)=nran(10000)
        do icpu=2,ncpu
 20      iseed=nran(10000)
         j=icpu
 10      if (iseed.lt.seed(j-1)) then
          j=j-1
          goto 10
         endif
         if (iseed.eq.seed(j-1)) goto 20
         do l=icpu-1,j,-1
          seed(l+1)=seed(l)
         enddo
         seed(j)=iseed
        enddo
        end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: dasum.f daxpy.f  ddot.f dgeco.f dgedi.f dgefa.f 
C                dscal.f dswap.f idamax.f
C       TYPE   : collection of subroutines 
C       PURPOSE: calculate inverse and determinant (look at 
C                subroutine dgedi.f) 
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
*DECK DASUM
      DOUBLE PRECISION FUNCTION DASUM (N, DX, INCX)
C***BEGIN PROLLOGUE  DASUM
C***PURPOSE  Compute the sum of the magnitudes of the elements of a
C            vector.
C***LLIBRARY   SLLATEC (BLLAS)
C***CATEGORY  D1A3A
C***TYPE      DOUBLE PRECISION (SASUM-S, DASUM-D, SCASUM-C)
C***KEYWORDS  BLLAS, LLINEAR ALLGEBRA, SUM OF MAGNITUDES OF A VECTOR
C***AUTHOR  LLawson, C. LL., (JPLL)
C           Hanson, R. J., (SNLLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPLL)
C***DESCRIPTION
C
C                B LL A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DASUM  double precision result (zero if N .LE. 0)
C
C     Returns sum of magnitudes of double precision DX.
C     DASUM = sum from 0 to N-1 of ABS(DX(IX+I*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. LL. LLawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  DASUM
      DOUBLE PRECISION DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  DASUM
      DASUM = 0.0D0
      IF (N .LE. 0) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DASUM = DASUM + ABS(DX(IX))
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 6.
C
   20 M = MOD(N,6)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DASUM = DASUM + ABS(DX(I))
   30 CONTINUE
      IF (N .LT. 6) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DASUM = DASUM + ABS(DX(I)) + ABS(DX(I+1)) + ABS(DX(I+2)) +
     1          ABS(DX(I+3)) + ABS(DX(I+4)) + ABS(DX(I+5))
   50 CONTINUE
      RETURN
      END
*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
C***BEGIN PROLLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LLIBRARY   SLLATEC (BLLAS)
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLLAS, LLINEAR ALLGEBRA, TRIAD, VECTOR
C***AUTHOR  LLawson, C. LL., (JPLL)
C           Hanson, R. J., (SNLLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPLL)
C***DESCRIPTION
C
C                B LL A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LLY+I*INCY) with DA*DX(LLX+I*INCX) +
C       DY(LLY+I*INCY),
C     where LLX = 1 if INCX .GE. 0, else LLX = 1+(1-N)*INCX, and LLY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. LL. LLawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LLX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
C***BEGIN PROLLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***LLIBRARY   SLLATEC (BLLAS)
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLLAS, INNER PRODUCT, LLINEAR ALLGEBRA, VECTOR
C***AUTHOR  LLawson, C. LL., (JPLL)
C           Hanson, R. J., (SNLLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPLL)
C***DESCRIPTION
C
C                B LL A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LLX+I*INCX) * DY(LLY+I*INCY),
C     where LLX = 1 if INCX .GE. 0, else LLX = 1+(1-N)*INCX, and LLY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. LL. LLawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LLX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DGECO
      SUBROUTINE DGECO (A, LLDA, N, IPVT, RCOND, Z)
C***BEGIN PROLLOGUE  DGECO
C***PURPOSE  Factor a matrix using Gaussian elimination and estimate
C            the condition number of the matrix.
C***LLIBRARY   SLLATEC (LLINPACK)
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGECO-S, DGECO-D, CGECO-C)
C***KEYWORDS  CONDITION NUMBER, GENERALL MATRIX, LLINEAR ALLGEBRA, LLINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGECO factors a double precision matrix by Gaussian elimination
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, DGEFA is slightly faster.
C     To solve  A*X = B , follow DGECO by DGESLL.
C     To compute  INVERSE(A)*C , follow DGECO by DGESLL.
C     To compute  DETERMINANT(A) , follow DGECO by DGEDI.
C     To compute  INVERSE(A) , follow DGECO by DGEDI.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LLDA, N)
C                the matrix to be factored.
C
C        LLDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = LL*U  where
C                LL  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an INTEGER vector of pivot indices.
C
C        RCOND   DOUBLE PRECISION
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILLON  may cause
C                relative perturbations in  X  of size  EPSILLON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       DOUBLE PRECISION(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LLINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALED  DASUM, DAXPY, DDOT, DGEFA, DSCALL
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  DGECO
      INTEGER LLDA,N,IPVT(*)
      DOUBLE PRECISION A(LLDA,*),Z(*)
      DOUBLE PRECISION RCOND
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,LL
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  DGECO
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = MAX(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LLDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LLOCALL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLLY RESCALED TO AVOID
C     OVERFLLOW.
C
C     SOLLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 30
            S = ABS(A(K,K))/ABS(EK-Z(K))
            CALL DSCALL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCALL(N,S,Z,1)
C
C     SOLLVE TRANS(LL)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/ABS(Z(K))
            CALL DSCALL(N,S,Z,1)
  110    CONTINUE
         LL = IPVT(K)
         T = Z(LL)
         Z(LL) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCALL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLLVE LL*V = Y
C
      DO 140 K = 1, N
         LL = IPVT(K)
         T = Z(LL)
         Z(LL) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/ABS(Z(K))
            CALL DSCALL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCALL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 150
            S = ABS(A(K,K))/ABS(Z(K))
            CALL DSCALL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCALL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
*DECK DGEDI
      SUBROUTINE DGEDI (A, LLDA, N, IPVT, DET, WORK, JOB)
C***BEGIN PROLLOGUE  DGEDI
C***PURPOSE  Compute the determinant and inverse of a matrix using the
C            factors computed by DGECO or DGEFA.
C***LLIBRARY   SLLATEC (LLINPACK)
C***CATEGORY  D3A1, D2A1
C***TYPE      DOUBLE PRECISION (SGEDI-S, DGEDI-D, CGEDI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LLINEAR ALLGEBRA, LLINPACK, MATRIX
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEDI computes the determinant and inverse of a matrix
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LLDA, N)
C                the output from DGECO or DGEFA.
C
C        LLDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       inverse of original matrix if requested.
C                Otherwise unchanged.
C
C        DET     DOUBLE PRECISION(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if DGECO has set RCOND .GT. 0.0 or DGEFA has set
C        INFO .EQ. 0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LLINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALED  DAXPY, DSCALL, DSWAP
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  DGEDI
      INTEGER LLDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LLDA,*),DET(2),WORK(*)
C
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,LL,NM1
C***FIRST EXECUTABLE STATEMENT  DGEDI
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCALL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(LL)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            LL = IPVT(K)
            IF (LL .NE. K) CALL DSWAP(N,A(1,K),1,A(1,LL),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
*DECK DGEFA
      SUBROUTINE DGEFA (A, LLDA, N, IPVT, INFO)
C***BEGIN PROLLOGUE  DGEFA
C***PURPOSE  Factor a matrix using Gaussian elimination.
C***LLIBRARY   SLLATEC (LLINPACK)
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
C***KEYWORDS  GENERALL MATRIX, LLINEAR ALLGEBRA, LLINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LLDA, N)
C                the matrix to be factored.
C
C        LLDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = LL*U  where
C                LL  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESLL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LLINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALED  DAXPY, DSCALL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  DGEFA
      INTEGER LLDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LLDA,*)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,LL,NM1
C
C     GAUSSIAN ELLIMINATION WITH PARTIALL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND LL = PIVOT INDEX
C
         LL = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = LL
C
C        ZERO PIVOT IMPLLIES THIS COLLUMN ALLREADY TRIANGULLARIZED
C
         IF (A(LL,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (LL .EQ. K) GO TO 10
               T = A(LL,K)
               A(LL,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCALL(N-K,T,A(K+1,K),1)
C
C           ROW ELLIMINATION WITH COLLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(LL,J)
               IF (LL .EQ. K) GO TO 20
                  A(LL,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DSCALL
      SUBROUTINE DSCALL (N, DA, DX, INCX)
C***BEGIN PROLLOGUE  DSCALL
C***PURPOSE  Multiply a vector by a constant.
C***LLIBRARY   SLLATEC (BLLAS)
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCALL-S, DSCALL-D, CSCALL-C)
C***KEYWORDS  BLLAS, LLINEAR ALLGEBRA, SCALE, VECTOR
C***AUTHOR  LLawson, C. LL., (JPLL)
C           Hanson, R. J., (SNLLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPLL)
C***DESCRIPTION
C
C                B LL A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. LL. LLawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  DSCALL
      DOUBLE PRECISION DA, DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  DSCALL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
*DECK DSWAP
      SUBROUTINE DSWAP (N, DX, INCX, DY, INCY)
C***BEGIN PROLLOGUE  DSWAP
C***PURPOSE  Interchange two vectors.
C***LLIBRARY   SLLATEC (BLLAS)
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
C***KEYWORDS  BLLAS, INTERCHANGE, LLINEAR ALLGEBRA, VECTOR
C***AUTHOR  LLawson, C. LL., (JPLL)
C           Hanson, R. J., (SNLLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPLL)
C***DESCRIPTION
C
C                B LL A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DX  input vector DY (unchanged if N .LE. 0)
C       DY  input vector DX (unchanged if N .LE. 0)
C
C     Interchange double precision DX and double precision DY.
C     For I = 0 to N-1, interchange  DX(LLX+I*INCX) and DY(LLY+I*INCY),
C     where LLX = 1 if INCX .GE. 0, else LLX = 1+(1-N)*INCX, and LLY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. LL. LLawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LLX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  DSWAP
      DOUBLE PRECISION DX(*), DY(*), DTEMP1, DTEMP2, DTEMP3
C***FIRST EXECUTABLE STATEMENT  DSWAP
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 3.
C
   20 M = MOD(N,3)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF (N .LT. 3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70 CONTINUE
      RETURN
      END
*DECK IDAMAX
      INTEGER FUNCTION IDAMAX (N, DX, INCX)
C***BEGIN PROLLOGUE  IDAMAX
C***PURPOSE  Find the smallest index of that component of a vector
C            having the maximum magnitude.
C***LLIBRARY   SLLATEC (BLLAS)
C***CATEGORY  D1A2
C***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLLAS, LLINEAR ALLGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  LLawson, C. LL., (JPLL)
C           Hanson, R. J., (SNLLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPLL)
C***DESCRIPTION
C
C                B LL A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. LL. LLawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLLOGUE  IDAMAX
      DOUBLE PRECISION DX(*), DMAX, XMAG
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF (N .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increments not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = ABS(DX(IX))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increments equal to 1.
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
        XMAG = ABS(DX(I))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END
