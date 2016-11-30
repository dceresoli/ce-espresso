        subroutine symmetry(g,xnn,logbin,fname)
        include 'lisaqmc.dat'

        character*80 linemc
        character*80 fname,fn
        logical logbin
        double precision g,gbin,xnn,sum
        integer i,j,ic,ione,itwo,l,n, ios
        dimension g(NORB,0:mtau-1),gbin(NORB,0:mtau-1)
        dimension xnn(Norb,Norb)
        integer n_class_one,n_class_two
        integer eq_one,eq_two
        dimension n_class_one(Norb),eq_one(Norb,Norb)
        dimension n_class_two(Norb*Norb),
     &        eq_two(2,Norb*Norb,Norb*Norb)


	open(3000,file='dmft.in',status='old',form='formatted')

 10    read(3000,201,end=100)linemc
!         write(*,201)linemc
!         write(*,*)'spinfac = ',spinfac
 201     format(a80)
       if(index(linemc(1:80),'Equivalence').ne.0) then
       write(*,201)linemc
       read(3000,*)ione
        do i=1,ione
         read(3000,202)  n_class_one(i),
     &(eq_one(j,i),j=1,n_class_one(i))
         write(6,202)n_class_one(i),(eq_one(j,i),j=1,n_class_one(i))
        enddo
 202     format(i3,2x,10i3)
          read(3000,*)itwo
           do i=1,itwo
            read(3000,203)n_class_two(i),(eq_two(1,j,i),eq_two(2,j,i),
     &                j=1,n_class_two(i))
            write(6,203)n_class_two(i),(eq_two(1,j,i),eq_two(2,j,i),
     &                j=1,n_class_two(i))
           enddo
 203     format(i3,2x,20(2i2,2x))
          else
           goto 10
          endif
          if (spinfac.eq.2) then
           do ic=1,ione
            do j=1,n_class_one(ic)
             eq_one(n_class_one(ic)+j,ic)=eq_one(j,ic)+Nint/2
            enddo
            n_class_one(ic)=2*n_class_one(ic)
           enddo
           do ic=1,itwo
            do j=1,n_class_two(ic)
             eq_two(1,j+n_class_two(ic),ic)=eq_two(1,j,ic)+Nint/2
             eq_two(2,j+n_class_two(ic),ic)=eq_two(2,j,ic)+Nint/2 
             eq_two(1,j,ic+itwo)=eq_two(1,j,ic)+Nint/2
             eq_two(2,j,ic+itwo)=eq_two(2,j,ic)
             eq_two(1,j+n_class_two(ic),ic+itwo)=eq_two(1,j,ic)
             eq_two(2,j+n_class_two(ic),ic+itwo)=eq_two(2,j,ic)+Nint/2
            enddo
            n_class_two(ic)=2*n_class_two(ic)
            n_class_two(ic+itwo)=n_class_two(ic)   
           enddo
           itwo=2*itwo
          endif
          write(6,*)'test sym',spinfac
          do ic=1,ione
           write(6,*)n_class_one(ic),':  ',
     &      (eq_one(j,ic),j=1,n_class_one(ic))
          enddo
          do ic=1,itwo
            write(6,*)n_class_two(ic),':'
            write(6,*)(eq_two(1,j,ic),eq_two(2,j,ic),' ',
     &      j=1,n_class_two(ic))
          enddo
          

          do ic=1,ione
           do l=0,LL-1
            sum=0.d0
            n=n_class_one(ic)
            do i=1,n
             sum=sum+g(eq_one(i,ic),l)
            enddo
            sum=sum/n
            do i=1,n
             g(eq_one(i,ic),l)=sum
            enddo
           enddo
          enddo

          do ic=1,itwo
           sum=0.d0
           n=n_class_two(ic)
           do i=1,n
            sum=sum+xnn(eq_two(1,i,ic),eq_two(2,i,ic))
           enddo
           sum=sum/n
           do i=1,n
            xnn(eq_two(1,i,ic),eq_two(2,i,ic))=sum
           enddo
          enddo

          if (logbin) then
          write(6,*)'GBIN SYMMETRIZED',ione
          do ic=1,ione
           call mknam(fn,fname,ic)
           write(*,*)'mknam',fn,fname,ic
           open(60+ic,file=fn,form='formatted',status='unknown')
          enddo
           rewind(50)
 20        read(50,201,end=110)linemc
           do ic=1,ione
            write(60+ic,201)linemc
           enddo
           do i=1,Nint
            read(50,567)(gbin(i,l),l=0,LL-1)
           enddo
            
           do ic=1,ione
            do l=0,LL-1
             sum=0.d0
             n=n_class_one(ic)
             do i=1,n
              sum=sum+gbin(eq_one(i,ic),l)
             enddo
             sum=sum/n
             do i=1,n
              gbin(eq_one(i,ic),l)=sum
             enddo
            enddo
           enddo
           do ic=1,ione
            write(60+ic,567)(gbin(eq_one(1,ic),l),l=0,LL-1)
           enddo
           goto 20
          endif
 567       format(200f12.8)

           write(*,*)'symmetrization finished'
           return
  100      write(*,*)'no symmetrization'       

	close(3000)

  110      return  
           
           end

      SUBROUTINE MKNAM(FNAME,OLDNAM,ILOOP)

! create a filename with running index ILOOP from
! a given 'parental' file name OLDNAM

      CHARACTER*80 FNAME,OLDNAM
      CHARACTER*4  ALOOP

      WRITE(ALOOP,'(I4)')ILOOP

      FNAME=OLDNAM

      DO I=LEN(ALOOP),1,-1
         IF(ALOOP(I:I).NE.' ') IFROM=I
      ENDDO

      ITO=LEN(ALOOP)

      DO 7777 I=LEN(FNAME),1,-1
         IF(FNAME(I:I).NE.' ') THEN
            FNAME(I+1:LEN(FNAME))='_' // ALOOP(IFROM:ITO)
            GOTO 7778
         ENDIF
 7777 CONTINUE
 7778 CONTINUE

      RETURN
      END

