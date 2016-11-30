      SUBROUTINE MAXENT_CALC(nq,nt1,idg,deltag,beta_in,
     .nw_in,dw,steps,alpha0,nr,seed_in)
      INCLUDE 'cont1.h'

      INTEGER i,iq,nq,nr,iw,st,iqx,iqy,iqz,steps,nw_in, seed_in
      REAL*8    a(-wmax:wmax),xt1(tmax),x2(qmax),ss(3,qmax)
      REAL*8    as(-wmax:wmax,qmax)
      REAL*8    s0,s1,s2,alpha0,dw,CHI2,RAN
      REAL*8    dummy, beta_in
      CHARACTER*8 string

      beta = beta_in
      nw = nw_in
      nt = nt1
      
      WRITE(*,*)'Number of momenta = ',nq
      WRITE(*,*)'Number of times = ',nt
      WRITE(*,*)'idg, delta-GEstimate Delta-G = ',idg,deltag    
      WRITE(*,*)'Beta = ',beta
      WRITE(*,*)'Number of frequencies = ',nw
      WRITE(*,*)'Delta frequency = ',dw
      WRITE(*,*)'Number of annealing steps = ',steps
      WRITE(*,*)'Starting Alpha = ',alpha0
      WRITE(*,*)'Number of smoothed runs = ',nr
      seed=53718551
      seed=seed_in
      WRITE(*,*)'Random number seed = ',seed
     
      CALL READIN(nq)

      DO 10 iw=-nw,nw
         w(iw)=iw*dw
 10   CONTINUE

      DO 50 iq=1,nq
         WRITE(*,*)'Starting Maxent ',iq         
         CALL GETDATA(iq,iqx,iqy,iqz)
         CALL CONTINUE(steps,alpha0,nr,a)
         OPEN(UNIT=10,FILE='dos',STATUS='unknown')
cc         DO 20 i=1,100000
c           READ(10,*,END=30)string
c20      CONTINUE
c30      CONTINUE
c30      WRITE(10,*)'#',iqx,iqy,iqz
         DO 40 iw=-nw,nw
            as(iw,iq)=a(iw)
 40      CONTINUE
         CALL CALCXT(a,xt1)
         x2(iq)=CHI2(xt1)
         CALL SUMRULES(a,s0,s1,s2)
         ss(1,iq)=s0
         ss(2,iq)=s1
         ss(3,iq)=s2
 50   CONTINUE
         DO 44 iw=-nw,nw
            WRITE(10,80)w(iw),(as(iw,iq),iq=1,nq)
 44      CONTINUE
 80   FORMAT(6F15.8)

      OPEN(UNIT=20,FILE='sums.dat',STATUS='unknown')
      DO 200 i=1,nq
         WRITE(20,310)i,ss(1,i),ss(2,i),ss(3,i)
 200  CONTINUE
      CLOSE(10)
 310  FORMAT(I4,'   Sums:  ',3F9.5)

      END SUBROUTINE MAXENT_CALC

c------------------------------------------------------------------------------
      SUBROUTINE CONTINUE(steps,alpha0,nr,a)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER i,iq,nr,st,steps
      REAL*8    a(-wmax:wmax),mm(-wmax:wmax),xt1(tmax)
      REAL*8    aa(-wmax:wmax,qmax),x2(qmax),ss(3,qmax)
      REAL*8    n0,rfac,temp,s0,s1,s2,alpha,alpha0,ent,tr,ls,adiv
      REAL*8    CHI2,ENTROPY,RAN

      n0=1.
      CALL INITKER
      CALL INITF0
      CALL INITDLDA
      CALL FLATMODEL(mm)
      CALL NORMALIZE(n0,mm)
      DO 10 i=-nw,nw
         a(i)=RAN(seed)
 10   CONTINUE
      CALL NORMALIZE(n0,a)
      temp=10.
      rfac=1.

      alpha=alpha0
 20   CALL MAXENT(mm,a,rfac,temp,alpha,steps)
      ent=ENTROPY(a,mm)
      CALL LAMBDA(a,alpha,tr)
      ls=-2.*ent*alpha
      WRITE(*,*)'-2*alpha*S = ',ls
      WRITE(*,*)'Trace      = ',tr
      temp=0.001
      rfac=0.05
      IF (ABS(ls/tr-1.).LT.0.05) THEN
          GOTO 30
       ELSE
          IF (tr/ls.LT.0.05) THEN
             alpha=alpha*0.05
          ELSE
            alpha=alpha*(tr/ls)*(1.+0.001*(RAN(seed)-0.5))
         ENDIF
         GOTO 20
      ENDIF

 30   DO 100 st=1,nr
         WRITE(*,*)'Starting smoothing ',st         
         CALL SMOOTH(a,3)
         CALL NORMALIZE(n0,a)
         temp=0.005
         rfac=0.005
         CALL MAXENT(mm,a,rfac,temp,alpha,steps)
 100  CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE SMOOTH(a,ns)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER i,j,i1,i2,ii,ns
      REAL*8    a(-wmax:wmax),a1(-wmax:wmax)
      
      DO 20 i=-nw,nw
         i1=i-ns
         IF (i1.LT.-nw) i1=-nw
         i2=i+ns
         IF (i2.GT.nw) i2=nw
         a1(i)=0.
         ii=0
         DO 10 j=i1,i2
            ii=ii+1
            a1(i)=a1(i)+a(j)
 10      CONTINUE
         a1(i)=a1(i)/FLOAT(ii)
 20   CONTINUE
      DO 30 i=-nw,nw
         a(i)=a1(i)
 30   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE READIN(nq)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER xyz(3,qmax)
      REAL*8  axt1(tmax,qmax),sxt1(tmax,qmax)
      COMMON/BQDATA/axt1,sxt1,xyz

      INTEGER it,iit,iq,nq,nt1

      OPEN(UNIT=10,FILE='gtau.out',STATUS='old')
      READ(10,*)
      DO 10 it=1,nt
c        READ(10,*)t(it),(axt1(it,iq),iq=1,nq),(sxt1(it,iq),iq=1,nq)
        READ(10,*)t(it),(axt1(it,iq),iq=1,nq)
        DO 20 iq=1,nq
          xyz(1,iq)=1.
          xyz(2,iq)=1.
          xyz(3,iq)=1.
c          Print*,t(it),axt1(it,iq)
          if( idg.eq.0 )then
c            sxt1(it,iq) = sxt1(it,iq) * deltag 
            sxt1(it,iq) = 0.d0
           else
            sxt1(it,iq) = deltag 
          end if
          axt1(it,iq)=-axt1(it,iq)
c          IF (ABS(sxt1(it,iq)).LT.0.00001) sxt1(it,iq)=0.00001
          IF (ABS(sxt1(it,iq)).LT.0.00001) sxt1(it,iq)=0.00001
          sxt1(it,iq)=1./sxt1(it,iq)**2
 20     CONTINUE
 10   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE GETDATA(iq,iqx,iqy,iqz)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER xyz(3,qmax)
      REAL*8    axt1(tmax,qmax),sxt1(tmax,qmax)
      COMMON/BQDATA/axt1,sxt1,xyz

      INTEGER i,iq,iqx,iqy,iqz
      
      iqx=xyz(1,iq)
      iqy=xyz(2,iq)
      iqz=xyz(3,iq)
      DO 10 i=1,nt
         axt(i)=axt1(i,iq)
         sxt(i)=sxt1(i,iq)
 10   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE MAXENT(mm,a,rfac,temp,alpha,steps)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER i,j,j1,j2,jj,it,steps,acc,try
      REAL*8    a(-wmax:wmax),da(-wmax:wmax),xt1(tmax),xt2(tmax)
      REAL*8    dj1,dj2,x1,x2,temp,alpha,wgt,eps,RAN
      REAL*8    CHI2,aa1,aa2,de,p,mm1,mm2,mm(-wmax:wmax),rfac,arat,am

      am=a(-nw)
      DO 10 i=-nw,nw         
         da(i)=rfac*(a(i)/10.+0.01)
         IF (a(i).GT.am) am=a(i)
 10   CONTINUE
      eps=1.E-12

      print*, ran(seed)
      DO 200 i=1,steps
         CALL CALCXT(a,xt1)
         x1=CHI2(xt1)
         DO 100 jj=1,2*nw+1

 35         j1=MIN(INT(RAN(seed)*(2*nw+1)),2*nw)-nw
      print*, ran(seed)
 40         dj1=da(j1)*(RAN(seed)-0.5)
            IF (a(j1)+dj1.LT.0.) GOTO 40
 50         j2=MIN(INT(RAN(seed)*(2*nw+1)),2*nw)-nw
            IF (j1.EQ.j2) GOTO 50
            dj2=-dj1+0.05*da(j2)*(RAN(seed)-0.5)            
            IF (a(j2)+dj2.LT.0.) GOTO 35
            IF (a(j1).GT.0.1*am) try=try+1

            DO 70 it=1,nt
               xt2(it)=xt1(it)+dj1*ker(j1,it)+dj2*ker(j2,it)
 70         CONTINUE
            x2=CHI2(xt2)
            aa1=a(j1)+dj1
            aa2=a(j2)+dj2
            mm1=mm(j1)
            mm2=mm(j2)
            de=0.
            IF (aa1.GT.eps) de=de-f0(j1)*aa1*LOG(aa1/mm1)
            IF (a(j1).GT.eps) de=de+f0(j1)*a(j1)*LOG(a(j1)/mm1)
            IF (aa2.GT.eps) de=de-f0(j2)*aa2*LOG(aa2/mm2)
            IF (a(j2).GT.eps) de=de+f0(j2)*a(j2)*LOG(a(j2)/mm2)            
            wgt=((x1-x2)/2.+alpha*de)/temp
            IF (wgt.GT.0.d0) THEN
               p=1.d0
c??         ELSEIF (wgt.LT.10.d0) THEN               
            ELSEIF (wgt.GT.-100.d0) THEN               
               p=EXP(wgt)
            ELSE
               p=0.
            ENDIF
            IF (RAN(seed).LT.p) THEN
               IF (a(j1).GT.0.1*am) acc=acc+1
               a(j1)=aa1
               a(j2)=aa2
               DO 80 it=1,nt
                  xt1(it)=xt2(it)
 80            CONTINUE
               x1=x2
             ENDIF

 100      CONTINUE

          IF (MOD(i,100).EQ.0) THEN
             arat=FLOAT(acc)/FLOAT(try)
             IF (arat.GT.0.1) THEN
                IF (rfac.LT.0.01) rfac=rfac*1.5
             ELSE
                IF (rfac.GT.0.001) rfac=rfac/1.5
             ENDIF
             am=a(1)
             DO 110 j=-nw,nw
                da(j)=rfac*(0.1*a(j)+0.01*am)
                IF (a(j).GT.am) am=a(j)

 110         CONTINUE
             IF (MOD(i,steps/10).EQ.0) THEN
                WRITE(*,*)i,x1,acc,rfac
             ENDIF
             temp=temp/1.5
             acc=0
             try=0
          ENDIF
 200   CONTINUE
 300   FORMAT(I8,F12.6,I8,F12.6)
       RETURN
       END

c------------------------------------------------------------------------------
      SUBROUTINE NORMALIZE(n0,a)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER iw
      REAL*8    a(-wmax:wmax),n0,n1,f

      n1=0.
      DO 10 iw=-nw,nw
         n1=n1+a(iw)*f0(iw)
 10   CONTINUE
      f=n0/n1
      DO 20 iw=-nw,nw
         a(iw)=f*a(iw)
 20   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE CALCXT(a,xt)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER it,iw
      REAL*8    a(-wmax:wmax),xt(tmax)

      DO 20 it=1,nt
         xt(it)=0.
         DO 10 iw=-nw,nw
            xt(it)=xt(it)+a(iw)*ker(iw,it)
 10      CONTINUE
 20   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      REAL*8 FUNCTION CHI2(xt)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER i
      REAL*8    x2,xt(tmax)

      x2=0.
      DO 10 i=1,nt
          x2=x2+sxt(i)*(axt(i)-xt(i))**2
 10   CONTINUE
      CHI2=x2
      RETURN
      END

c------------------------------------------------------------------------------
      REAL*8 FUNCTION ENTROPY(a,mm)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER i
      REAL*8    a(-wmax:wmax),mm(-wmax:wmax),ent,eps

      ent=0.
      eps=1.E-12
      DO 10 i=-nw,nw
         IF (a(i).GT.eps.AND.mm(i).GT.eps) THEN
            ent=ent-a(i)*LOG(a(i)/mm(i))*f0(i)
         ENDIF
 10   CONTINUE
      ENTROPY=ent
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE LAMBDA(a,alpha,tr)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      REAL*8 dlda(2*wmax+1,2*wmax+1)
      COMMON/BLAM/dlda

      INTEGER i,j,iw,jw,indx(2*wmax+1)
      REAL*8    lam(2*wmax+1,2*wmax+1),alam(2*wmax+1,2*wmax+1)
      REAL*8    alami(2*wmax+1,2*wmax+1),a(-wmax:wmax)
      REAL*8    alpha,tr,dw

      dw=w(2)-w(1)
      DO 20 j=1,2*nw+1
         DO 10 i=1,2*nw+1
            iw=i-nw-1
            jw=j-nw-1
            lam(i,j)=SQRT(a(iw))*dlda(i,j)*SQRT(a(jw))*dw
            alam(i,j)=lam(i,j)
 10      CONTINUE
         alam(j,j)=lam(j,j)+alpha
 20   CONTINUE
      CALL INVERT(alam,alami,indx,2*nw+1,2*wmax+1)      
      tr=0.
      DO 40 i=1,2*nw+1
         DO 30 j=1,2*nw+1
            tr=tr+lam(i,j)*alami(j,i)
 30      CONTINUE
 40   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE INITKER
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER it,iw
      REAL*8    bw,tw,dw

      dw=w(2)-w(1)
      DO 20 iw=-nw,nw
         DO 10 it=1,nt
            tw=EXP(-t(it)*w(iw))
            bw=1.+EXP(-beta*w(iw))
            ker(iw,it)=dw*tw/bw
 10      CONTINUE
 20   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE INITDLDA
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      REAL*8 dlda(2*wmax+1,2*wmax+1)
      COMMON/BLAM/dlda

      INTEGER i,j,iw,jw,it
      REAL*8    dw
      
      dw=w(2)-w(1)
      DO 20 j=1,2*nw+1
      DO 20 i=1,2*nw+1
         iw=i-nw-1
         jw=j-nw-1
         dlda(i,j)=0.
         DO 10 it=1,nt
            dlda(i,j)=dlda(i,j)+ker(iw,it)*sxt(it)*ker(jw,it)/dw**2       
 10      CONTINUE
 20   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE INITF0
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER iw
      REAL*8    dw

      dw=w(2)-w(1)
      DO 10 iw=-nw,nw
         f0(iw)=dw
         f1(iw)=dw*w(iw)
         f2(iw)=dw*w(iw)**2
 10   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE FLATMODEL(mm)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER i
      REAL*8    mm(-wmax:wmax)

      DO 10 i=-nw,nw
         mm(i)=1.
 10   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE SUMRULES(a,s0,s1,s2)
c------------------------------------------------------------------------------

      INCLUDE 'cont1.h'

      INTEGER i
      REAL*8    a(-wmax:wmax),s0,s1,s2

      s0=0.
      s1=0.
      s2=0.
      DO 10 i=-nw,nw
         s0=s0+a(i)*f0(i)
         s1=s1+a(i)*f1(i)
         s2=s2+a(i)*f2(i)
 10   CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE INVERT(a,y,indx,n,np)
c------------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 a(np,np),y(np,np)
      DIMENSION indx(np)
      DO 12 i=1,n
         DO 11 j=1,n
            y(i,j)=0.
 11      CONTINUE
         y(i,i)=1.
 12   CONTINUE
      CALL LUDCMP(a,n,np,indx,d)
      DO 13 j=1,n
         CALL LUBKSB(a,n,np,indx,y(1,j))
 13   CONTINUE
      RETURN
      END
c------------------------------------------------------------------------------
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
c------------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=2*500+1,TINY=1.0E-20)
      REAL*8 A(NP,NP),VV(NMAX)
      DIMENSION INDX(N)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) STOP 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
c------------------------------------------------------------------------------
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
c------------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(NP,NP),B(N)
      DIMENSION INDX(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END

c------------------------------------------------------------------------------

      REAL*8 FUNCTION RAN(IDUM)
      use random_numbers, only: randy
      IMPLICIT REAL*8(A-H,O-Z)
      ran=randy(idum)
      RETURN
      END

c       USE THESE COMMENTED OUT LINES IF REAL*8 DESIRED.
      REAL*8 FUNCTION RAN_old(IDUM)
      IMPLICIT REAL*8(A-H,O-Z)
c      FUNCTION RAN(IDUM)
      save
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
      DIMENSION IR(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1)STOP 'J.GT.97.OR.J.LT.1'
      IY=IR(J)
      RAN=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END





