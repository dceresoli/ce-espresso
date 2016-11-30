! ----------  ----------  ----------  ----------  ----------  ----------
! \__ \__  \__ \__ \__ \__ \__ \__ \__ \__ \__ \__ \__ \__ \__ \__ \__
      SUBROUTINE nfourier2(disc, xmom2, xmom3, RIndata)
      include 'lisaqmc.dat'
!
       double precision rindata(mtau), bt,
     &        a(mtau),b(mtau),c(mtau),d(mtau), xmw(0:mtau), xm(mtau+1),
     &        w(0:mtau),q(0:mtau),r(0:mtau),uu(0:mtau),
     &        du(0:mtau),dn(0:mtau),p(0:mtau)
       double precision xpi,delta,disc,xmom2,xmom3,fdsum,sdsum
       double precision yp0,ypn,ypn2,om,ww
       integer i,j,k,iw
       COMPLEX*16   coutdata(0:momega),cdummy,explus,ex,csum,ci
!
!       PRINT *, 'from Fourier new disc, xmom2,3 = ', disc, xmom2, xmom3
!
       xpi=acos(-1.d0)
       ci=(0.d0,1.d0)
       Delta = Beta/LL
       DO I = 1,LL
          w(I-1) = rindata(I)
       ENDDO
          w(LL) = disc-rindata(1)
                                                                                        
! ----------  ----------  ----------  ----------  ----------  ----------
!   Spline interpolation:
!   G(tau) = a(i) + b(i) (tau-tau_i) + c(i) ( )^2 + d(i) ( )^3
!   For given sum of the first and second derivatives for the END points
!   Written by V. Oudovenko
! ----------  ----------  ----------  ----------  ----------  ----------
                                                                                        
      FDsum = xmom2   ! Sum of the first  derivatives at the END points
      SDsum = xmom3   ! Sum of the second derivatives at the END points
!
        yp0 = w(1)-w(0)
        ypn = w(LL)-w(LL-1)
       ypn2 = w(LL)+w(LL-2)-2.*w(LL-1)
!
       q(0) = 4.D0
       q(1) = 2.D0-0.5D0/q(0)
       r(0) = 0.5D0/q(0)
       p(0) = 2.D0
       p(1) = 2.D0-0.5D0/q(0)
       uu(0) = (6.d0/delta)*( (yp0+ypn)/delta -FDsum)+2.*SDsum
      du(0) = uu(0)
      dn(1) = 3.D0*ypn2/delta**2-0.5D0*SDsum+0.5D0*uu(0)/q(0)
!
      DO k = 1,LL-3
           q(k+1) = ( 2.D0-0.25D0/q(k) )
             r(k) = 0.5D0*r(k-1)/q(k)
             uu(k) = 3.D0/delta**2*( W(k+1)+W(k-1)-2.*W(k) )
            du(k) = uu(k)-0.5D0*du(k-1)/q(k-1)
           p(k+1) = p(k)-(2.d0*r(k)*r(k-1))
          dn(k+1) = dn(k)+(-1)**(k)*du(k)*2.D0*r(k)
      END DO
           uu(LL-2) = 3.D0/delta**2*( W(LL-2+1)+W(LL-2-1)-2.*W(LL-2) )
          du(LL-2) = uu(LL-2)-0.5D0*du(LL-3)/q(LL-3)
           r(LL-2) = 0.5D0*r(LL-3)/q(LL-2)
!
         xmw(LL-1) = dn(LL-2)-du(LL-2)*(0.5D0-r(LL-3))/q(LL-2)
         xmw(LL-1) = xmw(LL-1)/( p(LL-2)-(0.5D0-r(LL-3))**2/q(LL-2) )
         xmw(LL-2) = ( du(LL-2)-(0.5D0-r(LL-3))*xmw(LL-1) )/ q(LL-2)
!
      DO k = LL-3,1,-1
         XMw(k) = (du(k)+((-1)**k)*r(k-1)*XMw(LL-1)-0.5D0*XMw(k+1))/q(k)
      END DO
         XMw(0) = ( uu(0)+XMw(LL-1)-XMw(1) )/ 4.D0
!
      XMw(LL) = SDsum-XMw(0)
      DO k = 1,LL+1
          XM(k) = XMw(k-1)
!          PRINT *, k, XM(k)
!
      END DO
! ----------  ----------  ----------  ----------  ----------  ----------
! The following formulas are taken from Stoer and Bulirsch p. 98
! ----------  ----------  ----------  ----------  ----------  ----------
                                                                                        
       DO j = 1, LL
          a(j) = W(j-1)
          c(j) = XM(j)/2.D0
          b(j) = (W(j)-W(j-1))/delta - (2.D0*XM(j)+XM(j+1))*delta/6.
          d(j) = ( XM(j+1)-XM(j) )/(6.*delta)
       END DO
                                                                                        
! ----------  ----------  ----------  ----------  ----------  ----------
!c The spline multiplied by the exponential can be exlicitely integrated.
!c The following formulas were obtained using MATHEMATICA
! ----------  ----------  ----------  ----------  ----------  ----------
                                                                                        
        DO i = 0,Im
           om = (2.D0*(i)+1.D0)*xpi/Beta
           coutdata(i) = 0.D0
           DO j = 1,LL
                                                                                        
              cdummy = ci*om*delta*j
              explus = exp(cdummy)
                                                                                        
              cdummy = ci*om*delta*(j-1)
              ex = exp(cdummy)
                                                                                        
              coutdata(i) = coutdata(i) + explus*(
     &         ( -6.D0* d(j) )/om**4 +
     &         ( 2.D0*ci*c(j) + 6.D0*delta*ci*d(j)  )/om**3 +
     &         ( b(j)+ 2.D0*delta*c(j)+ 3.D0*delta**2*d(j) )/om**2 +
     &         (- ci*a(j) - delta*ci*b(j) - delta**2*ci*c(j) -
     &         delta**3*ci*d(j))/om)
!
              coutdata(i) = coutdata(i) + ex*(
     &        6.D0*d(j)/om**4 - 2.D0*ci*c(j)/om**3
     &        -b(j)/om**2 + ci*a(j)/om)
           ENDDO
        END DO

        do iw=0,Im
         ww=(Two*iw+One)*xpi/Beta
         write(8,*)ww,dreal(coutdata(iw)),imag(coutdata(iw))
        enddo

!
        RETURN
        END
                                                                                        

