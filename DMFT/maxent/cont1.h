c     COMMON BLOCKS FOR fcont.f
c-----------------------------------------------------------------------

      IMPLICIT REAL*8(A-H,O-Z)
c     IMPLICIT NONE

      INTEGER tmax,wmax,qmax
      PARAMETER (tmax=257)
      PARAMETER (wmax=500)
      PARAMETER (qmax=20)

      INTEGER nt
      INTEGER nw
      INTEGER seed
      REAL*8    beta
      REAL*8    t(tmax)
      REAL*8    w(-wmax:wmax)
      REAL*8    axt(tmax)
      REAL*8    sxt(tmax)
      REAL*8    txt(tmax)
      REAL*8    ker(-wmax:wmax,tmax)
      REAL*8    f0(-wmax:wmax)
      REAL*8    f1(-wmax:wmax)
      REAL*8    f2(-wmax:wmax)
      COMMON/TWBLOCK/beta,t,w,axt,sxt,txt,ker,f0,f1,f2,nt,nw,seed

c-----------------------------------------------------------------------
