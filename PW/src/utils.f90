!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! Author: Ivan Carnimeo
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! General-purpose routines for scalar products, printing,
! linear-algebra operators for exact-exchange and localization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matcalc (label, DoE, PrtMat, ninner, n, m, U, V, mat, ee)
  !
  USE kinds,    ONLY : dp
  USE io_global,ONLY : stdout
  USE becmod,   ONLY : calbec
  USE wvfct,    ONLY : current_k, wg
  IMPLICIT NONE
  !
  ! compute the (n,n) matrix representation <U|V>
  ! and energy from V (m,n) and U(m,n)
  !
  CHARACTER(len=*), INTENT(IN) :: label
  LOGICAL, INTENT(IN) :: DoE
  INTEGER, INTENT(IN) :: PrtMat, ninner, n, m
  COMPLEX(DP), INTENT(IN) :: U(ninner,n), V(ninner,m)
  REAL(DP), INTENT(OUT) :: mat(n,m), ee
  !
  INTEGER :: i
  CHARACTER(len=2) :: string

  CALL start_clock('matcalc')

  string = 'M-'
  mat = 0.0_dp
  CALL calbec(ninner, U, V, mat, m)

  IF(DoE) THEN
     IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
     IF( PrtMat > 1 ) CALL matprt(string//label,n,m,mat)
     string = 'E-'
     ee = 0.0_dp
     DO i = 1,n
        ee = ee + wg(i,current_k)*mat(i,i)
     ENDDO
     IF ( PrtMat > 0 ) WRITE(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock('matcalc')

END SUBROUTINE matcalc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matcalc_k (label, DoE, PrtMat, ik, ninner, n, m, U, V, mat, ee)
  !
  USE kinds,                ONLY : dp
  USE io_global,ONLY : stdout
  USE wvfct,                ONLY : wg, npwx
  USE becmod,               ONLY : calbec
  USE noncollin_module,     ONLY : noncolin, npol
  IMPLICIT NONE
  !
  ! compute the (n,n) matrix representation <U|V>
  ! and energy from V (m,n) and U(m,n)
  !
  LOGICAL, INTENT(IN) :: DoE
  INTEGER, INTENT(IN) :: PrtMat, ik, ninner, n, m
  COMPLEX(dp), INTENT(IN) :: U(ninner,n), V(ninner,m)
  COMPLEX(dp), INTENT(OUT):: mat(n,m)
  REAL(DP), INTENT(OUT) :: ee
  CHARACTER(len=*), INTENT(IN) :: label

  INTEGER :: i
  CHARACTER(len=2) :: string

  CALL start_clock('matcalc')

  string = 'M-'
  mat = (0.0_dp, 0.0_dp)
  IF(noncolin) THEN
    noncolin = .false.
    CALL calbec(ninner, U, V, mat, m)
    noncolin = .true.
  ELSE
    CALL calbec(ninner, U, V, mat, m)
  ENDIF

  IF(DoE) THEN
    IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
    IF( PrtMat > 1 ) CALL matprt_k(string//label,n,m,mat)
    string = 'E-'
    ee = 0.0_dp
    DO i = 1,n
      ee = ee + wg(i,ik)*DBLE(mat(i,i))
    ENDDO
    IF ( PrtMat > 0 ) WRITE(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock('matcalc')

END SUBROUTINE matcalc_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matprt(label,n,m,A)
  USE kinds, ONLY : dp
  USE io_global,ONLY : stdout
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n,m
  REAL(dp), INTENT(IN):: A(n,m)
  INTEGER :: i
  CHARACTER(len=50) :: frmt
  CHARACTER(len=*) :: label

  WRITE(stdout,'(A)') label
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f16.10)'
  DO i = 1,n
    WRITE(stdout,frmt) A(i,:)
  ENDDO

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matprt_k(label,n,m,A)
  USE kinds, ONLY : dp
  USE io_global,ONLY : stdout
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n,m
  COMPLEX(dp), INTENT(IN):: A(n,m)
  INTEGER :: i
  CHARACTER(len=50) :: frmt
  CHARACTER(len=*) :: label

  WRITE(stdout,'(A)') label//'(real)'
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f12.6)'
  DO i = 1,n
     WRITE(stdout,frmt) dreal(A(i,:))
  ENDDO

  WRITE(stdout,'(A)') label//'(imag)'
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f12.6)'
  DO i = 1,n
     WRITE(stdout,frmt) aimag(A(i,:))
  ENDDO
END SUBROUTINE matprt_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invchol(n,A)
  USE kinds, ONLY : dp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(dp), INTENT(IN):: A(n,n)
!
! given a matrix A, returns the inverse of the Cholesky decomposition of A
! for real matrices
!
  integer :: INFO

  INFO = -1
  CALL DPOTRF( 'L', n, A, n, INFO )
  CALL errinfo('DPOTRF','Cholesky failed in invchol.',INFO)
  INFO = -1
  CALL DTRTRI( 'L', 'N', n, A, n, INFO )
  CALL errinfo('DTRTRI','inversion failed in invchol.',INFO)

END SUBROUTINE invchol

SUBROUTINE invchol_k(n,A)
  USE kinds, ONLY : dp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX(dp), INTENT(IN):: A(n,n)
  !
  ! given a matrix A, returns the inverse of the Cholesky decomposition of A
  ! for cholesky matrices
  !
  INTEGER :: INFO

  INFO = -1
  CALL ZPOTRF( 'L', n, A, n, INFO )
  CALL errinfo('ZPOTRF','Cholesky failed in invchol.',INFO)
  INFO = -1
  CALL ZTRTRI( 'L', 'N', n, A, n, INFO )
  CALL errinfo('ZTRTRI','inversion failed in invchol.',INFO)

END SUBROUTINE invchol_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE errinfo(routine,message,INFO)
  IMPLICIT NONE
  INTEGER :: INFO
  CHARACTER(len=*) :: routine,message

  IF(INFO/=0) THEN
     WRITE(*,*) routine,' exited with INFO= ',INFO
     CALL errore(routine,message,1)
  ENDIF

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
