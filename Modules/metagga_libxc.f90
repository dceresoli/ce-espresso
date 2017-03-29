!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=========================================================================
!
!                        META-GGA FUNCTIONALS
!
!  Available functionals: 
!  - TPSS (Tao, Perdew, Staroverov & Scuseria)
!  - M06L
!
!  Available functionals via LIBXC:
!  - TPSS
!  - TB09
!  - SCAN
!=========================================================================

  
#if defined(__LIBXC)
!======================================================================
! Generic meta-GGA functional via LibXC
!======================================================================
subroutine mgga_libxc(X, C, rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
  USE kinds,            ONLY : DP
  USE funct,            ONLY : libxc_major, libxc_minor, libxc_micro, get_libxc_version
  use xc_f90_types_m
  use xc_f90_lib_m
  implicit none  
  integer, intent(in) :: X, C
  real(dp), intent(in) :: rho, grho, tau
  real(dp), intent(out):: sx, sc, v1x, v2x, v3x, v1c, v2c, v3c

  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  integer, SAVE :: major=0, minor=0, micro=0
  integer :: size = 1
  integer :: flags
  real(dp) :: lapl_rho, vlapl_rho ! not used?

  if (libxc_major == 0) call get_libxc_version
  if (libxc_major < 3 .or. (libxc_major == 3 .and. libxc_minor /= -1)) & 
      call errore('mgga_libxc','please, recompile with LibXC trunk (i.e. >3.0.0))',1)
 
  lapl_rho = grho
  
  ! on input: tau has been already converted from Rydberg to Hartree
  ! on output: energy and potential will be converted from Hartree to Rydberg by
  !            the v_of_rho_meta routine. v3x and v3x will be converted from
  !            Hartree^{-1} to Rydberg^{-1} by v_of_rho_meta.

  ! exchange
  sx = 0.d0
  v2x = 0.d0
  v3x = 0.d0
  if (X == 0) goto 10
  call xc_f90_func_init(xc_func, xc_info, X, XC_UNPOLARIZED)
  flags = xc_f90_info_flags(xc_info)
  if (and(flags, 1) == 1) then
     call xc_f90_mgga_exc_vxc(xc_func, size, rho, grho, lapl_rho, 1.0d0*tau,&
                              sx, v1x, v2x, vlapl_rho, v3x)  
     sx = sx * rho
  else
     call xc_f90_mgga_vxc(xc_func, size, rho, grho, lapl_rho, 1.0d0*tau,&
                          v1x, v2x, vlapl_rho, v3x)  
     sx = 0.d0
  endif
  call xc_f90_func_end(xc_func)
  !!v2x = v2x*2.d0   ! Hartree to Rydberg
  !!v3x = v3x*0.5d0  ! Hartree^{-1} to Rydberg^{-1}
  10 continue

  ! correlation  
  sc = 0.d0
  v2c = 0.d0
  v3c = 0.d0
  if (C == 0) goto 20
  call xc_f90_func_init(xc_func, xc_info, C, XC_UNPOLARIZED)   
  flags = xc_f90_info_flags(xc_info)
  if (and(flags, 1) == 1) then
     call xc_f90_mgga_exc_vxc(xc_func, size, rho, grho, lapl_rho, 1.0d0*tau, &
                              sc, v1c, v2c, vlapl_rho, v3c)  
     sc = sc * rho
  else
     call xc_f90_mgga_vxc(xc_func, size, rho, grho, lapl_rho, 1.0d0*tau, &
                          v1c, v2c, vlapl_rho, v3c)
     sc = 0.d0
  endif
  call xc_f90_func_end(xc_func)
  !!v2c = v2c*2.d0   ! Hartree to Rydberg
  !!v3c = v3c*0.5d0  ! Hartree^{-1} to Rydberg^{-1}
  20 continue

end subroutine mgga_libxc


!======================================================================
! TPSS meta-GGA
!======================================================================
subroutine tpsscxc(rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
  USE kinds,            ONLY : DP
  implicit none  
  real(DP), intent(in) :: rho, grho, tau
  real(dp), intent(out):: sx, sc, v1x, v2x, v3x, v1c, v2c, v3c
  ! XC_MGGA_X_TPSS = 202
  ! XC_MGGA_C_TPSS = 231
  call mgga_libxc(202, 0, rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
end subroutine tpsscxc


!======================================================================
! M06L meta-GGA
!======================================================================
subroutine m06lxc(rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
  USE kinds,            ONLY : DP
  implicit none  
  real(DP), intent(in) :: rho, grho, tau
  real(dp), intent(out):: sx, sc, v1x, v2x, v3x, v1c, v2c, v3c
  ! XC_MGGA_X_M06_L = 203
  ! XC_MGGA_C_M06_L = 233
  call mgga_libxc(203, 233, rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
end subroutine m06lxc


!======================================================================
! TB09 meta-GGA
!======================================================================
subroutine tb09cxc(rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
  USE kinds,            ONLY : DP
  implicit none  
  real(DP), intent(in) :: rho, grho, tau
  real(dp), intent(out):: sx, sc, v1x, v2x, v3x, v1c, v2c, v3c
  ! XC_MGGA_X_TB09 = 208
  ! XC_MGGA_C_TPSS = 231
  call mgga_libxc(208, 0, rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
end subroutine tb09cxc


!======================================================================
! SCAN meta-GGA
!======================================================================
subroutine SCANcxc(rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
  USE kinds,            ONLY : DP
  implicit none  
  real(DP), intent(in) :: rho, grho, tau
  real(dp), intent(out):: sx, sc, v1x, v2x, v3x, v1c, v2c, v3c
  ! XC_MGGA_X_SCAN = 263
  ! XC_MGGA_C_SCAN = 267
  call mgga_libxc(263, 267, rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
end subroutine SCANcxc


#endif

