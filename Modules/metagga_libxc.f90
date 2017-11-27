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

module libxc_version
  implicit none
  save

  integer :: libxc_major=0, libxc_minor=0, libxc_micro=0
  public :: libxc_major, libxc_minor, libxc_micro, get_libxc_version

  contains

  subroutine get_libxc_version
     implicit none
     interface
        subroutine xc_version(major, minor, micro) bind(c)
           use iso_c_binding
           integer(c_int) :: major, minor, micro
        end subroutine xc_version
     end interface
     call xc_version(libxc_major, libxc_minor, libxc_micro)
  end subroutine get_libxc_version

end module libxc_version


!======================================================================
! Generic meta-GGA functional via LibXC
!======================================================================
subroutine mgga_libxc(func_id, has_exc, rho, grho, tau, s, v1, v2, v3)
  USE kinds,            ONLY : DP
  use xc_f90_types_m
  use xc_f90_lib_m
  use libxc_version
  implicit none  
  integer, intent(in) :: func_id
  logical, intent(in) :: has_exc
  real(dp), intent(in) :: rho, grho, tau
  real(dp), intent(out):: s, v1, v2, v3

  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  integer :: size = 1
  real(dp) :: lapl_rho, vlapl_rho ! not used?
  lapl_rho = grho

  ! check libxc version
  if (libxc_major == 0) call get_libxc_version
  if (libxc_major < 3 .or. (libxc_major == 3 .and. libxc_minor /= -1)) & 
      call errore('mgga_libxc','please, recompile with LibXC trunk (i.e. >3.0.0))',1)
  
  ! on input, tau has been already converted from Rydberg to Hartree
  s = 0.d0
  v1 = 0.d0
  v2 = 0.d0
  v3 = 0.d0
  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
  if (has_exc) then
     call xc_f90_mgga_exc_vxc(xc_func, size, rho, grho, lapl_rho, tau, &
                              s, v1, v2, vlapl_rho, v3)  
     s = s * rho
  else
     call xc_f90_mgga_vxc(xc_func, size, rho, grho, lapl_rho, tau, &
                          v1, v2, vlapl_rho, v3)  
     s = 0.d0
  endif
  call xc_f90_func_end(xc_func)
  v2 = v2*2.d0   ! Hartree to Rydberg
  v3 = v3*0.5d0  ! Hartree^{-1} to Rydberg^{-1}
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
  call mgga_libxc(202, .true., rho, grho, tau, sx, v1x, v2x, v3x)
  ! XC_MGGA_C_TPSS = 231
  call mgga_libxc(231, .true., rho, grho, tau, sc, v1c, v2c, v3c)
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
  call mgga_libxc(203, .true., rho, grho, tau, sx, v1x, v2x, v3x)
  ! XC_MGGA_C_M06_L = 233
  call mgga_libxc(233, .true., rho, grho, tau, sc, v1c, v2c, v3c)
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
  call mgga_libxc(208, .false., rho, grho, tau, sx, v1x, v2x, v3x)
  ! XC_MGGA_C_TPSS = 231
  call mgga_libxc(231, .true., rho, grho, tau, sc, v1c, v2c, v3c)
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
  call mgga_libxc(263, .true., rho, grho, tau, sx, v1x, v2x, v3x)
  ! XC_MGGA_C_SCAN = 267
  call mgga_libxc(267, .true., rho, grho, tau, sc, v1c, v2c, v3c)
end subroutine SCANcxc


#endif

