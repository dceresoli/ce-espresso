!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine set_rc_rv()
  !-----------------------------------------------------------------------
  !
  !      input : all-electron wavefunctions + valence states
  !      output: separated core and valence charges 
  !
  use kinds, only : dp
  use ld1_parameters, only : nwfx
  
  use ld1inc, only : grid, aeccharge, aevcharge, nwf, oc, isw, rel, psi, &
                     core_state, aectau, aevtau
  implicit none
  real(dp), allocatable :: gradpsi(:,:)
  integer :: n, ns, is
  !
  !      calculates core charge density
  !
  aevcharge=0.0_DP
  aeccharge=0.0_DP
  do n=1,grid%mesh
     do ns=1,nwf
        if (oc(ns)>0.0_DP) then
           is=isw(ns)
           if (rel==2) then
              if (core_state(ns)) then
                 aeccharge(n)=aeccharge(n) &
                              +oc(ns)*( psi(n,1,ns)**2 + psi(n,2,ns)**2 )
              else
                 aevcharge(n,is)=aevcharge(n,is)+oc(ns)*(psi(n,1,ns)**2 &
                                                       + psi(n,2,ns)**2)
              endif
           else
              if (core_state(ns)) then
                 aeccharge(n) = aeccharge(n) + oc(ns)*psi(n,1,ns)**2
              else
                 aevcharge(n,is) = aevcharge(n,is) + oc(ns)*psi(n,1,ns)**2
              endif
           endif
        endif
     enddo
  enddo
  !
  !      the same for tau
  !
  allocate( gradpsi(grid%mesh,2) )
  aevtau=0.0_DP
  aectau=0.0_DP
  do ns=1,nwf
     if (oc(ns)>0.0_DP) then
        is=isw(ns)
        call grad_log(psi(1,1,ns), gradpsi(1,1), grid%rm1, grid%dx, grid%mesh, 4)
        call grad_log(psi(1,2,ns), gradpsi(1,2), grid%rm1, grid%dx, grid%mesh, 4)
        do n=1,grid%mesh
           if (rel==2) then
              if (core_state(ns)) then
                 aectau(n)=aectau(n) &
                              +oc(ns)*( gradpsi(n,1)**2 + gradpsi(n,2)**2 )
              else
                 aevtau(n,is)=aevtau(n,is)+oc(ns)*(gradpsi(n,1)**2 &
                                                       + gradpsi(n,2)**2)
              endif
           else
              if (core_state(ns)) then
                 aectau(n) = aectau(n) + oc(ns)*gradpsi(n,1)**2
              else
                 aevtau(n,is) = aevtau(n,is) + oc(ns)*gradpsi(n,1)**2
              endif
           endif
        enddo
     endif
  enddo
  return
end subroutine set_rc_rv
