!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhoscf (drhoscf, weight, ik, dbecsum, dpsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE ions_base,            ONLY : nat
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE gvecs,                ONLY : nls
  USE wvfct,                ONLY : npwx, nbnd
  USE uspp_param,           ONLY : nhm
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : ngk,igk_k
  USE qpoint,               ONLY : ikks, ikqs
  USE control_lr,           ONLY : nbnd_occ
  USE mp_bands,             ONLY : me_bgrp, inter_bgrp_comm, ntask_groups
  USE mp,                   ONLY : mp_sum
  USE fft_helper_subroutines

  IMPLICIT NONE
  !
  ! I/O variables
  INTEGER, INTENT (IN) :: ik
  ! input: the k point
  REAL(DP), INTENT (IN) :: weight
  ! input: the weight of the k point
  COMPLEX(DP), INTENT (IN) :: dpsi (npwx,nbnd)
  ! input: the perturbed wfc for the given k point
  COMPLEX(DP), INTENT (INOUT) :: drhoscf (dffts%nnr), dbecsum (nhm*(nhm+1)/2,nat)
  ! input/output: the accumulated change to the charge density and dbecsum
  !
  !   here the local variables
  !
  REAL(DP) :: wgt
  ! the effective weight of the k point

  COMPLEX(DP), ALLOCATABLE :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space
  COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_dpsi(:), tg_drho(:)

  INTEGER :: npw, npwq, ikk, ikq
  INTEGER :: ibnd, ir, ir3, ig, incr, v_siz, idx, ioff, ioff_tg, nxyp
  INTEGER :: right_inc, ntgrp
  ! counters

  CALL start_clock ('incdrhoscf')
  !
  ALLOCATE(dpsic(dffts%nnr))
  ALLOCATE(psi(dffts%nnr))
  !
  wgt = 2.d0 * weight / omega
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  incr = 1
  !
  IF ( dffts%have_task_groups ) THEN
     !
     v_siz = dffts%nnr_tg
     !
     ALLOCATE( tg_psi( v_siz ) )
     ALLOCATE( tg_dpsi( v_siz ) )
     ALLOCATE( tg_drho( v_siz ) )
     !
     incr = fftx_ntgrp(dffts)
     !
  ENDIF
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  do ibnd = 1, nbnd_occ(ikk), incr
     !
     IF ( dffts%have_task_groups ) THEN
        !
        tg_drho=(0.0_DP, 0.0_DP)
        tg_psi=(0.0_DP, 0.0_DP)
        tg_dpsi=(0.0_DP, 0.0_DP)
        !
        ioff   = 0
        CALL tg_get_recip_inc( dffts, right_inc )
        ntgrp = fftx_ntgrp( dffts )
        !
        DO idx = 1, ntgrp
           !
           ! ... dtgs%nogrp ffts at the same time. We prepare both
           ! evc (at k) and dpsi (at k+q)
           !
           IF( idx + ibnd - 1 <= nbnd_occ(ikk) ) THEN
              !
              DO ig = 1, npw
                 tg_psi( nls( igk_k( ig,ikk ) ) + ioff ) = evc( ig, idx+ibnd-1 )
              END DO
              DO ig = 1, npwq
                 tg_dpsi( nls( igk_k( ig,ikq ) ) + ioff ) = dpsi( ig, idx+ibnd-1 )
              END DO
              !
           END IF
           !
           ioff = ioff + right_inc
           !
        END DO
        CALL invfft ('tgWave', tg_psi, dffts)
        CALL invfft ('tgWave', tg_dpsi, dffts)

        do ir = 1, dffts%nr1x * dffts%nr2x * dffts%my_nr3p
           tg_drho (ir) = tg_drho (ir) + wgt * CONJG(tg_psi (ir) ) *  tg_dpsi (ir)
        enddo
        !
        ! reduce the group charge (equivalent to sum over bands of 
        ! orbital group)
        !
        CALL tg_reduce_rho( drhoscf, tg_drho, dffts )
        !
     ELSE
        !
        ! Normal case: no task groups
        !
        ! FFT to R-space of the unperturbed wfct's evc
        !
        psi (:) = (0.d0, 0.d0)
        do ig = 1, npw
           psi (nls (igk_k(ig,ikk) ) ) = evc (ig, ibnd)
        enddo
        CALL invfft ('Wave', psi, dffts)
        !
        ! FFT to R-space of the perturbed wfct's dpsi
        !
        dpsic(:) = (0.d0, 0.d0)
        do ig = 1, npwq
           dpsic (nls (igk_k(ig,ikq) ) ) = dpsi (ig, ibnd)
        enddo
        CALL invfft ('Wave', dpsic, dffts)
        !
        ! Calculation of the response charge-density
        !
        do ir = 1, dffts%nnr
           drhoscf (ir) = drhoscf (ir) + wgt * CONJG(psi (ir) ) * dpsic (ir)
        enddo
        !
     ENDIF
     !
  enddo ! loop on bands
  !
  ! Ultrasoft contribution
  ! Calculate dbecsum = <evc|vkb><vkb|dpsi>
  ! 
  CALL addusdbec (ik, weight, dpsi, dbecsum)
  !
  DEALLOCATE(psi)
  DEALLOCATE(dpsic)
  !
  IF ( dffts%have_task_groups ) THEN
     DEALLOCATE(tg_psi)
     DEALLOCATE(tg_dpsi)
     DEALLOCATE(tg_drho)
  ENDIF
  !
  CALL stop_clock ('incdrhoscf')
  !
  RETURN
  !
end subroutine incdrhoscf
