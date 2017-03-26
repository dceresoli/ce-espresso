!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE cp_restart_new
  !-----------------------------------------------------------------------------
  !
  ! ... This module contains subroutines to write and read data required to
  ! ... restart a calculation from the disk
  !
#if !defined(__OLDXML)
  !
  USE iotk_module
  USE qes_module
  USE qexsd_input, ONLY: qexsd_init_k_points_ibz
  USE qexsd_module, ONLY: qexsd_init_schema, qexsd_openschema, qexsd_closeschema,      &
                          qexsd_init_convergence_info, qexsd_init_algorithmic_info,    & 
                          qexsd_init_atomic_species, qexsd_init_atomic_structure,      &
                          qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
                          qexsd_init_magnetization,qexsd_init_band_structure,          &
                          qexsd_init_dipole_info, qexsd_init_total_energy,             &
                          qexsd_init_forces,qexsd_init_stress,                         &
                          qexsd_init_outputElectricField, input_obj => qexsd_input_obj
  USE io_files,  ONLY : iunpun, xmlpun_schema, prefix, tmp_dir, qexsd_fmt,&
       qexsd_version
  USE io_base,   ONLY : write_wfc
  USE xml_io_base,     ONLY  : read_wfc, write_rho_xml,read_print_counter, create_directory
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE parser,    ONLY : version_compare
  USE matrix_inversion
  !
  IMPLICIT NONE
  !
  SAVE
  !
  !
  INTEGER, PRIVATE :: iunout
  !
  ! variables to describe qexml current version
  ! and back compatibility
  !
  LOGICAL, PRIVATE :: qexml_version_before_1_4_0 = .FALSE.
  !
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_writefile( ndw, ascii, nfi, simtime, acc, nk, xk,          &
                             wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,   &
                             taui, cdmi, stau0, svel0, staum, svelm, force,  &
                             vnhp, xnhp0, xnhpm, nhpcl, nhpdim, occ0, occm,  &
                             lambda0,lambdam, xnhe0, xnhem, vnhe, ekincm,    &
                             et, rho, c02, cm2, ctot, iupdwn, nupdwn,        &
                             iupdwn_tot, nupdwn_tot, wfc, mat_z ) ! BS added wfc
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : gamma_only, force_pairing, trhow, &
                                           tksw, twfcollect, do_makov_payne, &
                                           smallmem, llondon, lxdm, ts_vdw,  &
                                           tfor, tpre
      USE control_flags,            ONLY : lwfpbe0nscf, lwfnscf, lwf ! Lingzhu Kong
      USE constants,                ONLY : e2
      USE dener,                    ONLY : detot
      USE io_files,                 ONLY : psfile, pseudo_dir, iunwfc, &
                                           nwordwfc, tmp_dir, diropn
      USE mp_images,                ONLY : intra_image_comm, me_image, &
                                           nproc_image
      USE mp_pools,                 ONLY : nproc_pool, intra_pool_comm, root_pool, inter_pool_comm
      USE mp_bands,                 ONLY : me_bgrp, nproc_bgrp, &
                                           my_bgrp_id, intra_bgrp_comm, &
                                           inter_bgrp_comm, root_bgrp, &
                                           ntask_groups
      USE mp_diag,                  ONLY : nproc_ortho
      USE mp_world,                 ONLY : world_comm, nproc
      USE run_info,                 ONLY : title
      USE gvect,                    ONLY : ngm, ngm_g
      USE gvecs,                    ONLY : ngms_g, ecuts, dual
      USE gvecw,                    ONLY : ngw, ngw_g, ecutwfc
      USE gvect,                    ONLY : ig_l2g, mill
      USE electrons_base,           ONLY : nspin, nelt, nel, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, s_to_r, ainv ! BS added ainv
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           amass, iforce, ind_bck
      USE funct,                    ONLY : get_dft_name, get_inlc, &
           dft_is_hybrid, get_exx_fraction, get_screening_parameter, &
           dft_is_nonlocc, get_nonlocc_name
      USE ldaU_cp,                  ONLY : lda_plus_U, ns, ldmx,Hubbard_l, &
                                           Hubbard_lmax, Hubbard_U
      USE energies,                 ONLY : enthal, ekin, eht, esr, eself, &
                                           epseu, enl, exc, vave
      USE mp,                       ONLY : mp_sum, mp_barrier
      USE fft_base,                 ONLY : dfftp, dffts, dfftb
      USE uspp_param,               ONLY : n_atom_wfc, upf
      USE global_version,           ONLY : version_number
      USE cp_main_variables,        ONLY : descla
      USE cp_interfaces,            ONLY : collect_lambda, collect_zmat
      USE kernel_table,             ONLY : vdw_table_name, kernel_file_name
      USE london_module,            ONLY : scal6, lon_rcut, in_c6
      USE tsvdw_module,             ONLY : vdw_isolated, vdw_econv_thr
      USE wrappers,                 ONLY : f_copy
      USE uspp,                     ONLY : okvan
      USE input_parameters,         ONLY : vdw_corr, london, starting_ns_eigenvalue
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN) :: ndw          !
      LOGICAL,               INTENT(IN) :: ascii        !
      INTEGER,               INTENT(IN) :: nfi          ! index of the current step
      REAL(DP),              INTENT(IN) :: simtime      ! simulated time
      REAL(DP),              INTENT(IN) :: acc(:)       !  
      INTEGER,               INTENT(IN) :: nk           ! number of kpoints
      REAL(DP),              INTENT(IN) :: xk(:,:)      ! k-points coordinates 
      REAL(DP),              INTENT(IN) :: wk(:)        ! k-points weights
      REAL(DP),              INTENT(IN) :: ht(3,3)      ! 
      REAL(DP),              INTENT(IN) :: htm(3,3)     ! 
      REAL(DP),              INTENT(IN) :: htvel(3,3)   ! 
      REAL(DP),              INTENT(IN) :: gvel(3,3)    ! 
      REAL(DP),              INTENT(IN) :: xnhh0(3,3)   ! 
      REAL(DP),              INTENT(IN) :: xnhhm(3,3)   ! 
      REAL(DP),              INTENT(IN) :: vnhh(3,3)    ! 
      REAL(DP),              INTENT(IN) :: taui(:,:)    ! 
      REAL(DP),              INTENT(IN) :: cdmi(:)      ! 
      REAL(DP),              INTENT(IN) :: stau0(:,:)   ! 
      REAL(DP),              INTENT(IN) :: svel0(:,:)   ! 
      REAL(DP),              INTENT(IN) :: staum(:,:)   ! 
      REAL(DP),              INTENT(IN) :: svelm(:,:)   ! 
      REAL(DP),              INTENT(IN) :: force(:,:)   ! 
      REAL(DP),              INTENT(IN) :: xnhp0(:)     ! 
      REAL(DP),              INTENT(IN) :: xnhpm(:)     ! 
      REAL(DP),              INTENT(IN) :: vnhp(:)      ! 
      INTEGER,               INTENT(IN) :: nhpcl        ! 
      INTEGER,               INTENT(IN) :: nhpdim       ! 
      REAL(DP),              INTENT(IN) :: occ0(:)      !  occupations of electronic states
      REAL(DP),              INTENT(IN) :: occm(:)      ! 
      REAL(DP),              INTENT(IN) :: lambda0(:,:,:) ! 
      REAL(DP),              INTENT(IN) :: lambdam(:,:,:) ! 
      REAL(DP),              INTENT(IN) :: xnhe0        ! 
      REAL(DP),              INTENT(IN) :: xnhem        ! 
      REAL(DP),              INTENT(IN) :: vnhe         ! 
      REAL(DP),              INTENT(IN) :: ekincm       ! 
      REAL(DP),              INTENT(IN) :: et(:,:)      !  eigenvalues
      REAL(DP),              INTENT(IN) :: rho(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: c02(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: cm2(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: ctot(:,:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: nupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn_tot(:)! 
      INTEGER,               INTENT(IN) :: nupdwn_tot(:)! 
      REAL(DP),              INTENT(IN) :: wfc(:,:)     ! BS 
      REAL(DP),    OPTIONAL, INTENT(IN) :: mat_z(:,:,:) ! 
      !
      LOGICAL               :: write_charge_density
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname
      CHARACTER(LEN=320)    :: filename, sourcefile
      CHARACTER(LEN=4)      :: cspin
      INTEGER               :: kunit, ib, ik_eff
      INTEGER               :: k1, k2, k3
      INTEGER               :: nk1, nk2, nk3
      INTEGER               :: j, i, iss, ig, nspin_wfc, iss_wfc
      INTEGER               :: is, ia, isa, ik, ierr
      INTEGER,  ALLOCATABLE :: ftmp(:,:)
      INTEGER,  ALLOCATABLE :: ityp(:)
      REAL(DP), ALLOCATABLE :: tau(:,:)
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      REAL(DP)              :: omega, htm1(3,3), h(3,3)
      REAL(DP)              :: a1(3), a2(3), a3(3)
      REAL(DP)              :: b1(3), b2(3), b3(3)
      REAL(DP)              :: nelec
      REAL(DP)              :: scalef
      LOGICAL               :: lsda
      REAL(DP)              :: s0, s1, cclock
      INTEGER               :: nbnd_tot
      INTEGER               :: natomwfc, nbnd_
      REAL(DP), ALLOCATABLE :: mrepl(:,:)
      CHARACTER(LEN=256)    :: tmp_dir_save
      LOGICAL               :: exst
      INTEGER               :: inlc
      CHARACTER(iotk_attlenx)  :: attr
      REAL(DP), ALLOCATABLE :: temp_vec(:), wfc_temp(:,:) ! BS 
      TYPE(output_type) :: output
      LOGICAL :: is_hubbard(nsp)
      REAL(dp):: hubbard_dum(3,nsp)
      COMPLEX(dp), ALLOCATABLE :: ns_dum(:,:,:,:)
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      k1  = 0
      k2  = 0
      k3  = 0
      nk1 = 0
      nk2 = 0
      nk3 = 0
      !
      ! ... subroutine body
      !
      write_charge_density = trhow
      !
      IF( nspin > 1 .AND. .NOT. force_pairing ) THEN
         !
         !  check if the array storing wave functions is large enought
         !
         IF( SIZE( c02, 2 ) < ( iupdwn( 2 ) + nupdwn(1) - 1 ) ) &
            CALL errore('cp_writefile',' wrong wave functions dimension ', 1 )
         !
      END IF
      !
      IF(  nupdwn_tot(1) < nupdwn(1) ) &
         CALL errore( " writefile ", " wrong number of states ", 1 )
      !
      nbnd_    = nupdwn(1) 
      nbnd_tot = MAX( nupdwn(1), nupdwn_tot(1) )
      nelec = nelt
      !
      ! ... Cell related variables
      ! ... Dirty trick to avoid bogus complaints because ht in intent(in)
      !
      h = ht
      CALL invmat( 3, h, htm1, omega )
      h = TRANSPOSE( ht )
      !
      a1 = ht(1,:)/alat
      a2 = ht(2,:)/alat
      a3 = ht(3,:)/alat
      !
      CALL recips( a1, a2, a3, b1, b2, b3 )
      !
      ! ... Compute array ityp, and tau
      !
      ALLOCATE( ityp( nat ) )
      ALLOCATE( tau( 3, nat ) )
      !
      isa = 0
      !
      DO is = 1, nsp
         !
         DO ia = 1, na(is)
            !
            isa = isa + 1
            ityp(isa) = is
            !
         END DO
         !
      END DO
      !
      natomwfc =  n_atom_wfc ( nat, ityp ) 
      !
      CALL s_to_r( stau0, tau, na, nsp, h )
      !   
      lsda = ( nspin == 2 )
      !
      ALLOCATE( ftmp( nbnd_tot , nspin ) )
      !
      ftmp = 0.0d0
      !
      DO iss = 1, nspin
         !
         ftmp( 1:nupdwn(iss), iss ) = occ0( iupdwn(iss) : iupdwn(iss) + nupdwn(iss) - 1 )
         !
      END DO
      !
      ! XML descriptor
      ! 
      WRITE(dirname,'(A,A,"_",I2,".save/")') TRIM(tmp_dir), TRIM(prefix), ndw
      CALL create_directory( TRIM(dirname) )
      !
      CALL qexsd_init_schema( iunpun )
      !
      IF ( ionode ) THEN
         !
         ! ... here we init the variables and finally write them to file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL qexsd_openschema(TRIM( dirname ) // TRIM( xmlpun_schema ))
         output%tagname="output"
!-------------------------------------------------------------------------------
! ... CONVERGENCE_INFO - TO BE VERIFIED
!-------------------------------------------------------------------------------
!
         CALL qexsd_init_convergence_info(output%convergence_info, &
              n_scf_steps=0, scf_error=0.0_dp, &
              opt_conv_ispresent=.FALSE., &
              n_opt_steps=0, grad_norm=0.0_dp )
         !
!-------------------------------------------------------------------------------
! ... ALGORITHMIC_INFO
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_algorithmic_info(output%algorithmic_info, &
              real_space_q=.FALSE., uspp=okvan, paw=.FALSE.)
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_SPECIES
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm,&
                 psfile, amass)
!-------------------------------------------------------------------------------
! ... ATOMIC_STRUCTURE
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_atomic_structure(output%atomic_structure, nsp, atm, ityp, &
              nat, tau(:,ind_bck(:)), alat, alat*a1(:), alat*a2(:), alat*a3(:), ibrav)
         !
!-------------------------------------------------------------------------------
! ... SYMMETRIES
!-------------------------------------------------------------------------------
         output%symmetries%lwrite=.false.
!-------------------------------------------------------------------------------
! ... BASIS SET
!-------------------------------------------------------------------------------
         CALL qexsd_init_basis_set(output%basis_set,gamma_only, ecutwfc/e2, ecutwfc*dual/e2, &
              dfftp%nr1, dfftp%nr2, dfftp%nr3, dffts%nr1, dffts%nr2, dffts%nr3, &
              .FALSE., dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, ngms_g, ngw_g, &
              b1(:), b2(:), b3(:) )
!-------------------------------------------------------------------------------
! ... XC FUNCTIONAL
!-------------------------------------------------------------------------------
         dft_name = get_dft_name()
         IF ( london ) vdw_corr = 'grimme-d2' ! why this?
         is_hubbard(:) = (Hubbard_U(:) > 0.0_dp)
         hubbard_dum(:,:)= 0.0_dp
         ALLOCATE (ns_dum(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
         ns_dum= (0.0_dp, 0.0_dp)
         CALL qexsd_init_dft(output%dft, dft_name, .true., dft_is_hybrid(), &
              0, 0, 0, ecutwfc, get_exx_fraction(), get_screening_parameter(),&
              'none', .false., 0.0_dp, lda_plus_u, 0, 2*Hubbard_lmax+1,.false.,&
              nspin, nsp, 2*Hubbard_lmax+1, nat, atm, ityp, Hubbard_U,&
              Hubbard_dum(1,:), Hubbard_dum(2,:), Hubbard_dum(3,:),Hubbard_dum,&
              starting_ns_eigenvalue, ns, ns_dum, 'atomic', dft_is_nonlocc(), &
              TRIM(vdw_corr), TRIM ( get_nonlocc_name()), scal6, in_c6, &
              lon_rcut, 0.0_dp, 0.0_dp, vdw_econv_thr, vdw_isolated, &
              is_hubbard, upf(1:nsp)%psd)
         DEALLOCATE(ns_dum)
!-------------------------------------------------------------------------------
! ... MAGNETIZATION
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_magnetization(output%magnetization, lsda, .false.,&
              .false., 0.0_dp, [0.0_dp,0.0_dp, 0.0_dp], 0.0_dp, .false.)
         !
!-------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------
         CALL  qexsd_init_total_energy(output%total_energy,enthal, 0.0_dp, eht,&
              vave, exc, 0.0_dp, 0.0_dp, 0.0_dp)
!-------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------
         ! TEMP
         CALL qexsd_init_k_points_ibz( input_obj%k_points_ibz, 'Gamma', &
              'CP',nk1,nk2,nk3,k1,k2,k3,1,xk,wk,alat,a1,.false.) 
         input_obj%bands%occupations%tagname="occupations"
         input_obj%bands%occupations%lread=.false.
         input_obj%bands%occupations%lwrite=.true.
         input_obj%bands%occupations%spin_ispresent=.false.
         input_obj%bands%occupations%occupations="fixed"
         ! TEMP
         CALL  qexsd_init_band_structure(output%band_structure,lsda, .false., &
              .false., nbnd_, nelec, natomwfc, .true., 0.0_dp , .false., & 
              [0.0_dp,0.0_dp], et,  DBLE( ftmp ), nspin, xk, [ngw_g], wk,  &
              STARTING_KPOINTS = input_obj%k_points_IBZ, &
              OCCUPATION_KIND = input_obj%bands%occupations, &
              WF_COLLECTED = twfcollect)
!-------------------------------------------------------------------------------
! ... FORCES
!-------------------------------------------------------------------------------
         !
         output%forces_ispresent=tfor
         CALL qexsd_init_forces(output%forces,nat,force,tfor)
         !
!-------------------------------------------------------------------------------
! ... STRESS - TO BE VERIFIED
!-------------------------------------------------------------------------------
         output%stress_ispresent=tpre
         ! may be wrong or incomplete
         IF ( tpre) h = -MATMUL( detot, ht ) / omega
         CALL qexsd_init_stress(output%stress, h, tpre ) 
!-------------------------------------------------------------------------------
! ... ACTUAL WRITING
!-------------------------------------------------------------------------------
         !
         CALL qes_write_output(iunpun,output)
         CALL qes_reset_output(output)
         !
!-------------------------------------------------------------------------------
! ... CLOSING
!-------------------------------------------------------------------------------
         !
         CALL qexsd_closeschema()
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... WRITE WFC
!-------------------------------------------------------------------------------
      DO iss = 1, nspin
         !
         ik_eff = iss
         filename = TRIM(dirname) // 'wfc' // TRIM(int_to_char(ik_eff))
         CALL write_wfc( iunpun, ik_eff, nk, iss, nspin, &
              c02, ngw_g, gamma_only, nbnd_tot, ig_l2g, ngw,  &
              filename, scalef, ionode, root_pool, intra_pool_comm )
         !
      END DO
!-------------------------------------------------------------------------------
! ... WRITE PSEUDOPOTENTIALS
!-------------------------------------------------------------------------------
     !
     ! ... copy pseudopotential files into the .save directory
     !
     DO is = 1, nsp
        sourcefile= TRIM(pseudo_dir)//psfile(is)
        filename  = TRIM(dirname)//psfile(is)
        IF ( TRIM(sourcefile) /= TRIM(filename) ) &
             ierr = f_copy(sourcefile, filename)
     END DO
     inlc = get_inlc()
     IF ( inlc > 0 ) THEN 
        sourcefile= TRIM(kernel_file_name)
        filename = TRIM(dirname)//TRIM(vdw_table_name)
        IF ( TRIM(sourcefile) /= TRIM(filename) ) & 
           ierr = f_copy(sourcefile, filename)
     END IF  
     !
!-------------------------------------------------------------------------------
! ... CHARGE DENSITY
!-------------------------------------------------------------------------------
      !
      IF (write_charge_density) then
         !
         filename = TRIM( dirname ) // 'charge-density'
         !
         IF ( nspin == 1 ) THEN
            !
            CALL write_rho_xml( filename, rho(:,1), &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, &
                                dfftp%ipp, dfftp%npp, ionode, intra_bgrp_comm, inter_bgrp_comm )
            !
         ELSE IF ( nspin == 2 ) THEN
            !
            ALLOCATE( rhoaux( SIZE( rho, 1 ) ) )
            !
            rhoaux = rho(:,1) + rho(:,2) 
            !
            CALL write_rho_xml( filename, rhoaux, &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, &
                                dfftp%ipp, dfftp%npp, ionode, intra_bgrp_comm, inter_bgrp_comm )
            !
            filename = TRIM( dirname ) // 'spin-polarization'
            !
            rhoaux = rho(:,1) - rho(:,2) 
            !
            CALL write_rho_xml( filename, rhoaux, &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, &
                                dfftp%ipp, dfftp%npp, ionode, intra_bgrp_comm, inter_bgrp_comm )
            !
            DEALLOCATE( rhoaux )
            !
         END IF
         !
      END IF ! write_charge_density

!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( ftmp )
      DEALLOCATE( tau  )
      DEALLOCATE( ityp )
      !
      s1 = cclock() 
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, &
                '(3X,"restart file written in ",F8.3," sec.",/)' ) ( s1 - s0 )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE cp_writefile
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_readfile( ndr, ascii, nfi, simtime, acc, nk, xk,   &
                            wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,     &
                            taui, cdmi, stau0, svel0, staum, svelm, force,    &
                            vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ0, occm,      &
                            lambda0, lambdam, b1, b2, b3, xnhe0, xnhem, vnhe, &
                            ekincm, c02, cm2, wfc, mat_z ) ! added wfc
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : gamma_only, force_pairing, iverbosity, twfcollect, lwf ! BS added lwf
      USE io_files,                 ONLY : iunpun, xmlpun, iunwfc, nwordwfc, &
                                           tmp_dir, diropn
      USE run_info,                 ONLY : title
      USE gvect,                    ONLY : ngm
      USE gvecw,                    ONLY : ngw, ngw_g
      USE electrons_base,           ONLY : nspin, nbnd, nelt, nel, &
                                           nupdwn, iupdwn, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, s_to_r, r_to_s
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           sort_tau, ityp, ions_cofmass
      USE gvect,       ONLY : ig_l2g, mill
      USE cp_main_variables,        ONLY : nprint_nfi, descla
      USE cp_interfaces,            ONLY : distribute_lambda, distribute_zmat
      USE mp,                       ONLY : mp_sum, mp_bcast
      USE mp_global,                ONLY : nproc_file, nproc_pool_file, &
                                           nproc_image_file, ntask_groups_file,&
                                           nproc_bgrp_file, nproc_ortho_file
      USE parameters,               ONLY : ntypx
      USE constants,                ONLY : eps8, angstrom_au, pi
      USE qes_types_module,         ONLY : output_type, parallel_info_type, &
           general_info_type
      USE qexsd_reader_module,      ONLY : qexsd_get_general_info, &
           qexsd_get_parallel_info, qexsd_get_output
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)    :: ndr          !  I/O unit number
      LOGICAL,               INTENT(IN)    :: ascii        !
      INTEGER,               INTENT(INOUT) :: nfi          ! index of the current step
      REAL(DP),              INTENT(INOUT) :: simtime      ! simulated time
      REAL(DP),              INTENT(INOUT) :: acc(:)       !
      INTEGER,               INTENT(IN)    :: nk           ! number of kpoints
      REAL(DP),              INTENT(INOUT) :: xk(:,:)      ! k-points coordinates
      REAL(DP),              INTENT(INOUT) :: wk(:)        ! k-points weights
      REAL(DP),              INTENT(INOUT) :: ht(3,3)      !
      REAL(DP),              INTENT(INOUT) :: htm(3,3)     !
      REAL(DP),              INTENT(INOUT) :: htvel(3,3)   !
      REAL(DP),              INTENT(INOUT) :: gvel(3,3)    !
      REAL(DP),              INTENT(INOUT) :: xnhh0(3,3)   !
      REAL(DP),              INTENT(INOUT) :: xnhhm(3,3)   !
      REAL(DP),              INTENT(INOUT) :: vnhh(3,3)    !
      REAL(DP),              INTENT(INOUT) :: taui(:,:)    !
      REAL(DP),              INTENT(INOUT) :: cdmi(:)      !
      REAL(DP),              INTENT(INOUT) :: stau0(:,:)   !
      REAL(DP),              INTENT(INOUT) :: svel0(:,:)   !
      REAL(DP),              INTENT(INOUT) :: staum(:,:)   !
      REAL(DP),              INTENT(INOUT) :: svelm(:,:)   !
      REAL(DP),              INTENT(INOUT) :: force(:,:)   ! 
      REAL(DP),              INTENT(INOUT) :: xnhp0(:)     !      
      REAL(DP),              INTENT(INOUT) :: xnhpm(:)     ! 
      REAL(DP),              INTENT(INOUT) :: vnhp(:)      !  
      INTEGER,               INTENT(INOUT) :: nhpcl        !  
      INTEGER,               INTENT(INOUT) :: nhpdim       !  
      REAL(DP),              INTENT(INOUT) :: occ0(:)      ! occupations
      REAL(DP),              INTENT(INOUT) :: occm(:)      !
      REAL(DP),              INTENT(INOUT) :: lambda0(:,:,:) !
      REAL(DP),              INTENT(INOUT) :: lambdam(:,:,:) !
      REAL(DP),              INTENT(INOUT) :: b1(3)        !
      REAL(DP),              INTENT(INOUT) :: b2(3)        !
      REAL(DP),              INTENT(INOUT) :: b3(3)        !
      REAL(DP),              INTENT(INOUT) :: xnhe0        !
      REAL(DP),              INTENT(INOUT) :: xnhem        !
      REAL(DP),              INTENT(INOUT) :: vnhe         !  
      REAL(DP),              INTENT(INOUT) :: ekincm       !  
      COMPLEX(DP),           INTENT(INOUT) :: c02(:,:)     ! 
      COMPLEX(DP),           INTENT(INOUT) :: cm2(:,:)     ! 
      REAL(DP),              INTENT(INOUT) :: wfc(:,:)     ! BS 
      REAL(DP),    OPTIONAL, INTENT(INOUT) :: mat_z(:,:,:) ! 
      !
      CHARACTER(LEN=256)   :: dirname, kdirname, filename
      CHARACTER(LEN=5)     :: kindex
      CHARACTER(LEN=4)     :: cspin
      INTEGER              :: strlen
      INTEGER              :: kunit
      INTEGER              :: k1, k2, k3
      INTEGER              :: nk1, nk2, nk3
      INTEGER              :: i, j, iss, ig, nspin_wfc, ierr, ik
      REAL(DP)             :: omega, htm1(3,3), hinv(3,3), scalef
      LOGICAL              :: found
      !
      ! ... variables read for testing pourposes
      !
      INTEGER               :: ibrav_
      CHARACTER(LEN=3)      :: atm_(ntypx)
      CHARACTER(LEN=3), ALLOCATABLE :: symbols(:) 
      INTEGER               :: nat_, nsp_, na_
      INTEGER               :: nk_, ik_, nt_
      LOGICAL               :: gamma_only_ , lsda_
      REAL(DP)              :: alat_, a1_(3), a2_(3), a3_(3)
      REAL(DP)              :: zv_ 
      REAL(DP)              :: ecutwfc_, ecutrho_
      INTEGER               :: nr1,nr2,nr3,nr1s,nr2s,nr3s,nr1b,nr2b,nr3b
      INTEGER               :: ngm_g, ngms_g, npw_g 
      INTEGER               :: iss_, nspin_, ngwt_, nbnd_ , nbnd_tot
      INTEGER               :: nstates_up_ , nstates_dw_ , ntmp, nel_(2)
      REAL(DP)              :: nelec_ 
      REAL(DP)              :: scalef_
      REAL(DP)              :: wk_
      INTEGER               :: nhpcl_, nhpdim_ 
      INTEGER               :: ib, nb
      INTEGER               :: ik_eff
      REAL(DP)              :: amass_(ntypx)
      INTEGER,  ALLOCATABLE :: ityp_(:) 
      INTEGER,  ALLOCATABLE :: isrt_(:) 
      REAL(DP), ALLOCATABLE :: tau_(:,:) 
      REAL(DP), ALLOCATABLE :: occ_(:) 
      INTEGER,  ALLOCATABLE :: if_pos_(:,:) 
      CHARACTER(LEN=256)    :: psfile_(ntypx)
      CHARACTER(LEN=80)     :: pos_unit
      REAL(DP)              :: s1, s0, cclock
      REAL(DP), ALLOCATABLE :: mrepl(:,:) 
      LOGICAL               :: exst, exist_wfc 
      CHARACTER(LEN=256)    :: tmp_dir_save
      INTEGER               :: io_bgrp_id
      TYPE ( output_type)   :: output_obj 
      TYPE (parallel_info_type) :: parinfo_obj
      TYPE (general_info_type ) :: geninfo_obj 
      CHARACTER(iotk_attlenx)  :: attr
      !
      ! ... look for an empty unit
      !
      CALL iotk_free_unit( iunout, ierr )
      CALL errore( 'cp_readfile', &
                   'no free units to read wavefunctions', ierr )
      !
      !
      CALL qexsd_init_schema( iunpun )
      !
      WRITE(dirname,'(A,A,"_",I2,".save/")') TRIM(tmp_dir), TRIM(prefix), ndr
      filename = TRIM( dirname ) // TRIM( xmlpun_schema )
      INQUIRE ( file=filename, exist=found )
      IF (.NOT. found ) &
         CALL errore ('cp_readfile', 'xml data file not found', 1)
      !
      CALL iotk_open_read( iunpun, TRIM(filename) )
      !
      CALL qexsd_get_general_info ( iunpun, geninfo_obj, found)
      IF ( .NOT. found ) THEN
         ierr = ierr + 1
      ELSE
         CALL qexsd_copy_general_info (geninfo_obj, qexsd_fmt, qexsd_version) 
      END IF
      !
      CALL qexsd_get_parallel_info ( iunpun, parinfo_obj, found ) 
      IF ( .NOT. found ) THEN
         ierr = ierr + 10
      ELSE
         CALL  qexsd_copy_parallel_info (parinfo_obj, nproc_file, &
              nproc_pool_file, nproc_image_file, ntask_groups_file, &
              nproc_bgrp_file, nproc_ortho_file)
      END IF
      !
      CALL qexsd_get_output ( iunpun, output_obj, found ) 
      IF ( .NOT. found ) ierr = ierr + 101
      IF ( ierr > 100) CALL errore ('cp_readfile', 'missing data in file', ierr)
      !
      CALL qexsd_copy_atomic_species (output_obj%atomic_species, nsp, atm, &
           psfile_, amass_)
      ALLOCATE ( tau_(3,nat), symbols(nat) )
      CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nat_, &
           tau_, symbols, alat_, a1_, a2_, a3_, ibrav_ )
      CALL qexsd_copy_basis_set ( output_obj%basis_set, gamma_only_, ecutwfc_,&
           ecutrho_, nr1s, nr2s, nr3s, nr1, nr2, nr3, found, nr1b, nr2b, nr3b, &
           ngm_g, ngms_g, npw_g, b1, b2, b3 )
      DEALLOCATE ( tau_, symbols )
      !
      CALL errore('cp_readfile','XSD under development',1)
      !
      RETURN
      !
    END SUBROUTINE cp_readfile
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_general_info (geninfo_obj, qexsd_fmt, qexsd_version) 
    !-------------------------------------------------------------------------------
    ! 
    USE qes_types_module,    ONLY: general_info_type
    !
    IMPLICIT NONE 
    !
    CHARACTER(LEN=*), INTENT(OUT) :: qexsd_fmt, qexsd_version
    TYPE (general_info_type ),INTENT(IN)  :: geninfo_obj   
    !
    qexsd_fmt = TRIM (geninfo_obj%xml_format%NAME)
    qexsd_version = TRIM ( geninfo_obj%xml_format%VERSION)
    !
  END SUBROUTINE qexsd_copy_general_info
  !
  !---------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_parallel_info (parinfo_obj, nproc_file, &
       nproc_pool_file, nproc_image_file, ntask_groups_file, &
       nproc_bgrp_file, nproc_ortho_file)
    !---------------------------------------------------------------------------    !
    USE qes_types_module, ONLY : parallel_info_type
    !
    IMPLICIT NONE 
    !
    TYPE ( parallel_info_type ),INTENT(IN)     :: parinfo_obj
    INTEGER, INTENT(OUT) :: nproc_file, nproc_pool_file, &
                            nproc_image_file, ntask_groups_file, &
                            nproc_bgrp_file, nproc_ortho_file
    ! 
    nproc_file = parinfo_obj%nprocs
    nproc_pool_file = nproc_file/parinfo_obj%npool
    nproc_image_file = nproc_file 
    ntask_groups_file = parinfo_obj%ntasks
    nproc_bgrp_file = nproc_image_file / parinfo_obj%npool / parinfo_obj%nbgrp 
    nproc_ortho_file = parinfo_obj%ndiag
    !
  END SUBROUTINE qexsd_copy_parallel_info  
  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_atomic_species (atomic_species, nsp, atm, psfile, amass)
    !---------------------------------------------------------------------------    !
    USE qes_types_module, ONLY : atomic_species_type
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_species
    INTEGER, INTENT(out) :: nsp
    CHARACTER(LEN=*), INTENT(out) :: atm(:), psfile(:)
    REAL(dp), INTENT(out) :: amass(:)
    !
    INTEGER :: isp
    !
    nsp = atomic_species%ntyp
    DO isp = 1, nsp 
       amass(isp) = 0.d0 
       IF (atomic_species%species(isp)%mass_ispresent) &
            amass(isp) = atomic_species%species(isp)%mass
       atm(isp) = TRIM ( atomic_species%species(isp)%name )
       psfile(isp) = TRIM ( atomic_species%species(isp)%pseudo_file) 
    END DO
    !
  END SUBROUTINE qexsd_copy_atomic_species

  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_atomic_structure (atomic_structure, nat, tau, &
       symbols, alat, a1, a2, a3, ibrav )
  !--------------------------------------------------------------------------
    
    USE qes_types_module, ONLY : atomic_structure_type
    USE constants,        ONLY : pi
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    INTEGER, INTENT(out)  :: nat, ibrav
    REAL(dp), INTENT(out) :: alat, a1(:), a2(:), a3(:), tau(:,:)
    CHARACTER(LEN = 3), INTENT(out) :: symbols(:) 
    INTEGER :: iat, idx
    !
    nat = atomic_structure%nat 
    alat = atomic_structure%alat 
    IF ( atomic_structure%bravais_index_ispresent ) THEN 
       ibrav = atomic_structure%bravais_index 
    ELSE 
       ibrav = 0
    END IF
    loop_on_atoms:DO iat = 1, nat
       idx = atomic_structure%atomic_positions%atom(iat)%index
       tau(:,idx) = atomic_structure%atomic_positions%atom(iat)%atom 
       symbols(idx)  = TRIM ( atomic_structure%atomic_positions%atom(idx)%name)
    END DO loop_on_atoms
    
    IF ( atomic_structure%alat_ispresent ) alat = atomic_structure%alat 
    a1(:) = atomic_structure%cell%a1
    a2(:) = atomic_structure%cell%a2
    a3(:) = atomic_structure%cell%a3

  END SUBROUTINE qexsd_copy_atomic_structure

  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_basis_set ( basis_set, gamma_only, ecutwfc, ecutrho, &
       nr1s, nr2s, nr3s, nr1, nr2, nr3, fft_box_ispresent, nr1b, nr2b, nr3b, &
       ngm_g, ngms_g, npw_g, b1, b2, b3 )
  !--------------------------------------------------------------------------   
    !
    USE qes_types_module, ONLY : basis_set_type
    !
    IMPLICIT NONE 
    TYPE ( basis_set_type ),INTENT(IN)         :: basis_set
    LOGICAL, INTENT(in)   :: fft_box_ispresent
    LOGICAL, INTENT(out)  :: gamma_only
    INTEGER, INTENT(out)  :: ngm_g, ngms_g, npw_g
    INTEGER, INTENT(out)  :: nr1s, nr2s, nr3s, nr1, nr2, nr3
    INTEGER, INTENT(out)  :: nr1b, nr2b, nr3b
    REAL(dp), INTENT(out) :: ecutwfc, ecutrho, b1(:), b2(:), b3(:)
    ! 
    ecutwfc = basis_set%ecutwfc
    ecutrho = basis_set%ecutrho
    gamma_only= basis_set%gamma_only
    nr1 = basis_set%fft_grid%nr1
    nr2 = basis_set%fft_grid%nr2          
    nr3 = basis_set%fft_grid%nr3
    nr1s= basis_set%fft_smooth%nr1
    nr2s= basis_set%fft_smooth%nr2
    nr3s= basis_set%fft_smooth%nr3
    IF ( fft_box_ispresent ) THEN
       nr1b = basis_set%fft_box%nr1
       nr2b = basis_set%fft_box%nr2
       nr3b = basis_set%fft_box%nr3
    END IF
    ngm_g     = basis_set%ngm
    ngms_g    = basis_set%ngms
    npw_g     = basis_set%npwx
    !
    b1 =  basis_set%reciprocal_lattice%b1
    b2 =  basis_set%reciprocal_lattice%b2
    b3 =  basis_set%reciprocal_lattice%b3
    !
  END SUBROUTINE qexsd_copy_basis_set
#endif
    !
  END MODULE cp_restart_new
