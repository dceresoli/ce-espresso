!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! external potential
MODULE extpot
    USE kinds,                          ONLY: DP
    IMPLICIT NONE
    SAVE

    LOGICAL :: textpot = .false.
    CHARACTER(256) :: external_potential
    REAL(DP), ALLOCATABLE :: ext_pot(:)

END MODULE extpot


SUBROUTINE add_extpot(nrxx, vltot, extpot)
    USE kinds,                  ONLY: DP
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nrxx
    REAL(DP), INTENT(INOUT) :: vltot(nrxx)
    REAL(DP), INTENT(IN) :: extpot(nrxx)

    vltot(:) = vltot(:) + extpot(:)

END SUBROUTINE add_extpot


SUBROUTINE read_extpot
    USE kinds,                  ONLY: DP
    USE extpot,                 ONLY: textpot, ext_pot, filename => external_potential
    USE io_global,              ONLY: stdout, ionode
    USE fft_base,               ONLY: dfftp
    USE scatter_mod,            ONLY: scatter_grid
    IMPLICIT NONE
    INTEGER :: i1, i2, i3, n1, n2, n3, err, iupot
    INTEGER, EXTERNAL :: find_free_unit
    REAL(DP), ALLOCATABLE :: raux(:,:,:)

    textpot = .false.
    if (len(trim(filename)) == 0) return

    ! allocate buffer
    allocate( raux(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x) )

    ! try to open file
    if (ionode) then
        iupot = find_free_unit()
        open(unit=iupot, name=filename, status='old', iostat=err)
        if (err /= 0) call errore("read_extpot", "error opening file "//trim(filename), err)


        ! read dimensions
        read(iupot,*,iostat=err) n1, n2, n3
        if (err /= 0) call errore("read_extpot", "error reading dimensions", err)

        if (n1 /= dfftp%nr1 .or. n2 /= dfftp%nr2 .or. n3 /= dfftp%nr3) &
            call errore("read_extpot", "dimensions do not match FFT grid", -1)

        write(stdout,*)
        write(stdout,'(5X,''Reading external potential from file: '',A)') trim(filename)
        write(stdout,*)
        do i3 = 1, n3
            do i2 = 1, n2
                do i1 = 1, n1
                    read(iupot,*,iostat=err) raux(i1,i2,i3)
                    if (err /= 0) call errore("read_extpot", "error reading potential", err)
                enddo
            enddo
        enddo
        close(unit=iupot)
    endif ! ionode

    textpot = .true.

#if defined(__MPI)
    allocate( ext_pot(dfftp%nnr) )
    call scatter_grid(dfftp, reshape(raux,(/dfftp%nr1x*dfftp%nr2x*dfftp%nr3x/)), ext_pot)
#else
    allocate( ext_pot(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
    call DCOPY(dfftp%nnrx, raux(1,1,1), 1, ext_pot(1), 1)
#endif
    deallocate( raux )

END SUBROUTINE read_extpot

