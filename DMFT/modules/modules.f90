! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

module kinds
	implicit none
	INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
	INTEGER, PARAMETER :: i4b = selected_int_kind(9) ! integers up to 10^9
end module kinds

module constants
	use kinds
	implicit none
	integer, parameter :: nsigmamax = 10
	integer, parameter :: hdimmax = 20
	
	complex(DP), parameter :: CONE = dcmplx(1d0,0d0)
	complex(DP), parameter :: CZERO = dcmplx(0d0,0d0)
	complex(DP), parameter :: CI = dcmplx(0d0,1d0)
	
	real(DP), PARAMETER :: PI=3.14159265358979323846_DP
	
end module constants

module input_parameters
	use kinds
	use constants
	implicit none

	real :: ecut = 250.0
	real :: beta
	real :: nelec
	real :: start_mu
	logical :: fix_mu
	integer :: &
		time_slices =100, &
		hdim, &
		nspin = 1, &
		nat, & ! Number of atoms
		nsigma ! Number of sigmas
	integer :: block_start(nsigmamax), block_size(nsigmamax) ! First interaction element and size of int block for i-th atom j-th sigma
	real(DP) :: U(nsigmamax,2*hdimmax,2*hdimmax)
	
	integer(i4b) :: nsweeps = 10000
	integer(i4b) :: nwarmup = 500
	integer :: nbin = 0
	real(DP) :: aaa=3.d-2
	
	real :: mixing=0.5
	
	integer :: iqmc(1:nsigmamax) = 3129
	
	integer(i4b) :: sweeps ! for maxent
	integer(i4b) :: energy_steps ! for maxent
	real :: delta_e
	
end module input_parameters

module system
	use kinds
	implicit none
	
	complex(DP), allocatable :: H(:,:,:,:)
	complex(DP), allocatable :: sigma(:,:,:,:,:)
	
	complex(DP), allocatable :: sigma_one(:,:,:,:)
	complex(DP), allocatable :: sigma_inf(:,:,:,:)
	
	real(DP) :: mu
	real(DP) :: nel_total
	integer :: nkp
	real(DP), allocatable :: wk(:),nel(:)

	logical :: fix_mu

	real(DP) :: ecut
	real(DP) :: beta
	integer :: &
		nomega, &
		time_slices, &
		hdim, &
		nspin, &
		nat, & ! Number of atoms
		nsigma, &
		nbin, &
		nsweeps, &
		nwarmup, &
		energy_steps, &
		maxent_sweeps
		
	real(DP) :: mixing
		
	real(DP) :: aaa
	
	integer, allocatable :: block_start(:), block_size(:)
	
	integer, allocatable :: iqmc(:)
	
	real(DP),allocatable :: U(:,:,:)
	
	real(DP) :: delta_e
	
end module system

module io
	implicit none
	
	integer, parameter :: stdout = 6
	integer, parameter :: file_input = 3000
	integer, parameter :: file_hamilt = 3001
	integer, parameter :: file_sigma = 3002
	integer, parameter :: file_Gomega = 3003
	integer, parameter :: file_Gtau = 3004
	integer, parameter :: file_moments = 3005
	integer, parameter :: file_sigma_old = 3006
	integer, parameter :: file_Gomegaout = 3007
	integer, parameter :: file_umatrix = 3008


	integer :: ios
	logical :: chatter = .true.
	
	contains
	
	subroutine error(text, istatus)
		implicit none
		character(len=*), intent(in) :: text
		integer, intent(in) :: istatus
		if(istatus /= 0) call abort!(text)
	end subroutine error
	
	subroutine read_input
	! Reads input parameters from dmft.in and moves to internal variables.
		use kinds
		use input_parameters
		use constants
		use system, only : ecut_ => ecut, &
					beta_ => beta, &
					hdim_ => hdim, &
					nel_total, &
					nel, &
					mu, &
					fix_mu_ => fix_mu, &
					time_slices_ => time_slices, &
					nspin_ => nspin, &
					nat_ => nat, &
					nsigma_ => nsigma, &
					block_start_ => block_start, &
					block_size_ => block_size, &
					nomega, &
					iqmc_ => iqmc, &
					U_ => U, &
					nbin_ => nbin, &
					nsweeps_ => nsweeps, &
					nwarmup_ => nwarmup, &
					aaa_ => aaa, &
					mixing_ =>mixing, &
					energy_steps_ => energy_steps, &
					delta_e_ => delta_e, &
					maxent_sweeps_ => maxent_sweeps
		implicit none
		
		integer :: i,j
		
		namelist / dmftin / ecut, beta, hdim, nelec, start_mu, fix_mu, time_slices, nspin, nat, &
					nsigma, block_start, block_size, iqmc, nbin, nsweeps, nwarmup, aaa, mixing
					
		namelist / maxent / energy_steps, delta_e, sweeps
					
		open (unit = file_input, file = 'dmft.in', status = 'old', form = &
						'formatted', err = 100, iostat = ios)
		
		read(file_input, dmftin, err=102, iostat = ios)

		! Coping input values to internal code variables
		ecut_ = ecut
		hdim_ = hdim
		nel_total = nelec
		mu = start_mu
		fix_mu_ = fix_mu
		time_slices_ = time_slices
		nspin_ = nspin
		nat_ = nat
		nsigma_ = nsigma
		beta_ = beta
		nbin_ = nbin
		nsweeps_ = nsweeps
		nwarmup_ = nwarmup
		aaa_ = aaa
		mixing_ = mixing

		allocate(iqmc_(nsigma))
		iqmc_ = iqmc
		
		nomega = INT(ecut*beta/(2.d0*PI))
		
		allocate(nel(nspin))
		if(fix_mu) then
			nel(:) = 0.d0
		else 
			nel(:) = nel_total/dble(nspin)
		end if
		
		! Information about interaction blocks
		
		allocate(block_start_(nsigma_))
		allocate(block_size_(nsigma_))
		
		block_start_ = CZERO
		block_size_ = CZERO
		
		do i=1, nsigma_
				block_start_(i) = block_start(i)
				block_size_(i) = block_size(i)
		end do
		
		
		!Maxent input
		sweeps = 3000
		delta_e = 0.05
		energy_steps = 200
		
		read(file_input, maxent, err=102, iostat = ios)
		energy_steps_ = energy_steps
		maxent_sweeps_ = sweeps
		delta_e_ = delta_e

		close(file_input)
			
		return

	100 call error('Error opening dmft.in', ios)
	102	call error("Can't read input", ios)

	end subroutine read_input


end module io
