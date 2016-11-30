! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine chem_pot
	use kinds
	use io
	use system
	use constants
	implicit none

	real(DP) :: tol, upper,lower
	logical :: iterate
	integer, parameter :: maxiter = 50
	integer :: iter = 0

	if(chatter) write(stdout,*) 'Chemical potential calculation start'
	
	tol = 1.d-6*(ecut/(DBLE(2*nomega+1)*PI/beta))**3
    tol = MAX(tol,1.d-6)
    
    if(chatter) write(stdout,*) 'Tolerance', tol
	
	upper = mu+2.d1
	lower = mu-2.d1
	iterate = .true.
	
	do while(iterate.AND.iter<maxiter)

		call find_bounds(upper,lower,iterate,tol)
		iter = iter + 1
	end do
	
	mu = (upper+lower)/2.d0

	if(chatter) write(stdout,*) 'Chemical potential calculation end'

end subroutine chem_pot



subroutine find_bounds(upper,lower,result,tol)
	use kinds
	use io
	use system
	use constants
	implicit none

	real(DP) :: upper,lower
	real(DP) :: tol, center
	real(DP), allocatable ::  nel_now(:)
	integer :: i,j,ispin
	logical :: result
	
	allocate(nel_now(nspin))
	
	result = .false.
	nel_now(:) = 0.d0
	    
    center = (upper + lower)/2.d0
    
    do ispin=1,nspin    
    	call nelectrons(ispin,nel_now(ispin),center)
    end do
    
   	if(chatter) then
    		write(stdout,'(a8,f9.5,a8,f9.5,a10,f9.5,a22,f7.3)') 'Upper:', upper, ' Lower:', lower, &
                     'Trial mu:', center, 'Number of electrons:', SUM(nel_now)
	end if
    
    if(ABS(nel_total-SUM(nel_now)).gt.tol) then
    	if(nel_total.gt.SUM(nel_now)) then
    		lower = center
    	else
    		upper = center
    	end if
    	
    	result = .true.
    else
    	result = .false.
    end if
    
    deallocate(nel_now)

end subroutine find_bounds


