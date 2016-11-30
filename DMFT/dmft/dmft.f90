! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

program dmft
	use kinds
	use io
	use system
	
	real :: time_start, time_end
	integer :: ispin
		
	call CPU_TIME(time_start)
		
	call read_input()	
	
	!Hamiltonian and start Sigma
		
	call read_hamilt()
	
	call read_sigma()
	
	if(chatter) call print_input()
		
	call asymptotics()
	
	if(fix_mu) then
		write(stdout,*) 'Fixed chemical potential'
		do ispin = 1,nspin
			call nelectrons(ispin,nel(ispin),mu)
			write(stdout,*) 'Number of electrons for spin',ispin,' equal to', nel(ispin)
		end do
	else 
		call chem_pot()
		write(stdout,*) 'Chemical potential equal to', mu
	end if

	!Bath green function
	call gloc()

	call program_end()
	
	call CPU_TIME(time_end)
	
	if(chatter) write(stdout,'(x,a17,f7.4,a8)') 'Calculation time:', time_end-time_start, 'seconds'
	
end program dmft

