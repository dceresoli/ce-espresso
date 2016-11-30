! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! Jan Kunes jan.kunes@physik.uni-augsburg.de
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

program maxent
	use kinds
	use io
	use system
	
	implicit none

		
	call read_input()	
	
	call maxent_calc(nsigma,time_slices+1,1,2.d-3,beta,energy_steps,delta_e,maxent_sweeps,1200.d0,1,12345)
	
	call program_end()
end program maxent

