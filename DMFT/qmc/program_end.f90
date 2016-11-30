! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine program_end
! Subroutine deallocates everything
	use io
	use system
	implicit none
	
	deallocate(block_start)
	deallocate(block_size)
	deallocate(iqmc)
	deallocate(U)

end subroutine program_end
