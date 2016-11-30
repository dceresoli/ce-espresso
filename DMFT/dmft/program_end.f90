! Copyright (C) 2009 Dmitry Korotin dmitry@korotin.name,
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.

subroutine program_end
! Subroutine deallocates everything
	use io
	use system
	implicit none
	
	if(allocated(block_start)) deallocate(block_start)
	if(allocated(block_size)) deallocate(block_size)
	if(allocated(wk)) deallocate(wk)
	if(allocated(H)) deallocate(H)
	if(allocated(sigma)) deallocate(sigma)
	if(allocated(sigma_one)) deallocate(sigma_one)
	if(allocated(sigma_inf)) deallocate(sigma_inf)
	if(allocated(nel)) deallocate(nel)

end subroutine program_end
