subroutine avost(err_msg)

use mpi

implicit none

integer(kind=4) ierr
character(*) :: err_msg

print *, err_msg
call MPI_FINALIZE(ierr)
stop

end subroutine avost
