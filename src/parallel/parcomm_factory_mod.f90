module parcomm_factory_mod

use parcomm_mod, only : parcomm_t
use mpi

implicit none

contains

subroutine create_parcomm(comm_w, parcomm)

    integer(kind=4), intent(in)  :: comm_w
    type(parcomm_t), intent(out) :: parcomm

    integer(kind=4) :: ierr

    parcomm%comm_w = comm_w!mpi_comm_world
    parcomm%myid = parcomm%get_mpi_rank()
    parcomm%np   = parcomm%get_mpi_proc_number()

end subroutine create_parcomm

subroutine create_group_parcomm(group_parcomm, parent_parcomm, ranks)

    type(parcomm_t), intent(out) :: group_parcomm
    type(parcomm_t), intent(in)  :: parent_parcomm
    integer(kind=4), intent(in)  :: ranks(:)

    integer(kind=4) :: group_comm
    integer(kind=4) :: ierr, r
    logical  :: is_in_ranks
    integer(kind=4), parameter :: color=1234 !Any constant

    is_in_ranks = .false.
    do r = 1,size(ranks)
        if(parent_parcomm%myid == ranks(r)) then
            is_in_ranks = .true.
            exit
        end if
    end do

    if(is_in_ranks) then
        call mpi_comm_split(parent_parcomm%comm_w, color, parent_parcomm%myid, group_comm, ierr)
        call create_parcomm(group_comm,group_parcomm)
    else
        call mpi_comm_split(parent_parcomm%comm_w, MPI_UNDEFINED, parent_parcomm%myid, group_comm, ierr)
    end if

end subroutine create_group_parcomm

end module parcomm_factory_mod
