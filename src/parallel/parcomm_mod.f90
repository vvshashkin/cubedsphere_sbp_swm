module parcomm_mod

use mpi

implicit none

type, public :: parcomm_t
    integer(kind=4)                :: myid, np
    integer(kind=mpi_integer_kind) :: comm_w
contains
    procedure, public :: get_mpi_rank
    procedure, public :: get_mpi_proc_number
    procedure, public :: barrier
    procedure, public :: abort
    procedure, public :: print
end type parcomm_t


!single acces point for global parcomm singleton
!can be thought as an analogue of use mpi
type(parcomm_t) :: parcomm_global

contains

subroutine init_global_parallel_enviroment()

    integer(kind=4) :: ierr

    call MPI_init(ierr)

    parcomm_global%comm_w = mpi_comm_world
    parcomm_global%myid   = parcomm_global%get_mpi_rank()
    parcomm_global%np     = parcomm_global%get_mpi_proc_number()

end subroutine init_global_parallel_enviroment

subroutine deinit_global_parallel_enviroment()

    integer(kind=4) :: ierr

    call parcomm_global%barrier()
    call mpi_finalize(ierr)

end subroutine deinit_global_parallel_enviroment

function get_mpi_rank(this) result(myid)

    class(parcomm_t), intent(in) :: this
    integer(kind=4) :: myid, ierr

    call MPI_comm_rank(this%comm_w, myid, ierr)

end function get_mpi_rank

function get_mpi_proc_number(this) result(np)

    class(parcomm_t), intent(in) :: this
    integer(kind=4) :: np, ierr

    call MPI_comm_size(this%comm_w, np, ierr)

end function get_mpi_proc_number

subroutine barrier(this)

    class(parcomm_t), intent(in) :: this
    integer(kind=4) :: np, ierr

    call mpi_barrier(this%comm_w, ierr)

end subroutine barrier

subroutine abort(this, error_message)

    class(parcomm_t), intent(in) :: this
    character(len=*), intent(in) :: error_message

    integer(kind=4) :: code, ierr

    print*, error_message
    call mpi_abort(this%comm_w, code, ierr)

end subroutine abort

subroutine print(this, message)

    class(parcomm_t), intent(in) :: this
    character(len=*), intent(in) :: message

    if (this%myid == 0) print*, message

end subroutine print

end module parcomm_mod
