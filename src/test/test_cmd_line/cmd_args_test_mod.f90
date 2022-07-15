module cmd_args_test_mod

implicit none

contains

subroutine test_cmd_args()
    use mpi
    use cmd_args_mod, only: cmd_arg_t, get_cmd_args

    type(cmd_arg_t), allocatable :: cmd_args(:)
    integer(kind=4) nargs
    integer(kind=4) i
    integer(kind=4) myid, Np, ierr

    call get_cmd_args(cmd_args, nargs)

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if(myid == Np-1) then
        do i=1,nargs
            print *, cmd_args(i)%str
        end do
    end if
end subroutine test_cmd_args

end module cmd_args_test_mod
