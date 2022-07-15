module global_diag_mod

use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t

implicit none

!type integral_t
!    integer(kind=4) :: mpi_reduce_op !code of mpi-reduce operation
!    logical         :: l_only_master !if we need integral value only on master process
!    logical         :: calc_after_
!    procedure(calc_integral_local_i), nopass, pointer :: calc_integral_local
!    procedure(calc_after_reduce_i), nopass, pointer   :: calc_after_reduce => null()
!
!    contains
!    procedure calc_integral => calc_integral
!
!end type integral_t

abstract interface
    function calc_integral_local_i(stvec, model_params) result(f)
        import state_abstract_t, model_parameters_abstract_t
        real(kind=8), allocatable                      :: f(:)
        class(state_abstract_t),            intent(in) :: stvec
        class(model_parameters_abstract_t), intent(in) :: model_params
    end function calc_integral_local_i

    function calc_after_reduce_i(f, model_params) result(q)
        import model_parameters_abstract_t
        real(kind=8), allocatable                      :: q(:)
        real(kind=8),                       intent(in) :: f(:)
        class(model_parameters_abstract_t), intent(in) :: model_params
    end function calc_after_reduce_i
end interface

contains

function calc_global_diag(myid, master_id, lall_reduce, stvec, model_params, &
                       calc_local, mpi_reduce_op, calc_after_reduce) result(value)

    use mpi

    real(kind=8), allocatable                      :: value(:)
    integer(kind=4),                    intent(in) :: myid, master_id
    logical,                            intent(in) :: lall_reduce
    class(state_abstract_t),            intent(in) :: stvec
    class(model_parameters_abstract_t), intent(in) :: model_params
    procedure(calc_integral_local_i),   pointer, &
                                        intent(in) :: calc_local
    integer(kind=4),                    intent(in) :: mpi_reduce_op
    procedure(calc_after_reduce_i),     pointer, optional, &
                                        intent(in) :: calc_after_reduce

    real(kind=8), allocatable :: f(:), reduced_f(:)
    integer(kind=4)           :: val_vec_size, ierr

    f = calc_local(stvec, model_params)

    val_vec_size = size(f)
    allocate(reduced_f(val_vec_size))

    if(lall_reduce) then
        call mpi_allreduce(f,reduced_f,val_vec_size, MPI_DOUBLE, mpi_reduce_op, &
                           MPI_COMM_WORLD, ierr)
    else
        call mpi_reduce(f,reduced_f,val_vec_size, MPI_DOUBLE, mpi_reduce_op, &
                        master_id, MPI_COMM_WORLD, ierr)
    end if

    if(present(calc_after_reduce) .and. (myid == master_id .or. lall_reduce)) then
        value = calc_after_reduce(reduced_f, model_params)
    else
        value = reduced_f
    end if

end function calc_global_diag

end module global_diag_mod
