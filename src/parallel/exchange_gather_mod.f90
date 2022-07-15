module exchange_gather_mod

use grid_field_mod,        only : grid_field_t
use exchange_abstract_mod, only : exchange_t
use buffer_mod,            only : buffer_t, pack_to_buf, unpack_from_buf
use parcomm_mod,           only : parcomm_t
use mpi

implicit none

type, extends(exchange_t), public :: exchange_gather_t

    type(buffer_t), allocatable :: send_buff(:), recv_buff(:)

    integer(kind=4) :: mpi_message_type = mpi_real8

    integer(kind=4) :: master_id
    integer(kind=4) :: recv_number, send_number

    integer(kind=4), allocatable, dimension(:) :: mpi_send_req, mpi_recv_req

    integer(kind=4), allocatable, dimension(:) :: recv_points_num, send_points_num
    integer(kind=4), allocatable, dimension(:) :: recv_to_tile_ind
    integer(kind=4), allocatable, dimension(:) :: send_from_tile_ind
    integer(kind=4), allocatable, dimension(:) :: recv_from_proc_id
    integer(kind=4), allocatable, dimension(:) :: send_tag, recv_tag

    integer(kind=4), allocatable, dimension(:) :: send_is, send_ie, send_js, send_je, send_ks, send_ke
    integer(kind=4), allocatable, dimension(:) :: recv_is, recv_ie, recv_js, recv_je, recv_ks, recv_ke
contains
    procedure, public:: do     => do_gather_exchange
    procedure, public:: do_vec => do_gather_exchange_vec
end type exchange_gather_t

contains

subroutine do_gather_exchange(this, f, parcomm)

    class(exchange_gather_t), intent(inout) :: this
    type(parcomm_t),          intent(in)    :: parcomm
    type(grid_field_t),       intent(inout) :: f

    integer(kind=4) :: ierr
    integer(kind=4) :: i, ind, ind_recv

    if (parcomm%myid == this%master_id) then

        this%mpi_recv_req = MPI_REQUEST_NULL

        do i = 1, this%recv_number
            call MPI_irecv(this%recv_buff(i)%p,        &
                           this%recv_points_num(i),    &
                           this%mpi_message_type,      &
                           this%recv_from_proc_id(i),  &
                           this%recv_tag(i),           &
                           mpi_comm_world,             &
                           this%mpi_recv_req(i),       &
                           ierr)
        end do

        do ind = 1, this%recv_number
            call mpi_waitany(this%recv_number, this%mpi_recv_req, i, mpi_status_ignore, ierr)

            call unpack_from_buf(f%tile(this%recv_to_tile_ind(i)), &
                 this%recv_buff(i)%p,              &
                 this%recv_is(i), this%recv_ie(i), &
                 this%recv_js(i), this%recv_je(i), &
                 this%recv_ks(i), this%recv_ke(i), &
                 this%recv_points_num(i) )

        end do
    else

        this%mpi_send_req = MPI_REQUEST_NULL

        do i = 1, this%send_number
            call pack_to_buf(f%tile(this%send_from_tile_ind(i)), &
                 this%send_buff(i)%p,              &
                 this%send_is(i), this%send_ie(i), &
                 this%send_js(i), this%send_je(i), &
                 this%send_ks(i), this%send_ke(i), &
                 'i',                              &
                 1,                                &
                 1,                                &
                 this%send_points_num(i) )
        end do

        do i = 1, this%send_number
            call MPI_isend(this%send_buff(i)%p,      &
                           this%send_points_num(i),  &
                           this%mpi_message_type,    &
                           this%master_id,           &
                           this%send_tag(i),         &
                           mpi_comm_world,           &
                           this%mpi_send_req(i),     &
                           ierr )
        end do

        call mpi_waitall(this%send_number, this%mpi_send_req, mpi_statuses_ignore, ierr)
    end if

end subroutine do_gather_exchange
subroutine do_gather_exchange_vec(this, u, v, parcomm)

    class(exchange_gather_t), intent(inout) :: this
    type(parcomm_t),          intent(in)    :: parcomm
    type(grid_field_t),       intent(inout) :: u, v

    integer(kind=4) :: code, ierr

    call parcomm%print('Gather vec exchange is not implemented :(. Abort!')
    call mpi_abort(parcomm%comm_w, 111, ierr)

end subroutine do_gather_exchange_vec
end module exchange_gather_mod
