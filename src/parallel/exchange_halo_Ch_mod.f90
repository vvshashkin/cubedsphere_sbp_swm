module exchange_halo_Ch_mod

use grid_field_mod,        only : grid_field_t
use exchange_abstract_mod, only : exchange_t
use exchange_halo_mod,     only : exchange_2D_halo_t
use buffer_mod,            only : buffer_t, pack_to_buf, unpack_from_buf, &
                                  pack_to_buf_u_vec, unpack_from_buf_vec, pack_to_buf_v_vec
use tile_mod,              only : tile_t
use parcomm_mod,           only : parcomm_t
use mpi

implicit none

type, extends(exchange_t), public :: exchange_2D_halo_Ch_t

    type(exchange_2D_halo_t) :: exch_u, exch_v

contains

    procedure, public:: do     => do_halo_exchange
    procedure, public:: do_u   => do_halo_exchange_u
    procedure, public:: do_v   => do_halo_exchange_v
    procedure, public:: do_vec => do_halo_exchange_vec

end type exchange_2D_halo_Ch_t

contains
subroutine do_halo_exchange(this, f, parcomm)

    class(exchange_2D_halo_Ch_t), intent(inout) :: this
    type(parcomm_t),             intent(in)    :: parcomm
    type(grid_field_t),          intent(inout) :: f

    integer(kind=4) :: ierr, myid
    integer(kind=4) :: i, ind, ind_recv

    call parcomm%abort('do_halo_exchange is not implemented for exchange_2D_halo_Ch_t! Abort')

end subroutine do_halo_exchange
subroutine do_halo_exchange_u(this, u, v, parcomm)

    class(exchange_2D_halo_Ch_t), intent(inout) :: this
    type(parcomm_t),             intent(in)    :: parcomm
    type(grid_field_t),          intent(inout) :: u, v

    integer(kind=4) :: ierr, myid
    integer(kind=4) :: i, ind, ind_recv

    this%exch_u%mpi_send_req = MPI_REQUEST_NULL
    this%exch_u%mpi_recv_req = MPI_REQUEST_NULL

    do i = 1, this%exch_u%recv_number
        call MPI_irecv(this%exch_u%recv_buff(i)%p,        &
                       this%exch_u%recv_points_num(i),    &
                       this%exch_u%mpi_message_type,      &
                       this%exch_u%recv_from_proc_id(i),  &
                       this%exch_u%recv_tag(i),           &
                       parcomm%comm_w,             &
                       this%exch_u%mpi_recv_req(i),       &
                       ierr)
    end do

    do i = 1, this%exch_u%send_number
        call pack_to_buf_v_vec(                          &
             u%tile(this%exch_u%send_from_tile_ind(i)),         &
             v%tile(this%exch_u%send_from_tile_ind(i)),         &
             this%exch_u%send_buff(i)%p,                        &
             this%exch_u%send_tile(i)%is, this%exch_u%send_tile(i)%ie, &
             this%exch_u%send_tile(i)%js, this%exch_u%send_tile(i)%je, &
             this%exch_u%send_tile(i)%ks, this%exch_u%send_tile(i)%ke, &
             this%exch_u%first_dim_index(i),                    &
             this%exch_u%send_i_step(i),                        &
             this%exch_u%send_j_step(i),                        &
             this%exch_u%send_points_num(i) )
    end do

    do i = 1, this%exch_u%send_number
        call MPI_isend(this%exch_u%send_buff(i)%p,      &
                       this%exch_u%send_points_num(i),  &
                       this%exch_u%mpi_message_type,    &
                       this%exch_u%send_to_proc_id(i),  &
                       this%exch_u%send_tag(i),         &
                       parcomm%comm_w,           &
                       this%exch_u%mpi_send_req(i),     &
                       ierr )
    end do

    do ind = 1, this%exch_u%recv_number
        call mpi_waitany(this%exch_u%recv_number, this%exch_u%mpi_recv_req, i, mpi_status_ignore, ierr)

        call unpack_from_buf(                            &
             v%tile(this%exch_u%recv_to_tile_ind(i)),           &
             this%exch_u%recv_buff(i)%p,                        &
             this%exch_u%recv_tile(i)%is, this%exch_u%recv_tile(i)%ie, &
             this%exch_u%recv_tile(i)%js, this%exch_u%recv_tile(i)%je, &
             this%exch_u%recv_tile(i)%ks, this%exch_u%recv_tile(i)%ke, &
             this%exch_u%recv_points_num(i) )
    end do

    call mpi_waitall(this%exch_u%send_number, this%exch_u%mpi_send_req, mpi_statuses_ignore, ierr)

end subroutine do_halo_exchange_u
subroutine do_halo_exchange_v(this, u, v, parcomm)

    class(exchange_2D_halo_Ch_t), intent(inout) :: this
    type(parcomm_t),             intent(in)    :: parcomm
    type(grid_field_t),          intent(inout) :: u, v

    integer(kind=4) :: ierr, myid
    integer(kind=4) :: i, ind, ind_recv

    this%exch_v%mpi_send_req = MPI_REQUEST_NULL
    this%exch_v%mpi_recv_req = MPI_REQUEST_NULL

    do i = 1, this%exch_v%recv_number
        call MPI_irecv(this%exch_v%recv_buff(i)%p,        &
                       this%exch_v%recv_points_num(i),    &
                       this%exch_v%mpi_message_type,      &
                       this%exch_v%recv_from_proc_id(i),  &
                       this%exch_v%recv_tag(i),           &
                       parcomm%comm_w,             &
                       this%exch_v%mpi_recv_req(i),       &
                       ierr)
    end do

    do i = 1, this%exch_v%send_number
        call pack_to_buf_u_vec(                          &
             u%tile(this%exch_v%send_from_tile_ind(i)),         &
             v%tile(this%exch_v%send_from_tile_ind(i)),         &
             this%exch_v%send_buff(i)%p,                        &
             this%exch_v%send_tile(i)%is, this%exch_v%send_tile(i)%ie, &
             this%exch_v%send_tile(i)%js, this%exch_v%send_tile(i)%je, &
             this%exch_v%send_tile(i)%ks, this%exch_v%send_tile(i)%ke, &
             this%exch_v%first_dim_index(i),                    &
             this%exch_v%send_i_step(i),                        &
             this%exch_v%send_j_step(i),                        &
             this%exch_v%send_points_num(i) )
    end do

    do i = 1, this%exch_v%send_number
        call MPI_isend(this%exch_v%send_buff(i)%p,      &
                       this%exch_v%send_points_num(i),  &
                       this%exch_v%mpi_message_type,    &
                       this%exch_v%send_to_proc_id(i),  &
                       this%exch_v%send_tag(i),         &
                       parcomm%comm_w,           &
                       this%exch_v%mpi_send_req(i),     &
                       ierr )
    end do

    do ind = 1, this%exch_v%recv_number
        call mpi_waitany(this%exch_v%recv_number, this%exch_v%mpi_recv_req, i, mpi_status_ignore, ierr)

        call unpack_from_buf(                            &
             u%tile(this%exch_v%recv_to_tile_ind(i)),           &
             this%exch_v%recv_buff(i)%p,                        &
             this%exch_v%recv_tile(i)%is, this%exch_v%recv_tile(i)%ie, &
             this%exch_v%recv_tile(i)%js, this%exch_v%recv_tile(i)%je, &
             this%exch_v%recv_tile(i)%ks, this%exch_v%recv_tile(i)%ke, &
             this%exch_v%recv_points_num(i) )
    end do

    call mpi_waitall(this%exch_v%send_number, this%exch_v%mpi_send_req, mpi_statuses_ignore, ierr)

end subroutine do_halo_exchange_v
subroutine do_halo_exchange_vec(this, u, v, parcomm)

    class(exchange_2D_halo_Ch_t), intent(inout) :: this
    type(parcomm_t),             intent(in)    :: parcomm
    type(grid_field_t),          intent(inout) :: u, v

    integer(kind=4) :: ierr, myid
    integer(kind=4) :: i, ind, ind_recv

    call this%do_u(u, v, parcomm)
    call this%do_v(u, v, parcomm)

end subroutine do_halo_exchange_vec
end module exchange_halo_Ch_mod
