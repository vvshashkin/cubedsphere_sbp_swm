module exchange_halo_mod

use grid_field_mod,        only : grid_field_t
use exchange_abstract_mod, only : exchange_t
use buffer_mod,            only : buffer_t, pack_to_buf, unpack_from_buf, &
                                  pack_to_buf_vec, unpack_from_buf_vec
use tile_mod,              only : tile_t
use parcomm_mod,           only : parcomm_t
use mpi

implicit none

type, extends(exchange_t), public :: exchange_2D_halo_t

    type(buffer_t), allocatable :: send_buff(:), recv_buff(:)

    integer(kind=4) :: mpi_message_type = mpi_real8

    integer(kind=4) :: recv_number, send_number

    integer(kind=4), allocatable, dimension(:) :: mpi_send_req, mpi_recv_req

    integer(kind=4), allocatable, dimension(:) :: recv_points_num, send_points_num
    integer(kind=4), allocatable, dimension(:) :: recv_to_tile_ind, send_from_tile_ind
    integer(kind=4), allocatable, dimension(:) :: send_to_proc_id, recv_from_proc_id
    integer(kind=4), allocatable, dimension(:) :: send_tag, recv_tag

    type(tile_t), allocatable :: recv_tile(:), send_tile(:)

    integer(kind=4),  allocatable :: send_i_step(:), send_j_step(:)
    character(len=1), allocatable :: first_dim_index(:)

contains

    procedure, public:: do     => do_halo_exchange
    procedure, public:: do_vec => do_halo_exchange_vec

end type exchange_2D_halo_t

contains

subroutine do_halo_exchange(this, f, parcomm)

    class(exchange_2D_halo_t), intent(inout) :: this
    type(parcomm_t),           intent(in)    :: parcomm
    type(grid_field_t),        intent(inout) :: f

    integer(kind=4) :: ierr, myid
    integer(kind=4) :: i, ind, ind_recv

    this%mpi_send_req = MPI_REQUEST_NULL
    this%mpi_recv_req = MPI_REQUEST_NULL

    do i = 1, this%recv_number
        call MPI_irecv(this%recv_buff(i)%p,        &
                       this%recv_points_num(i),    &
                       this%mpi_message_type,      &
                       this%recv_from_proc_id(i),  &
                       this%recv_tag(i),           &
                       parcomm%comm_w,             &
                       this%mpi_recv_req(i),       &
                       ierr)
    end do

    do i = 1, this%send_number
        call pack_to_buf(f%tile(this%send_from_tile_ind(i)), &
             this%send_buff(i)%p,              &
             this%send_tile(i)%is, this%send_tile(i)%ie, &
             this%send_tile(i)%js, this%send_tile(i)%je, &
             this%send_tile(i)%ks, this%send_tile(i)%ke, &
             this%first_dim_index(i),          &
             this%send_i_step(i),              &
             this%send_j_step(i),              &
             this%send_points_num(i) )
    end do

    do i = 1, this%send_number
        call MPI_isend(this%send_buff(i)%p,      &
                       this%send_points_num(i),  &
                       this%mpi_message_type,    &
                       this%send_to_proc_id(i),  &
                       this%send_tag(i),         &
                       parcomm%comm_w,           &
                       this%mpi_send_req(i),     &
                       ierr )
    end do

    do ind = 1, this%recv_number
        call mpi_waitany(this%recv_number, this%mpi_recv_req, i, mpi_status_ignore, ierr)

        call unpack_from_buf(f%tile(this%recv_to_tile_ind(i)), &
             this%recv_buff(i)%p,              &
             this%recv_tile(i)%is, this%recv_tile(i)%ie, &
             this%recv_tile(i)%js, this%recv_tile(i)%je, &
             this%recv_tile(i)%ks, this%recv_tile(i)%ke, &
             this%recv_points_num(i) )

    end do

    call mpi_waitall(this%send_number, this%mpi_send_req, mpi_statuses_ignore, ierr)

end subroutine do_halo_exchange
subroutine do_halo_exchange_vec(this, u, v, parcomm)

    class(exchange_2D_halo_t), intent(inout) :: this
    type(parcomm_t),           intent(in)    :: parcomm
    type(grid_field_t),        intent(inout) :: u, v

    integer(kind=4) :: ierr
    integer(kind=4) :: i, ind, ind_recv

    this%mpi_send_req = MPI_REQUEST_NULL
    this%mpi_recv_req = MPI_REQUEST_NULL

    do i = 1, this%recv_number
        call MPI_irecv(this%recv_buff(i)%p,        &
                       2*this%recv_points_num(i),  &
                       this%mpi_message_type,      &
                       this%recv_from_proc_id(i),  &
                       this%recv_tag(i),           &
                       parcomm%comm_w,             &
                       this%mpi_recv_req(i),       &
                       ierr)
    end do

    do i = 1, this%send_number
        call pack_to_buf_vec(                     &
             u%tile(this%send_from_tile_ind(i)), &
             v%tile(this%send_from_tile_ind(i)), &
             this%send_buff(i)%p,              &
             this%send_tile(i)%is, this%send_tile(i)%ie, &
             this%send_tile(i)%js, this%send_tile(i)%je, &
             this%send_tile(i)%ks, this%send_tile(i)%ke, &
             this%first_dim_index(i),          &
             this%send_i_step(i),              &
             this%send_j_step(i),              &
             this%send_points_num(i) )
    end do

    do i = 1, this%send_number
        call MPI_isend(this%send_buff(i)%p,      &
                       2*this%send_points_num(i),  &
                       this%mpi_message_type,    &
                       this%send_to_proc_id(i),  &
                       this%send_tag(i),         &
                       parcomm%comm_w,           &
                       this%mpi_send_req(i),     &
                       ierr )
    end do

    do ind = 1, this%recv_number
        call mpi_waitany(this%recv_number, this%mpi_recv_req, i, mpi_status_ignore, ierr)

        call unpack_from_buf_vec(               &
             u%tile(this%recv_to_tile_ind(i)), &
             v%tile(this%recv_to_tile_ind(i)), &
             this%recv_buff(i)%p,               &
             this%recv_tile(i)%is, this%recv_tile(i)%ie, &
             this%recv_tile(i)%js, this%recv_tile(i)%je, &
             this%recv_tile(i)%ks, this%recv_tile(i)%ke, &
             this%recv_points_num(i) )

    end do

    call mpi_waitall(this%send_number, this%mpi_send_req, mpi_statuses_ignore, ierr)

end subroutine do_halo_exchange_vec
end module exchange_halo_mod
