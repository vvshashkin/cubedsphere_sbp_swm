module buffer_mod

use grid_field_mod, only : tile_field_t

implicit none

type, public :: buffer_t
    real(kind=8), allocatable :: p(:)
contains
    procedure, public :: init
end type buffer_t

contains

subroutine init(this, n)

    class(buffer_t) :: this
    integer(kind=4), intent(in) :: n

    allocate(this%p(1:n))

end subroutine init

subroutine unpack_from_buf(f, buf, is, ie, js, je, ks, ke, pts_num)

    type(tile_field_t), intent(inout) :: f
    integer(kind=4),    intent(in)    :: pts_num
    real(kind=8),       intent(inout) :: buf(pts_num)
    integer(kind=4),    intent(in)    :: is, ie, js, je, ks, ke

    integer(kind=4) :: ind, i, j ,k, idx

    idx = 0

    do k = ks, ke
        do j = js, je
            do i = is, ie
                idx = idx + 1
                f%p(i,j,k) = buf(idx)
            end do
        end do
    end do

    if (idx /= pts_num) print*, 'Error in unpacking!'

end subroutine unpack_from_buf

subroutine pack_to_buf(f, buf, is, ie, js, je, ks, ke, first_dim_index, send_i_step, send_j_step, pts_num)

    type(tile_field_t),    intent(in)    :: f
    integer(kind=4),       intent(in)    :: pts_num
    real(kind=8),          intent(inout) :: buf(pts_num)
    integer(kind=4),       intent(in)    :: is, ie, js, je, ks, ke
    integer(kind=4),       intent(in)    :: send_i_step, send_j_step
    character(len=1),      intent(in)    :: first_dim_index

    integer(kind=4) :: ind, i, j ,k, idx

    idx = 0

    if (first_dim_index == 'i') then
        if (send_j_step == 1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = js, je
                    do i = is, ie
                        idx = idx + 1
                        buf(idx) = f%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = js, je
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx) = f%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i = is, ie
                        idx = idx + 1
                        buf(idx) = f%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx) = f%p(i,j,k)
                    end do
                end do
            end do
        end if

    else if (first_dim_index == 'j') then

        if (send_j_step == 1 .and. send_i_step ==1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = js, je
                        idx = idx + 1
                        buf(idx) = f%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = js, je
                        idx = idx + 1
                        buf(idx) = f%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx) = f%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx) = f%p(i,j,k)
                    end do
                end do
            end do
        end if

        if ( pts_num /= idx ) print*, 'Error in packing!'

    else
        print*, 'Error in packing!!!'
    end if

end subroutine pack_to_buf
subroutine unpack_from_buf_vec(u, v, buf, is, ie, js, je, ks, ke, pts_num)

    type(tile_field_t),   intent(inout) :: u, v
    integer(kind=4),      intent(in)    :: pts_num
    real(kind=8),         intent(inout) :: buf(2*pts_num)
    integer(kind=4),      intent(in)    :: is, ie, js, je, ks, ke

    integer(kind=4) :: ind, i, j ,k, idx

    idx = 0

    do k = ks, ke
        do j = js, je
            do i = is, ie
                idx = idx + 1
                u%p(i,j,k) = buf(idx)
                v%p(i,j,k) = buf(idx+pts_num)
            end do
        end do
    end do

    if (idx /= pts_num) print*, 'Error in unpacking!'

end subroutine unpack_from_buf_vec
subroutine pack_to_buf_vec(u, v, buf, is, ie, js, je, ks, ke, first_dim_index, send_i_step, send_j_step, pts_num)

    type(tile_field_t),    intent(in)    :: u, v
    integer(kind=4),       intent(in)    :: pts_num
    real(kind=8),          intent(inout) :: buf(2*pts_num)
    integer(kind=4),       intent(in)    :: is, ie, js, je, ks, ke
    integer(kind=4),       intent(in)    :: send_i_step, send_j_step
    character(len=1),      intent(in)    :: first_dim_index

    integer(kind=4) :: ind, i, j ,k, idx

    idx = 0

    if (first_dim_index == 'i') then
        if (send_j_step == 1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = js, je
                    do i = is, ie
                        idx = idx + 1
                        buf(idx)         = u%p(i,j,k)
                        buf(idx+pts_num) = v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = js, je
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx)         = -u%p(i,j,k)
                        buf(idx+pts_num) =  v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i = is, ie
                        idx = idx + 1
                        buf(idx)         =  u%p(i,j,k)
                        buf(idx+pts_num) = -v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx)         = -u%p(i,j,k)
                        buf(idx+pts_num) = -v%p(i,j,k)
                    end do
                end do
            end do
        end if

    else if (first_dim_index == 'j') then

        if (send_j_step == 1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = js, je
                        idx = idx + 1
                        buf(idx)         = v%p(i,j,k)
                        buf(idx+pts_num) = u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = js, je
                        idx = idx + 1
                        buf(idx)         =  v%p(i,j,k)
                        buf(idx+pts_num) = -u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx)         = -v%p(i,j,k)
                        buf(idx+pts_num) =  u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx)         = -v%p(i,j,k)
                        buf(idx+pts_num) = -u%p(i,j,k)
                    end do
                end do
            end do
        end if

        if ( pts_num /= idx ) print*, 'Error in packing!'

    else
        print*, 'Error in packing!!!'
    end if

end subroutine pack_to_buf_vec
subroutine pack_to_buf_u_vec(u, v, buf, is, ie, js, je, ks, ke, first_dim_index, send_i_step, send_j_step, pts_num)

    type(tile_field_t),    intent(in)    :: u, v
    integer(kind=4),       intent(in)    :: pts_num
    real(kind=8),          intent(inout) :: buf(pts_num)
    integer(kind=4),       intent(in)    :: is, ie, js, je, ks, ke
    integer(kind=4),       intent(in)    :: send_i_step, send_j_step
    character(len=1),      intent(in)    :: first_dim_index

    integer(kind=4) :: ind, i, j ,k, idx

    idx = 0

    if (first_dim_index == 'i') then
        if (send_j_step == 1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = js, je
                    do i = is, ie
                        idx = idx + 1
                        buf(idx)         = u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = js, je
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx)         = -u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i = is, ie
                        idx = idx + 1
                        buf(idx)         =  u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx)         = -u%p(i,j,k)
                    end do
                end do
            end do
        end if

    else if (first_dim_index == 'j') then

        if (send_j_step == 1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = js, je
                        idx = idx + 1
                        buf(idx)         = v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = js, je
                        idx = idx + 1
                        buf(idx)         =  v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx)         = -v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx)         = -v%p(i,j,k)
                    end do
                end do
            end do
        end if

        if ( pts_num /= idx ) print*, 'Error in packing!'

    else
        print*, 'Error in packing!!!'
    end if

end subroutine pack_to_buf_u_vec
subroutine pack_to_buf_v_vec(u, v, buf, is, ie, js, je, ks, ke, first_dim_index, send_i_step, send_j_step, pts_num)

    type(tile_field_t),    intent(in)    :: u, v
    integer(kind=4),       intent(in)    :: pts_num
    real(kind=8),          intent(inout) :: buf(pts_num)
    integer(kind=4),       intent(in)    :: is, ie, js, je, ks, ke
    integer(kind=4),       intent(in)    :: send_i_step, send_j_step
    character(len=1),      intent(in)    :: first_dim_index

    integer(kind=4) :: ind, i, j ,k, idx

    idx = 0

    if (first_dim_index == 'i') then
        if (send_j_step == 1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = js, je
                    do i = is, ie
                        idx = idx + 1
                        buf(idx)         = v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = js, je
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx)         = v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i = is, ie
                        idx = idx + 1
                        buf(idx)         =  -v%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx)         = -v%p(i,j,k)
                    end do
                end do
            end do
        end if

    else if (first_dim_index == 'j') then

        if (send_j_step == 1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = js, je
                        idx = idx + 1
                        buf(idx)         = u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = js, je
                        idx = idx + 1
                        buf(idx)         =  -u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx)         = u%p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx)         = -u%p(i,j,k)
                    end do
                end do
            end do
        end if

        if ( pts_num /= idx ) print*, 'Error in packing!'

    else
        print*, 'Error in packing!!!'
    end if

end subroutine pack_to_buf_v_vec

end module buffer_mod
