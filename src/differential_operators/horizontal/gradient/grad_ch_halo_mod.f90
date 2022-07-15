module grad_ch_halo_mod

use domain_mod,             only : domain_t
use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use halo_mod,               only : halo_t
use parcomm_mod,            only : parcomm_global

implicit none

type, public, extends(grad_operator_t) ::  grad_ch_halo_t
    integer(kind=4)             :: order
    integer(kind=4)             :: input_halo_width
    integer(kind=4)             :: output_halo_width
    class(halo_t),  allocatable :: halo
contains
    procedure, public :: calc_grad => calc_grad_ch_halo
end type grad_ch_halo_t

contains

subroutine calc_grad_ch_halo(this, gx, gy, f, domain)
    class(grad_ch_halo_t),  intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: f
    !output:
    type(grid_field_t),     intent(inout) :: gx, gy

    integer(kind=4) :: t

    if(this%input_halo_width /= 0) &
        call this%halo%get_halo_scalar(f,domain,this%input_halo_width)

    do t = domain%partition%ts, domain%partition%te
        select case (this%order)
        case(2)
            call calc_grad_on_tile2(gx%tile(t), gy%tile(t), f%tile(t),             &
                                domain%mesh_x%tile(t),  domain%mesh_y%tile(t), &
                                domain%mesh_xy%tile(t), domain%mesh_xy%scale)
        case(4)
            call calc_grad_on_tile4(gx%tile(t), gy%tile(t), f%tile(t),             &
                                    domain%mesh_x%tile(t),  domain%mesh_y%tile(t), &
                                    domain%mesh_xy%tile(t), domain%mesh_xy%scale)
        case default
            call parcomm_global%abort("grad_ch_halo, unknown oreder")
        end select
    end do

end subroutine calc_grad_ch_halo

subroutine calc_grad_on_tile2(gx, gy, f, mesh_x, mesh_y, mesh_xy, scale)

    use mesh_mod,         only : tile_mesh_t
    use tile_mod,         only : tile_t
    use sbp_operator_mod, only : sbp_operator_t

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_xy
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: dx, dy
    integer(kind=4) :: i, j, k

    type(tile_t)    :: dx_tile, dy_tile

    do k=mesh_xy%ks, mesh_xy%ke
        do j=mesh_y%js,mesh_y%je
            do i = mesh_y%is, mesh_y%ie
                gx%p(i,j,k) = (f%p(i+1,j,k)-f%p(i,j,k)) / (scale*mesh_y%hx)
            end do
        end do
        do j=mesh_x%js,mesh_x%je
            do i = mesh_x%is, mesh_x%ie
                gy%p(i,j,k) = (f%p(i,j+1,k)-f%p(i,j,k)) / (scale*mesh_x%hy)
            end do
        end do
    end do

end subroutine calc_grad_on_tile2

subroutine calc_grad_on_tile4(gx, gy, f, mesh_x, mesh_y, mesh_xy, scale)

    use mesh_mod,         only : tile_mesh_t
    use tile_mod,         only : tile_t
    use sbp_operator_mod, only : sbp_operator_t

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_xy
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: dx, dy
    integer(kind=4) :: i, j, k

    type(tile_t)    :: dx_tile, dy_tile

    do k=mesh_xy%ks, mesh_xy%ke
        do j=mesh_y%js,mesh_y%je
            do i = mesh_y%is, mesh_y%ie
                gx%p(i,j,k) = (-f%p(i+2,j,k)+27.0_8*f%p(i+1,j,k)-  &
                                27.0_8*f%p(i,j,k)+f%p(i-1,j,k)) / &
                                                    (24.0_8*scale*mesh_y%hx)
            end do
        end do
        do j=mesh_x%js,mesh_x%je
            do i = mesh_x%is, mesh_x%ie
                gy%p(i,j,k) = (-f%p(i,j+2,k)+27.0_8*f%p(i,j+1,k)-  &
                                27.0_8*f%p(i,j,k)+f%p(i,j-1,k)) / &
                                                    (24.0_8*scale*mesh_x%hx)
            end do
        end do
    end do

end subroutine calc_grad_on_tile4

end module grad_ch_halo_mod
