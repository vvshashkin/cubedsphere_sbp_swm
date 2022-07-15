module grad_ch_sbp_mod

use domain_mod,             only : domain_t
use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global
use sbp_operator_mod,       only : sbp_operator_t

implicit none

type, public, extends(grad_operator_t) ::  grad_ch_sbp_t
    class(exchange_t),     allocatable :: exch_scalar_interior
    class(sbp_operator_t), allocatable :: sbp_op
contains
    procedure, public :: calc_grad => calc_grad_ch_sbp
end type grad_ch_sbp_t

contains

subroutine calc_grad_ch_sbp(this, gx, gy, f, domain)
    class(grad_ch_sbp_t),   intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: f
    !output:
    type(grid_field_t),     intent(inout) :: gx, gy

    integer(kind=4) :: t

    call this%exch_scalar_interior%do(f, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_on_tile(gx%tile(t), gy%tile(t), f%tile(t),  &
                                this%sbp_op, domain%mesh_xy%scale, &
        domain%mesh_x%tile(t), domain%mesh_y%tile(t), domain%mesh_xy%tile(t))
    end do

end subroutine calc_grad_ch_sbp

subroutine calc_grad_on_tile(gx, gy, f, sbp_op, scale, mesh_x, mesh_y, mesh_xy)

    use mesh_mod,         only : tile_mesh_t
    use tile_mod,         only : tile_t
    use sbp_operator_mod, only : sbp_operator_t

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_xy
    class(sbp_operator_t),  intent(in)    :: sbp_op
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: dx, dy
    integer(kind=4) :: i, j, k

    type(tile_t)    :: dx_tile, dy_tile

    dx_tile = tile_t(is=mesh_y%is, ie=mesh_y%ie, js=mesh_y%js, je=mesh_y%je, ks=1, ke=1)
    dy_tile = tile_t(is=mesh_x%is, ie=mesh_x%ie, js=mesh_x%js, je=mesh_x%je, ks=1, ke=1)

    dx = mesh_xy%hx*scale
    dy = mesh_xy%hy*scale

    do k = mesh_xy%ks, mesh_xy%ke

        dx_tile%ks = k; dx_tile%ke = k
        call sbp_op%apply(gx, dx_tile, mesh_y%nx, 'x', f)
        call sbp_op%add_penalty(gx,dx_tile,mesh_y%nx,'x','at_center',f)

        do j = dx_tile%js, dx_tile%je
            do i = dx_tile%is, dx_tile%ie
                gx%p(i,j,k) = gx%p(i,j,k) / dx
            end do
        end do

        dy_tile%ks = k; dy_tile%ke = k
        call sbp_op%apply(gy, dy_tile, mesh_x%ny, 'y', f)
        call sbp_op%add_penalty(gy,dy_tile,mesh_x%ny,'y','at_center',f)

        do j = dy_tile%js, dy_tile%je
            do i = dy_tile%is, dy_tile%ie
                gy%p(i,j,k) = gy%p(i,j,k) / dy
            end do
        end do
    end do

end subroutine calc_grad_on_tile

end module grad_ch_sbp_mod
