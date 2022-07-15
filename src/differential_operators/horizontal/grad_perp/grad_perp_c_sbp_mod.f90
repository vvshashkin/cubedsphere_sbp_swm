module grad_perp_c_sbp_mod

use grid_field_mod,         only : grid_field_t, tile_field_t
use domain_mod,             only : domain_t
use abstract_grad_perp_mod, only : grad_perp_operator_t
use sbp_operator_mod,       only : sbp_operator_t
use exchange_abstract_mod,  only : exchange_t

implicit none

type, public, extends(grad_perp_operator_t) :: grad_perp_c_sbp_t
    class(exchange_t),     allocatable :: exchange
    class(sbp_operator_t), allocatable :: sbp_diff
contains
    procedure, public :: calc_grad_perp
end type grad_perp_c_sbp_t

contains

subroutine calc_grad_perp(this, gu, gv, w, domain)

    class(grad_perp_c_sbp_t), intent(inout) :: this
    type(domain_t),           intent(in)    :: domain
    type(grid_field_t),       intent(inout) :: gu, gv
    type(grid_field_t),       intent(inout) :: w

    integer(kind=4) :: t

    call this%exchange%do(w,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_perp_c_sbp_tile(gu%tile(t), gv%tile(t), w%tile(t),        &
                                  this%sbp_diff, domain%mesh_x%tile(t),          &
                                  domain%mesh_y%tile(t), domain%mesh_xy%tile(t), &
                                  domain%mesh_xy%scale)
    end do


end subroutine calc_grad_perp

subroutine calc_grad_perp_c_sbp_tile(gu, gv, w, sbp_diff, mesh_x, mesh_y, mesh_xy, scale)
    use tile_mod,       only : tile_t
    use mesh_mod,       only : tile_mesh_t
    use grid_field_mod, only : grid_field_t

    type(tile_field_t),    intent(inout) :: gu, gv, w
    class(sbp_operator_t), intent(in)    :: sbp_diff
    type(tile_mesh_t),     intent(in)    :: mesh_x, mesh_y, mesh_xy
    real(kind=8),          intent(in)    :: scale

    integer(kind=4) :: i, j, k
    integer(kind=4) :: ks, ke
    type(tile_t)    :: work_u, work_v

    ks = mesh_xy%ks
    ke = mesh_xy%ke

    work_u = tile_t(is=mesh_x%is, ie = mesh_x%ie, js = mesh_x%js, je = mesh_x%je, ks = 1, ke = 1)
    work_v = tile_t(is=mesh_y%is, ie = mesh_y%ie, js = mesh_y%js, je = mesh_y%je, ks = 1, ke = 1)

    do k = ks, ke

        work_u%ks = k
        work_u%ke = k

        call sbp_diff%apply(gu, work_u, mesh_x%ny, 'y', w)
        call sbp_diff%add_penalty(gu,work_u,mesh_x%ny,'y','at_center',w)

        work_v%ks = k
        work_v%ke = k

        call sbp_diff%apply(gv, work_v, mesh_y%nx, 'x', w)
        call sbp_diff%add_penalty(gv,work_v,mesh_y%nx,'x','at_center',w)

        do j = mesh_x%js, mesh_x%je
            do i = mesh_x%is, mesh_x%ie
                gu%p(i,j,k) = -gu%p(i,j,k) / (scale*mesh_x%hy*mesh_x%J(i,j,k))
            end do
        end do

        do j = mesh_y%js, mesh_y%je
            do i = mesh_y%is, mesh_y%ie
                gv%p(i,j,k) = gv%p(i,j,k) / (scale*mesh_y%hx*mesh_y%J(i,j,k))
            end do
        end do

    end do

end subroutine calc_grad_perp_c_sbp_tile

end module grad_perp_c_sbp_mod
