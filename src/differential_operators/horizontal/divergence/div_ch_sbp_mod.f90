module div_ch_sbp_mod

use domain_mod,             only : domain_t
use abstract_div_mod,       only : div_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use halo_mod,               only : halo_vec_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global
use sbp_operator_mod,       only : sbp_operator_t
use halo_mod,               only : halo_t
use mesh_mod,               only : tile_mesh_t

implicit none

type, public, extends(div_operator_t) :: div_ch_sbp_t
    type(grid_field_t)                 :: Gu, Gv
    class(exchange_t), allocatable     :: exch_uv
    class(halo_t), allocatable         :: sync_edges
    class(sbp_operator_t), allocatable :: sbp_op
contains
    procedure, public :: calc_div => calc_div_ch_sbp
end type div_ch_sbp_t

contains

subroutine calc_div_ch_sbp(this, div, u, v, domain)
    class(div_ch_sbp_t),    intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    !output
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4), parameter :: halo_width = 1
    integer(kind=4) :: t

    do t = domain%partition%ts, domain%partition%te
        call multiply_uv_by_G_tile(this%Gu%tile(t), u%tile(t), domain%mesh_y%tile(t))
        call multiply_uv_by_G_tile(this%Gv%tile(t), v%tile(t), domain%mesh_x%tile(t))
    end do

    call this%exch_uv%do_vec(this%Gu,this%Gv,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_div_on_tile(div%tile(t), this%Gu%tile(t), this%Gv%tile(t),   &
                             this%sbp_op,domain%mesh_xy%scale,     &
                             domain%mesh_xy%tile(t))
    end do

    ! call this%sync_edges%get_halo_scalar(div,domain,1)

end subroutine calc_div_ch_sbp

subroutine multiply_uv_by_G_tile(Guv, uv, mesh_uv)

    type(tile_field_t), intent(inout) :: Guv, uv
    type(tile_mesh_t),  intent(in)    :: mesh_uv

    integer(kind=4) :: k, i, j

    do k = mesh_uv%ks, mesh_uv%ke
        do j = mesh_uv%js, mesh_uv%je
            do i = mesh_uv%is, mesh_uv%ie
                Guv%p(i,j,k) = uv%p(i,j,k)*mesh_uv%J(i,j,k)
            end do
        end do
    end do

end subroutine

subroutine calc_div_on_tile(div, Gu, Gv, sbp_op, scale, mesh_xy)

    use tile_mod, only : tile_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: Gu, Gv
    type(tile_mesh_t),      intent(in)    :: mesh_xy
    class(sbp_operator_t),  intent(in)    :: sbp_op
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: Dx(mesh_xy%is:mesh_xy%ie,mesh_xy%js:mesh_xy%je,1)
    real(kind=8)    :: Dy(mesh_xy%is:mesh_xy%ie,mesh_xy%js:mesh_xy%je,1)

    real(kind=8)    :: hx, hy
    integer(kind=4) :: i, j, k
    type(tile_t)    :: div_tile

    hx = mesh_xy%hx*scale
    hy = mesh_xy%hy*scale

    div_tile = tile_t(is=mesh_xy%is, ie=mesh_xy%ie, js=mesh_xy%js, je=mesh_xy%je,ks = 1, ke=1)

    do k = mesh_xy%ks, mesh_xy%ke

        div_tile%ks = k; div_tile%ke = k

        call sbp_op%apply(Dx, div_tile, div_tile, mesh_xy%nx, 'x', Gu)
        call sbp_op%add_penalty(Dx,div_tile,div_tile,mesh_xy%nx,'x','at_interface',Gu)

        call sbp_op%apply(Dy, div_tile, div_tile, mesh_xy%ny, 'y', Gv)
        call sbp_op%add_penalty(Dy,div_tile,div_tile,mesh_xy%ny,'y','at_interface',Gv)

        do j = mesh_xy%js, mesh_xy%je
            do i = mesh_xy%is, mesh_xy%ie
                div%p(i,j,k) = (Dx(i,j,1)/hx + Dy(i,j,1)/hy) / mesh_xy%J(i,j,k)
            end do
        end do
    end do

end subroutine calc_div_on_tile

end module div_ch_sbp_mod
