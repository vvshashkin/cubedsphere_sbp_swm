module div_ah_sbp_mod

use domain_mod,             only : domain_t
use abstract_div_mod,       only : div_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use halo_mod,               only : halo_vec_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global
use sbp_operator_mod,       only : sbp_operator_t
use halo_mod,               only : halo_t

implicit none

type, public, extends(div_operator_t) :: div_ah_sbp_t
    class(exchange_t), allocatable     :: exch_uv_interior
    class(halo_t), allocatable         :: sync_edges
    class(sbp_operator_t), allocatable :: sbp_op
contains
    procedure, public :: calc_div => calc_div_ah_sbp
end type div_ah_sbp_t

contains

subroutine calc_div_ah_sbp(this, div, u, v, domain)
    class(div_ah_sbp_t),    intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    !output
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4), parameter :: halo_width = 1
    integer(kind=4) :: t

    call this%exch_uv_interior%do_vec(u,v,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_div_on_tile(div%tile(t), u%tile(t), v%tile(t),   &
                              domain%mesh_xy%tile(t), this%sbp_op, &
                              domain%mesh_xy%scale)
    end do

    call this%sync_edges%get_halo_scalar(div,domain,1)

end subroutine calc_div_ah_sbp

subroutine calc_div_on_tile(div, u, v, mesh, sbp_op,scale)

    use mesh_mod, only : tile_mesh_t
    use tile_mod, only : tile_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: u, v
    type(tile_mesh_t),      intent(in)    :: mesh
    class(sbp_operator_t),  intent(in)    :: sbp_op
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k
    real(kind=8)    :: Gu(u%is:u%ie,u%js:u%je,1)
    real(kind=8)    :: Gv(v%is:v%ie,v%js:v%je,1)
    real(kind=8)    :: Dx(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    real(kind=8)    :: Dy(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    type(tile_t)    :: dxdy_tile, Gu_tile, Gv_tile

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx

    dxdy_tile = tile_t(is = is, ie=ie, js=js, je=je,ks = 1, ke=1)
    Gu_tile   = tile_t(is = u%is, ie=u%ie, js=u%js, je=u%je, ks = 1, ke=1)
    Gv_tile   = tile_t(is = v%is, ie=v%ie, js=v%js, je=v%je, ks = 1, ke=1)

    do k = ks, ke

        do j = max(js-3,1), min(je+3,mesh%ny)
            do i = max(is-3,1),min(ie+3,mesh%nx)
                Gu(i,j,1) = mesh%J(i,j,k)*u%p(i,j,k)
            end do
        end do
        call sbp_op%apply_array_to_array(Dx, dxdy_tile, dxdy_tile, mesh%nx, 'x', Gu, Gu_tile)

        do j = max(js-3,1), min(je+3,mesh%nx)
            do i = max(is-3,1),min(ie+3,mesh%nx)
                Gv(i,j,1) = mesh%J(i,j,k)*v%p(i,j,k)
            end do
        end do
        call sbp_op%apply_array_to_array(Dy, dxdy_tile, dxdy_tile, mesh%ny, 'y', Gv, Gv_tile)

        do j = js, je
            do i = is,ie
                div%p(i,j,k) = (Dx(i,j,1)+Dy(i,j,1)) / (mesh%J(i,j,k)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile

end module div_ah_sbp_mod
