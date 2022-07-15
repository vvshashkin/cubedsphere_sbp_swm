module grad_ah_sbp_mod

use domain_mod,             only : domain_t
use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global
use sbp_operator_mod,       only : sbp_operator_t
use halo_mod,               only : halo_vec_t

implicit none

private

type, public, extends(grad_operator_t) :: grad_ah_sbp_t
    class(exchange_t), allocatable     :: exch_scalar_interior
    !syncronize gradient components accross edges:
    class(sbp_operator_t), allocatable :: sbp_op
    class(halo_vec_t),     allocatable :: sync_edges
contains
    procedure, public :: calc_grad => calc_grad_ah_sbp
end type grad_ah_sbp_t

contains

subroutine calc_grad_ah_sbp(this, gx, gy, f, domain)
    class(grad_ah_sbp_t),   intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: f
    !output:
    type(grid_field_t),     intent(inout) :: gx, gy

    integer(kind=4) :: t

    call this%exch_scalar_interior%do(f,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_on_tile(gx%tile(t), gy%tile(t), f%tile(t),  &
                               domain%mesh_xy%tile(t), this%sbp_op,&
                               domain%mesh_xy%scale)
    end do

    call this%sync_edges%get_halo_vector(gx,gy,domain,0)

end subroutine calc_grad_ah_sbp

subroutine calc_grad_on_tile(gx, gy, f, mesh, sbp_op, scale)

    use mesh_mod,         only : tile_mesh_t
    use tile_mod,         only : tile_t
    use sbp_operator_mod, only : sbp_operator_t

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh
    class(sbp_operator_t),  intent(in)    :: sbp_op
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx, mult_loc
    real(kind=8)    :: fdx, fdy, fdx1, fdy1, fdx0, fdy0
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    real(kind=8)    :: Dx(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    real(kind=8)    :: Dy(mesh%is:mesh%ie,mesh%js:mesh%je,1)

    type(tile_t)    :: dxdy_tile

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    dxdy_tile = tile_t(is = is, ie=ie, js=js, je=je,ks = 1, ke=1)

    hx = mesh%hx

    do k = ks, ke

        dxdy_tile%ks = k; dxdy_tile%ke = k
        call sbp_op%apply(Dx, dxdy_tile, dxdy_tile, mesh%nx, 'x', f)
        call sbp_op%apply(Dy, dxdy_tile, dxdy_tile, mesh%ny, 'y', f)

        do j=js,je
            do i=is,ie
                gx%p(i,j,k) = Dx(i,j,1) / (hx*scale)
                gy%p(i,j,k) = Dy(i,j,1) / (hx*scale)
            end do
        end do
    end do

end subroutine calc_grad_on_tile

end module grad_ah_sbp_mod
