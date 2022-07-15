module v_nabla_sbp_mod

use abstract_v_nabla_mod, only : v_nabla_operator_t
use grid_field_mod,       only : grid_field_t, tile_field_t
use mesh_mod,             only : mesh_t, tile_mesh_t
use sbp_operator_mod,     only : sbp_operator_t

implicit none

type, public, extends(v_nabla_operator_t) :: v_nabla_sbp_operator_t
    class(sbp_operator_t), allocatable :: sbp_op
contains
    procedure, public  :: calc_v_nabla => calc_v_nabla
end type v_nabla_sbp_operator_t

contains

subroutine calc_v_nabla(this, f_tend, f, ut, vt, mesh)

    class(v_nabla_sbp_operator_t), intent(inout) :: this
    type(grid_field_t),            intent(in)    :: f      !advected field
    type(grid_field_t),            intent(in)    :: ut, vt !contravariant components
    type(mesh_t),                  intent(in)    :: mesh
    type(grid_field_t),            intent(inout) :: f_tend !advective tendency

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call calc_v_nabla_tile(f_tend%tile(t), f%tile(t), ut%tile(t), vt%tile(t), &
                                    this%sbp_op, mesh%tile(t), mesh%scale)
    end do
end subroutine calc_v_nabla

subroutine calc_v_nabla_tile(f_tend, f, ut, vt, sbp_op, mesh, scale)

    use tile_mod, only : tile_t

    type(tile_field_t),            intent(in)    :: f      !advected field
    type(tile_field_t),            intent(in)    :: ut, vt !contravariant components
    class(sbp_operator_t),         intent(in)    :: sbp_op
    type(tile_mesh_t),             intent(in)    :: mesh
    real(kind=8),                  intent(in)    :: scale
    type(tile_field_t),            intent(inout) :: f_tend !advective tendency

    real(kind=8)    :: Dx(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    real(kind=8)    :: Dy(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    type(tile_t)    :: dxdy_tile
    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8)    :: hx

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    dxdy_tile = tile_t(is = is, ie=ie, js=js, je=je,ks = 1, ke=1)

    hx = mesh%hx

    do k = ks, ke

        dxdy_tile%ks = k; dxdy_tile%ke = k

        call sbp_op%apply(Dx, dxdy_tile, dxdy_tile, mesh%nx, 'x', f)
        call sbp_op%apply(Dy, dxdy_tile, dxdy_tile, mesh%ny, 'y', f)

        do j = js, je
            do i = is, ie
                f_tend%p(i,j,k) =-(ut%p(i,j,k)*Dx(i,j,1)+vt%p(i,j,k)*Dy(i,j,1)) / (hx*scale)
            end do
        end do

    end do

end subroutine calc_v_nabla_tile

end module v_nabla_sbp_mod
