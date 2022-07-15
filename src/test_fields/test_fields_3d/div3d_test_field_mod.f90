module div3d_test_field_mod

!Implements cross-polar flow in horizontal + sine type vertical velocity (no account for orography!)

use test_fields_3d_mod,   only : scalar_field3d_t, vector_field3d_t
use grid_field_mod,       only : tile_field_t
use mesh_mod,             only : tile_mesh_t
use const_mod,            only : pi
use sph_coords_mod,       only : cart2sph, sph2cart_vec
use const_mod,            only : Earth_radii


implicit none

real(kind=8), parameter :: U0 = 20.0_8, W0 = 0.1_8

type, extends(vector_field3d_t) :: div3d_test_wind_t
    real(kind=8) :: h_top
    contains
    procedure :: get_vector_component_tile
end type div3d_test_wind_t

type, extends(scalar_field3d_t) :: div3d_test_field_t
    real(kind=8) :: h_top
    contains
    procedure :: get_scalar_field_tile
end type div3d_test_field_t

contains

subroutine get_vector_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(div3d_test_wind_t),     intent(in)    :: this
    type(tile_field_t),           intent(inout) :: v
    type(tile_mesh_t),            intent(in)    :: mesh
    integer(kind=4),              intent(in)    :: halo_width
    real(kind=8),                 intent(in)    :: base_vec(n_comp, &
                                                     mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                                     mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                                     mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: velocity(4)
    real(kind=8)    :: lam, phi, u_sph, v_sph

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                call cart2sph(mesh%rx(i,j,k), mesh%ry(i,j,k), mesh%rz(i,j,k), lam, phi)
                u_sph = U0*(sin(phi)*(sin(phi)**2-3.0_8*cos(phi)**2)*sin(lam)-0.5_8*cos(phi))
                v_sph = U0*sin(phi)**2*cos(lam)
                call sph2cart_vec(lam, phi, u_sph, v_sph, velocity(1), velocity(2), velocity(3))
                velocity(4) = W0*sin(4._8*pi*mesh%h(i,j,k) / this%h_top)
                v%p(i,j,k) = sum(velocity(1:n_comp)*base_vec(1:n_comp,i,j,k))
            end do
        end do
    end do
end subroutine get_vector_component_tile

subroutine get_scalar_field_tile(this,f,mesh,halo_width)
    class(div3d_test_field_t),    intent(in)    :: this
    type(tile_field_t),           intent(inout) :: f
    type(tile_mesh_t),            intent(in)    :: mesh
    integer(kind=4),              intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: lam, phi

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                call cart2sph(mesh%rx(i,j,k), mesh%ry(i,j,k), mesh%rz(i,j,k), lam, phi)
                f%p(i,j,k) = W0*4.0_8*pi/this%h_top * cos(4.0_8*pi/this%h_top* mesh%h(i,j,k))-&
                             U0*cos(lam)*sin(phi)*cos(phi) / Earth_radii
            end do
        end do
    end do

end subroutine get_scalar_field_tile

end module div3d_test_field_mod
