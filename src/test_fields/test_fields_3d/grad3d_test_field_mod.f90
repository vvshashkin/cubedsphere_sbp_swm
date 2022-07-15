module grad3d_test_field_mod

!Implements cross-polar flow in horizontal + sine type vertical velocity (no account for orography!)

use test_fields_3d_mod,   only : scalar_field3d_t, vector_field3d_t
use grid_field_mod,       only : tile_field_t
use mesh_mod,             only : tile_mesh_t
use const_mod,            only : pi
use sph_coords_mod,       only : cart2sph, sph2cart_vec
use const_mod,            only : Earth_radii, grav=>Earth_grav, kappa, Rd=>rgaz

implicit none

real(kind=8), parameter :: ps_amp = 0.1 !relative amplitude of ps variation

type, extends(scalar_field3d_t) :: grad3d_test_input_t
    real(kind=8) :: h_top
    real(kind=8) :: T0
    contains
    procedure :: get_scalar_field_tile
end type grad3d_test_input_t

type, extends(vector_field3d_t) :: grad3d_test_out_t
    real(kind=8) :: h_top
    real(kind=8) :: T0
    contains
    procedure :: get_vector_component_tile
end type grad3d_test_out_t

contains

subroutine get_scalar_field_tile(this,f,mesh,halo_width)
    class(grad3d_test_input_t),   intent(in)    :: this
    type(tile_field_t),           intent(inout) :: f
    type(tile_mesh_t),            intent(in)    :: mesh
    integer(kind=4),              intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: lam, phi, ps

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                call cart2sph(mesh%rx(i,j,k), mesh%ry(i,j,k), mesh%rz(i,j,k), lam, phi)
                ps = 1.0_8 + ps_amp*(sin(phi)+cos(phi)*sin(lam))
                f%p(i,j,k) = ps**kappa * exp(-mesh%h(i,j,k)*grav*kappa / (Rd*this%T0))
            end do
        end do
    end do

end subroutine get_scalar_field_tile

subroutine get_vector_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(grad3d_test_out_t),     intent(in)    :: this
    type(tile_field_t),           intent(inout) :: v
    type(tile_mesh_t),            intent(in)    :: mesh
    integer(kind=4),              intent(in)    :: halo_width
    real(kind=8),                 intent(in)    :: base_vec(n_comp, &
                                                     mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                                     mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                                     mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: grad(4)
    real(kind=8)    :: lam, phi, dp_dx, dp_dy, dp_dz, P_exp, ps

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                call cart2sph(mesh%rx(i,j,k), mesh%ry(i,j,k), mesh%rz(i,j,k), lam, phi)
                P_exp = exp(-mesh%h(i,j,k)*grav*kappa / (Rd*this%T0))
                ps    = 1.0_8 + ps_amp*(sin(phi)+cos(phi)*sin(lam))
                dp_dx = (kappa-1.0_8)*ps*ps_amp*cos(lam)*P_exp / Earth_radii
                dp_dy = (kappa-1.0_8)*ps*ps_amp*(cos(phi)-sin(phi)*sin(lam))*P_exp / Earth_radii
                dp_dz =-ps**kappa*P_exp*grav*kappa / (Rd*this%T0)
                call sph2cart_vec(lam, phi, dp_dx, dp_dy, grad(1), grad(2), grad(3))
                grad(4) = dp_dz
                v%p(i,j,k) = sum(grad(1:n_comp)*base_vec(1:n_comp,i,j,k))
            end do
        end do
    end do
end subroutine get_vector_component_tile

end module grad3d_test_field_mod
