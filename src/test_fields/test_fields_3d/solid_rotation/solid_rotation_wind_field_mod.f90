module solid_rotation_wind_field_mod

use test_fields_3d_mod, only : non_stationary_vector_field3d_t
use mesh_mod,           only : tile_mesh_t
use grid_field_mod,     only : tile_field_t
use const_mod,          only : pi, Day24h_sec

implicit none

type, extends(non_stationary_vector_field3d_t) :: solid_rotation_wind_field_t
    real(kind=8) :: w_max = 0.0_8
    real(kind=8) :: tau   = 12.0_8*Day24h_sec
    real(kind=8) :: Lz    = 10e3
    real(kind=8) :: u0    = 1.0_8
    real(kind=8) :: alpha = 0.0_8
    contains
    procedure :: get_vector_component_tile
end type solid_rotation_wind_field_t

contains

subroutine get_vector_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(solid_rotation_wind_field_t),  intent(in)    :: this
    type(tile_field_t),                  intent(inout) :: v
    type(tile_mesh_t),                   intent(in)    :: mesh
    integer(kind=4),                     intent(in)    :: halo_width
    real(kind=8), intent(in)    :: base_vec(n_comp, &
                                            mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                            mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                            mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: vs(1:4)
    real(kind=8)    :: rot_axis(3)

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    rot_axis(1:3) = [sin(this%alpha), 0.0_8, cos(this%alpha)]

    do k=ks,ke
        do j=js,je
            do i=is,ie
                vs(1) = this%u0*(rot_axis(2)*mesh%rz(i,j,k)-rot_axis(3)*mesh%ry(i,j,k))
                vs(2) = this%u0*(-rot_axis(1)*mesh%rz(i,j,k)+rot_axis(3)*mesh%rx(i,j,k))
                vs(3) = this%u0*(rot_axis(1)*mesh%ry(i,j,k)-rot_axis(2)*mesh%rx(i,j,k))
                vs(4) = this%w_max*cos(2.0_8*pi*this%t/this%tau)*sin(pi*mesh%h(i,j,k) / this%Lz)
                v%p(i,j,k) = sum(vs(1:n_comp)*base_vec(1:n_comp,i,j,k))
            end do
        end do
    end do
end subroutine get_vector_component_tile

end module solid_rotation_wind_field_mod
