module vertical_test_field_mod

use test_fields_3d_mod,            only : scalar_field3d_t, vector_field3d_t
use abstract_vertical_profile_mod, only : vertical_profile_t
use grid_field_mod,                only : tile_field_t
use mesh_mod,                      only : tile_mesh_t

implicit none

type, extends(scalar_field3d_t) :: vertical_ExnerP_t
    real(kind=8) :: t0, p0
    class(vertical_profile_t), allocatable :: vert_profile

    contains
    procedure :: get_scalar_field_tile
end type vertical_ExnerP_t

type, extends(vector_field3d_t) :: vertical_ExnerP_grad_t
    real(kind=8) :: t0, p0
    class(vertical_profile_t), allocatable :: vert_profile

    contains
    procedure :: get_vector_component_tile
end type vertical_ExnerP_grad_t

type, extends(scalar_field3d_t) :: vertical_theta_t
    real(kind=8) :: t0, p0
    class(vertical_profile_t), allocatable :: vert_profile

    contains
    procedure :: get_scalar_field_tile=>get_theta_field_tile
end type vertical_theta_t

type, extends(vector_field3d_t) :: vertical_theta_grad_t
    real(kind=8) :: t0, p0
    class(vertical_profile_t), allocatable :: vert_profile

    contains
    procedure :: get_vector_component_tile=>get_dthetadz_component_tile
end type vertical_theta_grad_t

contains

subroutine get_scalar_field_tile(this,f,mesh,halo_width)
    class(vertical_ExnerP_t),    intent(in)    :: this
    type(tile_field_t),          intent(inout) :: f
    type(tile_mesh_t),           intent(in)    :: mesh
    integer(kind=4),             intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                f%p(i,j,k) = this%vert_profile%calc_ExnerP(mesh%h(i,j,k),this%t0,this%p0)
            end do
        end do
    end do

end subroutine get_scalar_field_tile

subroutine get_vector_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(vertical_ExnerP_grad_t), intent(in)    :: this
    type(tile_field_t),            intent(inout) :: v
    type(tile_mesh_t),             intent(in)    :: mesh
    integer(kind=4),               intent(in)    :: halo_width
    real(kind=8),            intent(in)    :: base_vec(n_comp, &
                                                       mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                                       mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                                       mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: w

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    if(n_comp<4) then
        v%p(is:ie,js:je,ks:ke) = 0.0_8
        return
    end if

    do k=ks,ke
        do j=js,je
            do i=is,ie
                w = this%vert_profile%calc_dExnerP_dz(mesh%h(i,j,k),this%t0,this%p0)
                v%p(i,j,k) = w*base_vec(4,i,j,k)
            end do
        end do
    end do
end subroutine get_vector_component_tile

subroutine get_theta_field_tile(this,f,mesh,halo_width)
    class(vertical_theta_t),     intent(in)    :: this
    type(tile_field_t),          intent(inout) :: f
    type(tile_mesh_t),           intent(in)    :: mesh
    integer(kind=4),             intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                f%p(i,j,k) = this%vert_profile%calc_theta(mesh%h(i,j,k),this%t0,this%p0)
            end do
        end do
    end do

end subroutine get_theta_field_tile

subroutine get_dthetadz_component_tile(this,v,mesh,halo_width, &
                                     base_vec, n_comp)

    class(vertical_theta_grad_t),  intent(in)    :: this
    type(tile_field_t),            intent(inout) :: v
    type(tile_mesh_t),             intent(in)    :: mesh
    integer(kind=4),               intent(in)    :: halo_width
    real(kind=8),            intent(in)    :: base_vec(n_comp, &
                                                       mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                                       mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                                       mesh%ks:mesh%ke)
    integer(kind=4),         intent(in)    :: n_comp

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: w

    is = mesh%is-halo_width; ie = mesh%ie+halo_width
    js = mesh%js-halo_width; je = mesh%je+halo_width
    ks = mesh%ks; ke = mesh%ke

    if(n_comp<4) then
        v%p(is:ie,js:je,ks:ke) = 0.0_8
        return
    end if

    do k=ks,ke
        do j=js,je
            do i=is,ie
                w = this%vert_profile%calc_dtheta_dz(mesh%h(i,j,k),this%t0,this%p0)
                v%p(i,j,k) = w*base_vec(4,i,j,k)
            end do
        end do
    end do
end subroutine get_dthetadz_component_tile

end module vertical_test_field_mod
