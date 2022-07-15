module test_fields_3d_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use parcomm_mod,    only : parcomm_global

implicit none

type, abstract :: scalar_field3d_t
contains
    procedure, public :: get_scalar_field
    procedure(generate_scalar_field_tile), deferred :: get_scalar_field_tile
end type scalar_field3d_t

type, abstract :: vector_field3d_t
contains
    procedure, public :: get_vector_field
    procedure, public :: get_vertical_component
    procedure(generate_vector_component_tile), deferred :: get_vector_component_tile
end type vector_field3d_t

type, abstract, extends(vector_field3d_t) :: non_stationary_vector_field3d_t
    real(kind=8) :: t = 0.0_8
end type non_stationary_vector_field3d_t

abstract interface
    subroutine generate_scalar_field_tile(this,f,mesh,halo_width)
        import scalar_field3d_t, tile_field_t, tile_mesh_t
        class(scalar_field3d_t), intent(in)    :: this
        type(tile_field_t),      intent(inout) :: f
        type(tile_mesh_t),       intent(in)    :: mesh
        integer(kind=4),         intent(in)    :: halo_width
    end subroutine generate_scalar_field_tile
    subroutine generate_vector_component_tile(this,v,mesh,halo_width, &
                                              base_vec, n_comp)
        import vector_field3d_t, tile_field_t, tile_mesh_t
        class(vector_field3d_t), intent(in)    :: this
        type(tile_field_t),      intent(inout) :: v
        type(tile_mesh_t),       intent(in)    :: mesh
        integer(kind=4),         intent(in)    :: halo_width
        real(kind=8),            intent(in)    :: base_vec(n_comp, &
                                                           mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width,&
                                                           mesh%js-mesh%halo_width:mesh%je+mesh%halo_width,&
                                                           mesh%ks:mesh%ke)
        integer(kind=4),         intent(in)    :: n_comp
    end subroutine generate_vector_component_tile
end interface

contains

subroutine get_scalar_field(this,f,mesh,halo_width)
    class(scalar_field3d_t), intent(in)    :: this
    type(grid_field_t),      intent(inout) :: f
    type(mesh_t),            intent(in)    :: mesh
    integer(kind=4),         intent(in)    :: halo_width

    integer(kind=4) :: t

    do t=mesh%ts, mesh%te
        call this%get_scalar_field_tile(f%tile(t),mesh%tile(t),halo_width)
    end do
end subroutine get_scalar_field

subroutine get_vector_field(this,u,v,w,mesh_u,mesh_v,mesh_w,halo_width,components_type)
    class(vector_field3d_t), intent(in)    :: this
    type(grid_field_t),      intent(inout) :: u,v,w
    type(mesh_t),            intent(in)    :: mesh_u, mesh_v, mesh_w
    integer(kind=4),         intent(in)    :: halo_width
    character(len=*),        intent(in)    :: components_type

    integer(kind=4) :: t

    select case(components_type)
    case("covariant")
        do t=mesh_u%ts, mesh_u%te
            call this%get_vector_component_tile(u%tile(t),mesh_u%tile(t),halo_width,&
                                                mesh_u%tile(t)%a1,size(mesh_u%tile(t)%a1,1))
        end do
        do t=mesh_v%ts, mesh_v%te
            call this%get_vector_component_tile(v%tile(t),mesh_v%tile(t),halo_width,&
                                                mesh_v%tile(t)%a2,size(mesh_v%tile(t)%a2,1))
        end do
        do t=mesh_w%ts, mesh_w%te
            call this%get_vector_component_tile(w%tile(t),mesh_w%tile(t),halo_width,&
                                                mesh_w%tile(t)%a3,size(mesh_w%tile(t)%a3,1))
        end do
    case("contravariant")
        do t=mesh_u%ts, mesh_u%te
            call this%get_vector_component_tile(u%tile(t),mesh_u%tile(t),halo_width,&
                                                mesh_u%tile(t)%b1,size(mesh_u%tile(t)%b1,1))
        end do
        do t=mesh_v%ts, mesh_v%te
            call this%get_vector_component_tile(v%tile(t),mesh_v%tile(t),halo_width,&
                                                mesh_v%tile(t)%b2,size(mesh_v%tile(t)%b2,1))
        end do
        do t=mesh_w%ts, mesh_w%te
            call this%get_vector_component_tile(w%tile(t),mesh_w%tile(t),halo_width,&
                                                mesh_w%tile(t)%b3,size(mesh_w%tile(t)%b3,1))
        end do
    case default
        call parcomm_global%abort("get_vector_field, unknown components_type: " //&
                                  components_type)
    end select
end subroutine get_vector_field

subroutine get_vertical_component(this,w,mesh_w,halo_width,component_type)
    class(vector_field3d_t), intent(in)    :: this
    type(grid_field_t),      intent(inout) :: w
    type(mesh_t),            intent(in)    :: mesh_w
    integer(kind=4),         intent(in)    :: halo_width
    character(len=*),        intent(in)    :: component_type

    integer(kind=4) :: t

    do t=mesh_w%ts, mesh_w%te
        select case(component_type)
        case("covariant")
            call this%get_vector_component_tile(w%tile(t),mesh_w%tile(t),halo_width,&
                                                mesh_w%tile(t)%a3,size(mesh_w%tile(t)%a3,1))
        case("contravariant")
            call this%get_vector_component_tile(w%tile(t),mesh_w%tile(t),halo_width,&
                                                mesh_w%tile(t)%b3,size(mesh_w%tile(t)%b3,1))
        end select
    end do
end subroutine get_vertical_component

end module test_fields_3d_mod
