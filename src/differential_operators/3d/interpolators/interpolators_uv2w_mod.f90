module interpolators_uv2w_mod

use abstract_interpolators3d_mod,   only : interpolator_uv2w_t
use grid_field_mod,                 only : grid_field_t
use domain_mod,                     only : domain_t
use abstract_vertical_operator_mod, only : vertical_operator_t
use abstract_interpolators2d_mod,   only : interpolator2d_vec2vec_t

implicit none

type, extends(interpolator_uv2w_t) :: uv2w_colocated_t
contains
    procedure :: interp_uv2w => interp_uv2w_colocated
end type uv2w_colocated_t

type, extends(interpolator_uv2w_t) :: uv2w_hor_colocated_t
    class(vertical_operator_t), allocatable :: p2w_oper
contains
    procedure :: interp_uv2w => interp_uv2w_hor_colocated
end type uv2w_hor_colocated_t

type, extends(interpolator_uv2w_t) :: uv2w_staggered_t
    class(vertical_operator_t), allocatable      :: p2w_oper
    type(grid_field_t)                           :: up, vp
    class(interpolator2d_vec2vec_t), allocatable :: uv2p_oper
contains
    procedure :: interp_uv2w => interp_uv2w_staggered
end type uv2w_staggered_t

contains

subroutine interp_uv2w_colocated(this, uw, vw, u, v, domain)
    class(uv2w_colocated_t),      intent(inout) :: this
    type(grid_field_t),           intent(inout) :: u, v
    type(domain_t),               intent(in)    :: domain
    !output:
    type(grid_field_t),           intent(inout) :: uw, vw

    call uw%assign(1.0_8,u,domain%mesh_u)
    call vw%assign(1.0_8,v,domain%mesh_v)
end subroutine interp_uv2w_colocated

subroutine interp_uv2w_hor_colocated(this, uw, vw, u, v, domain)
    class(uv2w_hor_colocated_t),  intent(inout) :: this
    type(grid_field_t),           intent(inout) :: u, v
    type(domain_t),               intent(in)    :: domain
    !output:
    type(grid_field_t),           intent(inout) :: uw, vw

    call this%p2w_oper%apply(uw,u,domain)
    call this%p2w_oper%apply(vw,v,domain)
end subroutine interp_uv2w_hor_colocated

subroutine interp_uv2w_staggered(this, uw, vw, u, v, domain)
    class(uv2w_staggered_t),      intent(inout) :: this
    type(grid_field_t),           intent(inout) :: u, v
    type(domain_t),               intent(in)    :: domain
    !output:
    type(grid_field_t),           intent(inout) :: uw, vw

    call this%uv2p_oper%interp2d_vec2vec(this%up,this%vp,u,v,domain)
    call this%p2w_oper%apply(uw,this%up,domain)
    call this%p2w_oper%apply(vw,this%vp,domain)
end subroutine interp_uv2w_staggered

end module interpolators_uv2w_mod
