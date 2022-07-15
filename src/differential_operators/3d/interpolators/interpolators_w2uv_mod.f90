module interpolators_w2uv_mod

use abstract_interpolators3d_mod,   only : interpolator_w2uv_t
use grid_field_mod,                 only : grid_field_t
use domain_mod,                     only : domain_t
use abstract_vertical_operator_mod, only : vertical_operator_t
use abstract_interpolators2d_mod,   only : interpolator2d_scalar2vec_t

implicit none

type, extends(interpolator_w2uv_t) :: w2uv_colocated_t
contains
    procedure :: interp_w2uv => interp_w2uv_colocated
end type w2uv_colocated_t

type, extends(interpolator_w2uv_t) :: w2uv_hor_colocated_t
    class(vertical_operator_t), allocatable :: w2p_oper
contains
    procedure :: interp_w2uv => interp_w2uv_hor_colocated
end type w2uv_hor_colocated_t

type, extends(interpolator_w2uv_t) :: w2uv_staggered_t
    class(vertical_operator_t), allocatable         :: w2p_oper
    type(grid_field_t)                              :: wp
    class(interpolator2d_scalar2vec_t), allocatable :: h2v_oper
contains
    procedure :: interp_w2uv => interp_w2uv_staggered
end type w2uv_staggered_t

contains

subroutine interp_w2uv_colocated(this, wu, wv, w, domain)
    class(w2uv_colocated_t),    intent(inout) :: this
    type(grid_field_t),         intent(inout) :: w
    type(domain_t),             intent(in)    :: domain
    !output:
    type(grid_field_t),         intent(inout) :: wu, wv

    call wu%assign(1.0_8,w,domain%mesh_u)
    call wv%assign(1.0_8,w,domain%mesh_v)
end subroutine interp_w2uv_colocated

subroutine interp_w2uv_hor_colocated(this, wu, wv, w, domain)
    class(w2uv_hor_colocated_t),  intent(inout) :: this
    type(grid_field_t),           intent(inout) :: w
    type(domain_t),               intent(in)    :: domain
    !output:
    type(grid_field_t),           intent(inout) :: wu, wv

    call this%w2p_oper%apply(wu,w,domain)
    call wv%assign(1.0_8,wu,domain%mesh_v)
end subroutine interp_w2uv_hor_colocated

subroutine interp_w2uv_staggered(this, wu, wv, w, domain)
    class(w2uv_staggered_t),      intent(inout) :: this
    type(grid_field_t),           intent(inout) :: w
    type(domain_t),               intent(in)    :: domain
    !output:
    type(grid_field_t),           intent(inout) :: wu, wv

    call this%w2p_oper%apply(this%wp,w,domain)
    call this%h2v_oper%interp2d_scalar2vec(wu,wv,this%wp,domain)
end subroutine interp_w2uv_staggered

end module interpolators_w2uv_mod
