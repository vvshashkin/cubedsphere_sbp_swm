module div_3d_hor_vert_mod

use grid_field_mod,                 only : grid_field_t
use domain_mod,                     only : domain_t
use abstract_div_3d_mod,            only : div_3d_operator_t
use abstract_div_mod,               only : div_operator_t
use abstract_vertical_operator_mod, only : vertical_operator_t
use vec_math_mod,                   only : multiply_by_J, divide_by_J_self

implicit none

type, public, extends(div_3d_operator_t) :: div_3d_hor_vert_t
    class(div_operator_t),      allocatable :: div_uv_op
    class(vertical_operator_t), allocatable :: diff_eta_op
    type(grid_field_t) :: Jw, diff_eta
contains
    procedure, public :: calc_div
end type div_3d_hor_vert_t

contains

subroutine calc_div(this, div, u, v, w, domain)

    class(div_3d_hor_vert_t),  intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: u, v, w
    type(grid_field_t),        intent(inout) :: div

    call this%div_uv_op%calc_div(div, u, v, domain)

    call multiply_by_J(this%Jw, w, domain%mesh_w)

    call this%diff_eta_op%apply(this%diff_eta, this%Jw, domain)

    !use diff_eta to store diff_eta(Jw)/J
    call divide_by_J_self(this%diff_eta, domain%mesh_p)

    call div%update(1.0_8, this%diff_eta, domain%mesh_p)

end subroutine calc_div

end module div_3d_hor_vert_mod
