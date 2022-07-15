module grad_3d_hor_vert_mod

use grid_field_mod,                 only : grid_field_t
use domain_mod,                     only : domain_t
use abstract_grad_3d_mod,           only : grad_3d_operator_t
use abstract_grad_mod,              only : grad_operator_t
use abstract_vertical_operator_mod, only : vertical_operator_t

implicit none

type, public, extends(grad_3d_operator_t) :: grad_3d_hor_vert_t
    class(grad_operator_t),     allocatable :: grad_xy  !horizontal part of grad
    class(vertical_operator_t), allocatable :: diff_eta !vertical   part of grad
contains
    procedure, public:: calc_grad
end type grad_3d_hor_vert_t

contains

subroutine calc_grad(this, gx, gy, gz, f, domain)

    class(grad_3d_hor_vert_t), intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f
    type(grid_field_t),        intent(inout) :: gx, gy, gz

    call this%grad_xy%calc_grad(gx, gy, f, domain)
    call this%diff_eta%apply(gz, f, domain)

end subroutine calc_grad

end module grad_3d_hor_vert_mod
