module abstract_grad_3d_mod

use grid_field_mod,                 only : grid_field_t
use domain_mod,                     only : domain_t

implicit none

type, public :: grad_3d_operator_t
contains
    procedure, public:: calc_grad
end type grad_3d_operator_t

abstract interface
    subroutine calc_grad(this, gx, gy, gz, f, domain)
        import grad_3d_operator_t, grid_field_t, domain_t
        class(grad_3d_operator_t), intent(inout) :: this
        type(domain_t),            intent(in)    :: domain
        type(grid_field_t),        intent(inout) :: f
        !output:
        type(grid_field_t),     intent(inout) :: gx, gy, gz
    end subroutine calc_grad
end interface

contains

end module abstract_grad_3d_mod
