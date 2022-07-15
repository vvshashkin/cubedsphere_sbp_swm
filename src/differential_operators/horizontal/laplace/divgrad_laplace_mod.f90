module divgrad_laplace_mod

use abstract_laplace_mod,   only: laplace_operator_t
use grid_field_mod,         only : grid_field_t
use domain_mod,             only : domain_t
use abstract_grad_mod,      only : grad_operator_t
use abstract_div_mod,       only : div_operator_t
use abstract_co2contra_mod, only : co2contra_operator_t

implicit none

type, extends(laplace_operator_t) :: divgrad_laplace_t
    class(grad_operator_t),      allocatable :: grad_operator
    class(co2contra_operator_t), allocatable :: co2contra_operator
    class(div_operator_t),       allocatable :: div_operator
    type(grid_field_t)                       :: gx, gy, gxt, gyt
    contains
    procedure :: calc_laplace
end type divgrad_laplace_t

contains

subroutine calc_laplace(this, f1, f, domain)
    class(divgrad_laplace_t),  intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f
    !output:
    type(grid_field_t),        intent(inout) :: f1

    call this%grad_operator%calc_grad(this%gx, this%gy, f, domain)
    call this%co2contra_operator%transform(this%gxt, this%gyt, this%gx, this%gy, domain)
    call this%div_operator%calc_div(f1, this%gxt, this%gyt, domain)
end subroutine calc_laplace

end module divgrad_laplace_mod
