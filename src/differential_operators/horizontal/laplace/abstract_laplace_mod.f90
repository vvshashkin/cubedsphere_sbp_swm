module abstract_laplace_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

type, abstract, public :: laplace_operator_t
    contains
        procedure(laplace_calc_procedure), deferred :: calc_laplace
end type laplace_operator_t

abstract interface
    subroutine laplace_calc_procedure(this, f1, f, domain)
        import laplace_operator_t, grid_field_t, domain_t
        class(laplace_operator_t), intent(inout) :: this
        type(domain_t),            intent(in)    :: domain
        type(grid_field_t),        intent(inout) :: f
        !output:
        type(grid_field_t),        intent(inout) :: f1
    end subroutine laplace_calc_procedure
end interface

end module abstract_laplace_mod
