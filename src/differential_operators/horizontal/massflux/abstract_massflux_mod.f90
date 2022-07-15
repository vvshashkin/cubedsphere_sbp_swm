module abstract_massflux_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

type, abstract, public :: massflux_operator_t

contains

procedure(massflux_calc_procedure), deferred :: calc_massflux

end type massflux_operator_t

abstract interface
    subroutine massflux_calc_procedure(this, fx, fy, f, u, v, domain)
        import massflux_operator_t, grid_field_t, domain_t
        class(massflux_operator_t), intent(inout) :: this
        type(domain_t),             intent(in)    :: domain
        type(grid_field_t),         intent(inout) :: f, u, v
        !output:
        type(grid_field_t),         intent(inout) :: fx, fy
    end subroutine massflux_calc_procedure
end interface

    contains

end module abstract_massflux_mod
