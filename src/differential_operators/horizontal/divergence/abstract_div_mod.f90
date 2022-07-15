module abstract_div_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t

implicit none

type, abstract, public :: div_operator_t

contains

procedure(div_calc_procedure), deferred :: calc_div

end type div_operator_t

abstract interface
    subroutine div_calc_procedure(this, div, u, v, domain)
        import div_operator_t, grid_field_t, domain_t
        class(div_operator_t),  intent(inout) :: this
        type(domain_t),         intent(in)    :: domain
        type(grid_field_t),     intent(inout) :: u, v
        !out put
        type(grid_field_t),     intent(inout) :: div
    end subroutine div_calc_procedure
end interface

    contains

end module abstract_div_mod
