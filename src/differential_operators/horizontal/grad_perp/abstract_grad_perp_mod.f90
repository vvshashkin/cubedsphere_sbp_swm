module abstract_grad_perp_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t

implicit none

type, public, abstract :: grad_perp_operator_t

contains
    procedure(grad_perp_calc_i), deferred :: calc_grad_perp
end type grad_perp_operator_t

abstract interface
    subroutine grad_perp_calc_i(this, gu, gv, w, domain)
        import grad_perp_operator_t, grid_field_t, domain_t
        class(grad_perp_operator_t),  intent(inout) :: this
        type(domain_t),               intent(in)    :: domain
        type(grid_field_t),           intent(inout) :: w
        !output
        type(grid_field_t),     intent(inout) :: gu, gv
    end subroutine grad_perp_calc_i
end interface

contains

end module abstract_grad_perp_mod
