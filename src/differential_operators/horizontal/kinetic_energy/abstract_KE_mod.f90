module abstract_KE_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

type, public, abstract :: KE_operator_t
contains
    procedure(calc_KE_i), public, deferred :: calc_KE
end type

abstract interface
    subroutine calc_KE_i(this, KE, u, v, ut, vt, domain)
        import KE_operator_t, grid_field_t, domain_t
        class(KE_operator_t), intent(inout)  :: this
        type(domain_t),       intent(in)    :: domain
        type(grid_field_t),   intent(inout) :: u, v, ut, vt
        type(grid_field_t),   intent(inout) :: KE
    end subroutine calc_KE_i
end interface

contains

end module abstract_KE_mod
