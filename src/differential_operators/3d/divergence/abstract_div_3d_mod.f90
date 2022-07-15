module abstract_div_3d_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t

implicit none

type, abstract, public :: div_3d_operator_t
contains
    procedure(calc_div_i), deferred :: calc_div
end type div_3d_operator_t

abstract interface
    subroutine calc_div_i(this, div, u, v, w, domain)
        import div_3d_operator_t, grid_field_t, domain_t
        class(div_3d_operator_t),  intent(inout) :: this
        type(domain_t),            intent(in)    :: domain
        type(grid_field_t),        intent(inout) :: u, v, w
        !output
        type(grid_field_t),        intent(inout) :: div
    end subroutine calc_div_i
end interface

contains

end module abstract_div_3d_mod
