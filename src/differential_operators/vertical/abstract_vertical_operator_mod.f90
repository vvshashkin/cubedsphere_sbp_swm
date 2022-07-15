module abstract_vertical_operator_mod

use domain_mod,     only : domain_t
use grid_field_mod, only : grid_field_t

implicit none

type, abstract :: vertical_operator_t
    contains
    procedure(apply_vert_operator), deferred :: apply
end type vertical_operator_t

abstract interface
    subroutine apply_vert_operator(this, f_out, f_in, domain)
        import vertical_operator_t, grid_field_t, domain_t
        class(vertical_operator_t), intent(in)    :: this
        type(grid_field_t),         intent(inout) :: f_out
        type(grid_field_t),         intent(in)    :: f_in
        type(domain_t),             intent(in)    :: domain
    end subroutine apply_vert_operator
end interface

contains

end module abstract_vertical_operator_mod
