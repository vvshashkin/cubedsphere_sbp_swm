module identity_vertical_operator_mod

use domain_mod,                     only : domain_t
use grid_field_mod,                 only : grid_field_t
use abstract_vertical_operator_mod, only : vertical_operator_t
use parcomm_mod,                    only : parcomm_global


implicit none

type, extends(vertical_operator_t) :: identity_vertical_operator_t
    character(len=1) :: points_type
contains
    procedure, public :: apply
end type identity_vertical_operator_t

contains


subroutine apply(this, f_out, f_in, domain)
    class(identity_vertical_operator_t), intent(in)    :: this
    type(grid_field_t),                  intent(inout) :: f_out
    type(grid_field_t),                  intent(in)    :: f_in
    type(domain_t),                      intent(in)    :: domain


    select case(this%points_type)
    case("p")
        call f_out%assign(1.0_8, f_in, domain%mesh_p)
    case default
        call parcomm_global%abort("Unknown points type in identity_vertical_operator_t"//&
                                  this%points_type)
    end select

end subroutine apply

end module identity_vertical_operator_mod
