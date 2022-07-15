module explicit_Eul1_mod

use stvec_mod,      only : stvec_t
use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t
use domain_mod,     only : domain_t

implicit none

private

type, public, extends(timescheme_t) :: explicit_Eul1_t
    class(stvec_t), allocatable :: tendency
contains
    procedure, public :: step => step_explicit_Eul1
end type explicit_Eul1_t

contains

subroutine step_explicit_Eul1(this, v0, operator, domain, dt)

    class(explicit_Eul1_t), intent(inout) :: this
    class(stvec_t),         intent(inout) :: v0
    class(operator_t),      intent(inout) :: operator
    type(domain_t),         intent(in)    :: domain
    real(kind=8),           intent(in)    :: dt

    call operator%apply(this%tendency, v0, domain)
    call v0%update(dt, this%tendency, domain)
end subroutine step_explicit_Eul1

end module explicit_Eul1_mod
