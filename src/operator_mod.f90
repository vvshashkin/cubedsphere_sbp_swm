module operator_mod

use stvec_mod,   only : stvec_t
use domain_mod,  only : domain_t
use parcomm_mod, only : parcomm_global

implicit none

private

type, abstract, public :: operator_t
contains
    procedure(apply_i),  public, deferred :: apply !vout=A*vin
    procedure,           public           :: solve !vout=inverse(I-dt*A)*rhs
    procedure,           public           :: get_diagnostics
end type operator_t

abstract interface
    subroutine apply_i(this, vout, vin, domain)
        import stvec_t, operator_t, domain_t
        class(operator_t), intent(inout) :: this
        class(stvec_t),    intent(inout) :: vout !inout to enable preallocated vectors
        class(stvec_t),    intent(inout) :: vin
        type(domain_t),    intent(in)    :: domain
    end subroutine apply_i
end interface

contains

subroutine solve(this, vout, rhs, dt, domain)
    class(operator_t), intent(inout) :: this
    class(stvec_t),    intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),    intent(inout) :: rhs
    real(kind=8),      intent(in)    :: dt
    type(domain_t),    intent(in)    :: domain

    call parcomm_global%abort("Solve function not implemented for specific operator class")
end subroutine solve

function get_diagnostics(this, v, domain) result(diagnostics)

    use key_value_mod, only : key_value_r8_t

    class(operator_t),     intent(inout) :: this
    class(stvec_t),        intent(inout) :: v
    type(domain_t),        intent(in)    :: domain

    type(key_value_r8_t)  :: diagnostics

    diagnostics = key_value_r8_t()

    call parcomm_global%abort("get_diagnostics function not implemented for specific operator class")
end function get_diagnostics

function get_diagnostics_tend(this, v, vtend, domain) result(diagnostics)

    use key_value_mod, only : key_value_r8_t

    class(operator_t),     intent(inout) :: this
    class(stvec_t),        intent(inout) :: v, vtend
    type(domain_t),        intent(in)    :: domain

    type(key_value_r8_t)  :: diagnostics

    diagnostics = key_value_r8_t()

    call parcomm_global%abort("get_diagnostics_tend function not implemented for specific operator class")
end function get_diagnostics_tend

end module operator_mod
