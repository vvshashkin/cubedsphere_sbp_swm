module abstract_co2contra_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t
use parcomm_mod,    only : parcomm_global

implicit none

type, abstract, public :: co2contra_operator_t

contains

procedure(transform_procedure), deferred :: transform
procedure, public  :: transform2co

end type co2contra_operator_t

abstract interface
    subroutine transform_procedure(this, u_contra, v_contra, u_cov, v_cov, domain)
        import co2contra_operator_t, grid_field_t, domain_t
        class(co2contra_operator_t), intent(inout) :: this
        type(domain_t),              intent(in)    :: domain
        type(grid_field_t),          intent(inout) :: u_cov, v_cov
        !output:
        type(grid_field_t),          intent(inout) :: u_contra, v_contra
    end subroutine transform_procedure
end interface

contains

subroutine transform2co(this, u_cov, v_cov, u_contra, v_contra, domain)
    class(co2contra_operator_t), intent(inout) :: this
    type(domain_t),              intent(in)    :: domain
    type(grid_field_t),          intent(inout) :: u_contra, v_contra
    !output:
    type(grid_field_t),          intent(inout) :: u_cov, v_cov

    call parcomm_global%abort("contra2co not implemented")
end subroutine transform2co

end module abstract_co2contra_mod
