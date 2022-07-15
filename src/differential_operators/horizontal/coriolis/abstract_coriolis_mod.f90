module abstract_coriolis_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t
use parcomm_mod,       only : parcomm_global

implicit none

type, abstract, public :: coriolis_operator_t
contains
    procedure(calc_coriolis_i),         deferred :: calc_coriolis
    procedure :: calc_coriolis_contra
    procedure(calc_coriolis_vec_inv_i), deferred :: calc_coriolis_vec_inv
end type coriolis_operator_t

abstract interface
    subroutine calc_coriolis_i(this, cor_u, cor_v, ut, vt, domain)
        import coriolis_operator_t, grid_field_t, domain_t
        class(coriolis_operator_t), intent(inout) :: this
        type(domain_t),             intent(in)    :: domain
        type(grid_field_t),         intent(inout) :: ut, vt!contravariant components
        type(grid_field_t),         intent(inout) :: cor_u, cor_v!covariant components
    end subroutine calc_coriolis_i
    subroutine calc_coriolis_vec_inv_i(this, cor_u, cor_v, hu, hv, h, curl, domain)
        import coriolis_operator_t, grid_field_t, domain_t
        class(coriolis_operator_t), intent(inout) :: this
        type(domain_t),             intent(in)    :: domain
        type(grid_field_t),         intent(inout) :: hu, hv!massflux contravariant components
        type(grid_field_t),         intent(inout) :: h, curl
        type(grid_field_t),         intent(inout) :: cor_u, cor_v!covariant components
    end subroutine calc_coriolis_vec_inv_i
end interface

contains

subroutine calc_coriolis_contra(this, cor_ut, cor_vt, ut, vt, domain)
    class(coriolis_operator_t), intent(inout) :: this
    type(domain_t),             intent(in)    :: domain
    type(grid_field_t),         intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),         intent(inout) :: cor_ut, cor_vt!contravariant components

    call parcomm_global%abort("Calc_coriolis contra not implemented")
end subroutine calc_coriolis_contra

end module abstract_coriolis_mod
