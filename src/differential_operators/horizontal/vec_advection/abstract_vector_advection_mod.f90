module abstract_vector_advection_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t

implicit none

type, abstract, public :: vector_advection_operator_t
contains
    procedure(calc_vec_advection_i),         deferred :: calc_vec_advection
    procedure(calc_vec_advection_contra_i),  deferred :: calc_vec_advection_contra
end type vector_advection_operator_t

abstract interface
    subroutine calc_vec_advection_i(this, u_tend, v_tend, u, v, ut, vt, domain)
        import vector_advection_operator_t, grid_field_t, domain_t

        class(vector_advection_operator_t), intent(inout) :: this
        type(domain_t),                     intent(in)    :: domain
        type(grid_field_t),                 intent(inout) :: u,  v!covariant components
        type(grid_field_t),                 intent(inout) :: ut, vt!contravariant components
        type(grid_field_t),                 intent(inout) :: u_tend, v_tend!advective tendencies
    end subroutine calc_vec_advection_i
    subroutine calc_vec_advection_contra_i(this, u_tend, v_tend, ut, vt, domain)
        import vector_advection_operator_t, grid_field_t, domain_t

        class(vector_advection_operator_t), intent(inout) :: this
        type(domain_t),                     intent(in)    :: domain
        type(grid_field_t),                 intent(inout) :: ut, vt!contravariant components
        type(grid_field_t),                 intent(inout) :: u_tend, v_tend!advective tendencies
    end subroutine calc_vec_advection_contra_i
end interface

contains

end module abstract_vector_advection_mod
