module abstract_vector_advection3d_mod

use domain_mod,      only : domain_t
use grid_field_mod,  only : grid_field_t

implicit none

type, abstract :: vector_advection3d_t
contains
    procedure(calc_vec_adv3d), deferred :: calc_vec_adv3d
end type vector_advection3d_t

abstract interface
    subroutine calc_vec_adv3d(this,u_tend,v_tend,w_tend,u,v,w,eta_dot,domain)
        import vector_advection3d_t, grid_field_t, domain_t
        class(vector_advection3d_t), intent(inout) :: this
        type(grid_field_t),          intent(inout) :: u, v, w, eta_dot
        type(domain_t),              intent(in)    :: domain
        !output
        type(grid_field_t),          intent(inout) :: u_tend, v_tend, w_tend
    end subroutine calc_vec_adv3d
end interface

end module abstract_vector_advection3d_mod
