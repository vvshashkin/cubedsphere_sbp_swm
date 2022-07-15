module shallow_atm_vecadv_mod

use abstract_vector_advection3d_mod, only : vector_advection3d_t
use abstract_vector_advection_mod,   only : vector_advection_operator_t
use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use abstract_interpolators3d_mod,    only : interpolator_w2uv_t
use domain_mod,                      only : domain_t
use grid_field_mod,                  only : grid_field_t
use abstract_adv_z_mod,              only : adv_z_t

type, extends(vector_advection3d_t) :: shallow_atm_staggered_vecadv_t
    class(vector_advection_operator_t), allocatable :: uv_hor_advection_op
    class(adv_z_t),                     allocatable :: uv_z_advec_op
    class(interpolator_w2uv_t),         allocatable :: interp_w2uv
    class(scalar_advection3d_t),        allocatable :: w_advection3d_op

    type(grid_field_t) :: eta_dot_u, eta_dot_v, u_tend_z, v_tend_z
    contains
    procedure calc_vec_adv3d => calc_vec_adv3d_staggered
end type shallow_atm_staggered_vecadv_t

contains

subroutine calc_vec_adv3d_staggered(this,u_tend,v_tend,w_tend,u,v,w,eta_dot,domain)
    class(shallow_atm_staggered_vecadv_t), intent(inout) :: this
    type(grid_field_t),                    intent(inout) :: u, v, w, eta_dot
    type(domain_t),                        intent(in)    :: domain
    !output
    type(grid_field_t),                    intent(inout) :: u_tend, v_tend, w_tend

    ! u&v vectical advection
    call this%interp_w2uv%interp_w2uv(this%eta_dot_u, this%eta_dot_v, eta_dot, domain)
    call this%uv_z_advec_op%calc_z_adv(this%u_tend_z,u,this%eta_dot_u,domain%mesh_u)
    call this%uv_z_advec_op%calc_z_adv(this%v_tend_z,v,this%eta_dot_v,domain%mesh_v)

    !u&v horizontal advection including metric terms tends
    call this%uv_hor_advection_op%calc_vec_advection_contra(u_tend,v_tend,u,v,domain)

    !total u&v advective tendency
    call u_tend%update(1.0_8,this%u_tend_z,domain%mesh_u)
    call v_tend%update(1.0_8,this%v_tend_z,domain%mesh_v)

    !w-advection
    call this%w_advection3d_op%calc_adv3d(w_tend,w,u,v,eta_dot,domain)
end subroutine calc_vec_adv3d_staggered

end module shallow_atm_vecadv_mod
