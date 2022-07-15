module advection3d_oper_mod

use operator_mod,                    only : operator_t
use grid_field_mod,                  only : grid_field_t
use stvec_mod,                       only : stvec_t
use stvec_nh_mod,                    only : stvec_nh_t
use parcomm_mod,                     only : parcomm_global
use domain_mod,                      only : domain_t
use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use test_fields_3d_mod,              only : vector_field3d_t, non_stationary_vector_field3d_t

implicit none

type, extends(operator_t) :: advection3d_operator_t

    class(scalar_advection3d_t), allocatable :: p_adv_oper, theta_adv_oper
    class(vector_field3d_t),     allocatable :: wind_generator
    type(grid_field_t) :: u_adv, v_adv, eta_dot_adv
    contains
    procedure, public :: apply

end type advection3d_operator_t

contains

subroutine apply(this, vout, vin, domain)

    class(advection3d_operator_t), intent(inout) :: this
    class(stvec_t),                intent(inout) :: vout
    class(stvec_t),                intent(inout) :: vin
    type(domain_t),                intent(in)    :: domain

    select type(vout)
    class is (stvec_nh_t)
    select type(vin)
    class is (stvec_nh_t)
        select type(wind_generator=>this%wind_generator)
        class is (non_stationary_vector_field3d_t)
            wind_generator%t = vin%model_time
        end select
        call this%wind_generator%get_vector_field(this%u_adv, this%v_adv, this%eta_dot_adv, &
                                                  domain%mesh_u, domain%mesh_v, domain%mesh_w, &
                                                  0,"contravariant")

        call vout%u%assign(0.0_8, domain%mesh_u)
        call vout%v%assign(0.0_8, domain%mesh_v)
        call vout%eta_dot%assign(0.0_8, domain%mesh_w)
        !call vout%theta%assign(0.0_8, domain%mesh_w)
        call this%theta_adv_oper%calc_adv3d(vout%theta,vin%theta,this%u_adv,this%v_adv, &
                                            this%eta_dot_adv,domain)
        call this%p_adv_oper%calc_adv3d(vout%P,vin%P,this%u_adv,this%v_adv, &
                                        this%eta_dot_adv,domain)
        vout%model_time = 1.0_8 != dt / dt
    class default
        call parcomm_global%abort("Ptheta_linear_nh_oper_t, vin type error")
    end select
    class default
        call parcomm_global%abort("Ptheta_linear_nh_oper_t, vout type error")
    end select


end subroutine apply

end module advection3d_oper_mod
