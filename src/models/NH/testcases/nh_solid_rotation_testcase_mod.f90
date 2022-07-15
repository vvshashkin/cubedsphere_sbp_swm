module NH_solid_rotation_testcase_mod

use parcomm_mod,        only : parcomm_global
use domain_mod,         only : domain_t
use stvec_mod,          only : stvec_t
use stvec_nh_mod,       only : stvec_nh_t
use test_fields_3d_mod, only : scalar_field3d_t

use solid_rotation_fields_factory_mod, only : create_solid_rotation_field_generators

implicit none

contains

subroutine get_nh_solid_rotation_initial_conditions(stvec, domain)
    class(stvec_t),   intent(inout) :: stvec
    type(domain_t),   intent(in)    :: domain

    select type(stvec)
    type is (stvec_nh_t)
        call get_solid_rotation_initial_conditions_nh_stvec(stvec, domain)
    class default
        call parcomm_global%abort("NH solid rotation testcases does not support this type of stvec")
    end select

end subroutine get_nh_solid_rotation_initial_conditions

subroutine get_solid_rotation_initial_conditions_nh_stvec(stvec, domain)

    use grid_field_mod,                only: grid_field_t
    use grid_field_factory_mod,        only: create_grid_field
    use test_fields_3d_mod,            only: scalar_field3d_t, vector_field3d_t
    use const_mod,                     only: Earth_grav

    class(stvec_nh_t),   intent(inout) :: stvec
    type(domain_t),      intent(in)    :: domain

    real(kind=8), parameter :: Nb = 0.01_8, T0 = 300.0_8, pref = 1e5_8, U0 = 20.0_8

    class(scalar_field3d_t), allocatable :: theta0_generator
    class(scalar_field3d_t), allocatable :: P0_generator
    class(vector_field3d_t), allocatable :: wind_gen

    call create_solid_rotation_field_generators("isoterm",u0=U0,omega=domain%mesh_p%omega, &
                                                sphere_rad = domain%mesh_p%scale, Nb = Nb, &
                                                grav = Earth_grav, alpha=0.0_8,      &
                                                theta_gen  = theta0_generator, &
                                                PExner_gen = P0_generator,  &
                                                wind_gen   = wind_gen)

    call wind_gen%get_vector_field(stvec%u,stvec%v,stvec%eta_dot, &
                                   domain%mesh_u,domain%mesh_v,domain%mesh_w,0,"contravariant")

    call P0_generator%get_scalar_field(stvec%P,domain%mesh_p,0)
    call theta0_generator%get_scalar_field(stvec%theta,domain%mesh_w,0)

end subroutine get_solid_rotation_initial_conditions_nh_stvec

end module NH_solid_rotation_testcase_mod
