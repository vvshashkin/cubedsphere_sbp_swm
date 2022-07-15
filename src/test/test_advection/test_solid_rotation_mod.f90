module test_solid_rotation_mod

use domain_mod,                     only : domain_t
use domain_factory_mod,             only : create_domain
use stvec_mod,                      only : stvec_t
use stvec_advection_mod,            only : stvec_advection_t
use stvec_advection_factory_mod,    only : create_stvec_advection
use operator_mod,                   only : operator_t
use operator_advection_factory_mod, only : create_advection_operator
use timescheme_mod,                 only : timescheme_t
use timescheme_factory_mod,         only : create_timescheme
use outputer_abstract_mod,          only : outputer_t
use outputer_factory_mod,           only : create_master_paneled_outputer,&
                                           create_latlon_outputer
use parcomm_mod,                    only : parcomm_global

use test_fields_mod, only : solid_rotation_t, solid_rotated_scalar_field_t
implicit none

type(solid_rotation_t),             allocatable :: velocity_field
type(solid_rotated_scalar_field_t), allocatable :: advected_field

contains

subroutine test_solid_rotation()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    type(domain_t) :: domain
    integer(kind=4), parameter :: N = 32, nz = 1, halo_width = 10

    real(kind=8),    parameter :: rotation_period = 1.0_8, rotation_axis_angle = pi/8
    integer(kind=4), parameter :: N_periods=2

    real(kind=8), parameter    :: dt = 0.005_8
    real(kind=8), parameter    :: tau_write = 0.01_8
    integer(kind=4), parameter :: nstep_write = nint(tau_write/dt)

    class(stvec_t),      allocatable :: state, state_err
    class(operator_t),   allocatable :: operator
    class(timescheme_t), allocatable :: timescheme
    class(outputer_t),   allocatable :: outputer

    real(kind=8)    :: time, l2err
    integer(kind=4) :: it

    call create_domain(domain, "cube", "C", N, nz)

    !WORKAROUND
    if (parcomm_global%myid==0) print*, "CFL = ", real(dt*2*pi/rotation_period/(2*pi/4/N),4)

    call create_stvec_advection(state,     domain, halo_width, 0)
    call create_stvec_advection(state_err, domain,          0, 0)

    !Ah
    !call create_advection_operator(operator, "massflux_colocated", "divergence_ah42_sbp", domain)
    !C
    call create_advection_operator(operator, "massflux_c_up4", "divergence_c_sbp42", domain)
    !call create_advection_operator(operator, "massflux_c_sbp42", "divergence_c_sbp42", domain)

    call create_timescheme(timescheme, state, 'rk4')

    !call create_master_paneled_outputer(outputer, "p", domain)
    call create_latlon_outputer(outputer, 2*N+1, 4*N, "A", domain)

    call init_field_generators(rotation_period, rotation_axis_angle)

    call get_exact_solution(state, domain, 0.0_8)

    do it = 1, nint(N_periods*rotation_period/dt)

        call timescheme%step(state, operator, domain, dt)

        time = it*dt

        call get_exact_solution(state_err, domain, time)
        call state_err%update(-1.0_8, state, domain)

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_advection_t)
                call outputer%write(state%h, domain, 'h_adv.dat', int(it/nstep_write))
            end select
            select type(state_err)
            class is (stvec_advection_t)
                call outputer%write(state_err%h, domain, 'h_adv_err.dat', int(it/nstep_write))
                l2err = l2norm(state_err%h, domain%mesh_p, domain%parcomm)
                if (parcomm_global%myid==0) print*, "Periods = ", real(time ,4), &
                                                    "L2err =", real(l2err,4)
            end select
        end if

    end do


end subroutine test_solid_rotation

subroutine init_field_generators(period, rotation_axis_angle)

    use test_fields_mod, only : gaussian_hill_scalar_field_generator_t, &
                                scalar_field_generator_t
    use const_mod,       only : pi

    real(kind=8) :: period, rotation_axis_angle

    class(scalar_field_generator_t), allocatable :: tmp

    velocity_field = solid_rotation_t(alpha = rotation_axis_angle, u0 = 2*pi/period, &
                                      axis = [-sin(rotation_axis_angle), 0.0_8, cos(rotation_axis_angle)])

    tmp = gaussian_hill_scalar_field_generator_t()

    advected_field = solid_rotated_scalar_field_t( &
                   solid_rotation_vector_field = velocity_field, &
                   scalar_field = tmp)

end subroutine init_field_generators

subroutine get_exact_solution(state, domain, time)

    use test_fields_mod, only : set_vector_test_field, &
                                set_scalar_test_field

    class(stvec_t), intent(inout) :: state
    type(domain_t), intent(in)    :: domain
    real(kind=8),   intent(in)    :: time

    advected_field%time = time

    select type(state)
    class is (stvec_advection_t)

        call set_scalar_test_field(state%h, advected_field, domain%mesh_p, 0)
        call set_vector_test_field(state%u, state%v, velocity_field, &
              domain%mesh_u, domain%mesh_v, 0, "contravariant")

    end select


end subroutine get_exact_solution


end module test_solid_rotation_mod
