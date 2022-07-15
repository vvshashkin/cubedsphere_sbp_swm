module ts2_mod

use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use stvec_mod,                only : stvec_t
use stvec_swm_mod,            only : stvec_swm_t
use stvec_swm_factory_mod,    only : create_stvec_swm
use operator_mod,             only : operator_t
use operator_swm_factory_mod, only : create_swm_operator
use timescheme_mod,           only : timescheme_t
use timescheme_factory_mod,   only : create_timescheme
use outputer_abstract_mod,    only : outputer_t, outputer_vector_t
use outputer_factory_mod,     only : create_master_paneled_outputer,&
                                     create_latlon_outputer, create_latlon_vec_outputer
use parcomm_mod,              only : parcomm_global
use config_swm_mod,           only : config_swm_t

use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator

use config_ts2_mod, only : config_ts2_t

use test_fields_mod, only : solid_rotation_t, ts2_height_generator_t

use key_value_mod, only : key_value_r8_t

implicit none

type(config_ts2_t) :: config_ts2

type(ts2_height_generator_t) :: height_field
type(solid_rotation_t)       :: velocity_field

contains

subroutine run_ts2()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    type(config_swm_t) :: config
    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state, state_err, state_ex
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer
    class(outputer_vector_t), allocatable :: outputer_vec

    type(operator_swm_diff_t),   allocatable :: operator_diff
    class(timescheme_t),         allocatable :: timescheme_diff

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: namelist_string
    type(key_value_r8_t) :: diagnostics

    real(kind=8)     :: dt
    real(kind=8)     :: tau_write
    integer(kind=4)  :: nstep_write, nstep_diagnostics

    real(kind=8)    :: time
    real(kind=8)    :: l2_err_h, l2_ex_h, l2_err_u, l2_ex_u
    real(kind=8)    :: linf_err_h, linf_ex_h, linf_err_u, linf_ex_u
    integer(kind=4) :: it

    call read_namelist_as_str(namelist_string, "namelist_swm", parcomm_global%myid)

    call config%parse(namelist_string)

    config%config_domain%config_metric%omega = config_ts2%omega
    config%config_domain%config_metric%scale = config_ts2%a

    height_field = ts2_height_generator_t(h_mean = config_ts2%h_mean, &
                u0 = config_ts2%u0, omega = config_ts2%omega, a = config_ts2%a, &
                grav = config_ts2%grav, axis = config%config_domain%config_metric%rotation_axis)

    velocity_field = solid_rotation_t(u0 = config_ts2%u0, &
                                      axis = config%config_domain%config_metric%rotation_axis)

    dt = config%dt
    tau_write = config%tau_write
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(config%tau_diagnostics/dt)

    if (parcomm_global%myid==0) print*, "Advective CFL = ", &
            real(dt*config_ts2%u0/(2*pi/4/config%config_domain%N)/config_ts2%a,4)

    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(config_ts2%h_mean*config_ts2%grav)    &
            /(2*pi/4/config%config_domain%N)/config_ts2%a,4)

    call create_domain(domain, config%config_domain)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_ex,  domain, 0         , 0)
    call create_stvec_swm(state_err, domain, 0         , 0)

    call create_swm_operator(operator, config_ts2%grav, config, domain)
    call create_timescheme(timescheme, state, 'rk4')

    call create_swm_diff_operator(operator_diff, config, domain)
    call create_timescheme(timescheme_diff, state, config%diff_time_scheme)

    print*, 4*domain%partition%Nh, 2*domain%partition%Nh+1

    if(config%config_domain%staggering_type == "Ah") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", &
                                   config%v_components_type, domain)
    else if(config%config_domain%staggering_type == "C") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "A", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "C", &
                                       config%v_components_type, domain)
    else
        call parcomm_global%abort("This staggering is not implemented in"//&
                                  " swe test 2 output:"//&
                                  config%config_domain%staggering_type)
    end if

    call get_exact_solution(state,    domain, config%v_components_type)
    call get_exact_solution(state_ex, domain, config%v_components_type)

    select type(state_ex)
    class is (stvec_swm_t)
        l2_ex_h = l2norm(state_ex%h,  domain%mesh_p, domain%parcomm)
        l2_ex_u = sqrt(l2norm(state_ex%u,  domain%mesh_u, domain%parcomm)**2+&
                       l2norm(state_ex%v,  domain%mesh_v, domain%parcomm)**2)
        linf_ex_h = state_ex%h%maxabs(domain%mesh_p, domain%parcomm)
        linf_ex_u = max(state_ex%u%maxabs(domain%mesh_u, domain%parcomm), &
                        state_ex%v%maxabs(domain%mesh_v, domain%parcomm))
    end select

    select type(state)
    class is (stvec_swm_t)
        call outputer%write(state%h, domain, 'h.dat', 1)
        call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', 1)
    end select


    do it = 1, int(config%simulation_time/dt)

        call timescheme%step(state, operator, domain, dt)
        call timescheme_diff%step(state, operator_diff, domain, dt)

        time = it*dt

        call state_err%assign(1.0_8, state_ex, -1.0_8, state, domain)

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            if(parcomm_global%myid == 0) call diagnostics%print()

            select type(state_err)
            class is (stvec_swm_t)
                !call outputer%write(state_err%h, domain, 'h_err.dat', int(it/nstep_write))
                !call outputer_vec%write(state_err%u, state_err%v, domain, 'u_err.dat', 'v_err.dat', int(it/nstep_write))

                l2_err_h = l2norm(state_err%h,  domain%mesh_p, domain%parcomm) / l2_ex_h
                l2_err_u = sqrt(l2norm(state_err%u,  domain%mesh_u, domain%parcomm)**2+&
                                l2norm(state_err%v,  domain%mesh_v, domain%parcomm)**2) / l2_ex_u
                linf_err_h = state_err%h%maxabs(domain%mesh_p, domain%parcomm) / linf_ex_h
                linf_err_u = max(state_err%u%maxabs(domain%mesh_u, domain%parcomm), &
                                 state_err%v%maxabs(domain%mesh_v, domain%parcomm)) / linf_ex_u
                if (parcomm_global%myid==0) print '(A,F12.4,4(A,E15.7))', &
                                                  "Hours = ", real(time/3600 ,4), &
                                                  " l2_h =", real(l2_err_h,4), " linf_h = ", real(linf_err_h,4), &
                                                  " l2_u =", real(l2_err_u,4), " linf_u = ", real(linf_err_u,4)
            end select
        end if

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', int(it/nstep_write)+1)
            end select
        end if
    end do

end subroutine run_ts2

subroutine get_exact_solution(state, domain, v_components_type)

    use test_fields_mod, only : set_vector_test_field, &
                                set_scalar_test_field

    class(stvec_t), intent(inout) :: state
    type(domain_t), intent(in)    :: domain
    character(*),   intent(in)    :: v_components_type

    select type(state)
    class is (stvec_swm_t)

        call set_scalar_test_field(state%h, height_field, domain%mesh_p, 0)
        call set_vector_test_field(state%u, state%v, velocity_field, &
              domain%mesh_u, domain%mesh_v, 0, v_components_type)

    end select


end subroutine get_exact_solution

end module ts2_mod
