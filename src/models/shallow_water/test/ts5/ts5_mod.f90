module ts5_mod

use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use grid_field_mod,           only : grid_field_t
use grid_field_factory_mod,   only : create_grid_field
use stvec_mod,                only : stvec_t
use stvec_swm_mod,            only : stvec_swm_t
use stvec_swm_factory_mod,    only : create_stvec_swm
use operator_mod,             only : operator_t
use operator_swm_mod,         only : operator_swm_t
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

use config_ts5_mod, only : config_ts5_t

use test_fields_mod, only : solid_rotation_t, ts2_height_generator_t, &
                            ts5_orography_generator_t,                &
                            set_scalar_test_field, set_vector_test_field

use key_value_mod, only : key_value_r8_t

implicit none

type(config_ts5_t) :: config_ts5

type(ts2_height_generator_t)    :: height_field
type(ts5_orography_generator_t) :: orography_generator
type(solid_rotation_t)          :: velocity_field

type(grid_field_t)              :: h_orog, h_total

contains

subroutine run_ts5()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    type(config_swm_t) :: config
    type(domain_t)     :: domain
    type(grid_field_t) :: curl

    class(stvec_t),           allocatable :: state
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer, outputer_curl
    class(outputer_vector_t), allocatable :: outputer_vec

    type(operator_swm_diff_t),   allocatable :: operator_diff
    class(timescheme_t),         allocatable :: timescheme_diff

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: namelist_string
    type(key_value_r8_t) :: diagnostics

    real(kind=8)     :: dt
    real(kind=8)     :: tau_write
    integer(kind=4)  :: nstep_write, nstep_diagnostics

    real(kind=8)    :: time, l2err, l2_ex, l2u, l2v
    integer(kind=4) :: it

    call read_namelist_as_str(namelist_string, "namelist_swm", parcomm_global%myid)

    call config%parse(namelist_string)

    config%config_domain%config_metric%omega = config_ts5%omega
    config%config_domain%config_metric%scale = config_ts5%a

    height_field = ts2_height_generator_t(h_mean = config_ts5%h_mean, &
                u0 = config_ts5%u0, omega = config_ts5%omega, a = config_ts5%a, &
                grav = config_ts5%grav)

    orography_generator = ts5_orography_generator_t(r_mount=config_ts5%r_mount,   &
                          h_mount=config_ts5%h_mount, x_mount=config_ts5%x_mount, &
                          y_mount=config_ts5%y_mount, z_mount=config_ts5%z_mount)

    velocity_field = solid_rotation_t(u0 = config_ts5%u0)

    dt = config%dt
    tau_write = config%tau_write
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(config%tau_diagnostics/dt)

    if (parcomm_global%myid==0) print*, "Advective CFL = ", &
            real(dt*config_ts5%u0/(2*pi/4/config%config_domain%N)/config_ts5%a,4)

    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(config_ts5%h_mean*config_ts5%grav)    &
            /(2*pi/4/config%config_domain%N)/config_ts5%a,4)

    call create_domain(domain, config%config_domain)

    call create_grid_field(h_orog,0,0,domain%mesh_p)
    call create_grid_field(h_total,0,0,domain%mesh_p)

    call create_stvec_swm(state,     domain, halo_width, 0)

    call create_swm_operator(operator, config_ts5%grav, config, domain, orography_generator)
    call create_timescheme(timescheme, state, 'rk4')

    call create_swm_diff_operator(operator_diff, config, domain)
    call create_timescheme(timescheme_diff, state, config%diff_time_scheme)

    print*, 4*domain%partition%Nh, 2*domain%partition%Nh+1

    if(config%config_domain%staggering_type == "Ah") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_outputer(outputer_curl,     2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", &
                                   config%v_components_type, domain)
    else if(config%config_domain%staggering_type == "C") then
        call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "A", domain)
        call create_latlon_outputer(outputer_curl,     2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "C", &
                                       config%v_components_type, domain)
    else
        call parcomm_global%abort("This staggering is not implemented in"//&
                                  " swe test 2 output:"//&
                                  config%config_domain%staggering_type)
    end if

    call create_grid_field(curl, halo_width, 0, domain%mesh_q)

    select type(state)
    class is (stvec_swm_t)

        call set_scalar_test_field(state%h, height_field, domain%mesh_p, 0)
        call set_scalar_test_field(h_orog, orography_generator, domain%mesh_p, 0)
        call state%h%update(-1.0_8,h_orog,domain%mesh_p)

        call set_vector_test_field(state%u, state%v, velocity_field, &
                                   domain%mesh_u, domain%mesh_v, 0, config%v_components_type)

        call h_total%assign(1.0_8,state%h,1.0_8,h_orog,domain%mesh_p)
        call outputer%write(h_total, domain, 'h.dat', 1)
        call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', 1)
        select type(operator)
        class is (operator_swm_t)
            call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
            call outputer_curl%write(curl, domain, 'curl.dat', 1)
        end select
    end select


    do it = 1, int(config%simulation_time/dt)

        call timescheme%step(state, operator, domain, dt)
        select type(state)
        class is (stvec_swm_t)
            call state%h%assign(1.0_8,state%h,1.0_8,h_orog,domain%mesh_p)
            call timescheme_diff%step(state, operator_diff, domain, dt)
            call state%h%assign(1.0_8,state%h,-1.0_8,h_orog,domain%mesh_p)
        end select

        time = it*dt

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            if(parcomm_global%myid == 0) call diagnostics%print()
        end if

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)
                call h_total%assign(1.0_8,state%h,1.0_8,h_orog,domain%mesh_p)
                call outputer%write(h_total, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u, state%v, domain, 'u.dat', 'v.dat', int(it/nstep_write)+1)
                select type(operator)
                class is (operator_swm_t)
                    call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
                    call outputer_curl%write(curl, domain, 'curl.dat', int(it/nstep_write)+1)
                end select
            end select

            if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
                                                "irec=",int(it/nstep_write)+1

        end if
    end do

end subroutine run_ts5

end module ts5_mod
