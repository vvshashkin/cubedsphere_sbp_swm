module forcing_test_mod

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

use operator_swm_mod,              only : operator_swm_t
use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator

use const_mod,  only : Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

use random_friction_mod, only : random_friction_t, initialize_random_friction, &
                                random_scalar_t,   initialize_random_scalar
use key_value_mod,   only : key_value_r8_t

implicit none

real(kind=8) :: grav    = Earth_grav
real(kind=8) :: omega   = Earth_omega
real(kind=8) :: a       = Earth_radii
real(kind=8) :: h_mean  = 3.5e3_8
real(kind=8) :: tau_forcing_h = 86400._8*100._8
! real(kind=8) :: tau_fric_decorrel = 1e4_8

contains

subroutine run_forcing_test()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    type(config_swm_t) :: config
    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state, state_eq, state_buff
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer
    class(outputer_vector_t), allocatable :: outputer_vec

    type(operator_swm_diff_t),   allocatable :: operator_diff
    class(timescheme_t),         allocatable :: timescheme_diff
    type(random_friction_t) :: random_friction
    type(random_scalar_t)   :: random_scalar

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: namelist_string
    type(key_value_r8_t) :: diagnostics

    real(kind=8)      :: dt
    real(kind=8)      :: tau_write
    integer(kind=4)   :: nstep_write, nstep_diagnostics

    real(kind=8)    :: time, l2err, l2_ex
    integer(kind=4) :: it

    call read_namelist_as_str(namelist_string, "namelist_swm", parcomm_global%myid)

    call config%parse(namelist_string)

    config%config_domain%config_metric%omega = omega
    config%config_domain%config_metric%scale = a

    dt = config%dt
    tau_write = config%tau_write
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(config%tau_diagnostics/dt)

    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(h_mean*grav)    &
            /(2*pi/4/config%config_domain%N)/a,4)

    call create_domain(domain, config%config_domain)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_eq,  domain, halo_width, 0)
    call create_stvec_swm(state_buff,domain, halo_width, 0)

    call create_swm_operator(operator, grav, config, domain)
    call create_timescheme(timescheme, state, 'rk4')

    call create_swm_diff_operator(operator_diff, config, domain)
    call create_timescheme(timescheme_diff, state, config%diff_time_scheme)

    ! call initialize_random_friction(random_friction, l=6, tau=1e4_8, &
    !                                 sigma_amp=3e-5_8, domain=domain)

    call initialize_random_scalar(random_scalar, l=6, tau=1e5_8, &
                                  amp=500._8, domain=domain)

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
                                  " barotropic instability swe test output:"//&
                                  config%config_domain%staggering_type)
    end if

    select type(state_eq)
    class is (stvec_swm_t)
        call init_equilibrium(state_eq, domain)
        ! call set_scalar_test_field(state_eq%h,forcing_test_height_generator,domain%mesh_p,0)
        ! call set_vector_test_field(state_eq%u,state_eq%v,forcing_test_wind_generator, &
        !                            domain%mesh_u,domain%mesh_v,0,config%v_components_type)
    end select

    select type(state)
    class is (stvec_swm_t)
        call state%u%assign(0.0_8,domain%mesh_u)
        call state%v%assign(0.0_8,domain%mesh_v)
        call state%h%assign(h_mean,domain%mesh_p)
        ! call state%assign(1.0_8,state_eq,0.0_8,state_eq,domain)

        call outputer%write(state%h, domain, 'h.dat',1)
        call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",1)
    end select
    select type(state_buff)
    class is (stvec_swm_t)
        call state_buff%u%assign(0.0_8,domain%mesh_u)
        call state_buff%v%assign(0.0_8,domain%mesh_v)
        call state_buff%h%assign(h_mean,domain%mesh_p)
    end select

    do it = 1, int(config%simulation_time/dt)

        select type(operator)
        type is (operator_swm_t)
            call random_scalar%apply_update_forcing(operator%h_surf,domain,dt)
            if(mod(it,nstep_write) == 0) then
            	call outputer%write(operator%h_surf, domain, 'h_surf.dat', it/nstep_write)
            end if
        end select

        call timescheme%step(state, operator, domain, dt)
        call timescheme_diff%step(state, operator_diff, domain, dt)

        ! call state%update(dt/tau_forcing, state_eq, -dt/tau_forcing, state, domain)
        select type(state_eq)
        class is (stvec_swm_t)
        select type(state)
        class is(stvec_swm_t)
            ! call random_friction%apply_update_forcing(state_eq%u,state_eq%v, &
            !                                           state%u,state%v,domain,dt)
            ! call state%u%update(dt,state_eq%u,domain%mesh_u)
            ! call state%v%update(dt,state_eq%v,domain%mesh_v)
            call state%h%update(dt/tau_forcing_h, state_eq%h, &
                               -dt/tau_forcing_h, state%h, domain%mesh_p)
            ! call outputer_vec%write(state%u,state%v,domain,"ru.dat","rv.dat",it)

        end select
        end select

        time = it*dt

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            if(domain%parcomm%myid==0) call diagnostics%print()
        end if

        call state_buff%update(1.0_8/nstep_write,state,domain) 
        !if(mod(it,nstep_write) == 0) then
        !    select type(state_buff)
        !    class is (stvec_swm_t)
        !        call outputer%write(state_buff%h, domain, 'h.dat', int(it/nstep_write)+1)
        !        call outputer_vec%write(state_buff%u,state_buff%v,domain,"u.dat","v.dat",int(it/nstep_write)+1)
        !        if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
        !                                            "rec=",int(it/nstep_write)+1
        !        call state_buff%h%assign(0.0_8,domain%mesh_p)
        !        call state_buff%u%assign(0.0_8,domain%mesh_u)
        !        call state_buff%v%assign(0.0_8,domain%mesh_v)
        !    end select
        !end if
        if(mod(it,nstep_write) == 0) then
            select type(state)
            class is (stvec_swm_t)
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",int(it/nstep_write)+1)
                if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
                                                    "rec=",int(it/nstep_write)+1
            end select
        end if
    end do

end subroutine run_forcing_test

subroutine init_equilibrium(state,domain)

    use const_mod, only : pi

    type(stvec_swm_t), intent(inout) :: state
    type(domain_t),    intent(in)    :: domain

    integer(kind=4) :: i, j, k, t

    call state%u%assign(0.0_8,domain%mesh_u)
    call state%v%assign(0.0_8,domain%mesh_v)

    do t=domain%mesh_p%ts,domain%mesh_p%te
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js, domain%mesh_p%tile(t)%je
                do i = domain%mesh_p%tile(t)%is, domain%mesh_p%tile(t)%ie
                    if(abs(domain%mesh_p%tile(t)%rz(i,j,k))<sin(0.25_8*pi)) then
                        state%h%tile(t)%p(i,j,k) = 4000._8
                    else
                        state%h%tile(t)%p(i,j,k) = 1000._8
                    end if
                    !state%h%tile(t)%p(i,j,k) = 1000._8*(1._8+cos(asin(domain%mesh_p%tile(t)%rz(i,j,k)))**4)
                end do
            end do
        end do
    end do

end subroutine init_equilibrium

end module forcing_test_mod
