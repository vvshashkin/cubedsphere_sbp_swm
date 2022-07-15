module config_nh_model_mod

use config_mod,        only : config_t
use config_domain_mod, only : config_domain_t
use config_postnh_mod, only : config_postnh_t

implicit none

type, extends(config_t) :: config_nh_model_t

    type(config_domain_t) :: config_domain
    type(config_postnh_t) :: config_postprocessing
    class(config_t), allocatable :: config_operator

    character(:), allocatable :: timescheme_name
    character(:), allocatable :: testcase_name
    character(:), allocatable :: operator_type

    real(kind=8) :: dt
    real(kind=8) :: tau_write
    real(kind=8) :: tau_diagnostics
    real(kind=8) :: simulation_time

    contains
    procedure, public :: parse

end type config_nh_model_t

contains

subroutine parse(this, config_string)

    use const_mod,  only : Earth_radii
    use config_nh_operator_mod, only : get_nh_operator_config

    class(config_nh_model_t),  intent(inout) :: this
    character(len=*),          intent(in)    :: config_string

    character(len=255) :: timescheme_name
    character(len=255) :: testcase_name
    character(len=255) :: operator_type

    real(kind=8) :: dt, tau_write, tau_diagnostics
    real(kind=8) :: simulation_time_days=0._8, simulation_time_hours=0._8, &
                    simulation_time_min=0._8, simulation_time_sec=0._8

    namelist /nh_model/ timescheme_name, testcase_name, operator_type,   &
                        dt, tau_write, tau_diagnostics,                  &
                        simulation_time_sec, simulation_time_min,        &
                        simulation_time_hours, simulation_time_days

    read(config_string, nh_model)

    this%timescheme_name   = trim(timescheme_name)
    this%testcase_name     = trim(testcase_name)
    this%operator_type     = trim(operator_type)

    this%dt = dt
    this%tau_write = tau_write
    this%tau_diagnostics = tau_diagnostics
    this%simulation_time = 86400._8*simulation_time_days+3600._8*simulation_time_hours+&
                                 60._8*simulation_time_min+simulation_time_sec

    call this%config_domain%parse(config_string)
    this%config_domain%config_metric%vertical_scale = this%config_domain%h_top
    !this%config_domain%config_metric%scale = Earth_radii

    call this%config_postprocessing%parse(config_string)

    this%config_operator = get_nh_operator_config(this%operator_type)
    call this%config_operator%parse(config_string)
end subroutine parse


end module config_nh_model_mod
