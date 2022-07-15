module config_swm_mod

use config_mod,        only : config_t
use config_domain_mod, only : config_domain_t

implicit none

type, public, extends(config_t) :: config_swm_t

    type(config_domain_t) :: config_domain

    character(:), allocatable :: swm_op_type
    character(:), allocatable :: v_components_type
    character(:), allocatable :: div_op_name
    character(:), allocatable :: grad_op_name
    character(:), allocatable :: curl_op_name
    character(:), allocatable :: coriolis_op_name
    character(:), allocatable :: KE_op_name
    character(:), allocatable :: massflux_op_name
    character(:), allocatable :: co2contra_op_name
    character(:), allocatable :: vector_adv_op_name
    character(:), allocatable :: quadrature_name
    character(:), allocatable :: diff_time_scheme
    character(:), allocatable :: hordiff_uv_name
    character(:), allocatable :: hordiff_h_name
    real(kind=8)              :: uv_diff_coeff = 0.0_8
    real(kind=8)              :: h_diff_coeff = 0.0_8
    real(kind=8)              :: dt              = 180.0_8
    real(kind=8)              :: tau_write       = 180.0_8
    real(kind=8)              :: tau_diagnostics = 180.0_8
    real(kind=8)              :: simulation_time = 86400.0_8

contains
    procedure, public :: parse
end type config_swm_t

contains

subroutine parse(this, config_string)

    use parcomm_mod,       only : parcomm_global

    class(config_swm_t),  intent(inout) :: this
    character(len=*),     intent(in)    :: config_string

    character(len=255) :: swm_op_type
    character(len=255) :: v_components_type
    character(len=255) :: div_op_name
    character(len=255) :: grad_op_name
    character(len=255) :: curl_op_name
    character(len=255) :: coriolis_op_name
    character(len=255) :: KE_op_name
    character(len=255) :: massflux_op_name
    character(len=255) :: co2contra_op_name
    character(len=255) :: vector_adv_op_name
    character(len=255) :: quadrature_name
    character(len=255) :: diff_time_scheme
    character(len=255) :: hordiff_uv_name
    character(len=255) :: hordiff_h_name
    real(kind=8)       :: uv_diff_coeff = 0.0_8
    real(kind=8)       :: h_diff_coeff = 0.0_8
    real(kind=8)       :: dt = 180.0
    real(kind=8)       :: tau_write = 180.0
    real(kind=8)       :: tau_diagnostics = 180.0
    real(kind=8)       :: simulation_time_days  = 0.0_8
    real(kind=8)       :: simulation_time_hours = 0.0_8
    real(kind=8)       :: simulation_time_min   = 0.0_8
    real(kind=8)       :: simulation_time_sec   = 3600.0_8

    namelist /shallow_water_model/ swm_op_type, v_components_type, &
                                   div_op_name, grad_op_name, curl_op_name, &
                                   coriolis_op_name, KE_op_name, massflux_op_name, &
                                   co2contra_op_name, vector_adv_op_name, &
                                   quadrature_name, hordiff_uv_name, &
                                   diff_time_scheme, &
                                   uv_diff_coeff, hordiff_h_name, h_diff_coeff, &
                                   dt, tau_write, tau_diagnostics, &
                                   simulation_time_sec, simulation_time_min, &
                                   simulation_time_hours, simulation_time_days

    read(config_string, shallow_water_model)

    call this%config_domain%parse(config_string)

    this%swm_op_type        = trim(swm_op_type)
    this%v_components_type  = trim(v_components_type)
    this%div_op_name        = trim(div_op_name)
    this%grad_op_name       = trim(grad_op_name)
    this%curl_op_name       = trim(curl_op_name)
    this%coriolis_op_name   = trim(coriolis_op_name)
    this%KE_op_name         = trim(KE_op_name)
    this%massflux_op_name   = trim(massflux_op_name)
    this%co2contra_op_name  = trim(co2contra_op_name)
    this%vector_adv_op_name = trim(vector_adv_op_name)
    this%quadrature_name    = trim(quadrature_name)
    this%diff_time_scheme   = trim(diff_time_scheme)
    this%hordiff_uv_name    = trim(hordiff_uv_name)
    this%hordiff_h_name     = trim(hordiff_h_name)
    this%uv_diff_coeff      = uv_diff_coeff
    this%h_diff_coeff       = h_diff_coeff
    this%dt = dt
    this%tau_write = tau_write
    this%tau_diagnostics = tau_diagnostics
    this%simulation_time = 86400._8*simulation_time_days+3600._8*simulation_time_hours+&
                                 60._8*simulation_time_min+simulation_time_sec

end subroutine parse

end module config_swm_mod
