module config_ts2_mod

use config_mod, only : config_t
use const_mod,  only : Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

implicit none

type, public, extends(config_t) :: config_ts2_t

    integer(kind=4) :: N, Nz

    real(kind=8) :: grav   = Earth_grav
    real(kind=8) :: omega  = Earth_omega
    real(kind=8) :: h_mean = 29400/Earth_grav
    real(kind=8) :: a      = Earth_radii
    real(kind=8) :: u0     = pi*Earth_radii/6.0_8/Earth_sidereal_T

contains
    procedure, public :: parse
end type config_ts2_t

contains

subroutine parse(this, config_string)

    use parcomm_mod, only : parcomm_global

    class(config_ts2_t), intent(inout) :: this
    character(len=*),    intent(in)    :: config_string

    integer(kind=4) :: N, Nz

    real(kind=8) :: grav   = Earth_grav
    real(kind=8) :: omega  = Earth_omega
    real(kind=8) :: h_mean = 29400/Earth_grav
    real(kind=8) :: a      = Earth_radii
    real(kind=8) :: u0     = pi*Earth_radii/6.0_8/Earth_sidereal_T

    namelist /ts2/ grav, omega, h_mean, a, u0

    read(config_string, ts2)

    this%grav   = grav
    this%omega  = omega
    this%h_mean = h_mean
    this%a      = a
    this%u0     = u0

end subroutine parse

end module config_ts2_mod
