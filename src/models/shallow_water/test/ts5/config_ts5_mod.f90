module config_ts5_mod

use config_mod, only : config_t
use const_mod,  only : Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

implicit none

type, public, extends(config_t) :: config_ts5_t

    integer(kind=4) :: N, Nz

    real(kind=8) :: grav    = Earth_grav
    real(kind=8) :: omega   = Earth_omega
    real(kind=8) :: h_mean  = 5960.0_8
    real(kind=8) :: a       = Earth_radii
    real(kind=8) :: u0      = 20.0_8
    real(kind=8) :: r_mount = pi / 9.0_8
    real(kind=8) :: h_mount = 2000.0_8
    real(kind=8) :: x_mount = 0.0_8
    real(kind=8) :: y_mount = 0.5_8*sqrt(3.0_8)
    real(kind=8) :: z_mount = 0.5_8

contains
    procedure, public :: parse
end type config_ts5_t

contains

subroutine parse(this, config_string)

    use parcomm_mod, only : parcomm_global

    class(config_ts5_t), intent(inout) :: this
    character(len=*),    intent(in)    :: config_string

    integer(kind=4) :: N, Nz

    real(kind=8) :: grav   = Earth_grav
    real(kind=8) :: omega  = Earth_omega
    real(kind=8) :: h_mean = 29400/Earth_grav
    real(kind=8) :: a      = Earth_radii
    real(kind=8) :: u0      = 20.0_8
    real(kind=8) :: r_mount = pi / 9.0_8
    real(kind=8) :: h_mount = 2000.0_8
    real(kind=8) :: x_mount = 0.0_8
    real(kind=8) :: y_mount = 0.5_8*sqrt(3.0_8)
    real(kind=8) :: z_mount = 0.5_8

    namelist /ts5/ grav, omega, h_mean, a, u0, r_mount, h_mount, &
                   x_mount, y_mount, z_mount

    read(config_string, ts5)

    this%grav    = grav
    this%omega   = omega
    this%h_mean  = h_mean
    this%a       = a
    this%u0      = u0
    this%r_mount = r_mount
    this%h_mount = h_mount
    this%x_mount = x_mount
    this%y_mount = y_mount
    this%z_mount = z_mount

end subroutine parse

end module config_ts5_mod
