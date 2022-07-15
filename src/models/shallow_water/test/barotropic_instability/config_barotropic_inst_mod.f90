module config_barotropic_inst_mod

use config_mod, only : config_t
use const_mod,  only : Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

implicit none

type, public, extends(config_t) :: config_barotropic_inst_t

    integer(kind=4) :: N, Nz

    real(kind=8)    :: grav   = Earth_grav
    real(kind=8)    :: omega  = Earth_omega
    real(kind=8)    :: H0     = 11000._8
    real(kind=8)    :: a      = Earth_radii
    real(kind=8)    :: u0     = 80.0_8
    real(kind=8)    :: h_pert = 120.0_8
    integer(kind=4) :: Nq     = 120

contains
    procedure, public :: parse
end type config_barotropic_inst_t

contains

subroutine parse(this, config_string)

    use parcomm_mod, only : parcomm_global

    class(config_barotropic_inst_t), intent(inout) :: this
    character(len=*),                intent(in)    :: config_string

    integer(kind=4) :: N, Nz

    real(kind=8)    :: grav   = Earth_grav
    real(kind=8)    :: omega  = Earth_omega
    real(kind=8)    :: H0     = 29400/Earth_grav
    real(kind=8)    :: a      = Earth_radii
    real(kind=8)    :: u0     = pi*Earth_radii/6.0_8/Earth_sidereal_T
    real(kind=8)    :: h_pert = 120.0_8
    integer(kind=4) :: Nq     = 120

    namelist /barotropic_inst/ grav, omega, H0, a, u0, Nq

    read(config_string, barotropic_inst)

    this%grav   = grav
    this%omega  = omega
    this%H0     = H0
    this%a      = a
    this%u0     = u0
    this%h_pert = h_pert
    this%Nq     = Nq

end subroutine parse

end module config_barotropic_inst_mod
