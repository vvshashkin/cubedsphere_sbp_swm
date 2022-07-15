module isothermal_profile_mod

use abstract_vertical_profile_mod, only : vertical_profile_t
use const_mod,                     only : Earth_grav, rgaz, kappa

implicit none

type, extends(vertical_profile_t) :: isothermal_profile_t

    real(kind=8) :: grav = Earth_grav
    real(kind=8) :: rgaz = rgaz, kappa = kappa
    real(kind=8) :: p0 = 1e5_8 !for potential temperature and Exner function definitions

    contains
    procedure, public :: calc_pressure
    procedure, public :: calc_temp
    procedure, public :: calc_theta
    procedure, public :: calc_dtheta_dz
    procedure, public :: calc_ExnerP
    procedure, public :: calc_dExnerP_dz
end type isothermal_profile_t

contains

pure elemental function calc_temp(this, h, t_surf, p_surf) result(t)
    class(isothermal_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: t

    t = t_surf

end function calc_temp

pure elemental function calc_pressure(this, h, t_surf, p_surf) result(p)
    class(isothermal_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: p

    p = p_surf*exp(-h*this%grav / (this%rgaz*t_surf))

end function calc_pressure

pure elemental function calc_theta(this, h, t_surf, p_surf) result(t)
    class(isothermal_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: t

    t = t_surf*(this%p0/p_surf*exp(h*this%grav / (t_surf*this%rgaz)))**this%kappa

end function calc_theta

pure elemental function calc_dtheta_dz(this, h, t_surf, p_surf) result(t)
    class(isothermal_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: t

    t = t_surf*(this%p0/p_surf)**this%kappa*exp(this%kappa*h*this%grav / (t_surf*this%rgaz))*&
                this%kappa*this%grav / (t_surf*this%rgaz)

end function calc_dtheta_dz

pure elemental function calc_ExnerP(this, h, t_surf, p_surf) result(p)
    class(isothermal_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: p

    p = (p_surf/this%p0*exp(-h*this%grav / (this%rgaz*t_surf)))**this%kappa

end function calc_ExnerP

pure elemental function calc_dExnerP_dz(this, h, t_surf, p_surf) result(p)
    class(isothermal_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: p

    p = -(p_surf/this%p0)**this%kappa*exp(-this%kappa*h*this%grav / (this%rgaz*t_surf))*&
                                           this%kappa*this%grav / (this%rgaz*t_surf)

end function calc_dExnerP_dz

end module isothermal_profile_mod
