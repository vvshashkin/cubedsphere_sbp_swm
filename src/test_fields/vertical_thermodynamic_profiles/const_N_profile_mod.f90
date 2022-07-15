module const_N_profile_mod

use abstract_vertical_profile_mod, only : vertical_profile_t
use const_mod,                     only : Earth_grav, Cp, kappa

implicit none

type, extends(vertical_profile_t) :: const_N_profile_t

    real(kind=8) :: N
    real(kind=8) :: grav = Earth_grav
    real(kind=8) :: Cp = Cp, kappa = kappa

    contains
    procedure, public :: calc_pressure
    procedure, public :: calc_temp
    procedure, public :: calc_theta
    procedure, public :: calc_dtheta_dz
    procedure, public :: calc_ExnerP
    procedure, public :: calc_dExnerP_dz
end type const_N_profile_t

contains

pure elemental function calc_temp(this, h, t_surf, p_surf) result(t)
    class(const_N_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: t

    t = this%calc_theta(h,t_surf,p_surf) / this%calc_ExnerP(h,t_surf,p_surf)

end function calc_temp

pure elemental function calc_pressure(this, h, t_surf, p_surf) result(p)
    class(const_N_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: p

    p = this%calc_ExnerP(h,t_surf,p_surf) ** (1.0_8/this%kappa)*p_surf

end function calc_pressure

pure elemental function calc_theta(this, h, t_surf, p_surf) result(t)
    class(const_N_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: t

    t = t_surf*exp(this%N**2*h / this%grav)

end function calc_theta

pure elemental function calc_dtheta_dz(this, h, t_surf, p_surf) result(t)
    class(const_N_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: t

    t = this%N**2/this%grav*t_surf*exp(this%N**2*h / this%grav)

end function calc_dtheta_dz

pure elemental function calc_ExnerP(this, h, t_surf, p_surf) result(p)
    class(const_N_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: p

    p = 1.0_8+this%grav**2/(this%N**2*this%Cp*t_surf)*(exp(-this%N**2/this%grav*h)-1.0_8)

end function calc_ExnerP

pure elemental function calc_dExnerP_dz(this, h, t_surf, p_surf) result(p)
    class(const_N_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: p

    p = -this%grav/(this%Cp*t_surf)*exp(-this%N**2/this%grav*h)

end function calc_dExnerP_dz

end module const_N_profile_mod
