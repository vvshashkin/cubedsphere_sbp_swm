module abstract_vertical_profile_mod

implicit none

type, abstract :: vertical_profile_t
    contains
    procedure, public :: calc_pressure   => return_NaN
    procedure, public :: calc_temp       => return_NaN
    procedure, public :: calc_theta      => return_NaN
    procedure, public :: calc_dtheta_dz  => return_NaN
    procedure, public :: calc_ExnerP     => return_NaN
    procedure, public :: calc_dExnerP_dz => return_NaN
end type vertical_profile_t

contains

pure elemental function return_NaN(this, h, t_surf, p_surf) result(p)
    use, intrinsic :: ieee_arithmetic, only : ieee_value, IEEE_SIGNALING_NAN

    class(vertical_profile_t), intent(in) :: this
    real(kind=8), intent(in) :: h, t_surf, p_surf

    real(kind=8) :: p

    p = ieee_value(p,IEEE_SIGNALING_NAN)

end function return_NaN

end module abstract_vertical_profile_mod
