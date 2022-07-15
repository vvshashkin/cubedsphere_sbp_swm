module vertical_transform_mod

    use abstract_vertical_transform_mod, only : vertical_transform_t

    implicit none

    !basic case: eta = (z-h_surf)/(h_top-h_surf)
    !            z = h_surf+eta*(h_top-h_surf)
    type, public, extends(vertical_transform_t) :: vertical_transform_default_t
    contains
        procedure :: calc_z
        procedure :: calc_dz_deta
        procedure :: calc_dz_dh_surf
        procedure :: calc_dz_dh_top
    end type

contains

    pure function calc_z(this, h_surf, h_top, eta) result(z)
        class(vertical_transform_default_t), intent(in) :: this
        real(kind=8), intent(in) :: h_surf, h_top, eta
        real(kind=8) :: z

        z = h_surf+(h_top-h_surf)*eta
    end function calc_z

    pure function calc_dz_deta(this, h_surf, h_top, eta) result(dz_deta)
        class(vertical_transform_default_t), intent(in) :: this
        real(kind=8), intent(in) :: h_surf, h_top, eta
        real(kind=8) :: dz_deta

        dz_deta = (h_top-h_surf)
    end function calc_dz_deta

    pure function calc_dz_dh_surf(this,eta) result(dz_dh)
        class(vertical_transform_default_t), intent(in) :: this
        real(kind=8), intent(in) :: eta
        real(kind=8) :: dz_dh

        dz_dh = 1.0_8-eta
    end function calc_dz_dh_surf

    pure function calc_dz_dh_top(this, eta) result(dz_dh)
        class(vertical_transform_default_t), intent(in) :: this
        real(kind=8), intent(in) :: eta
        real(kind=8) :: dz_dh

        dz_dh = eta
    end function calc_dz_dh_top

end module vertical_transform_mod
