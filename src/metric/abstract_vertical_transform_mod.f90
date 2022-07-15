module abstract_vertical_transform_mod

    implicit none

    type, abstract :: vertical_transform_t
    contains
        procedure(calc_z_i),     public, deferred :: calc_z
        procedure(calc_z_i),     public, deferred :: calc_dz_deta
        procedure(calc_dz_dh_i), public, deferred :: calc_dz_dh_surf
        procedure(calc_dz_dh_i), public, deferred :: calc_dz_dh_top
    end type

    abstract interface
        pure function calc_z_i(this, h_surf, h_top, eta) result(z)
            import vertical_transform_t
            class(vertical_transform_t), intent(in) :: this
            real(kind=8),                intent(in) :: h_surf, h_top, eta
            real(kind=8)                            :: z
        end function calc_z_i

        pure function calc_dz_dh_i(this, eta) result(dz_dh)
            import vertical_transform_t
            class(vertical_transform_t), intent(in) :: this
            real(kind=8),                intent(in) :: eta
            real(kind=8)                            :: dz_dh
        end function calc_dz_dh_i
    end interface

contains

end module abstract_vertical_transform_mod
