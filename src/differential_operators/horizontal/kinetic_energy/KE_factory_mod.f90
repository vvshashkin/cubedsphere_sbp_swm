module KE_factory_mod

use abstract_KE_mod, only : KE_operator_t
use domain_mod,      only : domain_t
use parcomm_mod,     only : parcomm_global

implicit none

contains

subroutine create_KE_operator(ke_operator, ke_operator_name, domain)

    use ke_colocated_mod, only : ke_colocated_t

    class(KE_operator_t), allocatable, intent(out) :: ke_operator
    character(len=*),                  intent(in)  :: ke_operator_name
    type(domain_t),                    intent(in)  :: domain

    select case(ke_operator_name)

    case("KE_colocated")
        ke_operator = ke_colocated_t()
    case("KE_Cgrid_sbp42")
        call create_KE_Cgrid_sbp(ke_operator, "interp2d_uv2pvec_C_sbp42", domain)
    case("KE_Cgrid_sbp21")
        call create_KE_Cgrid_sbp(ke_operator, "interp2d_uv2pvec_C_sbp21", domain)
    case default
        call parcomm_global%abort("Unknown KE operator: "//ke_operator_name)
    end select

end subroutine create_KE_operator

subroutine create_KE_Cgrid_sbp(ke_operator, sbp_i2c_interp_name, domain)

    use grid_field_factory_mod,       only : create_grid_field
    use ke_Cgrid_mod,                 only : ke_Cgrid_t
    use interpolator2d_factory_mod,   only : create_vec2vec_interpolator2d

    class(KE_operator_t), allocatable, intent(out) :: ke_operator
    character(len=*),                  intent(in)  :: sbp_i2c_interp_name
    type(domain_t),                    intent(in)  :: domain

    type(ke_Cgrid_t), allocatable :: ke_cgrid_op
    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 4

    allocate(ke_cgrid_op)

    call create_grid_field(ke_cgrid_op%KE_u, halo_width, 0, domain%mesh_u)
    call create_grid_field(ke_cgrid_op%KE_v, halo_width, 0, domain%mesh_v)

    call create_grid_field(ke_cgrid_op%KE_uh, 0, 0, domain%mesh_u)
    call create_grid_field(ke_cgrid_op%KE_vh, 0, 0, domain%mesh_v)

    !WORKAROUND
    call create_vec2vec_interpolator2d(ke_cgrid_op%interp_op, sbp_i2c_interp_name , domain)

    call move_alloc(ke_cgrid_op, ke_operator)

end subroutine create_KE_Cgrid_sbp

end module KE_factory_mod
