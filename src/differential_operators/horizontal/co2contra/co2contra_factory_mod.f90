module co2contra_factory_mod

use domain_mod,              only : domain_t
use abstract_co2contra_mod,  only : co2contra_operator_t
use co2contra_colocated_mod, only : co2contra_colocated_t
use co2contra_Cgrid_mod,     only : co2contra_c_sbp_t, co2contra_c_sbp_new_t!, co2contra_c_halo_t
use co2contra_ch_mod,        only : co2contra_ch_sbp_t
use parcomm_mod,             only : parcomm_global
use exchange_factory_mod,    only : create_symmetric_halo_vec_exchange_C
use sbp_factory_mod,         only : create_sbp_operator

implicit none

contains

function create_co2contra_operator(domain, co2contra_operator_name) result(co2contra)
    type(domain_t),    intent(in)    :: domain
    character(len=*),  intent(in)    :: co2contra_operator_name

    class(co2contra_operator_t), allocatable :: co2contra

    select case(co2contra_operator_name)
    case("co2contra_colocated")
        co2contra = co2contra_colocated_t()
    case("co2contra_c_sbp21", "co2contra_c_sbp42")
        co2contra = create_co2contra_c_sbp_operator(domain, co2contra_operator_name)
    case("co2contra_c_sbp21_new", "co2contra_c_sbp42_new")
        co2contra = create_co2contra_c_sbp_new_operator(domain, co2contra_operator_name)
    case("co2contra_ch_sbp21")
        co2contra = create_co2contra_ch_sbp_operator(domain, co2contra_operator_name)
    case("co2contra_ch_sbp42")
        co2contra = create_co2contra_ch_sbp_operator(domain, co2contra_operator_name)
    case default
        call parcomm_global%abort("unknown co2contra operator: "//co2contra_operator_name)
    end select

end function create_co2contra_operator

function create_co2contra_c_sbp_operator(domain, co2contra_operator_name) result(co2contra)

    type(domain_t),          intent(in) :: domain
    character(len=*),        intent(in) :: co2contra_operator_name
    type(co2contra_c_sbp_t)             :: co2contra

    integer(kind=4) :: halo_width

    co2contra%operator_name = co2contra_operator_name

    select case(co2contra_operator_name)
    case("co2contra_c_sbp21")
        halo_width = 1
        co2contra%sbp_interp_h2v = create_sbp_operator("W21_stagered_interp_c2i")
        co2contra%sbp_interp_v2h = create_sbp_operator("W21_stagered_interp_i2c")
    case("co2contra_c_sbp42")
        halo_width = 3
        co2contra%sbp_interp_h2v = create_sbp_operator("W42_stagered_interp_c2i")
        co2contra%sbp_interp_v2h = create_sbp_operator("W42_stagered_interp_i2c")
    case default
        call parcomm_global%abort("unknown co2contra_c_sbp operator "// co2contra_operator_name)
    end select

    co2contra%exchange_inner = &
          create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                                     domain%topology, halo_width, 'full')
end function

function create_co2contra_c_sbp_new_operator(domain, co2contra_operator_name) result(co2contra)

    use interpolator2d_factory_mod,    only : create_vec2vec_interpolator2d
    use grid_field_factory_mod,        only : create_grid_field

    type(domain_t),          intent(in) :: domain
    character(len=*),        intent(in) :: co2contra_operator_name
    type(co2contra_c_sbp_new_t)         :: co2contra

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 4

    co2contra%operator_name = co2contra_operator_name

    select case(co2contra_operator_name)
    case("co2contra_c_sbp21_new")
        call create_vec2vec_interpolator2d(co2contra%interp_h2v_op, "interp2d_pvec2uv_C_sbp21", domain)
        call create_vec2vec_interpolator2d(co2contra%interp_v2h_op, "interp2d_uv2pvec_C_sbp21", domain)
    case("co2contra_c_sbp42_new")
        call create_vec2vec_interpolator2d(co2contra%interp_h2v_op, "interp2d_pvec2uv_C_sbp42", domain)
        call create_vec2vec_interpolator2d(co2contra%interp_v2h_op, "interp2d_uv2pvec_C_sbp42", domain)
    case default
        call parcomm_global%abort("unknown co2contra_c_sbp operator "// co2contra_operator_name)
    end select

    call create_grid_field(co2contra%uh, halo_width, 0, domain%mesh_u)
    call create_grid_field(co2contra%vh, halo_width, 0, domain%mesh_v)

end function create_co2contra_c_sbp_new_operator

function create_co2contra_ch_sbp_operator(domain, co2contra_operator_name) result(co2contra)

    use interpolator2d_factory_mod,   only : create_vec2vec_interpolator2d
    use grid_field_factory_mod,       only : create_grid_field

    type(domain_t),             intent(in) :: domain
    character(len=*),           intent(in) :: co2contra_operator_name
    type(co2contra_ch_sbp_t)               :: co2contra

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 4

    co2contra%operator_name = co2contra_operator_name

    select case(co2contra_operator_name)
    case("co2contra_ch_sbp21")
        call create_vec2vec_interpolator2d(co2contra%interp_q2uv_op, &
                                           "interp2d_qvec2uv_Ch_sbp21", domain)
        call create_vec2vec_interpolator2d(co2contra%interp_uv2q_op, &
                                           "interp2d_uv2qvec_Ch_sbp21", domain)
    case("co2contra_ch_sbp42")
        call create_vec2vec_interpolator2d(co2contra%interp_q2uv_op, &
                                           "interp2d_qvec2uv_Ch_sbp42", domain)
        call create_vec2vec_interpolator2d(co2contra%interp_uv2q_op, &
                                           "interp2d_uv2qvec_Ch_sbp42", domain)
    case default
        call parcomm_global%abort("unknown co2contra_c_sbp operator "// co2contra_operator_name)
    end select

    call create_grid_field(co2contra%uq, halo_width, 0, domain%mesh_o)
    call create_grid_field(co2contra%vq, halo_width, 0, domain%mesh_o)

end function create_co2contra_ch_sbp_operator

! function create_co2contra_c_halo_operator(domain) result(co2contra)
!
!     use halo_factory_mod,       only : create_vector_halo_procedure
!
!     type(domain_t),          intent(in) :: domain
!     type(co2contra_c_halo_t)            :: co2contra
!
!     integer(kind=4) :: halo_width
!
!     call create_vector_halo_procedure(co2contra%halo,domain,2,"ecs_C_vec")
! end function

end module co2contra_factory_mod
